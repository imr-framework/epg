import math
from collections import OrderedDict

import gpi
import h5py
import numpy as np
import copy

from pulseq.core import convert
from pulseq.core.Sequence.sequence import Sequence
from pulseq.core.calc_duration import calc_duration
from pulseq.core.make_adc import makeadc
from pulseq.core.makearbitrary_grad import makearbitrary_grad
from pulseq.core.make_block import make_block_pulse
from pulseq.core.make_delay import make_delay
from pulseq.core.make_sinc import make_sinc_pulse
from pulseq.core.make_trap import make_trapezoid


class ExternalNode(gpi.NodeAPI):
    def initUI(self):
        # Widgets
        self.addWidget('TextBox', 'Events you defined')
        self.addWidget('StringBox', 'Order of events')
        self.addWidget('PushButton', 'ComputeEvents', button_title="Compute events")
        self.addWidget('TextBox', 'Sequence details')

        # IO Ports
        self.addInPort(title='input', type='DICT')
        self.addInPort(title='epi_readout', type='DICT', obligation=gpi.OPTIONAL)
        self.addOutPort(title='seq_output', type='DICT')
        self.addOutPort(title='adc_output', type='NPYarray')
        self.addOutPort(title='rf_mag_output', type='NPYarray')
        self.addOutPort(title='rf_phase_output', type='NPYarray')
        self.addOutPort(title='trap_x_output', type='NPYarray')
        self.addOutPort(title='trap_y_output', type='NPYarray')
        self.addOutPort(title='trap_z_output', type='NPYarray')

        self.all_event_defs = OrderedDict()
        self.all_events_ordered = OrderedDict()

        return 0

    def validate(self):
        self.in_dict = self.getData('input')
        if 'sequence' in self.in_dict:
            self.seq = self.in_dict['sequence']
        else:
            self.all_event_defs = self.in_dict['all_event_defs']
            self.all_events_ordered = self.in_dict['all_events_ordered']

            events_added_text = str()
            for unique_name in self.all_events_ordered:
                events_added_text += unique_name + '\n'
            self.setAttr('Events you defined', **{'val': events_added_text})

    def compute(self):
        if 'ComputeEvents' in self.widgetEvents():
            user_ordered_events = self.getVal('Order of events').split(',')
            if len(user_ordered_events) == 0:
                raise ValueError('Enter [Unique Name] of Events in the order you want Events to be added.')

            self.make_event_holders()
            self.add_blocks_to_seq(user_ordered_events)
            self.set_plot_outputs()

            # To display details of Sequence object in TextBox
            self.setAttr('Sequence details', val=str(self.seq))

            # Setting output
            self.in_dict['seq'] = self.seq
            self.setData('seq_output', self.in_dict)

            return 0

    def make_event_holders(self):
        """Make appropriate Holder objects depending on the Event type."""

        self.system = self.in_dict['system']
        # arbgrad_file_path is only for arbitrary gradients
        arbgrad_file_path = self.all_event_defs['file_path'] if 'file_path' in self.all_event_defs else None
        self.all_event_holders = {}

        for event in self.all_event_defs:
            event_unique_name = event['event_unique_name']
            event_name = event['event_name']
            event_values = list(event['event_values'].values())
            include_in_loop = event['include_in_loop']

            if event_name == 'Delay':
                params = self.parse_config_params(event_values)
                delay = make_delay(params[0])
                self.all_event_holders[event_unique_name] = delay, include_in_loop
            elif event_name == 'SincRF':
                include_gz = event['include_gz']
                max_grad, max_slew, flip_angle, duration, freq_offset, phase_offset, time_bw_product, apodization, slice_thickness = self.parse_config_params(
                    event_values)
                flip_angle = math.radians(flip_angle)
                max_grad = convert.convert_from_to(max_grad, 'mT/m')
                max_slew = convert.convert_from_to(max_slew, 'mT/m/ms')
                max_grad = self.system.max_grad if max_grad == 0 else max_grad
                max_slew = self.system.max_slew if max_slew == 0 else max_slew
                kwargs_for_sinc = {"flip_angle": flip_angle, "system": self.system, "duration": duration,
                                   "freq_offset": freq_offset, "phase_offset": phase_offset,
                                   "time_bw_product": time_bw_product, "apodization": apodization,
                                   "max_grad": max_grad, "max_slew": max_slew, "slice_thickness": slice_thickness}
                if include_gz:
                    rf, gz = make_sinc_pulse(kwargs_for_sinc, 2)
                    self.all_event_holders[event_unique_name] = rf, include_in_loop
                    self.all_event_holders['gz_' + event_unique_name] = gz, include_in_loop
                else:
                    rf = make_sinc_pulse(kwargs_for_sinc)
                    self.all_event_holders[event_unique_name] = rf, include_in_loop
            elif event_name == 'BlockRF':
                include_gz = event['include_gz']
                max_grad, max_slew, flip_angle, duration, freq_offset, phase_offset, time_bw_product, bandwidth, slice_thickness = self.parse_config_params(
                    event_values)
                flip_angle = math.radians(flip_angle)
                max_grad = convert.convert_from_to(max_grad, 'mT/m')
                max_slew = convert.convert_from_to(max_slew, 'mT/m/ms')
                max_grad = self.system.max_grad if max_grad == 0 else max_grad
                max_slew = self.system.max_slew if max_slew == 0 else max_slew
                kwargs_for_block = {"flip_angle": flip_angle, "system": self.system, "duration": duration,
                                    "freq_offset": freq_offset, "phase_offset": phase_offset,
                                    "time_bw_product": time_bw_product, "bandwidth": bandwidth,
                                    "max_grad": max_grad, "max_slew": max_slew, "slice_thickness": slice_thickness}
                if include_gz:
                    rf, gz = make_block_pulse(kwargs_for_block, 2)
                    self.all_event_holders[event_unique_name] = rf, include_in_loop
                    self.all_event_holders['gz_' + event_unique_name] = gz, include_in_loop
                else:
                    rf = make_block_pulse(kwargs_for_block)
                    self.all_event_holders[event_unique_name] = rf, include_in_loop
            elif event_name == 'G':
                channel = event_values.pop(0)
                max_grad, max_slew, duration, area, flat_time, flat_area, amplitude, rise_time = self.parse_config_params(
                    event_values)

                # area, flat_area and amplitude should be reset to -1 if user does not input any values. This is
                # because the default values are -1 in maketrap method.
                area = area if area != 0 else -1
                flat_area = flat_area if flat_area != 0 else -1
                amplitude = amplitude if amplitude != 0 else -1

                max_grad = convert.convert_from_to(max_grad, 'mT/m')
                max_slew = convert.convert_from_to(max_slew, 'mT/m/ms')
                max_grad = self.system.max_grad if max_grad == 0 else max_grad
                max_slew = self.system.max_slew if max_slew == 0 else max_slew
                kwargs_for_trap = {"channel": channel, "system": self.system, "duration": duration, "area": area,
                                   "flat_time": flat_time, "flat_area": flat_area, "amplitude": amplitude,
                                   "max_grad": max_grad, "max_slew": max_slew, "rise_time": rise_time}
                trap = make_trapezoid(kwargs_for_trap)
                self.all_event_holders[event_unique_name] = trap, include_in_loop
            elif event_name == 'GyPre':
                duration, area = self.parse_config_params(event_values)
                Ny = self.system.Ny
                delta_k = 1 / self.system.fov
                gy_pre_list = []

                for i in range(int(Ny)):
                    kwargs_for_gy_pre = {"channel": 'y', "system": self.system, "area": (i - Ny / 2) * delta_k,
                                         "duration": duration}
                    if area != 0:
                        kwargs_for_gy_pre['area'] = area
                    gy_pre = make_trapezoid(kwargs_for_gy_pre)
                    gy_pre_list.append(gy_pre)

                self.all_event_holders[event_unique_name] = gy_pre_list, True
            elif event_name == 'ArbGrad':
                channel = event_values.pop(0)
                max_grad, max_slew = self.parse_config_params(event_values)
                file = h5py.File(gpi.TranslateFileURI(arbgrad_file_path), "r")
                self.dataset = str()

                def append_if_dataset(name, obj):
                    if isinstance(obj, h5py.Dataset):
                        self.dataset = name
                        return True

                file.visititems(append_if_dataset)

                waveform = file[self.dataset].value
                kwargs_for_arb_grad = {"channel": channel, "waveform": waveform, "max_grad": max_grad,
                                       "max_slew": max_slew, "system": self.system}
                arb_grad = makearbitrary_grad(kwargs_for_arb_grad)
                self.all_event_holders[event_unique_name] = arb_grad, include_in_loop
            elif event_name == 'ADC':
                num_samples, dwell, duration, delay, freq_offset, phase_offset = self.parse_config_params(
                    event_values)
                kwargs_for_adc = {"num_samples": num_samples, "system": self.system, "dwell": dwell,
                                  "duration": duration, "delay": delay, "freq_offset": freq_offset,
                                  "phase_offset": phase_offset}
                adc = makeadc(kwargs_for_adc)
                self.all_event_holders[event_unique_name] = adc, include_in_loop

    def parse_config_params(self, all_params):
        """
        Parse simple single term analytical expressions.
        Syntax: [-]event_unique_name.event_property[operator][operand]
        [] - Optional
        """

        parsed_params = []
        for param in all_params:
            try:
                parsed_params.append(float(param))
            except ValueError:
                # Syntax: [-]event_unique_name.event_property[operator][operand]
                # [] - Optional
                if param.count('.') == 1:
                    sign = -1 if param[0] == '-' else 1
                    param = param[1:] if sign == -1 else param
                    p_split = param.split('.')
                    event_unique_name = p_split[0]
                    event_property = p_split[1]
                    operator = '/' if '/' in event_property else '*' if '*' in event_property else ''
                    operand = float(event_property.split(operator)[1]) if operator != '' else 1
                    event_property = event_property.split(operator)[0] if operator != '' else event_property
                    if event_unique_name in self.all_event_holders:
                        event_holder = self.all_event_holders[event_unique_name][0]
                        value = sign * float(getattr(event_holder, event_property))
                        value = value / operand if operator == '/' else value * operand
                        parsed_params.append(float(value))
                    else:
                        raise ValueError(
                            param + '\nValue not found. The syntax for referring to other event parameters_params is: ['
                                    'sign][*event_unique_name].[*event_property][operator][operand]. * - '
                                    'required')
                else:
                    raise ValueError(
                        param + '\nValue not found. The syntax for referring to other event parameters_params is: ['
                                'sign][*event_unique_name].[*event_property][operator][operand]. * - '
                                'required')

        if len(parsed_params) == 0:
            raise ValueError('Please make sure you have input correct configuration parameters.')
        return parsed_params

    def add_blocks_to_seq(self, user_ordered_events):
        """
        This method creates a Sequence object and adds Events.

        - 'epi_events' is not a Node; so if unique_node_name == 'epi_events' that means the single-shot EPI readout
        Events have to be added to the Sequence object - The single-shot EPI readout Events are defined in
        self.ss_epi_event_unique_names_list
        - At the end of Every EPI readout, the amplitude has to be inverted. This is done by bookmarking the index of
        the Event whose amplitude has to be inverted. This index is used at the end of every for loop iteration while
        adding the Events to the Sequence object
        """
        self.seq = Sequence(self.system)
        Ny = int(self.system.Ny)
        epi_readout = self.getData('epi_readout')
        epi_nodes = epi_readout['epi_events'].split(',') if epi_readout is not None else None
        epi_event_amp_inv = epi_readout['epi_event_amplitude_inversion'] if epi_readout is not None else None
        epi_events_to_add = []

        # For example, in ss_epi_demo.py:
        # for i in range(self.Ny):
        # - seq.add_block(gx, adc)
        # - seq.add_block(gy)
        # Two Nodes encapsulate the Events (gx, adc) and (gy). These Events *together* have to be added Ny times.
        # A typical implementation would add the Events of each Node Ny times, and then the next set of Events belonging
        # to the next Node. However, we want to add all these Events Ny times, in that order.
        # First, check that user has mentioned EPI Events
        if epi_readout is not None:
            for i in range(Ny):
                for each_epi_node in epi_nodes:
                    _epi_list = []

                    for epi_event_name in self.all_events_ordered[each_epi_node]:
                        epi_event = self.all_event_holders[epi_event_name][0]
                        print(id(epi_event))
                        # Invert amplitude
                        if i != 0 and epi_event_name == epi_event_amp_inv:
                            epi_event.amplitude = -epi_event.amplitude

                        # Create a deepcopy of _epi_event; because Python points multiple objects to same memory
                        _temp_epi_event = copy.deepcopy(epi_event)
                        _epi_list.append(_temp_epi_event)

                    epi_events_to_add.append(_epi_list)
            # If EPI readout Events are present, the other Events have to only be added once. So, Ny = 1
            Ny = 1

        # If a particular Event does not need to be added Ny times, but only once, add the unique_event_name to a
        # blacklist
        _events_blacklist = []
        for i in range(Ny):
            for unique_node_name in user_ordered_events:
                if unique_node_name == 'epi_events':
                    if len(epi_events_to_add) == 0:
                        raise ValueError('Please specify the Events which are EPI readout Events')
                    for _epi_list in epi_events_to_add:
                        self.seq.add_block(*_epi_list)
                else:
                    events_to_add = []
                    for event_unique_name in self.all_events_ordered[unique_node_name]:
                        if event_unique_name not in _events_blacklist:
                            event_def = self.get_eventdef_for_unique_name(event_unique_name)
                            epi_event = self.all_event_holders[event_unique_name][0]

                            events_to_add.append(epi_event[i] if event_unique_name == 'gy_pre' else epi_event)

                            # Initialize gz as None. Since the user can assign any unique name to the Rf Event, there is
                            # no definite way to check if the Gz Event has to be added. So, append the Rf Event's unique
                            # name to gz (as 'gz_'), and check if this exists in all_event_holders.
                            # (Look at make_event_holders method to refer to the addition of Gz Event to alL_event_holders.)
                            gz = None
                            if ('gz_' + event_unique_name) in self.all_event_holders.keys():
                                gz = self.all_event_holders['gz_' + event_unique_name][0]
                            if gz is not None:
                                events_to_add.append(gz)

                            # Remove Event if it has to be added only once
                            if not event_def['include_in_loop']:
                                _events_blacklist.append(event_unique_name)
                                # TODO: Optimize code
                                # self.all_events_ordered[unique_node_name].remove(event_unique_name)
                                # _keys = list(self.all_event_holders.keys())
                                # self.log.node(_keys)
                                # _values = list(self.all_event_holders.values())
                                # _index = _keys.index(event_unique_name)
                                # _keys.pop(_index)
                                # _values.pop(_index)
                                # self.log.node(_keys)
                                # self.all_event_holders = OrderedDict(zip(_keys, _values))

                    self.seq.add_block(*events_to_add)

    def get_eventdef_for_unique_name(self, event_unique_name):
        """This method returns the Event definition from all_event_defs matching a given event_unique_name."""

        for event in self.all_event_defs:
            if event['event_unique_name'] == event_unique_name:
                return event

    def set_plot_outputs(self):
        t0, time_range = 0, [0, np.inf]
        adc_values = [[], []]
        rf_mag_values = [[], []]
        rf_phase_values = [[], []]
        t_x_values = [[], []]
        t_y_values = [[], []]
        t_z_values = [[], []]

        for iB in range(1, len(self.seq.block_events)):
            block = self.seq.get_block(iB)
            is_valid = time_range[0] <= t0 <= time_range[1]
            if is_valid:
                if block is not None:
                    if 'adc' in block:
                        adc = block['adc']
                        t = adc.delay + [(x * adc.dwell) for x in range(0, int(adc.num_samples))]
                        adc_values[0].extend(t0 + t)
                        adc_values[1].extend(np.ones(len(t)))
                    if 'rf' in block:
                        rf = block['rf']
                        t = rf.t
                        rf_mag_values[0].extend(t0 + t[0])
                        rf_mag_values[1].extend(abs(rf.signal))
                        rf_phase_values[0].extend(t0 + t[0])
                        rf_phase_values[1].extend(np.angle(rf.signal))
                    grad_channels = ['gx', 'gy', 'gz']
                    for x in range(0, len(grad_channels)):
                        if grad_channels[x] in block:
                            grad = block[grad_channels[x]]
                            if grad.type == 'grad':
                                t = grad.t
                                waveform = 1e-3 * grad.waveform
                            else:
                                t = np.cumsum([0, grad.rise_time, grad.flat_time, grad.fall_time])
                                waveform = [1e-3 * grad.amplitude * x for x in [0, 1, 1, 0]]
                                if grad.channel == 'x':
                                    t_x_values[0].extend(t0 + t)
                                    t_x_values[1].extend(waveform)
                                elif grad.channel == 'y':
                                    t_y_values[0].extend(t0 + t)
                                    t_y_values[1].extend(waveform)
                                elif grad.channel == 'z':
                                    t_z_values[0].extend(t0 + t)
                                    t_z_values[1].extend(waveform)
            t0 += calc_duration(block)

        # Setting outputs
        # ADC
        adc_output = np.array(adc_values)
        adc_output = adc_output.transpose()
        self.setData('adc_output', adc_output)

        # RF Mag
        rf_mag_output = np.array(rf_mag_values)
        rf_mag_output = rf_mag_output.transpose()
        self.setData('rf_mag_output', rf_mag_output)

        # RF Phase
        rf_ph_output = np.array(rf_phase_values)
        rf_ph_output = rf_ph_output.transpose()
        self.setData('rf_phase_output', rf_ph_output)

        # TrapX
        t_x_output = np.array(t_x_values)
        t_x_output = t_x_output.transpose()
        self.setData('trap_x_output', t_x_output)

        # TrapY
        t_y_output = np.array(t_y_values)
        t_y_output = t_y_output.transpose()
        self.setData('trap_y_output', t_y_output)

        # TrapZ
        t_z_output = np.array(t_z_values)
        t_z_output = t_z_output.transpose()
        self.setData('trap_z_output', t_z_output)
