from collections import OrderedDict

import gpi
from gpi import QtGui


class AddBlockWidgets(gpi.GenericWidgetGroup):
    """A unique widget that display a variable number of StringBoxes (or FileBrowser) depending on the Event being
    configured."""

    valueChanged = gpi.Signal()

    def __init__(self, title, parent=None):
        super(AddBlockWidgets, self).__init__(title, parent)
        self.button_names_list = ['Off', 'Delay', 'SincRF', 'BlockRF', 'G', 'GyPre', 'ArbGrad', 'ADC']
        self.clicked_button_name, self.clicked_button_index = '', 0
        self.buttons_list, self.string_box_list = [], []

        # Labels for StringBoxes to configure Events
        self.delay_labels = ['Unique Event name', 'Delay (s)']
        self.sinc_rf_labels = ['Unique Event name', 'Maximum Gradient (mT/m)', 'Maximum Slew (T/m/s)',
                               'Flip Angle (deg)', 'Duration (s)', 'Frequency Offset', 'Phase Offset',
                               'Time Bw Product', 'Apodization', 'Slice Thickness (m)']
        self.block_rf_labels = ['Unique Event name', 'Maximum Gradient (mT/m)', 'Maximum Slew (T/m/s)',
                                'Flip Angle (deg)', 'Duration (s)', 'Frequency Offset', 'Phase Offset',
                                'Time Bw Product', 'Bandwidth', 'Slice Thickness (m)']
        self.trap_labels = ['Unique Event name', 'Channel', 'Maximum Gradient (mT/m)', 'Maximum Slew (T/m/s)',
                            'Duration (s)', 'Area', 'Flat Time (s)', 'Flat Area', 'Amplitude (Hz)', 'Rise Time (s)']
        self.gy_pre_labels = ['Unique Event name', 'Duration (s)', 'Area']
        self.arb_grad_labels = ['Unique Event name', 'Channel', 'Maximum Gradient (mT/m)', 'Maximum Slew (T/m/s)']
        self.adc_labels = ['Unique Event name', 'Number of samples', 'Dwell (s)', 'Duration (s)', 'Delay (s)',
                           'Frequency Offset', 'Phase Offset']
        # Placeholders for StringBoxes to configure Events
        # The placeholders are also the keys to retrieve values from the StringBoxes
        self.delay_placeholders = ['event_unique_name', 'delay']
        self.sinc_rf_placeholders = ['event_unique_name', 'max_grad', 'max_slew', 'flip_angle', 'duration',
                                     'freq_offset', 'phase_offset', 'time_bw_prod', 'apodization', 'slice_thickness']
        self.block_rf_placeholders = ['event_unique_name', 'max_grad', 'max_slew', 'flip_angle', 'duration',
                                      'freq_offset', 'phase_offset', 'time_bw_prod', 'bandwidth', 'slice_thickness']
        self.trap_placeholders = ['event_unique_name', 'channel', 'max_grad', 'max_slew', 'duration', 'area',
                                  'flat_time', 'flat_area', 'amplitude', 'rise_time']
        self.gy_pre_placeholders = ['event_unique_name', 'duration', 'area']
        self.arb_grad_placeholders = ['event_unique_name', 'channel', 'max_grad', 'max_slew']
        self.adc_placeholders = ['event_unique_name', 'num_samples', 'dwell', 'duration', 'delay', 'freq_offset',
                                 'phase_offset']

        # Variable to denote the maximum number of StringBoxes to be added; obviously depends on the Event which has the
        # maximum number of configuration parameters
        self.num_string_boxes = max(len(self.delay_labels), len(self.sinc_rf_labels), len(self.block_rf_labels),
                                    len(self.trap_labels), len(self.gy_pre_labels), len(self.arb_grad_labels),
                                    len(self.adc_labels))

        # First index is None because the first button is 'Off'. Look into event_def['event_values'] in get_val()
        self.labels = [None, self.delay_labels, self.sinc_rf_labels, self.block_rf_labels, self.trap_labels,
                       self.gy_pre_labels, self.arb_grad_labels, self.adc_labels]
        self.placeholders = [None, self.delay_placeholders, self.sinc_rf_placeholders, self.block_rf_placeholders,
                             self.trap_placeholders, self.gy_pre_placeholders, self.arb_grad_placeholders,
                             self.adc_placeholders]

        self.wdg_layout = QtGui.QGridLayout()
        self.add_event_pushbuttons()
        self.add_config_stringboxes()
        self.add_include_in_loop_pushbutton()
        self.add_include_gz_pushbutton()
        self.add_file_browser()

        self.setLayout(self.wdg_layout)
        self.buttons_list[0].setChecked(True)

    def add_event_pushbuttons(self):
        """Adding PushButtons for the Events."""
        col_count = 0
        for name in self.button_names_list:
            new_button = QtGui.QPushButton(name)
            new_button.setCheckable(True)
            new_button.setAutoExclusive(True)
            new_button.clicked.connect(self.button_clicked)
            new_button.clicked.connect(self.valueChanged)
            # Syntax: addWidget(widget, row, col, rowSpan, colSpan)
            self.wdg_layout.addWidget(new_button, 0, col_count, 1, 1)
            self.buttons_list.append(new_button)
            col_count += 1

    def add_config_stringboxes(self):
        """Adding StringBoxes for configuring the Events."""
        for x in range(self.num_string_boxes):
            string_box = gpi.StringBox(str(x))
            string_box.set_visible(False)
            # Syntax: addWidget(widget, row, col, rowSpan, colSpan)
            self.wdg_layout.addWidget(string_box, x + 1, 1, 1, 6)
            self.string_box_list.append(string_box)

    def add_include_in_loop_pushbutton(self):
        """Adding PushButton to toggle Event being included/excluded in loop."""
        self.include_in_loop_pushbutton = QtGui.QPushButton('Add event in loop')
        self.include_in_loop_pushbutton.setCheckable(True)
        self.include_in_loop_pushbutton.setChecked(True)
        self.include_in_loop_pushbutton.setVisible(False)
        self.include_in_loop_pushbutton.clicked.connect(self.button_clicked)
        self.include_in_loop_pushbutton.clicked.connect(self.valueChanged)
        # Syntax: addWidget(widget, row, col, rowSpan, colSpan)
        self.wdg_layout.addWidget(self.include_in_loop_pushbutton, 11, 1, 1, 6)

    def add_include_gz_pushbutton(self):
        """Adding PushButton toggle for Gz along with Rf."""
        self.include_gz_pushbutton = QtGui.QPushButton('Add Gz event with Rf')
        self.include_gz_pushbutton.setCheckable(True)
        self.include_gz_pushbutton.setVisible(False)
        self.include_gz_pushbutton.clicked.connect(self.button_clicked)
        self.include_gz_pushbutton.clicked.connect(self.valueChanged)
        # Syntax: addWidget(widget, row, col, rowSpan, colSpan)
        self.wdg_layout.addWidget(self.include_gz_pushbutton, 12, 1, 1, 6)

    def add_file_browser(self):
        """Adding FileBrowser necessary for ArbGrad Event."""
        self.file_browser = gpi.OpenFileBrowser('Read .hdf5')
        self.file_browser.set_button_title('Read .hdf5')
        self.file_browser.set_visible(False)
        # Syntax: addWidget(widget, row, col, rowSpan, colSpan)
        self.wdg_layout.addWidget(self.file_browser, 5, 1, 1, 6)

    # Getter
    def get_val(self):
        if self.clicked_button_index == self.button_names_list.index('Off'):
            # 'Off' PushButton selected, return empty dict
            return {}
        else:
            """
            event_def contains:
            - event_name : str
                Event name, corresponds to Event button that is selected. Required to correctly construct the Event 
                in GenSeq_GPI Node
            - event_unique_name : str
                Unique Event name; user input
            - event_values : OrderedDict
                key-value pairs of Event parameters_params and values
            - include_in_loop : bool
                If Event should be added to Sequence Ny times
            - include_gz : bool
                If Gz Event should be added to Sequence with Rf Event
            - file_path : str
                Path to .hdf5 file required for arbitrary gradient Event
            """
            event_def = dict()
            keys = self.placeholders[self.clicked_button_index][1:]
            values = [x.get_val() for x in self.string_box_list[1:]]
            for x in values:
                # Take default 0 if user has not entered any value
                if x == '':
                    values[values.index(x)] = '0'
            event_def['event_values'] = OrderedDict(zip(keys, values))
            event_def['event_name'] = self.clicked_button_name
            event_def['event_unique_name'] = self.string_box_list[0].get_val()
            if event_def['event_unique_name'] == '':
                # Return None to raise an error in compute()
                return None

            event_def['include_in_loop'] = self.include_in_loop_pushbutton.isChecked()
            if self.clicked_button_index == 2 or self.clicked_button_index == 3:
                # For Rf event, check if Gz has to be included
                event_def['include_gz'] = self.include_gz_pushbutton.isChecked()
            elif self.clicked_button_index == 6:
                # For arbitrary gradient event, retrieve file path
                event_def['file_path'] = self.file_browser.get_val()
            return event_def

    # Setter
    def set_val(self, val):
        self.hide_config_widgets()
        if len(val) != 0:
            event_values = val['event_values']
            event_values['event_unique_name'] = val['event_unique_name']
            self.clicked_button_name = val['event_name']
            self.clicked_button_index = self.button_names_list.index(self.clicked_button_name)
            self.buttons_list[self.clicked_button_index].setChecked(True)
            self.show_config_widgets()
            self.include_in_loop_pushbutton.setChecked(bool(val['include_in_loop']))
            if 'include_gz' in val:
                self.include_gz_pushbutton.setChecked(val['include_gz'])
            labels = self.labels[self.clicked_button_index]
            placeholders = self.placeholders[self.clicked_button_index]
            for x in range(len(placeholders)):
                self.string_box_list[x].setTitle(labels[x])
                self.string_box_list[x].set_placeholder(placeholders[x])
                self.string_box_list[x].set_val(event_values[placeholders[x]])

    def button_clicked(self):
        """Identifies the button that was clicked and stores the name and ID of the button."""
        for button in self.buttons_list:
            if button.isChecked():
                self.clicked_button_index = self.buttons_list.index(button)
                self.clicked_button_name = self.button_names_list[self.clicked_button_index]
        self.show_config_widgets()

    def show_config_widgets(self):
        """Show appropriate number of StringBoxes and relevant Widgets based on the button that was clicked."""
        self.hide_config_widgets()

        if self.clicked_button_index != 0 and self.clicked_button_index != 5:
            self.include_in_loop_pushbutton.setVisible(True)
        if self.clicked_button_index == 1:
            # Delay
            [self.string_box_list[x].set_visible(True) for x in range(len(self.delay_placeholders))]
            [self.string_box_list[x].setTitle(self.delay_labels[x]) for x in range(len(self.delay_labels))]
            [self.string_box_list[x].set_placeholder(self.delay_placeholders[x]) for x in
             range(len(self.delay_placeholders))]
        elif self.clicked_button_index == 2:
            # SincRF
            [self.string_box_list[x].set_visible(True) for x in range(len(self.sinc_rf_placeholders))]
            [self.string_box_list[x].setTitle(self.sinc_rf_labels[x]) for x in range(len(self.sinc_rf_labels))]
            [self.string_box_list[x].set_placeholder(self.sinc_rf_placeholders[x]) for x in
             range(len(self.sinc_rf_placeholders))]
            self.include_gz_pushbutton.setVisible(True)
        elif self.clicked_button_index == 3:
            # BlockRF
            [self.string_box_list[x].set_visible(True) for x in range(len(self.block_rf_placeholders))]
            [self.string_box_list[x].setTitle(self.block_rf_labels[x]) for x in range(len(self.block_rf_labels))]
            [self.string_box_list[x].set_placeholder(self.block_rf_placeholders[x]) for x in
             range(len(self.block_rf_placeholders))]
            self.include_gz_pushbutton.setVisible(True)
        elif self.clicked_button_index == 3:
            # G
            [self.string_box_list[x].set_visible(True) for x in range(len(self.trap_placeholders))]
            [self.string_box_list[x].setTitle(self.trap_labels[x]) for x in range(len(self.trap_labels))]
            [self.string_box_list[x].set_placeholder(self.trap_placeholders[x]) for x in
             range(len(self.trap_placeholders))]
        elif self.clicked_button_index == 4:
            # G
            [self.string_box_list[x].set_visible(True) for x in range(len(self.trap_placeholders))]
            [self.string_box_list[x].setTitle(self.trap_labels[x]) for x in range(len(self.trap_labels))]
            [self.string_box_list[x].set_placeholder(self.trap_placeholders[x]) for x in
             range(len(self.trap_placeholders))]
        elif self.clicked_button_index == 5:
            # GyPre
            [self.string_box_list[x].set_visible(True) for x in range(len(self.gy_pre_placeholders))]
            [self.string_box_list[x].setTitle(self.gy_pre_labels[x]) for x in range(len(self.gy_pre_labels))]
            [self.string_box_list[x].set_placeholder(self.gy_pre_placeholders[x]) for x in
             range(len(self.gy_pre_placeholders))]
        elif self.clicked_button_index == 6:
            # Arbitrary Grad
            [self.string_box_list[x].set_visible(True) for x in range(len(self.arb_grad_placeholders))]
            [self.string_box_list[x].setTitle(self.arb_grad_labels[x]) for x in range(len(self.arb_grad_labels))]
            [self.string_box_list[x].set_placeholder(self.arb_grad_placeholders[x]) for x in
             range(len(self.arb_grad_placeholders))]
            self.file_browser.set_visible(True)
        elif self.clicked_button_index == 7:
            # ADC
            [self.string_box_list[x].set_visible(True) for x in range(len(self.adc_placeholders))]
            [self.string_box_list[x].setTitle(self.adc_labels[x]) for x in range(len(self.adc_labels))]
            [self.string_box_list[x].set_placeholder(self.adc_placeholders[x]) for x in
             range(len(self.adc_placeholders))]

    def hide_config_widgets(self):
        """Hide all Widgets."""
        [x.set_visible(False) for x in self.string_box_list]
        [x.set_val("") for x in self.string_box_list]
        self.include_in_loop_pushbutton.setVisible(False)
        self.include_gz_pushbutton.setVisible(False)
        self.file_browser.set_visible(False)


class ExternalNode(gpi.NodeAPI):
    """
    This Node provides options for setting up the event that needs to be added. Event parameters_params should be set
    to 0 if left unconfigured. Up to 6 simultaneous events can be added in one block. The 'ComputeEvents' button gathers
    the input data into a dict object. The output of this Node (or a chain of AddBlock Nodes) has to be supplied to a
    GenSeq Node.

     Units:
     - duration: s
     - flip_angle : deg
     - flatTime : s
     - riseTime : s
     - dwell : s
     - delay : s
     - sliceThickness : m
     - amplitude : Hz
    """

    def initUI(self):
        # Init constant(s)
        self.num_concurrent_events = 6

        # Widgets
        self.addWidget('StringBox', 'Unique Node Name')
        for x in range(self.num_concurrent_events):
            self.addWidget('AddBlockWidgets', 'Event ' + str(x + 1))
        self.addWidget('PushButton', 'ComputeEvents', button_title="Compute events")

        # IO Ports
        self.addInPort('input', 'DICT')
        self.addOutPort('output', 'DICT')

        return 0

    def validate(self):
        if 'ComputeEvents' in self.widgetEvents() or 'input' in self.portEvents():
            self.unique_node_name = self.getVal('Unique Node Name')
            self.setDetailLabel(self.unique_node_name)

    def compute(self):
        """
        1. all_events_ordered is an OrderedDict, to preserve the order of user-configured Events.
        2. For unique_node_name in all_events_ordered, retrieve ordered_events; this will be a sequence of
            simultaneously occurring Events.
        3. Iterate through ordered_events, and retrieve corresponding event definitions from all_event_defs

        Structure
        ---------
        in_dict:
        - all_event_defs: [key-value pairs of event definitions]
        - all_events_ordered: {unique_node_name: [event_unique_name1, event_unique_name2 ...], ...}

        Variables
        ---------
        all_event_defs : list
            List of dict objects that are key-value pairs of event definitions.
        all_events_ordered: OrderedDict()
            OrderedDict consisting of key-value pairs of unique_node_name and ordered_events (see below).
        ordered_events : list
            List containing event_unique_name of all Events configured by the user in this Node. The order in which the
            user configured the Events is preserved.
        """
        if 'ComputeEvents' in self.widgetEvents() or 'input' in self.portEvents():
            in_dict = self.getData('input')
            all_event_defs = in_dict['all_event_defs'] if 'all_event_defs' in in_dict else []
            all_events_ordered = in_dict['all_events_ordered'] if 'all_events_ordered' in in_dict else OrderedDict()
            ordered_events = []

            for x in range(self.num_concurrent_events):
                event_def = self.getVal('Event ' + str(x + 1))
                # If event_def is None, raise a flag prompting the user to enter Unique Event Name.
                if event_def is None:
                    self.log.critical('Enter Unique Event Name')
                    return 1
                if len(event_def) != 0:
                    all_event_defs.append(event_def)
                    ordered_events.append(event_def['event_unique_name'])

            all_events_ordered[self.unique_node_name] = ordered_events

            in_dict['all_event_defs'] = all_event_defs
            in_dict['all_events_ordered'] = all_events_ordered
            self.setData('output', in_dict)
            self.log.node(in_dict)

            return 0
