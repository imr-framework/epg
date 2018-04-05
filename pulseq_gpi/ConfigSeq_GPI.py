import gpi

from pulseq.core.opts import Opts


class ExternalNode(gpi.NodeAPI):
    """This node providers options for configuring the properties of the pulse sequence."""

    def initUI(self):
        # Widgets
        self.addWidget('StringBox', 'Maximum Gradient', placeholder="max_grad")
        self.addWidget('ComboBox', 'Maximum Gradient Unit', placeholder="max_grad_unit",
                       items=['Hz/m', 'mT/m', 'rad/ms/mm'], val='mT/m')
        self.addWidget('StringBox', 'Maximum Slew Rate', placeholder="max_slew")
        self.addWidget('ComboBox', 'Maximum Slew Rate Unit', placeholder="max_slew_unit",
                       items=['Hz/m/s', 'mT/m/ms', 'T/m/s', 'rad/ms/mm/ms'], val='T/m/s')
        self.addWidget('StringBox', 'Repetition Time (s)', placeholder="TR")
        self.addWidget('StringBox', 'Echo Time (s)', placeholder="TE")
        self.addWidget('StringBox', 'Field of View', placeholder="fov")
        self.addWidget('StringBox', 'Nx', placeholder="Nx")
        self.addWidget('StringBox', 'Ny', placeholder="Ny")
        self.addWidget('StringBox', 'Rise Time (s)', placeholder="rise_time")
        self.addWidget('StringBox', 'RF Dead Time (s)', placeholder="rf_dead_time")
        self.addWidget('StringBox', 'ADC Dead Time (s)', placeholder="adc_deadTime")
        self.addWidget('StringBox', 'RF Raster Time (s)', placeholder="rf_raster")
        self.addWidget('StringBox', 'Gradient Raster Time (s)', placeholder="grad_raster")
        self.addWidget('TextBox', 'System limits')
        self.addWidget('PushButton', 'ComputeEvents', button_title="Compute events")

        # IO Ports
        self.addOutPort('output', 'DICT')

        return 0

    def compute(self):
        if 'ComputeEvents' in self.widgetEvents() or 'input' in self.portEvents():
            out_dict = {}
            max_grad = self.getVal('Maximum Gradient')
            max_grad_unit = self.getVal('Maximum Gradient Unit')
            max_slew = self.getVal('Maximum Slew Rate')
            max_slew_unit = self.getVal('Maximum Slew Rate Unit')
            tr = self.getVal('Repetition Time (s)')
            te = self.getVal('Echo Time (s)')
            fov = self.getVal('Field of View')
            Nx = self.getVal('Nx')
            Ny = self.getVal('Ny')
            rise_time = self.getVal('Rise Time (s)')
            rf_dead_time = self.getVal('RF Dead Time (s)')
            adc_dead_time = self.getVal('ADC Dead Time (s)')
            rf_raster_time = self.getVal('RF Raster Time (s)')
            rf_raster_time = rf_raster_time if rf_raster_time != '' else 1e-6
            grad_raster_time = self.getVal('Gradient Raster Time (s)')
            grad_raster_time = grad_raster_time if grad_raster_time != '' else 10e-6

            kwargs_for_opts = {'grad_unit': max_grad_unit, 'slew_unit': max_slew_unit}
            keys_for_opts = ['max_grad', 'max_slew', 'tr', 'te', 'fov', 'Nx', 'Ny',
                             'rise_tme', 'rf_dead_time', 'adc_dead_time', 'rf_raster_time', 'grad_raster_time']
            values_for_opts = [max_grad, max_slew, tr, te, fov, Nx, Ny, rise_time,
                               rf_dead_time, adc_dead_time, rf_raster_time, grad_raster_time]
            for i in range(len(values_for_opts)):
                kwargs_for_opts[keys_for_opts[i]] = float(values_for_opts[i])

            system = Opts(kwargs_for_opts)
            out_dict['system'] = system

            self.setData('output', out_dict)

            # To display the computed info in the TextBox
            self.setAttr('System limits', val=str(system))

            return 0
