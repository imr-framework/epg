from collections import OrderedDict

import gpi
from gpi import QtGui


class ExternalNode(gpi.NodeAPI):
    def initUI(self):
        # Widgets
        self.addWidget('StringBox', 'Unique Node Names containing EPI Events')
        self.addWidget('StringBox', 'Unique Event Name of EPI Event undergoing amplitude inversion')
        self.addWidget('PushButton', 'ComputeEvents', button_title="Compute events")

        # IO Ports
        self.addOutPort('output', 'DICT')

        return 0

    def compute(self):
        if 'ComputeEvents' in self.widgetEvents():
            epi_events = self.getVal('Unique Node Names containing EPI Events')
            epi_event_amp_inv = self.getVal('Unique Event Name of EPI Event undergoing amplitude inversion')
            self.setData('output', {'epi_events': epi_events, 'epi_event_amplitude_inversion': epi_event_amp_inv})

            return 0
