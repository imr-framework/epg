import gpi

from pulseq.core.Sequence.sequence import Sequence
from pulseq.core.opts import Opts


class ExternalNode(gpi.NodeAPI):
    """This node lets the user specify a file name and save location to write open-source file format seq files.    """

    def initUI(self):
        # IO Ports
        self.addOutPort(title='output', type='DICT')

        # Widgets
        self.addWidget('OpenFileBrowser', 'File location', button_title='Browse')
        return 0

    def compute(self):
        if 'File location' in self.widgetEvents():
            file_location = self.getVal('File location')
            seq = Sequence(Opts())
            seq.read(file_location)
            self.setData('output', {"seq": seq})
