import gpi


class ExternalNode(gpi.NodeAPI):
    """This node lets the user specify a file name and save location to write open-source file format seq files."""

    def initUI(self):
        # IO Ports
        self.addInPort(title='input', type='DICT')

        # Widgets
        self.addWidget('SaveFileBrowser', 'File location', button_title='Browse')
        self.addWidget('PushButton', 'Write seq file', button_title='Save now')

        return 0

    def compute(self):
        if 'Write seq file' in self.widgetEvents():
            in_dict = self.getData('input')
            seq = in_dict['seq']
            file_location = self.getVal('File location')
            seq.write(file_location)

            return 0
