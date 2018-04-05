import os

import gpi
import matplotlib.pyplot as plt
import numpy as np
from gpi import QtGui

import pulseq.core


class ExternalNode(gpi.NodeAPI):
    def initUI(self):
        # Widgets
        self.addWidget('DisplayBox', 'Sequence Plot')

        # IO Ports
        self.addInPort('input', 'DICT')

    def compute(self):
        self.plot()
        return 0

    def plot(self):
        seq = self.getData('input')['seq']
        fig = plt.figure(figsize=(7, 10))
        f11, f12, f13 = fig.add_subplot(611), fig.add_subplot(612), fig.add_subplot(613)
        f2 = [fig.add_subplot(614), fig.add_subplot(615), fig.add_subplot(616)]
        t0, time_range = 0, [0, np.inf]
        for iB in range(1, len(seq.block_events)):
            block = seq.get_block(iB)
            is_valid = time_range[0] <= t0 <= time_range[1]
            if is_valid:
                if block is not None:
                    if 'adc' in block:
                        adc = block['adc']
                        t = adc.delay + [(x * adc.dwell) for x in range(0, int(adc.num_samples))]
                        f11.plot((t0 + t), np.zeros(len(t)))
                    if 'rf' in block:
                        rf = block['rf']
                        t = rf.t
                        f12.plot(np.squeeze(t0 + t), abs(rf.signal))
                        f13.plot(np.squeeze(t0 + t), np.angle(rf.signal))
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
                            f2[x].plot(np.squeeze(t0 + t), waveform)
            t0 += core.calcduration.calcduration(block)

        f11.set_ylabel('adc')
        f12.set_ylabel('rf mag hz')
        f13.set_ylabel('rf phase rad')
        [f2[x].set_ylabel(grad_channels[x]) for x in range(3)]

        # Save the plot temporarily as a PNG file to be displayed by the DisplayBox Widget
        plt.savefig('plot_temp.png')
        img = QtGui.QImage()
        img.load('plot_temp.png')
        self.setAttr('Sequence Plot', val=img)
        os.remove('plot_temp.png')

    def execType(self):
        """Could be GPI_THREAD, GPI_PROCESS, GPI_APPLOOP"""
        return gpi.GPI_APPLOOP
