import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LogNorm

from DebugUtils.GenerateRMFitView import makeWaveformArray
from DebugUtils.RecoMoreReader import readRecoMore, readWCWaveforms, RecoMoreEvent
from CReader import *


class RecoMoreFitExaminer:
    def __init__(self, rawDataPath: str, recoMoreDataPath: str):
        self.RMPEs = ReadRecoMoreOutput(recoMoreDataPath)
        self.rawWFs = ReadWCDataFile(rawDataPath)

        self.reducedChiSqs = []
        self.amps = []
        self.times = []
        for event_ in self.RMPEs:
            for channel_ in event_.SiPM:
                self.reducedChiSqs.append(channel_.redChiSq)
                if len(channel_.pes) > 0:
                    for PE in channel_.pes:
                        self.amps.append(PE.amplitude)
                        self.times.append(PE.time)

    def getEventPair(self, eventNumber: int, channelNumber: int):
        RMEvent = None
        rawEvent = None
        for event_ in self.RMPEs:
            if (eventNumber == event_.eventNumber) and (channelNumber == event_.channelNumber):
                RMEvent = event_

        for WF_ in self.rawWFs.getEvents():
            if WF_.eventID == eventNumber:
                for channelWF in WF_.chData:
                    if channelWF.channel == channelNumber:
                        rawEvent = channelWF

        return RMEvent, rawEvent

    def plotSingleEvent(self, eventNumber: int, channelNumber: int):
        RMEvent_, rawEvent_ = self.getEventPair(eventNumber, channelNumber)
        rawEvent_.xVals = [i * 0.3125 for i in range(len(rawEvent_.waveform))]
        xs, ys = makeWaveformArray(rawEvent_, RMEvent_)
        rawEvent_.xVals = np.array(rawEvent_.xVals)
        plt.plot(rawEvent_.xVals, rawEvent_.waveform)
        plt.title('Event: {}, Channel: {}'.format(RMEvent_.eventNumber,
                                                  RMEvent_.channelNumber))
        formattedPEList = ""
        if len(RMEvent_.PEData) > 0:
            formattedPEList += "\nAmp, time:"
        for idx, PE in enumerate(RMEvent_.PEData):
            formattedPEList += "\n{:.4f}, {:.3f}".format(PE.amplitude, PE.time)

        plt.plot(xs, ys, color='k', label='Fit:\nReduced ChiSq {:.6f}'
                                          '\nBaseline: {:.5f}'
                 .format(RMEvent_.redChiSq,
                         RMEvent_.baseline) + formattedPEList)

        plt.legend()
        plt.show()

    def plotAllEvents(self):
        eventChannelNumbers = [(event_.eventNumber, event_.channelNumber) for event_ in self.RMPEs]
        for event_ in eventChannelNumbers:
            self.plotSingleEvent(event_[0], event_[1])

    def plotAmps(self):
        plt.hist(self.amps, bins=1000)
        plt.xlabel("PE amplitudes (V)")
        plt.yscale('log')
        plt.show()

    def plotTimes(self):
        plt.hist(self.times, bins=1000)
        plt.xlabel("PE times (ns)")
        plt.show()

    def plotChiSq(self):
        plt.hist(self.reducedChiSqs, bins=1000)
        plt.xlabel(r"$\chi^2_r$")
        plt.show()

    def timeAmpCorrelation(self):
        times_ = []
        amps_ = []
        for event_ in self.RMPEs:
            for channel_ in event_.SiPM:
                self.reducedChiSqs.append(channel_.redChiSq)
                if len(channel_.pes) > 0:
                    for PE in channel_.pes:
                        amps_.append(PE.amplitude)
                        times_.append(PE.time - channel_.pes[0].time)

        # H, xEdges, yEdges = np.histogram2d(self.times, self.amps, bins=200)

        H, xEdges, yEdges = np.histogram2d(times_, amps_, bins=200)
        H = H.transpose()
        plt.imshow(H, origin='lower',
                   extent=[xEdges[0], xEdges[-1], yEdges[0], yEdges[-1]], aspect='auto', norm=LogNorm())
        plt.ylabel('Amplitude (V)')
        plt.xlabel('Time (ns)')
        plt.show()

    def plotSumAmps(self, channel: int):
        sumPES = []

        for event_ in self.RMPEs:
            for channel_ in event_.SiPM:
                if channel_.ch == channel:
                    runSum = 0
                    for PE in channel_.pes:
                        if PE.amplitude > 0.006:
                            runSum += PE.amplitude
                    if runSum > 0:
                        sumPES.append(runSum)
        plt.hist(sumPES, bins=300)
        plt.xlabel("Summed amplitude (V)")
        plt.show()


if __name__ == "__main__":
    # recoMoreFileName = "/Users/joshuaporter/OneDrive - University of Sussex/liquidOLab/data/WavecatcherRuns/Runs/R193/R193PES.dat"
    # rawFileName = "/Users/joshuaporter/OneDrive - University of Sussex/liquidOLab/data/WavecatcherRuns/Runs/R193/R193.bin"
    recoMoreFileName = "/Users/joshuaporter/OneDrive - University of Sussex/liquidOLab/data/WavecatcherRuns/Runs/R110/R110PES.dat"
    rawFileName = "/Users/joshuaporter/OneDrive - University of Sussex/liquidOLab/data/WavecatcherRuns/Runs/R110/R110.dat"

    examiner = RecoMoreFitExaminer(recoMoreDataPath=recoMoreFileName, rawDataPath=rawFileName)
    # examiner.plotAllEvents()
    examiner.plotSumAmps(3)
    examiner.timeAmpCorrelation()
    examiner.plotAmps()
    examiner.plotTimes()
    examiner.plotChiSq()
