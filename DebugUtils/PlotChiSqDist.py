import matplotlib.pyplot as plt
import numpy as np

from DebugUtils.GenerateRMFitView import makeWaveformArray
from DebugUtils.RecoMoreReader import readRecoMore, readWCWaveforms


class RecoMoreFitExaminer:
    def __init__(self, rawDataPath: str, recoMoreDataPath: str):
        with open(recoMoreDataPath, 'r') as RMFile:
            self.RMPEs = readRecoMore(RMFile)

        with open(rawDataPath, 'r') as rawFile:
            self.rawWFs = readWCWaveforms(rawFile)

        self.reducedChiSqs = []
        self.amps = []
        self.times = []
        for event_ in self.RMPEs:
            if len(event_.PEData) > 0:
                self.reducedChiSqs.append(event_.redChiSq)
                for PE in event_.PEData:
                    self.amps.append(PE.amplitude)
                    self.times.append(PE.time)

    def getEventPair(self, eventNumber: int, channelNumber: int):
        RMEvent = None
        rawEvent = None
        for event_ in self.RMPEs:
            if (eventNumber == event_.eventNumber) and (channelNumber == event_.channelNumber):
                RMEvent = event_

        for WF_ in self.rawWFs:
            if (WF_.eventNumber == eventNumber) and (WF_.channelNumber == channelNumber):
                rawEvent = WF_

        return RMEvent, rawEvent

    def plotSingleEvent(self, eventNumber: int, channelNumber: int):
        RMEvent_, rawEvent_ = self.getEventPair(eventNumber, channelNumber)
        xs, ys = makeWaveformArray(rawEvent_, RMEvent_)
        rawEvent_.xVals = np.array(rawEvent_.xVals)
        plt.plot(rawEvent_.xVals, rawEvent_.rawData)
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
        plt.show()

    def plotTimes(self):
        plt.hist(self.times, bins=1000)
        plt.xlabel("PE times (ns)")
        plt.show()

    def plotChiSq(self):
        plt.hist(self.reducedChiSqs, bins=1000)
        plt.xlabel(r"$\chi^2_r$")
        plt.show()


if __name__ == "__main__":
    recoMoreFileName = "/Users/joshuaporter/OneDrive - University of Sussex/liquidOLab/data/WavecatcherRuns/Runs/R110/R110PES.dat"
    rawFileName = "/Users/joshuaporter/OneDrive - University of Sussex/liquidOLab/data/WavecatcherRuns/Runs/R110/R110.dat"

    examiner = RecoMoreFitExaminer(recoMoreDataPath=recoMoreFileName, rawDataPath=rawFileName)
    # examiner.plotAllEvents()
    examiner.plotAmps()
    examiner.plotTimes()
    examiner.plotChiSq()
