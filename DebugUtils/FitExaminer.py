import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm

from DebugUtils.CalibFitter import fit
from DebugUtils.GenerateRMFitView import makeWaveformArray
from CReader import *


class RecoMoreFitExaminer:
    def __init__(self, rawDataPath: str, recoMoreDataPath: str):
        """
        Utility for linking together raw and RecoMore data.
        :param rawDataPath: Path to the run's raw data file.
        :param recoMoreDataPath: Path to the run's processed RecoMore file.
        """
        self.RMPEs = ReadRecoMoreOutput(recoMoreDataPath)
        self.rawWFs = ReadWCDataFile(rawDataPath)

        self.reducedChiSqs = []
        self.amps = []
        self.times = []
        for event_ in self.RMPEs.getEvents():
            for channel_ in event_.SiPM:
                self.reducedChiSqs.append(channel_.redChiSq)
                if len(channel_.pes) > 0:
                    for PE_ in channel_.pes:
                        self.amps.append(PE_.amplitude)
                        self.times.append(PE_.time)

    def getEventPair(self, eventID: int, channelNumber: int):
        """
        Gets the RecoMore fit data and raw waveform for an event in a given channel.
        :param eventID:
        :param channelNumber:
        :return:
        """
        RMEvent = self.RMPEs.getEventChannel(eventID, channelNumber)
        rawEvent = self.rawWFs.getEventChannel(eventID, channelNumber)

        return RMEvent, rawEvent

    def plotSingleEvent(self, eventID: int, channelNumber: int):
        RMEvent_, rawEvent_ = self.getEventPair(eventID, channelNumber)
        rawEvent_.xVals = [i * 0.3125 for i in range(len(rawEvent_.waveform))]
        xs, ys = makeWaveformArray(rawEvent_, RMEvent_)
        rawEvent_.xVals = np.array(rawEvent_.xVals)
        plt.plot(rawEvent_.xVals, rawEvent_.waveform)
        plt.title('Event: {}, Channel: {}'.format(eventID,
                                                  RMEvent_.ch))
        formattedPEList = ""
        if len(RMEvent_.pes) > 0:
            formattedPEList += "\n\n$PEs$\nAmp (V), time (ns):"
        for idx, PE_ in enumerate(RMEvent_.pes):
            formattedPEList += "\n {:.4f},   {:.3f}".format(PE_.amplitude, PE_.time)

        plt.plot(xs, ys, color='k', label='$\chi^2_r$: {:.2f}'
                                          '\nBaseline: {:.2e} V'
                 .format(RMEvent_.redChiSq,
                         RMEvent_.baseline) + formattedPEList)
        plt.xlabel("Time (ns)")
        plt.ylabel("Amplitude (V)")
        plt.grid()
        plt.legend()
        plt.tight_layout()
        plt.show()

    def plotAllEvents(self, numPEThresh: int = 0):
        for event_ in self.RMPEs.getEvents():
            for channel in event_.SiPM:
                if len(channel.pes) > numPEThresh:
                    self.plotSingleEvent(eventID=event_.eventID, channelNumber=channel.ch)

    def plotAmps(self):
        plt.hist(self.amps, bins=300)
        plt.xlabel("PE amplitudes (V)")
        # plt.yscale('log')
        plt.grid()
        plt.tight_layout()
        plt.show()

    def plotTimes(self):
        plt.hist(self.times, bins=300)
        plt.xlabel("PE times (ns)")
        plt.tight_layout()
        plt.grid()
        plt.show()

    def plotChiSq(self):
        plt.hist(self.reducedChiSqs, bins=300)
        plt.xlabel(r"$\chi^2_r$")
        plt.tight_layout()
        plt.grid()
        plt.show()

    def timeAmpCorrelation(self, channel=None):
        times_ = []
        amps_ = []
        for event_ in self.RMPEs.getEvents():
            for channel_ in event_.SiPM:
                if channel is not None:
                    if channel_.ch != channel:
                        continue
                self.reducedChiSqs.append(channel_.redChiSq)
                if len(channel_.pes) > 0:
                    for PE_ in channel_.pes:
                        amps_.append(PE_.amplitude)
                        times_.append(PE_.time - channel_.pes[0].time)

        H, xEdges, yEdges = np.histogram2d(times_, amps_, bins=200)
        H = H.transpose()
        plt.imshow(H, origin='lower',
                   extent=[xEdges[0], xEdges[-1], yEdges[0], yEdges[-1]], aspect='auto', norm=LogNorm())
        plt.ylabel('Amplitude (V)')
        plt.xlabel('Time since first PE in waveform (ns)')
        plt.tight_layout()
        plt.show()

    def plotSumAmps(self, channel: int = None, PEThresh: float = 0):
        sumPES = {}
        minRunSum = 1000
        maxRunSum = 0

        for event_ in self.RMPEs.getEvents():
            for channel_ in event_.SiPM:
                if channel is not None:
                    if channel_.ch != channel:
                        continue
                runSum = 0
                for PE_ in channel_.pes:
                    if PE_.amplitude > PEThresh:
                        runSum += PE_.amplitude
                if runSum > 0:
                    if channel_.ch not in sumPES:
                        sumPES[channel_.ch] = []
                    sumPES[channel_.ch].append(runSum)

                    if runSum > maxRunSum:
                        maxRunSum = runSum
                    if runSum < minRunSum:
                        minRunSum = runSum

        bins = np.linspace(minRunSum, maxRunSum, 600)
        for key, val in sumPES.items():
            fit(val)

        for key, val in sumPES.items():
            plt.hist(val, bins=bins, label='{}'.format(key), histtype='step')

        plt.xlabel("Summed amplitude (V)")
        plt.legend()
        plt.tight_layout()
        plt.show()


if __name__ == "__main__":
    # recoMoreFileName = "/Users/joshuaporter/Library/CloudStorage/OneDrive-UniversityofSussex/liquidOLab/dataSTOP_DO_NOT_WRITE_HERE/WavecatcherRuns/Runs/R193/R193PES.dat"
    # rawFileName = "/Users/joshuaporter/Library/CloudStorage/OneDrive-UniversityofSussex/liquidOLab/dataSTOP_DO_NOT_WRITE_HERE/WavecatcherRuns/Runs/R193/R193.bin"
    matplotlib.rcParams.update({'font.size': 16})
    matplotlib.use("TKAgg", force=True)
    recoMoreFileName = "/home/joshuap/Software/JoshSoftware/RecoMore/R185PES.dat"
    rawFileName = "/home/joshuap/Software/JoshSoftware/RecoMore/R185.bin"

    examiner = RecoMoreFitExaminer(recoMoreDataPath=recoMoreFileName, rawDataPath=rawFileName)
    # examiner.plotAllEvents()
    # examiner.plotSumAmps(PEThresh=0.008)
    examiner.timeAmpCorrelation()
    examiner.plotAmps()
    examiner.plotTimes()
    examiner.plotChiSq()
