import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm

from DebugUtils.CalibFitter import fit
from DebugUtils.GenerateRMFitView import makeWaveformArray
from CReader import *
import matplotlib.font_manager
import matplotlib as mpl


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
        fig, (ax1, ax2) = plt.subplots(nrows=2, sharex=True, figsize=(12, 9), height_ratios=[2, 1])

        ax1.plot(rawEvent_.xVals, rawEvent_.waveform, lw=3)
        plt.suptitle('Event: {}, Channel: {}'.format(eventID,
                                                     RMEvent_.ch))
        formattedPEList = ""
        if len(RMEvent_.pes) > 0:
            formattedPEList += "\n\n$PEs$\nAmp (V), time (ns):"
        for idx, PE_ in enumerate(RMEvent_.pes):
            formattedPEList += "\n {:.4f},   {:.3f}".format(PE_.amplitude, PE_.time)

        ax1.plot(xs, ys, color='k', label='$\chi^2_r$: {:.2f}'
                                          '\nBaseline: {:.2e} V'
                 .format(RMEvent_.redChiSq,
                         RMEvent_.baseline) + formattedPEList)
        ax2.set_xlabel("Time (ns)")
        ax1.set_ylabel("Voltage (V)")
        ax1.grid()
        ax1.legend(fontsize=16)

        ax2.plot(xs, ys - rawEvent_.waveform, color='k')
        ax2.set_ylabel("Residuals (V)")
        ax2.grid()

        plt.subplots_adjust(hspace=0.)
        plt.tight_layout()
        # plt.savefig("Event_{}_Channel_{}.png".format(eventID, RMEvent_.ch), dpi=300, bbox_inches='tight')
        # plt.close()
        plt.show()

    def plotAllEvents(self, numPEThresh: int = 0):
        for event_ in self.RMPEs.getEvents():
            for channel in event_.SiPM:
                if len(channel.pes) > numPEThresh:
                    self.plotSingleEvent(eventID=event_.eventID, channelNumber=channel.ch)

    def plotAmps(self):
        plt.figure(figsize=(12, 9))
        plt.hist(self.amps, bins=1300)
        plt.xlabel("PE amplitudes (V)")
        # plt.yscale('log')
        plt.grid()
        plt.tight_layout()
        plt.show()

    def plotTimes(self):
        plt.figure(figsize=(12, 9))
        plt.hist(self.times, bins=300)
        plt.xlabel("PE times (ns)")
        plt.tight_layout()
        plt.grid()
        plt.show()

    def plotChiSq(self):
        plt.figure(figsize=(12, 9))
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

        H, xEdges, yEdges = np.histogram2d(times_, amps_, bins=(500, np.linspace(0, 0.1, 200)))
        H = H.transpose()
        plt.figure(figsize=(12, 9))
        # plt.ylim([0, 0.1])
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

        plt.figure(figsize=(12, 9))
        for key, val in sumPES.items():
            plt.hist(val, bins=bins, label='{}'.format(key), histtype='step')

        plt.xlabel("Summed amplitude (V)")
        plt.legend()
        plt.grid()
        plt.tight_layout()

        plt.show()

    def plotSumAmpsRaw(self, channel: int = None, PEThresh: float = 0):
        sumPES = {}
        minRunSum = 1000
        maxRunSum = -1000

        for event_ in self.rawWFs.getEvents():
            for channel_ in event_.chData:
                if channel is not None:
                    if channel_.channel != channel:
                        continue
                runSum = np.sum(channel_.waveform)
                # if runSum > 0:
                if channel_.channel not in sumPES:
                    sumPES[channel_.channel] = []
                sumPES[channel_.channel].append(runSum)

                if runSum > maxRunSum:
                    maxRunSum = runSum
                if runSum < minRunSum:
                    minRunSum = runSum

        bins = np.linspace(minRunSum, maxRunSum, 200)
        # for key, val in sumPES.items():
        #     fit(val)

        plt.figure(figsize=(12, 9))
        for key, val in sumPES.items():
            plt.hist(val, bins=bins, label='{}'.format(key), histtype='step')

        plt.xlabel("Summed amplitude (V)")
        plt.legend()
        plt.grid()
        plt.tight_layout()

        plt.show()

    def recoMoreVsRawComparison(self, channel: int = None, PEThresh: float = 0, cutMode: str = 'time', runName: str = ''):
        colour = '#164ba1'
        binCount = 1600

        # Check that cutMode is either 'none', 'time', 'negativeAmp' or 'timeAndNegativeAmp'
        if cutMode not in ['none', 'time', 'negativeAmp', 'timeAndNegativeAmp']:
            print("Invalid cutMode. Please choose from 'none', 'time', 'negativeAmp' or 'timeAndNegativeAmp'.")
            return

        if (cutMode == 'time') or (cutMode == 'timeAndNegativeAmp'):
            lowerTimeCut = 30
            upperTimeCut = 70

        # RecoMore
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
                        timeDiff = PE_.time - channel_.pes[0].time

                        if (cutMode == 'time') or (cutMode == 'timeAndNegativeAmp'):
                            if (timeDiff > -lowerTimeCut) and (timeDiff < upperTimeCut):
                                runSum += PE_.amplitude
                        else:
                            runSum += PE_.amplitude

                if runSum > 0:
                    if channel_.ch not in sumPES:
                        sumPES[channel_.ch] = []
                    sumPES[channel_.ch].append(runSum)

                    if runSum > maxRunSum:
                        maxRunSum = runSum
                    if runSum < minRunSum:
                        minRunSum = runSum

        bins = np.linspace(minRunSum, maxRunSum, binCount)

        fig, (ax1, ax2) = plt.subplots(nrows=2, figsize=(12, 9), height_ratios=[1, 1])
        for key, val in sumPES.items():
            n, bins_ = np.histogram(val, bins=bins)
            bins_ = bins_[:-1]/np.max(bins_[:-1])
            # bins = bins - np.min(bins)
            ax1.fill_between(bins_, n, alpha=0.4, step='pre', color=colour)
            ax1.step(bins_, n, label='RecoMore', lw=2, color='k')
            # plt.hist(val, bins=bins, label='{} RM'.format(key), histtype='step', lw=2)

            # Set upper x limit to where there are no more bins with a height greater than 25.
            upperXLim = bins_[np.where(n > 25)[0][-1]]
            lowerXLim = bins_[np.where(n > 25)[0][0]]
            ax1.set_xlim([lowerXLim, upperXLim])

        # Raw integration
        sumPESRAW = {}
        minRunSumRAW = 1000
        maxRunSumRAW = 0

        for event_ in self.rawWFs.getEvents():
            for channel_ in event_.chData:
                if channel is not None:
                    if channel_.channel != channel:
                        continue
                rawWF = np.array(channel_.waveform)

                if cutMode == 'none':
                    runSum = -np.sum(rawWF)
                elif cutMode == 'negativeAmp':
                    negAmpMask = rawWF < 0
                    runSum = -np.sum(rawWF[negAmpMask])
                elif cutMode == 'time':
                    minAmpTimeIDX = np.argmin(rawWF)
                    lowerIndex = minAmpTimeIDX - (int(lowerTimeCut / 0.3125))
                    upperIndex = minAmpTimeIDX + (int(upperTimeCut / 0.3125))
                    timeCutWF = rawWF[lowerIndex:upperIndex]
                    runSum = -np.sum(timeCutWF)
                elif cutMode == 'timeAndNegativeAmp':
                    minAmpTimeIDX = np.argmin(rawWF)
                    lowerIndex = minAmpTimeIDX - (int(lowerTimeCut / 0.3125))
                    upperIndex = minAmpTimeIDX + (int(upperTimeCut / 0.3125))
                    timeCutWF = rawWF[lowerIndex:upperIndex]
                    negAmpMask = timeCutWF < 0
                    runSum = -np.sum(timeCutWF[negAmpMask])

                runSum = runSum + np.random.normal(runSum, 0.0005, 1)[0]
                if channel_.channel not in sumPESRAW:
                    sumPESRAW[channel_.channel] = []
                sumPESRAW[channel_.channel].append(runSum)

                if runSum > maxRunSumRAW:
                    maxRunSumRAW = runSum
                if runSum < minRunSumRAW:
                    minRunSumRAW = runSum

        binsRAW = np.linspace(minRunSumRAW, maxRunSumRAW, binCount)

        for key, val in sumPESRAW.items():
            n, binsRAW_ = np.histogram(val, bins=binsRAW)
            binsRAW_ = binsRAW_[:-1]/np.max(binsRAW_[:-1])
            # bins = bins - np.min(bins)
            ax2.fill_between(binsRAW_, n, alpha=0.4, step='pre', color=colour)
            ax2.step(binsRAW_, n, label='Integration', lw=2, color='k')

            # Set upper x limit to where there are no more bins with a height greater than 25.
            upperXLim = binsRAW_[np.where(n > 25)[0][-1]]
            lowerXLim = binsRAW_[np.where(n > 25)[0][0]]
            ax2.set_xlim([lowerXLim, upperXLim])

        ax1.set_ylim(bottom=0)
        ax2.set_ylim(bottom=0)
        leg1 = ax1.legend(fontsize=34, handlelength=0, handletextpad=0, frameon=False, alignment='left')
        leg2 = ax2.legend(fontsize=34, handlelength=0, handletextpad=0, frameon=False, alignment='left')
        for t in leg1.texts:
            t.set_alpha(0.8)
        for t in leg2.texts:
            t.set_alpha(0.8)

        plt.setp(ax1.spines.values(), linewidth=1.7)
        plt.setp(ax2.spines.values(), linewidth=1.7)

        ax1.grid()
        ax2.grid()
        # ax1.set_axis_off()
        # ax2.set_axis_off()
        ax1.set_xticks([])
        ax2.set_xticks([])

        ax1.set_yticks([])
        ax2.set_yticks([])

        ax1.grid()
        ax2.grid()

        plt.subplots_adjust(hspace=0.03)
        # ax1.set_yscale('log')
        # ax2.set_yscale('log')
        # plt.tight_layout()
        ax2.set_xlabel("Waveform amplitude", fontsize=26)
        fig.text(0.1, 0.5, 'Count', ha='center', va='center', rotation='vertical', fontsize=26)

        # plt.savefig("SumAmps_IntegrationTimeCut.png", dpi=300, bbox_inches='tight')
        plt.savefig("SumAmps_{}_{}.png".format(cutMode, runName), dpi=300, bbox_inches='tight')
        plt.close()

        # plt.show()


if __name__ == "__main__":
    print(matplotlib.font_manager.findSystemFonts(fontpaths=None, fontext='ttf')[:10])

    mpl.rcParams['font.family'] = 'Roboto'

    # recoMoreFileName = "/Users/joshuaporter/Library/CloudStorage/OneDrive-UniversityofSussex/liquidOLab/dataSTOP_DO_NOT_WRITE_HERE/WavecatcherRuns/Runs/R193/R193PES.dat"
    # rawFileName = "/Users/joshuaporter/Library/CloudStorage/OneDrive-UniversityofSussex/liquidOLab/dataSTOP_DO_NOT_WRITE_HERE/WavecatcherRuns/Runs/R193/R193.bin"
    matplotlib.rcParams.update({'font.size': 20})
    matplotlib.use("TKAgg", force=True)
    # recoMoreFileName = "/home/joshuap/Software/JoshSoftware/RecoMore/R185PES.dat"
    # rawFileName = "/home/joshuap/Software/JoshSoftware/RecoMore/R185.bin"

    # Laser run
    # recoMoreFileName = "/home/joshuap/Downloads/R15PES.dat"
    # rawFileName = "/home/joshuap/Downloads/R15.dat"

    # recoMoreFileName = "/home/joshuap/Downloads/R22121331PES.dat"
    # rawFileName = "/home/joshuap/Downloads/R22121331.bin"

    # recoMoreFileName = "/home/joshuap/Downloads/0-78PES.dat"
    # rawFileName = "/home/joshuap/Downloads/0-78.dat"

    # Scintillator run
    recoMoreFileName = "/home/joshuap/PycharmProjects/SiPMSimulations/testOutputPES.dat"
    rawFileName = "/home/joshuap/PycharmProjects/SiPMSimulations/testOutput.dat"

    examiner = RecoMoreFitExaminer(recoMoreDataPath=recoMoreFileName, rawDataPath=rawFileName)
    examiner.plotAllEvents()
    # examiner.plotSumAmps(PEThresh=0.00)
    # examiner.plotSumAmpsRaw(PEThresh=0.00)
    examiner.recoMoreVsRawComparison(channel=0, PEThresh=0.00, cutMode='time', runName='R118')
    # examiner.recoMoreVsRawComparison(channel=0, PEThresh=0.00, cutMode='none', runName='R118')
    # examiner.recoMoreVsRawComparison(channel=0, PEThresh=0.00, cutMode='negativeAmp', runName='R118')
    examiner.recoMoreVsRawComparison(channel=0, PEThresh=0.00, cutMode='timeAndNegativeAmp', runName='R118')

    # examiner.recoMoreVsRawComparison(channel=8, PEThresh=0.00, cutMode='time', runName='R15')
    # examiner.recoMoreVsRawComparison(channel=8, PEThresh=0.00, cutMode='none', runName='R15')
    # examiner.recoMoreVsRawComparison(channel=8, PEThresh=0.00, cutMode='negativeAmp', runName='R15')
    # examiner.recoMoreVsRawComparison(channel=8, PEThresh=0.00, cutMode='timeAndNegativeAmp', runName='R15')
    # examiner.timeAmpCorrelation()
    # examiner.plotAmps()
    # examiner.plotTimes()
    # examiner.plotChiSq()
