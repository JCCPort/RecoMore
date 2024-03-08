from typing import List

import numpy as np
import pandas as pd

from DebugUtils.PythonRecoMoreReader import PE


def specialRound(x, base):
    return base * np.round(float(x)/base, 0)


def NPEPDFFunc(x, PEs: List[PE], baseline, idealWF):
    value = baseline
    for PE_ in PEs:
        amp = PE_.amplitude
        time = PE_.time

        PDFBin = 3201 + int(np.floor(0.5 + ((x - time)/100) * (1 / (0.01 * 0.003125))))
        if (PDFBin >= 0) and (PDFBin < 105601):
            value += (amp * idealWF[PDFBin])
    return value


def makeWaveformArray(waveform, RMEvent):
    # pdf = pd.read_csv("../PEFinder/pdf/ch{}.txt".format(waveform.channel), sep="\t", names=['time', 'amp'])
    pdf = pd.read_csv("/home/joshuap/PycharmProjects/SiPMSimulations/OutputDataForRecoMore/pdf/ch{}.txt".format(waveform.channel), sep="\t", names=['time', 'amp'])
    # TODO(Josh): Should have this path as a variable so other SiPMs can be used
    amps = pdf['amp'].values
    times = pdf['time'].values

    interpedAmps = []
    interpedTimes = []

    prevAmp = amps[0]
    prevTime = times[0]
    for i in range(1, len(amps)):
        delta_vAmp = (amps[i] - prevAmp) / 10
        delta_vTime = (times[i] - prevTime) / 10
        for step in range(10):
            interpedAmps.append((prevAmp + step * delta_vAmp))
            interpedTimes.append((prevTime + step * delta_vTime))
        prevAmp = amps[i]
        prevTime = times[i]
    interpedAmps.append(amps[-1])
    interpedTimes.append(times[-1])

    # xS = np.array(interpedTimes)
    yS = np.array(interpedAmps)
    sumYs = np.zeros(len(waveform.waveform))

    for i in range(len(waveform.xVals)):
        x = waveform.xVals[i]
        sumYs[i] = NPEPDFFunc(x, RMEvent.pes, RMEvent.baseline, yS)

    return waveform.xVals, sumYs
