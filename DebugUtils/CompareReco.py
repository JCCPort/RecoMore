import subprocess

import matplotlib.pyplot as plt
import sys

import numpy as np

from RecoMoreReader import readRecoMore

if __name__ == "__main__":
    # recoMoreFileName = sys.argv[1]
    # recozorFileName = sys.argv[2]

    recoMoreFileName = "testData/RUN_940PES.dat"
    recozorFileName = "testData/RUN_940_RECOZ.root"

    with open(recoMoreFileName, 'r') as file_object:
        RMPEs = readRecoMore(file_object)

    # return_code = subprocess.call("root 'GetRecozorList.cc(\"{}\")'".format(recozorFileName), shell=True)

    recozorAmps = np.fromfile("amp.dat", dtype=[("amplitude", np.double)])
    recozorTimes = np.fromfile("time.dat", dtype=[("time", np.double)])

    recozorAmps = [amp[0] for amp in recozorAmps]
    recozorTimes = [time[0] for time in recozorTimes]

    recomoreAmps = []
    recomoreTimes = []
    for event in RMPEs:
        for PE in event.PEData:
            recomoreAmps.append(PE.amplitude*1000)
            recomoreTimes.append(PE.time)

    bins = np.linspace(-25, 120, 100)
    plt.hist(recozorAmps, label="Recozor", density=True, histtype='step', bins=bins)
    plt.hist(recomoreAmps, label="RecoMore", density=True, histtype='step', bins=bins)
    plt.xlabel("Reco amplitudes")
    plt.legend()
    # plt.savefig("compareAmps.png")
    plt.show()

    plt.hist(recozorTimes, label="Recozor", density=True, histtype='step')
    plt.hist(recomoreTimes, label="RecoMore", density=True, histtype='step')
    plt.xlabel("Reco times")
    plt.legend()
    # plt.savefig("compareTimes.png")
    plt.show()
            


