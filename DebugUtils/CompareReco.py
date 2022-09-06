import subprocess

import matplotlib.pyplot as plt
import sys

import numpy as np

from RecoMoreReader import readRecoMore

if __name__ == "__main__":
    recoMoreFileName = sys.argv[1]
    recozorFileName = sys.argv[2]

    with open(recoMoreFileName, 'r') as file_object:
        RMPEs = readRecoMore(file_object)

    return_code = subprocess.call("root 'GetRecozorList.cc(\"{}\")'".format(recozorFileName), shell=True)

    recozorAmps = np.fromfile("amp.dat", dtype=[("amplitude", np.longlong)])
    recozorTimes = np.fromfile("time.dat", dtype=[("time", np.longlong)])

    recomoreAmps = []
    recomoreTimes = []
    for event in RMPEs:
        for PE in event.PEData:
            recomoreAmps.append(PE.amplitude)
            recomoreTimes.append(PE.time)
            
    plt.hist(recozorAmps, label="Recozor")
    plt.hist(recomoreAmps, label="RecoMore")
    plt.xlabel("Reco amplitudes")
    plt.legend()
    plt.savefig("compareAmps.png")

    plt.hist(recozorTimes, label="Recozor")
    plt.hist(recomoreTimes, label="RecoMore")
    plt.xlabel("Reco times")
    plt.legend()
    plt.savefig("compareTimes.png")
            


