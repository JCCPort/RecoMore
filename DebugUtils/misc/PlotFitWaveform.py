import matplotlib.pyplot as plt
import pandas as pd

if __name__ == '__main__':
    raw = pd.read_csv("../cmake-build-debug/rawWaveform.csv", names=['Val'])
    fit = pd.read_csv("../cmake-build-debug/fit.csv", names=['Val'])
    res = pd.read_csv("../cmake-build-debug/residual.csv", names=['Val'])
    fullFit = pd.read_csv("../cmake-build-debug/fullFit.csv", names=['Val'])

    plt.plot(raw['Val'], label="raw waveform")
    plt.plot(fit['Val'], label="Fit")
    plt.plot(res['Val'], label="res")
    plt.plot(fullFit['Val'], label="fullFit", color='k')
    plt.legend()
    plt.show()
