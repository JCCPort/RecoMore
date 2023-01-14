import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt


def hitAmpDistribution(x, A, c, mu, Q1, sigma1, sigma2, ratio):
    n = 100
    val = 0
    fact = 1

    expMinusMu = np.exp(-mu)
    sigma1RootPi = sigma1*np.sqrt(2*np.pi)
    sigma2RootPi = sigma2*np.sqrt(2*np.pi)

    for i in range(1, int(n)):
        fact = fact * i
        poissonPart = ((mu**i) / fact) * expMinusMu
        gaussPart1 = (ratio / (i*sigma1RootPi)) * np.exp(-(((x-c) - i*Q1)**2) / (4*(i*sigma1)**2))
        gaussPart2 = ((1-ratio) / (i*sigma2RootPi)) * np.exp(-(((x-c) - i*Q1)**2) / (4*(i*sigma2)**2))
        val += poissonPart*(gaussPart1 + gaussPart2)
    return val*A


def fit(amps):
    n, bins = np.histogram(amps, bins=1000)
    bins = [(bins[i] + bins[i+1])/2 for i in range(len(bins)-1)]

    p0 = [np.max(n)/10.5, 0.54136343e-02, 9.6*2, 3.44230966e-02/2, 1.82453950e-04, 1.8453950e-04, 1]
    # p0 = [np.max(n), 0, 35, 30, np.mean(amps) / 350, 0.005, np.mean(amps) / 350, 0.005]

    try:
        pOpt, pCov = curve_fit(hitAmpDistribution, xdata=bins, ydata=n, p0=p0)
    except RuntimeError:
        pOpt = p0

    fitVals = np.linspace(np.min(bins), np.max(bins), 1000)

    print(p0)
    print(pOpt)

    plt.plot(fitVals, hitAmpDistribution(fitVals, *pOpt), label="Fit")
    plt.step(bins, n, label="Data")
    plt.xlabel("Amplitude (V)")
    plt.ylabel("Count")
    plt.legend()
    plt.tight_layout()
    plt.show()
