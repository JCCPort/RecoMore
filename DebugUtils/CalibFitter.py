import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
# from tabulate import tabulate


def hitAmpDistribution(x, A, c, mu, Q, sigma1, sigma2, ratio):
    n = 100
    val = 0
    fact = 1

    expMinusMu = np.exp(-mu)
    sigma1RootPi = sigma1*np.sqrt(2*np.pi)
    sigma2RootPi = sigma2*np.sqrt(2*np.pi)

    for i in range(1, int(n)):
        fact = fact * i
        poissonPart = ((mu**i) / fact) * expMinusMu
        gaussPart1 = (ratio / (i*sigma1RootPi)) * np.exp(-(((x-c) - i * Q) ** 2) / (4 * (i * sigma1) ** 2))
        gaussPart2 = ((1-ratio) / (i*sigma2RootPi)) * np.exp(-(((x-c) - i * Q) ** 2) / (4 * (i * sigma2) ** 2))
        val += poissonPart*(gaussPart1 + gaussPart2)
    return val*A


def fit(amps):
    n, bins = np.histogram(amps, bins=1600)
    bins = [(bins[i] + bins[i+1])/2 for i in range(len(bins)-1)]

    p0 = [np.max(n)/10.5, 0.54136343e-02, 9.6*2, 3.44230966e-02/2, 1.82453950e-04, 1.8453950e-04, 1]
    # p0 = [np.max(n), 0, 35, 30, np.mean(amps) / 350, 0.005, np.mean(amps) / 350, 0.005]

    bounds = [[1, -np.inf, 0, 0, 0, 0, -1],
              [np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, 1]]

    try:
        pOpt, pCov = curve_fit(hitAmpDistribution, xdata=bins, ydata=n, p0=p0, sigma=np.sqrt(n+1), bounds=bounds)
    except RuntimeError:
        pOpt = p0

    fitVals = np.linspace(np.min(bins), np.max(bins), 5000)

    params = [p0, bounds[0], bounds[1], pOpt]

    chiSq = [0, 0, 0, 0]
    for idx, chiSq_ in enumerate(chiSq):
        for i in range(0, len(bins)):
            chiSq[idx] += (n[i] - hitAmpDistribution(bins[i], *params[idx]))**2/(np.sqrt(n[i] + 1))**2

    chiSq = np.array(chiSq)
    reducedChiSq = chiSq / (len(bins) - 6)

    p0.append(reducedChiSq[0])
    bounds[0].append(reducedChiSq[1])
    bounds[1].append(reducedChiSq[2])
    displayPOpt = list(pOpt)
    displayPOpt.append(reducedChiSq[3])

    # TODO(josh): Add this to be displayed on the graph
    # print(tabulate([p0, bounds[0], bounds[1], displayPOpt],
    #                headers=['A', 'c', 'mu', 'Q', 'sigma1', 'sigma2', 'ratio', 'reducedChiSq'],
    #                showindex=["Initial guess", "Lower bound", "Upper bound", "Fit"],
    #                floatfmt=['.1f', '.2f', '.4f', '.2f', '.3f', '.5f', '.5f', '.2f', '.2f'],
    #                numalign='left'))

    plt.plot(fitVals, hitAmpDistribution(fitVals, *pOpt), label="Fit")
    plt.step(bins, n, label="Data")
    plt.xlabel("Amplitude (V)")
    plt.ylabel("Count")
    plt.legend()
    plt.tight_layout()
    plt.show()
