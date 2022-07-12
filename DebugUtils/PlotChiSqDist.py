import matplotlib.pyplot as plt
import numpy as np

if __name__ == "__main__":
    dat = np.genfromtxt("../cmake-build-release/RecoMoreRedChisqs.csv", delimiter="\n")
    plt.hist(dat, bins=1000)
    plt.show()
