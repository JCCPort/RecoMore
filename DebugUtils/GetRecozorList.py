import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import ROOT

# ROOT.gROOT.ProcessLine("struct Event { Int_t ID; Double_t Data; };")
# ROOT.gROOT.ProcessLine("typedef struct {float amplitude;float time;} pmt_data ;")
# ROOT.gROOT.ProcessLine("typedef struct {float amplitude;float amplitude_error;float time;float time_error;float found_amplitude;float found_time;} pe_data;")
# ROOT.gROOT.ProcessLine("typedef struct {unsigned short ch;unsigned short id;float chi2ndf;float baseline;std::vector<pe_data> pes;} sipm_data ;")
# ROOT.gROOT.ProcessLine("typedef struct {unsigned short event_id;pmt_data pmt;std::vector<sipm_data> sipm;unsigned int npe_intime;unsigned int npe_offtime;} event_data ;")


ROOT.gROOT.ProcessLine("typedef struct {float amplitude;float time;} pmt_data ;"
                       "typedef struct {float amplitude;float amplitude_error;float time;float time_error;float found_amplitude;float found_time;} pe_data;"
                       "typedef struct {unsigned short ch;unsigned short id;float chi2ndf;float baseline;std::vector<pe_data> pes;} sipm_data ;"
                       "typedef struct {unsigned short event_id;pmt_data pmt;std::vector<sipm_data> sipm;unsigned int npe_intime;unsigned int npe_offtime;} event_data ;")



if __name__ == "__main__":
    fileName: str = "/sps/dchooz/diana/IJC/LiquidO/mLiquidO-data/run2_recoz/RUN_910/RUN_910_01.root"
    inFile = ROOT.TFile.Open(fileName,"READ")
    recoTree = inFile.Get("reco_tree")

    minTimeWindow = 10
    maxTimeWindow = 50

    times = []
    sipmAmps = []

    for entryNum in range(0,recoTree.GetEntries()):
        recoTree.GetEntry(entryNum)
        sipm = getattr(recoTree, "sipm")
        pmt = getattr(recoTree, "pmt")
        pes = getattr(sipm, "pes")
        for pe in range(0,pes.GetEntries()):
            peTime = getattr(pe.time)
            if (peTime > minTimeWindow) & (peTime < maxTimeWindow):
                times.append(peTime)
                sipmAmps.append(getattr(pe.amplitude))

    timeDiffs = np.diff(np.array(times))
    plt.hist(timeDiffs)
    plt.show()


