#include "Globals.h"

std::vector<int>    skipChannels{32, 36, 40, 44, 48, 52, 56, 60};
int                 pdfNSamples        = 105601;
float               pdfSamplingRate    = 0.003125;
int                 pdfT0Sample        = 3201;
float               pdfResidualRMS     = 0.827/1000;
double              meanReducedChiSq   = 0;
std::vector<double> reducedChiSqs{};
float               samplingRate2Inv   = 1.0f / (0.01f * pdfSamplingRate);
float               pdfT0SampleConv    = pdfT0Sample;
double              WFSigThresh        = 0.005;
int                 maxPEs             = 100;
double              ampDiff            = 0;
double              timeDiff           = 0;
double              baselineDiff       = 0;
int                 sysProcPECount     = 0;
int                 sysProcWFCount     = 0;

bool         saveWaveforms = false;
unsigned int waveformCount = 0;