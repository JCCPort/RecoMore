#include "Globals.h"

std::vector<int>    skipChannels{32, 36, 40, 44, 48, 52, 56, 60};
const int           pdfNSamples        = 105601;
const float         pdfSamplingRate    = 0.003125;
const int           pdfT0Sample        = 3201;
const float         PEFinderTimeOffset = 0.015935;
const float         pdfResidualRMS     = 0.827;
double              meanReducedChiSq   = 0;
std::vector<double> reducedChiSqs{};
const float         samplingRate2Inv   = 1.0f / (0.01f * pdfSamplingRate);
const float         pdfT0SampleConv    = pdfT0Sample;
const double        WFSigThresh        = 0.007;
const int           maxPEs             = 100;
float               ampDiff            = 0;
float               timeDiff           = 0;
float               baselineDiff       = 0;
int                 sysProcPECount     = 0;
int                 sysProcWFCount     = 0;

bool         saveWaveforms = true;
unsigned int waveformCount = 0;