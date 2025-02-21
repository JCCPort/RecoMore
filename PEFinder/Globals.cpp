#include "Globals.h"

std::vector<int>    skipChannels{32, 36, 40, 44, 48, 52, 56, 60, 15};
unsigned int        pdfExternalInterpFactor = 10;
unsigned int        pdfInternalInterpFactor = 10;
unsigned int        totalInterpFactor  = pdfExternalInterpFactor * pdfInternalInterpFactor;
unsigned int        pdfNSamples        = (10560 * pdfInternalInterpFactor) + 1; // This is a count, not an index position
float               pdfSamplingRate    = 0.3125f / static_cast<float>(totalInterpFactor); // 0.3125 is true sampling rate
float               trueSamplingRate   = 0.3125f;
int                 pdfT0Sample        = static_cast<int>(320 * pdfInternalInterpFactor) + 1;
float               pdfResidualRMS     = 0.827/1000;
double              meanReducedChiSq   = 0;
std::vector<double> reducedChiSqs{};
float               samplingRate2Inv   = 1.0f / (pdfSamplingRate); //TODO(Josh): Change to double? Loss of precision isn't too significant
float               pdfT0SampleConv    = static_cast<float>(pdfT0Sample);
double              WFSigThresh        = 0.0075;
int                 maxPEs             = 100;
double              ampDiff            = 0;
double              timeDiff           = 0;
double              baselineDiff       = 0;
int                 sysProcPECount     = 0;
int                 sysProcWFCount     = 0;

bool         saveWaveforms = false;
unsigned int waveformCount = 0;

float parameterTolerance = 1e-8;