#include "Globals.h"

std::vector<int>    skipChannels{32, 36, 40, 44, 48, 52, 56, 60, 15};
// TODO(Josh): Clarify the pdf___ numbers, why they have those values
unsigned int pdfExternalInterpFactor = 10;
unsigned int pdfInternalInterpFactor = 10;
unsigned int totalInterpFactor = pdfExternalInterpFactor * pdfInternalInterpFactor;
int                 pdfNSamples        = 105601; // This is a count, not an index position
float               pdfSamplingRate    = 0.3125f / (float)totalInterpFactor; // 0.3125 is true sampling rate
float               trueSamplingRate   = 0.3125f;
int                 pdfT0Sample        = 3201; // TODO(Josh): Now this IS an index position?
float               pdfResidualRMS     = 0.827/1000;
double              meanReducedChiSq   = 0;
std::vector<double> reducedChiSqs{};
float               samplingRate2Inv   = 1.0f / (pdfSamplingRate); //TODO(Josh): Change to double? Loss of precision isn't too significant
float               pdfT0SampleConv    = pdfT0Sample;
double              WFSigThresh        = 0.0025;
int                 maxPEs             = 100;
double              ampDiff            = 0;
double              timeDiff           = 0;
double              baselineDiff       = 0;
int                 sysProcPECount     = 0;
int                 sysProcWFCount     = 0;

bool         saveWaveforms = false;
unsigned int waveformCount = 0;