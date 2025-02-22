#include "Globals.h"

std::vector<int>    skipChannels{32, 36, 40, 44, 48, 52, 56, 60, 15};
unsigned int        pdfExternalInterpFactor = 10;
unsigned int        pdfInternalInterpFactor = 10;
float               trueSamplingRate   = 0.3125f;
float               pdfResidualRMS     = 0.827/1000;
double              meanReducedChiSq   = 0;
std::vector<double> reducedChiSqs{};
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