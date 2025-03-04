#include "Globals.h"

unsigned int        templateInternalInterpFactor = 10;
float               templateResidualRMS     = 0.827/1000;
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

float parameterTolerance = 1e-6f;