#include "Globals.h"

unsigned int winS = 1024;
unsigned int numC = 64;
std::vector<unsigned int> skipChannels{32, 36, 40, 44, 48, 52, 56, 60};
const int baselineNSamples = 32;
const float baselineRange = 1.5;
const float PEThreshold = -7.5;
const float mv2pe = 25.0;
const float inTimeWindowTMin = 70.0; // ns
const float inTimeWindowTMax = 170.0; // ns
const int pdfNSamples = 105601;
const double pdfSamplingRate = 0.003125;
const int pdfT0Sample = 3201;
const float PEFinderTimeOffset = 0.015935;
const float pdfResidualRMS = 0.827;
double meanReducedChisq = 0;
std::vector<double> reducedChiSqs {};
const double samplingRate2Inv = 1 / (0.01 * pdfSamplingRate);
const double pdfT0SampleConv = (double) pdfT0Sample;
const double WFSigThresh = 0.007;
const int maxPEs = 100;
float ampDiff = 0;
float timeDiff = 0;
float baselineDiff = 0;
int sysProcPECount = 0;
int sysProcWFCount = 0;

bool saveWaveforms = true;
unsigned int waveformCount = 0;