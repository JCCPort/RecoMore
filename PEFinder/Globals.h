#ifndef RECOMORE_GLOBALS_H
#define RECOMORE_GLOBALS_H

#include <vector>

extern std::vector<int> skipChannels; // Channels to skip (main use is if a channel is a trigger).

extern const int pdfNSamples; // Expected number of rows in ideal waveform file.

extern const float pdfSamplingRate; // Sampling rate used to create ideal waveforms.

extern const int pdfT0Sample;
extern const float PEFinderTimeOffset;
extern const float pdfResidualRMS;

extern double              meanReducedChiSq;
extern std::vector<double> reducedChiSqs;

extern const float samplingRate2Inv; // Effectively the data sampling rate... should probably just have that
extern const float pdfT0SampleConv; // Conversion of sampling rate to a double

extern const double WFSigThresh; // Threshold for considering either PEs or after-pulses (V).
extern const int maxPEs;

extern float ampDiff;
extern float timeDiff;
extern float baselineDiff;

extern int sysProcPECount; // Count of number of PEs processed during entire application run. Used to calculate average error in parameter estimate.
extern int sysProcWFCount; // Count of number of waveforms processed during entire application run. Used for average reduced chiSq.

extern bool saveWaveforms;
extern unsigned int waveformCount;

#endif //RECOMORE_GLOBALS_H
