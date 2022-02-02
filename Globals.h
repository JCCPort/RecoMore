#ifndef RECOMORE_GLOBALS_H
#define RECOMORE_GLOBALS_H

#include <vector>

extern unsigned int winS; // Number of samples in the trigger window, x length of event data.

extern unsigned int numC; // Number of channels in readout.

extern std::vector<unsigned int> skipChannels; // Channels to skip (main use is if a channel is a trigger).

extern const int baselineNSamples; // Number of bins used for baseline computation (input value for the fit).

// If computed baseline is outside a given range, 0.0 will be used
// happening (rarely) with early pulses and with noisy waveform).
extern const float baselineRange; // 0+-X mV allowed

extern const float PEThreshold; // Amplitude threshold for PE finder algorithm (mV).

extern const float mv2pe; // PE amplitude calibration (mV/PE).

extern const float inTimeWindowTMin; // Time window for PEs considered in-time (ns).
extern const float inTimeWindowTMax; // Time window for PEs considered in-time (ns).

extern const int pdfNSamples; // Expected number of rows in ideal waveform file.

extern const double pdfSamplingRate; // Sampling rate used to create ideal waveforms.

extern const int pdfT0Sample;
extern const float PEFinderTimeOffset;
extern const float pdfResidualRMS;

extern double meanReducedChisq;

extern const double samplingRate2Inv; // Effectively the data sampling rate... should probably just have that
extern const double pdfT0SampleConv; // Conversion of sampling rate to a double

extern const double WFSigThresh; // Threshold for considering either PEs or after-pulses (V).
extern const int maxPEs;

extern float ampDiff;
extern float timeDiff;
extern float baselineDiff;

extern int sysProcPECount; // Count of number of PEs processed during entire application run. Used to calculate average error in parameter estimate.
extern int sysProcWFCount; // Count of number of waveforms processed during entire application run. Used for average reduced chisq.



#endif //RECOMORE_GLOBALS_H
