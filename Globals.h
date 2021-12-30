#ifndef RECOMORE_GLOBALS_H
#define RECOMORE_GLOBALS_H

#include <vector>

//  Number of samples in the trigger window, x length of event data.
extern unsigned int winS;

//  Number of channels in readout.
extern unsigned int numC;

// Channels to skip.
extern std::vector<unsigned int> skipChannels;

// Number of bins used for baseline computation (input value for the fit).
extern const int baselineNSamples;

// If computed baseline is outside a given range, 0.0 will be used
// happening (rarely) with early pulses and with noisy waveform).
extern const float baselineRange; // 0+-X mV allowed

// Amplitude threshold for PE finder algorithm.
extern const float PEThreshold; // mV

// PE amplitude calibration.
extern const float mv2pe; // mV/PE

// Time window for PEs considered in-time.
extern const float inTimeWindowTMin; // ns
extern const float inTimeWindowTMax; // ns

// Expected number of rows in ideal waveform file.
extern const int pdfNSamples;

// Sampling rate used to create ideal waveforms.
extern const double pdfSamplingRate;

extern const int pdfT0Sample;
extern const float PEFinderTimeOffset;
extern const float pdfResidualRMS;

extern double meanReducedChisq;

extern const double samplingRate2Inv; // Effectively the data sampling rate... should probably just have that
extern const double pdfT0SampleConv; // Conversion of sampling rate to a double

extern const double WFSigThresh; // Threshold for considering either PEs or after-pulses
extern const int maxPEs;



#endif //RECOMORE_GLOBALS_H
