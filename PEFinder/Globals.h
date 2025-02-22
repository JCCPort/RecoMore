#ifndef RECOMORE_GLOBALS_H
#define RECOMORE_GLOBALS_H

#include <vector>

extern std::vector<int> skipChannels; // Channels to skip (main use is if a channel is a trigger).

extern unsigned int templateExternalInterpFactor; // How much oversampling was done when making the ideal waveform templates.
extern unsigned int templateInternalInterpFactor; // How much more oversampling is done in RecoMore when the templates are read in.

extern float trueSamplingRate; // True sampling rate of the data.

extern float templateResidualRMS; // Average RMS in amplitude for a point on the waveform. Used as uncertainty in y value for chi-sq calculation.

extern double meanReducedChiSq; // Average reduced chi-sq, printed at end of run to quickly check things are working/a change improved things.
extern std::vector<double> reducedChiSqs; // Used for making CSV of reduced chi-sq for debugging. DEPRECATED.

extern double WFSigThresh; // Threshold for considering either PEs or after-pulses (V).
extern int maxPEs;  // Maximum number of PEs RecoMore will find before stopping looking and doing the fit.

extern double ampDiff; // Running average of the difference between initial amplitude guess and amplitude fit value.
extern double timeDiff; // Running average of the difference between initial time guess and time fit value.
extern double baselineDiff; // Running average of the difference between initial baseline guess and baseline fit value.

extern int sysProcPECount; // Count of number of PEs processed during entire application run. Used to calculate average error in parameter estimate.
extern int sysProcWFCount; // Count of number of waveforms processed during entire application run. Used for average reduced chiSq.

extern bool saveWaveforms; // Save debug CSVs or not. This is DEPRECATED and will be removed in the future.
extern unsigned int waveformCount; // Count of waveform number used for outputting debug CSVs.

extern float parameterTolerance; // Tolerance for the fit.

#endif //RECOMORE_GLOBALS_H
