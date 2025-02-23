#ifndef RECOMORE_PEFIT_H
#define RECOMORE_PEFIT_H

#include <iostream>
#include <atomic>
#include <mutex>
#include <thread>
#include <ceres/cubic_interpolation.h>
#include "DataStructures.h"
#include "DataWriting.h"
#include "PETemplate.h"

void fitEvent(const DigitiserEvent *event, const std::unordered_map<unsigned int, PETemplate>& PETemplates, std::shared_ptr<SyncFile> outputFile, std::mutex &lock, float sampleSpacing);

bool batchFitEvents(const std::vector<DigitiserEvent> &                 events,
                    std::atomic<unsigned long> &                        count,
                    std::mutex &                                        lock,
                    const std::unordered_map<unsigned int, PETemplate>& PETemplates,
                    const std::shared_ptr<SyncFile> &                   file,
                    float                                               sampleSpacing);

void updateGuessCorrector(const std::vector<double>& amps, const std::vector<double>& times,
                                 const std::vector<double>& initialAmps, const std::vector<double>& initialTimes,
                                 const double baseline, const double initBaseline, const unsigned int numPEs);

bool getNextPEGuess(const std::vector<float>& residualWF, Photoelectron *guessPE, float sampleSpacing);

void amplitudeCorrection(FitParams fitParams, const std::vector<float>& waveform, const ceres::CubicInterpolator<ceres::Grid1D<float>>* templateInterpolator, const PETemplate* ChPETemplate, float sampleSpacing);


#endif //RECOMORE_PEFIT_H
