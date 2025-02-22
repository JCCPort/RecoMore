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

float NPEPDFFunc(float X, const std::vector<float> &p, const std::vector<double> *idealWaveform);

void
fitEvent(const DigitiserEvent *event, const std::unordered_map<unsigned int, PETemplate>& PETemplates, std::shared_ptr<SyncFile> outputFile, std::mutex &lock);

bool batchFitEvents(const std::vector<DigitiserEvent> &events, std::atomic<unsigned long> &count, std::mutex &lock,
                    const std::unordered_map<unsigned int, PETemplate>& PETemplates, const std::shared_ptr<SyncFile> &file);

void updateGuessCorrector(const std::vector<double>& amps, const std::vector<double>& times,
                                 const std::vector<double>& initialAmps, const std::vector<double>& initialTimes,
                                 double baseline, double initBaseline, const std::vector<Photoelectron>& pesFound);

bool getNextPEGuess(DigitiserChannel *residualWF, Photoelectron *guessPE, const double baseline, std::vector<Photoelectron> pesFound_, const DigitiserChannel& channel_, const std::vector<double> *idealWF_, ceres::CubicInterpolator<ceres::Grid1D<float>>* PDFInterpolator, std
                           ::vector<float> xValues, const PETemplate* ChPETemplate);

void amplitudeCorrection(std::vector<Photoelectron> *pesFound, std::vector<float> *params, const std::vector<float>& waveform, const std::vector<double> *);

#endif //RECOMORE_PEFIT_H
