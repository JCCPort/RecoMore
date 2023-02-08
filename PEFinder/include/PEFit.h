#ifndef RECOMORE_PEFIT_H
#define RECOMORE_PEFIT_H

#include <iostream>
#include <atomic>
#include <mutex>
#include <thread>
#include <ceres/cubic_interpolation.h>
#include "DataStructures.h"
#include "DataWriting.h"

float NPEPDFFunc(float X, const std::vector<float> &p, const std::vector<double> *idealWaveform);

void
fitPE(const DigitiserEvent *event, const std::vector<std::vector<double>> *idealWaveforms, std::shared_ptr<SyncFile> outputFile, std::mutex &lock);

bool fitBatchPEs(const std::vector<DigitiserEvent> &events, std::atomic<unsigned long> &count, std::mutex &m,
                 const std::vector<std::vector<double>> *idealWaveforms, const std::shared_ptr<SyncFile> &file);

void updateGuessCorrector(const std::vector<double>& amps, const std::vector<double>& times,
						  const std::vector<double>& initialAmps, const std::vector<double>& initialTimes,
						  float baseline, float initialBaseline, const std::vector<Photoelectron>& pesFound);

bool getNextPEGuess(DigitiserChannel residualWF, Photoelectron *guessPE);

#endif //RECOMORE_PEFIT_H
