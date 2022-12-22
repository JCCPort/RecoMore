#ifndef RECOMORE_PEFIT_H
#define RECOMORE_PEFIT_H

#include <iostream>
#include <atomic>
#include <mutex>
#include <thread>
#include <ceres/cubic_interpolation.h>
#include "../include/DataStructures.h"
#include "DataWriting.h"

float NPEPDFFunc(float X, const std::vector<float> &p, const std::vector<double> *idealWaveform);

void
fitPE(const EventData *event, const std::vector<std::vector<double>> *idealWaveforms, std::shared_ptr<SyncFile> file, std::mutex &m);

bool fitBatchPEs(const std::vector<EventData> &events, std::atomic<unsigned long> &count, std::mutex &m,
                 const std::vector<std::vector<double>> *idealWaveforms, const std::shared_ptr<SyncFile> &file);

#endif //RECOMORE_PEFIT_H
