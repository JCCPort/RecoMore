#ifndef RECOMORE_PEFIT_H
#define RECOMORE_PEFIT_H

#include <iostream>
#include <atomic>
#include <mutex>
#include <thread>
#include "../include/DataStructures.h"
#include "DataWriting.h"

float npe_pdf_func(float X, const std::vector<float> &p, std::vector<double> *idealWaveform);

void
fitPE(const EventData *event, const std::vector<std::vector<double>> *idealWaveforms, std::shared_ptr<SyncFile> file);

#endif //RECOMORE_PEFIT_H
