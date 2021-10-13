#ifndef RECOMORE_PEFIT_H
#define RECOMORE_PEFIT_H

#include <iostream>
#include <atomic>
#include <mutex>
#include <thread>
#include "../include/DataStructures.h"

double npe_pdf_func(const double *x, const std::vector<double>& p, std::vector<float> idealWaveform);

void fitPE(const EventData& event, const std::shared_ptr<std::vector<EventFitData>>& PEList, const std::vector<std::vector<float>>&);

#endif //RECOMORE_PEFIT_H
