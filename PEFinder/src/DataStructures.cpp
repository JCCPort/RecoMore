#include <stdexcept>
#include "../include/DataStructures.h"
#include "../include/Utils.h"

void DigitiserRun::addEvent(const DigitiserEvent &wf) {
    events.emplace_back(wf);
}

DigitiserEvent DigitiserRun::getEvent(const unsigned int eventNumber) {
    for (auto &event: events) {
        if (event.ID == eventNumber) {
            return event;
        }
    }
    throw std::runtime_error("Event " + std::to_string(eventNumber) + " not found.");
}

DigitiserChannel DigitiserRun::getEventChannel(const unsigned int eventNumber, const unsigned int channelNumber) {
    for (auto &event: events) {
        if (event.ID == eventNumber) {
            for (auto &channelWF: event.channels) {
                if (channelWF.ID == channelNumber) {
                    return channelWF;
                }
            }
        }
    }
    throw std::runtime_error(
            "Event " + std::to_string(eventNumber) + ", Channel " + std::to_string(channelNumber) + " not found.");
}


void FitRun::addEvent(const FitEvent &eventFit) {
    events.emplace_back(eventFit);
}

FitEvent FitRun::getEvent(const unsigned int eventNumber) {
    for (auto &event: events) {
        if (event.ID == eventNumber) {
            return event;
        }
    }
    throw std::runtime_error("Event " + std::to_string(eventNumber) + " not found.");
}

FitChannel FitRun::getEventChannel(const unsigned int eventNumber, const unsigned int channelNumber) {
    for (auto &event: events) {
        if (event.ID == eventNumber) {
            for (auto &channelWF: event.channels) {
                if (channelWF.ID == channelNumber) {
                    return channelWF;
                }
            }
        }
    }
    throw std::runtime_error(
            "Event " + std::to_string(eventNumber) + ", Channel " + std::to_string(channelNumber) + " not found.");
}

void FitRun::setEvents(const std::vector<FitEvent> &setData) {
    events = setData;
}

std::vector<unsigned int> FitRun::getEventIDs() {
    std::vector<unsigned int> runIDs;
    runIDs.reserve(events.size());
    for (const auto &event: events) {
        runIDs.push_back(event.ID);
    }
    return runIDs;
}


// TODO(josh): Implement this at some point to reduce the number of different ways the parameters are stored and moved.
FitParams::FitParams(const unsigned int numPEs_, const double baseline_,
                                        const std::vector<double>& amplitudes,
                                        const std::vector<double>& times){
    if (amplitudes.size() != times.size()) {
        throw std::runtime_error("Amplitude vector and time vector must be the same size");
    }
    // In this constructor, numPEs is expected to equal amplitudes.size()
    numPEs = numPEs_;
    baseline = baseline_;
    // Create a Photoelectron for each amplitude/time pair.
    for (unsigned int i = 0; i < amplitudes.size(); i++) {
        Photoelectron pe{static_cast<float>(amplitudes[i]), 0.f,
                            static_cast<float>(times[i]), 0.f,
                            0.f, 0.f};
        PEs.push_back(pe);
    }
}

FitParams::FitParams(const double baseline_, const std::vector<Photoelectron> &PEs_){
    numPEs = PEs.size();
    baseline = baseline_;
    PEs = PEs_;
}

/**
 * This method fills the solverParams vector with pointers to the parameters that the solver will adjust.
 * @param solverParams
 * @param times
 * @param amplitudes
 * @param baseline_
 */
void FitParams::makeSolverParams(std::vector<double*>* solverParams, std::vector<double>* times, std::vector<double>* amplitudes, double* baseline_) {
    *baseline_ = baseline;
    for (const auto &pe: PEs) {
        amplitudes->push_back(pe.amplitude);
        times->push_back(pe.time);
    }

    // Reserve space for baseline plus two pointers per photoelectron.
    solverParams->reserve(1 + 2 * PEs.size());
    // First parameter: baseline
    solverParams->push_back(baseline_);
    // For each photoelectron, add pointers to amplitude and time.
    for (unsigned int i = 0; i < PEs.size(); i++) {
        solverParams->push_back(&amplitudes->at(i));
        solverParams->push_back(&times->at(i));
    }
}

int FitParams::getNumParams() const {
    // Here the total number of parameters is:
    // 1 (baseline) + 2 per photoelectron + 1 (for numPEs_) = 2 * numPEs_ + 2.
    return static_cast<int>(2 * PEs.size() + 2);
}

void FitParams::sortPEsByTime() {
    std::ranges::sort(PEs, comparePETime());
}
