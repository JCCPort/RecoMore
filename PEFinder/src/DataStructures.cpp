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
FitParams::FitParams(const unsigned int numPEs, const double baseline,
                                        const std::vector<double>& amplitudes,
                                        const std::vector<double>& times)
{
    if (amplitudes.size() != times.size()) {
        throw std::runtime_error("Amplitude vector and time vector must be the same size");
    }
    // In this constructor, numPEs is expected to equal amplitudes.size()
    numPEs_ = numPEs;
    baseline_ = baseline;
    // Create a Photoelectron for each amplitude/time pair.
    for (unsigned int i = 0; i < amplitudes.size(); i++) {
        Photoelectron pe{static_cast<float>(amplitudes[i]), 0.f,
                            static_cast<float>(times[i]), 0.f,
                            0.f, 0.f};
        PEParams_.push_back(pe);
    }
}

FitParams::FitParams(const double baseline, const std::vector<Photoelectron> &PEs)
{
    numPEs_ = PEs.size();
    baseline_ = baseline;
    PEParams_ = PEs;
}

void FitParams::makeSolverParams(std::vector<double*>* solverParams, std::vector<double>* times, std::vector<double>* amplitudes, double* baseline) {
    *baseline = baseline_;
    for (const auto &pe: PEParams_) {
        amplitudes->push_back(pe.amplitude);
        times->push_back(pe.time);
    }


    // Reserve space for baseline plus two pointers per photoelectron.
    solverParams->reserve(1 + 2 * PEParams_.size());
    // First parameter: baseline
    solverParams->push_back(&baseline_);
    // For each photoelectron, add pointers to amplitude and time.
    for (unsigned int i = 0; i < PEParams_.size(); i++) {
        solverParams->push_back(&amplitudes->at(i));
        solverParams->push_back(&times->at(i));
    }
}

std::vector<float> FitParams::makeGuesserParams() {
    // Create a vector with the first element as the number of photoelectrons,
    // followed by the baseline and then each photoelectron's amplitude and time.
    std::vector<float> temp_;
    temp_.push_back(static_cast<float>(numPEs_));
    temp_.push_back(static_cast<float>(baseline_));
    for (const auto &pe: PEParams_) {
        temp_.push_back(pe.amplitude);
        temp_.push_back(pe.time);
    }
    return temp_;
}

int FitParams::getNumParams() const {
    // Here the total number of parameters is:
    // 1 (baseline) + 2 per photoelectron + 1 (for numPEs_) = 2 * numPEs_ + 2.
    return static_cast<int>(2 * PEParams_.size() + 2);
}

void FitParams::sortPEsByTime() {
    std::ranges::sort(PEParams_, comparePETime());
}
