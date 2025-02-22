#include <stdexcept>
#include "../include/DataStructures.h"

void DigitiserRun::addEvent(const DigitiserEvent &wf) {
    events.emplace_back(wf);
}

DigitiserEvent DigitiserRun::getEvent(unsigned int eventNumber) {
    for (auto &event: events) {
        if (event.ID == eventNumber) {
            return event;
        }
    }
    throw std::runtime_error("Event " + std::to_string(eventNumber) + " not found.");
}

DigitiserChannel DigitiserRun::getEventChannel(unsigned int eventNumber, unsigned int channelNumber) {
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

FitEvent FitRun::getEvent(unsigned int eventNumber) {
    for (auto &event: events) {
        if (event.ID == eventNumber) {
            return event;
        }
    }
    throw std::runtime_error("Event " + std::to_string(eventNumber) + " not found.");
}

FitChannel FitRun::getEventChannel(unsigned int eventNumber, unsigned int channelNumber) {
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
FitParams::FitParams(const unsigned int numPEs, double baseline, const std::vector<double>& amplitudes,
                                      const std::vector<double>& times) {
    numPEs_ = numPEs;
    baseline_ = baseline;
    PEParams_.emplace_back(baseline);
    if (amplitudes.size() != times.size()) {
        throw std::runtime_error("Amplitude vector and time vector must be the same size");
    }
    for (unsigned int i = 0; i < amplitudes.size(); i++) {
        PEParams_.emplace_back(amplitudes[i]);
        PEParams_.emplace_back(times[i]);
    }
}


FitParams::FitParams(double baseline, const std::vector<Photoelectron> &PEs) {
    numPEs_ = PEs.size();
    baseline_ = baseline;
    for (const auto pe: PEs) {
        PEParams_.push_back(pe.amplitude);
        PEParams_.push_back(pe.time);
    }
}


void FitParams::makeSolverParams(std::vector<double*>* solverParams, std::vector<double>* times, std::vector<double>* amplitudes, double* baseline) {
    solverParams->reserve(PEParams_.size());
    solverParams->push_back(&baseline_);
    for (auto &param: PEParams_) {
        solverParams->push_back(&param);
    }
}

std::vector<float> FitParams::makeGuesserParams() {
    std::vector<float> temp_;
    temp_.push_back((float) numPEs_);
    for (auto param: PEParams_) {
        temp_.push_back((float) param);
    }
    return temp_;
}

int FitParams::getNumParams() const {
    return static_cast<int>(PEParams_.size() * 2) + 1;
}
