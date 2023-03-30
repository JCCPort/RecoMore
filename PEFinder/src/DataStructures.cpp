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
[[maybe_unused]] FitParams::FitParams(unsigned int numPEs, double baseline, std::vector<double> amplitudes,
                                      std::vector<double> times) {
    numPEs_ = numPEs;
    PEParams_.emplace_back(baseline);
    if (amplitudes.size() != times.size()) {
        throw std::runtime_error("Amplitude vector and time vector must be the same size");
    }
    for (unsigned int i = 0; i < amplitudes.size(); i++) {
        PEParams_.emplace_back(amplitudes[i]);
        PEParams_.emplace_back(times[i]);
    }
}


[[maybe_unused]] FitParams::FitParams(double baseline, const std::vector<Photoelectron> &PEs) {
    numPEs_ = PEs.size();
    baseline_ = baseline;
    for (const auto pe: PEs) {
        PEParams_.push_back(pe.amplitude);
        PEParams_.push_back(pe.time);
    }
}


[[maybe_unused]] void FitParams::makeFitterParams(std::vector<double *> temp_) {
    temp_.reserve(PEParams_.size());
    temp_.push_back(&baseline_);
    for (auto &param: PEParams_) {
        temp_.push_back(&param);
    }
}

[[maybe_unused]] std::vector<float> FitParams::makeGuesserParams() {
    std::vector<float> temp_;
    temp_.push_back((float) numPEs_);
    for (auto param: PEParams_) {
        temp_.push_back((float) param);
    }
    return temp_;
}

int FitParams::getNumParams() const {
    return (int) PEParams_.size() + 2;
}
