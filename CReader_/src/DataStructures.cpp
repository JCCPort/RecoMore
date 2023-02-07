#include <stdexcept>
#include "../include/DataStructures.h"

void WCData::addRow(const EventData &wf) {
	events_.emplace_back(wf);
}

EventData WCData::getEvent(int eventNumber) {
	for(auto & event : events_){
		if(event.eventID == eventNumber){
			return event;
		}
	}
	throw std::runtime_error("Event " + std::to_string(eventNumber) + " not found.");
}

ChannelData WCData::getChannelWaveform(int eventNumber, int channelNumber) {
	for(auto & event : events_){
		if(event.eventID == eventNumber){
			for(auto & channelWF : event.chData){
				if(channelWF.channel == channelNumber){
					return channelWF;
				}
			}
		}
	}
	throw std::runtime_error("Event " + std::to_string(eventNumber) + ", Channel " + std::to_string(channelNumber) + " not found.");
}


void FitData::addRow(const EventFitData &eventFit) {
	fitEvents_.emplace_back(eventFit);
}

EventFitData FitData::getEventFit(int eventNumber) {
	for(auto & event : fitEvents_){
		if(event.eventID == eventNumber){
			return event;
		}
	}
	throw std::runtime_error("Event " + std::to_string(eventNumber) + " not found.");
}

ChannelFitData FitData::getChannelFit(int eventNumber, int channelNumber) {
	for(auto & event : fitEvents_){
		if(event.eventID == eventNumber){
			for(auto & channelWF : event.SiPM){
				if(channelWF.ch == channelNumber){
					return channelWF;
				}
			}
		}
	}
	throw std::runtime_error("Event " + std::to_string(eventNumber) + ", Channel " + std::to_string(channelNumber) + " not found.");
}

void FitData::setRows(const std::vector<EventFitData> &setData) {
	fitEvents_ = setData;
}


// TODO(josh): Implement this at some point to reduce the number of different ways the parameters are stored and moved.
[[maybe_unused]] FitParams::FitParams(unsigned int numPEs, double baseline, std::vector<double> amplitudes, std::vector<double> times) {
	numPEs_ = numPEs;
	params.emplace_back(baseline);
	if (amplitudes.size() != times.size()) {
		throw std::runtime_error("Amplitude vector and time vector must be the same size");
	}
	for (unsigned int i = 0; i < amplitudes.size(); i++) {
		params.emplace_back(amplitudes[i]);
		params.emplace_back(times[i]);
	}
}


[[maybe_unused]] FitParams::FitParams(double baseline, const std::vector<PEData> &PEs) {
	numPEs_ = PEs.size();
	params.push_back(baseline);
	for (auto &pe: PEs) {
		params.push_back(pe.amplitude);
		params.push_back(pe.time);
	}
}


[[maybe_unused]] std::vector<double *> FitParams::makeFitterParams() {
	std::vector<double *> temp_;
	temp_.reserve(params.size());
	for (auto &param: params) {
		temp_.push_back(&param);
	}
	return temp_;
}

[[maybe_unused]] std::vector<float> FitParams::makeGuesserParams() {
	std::vector<float> temp_;
	temp_.push_back((float) numPEs_);
	for (auto param: params) {
		temp_.push_back((float)param);
	}
	return temp_;
}
