#include <stdexcept>
#include "../include/DataStructures.h"

void WCData::addRow(const EventData &wf) {
	events_.emplace_back(wf);
}

FitParams::FitParams(unsigned int numPEs, double baseline, std::vector<double> amplitudes, std::vector<double> times) {
	numPEs_ = numPEs;
	params.emplace_back(baseline);
	if(amplitudes.size() != times.size()){
		throw std::runtime_error("Amplitude vector and time vector must be the same size");
	}
	for(unsigned int i = 0; i < amplitudes.size(); i++){
		params.emplace_back(amplitudes[i]);
		params.emplace_back(times[i]);
	}
}


FitParams::FitParams(double baseline, const std::vector<PEData>& PEs) {
	numPEs_ = PEs.size();
	params.push_back(baseline);
	for(auto& pe: PEs){
		params.push_back(pe.amplitude);
		params.push_back(pe.time);
	}
}


std::vector<double *> FitParams::makeFitterParams(){
	std::vector<double *> temp_;
	for(auto &param: params){
		temp_.push_back(&param);
	}
	return temp_;
}

std::vector<float> FitParams::makeGuesserParams(){
	std::vector<float> temp_;
	temp_.push_back((float) numPEs_);
	for(auto param: params){
		temp_.push_back(param);
	}
	return temp_;
}