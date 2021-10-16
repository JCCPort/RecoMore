#include "../include/DataWriting.h"

// TODO(josh): the call to _sf->write must be done at the end of the event (more precisely the end of each unit of data
//  processed by a single call to the fitting function).
void Writer::writeFitPE(PEData PEVar) {
	std::string writeString;
	writeString += std::to_string(PEVar.amplitude);
	writeString += ',';
	writeString += std::to_string(PEVar.amplitudeError);
	writeString += ',';
	writeString += std::to_string(PEVar.time * 100);
	writeString += ',';
	writeString += std::to_string(PEVar.timeError);
	writeString += ',';
	writeString += std::to_string(PEVar.foundAmplitude);
	writeString += ',';
	writeString += std::to_string(PEVar.foundTime);
	writeString += '\n';
	_sf->write(writeString);
}

void Writer::writeWaveformInfo(const ChannelFitData& wfDat) {
	std::string writeString;
	writeString += "WAVEFORM META INFO\n";
	writeString += "Ch=";
	writeString += std::to_string(wfDat.ch);
	writeString += ',';
	writeString += "RedChisq=";
	writeString += std::to_string(wfDat.chi2ndf);
	writeString += ',';
	writeString += "Baseline=";
	writeString += std::to_string(wfDat.baseline);
	writeString += '\n';
	_sf->write(writeString);
	for(auto& pe : wfDat.pes){
		writeFitPE(pe);
	}
}


void Writer::writeEventInfo(const EventFitData& evData){
    std::string writeString;
    writeString += "EVENT=";
    writeString += std::to_string(evData.eventID);
    writeString += '\n';
    _sf->write(writeString);
	for(auto& chDat: evData.sipm){
		writeWaveformInfo(chDat);
	}
}