#include "../include/DataWriting.h"

// TODO(josh): the call to _sf->write must be done at the end of the event (more precisely the end of each unit of data
//  processed by a single call to the fitting function).
std::string Writer::writeFitPE(PEData PEVar) {
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
	return writeString;
}

std::string Writer::writeWaveformInfo(const ChannelFitData& wfDat) {
	std::string writeString;
	writeString += "Ch=";
	writeString += std::to_string(wfDat.ch);
	writeString += ',';
	writeString += "RedChisq=";
	writeString += std::to_string(wfDat.chi2ndf);
	writeString += ',';
	writeString += "Baseline=";
	writeString += std::to_string(wfDat.baseline);
	writeString += '\n';
	for (auto &pe: wfDat.pes) {
		writeString += writeFitPE(pe);
	}
	writeString += "ENDPES\n";
	return writeString;
}


void Writer::writeEventInfo(const EventFitData& evData){
    std::string writeString;
    writeString += "EVENT=";
    writeString += std::to_string(evData.eventID);
    writeString += '\n';
	for(auto& chDat: evData.sipm){
		writeString += writeWaveformInfo(chDat);
	}
	_sf->write(writeString);
}