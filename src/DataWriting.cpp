#include "../include/DataWriting.h"

void Writer::writeFitPE(PEData PEVar) {
	std::string writeString;
	writeString += std::to_string(PEVar.amplitude);
	writeString += ',';
	writeString += std::to_string(PEVar.amplitudeError);
	writeString += ',';
	writeString += std::to_string(PEVar.time);
	writeString += ',';
	writeString += std::to_string(PEVar.timeError);
	writeString += ',';
	writeString += std::to_string(PEVar.foundAmplitude);
	writeString += ',';
	writeString += std::to_string(PEVar.foundTime);
	writeString += '\n';
	_sf->write(writeString);
}

void Writer::writeWaveformInfo(const waveformData& wfDat) {
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
}


void Writer::writeEventInfo(const eventData& evData){
    std::string writeString;
    writeString += "EVENT=\n";
    writeString += std::to_string(evData.eventID);
    writeString += '\n';
    _sf->write(writeString);
}