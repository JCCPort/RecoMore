#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>
#include "../include/DataWriting.h"

// TODO(josh): the call to _sf->write must be done at the end of the event (more precisely the end of each unit of data
//  processed by a single call to the fitting function).
/**
 *
 * @param PEVar
 * @return
 */
std::string Writer::writeFitPE(PEData PEVar) {
	std::string writeString;
	writeString += std::to_string(PEVar.amplitude);
	writeString += ',';
//	writeString += std::to_string(PEVar.amplitudeError);
//	writeString += ',';
	writeString += std::to_string(PEVar.time * 100);
//	writeString += ',';
//	writeString += std:q:to_string(PEVar.timeError);
//	writeString += ',';
//	writeString += std::to_string(PEVar.foundAmplitude);
//	writeString += ',';
//	writeString += std::to_string(PEVar.foundTime);
	writeString += '\n';
	return writeString;
}

/**
 *
 * @param wfDat
 * @return
 */
std::string Writer::writeWaveformInfo(const ChannelFitData &wfDat) {
	std::string writeString;
	writeString += "Ch=";
	writeString += std::to_string(wfDat.ch);
	writeString += ", ";
	writeString += "RedChisq=";
	writeString += std::to_string(wfDat.chi2ndf);
	writeString += ", ";
	writeString += "Baseline=";
	writeString += std::to_string(wfDat.baseline);
	writeString += '\n';
	for (auto &pe: wfDat.pes) {
		writeString += writeFitPE(pe);
	}
	writeString += "\n";
	return writeString;
}

/**
 *
 * @param evData
 */
void Writer::writeEventInfo(const EventFitData &evData) {
	std::string writeString;
	writeString += "EVENT=";
	writeString += std::to_string(evData.eventID);
	writeString += ", ";
	writeString += "DATE=";
	writeString += evData.date;
	writeString += ", ";
	writeString += "TDCCorrTime=";
	std::vector<std::string> timeSplitString;
	boost::split(timeSplitString, evData.TDCCorrTime.substr(10, 11), boost::is_any_of("."));
	writeString += evData.TDCCorrTime.substr(0, 8) + "." + timeSplitString[0] + timeSplitString[1] + timeSplitString[2] + "s";
	writeString += '\n';
	for (auto &chDat: evData.sipm) {
		writeString += writeWaveformInfo(chDat);
	}
	writeString += "\n\n";
	_sf->write(writeString);
}