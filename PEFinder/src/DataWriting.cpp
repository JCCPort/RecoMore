#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>
#include "../include/DataWriting.h"

/**
 *
 * @param PEVar
 * @return
 */
std::string Writer::writeFitPE(Photoelectron PEVar) {
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
std::string Writer::writeWaveformInfo(const FitChannel &wfDat) {
	std::string writeString;
	writeString += "Ch=";
	writeString += std::to_string(wfDat.channel);
	writeString += ", ";
	writeString += "RedChiSq=";
	writeString += std::to_string(wfDat.redChiSq);
	writeString += ", ";
	writeString += "Baseline=";
	writeString += std::to_string(wfDat.baseline);
	writeString += '\n';
	for (auto &pe: wfDat.PEs) {
		writeString += writeFitPE(pe);
	}
	writeString += "\n";
	return writeString;
}

/**
 *
 * @param evData
 */
void Writer::writeEventInfo(const FitEvent &evData) {
	if(writeMode_ == text){
		std::string writeString;
		writeString += "EVENT=";
		writeString += std::to_string(evData.eventID);
		writeString += ", ";
		writeString += "DATE=";
		writeString += evData.date;
		writeString += ", ";
		writeString += "TDCCorrTime=";
		std::vector<std::string> timeSplitString;
		std::string tempString = evData.correctedTime.substr(10, 11);
		boost::split(timeSplitString, tempString, boost::is_any_of("."));
		writeString +=
				evData.correctedTime.substr(0, 8) + "." + timeSplitString[0] + timeSplitString[1] + timeSplitString[2] + "s";
		writeString += '\n';
		for (auto &chDat: evData.channels) {
			writeString += writeWaveformInfo(chDat);
		}
		writeString += "\n\n";
		_sf->write(writeString);
	} else if (writeMode_ == binary){
		_sf->binaryWrite(evData);
	}
}