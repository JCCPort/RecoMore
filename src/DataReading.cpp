#include <string>
#include <cstring>
#include <regex>
#include <boost/format.hpp>
#include <boost/fusion/adapted.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/phoenix.hpp>
#include <fstream>
#include <iomanip>
#include "../include/DataReading.h"
#include "../Utils.h"

WCData ReadWCDataFile(const std::string &fileName){
	// TODO(josh): Add error checking to if the data file is corrupted/invalid
	std::string ending = fileName.substr(fileName.length() - 4);
	WCData returnDat;
	if(ending == ".dat"){
		returnDat = ReadWCDataFileDat(fileName);
	}
	else if(ending == ".bin"){
		returnDat = ReadWCDataFileBinary(fileName);
	}
	else{
		throw std::runtime_error("Provided data file (" + fileName + ") is not one of the accepted formats .dat, .bin.");
	}
	std::cout << "Data file read" << std::endl;
	return returnDat;
}

// Parser shamelessly stolen from https://www.boost.org/doc/libs/1_68_0/libs/spirit/example/qi/num_list2.cpp
namespace qi = boost::spirit::qi;

namespace client {
	namespace qi = boost::spirit::qi;
	namespace ascii = boost::spirit::ascii;
	namespace phoenix = boost::phoenix;

	template<typename Iterator>
	bool parse_numbers(Iterator first, Iterator last, std::vector<float> &v) {
		using qi::float_;
		using qi::phrase_parse;
		using qi::_1;
		using phoenix::push_back;
		using qi::eol;
		bool r = phrase_parse(first, last,

				//  Begin grammar
				              (
				                      float_[push_back(phoenix::ref(v), _1)]
						                      >> *(' ' >> float_[push_back(phoenix::ref(v), _1)])
		                      ),
				//  End grammar

				              eol);

		if (first != last) // fail if we did not get a full match
			return false;
		return r;
	}
}


/**
 *
 * @param fileName Path to the raw data file to run RecoMore over.
 * @return Raw data events parsed into a WCData instance that contains a list of events to be split across multiple threads.
 */
WCData ReadWCDataFileDat(const std::string &fileName) {
	// Defining regular expression searches to be used for getting event and channel numbers.
	std::regex eventNumberRegex("=== EVENT (\\d*) ===\\r");
	std::regex channelNumberRegex(R"(=== CH: (\d*) EVENTID: (\d*) FCR: (\d*) ===\r)");
	std::regex timeRegex("(=== UnixTime = (\\d*\\.\\d*)"
	                     " date = (\\d*\\.\\d*\\.\\d*)"
	                     " time = (\\d*h\\.\\d*m\\.\\d*s\\.\\d*ms) == "
	                     "TDC = (\\d*) == "
	                     "TDC corrected time = (\\d*h\\d*m\\d*s,\\d*\\.\\d*\\.\\d*ns) == "
	                     "Nb of channels = (\\d*) ===\r\n)");

	WCData readData;
	WaveformData wf;

	FILE *fp = fopen(fileName.c_str(), "r");
	std::ifstream input_file(fileName.c_str(), std::ios::binary | std::ios::in);
	if (fp == nullptr) {
		throw std::runtime_error("WaveCatcher data file: " + fileName + " not found.");
	}

	char *line = nullptr;
	size_t len = 0;

	// Skipping lines of unneeded metadata.
	getline(&line, &len, fp);
	getline(&line, &len, fp);
	getline(&line, &len, fp);
	getline(&line, &len, fp);
	getline(&line, &len, fp);

	std::cmatch eventMatch;
	// This loop has should bring the variable `line` to the next line in the data file that gives the event number.
	//  when it doesn't it means that the end of the file has been reached.
	while (std::regex_search(&line[0], eventMatch, eventNumberRegex)) {
		EventData event;
		event.eventID = stoi(eventMatch[1].str());

		getline(&line, &len, fp);

		std::cmatch timeMatch;
		std::regex_search(&line[0], timeMatch, timeRegex);
		event.TDCCorrTime = timeMatch[6].str();
		event.date = timeMatch[3].str();

		getline(&line, &len, fp);

		std::cmatch channelMatch;
		// Loop over all the channels for a given event. Will stop being true if you've just parsed the last channel for
		//  the event.
		while (std::regex_search(&line[0], channelMatch, channelNumberRegex)) {
			wf.channel = stoi(channelMatch[1].str());
			getline(&line, &len, fp);

			std::vector<float> temp;
			temp.reserve(1024);
			std::string lineStr = std::string(line);
			client::parse_numbers(lineStr.begin(), lineStr.end(), temp);
			wf.waveform = temp;

			event.chData.push_back(wf);
			getline(&line, &len, fp);
		}
		readData.addRow(event);
	}
	fclose(fp);
	if (line)
		free(line);
	return readData;
}
#pragma pack(1)
struct WCBinaryEventData {
	int eventNumber;
	double UNIXTime;
	unsigned int year;
	unsigned int month;
	unsigned int day;
	unsigned int hour;
	unsigned int minute;
	unsigned int second;
	unsigned int millisecond;
	unsigned long long int TDCSAMIndex;
	int nChannelStored;
};

struct WCChannelDataNoMeasurement {
	int channel;
	int eventIDSAMIndex;
	int firstCellToPlotSAMIndex;
	short waveform[1024];
};
#pragma pack()

/**
 *
 * @param fileName Path to the raw data file to run RecoMore over.
 * @return Raw data events parsed into a WCData instance that contains a list of events to be split across multiple threads.
 */
WCData ReadWCDataFileBinary(const std::string &fileName) {
	WCData readData;
	WaveformData wf;

	std::ifstream input_file(fileName.c_str(), std::ios::binary | std::ios::in);
	if (!input_file) {
		throw std::runtime_error("WaveCatcher data file: " + fileName + " not found.");
	}

	std::string line;

	// Skipping lines of unneeded metadata.
	getline(input_file, line, '\n');
	getline(input_file, line, '\n');
	getline(input_file, line, '\n');
	getline(input_file, line, '\n');

	unsigned long long int previousEventTDC = 0;
	double previousEventTime = 0;
	unsigned long long int TDCOverflowCounter = 0;
	const double TDC2ns = 5.0E-9; // SAMLONG clock @ 200 MHz
	const double TDCMax = std::pow(2, 40);
	const float ADC2mV = (2500. / 4096.);
	WCBinaryEventData event_{};
	while (input_file.read((char *) (&event_), sizeof(event_))) {
		EventData event;
		event.eventID = event_.eventNumber;
		event.date = std::to_string(event_.year) + "." + std::to_string(event_.month) + "." + std::to_string(event_.day);

		if (event_.TDCSAMIndex < previousEventTDC)
			TDCOverflowCounter++;

		double outputTime = (double)TDCOverflowCounter * TDCMax * TDC2ns + (double)event_.TDCSAMIndex * TDC2ns;
		double outputDeltaTime = outputTime - previousEventTime;
		previousEventTime = outputTime;
		std::string subSec = to_string_with_precision(outputDeltaTime + (event_.millisecond / 1000.), 9);

		std::ostringstream ssh;
		ssh << std::setw(2) << std::setfill('0') << std::to_string(event_.hour);

		std::ostringstream ssm;
		ssm << std::setw(2) << std::setfill('0') << std::to_string(event_.minute);

		std::ostringstream sss;
		sss << std::setw(2) << std::setfill('0') << std::to_string(event_.second);

		event.TDCCorrTime = ssh.str() + "h" + ssm.str() + "m" + sss.str() + "s"
		                    + "," + subSec.substr(2, 3) + "." + subSec.substr(5, 3) + "." + subSec.substr(8, 3) + "ns";


		for(int j = 0; j < event_.nChannelStored; j++){
			WCChannelDataNoMeasurement waveform{};
			input_file.read((char *) (&waveform), sizeof(WCChannelDataNoMeasurement));
			wf.channel = waveform.channel;
			std::vector<float> temp;
			for(short i : waveform.waveform){
				temp.push_back(ADC2mV * (float)i / 10000);  // There was a factor of 10 originally in recozor, it became 10000 because we're using V not mV
			}
			wf.waveform = temp;
			event.chData.push_back(wf);
		}
		readData.addRow(event);
	}
	return readData;
}


/**
 * This function reads the ideal PDFs for each channel which is then used for finding PEs in data by overlapping the ideal PDF with
 * data PEs and subtracting them until none remain.
 * @param ch The channel to fill in.
 * @param interpFactor Number of points in interpolated waveform divided by number of points in original waveform.
 * @param idealWFDir Path to directory containing ideal PDFs for each channel.
 * @param expectedSize Check that the PDF is the expected length.
 * @return
 */
std::vector<double>
readIdealWFs(unsigned int ch, int interpFactor, const std::string &idealWFDir, unsigned int expectedSize) {
	std::string idealWFPath = idealWFDir + "ch" + std::to_string(ch) + ".txt";
	std::ifstream idealWFFile(idealWFPath, std::ifstream::in);

	if (!idealWFFile.is_open()) {
		throw std::runtime_error("Ideal PE PDF: " + idealWFPath + " not found.");
	}

	double idealWFTime;
	double idealWFAmp;
	double prevAmp;

	// Parse first line
	idealWFFile >> idealWFTime >> idealWFAmp;
	std::vector<double> waveform;
	waveform.emplace_back(idealWFAmp);
	prevAmp = idealWFAmp;

	// Parse next line
	while (idealWFFile >> idealWFTime >> idealWFAmp) {
		double delta_v = (idealWFAmp - prevAmp) / double(interpFactor);

		for (int step = 1; step < interpFactor; ++step)  // Add linearly interpolated points to ideal PDF
			waveform.emplace_back(prevAmp + double(step) * delta_v);

		waveform.emplace_back(idealWFAmp);
		prevAmp = idealWFAmp;
	}

	if (waveform.size() != expectedSize) {
		throw std::runtime_error("Unexpected number of samples in " + idealWFPath);
	}

	return waveform;
}