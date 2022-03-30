#include <string>
#include <cstring>
#include <regex>
#include <boost/fusion/adapted.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/phoenix.hpp>
#include <fstream>
#include "../include/DataReading.h"
#include "../Utils.h"

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


const char* BinaryStringToText(std::string binaryString) {
	std::string text = "";
	std::stringstream sstream(binaryString);
	while (sstream.good())
	{
		std::bitset<8> bits;
		sstream >> bits;
		text += char(bits.to_ulong());
	}
	return text.c_str();
}


/**
 *
 * @param fileName Path to the raw data file to run RecoMore over.
 * @return Raw data events parsed into a WCData instance that contains a list of events to be split across multiple threads.
 */
WCData ReadWCDataFile(const std::string &fileName) {
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
struct eventData {
	int EventNumber;
	double EpochTime;
	unsigned int Year;
	unsigned int Month;
	unsigned int Day;
	unsigned int Hour;
	unsigned int Minute;
	unsigned int Second;
	unsigned int Millisecond;
	unsigned long long int TDCsamIndex;
	int nchannelstored;
};

struct waveformNoMeas {
	int channel;
	int EventIDsamIndex;
	int FirstCellToPlotsamIndex;
	short waveform[1024];
};
#pragma pack()

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

	unsigned long long int previous_event_tdc = 0;
	double previous_event_time = 0;
	unsigned long long int tdc_overflow_counter = 0;
	const double tdc2ns = 5.0E-9; // SAMLONG clock @ 200 MHz
	const double tdcmax = pow(2, 40);
	const float adc2mv = (2500. / 4096.);
	eventData event_;
	while (input_file.read((char *) (&event_), sizeof(event_))) {
		EventData event;
		event.eventID = event_.EventNumber;
		event.date = std::to_string(event_.Year) + "." + std::to_string(event_.Month) + "." + std::to_string(event_.Day);

		if (event_.TDCsamIndex < previous_event_tdc)
			tdc_overflow_counter++;

		double output_time = tdc_overflow_counter * tdcmax * tdc2ns + event_.TDCsamIndex * tdc2ns;
		double output_delta_time = output_time - previous_event_time;
		previous_event_time = output_time;

		std::string subSec = to_string_with_precision(output_delta_time+(event_.Millisecond/1000.), 9);
		event.TDCCorrTime = std::to_string(event_.Hour) + "h" + std::to_string(event_.Minute) + "m" + std::to_string(event_.Second) + "s"
				+ "," + subSec.substr(2, 3) + "." + subSec.substr(5, 3) + "." + subSec.substr(8, 3) + "ns";

		for(int j = 0; j < event_.nchannelstored; j++){
			waveformNoMeas waveform;
			input_file.read((char *) (&waveform), sizeof(waveformNoMeas));
			wf.channel = waveform.channel;
			std::vector<float> temp;
			for(short i : waveform.waveform){
				temp.push_back(adc2mv * i / 10000);
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
// TODO(josh): Fix performance of this as it seems to be slower than the WC data opening!
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