#include <string>
#include <cstring>
#include <regex>
#include <boost/fusion/adapted.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/phoenix.hpp>
#include <fstream>
#include "../include/DataReading.h"

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

WCData ReadWCDataFile(const std::string &fileName) {
	// Defining regular expression searches to be used for getting event and channel numbers.
	std::regex eventNumberRegex("=== EVENT (\\d*) ===\\r");
	std::regex channelNumberRegex(R"(=== CH: (\d*) EVENTID: (\d*) FCR: (\d*) ===\r)");

	WCData readData;
	WaveformData wf;

	FILE *fp = fopen(fileName.c_str(), "r");
	if (fp == nullptr)
		exit(EXIT_FAILURE);

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
		getline(&line, &len, fp);

		std::cmatch channelMatch;
		// Loop over all the channels for a given event. Will stop being true if you've just parsed the last channel for
		//  the event.
		while (std::regex_search(&line[0], channelMatch, channelNumberRegex)) {
			wf.channel = stoi(channelMatch[1].str());
			getline(&line, &len, fp);

			std::vector<float> temp;
			temp.reserve(10000);

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


/**
 *
 * @param array Array of length equal to the number of channels being read in, where each entry is an array of the ideal
 * waveform for the corresponding channel.
 * @param ch The channel to fill in.
 * @param interpFactor Number of points in interpolated waveform divided by number of points in original waveform.
 */
std::vector<double> readIdealWFs(unsigned int ch, int interpFactor, const std::string& idealWFDir, unsigned int expectedSize){
	std::string idealWFPath = idealWFDir + "ch" + std::to_string(ch) + ".txt";
	std::ifstream idealWFFile(idealWFPath, std::ifstream::in);

	if (!idealWFFile.is_open()) {
		throw std::runtime_error(idealWFPath + " not found.");
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

		for (int step = 1; step < interpFactor; ++step)
			waveform.emplace_back(prevAmp + double(step) * delta_v);

		waveform.emplace_back(idealWFAmp);
		prevAmp = idealWFAmp;
	}

	if (waveform.size() != expectedSize){
		throw std::runtime_error("Unexpected number of samples in " + idealWFPath);
	}

	return waveform;
}