#include <string>
#include <cstring>
#include <regex>
#include <boost/fusion/adapted.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/phoenix.hpp>
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
	using qi::double_;
	using qi::phrase_parse;

	// Defining regular expression searches to be used for getting event and channel numbers.
	std::regex eventNumberRegex("=== EVENT (\\d*) ===\\r");
	std::regex channelNumberRegex(R"(=== CH: (\d*) EVENTID: (\d*) FCR: (\d*) ===\r)");

	WCData readData;
	Waveform wf;


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
		wf.event_ = stoi(eventMatch[1].str());
		getline(&line, &len, fp);
		getline(&line, &len, fp);

		std::cmatch channelMatch;
		// Loop over all the channels for a given event. Will stop being true if you've just parsed the last channel for
		//  the event.
		while (std::regex_search(&line[0], channelMatch, channelNumberRegex)) {
			wf.channel_ = stoi(channelMatch[1].str());
			getline(&line, &len, fp);

			std::vector<float> temp;
			temp.reserve(1000);

			std::string lineStr = std::string(line);
			client::parse_numbers(lineStr.begin(), lineStr.end(), temp);
			wf.waveform_ = temp;
			readData.addRow(wf);

			getline(&line, &len, fp);
		}
	}
	fclose(fp);
	if (line)
		free(line);
	return readData;
}