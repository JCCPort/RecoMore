#include <string>
#include <fstream>
#include <sstream>
#include <cstring>
#include <regex>
#include <c++/10/iostream>
#include "../include/DataReading.h"

bool startsWith(const char *a, const char *b) {
	if (strncmp(a, b, strlen(b)) == 0) return true;
	return false;
}


WCData ReadWCDataFile(const std::string &fileName) {
	std::regex eventNumberRegex("=== EVENT (\\d*) ===\\r");
	std::regex channelNumberRegex(R"(=== CH: (\d*) EVENTID: (\d*) FCR: (\d*) ===\r)");

	WCData readData;
	Waveform wf;

	std::ifstream file(fileName);
	if (file.is_open()) {
		std::string line;
		std::getline(file, line);
		std::getline(file, line);
		std::getline(file, line);
		std::getline(file, line);
		std::getline(file, line);

		std::smatch eventMatch;
		while (std::regex_search(line, eventMatch, eventNumberRegex)) {
			wf.event_ = stoi(eventMatch[1].str());
			std::getline(file, line);
			std::getline(file, line);

			std::smatch channelMatch;
			while (std::regex_search(line, channelMatch, channelNumberRegex)) {
				wf.channel_ = stoi(channelMatch[1].str());

				std::getline(file, line);
				std::istringstream waveformStream(line);

				std::vector<float> temp{std::istream_iterator<float>(waveformStream), std::istream_iterator<float>()};
				wf.waveform_ = temp;
				readData.addRow(wf);

				std::getline(file, line);
			}
		}
		file.close();
	}
	return readData;
}
