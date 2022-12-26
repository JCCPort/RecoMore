#ifndef RECOMORE_UTILS_H
#define RECOMORE_UTILS_H

#include <iostream>
#include <atomic>
#include <thread>
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>
#include <iomanip>
#include "../Globals.h"


/**
 *
 * @param count
 * @param m
 * @param dataLength
 */
void displayProgress(std::atomic<unsigned long> &count, std::mutex &m, unsigned int dataLength) {
	while (count != dataLength) {
		std::this_thread::sleep_for(std::chrono::seconds(1));
		m.lock();
		std::cout << "Processed: " << std::setw(20) << count << "/" << dataLength << " \r" << std::flush;
		m.unlock();
	}
}


std::string defaultOutputName(std::string inputName){
	std::string inputFile = std::string(std::move(inputName));
	std::vector<std::string> pathDirSplit;
	boost::split(pathDirSplit, inputFile, boost::is_any_of("/"));
	
	std::vector<std::string> fileExtSplit;
	boost::split(fileExtSplit, pathDirSplit.back(), boost::is_any_of("."));
	
	std::string directory;
	for(int i = 0; i < pathDirSplit.size() - 1; i++){
		directory += pathDirSplit[i];
		directory += "/";
	}
	
	std::string outputFile = directory + fileExtSplit[0] + "PES.dat";
	return outputFile;
}

#endif //RECOMORE_UTILS_H

