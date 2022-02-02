#ifndef RECOMORE_DATAWRITING_H
#define RECOMORE_DATAWRITING_H

#include <fstream>
#include <memory>
#include <mutex>
#include <string>
#include <utility>
#include <iostream>
#include "DataStructures.h"

class SyncFile {
public:
	explicit SyncFile(const std::string &path) : path_(path) { myFile_.open(path, std::ios_base::out); }

	void write(const std::string &dataToWrite) {
		std::lock_guard<std::mutex> lock(writerMutex_);
		writeCache_ += dataToWrite;
		cachedStates_++;
		if (cachedStates_ == 1e2) {
			myFile_ << writeCache_;
			writeCache_.clear();
			cachedStates_ = 0;
		}
	}

	void closeFile() {
		std::lock_guard<std::mutex> lock(writerMutex_);
		myFile_ << writeCache_;
		writeCache_.clear();
		cachedStates_ = 0;
		std::cout << "Finished writing to file " << path_ << std::endl;
		myFile_.close();

		std::ifstream file2(path_);
		if (file2.peek() == std::ifstream::traits_type::eof()) {
			std::cout << "Output file is empty" << std::endl;
		}
	};

private:
	unsigned int cachedStates_ = 0;
	std::string writeCache_;
	std::ofstream myFile_;
	std::mutex writerMutex_;
	std::string path_;
};

class Writer {
public:
	explicit Writer(std::shared_ptr<SyncFile> sf) : _sf(std::move(sf)) {}

	static std::string writeWaveformInfo(const ChannelFitData &wfDat);

	static std::string writeFitPE(PEData PEVar);

	void writeEventInfo(const EventFitData &evData);

private:
	std::shared_ptr<SyncFile> _sf;
};

#endif //RECOMORE_DATAWRITING_H
