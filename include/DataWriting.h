#ifndef RECOMORE_DATAWRITING_H
#define RECOMORE_DATAWRITING_H

#include <fstream>
#include <memory>
#include <mutex>
#include <string>
#include <utility>
#include "DataStructures.h"

class SyncFile {
public:
	explicit SyncFile(const std::string &path) : _path(path) { myFile.open(path, std::ios_base::out); }

	void write(const std::string &dataToWrite) {
		std::lock_guard<std::mutex> lock(_writerMutex);
		writeCache_ += dataToWrite;
		cachedStates_++;
		if(cachedStates_ == 1e6){
			myFile << writeCache_;
			writeCache_.clear();
			cachedStates_ = 0;
		}
	}

	void closeFile() {
		std::lock_guard<std::mutex> lock(_writerMutex);
		myFile << writeCache_;
		writeCache_.clear();
		cachedStates_ = 0;
		myFile.close();
	};

private:
	unsigned int cachedStates_ = 0;
	std::string writeCache_;
	std::ofstream myFile;
	std::string _path;
	std::mutex _writerMutex;
};

class Writer {
public:
	explicit Writer(std::shared_ptr<SyncFile> sf) : _sf(std::move(sf)) {}

	void write(std::string &dataToWrite) { _sf->write(dataToWrite); }

	void writeWaveformInfo(const ChannelFitData&);
	void writeFitPE(PEData);
	void writeEventInfo(const EventFitData& evData);

private:
	std::shared_ptr<SyncFile> _sf;
};

#endif //RECOMORE_DATAWRITING_H
