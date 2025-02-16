#ifndef RECOMORE_DATAWRITING_H
#define RECOMORE_DATAWRITING_H

#include <fstream>
#include <memory>
#include <mutex>
#include <string>
#include <utility>
#include <iostream>
#include <boost/archive/binary_oarchive.hpp>
#include "DataStructures.h"

enum WriteMode{
	text = 0,
	binary = 1
};


class SyncFile {
public:
	explicit SyncFile(const std::string &path, WriteMode writeMode) : path_(path), writeMode_(writeMode) { myFile_.open(path, std::ios_base::out); }

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
	
	void binaryWrite(const FitEvent &evData) {
		std::lock_guard<std::mutex> lock(writerMutex_);
		binaryWriteCache_.push_back(evData);
	}

	void closeFile() {
		if(writeMode_ == text){
			std::lock_guard<std::mutex> lock(writerMutex_);
			myFile_ << writeCache_;
			writeCache_.clear();
			cachedStates_ = 0;
		} else if(writeMode_ == binary){
			for(int i = 0; i < binaryWriteCache_.size(); i++){ // TODO(josh): Decide on a better way to handle this abomination.
				for(int j = 0; j < binaryWriteCache_[i].channels.size(); j++){
					for(int k = 0; k < binaryWriteCache_[i].channels[j].PEs.size(); k++){
						binaryWriteCache_[i].channels[j].PEs[k].time      = binaryWriteCache_[i].channels[j].PEs[k].time;
						binaryWriteCache_[i].channels[j].PEs[k].timeError   = binaryWriteCache_[i].channels[j].PEs[k].timeError;
						binaryWriteCache_[i].channels[j].PEs[k].initialTime = binaryWriteCache_[i].channels[j].PEs[k].timeError;
					}
				}
			}
			boost::archive::binary_oarchive oa(myFile_);
			oa << binaryWriteCache_;
		}

		std::cout << "Finished writing to file " << path_ << std::endl;
		myFile_.close();

		std::ifstream file2(path_);
		if (file2.peek() == std::ifstream::traits_type::eof()) {
			std::cout << "Output file is empty" << std::endl;
		}
	};
	
	WriteMode getWriteMode() const {
		return writeMode_;
	}

private:
	WriteMode writeMode_;
	unsigned int cachedStates_ = 0;
	std::string writeCache_;
	
	std::vector<FitEvent> binaryWriteCache_;
	
	std::ofstream myFile_;
	std::mutex writerMutex_;
	std::string path_;
};

class Writer {
public:
	explicit Writer(std::shared_ptr<SyncFile> sf) : _sf(std::move(sf)) {writeMode_ = _sf->getWriteMode();}

	static std::string writeWaveformInfo(const FitChannel &);

	static std::string writeFitPE(Photoelectron);

	void writeEventInfo(const FitEvent &);
	
private:
	std::shared_ptr<SyncFile> _sf;
	WriteMode writeMode_;
};

#endif //RECOMORE_DATAWRITING_H
