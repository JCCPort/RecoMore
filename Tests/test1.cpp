#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE TestModule1

#include <boost/test/unit_test.hpp>

#include <iostream>
#include <atomic>
#include <thread>
#include "../PEFinder/include/Utils.h"
#include "../PEFinder/include/PEFit.h"
#include "../PEFinder/include/ThreadPool.h"



class SystemTest1{
public:
	bool runOverTestData();
	void openNewRunFile(){newData = ReadRecoMoreOutput("../TestData/R185PES.dat");};
	void openOldRunFile(){oldData = ReadRecoMoreOutput("../TestData/R185PES_Reference.dat");};
	
	bool comparisons();
	
	SystemTest1() = default;
	
private:
	FitRun newData;
	FitRun oldData;
	
	double timeSimilarity = 1e-6;
	double ampSimilarity = 1e-6;
	double redChiSqSimilarity = 1e-3;
	
public:
	bool numEventsSame = true;
	bool numEventsSameRan = false;
	
	bool numPESInEventSame = true;
	bool numPESInEventSameRan = false;
	
	bool redChiSqSame = true;
	bool redChiSqSameRan = false;
	
	bool timesSame = true;
	bool timesSameRan = false;
	
	bool ampsSame = true;
	bool ampsSameRan = false;
};


bool SystemTest1::runOverTestData() {
	// Input arguments
	std::string inputFileName = "../TestData/R185.bin";
	std::string outputFileName;
	outputFileName = defaultOutputName(inputFileName);
	auto pdfDir = "../TestData/pdf/";
	saveWaveforms = false;
	DigitiserRun data       = ReadWCDataFile(inputFileName);
	int numThreads  = 8;
	int batchNumber = 8;
	
	// Global parameters
	skipChannels = {32, 36, 40, 44, 48, 52, 56, 60};
	pdfNSamples        = 105601;
	pdfSamplingRate    = 0.003125;
	pdfT0Sample        = 3201;
	pdfResidualRMS     = 0.827/1000;
	meanReducedChiSq   = 0;
	samplingRate2Inv   = 1.0f / (0.01f * pdfSamplingRate);
	pdfT0SampleConv    = pdfT0Sample;
	WFSigThresh        = 0.005;
	maxPEs             = 100;
	ampDiff            = 0;
	timeDiff           = 0;
	baselineDiff       = 0;
	sysProcPECount     = 0;
	sysProcWFCount     = 0;
	
	// Processing
	std::shared_ptr<SyncFile> file;
	file = std::make_shared<SyncFile>(outputFileName, text);
	
	Writer writer(file);
	
	static std::atomic<unsigned long> count{0};
	std::mutex m;
	
	std::vector<std::vector<double>> idealWaveforms{64};
	for (int ch = 0; ch < 64; ch++) {
		if (std::count(skipChannels.begin(), skipChannels.end(), ch)) {
			continue;
		}
		idealWaveforms.at(ch) = readIdealWFs(ch, 10, pdfDir, pdfNSamples);
	}
	
	std::thread progressThread(displayProgress, std::reference_wrapper(count), std::reference_wrapper(m),
	                           data.getEvents().size());
	
	if (numThreads == 1) {
		batchFitEvents(data.getEvents(), count, m, &idealWaveforms, file);
	} else {
		// Determining how many events each thread should run over.
		unsigned int threadRepeatCount[batchNumber];
		unsigned int threadsWithExtra = data.getEvents().size() % batchNumber;
		unsigned int minRepeatsPerThread = data.getEvents().size() / batchNumber;
		
		for (unsigned int i = 0; i < batchNumber; i++) {
			if (i < threadsWithExtra) {
				threadRepeatCount[i] = minRepeatsPerThread + 1;
			} else {
				threadRepeatCount[i] = minRepeatsPerThread;
			}
		}
		
		BS::thread_pool pool(numThreads);
		
		unsigned int eventPos = 0;
		for(int i = 0; i < batchNumber; i++){
			std::vector passData = slice(data.getEvents(), eventPos, eventPos + threadRepeatCount[i] - 1);
			pool.push_task(batchFitEvents, passData, std::reference_wrapper(count), std::reference_wrapper(m), &idealWaveforms,
			               file);
			eventPos += threadRepeatCount[i];
		}
	}
	
	progressThread.join();
	
	file->closeFile();
	
	return true;
}


bool SystemTest1::comparisons(){
	// Check number of events is the same
	if(newData.getEvents().size() != oldData.getEvents().size()){
		numEventsSame = false;
		std::cout << newData.getEvents().size() << "\n" << oldData.getEvents().size() << std::endl;
		return false;
	}
	numEventsSameRan = true;
	
	// Loop over events
	std::vector<unsigned int> newEventIDs = newData.getEventIDs();
	
	for(unsigned int i : newEventIDs){
		// Loop over channels in event
		FitEvent newEventData = newData.getEvent(i);
		FitEvent oldEventData = oldData.getEvent(i);
		
		std::vector<unsigned int> newEventChannels = newEventData.getChannelIDs();
		std::vector<unsigned int> oldEventChannels = oldEventData.getChannelIDs();
		for(unsigned int j : newEventChannels){
			
			if(not((std::count(newEventChannels.begin(), newEventChannels.end(), j)) and
					std::count(oldEventChannels.begin(), oldEventChannels.end(), j))) {
				continue;
			}
			
			auto newChannelData = newEventData.getChannel(j);
			auto oldChannelData = oldEventData.getChannel(j);
			
			// Check number of PEs in channel fit is the same
			numPESInEventSame = numPESInEventSame and (newChannelData.PEs.size() == oldChannelData.PEs.size());
			if((newChannelData.PEs.size() != oldChannelData.PEs.size())){
				std::cout << "Mis-match in number of PEs for event " << newEventData.ID << ", channel " << newChannelData.ID << std::endl;
				std::cout << newChannelData.PEs.size() << "\t" << oldChannelData.PEs.size() << std::endl;
				return false;
			}
			numPESInEventSameRan = true;
			
			// Check reduced ChiSqs are the same
			redChiSqSame = redChiSqSame and (std::abs(newChannelData.reducedChiSq - oldChannelData.reducedChiSq) <= redChiSqSimilarity);
			if(std::abs(newChannelData.reducedChiSq - oldChannelData.reducedChiSq) > redChiSqSimilarity){
				std::cout << "Mis-match in redChiSq for event " << newEventData.ID << ", channel " << newChannelData.ID << std::endl;
				std::cout << newChannelData.reducedChiSq << "\t" << oldChannelData.reducedChiSq << std::endl;
				return false;
			}
			redChiSqSameRan = true;

			
			// Checking if fit values for each PE is the same
			for(int k = 0; k < newChannelData.PEs.size(); k++){
				auto newPE = newChannelData.PEs[k];
				auto oldPE = oldChannelData.PEs[k];
				timesSame = timesSame and (std::abs(newPE.time - oldPE.time) <= timeSimilarity);
				if(std::abs(newPE.time - oldPE.time) > timeSimilarity){
					std::cout << "Mis-match in time for event " << newEventData.ID << ", channel " << newChannelData.ID << std::endl;
					std::cout << newPE.time << "\t" << oldPE.time << std::endl;
					return false;
				}
				timesSameRan = true;
				
				ampsSame = ampsSame and (std::abs(newPE.amplitude - oldPE.amplitude) <= ampSimilarity);
				if(std::abs(newPE.amplitude - oldPE.amplitude) > ampSimilarity){
					std::cout << "Mis-match in amplitude for event " << newEventData.ID << ", channel " << newChannelData.ID << std::endl;
					std::cout << newPE.amplitude << "\t" << oldPE.amplitude << std::endl;
					return false;
				}
				ampsSameRan = true;
			}
		}
	}
	return true;
}

SystemTest1 systemTest1;

BOOST_AUTO_TEST_SUITE(SystemTest)

BOOST_AUTO_TEST_CASE(CheckFittingRan)
{
	BOOST_CHECK(systemTest1.runOverTestData());
	std::cout << "Ran fitting over test data" << std::endl;
}

BOOST_AUTO_TEST_CASE(RunFilesOpened)
{
	systemTest1.openNewRunFile();
	systemTest1.openOldRunFile();
}

BOOST_AUTO_TEST_CASE(ComparisonsRan)
{
	systemTest1.comparisons();
}

BOOST_AUTO_TEST_CASE(CheckEventSize)
{
	BOOST_CHECK((systemTest1.numEventsSame and systemTest1.numEventsSameRan));
}

BOOST_AUTO_TEST_CASE(CheckNumPESInEvent)
{
	BOOST_CHECK((systemTest1.numPESInEventSame and systemTest1.numPESInEventSameRan));
}

BOOST_AUTO_TEST_CASE(CheckRedChiSqs)
{
	BOOST_CHECK((systemTest1.redChiSqSame and systemTest1.redChiSqSameRan));
}

BOOST_AUTO_TEST_CASE(CheckTimesSame)
{
	BOOST_CHECK((systemTest1.timesSame and systemTest1.timesSameRan));
}

BOOST_AUTO_TEST_CASE(CheckAmpsSame)
{
	BOOST_CHECK((systemTest1.ampsSame and systemTest1.ampsSameRan));
}


BOOST_AUTO_TEST_SUITE_END()
