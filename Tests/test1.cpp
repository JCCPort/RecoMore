#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE TestModule1

#include <boost/test/unit_test.hpp>

#include <iostream>
#include <atomic>
#include <thread>
#include "../PEFinder/include/Utils.h"
#include "../PEFinder/include/DataReading.h"
#include "../PEFinder/include/DataStructures.h"
#include "../PEFinder/include/PEFit.h"
#include "../PEFinder/Globals.h"
#include "../PEFinder/include/DataWriting.h"

#include "../PEFinder/include/ThreadPool.h"
#include "../PEFinder/include/argparse.h"



class SystemTest1{
public:
	bool runOverTestData();
	void openNewRunFile(){newData = ReadRecoMoreOutput("../TestData/R185PES.dat");};
	void openOldRunFile(){oldData = ReadRecoMoreOutput("../TestData/R185PES_Reference.dat");};
	
	bool comparisons();
	
	SystemTest1(){}
	
private:
	FitData newData;
	FitData oldData;
	
	double timeSimilarity = 1e-6;
	double ampSimilarity = 1e-6;
	double redChiSqSimilarity = 1e-3;
	int numPESSimilarity = 5;
	
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
	WCData data = ReadWCDataFile(inputFileName);
	int numThreads = 8;
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
		fitBatchPEs(data.getEvents(), count, m, &idealWaveforms, file);
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
		
		ThreadPool pool(numThreads);
		
		unsigned int eventPos = 0;
		for(int i = 0; i < batchNumber; i++){
			std::vector passData = slice(data.getEvents(), eventPos, eventPos + threadRepeatCount[i] - 1);
			pool.push_task(fitBatchPEs, passData, std::reference_wrapper(count), std::reference_wrapper(m), &idealWaveforms,
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
	if(newData.getFitEvents().size() != oldData.getFitEvents().size()){
		numEventsSame = false;
		std::cout << newData.getFitEvents().size() << "\n" << oldData.getFitEvents().size() << std::endl;
		return false;
	}
	numEventsSameRan = true;
	
	// Loop over events
	for(int i = 1; i <= newData.getFitEvents().size(); i++){
		// Loop over channels in event
		EventFitData newEventData = newData.getEventFit(i);
		EventFitData oldEventData = oldData.getEventFit(i);
		for(int j = 0; j < newData.getFitEvents()[i].SiPM.size(); j++){
			std::vector<unsigned short> newEventChannels = newEventData.getChannels();
			std::vector<unsigned short> oldEventChannels = oldEventData.getChannels();
			if(not((std::count(newEventChannels.begin(), newEventChannels.end(), j)) and
					std::count(oldEventChannels.begin(), oldEventChannels.end(), j))) {
				continue;
			}
			
			auto newChannelData = newEventData.getChannel(j);
			auto oldChannelData = oldEventData.getChannel(j);
			
			// Check number of PEs in channel fit is the same
			numPESInEventSame = numPESInEventSame and (newChannelData.pes.size() == oldChannelData.pes.size());
			if((newChannelData.pes.size() != oldChannelData.pes.size())){
				std::cout << "Mis-match in number of PEs for event " << newData.getFitEvents()[i].eventID << ", channel " << newChannelData.ch << std::endl;
				std::cout << newChannelData.pes.size() << "\t" << oldChannelData.pes.size() << std::endl;
				return false;
			}
			numPESInEventSameRan = true;
			
			// Check reduced ChiSqs are the same
			redChiSqSame = redChiSqSame and (std::abs(newChannelData.redChiSq - oldChannelData.redChiSq) <= redChiSqSimilarity);
			if(std::abs(newChannelData.redChiSq - oldChannelData.redChiSq) > redChiSqSimilarity){
				std::cout << "Mis-match in redChiSq for event " << newData.getFitEvents()[i].eventID << ", channel " << newChannelData.ch << std::endl;
				std::cout << newChannelData.redChiSq << "\t" << oldChannelData.redChiSq << std::endl;
				return false;
			}
			redChiSqSameRan = true;

			
			// Checking if fit values for each PE is the same
			for(int k = 0; k < newChannelData.pes.size(); k++){
				timesSame = timesSame and (std::abs(newChannelData.pes[k].time - oldChannelData.pes[k].time) <= timeSimilarity);
				if(std::abs(newChannelData.pes[k].time - oldChannelData.pes[k].time) > timeSimilarity){
					std::cout << "Mis-match in time for event " << newData.getFitEvents()[i].eventID << ", channel " << newChannelData.ch << std::endl;
					std::cout << newChannelData.pes[k].time << "\t" << oldChannelData.pes[k].time << std::endl;
					return false;
				}
				timesSameRan = true;
				
				ampsSame = ampsSame and (std::abs(newChannelData.pes[k].amplitude - oldChannelData.pes[k].amplitude) <= ampSimilarity);
				if(std::abs(newChannelData.pes[k].amplitude - oldChannelData.pes[k].amplitude) > ampSimilarity){
					std::cout << "Mis-match in amplitude for event " << newData.getFitEvents()[i].eventID << ", channel " << newChannelData.ch << std::endl;
					std::cout << newChannelData.pes[k].amplitude << "\t" << oldChannelData.pes[k].amplitude << std::endl;
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
	systemTest1.openNewRunFile();
	systemTest1.openOldRunFile();
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
