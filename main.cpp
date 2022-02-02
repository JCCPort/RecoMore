#include <iostream>
#include <atomic>
#include <mutex>
#include <thread>
#include "Utils.h"
#include "include/DataReading.h"
#include "include/DataStructures.h"
#include "include/PEFit.h"
#include "Globals.h"
#include "include/DataWriting.h"
#include <boost/algorithm/string.hpp>

#include "include/ThreadPool.h"



int main(int argc, char** argv) {
	// TODO(josh): I strongly suspect we're being slowed down due to competing access for ideal waveform data.
	//  somewhere there's an inefficiency as CPU usage isn't being maxed out

	std::string inputFile = std::string(argv[1]);
	std::vector<std::string> splitString;
	boost::split(splitString, inputFile, boost::is_any_of("/"));
	std::string pdfDir = std::string(argv[2]);

	std::vector<std::string> splitString2;
	boost::split(splitString2, splitString.back(), boost::is_any_of("."));

	std::string directory;
	for(int i = 0; i < splitString.size() - 1; i++){
		directory += splitString[i];
		directory += "/";
	}

	std::string outputFile = directory + splitString2[0] + "PES.dat";

	// TODO(josh): Way to exclude specific channels from being read
	WCData data = ReadWCDataFile(inputFile);
	// TODO(josh): Add error checking to if the data file is corrupted/invalid

	unsigned int numThreads = 16;
	unsigned int batchNumber = 200;
	unsigned int numChannels = 16;
	static std::atomic<unsigned long> count{0};
	std::mutex m;

	auto file = std::make_shared<SyncFile>(outputFile);
	Writer writer(file);

	// TODO(josh): Use info read in from wavecatcher data file to determine what channels ideal PDFs to load.
	std::vector<std::vector<double>> idealWaveforms{numChannels};
	for (int ch = 0; ch < numChannels; ch++) {
		if ((ch == 32) || (ch == 36) || (ch == 40) || (ch == 44) || (ch == 48) || (ch == 52) || (ch == 56) ||
		    (ch == 60)) {
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
	std::cout << "Mean reduced chisq:\t" << meanReducedChisq << std::endl;

	return 0;
}
