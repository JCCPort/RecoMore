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
#include <utility>

#include "include/ThreadPool.h"
#include "include/argparse.h"


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


int main(int argc, char** argv) {
	argparse::ArgumentParser program("RecoMore");
	program.add_argument("-i", "--input")
		.required()
		.help("Path to raw data file (.dat or .bin).");
	program.add_argument("-o", "--output")
		.default_value(defaultOutputName(program.get<std::string>("-i")))
		.help("Path for output reco file.");
	program.add_argument("--pdf_dir")
		.default_value("../pdf/")
		.help("Path for ideal PDFs to use for fitting.");
	program.add_argument("--n_threads")
		.default_value("1")
		.help("Number of threads to run RecoMore on.")
		.scan<'d', unsigned int>();
	program.add_argument("--batch_size")
		.default_value(200)
		.help("Number of batches to split run into. Smaller values will reduce overhead, "
			  "but increase the likelihood of one thread hanging.")
		.scan<'d', unsigned int>();
	program.add_argument("--skip_channels")
		.nargs(argparse::nargs_pattern::any)
		.default_value(std::vector<unsigned int>{})
		.help("Channels to skip.");
	
	std::string inputFileName = program.get<std::string>("-i");
	std::string outputFileName = program.get<std::string>("-o");
	std::string pdfDir = program.get<std::string>("--pdf_dir");
	unsigned int numThreads = program.get<unsigned int>("--n_threads");
	unsigned int batchNumber = program.get<unsigned int>("--batch_size");
	skipChannels = program.get<std::vector<unsigned int>>("--skip_channels");
	
	WCData data = ReadWCDataFile(inputFileName);
	auto file = std::make_shared<SyncFile>(outputFileName);
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
	std::cout << "Mean reduced ChiSq:\t" << meanReducedChisq << std::endl;

	return 0;
}
