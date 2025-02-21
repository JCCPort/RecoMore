#include <iostream>
#include <atomic>
#include <thread>
#include "include/Utils.h"
#include "include/DataReading.h"
#include "include/DataStructures.h"
#include "include/PETemplate.h"
#include "include/PEFit.h"
#include "Globals.h"
#include "include/DataWriting.h"

#include "include/ThreadPool.h"
#include "include/argparse.h"


int main(int argc, char** argv) {
	argparse::ArgumentParser program("\n\n"
	                                 "██████╗ ███████╗ ██████╗ ██████╗ ███╗   ███╗ ██████╗ ██████╗ ███████╗\n"
	                                 "██╔══██╗██╔════╝██╔════╝██╔═══██╗████╗ ████║██╔═══██╗██╔══██╗██╔════╝\n"
	                                 "██████╔╝█████╗  ██║     ██║   ██║██╔████╔██║██║   ██║██████╔╝█████╗  \n"
	                                 "██╔══██╗██╔══╝  ██║     ██║   ██║██║╚██╔╝██║██║   ██║██╔══██╗██╔══╝  \n"
	                                 "██║  ██║███████╗╚██████╗╚██████╔╝██║ ╚═╝ ██║╚██████╔╝██║  ██║███████╗\n"
	                                 "╚═╝  ╚═╝╚══════╝ ╚═════╝ ╚═════╝ ╚═╝     ╚═╝ ╚═════╝ ╚═╝  ╚═╝╚══════╝\n"
	                                 "                                                                     ");
	program.add_description("Multi-threaded implementation of the RecoZoR PE finding algorithm.");
	program.add_argument("-i", "--input")
		.required()
		.help("Path to raw data file (.dat or .bin).");
	program.add_argument("-o", "--output")
		.help("Path for output reco file. Defaults to input file name with 'PES' appended.");
	program.add_argument("--txt-output")
         	.default_value(false)
	        .implicit_value(true)
		.help("Output reco file saved as text. Binary is default.");
	program.add_argument("--pdf_dir")
		.default_value(std::string("../pdf/"))
		.help("Path for ideal PDFs to use for fitting.");
	program.add_argument("--n_threads")
		.default_value(1)
		.help("Number of threads to run RecoMore on.")
		.scan<'i', int>();
	program.add_argument("--num_batches")
		.default_value(8)
		.help("Number of batches to split run into. Smaller values will reduce overhead, "
			  "but increase the likelihood of one thread hanging.")
		.scan<'i', int>();
	program.add_argument("--skip_channels")
		.nargs(argparse::nargs_pattern::any)
		.default_value(skipChannels)
		.scan<'i', int>()
		.help("Channels to skip. Space separated.");
	program.add_argument("--save_waveforms")
		.default_value(false)
		.help("Save waveforms with initial and final fits to csv files.");
    program.add_argument("--parameter_tolerance")
        .default_value(1e-8f)
        .help("Tolerance for the fit.")
        .scan<'e', float>();
	program.add_argument("--positive_pulse")
		.default_value(false)
		.implicit_value(true)
		.help("Use flag if SiPM pulse is positive. Default is negative.");
	
	program.parse_args(argc, argv);
	
	std::string inputFileName = program.get<std::string>("-i");
	std::string outputFileName;
	if(program.is_used("-o")){
		outputFileName = program.get<std::string>("-o");
	} else {
		outputFileName = defaultOutputName(program.get<std::string>("-i"));
	}

    if(outputFileName == inputFileName){
        throw std::runtime_error("Input and output file names are the same. Please specify a different output file name.");
    }

	auto pdfDir = program.get<std::string>("--pdf_dir");
	unsigned int numThreads = program.get<int>("--n_threads");
	unsigned int batchNumber = program.get<int>("--num_batches");
	skipChannels = program.get<std::vector<int>>("--skip_channels");
	saveWaveforms = program.get<bool>("--save_waveforms");
    parameterTolerance = program.get<float>("--parameter_tolerance");
	const bool positivePulse = program.get<bool>("--positive_pulse");
	DigitiserRun data = ReadWCDataFile(inputFileName, positivePulse);

	if (batchNumber < numThreads)
	{
		std::cout << "WARNING: Number of batches is less than number of threads. Not all threads will be used." << std::endl;
	}
 
	std::shared_ptr<SyncFile> file;
	const bool textOutput = program.get<bool>("--txt-output");
	if(!textOutput) {
	  file = std::make_shared<SyncFile>(outputFileName, binary);
	} else {
	  file = std::make_shared<SyncFile>(outputFileName, text);
	}

	Writer writer(file);

	static std::atomic<unsigned long> count{0};
	std::mutex                        progressTrackerLock;
	std::mutex						  meanReducedChiSqLock;

	// TODO(josh): Implement checking of channels from the data file to ensure the correct number of ideal waveforms are read in.
	unsigned int numChannels = 16;
	std::unordered_map<unsigned int, PETemplate> idealWaveforms;
	idealWaveforms.reserve(numChannels);  // Optional but can help performance

	for (int ch = 0; ch < numChannels; ++ch) {
		// If ch is in skipChannels, we continue (skip it)
		if (std::ranges::count(skipChannels, ch) > 0) {
			continue;
		}

		// Otherwise, construct the PETemplate and store it under key = ch
		idealWaveforms.emplace(
			ch,
			PETemplate(ch, pdfInternalInterpFactor, pdfDir, positivePulse)
		);
	}

	std::thread progressThread(displayProgress, std::reference_wrapper(count), std::reference_wrapper(progressTrackerLock),
	                           data.getEvents().size());

	// Record current time for timing purposes
	auto start = std::chrono::high_resolution_clock::now();

	if (numThreads == 1) {
		std::cout << "Processing with one thread." << std::endl;
		batchFitEvents(data.getEvents(), std::reference_wrapper(count), std::reference_wrapper(meanReducedChiSqLock), idealWaveforms, file);
	} else {
		std::cout << "Processing with " << numThreads << " threads." << std::endl;
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
			std::vector<DigitiserEvent> passData = slice(data.getEvents(), eventPos, eventPos + threadRepeatCount[i] - 1);
			pool.push_task(batchFitEvents, passData, std::reference_wrapper(count), std::reference_wrapper(meanReducedChiSqLock), idealWaveforms,
			               file);
			eventPos += threadRepeatCount[i];
		}
	}

	progressThread.join();

	file->closeFile();
	std::cout << "Mean reduced ChiSq:\t" << meanReducedChiSq << std::endl;

	// Record end time for timing purposes
	auto end = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed = end - start;
	std::cout << "Time taken: " << elapsed.count() << "s" << std::endl;

	return 0;
}
