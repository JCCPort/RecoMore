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
	auto &arg= program.add_argument("--template_dir")
		.default_value(std::string("../templates/"))
		.help("Path for PE templates to use for fitting.");
	program.add_hidden_alias_for(arg, "--pdf_dir");
	program.add_argument("--n_threads")
		.default_value(1)
		.help("Number of threads to run RecoMore on.")
		.scan<'i', int>();
	program.add_argument("--num_batches")
		.default_value(8)
		.help("Number of batches to split run into. Smaller values will reduce overhead, "
			  "but increase the likelihood of one thread hanging.")
		.scan<'i', int>();
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
	program.add_argument("--sample_spacing")
		.default_value(0.3125f)
		.help("Sample spacing in data to process in ns.");

	// =========================================================================
	//    Mutually exclusive group for channels:
	//    Either provide a list of channels or a file containing them
	// =========================================================================
	auto &channelsGroup = program.add_mutually_exclusive_group(/* required= */ true);

	channelsGroup.add_argument("--channels-list")
				 .help("Space-separated list of channel numbers to process.")
				 .nargs(argparse::nargs_pattern::any)
				 .scan<'i', int>();

	channelsGroup.add_argument("--channels-file")
				 .help("Path to a file containing channel numbers.");
	// =========================================================================
	
	program.parse_args(argc, argv);
	
	auto inputFileName = program.get<std::string>("-i");
	std::string outputFileName;
	if(program.is_used("-o")){
		outputFileName = program.get<std::string>("-o");
	} else {
		outputFileName = defaultOutputName(program.get<std::string>("-i"));
	}

    if(outputFileName == inputFileName){
        throw std::runtime_error("Input and output file names are the same. Please specify a different output file name.");
    }

	auto templateDir = program.get<std::string>("--template_dir");
	unsigned int numThreads = program.get<int>("--n_threads");
	unsigned int batchNumber = program.get<int>("--num_batches");
	saveWaveforms = program.get<bool>("--save_waveforms");
    parameterTolerance = program.get<float>("--parameter_tolerance");
	const bool positivePulse = program.get<bool>("--positive_pulse");
	const auto sampleSpacing = program.get<float>("--sample_spacing");


	// =========================================================================
	//    Determine channels to use
	// =========================================================================
	std::vector<int> channels;
	if (program.is_used("--channels-list")) {
		channels = program.get<std::vector<int>>("--channels-list");
	} else {
		// Must be from --channels-file
		auto channelsPath = program.get<std::string>("--channels-file");
		channels = parseChannelsFromFile(channelsPath);
	}

	if (batchNumber < numThreads)
	{
		std::cout << "WARNING: Number of batches is less than number of threads. Not all threads will be used." << std::endl;
	}

	DigitiserRun data = ReadWCDataFile(inputFileName, positivePulse);

	std::shared_ptr<SyncFile> file;
	if(const bool textOutput = program.get<bool>("--txt-output"); !textOutput) {
	  file = std::make_shared<SyncFile>(outputFileName, binary);
	} else {
	  file = std::make_shared<SyncFile>(outputFileName, text);
	}

	Writer writer(file);

	static std::atomic<unsigned long> count{0};
	std::mutex                        progressTrackerLock;
	std::mutex						  meanReducedChiSqLock;

	// =========================================================================
	//    Construct PETemplates *only* for the specified channels
	// =========================================================================
	std::unordered_map<unsigned int, PETemplate> idealWaveforms;
	idealWaveforms.reserve(channels.size());

	// Build a PETemplate for each channel in "channels"
	for (int ch : channels) {
		idealWaveforms.emplace(
			ch,
			PETemplate(ch, templateInternalInterpFactor, templateDir, positivePulse)
		);
	}

	std::thread progressThread(displayProgress, std::reference_wrapper(count), std::reference_wrapper(progressTrackerLock),
	                           data.getEvents().size());

	// Record current time for timing purposes
	auto start = std::chrono::high_resolution_clock::now();

	if (numThreads == 1) {
		std::cout << "Processing with one thread." << std::endl;
		batchFitEvents(data.getEvents(), std::reference_wrapper(count), std::reference_wrapper(meanReducedChiSqLock), idealWaveforms, file, sampleSpacing);
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
			               file, sampleSpacing);
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
