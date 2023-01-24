#include <iostream>
#include <atomic>
#include <thread>
#include "include/Utils.h"
#include "include/DataReading.h"
#include "include/DataStructures.h"
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
	
	program.parse_args(argc, argv);
	
	std::string inputFileName = program.get<std::string>("-i");
	std::string outputFileName;
	if(program.is_used("-o")){
		outputFileName = program.get<std::string>("-o");
	} else {
		outputFileName = defaultOutputName(program.get<std::string>("-i"));
	}
	auto pdfDir = program.get<std::string>("--pdf_dir");
	unsigned int numThreads = program.get<int>("--n_threads");
	unsigned int batchNumber = program.get<int>("--num_batches");
	skipChannels = program.get<std::vector<int>>("--skip_channels");
	saveWaveforms = program.get<bool>("--save_waveforms");
	WCData data = ReadWCDataFile(inputFileName);
	std::shared_ptr<SyncFile> file;
	bool textOutput = program.get<bool>("--txt-output");
	std::cout<<textOutput<<std::endl;
	if(!textOutput) {
	  file = std::make_shared<SyncFile>(outputFileName, binary);
	} else {
	  file = std::make_shared<SyncFile>(outputFileName, text);
	}
	
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
	std::cout << "Mean reduced ChiSq:\t" << meanReducedChiSq << std::endl;

	return 0;
}
