#include <iostream>
#include <atomic>
#include <mutex>
#include <thread>
#include "include/DataReading.h"
#include "include/DataStructures.h"
#include "include/PEFit.h"
#include "Globals.h"
#include "include/DataWriting.h"
#include <boost/algorithm/string.hpp>

#include "include/ThreadPool.h"

/**
 *
 * @tparam T
 * @param v
 * @param m
 * @param n
 * @return
 */
template<typename T>
std::vector<T> slice(std::vector<T> const &v, unsigned int m, unsigned int n) {
	auto first = v.cbegin() + m;
	auto last = v.cbegin() + n + 1;

	std::vector<T> vec(first, last);
	return vec;
}

/**
 *
 * @param count
 * @param m
 * @param dataLength
 */
void displayProgress(std::atomic<unsigned long> &count, std::mutex &m, unsigned int dataLength) {
	while (count != dataLength) {
		std::this_thread::sleep_for(std::chrono::seconds(1));
		m.lock();
		std::cout << "Processed: " << count << "/" << dataLength << std::endl;
		m.unlock();

	}

}

/**
 *
 * @param events
 * @param count
 * @param m
 * @param idealWaveforms
 * @param file
 * @return
 */
bool fitBatchPEs(const std::vector<EventData> &events, std::atomic<unsigned long> &count, std::mutex &m,
                 const std::vector<std::vector<double>> * idealWaveforms, const std::shared_ptr<SyncFile> &file) {
	for (auto &event: events) {
		m.lock();
		++count;
		m.unlock();
		fitPE(event, idealWaveforms, file);
	}
	return true;
}


int main(int argc, char** argv) {
	// TODO(josh): I strongly suspect we're being slowed down due to competing access for ideal waveform data.
	//  somewhere there's an inefficiency as CPU usage isn't being maxed out

	std::string inputFile = std::string(argv[1]);
	std::vector<std::string> splitString;
	boost::split(splitString, inputFile, boost::is_any_of("."));
	std::string pdfDir = std::string(argv[2]);
	std::string outputFile = splitString[0] + "PES.dat";

	WCData data = ReadWCDataFile(inputFile);

	unsigned int numThreads = 8;
	unsigned int batchSize = 20;
	static std::atomic<unsigned long> count{0};
	std::mutex m;

	auto file = std::make_shared<SyncFile>(outputFile);
	Writer writer(file);

	std::vector<std::vector<double>> idealWaveforms{64};
	for (int ch = 0; ch < 64; ch++) {
		if ((ch == 32) || (ch == 36) || (ch == 40) || (ch == 44) || (ch == 48) || (ch == 52) || (ch == 56) ||
		    (ch == 60)) {
			continue;
		}
		idealWaveforms.at(ch) = readIdealWFs(ch, 10, pdfDir, pdfNSamples);
	}


	ThreadPool pool(numThreads);

	unsigned int eventPos = 0;
	std::cout << data.getEvents().size() << std::endl;
	std::thread progressThread(displayProgress, std::reference_wrapper(count), std::reference_wrapper(m),
	                           data.getEvents().size());

	for(int i = 0; (int)(i < data.getEvents().size()/batchSize); i++) {
		std::vector passData = slice(data.getEvents(), eventPos, eventPos + batchSize - 1);
		pool.push_task(fitBatchPEs, passData, std::reference_wrapper(count), std::reference_wrapper(m), &idealWaveforms,
		               file);
		eventPos += batchSize;
	}

	progressThread.join();
	file->closeFile();

	return 0;
}
