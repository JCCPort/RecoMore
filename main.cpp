#include <iostream>
#include <atomic>
#include <mutex>
#include <thread>
#include "include/DataReading.h"
#include "include/DataStructures.h"
#include "include/PEFit.h"
#include "Globals.h"
#include "include/DataWriting.h"

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
                 std::vector<std::vector<double>> *idealWaveforms, const std::shared_ptr<SyncFile> &file) {
	for (auto &event: events) {
		m.lock();
		++count;
		m.unlock();
		fitPE(event, idealWaveforms, file);
	}
	return true;
}


int main() {
	WCData data = ReadWCDataFile("/home/josh/CLionProjects/RecoMore/R32.dat");
	unsigned int numThreads = 6;
	unsigned int batchSize = 1;
	// TODO(josh): Investigate why this code ran so quickly last night 23:45ish 18/10/2021
	static std::atomic<unsigned long> count{0};
	std::mutex m;

	auto file = std::make_shared<SyncFile>("R32PES.dat");
	Writer writer(file);

	std::vector<std::vector<double>> idealWaveforms{64};
	for (int ch = 0; ch < 64; ch++) {
		if ((ch == 32) || (ch == 36) || (ch == 40) || (ch == 44) || (ch == 48) || (ch == 52) || (ch == 56) ||
		    (ch == 60)) {
			continue;
		}
		idealWaveforms.at(ch) = readIdealWFs(ch, 10, "/home/josh/CLionProjects/RecoMore/pdf/", pdfNSamples);
	}


	thread_pool pool(numThreads);

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
