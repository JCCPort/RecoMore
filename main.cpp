#include <iostream>
#include <atomic>
#include <mutex>
#include <thread>
#include "include/DataReading.h"
#include "include/DataStructures.h"
#include "include/PEFit.h"
#include "Globals.h"

template<typename T>
std::vector<T> slice(std::vector<T> const &v, unsigned int m, unsigned int n) {
	auto first = v.cbegin() + m;
	auto last = v.cbegin() + n + 1;

	std::vector<T> vec(first, last);
	return vec;
}

bool fitBatchPEs(const std::vector<EventData> &events, std::atomic<unsigned long> &count, std::mutex &m,
                 const std::shared_ptr<std::vector<EventFitData>> &PEList,
                 const std::vector<std::vector<double>> &idealWaveforms) {

	for (auto &event: events) {
		m.lock();
		++count;
		m.unlock();
		fitPE(event, PEList, idealWaveforms);
	}

	return true;
}


int main() {
	WCData data = ReadWCDataFile("/home/josh/CLionProjects/RecoMore/R45.dat");
	unsigned int numThreads = 1;

	static std::atomic<unsigned long> count{0};
	std::mutex m;

	auto PEList = std::make_shared<std::vector<EventFitData>>();

	std::vector<std::vector<double>> idealWaveforms{64};
	for (int ch = 0; ch < 64; ch++) {
		if((ch == 32) || (ch == 36) || (ch == 40) || (ch == 44) || (ch == 48) || (ch == 52) || (ch == 56) || (ch == 60)){
			continue;
		}
		idealWaveforms.at(ch) = readIdealWFs(ch, 10, "/home/josh/CLionProjects/RecoMore/pdf/", pdfNSamples);
	}


	// If only one thread is being used there's no need to use std::thread, just call batchRun conventionally.
	if (numThreads == 1) {
		fitBatchPEs(data.getEvents(), count, m, PEList, idealWaveforms);
		return 1;
	}

	// Determining how many rays each thread should simulate.
	unsigned int threadRepeatCount[numThreads];
	unsigned int threadsWithExtra = data.getEvents().size() % numThreads;
	unsigned int minRepeatsPerThread = data.getEvents().size() / numThreads;

	for (unsigned int i = 0; i < numThreads; i++) {
		if (i < threadsWithExtra) {
			threadRepeatCount[i] = minRepeatsPerThread + 1;
		} else {
			threadRepeatCount[i] = minRepeatsPerThread;
		}
	}

	std::vector<std::thread> threads;
	// Carrying out the multi-threaded simulations.
	unsigned int eventPos = 0;
	std::cout << data.getEvents().size() << std::endl;
	for (unsigned int i = 0; i < numThreads; i++) {
		std::vector passData = slice(data.getEvents(), eventPos, eventPos + threadRepeatCount[i] - 1);
		std::thread t(fitBatchPEs, passData, std::reference_wrapper(count), std::reference_wrapper(m),
		              std::reference_wrapper(PEList), idealWaveforms);
		threads.push_back(std::move(t));
		eventPos += threadRepeatCount[i];
	}

	// Waiting for all threads to complete before continuing code execution.
	for (auto &th: threads) {
		th.join();
	}

	return 1;
}
