#include <iostream>
#include <c++/10/atomic>
#include <c++/10/mutex>
#include <c++/10/memory>
#include <c++/10/thread>
#include "include/DataReading.h"
#include "include/DataStructures.h"

template<typename T>
std::vector<T> slice(std::vector<T> const &v, unsigned int m, unsigned int n)
{
	auto first = v.cbegin() + m;
	auto last = v.cbegin() + n + 1;

	std::vector<T> vec(first, last);
	return vec;
}

void fitPE(Waveform wf, const std::shared_ptr<std::vector<EventData>>& PEList){

}

bool fitBatchPEs(std::vector<Waveform> wfs, std::atomic<unsigned long> &count, std::mutex &m, const std::shared_ptr<std::vector<EventData>>& PEList) {

	for (auto& wf : wfs) {
		m.lock();
		++count;
		m.unlock();
		fitPE(wf, PEList);
	}

	return true;
}


int main() {
	WCData data = ReadWCDataFile("/home/josh/CLionProjects/RecoMore/R43.dat");
	unsigned int numThreads = 10;

	static std::atomic<unsigned long> count{0};
	std::mutex m;

	auto PEList = std::make_shared<std::vector<EventData>>();


	// If only one thread is being used there's no need to use std::thread, just call batchRun conventionally.
	if (numThreads == 1) {
		fitBatchPEs(data.getWaveforms(), count, m, PEList);
		return 1;
	}

	// Determining how many rays each thread should simulate.
	unsigned int threadRepeatCount[numThreads];
	unsigned int threadsWithExtra = data.getWaveforms().size() % numThreads;
	unsigned int minRepeatsPerThread = data.getWaveforms().size() / numThreads;

	for (unsigned int i = 0; i < numThreads; i++) {
		if (i < threadsWithExtra) {
			threadRepeatCount[i] = minRepeatsPerThread + 1;
		} else {
			threadRepeatCount[i] = minRepeatsPerThread;
		}
	}

	std::vector<std::thread> threads;
	// Carrying out the multi-threaded simulations.
	for (unsigned int i = 0; i < numThreads; i++) {
		std::vector passData = slice(data.getWaveforms(), threadRepeatCount[i], threadRepeatCount[i+1]);
		std::thread t(fitBatchPEs, passData, std::reference_wrapper(count), std::reference_wrapper(m), std::reference_wrapper(PEList));
		threads.push_back(std::move(t));
	}

	// Waiting for all threads to complete before continuing code execution.
	for (auto &th : threads) {
		th.join();
	}

	return 1;
}
