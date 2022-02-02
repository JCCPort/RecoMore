#ifndef RECOMORE_UTILS_H
#define RECOMORE_UTILS_H

#include <iostream>
#include <atomic>
#include <mutex>
#include <thread>
#include "Globals.h"


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

#endif //RECOMORE_UTILS_H

