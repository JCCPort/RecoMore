#ifndef RECOMORE_UTILS_H
#define RECOMORE_UTILS_H

#include <list>
#include <cmath>
#include <exception>
#include <iostream>
#include <atomic>
#include <mutex>
#include <thread>
#include <boost/algorithm/string.hpp>
#include "include/DataStructures.h"
#include "include/DataReading.h"
#include "Globals.h"
#include <fstream>
#include <sstream>

struct comparePETime {
	inline bool operator()(const PEData &PE1, const PEData &PE2) {
		return (PE1.time < PE2.time);
	}
};

template<typename T>
[[maybe_unused]] void nanInfChecker(const std::list<T> list) {
	for (auto entry: list) {
		if (std::isnan(entry)) {
			throw std::runtime_error("Number should not be NaN.");
		} else if (std::isinf(entry)) {
			throw std::runtime_error("Number should not be inf.");
		}
	}
}


template<typename T>
float averageVector(std::vector<T> vec, int start, int end, float cut) {
	float sum = 0;
	for (int i = start; i < end; i++) {
		float temp = vec[i];
		if (std::abs(temp) < cut) {
			sum += vec[i];
		}
	}
	return sum / ((float) end - (float) start);
}

template<typename T>
void writeVector(const std::string &fileName, std::vector<T> vector) {
	std::ofstream file;
	file.open(fileName, std::ofstream::trunc);
	for (float k: vector) {
		file << k << "\n";
	}
	file.close();
}

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


void displayProgress(std::atomic<unsigned long> &count, std::mutex &m, unsigned int dataLength);

std::string defaultOutputName(std::string inputName);

template <typename T>
std::string to_string_with_precision(const T a_value, const int n) {
	std::ostringstream out;
	out.precision(n);
	out << std::fixed << a_value;
	return out.str();
}
#endif //RECOMORE_UTILS_H
