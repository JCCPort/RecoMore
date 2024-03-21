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
#include "DataStructures.h"
#include "DataReading.h"
#include "../Globals.h"
#include <fstream>
#include <sstream>

struct comparePETime {
	inline bool operator()(const Photoelectron &PE1, const Photoelectron &PE2) {
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
[[maybe_unused]] float averageVector(std::vector<T> vec, int start, int end, float cut) {
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
float averageVector(std::vector<T> vec, float cut) {
	float sum = 0;
	for (auto val: vec) {
		float temp = val;
		if (std::abs(temp) < cut) {
			sum += val;
		}
	}
	return sum / vec.size();
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

template<typename T>
double calculateMean(const std::vector<T>& numbers) {
    T sum = T(); // Initializes to zero (for numeric types)
    for (const auto& num : numbers) {
        sum += num;
    }
    return sum / static_cast<double>(numbers.size());
}

template<typename T>
double calculateVariance(const std::vector<T>& numbers, double mean) {
    double variance = 0.0;
    for (const auto& num : numbers) {
        variance += std::pow(num - mean, 2);
    }
    return variance / numbers.size();
}

template<typename T>
double calculateStandardDeviation(const std::vector<T>& numbers) {
    double mean = calculateMean(numbers);
    double variance = calculateVariance(numbers, mean);
    return std::sqrt(variance);
}

#endif //RECOMORE_UTILS_H
