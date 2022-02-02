#ifndef RECOMORE_UTILS_H
#define RECOMORE_UTILS_H

#include <list>
#include <cmath>
#include <exception>
#include "include/DataStructures.h"

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

struct comparePETime {
	inline bool operator()(const PEData &PE1, const PEData &PE2) {
		return (PE1.time < PE2.time);
	}
};

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

#endif //RECOMORE_UTILS_H
