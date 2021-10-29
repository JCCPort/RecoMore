#ifndef RECOMORE_UTILS_H
#define RECOMORE_UTILS_H
#include <list>
#include <cmath>
#include <exception>
#include "include/DataStructures.h"

template<typename T>
void valueChecker(const std::list<T> list){
//	for(auto entry: list){
////	    assert(std::isnan(entry));
////	    assert(std::isinf(entry));
//		if(std::isnan(entry)){
//			throw std::runtime_error("Number should not be NaN");
//		}
//		else if(std::isinf(entry)){
//			throw std::runtime_error("Number should not be inf");
//		}
//	}
}

struct compareTime
{
	inline bool operator() (const PEData& PE1, const PEData& PE2)
	{
		return (PE1.time < PE2.time);
	}
};

#endif //RECOMORE_UTILS_H
