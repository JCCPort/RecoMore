#ifndef RECOMORE_UTILS_H
#define RECOMORE_UTILS_H
#include <list>
#include <cmath>
#include <c++/9/stdexcept>

template<typename T>
void valueChecker(const std::list<T> list){
	for(auto entry: list){
		if(std::isnan(entry)){
			throw std::runtime_error("Number should not be NaN");
		}
		else if(std::isinf(entry)){
			throw std::runtime_error("Number should not be inf");
		}
	}
}

#endif //RECOMORE_UTILS_H
