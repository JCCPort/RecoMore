#include <vector>
#include "utils.h"

unsigned int binarySearch(std::vector<float> vec, float x, int low, int high) {
	int mid;
	while (low < high) {
		mid = (high + low) / 2;
		if (vec[mid] == x) {
			break;
		} else if (vec[mid] > x) {
			high = mid - 1;
		} else {
			low = mid + 1;
		}
	}
	mid = (high + low) / 2;
	if (x <= vec[mid])
		return mid;
	else
		return mid + 1;
}