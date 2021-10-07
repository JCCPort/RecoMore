#include <iostream>
#include "include/DataReading.h"
#include "include/DataStructures.h"

int main() {
	WCData data = ReadWCDataFile("/home/josh/CLionProjects/RecoMore/R43.dat");
	std::cout << "Hello, World!" << std::endl;
	return 0;
}
