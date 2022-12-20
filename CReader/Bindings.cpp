#include <pybind11/pybind11.h>
#include "../include/DataStructures.h"
#include "../include/DataReading.h"

namespace py = pybind11;


PYBIND11_MODULE(CReader, m) {
	py::class_<PEData>(m, "PE")
			.def(py::init<const float&, const float&, const float&, const float&, const float&, const float&>());
	py::class_<ChannelFitData>(m, "ChannelFitData")
	        .def(py::init<const unsigned short&, const float&, const float&, const std::vector<PEData>&>());
	py::class_<EventFitData>(m, "EventFitData")
	        .def(py::init<const unsigned int&, const std::string&, const std::string&, const std::vector<const ChannelFitData>>());
	py::class_<WaveformData>(m, "WaveformData")
	        .def(py::init<const unsigned short&, const std::vector<float>&>());
	py::class_<EventData>(m, "EventData")
			.def(py::init<const unsigned int&, const std::string&, const std::string&, const std::vector<WaveformData>>());
	
	py::class_<WCData>(m, "WCData")
	        .def("addRow", &WCData::addRow, "Add entry row to WCData")
			.def("getEvents", &WCData::getEvents, "Get rows");
	
	m.def("ReadWCDataFileDat", &ReadWCDataFileDat, "Read plain text WaveCatcher data files.");
	m.def("ReadWCDataFileBinary", &ReadWCDataFileBinary, "Read binary WaveCatcher data files.");
}
