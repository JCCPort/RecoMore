#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "../../../include/DataStructures.h"
#include "../../../include/DataReading.h"

namespace py = pybind11;

PYBIND11_MODULE(_core, m) {
	m.doc() = "Wrapped C++ readers for data files.";
	py::class_<PEData>(m, "PE", py::dynamic_attr())
			.def(py::init<const float&, const float&, const float&, const float&, const float&, const float&>())
			.def_readwrite("amplitude", &PEData::amplitude)
			.def_readwrite("amplitudeError", &PEData::amplitudeError)
			.def_readwrite("time", &PEData::time)
			.def_readwrite("timeError", &PEData::timeError)
			.def_readwrite("foundAmplitude", &PEData::foundAmplitude)
			.def_readwrite("foundTime", &PEData::foundTime);
	py::class_<ChannelFitData>(m, "ChannelFitData", py::dynamic_attr())
	        .def(py::init<const unsigned short&, const float&, const float&, const std::vector<PEData>&>())
			.def_readwrite("ch", &ChannelFitData::ch)
			.def_readwrite("redChiSq", &ChannelFitData::redChiSq)
			.def_readwrite("baseline", &ChannelFitData::baseline)
			.def_readwrite("pes", &ChannelFitData::pes);
	py::class_<EventFitData>(m, "EventFitData", py::dynamic_attr())
	        .def(py::init<const unsigned int&, const std::string&, const std::string&, const std::vector<ChannelFitData>>())
			.def_readwrite("eventID", &EventFitData::eventID)
			.def_readwrite("TDCCorrTime", &EventFitData::TDCCorrTime)
			.def_readwrite("date", &EventFitData::date)
			.def_readwrite("SiPM", &EventFitData::SiPM);
	py::class_<ChannelData>(m, "ChannelData", py::dynamic_attr())
	        .def(py::init<const unsigned short&, const std::vector<float>&>())
			.def_readwrite("channel", &ChannelData::channel)
			.def_readwrite("waveform", &ChannelData::waveform);
	py::class_<EventData>(m, "EventData", py::dynamic_attr())
			.def(py::init<const unsigned int&, const std::string&, const std::string&, const std::vector<ChannelData>>())
			.def_readwrite("eventID", &EventData::eventID)
			.def_readwrite("TDCCorrTime", &EventData::TDCCorrTime)
			.def_readwrite("date", &EventData::date)
			.def_readwrite("chData", &EventData::chData);
	
	py::class_<WCData>(m, "WCData")
	        .def("addRow", &WCData::addRow, "Add entry row to WCData")
			.def("getEvents", &WCData::getEvents, "Get rows");
	
	m.def("ReadWCDataFileDat", &ReadWCDataFileDat, "Read plain text WaveCatcher data files.");
	m.def("ReadWCDataFileBinary", &ReadWCDataFileBinary, "Read binary WaveCatcher data files.");
}
