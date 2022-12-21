#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "../include/DataStructures.h"
#include "../include/DataReading.h"

namespace py = pybind11;

PYBIND11_MODULE(Bindings, m) {
	m.doc() = "Wrapped C++ readers for data files.";
	py::class_<PEData>(m, "PE")
			.def(py::init<const float&, const float&, const float&, const float&, const float&, const float&>())
			.def_readwrite("amplitude", &PEData::amplitude)
			.def_readwrite("amplitudeError", &PEData::amplitudeError)
			.def_readwrite("time", &PEData::time)
			.def_readwrite("timeError", &PEData::timeError)
			.def_readwrite("foundAmplitude", &PEData::foundAmplitude)
			.def_readwrite("foundTime", &PEData::foundTime);
	py::class_<ChannelFitData>(m, "ChannelFitData")
	        .def(py::init<const unsigned short&, const float&, const float&, const std::vector<PEData>&>())
			.def_readwrite("ch", &ChannelFitData::ch)
			.def_readwrite("redChiSq", &ChannelFitData::redChiSq)
			.def_readwrite("baseline", &ChannelFitData::baseline)
			.def_readwrite("pes", &ChannelFitData::pes);
	py::class_<EventFitData>(m, "EventFitData")
	        .def(py::init<const unsigned int&, const std::string&, const std::string&, const std::vector<ChannelFitData>>())
			.def_readwrite("eventID", &EventFitData::eventID)
			.def_readwrite("TDCCorrTime", &EventFitData::TDCCorrTime)
			.def_readwrite("date", &EventFitData::date)
			.def_readwrite("SiPM", &EventFitData::SiPM);
	py::class_<WaveformData>(m, "WaveformData")
	        .def(py::init<const unsigned short&, const std::vector<float>&>())
			.def_readwrite("channel", &WaveformData::channel)
			.def_readwrite("waveform", &WaveformData::waveform);
	py::class_<EventData>(m, "EventData")
			.def(py::init<const unsigned int&, const std::string&, const std::string&, const std::vector<WaveformData>>())
			.def_readwrite("eventID", &EventData::eventID)
			.def_readwrite("TDCCorrTime", &EventData::TDCCorrTime)
			.def_readwrite("date", &EventData::date)
			.def_readwrite("chData", &EventData::chData);
	
	py::class_<WCData>(m, "WCData")
	        .def("addRow", &WCData::addRow, "Add entry row to WCData")
			.def("getEvents", &WCData::getEvents, "Get rows");
	
	m.def("ReadWCDataFileDat", &ReadWCDataFileDat, "Read plain text WaveCatcher data files.");
	m.def("ReadWCDataFileBinary", &ReadWCDataFileBinary, "Read binary WaveCatcher data files.");

#ifdef VERSION_INFO
	m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
	m.attr("__version__") = "dev";
#endif
}
