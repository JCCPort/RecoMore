#include "../include/PEFit.h"
#include "../Globals.h"
#include <cmath>
#include <fstream>
#include "ceres/ceres.h"
#include "glog/logging.h"

struct npe_pdf_functor {
	npe_pdf_functor(double x, double y) : x_(x), y_(y) {};
	template <typename T>
	bool operator()(const T* &p, const T* idealWaveform, T* residual) const {
		const int NPE = (int) (p[0]);
		const double BASELINE = p[1];

		double value = BASELINE;

		// This is adding up the contribution of each fit PE to a specific bin
		for (int PE = 0; PE < NPE; ++PE) {
			const double PE_CHARGE = p[2 + PE * 2];
			const double PE_TIME = p[3 + PE * 2];
			const int PE_PDF_BIN =
					pdfT0Sample + std::floor(0.5 + (x_ - PE_TIME) / (0.01 * pdfSamplingRate)); // v4 use PE_PDF[PE_PDF_CH]
			if ((PE_PDF_BIN >= 0) && (PE_PDF_BIN < pdfNSamples)) {
				float thing = idealWaveform[PE_PDF_BIN];
				value += PE_CHARGE * thing;
			}
		}

		residual[0] = y_ - value;
		return true;
	}
private:
	double x_;
	double y_;
};


double npe_pdf_func(const double X, const std::vector<double> &p, std::vector<float> idealWaveform) {
	// This way of passing x as a list then choosing the 0th index is from ROOT's syntax for fitting where you can fit
	//  with an arbitrary number of dimensions (N input values -> single output, an N dimensional function)
//	const double X = x[0];

	// parameter of the function:
	// p[0] : number N of PE (fixed)
	// p[1] : baseline
	// p[2] : ampl of PE #0
	// p[3] : time of PE #0
	// [...]
	// p[2+N*2] : ampl of PE #N
	// p[3+N*2] : time of PE #N

	const int NPE = (int) (p[0]);
	const double BASELINE = p[1];

	double value = BASELINE;

	// This is adding up the contribution of each fit PE to a specific bin
	for (int PE = 0; PE < NPE; ++PE) {
		const double PE_CHARGE = p[2 + PE * 2];
		const double PE_TIME = p[3 + PE * 2];
		const int PE_PDF_BIN =
				pdfT0Sample + std::floor(0.5 + (X - PE_TIME) / (0.01 * pdfSamplingRate)); // v4 use PE_PDF[PE_PDF_CH]
		if ((PE_PDF_BIN >= 0) && (PE_PDF_BIN < pdfNSamples)) {
			float thing = idealWaveform[PE_PDF_BIN];
			value += PE_CHARGE * thing;
		}
	}

	return value;
}


void fitPE(const EventData &event, const std::shared_ptr<std::vector<EventFitData>> &PEList,
                        const std::vector<std::vector<float>> &idealWaveforms) {
	EventFitData evFitDat;
	evFitDat.eventID = event.eventID;
	for (unsigned int i = 0; i < event.chData.size(); ++i) {
		auto channelWaveform = event.chData[i]; // This will be the variable that is modified to be the residual distribution after each iteration

		ChannelFitData chFit;
		unsigned int ch = channelWaveform.channel;
		chFit.ch = ch;

		// Baseline calculation
		float baseline = 0;

		// Start loop that will break when no more PEs are present
		std::vector<PEData> pesFound;
		unsigned int numPEsFound = 0;

		// This is getting estimate PEs that will then be passed as initial guesses to the minimiser.
		while (true) {


			// Making the vector of parameters for the objective function, not present in ROOT version as it uses a static object for the params,
			//  which is worth considering
			// TODO(josh): Make a function that does this.
			std::vector<double> params;
			params.push_back(pesFound.size());
			params.push_back(baseline);
			for (auto pe2: pesFound) {
				params.push_back(pe2.amplitude);
				params.push_back(pe2.time);
			}
			// Amplitude adjustment: if the latest PE found is before other one(s),
			//  its tail is going to add some amplitude to the following one.
			for (unsigned int j = 0; j < pesFound.size(); ++j) {
				unsigned int peTimeBinPos = std::floor(pesFound[j].time / pdfSamplingRate); // This should use a variable
				// corresponding to the input sampling rate, however for not they're the same
				float fitVal = npe_pdf_func(pesFound[j].time, params, idealWaveforms[ch]);
				float extraAmplitude =
						fitVal - channelWaveform.waveform[peTimeBinPos];
				float newAmplitude = pesFound[j].amplitude + extraAmplitude;
				pesFound[j].amplitude = newAmplitude;

			}

			std::ofstream myfile;
			myfile.open ("rawWaveform.csv");
			for(float k : channelWaveform.waveform){
				myfile << k << "\n";
			}
			myfile.close();


			// Compute residual
			std::ofstream myfile2;
			myfile2.open ("fit.csv");
			for(unsigned int k = 1; k <= channelWaveform.waveform.size(); ++k){
				double val = npe_pdf_func((k - 0.5 + 1) * pdfSamplingRate, params, idealWaveforms[ch]);
				myfile2 << val << "\n";
				channelWaveform.waveform[k] = channelWaveform.waveform[k] - val;
			}
			myfile2.close();


			std::ofstream myfile3;
			myfile3.open ("residual.csv");
			for(float k : channelWaveform.waveform){
				myfile3 << k << "\n";
			}
			myfile3.close();



			// Get initial guesses for the next PE
			auto minPosIt = std::min_element(channelWaveform.waveform.begin(), channelWaveform.waveform.end());
			unsigned int minTimePos = std::distance(channelWaveform.waveform.begin(), minPosIt);

			if(-channelWaveform.waveform[minTimePos] < 0.015){
				break;
			}

			PEData guessPE;
			guessPE.amplitude = -channelWaveform.waveform[minTimePos];
			guessPE.time = minTimePos * pdfSamplingRate;

			if ((minTimePos > 1) && (minTimePos < 1024)) {
				// improve initial time for a new PE based on average time
				// over 3 consecutive sample ponderated by the amplitude
				// of each sample... help a lot to resolve PEs very closed!
				float time_sum = 0;
				float ponderation_sum = 0;
				for (unsigned int b = minTimePos - 1; b <= minTimePos + 1; ++b) {
					float binCenter = (minTimePos - 0.5) * pdfSamplingRate;
					float binVal = channelWaveform.waveform[minTimePos];
					time_sum += binCenter * binVal;
					ponderation_sum += binVal;
				}

				guessPE.time = (PEFinderTimeOffset * 0.15) + time_sum / ponderation_sum;
			}

			pesFound.push_back(guessPE);
			channelWaveform = event.chData[i];

		}

		// The variable to solve for with its initial value.
		double initial_x = 5.0;
		double x = initial_x;

		using namespace ceres;
		// Build the problem.
		Problem problem;

		// Set up the only cost function (also known as residual). This uses
		// auto-differentiation to obtain the derivative (jacobian).
		CostFunction* cost_function =
				new AutoDiffCostFunction<npe_pdf_functor, 1, 1>(new npe_pdf_functor);
		problem.AddResidualBlock(cost_function, nullptr, &x);

		// Run the solver!
		Solver::Options options;
		options.linear_solver_type = ceres::DENSE_QR;
		options.minimizer_progress_to_stdout = true;
		Solver::Summary summary;
		Solve(options, &problem, &summary);


		std::cout << "hey ho" << std::endl;
	}

}