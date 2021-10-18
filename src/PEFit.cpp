#include "../include/PEFit.h"
#include "../Globals.h"
#include <cmath>
#include "ceres/ceres.h"
#include "ceres/cubic_interpolation.h"
#include <memory>
#include <utility>


struct npe_pdf_functor {
	npe_pdf_functor(double x, double y, ceres::CubicInterpolator<ceres::Grid1D<double, true>> *PDFInterpolator, int numPES) : x_(
			x), y_(y), PDFInterpolator_(PDFInterpolator), numPES_(numPES) {};

	template<typename T>
	bool operator()(T const* const* params, T *residual) const {
		// TODO(josh): Modify so that starting value of residual[0] = baseline
		T f;
		T X_(x_);
		residual[0] = params[0][0];
		for(int i = 0; i < numPES_; i++){
			auto evalVal = ((double) pdfT0Sample + ((X_ - params[(2*i)+2][0]) / (0.01 * pdfSamplingRate)));
			PDFInterpolator_->Evaluate(evalVal, &f);
			auto thing = (params[(2*i)+1][0] * f);
			residual[0] -= thing;
		}
		residual[0] += y_;
		return true;
	}

private:
	const double x_;
	const double y_;
	ceres::CubicInterpolator<ceres::Grid1D<double> > *PDFInterpolator_;
	int numPES_;
};


double npe_pdf_func(double X, const std::vector<double> &p, std::vector<double> *idealWaveform) {
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
		// TODO(josh): Start using interpolation for this too instead of using floor
		const int PE_PDF_BIN =
				pdfT0Sample + (int) std::floor(0.5 + (X - PE_TIME) / (0.01 * pdfSamplingRate));
		if ((PE_PDF_BIN >= 0) && (PE_PDF_BIN < pdfNSamples)) {
			double thing = idealWaveform->at(PE_PDF_BIN);
			value += PE_CHARGE * thing;
		}
	}

	return value;
}


void
fitPE(const EventData &event, const std::vector<std::vector<double>> &idealWaveforms, std::shared_ptr<SyncFile> file) {
	EventFitData evFitDat;
	evFitDat.eventID = event.eventID;
	for (const auto waveformData: event.chData) {
		auto residualWaveform = waveformData; // This will be the variable that is modified to be the residual distribution after each iteration

		ChannelFitData chFit;
		unsigned int ch = waveformData.channel;
		chFit.ch = ch;

		// Making a pointer to the ideal waveform for this channel to improve speed of passing.
		std::vector<double> tmp = idealWaveforms[ch];
		std::vector<double> *chIdealWaveform = &tmp;

		// Baseline calculation
		double baseline = 0;
		chFit.baseline = baseline; // Will want to replace this with the fit baseline

		// Start loop that will break when no more PEs are present
		std::vector<PEData> pesFound;
		unsigned int numPEsFound = 0;

		// This is getting estimate PEs that will then be passed as initial guesses to the minimiser.
		PEData guessPE;
		while (true) {

			// Making the vector of parameters for the objective function, not present in ROOT version as it uses a static object for the params,
			//  which is worth considering
			// TODO(josh): Make a function that does this.
			std::vector<double> params;
			params.push_back((double) pesFound.size());
			params.push_back(baseline);
			for (auto pe: pesFound) {
				params.push_back(pe.amplitude);
				params.push_back(pe.time);
			}

			// Amplitude adjustment: if the latest PE found is before other one(s),
			//  its tail is going to add some amplitude to the following one.
			for (auto &pe: pesFound) {
				unsigned int peTimeBinPos = std::floor(pe.time / pdfSamplingRate); // This should use a variable
				// corresponding to the input sampling rate, however for now they're the same
				double fitVal = npe_pdf_func(pe.time, params, chIdealWaveform);
				double extraAmplitude = fitVal - residualWaveform.waveform[peTimeBinPos];
				double newAmplitude = pe.amplitude + extraAmplitude;
				pe.amplitude = newAmplitude;

			}

			std::ofstream myfile;
			myfile.open("rawWaveform.csv");
			for (float k: waveformData.waveform) {
				myfile << k << "\n";
			}
			myfile.close();


			// Compute residual
			std::ofstream myfile2;
			myfile2.open("fit.csv");
			for (unsigned int k = 0; k < residualWaveform.waveform.size(); ++k) {
				double val = npe_pdf_func((k + 0.5) * pdfSamplingRate, params, chIdealWaveform);
				myfile2 << val << "\n";
				residualWaveform.waveform[k] = residualWaveform.waveform[k] - val;
			}
			myfile2.close();


			std::ofstream myfile3;
			myfile3.open("residual.csv");
			for (float k: residualWaveform.waveform) {
				myfile3 << k << "\n";
			}
			myfile3.close();



			// Get initial guesses for the next PE
			auto minPosIt = std::min_element(residualWaveform.waveform.begin(), residualWaveform.waveform.end());
			unsigned int minTimePos = std::distance(residualWaveform.waveform.begin(), minPosIt);

			if (-residualWaveform.waveform[minTimePos] < 0.013) {
				break;
			}

			guessPE.amplitude = -residualWaveform.waveform[minTimePos];
			guessPE.time = minTimePos * pdfSamplingRate;

			if ((minTimePos > 1) && (minTimePos < 1024)) {
				// improve initial time for a new PE based on average time
				// over 3 consecutive sample ponderated by the amplitude
				// of each sample... help a lot to resolve PEs very closed!
				double timeSum = 0;
				double ponderationSum = 0;
				for (unsigned int b = minTimePos - 1; b <= minTimePos + 1; ++b) {
					double binCenter = (minTimePos - 0.5) * (pdfSamplingRate);
					double binVal = residualWaveform.waveform[minTimePos];
					timeSum += binCenter * binVal;
					ponderationSum += binVal;
				}

				guessPE.time = (PEFinderTimeOffset * 0.10) + timeSum / ponderationSum;
			}
			numPEsFound += 1;
			pesFound.push_back(guessPE);
			residualWaveform = waveformData;

			if (numPEsFound > 100) {
				break;
			}
		}

		std::vector<double> xValues;
		for (unsigned int j = 0; j <= 1024; j++) {
			xValues.push_back((double) j * pdfSamplingRate);
		}


		using namespace ceres;
		// Build the problem.
		Problem problem;

		std::vector<double> initialTimes;
		std::vector<double> initialAmplitudes;
		for (auto &k: pesFound) {
			initialAmplitudes.push_back(k.amplitude);
			initialTimes.push_back(k.time);
		}

		std::vector<double> times = initialTimes;
		std::vector<double> amplitudes = initialAmplitudes;


		ceres::Grid1D grid = ceres::Grid1D<double>(chIdealWaveform->data(), 0, (int) chIdealWaveform->size());
		auto idealPDFInterpolator = new ceres::CubicInterpolator<ceres::Grid1D<double> >(grid);

		// Set up the only cost function (also known as residual). This uses
		// auto-differentiation to obtain the derivative (Jacobian).

		double baseline2 = 0;
		std::vector<double*> params;
		params.push_back(&baseline2);
		for (int i = 0; i < pesFound.size(); i++) {
			params.push_back(&amplitudes[i]);
			params.push_back(&times[i]);
		}

		for (unsigned int j = 0; j < waveformData.waveform.size(); ++j) {
//			auto costFunction = new AutoDiffCostFunction<npe_pdf_functor, 1, 1, 1, 1>(new npe_pdf_functor(xValues[j],
//			                                                                                              waveformData.waveform[j],
//			                                                                                              idealPDFInterpolator));
			auto costFunction = new DynamicAutoDiffCostFunction<npe_pdf_functor, 4>(new npe_pdf_functor(xValues[j],
			                                                                                              waveformData.waveform[j],
			                                                                                              idealPDFInterpolator,
			                                                                                            pesFound.size()));

			costFunction->AddParameterBlock(1);
			for (auto pe: pesFound){
				costFunction->AddParameterBlock(1);
				costFunction->AddParameterBlock(1);
			}
			costFunction->SetNumResiduals(1);

			problem.AddResidualBlock(costFunction, nullptr, params);

//			for (int k = 0; k < pesFound.size(); k++) {
//				auto costFunction = new AutoDiffCostFunction<npe_pdf_functor, 1, 1, 1, 1>(new npe_pdf_functor(xValues[j],
//				                                                                                              waveformData.waveform[j],
//				                                                                                              idealPDFInterpolator));
//				problem.AddResidualBlock(costFunction, nullptr, &times[k], &amplitudes[k], &baseline2);
//			}
			// TODO(josh): The minimum amplitude for the actual data is always 63-64 entries ahead of the minimum amplitude for the ideal waveform.
			//  This is likely relevant to the length of the waveform I'm using being 960 entries, but the ideal waveform being made for 1024 entries
		}

		for (int k = 0; k < pesFound.size(); k++) {
//			problem.SetParameterLowerBound(&times[k], 0, initialTimes[k] * 0.95);
//			problem.SetParameterLowerBound(&amplitudes[k], 0, initialAmplitudes[k] * 1.8);
//			problem.SetParameterUpperBound(&times[k], 0, initialTimes[k] * 1.05);
//			problem.SetParameterUpperBound(&amplitudes[k], 0, initialAmplitudes[k] * 0.3);
		}
//		problem.SetParameterLowerBound(params[0], 0, -0.001);
//		problem.SetParameterLowerBound(params[0], 0, 0.001);

		// Run the solver!
		Solver::Options options;
//		options.function_tolerance = 1e-12;
//		options.gradient_tolerance = 1e-12;
//		options.parameter_tolerance = 1e-12;
//		options.linear_solver_type = ceres::DENSE_QR;
//		options.linear_solver_type = ceres::DENSE_SCHUR;
//		options.linear_solver_type = ceres::SPARSE_NORMAL_CHOLESKY;
//		options.trust_region_strategy_type = ceres::LEVENBERG_MARQUARDT;
		options.minimizer_progress_to_stdout = true;
		Solver::Summary summary;
		Solve(options, &problem, &summary);

		delete idealPDFInterpolator;

		std::cout << summary.FullReport() << "\n";



		std::vector<double> params2;
		params2.push_back(pesFound.size());
		params2.push_back(baseline);
		for (int k = 0; k < pesFound.size(); k++) {
			params2.push_back(amplitudes[k]);
			params2.push_back(times[k]);
		}

		std::ofstream myfile4;
		myfile4.open("fullFit.csv");
		for (unsigned int k = 0; k < waveformData.waveform.size(); k++) {
			auto thing = npe_pdf_func(k * pdfSamplingRate, params2, chIdealWaveform);
			myfile4 << thing << "\n";
		}
		myfile4.close();

//		for (int k = 0; k < pesFound.size(); k++) {
//			std::cout << "Amplitude:\t" << initialAmplitudes[k] << " -> " << amplitudes[k] << std::endl;
//			std::cout << "Time:\t\t" << initialTimes[k] << " -> " << times[k] << std::endl;
//			std::cout << std::endl;
//		}

		std::vector<PEData> FitPEs;

		for (int k = 0; k < pesFound.size(); k++) {
			PEData pe;
			pe.amplitude = amplitudes[k];
			pe.time = times[k];
			FitPEs.emplace_back(pe);
		}

		chFit.pes = FitPEs;
		evFitDat.sipm.push_back(chFit);


//		std::cout << "hey ho" << std::endl;
		// TODO(josh): Check if chisq is better with the fit, if not, keep initial params.
		//  BUT THIS IS VERY HACKY AND WE SHOULD FIGURE OUT WHY IT ISN'T WORKING

		// TODO(josh): No magic numbers for different lengths of data

		// TODO(josh): It's worsening the fit when there are two PEs very close to each other

		// TODO(josh): We must consider after pulses giving a weird shape?

		// Note: Bounds slow the fit, buuut when running in release mode it still is very fast
		// Note: DENSE_QR is faster than LEVENBERG_MARQUARDT

		// The fundamental question is: why does the fit seem to give worse parameters than the initial guesses for some fits?!?!?
	}
	Writer writer(std::move(file));
	writer.writeEventInfo(evFitDat);

}