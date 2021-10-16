#include "../include/PEFit.h"
#include "../Globals.h"
#include <cmath>
#include <fstream>
#include "ceres/ceres.h"
#include "ceres/cubic_interpolation.h"
#include "ceres/loss_function.h"
#include "glog/logging.h"
#include <eigen3/Eigen/Core>
#include <memory>
#include <utility>
#include <memory>

double samplingRate2 = 0.01 * pdfSamplingRate;

struct npe_pdf_functor {
	npe_pdf_functor(double x, double y, const ceres::CubicInterpolator<ceres::Grid1D<double, true>>& compute_distortion) : x_(x), y_(y), compute_distortion_(compute_distortion) {};

	template<typename T>
	bool operator()(const T *const time, const T *const charge, const T *const baseline, T *residual) const {
			// TODO(josh): Modify so that starting value of residual[0] = baseline
			T f;
			auto evalVal = ((double)pdfT0Sample + ((x_ - time[0]) / (samplingRate2)));
			compute_distortion_.Evaluate(evalVal, &f);
//			if(evalVal <= (double)0){
//				T k;
//				compute_distortion_.Evaluate(-((double)pdfT0Sample + (0.5 + (x_ - time[0]) / (samplingRate2)) - (double)100000), &f);
//				residual[0] = k;
//				// TODO(josh): This isn't correct it's just I need to set residual to some value before
//				//  exiting the function, I believe it should be set to baseline - y_  or something like that...
//				return true;
//			}
//			if((x_ >= ))
//			std::cout << "x val:\t" << x_ << std::endl;
//			std::cout << "PE test time:\t" << time[0] << std::endl;
//			std::cout << "x Point for interpolation:\t" << evalVal << std::endl;
//			residual[0] = -(charge[0] * f);
//			std::cout << "Interpolated y value:\t" << residual[0] << std::endl;
//			std::cout << "Data y value:\t" << y_ << std::endl;
//		}

		auto thing = (charge[0] * f) - baseline[0];
		residual[0] = y_ - thing;
		return true;
	}

private:
	const double x_;
	const double y_;
	const ceres::CubicInterpolator<ceres::Grid1D<double> >&  compute_distortion_;
};

// y_obs,i -  (baseline + sum_i=0^NPE y_func,i)


double npe_pdf_func(double X, const std::vector<double> &p, std::vector<double> idealWaveform) {
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
			double thing = idealWaveform[PE_PDF_BIN];
			value += PE_CHARGE * thing;
		}
	}

	return value;
}


void fitPE(const EventData &event, const std::shared_ptr<std::vector<EventFitData>> &FitList,
           const std::vector<std::vector<double>> &idealWaveforms) {
	EventFitData evFitDat;
	evFitDat.eventID = event.eventID;
	for (unsigned int i = 0; i < event.chData.size(); ++i) {
		auto channelWaveform = event.chData[i]; // This will be the variable that is modified to be the residual distribution after each iteration

		ChannelFitData chFit;
		unsigned int ch = channelWaveform.channel;
		chFit.ch = ch;

		// Baseline calculation
		double baseline = 0;
		chFit.baseline = baseline; // Will want to replace this with the fit baseline

		// Start loop that will break when no more PEs are present
		std::vector<PEData> pesFound;
		unsigned int numPEsFound = 0;

		// This is getting estimate PEs that will then be passed as initial guesses to the minimiser.
		while (true) {


			// Making the vector of parameters for the objective function, not present in ROOT version as it uses a static object for the params,
			//  which is worth considering
			// TODO(josh): Make a function that does this.
			std::vector<double> params;
			params.push_back((double)pesFound.size());
			params.push_back(baseline);
			for (auto pe2: pesFound) {
				params.push_back(pe2.amplitude);
				params.push_back(pe2.time);
			}
			// Amplitude adjustment: if the latest PE found is before other one(s),
			//  its tail is going to add some amplitude to the following one.
			for (unsigned int j = 0; j < pesFound.size(); ++j) {
				unsigned int peTimeBinPos = std::floor(
						pesFound[j].time / pdfSamplingRate); // This should use a variable
				// corresponding to the input sampling rate, however for not they're the same
				double fitVal = npe_pdf_func(pesFound[j].time, params, idealWaveforms[ch]);
				double extraAmplitude =
						fitVal - channelWaveform.waveform[peTimeBinPos];
				double newAmplitude = pesFound[j].amplitude + extraAmplitude;
				pesFound[j].amplitude = newAmplitude;

			}

			std::ofstream myfile;
			myfile.open("rawWaveform.csv");
			for (float k: channelWaveform.waveform) {
				myfile << k << "\n";
			}
			myfile.close();


			// Compute residual
			std::ofstream myfile2;
			myfile2.open("fit.csv");
			for (unsigned int k = 1; k <= channelWaveform.waveform.size(); ++k) {
				double val = npe_pdf_func((k - 0.5 + 1) * pdfSamplingRate, params, idealWaveforms[ch]);
				myfile2 << val << "\n";
				channelWaveform.waveform[k] = channelWaveform.waveform[k] - val;
			}
			myfile2.close();


			std::ofstream myfile3;
			myfile3.open("residual.csv");
			for (float k: channelWaveform.waveform) {
				myfile3 << k << "\n";
			}
			myfile3.close();



			// Get initial guesses for the next PE
			auto minPosIt = std::min_element(channelWaveform.waveform.begin(), channelWaveform.waveform.end());
			unsigned int minTimePos = std::distance(channelWaveform.waveform.begin(), minPosIt);

			if (-channelWaveform.waveform[minTimePos] < 0.015) {
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

				guessPE.time = (PEFinderTimeOffset * 0.25) + time_sum / ponderation_sum;
			}

			pesFound.push_back(guessPE);
			channelWaveform = event.chData[i];

		}

		channelWaveform = event.chData[i];

		std::vector<double> params;
		params.push_back(baseline);
		for (auto pe2: pesFound) {
			params.push_back(pe2.amplitude);
			params.push_back(pe2.time);
		}

		std::vector<float> xValues;
		for (unsigned int j = 0; j <= 1024; j++) {
			xValues.push_back(j * pdfSamplingRate);
		}


		using namespace ceres;
		// Build the problem.
		Problem problem;

		std::vector<double> initialTimes;
		std::vector<double> initialAmplitudes;
		for (auto & k : pesFound){
			initialAmplitudes.push_back(k.amplitude);
			initialTimes.push_back(k.time);
		}

		std::vector<double> times = initialTimes;
		std::vector<double> amplitudes = initialAmplitudes;


		ceres::Grid1D grid = ceres::Grid1D<double>(idealWaveforms[ch].data(), 0, idealWaveforms[ch].size());
		auto compute_distortion = ceres::CubicInterpolator<ceres::Grid1D<double> >(grid);

		// Set up the only cost function (also known as residual). This uses
		// auto-differentiation to obtain the derivative (Jacobian).
		for (unsigned int j = 0; j < channelWaveform.waveform.size(); ++j) {
			CostFunction *cost_function =
					new AutoDiffCostFunction<npe_pdf_functor, 1, 1, 1, 1>(new npe_pdf_functor(xValues[j],
					                                                                       channelWaveform.waveform[j],
					                                                                       compute_distortion));


			double sumThing = 0;
			for (int k = 0; k < pesFound.size(); k++) {

				double f2;
				auto evalVal = ((double)pdfT0Sample + ((xValues[j] - times[k]) / (samplingRate2)));
				compute_distortion.Evaluate(evalVal, &f2);
				auto thing = (amplitudes[k] * f2);
				sumThing += thing;

				double baseline2 = 0;
				problem.AddResidualBlock(cost_function, nullptr, &times[k], &amplitudes[k], &baseline2);
				problem.SetParameterLowerBound(&times[k], 0,  times[k]*0.96);
				problem.SetParameterLowerBound(&amplitudes[k], 0, amplitudes[k]*0.96);
				problem.SetParameterUpperBound(&times[k], 0,  times[k]*1.04);
				problem.SetParameterUpperBound(&amplitudes[k], 0, amplitudes[k]*1.04);
			}
			// TODO(josh): The minimum amplitude for the actual data is always 63-64 entries ahead of the minimum amplitude for the ideal waveform.
			//  This is likely relevant to the length of the waveform I'm using being 960 entries, but the ideal waveform being made for 1024 entries
//			std::cout << j << "\t" << "Actual waveform height:\t" << channelWaveform.waveform[j] << "\t\tFunctor calc waveform height:\t" << sumThing << std::endl;
//			std::cout << "hello" << std::endl;
		}


		// Run the solver!
		Solver::Options options;
//		options.function_tolerance = 1e-12;
//		options.gradient_tolerance = 1e-12;
//		options.parameter_tolerance = 1e-12;
		options.linear_solver_type = ceres::DENSE_QR;
//		options.linear_solver_type = ceres::DENSE_SCHUR;
//		options.linear_solver_type = ceres::SPARSE_NORMAL_CHOLESKY;
//		options.trust_region_strategy_type = ceres::LEVENBERG_MARQUARDT;
		options.minimizer_progress_to_stdout = false;
		Solver::Summary summary;
		Solve(options, &problem, &summary);

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
		for (unsigned int k = 0; k < channelWaveform.waveform.size(); k++) {
			auto thing = npe_pdf_func(k * pdfSamplingRate, params2, idealWaveforms[ch]);
			myfile4 << thing << "\n";
		}
		myfile4.close();

		for (int k = 0; k < pesFound.size(); k++) {
			std::cout << "Amplitude:\t" << initialAmplitudes[k] << " -> " << amplitudes[k] << std::endl;
			std::cout << "Time:\t\t" << initialTimes[k] << " -> " << times[k] << std::endl;
			std::cout << std::endl;
		}

		std::vector<PEData> FitPEs;

		for (int k = 0; k < pesFound.size(); k++) {
			PEData pe;
			pe.amplitude = amplitudes[k];
			pe.time = times[k];
			FitPEs.emplace_back(pe);
		}

		chFit.pes = FitPEs;
		evFitDat.sipm.push_back(chFit);


		std::cout << "hey ho" << std::endl;
		// TODO(josh): Check if chisq is better with the fit, if not, keep initial params.
		//  BUT THIS IS VERY HACKY AND WE SHOULD FIGURE OUT WHY IT ISN'T WORKING

		// TODO(josh): No magic numbers for different lengths of data

		// TODO(josh): It's worsening the fit when there are two PEs very close to each other

		// TODO(josh): We must consider after pulses giving a weird shape?

		// Note: Bounds slow the fit, buuut when running in release mode it still is very fast
		// Note: DENSE_QR is faster than LEVENBERG_MARQUARDT

		// The fundamental question is: why does the fit seem to give worse parameters than the initial guesses for some fits?!?!?
	}
	FitList->push_back(evFitDat);

}