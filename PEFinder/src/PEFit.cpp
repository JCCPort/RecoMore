#include "../include/PEFit.h"
#include "../include/Utils.h"
#include "../include/PETemplate.h"
#include <mutex>
#include <atomic>
#include <cmath>
#include "ceres/ceres.h"
#include "ceres/cubic_interpolation.h"
#include <memory>
#include <utility>


template<typename T>
inline T diffAtTime(T const& baseline, const unsigned int numPEs, T const *const *params, T const& X, ceres::CubicInterpolator<ceres::Grid1D<float> > *PDFInterpolator_) {
	T accum = baseline;
	for (unsigned int i = 0; i < numPEs; i++) {
		const unsigned int i2 = 2 * i;
		const T amplitude = params[i2 + 1][0];
		const T timeShift = params[i2 + 2][0];

		const T pos = X - timeShift;
		if (pos >= T(0) && pos < T(pdfNSamples)) {
			T f;
			PDFInterpolator_->Evaluate(pos, &f);
			accum += amplitude * f;
		}
	}
	return accum;
}


struct NPEPDFFunctor {
	/**
	 * 
	 * @param waveformTimes These need to be the indices of the PETemplate corresponding to the samples of the waveform being processed,
	 * @param waveformAmplitudes 
	 * @param PDFInterpolator 
	 * @param numPES 
	 */
	NPEPDFFunctor(std::vector<float> waveformTimes, std::vector<float> waveformAmplitudes,
	              ceres::CubicInterpolator<ceres::Grid1D<float, true>> *PDFInterpolator, const unsigned int numPES)
			: waveformTimes_(std::move(waveformTimes)), waveformAmplitudes_(std::move(waveformAmplitudes)), PDFInterpolator_(PDFInterpolator), numPES_(numPES) {};

	/**
	 *
	 * @tparam T
	 * @param params Array of parameters to fit. The first two are the amplitude and time of the first PE, the next two are the amplitude and time of the second PE, and so on.
	 * @param residual Sum difference between the fit and the data.
	 * @return
	 */
	template<typename T>
	inline bool operator()(T const *const *params, T *__restrict__ residual) const {
		const T baseline = params[0][0];  // The "DC offset" or baseline
		for (unsigned int j = 0; j < waveformTimes_.size(); j++) {
			residual[j] = diffAtTime(baseline, numPES_, params, T(waveformTimes_[j]), PDFInterpolator_) - T(waveformAmplitudes_[j]);
		}

		return true;
	}

	void updateNumPEs(const unsigned int numPEs) {
		if (numPEs != numPES_) {
			numPES_ = numPEs;
		}
	}

private:
	const std::vector<float>                        waveformTimes_;
	const std::vector<float>                        waveformAmplitudes_;
	ceres::CubicInterpolator<ceres::Grid1D<float> > *PDFInterpolator_;
	unsigned int									numPES_;
};

/**
 * Function that returns the fit y-value of a waveform at a given time (in ns) X, given a set of parameters p and an ideal waveform.
 * 	 Parameters are structured as follows:
 *	 p[0] : number N of PE (fixed)
 *	 p[1] : baseline
 *	 p[2] : amplitude of PE #0
 *	 p[3] : time of PE #0
 *	 [...]
 *	 p[2+N*2] : amplitude of PE #N
 *   p[3+N*2] : time of PE #N
 *
 * @param X Time at which y-value is evaluated (in ns).
 * @param p Vector of parameters containing the number of PE, the baseline, and the amplitude and time of each PE.
 * @param idealWaveform Pointer to the ideal waveform used to fit the data.
 * @return
 */
inline float NPEPDFFunc(const float X, const std::vector<float> &p, const std::vector<double> *idealWaveform) {
	// TODO(josh): Replace p with an instance of fit parameters type
	
	// This way of passing x as a list then choosing the 0th index is from ROOT's syntax for fitting where you can fit
	//  with an arbitrary number of dimensions (N input values -> single output, an N dimensional function)
	const int   NPE      = static_cast<int>(p[0]);
	const float BASELINE = p[1];
	
	float value = BASELINE;
	
	// This is adding up the contribution of each fit PE to a specific bin
	for (int PE = 0; PE < NPE; ++PE) {
		const float PE_CHARGE  = p[2 + PE * 2];
		const float PE_TIME    = p[3 + PE * 2];
		// TODO(josh): Start using interpolation for this too instead of using floor
        const int distCurrPECurrXPos = static_cast<int>(std::floor(0.5 + (X - PE_TIME) * samplingRate2Inv)); // What is the time difference (in terms of bins) between the current PE being considered, and the time bin being considered.
		if (const int PE_PDF_BIN = pdfT0Sample + distCurrPECurrXPos; (PE_PDF_BIN >= 0) && (PE_PDF_BIN < pdfNSamples)) {
			value += PE_CHARGE * static_cast<float>(idealWaveform->at(PE_PDF_BIN));
		}
	}
	return value;
}

inline float NPEPDFFuncCubic(
	const float X,
	const std::vector<float>& p,
	const ceres::CubicInterpolator<ceres::Grid1D<float, true>>* PDFInterpolator,
	const unsigned int pdfNSamples, const PETemplate* ChPETemplate)
{
	// p[0] = NPE, p[1] = baseline, etc.
	const int   NPE      = static_cast<int>(p[0]);
	const float BASELINE = p[1];

	float value = BASELINE;

	for (int PE = 0; PE < NPE; ++PE) {
		const float PE_CHARGE = p[2 + PE * 2];

		const float PE_TIME   = p[3 + PE * 2]  / ChPETemplate->getTimeSpacing();

		// Reproduce your original approach to computing bin/time difference
		// but as a floating-point "pos".
		// Example: floor(0.5 + (X - PE_TIME) * samplingRate2Inv)
		//
		// Instead of forcibly flooring, we let the CubicInterpolator do
		// true interpolation. For a minimal change, we keep the same
		// center shift of +0.5, but that’s optional.
		const auto pos = static_cast<float>(X - PE_TIME);

		// Now we want to interpolate if we're in range
		if (pos >= 0.0f && pos < static_cast<float>(pdfNSamples)) {
			// Call the double overload
			auto posDouble = static_cast<double>(pos);
			double pdfValDouble = 0.0;
			PDFInterpolator->Evaluate(posDouble, &pdfValDouble);

			// Accumulate
			const auto pdfVal = static_cast<float>(pdfValDouble);
			value += PE_CHARGE * pdfVal;
		}
	}

	return value;
}



inline void updateGuessCorrector(const std::vector<double>& amps, const std::vector<double>& times,
                                 const std::vector<double>& initialAmps, const std::vector<double>& initialTimes,
                                 const double baseline, const double initBaseline, const std::vector<Photoelectron>& pesFound){
	for (int k   = 0; k < pesFound.size(); k++) {
		sysProcPECount++;
		ampDiff  = ampDiff + ((amps[k] - initialAmps[k]) - ampDiff) / sysProcPECount;
		timeDiff = timeDiff + ((times[k] - initialTimes[k]) - timeDiff) / sysProcPECount;
		// TODO(josh): Keep the below but only for debug builds
#ifdef IS_DEBUG
        if(ampDiff != ampDiff){
            throw std::runtime_error("AmpDiff is NaN");
        }
        if(timeDiff != timeDiff){
            throw std::runtime_error("TimeDiff is NaN");
        }
#endif
	}
	baselineDiff = baselineDiff + ((baseline - initBaseline) - baselineDiff) / sysProcPECount;
}


inline bool getNextPEGuess(DigitiserChannel *residualWF, Photoelectron *guessPE, const double baseline, std::vector<Photoelectron> pesFound_, const DigitiserChannel& channel_, const std::vector<double> *idealWF_, ceres::CubicInterpolator<ceres::Grid1D<float>>* PDFInterpolator, std
                           ::vector<float> xValues, const PETemplate* ChPETemplate) {
	// Get initial guesses for the next PE
	const auto         minPosIt   = std::ranges::min_element(residualWF->waveform);
	const unsigned int minTimePos = std::distance(residualWF->waveform.begin(), minPosIt);
	
	// If lowest point in waveform isn't below threshold there are no more PEs
	if (-residualWF->waveform[minTimePos] < WFSigThresh) {
		return false;
	}
	
	guessPE->amplitude = -residualWF->waveform[minTimePos];
	guessPE->time      = static_cast<float>(minTimePos) * trueSamplingRate;

    std::vector<Photoelectron> pesFoundLocal = std::move(pesFound_);
    pesFoundLocal.push_back(*guessPE);

    std::vector<float> paramsLocal;
    paramsLocal.reserve(2 * pesFoundLocal.size() + 2);
    paramsLocal.push_back(static_cast<float>(paramsLocal.size()));
    paramsLocal.push_back(static_cast<float>(baseline));
    for (const auto &pe: pesFoundLocal) {
        paramsLocal.push_back(pe.amplitude);
        paramsLocal.push_back(pe.time);
    }

	std::vector<float> tempResidual = channel_.waveform;
	for (unsigned int  k = 0; k < tempResidual.size(); ++k) {
		const float X = xValues[k];
		const float fitVal = NPEPDFFuncCubic(
			X,                    // time in index position
			paramsLocal,              // vector of parameters
			PDFInterpolator, // your cubic interpolator
			pdfNSamples,
			ChPETemplate
		);
		tempResidual[k] = tempResidual[k] - fitVal + static_cast<float>(baseline);
	}

//     // This is effectively checking in what direction the residual is skewed.
//     // (t1*a1)/(a1*a2*a3) + (t2*a2)/(a1*a2*a3) + (t3*a3)/(a1*a2*a3)
//     // If the estimated PE time is larger than truth the residual will be negative on the left,
//     // and positive on the right (of the PE time), this means a1/(a1*a2*a3) will be less than one,
//     // and a3/(a1+a2+a3) will be greater than one, shifting the time to the right...
//     if ((minTimePos > 1) && (minTimePos < tempResidual.size() - 1)) {
//         // improve initial time for a new PE based on average time
//         // over 3 consecutive sample ponderated by the amplitude
//         // of each sample... help a lot to resolve PEs very close!
//         float            timeSum        = 0;
//         float            ponderationSum = 0;
//         for (unsigned int b = minTimePos - 1; b <= minTimePos + 1; ++b) {
//             const float binCenter = (static_cast<float>(b) - 0.5f) * (trueSamplingRate);
//             const float binVal    = tempResidual[b]*1000;
//             timeSum += binCenter * binVal;
//             ponderationSum += binVal;
//         }
//         if((std::abs(timeSum) > 1e-10) & (std::abs(ponderationSum) > 1e-10)){
//             const double newVal = 0.015935 + timeSum/ponderationSum;
//             if(newVal != newVal){
//                 std::cout << "NewVal is NaN" << std::endl;
//             }
//             if(((newVal / trueSamplingRate) > 0) & ((newVal / trueSamplingRate) < static_cast<double>(tempResidual.size()))){
//                 guessPE->time = static_cast<float>(newVal);
//             }
// //            guessPE->time = newVal;
//         }
//     }
	return true;
}


inline void amplitudeCorrection(std::vector<Photoelectron> *pesFound, std::vector<float> *params, const std::vector<float>& waveform, const std::vector<double> *chIdealWF, ceres::CubicInterpolator<ceres::Grid1D<float>>* PDFInterpolator, std
						   ::vector<float> xValues, const PETemplate* ChPETemplate){
	for (int i = 0; i < pesFound->size(); i++) {
		// TODO(josh): Improve the adjustment by averaging the shift based off of a few bins around the PE time
		const unsigned int peTimeBinPos   = std::floor(pesFound->at(i).time / trueSamplingRate);
		const float fitVal = NPEPDFFuncCubic(
			ChPETemplate->getFractionalIndex(pesFound->at(i).time),                    // time in index position
			*params,              // vector of parameters
			PDFInterpolator, // your cubic interpolator
			pdfNSamples,
			ChPETemplate
		);
		const float        extraAmplitude = fitVal - waveform[peTimeBinPos];
		if (const float newAmplitude   = pesFound->at(i).amplitude + extraAmplitude; newAmplitude > WFSigThresh) { // TODO(josh): We need to consider the situations that this would ever be true
			params->at(2 + (2 * i)) = newAmplitude;
			pesFound->at(i).amplitude = newAmplitude;
		}
	}
}


void
fitEvent(const DigitiserEvent *event, const std::unordered_map<unsigned int, PETemplate>& PETemplates, std::shared_ptr<SyncFile> outputFile, std::mutex &lock) {
	std::vector<FitChannel> chFits;
	for (const auto &channel: event->channels) { // Looping through all channels for a given event
		auto residualWF = channel; // This will be the variable that is modified to be the residual distribution after each iteration

		// TODO(josh): Use time instead of index position for fitting procedure. Should just be that I need to create the interpolator using times and then remove the time->index conversion. Okay Ceres' interpolator cannot handle that. So will need
		//  to make some tool for making the interpolator. The requirement will be even spacing of the time values and passing the t0.
		FitChannel         chFit{};
		const unsigned int ch = channel.ID;
		if (std::ranges::count(skipChannels, ch)) {
			continue;
		}
		chFit.ID = ch;

		constexpr float preSigWindowTime = 18.f;  // ns. This is the time window before the signal that is used to calculate the noise level.
        std::vector<float> preSignalWF(residualWF.waveform.begin(), residualWF.waveform.begin() + static_cast<int>(preSigWindowTime / trueSamplingRate));
//        auto stdDevNoise = calculateVariance(preSignalWF, calculateMean(preSignalWF));
        auto stdDevNoise = calculateStandardDeviation(preSignalWF);
        pdfResidualRMS = static_cast<float>(stdDevNoise);
		
		// Making a pointer to the ideal waveform for this channel to improve speed of passing.
		if (ch > PETemplates.size() - 1) {
			throw std::runtime_error("Channel number exceeds the number of ideal waveforms");
		}

		const PETemplate* ChPETemplate = &PETemplates.at(ch);
		const std::vector<double>* chIdealWF = ChPETemplate->getVoltagesRef();
		
		// Baseline calculation
		float initBaseline = averageVector(residualWF.waveform, 0.01);
		chFit.baseline = initBaseline;

		// Start loop that will break when no more PEs are present
		std::vector<Photoelectron> pesFound;
		unsigned int               numPEsFound = 0;

		// This is getting estimate PEs that will then be passed as initial guesses to the minimiser.
		Photoelectron guessPE{};

		// Creating the x values that the solver will use. These are the index positions on the ideal WF for the positions on the real WF.
		std::vector<float> xValues;
		// std::vector<double> xValues2;
		for (unsigned int  j = 0; j < channel.waveform.size(); j++) {
			// xValues.push_back((static_cast<float>(j) * static_cast<float>(totalInterpFactor)) + pdfT0SampleConv);  // Multiplying index to match position on ideal PDF
			xValues.push_back(ChPETemplate->getFractionalIndex(static_cast<float>(j) * trueSamplingRate) + 1.f);
		}
		// TODO(josh): So what Ceres needs is corresponding index positions of waveform times, given that it uses an interpolator, and xValues are floats anyway, shouldn't need to have positions in the template that perfectly
		//  match those in the waveform. So we can just use the times of the waveform and convert them to the index positions on the ideal waveform.

		// Creating interpolator to allow for the ideal waveform to be used as a continuous (and differentiable) function.
		std::vector<float> temp(chIdealWF->begin(), chIdealWF->end());
		auto grid                 = ceres::Grid1D<float>(temp.data(), 0, static_cast<int>(chIdealWF->size()));
		auto               idealPDFInterpolator = new ceres::CubicInterpolator<ceres::Grid1D<float> >(grid);

		// Set up the only cost function (also known as residual). This uses
		// auto-differentiation to obtain the derivative (Jacobian).
		auto functor = new NPEPDFFunctor(xValues,
										 channel.waveform,
										 idealPDFInterpolator,
										 pesFound.size());
		auto costFunction = new ceres::DynamicAutoDiffCostFunction<NPEPDFFunctor>(functor);

		// Initial baseline removal.
		for (float &k: residualWF.waveform) {
			k = k - initBaseline;
		}

		while (true) {

			std::ranges::sort(pesFound, comparePETime()); // Need to sort for amplitude adjustment.

			std::vector<float> params;
			params.push_back(static_cast<float>(pesFound.size()));
			params.push_back(initBaseline);
			for (const auto &pe: pesFound) {
				params.push_back(pe.amplitude);
				params.push_back(pe.time);
			}

			// Amplitude adjustment: if the latest PE found is before other one(s),
			//  its tail is going to add some amplitude to the following one. Compares
			//  real and fit amplitude at the time bin corresponding to the PE time.
			amplitudeCorrection(&pesFound, &params, channel.waveform, chIdealWF, idealPDFInterpolator, xValues, ChPETemplate);

			// Compute residual
			for (unsigned int k = 0; k < xValues.size(); ++k) {
				float X = xValues[k];
				float fitVal = NPEPDFFuncCubic(
					X,                    // time in index position
					params,              // vector of parameters
					idealPDFInterpolator, // your cubic interpolator
					pdfNSamples,
					ChPETemplate
				);
				residualWF.waveform[k] = residualWF.waveform[k] - fitVal + initBaseline;
			}

			// Keep correcting baseline as new PEs are found.
			if(!pesFound.empty()){
				initBaseline = averageVector(residualWF.waveform, 0.01);
				for (float &k: residualWF.waveform) {
					k = k - initBaseline;
				}
			}

			if(!getNextPEGuess(&residualWF, &guessPE, initBaseline, pesFound, channel, chIdealWF, idealPDFInterpolator, xValues, ChPETemplate)){
				break;
			}

			numPEsFound += 1;
			pesFound.push_back(guessPE);
			residualWF = channel;

			if (numPEsFound >= maxPEs) {  // To handle the possibility of the algorithm being overly keen.
				break;
			}
		} // End of PE find loop


        // Carry out one final adjustment to the amplitude of the PEs found and recalculate baseline, residual.
        std::vector<float> params2;
        params2.push_back(static_cast<float>(pesFound.size()));
        params2.push_back(initBaseline);
        for (const auto &pe: pesFound) {
            params2.push_back(pe.amplitude);
            params2.push_back(pe.time);
        }
        amplitudeCorrection(&pesFound, &params2, channel.waveform, chIdealWF, idealPDFInterpolator, xValues, ChPETemplate);

        initBaseline = averageVector(residualWF.waveform, 0.01);
        // =========================================

        // TODO(JOSH): ADD A CHECKER TO SEE IF FIT IS GOOD ENOUGH AT THIS POINT AND THEN SKIP USING CERES IF IT IS.
		
		if(numPEsFound == 0){
			continue;
		}
		
		ceres::Problem problem{};
		
		std::vector<double> initialTimes;
		std::vector<double> initialAmplitudes;
		for (auto &k: pesFound) {
			initialAmplitudes.push_back(k.amplitude);
			initialTimes.push_back(k.time);
		}
		
		// Want to create a copy of the initial estimates to modify in below running correction.
		std::vector<double> times      = initialTimes;
		std::vector<double> amplitudes = initialAmplitudes;
		for (int i = 0; i < times.size(); i++) {
			times[i] += timeDiff;
			amplitudes[i] += ampDiff;
		}
		
		// Converting the time into an ideal waveform PDF index to simplify NPEPDFFunctor method call
		for (double &time: times) {
			time = time * samplingRate2Inv;
		}
		
		
		// Creating vector of references to parameter values that the fitter will use and modify. Note this means that
		// the references are the initial values before the fit and the final values after the fit.
		double              baseline = initBaseline;
		std::vector<double*> params   = {};
		params.push_back(&baseline);
		for (int i = 0; i < pesFound.size(); i++) {
			params.push_back(&amplitudes[i]);
			params.push_back(&times[i]);
		}

		functor->updateNumPEs(pesFound.size());

		costFunction->SetNumResiduals(static_cast<int>(channel.waveform.size()));

		costFunction->AddParameterBlock(1); // Baseline param
		for ([[maybe_unused]]const auto &pe: pesFound) {
			costFunction->AddParameterBlock(1); // Params for one PE amplitude
			costFunction->AddParameterBlock(1); // Params for one PE time
		}

		problem.AddResidualBlock(costFunction, nullptr, params);

//		for (auto & param : params) {
//			problem.SetParameterLowerBound(param, 0, 0);
//		}

        for (unsigned int i = 1; i < params.size(); i++) {
            problem.SetParameterLowerBound(params[i], 0, 0);
        }

		// Run the solver!
        if(std::abs(parameterTolerance) >= 1e-10){
            ceres::Solver::Options options;
            options.minimizer_progress_to_stdout = false;
            options.parameter_tolerance = parameterTolerance;

            std::unique_ptr<ceres::ParameterBlockOrdering> ordering(new ceres::ParameterBlockOrdering);

			// Add parameter blocks to the ordering object
            ordering->AddElementToGroup(params[0], 0);
            for (int i = 0; i < pesFound.size(); i++) {
                ordering->AddElementToGroup(params[2*i + 1], i);
                ordering->AddElementToGroup(params[2*i + 2], i);
            }
			// Continue for other parameters, possibly using more sophisticated grouping based on the analysis

			// Use this ordering in the solver options
            options.linear_solver_ordering = std::move(ordering);

            ceres::Solver::Summary summary;
            Solve(options, &problem, &summary);

//        std::cout << summary.FullReport() << "\n";
        }

		std::vector<double> postFitParams{};
		for (auto &param : params) {
			postFitParams.push_back(*param);
		}
		
		// Updating the amplitudes and times using params which is a vector of pointers to doubles
//		std::cout << "\n\n=====================" << std::endl;
//		std::cout << "Fit values" << std::endl;
//		std::cout << "Baseline : " << (postfitParams[0]) << std::endl;
//		baseline = (postfitParams[0]);
//		for (int i = 0; i < pesFound.size(); i++) {
//			std::cout << "Amplitude " << i << " : " << (postfitParams[2*i + 1]) << std::endl;
//			std::cout << "Time " << i << " : " << (postfitParams[2*i + 2]) << std::endl;
//			amplitudes[i] = (postfitParams[2*i + 1]);
//			times[i] = (postfitParams[2*i + 2]);
//		}
//		std::cout << "=====================" << std::endl;
		
		// Going back from ideal waveform PDF index to time
		for (double &time: times) {
			time = time / samplingRate2Inv;
		}
		
		std::vector<Photoelectron> FitPEs;
		for (int            k = 0; k < pesFound.size(); k++) {
			Photoelectron pe{};
			pe.amplitude = static_cast<float>(amplitudes[k]);
			pe.time      = static_cast<float>(times[k]);
			FitPEs.push_back(pe);
		}
		chFit.PEs      = FitPEs;
		chFit.baseline = static_cast<float>(baseline);
		
		std::vector<float> finalParams;
		finalParams.push_back(static_cast<float>(FitPEs.size()));
		finalParams.push_back(static_cast<float>(baseline));
		for (const auto &PE: FitPEs) {
			finalParams.push_back(static_cast<float>(PE.amplitude));
			finalParams.push_back(static_cast<float>(PE.time));
		}

        std::vector<double> observedValues;
        std::vector<double> predictedValues;

		float             chiSq    = 0;
		for (unsigned int j        = 0; j < channel.waveform.size(); j++) {
			const float observed = channel.waveform[j];
			const float expected = NPEPDFFunc(static_cast<float>(j) * trueSamplingRate, finalParams, chIdealWF);
            observedValues.push_back(observed);
            predictedValues.push_back(expected);
			chiSq += static_cast<float>(std::pow(observed - expected, 2)) / (pdfResidualRMS * pdfResidualRMS);
            //TODO(josh): Where does the 1000 come from? Is it to convert from mV to V?
			//TODO(josh): Need to calculate the residual RMS on a per waveform basis?
		}
		float            redChiSq = chiSq / (static_cast<float>(channel.waveform.size()) - (static_cast<float>(finalParams.size()) - 1));
		
		sysProcWFCount++;
		meanReducedChiSq = meanReducedChiSq + (redChiSq - meanReducedChiSq) / sysProcWFCount;
		
		chFit.reducedChiSq = redChiSq;
		
		lock.lock();
		reducedChiSqs.emplace_back(redChiSq);
		lock.unlock();
		
		chFits.push_back(chFit);
		
		// Update correction values for initial guesses
		updateGuessCorrector(amplitudes, times, initialAmplitudes, initialTimes, baseline, initBaseline, pesFound);
		
		delete idealPDFInterpolator;
		
//		delete[] x1;
//		for (int i = 0; i < pesFound.size(); i++) {
//			delete[] x2[i];
//		}
//		delete[] x2;
	}
	
	FitEvent evFitDat{event->ID, event->correctedTime, event->date, chFits};
	
	Writer writer(std::move(outputFile));
	writer.writeEventInfo(evFitDat);
}

/**
 *
 * @param events
 * @param count
 * @param lock
 * @param PETemplates
 * @param file
 * @return
 */
bool batchFitEvents(const std::vector<DigitiserEvent> &events, std::atomic<unsigned long> &count, std::mutex &lock,
                    const std::unordered_map<unsigned int, PETemplate>& PETemplates, const std::shared_ptr<SyncFile> &file) {
	for (const auto &event: events) {
		lock.lock();
		++count;
		lock.unlock();
		fitEvent(&event, PETemplates, file, lock);
	}
	return true;
}