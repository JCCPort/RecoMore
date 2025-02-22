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


struct NPEPDFFunctor {
	/**
	 * 
	 * @param waveformTimes These need to be the indices of the PETemplate corresponding to the samples of the waveform being processed,
	 * @param waveformAmplitudes 
	 * @param PDFInterpolator
	 * @param numPDFSamples
	 * @param numPES 
	 */
	NPEPDFFunctor(std::vector<float> waveformTimes, std::vector<float> waveformAmplitudes,
	              ceres::CubicInterpolator<ceres::Grid1D<float, true>> *PDFInterpolator, const unsigned int numPDFSamples, const unsigned int numPES)
			: waveformTimes_(std::move(waveformTimes)), waveformAmplitudes_(std::move(waveformAmplitudes)), PDFInterpolator_(PDFInterpolator), numPDFSamples_(numPDFSamples), numPES_(numPES) {};

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
			T accum = baseline;
			const T X = T(waveformTimes_[j]);
			for (unsigned int i = 0; i < numPES_; i++) {
				const unsigned int i2 = 2 * i;
				const T amplitude = params[i2 + 1][0];
				const T timeShift = params[i2 + 2][0];

				const T pos = X - timeShift;
				if (pos >= T(0) && pos < T(numPDFSamples_)) {
					T f;
					PDFInterpolator_->Evaluate(pos, &f);
					accum += amplitude * f;
				}
			}
			residual[j] = accum - T(waveformAmplitudes_[j]);
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
	unsigned int 									numPDFSamples_;
	unsigned int									numPES_;
};


inline float NPEPDFFuncCubic(
	const float X,
	const FitParams& fitParams,
	const ceres::CubicInterpolator<ceres::Grid1D<float, true>>* PDFInterpolator,
	const PETemplate* ChPETemplate)
{
	float value = fitParams.baseline_;

	for (int PE = 0; PE < fitParams.numPEs_; ++PE) {
		const float PE_CHARGE = fitParams.PEParams_[PE].amplitude;
		const float PE_TIME   = fitParams.PEParams_[PE].time  / ChPETemplate->getTimeSpacing();

		// pos is the position on the PE template waveform relative to the PE time in terms of index position
		const float pos = X - PE_TIME;

		// Now we want to interpolate if we're in range
		if (pos >= 0.0f && pos < static_cast<float>(ChPETemplate->getNumberOfSamples())) {
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
                                 const double baseline, const double initBaseline, const unsigned int numPEs){
	for (int k   = 0; k < numPEs; k++) {
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


inline bool getNextPEGuess(DigitiserChannel *residualWF, Photoelectron *guessPE) {
	// Get initial guesses for the next PE
	const auto         minPosIt   = std::ranges::min_element(residualWF->waveform);
	const unsigned int minTimePos = std::distance(residualWF->waveform.begin(), minPosIt);

	// If lowest point in waveform isn't below threshold there are no more PEs
	if (-residualWF->waveform[minTimePos] < WFSigThresh) {
		return false;
	}

	guessPE->amplitude = -residualWF->waveform[minTimePos];
	guessPE->time      = static_cast<float>(minTimePos) * trueSamplingRate;

	return true;
}


inline void amplitudeCorrection(FitParams fitParams, const std::vector<float>& waveform, const ceres::CubicInterpolator<ceres::Grid1D<float>>* PDFInterpolator, const PETemplate* ChPETemplate){
	for (int i = 0; i < fitParams.numPEs_; ++i) {
		// TODO(josh): Improve the adjustment by averaging the shift based off of a few bins around the PE time
		const unsigned int peTimeBinPos   = std::floor(fitParams.PEParams_[i].time / trueSamplingRate);

		const float fitVal = NPEPDFFuncCubic(
			ChPETemplate->getFractionalIndex(fitParams.PEParams_[i].time),                    // time in index position
			fitParams,              // vector of parameters
			PDFInterpolator, // your cubic interpolator
			ChPETemplate
		);
		const float        extraAmplitude = fitVal - waveform[peTimeBinPos];
		if (const float newAmplitude   = fitParams.PEParams_[i].amplitude + extraAmplitude; newAmplitude > WFSigThresh) { // TODO(josh): We need to consider the situations that this would ever be true
			fitParams.PEParams_[i].amplitude = newAmplitude;
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
		unsigned int               numPEsFound = 0;

		// This is getting estimate PEs that will then be passed as initial guesses to the minimiser.
		Photoelectron guessPE{};

		// Creating the x values that the solver will use. These are the index positions on the ideal WF for the positions on the real WF.
		std::vector<float> xValues;
		for (unsigned int  j = 0; j < channel.waveform.size(); j++) {
			xValues.push_back(ChPETemplate->getFractionalIndex(static_cast<float>(j) * trueSamplingRate) + 1.f);
		}

		// Creating interpolator to allow for the ideal waveform to be used as a continuous (and differentiable) function.
		std::vector<float> temp(chIdealWF->begin(), chIdealWF->end());
		auto grid                 = ceres::Grid1D<float>(temp.data(), 0, static_cast<int>(chIdealWF->size()));
		auto               idealPDFInterpolator = new ceres::CubicInterpolator<ceres::Grid1D<float> >(grid);

		// Set up the only cost function (also known as residual). This uses
		// auto-differentiation to obtain the derivative (Jacobian).
		auto functor = new NPEPDFFunctor(xValues,
										 channel.waveform,
										 idealPDFInterpolator,
										 ChPETemplate->getNumberOfSamples(),
										 0);
		auto costFunction = new ceres::DynamicAutoDiffCostFunction<NPEPDFFunctor>(functor, ceres::Ownership::DO_NOT_TAKE_OWNERSHIP);

		// Initial baseline removal.
		for (float &k: residualWF.waveform) {
			k = k - initBaseline;
		}

		FitParams fitParams;

		while (true) {
			fitParams.sortPEsByTime();
			fitParams.baseline_ = initBaseline;

			// Amplitude adjustment: if the latest PE found is before other one(s),
			//  its tail is going to add some amplitude to the following one. Compares
			//  real and fit amplitude at the time bin corresponding to the PE time.
			amplitudeCorrection(fitParams, channel.waveform, idealPDFInterpolator, ChPETemplate);

			// Compute residual
			for (unsigned int k = 0; k < xValues.size(); ++k) {
				float fitVal = NPEPDFFuncCubic(
					xValues[k],                    // time in index position
					fitParams,              // vector of parameters
					idealPDFInterpolator, // your cubic interpolator
					ChPETemplate
				);
				residualWF.waveform[k] = residualWF.waveform[k] - fitVal + initBaseline;
			}

			// Keep correcting baseline as new PEs are found.
			if(fitParams.numPEs_ != 0){
				initBaseline = averageVector(residualWF.waveform, 0.01);
				for (float &k: residualWF.waveform) {
					k = k - initBaseline;
				}
			}

			if(!getNextPEGuess(&residualWF, &guessPE)){
				break;
			}

			numPEsFound += 1;
			fitParams.addPE(guessPE);

			residualWF = channel;

			if (numPEsFound >= maxPEs) {  // To handle the possibility of the algorithm being overly keen.
				break;
			}
		} // End of PE find loop


        // Carry out one final adjustment to the amplitude of the PEs found and recalculate baseline, residual.
		amplitudeCorrection(fitParams, channel.waveform, idealPDFInterpolator, ChPETemplate);

        initBaseline = averageVector(residualWF.waveform, 0.01);
        // =========================================

        // TODO(JOSH): ADD A CHECKER TO SEE IF FIT IS GOOD ENOUGH AT THIS POINT AND THEN SKIP USING CERES IF IT IS.
		
		if(numPEsFound == 0){
			continue;
		}

		// Want to create a copy of the initial estimates for later improvement of correction factors.
		std::vector<double> initialTimes;
		std::vector<double> initialAmplitudes;
		for (auto &k: fitParams.PEParams_) {
			initialAmplitudes.push_back(k.amplitude);
			initialTimes.push_back(k.time);
		}
		

		std::vector<double> times;
		std::vector<double> amplitudes;
		double              baseline;

		// Creating vector of references to parameter values that the fitter will use and modify. Note this means that
		// the references are the initial values before the fit and the final values after the fit.
		std::vector<double*> params   = {};
		fitParams.makeSolverParams(&params, &times, &amplitudes, &baseline);

		// Make adjustment with running correction factors.
		for (int i = 0; i < times.size(); i++) {
			times[i] += timeDiff;
			amplitudes[i] += ampDiff;
		}
		
		// Converting the time into an ideal waveform PDF index to simplify NPEPDFFunctor method call
		for (double &time: times) {
			time = time / ChPETemplate->getTimeSpacing();
		}

		costFunction->SetNumResiduals(static_cast<int>(channel.waveform.size()));

		costFunction->AddParameterBlock(1); // Baseline param
		for ([[maybe_unused]]const auto &pe: fitParams.PEParams_) {
			costFunction->AddParameterBlock(1); // Params for one PE amplitude
			costFunction->AddParameterBlock(1); // Params for one PE time
		}

		functor->updateNumPEs(fitParams.numPEs_);

		ceres::Problem::Options option;
		option.loss_function_ownership = ceres::Ownership::DO_NOT_TAKE_OWNERSHIP;
		option.cost_function_ownership = ceres::Ownership::DO_NOT_TAKE_OWNERSHIP;
		option.manifold_ownership = ceres::Ownership::DO_NOT_TAKE_OWNERSHIP;
		ceres::Problem problem(option);
		problem.AddResidualBlock(costFunction, nullptr, params);

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
            for (int i = 0; i < fitParams.numPEs_; i++) {
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
			time = time * ChPETemplate->getTimeSpacing();
		}

		FitParams finalFitParams{fitParams.numPEs_, baseline, amplitudes, times};

		chFit.PEs      = finalFitParams.PEParams_;
		chFit.baseline = static_cast<float>(baseline);

        std::vector<double> observedValues;
        std::vector<double> predictedValues;


		float chiSq = 0;
		for (unsigned int j = 0; j < channel.waveform.size(); j++) {
			const float observed = channel.waveform[j];
			float expected = NPEPDFFuncCubic(
				xValues[j],                    // time in index position
				finalFitParams,              // vector of parameters
				idealPDFInterpolator, // your cubic interpolator
				ChPETemplate
			);
            observedValues.push_back(observed);
            predictedValues.push_back(expected);
			chiSq += static_cast<float>(std::pow(observed - expected, 2)) / (pdfResidualRMS * pdfResidualRMS);
		}
		float            redChiSq = chiSq / (static_cast<float>(channel.waveform.size()) - (static_cast<float>(finalFitParams.getNumParams()) - 1));
		
		sysProcWFCount++;
		meanReducedChiSq = meanReducedChiSq + (redChiSq - meanReducedChiSq) / sysProcWFCount;
		
		chFit.reducedChiSq = redChiSq;
		
		lock.lock();
		reducedChiSqs.emplace_back(redChiSq);
		lock.unlock();
		
		chFits.push_back(chFit);
		
		// Update correction values for initial guesses
		updateGuessCorrector(amplitudes, times, initialAmplitudes, initialTimes, baseline, initBaseline, finalFitParams.numPEs_);
		
		delete idealPDFInterpolator;
		delete functor;
		delete costFunction;
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