Installation (one path that's working for me currently).

 - Activate your conda environment in console.
 - Navigate to this directory.
 - `python3 setup.py build_ext --inplace install`

Look at `DebugUtils/PlotChiSqDist.py`, `self.rawWFs = ReadWCDataFile(rawDataPath)` for an example of how to use CReader.