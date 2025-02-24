This is an alternative implementation of the RecoZoR algorithm currently used in LiquidO that doesn't require ROOT as a dependency, exists as a standalone tool separate to other LiquidO analysis software, and also can be run on multiple cores. Combining the improved single thread performance with multi-core processing, RecoMore can very easily be 1-2 orders of magnitude faster than the original software.

The dependencies are all easily installable on OSX and Linux systems, however Ceres-solver is not as easy to install for Windows. The current workaround I'm using is to install WSL2 and then run RecoMore on there - note that if try to run RecoMore on data that is outside of your WSL virtual disk it will be very slow.

**Dependencies**:

- Boost
- Ceres-solver
- CMake
- pybind11 (for CReader)
- pybind11-stubgen (for CReader type hints)

For Debian systems, these packages are easily available using apt install:
- libboost-all-dev
- libceres-dev
- cmake

## RecoMore

To build RecoMore:
- Make a build directory
- `cd` into build directory
- `cmake -DCMAKE_BUILD_TYPE=Release <pathToProjectDir>`
- `make`

To run RecoMore:

```./PEFinder -i <pathToDataFile> --pdf_dir <pathToIdealPDFsDir> --n_threads <numberOFThreads> --n_batches <numberOfBatches>```

For more info run `./PEFinder -h` to get the help output.

Note: Make sure to include the trailing slash of the path to the PDF dir

## CReader

CReader is a python wrapper around C++ file parsers to combine high performance with easy usage. To initialise
pybind11 simply run `git submodule update --init` from the repo root directory.

To build and install CReader to be a usable python module:
- Navigate to the `CReader_` directory.
- Make sure you are in the same python environment that you intend to run CReader in.
- `python3 setup.py build_ext --inplace install`
- To rebuild the stubs (used for IDE type hints etc):
  ```
  pybind11-stubgen CReader
     --output-dir="./stubs/generated"
     --root-module-suffix=""
     --ignore-invalid=all
     --no-setup-py
  ```

Import with `from CReader import *`


## DebugUtils

DebugUtils contains some utilities for debugging/analysing RecoMore output along with some examples of using
CReader. 
