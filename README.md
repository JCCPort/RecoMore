This is an alternative implementation of the RecoZoR algorithm currently used in LiquidO that doesn't require ROOT as a dependency, exists as a standalone tool separate to other LiquidO analysis software, and also can be run on multiple cores. Combining the improved single thread performance with multi-core processing, RecoMore can very easily be 1-2 orders of magnitude faster than the original software.

The dependencies are all easily installable on OSX and Linux systems, however Ceres-solver is not as easy to install for Windows. The current workaround I'm using is install WSL2 and then run RecoMore on there - note that if try to run RecoMore on data that is outside of your WSL virtual disk it will be very slow.

Dependencies:

- Boost
- Ceres-solver
- CMake


To build:
- Make a build directory
- `cd` into build directory
- `cmake <pathToProjectDir>`
- `make`


To run:
- `./RecoMore <pathToDataFile> <pathToIdealPDFs>`

Note: Make sure to include the trailing slash of the path to the PDF dir
