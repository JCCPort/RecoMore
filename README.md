This is an alternative to the current LiquidO RecoZoR implementation that doesn't require ROOT as a dependency, exists as a standalone tool separate to other LiquidO analysis software, and also can be run on multiple cores. Combining the improved single thread performance with multi-core processing, RecoMore can very easily be 1-2 orders of magnitude faster than the original software.

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
