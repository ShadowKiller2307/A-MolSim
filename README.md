# MolSim Group A

## Members

- Alexander Wachenfeld
- Erick Andreas Lazar
- Felix Huber

## Code

- Available on [GitHub](https://github.com/ShadowKiller2307/A-MolSim)
- Branch: Refactor
- Commit:

## Tools and requirements

- Compiler: GCC 13.1.0
- CMake: 3.22.1
- Make: 4.3

Other compilers, such as Clang, and to a reasonable extent, older versions of these, might - and most probably will - work as well. This wasn't verified by us however. Other requirements include:

- Libxerces via `sudo apt install libxerces-c-dev`
- Doxygen via `sudo apt install doxygen`

All the following dependencies the program relies on will be downloaded automatically by cmake:

- spdlog
- gtest
- [json](https://github.com/nlohmann/json)

## Build

run `mkdir build && cd build && cmake .. && make` in the root folder of the project to

1. create a folder and
2. switch to it
3. create the build files
4. make the executeable

## Run

The main executeable **MolSim** is located in `build/src` and can be run with `./src/MolSim <CL arguments>`  
Possible arguments and other recommended tips can be obtained by passing in `--help` or `-h`  
The output will be generated in a folder `output` created by cmake in the root directory
