#include "FileReader.h"
#include "MolSim.h"
#include "outputWriter/XYZWriter.h"
#include "outputWriter/VTKWriter.h"
#include "utils/ArrayUtils.h"
#include "ParticleContainer.h"
#include "ForceV1.h"
#include "spdlog/spdlog.h"
#include "spdlog/sinks/rotating_file_sink.h"
#include "ParticleGenerator.h"
#include "LogLevel.h"

#include <iostream>
#include <fstream>
#include <getopt.h>
#include <vector>
#include <chrono>
#include <nlohmann/json.hpp>

using json = nlohmann::json;

constexpr double startTime{0};
/// @brief how long the simulation should run
double endTime{5};

/// @brief timestep between each iteration
double deltaT{0.0002};

/// @brief possible ways to obtain a set of particles as well as their starting positions and velocities
enum particleSources
{
    generator,
    picture,
    txtFile
};

/// @brief the current way the program is using to get the particles, for worksheet 2 default to using the particle generator
particleSources pSource = generator;

LogLevel logLevel = standard;

auto maxSize = 5 * 1024 * 1024;
auto maxFiles = 4;

auto fileLogger = spdlog::rotating_logger_mt("fileLogger", "/logs/log.txt", maxSize, maxFiles);

/// @brief how much output the program should generate to stdout and (log-)files
enum LogLevel
{
    debug,           // debug, print all, slowest
    standard,        // default, no debug
    noCOut,          // disable std::cout and logging but still write Files, minor performance improvement
    noFiles,         // don't write output to files but still std::cout and logging, big performance improvement but no data :(
    onlyCalculations // combination of noCOut and noFiles, no output whatsoever, fastest
};
/// @brief the current logLevel
LogLevel logLevel = standard;

/// @brief which output filetype to use, default to VTK
bool outputModeVTK = true;

/// @brief prefix for the name of each outputfile
std::string outName{"MD_vtk"};

// default values for the particle generator
ParticleGenerator particleGenerator;
double h = 0.5;
double mass = 1.0; //<- vllt noch was sinnvolleres hierhin

ParticleContainer particleContainer{};

int main(int argc, char *argsv[]) {

    fileLogger->set_level(toSpdLevel(logLevel));
    fileLogger->set_pattern("[%Y-%m-%d %H:%M] [%l] %v");


    fileLogger->log(fileLogger->level(), "First message of the logger.\n");
    fileLogger->log(fileLogger->level(), "Hello from MolSim for PSE\n");

    if (argc < 2) {
        fileLogger->log(fileLogger->level(),"Erroneous program call\n");

        printHelp();
        return EXIT_FAILURE;
    }

    option longOpts[] = {
        {"deltaT", required_argument, nullptr, 'd'},
        {"endTime", required_argument, nullptr, 'e'},
        {"help", no_argument, nullptr, 'h'}, // print the help.txt file to stdout
        {"logLevel", required_argument, nullptr, 'l'},
        // TODO: update help.txt with these
        {"inputGenerator", no_argument, nullptr, 'g'}, // funny
        {"inputPicture", no_argument, nullptr, 'p'},   // because
        {"inputText", no_argument, nullptr, 't'},      // GPT
        {"outName", required_argument, nullptr, 'o'},
        {nullptr, 0, nullptr, 0}};

    int longOptsIndex = 0;
    // TODO: Extend with the command line arguments for the generator
    while (true) {
        int c = getopt_long(argc, argsv, "d:e:hl:xgpt", longOpts, &longOptsIndex);
        if (c == -1) {
            break;
        }
        switch (c) {
            case 'd':
                deltaT = std::stod(optarg);
                break;
            case 'e':
                endTime = std::stod(optarg);
                break;
            case 'h':
                printHelp();
                break;
            case 'l': {
                int temp = std::stoi(optarg);
                if (temp >= 0 && temp <= 4) {
                    logLevel = static_cast<LogLevel>(temp);
                    spdlog::level::level_enum level = toSpdLevel(logLevel);
                    fileLogger->set_level(level);


                }
                break;
            }
            case 'x':
                outputModeVTK = false;
                outName = "MD_xyz";
                break;
            case 'g':
                pSource = generator;
                break;
            case 'p':
                pSource = picture;
                break;
            case 't':
                pSource = txtFile;
                break;
            case 'o':
                outName = optarg;
                break;
            default:
                fileLogger->log(fileLogger->level(),
                                "Error parsing arguments. Maybe you gave one that isn't recognized.\n");

                break;
        }
    }

    if (optind >= argc) {
        fileLogger->log(fileLogger->level(), "Input missing as an argument, aborting\n");

        printHelp();
        return EXIT_FAILURE;
    }

    switch (pSource) {
        case generator: {
            std::ifstream f(std::string("_").compare(argsv[optind]) == 0 ? "../input/collision.json" : argsv[optind]);
            json jsonParticles = json::parse(f);

            for (size_t i = 0; i < jsonParticles[0]; i++) {
                std::array<double, 3UL> llfc = jsonParticles[i + 1]["x"];
                std::array<double, 3UL> particleVelocity = jsonParticles[i + 1]["v"];
                std::array<unsigned int, 3UL> particlePerDimension = jsonParticles[i + 1]["N"];
                // Here the particles will be generated by the ParticleGenerator
                particleGenerator.instantiateCuboid(particleContainer, llfc, particlePerDimension, particleVelocity, h,
                                                    mass, i);
            }
            break;
        }
    }
    case picture:
        break;
    case txtFile:
        FileReader fileReader;
        std::vector<Particle> particles;
        fileReader.readFile(particles, (std::string("_").compare(argsv[optind]) == 0 ? (char *)"../input/eingabe-sonne.txt" : argsv[optind]));
        // Initialising the ParticleContainer with particles
        particleContainer.setParticles(particles);
        break;
    }

    particleContainer.setForceCalculator(0);
    /*
     * Mode 0 GravitationalForce
     * Mode 1 LennardJonesForce
     * Mode 2 infinity: restliche Forces
     */
    particleContainer.setDeltaTTwo(deltaT);

    double currentTime = startTime;
    int iteration = 0;
    if (logLevel < noCOut) {
        fileLogger->log(fileLogger->level(),
                        outputModeVTK ? "Plotting particles with VTK." : "Plotting particles with XYZ\n");

    }
    fileLogger->log(fileLogger->level(), "This might take a while\n");


    // Start measuring time
    auto begin = std::chrono::high_resolution_clock::now();

    // for this loop, we assume: current x, current f and current v are known
    while (currentTime < endTime) {
        // calculate new x
        particleContainer.calculatePosition();
        // calculate new f
        particleContainer.calculateForces();
        // calculate new v
        particleContainer.calculateVelocity();

        iteration++;
        if (iteration % 10 == 0) {
            if (logLevel < noFiles) {
                plotParticles(iteration);
            }
            if ((logLevel < noCOut || logLevel == noFiles) && iteration % 100 == 0) {
                fileLogger->log(fileLogger->level(), "Iteration {} finished. ({}%)", iteration, std::round(iteration * 10000 / (endTime / deltaT)) / 100);
            }
        }
        currentTime += deltaT;
    }

    // Stop measuring time and calculate the elapsed time
    auto end = std::chrono::high_resolution_clock::now();
    int64_t diff = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();

    fileLogger->log(fileLogger->level(), "Output written, took {} milliseconds. (about {}) Terminating...\n", diff,
                    (iteration > diff ? std::to_string(iteration / diff) + " iter/ms" :
                     std::to_string(diff / iteration) +
                     " ms/iter"));

    return 0;
}

void plotParticles(int iteration) {
    if (outputModeVTK) {
        outputWriter::VTKWriter writer;
        writer.initializeOutput(particleContainer.getParticles().size());
        for (auto &p: particleContainer.getParticles()) {
            writer.plotParticle(p);
        }
        writer.writeFile("../output/" + outName, iteration);
    } else {
        outputWriter::XYZWriter writer;
        writer.plotParticles(particleContainer.getParticles(), "../output/" + outName, iteration);
    }
}

void printHelp() {
    std::ifstream file("../help.txt");
    if (file.is_open()) {
        std::cout << file.rdbuf();
    }
}

spdlog::level::level_enum toSpdLevel(LogLevel level) {
    switch (level) {
//TODO
        case debug:
            return spdlog::level::debug;
        case standard:
            return spdlog::level::info;
        case noCOut:
            return spdlog::level::off;
        case noFiles:
            return spdlog::level::info;
        case onlyCalculations:
            return spdlog::level::off;
    }
}
