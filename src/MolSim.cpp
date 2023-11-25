#include "FileReader.h"
#include "MolSim.h"
#include "outputWriter/XYZWriter.h"
#include "outputWriter/VTKWriter.h"
#include "utils/ArrayUtils.h"
#include "ParticleContainer.h"
#include "ParticleContainerDS.h"

#include "spdlog/sinks/rotating_file_sink.h"
#include "ParticleGenerator.h"
#include "LogLevel.h"
#include "LogManager.h"
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

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
    txtFile,
    xml
};

/// @brief the current way the program is using to get the particles, for worksheet 2 default to using the particle generator
particleSources pSource = generator;

auto fileLogger = LogManager::getInstance().getLogger();
LogManager &logManager = LogManager::getInstance();

/// @brief the current logLevel
LogLevel logLevel = standard;

/// @brief which output filetype to use, default to VTK
bool outputModeVTK = true;

bool writeToJSON = false;

/// @brief prefix for the name of each outputfile
std::string outName{"MD_vtk"};

// default values for the particle generator
ParticleGenerator particleGenerator;
double h = 1.122462048;
double mass = 1.0; //<- vllt noch was sinnvolleres hierhin

ParticleContainer* particleContainer;

int main(int argc, char *argsv[])
{
    ParticleContainerDS containerDs{};
    particleContainer = &containerDs;
    logManager.setLogLevel(toSpdLevel(standard));

    fileLogger->log(fileLogger->level(), "First message of the logger.\n");
    fileLogger->log(fileLogger->level(), "Hello from MolSim for PSE\n");

    if (argc < 2)
    {
        fileLogger->log(logManager.getLevel(), "Erroneous program call\n");

        printHelp();
        return EXIT_FAILURE;
    }

    option longOpts[] = {
        {"deltaT", required_argument, nullptr, 'd'},
        {"endTime", required_argument, nullptr, 'e'},
        {"help", no_argument, nullptr, 'h'},
        {"logLevel", required_argument, nullptr, 'l'},
        // TODO: update help.txt with these
        {"inputGenerator", no_argument, nullptr, 'g'}, // funny
        {"inputPicture", no_argument, nullptr, 'p'},   // because
        {"inputText", no_argument, nullptr, 't'},      // GPT
        {"outName", required_argument, nullptr, 'o'},
        {nullptr, 0, nullptr, 0}};

    int longOptsIndex = 0;
    while (true)
    {
        int c = getopt_long(argc, argsv, "d:e:hl:xgptw", longOpts, &longOptsIndex);
        if (c == -1)
        {
            break;
        }
        switch (c)
        {
        case 'd':
            deltaT = std::stod(optarg);
            break;
        case 'e':
            endTime = std::stod(optarg);
            break;
        case 'h':
            printHelp();
            break;
        case 'l':
        {
            int temp = std::stoi(optarg);
            if (temp >= 0 && temp <= 4)
            {
                logLevel = static_cast<LogLevel>(temp);
                spdlog::level::level_enum level = toSpdLevel(logLevel);
                logManager.setLogLevel(level);
            }
            break;
        }
        case 'g':
            pSource = generator;
            break;
        case 'p':
            pSource = picture;
            break;
        case 't':
            pSource = txtFile;
            break;
        case 'x':
            pSource = xml;
            break;
        case 'o':
            outName = optarg;
            break;
        case 'w':
            writeToJSON = true;
            break;
        default:
            fileLogger->log(logManager.getLevel(),
                            "Error parsing arguments. Maybe you gave one that isn't recognized.\n");

            break;
        }
    }

    if (optind >= argc)
    {
        fileLogger->log(logManager.getLevel(), "Input missing as an argument, aborting\n");

        printHelp();
        return EXIT_FAILURE;
    }

    switch (pSource)
    {
    case generator:
    {
        std::ifstream f(std::string("_").compare(argsv[optind]) == 0 ? "../input/collision.json" : argsv[optind]);
        json jsonParticles = json::parse(f);

        for (size_t i = 0; i < jsonParticles[0]; i++)
        {
            std::array<double, 3UL> llfc = jsonParticles[i + 1]["x"];
            std::array<double, 3UL> particleVelocity = jsonParticles[i + 1]["v"];
            std::array<unsigned int, 3UL> particlePerDimension = jsonParticles[i + 1]["N"];
            // Here the particles will be generated by the ParticleGenerator
            particleGenerator.instantiateCuboid(*particleContainer, llfc, particlePerDimension, particleVelocity, h,
                                                mass, i);
        }
        break;
    }
    case picture:
    {
        // gross C code, brace yourself

        int width, height, bpp;
        uint8_t *rgb_image = stbi_load(std::string("_").compare(argsv[optind]) == 0 ? "../input/Cool MolSim.png" : argsv[optind], &width, &height, &bpp, 3);
        std::vector<Particle> particles;
        for (int i = 0; i < height; ++i)
        {
            for (int j = 0; j < width; j++)
            {
                int index = i * width * 3 + j * 3;
                uint8_t r = rgb_image[index + 0];
                uint8_t g = rgb_image[index + 1];
                uint8_t b = rgb_image[index + 2];
                if (r != 255 || g != 255 || b != 255)
                {
                    std::array<double, 3> x_arg{j * h, -i * h, 0};
                    std::array<double, 3> v_arg{0.0};
                    if (r == 255)
                    {
                        v_arg = {0.0, -20.0, 0.0};
                    }
                    particles.emplace_back(x_arg, v_arg, 1, r == 255 ? 1 : 0);
                }
            }
        }
        particleContainer->setParticles(particles);
        stbi_image_free(rgb_image);
        break;
    }
    case txtFile:
    {
        FileReader fileReader;
        std::vector<Particle> particles;
        fileReader.readFile(particles, (std::string("_").compare(argsv[optind]) == 0 ? (char *)"../input/eingabe-sonne.txt" : argsv[optind]));
        // Initialising the ParticleContainerDS with particles
        particleContainer->setParticles(particles);
        break;
    }
    case xml:
        // TODO
        break;
    }

    if (writeToJSON)
    {
        // TODO
    }

    particleContainer->setForceCalculator(1);
    /*
     * Mode 0 GravitationalForce
     * Mode 1 LennardJonesForce
     * Mode 2 infinity: restliche Forces
     */
    particleContainer->setDeltaTTwo(deltaT);

    double currentTime = startTime;
    int iteration = 0;
    if (logLevel < noCOut)
    {
        fileLogger->log(logManager.getLevel(),
                        outputModeVTK ? "Plotting particles with VTK." : "Plotting particles with XYZ\n");
    }
    fileLogger->log(logManager.getLevel(), "This might take a while\n");

    // Start measuring time
    auto begin = std::chrono::high_resolution_clock::now();

    // for this loop, we assume: current x, current f and current v are known
    while (currentTime < endTime)
    {
        if (iteration % 10 == 0)
        {
            if (logLevel < noFiles)
            {
                plotParticles(iteration);
            }
            if (iteration % 100 == 0)
            {
                fileLogger->log(logManager.getLevel(), "Iteration {} finished. ({}%)", iteration, std::round(iteration * 10000 / (endTime / deltaT)) / 100);
            }
        }

        // TODO maybe some sort of iteratorOverTheHaloParticles to delete the particles in the
        // Halo area
        // calculate new x
        particleContainer->calculatePosition();
        // calculate new f
        particleContainer->calculateForces();
        // calculate new v
        particleContainer->calculateVelocity();

        iteration++;
        currentTime += deltaT;
    }

    // Stop measuring time and calculate the elapsed time
    auto end = std::chrono::high_resolution_clock::now();
    int64_t diff = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();

    fileLogger->log(logManager.getLevel(), "Output written, took {} milliseconds. (about {}) Terminating...\n", diff,
                    (iteration > diff ? std::to_string(iteration / diff) + " iter/ms" : std::to_string(diff / iteration) + " ms/iter"));

    return 0;
}

void plotParticles(int iteration)
{
    if (outputModeVTK)
    {
        outputWriter::VTKWriter writer;
        writer.initializeOutput(particleContainer->getParticles().size());
        for (auto &p : particleContainer->getParticles())
        {
            writer.plotParticle(p);
        }
        writer.writeFile("../output/" + outName, iteration);
    }
    else
    {
        outputWriter::XYZWriter writer;
        writer.plotParticles(particleContainer->getParticles(), "../output/" + outName, iteration);
    }
}

void printHelp()
{
    std::ifstream file("../help.txt");
    if (file.is_open())
    {
        std::cout << file.rdbuf();
    }
}

spdlog::level::level_enum toSpdLevel(LogLevel level)
{
    switch (level)
    {
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
    default:
        return spdlog::level::info;
    }
}
