#include "FileReader.h"
#include "MolSim.h"
#include "outputWriter/XYZWriter.h"
#include "outputWriter/VTKWriter.h"
#include "utils/ArrayUtils.h"
#include "ParticleContainer.h"
#include "ForceV1.h"
#include "spdlog/spdlog.h"
#include "ParticleGenerator.h"


#include <iostream>
#include <fstream>
#include <getopt.h>
#include <vector>
#include <nlohmann/json.hpp>
using json = nlohmann::json;

constexpr double startTime{0};
double endTime{1000};
double deltaT{0.014};
bool generate = false; //bool which indicates whether the particles are read from a file or generated
ParticleGenerator particleGenerator;
bool outputModeVTK = true; // default to VTK
//default values for the particle generator
std::array<double, 3> llfc{0.0, 0.0, 0.0}; // lower left frontside corner
std::array<unsigned int, 3> particlePerDimension{10, 10, 10};
double h = 0.5;
double mass = 1.0; //<- vllt noch was sinnvolleres hierhin
std::array<double, 3> particleVelocity = {1.0, 1.0, 1.0};




std::vector<Particle> particles;
// ParticleContainer particles
ParticleContainer particleContainer{};

//TODO: Implement mode for runtime measurement which disables all I/O
int main(int argc, char *argsv[]) {
    spdlog::info("Erste Nachricht durch den Logger");
    std::cout << "Hello from MolSim for PSE!" << std::endl;
    if (argc < 2) {
        std::cout << "Erroneous programme call! " << std::endl;
        printHelp();
        return EXIT_FAILURE;
    }


    option longOpts[] = {
            {"deltaT",  required_argument, nullptr, 'd'},
            {"endTime", required_argument, nullptr, 'e'},
            {"help",    no_argument,       nullptr, 'h'},
            {nullptr, 0,                   nullptr, 0}
    };

    int longOptsIndex = 0;
    //TODO: Extend with the command line arguments for the generator
    while (int c = getopt_long(argc, argsv, "e:d:hx", longOpts, &longOptsIndex) != -1) {
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
            case 'x':
                outputModeVTK = false;
                break;
            default:
                break;
        }
    }

    if (optind >= argc) {
        std::cout << "Input file missing as an argument, aborting" << std::endl;
        printHelp();
        return EXIT_FAILURE;
    }

    if (!generate) {
        FileReader fileReader;
        //ParticleContainer particleContainer;
        fileReader.readFile(particles, argsv[optind]);
        //Initialising the ParticleContainer with particles
        particleContainer.setParticles(particles);
    } else{
        // Here the particles will be generated by the ParticleGenerator
        particleGenerator.instantiateCuboid(particleContainer, llfc, particlePerDimension, h, mass, particleVelocity, 0);
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

    std::time_t start = std::time(nullptr);
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
            plotParticles(iteration);
        }
        //TODO: das durch logging Nachricht ersetzen
        std::cout << "Iteration " << iteration << " finished." << std::endl;

        currentTime += deltaT;
    }
    std::time_t end = std::time(nullptr);
    auto diff = end - start;

    std::cout << "Output written, took " << diff << " seconds. Terminating..." << std::endl;
    return 0;
}

void plotParticles(int iteration) {
    std::string outName("../output/MD_vtk");

    if (outputModeVTK) {
        outputWriter::VTKWriter writer;
        writer.initializeOutput(particleContainer.getParticles()->size());
        for (auto &p: *particleContainer.getParticles()) {
            writer.plotParticle(p);
        }
        writer.writeFile(outName, iteration);
    } else {
        outputWriter::XYZWriter writer;
        writer.plotParticles(*particleContainer.getParticles(), outName, iteration);
    }
}

void printHelp() {
    std::ifstream file("../help.txt");

    if (file.is_open()) {
        std::cout << file.rdbuf();
    }
}

