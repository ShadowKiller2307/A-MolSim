#include "FileReader.h"
#include "MolSim.h"
#include "outputWriter/XYZWriter.h"
#include "outputWriter/VTKWriter.h"
#include "utils/ArrayUtils.h"
#include "ParticleContainer.h"
#include "ForceV1.h"


#include <iostream>
#include <fstream>
#include <getopt.h>
#include <vector>

constexpr double startTime{0};
double endTime{1000};
double deltaT{0.014};
bool outputModeVTK = true; // default to VTK

std::vector<Particle> particles;
// ParticleContainer particles
ParticleContainer particleContainer{};

//TODO: Implement mode for runtime measurement which disables all I/O
int main(int argc, char *argsv[]) {

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

    FileReader fileReader;
    //ParticleContainer particleContainer;
    fileReader.readFile(particles, argsv[optind]);
    //Initialising the ParticleContainer with particles
    particleContainer.setParticles(particles);
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

