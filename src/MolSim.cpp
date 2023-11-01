
#include "FileReader.h"
#include "MolSim.h"
#include "outputWriter/XYZWriter.h"
#include "outputWriter/VTKWriter.h"
#include "utils/ArrayUtils.h"

#include <iostream>
#include <fstream>
#include <getopt.h>
#include <vector>

constexpr double startTime{0};
double endTime{1000};
double deltaT{0.014};
bool outputModeVTK = true; // default to VTK
std::size_t numParticels{0};
std::vector<Particle> particles;

int main(int argc, char *argsv[]) {

  std::cout << "Hello from MolSim for PSE!" << std::endl;
  if (argc < 2) {
    std::cout << "Erroneous programme call! " << std::endl;
    printHelp();
    return EXIT_FAILURE;
  }

  option longOpts[] = {
    {"deltaT", required_argument, NULL, 'd'},
    {"endTime", required_argument, NULL, 'e'},
    {"help", no_argument, NULL, 'h'},
    {nullptr, 0, nullptr, 0}
  };

  int longOptsIndex = 0;
  while (int c = getopt_long(argc, argsv, "e:d:hx", longOpts, &longOptsIndex) != -1)
  {
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
    case 'x':
      outputModeVTK = false;
      break;
    default:
      break;
    }
  }

  if(optind >= argc){
    std::cout << "Inputfile missing as an argument, aborting" << std::endl;
    printHelp();
    return EXIT_FAILURE;
  }

  FileReader fileReader;
  fileReader.readFile(particles, argsv[optind]);


  double currentTime = startTime;
  int iteration = 0;

  std::time_t start = std::time(nullptr);
  // for this loop, we assume: current x, current f and current v are known
  while (currentTime < endTime) {
    // calculate new x
    calculateX();
    // calculate new f
    calculateF();
    // calculate new v
    calculateV();

    iteration++;
    if (iteration % 10 == 0) {
      plotParticles(iteration);
    }
    std::cout << "Iteration " << iteration << " finished." << std::endl;

    currentTime += deltaT;
  }
  std::time_t end = std::time(nullptr);
  auto diff = end - start;

  std::cout << "Output written, took " << diff <<" seconds. Terminating..." << std::endl;
  return 0;
}

/**
 * This function is used in the velocity calculation.
 * @param array The array to divide by the scalar. The values are changed inplace.
 * @param scalar The double to divide by.
 * @param isDivision The mode of the operation, where false equals multiplication and true equals division.
 */

void scalarOperations(std::array<double,3> &array, double scalar, bool isDivision){
  if (isDivision) {
    for (double & i : array) {
      i /= scalar;
    }
  }
  else {
    for (double & i : array) {
      i *= scalar;
    }
  }
}

/**
 * @brief Helper function to quickly calculate the euclidean norm of an array
 * 
 * @param arr The array in question to calculate the euclidean norm
 * @return The euclidean norm of the array as a double 
 */

double euclideanNorm(const std::array<double, 3>& arr) {
  double sum = 0.0;
  for (const auto& element : arr) {
    sum += element * element;
  }
  return std::sqrt(sum);
}

/// @brief Calculates the force acting on each particle by looping over them pairwise, calculating the force for each pair and adding it to them respectively

void calculateF() {
  for (auto & p : particles)
  {
    auto oldForce = p.getF();
    std::array<double, 3> zero = {0.0, 0.0, 0.0};
    p.setF(zero);
    p.setOldF(oldForce);
  }

  for (size_t i = 0; i < particles.size() - 1; ++i)
  {
    Particle &pi = particles.at(i);
    for (size_t j = i + 1; j < particles.size(); ++j)
    {
      Particle &pj = particles.at(j);
      double scalar = pi.getM() * pj.getM() / std::pow(euclideanNorm(pi.getX() - pj.getX()), 3);
      std::array<double, 3> force = scalar * (pj.getX() - pi.getX());
      std::array<double, 3> resultingForce = pi.getF() + force;
      pi.setF(resultingForce);

      scalarOperations(force, -1.0, false);
        resultingForce = pj.getF() + force;
      pj.setF(resultingForce);
    }
  }
}

///Calculates the new position for every particle according to Velocity-Störmer-Verlet

void calculateX() {
  for (auto &p : particles) {
    std::array<double, 3> force = p.getF();
    scalarOperations(force,  2*p.getM(), true);
    scalarOperations(force, std::pow(deltaT, 2), false);
    std::array<double, 3> newPosition = p.getX() + deltaT * p.getV() + force;
    p.setX(newPosition);
  }
}

/// Calculates the new velocity according to Velocity-Störmer-Verlet.

void calculateV() {
  for (auto &p : particles) {
    double twoTimesMass = 2 * p.getM();
    std::array<double,3> sumOfForces = p.getOldF() + p.getF();
    scalarOperations(sumOfForces, twoTimesMass, true);

    scalarOperations(sumOfForces,deltaT,false);
    std::array<double, 3> newVelocity = p.getV() + sumOfForces;
    p.setV(newVelocity);
  }
}

void plotParticles(int iteration) {
  std::string outName("output/MD_vtk");
  if(outputModeVTK){
    outputWriter::VTKWriter writer;
    writer.initializeOutput(particles.size());
    for (auto &p : particles)
    {
      writer.plotParticle(p);
    }
    writer.writeFile(outName, iteration);
  }
  else
  {
    outputWriter::XYZWriter writer;
    writer.plotParticles(particles, outName, iteration);
  }
}

void printHelp() {
  std::ifstream file("help.txt");

  if (file.is_open()){
    std::cout << file.rdbuf();
  }
}
