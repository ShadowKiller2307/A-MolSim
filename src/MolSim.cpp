
#include "FileReader.h"
#include "outputWriter/XYZWriter.h"
#include "utils/ArrayUtils.h"

#include <iostream>
#include <fstream>
#include <getopt.h>
#include <list>

/**** forward declaration of the calculation functions ****/

/**
 * calculate the force for all particles
 */
void calculateF();

/**
 * calculate the position for all particles
 */
void calculateX();

/**
 * calculate the position for all particles
 */
void calculateV();

/**
 * plot the particles to a xyz-file
 */
void plotParticles(int iteration);

void printHelp();

constexpr double start_time{0};
double end_time{1000};
double delta_t{0.014};

// TODO: what data structure to pick?
std::list<Particle> particles;

int main(int argc, char *argsv[]) {

  std::cout << "Hello from MolSim for PSE!" << std::endl;
  if (argc < 2) {
    std::cout << "Erroneous programme call! " << std::endl;
    printHelp();
    return EXIT_FAILURE;
  }

  option long_opts[] = {
    {"deltaT", required_argument, NULL, 'd'},
    {"endTime", required_argument, NULL, 'e'},
    {"help", no_argument, NULL, 'h'},
    {NULL, 0, NULL, 0}
  };

  int long_opts_index = 0;
  while (int c = getopt_long(argc, argsv, "e:d:h", long_opts, &long_opts_index) != -1)
  {
    switch (c)
    {
    case 'd':
      delta_t = std::stod(optarg);
      break;
    case 'e':
      end_time = std::stod(optarg);
      break;
    case 'h':
      printHelp();
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

  double current_time = start_time;

  int iteration = 0;

  // for this loop, we assume: current x, current f and current v are known
  while (current_time < end_time) {
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

    current_time += delta_t;
  }

  std::cout << "output written. Terminating..." << std::endl;
  return 0;
}

/**
 * This function is used in the velocity calculation.
 * @param array The array I want to divide by the scalar.
 * @param scalar The double I want to divide by.
 * @param mode The mode of the operation, where false equals division and true equals multiplication.
 * @return A new std::array<double,3> where every value is the value of the parameter array divided by the parameter scalar.
 */

std::array<double,3> scalar_Operations(std::array<double,3> &array, double scalar, bool mode){
  if (mode) {
    std::array<double, 3> new_array{}; //TODO: Do we need a new array or should we just modify it
    for (size_t i = 0; i < array.size(); ++i) {
      new_array[i] = array[i] / scalar;
    }
    return new_array;
  }
  else {
    std::array<double, 3> new_array{};
    for (size_t i = 0; i < array.size(); ++i) {
      new_array[i] = array[i] * scalar;
    }
    return new_array;
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

/**
 * @brief Calculates the force acting on each particle by looping over them pairwise, calculating the force for each pair and adding it to them respectively
 * 
 */

void calculateF() {
  std::list<Particle>::iterator it1;
  std::list<Particle>::iterator it2;
  it1 = particles.begin();

  for (; it1 != (particles.end()--); it1++)
  {
    it2 = it1;
    it2++;
    for (; it2 != particles.end(); it2++)
    {
      double scalar = it1->getM() * it2->getM() / std::pow(euclideanNorm(it1->getX() - it1->getX()), 3);
      std::array<double, 3> force = scalar * (it2->getX() - it1->getX());
      std::array<double, 3> resultingforce = it1->getF() + force;
      it1->setF(resultingforce);

      auto inverseForce = scalar_Operations(force, -1.0, false);
      resultingforce = it2->getF() + inverseForce;
      it2->setF(resultingforce);
    }
  }
}

///Calculates the new position for every particle according to Velocity-Störmer-Verlet

void calculateX() {
  for (auto &p : particles) {
    // TODO: insert calculation of position updates here!
    std::array<double, 3> force = p.getF();
    std::array<double, 3> temp_array = scalar_Operations(force,  2*p.getM(), false);
    std::array<double,3> new_position = p.getX() + delta_t * p.getV() + scalar_Operations(temp_array, std::pow(delta_t, 2), true);
    p.setF(new_position);
  }
}


/// Calculates the new velocity according to Velocity-Störmer-Verlet.

void calculateV() {
  for (auto &p : particles) {
    // TODO: insert calculation of veclocity updates here!
    double two_times_mass = 2 * p.getM();
    std::array<double,3> sum_of_forces = p.getOldF() + p.getF();
    std::array<double,3> new_array = scalar_Operations(sum_of_forces, two_times_mass, true);

    std::array<double,3> new_velocity = p.getV()+(scalar_Operations(new_array,delta_t,false));
    p.setV(new_velocity);
  }
}

void plotParticles(int iteration) {

  std::string out_name("MD_vtk");

  outputWriter::XYZWriter writer;
  writer.plotParticles(particles, out_name, iteration);
}

void printHelp()
{
  std::ifstream file("help.txt");

  if (file.is_open()){
    std::cout << file.rdbuf();
  }
}
