
#include "FileReader.h"
#include "outputWriter/XYZWriter.h"
#include "utils/ArrayUtils.h"

#include <iostream>
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

constexpr double start_time = 0;
constexpr double end_time = 1000;
constexpr double delta_t = 0.014;

// TODO: what data structure to pick?
std::list<Particle> particles;

int main(int argc, char *argsv[]) {

  std::cout << "Hello from MolSim for PSE!" << std::endl;
  if (argc != 2) {
    std::cout << "Erroneous programme call! " << std::endl;
    std::cout << "./molsym filename" << std::endl;
  }

  FileReader fileReader;
  fileReader.readFile(particles, argsv[1]);

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

//TODO add documentation
double euclideanNorm(const std::array<double, 3>& arr) {
    double sum = 0.0;
    for (const auto& element : arr) {
        sum += element * element;
    }
    return std::sqrt(sum);
}

//TODO add documentation
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
      std::array<double, 3> resultingforce = it2->getF() + inverseForce;
      it2->setF(resultingforce);
    }
  }
}

void calculateX() {
  for (auto &p : particles) {
    // TODO: insert calculation of position updates here!
  }
}


 /// Calculates the new velocity according to Velocity-Stormer-Verlet.

void calculateV() {
  for (auto &p : particles) {

    // TODO: insert calculation of veclocity updates here!
    double two_times_mass = 2 * p.getM();
    std::array<double,3> sum_of_forces = p.getOldF() + p.getF();
    std::array<double,3> new_array = scalar_Operations(sum_of_forces,two_times_mass, false);

    std::array<double,3> new_velocity = p.getV()+(delta_t*new_array);
    p.setV(new_velocity);
  }
}

void plotParticles(int iteration) {

  std::string out_name("MD_vtk");

  outputWriter::XYZWriter writer;
  writer.plotParticles(particles, out_name, iteration);
}
