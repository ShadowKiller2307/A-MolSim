#pragma once
#include <array>
#include <vector>
#include <stdint.h>
#include "Particle.h"

class ParticleContainer;

class particleGenerator
{
private:
	static int generateNumber_;

public:
	static void instantiateCuboid(ParticleContainer &container, const std::array<double, 3> &llfc, const std::array<unsigned int, 3> &particlePerDimension,
								  std::array<double, 3> &particleVelocity, double h, double mass, int generateNumber, double meanVelocity);
	static void instantiateSphere(ParticleContainer &container, const std::array<double, 3> &center, const uint32_t &sphereRadius, std::array<double, 3> &particleVelocity,
								  double h, double mass, int generateNumber, double meanVelocity);
	static void instantiateJSON(ParticleContainer &container, const std::string &path, int generateNumber, double meanVelocity);
	static void instantiateXML();
};
