#pragma once
#include <array>
#include <vector>
#include <stdint.h>
#include "Particle.h"

class ParticleContainer;

/*
	to use, pass for example to a function
	optinals({.deltaT = 0.00069, .endTime = 7})
	to set h and mass but use te default values
	for meanVelocity and generateNumber
*/
struct SimParams
{
public:
	double deltaT = 0.0002;
	double endTime = 5;
};

class particleGenerator
{
private:
	static int generateNumber_;
	static constexpr const double meanVelocity_ = 0.1;
	static constexpr const int h_ = 1.122462048;
	static constexpr const int m_ = 1.0;

	static ParticleContainer *createContainer(double deltaT, double endTime);
	static ParticleContainer *createContainer(double deltaT, double endTime, std::array<double, 3>);

public:
	static void instantiateCuboid(ParticleContainer **container, const std::array<double, 3> &llfc,
								  const std::array<unsigned int, 3> &particlesPerDimension, std::array<double, 3> &particleVelocity, double h, double m, int type);

	static void instantiateSphere(ParticleContainer **container, const std::array<double, 3> &center,
								  const uint32_t &sphereRadius, std::array<double, 3> &particleVelocity, double h, double m, bool is2D, int type);

	static void instantiateJSON(ParticleContainer **container, const std::string &path, SimParams params);

	static void instantiatePicture(ParticleContainer **container, const std::string &path, SimParams params);

	static void instantiateTxt(ParticleContainer **container, const std::string &path, SimParams params);

	static void instantiateXML(ParticleContainer **container, const std::string &path, SimParams params);
};
