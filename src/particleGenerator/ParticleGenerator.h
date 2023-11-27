#pragma once
#include <array>
#include <vector>
#include <stdint.h>
#include "Particle.h"

class ParticleContainer;

/*
	to use, pass for example to a function
	optinals({.h = 0.9, .mass = 4.0})
	to set h and mass but use te default values
	for meanVelocity and generateNumber
*/
struct optionals
{
public:
	double h = 1.122462048;
	double meanVelocity = 0.1;
	int generateNumber = -1;
	double mass = 1.0;
	double deltaT = -1.0;
	double endTime = -1.0;
};

class particleGenerator
{
private:
	static int generateNumber_;

public:
	static void instantiateCuboid(ParticleContainer **container, const std::array<double, 3> &llfc, const std::array<unsigned int, 3> &particlesPerDimension,
								  std::array<double, 3> &particleVelocity, optionals optionalArguments);

	static void instantiateSphere(ParticleContainer **container, const std::array<double, 3> &center, const uint32_t &sphereRadius, std::array<double, 3> &particleVelocity, optionals optionalArguments);

	static void instantiateJSON(ParticleContainer **container, const std::string &path, optionals optionalArguments);

	static void instantiatePicture(ParticleContainer **container, const std::string &path, optionals optionalArguments);

	static void instantiateTxt(ParticleContainer **container, const std::string &path, optionals optionalArguments);

	static void instantiateXML();

	static ParticleContainer *createContainer(double deltaT, double endTime);
	static ParticleContainer *createContainer(double deltaT, double endTime, std::array<double, 3>);
};
