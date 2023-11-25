#include "particleContainers/ParticleContainer.h"

void ParticleContainer::calculatePosition()
{
	for (size_t i = 0; i < particles.size(); i++)
	{
		Particle &p = particles.at(i);
		std::array<double, 3> force = p.getF();
		double factor = std::pow(deltaT, 2) / (2 * p.getM());
		force = factor * force;
		std::array<double, 3> newPosition = p.getX() + deltaT * p.getV() + force;
		// TODO (ADD): Log
		// ParticleContainer::debugLog("The new position for particle {} is {}.\n", i, ArrayUtils::to_string(newPosition));

		p.setX(newPosition);
	}
}

void ParticleContainer::calculateVelocity()
{
	for (size_t i = 0; i < particles.size(); i++)
	{
		Particle &p = particles.at(i);
		double factor = deltaT / (2 * p.getM());
		std::array<double, 3> sumOfForces = p.getOldF() + p.getF();
		sumOfForces = factor * sumOfForces;
		std::array<double, 3> newVelocity = p.getV() + sumOfForces;

		// TODO (ADD): Log
		// ParticleContainer::debugLog("The new velocity for particle {} is {}.\n", i, ArrayUtils::to_string(newVelocity));
		p.setV(newVelocity);
	}
}

void iterOverPairs(const std::function<void(Particle &a, Particle &b)> &f)
{
}
void add(const std::array<double, 3> &x_arg, const std::array<double, 3> &v_arg, double mass, int type)
{
}