#include "particleContainers/ParticleContainerDirSum.h"
#include "ParticleContainerDirSum.h"
#include <iostream>

ParticleContainerDirSum::ParticleContainerDirSum(double deltaT, double endTime, int writeFrequency, std::function<void(Particle &a, Particle &b)> f) : ParticleContainer(deltaT, endTime, writeFrequency, f)
{
}

void ParticleContainerDirSum::add(const std::array<double, 3> &x_arg, const std::array<double, 3> &v_arg, double mass, int type)
{
	particles_.push_back(new Particle(x_arg, v_arg, mass, type));
}

void ParticleContainerDirSum::iterOverInnerPairs(const std::function<void(Particle &a, Particle &b)> &f)
{
	for (size_t i = 0; i < particles_.size() - 1; ++i)
	{
		Particle &pi = *particles_.at(i);
		for (size_t j = i + 1; j < particles_.size(); ++j)
		{
			Particle &pj = *particles_.at(j);
			f(pi, pj);
		}
	}
}

void ParticleContainerDirSum::calculateForces()
{
	for (auto &p : particles_)
	{
		auto oldForce = p->getF();
		std::array<double, 3> zero = {0.0, 0.0, 0.0};
		p->setF(zero);
		p->setOldF(oldForce);
	}
	iterOverInnerPairs(force_);
}
