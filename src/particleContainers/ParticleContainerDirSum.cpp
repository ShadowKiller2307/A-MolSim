#include "particleContainers/ParticleContainerDirSum.h"
#include "ParticleContainerDirSum.h"
#include <iostream>

ParticleContainerDirSum::ParticleContainerDirSum(double deltaT, double endTime, int writeFrequency) : ParticleContainer(deltaT, endTime, writeFrequency)
{
}

void ParticleContainerDirSum::add(const std::array<double, 3> &x_arg, const std::array<double, 3> &v_arg, double mass, int type, double epsilon, double sigma)
{
	particles_.emplace_back(x_arg, v_arg, mass, type);
}

void ParticleContainerDirSum::addCompleteParticle(Particle &p)
{
	particles_.push_back(p);
}

void ParticleContainerDirSum::iterOverInnerPairs(const std::function<void(Particle &a, Particle &b)> &f)
{
	for (size_t i = 0; i < particles_.size() - 1; ++i)
	{
		Particle &pi = particles_.at(i);
		for (size_t j = i + 1; j < particles_.size(); ++j)
		{
			Particle &pj = particles_.at(j);
			calcF(pi, pj);
			std::cout << "force pi" << pi.getF() << std::endl;
			std::cout << "force pj" << pj.getF() << std::endl;
		}
	}
}

void ParticleContainerDirSum::iterOverAllParticles(const std::function<void(std::vector<Particle>::iterator)> &f)
{
	for (auto it = particles_.begin();;)
	{
		if (it == particles_.end())
		{
			break;
		}
		f(it);
	}
}

size_t ParticleContainerDirSum::getAmountOfParticles() const
{
	return particles_.size();
}

void ParticleContainerDirSum::calculateForces()
{
	for (auto &p : particles_)
	{
		auto oldForce = p.getF();
		std::array<double, 3> zero = {0.0, 0.0, 0.0};
		p.setF(zero);
		p.setOldF(oldForce);
	}
	iterOverInnerPairs(force_);
}
