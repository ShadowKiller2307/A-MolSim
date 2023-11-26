#include "particleContainers/ParticleContainerDirSum.h"

ParticleContainerDirSum::ParticleContainerDirSum(double deltaT, double endTime) : ParticleContainer(deltaT, endTime)
{
}

void ParticleContainerDirSum::add(const std::array<double, 3> &x_arg, const std::array<double, 3> &v_arg, double mass, int type)
{
	particles_.emplace_back(Particle(x_arg, v_arg, mass, type));
}

void ParticleContainerDirSum::iterOverInnerPairs(const std::function<void(Particle &a, Particle &b)> &f)
{
	for (auto &p : particles_)
	{
		auto oldForce = p.getF();
		std::array<double, 3> zero = {0.0, 0.0, 0.0};
		p.setF(zero);
		p.setOldF(oldForce);
	}
	for (size_t i = 0; i < particles_.size() - 1; ++i)
	{
		Particle &pi = particles_.at(i);
		for (size_t j = i + 1; j < particles_.size(); ++j)
		{
			Particle &pj = particles_.at(j);
			f(pi, pj);
		}
	}
}