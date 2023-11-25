#include "particleContainers/ParticleContainerDirSum.h"

void ParticleContainerDirSum::add(const std::array<double, 3> &x_arg, const std::array<double, 3> &v_arg, double mass, int type)
{
	particles.emplace_back(Particle(x_arg, v_arg, mass, type));
}

void ParticleContainerDirSum::iterOverPairs(const std::function<void(Particle &a, Particle &b)> &f)
{
	for (auto &p : particles)
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
			f(pi, pj);
		}
	}
}