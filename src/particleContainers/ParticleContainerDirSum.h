#pragma once
#include "particleContainers/ParticleContainer.h"

class ParticleContainerDirSum : public ParticleContainer
{
private:
public:
	ParticleContainerDirSum(double deltaT, double endTime, int writeFrequency, std::function<void(Particle &a, Particle &b)> f);
	~ParticleContainerDirSum() = default;
	void add(const std::array<double, 3> &x_arg, const std::array<double, 3> &v_arg, double mass, int type) override;
	void iterOverInnerPairs(const std::function<void(Particle &a, Particle &b)> &f) override;
};
