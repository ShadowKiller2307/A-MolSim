#include "particleContainers/ParticleContainer.h"

class ParticleContainerDirSum : public ParticleContainer
{
private:
public:
	ParticleContainerDirSum();
	~ParticleContainerDirSum();
	void add(const std::array<double, 3> &x_arg, const std::array<double, 3> &v_arg, double mass, int type) override;
	void iterOverPairs(const std::function<void(Particle &a, Particle &b)> &f) override;
};
