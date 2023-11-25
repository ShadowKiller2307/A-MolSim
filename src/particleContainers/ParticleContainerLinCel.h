#include "particleContainers/ParticleContainer.h"
using cell = std::vector<Particle>;

class ParticleContainerLinCel : ParticleContainer
{
private:
	std::vector<cell> cells;

public:
	ParticleContainerLinCel(/* args */);
	~ParticleContainerLinCel();
};
