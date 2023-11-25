#include "particleGenerator/ParticleGenerator.h"
#include "particleContainers/ParticleContainer.h"
#include "utils/MaxwellBoltzmannDistribution.h"
#include "utils/ArrayUtils.h"
#include "Particle.h"
#include <vector>

int particleGenerator::generateNumber_ = 0;

void particleGenerator::instantiateCuboid(ParticleContainer &container, const std::array<double, 3> &llfc, const std::array<unsigned int, 3> &particlePerDimension,
										  std::array<double, 3> &particleVelocity, double h, double mass, int generateNumber = generateNumber_, double meanVelocity = 0.1)
{
	if (generateNumber < 0)
	{
		generateNumber = generateNumber_;
	}
	for (size_t i = 0; i < particlePerDimension[0]; ++i)
	{
		for (size_t j = 0; j < particlePerDimension[1]; ++j)
		{
			for (size_t k = 0; k < particlePerDimension[2]; ++k)
			{
				// TODO (ASK): should the mbVelocity be calculated for every particle instead?
				std::array<double, 3> mbVelocity = MaxwellBoltzmannDistribution::maxwellBoltzmannDistributedVelocity(meanVelocity, 2); // TODO (ASK): Is the mean here the same as the average?
				std::array<double, 3> x_arg{i * h + llfc[0], j * h + llfc[1], k * h + llfc[2]};
				std::array<double, 3> v_arg{particleVelocity + mbVelocity}; // Calculate via the Brownian motion
				container.add(x_arg, v_arg, mass, generateNumber);			// using the add function in order to be able to add elements to both container implementations
			}
		}
	}
}
