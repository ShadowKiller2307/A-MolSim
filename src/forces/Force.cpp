#include "forces/Force.h"
#include "Force.h"

std::function<void(Particle &a, Particle &b)> Force::innerPairs()
{
	return std::function<void(Particle & a, Particle & b)>();
}
std::function<void(Particle &a, Particle &b)> Force::boundaryPairs()
{
	return std::function<void(Particle & a, Particle & b)>();
}
