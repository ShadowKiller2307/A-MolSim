#include "GravPot.h"

std::function<void(Particle &a, Particle &b)> GravPot::innerPairs()
{
	return ([](Particle &a, Particle &b) {

	});
}

std::function<void(Particle &a, Particle &b)> GravPot::boundaryPairs()
{
	return ([](Particle &a, Particle &b) {

	});
}