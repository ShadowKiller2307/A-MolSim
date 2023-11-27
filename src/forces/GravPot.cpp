#include "GravPot.h"

std::function<void(Particle &a, Particle &b)> GravPot::innerPairs()
{
	return ([](Particle &a, Particle &b)
			{
			double scalar = a.getM() * b.getM() / std::pow(ArrayUtils::L2Norm(a.getX()-b.getX()), 3);
        	std::array<double, 3> force = scalar * (b.getX() - a.getX());
        	std::array<double, 3> resultingForce = a.getF() + force;
        	a.setF(resultingForce);

			force = -1 * force;
			resultingForce = b.getF() + force;
        	b.setF(resultingForce); });
}

std::function<void(Particle &a, Particle &b)> GravPot::boundaryPairs()
{
	return ([](Particle &a, Particle &b) {

	});
}