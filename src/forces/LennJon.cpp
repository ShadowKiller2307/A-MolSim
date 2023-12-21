#include "LennJon.h"
#include <iostream>

std::function<void(Particle &a, Particle &b)> LennJon::innerPairs()
{
	return ([this](Particle &a, Particle &b)
			{
			auto diff = a.getX() - b.getX();
			double norm = ArrayUtils::L2Norm(diff);
			double first = this->ntfEpsilon_ / (norm * norm);
			double frac = this->sigma_ / norm;
			double pow6 = std::pow(frac, 6);
			double pow12 = 2 * std::pow(pow6, 2);
			double middle = (pow6 - pow12);
			auto force = (first * middle) * diff;
			a.addF(force);
			force = -1 * force;
			b.addF(force); });
}

std::function<void(Particle &a, Particle &b)> LennJon::boundaryPairs()
{
	return ([this](Particle &a, Particle &b)
			{
			auto diff = a.getX() - b.getX();
			double norm = ArrayUtils::L2Norm(diff);
			double first = this->ntfEpsilon_ / (norm * norm);
			double frac = this->sigma_ / norm;
			double pow6 = std::pow(frac, 6);
			double pow12 = 2 * std::pow(pow6, 2);
			double middle = (pow6 - pow12);
			auto force = (first * middle) * diff;
			a.addF(force); });
}

LennJon::LennJon(double epsilon, double sigma)
{
	ntfEpsilon_ = -24 * epsilon;
	sigma_ = sigma;
}
