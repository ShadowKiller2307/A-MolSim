#pragma once
#include "forces/Force.h"

class LennJon : public Force
{
private:
	double ntfEpsilon_ /*(Negative TwentyFour Epsilon)*/, sigma_;

public:
	std::function<void(Particle &a, Particle &b)> innerPairs() override;
	std::function<void(Particle &a, Particle &b)> boundaryPairs() override;
	LennJon(double epsilon, double sigma);
	~LennJon() = default;
};
