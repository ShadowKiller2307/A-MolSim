#pragma once
#include "forces/Force.h"

class GravPot : public Force
{
public:
	std::function<void(Particle &a, Particle &b)> innerPairs() override;
	std::function<void(Particle &a, Particle &b)> boundaryPairs() override;
};
