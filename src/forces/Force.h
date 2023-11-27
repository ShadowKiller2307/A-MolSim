#pragma once
#include "Particle.h"
#include "utils/ArrayUtils.h"
#include <functional>

class Force
{
public:
	virtual std::function<void(Particle &a, Particle &b)> innerPairs();
	virtual std::function<void(Particle &a, Particle &b)> boundaryPairs();
};
