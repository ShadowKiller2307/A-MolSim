#pragma once
#include "utils/ArrayUtils.h"
#include "Particle.h"
#include <functional>

class Force
{
public:
	virtual std::function<void(Particle &a, Particle &b)> innerPairs();
	virtual std::function<void(Particle &a, Particle &b)> boundaryPairs();
};
