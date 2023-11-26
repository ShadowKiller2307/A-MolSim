#pragma once
#include "utils/ArrayUtils.h"
#include "logOutputManager/outputManager.h"
#include "Particle.h"
#include <functional>
#include <vector>

class ParticleContainer
{
protected:
	double deltaT_, startTime_, endTime_;
	size_t iteration_;
	std::vector<Particle> particles_;
	ParticleContainer(const double deltaT, const double endTime);
	outputManager outManager_;

public:
	/**
	 * @brief iterate over every particle pair in the container and apply the lambda function
	 * @param f the lambda function
	 */
	virtual void iterOverInnerPairs(const std::function<void(Particle &a, Particle &b)> &f);

	/**
	 * @brief calculate the forces for all particles
	 * @param None
	 * @return void
	 */
	virtual void calculateForces();

	/**
	 * @brief calculate the velocity for all particles
	 * @param None
	 * @return void
	 */
	virtual void calculateVelocity();

	/**
	 * @brief calculate the position for all particles
	 * @param None
	 * @return void
	 */
	virtual void calculatePosition();

	/**
	 * @brief adds a particle to the container
	 *	@param x_arg position of the particle
	 *	@param v_arg velocity of the particle
	 *	@param mass mass of the particle
	 *	@param type typenumber of the particle
	 */
	virtual void add(const std::array<double, 3> &x_arg, const std::array<double, 3> &v_arg, double mass, int type);

	void simulateParticles();

	std::vector<Particle> &getParticles();

	void setParticles(std::vector<Particle> &particles);
};