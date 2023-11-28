#pragma once
#include "logOutputManager/OutputManager.h"
#include "utils/ArrayUtils.h"
#include "Particle.h"

#include <functional>
#include <vector>

class OutputManager;

class ParticleContainer
{
protected:
	double deltaT_, startTime_, endTime_;
	size_t iteration_;
	std::vector<Particle> particles_;
	ParticleContainer(const double deltaT, const double endTime);
	OutputManager *outManager_;
	std::function<void(Particle &a, Particle &b)> force;
	int outputEveryNIterations_;

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

	void writeJSON(std::string &name);

	std::vector<Particle> &getParticles();

	void setParticles(std::vector<Particle> &particles);

	void setForce(const std::function<void(Particle &a, Particle &b)> f);

	const double getDeltaT() const;

	const double getEndTime() const;

	~ParticleContainer();
};