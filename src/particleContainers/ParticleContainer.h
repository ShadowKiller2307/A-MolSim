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
	int outputEveryNIterations_;
	size_t iteration_;
	std::vector<Particle> particles_;
	OutputManager *outManager_;
	std::function<void(Particle &a, Particle &b)> force_;
	ParticleContainer(double deltaT, double endTime, int writeFrequency, std::function<void(Particle &a, Particle &b)> f);

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

	/**
	 * @brief runs the simulation
	 * @param None
	 * @return void
	 */
	void simulateParticles();

	void writeJSON(std::string &name);

	std::vector<Particle> &getParticles();

	void setParticles(std::vector<Particle> &particles);

	const double getDeltaT() const;

	const double getEndTime() const;

	~ParticleContainer();
};