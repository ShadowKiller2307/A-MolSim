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
public:
	virtual ~ParticleContainer() = 0;

	/**
	 * @brief iterate over every particle pair in the container and apply the lambda function
	 * @param f the lambda function
	 */
	virtual void iterOverInnerPairs(const std::function<void(Particle &a, Particle &b)> &f) = 0;

	virtual void iterOverAllParticles(const std::function<void(std::vector<Particle>::iterator)> &f) = 0;

	/**
	 * @brief calculate the forces for all particles
	 * @param None
	 * @return void
	 */
	virtual void calculateForces() = 0;

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
	virtual Particle *add(const std::array<double, 3> &x_arg, const std::array<double, 3> &v_arg, double mass, int type) = 0;

	virtual void addCompleteParticle(Particle &p) = 0;

	/**
	 * @brief runs the simulation
	 * @param None
	 * @return void
	 */
	virtual void simulateParticles();

	void writeJSON(std::string &name);

	std::vector<Particle> &getParticles();

	void setParticles(std::vector<Particle> &particles);

	const double getDeltaT() const;

	const double getEndTime() const;

	void calcF(Particle &a, Particle &b);

	virtual size_t getAmountOfParticles() const = 0;
};