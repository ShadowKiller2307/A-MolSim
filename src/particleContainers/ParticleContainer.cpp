#include "particleContainers/ParticleContainer.h"
#include "logOutputManager/LogManager.h"

ParticleContainer::ParticleContainer(const double deltaT, const double endTime) : deltaT_(deltaT), endTime_(endTime)
{
	outManager_ = outputManager();
	startTime_ = 0.0;
}

void ParticleContainer::iterOverInnerPairs(const std::function<void(Particle &a, Particle &b)> &f)
{
}

void ParticleContainer::calculateForces()
{
}

void ParticleContainer::calculateVelocity()
{
	for (size_t i = 0; i < particles_.size(); i++)
	{
		Particle &p = particles_.at(i);
		double factor = deltaT_ / (2 * p.getM());
		std::array<double, 3> sumOfForces = p.getOldF() + p.getF();
		sumOfForces = factor * sumOfForces;
		std::array<double, 3> newVelocity = p.getV() + sumOfForces;

		// TODO (ADD): Log
		// ParticleContainer::debugLog("The new velocity for particle {} is {}.\n", i, ArrayUtils::to_string(newVelocity));
		p.setV(newVelocity);
	}
}

void ParticleContainer::calculatePosition()
{
	for (size_t i = 0; i < particles_.size(); i++)
	{
		Particle &p = particles_.at(i);
		std::array<double, 3> force = p.getF();
		double factor = std::pow(deltaT_, 2) / (2 * p.getM());
		force = factor * force;
		std::array<double, 3> newPosition = p.getX() + deltaT_ * p.getV() + force;
		// TODO (ADD): Log
		// ParticleContainer::debugLog("The new position for particle {} is {}.\n", i, ArrayUtils::to_string(newPosition));

		p.setX(newPosition);
	}
}

void ParticleContainer::add(const std::array<double, 3> &x_arg, const std::array<double, 3> &v_arg, double mass, int type)
{
}

void ParticleContainer::simulateParticles()
{
	while (startTime_ < endTime_)
	{
		// if (LogManager::getInstance().getOutputFiles() && iteration % outputEveryNIterations == 0)
		// {
		// 	plotParticles(iteration);
		// }
		outManager_.plotParticles(particles_, iteration_);
		if (iteration_ % 100 == 0)
		{
			// TODO (ADD): Log
			// fileLogger->log(logManager.getLevel(), "Iteration {} finished. ({}%)", iteration, std::round(iteration * 10000 / (endTime / deltaT)) / 100);
		}

		// calculate new x
		calculatePosition();
		// calculate new f
		calculateForces();
		// calculate new v
		calculateVelocity();

		iteration_++;
		startTime_ += deltaT_;
	}
}

std::vector<Particle> &ParticleContainer::getParticles()
{
	return particles_;
}

void ParticleContainer::setParticles(std::vector<Particle> &particles)
{
	this->particles_ = particles;
}