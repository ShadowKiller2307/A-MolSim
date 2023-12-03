#include "particleContainers/ParticleContainer.h"
#include "boundaryConditions/Reflecting.h"
#include "logOutputManager/LogManager.h"
#include "boundaryConditions/Outflow.h"
#include <iostream>
#include <iomanip>

ParticleContainer::ParticleContainer(double deltaT, double endTime, int writeFrequency, std::function<void(Particle &a, Particle &b)> f) : deltaT_(deltaT), endTime_(endTime)
{
	outManager_ = new OutputManager();
	startTime_ = 0.0;
	if (writeFrequency <= 0)
	{
		outManager_->outputFiles = false;
	}
	outputEveryNIterations_ = writeFrequency;
	force_ = f;
}

void ParticleContainer::iterOverInnerPairs(const std::function<void(Particle &a, Particle &b)> &f)
{
}

void ParticleContainer::calculateForces()
{
	// TODO REMOVE
	LogManager::errorLog("dont want to be here");
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
	auto begin = std::chrono::high_resolution_clock::now();

	while (startTime_ < endTime_)
	{
		if (outManager_->outputFiles && iteration_ % outputEveryNIterations_ == 0)
		{
			outManager_->plotParticles(particles_, iteration_);
		}
		if (iteration_ % 100 == 0)
		{
			LogManager::infoLog("Iteration {} finished. ({}%)", iteration_, std::round(iteration_ * 10000 / (endTime_ / deltaT_)) / 100);
			if (iteration_)
			{
				auto end = std::chrono::high_resolution_clock::now();
				size_t diff = std::chrono::duration_cast<std::chrono::seconds>(end - begin).count();
				auto remainig = (static_cast<double>(diff) / iteration_) * (endTime_ - startTime_) / deltaT_;
				std::string min, sec;
				if (remainig > 60)
				{
					min = std::to_string(static_cast<int>(remainig / 60)) + "m, ";
					remainig = std::fmod(remainig, 60);
				}
				else
				{
					min = "0m, ";
				}
				sec = std::to_string(static_cast<int>(remainig)) + "s   \r";
				std::cout << "ETA: " << min << sec << std::flush;
			}
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
	//  TODO (ADD): Log
	auto end = std::chrono::high_resolution_clock::now();
	size_t diff = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();
	std::cout << "Output written, took " + std::to_string(diff) + " milliseconds. (about " + (iteration_ > diff ? std::to_string(iteration_ / diff) + " iter/ms" : std::to_string(diff / iteration_) + " ms/iter") + ") Terminating...\n";
}

void ParticleContainer::writeJSON(std::string &name)
{
	outManager_->writeJSON(name, *this);
}

const double ParticleContainer::getDeltaT() const
{
	return deltaT_;
}

const double ParticleContainer::getEndTime() const
{
	return endTime_;
}

ParticleContainer::~ParticleContainer()
{
	delete outManager_;
}

std::vector<Particle> &ParticleContainer::getParticles()
{
	return particles_;
}

void ParticleContainer::setParticles(std::vector<Particle> &particles)
{
	this->particles_ = particles;
}