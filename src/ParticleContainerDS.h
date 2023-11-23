/**
 * @brief Class to encapsulate a vector of particles, it has an attribute ForceCalulator
 * to define which calculation method should be used on the particles
 */
#pragma once

#include "Particle.h"
#include <vector>
#include <functional>
#include "ForceCalculator.h"
#include "LogManager.h"
#include "spdlog/spdlog.h"
#include "ParticleContainer.h"

class ForceCalculator; //needs to be defined so that the compiler doesn't throw an error



class ParticleContainerDS : public ParticleContainer
{

public:

    /**
//     * @brief calculate the velocity for all particles
//     * @param None
//     * @return void
//     */
//    void calculateVelocity();
//    /**
//     * @brief calculate the position for all particles
//     * @param None
//     * @return void
//     */
//    void calculatePosition();

    explicit ParticleContainerDS();

    /**
     * @brief sets the particles for the container
     * @param particles1 the particle vector to be set for the container
     * @return void
     *//*
    void setParticles(const std::vector<Particle> &particles1);

    *//***
     * @brief Checks if the level of the Log Manager is debug
     * @return Returns debug mode if true, off mode if false
     */

    void add(Particle &a) override;
/*
    *//**
    * @brief calculate the force for all particles
    * @param None
    * @return void
    *//*
    void calculateForces();*/
    void iterOverPairs(const std::function<void(Particle &a, Particle &b)> &forceLambda) override;
};
