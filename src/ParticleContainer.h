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

class ForceCalculator; //needs to be defined so that the compiler doesn't throw an error



class ParticleContainer
{

private:
    std::vector<Particle> particles;
    ForceCalculator *forceCalculator;
    double deltaTTwo;


public:
    std::vector<Particle> &getParticles();

    /**
     * @brief calculate the force for all particles
     * @param None
     * @return void
     */
    void calculateForces();

    /**
     * @brief iterate over every particle pair in the container and apply the lambda function
     * @param f the lambda function
     */
    void iterOverPairs(const std::function<void(Particle &a, Particle &b)> &f);

    /**
     * @brief calculate the velocity for all particles
     * @param None
     * @return void
     */
    void calculateVelocity();
    /**
     * @brief calculate the position for all particles
     * @param None
     * @return void
     */
    void calculatePosition();

    explicit ParticleContainer();

    /**
     * @brief sets the particles for the container
     * @param particles1 the particle vector to be set for the container
     * @return void
     */
    void setParticles(const std::vector<Particle> &particles1);
    /**
     * @brief sets the ForceCalculator for the container
     * @param mode represent the different approaches for the force calculation
     * @return void
     */
    void setForceCalculator(int mode);

    /**
     * @brief sets the deltaT for the container
     * @param deltaT the deltaT passed by the user or the default value
     * @return void
     */
    void setDeltaTTwo(double deltaT);

    /***
     * @brief Checks if the level of the Log Manager is debug
     * @return Returns debug mode if true, off mode if false
     */
    static spdlog::level::level_enum isDebug();

    double getDeltaTwo();

};
