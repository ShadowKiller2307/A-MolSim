/**
 * @brief Class which reads the input for the calculation, plots the particles
 * and calls the calculation functions on its attribute particleContainer
 */
#include "ForceCalculator.h"

#pragma once

/**
 * @brief calculate the force for all particles
 * @param None
 * @return void
 */
void calculateF();

/**
 * @brief calculate the position for all particles
 * @param None
 * @return void
 */
void calculateX();

/**
 * @brief calculate the velocity for all particles
 * @param None
 * @return void
 */
void calculateV();

/**
 * @brief plot the particles to a xyz-file
 * @param iteration amount of iteration over the particles
 * @return void
 */
void plotParticles(int iteration);

/**
 * @brief print the usage of the program to stdout
 * @param None
 * @return void
 */
void printHelp();

void SetForceCalculator(int mode);
