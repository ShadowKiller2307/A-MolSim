#pragma once

/// calculate the force for all particles
void calculateF();

/// calculate the position for all particles
void calculateX();

/// calculate the position for all particles
void calculateV();

/// plot the particles to a xyz-file
void plotParticles(int iteration);

/// print the usage of the program to stdout
void printHelp();