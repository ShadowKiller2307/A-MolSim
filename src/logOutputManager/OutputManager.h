#pragma once
#include "particleContainers/ParticleContainer.h"
#include "outputWriter/VTKWriter.h"
#include "Particle.h"

#include <vector>

class ParticleContainer;

class OutputManager
{
private:
	outputWriter::VTKWriter writer;

public:
	bool outputFiles;
	bool outputBaseName;
	void plotParticles(const std::vector<Particle> &particles, const size_t iteration);
	void writeJSON(std::string &name, ParticleContainer &container);
	OutputManager();
};
