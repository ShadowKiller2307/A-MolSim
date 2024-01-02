#pragma once
#include "particleContainers/ParticleContainerLinCel.h"
#include "outputWriter/VTKWriter.h"
#include "Particle.h"

#include <vector>

class ParticleContainerLinCel;

class OutputManager
{
private:
	outputWriter::VTKWriter writer;

public:
	bool outputFiles;
	bool outputBaseName;
	void plotParticles(const std::vector<Particle> &particles, const size_t iteration);
	void writeJSON(std::string &name, ParticleContainerLinCel *container);
	OutputManager();
};
