#include "logOutputManager/OutputManager.h"
#include "particleContainers/ParticleContainer.h"
#include "Particle.h"

#include <nlohmann/json.hpp>
#include <iostream>
#include <fstream>
#include <vector>

using json = nlohmann::json;

OutputManager::OutputManager()
{
	writer = outputWriter::VTKWriter();
	outputFiles = true;
}

void OutputManager::plotParticles(const std::vector<Particle *> &particles, const size_t iteration)
{
	auto w = outputWriter::VTKWriter();
	w.initializeOutput(particles.size());
	for (auto &p : particles)
	{
		w.plotParticle(*p);
	}
	w.writeFile("../output/MD_vtk", iteration);
}

void OutputManager::writeJSON(std::string &name, ParticleContainer &container)
{
	json output = {
		{"params", {{"numParticles", container.getParticles().size()}, {"deltaT", container.getDeltaT()}, {"endTime", container.getEndTime()}}},
		{"particles", json::array()}};
	// TODO: Add Particles to JSON
	std::ofstream o(name);
	o << std::setw(4) << output << std::endl;
}
