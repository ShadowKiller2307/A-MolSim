#include "logOutputManager/outputManager.h"
#include "Particle.h"
#include <vector>

outputManager::outputManager()
{
	writer = outputWriter::VTKWriter();
}

void outputManager::plotParticles(const std::vector<Particle> &particles, const size_t iteration)
{
	auto w = outputWriter::VTKWriter();
	w.initializeOutput(particles.size());
	for (auto &p : particles)
	{
		w.plotParticle(p);
	}
	w.writeFile("../output/MD_vtk", iteration);
}

bool outputManager::getOutputFiles()
{
	return outputFiles;
}