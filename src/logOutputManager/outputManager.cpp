#include "logOutputManager/outputManager.h"
#include "Particle.h"
#include <vector>

outputManager::outputManager()
{
	writer = outputWriter::VTKWriter();
}

void outputManager::plotParticles(const std::vector<Particle> &particles, const size_t iteration)
{
	writer.initializeOutput(particles.size());
	for (auto &p : particles)
	{
		writer.plotParticle(p);
	}
	writer.writeFile("../output/MD_vtk", iteration);
}

bool outputManager::getOutputFiles()
{
	return outputFiles;
}