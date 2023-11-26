#include "outputWriter/VTKWriter.h"
#include "Particle.h"
#include <vector>

class outputManager
{
private:
	outputWriter::VTKWriter writer;
	bool outputFiles = true;

public:
	void plotParticles(const std::vector<Particle> &particles, const size_t iteration);
	outputManager();
	bool getOutputFiles();
};
