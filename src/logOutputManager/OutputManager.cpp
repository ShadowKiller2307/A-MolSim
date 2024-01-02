#include "logOutputManager/OutputManager.h"
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

void OutputManager::plotParticles(const std::vector<Particle> &particles, const size_t iteration)
{
	auto w = outputWriter::VTKWriter();
	w.initializeOutput(particles.size());
	for (auto &p : particles)
	{
		w.plotParticle(p);
	}
	w.writeFile("../output/MD_vtk", iteration);
}

void OutputManager::writeJSON(std::string &name, ParticleContainerLinCel *container)
{
	std::string containerType;
	if (dynamic_cast<ParticleContainerLinCel *>(container) != nullptr)
	{
		containerType = "LinCel";
	}
	else
	{
		containerType = "DirSum";
	}
	json output = {{"params", {{"numParticles", container->getAmountOfParticles()}, {"deltaT", container->getDeltaT()}, {"endTime", container->getEndTime()}, {"containerType", containerType}}}, {"particles", json::array()}};
	if (containerType == "LinCel")
	{
		std::string conditionsString = "";
		auto newPointer = dynamic_cast<ParticleContainerLinCel *>(container);
		auto conditions = newPointer->getConditions();
		for (auto &c : conditions)
		{
			switch (c)
			{
			case BoundaryCondition::Outflow:
				conditionsString += "o";
				break;
			case BoundaryCondition::Reflecting:
				conditionsString += "r";
				break;
			case BoundaryCondition::Periodic:
				conditionsString += "p";
				break;
			}
		}
		output["params"]["boundaries"] = conditionsString;
		output["params"]["domainSize"] = newPointer->getDomainSize();
		output["params"]["cutoffRadius"] = newPointer->getCutOffRadius();
	}
	container->iterOverAllParticles([&](std::vector<Particle>::iterator it)
									{ output["particles"].push_back(it->toJSON()); });
	std::ofstream o(name);
	o << std::setw(4) << output << std::endl;
}
