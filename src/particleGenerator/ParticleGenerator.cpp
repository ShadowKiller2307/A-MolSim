#include "particleGenerator/ParticleGenerator.h"
#include "particleContainers/ParticleContainerDirSum.h"
#include "utils/MaxwellBoltzmannDistribution.h"
#include "utils/ArrayUtils.h"
#include "fileReader/FileReader.h"
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#include "Particle.h"
#include <nlohmann/json.hpp>
#include <fstream>
#include <vector>
#include "ParticleGenerator.h"
#include <iostream>

using json = nlohmann::json;

int particleGenerator::generateNumber_ = 0; // to make the compiler is happy, initialize the member here

void particleGenerator::instantiateCuboid(ParticleContainer **container, const std::array<double, 3> &llfc,
										  const std::array<unsigned int, 3> &particlesPerDimension, std::array<double, 3> &particleVelocity,
										  optionals optArgs)
{
	if (optArgs.generateNumber < 0)
	{
		optArgs.generateNumber = generateNumber_++;
	}
	for (size_t i = 0; i < particlesPerDimension[0]; ++i)
	{
		for (size_t j = 0; j < particlesPerDimension[1]; ++j)
		{
			for (size_t k = 0; k < particlesPerDimension[2]; ++k)
			{
				// TODO (ASK): should the mbVelocity be calculated for every particle instead?
				std::array<double, 3> mbVelocity = MaxwellBoltzmannDistribution::maxwellBoltzmannDistributedVelocity(optArgs.meanVelocity, 2); // TODO (ASK): Is the mean here the same as the average?
				std::array<double, 3> x_arg{i * optArgs.h + llfc[0], j * optArgs.h + llfc[1], k * optArgs.h + llfc[2]};
				std::array<double, 3> v_arg{particleVelocity + mbVelocity};			   // Calculate via the Brownian motion
				(*container)->add(x_arg, v_arg, optArgs.mass, optArgs.generateNumber); // using the add function in order to be able to add elements to all container implementations
			}
		}
	}
}

void particleGenerator::instantiateSphere(ParticleContainer **container, const std::array<double, 3> &center, const uint32_t &sphereRadius,
										  std::array<double, 3> &particleVelocity, optionals optArgs)
{
	if (optArgs.generateNumber < 0)
	{
		optArgs.generateNumber = generateNumber_++;
	}
	// TODO (ADD): Sphere instantiation
}

void particleGenerator::instantiateJSON(ParticleContainer **container, const std::string &path, optionals optArgs)
{
	std::ifstream f(path);
	json jsonContent = json::parse(f);
	json params = jsonContent["params"];
	if (!(*container))
	{
		double deltaT = optArgs.deltaT > 0 ? optArgs.deltaT : static_cast<double>(params["deltaT"]);
		double endTime = optArgs.endTime > 0 ? optArgs.endTime : static_cast<double>(params["endTime"]);
		(*container) = createContainer(deltaT, endTime);
	}
	// jsonContent["lol"] // returns null

	for (size_t i = 0; i < jsonContent["params"]["numParticles"]; i++)
	{
		auto j = jsonContent["particles"][i];
		std::array<double, 3UL> llfc = j["x"];
		std::array<double, 3UL> particleVelocity = j["v"];
		std::array<unsigned int, 3UL> particlePerDimension = j["N"];
		instantiateCuboid(container, llfc, particlePerDimension, particleVelocity, optArgs);
	}
}

void particleGenerator::instantiatePicture(ParticleContainer **container, const std::string &path, optionals optArgs)
{
	// gross C code, brace yourself
	if (optArgs.generateNumber < 0)
	{
		optArgs.generateNumber = generateNumber_++;
	}
	int width, height, bpp;
	char *charPath = new char[path.size()];
	strcpy(charPath, path.c_str());										// convert std::string to char*
	charPath[path.size() - 1] = '\0';									// explicitly set null terminator
	uint8_t *rgb_image = stbi_load(charPath, &width, &height, &bpp, 3); // load the image
	if (!rgb_image)
	{
		std::cout << "Error, could not load image!" << std::endl;
		exit(0);
	}
	std::vector<Particle> &particles = (*container)->getParticles();
	for (int i = 0; i < height; ++i)
	{
		for (int j = 0; j < width; j++)
		{
			int index = i * width * 3 + j * 3;
			uint8_t r = rgb_image[index + 0];
			uint8_t g = rgb_image[index + 1];
			uint8_t b = rgb_image[index + 2];
			if (r != 255 || g != 255 || b != 255)
			{
				std::array<double, 3> x_arg{j * optArgs.h, -i * optArgs.h, 0};
				std::array<double, 3> v_arg{0.0};
				if (r == 255)
				{
					v_arg = {0.0, r / 255.0 * -30.0, 0.0};
				}
				particles.emplace_back(x_arg, v_arg, optArgs.mass, optArgs.generateNumber);
			}
		}
	}
	stbi_image_free(rgb_image);
}

void particleGenerator::instantiateTxt(ParticleContainer **container, const std::string &path, optionals optArgs)
{
	if (!(*container))
	{
		double deltaT = optArgs.deltaT;
		double endTime = optArgs.endTime;
		(*container) = createContainer(deltaT, endTime);
	}
	FileReader fr = FileReader();
	auto actualPath = std::string("_.txt").compare(path) == 0 ? "../input/eingabe-sonne.txt" : path;
	char *charPath = new char[actualPath.size() + 1];
	strcpy(charPath, actualPath.c_str()); // covert std::string to char*
	charPath[actualPath.size()] = '\0';	  // explicitly set null terminator
	fr.readFile(container, charPath);
}

void particleGenerator::instantiateXML()
{
	// TODO (ADD): XML instantiation
}

ParticleContainer *particleGenerator::createContainer(double deltaT, double endTime)
{
	return new ParticleContainerDirSum{deltaT, endTime};
}
