#include "particleContainers/ParticleContainerDirSum.h"
#include "particleGenerator/ParticleGenerator.h"
#include "utils/MaxwellBoltzmannDistribution.h"
#include "logOutputManager/LogManager.h"
#include "fileReader/FileReader.h"
#include "ParticleGenerator.h"
#include "utils/ArrayUtils.h"
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#include "Particle.h"
#include <nlohmann/json.hpp>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <vector>

using json = nlohmann::json;

int particleGenerator::generateNumber_ = 0; // to make the compiler is happy, initialize the member here

void particleGenerator::instantiateCuboid(ParticleContainer **container, const std::array<double, 3> &llfc, const std::array<unsigned int, 3> &particlesPerDimension, std::array<double, 3> &particleVelocity, double h, double m, int type = -1)
{
	if (type < 0)
	{
		type = generateNumber_++;
	}
	for (size_t i = 0; i < particlesPerDimension[0]; ++i)
	{
		for (size_t j = 0; j < particlesPerDimension[1]; ++j)
		{
			for (size_t k = 0; k < particlesPerDimension[2]; ++k)
			{
				// TODO (ASK): should the mbVelocity be calculated for every particle instead?
				std::array<double, 3> mbVelocity = MaxwellBoltzmannDistribution::maxwellBoltzmannDistributedVelocity(meanVelocity_, 2); // TODO (ASK): Is the mean here the same as the average?
				std::array<double, 3> x_arg{i * h + llfc[0], j * h + llfc[1], k * h + llfc[2]};
				std::array<double, 3> v_arg{particleVelocity + mbVelocity}; // Calculate via the Brownian motion
				(*container)->add(x_arg, v_arg, m, type);					// using the add function in order to be able to add elements to all container implementations
			}
		}
	}
}

void particleGenerator::instantiateSphere(ParticleContainer **container, const std::array<double, 3> &center, const uint32_t &sphereRadius, std::array<double, 3> &particleVelocity, double h, double m, bool is2D, int type = -1)
{
	if (type < 0)
	{
		type = generateNumber_++;
	}
	//(*container)->add();

	// TODO (ADD): Sphere instantiation
}

void particleGenerator::instantiateJSON(ParticleContainer **container, const std::string &path, SimParams params)
{
	std::ifstream f(path);
	json jsonContent = json::parse(f);
	json JSONparams = jsonContent["params"];
	if (!(*container))
	{
		double deltaT = params.deltaT > 0 ? params.deltaT : static_cast<double>(JSONparams["deltaT"]);
		double endTime = params.endTime > 0 ? params.endTime : static_cast<double>(JSONparams["endTime"]);
		(*container) = createContainer(deltaT, endTime);
	}
	// jsonContent["lol"] // returns null

	for (size_t i = 0; i < jsonContent["params"]["numParticles"]; i++)
	{
		auto j = jsonContent["particles"][i];

		std::array<double, 3UL> pos = j["x"];
		std::array<double, 3UL> particleVelocity = j["v"];
		double h = j.contains("h") ? static_cast<double>(j["h"]) : h_;
		double m = j.contains("m") ? static_cast<double>(j["m"]) : m_;
		double type = j.contains("type") ? static_cast<double>(j["type"]) : generateNumber_++;
		if (j["shape"] == "cuboid")
		{
			std::array<unsigned int, 3UL> particlePerDimension = j["N"];
			instantiateCuboid(container, pos, particlePerDimension, particleVelocity, h, m, type);
		}
		else if (j["shape"] == "sphere")
		{
			uint32_t radius = j["R"];
			instantiateSphere(container, pos, radius, particleVelocity, h, m, type);
		}
		else
		{
			LogManager::errorLog("Type {} is not a valid type!", j["type"]);
			continue;
		}
	}
}

void particleGenerator::instantiatePicture(ParticleContainer **container, const std::string &path, SimParams params)
{
	// gross C code, brace yourself
	if (!(*container))
	{
		(*container) = createContainer(params.deltaT, params.endTime);
	}
	int width, height, bpp;
	char *charPath = new char[path.size()];
	strcpy(charPath, path.c_str());										// convert std::string to char*
	charPath[path.size() - 1] = '\0';									// explicitly set null terminator
	uint8_t *rgb_image = stbi_load(charPath, &width, &height, &bpp, 3); // load the image
	if (!rgb_image)
	{
		LogManager::errorLog("Error, could not load image!");
		exit(1);
	}
	std::vector<std::array<uint8_t, 3>> lookup;
	for (int i = 0; i < height; ++i)
	{
		for (int j = 0; j < width; j++)
		{
			int index = i * width * 3 + j * 3;
			uint8_t r = rgb_image[index + 0]; // x-movement
			uint8_t g = rgb_image[index + 1]; // y-movement
			uint8_t b = rgb_image[index + 2]; // mass
			if (!(r == 255 && g == 255 && b == 255))
			{
				int8_t sr = r <= INT8_MAX ? static_cast<int8_t>(r) : static_cast<int>(r - INT8_MIN) + INT8_MIN;
				int8_t sg = g <= INT8_MAX ? static_cast<int8_t>(g) : static_cast<int>(g - INT8_MIN) + INT8_MIN;
				std::array<double, 3> x_arg{j * h_, -i * h_, 0};
				std::array<double, 3> v_arg = {sr / 127.0 * 40.0, sg / 127.0 * 40.0, 0.0};
				std::array<uint8_t, 3> val = {r, g, b};
				auto index = std::find(lookup.begin(), lookup.end(), val);
				int type;
				if (index == lookup.end())
				{
					lookup.push_back(val);
					type = lookup.size() - 1;
				}
				else
				{
					type = std::distance(lookup.begin(), index);
				}
				(*container)->add(x_arg, v_arg, b, type);
			}
		}
	}
	stbi_image_free(rgb_image);
}

void particleGenerator::instantiateTxt(ParticleContainer **container, const std::string &path, SimParams params)
{
	if (!(*container))
	{
		(*container) = createContainer(params.deltaT, params.endTime);
	}
	FileReader fr = FileReader();
	auto actualPath = std::string("_.txt").compare(path) == 0 ? "../input/eingabe-sonne.txt" : path;
	char *charPath = new char[actualPath.size() + 1];
	strcpy(charPath, actualPath.c_str()); // covert std::string to char*
	charPath[actualPath.size()] = '\0';	  // explicitly set null terminator
	fr.readFile(container, charPath);
}

void particleGenerator::instantiateXML(ParticleContainer **container, const std::string &path, SimParams params)
{
	// TODO (ADD): XML instantiation
}

ParticleContainer *particleGenerator::createContainer(double deltaT, double endTime)
{
	return new ParticleContainerDirSum{deltaT, endTime};
}
