#include "particleContainers/ParticleContainerLinCel.h"
#include "particleGenerator/ParticleGenerator.h"
#include "utils/MaxwellBoltzmannDistribution.h"
#include "logOutputManager/LogManager.h"
#include "xmlSchema/XMLReader.h"
#include "ParticleGenerator.h"
#include "utils/ArrayUtils.h"
#include "forces/LennJon.h"

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

void particleGenerator::instantiateCuboid(ParticleContainerLinCel **container, const std::array<double, 3> &llfc,
										  const std::array<unsigned int, 3> &particlesPerDimension,
										  std::array<double, 3> &particleVelocity, double h, double m, int type, double epsilon, double sigma, double initT = 0.0)
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
				initT = fabs(initT);
				std::array<double, 3> mbVelocity;
				if (initT >= 0.0001)
				{
					double initTVelocity = std::sqrt((initT / m));
					mbVelocity = MaxwellBoltzmannDistribution::maxwellBoltzmannDistributedVelocity(
						initTVelocity, 2);
				}
				else
				{
					mbVelocity = MaxwellBoltzmannDistribution::maxwellBoltzmannDistributedVelocity(
						meanVelocity_, 2); // TODO (ASK): Is the mean here the same as the average?
				}
				std::array<double, 3> x_arg{i * h + llfc[0], j * h + llfc[1], k * h + llfc[2]};
				std::array<double, 3> v_arg{particleVelocity + mbVelocity}; // Calculate via the Brownian motion
				auto p = (*container)->add(x_arg, v_arg, m, type);			// using the add function in order to be able to add elements to all container implementations
				if (p)
				{
					p->setEpsilon(epsilon);
					p->setEpsilon(sigma);
				}
			}
		}
	}
}

void particleGenerator::instantiateSphere(ParticleContainerLinCel **container, const std::array<double, 3> &center,
										  const int32_t &sphereRadius, const std::array<double, 3> &particleVelocity,
										  double h, double m, bool is2D, int type, double epsilon, double sigma, double initT = 0.0)
{
	if (type < 0)
	{
		type = generateNumber_++;
	}
	if (is2D)
	{
		for (int32_t i = -sphereRadius + 1; i < sphereRadius; ++i)
		{
			for (int32_t j = -sphereRadius + 1; j < sphereRadius; ++j)
			{
				if (std::sqrt(i * i + j * j) <= sphereRadius)
				{
					initT = fabs(initT);
					std::array<double, 3> mbVelocity;
					if (initT >= 0.0001)
					{
						double initTVelocity = std::sqrt((initT / m));
						mbVelocity = MaxwellBoltzmannDistribution::maxwellBoltzmannDistributedVelocity(
							initTVelocity, 2);
					}
					else
					{
						mbVelocity = MaxwellBoltzmannDistribution::maxwellBoltzmannDistributedVelocity(
							meanVelocity_, 2); // TODO (ASK): Is the mean here the same as the average?
					}
					std::array<double, 3> x_arg{i * h + center[0], j * h + center[1], 0};
					std::array<double, 3> v_arg{particleVelocity + mbVelocity}; // Calculate via the Brownian motion
					auto p = (*container)->add(x_arg, v_arg, m, type);			// using the add function in order to be able to add elements to all container implementations
					if (p)
					{
						p->setEpsilon(epsilon);
						p->setEpsilon(sigma);
					}
				}
			}
		}
	}
	else
	{
		for (int32_t i = -sphereRadius + 1; i < sphereRadius; ++i)
		{
			for (int32_t j = -sphereRadius + 1; j < sphereRadius; ++j)
			{
				for (size_t k = -sphereRadius + 1; i < sphereRadius; ++k)
				{
					if (std::sqrt(i * i + j * j + k * k) <= sphereRadius)
					{
						std::array<double, 3> mbVelocity = MaxwellBoltzmannDistribution::maxwellBoltzmannDistributedVelocity(
							meanVelocity_, 3); // TODO (ASK): Is the mean here the same as the average?
						std::array<double, 3> x_arg{i * h + center[0], j * h + center[1], k * h + center[2]};
						std::array<double, 3> v_arg{particleVelocity + mbVelocity}; // Calculate via the Brownian motion
						auto p = (*container)->add(x_arg, v_arg, m, type);			// using the add function in order to be able to add elements to all container implementations
						if (p)
						{
							p->setEpsilon(epsilon);
							p->setEpsilon(sigma);
						}
					}
				}
			}
		}
	}
}

void particleGenerator::instantiateJSON(ParticleContainerLinCel **container, const std::string &path, Force &force,
										SimParams params)
{
	std::ifstream f(path);
	json jsonContent = json::parse(f);
	json JSONparams = jsonContent["params"];
	if (!(*container))
	{
		double deltaT = params.deltaT > 0 ? params.deltaT : static_cast<double>(JSONparams["deltaT"]);
		double endTime = params.endTime > 0 ? params.endTime : static_cast<double>(JSONparams["endTime"]);
		std::string containerType = params.containerType != "" ? params.containerType : static_cast<std::string>(JSONparams["containerType"]);

		int writeFrequency;
		if (params.writeFrequency > 0)
		{
			writeFrequency = params.writeFrequency;
		}
		else if (JSONparams.contains("writeFrequency"))
		{
			writeFrequency = static_cast<int>(JSONparams["writeFrequency"]);
		}
		else
		{
			writeFrequency = 100;
		}

		if (containerType == "LinCel")
		{
			std::string bounds =
				params.boundaries != "" ? params.boundaries : static_cast<std::string>(JSONparams["boundaries"]);
			std::array<double, 3> domainSize;
			for (int i = 0; i < 3; i++)
			{
				domainSize[i] = params.domainSize[i] > 0 ? params.domainSize[i]
														 : static_cast<std::array<double, 3>>(JSONparams["domainSize"])[i];
			}
			double cutoffRadius = params.cutoffRadius > 0 ? params.cutoffRadius : static_cast<double>(JSONparams["cutoffRadius"]);
			if (!JSONparams.contains("nThermo") || !JSONparams.contains("tInit"))
			{
				(*container) = new ParticleContainerLinCel(deltaT, endTime, writeFrequency, domainSize, bounds,
														   cutoffRadius, false);
			}
			else
			{
				int nThermostat = static_cast<int>(JSONparams["nThermo"]);
				double tInit = static_cast<double>(JSONparams["tInit"]);
				double tempTarget = JSONparams.contains("tempTarget") ? static_cast<double>(JSONparams["tempTarget"]) : tInit;

				if (JSONparams.contains("gGrav"))
				{
					(*container) = new ParticleContainerLinCel(deltaT, endTime, writeFrequency, domainSize, bounds,
															   cutoffRadius, true, nThermostat, true, tInit, tempTarget, 0.5, static_cast<double>(JSONparams["gGrav"]));
				}
				else
				{
					(*container) = new ParticleContainerLinCel(deltaT, endTime, writeFrequency, domainSize, bounds,
															   cutoffRadius, true, nThermostat, true, tInit, tempTarget, 0.5);
				}
			}
		}
		else
		{
			LogManager::errorLog("Contianer type \"{}\" is unknown!", containerType);
			exit(1);
		}
	}

	for (size_t i = 0; i < jsonContent["params"]["numParticles"]; i++)
	{
		auto j = jsonContent["particles"][i];

		std::array<double, 3UL> pos = j["x"];
		std::array<double, 3UL> particleVelocity = j["v"];
		double h = j.contains("h") ? static_cast<double>(j["h"]) : h_;
		double m = j.contains("m") ? static_cast<double>(j["m"]) : m_;
		int type = j.contains("type") ? static_cast<int>(j["type"]) : generateNumber_++;
		double sigm = j.contains("sigma") ? static_cast<double>(j["sigma"]) : 1;
		double epsi = j.contains("epsilon") ? static_cast<double>(j["epsilon"]) : 5;
		if (j["shape"] == "cuboid")
		{
			std::array<unsigned int, 3UL> particlePerDimension = j["N"];
			instantiateCuboid(container, pos, particlePerDimension, particleVelocity, h, m, type, epsi, sigm);
		}
		else if (j["shape"] == "sphere")
		{
			uint32_t radius = j["R"];
			instantiateSphere(container, pos, radius, particleVelocity, h, m, true, type, epsi, sigm);
		}
		else if (j["shape"] == "particle")
		{
			auto particle = Particle(j["x"], j["v"], j["m"], j["type"]);
			auto f = static_cast<std::array<double, 3>>(j["old_f"]);
			particle.setOldF(f);
			auto old_f = static_cast<std::array<double, 3>>(j["f"]);
			particle.setF(old_f);
			particle.setEpsilon(epsi);
			particle.setSigma(sigm);
			(*container)->addCompleteParticle(particle);
		}
		else
		{
			if (!j["type"])
			{
				LogManager::errorLog("You forgot to specifie a type!");
			}
			else
			{
				LogManager::errorLog("Type {} is not a valid type!", j["type"]);
			}
			continue;
		}
	}
}

void particleGenerator::instantiateXML(ParticleContainerLinCel **container, std::string &path, Force &force, SimParams clArgs)
{

	XMLReader xmlReader(path);
	xmlReader.extractSimulationParameters();
	xmlReader.extractCuboid();
	xmlReader.extractSphere();

	SimulationConstructor simConst = xmlReader.getSimulationConstructor();
	auto cuboidConst = xmlReader.getCuboidConstructors();
	auto sphereConst = xmlReader.getSphereConstructors();

	LogManager::warnLog("Now reading simulation parameters. Commandline arguments have a higher priority if they exist.\n");

	double delta_t = clArgs.deltaT > 0 ? clArgs.deltaT : simConst.getDelta_t();
	double t_end = clArgs.endTime > 0 ? clArgs.endTime : simConst.getT_end();
	std::string containerT = !clArgs.containerType.empty() ? clArgs.containerType : simConst.getContainerType();
	std::string boundaries = !clArgs.boundaries.empty() ? clArgs.boundaries : simConst.getBoundaries();
	int writeFrequency = clArgs.writeFrequency > 0 ? clArgs.writeFrequency : simConst.getWriteFrequency();
	double cutOffRadius = clArgs.cutoffRadius > 0 ? clArgs.cutoffRadius : simConst.getCutOffRadius();

	double gGrav = simConst.getGrav();
	bool useThermostat = simConst.isUseThermostat();
	double initT = simConst.getInitialTemp();
	unsigned int nT = simConst.getNThermostat();
	bool isGradual = simConst.getIsGradual();
	double tempTarget = simConst.getTempTarget();
	double maxDiff = simConst.getMaxDiff();

	std::array<double, 3> domainSize{};
	for (int i = 0; i < 3; i++)
	{
		domainSize[i] = clArgs.domainSize.at(i) > 0 ? clArgs.domainSize.at(i) : simConst.getDomainSize().at(i);
	}

	if (containerT == "LinCel")
	{
		if (xmlReader.isThermostatPresent())
		{
			(*container) = new ParticleContainerLinCel(delta_t, t_end, writeFrequency, domainSize, boundaries,
													   cutOffRadius, useThermostat, nT, isGradual, initT,
													   tempTarget, maxDiff, gGrav);
		}
		else
		{
			(*container) = new ParticleContainerLinCel(delta_t, t_end, writeFrequency, domainSize, boundaries,
													   cutOffRadius);
		}

		for (auto &cuboid : cuboidConst)
		{
			double h = h_;
			if (cuboid.getH() != -1)
			{
				h = cuboid.getH();
			}
			double epsilon = 5;
			double sigma = 1;
			if (cuboid.getSigma() != -1)
			{
				sigma = cuboid.getSigma();
			}
			if (cuboid.getEpsilon() != -1)
			{
				epsilon = cuboid.getEpsilon();
			}
			instantiateCuboid(container, cuboid.getLlfc(), cuboid.getParticlesPerDimension(),
							  const_cast<std::array<double, 3> &>(cuboid.getParticleVelocity()),
							  h, cuboid.getMass(), cuboid.getType(), epsilon, sigma);
			LogManager::debugLog("Instantiated a cuboid from xml\n");
		}
		for (auto &sphere : sphereConst)
		{
			double h = h_;
			if (sphere.getDistance() != -1)
			{
				h = sphere.getDistance();
			}
			double epsilon = 5;
			double sigma = 1;
			if (sphere.getSigma() != -1)
			{
				sigma = sphere.getSigma();
			}
			if (sphere.getEpsilon() != -1)
			{
				epsilon = sphere.getEpsilon();
			}
			int type = sphere.getType();
			instantiateSphere(container, sphere.getCenterCoordinates(), sphere.getRadius(), sphere.getInitialVelocity(),
							  h, sphere.getMass(), true, type, epsilon, sigma);
			LogManager::debugLog("Instantiated a sphere from xml\n");
		}
	}
	else
	{
		LogManager::errorLog("Type {} is not a valid type!", containerT);
	}
}
