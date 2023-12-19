#include "particleContainers/ParticleContainer.h"
#include "particleGenerator/ParticleGenerator.h"
#include "forces/GravPot.h"
#include "forces/LennJon.h"
#include "MolSim.h"
#include <iostream>
#include <fstream>
#include <getopt.h>
#include <chrono>
#include "logOutputManager/LogManager.h"
#include "particleContainers/ParticleContainerLinCel.h"

int main(int argc, char *const argv[])
{
	// std::cout << "bis hier ok1\n";
	LogManager &logManager = LogManager::getInstance();
	logManager.setLogLevel(spdlog::level::info); // standard default level

	auto params = SimParams();
	bool writeToJSON = false;
	std::string outName;
	ParticleContainer *container = nullptr;
	option longOpts[] = {
		{"deltaT", required_argument, nullptr, 'd'},
		{"endTime", required_argument, nullptr, 'e'},
		{"help", no_argument, nullptr, 'h'},
		{"logLevel", required_argument, nullptr, 'l'},
		{"outName", required_argument, nullptr, 'o'},
		{nullptr, 0, nullptr, 0}};

	int longOptsIndex = 0;
	while (true)
	{
		int c = getopt_long(argc, argv, "d:e:hl:w", longOpts, &longOptsIndex);
		if (c == -1)
		{
			break;
		}
		switch (c)
		{
		case 'd':
			params.deltaT = std::stod(optarg);
			break;
		case 'e':
			params.endTime = std::stod(optarg);
			break;
		case 'h':
			printHelp();
			exit(0);
		case 'l':
		{
			int temp = std::stoi(optarg);
			if (temp >= 0 && temp < 5)
			{
				logManager.setLogLevel(mapIntToLevel(temp));
			}
			else
			{
				LogManager::errorLog("Invalid Log Level <{}>!", temp);
			}
			break;
		}
		case 'o':
			outName = optarg;
			break;
		case 'w':
			writeToJSON = true;
			break;
		default:
			LogManager::errorLog("Error parsing arguments. Maybe you gave one that isn't recognized.\n");
			break;
		}
	}

	if (optind >= argc)
	{
		LogManager::errorLog("Input missing as an argument, aborting\n");
		printHelp();
		return EXIT_FAILURE;
	}

	// read the inputfile depending on the file ending
	auto str = std::string(argv[optind]);
	int dot = str.find_last_of(".");
	auto ending = str.substr(dot + 1);
	std::string path;

	// C++ can't do switch statements on strings😔, this does the job, let's not overcomplicate things
	if (ending == "xml")
	{
		path = std::string("_.xml").compare(argv[optind]) == 0 ? "../input/default.xml" : argv[optind];
		auto force = LennJon(5.0, 1.0);
		particleGenerator::instantiateXML(&container, path, force, params);
	}
	else if (ending == "json")
	{
		path = std::string("_.json").compare(argv[optind]) == 0 ? "../input/collision.json" : argv[optind];
		auto force = LennJon(5.0, 1.0);
		particleGenerator::instantiateJSON(&container, path, force, params);
	}
	/*else if (ending == "png")
	{
		path = std::string("_.png").compare(argv[optind]) == 0 ? "../input/Cool MolSim.png" : argv[optind];
		auto force = LennJon(5.0, 1.0);
		particleGenerator::instantiatePicture(&container, path, force, SimParams{.deltaT = 0.0002, .endTime = 5});
	}*/
	else if (ending == "txt")
	{
		auto force = GravPot();
		particleGenerator::instantiateTxt(&container, argv[optind], force, SimParams{.deltaT = 0.014, .endTime = 1000, .containerType = "DirSum"});
	}
	if (writeToJSON)
	{
		int slash = str.find_last_of("/");
		auto path = str.substr(0, slash + 1);
		auto newName = path + "generated_" + str.substr(slash + 1, str.size() - slash - 6) + ".json";
		container->writeJSON(newName);
		return 0;
	}
	if (true)
	{
		ParticleContainerLinCel *lincelContainer = dynamic_cast<ParticleContainerLinCel *>(container);
		lincelContainer->simulateParticles();
	}
	// std::cout << "bis hier ok3\n";
	delete container;
	return 0;
}

void printHelp()
{
	std::ifstream file("../help.txt");
	if (file.is_open())
	{
		std::cout << file.rdbuf();
	}
}
spdlog::level::level_enum mapIntToLevel(int programArgument)
{
	switch (programArgument)
	{
	case 0:
		return spdlog::level::level_enum::err;
	case 1:
		return spdlog::level::level_enum::warn;
	case 2:
		return spdlog::level::level_enum::info;
	case 3:
		return spdlog::level::level_enum::debug;
	default:
		return spdlog::level::level_enum::off;
	}
}
