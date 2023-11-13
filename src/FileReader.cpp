/*
 * FileReader.cpp
 *
 *  Created on: 23.02.2010
 *      Author: eckhardw
 */

#include "FileReader.h"

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>

FileReader::FileReader() = default;

FileReader::~FileReader() = default;

void FileReader::readFile(std::vector<Particle> &particles, char *filename)
{
  std::array<double, 3> x;
  std::array<double, 3> v;
  double m;
  size_t numParticles = 0;

  std::ifstream inputFile(filename);
  std::string tmpString;

  if (inputFile.is_open())
  {

    getline(inputFile, tmpString);
    std::cout << "Read line: " << tmpString << std::endl;

    while (tmpString.empty() or tmpString[0] == '#')
    {
      getline(inputFile, tmpString);
      std::cout << "Read line: " << tmpString << std::endl;
    }

    std::istringstream numstream(tmpString);
    numstream >> numParticles;
    particles.resize(numParticles); // reserve enough space for all particles
    std::cout << "Reading " << numParticles << " particles." << std::endl;
    for (size_t i = 0; i < numParticles; i++)
    {
      getline(inputFile, tmpString);
      std::cout << "Read line: " << tmpString << std::endl;
      std::istringstream datastream(tmpString);

      for (auto &xj : x)
      {
        datastream >> xj;
      }
      if (datastream.eof())
      {
        std::cout
            << "Error reading file: eof reached unexpectedly reading from line "
            << i << std::endl;
        exit(-1);
      }
      for (auto &vj : v)
      {
        datastream >> vj;
      }
      if (datastream.eof())
      {
        std::cout
            << "Error reading file: eof reached unexpectedly reading from line "
            << i << std::endl;
        exit(-1);
      }
      datastream >> m;
      particles.at(i) = Particle(x, v, m, 0);
    }
  }
  else
  {
    std::cout << "Error: could not open file " << filename << std::endl;
    exit(-1);
  }
}
