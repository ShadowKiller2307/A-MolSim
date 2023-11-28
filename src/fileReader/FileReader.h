/*
 * FileReader.h
 *
 *  Created on: 23.02.2010
 *      Author: eckhardw
 */

#pragma once

#include "particleContainers/ParticleContainer.h"
#include "Particle.h"

#include <vector>

class FileReader
{

public:
  FileReader();
  virtual ~FileReader();

  void readFile(ParticleContainer **container, char *filename);
};
