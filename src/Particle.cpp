/*
 * Particle.cpp
 *
 *  Created on: 23.02.2010
 *      Author: eckhardw
 */

#include "Particle.h"

#include <iostream>
#include "utils/ArrayUtils.h"

Particle::Particle(int type_arg)
{
  type = type_arg;
  // std::cout << "Particle generated!" << std::endl;
  f = {0., 0., 0.};
  old_f = {0., 0., 0.};
}

Particle::Particle(std::array<double, 3> x_arg, std::array<double, 3> v_arg, double m_arg, int type_arg) : x(x_arg), v(v_arg), m(m_arg), type(type_arg)
{
  f = {0., 0., 0.};
  old_f = {0., 0., 0.};
}

Particle::~Particle() {}

const std::array<double, 3> &Particle::getX() const { return x; }

const std::array<double, 3> &Particle::getV() const { return v; }

const std::array<double, 3> &Particle::getF() const { return f; }

const std::array<double, 3> &Particle::getOldF() const { return old_f; }

void Particle::setF(std::array<double, 3> &force)
{
  f = force;
}

void Particle::addF(std::array<double, 3> &force)
{
  f = f + force;
}

void Particle::setX(std::array<double, 3> &xPosition)
{
  x = xPosition;
}

void Particle::setV(std::array<double, 3> &velocity)
{
  v = velocity;
}

void Particle::setOldF(std::array<double, 3> &force)
{
  old_f = force;
}

double Particle::getM() const { return m; }

int Particle::getType() const { return type; }

std::string Particle::toString() const
{
  std::stringstream stream;
  stream << "Particle: X:" << x << " v: " << v << " f: " << f
         << " old_f: " << old_f << " type: " << type;
  return stream.str();
}

json Particle::toJSON() const
{
  nlohmann::json j;
  j["shape"] = "particle";
  j["x"] = x;
  j["v"] = v;
  j["f"] = f;
  j["old_f"] = old_f;
  j["type"] = type;
  j["m"] = m;
  return j;
}

bool Particle::operator==(Particle &other)
{
  return (x == other.x) and (v == other.v) and (f == other.f) and
         (type == other.type) and (m == other.m) and (old_f == other.old_f);
}

std::ostream &operator<<(std::ostream &stream, Particle &p)
{
  stream << p.toString();
  return stream;
}
