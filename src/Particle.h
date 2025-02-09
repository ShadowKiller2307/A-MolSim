/*
 * Particle.h
 *
 *  Created on: 23.02.2010
 *      Author: eckhardw
 */

#pragma once
#include "nlohmann/json.hpp"
#include <array>
#include <string>

using json = nlohmann::json;
class Particle
{

private:
  /**
   * Position of the particle
   */
  std::array<double, 3> x;

  /**
   * Velocity of the particle
   */
  std::array<double, 3> v;

  /**
   * Force effective on this particle
   */
  std::array<double, 3> f;

  /**
   * Force which was effective on this particle
   */
  std::array<double, 3> old_f;

  /**
   * Mass of this particle
   */
  double m;

  /**
   * Type of the particle. Use it for whatever you want (e.g. to separate
   * molecules belonging to different bodies, matters, and so on)
   */
  int type;

  double sigma = 1;

  double epsilon = 5;

  /**
   * @brief stores pointers to the 8 neighbours of each particle in the following order:
   * 1 2 3
   * 4 p 5
   * 6 7 8
   * where p is the particle itself. If no neighbours is present at that position, a nullpointer is stored instead
   */
  std::array<std::shared_ptr<Particle>, 8> neighbours;

public:
  explicit Particle(int type = 0);

  Particle(
      // for visualization, we need always 3 coordinates
      // -> in case of 2d, we use only the first and the second
      std::array<double, 3> x_arg, std::array<double, 3> v_arg, double m_arg,
      int type = 0, int epsilon_arg = 5, int sigma_arg = 1);

  virtual ~Particle();

  const std::array<double, 3> &getX() const;

  const std::array<double, 3> &getV() const;

  const std::array<double, 3> &getF() const;

  const std::array<double, 3> &getOldF() const;

  double getSigma();

  double getEpsilon();

  void setF(std::array<double, 3> &force);

  void addF(std::array<double, 3> &force);

  void setX(std::array<double, 3> &xPosition);

  void addX(std::array<double, 3> &xPosition);

  void setV(std::array<double, 3> &velocity);

  void setType(int type);

  void setSigma(double sigma_arg);

  void setEpsilon(double epsilon_arg);

  void setOldF(std::array<double, 3> &force);

  double getM() const;

  int getType() const;

  bool operator==(Particle &other);

  std::string toString() const;

  json toJSON() const;
};

std::ostream &operator<<(std::ostream &stream, Particle &p);
