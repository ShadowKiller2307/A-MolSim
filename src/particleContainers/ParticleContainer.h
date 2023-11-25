#include "utils/ArrayUtils.h"
#include "Particle.h"
#include <functional>
#include <vector>

class ParticleContainer
{
protected:
	double deltaT;
	std::vector<Particle> particles;

public:
	/**
	 * @brief iterate over every particle pair in the container and apply the lambda function
	 * @param f the lambda function
	 */
	virtual void iterOverPairs(const std::function<void(Particle &a, Particle &b)> &f);

	/**
	 * @brief calculate the velocity for all particles
	 * @param None
	 * @return void
	 */
	void calculateVelocity();

	/**
	 * @brief calculate the position for all particles
	 * @param None
	 * @return void
	 */
	void calculatePosition();

	/**
	 * @brief adds a particle to the container
	 *	@param x_arg position of the particle
	 *	@param v_arg velocity of the particle
	 *	@param mass mass of the particle
	 *	@param type typenumber of the particle
	 */
	virtual void add(const std::array<double, 3> &x_arg, const std::array<double, 3> &v_arg, double mass, int type);
};