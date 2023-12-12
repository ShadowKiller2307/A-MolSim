#include "Periodic.h"


//Periodic::Periodic(double position, int direction, std::function<void(Particle &, Particle &)> &forceLambda) : BoundaryCondition(position, direction, forceLambda) {}

bool Periodic::affectsForce() {
    return true;
}

/*
void Periodic::applyHaloCondition(Particle &a) {

    std::array<double,3> velocity = a.getV();
    double velocityInX = velocity.at(0);
    double velocityInY = velocity.at(1);

    if(velocityInX > 0 && velocityInY == 0){ //the particle moves right

    }
    else if(velocityInX < 0 && velocityInY == 0){ //the particle moves left

    }
    else if(velocityInX == 0 && velocityInY >0){ //the particle moves up

    }
    else if(velocityInX == 0 && velocityInY < 0){ //the particle moves down

    }


    else if(velocityInX > 0 && velocityInY>0){ // the particle moves up right

    }
    else if(velocityInX<0 && velocityInY > 0 ) { //the particle moves up left

    }
    else if(velocityInX>0 &&velocityInY < 0){ //the particle moves down right

    }
    else if(velocityInX < 0&&velocityInY < 0){ //the particle moves down left

    }

}*/
/**
 * order:
 * left plane
 * right plane
 * lower plane
 * upper plane
 * back plane
 * front plane
 */
 /*
  * first mirror the x-Position
  * then the y position if it is still in the Halo zone
  * then the z position if it is still in the Halo zone
  */
void Periodic::applyHaloCondition(Particle &a) {
    //direction and position from the parent class can also be used for the halo cells
    // calculate the new position for particle {a} according to the periodic boundary
    //directions: x == 0
    // y == 1
    // z == 2
    double diffABound; // difference between the position of particle a and the position of the domain boundary
     if (direction_ == 0) { //
         if (position_ == 0) { // in the halo cell next to the left domain boundary
             double newX = a.getX().at(0) + domainSize_[0];
             std::array<double, 3> newPos {newX,a.getX().at(1),a.getX().at(2)};
             a.setX(newPos);
         }
         else { // in the halo cell next to the right domain boundary
             double newX = a.getX().at(0) - domainSize_[0];
             std::array<double, 3> newPos {newX,a.getX().at(1),a.getX().at(2)};
             a.setX(newPos);
         }
     }
     if (direction_ == 1) { //
         if (position_ == 0) {  // in the halo cell below the lower domain boundary
             double newY = a.getX().at(1) + domainSize_[1];
             std::array<double, 3> newPos {a.getX().at(0), newY, a.getX().at(2)};
             a.setX(newPos);
         }
         else { // in the halo cell over the top domain boundary
             double newY = a.getX().at(1) - domainSize_[1];
             std::array<double, 3> newPos {a.getX().at(0), newY, a.getX().at(2)};
             a.setX(newPos);
         }
     }
     if (direction_ == 2) { //
         if (position_ == 0) {  // in the halo cell below the lower domain boundary
             double newZ = a.getX().at(2) + domainSize_[2];
             std::array<double, 3> newPos {a.getX().at(0), a.getX().at(1), newZ};
             a.setX(newPos);
         }
         else { // in the halo cell over the top domain boundary
             double newZ = a.getX().at(2) - domainSize_[2];
             std::array<double, 3> newPos {a.getX().at(0), a.getX().at(1), newZ};
             a.setX(newPos);
         }
     }
}



Periodic::Periodic(double position, int direction, std::function<void(Particle &, Particle &)> &forceLambda, std::array<double, 3> domainSize)
        : BoundaryCondition(position, direction, forceLambda, domainSize) {

}

/**
 * @brief forces from the particles on the other side of the boundary should ba applied
 * @param a
 */
void Periodic::applyBoundCondition(Particle &a, std::vector<cell> &particlesOtherSide) {
    for (BoundaryCondition::cell &current : particlesOtherSide) {
        for (Particle &b : current) {
            forceLambda_(a, b);
        }
    }
}
