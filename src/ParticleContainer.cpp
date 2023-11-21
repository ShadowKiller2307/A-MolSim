#include "ParticleContainer.h"
#include "Particle.h"



spdlog::level::level_enum ParticleContainer::isDebug(){
    if(LogManager::getInstance().getLevel() == spdlog::level::debug){
        return spdlog::level::debug;
    }
    return spdlog::level::off;
}
