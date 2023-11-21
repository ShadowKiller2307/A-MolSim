#include "ParticleContainer.h"
#include "Particle.h"



spdlog::level::level_enum ParticleContainer::isDebug(){
    if(LogManager::getInstance().getLevel() == spdlog::level::debug){
        return spdlog::level::debug;
    }
    return spdlog::level::off;
}

spdlog::level::level_enum ParticleContainer::isInfo() {
    if(LogManager::getInstance().getLevel() == spdlog::level::info){
        return spdlog::level::info;
    }
    return spdlog::level::off;
}