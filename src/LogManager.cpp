
#include "LogManager.h"

LogManager &LogManager::getInstance()
{
    static LogManager instance;
    return instance;
}

LogManager::LogManager()
{

    logLevel = spdlog::level::info;
    spdlog::set_pattern("[%Y-%m-%d %H:%M] [%l] %v");
    logger = spdlog::rotating_logger_mt("fileLogger", "../logs/log.txt",
                                        maxSize, maxFiles, true);
    logger->set_level(this->logLevel);
}

std::shared_ptr<spdlog::logger> LogManager::getLogger()
{
    return logger;
}

spdlog::level::level_enum LogManager::getLevel()
{
    return this->logLevel;
}

void LogManager::setLogLevel(spdlog::level::level_enum level)
{
    this->logLevel = level;
    this->logger->set_level(level);
}
