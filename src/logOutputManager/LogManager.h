#pragma once
#include "spdlog/spdlog.h"
#include <spdlog/sinks/rotating_file_sink.h>

class LogManager
{
public:
    /// Returns the Singleton instance of the LogManager
    static LogManager &getInstance();

    /// Returns the logger of the LogManager
    std::shared_ptr<spdlog::logger> getLogger();

    /// Sets the log level of the LogManager
    void setLogLevel(spdlog::level::level_enum level);

    /// Returns the log level of the LogManager
    spdlog::level::level_enum getLevel();

    /// Deletes the copy constructor of the class
    LogManager(const LogManager &) = delete;

    /// Deletes the assignment operator
    LogManager &operator=(const LogManager &) = delete;

private:
    std::shared_ptr<spdlog::logger> logger;
    const int maxSize = 5 * 1024 * 1024;
    const int maxFiles = 4;
    spdlog::level::level_enum logLevel;

    LogManager();
    ~LogManager() = default;
};
