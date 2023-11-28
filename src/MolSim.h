#pragma once
#include "spdlog/spdlog.h"





/**
 * @brief print the usage of the program to stdout
 * @param None
 * @return void
 */
void printHelp();
spdlog::level::level_enum mapIntToLevel(int programArgument);