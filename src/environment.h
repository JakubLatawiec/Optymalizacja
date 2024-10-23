#pragma once

#include <string>
#include "matrix.h"

void create_environment(std::string env_name);
void save_to_file(std::string filename, std::string& data);
void save_to_file(std::string filename, matrix& data);