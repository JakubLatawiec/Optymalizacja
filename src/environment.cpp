#include "environment.h"
#include "configuration.h"

#include <filesystem>
#include <sstream>
#include <fstream>

void create_environment(std::string env_name)
{
	//Tworzenie folderu jeœli ten nie istnieje
	if (!std::filesystem::exists(DATA_PATH + env_name))
		std::filesystem::create_directory(DATA_PATH + env_name);

	//Tworzenie stringa który jest czasem i dat¹ wykonywania programu
	auto now = std::chrono::system_clock::now();
	std::time_t now_time = std::chrono::system_clock::to_time_t(now);
	std::tm* local_time = std::localtime(&now_time);
	std::ostringstream date;
	date << std::put_time(local_time, "%Y-%m-%d_%H-%M-%S");

	//Tworzenie folderu z dat¹
	std::filesystem::create_directory(DATA_PATH + env_name + "/" + date.str());

	//Przypisanie folderu z plikami dla stworzonego œrodowiska
	FILE_PATH = DATA_PATH + env_name + "/" + date.str() + "/";
}

void save_to_file(std::string filename, std::string& data)
{
	std::ofstream file(FILE_PATH + filename);
	file << data;
	file.close();
}

void save_to_file(std::string filename, matrix& data)
{
	std::ofstream file(FILE_PATH + filename);
	file << data;
	file.close();
}
