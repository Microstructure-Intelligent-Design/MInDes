#pragma once
#include <filesystem>
#include <vector>
#include <iostream>
#include "../base/whereami.h"
namespace pf {
	// get program's path, powered by "whereami.h", by Gregory Pakosz, https://github.com/gpakosz/whereami
	// use lambda to simplify expression.
	const std::filesystem::path program_path = std::filesystem::canonical([]() {
		int length = wai_getExecutablePath(nullptr, 0, nullptr);
		if (length <= 0) { 
			std::cout << "Failed to get executable path length";
			std::exit(0);
		};

		std::vector<char> wai_path(length + 1);
		wai_getExecutablePath(wai_path.data(), length, nullptr);

		return std::filesystem::path(wai_path.data());
		}()).remove_filename();
}