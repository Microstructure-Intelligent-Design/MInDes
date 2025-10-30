#pragma once
#include "MainIterator.h"
int main(int argc, char* argv[]) {
	// run application
	pf::main_iterator::init_modules(argc, argv);
	pf::main_iterator::run();
}