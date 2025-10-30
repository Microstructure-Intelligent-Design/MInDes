#pragma once
#include <iostream>
#include <string>
#include "SimuInfo.h"
#ifdef _WIN32
#include <windows.h>
namespace pf {
	inline void selectPMFile(SimuInfo& _simu_info) {
		char szBuffer[MAX_PATH] = { 0 };
		char myDocumentsPath[MAX_PATH] = { 0 };

		OPENFILENAMEA file = { 0 };
		file.lStructSize = sizeof(file);
		file.lpstrFilter = "MInDes Input File(*.mindes)\0*.mindes\0All Documents(*.*)\0*.*\0";
		file.lpstrTitle = "Select a MInDes Input File";
		file.lpstrFile = szBuffer;
		file.nMaxFile = MAX_PATH;
		file.Flags = OFN_PATHMUSTEXIST | OFN_FILEMUSTEXIST | OFN_EXPLORER | OFN_NOCHANGEDIR;

		if (GetOpenFileNameA(&file)) {
			_simu_info.simu_path = file.lpstrFile;
			_simu_info.is_simu_ready = true;
		}
		else {
			DWORD dwError = CommDlgExtendedError();
			if (dwError != 0) {
				pf::printf_color_on_control("Error opening file dialog, error code: " + dwError, 31);
				printf("\n");
				exit(dwError);
			}
			else {
				pf::printf_color_on_control("File selection cancelled by user.\n", 31);
				_simu_info.is_simu_ready = false;
			}
		}
	}
	inline void selectDataFile(std::string& file_path) {
		char szBuffer[MAX_PATH] = { 0 };
		char myDocumentsPath[MAX_PATH] = { 0 };

		OPENFILENAMEA file = { 0 };
		file.lStructSize = sizeof(file);
		file.lpstrFilter = "DataFile (*.dat)\0*.dat\0All Documents(*.*)\0*.*\0";
		file.lpstrTitle = "Select a MInDes Input File";
		file.lpstrFile = szBuffer;
		file.nMaxFile = MAX_PATH;
		file.Flags = OFN_PATHMUSTEXIST | OFN_FILEMUSTEXIST | OFN_EXPLORER | OFN_NOCHANGEDIR;

		if (GetOpenFileNameA(&file)) {
			file_path = file.lpstrFile;
		}
		else {
			DWORD dwError = CommDlgExtendedError();
			if (dwError != 0) {
				pf::printf_color_on_control("Error opening file dialog, error code: " + dwError, 31);
				printf("\n");
				exit(dwError);
			}
			else {
				pf::printf_color_on_control("File selection cancelled by user.\n", 31);
			}
		}
	}
	inline void selectAllFile(std::string& file_path) {
		char szBuffer[MAX_PATH] = { 0 };
		char myDocumentsPath[MAX_PATH] = { 0 };

		OPENFILENAMEA file = { 0 };
		file.lStructSize = sizeof(file);
		file.lpstrFilter = "All Documents(*.*)\0*.*\0";
		file.lpstrTitle = "Select a MInDes Input File";
		file.lpstrFile = szBuffer;
		file.nMaxFile = MAX_PATH;
		file.Flags = OFN_PATHMUSTEXIST | OFN_FILEMUSTEXIST | OFN_EXPLORER | OFN_NOCHANGEDIR;

		if (GetOpenFileNameA(&file)) {
			file_path = file.lpstrFile;
		}
		else {
			DWORD dwError = CommDlgExtendedError();
			if (dwError != 0) {
				pf::printf_color_on_control("Error opening file dialog, error code: " + dwError, 31);
				printf("\n");
				exit(dwError);
			}
			else {
				pf::printf_color_on_control("File selection cancelled by user.\n", 31);
			}
		}
	}
}
#endif // _WIN32
