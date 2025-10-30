#pragma once
#include "../Module.h"
#include "inputfiles/QuickStartUp.h"
#include "inputfiles/infile_selector.h"
#include "ioFiles_Params.h"
#include "inputfiles/InputFileReader.h"
#include "inputfiles/UserStartUp.h"
namespace pf {
	static void init_input_modules() {
		stringstream out; 
		out << "============================== Parameters Definition Format =============================" << endl;
		out << "                             Parameter_Name = Parameter_Value" << endl;
		out << "======================================= M a c r o =======================================" << endl;
		out << "1.tube               " << ", $TUBE[1,10,2]$ = 1,3,5,7,9" << endl;
		out << "2.rand               " << ", $RAND_INT[1,10]$ = 1 - 10" << endl;
		out << "                     " << ", $RAND_REAL[1,10]$ = 1.000000 - 10.000000" << endl;
		out << "========================= Define Custom Variables and Functions =========================" << endl;
		out << "Define.Var  = VarName,0.1" << endl;
		out << "Define.Func = FuncName@{[VarName*pow(VarName,2)]}@" << endl;
		out << "default functions      : \"pow(val, ord)\", \"sqrt(val)\", \"abs(val)\", \"exp(val)\", \"ln(val)\"," << endl;
		out << "                         \"log(base_val, val)\", \"sin(val)\", \"cos(val)\", \"tan(val)\"," << endl;
		out << "                         \"asin(val)\", \"acos(val)\", \"atan(val)\", \"cos(val)\", " << endl;
		out << "                         \"tan(val)\", \"asin(val)\", \"acos(val)\", \"atan(val)\"" << endl;
		out << "=========================================================================================" << endl;
		add_string_to_file(out.str(), input_output_files_parameters::DebugFile_Path);
		InputFileReader::get_instance()->debug_infile_and_valid_words();
		InputFileReader::get_instance()->debug_custom_variavle_and_funcs();
		out.str("");
		out << "=========================================================================================" << endl;
		out << "================================= Parameters Definition =================================" << endl;
		out << "=========================================================================================" << endl;
		add_string_to_file(out.str(), input_output_files_parameters::DebugFile_Path);
		WriteLog("> MODULE INIT : Input File (.mindes) Ready !\n");
	}
}