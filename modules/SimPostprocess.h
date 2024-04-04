/*
This file is a part of the microstructure intelligent design software project.

Created:     Qi Huang 2023.04

Modified:    Qi Huang 2023.04;

Copyright (c) 2019-2023 Science center for phase diagram, phase transition, material intelligent design and manufacture, Central South University, China

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free
	Software Foundation, either version 3 of the License, or (at your option) any later version.

	This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
	FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

	You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/


#pragma once
#include "Base.h"
#include "sim_postprocess/ElectricField.h"
#include "sim_postprocess/FluidField.h"
#include "sim_postprocess/MechanicalField.h"
#include "sim_postprocess/WritePhiCT.h"
#include "sim_postprocess/CpuMemoryUsage.h"
#include "sim_postprocess/Statistics.h"
#include "sim_postprocess/Noise.h"
#include "sim_models/Source/Reaction/ElectrodeReaction.h"
#include "sim_postprocess/FieldOptimization.h"

namespace pf {
	namespace sim_postprocess {
		bool is_electric_field_on = false;
		bool is_fluid_field_on = false;
		bool is_mechanical_field_on = false;
		static void init(FieldStorage_forPhaseNode& phaseMesh) {
			bool infile_debug = false;
			InputFileReader::get_instance()->read_bool_value("InputFile.debug", infile_debug, false);

			if (infile_debug)
				InputFileReader::get_instance()->debug_writer->add_string_to_txt("# Postprocess.PhysicalFields.xxx = true/false (on/off) \n", InputFileReader::get_instance()->debug_file);

			InputFileReader::get_instance()->read_bool_value("Postprocess.PhysicalFields.mechanics", is_mechanical_field_on, infile_debug);

			InputFileReader::get_instance()->read_bool_value("Postprocess.PhysicalFields.fluid", is_fluid_field_on, infile_debug);

			InputFileReader::get_instance()->read_bool_value("Postprocess.PhysicalFields.electric", is_electric_field_on, infile_debug);

			if (is_electric_field_on)
				pf::electric_field::init(phaseMesh);
			if (is_fluid_field_on)
				pf::fluid_field::init(phaseMesh);
			if (is_mechanical_field_on)
				pf::mechanical_field::init(phaseMesh);
			noise::init(phaseMesh);
			field_optimization::init(phaseMesh);
			statistics::init(phaseMesh, is_mechanical_field_on);

			pf::WritePhiCT::init(phaseMesh);
			pf::CpuMemoryUsage::init(phaseMesh);
		}
		static void exec_pre(FieldStorage_forPhaseNode& phaseMesh) {
			if (is_electric_field_on) {

				//for (int x = 0; x < phaseMesh.limit_x; x++)
				//	for (int y = 0; y < phaseMesh.limit_y; y++)
				//		for (int z = 0; z < phaseMesh.limit_z; z++) {
				//			PhaseNode& node = phaseMesh(x, y, z);
				//			node.customValues[101] = 1.0;
				//		}
				pf::electric_field::exec_pre(phaseMesh);
			}
			if (is_fluid_field_on)
				pf::fluid_field::exec_pre(phaseMesh);
			if (is_mechanical_field_on)
				pf::mechanical_field::exec_pre(phaseMesh);
			noise::exec_pre(phaseMesh);
			field_optimization::exec_pre(phaseMesh);
			statistics::exec_pre(phaseMesh);

			pf::CpuMemoryUsage::exec_pre(phaseMesh);
		}
		static string exec_loop(FieldStorage_forPhaseNode& phaseMesh) {
			string report = "";
			if (is_electric_field_on)
				report += pf::electric_field::exec_loop(phaseMesh);
			if (is_fluid_field_on)
				report += pf::fluid_field::exec_loop(phaseMesh);
			if (is_mechanical_field_on)
				report += pf::mechanical_field::exec_loop(phaseMesh);
			report += noise::exec_loop(phaseMesh);
			report += field_optimization::exec_loop(phaseMesh);
			report += statistics::exec_loop(phaseMesh);

			return report;
		}
		static void deinit(FieldStorage_forPhaseNode& phaseMesh) {
			if (is_electric_field_on)
				pf::electric_field::deinit(phaseMesh);
			if (is_fluid_field_on)
				pf::fluid_field::deinit(phaseMesh);
			if (is_mechanical_field_on)
				pf::mechanical_field::deinit(phaseMesh);
			statistics::deinit(phaseMesh);
			noise::deinit(phaseMesh);
			field_optimization::deinit(phaseMesh);

			pf::CpuMemoryUsage::deinit(phaseMesh);
		}
		static void write_scalar(ofstream& fout, FieldStorage_forPhaseNode& phaseMesh) {
			pf::WritePhiCT::write_scalar(fout, phaseMesh);
			if (is_electric_field_on)
				pf::electric_field::write_scalar(fout, phaseMesh);
			if (is_fluid_field_on)
				pf::fluid_field::write_scalar(fout, phaseMesh);
			if (is_mechanical_field_on)
				pf::mechanical_field::write_scalar(fout, phaseMesh);
			electrode_reaction::write_scalar(fout, phaseMesh);
			noise::write_scalar(fout, phaseMesh);
			field_optimization::write_scalar(fout, phaseMesh);
		}
		static void write_vec3(ofstream& fout, FieldStorage_forPhaseNode& phaseMesh) {
			pf::WritePhiCT::write_vec3(fout, phaseMesh);
			if (is_electric_field_on)
				pf::electric_field::write_vec3(fout, phaseMesh);
			if (is_fluid_field_on)
				pf::fluid_field::write_vec3(fout, phaseMesh);
			if (is_mechanical_field_on)
				pf::mechanical_field::write_vec3(fout, phaseMesh);
			noise::write_vec3(fout, phaseMesh);
			field_optimization::write_vec3(fout, phaseMesh);
		}
		static void load_module() {
			Solvers::get_instance()->create_a_new_module(init, exec_pre, exec_loop, deinit, write_scalar, write_vec3);
		};
	}
}