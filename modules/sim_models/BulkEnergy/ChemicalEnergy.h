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
#include "../../Base.h"


namespace pf {
	enum dfdphiType { // phi , phase con , grand potential , temperature
		dfdphi_Const,
		dfdphi_DOUBEL_WELL,
		dfdphi_LQ_CHEN,
		dfdcon_1_CRACK_WELL,
		dfdcon_1_CRACK_OBSTACLE,
		dfdcon_N_CRACK_WELL,
		dfdcon_N_CRACK_OBSTACLE,
		dfdphi_HighOrder,
	};
	enum dfdconType { // phi , total con , temperature
		dfdcon_Const,
		dfdcon_HighOrder,
	};
	namespace preset_function {
		struct hlf_entry {
			hlf_entry() {
				index = 0;
			}
			int index;
			HLFunc _hlf;
		};
		class hlf_box
		{
		public:
			hlf_box() { index = 0; }
			~hlf_box() { _func_box.clear(); }
			vector<hlf_entry> _func_box;
			int index;
			typedef std::vector<hlf_entry>::iterator iterator;
			typedef std::vector<hlf_entry>::const_iterator citerator;
			iterator  begin() { return _func_box.begin(); };
			iterator  end() { return _func_box.end(); };
			hlf_entry& operator[](const int index) {
				for (auto i = _func_box.begin(); i < _func_box.end(); ++i) {
					if (i->index == index) return *i;
				}
				cout << "_func_box error, can't find the value index : " << index << endl;
				SYS_PROGRAM_STOP;
			}
			hlf_box& operator=(const hlf_box& n) {
				_func_box = n._func_box;
				index = n.index;
				return *this;
			}
			void add_hlf(int index, HLFunc _func) {
				for (auto i = _func_box.begin(); i < _func_box.end(); ++i)
					if (i->index == index) {
						i->_hlf = _func;
						return;
					}
				hlf_entry elem;
				elem.index = index;
				elem._hlf = _func;
				_func_box.push_back(elem);
			}
			void erase(int index) {
				for (auto i = _func_box.begin(); i < _func_box.end();) {
					if (i->index == index) {
						i = _func_box.erase(i);
					}
					else
						++i;
				}
			}
			void clear() {
				_func_box.clear();
			}
			int size() {
				return int(_func_box.size());
			}

		private:

		};
		class hlf_box2
		{
		public:
			hlf_box2() {  }
			~hlf_box2() { _func_box.clear(); }
			vector<hlf_box> _func_box;
			typedef std::vector<hlf_box>::iterator iterator;
			typedef std::vector<hlf_box>::const_iterator citerator;
			iterator  begin() { return _func_box.begin(); };
			iterator  end() { return _func_box.end(); };
			hlf_box& operator[](const int index) {
				for (auto i = _func_box.begin(); i < _func_box.end(); ++i) {
					if (i->index == index) return *i;
				}
				cout << "_func_box error, can't find the value index : " << index << endl;
				SYS_PROGRAM_STOP;
			}
			hlf_box2& operator=(const hlf_box2& n) {
				_func_box = n._func_box;
				return *this;
			}
			void add_hlf(int index1, int index2, HLFunc _func) {
				for (auto i = _func_box.begin(); i < _func_box.end(); ++i) {
					if (i->index == index1) {
						for (auto j = i->begin(); j < i->end(); ++j)
							if (j->index == index2) {
								j->_hlf = _func;
								return;
							}
						hlf_entry elem2;
						elem2.index = index2;
						elem2._hlf = _func;
						i->_func_box.push_back(elem2);
					}
				}
				hlf_entry elem2;
				elem2.index = index2;
				elem2._hlf = _func;
				hlf_box elem1;
				elem1.index = index1;
				elem1._func_box.push_back(elem2);
				_func_box.push_back(elem1);
			}
			void erase(int index) {
				for (auto i = _func_box.begin(); i < _func_box.end();) {
					if (i->index == index) {
						i = _func_box.erase(i);
					}
					else
						++i;
				}
			}
			void erase(int index1, int index2) {
				for (auto i = _func_box.begin(); i < _func_box.end();i++) {
					for (auto j = i->begin(); j < i->end();) {
						if (i->index == index1 && j->index == index2) {
							j = i->_func_box.erase(j);
						}
						else
							++j;
					}
				}
			}
			void clear() {
				_func_box.clear();
			}
			int size() {
				return int(_func_box.size());
			}
		};
		vector<double> comp2hlf(pf::ConNode x) {
			vector<double> _x;
			for (auto c = x.begin(); c < x.end(); c++)
				_x.push_back(c->value);
			return _x;
		}
	}
	namespace chemical_energy {
		pf::ConEquationDomain _domain = pf::ConEquationDomain::CEDomain_Standard;
		vector<int> smooth_phases;
		static preset_function::hlf_box  fchem;
		static preset_function::hlf_box2 grand_con;
		static preset_function::hlf_box2 dgrand_con_du;
		//-----------------------------------------------------------------------------------------------
		// get chemical energy density for each phase
		static double fchem_hlfuncs(pf::PhaseNode& node, pf::PhaseEntry& phase) {
			return fchem[phase.property]._hlf.cal_func(preset_function::comp2hlf(phase.x));
		}
		// get diffusion potential of component
		static double dfchem_dcon_hlfuncs(pf::PhaseEntry& phase, vector<double> con, int con_i) {
			return fchem[phase.property]._hlf.cal_dfunc_dxi(con, con_i);
		}
		//-----------------------------------------------------------------------------------------------
		// fchem for each phase
		static double dfchem_dphi_const(pf::PhaseNode& node, pf::PhaseEntry& phase){
			return 0.0;
		}
		static double dfchem_dphi_con_hlfuncs(pf::PhaseNode& node, pf::PhaseEntry& phase) {
			return fchem[phase.property]._hlf.cal_func(preset_function::comp2hlf(node.x));  //- ?
		}
		static double dfchem_dphi_phase_con_hlfuncs(pf::PhaseNode& node, pf::PhaseEntry& phase) {
			double tail = 0.0;
			for (auto con = phase.x.begin(); con < phase.x.end(); con++)
				tail += con->value * phase.potential[con->index].value;
			return fchem[phase.property]._hlf.cal_func(preset_function::comp2hlf(phase.x)) - tail;
		}
		// double well : V * phi_a * phi_a * (1.0 - phi_a) * (1.0 - phi_a) + W * sum_a{ sum_b!=a{ phi_a * phi_a * phi_b * phi_b } }
		static double DOUBLE_WELL_A = 0.0, DOUBLE_WELL_B = 0.0;
		static double dfchem_dphi_double_well(pf::PhaseNode& node, pf::PhaseEntry& phase) {
			double phi = phase.phi, sum_phi_j = 0.0;
			for (auto p = node.begin(); p < node.end(); p++)
				if (p->index != phase.index)
					sum_phi_j += p->phi * p->phi;
			return 2.0 * DOUBLE_WELL_A * phi * (1.0 - phi) * (1.0 - 2.0 * phi) + 2.0 * DOUBLE_WELL_B * phi * sum_phi_j;
		}
		// LQ. Chen : sum_a{ - A/2.0 * phi_a * phi_a + B/4.0 * phi_a * phi_a * phi_a * phi_a } + C * sum_a{ sum_b!=a{ phi_a * phi_a * phi_b * phi_b } }
		static double LQ_Chen_A = 0.0, LQ_Chen_B = 0.0, LQ_Chen_C = 0.0;
		static double dfchem_dphi_LQ_Chen(pf::PhaseNode& node, pf::PhaseEntry& phase) {
			double phi = phase.phi, sum_phi_j = 0.0;
			for (auto p = node.begin(); p < node.end(); p++)
				if (p->index != phase.index)
					sum_phi_j += p->phi * p->phi;
			return -LQ_Chen_A * phi + LQ_Chen_B * phi * phi * phi + 2.0 * LQ_Chen_C * phi * sum_phi_j;
		}
		//-----------------------------------------------------------------------------------------------
		// const : potential = con
		static void dfchem_dphase_con_const(pf::PhaseNode& node, pf::PhaseEntry& phase) {
			for (auto p = phase.potential.begin(); p < phase.potential.end(); p++)
				p->value += phase.x[p->index].value;
		}
		static void dfchem_dphase_con_hlfuncs(pf::PhaseNode& node, pf::PhaseEntry& phase) {
			vector<double> hlf_x = preset_function::comp2hlf(phase.x);
			for (auto p = phase.potential.begin(); p < phase.potential.end(); p++)
				p->value += fchem[phase.property]._hlf.cal_dfunc_dxi(hlf_x, p->index);
		}
		//-----------------------------------------------------------------------------------------------
		// fchem for total con
		static double dfchem_dcon_i_const(pf::PhaseNode& node, int con_i) {
			return node.x[con_i].value;
		}
		//-----------------------------------------------------------------------------------------------
		// dx_i^a / du_i
		static double dphase_con_i_du_i_const(pf::PhaseNode& node, pf::PhaseEntry& phase, int con_i) {
			if (_domain == ConEquationDomain::CEDomain_Standard) {
				for (auto index = smooth_phases.begin(); index < smooth_phases.end(); index++)
					if (phase.index == *index)
						return 1.0;
				return 0.0;
			}
			else if (_domain == ConEquationDomain::CEDomain_Reverse) {
				for (auto index = smooth_phases.begin(); index < smooth_phases.end(); index++)
					if (phase.index == *index)
						return 0.0;
				return 1.0;
			}
			return 0.0;
		}
		static double dphase_con_i_du_i_hlfuncs(pf::PhaseNode& node, pf::PhaseEntry& phase, int con_i) {
			vector<double> c_i; c_i.push_back(phase.x[con_i].value);
			return dgrand_con_du[phase.property][con_i]._hlf.cal_func(c_i);
		}
		//-----------------------------------------------------------------------------------------------
		// x_i^a = func( u_i )
		static void phase_con_const(pf::PhaseNode& node, pf::PhaseEntry& phase) {
			if (_domain == ConEquationDomain::CEDomain_Standard) {
				for (auto index = smooth_phases.begin(); index < smooth_phases.end(); index++) {
					if (phase.index == *index) {
						for (auto comp = phase.x.begin(); comp < phase.x.end(); comp++)
							comp->value = node.potential[comp->index].value;
						return;
					}
				}
				for (auto comp = phase.x.begin(); comp < phase.x.end(); comp++)
					comp->value = 0.0;
			}
			else if (_domain == ConEquationDomain::CEDomain_Reverse) {
				for (auto comp = phase.x.begin(); comp < phase.x.end(); comp++)
					comp->value = node.potential[comp->index].value;
				for (auto index = smooth_phases.begin(); index < smooth_phases.end(); index++)
					if (phase.index == *index)
						for (auto comp = phase.x.begin(); comp < phase.x.end(); comp++)
							comp->value = 0.0;
			}
		}
		static void phase_con_hlfuncs(pf::PhaseNode& node, pf::PhaseEntry& phase) {
			bool is_cal = false;
			if (_domain == ConEquationDomain::CEDomain_Standard) {
				is_cal = false;
				for (auto index = smooth_phases.begin(); index < smooth_phases.end(); index++)
					if (phase.index == *index)
						is_cal = true;
			}
			else if (_domain == ConEquationDomain::CEDomain_Reverse) {
				is_cal = true;
				for (auto index = smooth_phases.begin(); index < smooth_phases.end(); index++)
					if (phase.index == *index)
						is_cal = false;
			}
			if (is_cal) {
				for (auto comp = phase.x.begin(); comp < phase.x.end(); comp++) {
					vector<double> hlf_x; hlf_x.push_back(node.potential[comp->index].value);
					comp->value = grand_con[phase.property][comp->index]._hlf.cal_func(hlf_x);
				}
			}
			else {
				for (auto comp = phase.x.begin(); comp < phase.x.end(); comp++)
					comp->value = 0.0;
			}
		}

		static double (*fchem_density)(pf::PhaseNode& node, pf::PhaseEntry& phase);

		static double (*dfchem_dphi)(pf::PhaseNode& node, pf::PhaseEntry& phase);  // main function of this model

		static void (*dfchem_dphase_con)(pf::PhaseNode& node, pf::PhaseEntry& phase);  // main function of phase con

		static double (*dfchem_dcon_i)(pf::PhaseNode& node, int con_i);  // main function of total con & grand potential

		static double (*dphase_con_i_du_i)(pf::PhaseNode& node, pf::PhaseEntry& phase, int con_i);  // grand potential

		static void (*phase_con)(pf::PhaseNode& node, pf::PhaseEntry& phase);  // grand potential

		static void load_high_order_linear_functions(bool infile_debug) {
			if (infile_debug)
				InputFileReader::get_instance()->debug_writer->add_string_to_txt("# High Order Multivariate Linear Equation : f(xi) = sum_i{ sum_j[ Aij * xi^j ] }, i: component, j: order\n", InputFileReader::get_instance()->debug_file);
			if (Solvers::get_instance()->parameters.ConEType == ConEquationType::CEType_PhaseX) {
				dfchem_dphase_con = dfchem_dphase_con_hlfuncs;
				for (auto phase = Solvers::get_instance()->parameters.Phases.begin(); phase < Solvers::get_instance()->parameters.Phases.end(); phase++) {
					string energy_key = "ModelsManager.Phi.BulkEnergy." + phase->phi_name, energy_input = "[" + to_string(phase->x.size()) + "*N]";
					if (InputFileReader::get_instance()->read_string_value(energy_key, energy_input, infile_debug)) {
						HLFunc hlf; int comps = phase->x.size(), orders = 1;
						vector<vector<input_value>> energy_value = InputFileReader::get_instance()->trans_matrix_2d_const_const_to_input_value(InputValueType::IVType_DOUBLE, energy_key, energy_input, infile_debug);
						for (int cIndex = 0; cIndex < energy_value.size(); cIndex++)
							if (energy_value[cIndex].size() > orders)
								orders = int(energy_value[cIndex].size());
						hlf.init(comps, orders);
						for (int cIndex = 0; cIndex < energy_value.size(); cIndex++)
							for (int oIndex = 0; oIndex < energy_value[cIndex].size(); oIndex++)
								if (cIndex < comps && oIndex < orders)
									hlf.parameter_ij(cIndex, oIndex, energy_value[cIndex][oIndex].double_value);
						fchem.add_hlf(phase->phi_property, hlf);
					}
					else {
						HLFunc hlf;
						hlf.init(phase->x.size(), 1);
						fchem.add_hlf(phase->phi_property, hlf);
					}
				}
			}
			else if (Solvers::get_instance()->parameters.ConEType == ConEquationType::CEType_GrandP) {
				dfchem_dphase_con = dfchem_dphase_con_hlfuncs;
				dphase_con_i_du_i = dphase_con_i_du_i_hlfuncs;
				phase_con = phase_con_hlfuncs;
				for (auto phase = Solvers::get_instance()->parameters.Phases.begin(); phase < Solvers::get_instance()->parameters.Phases.end(); phase++) {
					string energy_key = "ModelsManager.Phi.BulkEnergy.Fchem." + phase->phi_name, energy_input = "[" + to_string(phase->x.size()) + "*N]";
					if (InputFileReader::get_instance()->read_string_value(energy_key, energy_input, infile_debug)) {
						HLFunc hlf; int comps = phase->x.size(), orders = 1;
						vector<vector<input_value>> energy_value = InputFileReader::get_instance()->trans_matrix_2d_const_const_to_input_value(InputValueType::IVType_DOUBLE, energy_key, energy_input, infile_debug);
						for (int cIndex = 0; cIndex < energy_value.size(); cIndex++)
							if (energy_value[cIndex].size() > orders)
								orders = int(energy_value[cIndex].size());
						hlf.init(comps, orders);
						for (int cIndex = 0; cIndex < energy_value.size(); cIndex++)
							for (int oIndex = 0; oIndex < energy_value[cIndex].size(); oIndex++)
								if (cIndex < comps && oIndex < orders)
									hlf.parameter_ij(cIndex, oIndex, energy_value[cIndex][oIndex].double_value);
						fchem.add_hlf(phase->phi_property, hlf);
					}
					else {
						HLFunc hlf;
						hlf.init(phase->x.size(), 1);
						fchem.add_hlf(phase->phi_property, hlf);
					}
					for (auto con = phase->x.begin(); con < phase->x.end(); con++) {
						energy_key = "ModelsManager.Phi.BulkEnergy.GrandCon." + phase->phi_name + "_" + con->name, energy_input = "()";
						if (InputFileReader::get_instance()->read_string_value(energy_key, energy_input, infile_debug)) {
							HLFunc hlf; int orders = 0;
							vector<input_value> con_value = InputFileReader::get_instance()->trans_matrix_1d_const_to_input_value(InputValueType::IVType_DOUBLE, energy_key, energy_input, infile_debug);
							if (con_value.size() > orders)
								orders = int(con_value.size());
							hlf.init(1, orders);
							for (int oIndex = 0; oIndex < con_value.size(); oIndex++)
								if (oIndex < orders)
									hlf.parameter_ij(0, oIndex, con_value[oIndex].double_value);
							grand_con.add_hlf(phase->phi_property, con->index, hlf);
						}
						else {
							HLFunc hlf;
							hlf.init(1, 1);
							grand_con.add_hlf(phase->phi_property, con->index, hlf);
						}
						energy_key = "ModelsManager.Phi.BulkEnergy.dPhiCon_dU." + phase->phi_name + "_" + con->name, energy_input = "()";
						if (InputFileReader::get_instance()->read_string_value(energy_key, energy_input, infile_debug)) {
							HLFunc hlf; int orders = 0;
							vector<input_value> con_value = InputFileReader::get_instance()->trans_matrix_1d_const_to_input_value(InputValueType::IVType_DOUBLE, energy_key, energy_input, infile_debug);
							if (con_value.size() > orders)
								orders = int(con_value.size());
							hlf.init(1, orders);
							for (int oIndex = 0; oIndex < con_value.size(); oIndex++)
								if (oIndex < orders)
									hlf.parameter_ij(0, oIndex, con_value[oIndex].double_value);
							dgrand_con_du.add_hlf(phase->phi_property, con->index, hlf);
						}
						else {
							HLFunc hlf;
							hlf.init(1, 1);
							dgrand_con_du.add_hlf(phase->phi_property, con->index, hlf);
						}
					}
				}
			}
		}

		static void load_chemical_energy_model(bool infile_debug) {
			if (Solvers::get_instance()->parameters.PhiEType == PhiEquationType::PEType_AC_Standard || Solvers::get_instance()->parameters.PhiEType == PhiEquationType::PEType_CH_Standard) {
				int model_type = 0;
				if (infile_debug)
				InputFileReader::get_instance()->debug_writer->add_string_to_txt("# ModelsManager.PhiCon.BulkEnergy.type : 1 - DoubleWell, 2 - LQ_Chen, 3 - H_Liang , 7 - HighOrder\n", InputFileReader::get_instance()->debug_file);
				InputFileReader::get_instance()->read_int_value("ModelsManager.PhiCon.BulkEnergy.type", model_type, infile_debug);
				switch (model_type)
				{
				case pf::dfdphi_DOUBEL_WELL:
					dfchem_dphi = dfchem_dphi_double_well;
					InputFileReader::get_instance()->read_double_value("ModelsManager.PhiCon.BulkEnergy.DoubleWell.A", DOUBLE_WELL_A, infile_debug);
					InputFileReader::get_instance()->read_double_value("ModelsManager.PhiCon.BulkEnergy.DoubleWell.B", DOUBLE_WELL_B, infile_debug);
					break;
				case pf::dfdphi_LQ_CHEN:
					dfchem_dphi = dfchem_dphi_LQ_Chen;
					InputFileReader::get_instance()->read_double_value("ModelsManager.PhiCon.BulkEnergy.LQ_Chen.A", LQ_Chen_A, infile_debug);
					InputFileReader::get_instance()->read_double_value("ModelsManager.PhiCon.BulkEnergy.LQ_Chen.B", LQ_Chen_B, infile_debug);
					InputFileReader::get_instance()->read_double_value("ModelsManager.PhiCon.BulkEnergy.LQ_Chen.C", LQ_Chen_C, infile_debug);
					break;
				case pf::dfdphi_HighOrder:
					InputFileReader::get_instance()->debug_writer->add_string_to_txt("# No models for HighOrder with standard AC and CH equations\n", InputFileReader::get_instance()->debug_file);
					std::exit(0);
					dfchem_dphi = dfchem_dphi_phase_con_hlfuncs;
					load_high_order_linear_functions(infile_debug);
					break;
				default:
					break;
				}
			}
			else if (Solvers::get_instance()->parameters.PhiEType == PhiEquationType::PEType_AC_Pairwise || Solvers::get_instance()->parameters.ConEType == ConEquationType::CEType_GrandP) {
				int model_type = 0;
				if (infile_debug)
				InputFileReader::get_instance()->debug_writer->add_string_to_txt("# ModelsManager.PhiCon.BulkEnergy.type : 7 - HighOrder\n", InputFileReader::get_instance()->debug_file);
				InputFileReader::get_instance()->read_int_value("ModelsManager.PhiCon.BulkEnergy.type", model_type, infile_debug);
				switch (model_type)
				{
				case pf::PEType_Const:
					dfchem_dphi = dfchem_dphi_const;
					break;
				case pf::dfdcon_1_CRACK_WELL:
					break;
				case pf::dfdcon_1_CRACK_OBSTACLE:
					break;
				case pf::dfdcon_N_CRACK_WELL:
					break;
				case pf::dfdcon_N_CRACK_OBSTACLE:
					break;
				case pf::dfdphi_HighOrder:
					fchem_density = fchem_hlfuncs;
					dfchem_dphi = dfchem_dphi_phase_con_hlfuncs;
					load_high_order_linear_functions(infile_debug);
					break;
				default:
					break;
				}
			}
		}

		static void init(FieldStorage_forPhaseNode& phaseMesh) {
			bool infile_debug = false;
			InputFileReader::get_instance()->read_bool_value("InputFile.debug", infile_debug, false);
			fchem_density = dfchem_dphi_const;
			dfchem_dphi = dfchem_dphi_const;
			dfchem_dphase_con = dfchem_dphase_con_const;
			dfchem_dcon_i = dfchem_dcon_i_const;
			dphase_con_i_du_i = dphase_con_i_du_i_const;
			phase_con = phase_con_const;
			smooth_phases = Solvers::get_instance()->C_Solver.phase_indexes;
			_domain = Solvers::get_instance()->parameters.ConEDomain;
			if (Solvers::get_instance()->parameters.PhiEType == PhiEquationType::PEType_AC_Standard || Solvers::get_instance()->parameters.PhiEType == PhiEquationType::PEType_CH_Standard)
				dfchem_dphi = dfchem_dphi_double_well;

			load_chemical_energy_model(infile_debug);

		}

		static void deinit(FieldStorage_forPhaseNode& phaseMesh) {

		}
	}
}