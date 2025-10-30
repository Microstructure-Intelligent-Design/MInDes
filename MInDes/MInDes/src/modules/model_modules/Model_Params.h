#pragma once
#include "../base/Mesh_0.h"
#include "../base/RotationMatrix.h"
#include "../../MainIterator_Params.h"
#include "../input_modules/ioFiles_Params.h"
namespace pf {
	enum Dimension { One_Dimension, Two_Dimension, Three_Dimension };
	namespace mesh_parameters {
		// main mesh size
		inline size_t MESH_NX = 1;
		inline size_t MESH_NY = 1;
		inline size_t MESH_NZ = 1;
		inline Dimension dimention = Dimension::One_Dimension;
		// main mesh boundary condition
		inline BoundaryCondition x_down = BoundaryCondition::PERIODIC;
		inline BoundaryCondition y_down = BoundaryCondition::PERIODIC;
		inline BoundaryCondition z_down = BoundaryCondition::PERIODIC;
		inline BoundaryCondition x_up = BoundaryCondition::PERIODIC;
		inline BoundaryCondition y_up = BoundaryCondition::PERIODIC;
		inline BoundaryCondition z_up = BoundaryCondition::PERIODIC;
		// grid size
		inline REAL delt_r = 1.0;
	}
	namespace time_parameters {
		// real time
		inline REAL Real_Time = 0.0;
		// time interval for each simulation step
		inline REAL delt_t = 1.0;
	}
	namespace main_field {
		// main mesh for phase field
		inline bool is_phi_field_on = false;
		inline size_t phi_number = 0;
		inline REAL PHI_MAX_VARIATION = 0;
		inline Mesh_Boundry<std::vector<REAL>> phase_field;
		inline void init_phase_field() {
			phase_field.init(mesh_parameters::MESH_NX, mesh_parameters::MESH_NY, mesh_parameters::MESH_NZ, mesh_parameters::delt_r,
				mesh_parameters::x_down, mesh_parameters::x_up, mesh_parameters::y_down, mesh_parameters::y_up, mesh_parameters::z_down, mesh_parameters::z_up);
			for (long long z = 0; z < phase_field.Nz(); z++)
				for (long long y = 0; y < phase_field.Ny(); y++)
					for (long long x = 0; x < phase_field.Nx(); x++)
						phase_field(x, y, z).resize(phi_number, 0);
		}

		// main mesh for concentration field
		inline bool is_con_field_on = false;
		inline size_t con_number = 0;
		inline REAL CON_MAX_VARIATION = 0;
		inline Mesh_Boundry<std::vector<REAL>> concentration_field;
		inline void init_concentration_field() {
			concentration_field.init(mesh_parameters::MESH_NX, mesh_parameters::MESH_NY, mesh_parameters::MESH_NZ, mesh_parameters::delt_r,
				mesh_parameters::x_down, mesh_parameters::x_up, mesh_parameters::y_down, mesh_parameters::y_up, mesh_parameters::z_down, mesh_parameters::z_up);
			for (long long z = 0; z < concentration_field.Nz(); z++)
				for (long long y = 0; y < concentration_field.Ny(); y++)
					for (long long x = 0; x < concentration_field.Nx(); x++)
						concentration_field(x, y, z).resize(con_number, 0);
		}

		// main mesh for temperature field
		inline bool is_temp_field_on = false;
		inline REAL TEMP_MAX_VARIATION = 0;
		inline Mesh_Boundry<REAL> temperature_field;
		inline void init_temperature_field() {
			temperature_field.init(mesh_parameters::MESH_NX, mesh_parameters::MESH_NY, mesh_parameters::MESH_NZ, mesh_parameters::delt_r,
				mesh_parameters::x_down, mesh_parameters::x_up, mesh_parameters::y_down, mesh_parameters::y_up, mesh_parameters::z_down, mesh_parameters::z_up);
			for (long long z = 0; z < temperature_field.Nz(); z++)
				for (long long y = 0; y < temperature_field.Ny(); y++)
					for (long long x = 0; x < temperature_field.Nx(); x++)
						temperature_field(x, y, z) = 0;
		}
	}
}