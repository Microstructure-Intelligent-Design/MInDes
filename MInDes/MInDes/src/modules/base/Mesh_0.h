#pragma once
#include "MACRO_DEF.h"
#include <vector>
#include <iostream>
namespace pf {
	// - important
	inline size_t MESH_INDEX(size_t x, size_t y, size_t z, size_t Nx, size_t Ny) { return x + y * Nx + z * Nx * Ny; }
	inline size_t MESH_INDEX(long long x, long long y, long long z, size_t Nx, size_t Ny) { return x + y * Nx + z * Nx * Ny; }
	inline size_t MESH_INDEX(int x, int y, int z, size_t Nx, size_t Ny) { return x + y * Nx + z * Nx * Ny; }

	template <typename T>
	class Mesh {
	public:
		// - init
		Mesh(size_t NX, size_t NY, size_t NZ, REAL DeltR)
			: nx(NX), ny(NY), nz(NZ), dr(DeltR), data(NX* NY* NZ) {
			if (NX <= 0 || NY <= 0 || NZ <= 0) {
				std::cout << "ERROR ! size of mesh cannot be set as 0 !" << std::endl;
				SYS_PROGRAM_STOP;
			}
			if (DeltR <= 0) {
				std::cout << "ERROR ! DeltR is too small or smaller than 0 !" << std::endl;
				SYS_PROGRAM_STOP;
			}
		}
		Mesh() = default;
		void init(size_t NX, size_t NY, size_t NZ, REAL DeltR) {
			if (NX <= 0 || NY <= 0 || NZ <= 0) {
				std::cout << "ERROR ! size of mesh cannot be set as 0 !" << std::endl;
				SYS_PROGRAM_STOP;
			}
			if (DeltR <= 0) {
				std::cout << "ERROR ! DeltR is too small or smaller than 0 !" << std::endl;
				SYS_PROGRAM_STOP;
			}
			nx = NX;
			ny = NY;
			nz = NZ;
			dr = DeltR;
			data.resize(NX * NY * NZ);
		}
		~Mesh() {
			data.clear();
		}
		void clear() {
			nx = 0;
			ny = 0;
			nz = 0;
			dr = 1;
			data.clear();
			std::vector<T>().swap(data);
		}

		long long Nx() const { return nx; }
		long long Ny() const { return ny; }
		long long Nz() const { return nz; }
		REAL DeltR() const { return dr; }

		// - find element directly
		T& operator()(size_t x, size_t y, size_t z) {
			return data[MESH_INDEX(x, y, z, nx, ny)];
		}
		T& at(size_t x, size_t y, size_t z) {
			return data[MESH_INDEX(x, y, z, nx, ny)];
		}
		T& operator()(long long x, long long y, long long z) {
			return data[MESH_INDEX(x, y, z, nx, ny)];
		}
		T& at(long long x, long long y, long long z) {
			return data[MESH_INDEX(x, y, z, nx, ny)];
		}
		T& operator()(int x, int y, int z) {
			return data[MESH_INDEX(x, y, z, nx, ny)];
		}
		T& at(int x, int y, int z) {
			return data[MESH_INDEX(x, y, z, nx, ny)];
		}

	protected:
		size_t nx = 1;
		size_t ny = 1;
		size_t nz = 1;
		REAL dr = 1.0;
		std::vector<T> data;
	};

	enum class BoundaryCondition { FIXED, PERIODIC, ZEROFLUX};

	template <typename T>
	class Mesh_Boundry : public Mesh<T> {
	public:
		// - init
		Mesh_Boundry(size_t NX, size_t NY, size_t NZ, REAL DeltR, BoundaryCondition x_down, BoundaryCondition x_up
			, BoundaryCondition y_down, BoundaryCondition y_up, BoundaryCondition z_down, BoundaryCondition z_up) {
			if (NX <= 0 || NY <= 0 || NZ <= 0) {
				std::cout << "ERROR ! size of mesh cannot be set as 0 !" << std::endl;
				SYS_PROGRAM_STOP;
			}
			if (DeltR <= 0) {
				std::cout << "ERROR ! DeltR is too small or smaller than 0 !" << std::endl;
				SYS_PROGRAM_STOP;
			}
			Mesh<T>::nx = NX + 2;
			Mesh<T>::ny = NY + 2;
			Mesh<T>::nz = NZ + 2;
			Mesh<T>::dr = DeltR;
			bc_x_down = x_down;
			bc_x_up = x_up;
			bc_y_down = y_down;
			bc_y_up = y_up;
			bc_z_down = z_down;
			bc_z_up = z_up;
			Mesh<T>::data.resize((NX + 2) * (NY + 2) * (NZ + 2));
		}
		Mesh_Boundry() {
			Mesh<T>::nx = 2;
			Mesh<T>::ny = 2;
			Mesh<T>::nz = 2;
			Mesh<T>::dr = 1.0;
			bc_x_down = BoundaryCondition::PERIODIC;
			bc_x_up = BoundaryCondition::PERIODIC;
			bc_y_down = BoundaryCondition::PERIODIC;
			bc_y_up = BoundaryCondition::PERIODIC;
			bc_z_down = BoundaryCondition::PERIODIC;
			bc_z_up = BoundaryCondition::PERIODIC;
			Mesh<T>::data.resize(Mesh<T>::nx * Mesh<T>::ny * Mesh<T>::nz);
		};

		void init(size_t NX, size_t NY, size_t NZ, REAL DeltR, BoundaryCondition x_down, BoundaryCondition x_up
			, BoundaryCondition y_down, BoundaryCondition y_up, BoundaryCondition z_down, BoundaryCondition z_up) {
			if (NX <= 0 || NY <= 0 || NZ <= 0) {
				std::cout << "ERROR ! size of mesh cannot be set as 0 !" << std::endl;
				SYS_PROGRAM_STOP;
			}
			if (DeltR <= 0) {
				std::cout << "ERROR ! DeltR is too small or smaller than 0 !" << std::endl;
				SYS_PROGRAM_STOP;
			}
			Mesh<T>::nx = NX + 2;
			Mesh<T>::ny = NY + 2;
			Mesh<T>::nz = NZ + 2;
			Mesh<T>::dr = DeltR;
			bc_x_down = x_down;
			bc_x_up = x_up;
			bc_y_down = y_down;
			bc_y_up = y_up;
			bc_z_down = z_down;
			bc_z_up = z_up;
			Mesh<T>::data.resize((NX + 2) * (NY + 2) * (NZ + 2));
		}
		~Mesh_Boundry() {
			Mesh<T>::data.clear();
		}

		long long COMP_X_BGN() { return 1; };
		long long COMP_Y_BGN() { return 1; };
		long long COMP_Z_BGN() { return 1; };
		long long COMP_X_END() { return static_cast<long long>(Mesh<T>::nx) - 2; };
		long long COMP_Y_END() { return static_cast<long long>(Mesh<T>::ny) - 2; };
		long long COMP_Z_END() { return static_cast<long long>(Mesh<T>::nz) - 2; };

		T& at_boundary_x_down(long long y, long long z) {
			return Mesh<T>::data[MESH_INDEX(0LL, y, z, Mesh<T>::nx, Mesh<T>::ny)];
		}

		T& at_boundary_x_up(long long y, long long z) {
			return Mesh<T>::data[MESH_INDEX(static_cast<long long>(Mesh<T>::nx) - 1, y, z, Mesh<T>::nx, Mesh<T>::ny)];
		}

		T& at_boundary_y_down(long long x, long long z) {
			return Mesh<T>::data[MESH_INDEX(x, 0LL, z, Mesh<T>::nx, Mesh<T>::ny)];
		}

		T& at_boundary_y_up(long long x, long long z) {
			return Mesh<T>::data[MESH_INDEX(x, static_cast<long long>(Mesh<T>::ny) - 1, z, Mesh<T>::nx, Mesh<T>::ny)];
		}

		T& at_boundary_z_down(long long x, long long y) {
			return Mesh<T>::data[MESH_INDEX(x, y, 0LL, Mesh<T>::nx, Mesh<T>::ny)];
		}

		T& at_boundary_z_up(long long x, long long y) {
			return Mesh<T>::data[MESH_INDEX(x, y, static_cast<long long>(Mesh<T>::nz) - 1, Mesh<T>::nx, Mesh<T>::ny)];
		}

		BoundaryCondition BC_X_DOWN() { return bc_x_down; };
		BoundaryCondition BC_Y_DOWN() { return bc_y_down; };
		BoundaryCondition BC_Z_DOWN() { return bc_z_down; };
		BoundaryCondition BC_X_UP() { return bc_x_up; };
		BoundaryCondition BC_Y_UP() { return bc_y_up; };
		BoundaryCondition BC_Z_UP() { return bc_z_up; };

		void init_boundary_condition() {
			if (bc_x_down == BoundaryCondition::FIXED) {
#pragma omp parallel for
				for (int y = 0; y < Mesh<T>::ny; y++)
					for (int z = 0; z < Mesh<T>::nz; z++) {
						at_boundary_x_down(y, z) = Mesh<T>::at(1, y, z);
					}
			}
			if (bc_y_down == BoundaryCondition::FIXED) {
#pragma omp parallel for
				for (int x = 0; x < Mesh<T>::nx; x++)
					for (int z = 0; z < Mesh<T>::nz; z++) {
						at_boundary_y_down(x, z) = Mesh<T>::at(x, 1, z);
					}
			}
			if (bc_z_down == BoundaryCondition::FIXED) {
#pragma omp parallel for
				for (int x = 0; x < Mesh<T>::nx; x++)
					for (int y = 0; y < Mesh<T>::ny; y++) {
						at_boundary_z_down(x, y) = Mesh<T>::at(x, y, 1);
					}
			}
			if (bc_x_up == BoundaryCondition::FIXED) {
#pragma omp parallel for
				for (int y = 0; y < Mesh<T>::ny; y++)
					for (int z = 0; z < Mesh<T>::nz; z++) {
						at_boundary_x_up(y, z) = Mesh<T>::at(int(Mesh<T>::nx) - 2, y, z);
					}
			}
			if (bc_y_up == BoundaryCondition::FIXED) {
#pragma omp parallel for
				for (int x = 0; x < Mesh<T>::nx; x++)
					for (int z = 0; z < Mesh<T>::nz; z++) {
						at_boundary_y_up(x, z) = Mesh<T>::at(x, int(Mesh<T>::ny) - 2, z);
					}
			}
			if (bc_z_up == BoundaryCondition::FIXED) {
#pragma omp parallel for
				for (int x = 0; x < Mesh<T>::nx; x++)
					for (int y = 0; y < Mesh<T>::ny; y++) {
						at_boundary_z_up(x, y) = Mesh<T>::at(x, y, int(Mesh<T>::nz) - 2);
					}
			}
		}
		void do_boundary_condition() {
			if (bc_x_down == BoundaryCondition::PERIODIC) {
#pragma omp parallel for
				for (int y = 0; y < Mesh<T>::ny; y++)
					for (int z = 0; z < Mesh<T>::nz; z++) {
						at_boundary_x_down(y, z) = Mesh<T>::at(int(Mesh<T>::nx) - 2, y, z);
					}
			}
			else if (bc_x_down == BoundaryCondition::ZEROFLUX) {
#pragma omp parallel for
				for (int y = 0; y < Mesh<T>::ny; y++)
					for (int z = 0; z < Mesh<T>::nz; z++) {
						at_boundary_x_down(y, z) = Mesh<T>::at(1, y, z);
					}
			}
			if (bc_y_down == BoundaryCondition::PERIODIC) {
#pragma omp parallel for
				for (int x = 0; x < Mesh<T>::nx; x++)
					for (int z = 0; z < Mesh<T>::nz; z++) {
						at_boundary_y_down(x, z) = Mesh<T>::at(x, int(Mesh<T>::ny) - 2, z);
					}
			}
			else if (bc_y_down == BoundaryCondition::ZEROFLUX) {
#pragma omp parallel for
				for (int x = 0; x < Mesh<T>::nx; x++)
					for (int z = 0; z < Mesh<T>::nz; z++) {
						at_boundary_y_down(x, z) = Mesh<T>::at(x, 1, z);
					}
			}
			if (bc_z_down == BoundaryCondition::PERIODIC) {
#pragma omp parallel for
				for (int x = 0; x < Mesh<T>::nx; x++)
					for (int y = 0; y < Mesh<T>::ny; y++) {
						at_boundary_z_down(x, y) = Mesh<T>::at(x, y, int(Mesh<T>::nz) - 2);
					}
			}
			else if (bc_z_down == BoundaryCondition::ZEROFLUX) {
#pragma omp parallel for
				for (int x = 0; x < Mesh<T>::nx; x++)
					for (int y = 0; y < Mesh<T>::ny; y++) {
						at_boundary_z_down(x, y) = Mesh<T>::at(x, y, 1);
					}
			}
			if (bc_x_up == BoundaryCondition::PERIODIC) {
#pragma omp parallel for
				for (int y = 0; y < Mesh<T>::ny; y++)
					for (int z = 0; z < Mesh<T>::nz; z++) {
						at_boundary_x_up(y, z) = Mesh<T>::at(1, y, z);
					}
			}
			else if (bc_x_up == BoundaryCondition::ZEROFLUX) {
#pragma omp parallel for
				for (int y = 0; y < Mesh<T>::ny; y++)
					for (int z = 0; z < Mesh<T>::nz; z++) {
						at_boundary_x_up(y, z) = Mesh<T>::at(int(Mesh<T>::nx) - 2, y, z);
					}
			}
			if (bc_y_up == BoundaryCondition::PERIODIC) {
#pragma omp parallel for
				for (int x = 0; x < Mesh<T>::nx; x++)
					for (int z = 0; z < Mesh<T>::nz; z++) {
						at_boundary_y_up(x, z) = Mesh<T>::at(x, 1, z);
					}
			}
			else if (bc_y_up == BoundaryCondition::ZEROFLUX) {
#pragma omp parallel for
				for (int x = 0; x < Mesh<T>::nx; x++)
					for (int z = 0; z < Mesh<T>::nz; z++) {
						at_boundary_y_up(x, z) = Mesh<T>::at(x, int(Mesh<T>::ny) - 2, z);
					}
			}
			if (bc_z_up == BoundaryCondition::PERIODIC) {
#pragma omp parallel for
				for (int x = 0; x < Mesh<T>::nx; x++)
					for (int y = 0; y < Mesh<T>::ny; y++) {
						at_boundary_z_up(x, y) = Mesh<T>::at(x, y, 1);
					}
			}
			else if (bc_z_up == BoundaryCondition::ZEROFLUX) {
#pragma omp parallel for
				for (int x = 0; x < Mesh<T>::nx; x++)
					for (int y = 0; y < Mesh<T>::ny; y++) {
						at_boundary_z_up(x, y) = Mesh<T>::at(x, y, int(Mesh<T>::nz) - 2);
					}
			}
		}

	protected:
		BoundaryCondition bc_x_down = BoundaryCondition::PERIODIC;
		BoundaryCondition bc_x_up = BoundaryCondition::PERIODIC;
		BoundaryCondition bc_y_down = BoundaryCondition::PERIODIC;
		BoundaryCondition bc_y_up = BoundaryCondition::PERIODIC;
		BoundaryCondition bc_z_down = BoundaryCondition::PERIODIC;
		BoundaryCondition bc_z_up = BoundaryCondition::PERIODIC;
	};

}