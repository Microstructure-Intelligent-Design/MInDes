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
#include "../../solvers/base/sysTool.h"
#include "../../solvers/base/NormalNode.h"
#include "../../solvers/base/RotationMatrix.h"

#define NP_AnyPhase int(-1)
#define NP_x_DownLimit value
#define NP_x_UpLimit value2

namespace pf {
	using namespace std;

	enum NucleationProperty { DefiniteNucleation, ConditionalNucleation, UserDefined_Nucleated};
	enum Geometry{ Geo_None, Geo_Ellipsoid, Geo_Polyhedron};
	enum NucleationPosition { NP_Anywhere, NP_Bulk, NP_Interface };
	enum NucleationConcentration { NC_Default, NC_Defined};

	struct line_func_2D {
		// temp: a * x + b * y + c = 0
		line_func_2D() {
			a = 0.0;
			b = 0.0;
			c = 0.0;
		}
		line_func_2D(Point2D p1, Point2D p2) {
			init(p1, p2);
		}
		line_func_2D(double x1, double y1, double x2, double y2) {
			init(x1, y1, x2, y2);
		}
		void init(Point2D p1, Point2D p2) {
			a = p2.y - p1.y;
			b = p1.x - p2.x;
			c = p2.x * p1.y - p2.y * p1.x;
		}
		void init(double x1, double y1, double x2, double y2) {
			a = y2 - y1;
			b = x1 - x2;
			c = x2 * y1 - y2 * x1;
		}
		bool is_point_on_line(double x, double y) {
			if (a == 0 && b == 0)
				return false;
			double val = a * x + b * y + c;
			if (Is_Equality(val, 0.0))
				return true;
			else
				return false;
		}
		bool is_point_above_line(double x, double y) {
			if (a == 0 && b == 0)
				return false;
			if ((a * x + b * y + c) > 0.0)
				return true;
			else
				return false;
		}
		double a, b, c;
	};

	struct surf_func_3D {
		// temp: a * x + b * y + c * z + d = 0
		surf_func_3D() {
			a = 0.0;
			b = 0.0;
			c = 0.0;
			d = 0.0;
		}
		surf_func_3D(Point p1, Point p2, Point p3) {
			init(p1.x, p1.y, p1.z, p2.x, p2.y, p2.z, p3.x, p3.y, p3.z);
		}
		surf_func_3D(Point normal, Point point) {
			a = normal.x;
			b = normal.y;
			c = normal.z;
			d = -(a * point.x + b * point.y + c * point.z);
		}
		void init(double x1, double y1, double z1, double x2, double y2, double z2, double x3, double y3, double z3) {
			Vector3 t1(x2 - x1, y2 - y1, z2 - z1);
			Vector3 t2(x3 - x1, y3 - y1, z3 - z1);
			t1 = t1.cross(t2);
			a = t1[0];
			b = t1[1];
			c = t1[2];
			d = -(a * x1 + b * y1 + c * z1);
		}
		bool is_point_on_surf(double x, double y, double z) {
			double val = a * x + b * y + c * z + d;
			if (Is_Equality(val, 0.0))
				return true;
			else
				return false;
		}
		bool is_point_on_surf(Point p) {
			double val = a * p.x + b * p.y + c * p.z + d;
			if (Is_Equality(val, 0.0))
				return true;
			else
				return false;
		}
		bool is_point_above_surf(double x, double y, double z) {
			if ((a * x + b * y + c * z + d) > 0.0)
				return true;
			else
				return false;
		}
		bool is_point_above_surf(Point p) {
			if ((a * p.x + b * p.y + c * p.z + d) > 0.0)
				return true;
			else
				return false;
		}
		surf_func_3D& operator=(const surf_func_3D& n) {
			a = n.a;
			b = n.b;
			c = n.c;
			d = n.d;
			return *this;
		}
		double a, b, c, d;
	};
	struct Ellipsoid {
		//> (x - core.x) * (x - core.x) / radius_x / radius_x + (y - core.y) * (y - core.y) / radius_y / radius_y + (z - core.z) * (z - core.z) / radius_z / radius_z = 1.0
		Ellipsoid(){
			radius_x = 0.0;
			radius_y = 0.0;
			radius_z = 0.0;
			radian_x = 0.0;
			radian_y = 0.0;
			radian_z = 0.0;
			rotationGauge = RotationGauge::RG_XYX;
		}
		void set_core(double x, double y, double z) {
			core.x = x;
			core.y = y;
			core.z = z;
		}
		void set_radius(double x_radius, double y_radius, double z_radius) {
			radius_x = abs(x_radius) + SYS_EPSILON;
			radius_y = abs(y_radius) + SYS_EPSILON;
			radius_z = abs(z_radius) + SYS_EPSILON;
		}
		void set_rotation_radian_and_rotation_gauge(double _radian[3], RotationGauge _rotationGauge) {
			radian_x = -_radian[0];
			radian_y = -_radian[1];
			radian_z = -_radian[2];
			rotationGauge = _rotationGauge;
		}
		bool check_point_inside_ellipsoid(double x, double y, double z) {
			double test = (x - core.x) * (x - core.x) / radius_x / radius_x + (y - core.y) * (y - core.y) / radius_y / radius_y + (z - core.z) * (z - core.z) / radius_z / radius_z;
			if (test <= 1.0)
				return true;
			else
				return false;
		}
		bool check_point_inside_ellipsoid(Point p) {
			double test = (p.x - core.x) * (p.x - core.x) / radius_x / radius_x + (p.y - core.y) * (p.y - core.y) / radius_y / radius_y + (p.z - core.z) * (p.z - core.z) / radius_z / radius_z;
			if (test <= 1.0)
				return true;
			else
				return false;
		}
		Ellipsoid& operator=(const Ellipsoid& n) {
			this->core = n.core;
			this->radius_x = n.radius_x;
			this->radius_y = n.radius_y;
			this->radius_z = n.radius_z;
			this->radian_x = n.radian_x;
			this->radian_y = n.radian_y;
			this->radian_z = n.radian_z;
			rotationGauge = n.rotationGauge;
			return *this;
		}
		Point core;
		double radius_x;
		double radius_y;
		double radius_z;
		double radian_x;
		double radian_y;
		double radian_z;
		RotationGauge rotationGauge;
	};

	struct Polyhedron {
		Polyhedron(int inside_x = 0, int inside_y = 0, int inside_z = 0) {
			point_inside_polyhedron.set(inside_x, inside_y, inside_z);
			radian_x = 0.0;
			radian_y = 0.0;
			radian_z = 0.0;
			rotationGauge = RotationGauge::RG_XYX;
		}
		Polyhedron(Point inside_point) {
			point_inside_polyhedron.set(inside_point.x, inside_point.y, inside_point.z);
			radian_x = 0.0;
			radian_y = 0.0;
			radian_z = 0.0;
			rotationGauge = RotationGauge::RG_XYX;
		}
		void set_rotation_radian_and_rotation_gauge(double _radian[3], RotationGauge _rotationGauge) {
			radian_x = -_radian[0];
			radian_y = -_radian[1];
			radian_z = -_radian[2];
			rotationGauge = _rotationGauge;
		}
		void add_surf(Point p1, Point p2, Point p3) {
			surf_func_3D surf(p1, p2, p3);
			surfaces.push_back(surf);
		}
		void add_surf(Point norm, Point p) {
			surf_func_3D surf(norm, p);
			surfaces.push_back(surf);
		}
		void set_a_point_inside_polyhedron(int x, int y, int z) {
			point_inside_polyhedron.set(x, y, z);
		}
		void set_a_point_inside_polyhedron(Point p) {
			point_inside_polyhedron.set(p.x, p.y, p.z);
		}
		bool check_point_inside_polyhedron(int x, int y, int z) {
			bool is_inside = true;
			for (auto surf = surfaces.begin(); surf < surfaces.end(); surf++)
				if (surf->is_point_above_surf(x, y, z) != surf->is_point_above_surf(point_inside_polyhedron) && !surf->is_point_on_surf(x, y, z))
					is_inside = false;
			return is_inside;
		}
		bool check_point_inside_polyhedron(Point p) {
			bool is_inside = true;
			for (auto surf = surfaces.begin(); surf < surfaces.end(); surf++)
				if (surf->is_point_above_surf(p) != surf->is_point_above_surf(point_inside_polyhedron) && !surf->is_point_on_surf(p))
					is_inside = false;
			return is_inside;
		}
		Polyhedron& operator=(const Polyhedron& n) {
			this->surfaces = n.surfaces;
			this->point_inside_polyhedron = n.point_inside_polyhedron;
			this->radian_x = n.radian_x;
			this->radian_y = n.radian_y;
			this->radian_z = n.radian_z;
			rotationGauge = n.rotationGauge;
			return *this;
		}
		vector<surf_func_3D> surfaces;
		Point point_inside_polyhedron;
		double radian_x;
		double radian_y;
		double radian_z;
		RotationGauge rotationGauge;
	};

	struct GeometricRegion {
		Geometry geometryProperty;
		Ellipsoid ellipSolid;
		Polyhedron polyhedron;
		int generate_step; 
		int phaseProperty;  
		int phaseIndex;   
		XNode x;
		double temperature;
		double phi;
		bool isReverseRegion;
		bool isNormalized;

		double_box customValues;
		int_box customFlags;
		vec3_box customVec3s;
		vec6_box customVec6s;

		GeometricRegion(pf::Geometry _geometry_property = Geo_None, int _generate_step = 0, int _phase_index = 0, int _phase_property = 0, bool _isReverseRegion = false) {
			init(_geometry_property, _generate_step, _phase_index, _phase_property, _isReverseRegion);
			temperature = 0.0; 
			phi = 1.0;
			isNormalized = true;
		}
		void init(pf::Geometry _geometry_property = Geo_None, int _generate_step = 0, int _phase_index = 0, int _phase_property = 0, bool _isReverseRegion = false) {
			geometryProperty = _geometry_property;
			generate_step = _generate_step;
			phaseProperty = _phase_property;
			phaseIndex = _phase_index;
			isReverseRegion = _isReverseRegion;
		}
		GeometricRegion& operator=(const GeometricRegion& n) {
			this->generate_step = n.generate_step;
			this->geometryProperty = n.geometryProperty;
			this->temperature = n.temperature;
			this->phaseProperty = n.phaseProperty;
			this->phaseIndex = n.phaseIndex;
			this->ellipSolid = n.ellipSolid;
			this->polyhedron = n.polyhedron;
			this->x = n.x;
			this->phi = n.phi;
			this->customValues = n.customValues;
			this->customFlags = n.customFlags;
			this->customVec3s = n.customVec3s;
			this->customVec6s = n.customVec6s;
			this->isReverseRegion = n.isReverseRegion;
			this->isNormalized = n.isNormalized;
			return *this;
		}

	};

	class PointSet {
	public:
		PointSet(int _generate_step = 0, int _phase_index = 0, int _phase_property = 0, double _temperature = 0.0) {
			init(_generate_step, _phase_index, _phase_property, _temperature);
			is_normalized = true;
		}
		void init(int _generate_step = 0, int _phase_index = 0, int _phase_property = 0, double _temperature = 0.0) {
			generate_step = _generate_step;
			phaseProperty = _phase_property;
			phaseIndex = _phase_index;
			temperature = _temperature;
		}
		~PointSet() {
			points.clear();
		}
		vector<Point> points;
		vector<double> points_phi;
		int generate_step;
		int phaseProperty; 
		int phaseIndex;
		XNode x;
		double temperature;
		bool is_normalized;
		bool is_point_in_set(Point point, int x_limit, BoundaryCondition x_down_bc, BoundaryCondition x_up_bc, int y_limit = 1, BoundaryCondition y_down_bc = BoundaryCondition::ADIABATIC, BoundaryCondition y_up_bc = BoundaryCondition::ADIABATIC,
			int z_limit = 1, BoundaryCondition z_down_bc = BoundaryCondition::ADIABATIC, BoundaryCondition z_up_bc = BoundaryCondition::ADIABATIC) {
			int x = double_to_int(point.x), y = double_to_int(point.y), z = double_to_int(point.z);
			if (x < 0) {
				if (x_down_bc == BoundaryCondition::PERIODIC)
					x += x_limit;
				else
					x = 0;
			}
			else if (x >= x_limit) {
				if (x_up_bc == BoundaryCondition::PERIODIC)
					x -= x_limit;
				else
					x = x_limit - 1;
			}
			if (y < 0) {
				if (y_down_bc == BoundaryCondition::PERIODIC)
					y += y_limit;
				else
					y = 0;
			}
			else if (y >= y_limit) {
				if (y_up_bc == BoundaryCondition::PERIODIC)
					y -= y_limit;
				else
					y = y_limit - 1;
			}
			if (z < 0) {
				if (z_down_bc == BoundaryCondition::PERIODIC)
					z += z_limit;
				else
					z = 0;
			}
			else if (z >= z_limit) {
				if (z_up_bc == BoundaryCondition::PERIODIC)
					z -= z_limit;
				else
					z = z_limit - 1;
			}

			for (auto p = points.begin(); p < points.end(); p++)
				if (p->x == x && p->y == y && p->z == z)
					return true;
			return false;
		}
		void add_point(int _x, int _y, int _z, double _phi) {
			points.push_back(Point(_x, _y, _z));
			points_phi.push_back(_phi);
		}
		void add_point(Point _point, double _phi) {
			points.push_back(_point);
			points_phi.push_back(_phi);
		}
		PointSet& operator=(const PointSet& n) {
			this->generate_step = n.generate_step;
			this->points = n.points;
			this->temperature = n.temperature;
			this->phaseProperty = n.phaseProperty;
			this->phaseIndex = n.phaseIndex;
			this->x = n.x;
			this->points_phi = n.points_phi;
			this->is_normalized = n.is_normalized;
			return *this;
		}
	};

	struct Nucleus {
		int generate_step;    //range = [1, nstep)
		Point core;
		double phasefraction;   //range = (0, RADIUS_LIMIT)
		int phaseProperty;    //range = [0, phaseNum)
		int phaseIndex;     // 0 = basePhase, others[0, ~ ) = newPhase
		XNode x;
		Nucleus() {
			generate_step = 0;
			phasefraction = 1.0;
			phaseProperty = 0;
			phaseIndex = 0;
		}
		Nucleus& operator=(const Nucleus& n) {
			this->generate_step = n.generate_step;
			this->core = n.core;
			this->phasefraction = n.phasefraction;
			this->phaseProperty = n.phaseProperty;
			this->phaseIndex = n.phaseIndex;
			this->x = n.x;
			return *this;
		}
	};

	struct ConditionalPhase {
		int phaseProperty;    ///< property
		int phaseIndex;
		int nucleationPosition;
		vector<int> contactedPhases;   ///< property
		double nucleation_rate;   ///< 0.0 ~ 1.0
		int nucleation_con;
		XNode x;
		ConditionalPhase() {
			phaseProperty = 0;
			phaseIndex = 0;
			nucleationPosition = NucleationPosition::NP_Anywhere;
			nucleation_rate = 0.0;
			nucleation_con = NucleationConcentration::NC_Default;
		}
	};

	struct NucleationBox {
	public:
		std::vector<pf::Nucleus> nucleus_box;
		std::vector<pf::GeometricRegion> geometry_box;
		std::vector<pf::ConditionalPhase> condition_phi_box;
		std::vector<pf::PointSet> point_set_box;
		///< Property
		int nucleation_property;
		///< ConditionalNucleation
		int nucleation_step;
		NucleationBox() {
			nucleation_property = NucleationProperty::DefiniteNucleation;
			nucleation_step = 10000;
		};
		NucleationBox(NucleationProperty np) {
			nucleation_property = np;
			nucleation_step = 10000;
		}
		NucleationBox& operator=(const NucleationBox& n) {
			nucleus_box = n.nucleus_box;
			geometry_box = n.geometry_box;
			condition_phi_box = n.condition_phi_box;
			point_set_box = n.point_set_box;
			nucleation_property = n.nucleation_property;
			nucleation_step = n.nucleation_step;
			return *this;
		}
	};

	struct point_in_region_index {
		point_in_region_index(int re, int gi) {
			region = re;
			grain_index = gi;
		}
		int region;
		int grain_index;
	};
}
