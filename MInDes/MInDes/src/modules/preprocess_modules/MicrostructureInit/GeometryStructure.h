#pragma once
#include "../../base/VectorMatrix.h"
#include <cmath>
#include <iostream>
#include <random>
#include <sstream>
namespace pf {
	namespace geometry_structure {
		enum Geometry { Geo_None, Geo_Ellipsoid, Geo_Polyhedron, Geo_SegmentedCylinder, Geo_RectangularCuboid };
		enum class BoundaryCondition { FIXED, PERIODIC, ADIABATIC };
		const double GEO_EPSILON = 1e-6;
		class Vector3
		{
		public:
			Vector3() {
				storage[0] = 0.0;
				storage[1] = 0.0;
				storage[2] = 0.0;
			};
			Vector3(double i, double j, double k) {
				storage[0] = i;
				storage[1] = j;
				storage[2] = k;
			}

			Vector3(const Vector3& rhs)
			{
				for (int i = 0; i < 3; i++)
					storage[i] = rhs[i];
			};

			Vector3(const double vecinit[3])
			{
#ifdef DEBUG
				if (vecinit.size() != 3)
				{
					std::cout << "Error in dVector3::constructor()\n"
						<< "Initialization list size beyond storage range."
						<< "\nTerminating!!!" << std::endl;
					exit(13);
				}
#endif
				storage[0] = vecinit[0];
				storage[1] = vecinit[1];
				storage[2] = vecinit[2];
			}

			double& operator[](const int i)
			{
#ifdef DEBUG
				if (i > 2)
				{
					std::cout << "Error in dVector3::operator[]\n"
						<< "Access beyond storage range. i = "
						<< i << " > 2"
						<< "\nTerminating!!!" << std::endl;
					exit(13);
				}
#endif
				return storage[i];
			};
			double const& operator[](const int i) const
			{
#ifdef DEBUG
				if (i > 2)
				{
					std::cout << "Error in dVector3::operator[]\n"
						<< "Access beyond storage range. i = "
						<< i << " > 2"
						<< "\nTerminating!!!" << std::endl;
					exit(13);
				}
#endif
				return storage[i];
			};

			double getX(void) const
			{
				return storage[0];
			};

			void setX(const double newX)
			{
				storage[0] = newX;
			};

			double getY(void) const
			{
				return storage[1];
			};

			void setY(const double newY)
			{
				storage[1] = newY;
			};

			double getZ(void) const
			{
				return storage[2];
			};

			void setZ(const double newX)
			{
				storage[2] = newX;
			};

			void set_to_zero(void)
			{
				storage[0] = 0.0;
				storage[1] = 0.0;
				storage[2] = 0.0;
			};
			void set_to_unitX(void)
			{
				storage[0] = 1.0;
				storage[1] = 0.0;
				storage[2] = 0.0;
			};
			void set_to_unitY(void)
			{
				storage[0] = 0.0;
				storage[1] = 1.0;
				storage[2] = 0.0;
			};
			void set_to_unitZ(void)
			{
				storage[0] = 0.0;
				storage[1] = 0.0;
				storage[2] = 1.0;
			};
			bool operator==(const Vector3& rhs) const
			{
				for (int i = 0; i < 3; i++) {
					if (storage[i] != rhs[i])
						return false;
				}
				return true;
			};
			bool operator!=(const Vector3& rhs) const
			{
				for (int i = 0; i < 3; i++) {
					if (storage[i] != rhs[i])
						return true;
				}
				return false;
			};
			Vector3 operator*(const double m) const
			{
				Vector3 tmp;
				tmp[0] = storage[0] * m;
				tmp[1] = storage[1] * m;
				tmp[2] = storage[2] * m;
				return tmp;
			};
			Vector3 operator/(const double m) const
			{
				Vector3 tmp;
				tmp[0] = storage[0] / m;
				tmp[1] = storage[1] / m;
				tmp[2] = storage[2] / m;
				return tmp;
			};
			double operator*(const Vector3& rhs) const
			{
				return storage[0] * rhs[0] + storage[1] * rhs[1] + storage[2] * rhs[2];
			};
			double abs() const
			{
				return std::sqrt(storage[0] * storage[0] +
					storage[1] * storage[1] +
					storage[2] * storage[2]);
			};
			Vector3 cross(const Vector3& rhs) const
			{
				Vector3 tmp;
				tmp[0] = storage[1] * rhs[2] - storage[2] * rhs[1];
				tmp[1] = storage[2] * rhs[0] - storage[0] * rhs[2];
				tmp[2] = storage[0] * rhs[1] - storage[1] * rhs[0];
				return tmp;
			};
			Vector3 operator+(const Vector3& rhs) const
			{
				Vector3 tmp;
				tmp[0] = storage[0] + rhs[0];
				tmp[1] = storage[1] + rhs[1];
				tmp[2] = storage[2] + rhs[2];
				return tmp;
			};
			Vector3 operator-(const Vector3& rhs) const
			{
				Vector3 tmp;
				tmp[0] = storage[0] - rhs[0];
				tmp[1] = storage[1] - rhs[1];
				tmp[2] = storage[2] - rhs[2];
				return tmp;
			};
			Vector3& operator*=(const double m)
			{
				storage[0] *= m;
				storage[1] *= m;
				storage[2] *= m;
				return *this;
			};
			Vector3& operator/=(const double m)
			{
				storage[0] /= m;
				storage[1] /= m;
				storage[2] /= m;
				return *this;
			};
			Vector3& operator-=(const Vector3& rhs)
			{
				storage[0] = storage[0] - rhs[0];
				storage[1] = storage[1] - rhs[1];
				storage[2] = storage[2] - rhs[2];
				return *this;
			};
			Vector3& operator+=(const Vector3& rhs)
			{
				storage[0] = storage[0] + rhs[0];
				storage[1] = storage[1] + rhs[1];
				storage[2] = storage[2] + rhs[2];
				return *this;
			};
			Vector3& operator=(const Vector3& rhs)
			{
				storage[0] = rhs[0];
				storage[1] = rhs[1];
				storage[2] = rhs[2];
				return *this;
			};
			Vector3& operator=(const double rhs[3])
			{
				storage[0] = rhs[0];
				storage[1] = rhs[1];
				storage[2] = rhs[2];
				return *this;
			};
			double length(void) const
			{
				return std::sqrt(storage[0] * storage[0] +
					storage[1] * storage[1] +
					storage[2] * storage[2]);
			};
			Vector3& normalize(void)
			{
				double norm = std::sqrt(storage[0] * storage[0] +
					storage[1] * storage[1] +
					storage[2] * storage[2]);
				if (norm != 0.0)
				{
					storage[0] /= norm;
					storage[1] /= norm;
					storage[2] /= norm;
				}
				return *this;
			};
			Vector3 normalized(void) const
			{
				double norm = std::sqrt(storage[0] * storage[0] +
					storage[1] * storage[1] +
					storage[2] * storage[2]);
				Vector3 tmp;
				if (norm != 0.0)
				{
					tmp[0] = storage[0] / norm;
					tmp[1] = storage[1] / norm;
					tmp[2] = storage[2] / norm;
				}
				else
				{
					tmp.set_to_zero();
				}
				return tmp;
			};
			Vector3 Xreflected(void) const
			{
				Vector3 Out;
				for (int i = 0; i < 3; ++i)
				{
					Out[i] = storage[i];
				}
				Out[0] *= -1.0;
				return Out;
			};
			Vector3 Yreflected(void) const
			{
				Vector3 Out;
				for (int i = 0; i < 3; ++i)
				{
					Out[i] = storage[i];
				}
				Out[1] *= -1.0;
				return Out;
			};
			Vector3 Zreflected(void) const
			{
				Vector3 Out;
				for (int i = 0; i < 3; ++i)
				{
					Out[i] = storage[i];
				}
				Out[2] *= -1.0;
				return Out;
			};
			std::string print(void) const
			{
				std::stringstream out;

				out << "(" << storage[0] << ", "
					<< storage[1] << ", "
					<< storage[2] << ")";
				return out.str();
			};
			double* data(void)
			{
				return storage;
			};
			const double* const_data(void) const
			{
				return storage;
			};
		protected:
		private:
			double storage[3];
		};
		struct point_in_region_index {
			point_in_region_index(size_t re, size_t gi) {
				region = re;
				grain_index = gi;
			}
			size_t region;
			size_t grain_index;
		};
		inline bool isTwodoubleEquality(double a, double b) {
			if ((a - b) < (0.000001f) && (a - b) > -(0.000001f))
				return true;
			else
				return false;
		}
		inline int double_to_int(double a) {
			if ((a - int(a)) > 0.5f)
				return int(a) + 1;
			else
				return int(a);
		}
		struct Point {
			double x;
			double y;
			double z;
			Point() {
				x = 0;
				y = 0;
				z = 0;
			}
			Point(double _x, double _y, double _z) {
				x = _x;
				y = _y;
				z = _z;
			}
			Point(int _x, int _y, int _z) {
				x = double(_x);
				y = double(_y);
				z = double(_z);
			}
			Point(size_t _x, size_t _y, size_t _z) {
				x = double(_x);
				y = double(_y);
				z = double(_z);
			}
			void set(double _x, double _y, double _z) {
				x = _x;
				y = _y;
				z = _z;
			}
			double to_length(double _x, double _y, double _z) {
				return std::sqrt((_x - x) * (_x - x) + (_y - y) * (_y - y) + (_z - z) * (_z - z));
			}
			double to_length(Point p) {
				return std::sqrt((p.x - x) * (p.x - x) + (p.y - y) * (p.y - y) + (p.z - z) * (p.z - z));
			}
			double to_cosine(Point direction1, Point direction2) {
				Vector3 vec1(direction1.x - x, direction1.y - y, direction1.z - z), vec2(direction2.x - x, direction2.y - y, direction2.z - z);
				return vec1 * vec2 / vec1.length() / vec2.length();
			}
			void do_boundary(BoundaryCondition x_up_bc, BoundaryCondition y_up_bc, BoundaryCondition z_up_bc, BoundaryCondition x_down_bc, BoundaryCondition y_down_bc, BoundaryCondition z_down_bc, double x_limit, double y_limit, double z_limit) {
				if (x_down_bc == BoundaryCondition::ADIABATIC && x < 0)
					x = 0;
				else if (x_up_bc == BoundaryCondition::ADIABATIC && x >= x_limit)
					x = x_limit - 1;
				else if (x_down_bc == BoundaryCondition::FIXED && x < 0)
					x = 0;
				else if (x_up_bc == BoundaryCondition::FIXED && x >= x_limit)
					x = x_limit - 1;
				else if (x_down_bc == BoundaryCondition::PERIODIC && x < 0)
					x += x_limit;
				else if (x_up_bc == BoundaryCondition::PERIODIC && x >= x_limit)
					x -= x_limit;
				else if (y_down_bc == BoundaryCondition::ADIABATIC && y < 0)
					y = 0;
				else if (y_up_bc == BoundaryCondition::ADIABATIC && y >= y_limit)
					y = y_limit - 1;
				else if (y_down_bc == BoundaryCondition::FIXED && y < 0)
					y = 0;
				else if (y_up_bc == BoundaryCondition::FIXED && y >= y_limit)
					y = y_limit - 1;
				else if (y_down_bc == BoundaryCondition::PERIODIC && y < 0)
					y += y_limit;
				else if (y_up_bc == BoundaryCondition::PERIODIC && y >= y_limit)
					y -= y_limit;
				else if (z_down_bc == BoundaryCondition::ADIABATIC && z < 0)
					z = 0;
				else if (z_up_bc == BoundaryCondition::ADIABATIC && z >= z_limit)
					z = z_limit - 1;
				else if (z_down_bc == BoundaryCondition::FIXED && z < 0)
					z = 0;
				else if (z_up_bc == BoundaryCondition::FIXED && z >= z_limit)
					z = z_limit - 1;
				else if (z_down_bc == BoundaryCondition::PERIODIC && z < 0)
					z += z_limit;
				else if (z_up_bc == BoundaryCondition::PERIODIC && z >= z_limit)
					z -= z_limit;
				else
					return;
				return do_boundary(x_up_bc, y_up_bc, z_up_bc, x_down_bc, y_down_bc, z_down_bc, x_limit, y_limit, z_limit);
			}
			Point& operator=(const Point& n) {
				this->x = n.x;
				this->y = n.y;
				this->z = n.z;
				return *this;
			}
			Point operator+(const Point& n) {
				Point re;
				re.x = this->x + n.x;
				re.y = this->y + n.y;
				re.z = this->z + n.z;
				return re;
			}
			Point operator-(const Point& n) {
				Point re;
				re.x = this->x - n.x;
				re.y = this->y - n.y;
				re.z = this->z - n.z;
				return re;
			}
			Point operator*(const double n) {
				Point re;
				re.x = this->x * n;
				re.y = this->y * n;
				re.z = this->z * n;
				return re;
			}
			Point operator/(const double n) {
				Point re;
				re.x = this->x / n;
				re.y = this->y / n;
				re.z = this->z / n;
				return re;
			}
			Point& operator+=(const Point& n) {
				this->x += n.x;
				this->y += n.y;
				this->z += n.z;
				return *this;
			}
			Point& operator-=(const Point& n) {
				this->x -= n.x;
				this->y -= n.y;
				this->z -= n.z;
				return *this;
			}
			Point& operator*=(const double n) {
				this->x *= n;
				this->y *= n;
				this->z *= n;
				return *this;
			}
			Point& operator/=(const double n) {
				this->x /= n;
				this->y /= n;
				this->z /= n;
				return *this;
			}
			Point operator+(const Vector3& n) {
				Point re;
				re.x = this->x + n[0];
				re.y = this->y + n[1];
				re.z = this->z + n[2];
				return re;
			}
			Point operator-(const Vector3& n) {
				Point re;
				re.x = this->x - n[0];
				re.y = this->y - n[1];
				re.z = this->z - n[2];
				return re;
			}
			Point& operator+=(const Vector3& n) {
				this->x += n[0];
				this->y += n[1];
				this->z += n[2];
				return *this;
			}
			Point& operator-=(const Vector3& n) {
				this->x -= n[0];
				this->y -= n[1];
				this->z -= n[2];
				return *this;
			}
		};
		inline Vector3 get_vector(Point Tail, Point Head) {
			return Vector3(Head.x - Tail.x, Head.y - Tail.y, Head.z - Tail.z);
		}
		struct surf_func_3D {
			// temp: a * x + b * y + c * z + d = 0
			surf_func_3D() {
				a = 0;
				b = 0;
				c = 0;
				d = 0;
			}
			surf_func_3D(Point p1, Point p2, Point p3) {
				init(p1.x, p1.y, p1.z, p2.x, p2.y, p2.z, p3.x, p3.y, p3.z);
			}
			surf_func_3D(Vector3 normal, Point point) {
				a = normal[0];
				b = normal[1];
				c = normal[2];
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
			void move(double x, double y, double z) {
				d += -(a * x + b * y + c * z);
			}
			double distance(double _x, double _y, double _z) {
				return std::abs(a * _x + b * _y + c * _z + d) / std::sqrt(a * a + b * b + c * c);
			}
			bool is_point_on_surf(double x, double y, double z) {
				double val = a * x + b * y + c * z + d;
				if (isTwodoubleEquality(val, 0.0f))
					return true;
				else
					return false;
			}
			bool is_point_on_surf(Point p) {
				double val = a * p.x + b * p.y + c * p.z + d;
				if (isTwodoubleEquality(val, 0.0f))
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
			Ellipsoid() {
				radius_x = 0;
				radius_y = 0;
				radius_z = 0;
				radian_x = 0;
				radian_y = 0;
				radian_z = 0;
				rotationGauge = 0;
			}
			void set_core(double x, double y, double z) {
				core.x = x;
				core.y = y;
				core.z = z;
			}
			void set_radius(double x_radius, double y_radius, double z_radius) {
				radius_x = abs(x_radius) + GEO_EPSILON;
				radius_y = abs(y_radius) + GEO_EPSILON;
				radius_z = abs(z_radius) + GEO_EPSILON;
			}
			void move(double x, double y, double z) {
				core.x += x;
				core.y += y;
				core.z += z;
			}
			void add_rotation_radian(double radian1, double radian2, double radian3) {
				radian_x += -radian1;
				radian_y += -radian2;
				radian_z += -radian3;
			}
			void set_rotation_radian_and_rotation_gauge(double _radian[3], int _rotationGauge) {
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
			int rotationGauge;
		};
		struct SegmentedCylinder {
			SegmentedCylinder() {
				radius = 0;
				geometric_center = Point(0, 0, 0);
				radian_x = 0;
				radian_y = 0;
				radian_z = 0;
				rotationGauge = 0;
			}
			void add_point(Point p) {
				central_axis_line.push_back(p);
				geometric_center = Point(0, 0, 0);
				for (int index = 0; index < central_axis_line.size(); index++)
					geometric_center += central_axis_line[index];
				geometric_center /= double(central_axis_line.size());
			}
			void set_radius(double _radius) {
				radius = abs(_radius) + (0.000001f);
			}
			void move(double x, double y, double z) {
				geometric_center += Point(x, y, z);
				for (int index = 0; index < central_axis_line.size(); index++)
					central_axis_line[index] += Point(x, y, z);
			}
			void add_rotation_radian(double radian1, double radian2, double radian3) {
				radian_x += -radian1;
				radian_y += -radian2;
				radian_z += -radian3;
			}
			void set_rotation_radian_and_rotation_gauge(double _radian[3], int _rotationGauge) {
				radian_x = -_radian[0];
				radian_y = -_radian[1];
				radian_z = -_radian[2];
				rotationGauge = _rotationGauge;
			}
			bool check_point_inside_segmented_cylinder(double x, double y, double z) {
				Point p(x, y, z);
				for (int index = 1; index < central_axis_line.size(); index++) {
					double cosin1 = central_axis_line[index - 1].to_cosine(p, central_axis_line[index]),
						cosin2 = central_axis_line[index].to_cosine(p, central_axis_line[index - 1]);
					if (cosin1 > 0.0 && cosin2 > 0.0) {
						double length = central_axis_line[index].to_length(p), d2 = cosin2 * length,
							p2l_length_2 = length * length - d2 * d2;
						if (p2l_length_2 <= radius * radius)
							return true;
					}
					if (index < central_axis_line.size() - 1) {
						double length2 = (p.x - central_axis_line[index].x) * (p.x - central_axis_line[index].x)
							+ (p.y - central_axis_line[index].y) * (p.y - central_axis_line[index].y)
							+ (p.z - central_axis_line[index].z) * (p.z - central_axis_line[index].z);
						if (length2 <= radius * radius)
							return true;
					}
					if (index == central_axis_line.size() - 1) {
						if (isTwodoubleEquality(central_axis_line[0].x, central_axis_line[index].x)
							&& isTwodoubleEquality(central_axis_line[0].y, central_axis_line[index].y)
							&& isTwodoubleEquality(central_axis_line[0].z, central_axis_line[index].z)) {
							double length2 = (p.x - central_axis_line[index].x) * (p.x - central_axis_line[index].x)
								+ (p.y - central_axis_line[index].y) * (p.y - central_axis_line[index].y)
								+ (p.z - central_axis_line[index].z) * (p.z - central_axis_line[index].z);
							if (length2 <= radius * radius)
								return true;
						}
					}
				}
				return false;
			}
			bool check_point_inside_segmented_cylinder(Point p) {
				for (size_t index = 1; index < central_axis_line.size(); index++) {
					double cosin1 = central_axis_line[index - 1].to_cosine(p, central_axis_line[index]),
						cosin2 = central_axis_line[index].to_cosine(p, central_axis_line[index - 1]);
					if (cosin1 > 0.0 && cosin2 > 0.0) {
						double length = central_axis_line[index].to_length(p), d2 = cosin2 * length,
							p2l_length_2 = length * length - d2 * d2;
						if (p2l_length_2 < radius * radius)
							return true;
					}
					if (index < central_axis_line.size() - 1) {
						double length2 = (p.x - central_axis_line[index].x) * (p.x - central_axis_line[index].x)
							+ (p.y - central_axis_line[index].y) * (p.y - central_axis_line[index].y)
							+ (p.z - central_axis_line[index].z) * (p.z - central_axis_line[index].z);
						if (length2 <= radius * radius)
							return true;
					}
					if (index == central_axis_line.size() - 1) {
						if (isTwodoubleEquality(central_axis_line[0].x, central_axis_line[index].x)
							&& isTwodoubleEquality(central_axis_line[0].y, central_axis_line[index].y)
							&& isTwodoubleEquality(central_axis_line[0].z, central_axis_line[index].z)) {
							double length2 = (p.x - central_axis_line[index].x) * (p.x - central_axis_line[index].x)
								+ (p.y - central_axis_line[index].y) * (p.y - central_axis_line[index].y)
								+ (p.z - central_axis_line[index].z) * (p.z - central_axis_line[index].z);
							if (length2 <= radius * radius)
								return true;
						}
					}
				}
				return false;
			}
			SegmentedCylinder& operator=(const SegmentedCylinder& n) {
				this->central_axis_line = n.central_axis_line;
				this->radius = n.radius;
				this->geometric_center = n.geometric_center;
				this->radian_x = n.radian_x;
				this->radian_y = n.radian_y;
				this->radian_z = n.radian_z;
				rotationGauge = n.rotationGauge;
				return *this;
			}
			std::vector<Point> central_axis_line;
			double radius;
			Point geometric_center;
			double radian_x;
			double radian_y;
			double radian_z;
			int rotationGauge;
		};
		struct Polyhedron {
			Polyhedron(double inside_x = 0, double inside_y = 0, double inside_z = 0) {
				point_inside_polyhedron.set(inside_x, inside_y, inside_z);
				radian_x = 0;
				radian_y = 0;
				radian_z = 0;
				rotationGauge = 0;
			}
			Polyhedron(Point inside_point) {
				point_inside_polyhedron.set(inside_point.x, inside_point.y, inside_point.z);
				radian_x = 0;
				radian_y = 0;
				radian_z = 0;
				rotationGauge = 0;
			}
			void move(double x, double y, double z) {
				point_inside_polyhedron += Point(x, y, z);
				for (size_t index = 0; index < surfaces.size(); index++)
					surfaces[index].move(x, y, z);
			}
			void add_rotation_radian(double radian1, double radian2, double radian3) {
				radian_x += -radian1;
				radian_y += -radian2;
				radian_z += -radian3;
			}
			void set_rotation_radian_and_rotation_gauge(double _radian[3], int _rotationGauge) {
				radian_x = -_radian[0];
				radian_y = -_radian[1];
				radian_z = -_radian[2];
				rotationGauge = _rotationGauge;
			}
			void add_surf(Point p1, Point p2, Point p3) {
				surf_func_3D surf(p1, p2, p3);
				surfaces.push_back(surf);
			}
			void add_surf(Vector3 norm, Point p) {
				surf_func_3D surf(norm, p);
				surfaces.push_back(surf);
			}
			void set_a_point_inside_polyhedron(double x, double y, double z) {
				point_inside_polyhedron.set(x, y, z);
			}
			void set_a_point_inside_polyhedron(Point p) {
				point_inside_polyhedron.set(p.x, p.y, p.z);
			}
			bool check_point_inside_polyhedron(double x, double y, double z) {
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
			std::vector<surf_func_3D> surfaces;
			Point point_inside_polyhedron;
			double radian_x;
			double radian_y;
			double radian_z;
			int rotationGauge;
		};
		struct GeometricRegion {
			Geometry geometryProperty;
			Ellipsoid ellipSolid;
			Polyhedron polyhedron;
			SegmentedCylinder cylinder;
			size_t generate_step;
			size_t phaseIndex;
			Matrix1D<double> con;
			double temperature;
			double phi;
			bool is_phi = false;
			bool is_con = false;
			bool is_temp = false;
			bool isReverseRegion;
			bool isNormalized;

			GeometricRegion(Geometry _geometry_property = Geometry::Geo_None, int _generate_step = 0, int _phase_index = 0, bool _isReverseRegion = false) {
				init(_geometry_property, _generate_step, _phase_index, _isReverseRegion);
				temperature = 0;
				phi = 0;
				isNormalized = false;
			}
			void init(Geometry _geometry_property = Geometry::Geo_None, int _generate_step = 0, int _phase_index = 0, bool _isReverseRegion = false) {
				geometryProperty = _geometry_property;
				generate_step = _generate_step;
				phaseIndex = _phase_index;
				isReverseRegion = _isReverseRegion;
			}
			GeometricRegion& operator=(const GeometricRegion& n) {
				this->generate_step = n.generate_step;
				this->geometryProperty = n.geometryProperty;
				this->temperature = n.temperature;
				this->phaseIndex = n.phaseIndex;
				this->ellipSolid = n.ellipSolid;
				this->polyhedron = n.polyhedron;
				this->cylinder = n.cylinder;
				this->con = n.con;
				this->phi = n.phi;
				this->is_phi = n.is_phi;
				this->is_con = n.is_con;
				this->is_temp = n.is_temp;
				this->isReverseRegion = n.isReverseRegion;
				this->isNormalized = n.isNormalized;
				return *this;
			}
		};
		class PointSet {
		public:
			PointSet(size_t _generate_step = 0, size_t _phase_index = 0, double _temperature = 0) {
				init(_generate_step, _phase_index, _temperature);
				is_normalized = true;
			}
			void init(size_t _generate_step = 0, size_t _phase_index = 0, double _temperature = 0) {
				generate_step = _generate_step;
				phaseIndex = _phase_index;
				temperature = _temperature;
			}
			~PointSet() {
				points.clear();
			}
			std::vector<Point> points;
			std::vector<double> points_phi;
			size_t generate_step;
			size_t phaseIndex;
			Matrix1D<double> con;
			double temperature;
			bool is_normalized;
			bool is_phi = false;
			bool is_con = false;
			bool is_temp = false;
			bool is_point_in_set(Point point, int x_limit, BoundaryCondition x_down_bc, BoundaryCondition x_up_bc, int y_limit = 1,
				BoundaryCondition y_down_bc = BoundaryCondition::ADIABATIC, BoundaryCondition y_up_bc = BoundaryCondition::ADIABATIC,
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
			void add_point(double _x, double _y, double _z, double _phi) {
				points.push_back(Point(_x, _y, _z));
				points_phi.push_back(_phi);
			}
			void add_point(int _x, int _y, int _z, double _phi) {
				points.push_back(Point(_x, _y, _z));
				points_phi.push_back(_phi);
			}
			void add_point(size_t _x, size_t _y, size_t _z, double _phi) {
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
				this->phaseIndex = n.phaseIndex;
				this->con = n.con;
				this->points_phi = n.points_phi;
				this->is_normalized = n.is_normalized;
				this->is_phi = n.is_phi;
				this->is_con = n.is_con;
				this->is_temp = n.is_temp;
				return *this;
			}
		};
		struct NucleationBox {
		public:
			std::vector<GeometricRegion> geometry_box;
			std::vector<PointSet> point_set_box;
			NucleationBox() {

			};
			NucleationBox& operator=(const NucleationBox& n) {
				geometry_box = n.geometry_box;
				point_set_box = n.point_set_box;
				return *this;
			}
		};
		// - main
		inline NucleationBox nucleation_box;
		// - matrix
		inline size_t matrix_phi_index = 0;
		inline double matrix_phi_value = 0;
		inline bool is_matrix_phi = false;
		inline Matrix1D<REAL> matrix_con;
		inline bool is_matrix_con = false;
		inline double matrix_temperature = 0.0;
		inline bool is_matrix_temp = false;
	}
}