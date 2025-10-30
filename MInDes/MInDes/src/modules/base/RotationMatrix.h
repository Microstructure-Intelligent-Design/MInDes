#pragma once
#include "VectorMatrix.h"

namespace pf{
	enum RotationGauge { RG_XYX = 0, RG_XZX, RG_YXY, RG_YZY, RG_ZXZ, RG_ZYZ, RG_XYZ, RG_XZY, RG_YXZ, RG_YZX, RG_ZXY, RG_ZYX };

	namespace RotationMatrix {
		inline Matrix3x3 rotate_zero() {
			Matrix3x3 Rx;
			Rx.set_to_unity();
			return Rx;
		}
		inline Matrix3x3 rotate_x(REAL radian) {
			Matrix3x3 Rx;
			Rx.set_to_zero();
			Rx(0, 0) = 1;
			Rx(1, 1) = cos(radian);
			Rx(1, 2) = sin(radian);
			Rx(2, 1) = -sin(radian);
			Rx(2, 2) = cos(radian);
			return Rx;
		}
		inline Matrix3x3 rotate_y(REAL radian) {
			Matrix3x3 Ry;
			Ry.set_to_zero();
			Ry(1, 1) = 1;
			Ry(0, 0) = cos(radian);
			Ry(0, 2) = -sin(radian);
			Ry(2, 0) = sin(radian);
			Ry(2, 2) = cos(radian);
			return Ry;
		}
		inline Matrix3x3 rotate_z(REAL radian) {
			Matrix3x3 Rz;
			Rz.set_to_zero();
			Rz(2, 2) = 1;
			Rz(0, 0) = cos(radian);
			Rz(0, 1) = sin(radian);
			Rz(1, 0) = -sin(radian);
			Rz(1, 1) = cos(radian);
			return Rz;
		}
		inline Matrix3x3 rotationMatrix(Vector3 radian, RotationGauge rotationGauge = RG_XYX) {
			Matrix3x3 RM;
			RM.set_to_unity();
			switch (rotationGauge)
			{
			case pf::RG_XYX:
				RM = rotate_x(radian[2]) * (rotate_y(radian[1]) * rotate_x(radian[0]));
				return RM;
			case pf::RG_XZX:
				RM = rotate_x(radian[2]) * (rotate_z(radian[1]) * rotate_x(radian[0]));
				return RM;
			case pf::RG_YXY:
				RM = rotate_y(radian[2]) * (rotate_x(radian[1]) * rotate_y(radian[0]));
				return RM;
			case pf::RG_YZY:
				RM = rotate_y(radian[2]) * (rotate_z(radian[1]) * rotate_y(radian[0]));
				return RM;
			case pf::RG_ZXZ:
				RM = rotate_z(radian[2]) * (rotate_x(radian[1]) * rotate_z(radian[0]));
				return RM;
			case pf::RG_ZYZ:
				RM = rotate_z(radian[2]) * (rotate_y(radian[1]) * rotate_z(radian[0]));
				return RM;
			case pf::RG_XYZ:
				RM = rotate_z(radian[2]) * (rotate_y(radian[1]) * rotate_x(radian[0]));
				return RM;
			case pf::RG_XZY:
				RM = rotate_y(radian[2]) * (rotate_z(radian[1]) * rotate_x(radian[0]));
				return RM;
			case pf::RG_YXZ:
				RM = rotate_z(radian[2]) * (rotate_x(radian[1]) * rotate_y(radian[0]));
				return RM;
			case pf::RG_YZX:
				RM = rotate_x(radian[2]) * (rotate_z(radian[1]) * rotate_y(radian[0]));
				return RM;
			case pf::RG_ZXY:
				RM = rotate_y(radian[2]) * (rotate_x(radian[1]) * rotate_z(radian[0]));
				return RM;
			case pf::RG_ZYX:
				RM = rotate_x(radian[2]) * (rotate_y(radian[1]) * rotate_z(radian[0]));
				return RM;
			default:
				return RM;
			}
		}
		inline Matrix3x3 rotationMatrix_XYX(Vector3 radian) {
			Matrix3x3 RM;
			RM = rotate_x(radian[2]) * (rotate_y(radian[1]) * rotate_x(radian[0]));
			return RM;
		}
		inline Matrix3x3 rotationMatrix_XZX(Vector3 radian) {
			Matrix3x3 RM;
			RM = rotate_x(radian[2]) * (rotate_z(radian[1]) * rotate_x(radian[0]));
			return RM;
		}
		inline Matrix3x3 rotationMatrix_YXY(Vector3 radian) {
			Matrix3x3 RM;
			RM = rotate_y(radian[2]) * (rotate_x(radian[1]) * rotate_y(radian[0]));
			return RM;
		}
		inline Matrix3x3 rotationMatrix_YZY(Vector3 radian) {
			Matrix3x3 RM;
			RM = rotate_y(radian[2]) * (rotate_z(radian[1]) * rotate_y(radian[0]));
			return RM;
		}
		inline Matrix3x3 rotationMatrix_ZXZ(Vector3 radian) {
			Matrix3x3 RM;
			RM = rotate_z(radian[2]) * (rotate_x(radian[1]) * rotate_z(radian[0]));
			return RM;
		}
		inline Matrix3x3 rotationMatrix_ZYZ(Vector3 radian) {
			Matrix3x3 RM;
			RM = rotate_z(radian[2]) * (rotate_y(radian[1]) * rotate_z(radian[0]));
			return RM;
		}
	}
}