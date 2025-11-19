#pragma once
#include "VectorMatrix.h"

namespace pf{
	enum RotationGauge { RG_XYX = 0, RG_XZX, RG_YXY, RG_YZY, RG_ZXZ, RG_ZYZ, RG_XYZ, RG_XZY, RG_YXZ, RG_YZX, RG_ZXY, RG_ZYX };

	namespace RotationMatrix {
		Matrix3x3 rotate_zero();
		Matrix3x3 rotate_x(REAL radian);
		Matrix3x3 rotate_y(REAL radian);
		Matrix3x3 rotate_z(REAL radian);
		Matrix3x3 rotationMatrix(Vector3 radian, RotationGauge rotationGauge = RG_XYX);
		Matrix3x3 rotationMatrix_XYX(Vector3 radian);
		Matrix3x3 rotationMatrix_XZX(Vector3 radian);
		Matrix3x3 rotationMatrix_YXY(Vector3 radian);
		Matrix3x3 rotationMatrix_YZY(Vector3 radian);
		Matrix3x3 rotationMatrix_ZXZ(Vector3 radian);
		Matrix3x3 rotationMatrix_ZYZ(Vector3 radian);
	}
}