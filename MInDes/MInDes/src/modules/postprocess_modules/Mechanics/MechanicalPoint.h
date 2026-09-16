#pragma once
#include "../../base/VectorMatrix.h"
namespace pf {
	class ElasticPoint {
	public:
		vStress       Stress;                                        ///< Storage for stress
		vStrain       Strain;                                        ///< Storage for strain
		vStrain       StrainIncrement;                               ///< Storage for strain increment
		vStrain	      VirtualEigenStrain;                            ///< Storage for virtual eigenstrain
		vStrain       EffectiveEigenStrain;                         ///< Storage for effective eigenstrain
		Matrix6x6     EffectiveElasticConstant;                     ///< Storage for effective elastic constant
		void operator=(const ElasticPoint& n) {
			Stress = n.Stress;
			Strain = n.Strain;
			VirtualEigenStrain = n.VirtualEigenStrain;
			StrainIncrement = n.StrainIncrement;
			EffectiveEigenStrain = n.EffectiveEigenStrain;
			EffectiveElasticConstant = n.EffectiveElasticConstant;
		}
		//  (+)
		ElasticPoint operator+(const ElasticPoint& other) const {
			ElasticPoint result = *this;
			result += other;
			return result;
		}
		//  (-)
		ElasticPoint operator-(const ElasticPoint& other) const {
			ElasticPoint result = *this;
			result -= other;
			return result;
		}
		//  (*)
		ElasticPoint operator*(const REAL& other) const {
			ElasticPoint result = *this;
			result *= other;
			return result;
		}
		//  (/)
		ElasticPoint operator/(const REAL& other) const {
			ElasticPoint result = *this;
			result /= other;
			return result;
		}
		//  (+=)
		ElasticPoint& operator+=(const ElasticPoint& other) {
			Stress += other.Stress;
			Strain += other.Strain;
			StrainIncrement += other.StrainIncrement;
			VirtualEigenStrain += other.VirtualEigenStrain;
			EffectiveEigenStrain += other.EffectiveEigenStrain;
			EffectiveElasticConstant += other.EffectiveElasticConstant;
			return *this;
		}
		//  (-=)
		ElasticPoint& operator-=(const ElasticPoint& other) {
			Stress -= other.Stress;
			Strain -= other.Strain;
			StrainIncrement -= other.StrainIncrement;
			VirtualEigenStrain -= other.VirtualEigenStrain;
			EffectiveEigenStrain -= other.EffectiveEigenStrain;
			EffectiveElasticConstant -= other.EffectiveElasticConstant;
			return *this;
		}
		//  (*=)
		ElasticPoint& operator*=(const REAL& other) {
			Stress *= other;
			Strain *= other;
			StrainIncrement *= other;
			VirtualEigenStrain *= other;
			EffectiveEigenStrain *= other;
			EffectiveElasticConstant *= other;
			return *this;
		}
		//  (/=)
		ElasticPoint& operator/=(const REAL& other) {
			Stress /= other;
			Strain /= other;
			StrainIncrement /= other;
			VirtualEigenStrain /= other;
			EffectiveEigenStrain /= other;
			EffectiveElasticConstant /= other;
			return *this;
		}
		ElasticPoint() {
			
		};
	};
	class PlasticPoint {
	public:
		vStrain       PlasticStrain;                         ///< Storage for plastic strain
		REAL          AvePlasticStrain;                      ///< Storage for average plastic strain
		void operator=(const PlasticPoint& n) {
			PlasticStrain = n.PlasticStrain;
			AvePlasticStrain = n.AvePlasticStrain;
		}
		//  (+)
		PlasticPoint operator+(const PlasticPoint& other) const {
			PlasticPoint result = *this;
			result += other;
			return result;
		}
		//  (-)
		PlasticPoint operator-(const PlasticPoint& other) const {
			PlasticPoint result = *this;
			result -= other;
			return result;
		}
		//  (*) 
		PlasticPoint operator*(const REAL& other) const {
			PlasticPoint result = *this;
			result *= other;
			return result;
		}
		//  (/) 
		PlasticPoint operator/(const REAL& other) const {
			PlasticPoint result = *this;
			result /= other;
			return result;
		}
		//  (+=)
		PlasticPoint& operator+=(const PlasticPoint& other) {
			PlasticStrain += other.PlasticStrain;
			AvePlasticStrain += other.AvePlasticStrain;
			return *this;
		}
		//  (-=)
		PlasticPoint& operator-=(const PlasticPoint& other) {
			PlasticStrain -= other.PlasticStrain;
			AvePlasticStrain -= other.AvePlasticStrain;
			return *this;
		}
		//  (*=)
		PlasticPoint& operator*=(const REAL& other) {
			PlasticStrain *= other;
			AvePlasticStrain *= other;
			return *this;
		}
		//  (/=)
		PlasticPoint& operator/=(const REAL& other) {
			PlasticStrain /= other;
			AvePlasticStrain /= other;
			return *this;
		}
		PlasticPoint() {
			AvePlasticStrain = 0;
		};
	};
}
