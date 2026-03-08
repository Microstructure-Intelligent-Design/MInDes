#pragma once
#include "VectorMatrix.h"
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
		PlasticPoint() {
			AvePlasticStrain = 0;
		};
	};
}
