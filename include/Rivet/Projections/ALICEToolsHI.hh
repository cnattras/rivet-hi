// -*- C++ -*-
#ifndef RIVET_ALICEToolsHI_HH
#define RIVET_ALICEToolsHI_HH

#include "fastjet/ClusterSequenceArea.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/ParticleVn.hh"
#include "Rivet/Projections/EventPlane.hh"
#include "Rivet/Particle.hh"
#include "Rivet/Event.hh"	
#include <algorithm>
#include <cmath>

namespace Rivet {


	/// @brief Project Out the Background of an event using ALICE methodology
	///
	///
	/// @todo 
	class ALICEToolsHI : public ParticleVn {
		public:


		/// @name Constructors and destructors.
		//@{

		/// Constructor
		///
		ALICEToolsHI( const ChargedFinalState& cfs, ///< Input particles
                      const FastJets& antiKtJets, ///< FastJets Object for your anti-kt Jets
                      double maxAbsEta=0.9, ///< Highest |eta| for your particles
					  std::vector<int> vn = {2,3} ///< Vn of interest to remove flow
					  );

		/// Clone on the heap.
		DEFAULT_RIVET_PROJ_CLONE(ALICEToolsHI);

		//@}


		/// Returns the median energy (minus top two) jet energy density
		const double MedianRho() const { return _backgroundMedianDensity;} 
		
		/// Returns the average energy density
		const double AvgRho() const { return _backgroundRho0;} 
		
		
		/// Returns the underlying event energy den at a given phi
		const double RhoLocal(double phi) const ;

		/// Get the nth order event plane
		const double EventPlaneN(int N = 2) const { 
			for (unsigned int i = 0; i < _nOfInterest.size(); i++){
				if (_nOfInterest[i] == N) return _epN[i]; 
			}
			return -100.; 
			
		}
		
		/// Get a vector of all calculated event planes
		const std::vector<double> EventPlanes() const { return _epN; }


		/// Get the reaction plane resolution for v2
		const double ReactionPlaneResolution2() const {
				return sqrt((_r2PosEta * _r2NegEta)/(_nEvent * _r2DifEta)); 
		}
		
		/// Get the reaction plane resolution for v3
		const double ReactionPlaneResolution3() const {
				return sqrt((_r3PosEta * _r3NegEta)/(_nEvent * _r3DifEta)); 
		}

	protected:
		
		/// Apply the projection to the event.
		void project(const Event& e);
		
		/// Compare projections
		int compare(const Projection& p) const;
		

		
		Particles _eventPart;
		
		/// ALICE specific
		double _backgroundMedianDensity;
		double _backgroundRho0;
		double _maxAbsEta;
		double _sin3R;
		double _sin2R;
		double _jetR;
		std::vector<double> _epN;
		
		
		// Resolution terms 
		double _r2PosEta = 0.0;
		double _r2NegEta = 0.0;
		double _r2DifEta = 0.0;
		
		double _r3PosEta = 0.0;
		double _r3NegEta = 0.0;
		double _r3DifEta = 0.0;
		double _nEvent = 0.0;
		
	};

			
}


#endif
