// -*- C++ -*-
#ifndef RIVET_DetectorPixels_HH
#define RIVET_DetectorPixels_HH

#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Particle.hh"
#include "Rivet/Event.hh"

namespace Rivet
{


	/// @brief Project Out the Background of an event
	///
	///
	/// @todo
	class DetectorPixels : public FinalState
	{
		public:


			enum Detector
			{
			    ATLAS_MidRapidity,
			    ATLAS_FCal,
			    CMS_HCALMid,
			    CMS_HCALEndCap,
			    ALICE_V0Detectors,
			    ALICE_MidRapidity
			};

			/// @name Constructors and destructors.
			//@{

			/// Constructor
			DetectorPixels() : DetectorPixels( Cuts::abseta < 2.1, 0.,  2.1, 0.1,  0.1, false) {}

			DetectorPixels( const Cut& c, ///< Cut for input particles
			                double minAbsRap = 0, ///< Minimum absolute rapidiy of detector
			                double maxAbsRap = 2.1, ///< Minimum absolute rapidiy of detector
			                double towerWidthRap = 0.1, ///< Tower pixel width in rapidity
			                double towerWidthPhi = 0.1, ///< Tower pixel width in phi
			                bool includeNeutrals = false   ///< Include neutral particle energy
			              )
			{
				_init(c,minAbsRap,maxAbsRap,towerWidthRap,towerWidthPhi, includeNeutrals);
			}

			/// Constructor from specific detector
			DetectorPixels( Detector det ); ///< Detector to mock up



			/// Clone on the heap.
			DEFAULT_RIVET_PROJ_CLONE(DetectorPixels);

			//@}

			/// Detector pixel width in \f$\phi\f$
			const double PhiWidth() const
			{
				return _wPhi;
			}
			/// Detector pixel width in \f$\eta\f$
			const double RapWidth() const
			{
				return _wRap;
			}
			/// Number of pixels in \f$\phi\f$ direction
			const int PhiTowers() const
			{
				return _nPhi;
			}
			/// Number of pixels in \f$\eta\f$ direction
			const int RapTowers() const
			{
				return _nRap;
			}
			/// Minimum \f$|\eta|\f$ of the detector
			const double RapidityMin() const
			{
				return _minRap;
			}
			/// Maximum \f$|\eta|\f$ of the detector
			const double RapidityMax() const
			{
				return _maxRap;
			}
			/// Pixel area in \f$\phi-eta\f$ space
			const double TowerArea() const
			{
				return _wPhi * _wRap;
			}
			/// Get applied cuts to input ChargedFinalState
			const Cut Cuts() const
			{
				return _cuts;
			}
			/// Get the correct \f$\phi\f$ for all particles, including those less than 1 eV
			static double ghostPhi(const Particle& p)
			{
				return mapAngle0To2Pi(atan2( p.py(), p.px() ));
			}


		protected:

			/// Initialize the class
			void _init(const Cut& c , double minAbsRap, double maxAbsRap, double towerWidthRap, double towerWidthPhi, bool includeNeutrals);

			/// Apply the projection to the event.
			void project(const Event& e);

			/// Compare projections
			int compare(const Projection& p) const;


			std::vector<FourMomentum> _towerMom0;
			double _wPhi;
			int _nPhi;
			double _wRap;
			int _nRap;
			int _nRap2;
			double _minRap;
			double _maxRap;
			Cut _cuts;
			bool _inclNeut;


	};



}


#endif
