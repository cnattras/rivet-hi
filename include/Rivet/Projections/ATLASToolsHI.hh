// -*- C++ -*-
#ifndef RIVET_ATLASToolsHI_HH
#define RIVET_ATLASToolsHI_HH

#include "fastjet/ClusterSequenceArea.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/ParticleVn.hh"
#include "Rivet/Projections/EventPlane.hh"
#include "Rivet/Projections/DetectorPixels.hh"
#include "Rivet/Particle.hh"
#include "Rivet/Event.hh"
#include <algorithm>
#include <cmath>

namespace Rivet
{


	/// @brief Project Out the Background of an event using ATLAS method
	///
	///
	/// @todo
	class ATLASToolsHI : public FinalState
	{
		public:

			ATLASToolsHI() : ATLASToolsHI(7, 8) {}

			ATLASToolsHI( double trackJetPtQualifier,
			              double EMTowerQualifier
			            ) ;

			Particles calorimeterTowers(const Cut& c = Cuts::open()) const
			{
				return filter_select(_calorimeterTowers,c);
			}
			Particles calorimeterTowersRaw(const Cut& c = Cuts::open()) const
			{
				return filter_select(_calorimeterTowersRaw,c);
			}
			Particles calorimeterTowersInclGhosts(const Cut& c = Cuts::open()) const
			{
				return filter_select(_calorimeterTowersInclGhosts,c);
			}
			Jets calorimeterJets(const Cut& c = Cuts::open() ) const
			{
				return filter_select(_calorimeterJets, c);
			}
			const Jets trackJets() const
			{
				return _trackJets;
			}
			static double ghostPhi(const Particle& p)
			{
				return mapAngle0To2Pi(atan2( p.py(), p.px() ));
			}

			const double eventPlane2() const
			{
				return _psi2;
			}
			const double v2() const
			{
				return _v2;
			}
			/// Clone on the heap.
			DEFAULT_RIVET_PROJ_CLONE(ATLASToolsHI);


		protected:

			/// Apply the projection to the even_t.
			void project(const Event& e);


			/// Compare projections
			int compare(const Projection& p) const;

			double _minRap;
			double _maxRap;
			int _nRap;
			int _nRap2;
			int _nPhi;
			double _wRap;
			double _wPhi;


			std::vector<FourMomentum> _towerMom0;
			double _towerArea;
			std::vector<double> _rhoi; // rho i as a function of eta
			std::vector<int> _rhoiEntries;
			double _v2;
			double _psi2;

			std::vector<bool> _isSeed;
			Jets _calorimeterJets;
			Jets _trackJets;

			Particles _calorimeterTowers;
			Particles _calorimeterTowersInclGhosts;
			Particles _calorimeterTowersRaw;


			double _tJetPt;
			double _towEt;

		private:
			void _tagSeedJets1(const Jets& js);

			void _calculateRhoi1(const Jets& js);

			void _subtractInterimUE(Jets& js);
			
			FourMomentum _adjMom(const Particle& p);

			void _getSeedJets2(const Jets& caloJets, const Jets& trackJets, Jets& cSeeds, Jets& tSeeds);
			
			void _calculateRhoi2(const Particles& towers, const Jets& cSeeds, const Jets& tSeeds);
			
			void _calcV2Psi2(const Particles& towers, const Jets& cSeeds, const Jets& tSeeds);
			
			void _subtractUE();
			
			int _rapBin(double eta);

			void _finalizeJets(const Jets& js);
	};

}


#endif
