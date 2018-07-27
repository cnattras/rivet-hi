// -*- C++ -*-
#ifndef RIVET_CMSToolsHI_HH
#define RIVET_CMSToolsHI_HH


#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Particle.hh"
#include "Rivet/Event.hh"
#include "Rivet/Projections/DetectorPixels.hh"

namespace Rivet
{


	/// @brief Project Out the CMS method background of an event
	///
	///
	/// @todo
	class CMSToolsHI : public FinalState
	{
		public:

			enum BackgroundMethod
			{
			    CMS_Reflection,
			    CMS_NoisePedestal // dx.doi.org/10.1088/1126-6708/2006/05/026
			};

			/// @name Constructors and destructors.
			//@{

			/// Constructor
			///
			CMSToolsHI(const FinalState& fs, BackgroundMethod bgm, double kMultiplier = 1);


			/// Clone on the heap.
			DEFAULT_RIVET_PROJ_CLONE(CMSToolsHI);

			/// Get particles in cone in reflected eta from jetOfInterest
			const Particles CMSReflectionParticles(const Jet& jetOfInterest, double R_Param) const
			{
				Particles particles;
				particles = getConeParticles(-jetOfInterest.eta(),jetOfInterest.phi(), R_Param);
				return particles;
			}
			
			/// Get the summed FourMomentum of the reflected particles in the cone of reflected eta for jetOfInterest
			const FourMomentum CMSReflectionFourMomentum(const Jet& jetOfInterest, double R_Param) const
			{
				FourMomentum vec(0,0,0,0);
				Particles particles = CMSReflectionParticles(jetOfInterest, R_Param);
				foreach (const Particle &p, particles)
				{
					vec += p.momentum();
				}
				return vec;
			}

			/// Get all particles within R of jet's centroid
			const Particles getConeParticles(const Jet& jet,double R) const
			{
				return getConeParticles(jet.eta(), jet.phi(), R);
			}
			
			/// Get all particles within R of jet's centroid			
			const Particles getConeParticles(double eta, double phi,double R) const
			{
				Particles conePart;
				foreach (const Particle &p, _eventPart)
				{
					if (deltaR(p,eta,phi) < R) conePart.push_back(p);
				}
				return conePart;
			}
			static double ghostPhi(const Particle& p)
			{
				return mapAngle0To2Pi(atan2( p.py(), p.px() ));
			}


		protected:

			/// Apply the projection to the event.
			void project(const Event& e);


			int compare(const Projection& p) const;

			BackgroundMethod _bgm;
			double _kMultiplier;

		private:

			int _rapBinMid(double eta);
			int _rapBinEnd(double eta);
			void _calculateEAvgRho();
			void _adjustEnergies();
			void _findJets();
			void _calculateNewEAvgRho();
			void _finalizeOutputParticles();

			Jets _iConeJets;

			std::vector<double> _eAvgMid;
			std::vector<double> _eAvgEnd;
			std::vector<double> _eRhoMid;
			std::vector<double> _eRhoEnd;

			Particles _hcalMid;
			Particles _hcalEnd;
			Particles _hcalMidAdj;
			Particles _hcalEndAdj;

			double _minRapMid;
			double _maxRapMid;
			int _nRapMid;
			int _nRap2Mid;
			int _nPhiMid;
			double _wRapMid;
			double _wPhiMid;
			double _towerAreaMid;

			double _minRapEnd;
			double _maxRapEnd;
			int _nRapEnd;
			int _nRap2End;
			int _nPhiEnd;
			double _wRapEnd;
			double _wPhiEnd;
			double _towerAreaEnd;

			Particles _eventPart;

			DetectorPixels _dpMid;
	};



}


#endif
