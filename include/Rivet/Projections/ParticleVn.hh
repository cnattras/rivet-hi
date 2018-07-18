// -*- C++ -*-
#ifndef RIVET_ParticleVn_HH
#define RIVET_ParticleVn_HH

#include "Rivet/Projections/EventPlane.hh"
#include "Rivet/Particle.hh"
#include "Rivet/Event.hh"
#include "Rivet/Math/MathUtils.hh"

namespace Rivet
{


	/// @brief Project Out event Vn
	///
	///
	/// @todo
	class ParticleVn : public Projection
	{
		public:


			/// @name Constructors and destructors.
			//@{
			static const int SQRT_N = -1;
			/// Constructor
			///
			ParticleVn(const ChargedFinalState& cfs,const std::vector<int>& nOfInterest, int nBins = ParticleVn::SQRT_N)
			{

				addProjection(cfs,"CFS");
				addProjection(EventPlane(cfs,nOfInterest),"EventPlane");

				setName("ParticleVn");
				_nOfInterest = nOfInterest;
				_vn.resize(_nOfInterest.size());
				_nBins = nBins;
			}


			ParticleVn(const ChargedFinalState& cfs,int nOfInterest, int nBins = ParticleVn::SQRT_N)
			{

				addProjection(cfs,"CFS");
				addProjection(EventPlane(cfs,nOfInterest),"EventPlane");

				setName("ParticleVn");
				_nOfInterest.push_back(nOfInterest);
				_vn.resize(_nOfInterest.size());
				_nBins = nBins;
			}

			/// Clone on the heap.
			DEFAULT_RIVET_PROJ_CLONE(ParticleVn);

			/// Accessors
			const double GetVnOrdinal(int n = 0) const
			{

				if (n < (int)_vn.size() ) return _vn[n];
				else return -1.;
			}


			const double GetVn(int n = 2) const
			{
				for (unsigned int i = 0; i < _nOfInterest.size(); i++)
					{
						if (_nOfInterest[i] == n) return _vn[i];
					}
				return -1.;
			}

			const std::vector<double> GetVnList() const
			{
				return _vn;
			}

			const double GetEventPlaneOrdinal(int n = 0) const { return _epNvn[n]; }
		
			const double EventPlaneN(int N = 2) const { 
				for (unsigned int i = 0; i < _nOfInterest.size(); i++){
					if (_nOfInterest[i] == N) return _epNvn[i]; 
				}
				return -100.; 
				
			}
			const std::vector<double> EventPlanes() const { return _epNvn; }

		





		protected:

			/// Apply the projection to the event.
			void project(const Event& e)
			{
				_epNvn = apply<EventPlane>(e, "EventPlane").EventPlanes();
				const Particles& part = apply<ChargedFinalState>(e,"CFS").particles();
				for (unsigned int i = 0; i < _nOfInterest.size(); i++) {
					_vn[i] = _fitToCos(part,_epNvn[i],_nOfInterest[i]);
				}


			}

			/// Compare projections
			int compare(const Projection& p) const
			{
				return mkNamedPCmp(p, "CFS");
			}


			double deltaPhi2(double phi1, double phi2)
			{
				return mapAngle0To2Pi(phi1 - phi2);
			}


			/// cosine fit function
			double _fitToCos(const Particles& part, double EPn, int n)
			{
				int nBins;

				if (_nBins == -1) nBins = sqrt(part.size());
				else nBins = _nBins;

				std::vector<double> phiBins(nBins,0.0);

				double binWidth = 2*pi/nBins;


				foreach (const Particle &p, part)
				{
					if (p.pT() < 5. * GeV)	{
							phiBins[static_cast<int>(deltaPhi2(p.phi(),EPn)/binWidth)] += p.pT()*GeV;
						}
				}

				double integral = 0.0;
				double midOffset = binWidth * .5;
				double Vn = 0.0;

				foreach (const double &d, phiBins) integral += d;


				for (int i =0; i < nBins; i++)
					Vn += (phiBins[i]) * cos(n * (i * binWidth + midOffset ));


				Vn /= integral;
				return Vn;
			}

			std::vector<int> _nOfInterest;
			std::vector<double> _vn;
			std::vector<double> _epNvn;
			int _nBins;
	};


}


#endif
