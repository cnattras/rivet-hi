// -*- C++ -*-
#include "Rivet/Projections/EventPlane.hh"

namespace Rivet {

	/////////////////////
	/// Constructor
	/////////////////////
	EventPlane::EventPlane( const FinalState& fs, ///< Particles from which to calculate Event Planes
	                        const std::vector<int>& nOfInterest, ///< nth order event planes to calculate
	                        const Cut& c  ///< Cuts to apply to \code{fs} before calculating Event Planes
	                      )
	{
		setName("EventPlane");
		addProjection(fs, "FS");
		_nOfInterest = nOfInterest;
		_epN.resize(_nOfInterest.size());
		_pCut = c;
	}


	/////////////////////
	/// Constructor
	/////////////////////
	EventPlane::EventPlane( const FinalState& fs, ///< Particles from which to calculate Event Planes
	                        int nOfInterest, ///< nth order event plane to calculate
	                        const Cut& c ///< Cuts to apply to \code{fs} before calculating Event Planes
	                      )
	{
		setName("EventPlane");
		addProjection(fs, "FS");
		_nOfInterest.push_back(nOfInterest);
		_epN.resize(_nOfInterest.size());
		_pCut = c;
	}


	void EventPlane::project(const Event& e)
	{

		const Particles& parts = apply<FinalState>(e, "FS").particles(_pCut);

		for (unsigned int i = 0; i < _nOfInterest.size(); i++)
		{
			double n = static_cast<double>(_nOfInterest[i]);
			double Qx = 0.0;
			double Qy = 0.0;
			foreach (const Particle& p, parts)
			{
				Qx += p.pT()*cos(n*p.phi());
				Qy += p.pT()*sin(n*p.phi());
			}
			_epN[i] = atan2(Qy,Qx)/ n;
		}

	}

	int EventPlane::compare(const Projection& p) const
	{
		const EventPlane& other = dynamic_cast<const EventPlane&>(p);
		if (other._nOfInterest.size() != _nOfInterest.size()) return UNDEFINED;
		for (int i = 0; i < (signed)_nOfInterest.size(); ++i)
		{
			if (other._nOfInterest[i] != _nOfInterest[i]) return UNDEFINED;
		}
		return _pCut == other._pCut ? EQUIVALENT : UNDEFINED;
	}




}
