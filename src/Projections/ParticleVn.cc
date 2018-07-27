// -*- C++ -*-
#include "Rivet/Projections/ParticleVn.hh"

namespace Rivet {


	ParticleVn::ParticleVn( const FinalState& fs, ///< Particles to use to calculate Vn
	                        const std::vector<int>& nOfInterest, ///< nth order coefficients to calculate
	                        int nBins ///< Number of bins to use for granularity
	                      )
	{

		addProjection(fs,"FS");
		addProjection(EventPlane(fs,nOfInterest),"EventPlane");

		setName("ParticleVn");
		_nOfInterest = nOfInterest;
		_vn.resize(_nOfInterest.size());
		_nBins = nBins;
	}


	ParticleVn::ParticleVn( const FinalState& fs, ///< Particles to use to calculate Vn
	                        int nOfInterest, ///< nth order coefficient to calculate
	                        int nBins ///< Number of bins to use for granularity
	                      )
	{

		addProjection(fs,"FS");
		addProjection(EventPlane(fs,nOfInterest),"EventPlane");

		setName("ParticleVn");
		_nOfInterest.push_back(nOfInterest);
		_vn.resize(_nOfInterest.size());
		_nBins = nBins;
	}

	void ParticleVn::project(const Event& e)
	{
		_epNvn = apply<EventPlane>(e, "EventPlane").EventPlanes();
		const Particles& part = apply<FinalState>(e,"FS").particles();
		for (unsigned int i = 0; i < _nOfInterest.size(); i++)
		{
			_vn[i] = _fitToCos(part,_epNvn[i],_nOfInterest[i]);
		}


	}

	int ParticleVn::compare(const Projection& p) const
	{
		const ParticleVn& other = dynamic_cast<const ParticleVn&>(p);
		if (other._nOfInterest.size() != _nOfInterest.size()) return UNDEFINED;
		for (int i = 0; i < (signed)_nOfInterest.size(); ++i)
		{
			if (other._nOfInterest[i] != _nOfInterest[i]) return UNDEFINED;
		}
		if (_nBins != other._nBins) return UNDEFINED;
		return mkNamedPCmp(p, "FS");
	}
	double ParticleVn::_fitToCos(const Particles& part, double EPn, int n)
	{
		int nBins;

		if (_nBins == -1) nBins = sqrt(part.size());
		else nBins = _nBins;

		std::vector<double> phiBins(nBins,0.0);

		double binWidth = 2*pi/nBins;


		foreach (const Particle &p, part)
		{
			if (p.pT() < 5. * GeV)
			{
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




}
