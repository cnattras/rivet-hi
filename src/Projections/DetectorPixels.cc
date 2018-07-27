// -*- C++ -*-
#include "Rivet/Projections/DetectorPixels.hh"

namespace Rivet {

	void DetectorPixels::_init(const Cut& c , double minAbsRap, double maxAbsRap, double towerWidthRap, double towerWidthPhi, bool includeNeutrals)
	{
		setName("DetectorPixels");
		_inclNeut = includeNeutrals;
		if (_inclNeut)
		{
			FinalState fs(Cuts::pT > 0.15 * GeV && c);
			addProjection(fs, "FS");
		}
		else
		{
			ChargedFinalState cfs(Cuts::pT > 0.15 * GeV && c);
			addProjection(cfs, "CFS");
		}

		_nPhi = static_cast<int>(2*pi / towerWidthPhi);
		_wPhi = 2.*pi / _nPhi;
		_nRap = 2*static_cast<int>((maxAbsRap-minAbsRap)/towerWidthRap);
		_wRap = (2*(maxAbsRap-minAbsRap))/_nRap;

		_minRap = minAbsRap;
		_maxRap = maxAbsRap;
		_nRap2 = _nRap /2;
		_towerMom0 = std::vector<FourMomentum>(_nPhi*_nRap,FourMomentum(0,0,0,0));
		_cuts = c;
		_theParticles = Particles(_nPhi * _nRap);
	}

	void DetectorPixels::project(const Event& e)
	{
		Particles parts;
		if (_inclNeut)
		{
			parts = apply<FinalState>(e,"FS").particles();
		}
		else
		{
			parts = apply<ChargedFinalState>(e,"CFS").particles();

		}
		std::vector<FourMomentum> towerMom = _towerMom0;

		int binRap2 =  _nRap / 2;

		foreach (const Particle& p, parts)
		{
			int binPhi = static_cast<int>(ghostPhi(p)/_wPhi);

			int binRap;
			double pEta = p.eta();
			if (pEta < 0) binRap = static_cast<int>((p.eta()+_maxRap)/_wRap);
			else binRap = binRap2 + static_cast<int>((p.eta()-_minRap)/_wRap);
			towerMom[binPhi + (_nPhi * binRap)] += p.mom();
		}



		for (int i = 0; i < _nPhi * _nRap ; ++i)
		{
			int binRap = i / _nPhi;
			int binPhi = i % _nPhi;

			double rap;

			if (binRap < binRap2) rap = -1 * _maxRap + ((binRap+0.5) * _wRap) ;
			else rap =  _minRap + ((binRap+0.5-binRap2) * _wRap);



			double phi = (binPhi+0.5) * _wPhi;
			double p = towerMom[i].p();
			double E = towerMom[i].E();
			if (E < 1e-20) E =1e-20; // for non-zero towers, make into ghosts with ~0 Et
			if (p < 1e-20) p =1e-20; // ghost 0 < .1 MeV leads to miscalculated phi
			double theta = 2*atan(exp(-rap));

			FourMomentum mom(E,p*cos(phi)*sin(theta),p*sin(phi)*sin(theta),p*cos(theta));

			_theParticles[i] = Particle(i, mom, FourVector(0,0,0,0));

		}


	}
	DetectorPixels::DetectorPixels( Detector det ///< Detector to mock up
	                              )
	{
		switch (det)
		{
		case ATLAS_MidRapidity:
			_init(Cuts::abseta < 2.5, 0, 2.5, 0.1, 0.1, true);
			break;
		case ATLAS_FCal:
			_init(Cuts::abseta < 4.9 && Cuts::abseta > 3.2, 3.2, 4.9, 0.2, 0.2, true);
			break;
		case CMS_HCALMid:
			_init(Cuts::abseta < 1.74, 0, 1.74, 0.087, 0.087, true);
			break;
		case CMS_HCALEndCap:
			_init(Cuts::abseta < 3 && Cuts::abseta >= 1.74, 1.74, 3, 0.174, 0.174, true);
			break;
		case ALICE_MidRapidity:
			_init(Cuts::abseta < 0.9, 0, 0.9, 0.1, 0.1, false); // TODO look up true granularity
			break;
		case ALICE_V0Detectors:
			_init(Cuts::abseta < 5.1 && Cuts::abseta > 2.8, 2.8, 5.1, 0.2, 0.2, false); // TODO look up true granularity
			// and fix for assymetry of V0A and V0C
			break;
		}
	}

	int DetectorPixels::compare(const Projection& p) const
	{
		const DetectorPixels& other = dynamic_cast<const DetectorPixels&>(p);
		if (other._wPhi != _wPhi) return UNDEFINED;
		if (other._nPhi != _nPhi) return UNDEFINED;
		if (other._wRap != _wRap) return UNDEFINED;
		if (other._nRap != _nRap) return UNDEFINED;
		return EQUIVALENT;
	}


}
