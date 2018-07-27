// -*- C++ -*-
#include "Rivet/Projections/CMSToolsHI.hh"

namespace Rivet {

	CMSToolsHI::CMSToolsHI(const FinalState& fs, BackgroundMethod bgm, double kMultiplier) : _kMultiplier(kMultiplier)
	{
		setName("CMSToolsHI");
		addProjection(fs, "FS");
		_bgm = bgm;
		if (_bgm == CMS_NoisePedestal)
		{
			DetectorPixels HCALmid(DetectorPixels::CMS_HCALMid);
			DetectorPixels HCALend(DetectorPixels::CMS_HCALEndCap);
			declare(HCALmid, "HCALMID");
			declare(HCALend, "HCALEND");


			_minRapMid = HCALmid.RapidityMin();
			_maxRapMid = HCALmid.RapidityMax();
			_nPhiMid = HCALmid.PhiTowers();
			_wPhiMid = HCALmid.PhiWidth();
			_nRapMid = HCALmid.RapTowers();
			_wRapMid = HCALmid.RapWidth();
			_nRap2Mid = _nRapMid /2;
			_towerAreaMid = _wPhiMid * _wRapMid;

			_minRapEnd = HCALend.RapidityMin();
			_maxRapEnd = HCALend.RapidityMax();
			_nPhiEnd = HCALend.PhiTowers();
			_wPhiEnd = HCALend.PhiWidth();
			_nRapEnd = HCALend.RapTowers();
			_wRapEnd = HCALend.RapWidth();
			_nRap2End = _nRapEnd /2;
			_towerAreaEnd = _wPhiEnd * _wRapEnd;
		}
	}

	void CMSToolsHI::project(const Event& e)
	{
		_eventPart = applyProjection<FinalState>(e, "FS").particles();
		if (_bgm == CMS_NoisePedestal)
		{
			_dpMid = apply<DetectorPixels>(e,"HCALMID");
			_hcalMid = _dpMid.particles();
			_hcalEnd = apply<DetectorPixels>(e,"HCALEND").particles();

			_calculateEAvgRho();
			_adjustEnergies();
			_findJets();
			_calculateNewEAvgRho();
			_adjustEnergies();
			_finalizeOutputParticles();
		}
	}


	int CMSToolsHI::compare(const Projection& p) const
	{
		const CMSToolsHI& other = dynamic_cast<const CMSToolsHI&>(p);
		if (_kMultiplier != other._kMultiplier) return UNDEFINED;
		if (_bgm != other._bgm) return UNDEFINED;
		return mkNamedPCmp(p, "CFS");
	}

	int CMSToolsHI::_rapBinMid(double eta)
	{
		if (eta < 0) return static_cast<int>((eta+_maxRapMid)/_wRapMid);
		else return _nRap2Mid + static_cast<int>((eta-_minRapMid)/_wRapMid);
	}
	int CMSToolsHI::_rapBinEnd(double eta)
	{
		if (eta < 0) return static_cast<int>((eta+_maxRapEnd)/_wRapEnd);
		else return _nRap2End + static_cast<int>((eta-_minRapEnd)/_wRapEnd);
	}

	void CMSToolsHI::_calculateEAvgRho()
	{
		_eAvgMid = std::vector<double>(_nRapMid,0.0);
		_eRhoMid = std::vector<double>(_nRapMid,0.0);
		_eAvgEnd = std::vector<double>(_nRapEnd,0.0);
		_eRhoEnd = std::vector<double>(_nRapEnd,0.0);

		foreach (const Particle& p, _hcalMid)
		{
			_eAvgMid[_rapBinMid(p.eta())] +=  p.Et();
			_eRhoMid[_rapBinMid(p.eta())] +=  p.Et2();
			//cout << p.Et() << "\t";
		}
		cout << endl;
		for (int i = 0; i < _nRapMid; ++i)
		{
			_eAvgMid[i] /= _nPhiMid;
			_eRhoMid[i] /= _nPhiMid;
			_eRhoMid[i] -= _eAvgMid[i] * _eAvgMid[i];
			_eRhoMid[i] = sqrt(_eRhoMid[i]);
		}

		foreach (const Particle& p, _hcalEnd)
		{
			_eAvgEnd[_rapBinEnd(p.eta())] +=  p.Et();
			_eRhoEnd[_rapBinEnd(p.eta())] +=  p.Et2();
		}
		for (int i = 0; i < _nRapEnd; ++i)
		{
			_eAvgEnd[i] /= _nPhiEnd;
			_eRhoEnd[i] /= _nPhiEnd;
			_eRhoEnd[i] -= _eAvgEnd[i] * _eAvgEnd[i];
			_eRhoEnd[i] = sqrt(_eRhoEnd[i]);
		}

	}

	void CMSToolsHI::_adjustEnergies()
	{
		_hcalMidAdj.clear();
		_hcalEndAdj.clear();
		foreach (const Particle& p, _hcalMid)
		{
			double theta = p.theta();
			double phi = ghostPhi(p);
			double Et = p.Et();
			double Eadj = 1e-20; // ghost energy
			double EtUE = _eAvgMid[_rapBinMid(p.eta())] + _kMultiplier * _eRhoMid[_rapBinMid(p.eta())];
			if (Et-EtUE > 0.0) Eadj = (Et - EtUE) / sin(theta);

			Particle pNew(p.pid(), FourMomentum(Eadj,Eadj*cos(phi)*sin(theta),Eadj*sin(phi)*sin(theta),Eadj*cos(theta)));
			_hcalMidAdj.push_back(pNew);
		}
		foreach (const Particle& p, _hcalEnd)
		{
			double theta = p.theta();
			double phi = ghostPhi(p);
			double Et = p.Et();
			double Eadj = 1e-20; // ghost energy
			double EtUE = _eAvgEnd[_rapBinEnd(p.eta())] + _kMultiplier * _eRhoEnd[_rapBinEnd(p.eta())];
			if (Et-EtUE > 0.0) Eadj = (Et - EtUE) / sin(theta);

			Particle pNew(p.pid() + _hcalMid.size(), FourMomentum(Eadj,Eadj*cos(phi)*sin(theta),Eadj*sin(phi)*sin(theta),Eadj*cos(theta)));
			_hcalEndAdj.push_back(pNew);
		}
	}

	void CMSToolsHI::_findJets()
	{
		FastJets fj(_dpMid, FastJets::CMSCONE, 0.5, JetAlg::ALL_MUONS, JetAlg::ALL_INVISIBLES);
		Particles bothPixels(_hcalMidAdj);
		bothPixels.insert(bothPixels.end(), _hcalEndAdj.begin(), _hcalEndAdj.end());
		fj.calc(bothPixels);
		_iConeJets = fj.jets(Cuts::pt > 10 * GeV);

	}


	void CMSToolsHI::_calculateNewEAvgRho()
	{
		_eAvgMid = std::vector<double>(_nRapMid,0.0);
		_eRhoMid = std::vector<double>(_nRapMid,0.0);
		_eAvgEnd = std::vector<double>(_nRapEnd,0.0);
		_eRhoEnd = std::vector<double>(_nRapEnd,0.0);

		std::vector<int> nBinsMid(_nRapMid, _nPhiMid);
		std::vector<int> nBinsEnd(_nRapEnd, _nPhiEnd);


		foreach (const Particle& p, _hcalMid)
		{
			_eAvgMid[_rapBinMid(p.eta())] +=  p.Et();
			_eRhoMid[_rapBinMid(p.eta())] +=  p.Et2();
		}

		foreach (const Particle& p, _hcalEnd)
		{
			_eAvgEnd[_rapBinEnd(p.eta())] +=  p.Et();
			_eRhoEnd[_rapBinEnd(p.eta())] +=  p.Et2();
		}

		int bin;
		int hcalSize = _hcalMid.size();
		foreach (const Jet &j, _iConeJets)
		{
			foreach (const Particle& p, j.particles())
			{
				if (p.pid() < hcalSize)
				{
					bin = _rapBinMid(p.eta());
					_eAvgMid[bin] -=  _hcalMid[p.pid()].Et();
					_eRhoMid[bin] -=  _hcalMid[p.pid()].Et2();
		
					--nBinsMid[bin];
				}
				else
				{
					bin = _rapBinEnd(p.eta());
					_eAvgEnd[bin] -=  _hcalEnd[p.pid()-hcalSize].Et();
					_eRhoEnd[bin] -=  _hcalEnd[p.pid()-hcalSize].Et2();
					--nBinsEnd[bin];
				}

			}
		}

		for (int i = 0; i < _nRapMid; ++i)
		{
			_eAvgMid[i] /= nBinsMid[i];
			_eRhoMid[i] /= nBinsMid[i];
			_eRhoMid[i] -= _eAvgMid[i] * _eAvgMid[i];
			_eRhoMid[i] = sqrt(_eRhoMid[i]);

		}


		for (int i = 0; i < _nRapEnd; ++i)
		{
			_eAvgEnd[i] /= nBinsEnd[i];
			_eRhoEnd[i] /= nBinsEnd[i];
			_eRhoEnd[i] -= _eAvgEnd[i] * _eAvgEnd[i];
			_eRhoEnd[i] = sqrt(_eRhoEnd[i]);

		}

	}

	void CMSToolsHI::_finalizeOutputParticles()
	{
		_theParticles.clear();

		// output each tower with a PID of 2212, proton
		// prevents issues with fastjet making cuts
		foreach (const Particle &p, _hcalMidAdj)
		{
			if (p.Et() > 1e-10)
			{
				_theParticles.push_back(Particle(2212, p.momentum()));
			}
		}
		foreach (const Particle &p, _hcalEndAdj)
		{
			if (p.Et() > 1e-10)
			{
				_theParticles.push_back(Particle(2212, p.momentum()));
			}
		}
	}


}
