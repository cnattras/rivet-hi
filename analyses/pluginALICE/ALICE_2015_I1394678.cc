// -*- C++ -*-a
#include "fastjet/ClusterSequenceArea.hh"
#include "Rivet/HeavyIonAnalysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/ALICEToolsHI.hh"
#include <iostream>

namespace Rivet {


  /// @brief Add a short analysis description here
  class ALICE_2015_I1394678 : public HeavyIonAnalysis {
  public:

    /// Constructor
    ALICE_2015_I1394678(): HeavyIonAnalysis("ALICE_2015_I1394678"){}


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      HeavyIonAnalysis::init();
      // Initialise and register projections
      addCentralityMethod(HeavyIonAnalysis::ImpactParameter, 300,"ImpactParameterMethod");

      ChargedFinalState cfs(Cuts::abseta < 0.9);
      declare(cfs, "Tracks");

      FastJets fj(cfs, FastJets::ANTIKT, 0.2);
      fj.useJetArea(new fastjet::AreaDefinition(fastjet::active_area,
        fastjet::GhostedAreaSpec(0.7, 1, 1./200.)));
      declare(fj, "Jets");

      ALICEToolsHI ath(cfs, fj);
      declare(ath, "ATH");


      // Book histograms
      _histV2ChJet1   = bookScatter2D("d01-x01-y01");
      _histV2ChJet2   = bookScatter2D("d02-x01-y01");

      _histInPlane1   = bookHisto1D("0-5_centrality_in", refData(1,1,1));
      _histOutPlane1  = bookHisto1D("0-5_centrality_out", refData(1,1,1));

      _histInPlane2   = bookHisto1D("30-50_centrality_in", refData(2,1,1));
      _histOutPlane2  = bookHisto1D("30-50_centrality_out", refData(2,1,1));

      _histadd1       = bookHisto1D("add1", refData(1,1,1));
      _histsub1       = bookHisto1D("sub1", refData(1,1,1));

      _histadd2       = bookHisto1D("add2", refData(2,1,1));
      _histsub2       = bookHisto1D("sub2", refData(2,1,1));

      cout << "init\n";
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {

      const double c = centrality(event, "ImpactParameterMethod");

      if((c < 0.) || (c > 100.)){
        vetoEvent;
      }
      /// @todo Do the event by event analysis here
      FastJets fastjets = apply<FastJets>(event, "Jets");
      const Jets& jets = fastjets.jetsByPt(Cuts::pT > 0.15*GeV && Cuts::pT < 100*GeV && Cuts::abseta < 0.7);
      MSG_DEBUG("Jet multiplicity = " << jets.size());

      const auto clust_seq_area = fastjets.clusterSeqArea();

      const Particles& trks = apply<ChargedFinalState>(event, "Tracks").particles();
      MSG_DEBUG("Track multiplicity = " << trks.size());
      const ALICEToolsHI &ath = apply<ALICEToolsHI>(event,"ATH");
      //Find Event Plane 2 and 3
      const double EP2 = ath.EventPlaneN(2);

     foreach (const Jet& j, jets){

      const double area = clust_seq_area->area(j);

      double ptJet = j.pT() - (ath.RhoLocal(j.phi())*area);


      // plane determination for 0-5% centrality
      if (c >= 0 && c <= 5){

        // In-Plane criteria
        if ((deltaPhi(j,EP2) < (pi/4)) || (deltaPhi(j,EP2) > (3*pi/4))) {
          _histInPlane1->fill(ptJet, event.weight());
        }
        // Out-of-Plane criteria
        else {
          _histOutPlane1->fill(ptJet, event.weight());
        }
      }

      // plane determination for 30-50% centrality
      if (c >= 30 && c <= 50){


        // In-Plane criteria
        if ((deltaPhi(j,EP2) < (pi/4)) || (deltaPhi(j,EP2) > (3*pi/4))) {
          _histInPlane2->fill(ptJet, event.weight());
        }
        // Out-of-Plane criteria
        else {
          _histOutPlane2->fill(ptJet, event.weight());
        }
      }
    }

    _resolution2 = ath.ReactionPlaneResolution2();
    _resolution3 = ath.ReactionPlaneResolution3();



  }

    /// Normalise histograms etc., after the run
    void finalize() {
      //double EPRes = 1;

      //0-5% centrality
      add(_histInPlane1, _histOutPlane1, _histadd1);
      subtract(_histInPlane1, _histOutPlane1, _histsub1);

      scale(_histsub1, (pi/4)*(1/_resolution2));
      divide(_histsub1, _histadd1, _histV2ChJet1);

      //30-50% centrality
      add(_histInPlane2, _histOutPlane2, _histadd2);
      subtract(_histInPlane2, _histOutPlane2, _histsub2);

      scale(_histsub2, (pi/4)*(1/_resolution3));
      divide(_histsub2, _histadd2, _histV2ChJet2);

      cout << "finalize\n";
    }


    void add(Histo1DPtr h1, Histo1DPtr h2, Histo1DPtr hresult) const {
      const string path = hresult->path();
      *hresult = *h1 + *h2;
      hresult->setPath(path);
    }

    void subtract(Histo1DPtr h1, Histo1DPtr h2, Histo1DPtr hresult) const {
      const string path = hresult->path();
      *hresult = *h1 - *h2;
      hresult->setPath(path);
    }
    void divide(Histo1DPtr h1, Histo1DPtr h2, Scatter2DPtr hresult) const {
      const string path = hresult->path();
      *hresult = *h1 / *h2;
      hresult->setPath(path);
    }
    void integrate(Histo1DPtr h, Scatter2DPtr s) const {
      const string path = s->path();
      *s = toIntegralHisto(*h);
      s->setPath(path);
    }
    double deltaPhi2(double phi1, double phi2) {
      return mapAngle0To2Pi(phi1-phi2);
    }

    //@}

    /// @name Histogramst
    //@{

    Scatter2DPtr   _histV2ChJet1;   //0-5% centrality
    Scatter2DPtr   _histV2ChJet2;   //30-50% centrality

    Histo1DPtr     _histInPlane1;     //0-5% centrality
    Histo1DPtr     _histOutPlane1;    //0-5% centrality

    Histo1DPtr     _histInPlane2;     //30-50% centrality
    Histo1DPtr     _histOutPlane2;    //30-50% centrality

    Histo1DPtr     _histadd1;         //0-5% centrality
    Histo1DPtr     _histsub1;         //0-5% centrality

    Histo1DPtr     _histadd2;         //30-50% centrality
    Histo1DPtr     _histsub2;         //30-50% centrality

    double _resolution2;
    double _resolution3;
  //  double _Nevents0 = 0;
  //  double _Nevents30 = 0;

    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ALICE_2015_I1394678);


}
