// -*- C++ -*-
//
// Package:    JPsiKKKPATRunII
// Class:      JPsiKKKPATRunII
// 
/**\class JPsiKKKPATRunII JPsiKKKPATRunII.cc myAnalyzers/JPsiKKKPATRunII/src/JPsiKKKPATRunII.cc

 Description: <one line class summary>
Make rootTuple for JPsiKKK reconstruction

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  
//
//

#ifndef _JPsiKKKPATRunII_h
#define _JPsiKKKPATRunII_h

// system include files
#include <memory>

// user include files
#include "../interface/VertexReProducer.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/Common/interface/Handle.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/TransientTrackKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackFromFTSFactory.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/V0Candidate/interface/V0Candidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"

#include "CondFormats/L1TObjects/interface/L1GtTriggerMenu.h"
#include "CondFormats/DataRecord/interface/L1GtTriggerMenuRcd.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetupFwd.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerObjectMapRecord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetup.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/TrackReco/interface/DeDxData.h"

#include "TFile.h"
#include "TTree.h"

#include <vector>
#include <string>
//
// class decleration
//

using std::vector;
using namespace edm;
using namespace reco;
using namespace std;

class JPsiKKKPATRunII : public EDAnalyzer {
public:
  explicit JPsiKKKPATRunII(const ParameterSet&);
  ~JPsiKKKPATRunII();
  
private:
  virtual void beginJob() ;
  virtual void beginRun(Run const & iRun, EventSetup const& iSetup);
  virtual void analyze(const Event&, const EventSetup&);
  virtual void endJob() ;

/*	
  const reco::DeDxDataValueMap * energyLoss;
  Int_t iexception_dedx;

  virtual double getSigmaOfLogdEdx(double logde);
  virtual float getEnergyLoss(const reco::TrackRef & track);
  virtual double nsigmaofdedx(const reco::TrackRef & track, double & theo, double & sigma);
  virtual double getLogdEdx(double bg);
  virtual double GetMass(const reco::TrackRef & track);
  */

  // ----------member data ---------------------------
  string proccessName_;
  HLTConfigProvider hltConfig_;
  edm::EDGetTokenT<edm::TriggerResults> hlTriggerResults_;
  edm::EDGetTokenT<reco::BeamSpot> thebeamspot_;
  edm::EDGetTokenT<reco::VertexCollection> vtxSample;
  edm::EDGetTokenT<vector < pat::GenericParticle > > tracks_;
  edm::EDGetTokenT<L1GlobalTriggerReadoutRecord> L1Token;
  edm::EDGetTokenT< vector < pat::Muon >> muons_;
  edm::EDGetTokenT<reco::TrackCollection> revtxtrks_;
//  edm::EDGetTokenT<reco::GenParticleCollection>genParticles_;
  InputTag inputGEN_;
  bool doMC_;
  //bool doMC;
  int MCParticle;
  bool doJPsiMassCost;
  vector<std::string> qual;
  int MuPixHits_c;
  int MuSiHits_c;
  double MuNormChi_c;
  double MuD0_c;

  double JMaxM_c;
  double JMinM_c;
  int PiSiHits_c;
  double PiPt_c;
  double JPiPiDR_c;
  double XPiPiDR_c;
  bool UseXDr_c;	
  double JPiPiMax_c;
  double JPiPiMin_c;
	
  bool resolveAmbiguity_; 
  bool addXlessPrimaryVertex_;
  vector<string>      TriggersForMatching_;
  vector<string>      FiltersForMatching_;
  int  MatchingTriggerResult[50];
  bool Debug_;
  double Chi_Track_;

  vector<reco::TrackBase::TrackQuality> qualities;


  TTree* X_One_Tree_;
  TTree* X_One_TreeMC_;

  unsigned int        runNum, evtNum, lumiNum;
  
  vector<unsigned int>* trigRes;
  vector<std::string>* trigNames;
  vector<unsigned int>* L1TT;
  vector<std::string>* MatchTriggerNames;

  unsigned int        nX, nK, nJPsi, MCnX, nJ, nMu, nMC;
  float               priVtxX, priVtxY, priVtxZ, priVtxXE, priVtxYE, priVtxZE, priVtxChiNorm, priVtxChi, priVtxCL;
//  vector<float>       *MyPriVtxX, *MyPriVtxY, *MyPriVtxZ, *MyPriVtxXE, *MyPriVtxYE, *MyPriVtxZE, *MyPriVtxChiNorm, *MyPriVtxChi, *MyPriVtxCL;

  vector<float>       *xMass, *xVtxCL, *xVtxC2, *PsiTwoSMass, *SecondPairMass;
  vector<float>       *BsMass,*BsVtxCL, *BsVtxC2;
  vector<double>      *xPx, *xPy, *xPz, *xPxE, *xPyE, *xPzE;
  vector<float>       *xDecayVtxX, *xDecayVtxY, *xDecayVtxZ;
  vector<double>      *xDecayVtxXE, *xDecayVtxYE, *xDecayVtxZE;
  vector<float>       *MyPriVtxX, *MyPriVtxY, *MyPriVtxZ, *MyPriVtxXE, *MyPriVtxYE, *MyPriVtxZE, *MyPriVtxChiNorm, *MyPriVtxChi, *MyPriVtxCL;
  vector<float>       *PVtrackptXBLess, *PVtrackptYBLess, *PVtrackptZBLess, *PVCosThetaXBLess, *PVCosThetaYBLess, *PVCosThetaZBLess, *PVtrackptX, *PVtrackptY, *PVtrackptZ, *PVCosThetaX, *PVCosThetaY, *PVCosThetaZ,*PriVtxXCorrX, *PriVtxXCorrY, *PriVtxXCorrZ;
  vector<double>      *PriVtxXCorrEX, *PriVtxXCorrEY, *PriVtxXCorrEZ;
  vector<float>	      *PriVtxXCorrC2, *PriVtxXCorrCL;

  vector<double>      *xLxyPV, *xLxyPVE, *xCosAlpha, *xCTauPV, *xLxyBS, *xLxyBSE, *xExtraTrkabsDxywrtXvtxmin; 
//  vector<double>      *xCTauPVE, *xCTauBSE, *xExtraTrkabsDxywrtXvtxmin;
  vector<double>      *xExtraTrkabsDxywrtXvtxminMaxDz0dot5cm,*xExtraTrkabsDxywrtXvtxminMaxDz1dot0cm, *xExtraTrkabsDxywrtXvtxminMaxDz1dot5cm, *xExtraTrkabsDxywrtXvtxminMaxDz2dot0cm,*xExtraTrkabsDxywrtXvtxminMaxDz3dot0cm, *xExtraTrkabsDxywrtXvtxminMaxDz4dot0cm,*xExtraTrkabsDxywrtXvtxminMaxDz5dot0cm;
  vector<double>      *xExtraTrkabsDxywrtXvtxminMaxDz0dot1cm, *xExtraTrkabsDxywrtXvtxminMaxDz0dot05cm, *xExtraTrkabsDxywrtXvtxminMaxDz0dot04cm,*xExtraTrkabsDxywrtXvtxminMaxDz0dot03cm,*xExtraTrkabsDxywrtXvtxminMaxDz0dot02cm,*xExtraTrkabsDxywrtXvtxminMaxDz0dot01cm, *xCosAlphaBS, *xCTauBS, *xCTauPVE, *xCTauBSE;
  vector<int>         *JIndex;
  vector<int>         *pipIdx, *pimIdx, *pi3rdIdx;  //added pi3rdIdx by yik

  vector<float>       *JMass, *JVtxCL, *JVtxC2, *JPx, *JPy, *JPz;
  vector<float>       *JDecayVtxX, *JDecayVtxY, *JDecayVtxZ;
  vector<float>       *JDecayVtxXE, *JDecayVtxYE, *JDecayVtxZE;
  vector<int>         *mupIdx, *mumIdx;
  vector<bool>        *mupTrigMatch, *mumTrigMatch;
  vector<bool>        *JPsiMuonTrigMatch;

  vector<float>       *XEtawoFit, *XMasswoFit, *XPtwoFit, *mumPx, *mumPy, *mumPz, *mumfChi2;
  vector<float>       *mupPx, *mupPy, *mupPz, *mupfChi2;
  vector<int>         *mumfNDF, *mupfNDF;
  vector<double>      *muPx, *muPy, *muPz, *muD0;
  vector<float>       *cosAlpha, *mumuVtxCL, *mumuFLSig, *mumurVtxMag2D, *mumusigmaRvtxMag2D, *muD0E, *muDz, *muChi2, *muGlChi2, *mufHits;
  vector<bool>        *muFirstBarrel, *muFirstEndCap;
  vector<float>       *muDzVtx, *muDxyVtx;
  vector<int>         *muNDF, *muGlNDF, *muPhits, *muShits, *muGlMuHits, *muType;
  vector<float>       *muQual;
 vector<float>        *muTrack, *muCharge;
//  vector<float>       *MC_muPx, *MC_muPy, *MC_muPz, *MC_muCharge, *MC_muPdgId;
 // vector<int>	      *MC_muMother1, *MC_muMother2;
 

  vector<float>       *MCPx, *MCPy, *MCPz, *MCE, *MCJidx, *MCPhidx, *MCKidx, *MCBIdx, *MCStatus, *BIndex, *MCPt;
  vector<int>         *trQuality, *ThreeGoodTracks;
  vector<float>       *trPx, *trPy, *trPz, *trE;
 // vector<int>         *MCPdgId, *MCPdgIdAll, *MCnDecay, *MCNDaughters, *MCParent, *MCStatus, *MCDauJpsi_N, *MCDauK_N, *MCDauphi_N;
  vector<int>         *trNDF, *trPhits, *trShits;
  vector<float>     *trChi2;
  vector<float>     *trfHits;
  vector<bool>      *trFirstBarrel, *trFirstEndCap;
  vector<float>     *trDzVtx, *trDxyVtx, *trVx, *trVy;
  vector<double>     *tr_nsigdedx;
  vector<float>      *tr_dedx, *tr_dedxMass, *tr_theo, *tr_sigma;
	
  vector<float>       *MC_trPx, *MC_trPy, *MC_trPz,*MC_trphi, *MC_treta, *MC_trPdgId ,*MC_trE;
  vector<int>	      *MC_trMother1, *MC_trMother2;
  vector<float>       *MCmass, *MCPdgId, *trD0, *trD0E, *trCharge;
  vector<float>       *MCDauphi_Pt, *phiind, *psiind, *kaonind, *MCDauphi_Px, *MCDauphi_Py, *MCDauphi_Pz, *MCDauphi_mass, *MCDauphi_E, *MCmupDau_M, *MCmupDau_px, *MCmupDau_py, *MCmupDau_pz, *MCmupDau_E, *MCmumDau_M, *MCmumDau_px, *MCmumDau_py, *MCmumDau_pz, *MCmumDau_E, *MCkpDau_M, *MCkpDau_px, *MCkpDau_py, *MCkpDau_pz, *MCkpDau_E, *MCkmDau_M, *MCkmDau_px, *MCkmDau_py, *MCkmDau_pz, *MCkmDau_E, *BDauIdx, *MCDauJpsi_Pt, *MCDauJpsi_Px, *MCDauJpsi_Py, *MCDauJpsi_Pz, *MCDauJpsi_mass, *MCDauJpsi_E, *MCDauK_Pt, *MCDauK_Px, *MCDauK_Py, *MCDauK_Pz, *MCDauK_mass, *MCDauK_E;

  vector<float>       *fpi1Px, *fpi1Py, *fpi1Pz, *fpi1E;
  vector<float>       *fpi2Px, *fpi2Py, *fpi2Pz, *fpi2E;
  vector<float>       *fpi3Px, *fpi3Py, *fpi3Pz, *fpi3E; //added by yik
  vector<float>       *fmu1Px, *fmu1Py, *fmu1Pz, *fmu1E;
  vector<float>       *fmu2Px, *fmu2Py, *fmu2Pz, *fmu2E;
};

#endif
