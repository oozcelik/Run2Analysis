#include <memory>
#include "../interface/JPsiKKKPATRunII.h"
#include "../interface/VertexReProducer.h"
#include "TLorentzVector.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TrackingTools/PatternTools/interface/TwoTrackMinimumDistance.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GtFdlWord.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/Common/interface/RefToBase.h"
#include "DataFormats/Candidate/interface/ShallowCloneCandidate.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h"
#include "DataFormats/Math/interface/Error.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/Math/interface/Vector3D.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/TrackReco/interface/DeDxData.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/DeepCopyPointer.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/CLHEP/interface/AlgebraicObjects.h"
#include "DataFormats/CLHEP/interface/Migration.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"
#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"

#include "MagneticField/Engine/interface/MagneticField.h"

#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "RecoMuon/Records/interface/MuonRecoGeometryRecord.h"
#include "RecoMuon/DetLayers/interface/MuonDetLayerGeometry.h"
#include "Geometry/Records/interface/GlobalTrackingGeometryRecord.h"
#include "Geometry/CSCGeometry/interface/CSCGeometry.h"

#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"

#include <vector>
#include <utility>

#include "CLHEP/Matrix/Vector.h"
#include "CLHEP/Matrix/Matrix.h"
#include "CLHEP/Matrix/SymMatrix.h"
#include "CommonTools/UtilAlgos/interface/PhysObjectMatcher.h"
#include "CommonTools/UtilAlgos/interface/MCMatchSelector.h"
#include "CommonTools/UtilAlgos/interface/MatchByDRDPt.h"
#include "CommonTools/UtilAlgos/interface/MatchLessByDPt.h"

#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/JetReco/interface/GenJet.h"


// constants, enums and typedefs
typedef math::Error < 3 >::type CovarianceMatrix;
typedef ROOT::Math::SVector < double,3 > SVector3;
typedef ROOT::Math::SMatrix < double,3,3,ROOT::Math::MatRepSym < double,3 > >SMatrixSym3D;
typedef edm::Association<reco::GenParticleCollection> GenParticleMatch;
// constructors and destructor
JPsiKKKPATRunII::JPsiKKKPATRunII(const edm::ParameterSet & iConfig)
    :hlTriggerResults_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("HLTriggerResults"))),
    thebeamspot_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpotTag"))),  
     vtxSample(consumes< reco::VertexCollection > (iConfig.getParameter<edm::InputTag>("VtxSample"))),
     tracks_(consumes< vector < pat::GenericParticle > > (iConfig.getParameter<edm::InputTag>("Trak"))),
     muons_(consumes< vector <pat::Muon >> (iConfig.getParameter<edm::InputTag>("muons"))),
     revtxtrks_(consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("revtxtrks"))),
     doMC_(iConfig.getUntrackedParameter < bool > ("DoMonteCarloTree", false)),
     MCParticle(iConfig.getUntrackedParameter < int >("MonteCarloParticleId", 20443)),   // 20443 X, 100443 Psi(2S), 9120443
     doJPsiMassCost(iConfig.getUntrackedParameter < bool > ("DoJPsiMassConstraint", false)),
     qual(iConfig.getUntrackedParameter < std::vector < std::string > >("trackQualities")),
     MuPixHits_c(iConfig.getUntrackedParameter < int >("MinNumMuPixHits", 0)),
     MuSiHits_c(iConfig.getUntrackedParameter < int >("MinNumMuSiHits", 0)),
     MuNormChi_c(iConfig.getUntrackedParameter < double >("MaxMuNormChi2", 1000)),
     MuD0_c(iConfig.getUntrackedParameter < double >("MaxMuD0", 1000)),
     JMaxM_c(iConfig.getUntrackedParameter < double >("MaxJPsiMass", 4)),
     JMinM_c(iConfig.getUntrackedParameter < double >("MinJPsiMass", 2.2)),
     PiSiHits_c(iConfig.getUntrackedParameter < int >("MinNumTrSiHits", 0)),
     PiPt_c(iConfig.getUntrackedParameter < double >("MinTrPt", 0)),
     JPiPiDR_c(iConfig.getUntrackedParameter < double >("JPsiKKKMaxDR", 1)),
     XPiPiDR_c(iConfig.getUntrackedParameter < double >("XCandPiPiMaxDR", 1.1)),
     UseXDr_c(iConfig.getUntrackedParameter < bool > ("UseXDr", false)),
     JPiPiMax_c(iConfig.getUntrackedParameter < double >("JPsiKKKMaxMass", 50)),
     JPiPiMin_c(iConfig.getUntrackedParameter < double >("JPsiKKKMinMass", 0)),
     resolveAmbiguity_(iConfig.getUntrackedParameter < bool > ("resolvePileUpAmbiguity", true)),
     addXlessPrimaryVertex_(iConfig.getUntrackedParameter < bool > ("addXlessPrimaryVertex", true)),
     TriggersForMatching_(iConfig.getUntrackedParameter < std::vector < std::string > >("TriggersForMatching")),
     FiltersForMatching_(iConfig.getUntrackedParameter < std::vector < std::string > >("FiltersForMatching")),
     Debug_(iConfig.getUntrackedParameter < bool > ("Debug_Output", false)),
     Chi_Track_(iConfig.getUntrackedParameter < double >("Chi2NDF_Track", 10)),
     X_One_Tree_(0),
     runNum(0),
     evtNum(0),
     lumiNum(0),
     trigRes(0),
     trigNames(0),
     L1TT(0),
     MatchTriggerNames(0),
     nX(0),
     nK(0),
     nJPsi(0),
     MCnX(0),
     nJ(0),
     nMu(0),
     nMC(0),
     priVtxX(0),
     priVtxY(0),
     priVtxZ(0),
     priVtxXE(0),
     priVtxYE(0),
     priVtxZE(0),
     priVtxChiNorm(0),
     priVtxChi(0),
     priVtxCL(0),
     xMass(0),
     xVtxCL(0),
     xVtxC2(0),
     PsiTwoSMass(0),
     SecondPairMass(0),
     BsMass(0),
     BsVtxCL(0),
     BsVtxC2(0),
     xPx(0),
     xPy(0),
     xPz(0),
     xPxE(0),
     xPyE(0),
     xPzE(0),
     xDecayVtxX(0),
     xDecayVtxY(0),
     xDecayVtxZ(0),
     xDecayVtxXE(0),
     xDecayVtxYE(0),
     xDecayVtxZE(0),
     MyPriVtxX(0),
     MyPriVtxY(0),
     MyPriVtxZ(0),
     MyPriVtxXE(0),
     MyPriVtxYE(0),
     MyPriVtxZE(0),
     MyPriVtxChiNorm(0),
     MyPriVtxChi(0),
     MyPriVtxCL(0),
     PriVtxXCorrX(0),
     PriVtxXCorrY(0),
     PriVtxXCorrZ(0),
     PriVtxXCorrEX(0),
     PriVtxXCorrEY(0),
     PriVtxXCorrEZ(0),
     PriVtxXCorrC2(0),
     PriVtxXCorrCL(0),
     xLxyPV(0),
     xLxyPVE(0),
     xCosAlpha(0),
     xCTauPV(0),
     xLxyBS(0),
     xLxyBSE(0),
     xExtraTrkabsDxywrtXvtxmin(0),
     xExtraTrkabsDxywrtXvtxminMaxDz0dot5cm(0),xExtraTrkabsDxywrtXvtxminMaxDz1dot0cm(0), xExtraTrkabsDxywrtXvtxminMaxDz1dot5cm(0), xExtraTrkabsDxywrtXvtxminMaxDz2dot0cm(0),xExtraTrkabsDxywrtXvtxminMaxDz3dot0cm(0), xExtraTrkabsDxywrtXvtxminMaxDz4dot0cm(0),xExtraTrkabsDxywrtXvtxminMaxDz5dot0cm(0),
     xExtraTrkabsDxywrtXvtxminMaxDz0dot1cm(0), xExtraTrkabsDxywrtXvtxminMaxDz0dot05cm(0), xExtraTrkabsDxywrtXvtxminMaxDz0dot04cm(0),xExtraTrkabsDxywrtXvtxminMaxDz0dot03cm(0),xExtraTrkabsDxywrtXvtxminMaxDz0dot02cm(0),xExtraTrkabsDxywrtXvtxminMaxDz0dot01cm(0),
     xCosAlphaBS(0),
     xCTauBS(0),
     xCTauPVE(0),
     xCTauBSE(0),
     JIndex(0),
     pipIdx(0),
     pimIdx(0),
     pi3rdIdx(0),
     JMass(0),
     JVtxCL(0),
     JVtxC2(0),
     JPx(0),
     JPy(0),
     JPz(0),
     JDecayVtxX(0),
     JDecayVtxY(0),
     JDecayVtxZ(0),
     JDecayVtxXE(0),
     JDecayVtxYE(0),
     JDecayVtxZE(0),
     mupIdx(0),
     mumIdx(0),
     mupTrigMatch(0),
     mumTrigMatch(0),
     JPsiMuonTrigMatch(0),
     mumPx(0),
     mumPy(0),
     mumPz(0),
     mumfChi2(0),
     mupPx(0),
     mupPy(0),
     mupPz(0),
     mupfChi2(0),
     mumfNDF(0),
     mupfNDF(0),
     muPx(0),
     muPy(0),
     muPz(0),
     muD0(0),
     cosAlpha(0),
     mumuVtxCL(0),
     mumuFLSig(0),
     mumurVtxMag2D(0),
     mumusigmaRvtxMag2D(0),
     muD0E(0),
     muDz(0),
     muChi2(0),
     muGlChi2(0),
     mufHits(0),
     muFirstBarrel(0),
     muFirstEndCap(0),
     muDzVtx(0),
     muDxyVtx(0),
     muNDF(0),
     muGlNDF(0),
     muPhits(0),
     muShits(0),
     muGlMuHits(0),
     muType(0),
     muQual(0),
     muTrack(0),
     muCharge(0),
     MCPx(0),
     MCPy(0),
     MCPz(0),
     MCE(0),
MCJidx(0),
MCPhidx(0),
MCKidx(0),
     MCBIdx(0),
MCStatus(0),
     MCPt(0),
     trQuality(0),
     ThreeGoodTracks(0),
     trPx(0),
     trPy(0),
     trPz(0),
     trE(0),
     trNDF(0),
     trPhits(0),
     trShits(0),
     trChi2(0),
     trfHits(0),
     trFirstBarrel(0),
     trFirstEndCap(0),
     trDzVtx(0),
     trDxyVtx(0),
     trVx(0),
     trVy(0),
     tr_nsigdedx(0),
     tr_dedx(0),
     tr_dedxMass(0),
     tr_theo(0),
     tr_sigma(0),
     MC_trPx(0),
     MC_trPy(0),
     MC_trPz(0),
     MC_trphi(0),
     MC_treta(0),
     MC_trPdgId(0),
     MC_trE(0),
     MC_trMother1(0),
     MC_trMother2(0),
     MCmass(0),
     MCPdgId(0),
     trD0(0),
     trD0E(0),
     trCharge(0),
     MCDauphi_Pt(0),
     MCDauphi_Px(0),
     MCDauphi_Py(0),
     MCDauphi_Pz(0),
     MCDauphi_mass(0),
     MCDauphi_E(0),
     MCmupDau_M(0),
     MCmupDau_px(0),
     MCmupDau_py(0),
     MCmupDau_pz(0),
     MCmupDau_E(0),
     MCmumDau_M(0),
     MCmumDau_px(0),
     MCmumDau_py(0),
     MCmumDau_pz(0),
     MCmumDau_E(0),
     MCkpDau_M(0),
     MCkpDau_px(0),
     MCkpDau_py(0),
     MCkpDau_pz(0),
     MCkpDau_E(0),
     MCkmDau_M(0),
     MCkmDau_px(0),
     MCkmDau_py(0),
     MCkmDau_pz(0),
     MCkmDau_E(0),
     BDauIdx(0),
     MCDauJpsi_Pt(0),
     MCDauJpsi_Px(0),
     MCDauJpsi_Py(0),
     MCDauJpsi_Pz(0),
     MCDauJpsi_mass(0),
     MCDauJpsi_E(0),
     MCDauK_Pt(0),
     MCDauK_Px(0),
     MCDauK_Py(0),
     MCDauK_Pz(0),
     MCDauK_mass(0),
     MCDauK_E(0),
     fpi1Px(0),
     fpi1Py(0),
     fpi1Pz(0),
     fpi1E(0),
     fpi2Px(0),
     fpi2Py(0),
     fpi2Pz(0),
     fpi2E(0),
     fpi3Px(0),
     fpi3Py(0),
     fpi3Pz(0),
     fpi3E(0),
     fmu1Px(0),
     fmu1Py(0),
     fmu1Pz(0),
     fmu1E(0),
     fmu2Px(0),
     fmu2Py(0),
     fmu2Pz(0),
     fmu2E(0)
{
  L1Token = consumes<L1GlobalTriggerReadoutRecord>((edm::InputTag)"gtDigis");
}

JPsiKKKPATRunII::~JPsiKKKPATRunII()
{
    // do anything here that needs to be done at desctruction time
    // (e.g. close files, deallocate resources etc.)
}

// ------------ method called to for each event  ------------
void JPsiKKKPATRunII::analyze(const edm::Event & iEvent, const edm::EventSetup & iSetup)
{
    // get event content information

    using std::vector;
    using namespace edm;
    using namespace reco;
    using namespace std;

    runNum = iEvent.id().run();
    evtNum = iEvent.id().event();
    lumiNum = iEvent.id().luminosityBlock();

    ESHandle < MagneticField > bFieldHandle;
    iSetup.get < IdealMagneticFieldRecord > ().get(bFieldHandle);


    // first get HLT results
    edm::Handle < edm::TriggerResults > hltresults;
    try
    {
        iEvent.getByToken(hlTriggerResults_, hltresults);
    }
    catch(...)
    {
        cout << "Couldn't get handle on HLT Trigger!" << endl;
    }
    if (!hltresults.isValid())
    {
        cout << "No Trigger Results!" << endl;
    }
    else
    {
        int ntrigs = hltresults->size();
        if (ntrigs == 0)
        {
            cout << "No trigger name given in TriggerResults of the input " << endl;
        }

        // get hold of trigger names - based on TriggerResults object!
        edm::TriggerNames triggerNames_;
        triggerNames_ = iEvent.triggerNames(*hltresults);

        int ntriggers = TriggersForMatching_.size();
        for (int MatchTrig = 0; MatchTrig < ntriggers; MatchTrig++)
        {
            MatchingTriggerResult[MatchTrig] = 0;
            
        }

        for (int itrig = 0; itrig < ntrigs; itrig++)
        {
            string trigName = triggerNames_.triggerName(itrig);
   
            // check if trigger fired //  
            int hltflag = (*hltresults)[itrig].accept();
            ////// hltflag 0 or 1 ///////////
            trigRes->push_back(hltflag);
            trigNames->push_back(trigName);
    
            /// loop for selected triggers ///
            int ntriggers = TriggersForMatching_.size();
            for (int MatchTrig = 0; MatchTrig < ntriggers; MatchTrig++)
            {
                /// check if hlt triggers match with the selected ones. ////
                if (TriggersForMatching_[MatchTrig] == triggerNames_.triggerName(itrig))
                {
                    /// will leave the loop once it has matched ///
                    MatchingTriggerResult[MatchTrig] = hltflag;
                    break;
               }
             } 
           }
           for (int MatchTrig = 0; MatchTrig < ntriggers; MatchTrig++)
            {
              MatchTriggerNames->push_back(TriggersForMatching_[MatchTrig]);
            }
        }
    // get L1 trigger info
    edm::ESHandle < L1GtTriggerMenu > menuRcd;
    iSetup.get < L1GtTriggerMenuRcd > ().get(menuRcd);

    edm::Handle < L1GlobalTriggerReadoutRecord > gtRecord;
    iEvent.getByToken(L1Token, gtRecord);
    const DecisionWord dWord = gtRecord->decisionWord();

    const TechnicalTriggerWord ttWord = gtRecord->technicalTriggerWord();
    for (unsigned int l1i = 0; l1i != ttWord.size(); ++l1i)
    {
        L1TT->push_back(ttWord.at(l1i));
    }
    //////////////////////////// 
    
    Vertex thePrimaryV;
    Vertex theRecoVtx;
    Vertex theBeamSpotV;
    BeamSpot beamSpot;
    math::XYZPoint RefVtx;

    // get BeamSpot
    edm::Handle<reco::BeamSpot> beamSpotHandle;
    iEvent.getByToken(thebeamspot_, beamSpotHandle);

    if (beamSpotHandle.isValid())
    {
        beamSpot = *beamSpotHandle;
        theBeamSpotV = Vertex(beamSpot.position(), beamSpot.covariance3D());
    }
    else
        cout << "No beam spot available from EventSetup" << endl;


    Handle < VertexCollection > recVtxs;
    iEvent.getByToken(vtxSample, recVtxs);
    unsigned int nVtxTrks = 0;

    if (recVtxs->begin() != recVtxs->end())
    {
        if (addXlessPrimaryVertex_ || resolveAmbiguity_)
        {
            thePrimaryV = Vertex(*(recVtxs->begin()));
        }
        else
        {
            for (reco::VertexCollection::const_iterator vtx = recVtxs->begin(); vtx != recVtxs->end(); ++vtx)
            {
                if (nVtxTrks < vtx->tracksSize())
                {
                    nVtxTrks = vtx->tracksSize();
                    thePrimaryV = Vertex(*vtx);
                }
            }
        }
    }
    else
    {
        thePrimaryV = Vertex(beamSpot.position(), beamSpot.covariance3D());
    }

    RefVtx = thePrimaryV.position();
    priVtxX = (thePrimaryV.position().x());
    priVtxY = (thePrimaryV.position().y());
    priVtxZ = (thePrimaryV.position().z());
    priVtxXE = (thePrimaryV.xError());
    priVtxYE = (thePrimaryV.yError());
    priVtxZE = (thePrimaryV.zError());
    priVtxChiNorm = (thePrimaryV.normalizedChi2());
    priVtxChi = thePrimaryV.chi2();
    priVtxCL = ChiSquaredProbability((double) (thePrimaryV.chi2()), (double) (thePrimaryV.ndof()));

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // try reconstruction without fitting modules
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    Handle < vector < pat::GenericParticle > >thePATTrackHandle;
    iEvent.getByToken(tracks_, thePATTrackHandle);
    Handle <  vector <pat::Muon  >>thePATMuonHandle;
    iEvent.getByToken(muons_, thePATMuonHandle);

    if (Debug_)
    {
//        cout << "starting event with " << thePATTrackHandle->size() << " tracks, and " << thePATMuonHandle->size() << " muons" << endl;
    }

    for (unsigned int ndx = 0; ndx < qual.size(); ndx++)
    {
        qualities.push_back(reco::TrackBase::qualityByName(qual[ndx]));
    }

    if (thePATTrackHandle->size() >= 2 && thePATMuonHandle->size() >= 2)    // at least 2 tracks and two muons
    {

        for (vector < pat::GenericParticle >::const_iterator iTr = thePATTrackHandle->begin();
                iTr != thePATTrackHandle->end(); ++iTr)
        {

            pat::GenericParticle tr = *iTr;
            TrackRef tmpRef = tr.track();

            if (tmpRef->quality(qualities[2]))
            {
                trQuality->push_back(2);
            }
            else if (tmpRef->quality(qualities[1]))
            {
                trQuality->push_back(1);
            }
            else if (tmpRef->quality(qualities[0]))
            {
                trQuality->push_back(0);
            }

            trPx->push_back(tr.px());
            trPy->push_back(tr.py());
            trPz->push_back(tr.pz());
            trE->push_back(tr.energy());
            trNDF->push_back(tr.track()->ndof());
            trPhits->push_back(tr.track()->hitPattern().numberOfValidPixelHits());
            trShits->push_back(tr.track()->hitPattern().numberOfValidStripHits());
            trChi2->push_back(tr.track()->chi2());
            trD0->push_back(tr.track()->d0());
            trD0E->push_back(tr.track()->d0Error());
            trCharge->push_back(tr.charge());

/*            float hits =
                (1.0 * tr.track()->found()) / (tr.track()->found() + tr.track()->lost() +
                                               tr.track()->trackerExpectedHitsInner().numberOfHits() +
                                               tr.track()->trackerExpectedHitsOuter().numberOfHits());*/
//            trfHits->push_back(hits);
  /*          trFirstBarrel->push_back(tr.track()->hitPattern().hasValidHitInFirstPixelBarrel());
            trFirstEndCap->push_back(tr.track()->hitPattern().hasValidHitInFirstPixelEndcap());*/
            trDzVtx->push_back(tr.track()->dz(RefVtx));
            trDxyVtx->push_back(tr.track()->dxy(RefVtx));
            trVx->push_back(tr.track()->vx());
            trVy->push_back(tr.track()->vy());

            tr_nsigdedx->push_back(-9999.0);
            tr_dedx->push_back(-9999.0);
            tr_dedxMass->push_back(-9999.0);
            tr_theo->push_back(-9999.0);
            tr_sigma->push_back(-9999.0);

        }                      // end of foor loop for iTr = thePATTrackHandle->begin()


        // get X and J cands
        for (vector <pat::Muon >::const_iterator iMuonP = thePATMuonHandle->begin(); iMuonP != thePATMuonHandle->end(); ++iMuonP)
        {
            ++nMu;
            const reco::Muon * rmu = dynamic_cast < const reco::Muon * >(iMuonP->originalObject()); // first muon

            muPx->push_back(rmu->px());
            muPy->push_back(rmu->py());
            muPz->push_back(rmu->pz());
            muCharge->push_back(rmu->charge());


            if (rmu->track().isNull())      // rmu->track() returns innerTrack();
            {
                muD0->push_back(0);
                muD0E->push_back(0);
                muDz->push_back(0);
                muChi2->push_back(0);
                muNDF->push_back(-1);
                muPhits->push_back(0);
                muShits->push_back(0);

                muDzVtx->push_back(0);
                muDxyVtx->push_back(0);
                mufHits->push_back(0);
                muFirstBarrel->push_back(0);
                muFirstEndCap->push_back(0);
            }
            else
            {
                muD0->push_back(rmu->track()->d0());
                muD0E->push_back(rmu->track()->d0Error());
                muDz->push_back(rmu->track()->dz());
                muChi2->push_back(rmu->track()->chi2());
                muNDF->push_back(rmu->track()->ndof());
                muPhits->push_back(rmu->track()->hitPattern().numberOfValidPixelHits());
                muShits->push_back(rmu->track()->hitPattern().numberOfValidStripHits());
                muDzVtx->push_back(rmu->track()->dz(RefVtx));
                muDxyVtx->push_back(rmu->track()->dxy(RefVtx));
             /*   mufHits->push_back((1.0 * rmu->track()->found()) / (rmu->track()->found() + rmu->track()->lost() +
                                   rmu->track()->trackerExpectedHitsInner().numberOfHits() +
                                   rmu->track()->trackerExpectedHitsOuter().numberOfHits()));*/
/*                muFirstBarrel->push_back(rmu->track()->hitPattern().hasValidHitInFirstPixelBarrel());
                muFirstEndCap->push_back(rmu->track()->hitPattern().hasValidHitInFirstPixelEndcap());*/
            }                   // filling muon information for all the muon tracks

            if (rmu->globalTrack().isNull())    // rmu->globalTrack() returns globalTrack();
            {
                muGlChi2->push_back(0);
                muGlNDF->push_back(-1);
                muGlMuHits->push_back(0);
            }
            else
            {
                muGlChi2->push_back(rmu->globalTrack()->chi2());
                muGlNDF->push_back(rmu->globalTrack()->ndof());
                muGlMuHits->push_back(rmu->globalTrack()->hitPattern().numberOfValidMuonHits());
            }                   // filling global track information

            muType->push_back(rmu->type());
            int qm = 0;
            for (int qi = 1; qi != 24; ++qi)
            {
                if (muon::isGoodMuon(*rmu, muon::SelectionType(qi)))
                {
                    qm += 1 << qi;
                }
            }
            // 11-> TMOneStationLoose 12-> TMOneStationTight
            muQual->push_back(qm);
            muTrack->push_back(-1); // not implemented yet

            TrackRef muTrackP = iMuonP->track();    // first track is accepted as plus charged
            if (muTrackP.isNull()) continue;

            // 2012 Soft Muon
            bool SoftMuon1= false;
            if (muon::isGoodMuon(*rmu,muon::TMOneStationTight)
                    && rmu->track()->hitPattern().trackerLayersWithMeasurement() > 5
                    && rmu->innerTrack()->hitPattern().pixelLayersWithMeasurement() > 1
                    && rmu->innerTrack()->normalizedChi2() < 1.8
                    && fabs(rmu->innerTrack()->dxy(RefVtx)) < 3 && fabs(rmu->innerTrack()->dz(RefVtx)) < 30
               ) SoftMuon1 = true;

            if (!SoftMuon1) continue;

            // next check for mu-
            for ( vector <pat::Muon >::const_iterator iMuonM = iMuonP + 1; iMuonM != thePATMuonHandle->end(); ++iMuonM)
            {

               if (iMuonM->charge() * iMuonP->charge() > 0) continue;   // second track should have different charge than the first one
                TrackRef muTrackM = iMuonM->track();    // second muon track from PATmuons
                if (muTrackM.isNull()) continue;


                const reco::Muon * rmu2 = dynamic_cast < const reco::Muon * >(iMuonM->originalObject());    // second
                if (muon::overlap(*rmu, *rmu2)) continue;   // no overlap betweed the 1st and 2nd reconstrcuted muons

                bool SoftMuon2= false;
                if (muon::isGoodMuon(*rmu2,muon::TMOneStationTight)
                        && rmu2->track()->hitPattern().trackerLayersWithMeasurement() > 5
                        && rmu2->innerTrack()->hitPattern().pixelLayersWithMeasurement() > 1
                        && rmu2->innerTrack()->normalizedChi2() < 1.8
                        && fabs(rmu2->innerTrack()->dxy(RefVtx)) < 3 && fabs(rmu2->innerTrack()->dz(RefVtx)) < 30
                   ) SoftMuon2 = true;

                if (!SoftMuon2) continue;

                // Get The J/Psi information from the tracks
                TransientTrack muonPTT(muTrackP, &(*bFieldHandle));
                TransientTrack muonMTT(muTrackM, &(*bFieldHandle));

                // The mass of a muon and the insignificant mass sigma to avoid singularities in the covariance matrix.
                ParticleMass muon_mass = 0.10565837;    // pdg mass
                float muon_sigma = muon_mass * 1.e-6;

                // initial chi2 and ndf before kinematic fits.
                float chi = 0.;
                float ndf = 0.;

                vector < RefCountedKinematicParticle > muonParticles;
                KinematicParticleFactoryFromTransientTrack pFactory;
                muonParticles.push_back(pFactory.particle(muonPTT, muon_mass, chi, ndf, muon_sigma));
                muonParticles.push_back(pFactory.particle(muonMTT, muon_mass, chi, ndf, muon_sigma));
                KinematicParticleVertexFitter fitter;
                RefCountedKinematicTree psiVertexFitTree;
                psiVertexFitTree = fitter.fit(muonParticles);   // fit to the muon pair
                if (!psiVertexFitTree->isValid()) continue;
                psiVertexFitTree->movePointerToTheTop();
                RefCountedKinematicParticle psi_vFit_noMC = psiVertexFitTree->currentParticle();
                RefCountedKinematicVertex psi_vFit_vertex_noMC = psiVertexFitTree->currentDecayVertex();

                // added by asli
                std::vector < double >mumuVtxEVec;
                mumuVtxEVec.push_back(psi_vFit_vertex_noMC->error().cxx());
                mumuVtxEVec.push_back(psi_vFit_vertex_noMC->error().cyx());
                mumuVtxEVec.push_back(psi_vFit_vertex_noMC->error().cyy());
                mumuVtxEVec.push_back(psi_vFit_vertex_noMC->error().czx());
                mumuVtxEVec.push_back(psi_vFit_vertex_noMC->error().czy());
                mumuVtxEVec.push_back(psi_vFit_vertex_noMC->error().czz());
                SMatrixSym3D mumuVtxCov(mumuVtxEVec.begin(), mumuVtxEVec.end());
                SMatrixSym3D totalCovMuMu = mumuVtxCov + beamSpot.covariance3D();
                SVector3 mumudistanceVector2DD(psi_vFit_vertex_noMC->position().x()-beamSpot.x0(),psi_vFit_vertex_noMC->position().y()-beamSpot.y0(),0.);
                double mumurVtxMag2DD = ROOT::Math::Mag(mumudistanceVector2DD);
                mumurVtxMag2D->push_back(mumurVtxMag2DD);
                double mumusigmaRvtxMag2DD = sqrt(ROOT::Math::Similarity(totalCovMuMu, mumudistanceVector2DD)) / mumurVtxMag2DD;
                mumusigmaRvtxMag2D->push_back(mumusigmaRvtxMag2DD);
                float mumuFLSigFloat = (mumurVtxMag2DD / mumusigmaRvtxMag2DD);
                mumuVtxCL->push_back(ChiSquaredProbability((double) (psi_vFit_vertex_noMC->chiSquared()), (double) (psi_vFit_vertex_noMC->degreesOfFreedom())));
                mumuFLSig->push_back(mumuFLSigFloat);
                // finish added by asli

                psiVertexFitTree->movePointerToTheFirstChild();
                RefCountedKinematicParticle muPCandMC = psiVertexFitTree->currentParticle();
                psiVertexFitTree->movePointerToTheNextChild();
                RefCountedKinematicParticle muMCandMC = psiVertexFitTree->currentParticle();
                KinematicParameters psiMupKP = muPCandMC->currentState().kinematicParameters();
                KinematicParameters psiMumKP = muMCandMC->currentState().kinematicParameters();

                // added by asli
                math::XYZVector pperp(psiMumKP.momentum().x() + psiMupKP.momentum().x(), psiMumKP.momentum().y() + psiMupKP.momentum().y(), 0.);
                // vertex direction
                GlobalPoint displacementFromBeamspot(-1 * ((beamSpot.x0() - psi_vFit_vertex_noMC->position().x()) +
                                                     (psi_vFit_vertex_noMC->position().z() - beamSpot.z0()) * beamSpot.dxdz()),
                                                     -1 * ((beamSpot.y0() - psi_vFit_vertex_noMC->position().y()) +
                                                             (psi_vFit_vertex_noMC->position().z() - beamSpot.z0()) * beamSpot.dydz()), 0);
                Vertex::Point vperp(displacementFromBeamspot.x(), displacementFromBeamspot.y(), 0.);
                float cosAlphaFloat = vperp.Dot(pperp) / (vperp.R() * pperp.R());
                cosAlpha->push_back(cosAlphaFloat);
                // finish added by asli

                // Fill the jtree vectors
                JMass->push_back(psi_vFit_noMC->currentState().mass());
                JDecayVtxX->push_back(psi_vFit_vertex_noMC->position().x());
                JDecayVtxY->push_back(psi_vFit_vertex_noMC->position().y());
                JDecayVtxZ->push_back(psi_vFit_vertex_noMC->position().z());
                JDecayVtxXE->push_back(sqrt(psi_vFit_vertex_noMC->error().cxx()));
                JDecayVtxYE->push_back(sqrt(psi_vFit_vertex_noMC->error().cyy()));
                JDecayVtxZE->push_back(sqrt(psi_vFit_vertex_noMC->error().czz()));
                JVtxCL->push_back(ChiSquaredProbability((double) (psi_vFit_vertex_noMC->chiSquared()),(double) (psi_vFit_vertex_noMC->degreesOfFreedom())));
                JVtxC2->push_back(psi_vFit_vertex_noMC->chiSquared());
                JPx->push_back(psiMumKP.momentum().x() + psiMupKP.momentum().x());
                JPy->push_back(psiMumKP.momentum().y() + psiMupKP.momentum().y());
                JPz->push_back(psiMumKP.momentum().z() + psiMupKP.momentum().z());
                mumPx->push_back(psiMumKP.momentum().x());
                mumPy->push_back(psiMumKP.momentum().y());
                mumPz->push_back(psiMumKP.momentum().z());
                mupPx->push_back(psiMupKP.momentum().x());
                mupPy->push_back(psiMupKP.momentum().y());
                mupPz->push_back(psiMupKP.momentum().z());
                mupIdx->push_back(std::distance(thePATMuonHandle->begin(), iMuonP));
                mumIdx->push_back(std::distance(thePATMuonHandle->begin(), iMuonM));
                mumfChi2->push_back(muMCandMC->chiSquared());
                mumfNDF->push_back(muMCandMC->degreesOfFreedom());
                mupfChi2->push_back(muPCandMC->chiSquared());
                mupfNDF->push_back(muPCandMC->degreesOfFreedom());

                int ntriggers = TriggersForMatching_.size();
                for (int MatchTrig = 0; MatchTrig < ntriggers; MatchTrig++)
                {
                    if (MatchingTriggerResult[MatchTrig] != 0)
                    {
                        const pat::TriggerObjectStandAloneCollection mu1HLTMatches = iMuonP->triggerObjectMatchesByFilter(FiltersForMatching_[MatchTrig]);
                        const pat::TriggerObjectStandAloneCollection mu2HLTMatches = iMuonM->triggerObjectMatchesByFilter(FiltersForMatching_[MatchTrig]);
                        bool pass1 = mu1HLTMatches.size() > 0;
                        bool pass2 = mu2HLTMatches.size() > 0;
                        mupTrigMatch->push_back(pass1);
                        mumTrigMatch->push_back(pass2);
                        if ((pass1) && (pass2)) JPsiMuonTrigMatch->push_back(true);
                        else JPsiMuonTrigMatch->push_back(false);
//               cout << " pass-1 " << pass1 << " pass-2 "<< pass2 << " " << MatchingTriggerResult[MatchTrig]  << endl;
                    }
                    else
                        JPsiMuonTrigMatch->push_back(false);
                }

                ++nJ;
                muonParticles.clear();

                double Num1 = psi_vFit_vertex_noMC->chiSquared();
                double Num2 = psi_vFit_vertex_noMC->degreesOfFreedom();
                if (ChiSquaredProbability(Num1, Num2) < 0.001) continue;
                if (JMass->at(nJ - 1) < JMinM_c || JMass->at(nJ - 1) > JMaxM_c) continue;


                float mytrackDZmax = 2.0;
                bool TrackP_Pure;
                bool TrackM_Pure;
                bool Track3rd_Pure;

                for (vector < pat::GenericParticle >::const_iterator iTrackP = thePATTrackHandle->begin(); iTrackP != thePATTrackHandle->end(); ++iTrackP)   // 1st
                {
                    pat::GenericParticle MyTrackP = *iTrackP;
                    TrackP_Pure = MyTrackP.track()->quality(qualities[2]);

                    if (!TrackP_Pure) continue;
                    if (iTrackP->track().key() == rmu->track().key() || iTrackP->track().key() == rmu2->track().key()) continue;
                    if ((iTrackP->track()->chi2() / iTrackP->track()->ndof() > Chi_Track_) || iTrackP->pt() < PiPt_c)  continue;
                    if (abs((*iTrackP).track()->dz((math::XYZPoint) psi_vFit_vertex_noMC->position())) > mytrackDZmax) continue;

                    for (vector < pat::GenericParticle >::const_iterator iTrackM = iTrackP + 1; iTrackM != thePATTrackHandle->end(); ++iTrackM)   // 2nd
                    {
                        pat::GenericParticle MyTrackM = *iTrackM;
                        TrackM_Pure = MyTrackM.track()->quality(qualities[2]);

                        if (!TrackM_Pure) continue;
                        if (iTrackM->track().key() == rmu->track().key() || iTrackM->track().key() == rmu2->track().key()) continue;
                        if ((iTrackM->track()->chi2() / iTrackM->track()->ndof() > Chi_Track_) || iTrackM->pt() < PiPt_c)  continue;
                        if (abs((*iTrackM).track()->dz((math::XYZPoint) psi_vFit_vertex_noMC->position())) > mytrackDZmax) continue;

                        for (vector < pat::GenericParticle >::const_iterator iTrack3rd = iTrackM + 1; iTrack3rd != thePATTrackHandle->end(); ++iTrack3rd)   // 3rd
                        {
                            pat::GenericParticle MyTrack3rd = *iTrack3rd;
                            Track3rd_Pure = MyTrack3rd.track()->quality(qualities[2]);
                            if(!Track3rd_Pure) continue;
                            if (abs((iTrackP->charge() + iTrackM->charge() + iTrack3rd->charge())) != 1) continue;
                            if (iTrack3rd->track().key() == rmu->track().key() || iTrack3rd->track().key() == rmu2->track().key())  continue;
                            if ((iTrack3rd->track()->chi2() / iTrack3rd->track()->ndof() > Chi_Track_) || iTrack3rd->pt() < PiPt_c) continue;
                            if (abs((*iTrack3rd).track()->dz((math::XYZPoint) psi_vFit_vertex_noMC->position())) > mytrackDZmax)    continue;

                            int ChargeSingleKaon = 0;
                            math::XYZTLorentzVector kpkmpairlowp4(0, 0, 0, 9999.0);
                            math::XYZTLorentzVector SingleKaonp4(0, 0, 0, 9999.0);
                            int FindThePair = -999;
                            if ((iTrackP->charge() + iTrackM->charge()) == 0 && (iTrackP->p4() + iTrackM->p4()).M() < kpkmpairlowp4.M())
                            {
                                kpkmpairlowp4 = iTrackP->p4() + iTrackM->p4();
                                SingleKaonp4 = iTrack3rd->p4();
                            }
                            if ((iTrackP->charge() + iTrack3rd->charge()) == 0 && (iTrackP->p4() + iTrack3rd->p4()).M() < kpkmpairlowp4.M())
                            {
                                kpkmpairlowp4 = iTrackP->p4() + iTrack3rd->p4();
                                SingleKaonp4 = iTrackM->p4();
                            }
                            if ((iTrack3rd->charge() + iTrackM->charge()) == 0 && (iTrack3rd->p4() + iTrackM->p4()).M() < kpkmpairlowp4.M())
                            {
                                kpkmpairlowp4 = iTrack3rd->p4() + iTrackM->p4();
                                SingleKaonp4 = iTrackP->p4();
                            }
                            if (kpkmpairlowp4.M() > 1.06) continue;

                            TransientTrack kaonPTT(iTrackP->track(), &(*bFieldHandle));
                            TransientTrack kaonMTT(iTrackM->track(), &(*bFieldHandle));
                            TransientTrack kaon3rdTT(iTrack3rd->track(), &(*bFieldHandle));
                            const ParticleMass kaon_mass = 0.493677;    // changed
                            float kaon_sigma = kaon_mass * 1.e-6;
                            // Do mass constraint for JPsi cand and do mass constrained vertex fit creating the
                            // constraint with a small sigma to put in the resulting covariance matrix in order to
                            // avoid singularitieskaon_mass

                            vector < RefCountedKinematicParticle > vFitMCParticles;
                            vFitMCParticles.push_back(pFactory.particle(muonPTT, muon_mass, chi, ndf, muon_sigma));
                            vFitMCParticles.push_back(pFactory.particle(muonMTT, muon_mass, chi, ndf, muon_sigma));
                            vFitMCParticles.push_back(pFactory.particle(kaonPTT, kaon_mass, chi, ndf, kaon_sigma));
                            vFitMCParticles.push_back(pFactory.particle(kaonMTT, kaon_mass, chi, ndf, kaon_sigma));
                            vFitMCParticles.push_back(pFactory.particle(kaon3rdTT, kaon_mass, chi, ndf, kaon_sigma));

                            RefCountedKinematicParticle xCandMC;
                            RefCountedKinematicVertex xDecayVertexMC;
                            RefCountedKinematicTree vertexFitTree;
                            if (doJPsiMassCost)     // cout
                            {

                                ParticleMass psi_mass = 3.096916;
                                MultiTrackKinematicConstraint *j_psi_c = new TwoTrackMassKinematicConstraint(psi_mass);
                                KinematicConstrainedVertexFitter kcvFitter;
                                vertexFitTree = kcvFitter.fit(vFitMCParticles, j_psi_c);
                                if (!vertexFitTree->isValid())  continue;
                                vertexFitTree->movePointerToTheTop();
                                xCandMC = vertexFitTree->currentParticle();
                                xDecayVertexMC = vertexFitTree->currentDecayVertex();
                            }
                            else        // fit w/o mass constraint
                            {
                                KinematicParticleVertexFitter kcvFitter;
                                vertexFitTree = kcvFitter.fit(vFitMCParticles);
                                if (!vertexFitTree->isValid()) continue;
                                vertexFitTree->movePointerToTheTop();
                                xCandMC = vertexFitTree->currentParticle();
                                xDecayVertexMC = vertexFitTree->currentDecayVertex();
                            }
                             if (!xDecayVertexMC->vertexIsValid()) continue;
                             if (xDecayVertexMC->chiSquared() < 0 || xDecayVertexMC->chiSquared() > 10000) continue;

                            // NEW == B vtx chi2 probability cut from Kai's code
                            double Num11 = xDecayVertexMC->chiSquared();
                            double Num22 = xDecayVertexMC->degreesOfFreedom();
                            double xVtx_CL = ChiSquaredProbability(Num11, Num22);
                            if (xVtx_CL < 0.000001) continue;

                            if (xCandMC->currentState().mass() > 5.6 || xCandMC->currentState().mass() < 5.0) continue;


			    math::XYZPoint myXpoint(xDecayVertexMC->position().x(),xDecayVertexMC->position().y(),xDecayVertexMC->position().z());
                            xMass->push_back(xCandMC->currentState().mass());
                            xPx->push_back(xCandMC->currentState().globalMomentum().x());
                            xPy->push_back(xCandMC->currentState().globalMomentum().y());
                            xPz->push_back(xCandMC->currentState().globalMomentum().z());
                            xPxE->push_back(sqrt(xCandMC->currentState().kinematicParametersError().matrix()(3, 3)));
                            xPyE->push_back(sqrt(xCandMC->currentState().kinematicParametersError().matrix()(4, 4)));
                            xPzE->push_back(sqrt(xCandMC->currentState().kinematicParametersError().matrix()(5, 5)));
                            xVtxCL->push_back(ChiSquaredProbability((double) (xDecayVertexMC->chiSquared()),(double) (xDecayVertexMC->degreesOfFreedom())));
                            xVtxC2->push_back(xDecayVertexMC->chiSquared());
                            xDecayVtxX->push_back((*xDecayVertexMC).position().x());
                            xDecayVtxY->push_back((*xDecayVertexMC).position().y());
                            xDecayVtxZ->push_back((*xDecayVertexMC).position().z());
                            xDecayVtxXE->push_back(sqrt((*xDecayVertexMC).error().cxx()));
                            xDecayVtxYE->push_back(sqrt((*xDecayVertexMC).error().cyy()));
                            xDecayVtxZE->push_back(sqrt((*xDecayVertexMC).error().czz()));
                            vertexFitTree->movePointerToTheFirstChild();
                            RefCountedKinematicParticle mu1 = vertexFitTree->currentParticle();
                            vertexFitTree->movePointerToTheNextChild();
                            RefCountedKinematicParticle mu2 = vertexFitTree->currentParticle();
                            vertexFitTree->movePointerToTheNextChild();
                            RefCountedKinematicParticle kaon1 = vertexFitTree->currentParticle();
                            vertexFitTree->movePointerToTheNextChild();
                            RefCountedKinematicParticle kaon2 = vertexFitTree->currentParticle();
                            vertexFitTree->movePointerToTheNextChild();
                            RefCountedKinematicParticle kaon3 = vertexFitTree->currentParticle();
                            // getting parameters from B fit
                            fpi1Px->push_back(kaon1->currentState().globalMomentum().x());
                            fpi1Py->push_back(kaon1->currentState().globalMomentum().y());
                            fpi1Pz->push_back(kaon1->currentState().globalMomentum().z());
                            fpi1E->push_back(kaon1->currentState().kinematicParameters().energy());
                            fpi2Px->push_back(kaon2->currentState().globalMomentum().x());
                            fpi2Py->push_back(kaon2->currentState().globalMomentum().y());
                            fpi2Pz->push_back(kaon2->currentState().globalMomentum().z());
                            fpi2E->push_back(kaon2->currentState().kinematicParameters().energy());
                            fpi3Px->push_back(kaon3->currentState().globalMomentum().x());
                            fpi3Py->push_back(kaon3->currentState().globalMomentum().y());
                            fpi3Pz->push_back(kaon3->currentState().globalMomentum().z());
                            fpi3E->push_back(kaon3->currentState().kinematicParameters().energy());
                            fmu1Px->push_back(mu1->currentState().globalMomentum().x());
                            fmu1Py->push_back(mu1->currentState().globalMomentum().y());
                            fmu1Pz->push_back(mu1->currentState().globalMomentum().z());
                            fmu1E->push_back(mu1->currentState().kinematicParameters().energy());
                            fmu2Px->push_back(mu2->currentState().globalMomentum().x());
                            fmu2Py->push_back(mu2->currentState().globalMomentum().y());
                            fmu2Pz->push_back(mu2->currentState().globalMomentum().z());
                            fmu2E->push_back(mu2->currentState().kinematicParameters().energy());

                            // Close Bs Reconstruction Part // changed by ASLI

                            vector < TransientVertex > pvs;
                            if (addXlessPrimaryVertex_) ///// Xless revertexing //////
                            {
                                VertexReProducer revertex(recVtxs, iEvent);
                                Handle < TrackCollection > pvtracks;
                                iEvent.getByToken(revtxtrks_, pvtracks);
                                Handle < reco::BeamSpot > pvbeamspot;
                                iEvent.getByToken(thebeamspot_, pvbeamspot);
                                if (pvbeamspot.id() != beamSpotHandle.id())
                                    edm::
                                    LogWarning("Inconsistency") <<
                                                                "The BeamSpot used for PV reco is not the same used in this analyzer.";
                                const reco::Muon * rmu_1 = dynamic_cast < const reco::Muon * >(iMuonP->originalObject());
                                const reco::Muon * rmu_2 = dynamic_cast < const reco::Muon * >(iMuonM->originalObject());
                                if (rmu_1 != 0 && rmu_2 != 0 && rmu_1->track().id() == pvtracks.id()
                                        && rmu_2->track().id() == pvtracks.id() && iTrackP->track().id() == pvtracks.id()
                                        && iTrackM->track().id() == pvtracks.id()
                                   )
                                {
                                    TrackCollection XLess;
                                    XLess.reserve(pvtracks->size());
                                    for (size_t i = 0, n = pvtracks->size(); i < n; ++i)
                                    {
                                        if (i == rmu_1->track().key())   continue;
                                        if (i == rmu_2->track().key())   continue;
                                        if (i == iTrackP->track().key()) continue;
                                        if (i == iTrackM->track().key()) continue;
                                        XLess.push_back((*pvtracks)[i]);
                                   }
                                    pvs = revertex.makeVertices(XLess, *pvbeamspot, iSetup);
                                    if (!pvs.empty())
                                    {
                                        Vertex XLessPV = Vertex(pvs.front());
                                        thePrimaryV = XLessPV;
                                    }
                                }
                            }   // if (addXlessPrimaryVertex_)

                            if (resolveAmbiguity_)
                            {

                                float minDZ = 999999.;
                                for (VertexCollection::const_iterator itv = recVtxs->begin(), itvend =
                                            recVtxs->end(); itv != itvend; ++itv)
                                {
                                    float deltaZ = fabs((*xDecayVertexMC).position().z() - itv->position().z());
                                    if (deltaZ < minDZ)   // Vertex to be chosen the closest one along the z axis. 
                                    {
                                        minDZ = deltaZ;
                                        theRecoVtx = Vertex(*itv);
                                    }
                                }

                               ///// Primary Vertex with X /////
                                MyPriVtxX->push_back(theRecoVtx.position().x());
                                MyPriVtxY->push_back(theRecoVtx.position().y());
                                MyPriVtxZ->push_back(theRecoVtx.position().z());
                                MyPriVtxXE->push_back(theRecoVtx.xError());
                                MyPriVtxYE->push_back(theRecoVtx.yError());
                                MyPriVtxZE->push_back(theRecoVtx.zError());
                                MyPriVtxChiNorm->push_back(theRecoVtx.normalizedChi2());
                                MyPriVtxChi->push_back(theRecoVtx.chi2());
                                MyPriVtxCL->push_back(ChiSquaredProbability((double) (theRecoVtx.chi2()), (double) (theRecoVtx.ndof())));
                                float minDz = 999999.;
                                if (!addXlessPrimaryVertex_)
                                {
                                    for (VertexCollection::const_iterator itv = recVtxs->begin(), itvend =
                                                recVtxs->end(); itv != itvend; ++itv)
                                    {
                                        float deltaZ = fabs((*xDecayVertexMC).position().z() - itv->position().z());
                                        if (deltaZ < minDz)
                                        {
                                            minDz = deltaZ;
                                            thePrimaryV = Vertex(*itv);
                                        }
                                    }
                                }
                                else
                                {
                                    for (vector < TransientVertex >::iterator itv2 = pvs.begin(), itvend2 =
                                                pvs.end(); itv2 != itvend2; ++itv2)
                                    {
                                        float deltaZ = fabs((*xDecayVertexMC).position().z() - itv2->position().z());
                                        if (deltaZ < minDz)
                                        {
                                            minDz = deltaZ;
                                            Vertex XLessPV = Vertex(*itv2);
                                            thePrimaryV = XLessPV;
                                        }
                                    }
                                }
                            }   // if (resolveAmbiguity_)

                            ///////// PV without X /////////////////////////
                            PriVtxXCorrX->push_back(thePrimaryV.position().x());
                            PriVtxXCorrY->push_back(thePrimaryV.position().y());
                            PriVtxXCorrZ->push_back(thePrimaryV.position().z());
                            PriVtxXCorrEX->push_back(thePrimaryV.xError());
                            PriVtxXCorrEY->push_back(thePrimaryV.yError());
                            PriVtxXCorrEZ->push_back(thePrimaryV.zError());
                            PriVtxXCorrCL->push_back(ChiSquaredProbability((double) (thePrimaryV.chi2()),(double) (thePrimaryV.ndof())));
                            PriVtxXCorrC2->push_back(thePrimaryV.chi2());
                            TVector3 vtx;
                            vtx.SetXYZ((*xDecayVertexMC).position().x(), (*xDecayVertexMC).position().y(), 0);
                            TVector3 pvtx;
                            pvtx.SetXYZ(thePrimaryV.position().x(), thePrimaryV.position().y(), 0);
                            VertexDistanceXY vdistXY;
                            TVector3 pperp(xCandMC->currentState().globalMomentum().x(), xCandMC->currentState().globalMomentum().y(), 0);
                            AlgebraicVector vpperp(3);
                            vpperp[0] = pperp.x();
                            vpperp[1] = pperp.y();
                            vpperp[2] = 0.;
                            TVector3 vdiff = vtx - pvtx;
                            double cosAlpha = vdiff.Dot(pperp) / (vdiff.Perp() * pperp.Perp());
                            Measurement1D distXY = vdistXY.distance(Vertex(*xDecayVertexMC), Vertex(thePrimaryV));
                            double ctauPV = distXY.value() * cosAlpha * xCandMC->currentState().mass() / pperp.Perp();
                            GlobalError v1e = (Vertex(*xDecayVertexMC)).error();
                            GlobalError v2e = thePrimaryV.error();
                            AlgebraicSymMatrix vXYe = asHepMatrix(v1e.matrix()) + asHepMatrix(v2e.matrix());
                            double ctauErrPV = sqrt(vXYe.similarity(vpperp)) * xCandMC->currentState().mass() / (pperp.Perp2());
                            float lxyPV = vdiff.Dot(pperp) / pperp.Mag();
                            double lxyPVE = sqrt(vXYe.similarity(vpperp));  
                            xCosAlpha->push_back(cosAlpha);
                            xCTauPV->push_back(ctauPV);
                            xCTauPVE->push_back(ctauErrPV);
                            xLxyPV->push_back(lxyPV);
                            xLxyPVE->push_back(lxyPVE);
                            // Lifetime BS
                            pvtx.SetXYZ(theBeamSpotV.position().x(), theBeamSpotV.position().y(), 0);
                            vdiff = vtx - pvtx;
                            cosAlpha = vdiff.Dot(pperp) / (vdiff.Perp() * pperp.Perp());
                            distXY = vdistXY.distance(Vertex(*xDecayVertexMC), Vertex(theBeamSpotV));
                            double ctauBS = distXY.value() * cosAlpha * xCandMC->currentState().mass() / pperp.Perp();
                            GlobalError v1eB = (Vertex(*xDecayVertexMC)).error();
                            GlobalError v2eB = theBeamSpotV.error();
                            AlgebraicSymMatrix vXYeB = asHepMatrix(v1eB.matrix()) + asHepMatrix(v2eB.matrix());
                            double ctauErrBS = sqrt(vXYeB.similarity(vpperp)) * xCandMC->currentState().mass() / (pperp.Perp2());
                            float lxyBS = vdiff.Dot(pperp) / pperp.Mag();
                            double lxyBSE = sqrt(vXYeB.similarity(vpperp)); 
                            xCosAlphaBS->push_back(cosAlpha);
                            xCTauBS->push_back(ctauBS);
                            xCTauBSE->push_back(ctauErrBS);
                            xLxyBS->push_back(lxyBS);
                            xLxyBSE->push_back(lxyBSE);
                            JIndex->push_back(nJ - 1);
                            pipIdx->push_back(std::distance(thePATTrackHandle->begin(), iTrackP));
                            pimIdx->push_back(std::distance(thePATTrackHandle->begin(), iTrackM));
                            pi3rdIdx->push_back(std::distance(thePATTrackHandle->begin(), iTrack3rd));

                            nX++;
                            vFitMCParticles.clear();
                       }        // 3rd track
                    }           // 2nd loop over track (look for pi-)
                }               // 1st loop over track (look for pi+)
            }                   // 2nd loop over muons (look for mu-)
        }                       // first loop over muons (look for mu+)
    }                           // if two muons


// fill the tree and clear the vectors
    if (nJ > 0 || MCnX > 0)
    {
        X_One_Tree_->Fill();
    }

    /*if ( MCnX > 0 ) {
                X_One_Tree_->Fill();
    } */
/*    if (Debug_)
    {
        cout << "Resetting branches, had " << nX << " X cands and " << nJ << " J cands." << endl;
    }*/

    trigRes->clear();
    trigNames->clear();
    L1TT->clear();
    MatchTriggerNames->clear(); // trigPreScl->clear();MatchTriggerNames
    runNum = 0;
    evtNum = 0;
    lumiNum = 0;
    nJPsi = 0;
    nJ = 0;
    nX = 0;
    nMu = 0;
    nMC = 0;
    MCnX = 0;
    nK = 0;
    priVtxX = 0;
    priVtxY = 0;
    priVtxZ = 0;
    priVtxXE = 0;
    priVtxYE = 0;
    priVtxZE = 0;
    priVtxChiNorm = 0;
    priVtxChi = 0;
    priVtxCL = 0;
    xMass->clear();
    xVtxCL->clear();
    xVtxC2->clear();
    PsiTwoSMass->clear();
    SecondPairMass->clear();
    BsMass->clear();
    BsVtxCL->clear();
    BsVtxC2->clear();
    xPx->clear();
    xPy->clear();
    xPz->clear();
    xPxE->clear();
    xPyE->clear();
    xPzE->clear();
    xDecayVtxX->clear();
    xDecayVtxY->clear();
    xDecayVtxZ->clear();
    xDecayVtxXE->clear();
    xDecayVtxYE->clear();
    xDecayVtxZE->clear();
    MyPriVtxX->clear();
    MyPriVtxY->clear();
    MyPriVtxZ->clear();
    MyPriVtxXE->clear();
    MyPriVtxYE->clear();
    MyPriVtxZE->clear();
    MyPriVtxChiNorm->clear();
    MyPriVtxChi->clear();
    MyPriVtxCL->clear();
    PriVtxXCorrX->clear();
    PriVtxXCorrY->clear();
    PriVtxXCorrZ->clear();
    PriVtxXCorrEX->clear();
    PriVtxXCorrEY->clear();
    PriVtxXCorrEZ->clear();
    PriVtxXCorrC2->clear();
    PriVtxXCorrCL->clear();
    xLxyPV->clear();
    xLxyPVE->clear();
    xCosAlpha->clear();
    xCTauPV->clear();
    xLxyBS->clear();
    xLxyBSE->clear();
    xExtraTrkabsDxywrtXvtxmin->clear();
    xExtraTrkabsDxywrtXvtxminMaxDz0dot01cm->clear();
    xExtraTrkabsDxywrtXvtxminMaxDz0dot02cm->clear();
    xExtraTrkabsDxywrtXvtxminMaxDz0dot03cm->clear();
    xExtraTrkabsDxywrtXvtxminMaxDz0dot04cm->clear();
    xExtraTrkabsDxywrtXvtxminMaxDz0dot05cm->clear();
    xExtraTrkabsDxywrtXvtxminMaxDz0dot1cm->clear();
    xExtraTrkabsDxywrtXvtxminMaxDz0dot5cm->clear();
    xExtraTrkabsDxywrtXvtxminMaxDz1dot0cm->clear();
    xExtraTrkabsDxywrtXvtxminMaxDz1dot5cm->clear();
    xExtraTrkabsDxywrtXvtxminMaxDz2dot0cm->clear();
    xExtraTrkabsDxywrtXvtxminMaxDz3dot0cm->clear();
    xExtraTrkabsDxywrtXvtxminMaxDz4dot0cm->clear();
    xExtraTrkabsDxywrtXvtxminMaxDz5dot0cm->clear();
    xCosAlphaBS->clear();
    xCTauBS->clear();
    xCTauPVE->clear();
    xCTauBSE->clear();
    JIndex->clear();
    pipIdx->clear();
    pimIdx->clear();
    pi3rdIdx->clear();
    JMass->clear();
    JVtxCL->clear();
    JVtxC2->clear();
    JPx->clear();
    JPy->clear();
    JPz->clear();
    JDecayVtxX->clear();
    JDecayVtxY->clear();
    JDecayVtxZ->clear();
    JDecayVtxXE->clear();
    JDecayVtxYE->clear();
    JDecayVtxZE->clear();
    mupIdx->clear();
    mumIdx->clear();
    mupTrigMatch->clear();
    mumTrigMatch->clear();
    JPsiMuonTrigMatch->clear();
    mumPx->clear();
    mumPy->clear();
    mumPz->clear();
    mupPx->clear();
    mupPy->clear();
    mupPz->clear();
    mumfChi2->clear();
    mumfNDF->clear();
    mupfChi2->clear();
    mupfNDF->clear();
    muPx->clear();
    muPy->clear();
    muPz->clear();
    muD0->clear();
    muD0E->clear();
    cosAlpha->clear();
    mumuVtxCL->clear();
    mumuFLSig->clear();
    mumurVtxMag2D->clear();
    mumusigmaRvtxMag2D->clear();
    muDz->clear();
    muChi2->clear();
    muGlChi2->clear();
    mufHits->clear();
    muFirstBarrel->clear();
    muFirstEndCap->clear();
    muDzVtx->clear();
    muDxyVtx->clear();
    muNDF->clear();
    muGlNDF->clear();
    muPhits->clear();
    muShits->clear();
    muGlMuHits->clear();
    muType->clear();
    muQual->clear();
    muTrack->clear();
    muCharge->clear();
    ThreeGoodTracks->clear();
    trQuality->clear();
    trPx->clear();
    trPy->clear();
    trPz->clear();
    trE->clear();
    trNDF->clear();
    trPhits->clear();
    trShits->clear();
    trChi2->clear();
    trD0->clear();
    trD0E->clear();
    trCharge->clear();
    trfHits->clear();
    trFirstBarrel->clear();
    trFirstEndCap->clear();
    trDzVtx->clear();
    trDxyVtx->clear();
    trVx->clear();
    trVy->clear();
    tr_nsigdedx->clear();
    tr_dedx->clear();
    tr_dedxMass->clear();
    tr_theo->clear();
    tr_sigma->clear();
    fpi1Px->clear();
    fpi1Py->clear();
    fpi1Pz->clear();
    fpi1E->clear();
    fpi2Px->clear();
    fpi2Py->clear();
    fpi2Pz->clear();
    fpi2E->clear();
    fpi3Px->clear();
    fpi3Py->clear();
    fpi3Pz->clear();
    fpi3E->clear();
    fmu1Px->clear();
    fmu1Py->clear();
    fmu1Pz->clear();
    fmu1E->clear();
    fmu2Px->clear();
    fmu2Py->clear();
    fmu2Pz->clear();
    fmu2E->clear();
    MCPdgId->clear();
    MCDauphi_Pt->clear();
    MCDauphi_mass->clear();
    MCDauJpsi_Pt->clear();
    MCDauphi_Px->clear();
    MCDauphi_Py->clear();
    MCDauphi_Pz->clear();
    MCDauphi_E->clear();
    BDauIdx->clear();
    MCDauJpsi_Px->clear();
    MCDauJpsi_Py->clear();
    MCDauJpsi_Pz->clear();
    MCDauJpsi_mass->clear();
    MCDauJpsi_E->clear();
    MCmupDau_M->clear();
    MCmupDau_px->clear();
    MCmupDau_py->clear();
    MCmupDau_pz->clear();
    MCmupDau_E->clear();
    MCmumDau_M->clear();
    MCmumDau_px->clear();
    MCmumDau_py->clear();
    MCmumDau_pz->clear();
    MCmumDau_E->clear();
    MCkpDau_M->clear();
    MCkpDau_px->clear();
    MCkpDau_py->clear();
    MCkpDau_pz->clear();
    MCkpDau_E->clear();
    MCkmDau_M->clear();
    MCkmDau_px->clear();
    MCkmDau_py->clear();
    MCkmDau_pz->clear();
    MCkmDau_E->clear();
    MCDauK_Pt->clear();
    MCDauK_Px->clear();
    MCDauK_Py->clear();
    MCDauK_Pz->clear();
    MCDauK_mass->clear();
    MCDauK_E->clear();
    MCmass->clear();
    MCPx->clear();
    MCPy->clear();
    MCPz->clear();
    MCE->clear();
    MCJidx->clear();
    MCPhidx->clear();
    MCKidx->clear();
    MCBIdx->clear();
MCStatus->clear();
    MCPt->clear();
    // }                           // if(doMC)
    // clear*/
}                               // analyze


// ------------ method called once each job just before starting event loop  ------------
void JPsiKKKPATRunII::beginRun(edm::Run const &iRun, edm::EventSetup const &iSetup)
{
    // bool changed = true;
    // proccessName_="HLT";
    // hltConfig_.init(iRun,iSetup,proccessName_,changed);
}

void JPsiKKKPATRunII::beginJob()
{
    edm::Service < TFileService > fs;
    X_One_Tree_ = fs->make < TTree > ("X_data", "X(3872) Data");
    X_One_Tree_->Branch("TrigRes", &trigRes);
    X_One_Tree_->Branch("TrigNames", &trigNames);
    X_One_Tree_->Branch("MatchTriggerNames", &MatchTriggerNames);
    X_One_Tree_->Branch("L1TrigRes", &L1TT);
    X_One_Tree_->Branch("evtNum", &evtNum, "evtNum/i");
    X_One_Tree_->Branch("runNum", &runNum, "runNum/i");
    X_One_Tree_->Branch("lumiNum", &lumiNum, "lumiNum/i");
    X_One_Tree_->Branch("priVtxX", &priVtxX, "priVtxX/f");
    X_One_Tree_->Branch("priVtxY", &priVtxY, "priVtxY/f");
    X_One_Tree_->Branch("priVtxZ", &priVtxZ, "priVtxZ/f");
    X_One_Tree_->Branch("priVtxXE", &priVtxXE, "priVtxXE/f");
    X_One_Tree_->Branch("priVtxYE", &priVtxYE, "priVtxYE/f");
    X_One_Tree_->Branch("priVtxZE", &priVtxZE, "priVtxZE/f");
    X_One_Tree_->Branch("priVtxChiNorm", &priVtxChiNorm, "priVtxChiNorm/f");
    X_One_Tree_->Branch("priVtxChi", &priVtxChi, "priVtxChi/f");
    X_One_Tree_->Branch("priVtxCL", &priVtxCL, "priVtxCL/f");
    X_One_Tree_->Branch("nMu", &nMu, "nMu/i");
    X_One_Tree_->Branch("muPx", &muPx);
    X_One_Tree_->Branch("muPy", &muPy);
    X_One_Tree_->Branch("muPz", &muPz);
    X_One_Tree_->Branch("muD0", &muD0);
    X_One_Tree_->Branch("muD0E", &muD0E);
    X_One_Tree_->Branch("cosAlpha", &cosAlpha);
    X_One_Tree_->Branch("mumuVtxCL", &mumuVtxCL);
    X_One_Tree_->Branch("mumuFLSig", &mumuFLSig);
    X_One_Tree_->Branch("mumurVtxMag2D", &mumurVtxMag2D);
    X_One_Tree_->Branch("mumusigmaRvtxMag2D", &mumusigmaRvtxMag2D);
    X_One_Tree_->Branch("muDz", &muDz);
    X_One_Tree_->Branch("muChi2", &muChi2);
    X_One_Tree_->Branch("muNDF", &muNDF);
    X_One_Tree_->Branch("muPhits", &muPhits);
    X_One_Tree_->Branch("muShits", &muShits);
    X_One_Tree_->Branch("muGlChi2", &muGlChi2);
    X_One_Tree_->Branch("muGlNDF", &muGlNDF);
    X_One_Tree_->Branch("muGlMuHits", &muGlMuHits);
    X_One_Tree_->Branch("muType", &muType);
    X_One_Tree_->Branch("muQual", &muQual);
    X_One_Tree_->Branch("muTrack", &muTrack);
    X_One_Tree_->Branch("muCharge", &muCharge);
    X_One_Tree_->Branch("mufHits", &mufHits);
    X_One_Tree_->Branch("muFirstBarrel", &muFirstBarrel);
    X_One_Tree_->Branch("muFirstEndCap", &muFirstEndCap);
    X_One_Tree_->Branch("muDzVtx", &muDzVtx);
    X_One_Tree_->Branch("muDxyVtx", &muDxyVtx);
    X_One_Tree_->Branch("ThreeGoodTracks", &ThreeGoodTracks);
    X_One_Tree_->Branch("trackQuality", &trQuality);
    X_One_Tree_->Branch("TrackPx", &trPx);
    X_One_Tree_->Branch("TrackPy", &trPy);
    X_One_Tree_->Branch("TrackPz", &trPz);
    X_One_Tree_->Branch("TrackEnergy", &trE);
    X_One_Tree_->Branch("TrackNDF", &trNDF);
    X_One_Tree_->Branch("TrackPhits", &trPhits);
    X_One_Tree_->Branch("TrackShits", &trShits);
    X_One_Tree_->Branch("TrackChi2", &trChi2);
    X_One_Tree_->Branch("TrackD0", &trD0);
    X_One_Tree_->Branch("TrackD0Err", &trD0E);
    X_One_Tree_->Branch("TrackCharge", &trCharge);
    X_One_Tree_->Branch("trfHits", &trfHits);
    X_One_Tree_->Branch("trFirstBarrel", &trFirstBarrel);
    X_One_Tree_->Branch("trFirstEndCap", &trFirstEndCap);
    X_One_Tree_->Branch("trDzVtx", &trDzVtx);
    X_One_Tree_->Branch("trDxyVtx", &trDxyVtx);
    X_One_Tree_->Branch("trVx", &trVx);
    X_One_Tree_->Branch("trVy", &trVy);
    X_One_Tree_->Branch("tr_nsigdedx", &tr_nsigdedx);
    X_One_Tree_->Branch("tr_dedx", &tr_dedx);
    X_One_Tree_->Branch("tr_dedxMass", &tr_dedxMass);
    X_One_Tree_->Branch("tr_theo", &tr_theo);
    X_One_Tree_->Branch("tr_sigma", &tr_sigma);
    X_One_Tree_->Branch("nX", &nX, "nX/i");
    X_One_Tree_->Branch("nK", &nK, "nK/i");
    X_One_Tree_->Branch("nMC", &nMC, "nMC/i");
    X_One_Tree_->Branch("nJPsi", &nJPsi, "nJPsi/i");
    X_One_Tree_->Branch("MCnX", &MCnX, "MCnX/i");
    X_One_Tree_->Branch("xMass", &xMass);
    X_One_Tree_->Branch("xVtxCL", &xVtxCL);
    X_One_Tree_->Branch("xVtxC2", &xVtxC2);
    X_One_Tree_->Branch("PsiTwoSMass", &PsiTwoSMass);
    X_One_Tree_->Branch("SecondPairMass", &SecondPairMass);
    X_One_Tree_->Branch("BsMass", &BsMass);
    X_One_Tree_->Branch("BsVtxCL", &BsVtxCL);
    X_One_Tree_->Branch("BsVtxC2", &BsVtxC2);
    X_One_Tree_->Branch("xPx", &xPx);
    X_One_Tree_->Branch("xPy", &xPy);
    X_One_Tree_->Branch("xPz", &xPz);
    X_One_Tree_->Branch("xPxE", &xPxE);
    X_One_Tree_->Branch("xPyE", &xPyE);
    X_One_Tree_->Branch("xPzE", &xPzE);
    X_One_Tree_->Branch("xDecayVtxX", &xDecayVtxX);
    X_One_Tree_->Branch("xDecayVtxY", &xDecayVtxY);
    X_One_Tree_->Branch("xDecayVtxZ", &xDecayVtxZ);
    X_One_Tree_->Branch("xDecayVtxXE", &xDecayVtxXE);
    X_One_Tree_->Branch("xDecayVtxYE", &xDecayVtxYE);
    X_One_Tree_->Branch("xDecayVtxZE", &xDecayVtxZE);
    X_One_Tree_->Branch("PriVtxXCorrX", &PriVtxXCorrX);
    X_One_Tree_->Branch("PriVtxXCorrY", &PriVtxXCorrY);
    X_One_Tree_->Branch("PriVtxXCorrZ", &PriVtxXCorrZ);
    X_One_Tree_->Branch("PriVtxXCorrEX", &PriVtxXCorrEX);
    X_One_Tree_->Branch("PriVtxXCorrEY", &PriVtxXCorrEY);
    X_One_Tree_->Branch("PriVtxXCorrEZ", &PriVtxXCorrEZ);
    X_One_Tree_->Branch("PriVtxXCorrC2", &PriVtxXCorrC2);
    X_One_Tree_->Branch("PriVtxXCorrCL", &PriVtxXCorrCL);
    X_One_Tree_->Branch("MyPriVtxX", &MyPriVtxX);
    X_One_Tree_->Branch("MyPriVtxY", &MyPriVtxY);
    X_One_Tree_->Branch("MyPriVtxZ", &MyPriVtxZ);
    X_One_Tree_->Branch("MyPriVtxXE", &MyPriVtxXE);
    X_One_Tree_->Branch("MyPriVtxYE", &MyPriVtxYE);
    X_One_Tree_->Branch("MyPriVtxZE", &MyPriVtxZE);
    X_One_Tree_->Branch("MyPriVtxChiNorm", &MyPriVtxChiNorm);
    X_One_Tree_->Branch("MyPriVtxChi", &MyPriVtxChi);
    X_One_Tree_->Branch("MyPriVtxCL", &MyPriVtxCL);
    X_One_Tree_->Branch("xLxyPV", &xLxyPV);
    X_One_Tree_->Branch("xLxyPVE", &xLxyPVE);
    X_One_Tree_->Branch("xCosAlpha", &xCosAlpha);
    X_One_Tree_->Branch("xCTauPV", &xCTauPV);
    X_One_Tree_->Branch("xCTauPVE", &xCTauPVE);
    X_One_Tree_->Branch("xLxyBS", &xLxyBS);
    X_One_Tree_->Branch("xLxyBSE", &xLxyBSE);
    X_One_Tree_->Branch("xExtraTrkabsDxywrtXvtxmin", &xExtraTrkabsDxywrtXvtxmin);
    X_One_Tree_->Branch("xExtraTrkabsDxywrtXvtxminMaxDz0dot5cm", &xExtraTrkabsDxywrtXvtxminMaxDz0dot5cm);
    X_One_Tree_->Branch("xExtraTrkabsDxywrtXvtxminMaxDz1dot0cm", &xExtraTrkabsDxywrtXvtxminMaxDz1dot0cm);
    X_One_Tree_->Branch("xExtraTrkabsDxywrtXvtxminMaxDz1dot5cm", &xExtraTrkabsDxywrtXvtxminMaxDz1dot5cm);
    X_One_Tree_->Branch("xExtraTrkabsDxywrtXvtxminMaxDz2dot0cm", &xExtraTrkabsDxywrtXvtxminMaxDz2dot0cm);
    X_One_Tree_->Branch("xExtraTrkabsDxywrtXvtxminMaxDz3dot0cm", &xExtraTrkabsDxywrtXvtxminMaxDz3dot0cm);
    X_One_Tree_->Branch("xExtraTrkabsDxywrtXvtxminMaxDz4dot0cm", &xExtraTrkabsDxywrtXvtxminMaxDz4dot0cm);
    X_One_Tree_->Branch("xExtraTrkabsDxywrtXvtxminMaxDz5dot0cm", &xExtraTrkabsDxywrtXvtxminMaxDz5dot0cm);
    X_One_Tree_->Branch("xExtraTrkabsDxywrtXvtxminMaxDz0dot1cm", &xExtraTrkabsDxywrtXvtxminMaxDz0dot1cm);
    X_One_Tree_->Branch("xExtraTrkabsDxywrtXvtxminMaxDz0dot05cm", &xExtraTrkabsDxywrtXvtxminMaxDz0dot05cm);
    X_One_Tree_->Branch("xExtraTrkabsDxywrtXvtxminMaxDz0dot04cm", &xExtraTrkabsDxywrtXvtxminMaxDz0dot04cm);
    X_One_Tree_->Branch("xExtraTrkabsDxywrtXvtxminMaxDz0dot03cm", &xExtraTrkabsDxywrtXvtxminMaxDz0dot03cm);
    X_One_Tree_->Branch("xExtraTrkabsDxywrtXvtxminMaxDz0dot02cm", &xExtraTrkabsDxywrtXvtxminMaxDz0dot02cm);
    X_One_Tree_->Branch("xExtraTrkabsDxywrtXvtxminMaxDz0dot01cm", &xExtraTrkabsDxywrtXvtxminMaxDz0dot01cm);
    X_One_Tree_->Branch("xCosAlphaBS", &xCosAlphaBS);
    X_One_Tree_->Branch("xCTauBS", &xCTauBS);
    X_One_Tree_->Branch("xCTauBSE", &xCTauBSE);
    X_One_Tree_->Branch("JPsiIndex", &JIndex);
    X_One_Tree_->Branch("pipIndex", &pipIdx);
    X_One_Tree_->Branch("pimIndex", &pimIdx);
    X_One_Tree_->Branch("pi3rdIndex", &pi3rdIdx);
    X_One_Tree_->Branch("fitpi1Px", &fpi1Px);
    X_One_Tree_->Branch("fitpi1Py", &fpi1Py);
    X_One_Tree_->Branch("fitpi1Pz", &fpi1Pz);
    X_One_Tree_->Branch("fitpi1E", &fpi1E);
    X_One_Tree_->Branch("fitpi2Px", &fpi2Px);
    X_One_Tree_->Branch("fitpi2Py", &fpi2Py);
    X_One_Tree_->Branch("fitpi2Pz", &fpi2Pz);
    X_One_Tree_->Branch("fitpi2E", &fpi2E);
    X_One_Tree_->Branch("fitpi3Px", &fpi3Px);
    X_One_Tree_->Branch("fitpi3Py", &fpi3Py);
    X_One_Tree_->Branch("fitpi3Pz", &fpi3Pz);
    X_One_Tree_->Branch("fitpi3E", &fpi3E);
    X_One_Tree_->Branch("fitmu1Px", &fmu1Px);
    X_One_Tree_->Branch("fitmu1Py", &fmu1Py);
    X_One_Tree_->Branch("fitmu1Pz", &fmu1Pz);
    X_One_Tree_->Branch("fitmu1E", &fmu1E);
    X_One_Tree_->Branch("fitmu2Px", &fmu2Px);
    X_One_Tree_->Branch("fitmu2Py", &fmu2Py);
    X_One_Tree_->Branch("fitmu2Pz", &fmu2Pz);
    X_One_Tree_->Branch("fitmu2E", &fmu2E);
    X_One_Tree_->Branch("nJ", &nJ, "nJ/i");
    X_One_Tree_->Branch("JMass", &JMass);
    X_One_Tree_->Branch("JVtxCL", &JVtxCL);
    X_One_Tree_->Branch("JVtxC2", &JVtxC2);
    X_One_Tree_->Branch("JPx", &JPx);
    X_One_Tree_->Branch("JPy", &JPy);
    X_One_Tree_->Branch("JPz", &JPz);
    X_One_Tree_->Branch("JDecayVtxX", &JDecayVtxX);
    X_One_Tree_->Branch("JDecayVtxY", &JDecayVtxY);
    X_One_Tree_->Branch("JDecayVtxZ", &JDecayVtxZ);
    X_One_Tree_->Branch("JDecayVtxXE", &JDecayVtxXE);
    X_One_Tree_->Branch("JDecayVtxYE", &JDecayVtxYE);
    X_One_Tree_->Branch("JDecayVtxZE", &JDecayVtxZE);
    X_One_Tree_->Branch("mupIdx", &mupIdx);
    X_One_Tree_->Branch("mumIdx", &mumIdx);
    X_One_Tree_->Branch("mumPx", &mumPx);
    X_One_Tree_->Branch("mumPy", &mumPy);
    X_One_Tree_->Branch("mumPz", &mumPz);
    X_One_Tree_->Branch("mupPx", &mupPx);
    X_One_Tree_->Branch("mupPy", &mupPy);
    X_One_Tree_->Branch("mupPz", &mupPz);
    X_One_Tree_->Branch("mumfChi2", &mumfChi2);
    X_One_Tree_->Branch("mumfNDF", &mumfNDF);
    X_One_Tree_->Branch("mupfChi2", &mupfChi2);
    X_One_Tree_->Branch("mupfNDF", &mupfNDF);
    X_One_Tree_->Branch("JPsiMuonTrigMatch", &JPsiMuonTrigMatch);
    X_One_Tree_->Branch("mupTrigMatch", &mupTrigMatch);
    X_One_Tree_->Branch("mumTrigMatch", &mumTrigMatch);
    X_One_Tree_->Branch("MCPdgId", &MCPdgId);
    X_One_Tree_->Branch("MCDauphi_Pt", &MCDauphi_Pt);
    X_One_Tree_->Branch("MCDauphi_mass",&MCDauphi_mass);
    X_One_Tree_->Branch("MCDauphi_Px", & MCDauphi_Px);
    X_One_Tree_->Branch("MCDauphi_Py", & MCDauphi_Py);
    X_One_Tree_->Branch("MCDauphi_Pz", & MCDauphi_Pz);
    X_One_Tree_->Branch("MCDauphi_E", &MCDauphi_E);
    X_One_Tree_->Branch("MCmupDau_M",&MCmupDau_M);
    X_One_Tree_->Branch("MCmupDau_px",&MCmupDau_px);
    X_One_Tree_->Branch("MCmupDau_py",&MCmupDau_py);
    X_One_Tree_->Branch("MCmupDau_pz",&MCmupDau_pz);
    X_One_Tree_->Branch("MCmupDau_E",&MCmupDau_E);
    X_One_Tree_->Branch("MCmumDau_M",&MCmumDau_M);
    X_One_Tree_->Branch("MCmumDau_px",&MCmumDau_px);
    X_One_Tree_->Branch("MCmumDau_py",&MCmumDau_py);
    X_One_Tree_->Branch("MCmumDau_pz",&MCmumDau_pz);
    X_One_Tree_->Branch("MCmumDau_E",&MCmumDau_E);
    X_One_Tree_->Branch("MCkpDau_M", &MCkpDau_M);
    X_One_Tree_->Branch("MCkpDau_px", &MCkpDau_px);
    X_One_Tree_->Branch("MCkpDau_py", &MCkpDau_py);
    X_One_Tree_->Branch("MCkpDau_pz", &MCkpDau_pz);
    X_One_Tree_->Branch("MCkpDau_E", &MCkpDau_E);
    X_One_Tree_->Branch("MCkmDau_M", &MCkmDau_M);
    X_One_Tree_->Branch("MCkmDau_px", &MCkmDau_px);
    X_One_Tree_->Branch("MCkmDau_py", &MCkmDau_py);
    X_One_Tree_->Branch("MCkmDau_pz", &MCkmDau_pz);
    X_One_Tree_->Branch("MCkmDau_E", &MCkmDau_E);
    X_One_Tree_->Branch("BDauIdx", &BDauIdx);
    X_One_Tree_->Branch("MCDauJpsi_Pt", &MCDauJpsi_Pt);
    X_One_Tree_->Branch("MCDauJpsi_Px", &MCDauJpsi_Px);
    X_One_Tree_->Branch("MCDauJpsi_Py", &MCDauJpsi_Py);
    X_One_Tree_->Branch("MCDauJpsi_Pz", &MCDauJpsi_Pz);
    X_One_Tree_->Branch("MCDauJpsi_mass", &MCDauJpsi_mass);
    X_One_Tree_->Branch("MCDauJpsi_E", &MCDauJpsi_E);
    X_One_Tree_->Branch("MCDauK_Pt", &MCDauK_Pt);
    X_One_Tree_->Branch("MCDauK_Px", &MCDauK_Px);
    X_One_Tree_->Branch("MCDauK_Py", &MCDauK_Py);
    X_One_Tree_->Branch("MCDauK_Pz", &MCDauK_Pz);
    X_One_Tree_->Branch("MCDauK_mass", &MCDauK_mass);
    X_One_Tree_->Branch("MCDauK_E", &MCDauK_E);
    X_One_Tree_->Branch("MCmass", &MCmass);
    X_One_Tree_->Branch("MCPx", &MCPx);
    X_One_Tree_->Branch("MCPy", &MCPy);
    X_One_Tree_->Branch("MCPz", &MCPz);
    X_One_Tree_->Branch("MCPt", &MCPt);
    X_One_Tree_->Branch("MCE", &MCE);
    X_One_Tree_->Branch("MCJidx", &MCJidx);
    X_One_Tree_->Branch("MCPhidx", &MCPhidx);
    X_One_Tree_->Branch("MCKidx", &MCKidx);
    X_One_Tree_->Branch("MCBIdx", &MCBIdx);
    X_One_Tree_->Branch("MCStatus" , &MCStatus);
}                               // begin Job

// ------------ method called once each job just after ending the event loop  ------------
void JPsiKKKPATRunII::endJob()
{
    X_One_Tree_->GetDirectory()->cd();
    X_One_Tree_->Write();
}                               // endjob
//define this as a plug-in
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(JPsiKKKPATRunII);
