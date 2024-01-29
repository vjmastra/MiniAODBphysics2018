#ifndef _JPsiPi_h
#define _JPsiPi_h

// system include files
#include <memory>

// user include files
//#include "myAnalyzers/JPsiKsPAT/interface/JPsif0PAT.h"

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
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"

#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackFromFTSFactory.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/V0Candidate/interface/V0Candidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"

#include "DataFormats/PatCandidates/interface/PackedCandidate.h" // muy importante para MiniAOD
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"

#include "CondFormats/L1TObjects/interface/L1GtTriggerMenu.h"
#include "CondFormats/DataRecord/interface/L1GtTriggerMenuRcd.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetupFwd.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerObjectMapRecord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetup.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"

#include "RecoVertex/VertexPrimitives/interface/BasicSingleVertexState.h"
#include "RecoVertex/VertexPrimitives/interface/VertexState.h"

#include "TFile.h"
#include "TTree.h"

#include "TLorentzVector.h"
#include "TVector3.h"

//
// class decleration
//

class JPsiPi : public edm::EDAnalyzer {
public:
  explicit JPsiPi(const edm::ParameterSet&);
  ~JPsiPi();
  void fillPsi(const reco::Candidate& genpsi);
  int const getMuCat(reco::Muon const& muon) const;
//  bool IsTheSame(const reco::Track& tk, const pat::Muon& mu);
  bool IsTheSame(const pat::PackedCandidate& tk, const pat::Muon& mu);
  bool IsTheSame2(const pat::PackedCandidate& tk1, const pat::PackedCandidate& tk2); 
  Double_t Distance(const Double_t p1[], const Double_t p2[]);
  Double_t DistanceError(const Double_t p1[], const Double_t err1[], const Double_t p2[], const Double_t err2[]);
  double ProperTime(const double len[], const double mom[]);
  double ProperTimeErr(const double len[], const double errlen[], const double mom[], const double errmom[]);
  float DeltaR(const pat::PackedCandidate t1, const pat::TriggerObjectStandAlone t2);
  float DeltaPt(const pat::PackedCandidate t1, const pat::TriggerObjectStandAlone t2);
  bool MatchByDRDPt(const pat::PackedCandidate t1, const pat::TriggerObjectStandAlone t2);
  float DeltaR(const pat::Muon t1, const pat::TriggerObjectStandAlone t2);
  float DeltaPt(const pat::Muon t1, const pat::TriggerObjectStandAlone t2);
  bool MatchByDRDPt(const pat::Muon t1, const pat::TriggerObjectStandAlone t2);
  bool   isAncestor(const reco::Candidate*, const reco::Candidate*);
  double GetLifetime(TLorentzVector, TVector3, TVector3);

private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  void printout(const RefCountedKinematicVertex& myVertex) const;
  void printout(const RefCountedKinematicParticle& myParticle) const;
  void printout(const RefCountedKinematicTree& myTree) const;
 
  // ----------member data ---------------------------
  edm::EDGetTokenT<edm::View<pat::Muon>> dimuon_Label;
  edm::EDGetTokenT<edm::View<pat::PackedCandidate>> trakCollection_label;
  edm::EDGetTokenT<reco::VertexCollection> primaryVertices_Label;
  edm::EDGetTokenT<reco::VertexCollection> primaryVerticesWithBS_Label;
  edm::EDGetTokenT<reco::BeamSpot> BSLabel_;
  edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> builderToken_;
  edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerCollection_;
  edm::EDGetTokenT<edm::TriggerResults> triggerResults_Label;

  //std::string genParticles_;
  edm::EDGetTokenT<reco::GenParticleCollection> genParticles_;
  edm::EDGetTokenT<pat::PackedGenParticleCollection> packedGenToken_;

  bool OnlyBest_;
  bool isMC_;
  bool OnlyGen_;
  bool doMC_;

  TTree      *tree_;
//  TTree      *treeTest_;
  int mupCategory;
  int mumCategory;
  int mupME1Clean;
  int mumME1Clean;
  
  float       mumC2;
  int         mumNHits, mumNPHits; 
  float       mupC2;
  int         mupNHits, mupNPHits;
  float       mumdxy, mupdxy, mumdz, mupdz;
  float       muon_dca;

  bool        mu1soft, mu2soft, mu1tight, mu2tight;  
  bool        mu1PF, mu2PF, mu1loose, mu2loose;  
 
  int         muAcc, muTrig, weight;
 
  unsigned int nVtx;

  double       priVtxX, priVtxY, priVtxZ, priVtxXE, priVtxYE, priVtxZE, priVtxCL;
  double       priVtxXYE, priVtxXZE, priVtxYZE;

  double       bestVtxX, bestVtxY, bestVtxZ, bestVtxXE, bestVtxYE, bestVtxZE, bestVtxCL;
  double       bestVtxXYE, bestVtxXZE, bestVtxYZE;

  double       bestWBSVtxX, bestWBSVtxY, bestWBSVtxZ, bestWBSVtxXE, bestWBSVtxYE, bestWBSVtxZE, bestWBSVtxCL;
  double       bestWBSVtxXYE, bestWBSVtxXZE, bestWBSVtxYZE;

  double       trkVtxX, trkVtxY, trkVtxZ, trkVtxXE, trkVtxYE, trkVtxZE, trkVtxCL;
  double       trkVtxXYE, trkVtxXZE, trkVtxYZE;

  int         indexVtx;
  int         vRefTrk; 

  double       JDecayVtxX, JDecayVtxY, JDecayVtxZ;
  double      JDecayVtxXE, JDecayVtxYE, JDecayVtxZE;
  double      JDecayVtxXYE, JDecayVtxXZE, JDecayVtxYZE;

  double       BDecayVtxX, BDecayVtxY, BDecayVtxZ;
  double      BDecayVtxXE, BDecayVtxYE, BDecayVtxZE;
  double      BDecayVtxXYE, BDecayVtxXZE, BDecayVtxYZE;

  // *************************************

  unsigned int             nB;
  unsigned int             nMu;
  unsigned int             nJpsi;

  std::array<bool, 4>        triggerFlags;
  bool        mu1HLTmatched, mu2HLTmatched, trkHLT0matched, trkHLT1matched;

  float       B_mass, B_px, B_py, B_pz;
  float       B_pt, B_eta, B_phi;
  float       jPointingAngle, bPointingAngle;
  float       jPointingAngle2D, bPointingAngle2D;

  float       pi_px, pi_py, pi_pz, pi_charge;
  float       pi_pt, pi_eta, pi_phi;

  float       d0ValPi, d0ErrPi, d0SigPi;

  float       J_mass, J_px, J_py, J_pz, deltaRmumu;
  float       J_pt, J_eta, J_phi;

  float       J_px1, J_py1, J_pz1;
  float       J_pt1, J_eta1, J_phi1;

  float       J_px2, J_py2, J_pz2;
  float       J_pt2, J_eta2, J_phi2;
  int         J_charge1, J_charge2;
 
  double      jFlightLen, jFlightLenErr, jFlightLenSig;
  double      jFlightLen2D, jFlightLenErr2D, jFlightLenSig2D;
  double      bFlightLen, bFlightLenErr, bFlightLenSig;
  double      bFlightLen2D, bFlightLenErr2D, bFlightLenSig2D;
  double      bFlightTime, bFlightTimeErr;
  double      bFlightTime3D, bFlightTime3DErr;
  double      bFlightTimeM, bFlightTimeMErr;
  double      bFlightTimeBest, bFlightTimeBestErr;
  double      bFlightTimeBestWBS, bFlightTimeBestWBSErr;
  double      bFlightTimeTrk, bFlightTimeTrkErr;
  double      bFlightTimeOld, bFlightTimeOldErr;
  double      deltaR_J_pi, cosAngle_J_pi;

  float       J_chi2, B_chi2; 
  float       B_Prob, J_Prob; 

  int  run, event;
  int  lumiblock;

  TLorentzVector gen_bc_p4,gen_jpsi_p4,gen_pion3_p4,gen_muon1_p4,gen_muon2_p4;
  TVector3       gen_bc_vtx,gen_jpsi_vtx;
  float          gen_bc_ct;

  bool jpsiGenMatched;
  bool candGenMatched;
};

#endif
