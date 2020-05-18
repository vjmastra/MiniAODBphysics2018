#ifndef _Psi2Spi_h
#define _Psi2Spi_h

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

//
// class decleration
//

class Psi2Spi : public edm::EDAnalyzer {
public:
  explicit Psi2Spi(const edm::ParameterSet&);
  ~Psi2Spi();
  void fillPsi(const reco::Candidate& genpsi);
  int const getMuCat(reco::Muon const& muon) const;
//  bool IsTheSame(const reco::Track& tk, const pat::Muon& mu);
  bool IsTheSame(const pat::PackedCandidate& tk, const pat::Muon& mu);
  bool IsTheSame2(const pat::PackedCandidate& tk1, const pat::PackedCandidate& tk2); 
  Double_t Distance(const Double_t p1[], const Double_t p2[]);
  Double_t DistanceError(const Double_t p1[], const Double_t err1[], const Double_t p2[], const Double_t err2[]);
  float DeltaR(const pat::PackedCandidate t1, const pat::TriggerObjectStandAlone t2);
  float DeltaPt(const pat::PackedCandidate t1, const pat::TriggerObjectStandAlone t2);
  bool MatchByDRDPt(const pat::PackedCandidate t1, const pat::TriggerObjectStandAlone t2);

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
  edm::EDGetTokenT<reco::BeamSpot> BSLabel_;
  edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerCollection_;
  edm::EDGetTokenT<edm::TriggerResults> triggerResults_Label;

  std::string genParticles_;
  bool OnlyBest_;
  bool isMC_;
  bool OnlyGen_;
  bool doMC_;

  TTree      *tree_, *treeTest_;
  int mupCategory;
  int mumCategory;
  int mupME1Clean;
  int mumME1Clean;
  
  std::vector<float>       *mumC2;
  std::vector<int>         *mumNHits, *mumNPHits; 
  std::vector<float>       *mupC2;
  std::vector<int>         *mupNHits, *mupNPHits;
  std::vector<float>       *mumdxy, *mupdxy, *mumdz, *mupdz;
  std::vector<float>       *muon_dca;

  std::vector<int>         *tri_Dim25, *tri_JpsiTk, *tri_JpsiTkTk;
 
  std::vector<bool>        *mu1soft, *mu2soft, *mu1tight, *mu2tight;  
  std::vector<bool>        *mu1PF, *mu2PF, *mu1loose, *mu2loose;  
 
  int                      muAcc, muTrig, weight;
 
  // vertice primario CON mayor Pt
  unsigned int             nVtx;
  std::vector<float>       *priVtxX, *priVtxY, *priVtxZ, *priVtxXE, *priVtxYE, *priVtxZE, *priVtxCL;
  std::vector<float>       *priVtxXYE, *priVtxXZE, *priVtxYZE;

  std::vector<int>         *indexVtx, *nTracksFromPV;

  std::vector<int>         *vRefMuP, *vRefMuM, *vRefPi1, *vRefPi2, *vRefPi3;
 
  // ********************************** ************************************************************************

  std::vector<float>       *BDecayVtxX, *BDecayVtxY, *BDecayVtxZ;
  std::vector<double>      *BDecayVtxXE, *BDecayVtxYE, *BDecayVtxZE;
  std::vector<double>      *BDecayVtxXYE, *BDecayVtxXZE, *BDecayVtxYZE;

  // *************************************

  unsigned int             nB;
  unsigned int             nMu;
  unsigned int             nJpsi;
  unsigned int             nPsi2S;
  unsigned int             nJpsi_test;

  std::vector<float>       *B_mass, *B_px, *B_py, *B_pz;

  std::vector<float>       *piPi_mass, *psiPiPi_mass;
  std::vector<float>       *deltaR1, *deltaR2;
  std::vector<float>       *pointingAngle;

  std::vector<float>       *pi1_px, *pi1_py, *pi1_pz,  *pi1_charge;
  std::vector<float>       *pi1_px_track, *pi1_py_track, *pi1_pz_track;
  std::vector<float>       *pi2_px, *pi2_py, *pi2_pz,  *pi2_charge;
  std::vector<float>       *pi2_px_track, *pi2_py_track, *pi2_pz_track;
  std::vector<float>       *pi3_px, *pi3_py, *pi3_pz,  *pi3_charge;
  std::vector<float>       *pi3_px_track, *pi3_py_track, *pi3_pz_track;

  std::vector<float>       *J_mass, *J_px, *J_py, *J_pz;
  std::vector<float>       *J_pt1, *J_px1, *J_py1, *J_pz1;
  std::vector<float>       *J_pt2, *J_px2, *J_py2, *J_pz2;
  std::vector<int>         *J_charge1, *J_charge2;
 
  std::vector<float>       *flightLen, *flightLenErr, *flightLenSig;
/*
  std::vector<float>	   *mu1_px_test, *mu1_py_test, *mu1_pz_test, *mu1_charge_test;
  std::vector<float>       *mu2_px_test, *mu2_py_test, *mu2_pz_test, *mu2_charge_test;
  std::vector<float>       *Jpsi_dca_test;
  std::vector<float>       *Jpsi_vx_test, *Jpsi_vy_test, *Jpsi_vz_test;
 
  std::vector<float>       *Jpsi_mass_test, *Jpsi_prob_test, *Jpsi_chi2_test;
*/
  std::vector<float>       *J_chi2, *psi2S_chi2, *B_chi2; 
  std::vector<float>       *B_Prob, *J_Prob, *psi2S_Prob; 

  int  run, event;
  int  lumiblock;

};

#endif
