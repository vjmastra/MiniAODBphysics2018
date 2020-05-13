#ifndef _Psi2SLambda_h
#define _Psi2SLambda_h

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

class Psi2SLambda : public edm::EDAnalyzer {
public:
  explicit Psi2SLambda(const edm::ParameterSet&);
  ~Psi2SLambda();
  void fillPsi(const reco::Candidate& genpsi);
  void fillV0(const reco::Candidate& genv0);
  int const getMuCat(reco::Muon const& muon) const;
//  bool IsTheSame(const reco::Track& tk, const pat::Muon& mu);
  bool IsTheSame(const pat::PackedCandidate& tk, const pat::Muon& mu);
  bool IsTheSame2(const pat::PackedCandidate& tk, const pat::PackedCandidate& tk2); 
  Double_t Distance(const Double_t p1[], const Double_t p2[]);
  Double_t DistanceError(const Double_t p1[], const Double_t err1[], const Double_t p2[], const Double_t err2[]);

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
  edm::EDGetTokenT<edm::TriggerResults> triggerResults_Label;
  edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_;
  edm::EDGetTokenT<reco::VertexCompositePtrCandidateCollection> v0PtrCollection_;

  std::string genParticles_;
  bool OnlyBest_;
  bool isMC_;
  bool OnlyGen_;
  bool doMC_;

  TTree      *tree_; 
  TTree      *treeTest_;

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

  std::vector<int>         *vRefMuP, *vRefMuM, *vRefPi1, *vRefPi2, *vRefDau1, *vRefDau2;
  // ********************************** ************************************************************************

  std::vector<float>       *lBDecayVtxX, *lBDecayVtxY, *lBDecayVtxZ;
  std::vector<double>      *lBDecayVtxXE, *lBDecayVtxYE, *lBDecayVtxZE;
  std::vector<double>      *lBDecayVtxXYE, *lBDecayVtxXZE, *lBDecayVtxYZE;

  std::vector<float>       *VDecayVtxX, *VDecayVtxY, *VDecayVtxZ;
  std::vector<float>       *VDecayVtxXE, *VDecayVtxYE, *VDecayVtxZE;
  std::vector<float>       *VDecayVtxXYE, *VDecayVtxXZE, *VDecayVtxYZE;

  // *************************************

  unsigned int             nlB;
  unsigned int             nMu;
  unsigned int             nJpsi; 
  unsigned int             nPsi2S;

  std::vector<float>       *lB_mass, *lB_px, *lB_py, *lB_pz;

  std::vector<float>       *piPi_mass, *psiPiPi_mass;
  std::vector<float>       *deltaR1, *deltaR2;
  std::vector<float>       *pointingAngle;

  std::vector<float>       *lambda_mass, *lambda_px, *lambda_py, *lambda_pz;
  std::vector<float>       *lambda_pt1, *lambda_px1, *lambda_py1, *lambda_pz1;
  std::vector<float>       *lambda_pt2, *lambda_px2, *lambda_py2, *lambda_pz2;
  
  std::vector<float>       *lambda_px1_track, *lambda_py1_track, *lambda_pz1_track;
  std::vector<float>       *lambda_px2_track, *lambda_py2_track, *lambda_pz2_track;

  std::vector<float>       *pi1dxy, *pi2dxy, *pi1dz, *pi2dz;
  std::vector<float>       *pi1dxy_e, *pi2dxy_e, *pi1dz_e, *pi2dz_e;
  std::vector<int>         *lambda_charge1, *lambda_charge2;
  
  std::vector<float>       *J_mass, *J_px, *J_py, *J_pz;
  std::vector<float>       *J_pt1, *J_px1, *J_py1, *J_pz1;
  std::vector<float>       *J_pt2, *J_px2, *J_py2, *J_pz2;
  std::vector<int>         *J_charge1, *J_charge2;
  
  std::vector<float>       *flightLen, *flightLenErr, *flightLenSig;

  std::vector<float>       *Jtest_mass, *Jtest_prob;

  std::vector<float>       *lambda_chi2, *J_chi2, *psi2S_chi2, *lB_chi2;
  std::vector<float>       *lB_Prob, *J_Prob, *lambda_Prob, *psi2S_Prob;

  int  run, event;
  int  lumiblock;

};

#endif
