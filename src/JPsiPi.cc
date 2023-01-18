// -*- C++ -*-
//
// Package:    JPsiPi
// Class:      JPsiPi
// 
//=================================================
// original author:  Jhovanny Andres Mejia        |
//         created:  Fryday Sep 23                |
//         <jhovanny.andres.mejia.guisao@cern.ch> | 
//=================================================


// system include files
#include <memory>

// user include files
#include "myAnalyzers/JPsiKsPAT/src/JPsiPi.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"

#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "DataFormats/MuonReco/interface/MuonChamberMatch.h"

#include "DataFormats/Common/interface/RefToBase.h"
#include "DataFormats/Candidate/interface/ShallowCloneCandidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h"
#include "DataFormats/Math/interface/Error.h"
#include "RecoVertex/VertexTools/interface/VertexDistance.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/CLHEP/interface/Migration.h"

#include "RecoVertex/VertexPrimitives/interface/BasicSingleVertexState.h"
#include "RecoVertex/VertexPrimitives/interface/VertexState.h"

#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"

#include "TFile.h"
#include "TTree.h"

#include "TLorentzVector.h"

#include <vector>
#include <utility>
#include <string>
//
// constants, enums and typedefs
//

  typedef math::Error<3>::type CovarianceMatrix;

//
// static data member definitions
//

//
// constructors and destructor
//
JPsiPi::JPsiPi(const edm::ParameterSet& iConfig)
  :
  dimuon_Label(consumes<edm::View<pat::Muon>>(iConfig.getParameter<edm::InputTag>("dimuons"))),
  trakCollection_label(consumes<edm::View<pat::PackedCandidate>>(iConfig.getParameter<edm::InputTag>("Trak"))),
  primaryVertices_Label(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("primaryVertices"))),
  BSLabel_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("bslabel"))),
  builderToken_(esConsumes<TransientTrackBuilder, TransientTrackRecord>(edm::ESInputTag("", "TransientTrackBuilder"))),
  triggerCollection_(consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("TriggerInput"))),
  triggerResults_Label(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("TriggerResults"))),
  genParticles_ ( iConfig.getUntrackedParameter<std::string>("GenParticles",std::string("genParticles")) ),
  OnlyBest_(iConfig.getParameter<bool>("OnlyBest")),
  isMC_(iConfig.getParameter<bool>("isMC")),
  OnlyGen_(iConfig.getParameter<bool>("OnlyGen")),
  doMC_ ( iConfig.getUntrackedParameter<bool>("doMC",false) ),
  tree_(0), 

  mumC2(0), mumNHits(0), mumNPHits(0),
  mupC2(0), mupNHits(0), mupNPHits(0),
  mumdxy(0), mupdxy(0), mumdz(0), mupdz(0),
  muon_dca(0),

  mu1soft(0), mu2soft(0), mu1tight(0), mu2tight(0), 
  mu1PF(0), mu2PF(0), mu1loose(0), mu2loose(0),

  nVtx(0),
  priVtxX(0), priVtxY(0), priVtxZ(0), priVtxXE(0), priVtxYE(0), priVtxZE(0), priVtxCL(0),
  priVtxXYE(0), priVtxXZE(0), priVtxYZE(0),

  indexVtx(0), vRefTrk(0),

  JDecayVtxX(0), JDecayVtxY(0), JDecayVtxZ(0), 
  JDecayVtxXE(0), JDecayVtxYE(0), JDecayVtxZE(0),
  JDecayVtxXYE(0), JDecayVtxXZE(0), JDecayVtxYZE(0),

  BDecayVtxX(0), BDecayVtxY(0), BDecayVtxZ(0), 
  BDecayVtxXE(0), BDecayVtxYE(0), BDecayVtxZE(0),
  BDecayVtxXYE(0), BDecayVtxXZE(0), BDecayVtxYZE(0),

  // *******************************************************
  nB(0), nMu(0), nJpsi(0), 

  triggerFlags({}),
  mu1HLTmatched(0), mu2HLTmatched(0), trkHLT0matched(0), trkHLT1matched(0),

  B_mass(0), B_px(0), B_py(0), B_pz(0),
  B_pt(0), B_eta(0), B_phi(0),

  jPointingAngle(0), bPointingAngle(0),
  jPointingAngle2D(0), bPointingAngle2D(0),

  pi_px(0), pi_py(0), pi_pz(0), pi_charge(0),
  pi_pt(0), pi_eta(0), pi_phi(0),

  d0ValPi(0), d0ErrPi(0), d0SigPi(0),

  J_mass(0), J_px(0), J_py(0), J_pz(0), deltaRmumu(0),
  J_pt(0), J_eta(0), J_phi(0),

  J_px1(0), J_py1(0), J_pz1(0), 
  J_pt1(0), J_eta1(0), J_phi1(0),

  J_px2(0), J_py2(0), J_pz2(0), 
  J_pt2(0), J_eta2(0), J_phi2(0),

  J_charge1(0), J_charge2(0),

  jFlightLen(0), jFlightLenErr(0), jFlightLenSig(0),
  jFlightLen2D(0), jFlightLenErr2D(0), jFlightLenSig2D(0),
  bFlightLen(0), bFlightLenErr(0), bFlightLenSig(0),
  bFlightLen2D(0), bFlightLenErr2D(0), bFlightLenSig2D(0),
  bFlightTime(0), bFlightTimeErr(0),
  deltaR_J_pi(0), cosAngle_J_pi(0),
 
  J_chi2(0), B_chi2(0),
  B_Prob(0), J_Prob(0),

  run(0), event(0),
  lumiblock(0)

{
   //now do what ever initialization is needed
}

JPsiPi::~JPsiPi()
{

}

// ------------ method called to for each event  ------------
void JPsiPi::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using std::vector;
  using namespace edm;
  using namespace reco;
  using namespace std;
  
  // *********************************
  // Get event content information
  // *********************************  

  auto const &theB = iSetup.getData(builderToken_);

  edm::Handle< View<pat::PackedCandidate> > thePATTrackHandle;
  iEvent.getByToken(trakCollection_label,thePATTrackHandle);

  edm::Handle< View<pat::Muon> > thePATMuonHandle;
  iEvent.getByToken(dimuon_Label,thePATMuonHandle);

  edm::Handle<edm::TriggerResults> triggerBits;
  iEvent.getByToken(triggerResults_Label, triggerBits);

  edm::Handle<std::vector<pat::TriggerObjectStandAlone>> triggerCollection;
  iEvent.getByToken(triggerCollection_, triggerCollection);

  //Trigger Collections

  const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);
  const pat::TriggerObjectStandAloneCollection unPackedCollection;

  std::array<std::string, 4> pathNames = {{
"HLT_DoubleMu4_JpsiTrk_Bc_v", "HLT_DoubleMu4_MuMuTrk_Displaced_v", "HLT_DoubleMu4_3_LowMass_v", "HLT_DoubleMu4_LowMass_Displaced_v"}};
  std::string hltMuColl = "hltIterL3MuonCandidates";
  std::string hltTkColl0 = "hltBcJpsiTkAllConeTracksIter"; 
  std::string hltTkColl1 = "hltMuMuTkAllConeTracksIter";

  for (unsigned int i = 0; i  < triggerBits->size(); i++) {
    std::string iName = names.triggerName(i);
    bool flag = triggerBits->accept(i);
    if (flag) {
      for(unsigned int j = 0; j < pathNames.size(); j++) {
        if (iName.find(pathNames[j]) != std::string::npos) {
          triggerFlags[j] = flag;
        }
      }
    }
    //if (triggerBits->accept(i)) std::cout << "Trigger " << names.triggerName(i) << "\n";
  }

//  for (unsigned int j = 0; j < triggerFlags.size(); j++) std::cout << triggerFlags[j] << " ";
//  std::cout << std::endl;

//  if (triggerFlag) std::cout << "Trigger found!" << std::endl;  

  for (pat::TriggerObjectStandAlone trig : *triggerCollection) {
    trig.unpackPathNames(names);
    trig.unpackFilterLabels(iEvent, *triggerBits);  
    //std::cout << "\tTrigger object:  pt " << trig.pt() << ", eta " << trig.eta() << ", phi " << trig.phi() << '\n';
    //std::cout << "\t   Collection: " << trig.collection() << '\n';
    //std::cout << "\t   Type IDs:   ";
    //for (unsigned h = 0; h < trig.filterIds().size(); ++h) std::cout << " " << trig.filterIds()[h] ;
    //std::cout << std::endl;
  }

  // *********************************
  //Now we get the primary vertex 
  // *********************************

  reco::Vertex highestptVtx;

  // get primary vertex
  edm::Handle<std::vector<reco::Vertex> > recVtxs;
  iEvent.getByToken(primaryVertices_Label, recVtxs);

  highestptVtx = *(recVtxs->begin());
  nVtx = recVtxs->size();

  lumiblock = iEvent.id().luminosityBlock();
  run = iEvent.id().run();
  event = iEvent.id().event(); 

  //*****************************************
  //Let's begin by looking for J/psi->mu+mu-

  unsigned int nMu_tmp = thePATMuonHandle->size();
 
  for(View<pat::Muon>::const_iterator iMuon1 = thePATMuonHandle->begin(); iMuon1 != thePATMuonHandle->end(); ++iMuon1) { 

    if((iMuon1->track()).isNull()) continue;
    if(iMuon1->track()->pt()<4.0) continue;
 
    for(View<pat::Muon>::const_iterator iMuon2 = iMuon1+1; iMuon2 != thePATMuonHandle->end(); ++iMuon2) {  

      if((iMuon2->track()).isNull()) continue;      
      if(iMuon1==iMuon2) continue;
      if( (iMuon1->charge())*(iMuon2->charge()) == 1) continue;
      if(iMuon2->track()->pt()<4.0) continue;

      TrackRef glbTrackP;	  
      TrackRef glbTrackM;	  
  
      if(iMuon1->charge() == 1){glbTrackP = iMuon1->track();}
      if(iMuon1->charge() == -1){glbTrackM = iMuon1->track();}
  
      if(iMuon2->charge() == 1) {glbTrackP = iMuon2->track();}
      if(iMuon2->charge() == -1){glbTrackM = iMuon2->track();}

      if(!(glbTrackM->quality(reco::TrackBase::highPurity))) continue;
      if(!(glbTrackP->quality(reco::TrackBase::highPurity))) continue;
	  
      //Let's check the vertex and mass
      reco::TransientTrack muon1TT((theB).build(glbTrackP));
      reco::TransientTrack muon2TT((theB).build(glbTrackM));

      // *****  Trajectory states to calculate DCA for the 2 muons *********************
      FreeTrajectoryState mu1State = muon1TT.impactPointTSCP().theState();
      FreeTrajectoryState mu2State = muon2TT.impactPointTSCP().theState();
 
      if( !muon1TT.impactPointTSCP().isValid() || !muon2TT.impactPointTSCP().isValid() ) continue;

      // Measure distance between tracks at their closest approach
      ClosestApproachInRPhi cApp;
      cApp.calculate(mu1State, mu2State);
      if( !cApp.status() ) continue;
      muon_dca = fabs( cApp.distance() );	  
      if (muon_dca < 0. || muon_dca > 0.5) continue;
      //cout<<" closest approach  "<<dca<<endl;

      // *****  end DCA for the 2 muons *********************

      //The mass of a muon and the insignificant mass sigma 
      //to avoid singularities in the covariance matrix.
      ParticleMass muon_mass = 0.10565837; //pdg mass
      ParticleMass Jpsi_mass = 3.096916;
      float muon_sigma = muon_mass*1.e-6;

      //Creating a KinematicParticleFactory
      KinematicParticleFactoryFromTransientTrack pFactory;
  
      //initial chi2 and ndf before kinematic fits.
      float chi = 0.;
      float ndf = 0.;
      vector<RefCountedKinematicParticle> muonParticles;
      try {
        muonParticles.push_back(pFactory.particle(muon1TT,muon_mass,chi,ndf,muon_sigma));
        muonParticles.push_back(pFactory.particle(muon2TT,muon_mass,chi,ndf,muon_sigma));
      }
      catch(...) { 
        std::cout<<" Exception caught ... continuing 1 "<<std::endl; 
        continue;
      }

      KinematicParticleVertexFitter fitter;   

      RefCountedKinematicTree JVertexFitTree;
      try {
        JVertexFitTree = fitter.fit(muonParticles); 
      }
      catch (...) { 
        std::cout<<" Exception caught ... continuing 2 "<<std::endl; 
        continue;
      }

      if (!JVertexFitTree->isValid()) {
        //std::cout << "caught an exception in the psi vertex fit" << std::endl;
        continue; 
      }

      JVertexFitTree->movePointerToTheTop();
	  
      RefCountedKinematicParticle J_vFit_noMC = JVertexFitTree->currentParticle(); 
      RefCountedKinematicVertex J_vFit_vertex_noMC = JVertexFitTree->currentDecayVertex();
	  
      if(J_vFit_vertex_noMC->chiSquared() < 0 ) {
        std::cout << "negative chisq from psi fit" << endl;
        continue;
      }

      J_Prob = TMath::Prob(J_vFit_vertex_noMC->chiSquared(),(int)J_vFit_vertex_noMC->degreesOfFreedom());
      if(J_Prob<0.01) {
        continue;
      }

      if(J_vFit_noMC->currentState().mass()<2.9 || J_vFit_noMC->currentState().mass()>3.3) continue;

      TLorentzVector muon14V, muon24V;
      muon14V.SetXYZM(iMuon1->track()->px(), iMuon1->track()->py(), iMuon1->track()->pz(), muon_mass);
      muon24V.SetXYZM(iMuon2->track()->px(), iMuon2->track()->py(), iMuon2->track()->pz(), muon_mass);
      deltaRmumu = muon14V.DeltaR(muon24V);
 
      //jpsipi
      for (View<pat::PackedCandidate>::const_iterator iTrack = thePATTrackHandle->begin(); iTrack != thePATTrackHandle->end(); ++iTrack ) {

        if (!(iTrack->hasTrackDetails())) continue;
        if (iTrack->pt()<0.5) continue; //min value 0.5 for 2017-2018, 0.95 for 2015-2016
        if (!(iTrack->trackHighPurity())) continue;

        if ( IsTheSame(*iTrack,*iMuon1) || IsTheSame(*iTrack,*iMuon2) ) continue;

        reco::TransientTrack pionTT((theB).build(iTrack->pseudoTrack()));

        ParticleMass pion_mass = 0.13957018;
        //pion_mass = 0.493677; ///////////////////////////////kaon
        float pion_sigma = pion_mass*1.e-6;
        float chi = 0.;
        float ndf = 0.;

        // JPsi mass constraint is applied

        vector<RefCountedKinematicParticle> vFitMCParticles;
        vFitMCParticles.push_back(pFactory.particle(muon1TT,muon_mass,chi,ndf,muon_sigma));
        vFitMCParticles.push_back(pFactory.particle(muon2TT,muon_mass,chi,ndf,muon_sigma));
        vFitMCParticles.push_back(pFactory.particle(pionTT,pion_mass,chi,ndf,pion_sigma));

        MultiTrackKinematicConstraint *  Jpsi_c = new  TwoTrackMassKinematicConstraint(Jpsi_mass);
        KinematicConstrainedVertexFitter kcvFitter;
        RefCountedKinematicTree bVertexFitTree = kcvFitter.fit(vFitMCParticles, Jpsi_c);
        if (!bVertexFitTree->isValid()) continue;

        bVertexFitTree->movePointerToTheTop();
        RefCountedKinematicParticle bCandMC = bVertexFitTree->currentParticle();          
        RefCountedKinematicVertex bDecayVertexMC = bVertexFitTree->currentDecayVertex();
        if (!bDecayVertexMC->vertexIsValid()) continue;
        if ((bCandMC->currentState().mass() < 5.95) || (bCandMC->currentState().mass() > 6.55) ) continue; //5.95, 6.55 per bc 

        if(bDecayVertexMC->chiSquared()<0) {
           //std::cout << " continue from negative chi2 = " << bDecayVertexMC->chiSquared() << endl;
           continue;
        }

        B_Prob = TMath::Prob(bDecayVertexMC->chiSquared(),(int)bDecayVertexMC->degreesOfFreedom());
        if (B_Prob<0.01) continue;

        //Select PV that minimized the pointing angle

        reco::Vertex bestVtxBSIP;
        reco::Vertex vtxBS;

        jPointingAngle = -10.;
        bPointingAngle = -10.;
        jPointingAngle2D = -10.;
        bPointingAngle2D = -10.;

        Double_t jVtx[3] = {(*J_vFit_vertex_noMC).position().x(), (*J_vFit_vertex_noMC).position().y(), (*J_vFit_vertex_noMC).position().z()};
        Double_t sVtx[3] = {(*bDecayVertexMC).position().x(), (*bDecayVertexMC).position().y(), (*bDecayVertexMC).position().z()};
        TVector3 Jmomentum(J_vFit_noMC->currentState().globalMomentum().x(), J_vFit_noMC->currentState().globalMomentum().y(), J_vFit_noMC->currentState().globalMomentum().z());
        TVector3 Bmomentum(bCandMC->currentState().globalMomentum().x(), bCandMC->currentState().globalMomentum().y(), bCandMC->currentState().globalMomentum().z());

        Double_t jVtxT[3] = {jVtx[0], jVtx[1], 0.};
        Double_t sVtxT[3] = {sVtx[0], sVtx[1], 0.};
        TVector3 JmomentumT(J_vFit_noMC->currentState().globalMomentum().x(), J_vFit_noMC->currentState().globalMomentum().y(), 0.);
        TVector3 BmomentumT(bCandMC->currentState().globalMomentum().x(), bCandMC->currentState().globalMomentum().y(), 0.);

        for (size_t i = 0; i < recVtxs->size(); i++) {
          vtxBS = (*recVtxs)[i];
          //pointing Angle computation
          Double_t primaryVertex[3] = {vtxBS.x(), vtxBS.y(), vtxBS.z()};
          Double_t jFlightVec[3];
          Double_t bFlightVec[3];
          for (int i = 0; i < 3; i++) {
            jFlightVec[i] = jVtx[i] - primaryVertex[i];
            bFlightVec[i] = sVtx[i] - primaryVertex[i];
          }
          TVector3 jFlightDir(jFlightVec[0], jFlightVec[1], jFlightVec[2]);
          TVector3 bFlightDir(bFlightVec[0], bFlightVec[1], bFlightVec[2]);
          TVector3 jFlightDirT(jFlightVec[0], jFlightVec[1], 0.);
          TVector3 bFlightDirT(bFlightVec[0], bFlightVec[1], 0.);
          //best PV selection
          double cosAlphaXYb = TMath::Cos(bFlightDir.Angle(Bmomentum));
          if (cosAlphaXYb > bPointingAngle) {
            bPointingAngle = cosAlphaXYb;
            jPointingAngle = TMath::Cos(jFlightDir.Angle(Jmomentum));
            bPointingAngle2D = TMath::Cos(bFlightDirT.Angle(BmomentumT));
            jPointingAngle2D = TMath::Cos(jFlightDirT.Angle(JmomentumT));
            indexVtx = i;
            bestVtxBSIP = vtxBS;
          }
        }

        //if (bPointingAngle < 0.9) continue; //Cut from JpsiTk trigger            

        vRefTrk = (int)iTrack->vertexRef().key();

        priVtxX = bestVtxBSIP.x();
        priVtxY = bestVtxBSIP.y();
        priVtxZ = bestVtxBSIP.z();
        priVtxXE = bestVtxBSIP.covariance(0, 0);
        priVtxYE = bestVtxBSIP.covariance(1, 1);
        priVtxZE = bestVtxBSIP.covariance(2, 2);
        priVtxXYE = bestVtxBSIP.covariance(0, 1);
        priVtxXZE = bestVtxBSIP.covariance(0, 2);
        priVtxYZE = bestVtxBSIP.covariance(1, 2);
        priVtxCL = (TMath::Prob(bestVtxBSIP.chi2(), (int)bestVtxBSIP.ndof()));

        //Flight distance
            
        Double_t pVtx[3] = {bestVtxBSIP.x(), bestVtxBSIP.y(), bestVtxBSIP.z()};
        Double_t pVtxCov[6] = {bestVtxBSIP.covariance(0, 0), bestVtxBSIP.covariance(1, 1), bestVtxBSIP.covariance(2, 2), bestVtxBSIP.covariance(0, 1), bestVtxBSIP.covariance(0, 2), bestVtxBSIP.covariance(1, 2)}; //xx yy zz xy xz yz
        Double_t jVtxCov[6] = {J_vFit_vertex_noMC->error().cxx(), J_vFit_vertex_noMC->error().cyy(), J_vFit_vertex_noMC->error().czz(), J_vFit_vertex_noMC->error().cyx(), J_vFit_vertex_noMC->error().czx(), J_vFit_vertex_noMC->error().czy()}; //xx yy zz xy xz yz
        Double_t sVtxCov[6] = {bDecayVertexMC->error().cxx(), bDecayVertexMC->error().cyy(), bDecayVertexMC->error().czz(), bDecayVertexMC->error().cyx(), bDecayVertexMC->error().czx(), bDecayVertexMC->error().czy()}; //xx yy zz xy xz yz

        Double_t pVtxT[3] = {pVtx[0], pVtx[1], 0.};

        jFlightLen = Distance(pVtx, jVtx);
        jFlightLenErr = DistanceError(pVtx, pVtxCov, jVtx, jVtxCov);
        jFlightLenSig = jFlightLen/jFlightLenErr;

        bFlightLen = Distance(pVtx, sVtx);
        bFlightLenErr = DistanceError(pVtx, pVtxCov, sVtx, sVtxCov);
        bFlightLenSig = bFlightLen/bFlightLenErr;

        jFlightLen2D = Distance(pVtxT, jVtxT);
        jFlightLenErr2D = DistanceError(pVtxT, pVtxCov, jVtxT, jVtxCov);
        jFlightLenSig2D = jFlightLen2D/jFlightLenErr2D;

        bFlightLen2D = Distance(pVtxT, sVtxT);
        bFlightLenErr2D = DistanceError(pVtxT, pVtxCov, sVtxT, sVtxCov);
        bFlightLenSig2D = bFlightLen2D/bFlightLenErr2D;

        pi_pt = iTrack->pt();
        pi_eta = iTrack->eta();
        pi_phi = iTrack->phi();

        d0ValPi = iTrack->dxy();
        d0ErrPi = iTrack->dxyError();
        d0SigPi = (d0ErrPi == 0) ? 0 : fabs(d0ValPi/d0ErrPi);

        TLorentzVector jpsi4V, pi4V;
        jpsi4V.SetXYZM(J_vFit_noMC->currentState().globalMomentum().x(), J_vFit_noMC->currentState().globalMomentum().y(), J_vFit_noMC->currentState().globalMomentum().z(), Jpsi_mass);
        pi4V.SetXYZM(iTrack->px(), iTrack->py(), iTrack->pz(), pion_mass);
        deltaR_J_pi = jpsi4V.DeltaR(pi4V);
        cosAngle_J_pi = TMath::Cos(jpsi4V.Angle(pi4V.Vect()));

        //Trigger matching

        int sum = 0;
        for (unsigned int i = 0; i < triggerFlags.size(); i++) sum += triggerFlags[i];

        if (sum) {
          for (pat::TriggerObjectStandAlone obj : *triggerCollection) {
            if(MatchByDRDPt(*iMuon1, obj)) mu1HLTmatched = true;
            if(MatchByDRDPt(*iMuon2, obj)) mu2HLTmatched = true;
            if(obj.collection().find(hltTkColl0) != std::string::npos && MatchByDRDPt(*iTrack, obj)) trkHLT0matched = true;
            if(obj.collection().find(hltTkColl1) != std::string::npos && MatchByDRDPt(*iTrack, obj)) trkHLT1matched = true;
          }
        }
    
        //std::cout << mu1HLTmatched << " " << mu2HLTmatched << " " << trkHLT0matched << " " << trkHLT1matched << std::endl;

        // get children from final B fit
	bVertexFitTree->movePointerToTheFirstChild();
	RefCountedKinematicParticle mu1CandMC = bVertexFitTree->currentParticle();
	bVertexFitTree->movePointerToTheNextChild();
	RefCountedKinematicParticle mu2CandMC = bVertexFitTree->currentParticle();
        bVertexFitTree->movePointerToTheNextChild();
        RefCountedKinematicParticle piCandMC = bVertexFitTree->currentParticle();
		   
	KinematicParameters psiMu1KP = mu1CandMC->currentState().kinematicParameters();
	KinematicParameters psiMu2KP = mu2CandMC->currentState().kinematicParameters();
	KinematicParameters psiMupKP;
	KinematicParameters psiMumKP;
	       
	if ( mu1CandMC->currentState().particleCharge() > 0 ) psiMupKP = psiMu1KP;
	if ( mu1CandMC->currentState().particleCharge() < 0 ) psiMumKP = psiMu1KP;
	if ( mu2CandMC->currentState().particleCharge() > 0 ) psiMupKP = psiMu2KP;
        if ( mu2CandMC->currentState().particleCharge() < 0 ) psiMumKP = psiMu2KP;	 

        KinematicParameters pionKP = piCandMC->currentState().kinematicParameters();

 	GlobalVector Jp1vec(mu1CandMC->currentState().globalMomentum().x(),
				mu1CandMC->currentState().globalMomentum().y(),
				mu1CandMC->currentState().globalMomentum().z());

        GlobalVector Jp2vec(mu2CandMC->currentState().globalMomentum().x(),
				mu2CandMC->currentState().globalMomentum().y(),
 	      			mu2CandMC->currentState().globalMomentum().z());

 	GlobalVector pivec( piCandMC->currentState().globalMomentum().x(),
	      	  		piCandMC->currentState().globalMomentum().y(),
 	  			piCandMC->currentState().globalMomentum().z());

        // fill candidate variables now
		   
        if(nB==0){		    
	  nMu  = nMu_tmp;
	  // cout<< "*Number of Muons : " << nMu_tmp << endl;
	}		     

	B_mass = bCandMC->currentState().mass();
	B_px = bCandMC->currentState().globalMomentum().x();
	B_py = bCandMC->currentState().globalMomentum().y();
	B_pz = bCandMC->currentState().globalMomentum().z();
        B_pt = bCandMC->currentState().globalMomentum().perp();
        B_eta = bCandMC->currentState().globalMomentum().eta();
        B_phi = bCandMC->currentState().globalMomentum().phi();

        J_mass = J_vFit_noMC->currentState().mass();
        J_px = J_vFit_noMC->currentState().globalMomentum().x();
        J_py = J_vFit_noMC->currentState().globalMomentum().y();
        J_pz = J_vFit_noMC->currentState().globalMomentum().z();
        J_pt = J_vFit_noMC->currentState().globalMomentum().perp();
        J_eta = J_vFit_noMC->currentState().globalMomentum().eta();
        J_phi = J_vFit_noMC->currentState().globalMomentum().phi();

        J_px1 = psiMu1KP.momentum().x();
	J_py1 = psiMu1KP.momentum().y();
	J_pz1 = psiMu1KP.momentum().z();
        J_pt1 = Jp1vec.perp();
        J_eta1 = Jp1vec.eta();
        J_phi1 = Jp1vec.phi();
	J_charge1 = mu1CandMC->currentState().particleCharge();
 
	J_px2 = psiMu2KP.momentum().x();
	J_py2 = psiMu2KP.momentum().y();
	J_pz2 = psiMu2KP.momentum().z();
        J_pt2 = Jp2vec.perp();
        J_eta2 = Jp2vec.eta();
        J_phi2 = Jp2vec.phi();
	J_charge2 = mu2CandMC->currentState().particleCharge();

        pi_px = pionKP.momentum().x();
        pi_py = pionKP.momentum().y();
        pi_pz = pionKP.momentum().z();
        pi_charge = piCandMC->currentState().particleCharge();

	J_chi2 = J_vFit_vertex_noMC->chiSquared();
	B_chi2 = bDecayVertexMC->chiSquared();

        JDecayVtxX = (*J_vFit_vertex_noMC).position().x();
        JDecayVtxY = (*J_vFit_vertex_noMC).position().y();
        JDecayVtxZ = (*J_vFit_vertex_noMC).position().z();
        JDecayVtxXE = J_vFit_vertex_noMC->error().cxx();
        JDecayVtxYE = J_vFit_vertex_noMC->error().cyy();
        JDecayVtxZE = J_vFit_vertex_noMC->error().czz();
        JDecayVtxXYE = J_vFit_vertex_noMC->error().cyx();
        JDecayVtxXZE = J_vFit_vertex_noMC->error().czx();
        JDecayVtxYZE = J_vFit_vertex_noMC->error().czy();

	BDecayVtxX = (*bDecayVertexMC).position().x();
	BDecayVtxY = (*bDecayVertexMC).position().y();
	BDecayVtxZ = (*bDecayVertexMC).position().z();
	BDecayVtxXE = bDecayVertexMC->error().cxx();
	BDecayVtxYE = bDecayVertexMC->error().cyy();
	BDecayVtxZE = bDecayVertexMC->error().czz();
	BDecayVtxXYE = bDecayVertexMC->error().cyx();
	BDecayVtxXZE = bDecayVertexMC->error().czx();
	BDecayVtxYZE = bDecayVertexMC->error().czy();
		  
	mu1soft = iMuon1->isSoftMuon(bestVtxBSIP);
	mu2soft = iMuon2->isSoftMuon(bestVtxBSIP);
	mu1tight = iMuon1->isTightMuon(bestVtxBSIP);
	mu2tight = iMuon2->isTightMuon(bestVtxBSIP);
	mu1PF = iMuon1->isPFMuon();
	mu2PF = iMuon2->isPFMuon();
	mu1loose = muon::isLooseMuon(*iMuon1);
	mu2loose = muon::isLooseMuon(*iMuon2);

        mumC2 = glbTrackM->normalizedChi2();
	mumNHits = glbTrackM->numberOfValidHits();
	mumNPHits = glbTrackM->hitPattern().numberOfValidPixelHits();	       
	mupC2 = glbTrackP->normalizedChi2();
	mupNHits = glbTrackP->numberOfValidHits();
	mupNPHits = glbTrackP->hitPattern().numberOfValidPixelHits();
        mumdxy = glbTrackM->dxy(bestVtxBSIP.position());
	mupdxy = glbTrackP->dxy(bestVtxBSIP.position());
	mumdz = glbTrackM->dz(bestVtxBSIP.position());
	mupdz = glbTrackP->dz(bestVtxBSIP.position());
	
        //Proper time calculation
        
        double len[2] = {sVtx[0]-pVtx[0], sVtx[1]-pVtx[1]};

        double factor = (len[0]*len[1])/(len[0]*len[0] + len[1]*len[1]);
        double elen[3] = {	bestVtxBSIP.covariance(0, 0) + bDecayVertexMC->error().cxx(), 
				bestVtxBSIP.covariance(1, 1) + bDecayVertexMC->error().cyy(), 
				factor*(bestVtxBSIP.covariance(0, 1) + bDecayVertexMC->error().cyx())};

        double mom[2] = {B_px, B_py};
        double emom[3] = {	bCandMC->currentState().kinematicParametersError().matrix()(3, 3), 
				bCandMC->currentState().kinematicParametersError().matrix()(4, 4), 
				bCandMC->currentState().kinematicParametersError().matrix()(3, 4)};

        
        double tDecLen = ProperTime(len, mom);
        double tDecLenErr = ProperTimeErr(len, elen, mom, emom);

        double a = 0.;
        double b = 0.;
        double c = 0.;

        a = bCandMC->currentState().kinematicParametersError().matrix()(6, 6);
        a /= B_mass;

        abs(tDecLen) > 0 ? b = tDecLenErr/tDecLen : b = 0.;

        double mom2[3] = {mom[0]*mom[0], mom[1]*mom[1], fabs(mom[0]*mom[1])};
        for (int i = 0; i < 3; i++) c += mom2[i]*emom[i];
        B_pt > 0 ? c = TMath::Sqrt(c)/B_pt : c = 0.;
      
        B_pt > 0 ? bFlightTime = tDecLen * B_mass / B_pt : bFlightTime = -10.;
        const double lightspeed = 2.99792458;
        bFlightTime *= (100./lightspeed);
        bFlightTimeErr = bFlightTime * TMath::Sqrt( a*a + b*b + c*c);

        //std::cout << bFlightTime << " " << bFlightTimeErr << " " << bFlightTimeErr/bFlightTime << std::endl;

        //Fill and count
  
        tree_->Fill();    
        nB++;

        //Clear
        
        mu1HLTmatched = false; mu2HLTmatched = false; trkHLT0matched = false; trkHLT1matched = false;

        B_mass = 0.; B_px = 0.; B_py = 0.; B_pz = 0.;
        B_pt = 0.; B_eta = 0.; B_phi = 0.;

        J_mass = 0.; J_px = 0.; J_py = 0.; J_pz = 0.;
        J_pt = 0.; J_eta = 0.; J_phi = 0.;

        J_px1 = 0.; J_py1 = 0.; J_pz1 = 0.; J_charge1 = 0;
        J_pt1 = 0.; J_eta1 = 0.; J_phi1 = 0.;
        J_px2 = 0.; J_py2 = 0.; J_pz2 = 0.; J_charge2 = 0;
        J_pt2 = 0.; J_eta2 = 0.; J_phi2 = 0.;
        pi_px = 0.; pi_py = 0.; pi_pz = 0.; pi_charge = 0;      
        pi_pt = 0.; pi_eta = 0.; pi_phi = 0.;

        d0ValPi = 0.; d0ErrPi = 0.; d0SigPi = 0.;
        B_chi2 = 0.; B_Prob = 0.;

        priVtxX = 0.; priVtxY = 0.; priVtxZ = 0.;
        priVtxXE = 0.; priVtxYE = 0.; priVtxZE = 0.;
        priVtxXYE = 0.; priVtxXZE = 0.; priVtxYZE = 0.;
        priVtxCL = 0;

        indexVtx = 0; vRefTrk = 0;

        jFlightLen = 0.; jFlightLenErr = 0.; jFlightLenSig = 0.;
        jFlightLen2D = 0.; jFlightLenErr2D = 0.; jFlightLenSig2D = 0.;
        bFlightLen = 0.; bFlightLenErr = 0.; bFlightLenSig = 0.;
        bFlightLen2D = 0.; bFlightLenErr2D = 0.; bFlightLenSig2D = 0.;
        jPointingAngle = 0.; bPointingAngle = 0.;
        jPointingAngle2D = 0.; bPointingAngle2D = 0.;
        deltaR_J_pi = 0.; cosAngle_J_pi = 0.;

        bFlightTime = 0.; bFlightTimeErr = 0.;

        JDecayVtxX = 0.; JDecayVtxY = 0.; JDecayVtxZ = 0.;
        JDecayVtxXE = 0.; JDecayVtxYE = 0.; JDecayVtxZE = 0.;
        JDecayVtxXYE = 0.; JDecayVtxXZE = 0.; JDecayVtxYZE = 0.;

        BDecayVtxX = 0.; BDecayVtxY = 0.; BDecayVtxZ = 0.;
        BDecayVtxXE = 0.; BDecayVtxYE = 0.; BDecayVtxZE = 0.;
        BDecayVtxXYE = 0.; BDecayVtxXZE = 0.; BDecayVtxYZE = 0.;

        mumC2 = 0; mumNHits = 0; mumNPHits = 0;
        mumdxy = 0.; mumdz = 0.;
        mupC2 = 0; mupNHits = 0; mupNPHits = 0; 
        mupdxy = 0.; mupdz = 0.;

        mu1soft = 0; mu1PF = 0;
        mu1loose = 0; mu1tight = 0;
        mu2soft = 0; mu2PF = 0;
        mu2loose = 0; mu2tight = 0;
 
        muonParticles.clear();
	vFitMCParticles.clear();
	   
      }//pi

      nJpsi++;
      nB = 0;
 
      J_chi2 = 0.; J_Prob = 0.;
      deltaRmumu = 0.;
      muon_dca = 0.;

    }//muon2
  }//muon1
    
  nB = 0; nMu = 0; nJpsi = 0;
  triggerFlags = {{0, 0, 0, 0}};

}//end analyzer

bool JPsiPi::IsTheSame(const pat::PackedCandidate& tk, const pat::Muon& mu){
  double DeltaEta = fabs(mu.eta()-tk.eta());
  double DeltaP   = fabs(mu.p()-tk.p());
  if (DeltaP > double(M_PI)) {
    DeltaP -= double(2*M_PI);
    DeltaP  = fabs(DeltaP);
  }
  if (DeltaEta < 0.02 && DeltaP < 0.02) return true;
  return false;
}

bool JPsiPi::IsTheSame2(const pat::PackedCandidate& tk1, const pat::PackedCandidate& tk2){
  double DeltaEta = fabs(tk1.eta()-tk2.eta());
  double DeltaP   = fabs(tk1.p()-tk2.p());
  if (DeltaP > double(M_PI)) {
    DeltaP -= double(2*M_PI);
    DeltaP  = fabs(DeltaP);
  }
  if (DeltaEta < 0.02 && DeltaP < 0.02) return true;
  return false;
}

Double_t JPsiPi::Distance(const Double_t p1[], const Double_t p2[]){
  Double_t diff[3], diff2[3];
  Double_t dist = 0;
  for (int i = 0; i < 3; i++) {
    diff[i] = p2[i]-p1[i];
    diff2[i] = diff[i]*diff[i];
    dist += diff2[i];
  }
  return TMath::Sqrt(dist);
}

Double_t JPsiPi::DistanceError(const Double_t p1[], const Double_t err1[], const Double_t p2[], const Double_t err2[]){
  Double_t diff[3];
  for (int i = 0; i < 3; i++) {
    diff[i] = fabs(p2[i]-p1[i]);
  }
  Double_t deriv2[6] = {diff[0]*diff[0], diff[1]*diff[1], diff[2]*diff[2], diff[0]*diff[1], diff[0]*diff[2], diff[1]*diff[2]};
  for (int i = 0; i < 6; i++) {
    deriv2[i] = fabs(deriv2[i]);
  }
  Double_t covar[6];
  for (int i = 0; i < 6; i++) {
    covar[i] = err1[i] + err2[i];
  }
  Double_t dist = 0;
  Double_t distErr = 0;
  for (int i = 0; i < 6; i++) {
    distErr += deriv2[i]*covar[i];
  }
  distErr = TMath::Sqrt(distErr);
  dist = Distance(p1, p2);
  dist == 0 ? distErr = -1 : distErr = distErr/dist;
  return distErr;
}

double JPsiPi::ProperTime(const double len[], const double mom[]){
  double den = TMath::Sqrt(mom[0]*mom[0] + mom[1]*mom[1]);
  double ret = len[0]*mom[0] + len[1]*mom[1];
  if (den > 0) ret /= den;
    else ret = 0.;
  return ret;
}

double JPsiPi::ProperTimeErr(const double len[], const double errlen[], const double mom[], const double errmom[]){
  double pt2 = mom[0]*mom[0] + mom[1]*mom[1];
  if (!(pt2 > 0)) return 0.;

  double mom2[3] = {mom[0]*mom[0], mom[1]*mom[1], fabs(mom[0]*mom[1])};
  double term1 = 0.;
  for (int i = 0; i < 3; i++) term1 += mom2[i]*errlen[i];
  term1 /= pt2;

  double len2[3] = {len[0]*len[0], len[1]*len[1], fabs(len[0]*len[1])};
  double mom4[3] = {mom2[1]*mom2[1], mom2[0]*mom2[0], mom2[0]*mom2[1]};
  double term2 = 0.;
  for (int i = 0; i < 3; i++) term2 += mom4[i]*len2[i]*errmom[i];
  term2 /= pt2*pt2*pt2;

  return TMath::Sqrt(term1 + term2);
}

float JPsiPi::DeltaR(const pat::PackedCandidate t1, const pat::TriggerObjectStandAlone t2){
  float p1 = t1.phi();
  float p2 = t2.phi();
  float e1 = t1.eta();
  float e2 = t2.eta();
  float de = e1-e2;
  auto dp=std::abs(p1-p2);
  if (dp>float(M_PI)) dp-=float(2*M_PI);
  return sqrt(de*de + dp*dp);
}

float JPsiPi::DeltaPt(const pat::PackedCandidate t1, const pat::TriggerObjectStandAlone t2){
  return (fabs(t1.pt()-t2.pt())/t2.pt());
}

bool JPsiPi::MatchByDRDPt(const pat::PackedCandidate t1, const pat::TriggerObjectStandAlone t2){
  bool ptFlag = DeltaPt(t1, t2) < 2.0;
  bool dRFlag = DeltaR(t1, t2)  < 0.01; //from Adriano's DiMuonDiTrackProducer.cc code
  return ptFlag && dRFlag;
}

float JPsiPi::DeltaR(const pat::Muon t1, const pat::TriggerObjectStandAlone t2){
  float p1 = t1.phi();
  float p2 = t2.phi();
  float e1 = t1.eta();
  float e2 = t2.eta();
  float de = e1-e2;
  auto dp=std::abs(p1-p2);
  if (dp>float(M_PI)) dp-=float(2*M_PI);
  return sqrt(de*de + dp*dp);
}

float JPsiPi::DeltaPt(const pat::Muon t1, const pat::TriggerObjectStandAlone t2){
  return (fabs(t1.pt()-t2.pt())/t2.pt());
}

bool JPsiPi::MatchByDRDPt(const pat::Muon t1, const pat::TriggerObjectStandAlone t2){
  bool ptFlag = DeltaPt(t1, t2) < 2.0;
  bool dRFlag = DeltaR(t1, t2)  < 0.01; //from Adriano's DiMuonDiTrackProducer.cc code
  return ptFlag && dRFlag;
}

// ------------ method called once each job just before starting event loop  ------------

void JPsiPi::beginJob()
{

  std::cout << "Beginning analyzer job with value of doMC = " << doMC_ << std::endl;

  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>("ntuple","Bc->JpsiPi ntuple");

  tree_->Branch("idB",&nB,"nB/i");
  tree_->Branch("nMu",&nMu,"nMu/i");
  tree_->Branch("idJpsi",&nJpsi,"nJpsi/i");

  tree_->Branch("triggerFlags", &triggerFlags);
  tree_->Branch("mu1HLTmatched", &mu1HLTmatched);
  tree_->Branch("mu2HLTmatched", &mu2HLTmatched);
  tree_->Branch("trkHLT0matched", &trkHLT0matched);
  tree_->Branch("trkHLT1matched", &trkHLT1matched);

  tree_->Branch("B_mass", &B_mass);
  tree_->Branch("B_px", &B_px);
  tree_->Branch("B_py", &B_py);
  tree_->Branch("B_pz", &B_pz);
  tree_->Branch("B_pt", &B_pt);
  tree_->Branch("B_eta", &B_eta);
  tree_->Branch("B_phi", &B_phi);

  tree_->Branch("jCosAlpha", &jPointingAngle);
  tree_->Branch("bCosAlpha", &bPointingAngle);
  tree_->Branch("jCosAlpha2D", &jPointingAngle2D);
  tree_->Branch("bCosAlpha2D", &bPointingAngle2D);
  tree_->Branch("dR_J_pi", &deltaR_J_pi);
  tree_->Branch("cosAngle_J_pi", &cosAngle_J_pi);
 
  tree_->Branch("J_mass", &J_mass);
  tree_->Branch("J_px", &J_px);
  tree_->Branch("J_py", &J_py);
  tree_->Branch("J_pz", &J_pz);
  tree_->Branch("J_pt", &J_pt);
  tree_->Branch("J_eta", &J_eta);
  tree_->Branch("J_phi", &J_phi);
  tree_->Branch("deltaRmumu", &deltaRmumu);

  tree_->Branch("J_px1", &J_px1);
  tree_->Branch("J_py1", &J_py1);
  tree_->Branch("J_pz1", &J_pz1);
  tree_->Branch("J_pt1", &J_pt1);
  tree_->Branch("J_eta1", &J_eta1);
  tree_->Branch("J_phi1", &J_phi1);
  tree_->Branch("J_charge1", &J_charge1);

  tree_->Branch("J_px2", &J_px2);
  tree_->Branch("J_py2", &J_py2);
  tree_->Branch("J_pz2", &J_pz2);
  tree_->Branch("J_pt2", &J_pt2);
  tree_->Branch("J_eta2", &J_eta2);
  tree_->Branch("J_phi2", &J_phi2);
  tree_->Branch("J_charge2", &J_charge2);

  tree_->Branch("pi_pt", &pi_pt);
  tree_->Branch("pi_eta", &pi_eta);
  tree_->Branch("pi_phi", &pi_phi);
  tree_->Branch("pi_px", &pi_px);
  tree_->Branch("pi_py", &pi_py);
  tree_->Branch("pi_pz", &pi_pz);
  tree_->Branch("pi_charge", &pi_charge);

  tree_->Branch("d0ValPi", &d0ValPi);
  tree_->Branch("d0ErrPi", &d0ErrPi);
  tree_->Branch("d0SigPi", &d0SigPi);

  tree_->Branch("B_chi2", &B_chi2);
  tree_->Branch("J_chi2", &J_chi2);

  tree_->Branch("B_Prob",  &B_Prob);
  tree_->Branch("J_Prob",  &J_Prob);     

  // *************************

  tree_->Branch("priVtxX",&priVtxX);
  tree_->Branch("priVtxY",&priVtxY);
  tree_->Branch("priVtxZ",&priVtxZ);
  tree_->Branch("priVtxXE",&priVtxXE);
  tree_->Branch("priVtxYE",&priVtxYE);
  tree_->Branch("priVtxZE",&priVtxZE);
  tree_->Branch("priVtxXYE",&priVtxXYE);
  tree_->Branch("priVtxXZE",&priVtxXZE);
  tree_->Branch("priVtxYZE",&priVtxYZE);
  tree_->Branch("priVtxCL",&priVtxCL);

  tree_->Branch("indexVtx", &indexVtx);
  tree_->Branch("vRefTrk", &vRefTrk);

  tree_->Branch("nVtx",       &nVtx);
  tree_->Branch("run",        &run,       "run/I");
  tree_->Branch("event",        &event,     "event/I");
  tree_->Branch("lumiblock",&lumiblock,"lumiblock/I");

  tree_->Branch("jFlightLen", &jFlightLen);
  tree_->Branch("jFlightLenErr", &jFlightLenErr);
  tree_->Branch("jFlightLenSig", &jFlightLenSig);
  tree_->Branch("bFlightLen", &bFlightLen);
  tree_->Branch("bFlightLenErr", &bFlightLenErr);
  tree_->Branch("bFlightLenSig", &bFlightLenSig);

  tree_->Branch("jFlightLen2D", &jFlightLen2D);
  tree_->Branch("jFlightLenErr2D", &jFlightLenErr2D);
  tree_->Branch("jFlightLenSig2D", &jFlightLenSig2D);
  tree_->Branch("bFlightLen2D", &bFlightLen2D);
  tree_->Branch("bFlightLenErr2D", &bFlightLenErr2D);
  tree_->Branch("bFlightLenSig2D", &bFlightLenSig2D);
  tree_->Branch("bFlightTime", &bFlightTime);
  tree_->Branch("bFlightTimeErr", &bFlightTimeErr);

  tree_->Branch("JDecayVtxX",&JDecayVtxX);
  tree_->Branch("JDecayVtxY",&JDecayVtxY);
  tree_->Branch("JDecayVtxZ",&JDecayVtxZ);
  tree_->Branch("JDecayVtxXE",&JDecayVtxXE);
  tree_->Branch("JDecayVtxYE",&JDecayVtxYE);
  tree_->Branch("JDecayVtxZE",&JDecayVtxZE);
  tree_->Branch("JDecayVtxXYE",&JDecayVtxXYE);
  tree_->Branch("JDecayVtxXZE",&JDecayVtxXZE);
  tree_->Branch("JDecayVtxYZE",&JDecayVtxYZE);

  tree_->Branch("BDecayVtxX",&BDecayVtxX);
  tree_->Branch("BDecayVtxY",&BDecayVtxY);
  tree_->Branch("BDecayVtxZ",&BDecayVtxZ);
  tree_->Branch("BDecayVtxXE",&BDecayVtxXE);
  tree_->Branch("BDecayVtxYE",&BDecayVtxYE);
  tree_->Branch("BDecayVtxZE",&BDecayVtxZE);
  tree_->Branch("BDecayVtxXYE",&BDecayVtxXYE);
  tree_->Branch("BDecayVtxXZE",&BDecayVtxXZE);
  tree_->Branch("BDecayVtxYZE",&BDecayVtxYZE);

  tree_->Branch("mumC2",&mumC2);  
  tree_->Branch("mumNHits",&mumNHits);
  tree_->Branch("mumNPHits",&mumNPHits);
  tree_->Branch("mupC2",&mupC2);  
  tree_->Branch("mupNHits",&mupNHits);
  tree_->Branch("mupNPHits",&mupNPHits);
  tree_->Branch("mumdxy",&mumdxy);
  tree_->Branch("mupdxy",&mupdxy);
  tree_->Branch("muon_dca",&muon_dca);

  tree_->Branch("mumdz",&mumdz);
  tree_->Branch("mupdz",&mupdz);
  tree_->Branch("mu1soft",&mu1soft);
  tree_->Branch("mu2soft",&mu2soft);
  tree_->Branch("mu1tight",&mu1tight);
  tree_->Branch("mu2tight",&mu2tight);
  tree_->Branch("mu1PF",&mu1PF);
  tree_->Branch("mu2PF",&mu2PF);
  tree_->Branch("mu1loose",&mu1loose);
  tree_->Branch("mu2loose",&mu2loose);

}

// ------------ method called once each job just after ending the event loop  ------------
void JPsiPi::endJob() {
  tree_->GetDirectory()->cd();
  tree_->Write("", TObject::kOverwrite);
}

//define this as a plug-in
DEFINE_FWK_MODULE(JPsiPi);

