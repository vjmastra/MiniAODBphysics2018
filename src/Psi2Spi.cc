// -*- C++ -*-
//
// Package:    Psi2Spi
// Class:      Psi2Spi
// 
//=================================================
// original author:  Jhovanny Andres Mejia        |
//         created:  Fryday Sep 23                |
//         <jhovanny.andres.mejia.guisao@cern.ch> | 
//=================================================


// system include files
#include <memory>

// user include files
#include "myAnalyzers/JPsiKsPAT/src/Psi2Spi.h"

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
Psi2Spi::Psi2Spi(const edm::ParameterSet& iConfig)
  :
  dimuon_Label(consumes<edm::View<pat::Muon>>(iConfig.getParameter<edm::InputTag>("dimuons"))),
  trakCollection_label(consumes<edm::View<pat::PackedCandidate>>(iConfig.getParameter<edm::InputTag>("Trak"))),
  primaryVertices_Label(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("primaryVertices"))),
  BSLabel_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("bslabel"))),
  triggerResults_Label(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("TriggerResults"))),
//triggerObjects_(consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("objects"))),   

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

  tri_Dim25(0), tri_JpsiTk(0), tri_JpsiTkTk(0),

  mu1soft(0), mu2soft(0), mu1tight(0), mu2tight(0), 
  mu1PF(0), mu2PF(0), mu1loose(0), mu2loose(0),

  nVtx(0),
  priVtxX(0), priVtxY(0), priVtxZ(0), priVtxXE(0), priVtxYE(0), priVtxZE(0), priVtxCL(0),
  priVtxXYE(0), priVtxXZE(0), priVtxYZE(0),

  indexVtx(0), nTracksFromPV(0),
  
  vRefMuP(0), vRefMuM(0), vRefPi1(0), vRefPi2(0), vRefPi3(0),
 
  // ************************ ****************************************************

  BDecayVtxX(0), BDecayVtxY(0), BDecayVtxZ(0), BDecayVtxXE(0), BDecayVtxYE(0), BDecayVtxZE(0),
  BDecayVtxXYE(0), BDecayVtxXZE(0), BDecayVtxYZE(0),

  // *******************************************************
  nB(0), nMu(0), nJpsi(0), nPsi2S(0), nJpsi_test(0),
  B_mass(0), B_px(0), B_py(0), B_pz(0),

  piPi_mass(0), psiPiPi_mass(0),
  deltaR1(0), deltaR2(0), pointingAngle(0),

  pi1_px(0), pi1_py(0), pi1_pz(0), pi1_charge(0),
  pi1_px_track(0), pi1_py_track(0), pi1_pz_track(0),
  pi2_px(0), pi2_py(0), pi2_pz(0), pi2_charge(0),
  pi2_px_track(0), pi2_py_track(0), pi2_pz_track(0),
  pi3_px(0), pi3_py(0), pi3_pz(0), pi3_charge(0),
  pi3_px_track(0), pi3_py_track(0), pi3_pz_track(0),

  J_mass(0), J_px(0), J_py(0), J_pz(0),
  J_pt1(0), J_px1(0), J_py1(0), J_pz1(0), 
  J_pt2(0), J_px2(0), J_py2(0), J_pz2(0), 
  J_charge1(0), J_charge2(0),

  flightLen(0), flightLenErr(0), flightLenSig(0),
/*  
  mu1_px_test(0), mu1_py_test(0), mu1_pz_test(0), mu1_charge_test(0),
  mu2_px_test(0), mu2_py_test(0), mu2_pz_test(0), mu2_charge_test(0),
  Jpsi_dca_test(0),
  Jpsi_vx_test(0), Jpsi_vy_test(0), Jpsi_vz_test(0),

  Jpsi_mass_test(0), Jpsi_prob_test(0), Jpsi_chi2_test(0),
*/
  J_chi2(0), psi2S_chi2(0), B_chi2(0),
  B_Prob(0), J_Prob(0), psi2S_Prob(0),

  run(0), event(0),
  lumiblock(0)

{
   //now do what ever initialization is needed
}

Psi2Spi::~Psi2Spi()
{

}

// ------------ method called to for each event  ------------
void Psi2Spi::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using std::vector;
  using namespace edm;
  using namespace reco;
  using namespace std;
  
  // *********************************
  // Get event content information
  // *********************************  

  edm::ESHandle<TransientTrackBuilder> theB; 
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB); 

  edm::Handle< View<pat::PackedCandidate> > thePATTrackHandle;
  iEvent.getByToken(trakCollection_label,thePATTrackHandle);

  edm::Handle< View<pat::Muon> > thePATMuonHandle;
  iEvent.getByToken(dimuon_Label,thePATMuonHandle);

  edm::Handle<edm::TriggerResults> triggerBits;
//  edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;

  iEvent.getByToken(triggerResults_Label, triggerBits);
//  iEvent.getByToken(triggerObjects_, triggerObjects); 
//  const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);

//  for ( size_t iTrigObj = 0; iTrigObj < triggerObjects->size(); ++iTrigObj ) {
//    pat::TriggerObjectStandAlone obj( triggerObjects->at( iTrigObj ) );
//    obj.unpackPathNames(names);
//    obj.unpackFilterLabels(iEvent, *triggerBits);
//  }

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
 
  //Control loop on J/psi
/*
  for(View<pat::Muon>::const_iterator iMuonA = thePATMuonHandle->begin(); iMuonA != thePATMuonHandle->end(); ++iMuonA) {
    for(View<pat::Muon>::const_iterator iMuonB = iMuonA+1; iMuonB != thePATMuonHandle->end(); ++iMuonB) {
   
      if (iMuonA == iMuonB) continue;
      if ( (iMuonA->charge())*(iMuonB->charge()) == 1 ) continue;
     
      TrackRef muonTrackA = iMuonA->track();
      TrackRef muonTrackB = iMuonB->track();

      if (muonTrackA.isNull() || muonTrackB.isNull()) continue;
 
      if (muonTrackA->pt()<4.0) continue;
      if (muonTrackB->pt()<4.0) continue;

      if ( !(muonTrackA->quality(reco::TrackBase::highPurity)) ) continue;
      if ( !(muonTrackB->quality(reco::TrackBase::highPurity)) ) continue;     

      reco::TransientTrack muonATT((*theB).build(muonTrackA));
      reco::TransientTrack muonBTT((*theB).build(muonTrackB));

      FreeTrajectoryState muAState = muonATT.impactPointTSCP().theState();
      FreeTrajectoryState muBState = muonBTT.impactPointTSCP().theState();

      if( !muonATT.impactPointTSCP().isValid() || !muonBTT.impactPointTSCP().isValid() ) continue;

      ClosestApproachInRPhi cApp;
      cApp.calculate(muAState, muBState);
      if( !cApp.status() ) continue;
      float dca = fabs( cApp.distance() );
      if (dca < 0. || dca > 0.5) continue;

      ParticleMass muon_mass = 0.10565837;
      float muon_sigma = muon_mass*1.e-6;

      KinematicParticleFactoryFromTransientTrack pFactory;
 
      float chi = 0.;
      float ndf = 0.;
      vector<RefCountedKinematicParticle> muonParticles_test;
      try {
        muonParticles_test.push_back(pFactory.particle(muonATT,muon_mass,chi,ndf,muon_sigma));
        muonParticles_test.push_back(pFactory.particle(muonBTT,muon_mass,chi,ndf,muon_sigma));
      }
      catch(...) {
        std::cout<<" Exception caught ... continuing 1 "<<std::endl;
        continue;
      }

      KinematicParticleVertexFitter fitter;

      RefCountedKinematicTree Jpsi_VertexFitTree_test;
      try {
        Jpsi_VertexFitTree_test = fitter.fit(muonParticles_test);
      }
      catch (...) {
        std::cout<<" Exception caught ... continuing 2 "<<std::endl;
        continue;
      }

      if (!Jpsi_VertexFitTree_test->isValid()) {
        //std::cout << "caught an exception in the psi vertex fit" << std::endl;
        continue;
      }

      Jpsi_VertexFitTree_test->movePointerToTheTop();

      RefCountedKinematicParticle Jpsi_vFit_test = Jpsi_VertexFitTree_test->currentParticle();
      RefCountedKinematicVertex Jpsi_vFit_vertex_test = Jpsi_VertexFitTree_test->currentDecayVertex();

      if(Jpsi_vFit_vertex_test->chiSquared() < 0 ) {
        std::cout << "negative chisq from psi fit" << endl;
        continue;
      }

      double Jpsi_prob_tmp_test   = TMath::Prob(Jpsi_vFit_vertex_test->chiSquared(),(int)Jpsi_vFit_vertex_test->degreesOfFreedom());
      if(Jpsi_prob_tmp_test<0.01) {
        continue;
      }

      //some loose cuts go here
      double Jpsi_mass_tmp_test = Jpsi_vFit_test->currentState().mass();
      if( Jpsi_mass_tmp_test < 3.0 || Jpsi_mass_tmp_test > 3.2) continue;
   
      //Write

      Jpsi_VertexFitTree_test->movePointerToTheFirstChild();
      RefCountedKinematicParticle muACand_test = Jpsi_VertexFitTree_test->currentParticle();
      Jpsi_VertexFitTree_test->movePointerToTheNextChild();
      RefCountedKinematicParticle muBCand_test = Jpsi_VertexFitTree_test->currentParticle();

      KinematicParameters muAKP_test = muACand_test->currentState().kinematicParameters();
      KinematicParameters muBKP_test = muBCand_test->currentState().kinematicParameters();

      nJpsi_test++;     
      mu1_px_test->push_back(muAKP_test.momentum().x());
      mu1_py_test->push_back(muAKP_test.momentum().y());
      mu1_pz_test->push_back(muAKP_test.momentum().z());
      mu1_charge_test->push_back(muACand_test->currentState().particleCharge());
      mu2_px_test->push_back(muBKP_test.momentum().x());
      mu2_py_test->push_back(muBKP_test.momentum().x());
      mu2_pz_test->push_back(muBKP_test.momentum().x());
      mu2_charge_test->push_back(muBCand_test->currentState().particleCharge());
      Jpsi_vx_test->push_back((*Jpsi_vFit_vertex_test).position().x());
      Jpsi_vy_test->push_back((*Jpsi_vFit_vertex_test).position().y());
      Jpsi_vz_test->push_back((*Jpsi_vFit_vertex_test).position().z());
      Jpsi_dca_test->push_back(dca);
      Jpsi_mass_test->push_back(Jpsi_mass_tmp_test);
      Jpsi_prob_test->push_back(Jpsi_prob_tmp_test);   
      Jpsi_chi2_test->push_back(Jpsi_vFit_vertex_test->chiSquared());

      //Clean

      muonParticles_test.clear();
   
    }
  }
*/
  //*****************************************
  //Let's begin by looking for J/psi->mu+mu-

  unsigned int nMu_tmp = thePATMuonHandle->size();
 
  for(View<pat::Muon>::const_iterator iMuon1 = thePATMuonHandle->begin(); iMuon1 != thePATMuonHandle->end(); ++iMuon1) { 
    for(View<pat::Muon>::const_iterator iMuon2 = iMuon1+1; iMuon2 != thePATMuonHandle->end(); ++iMuon2) {  
      
      if(iMuon1==iMuon2) continue;
      if( (iMuon1->charge())*(iMuon2->charge()) == 1) continue;

      TrackRef glbTrackP;	  
      TrackRef glbTrackM;	  
  
      if(iMuon1->charge() == 1){glbTrackP = iMuon1->track();}
      if(iMuon1->charge() == -1){glbTrackM = iMuon1->track();}
  
      if(iMuon2->charge() == 1) {glbTrackP = iMuon2->track();}
      if(iMuon2->charge() == -1){glbTrackM = iMuon2->track();}
	  
      if( glbTrackP.isNull() || glbTrackM.isNull() ) {
        //std::cout << "continue due to no track ref" << endl;
        continue;
      }

      if(iMuon1->track()->pt()<4.0) continue;
      if(iMuon2->track()->pt()<4.0) continue;

      if(!(glbTrackM->quality(reco::TrackBase::highPurity))) continue;
      if(!(glbTrackP->quality(reco::TrackBase::highPurity))) continue;
	  
      //Let's check the vertex and mass
      reco::TransientTrack muon1TT((*theB).build(glbTrackP));
      reco::TransientTrack muon2TT((*theB).build(glbTrackM));

      // *****  Trajectory states to calculate DCA for the 2 muons *********************
      FreeTrajectoryState mu1State = muon1TT.impactPointTSCP().theState();
      FreeTrajectoryState mu2State = muon2TT.impactPointTSCP().theState();
 
      if( !muon1TT.impactPointTSCP().isValid() || !muon2TT.impactPointTSCP().isValid() ) continue;

      // Measure distance between tracks at their closest approach
      ClosestApproachInRPhi cApp;
      cApp.calculate(mu1State, mu2State);
      if( !cApp.status() ) continue;
      float dca = fabs( cApp.distance() );	  
      if (dca < 0. || dca > 0.5) continue;
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

      double J_Prob_tmp   = TMath::Prob(J_vFit_vertex_noMC->chiSquared(),(int)J_vFit_vertex_noMC->degreesOfFreedom());
      if(J_Prob_tmp<0.01) {
        continue;
      }
	  
      //some loose cuts go here

      if(J_vFit_noMC->currentState().mass()<3.0 || J_vFit_noMC->currentState().mass()>3.2) continue;
 
      nJpsi++;     
      //jpsipipi

      for (View<pat::PackedCandidate>::const_iterator iTrack1 = thePATTrackHandle->begin(); iTrack1 != thePATTrackHandle->end(); ++iTrack1 ) {
        for (View<pat::PackedCandidate>::const_iterator iTrack2 = iTrack1+1; iTrack2 != thePATTrackHandle->end(); ++iTrack2 ) {

          if (iTrack1==iTrack2) continue;
          if ((iTrack1->charge())*(iTrack2->charge()) != -1) continue;  //opposite charge
          if (iTrack1->pt()<0.95) continue;
          if (iTrack2->pt()<0.95) continue;
          if (!(iTrack1->trackHighPurity())) continue;
          if (!(iTrack2->trackHighPurity())) continue;
          if (IsTheSame(*iTrack1,*iMuon1) || IsTheSame(*iTrack1,*iMuon2)) continue;
          if (IsTheSame(*iTrack2,*iMuon1) || IsTheSame(*iTrack2,*iMuon2)) continue;

          reco::TransientTrack pion1TT((*theB).build(iTrack1->pseudoTrack()));
          reco::TransientTrack pion2TT((*theB).build(iTrack2->pseudoTrack())); 

          ParticleMass pion_mass = 0.13957018;
          float pion_sigma = pion_mass*1.e-6;
          float chi = 0.;
          float ndf = 0.;

          // Jpsipionpion invariant mass (before kinematic vertex fit)
          TLorentzVector pion14V, pion24V, Jpsi4V, psi2S4V;
          pion14V.SetXYZM(iTrack1->px(),iTrack1->py(),iTrack1->pz(),pion_mass);
          pion24V.SetXYZM(iTrack2->px(),iTrack2->py(),iTrack2->pz(),pion_mass);
          Jpsi4V.SetXYZM(J_vFit_noMC->currentState().globalMomentum().x(), J_vFit_noMC->currentState().globalMomentum().y(), J_vFit_noMC->currentState().globalMomentum().z(),J_vFit_noMC->currentState().mass());          
          psi2S4V = pion14V + pion24V + Jpsi4V;
          double piPiMass = (pion14V + pion24V).M();
          double psiPiPiMass = (psi2S4V).M();       

          if (piPiMass < 0.45) continue;
          if (psiPiPiMass < 3.4 || psiPiPiMass > 4.2) continue;

          float dR1_tmp;
          float dR2_tmp;

          dR1_tmp = Jpsi4V.DeltaR(pion14V);
          dR2_tmp = Jpsi4V.DeltaR(pion24V);

          if (dR1_tmp > 0.3) continue;
          if (dR2_tmp > 0.3) continue; 

          // JPsi mass constraint is applied

          vector<RefCountedKinematicParticle> vFitMCParticles;
          vFitMCParticles.push_back(pFactory.particle(muon1TT,muon_mass,chi,ndf,muon_sigma));
          vFitMCParticles.push_back(pFactory.particle(muon2TT,muon_mass,chi,ndf,muon_sigma));
          vFitMCParticles.push_back(pFactory.particle(pion1TT,pion_mass,chi,ndf,pion_sigma));
          vFitMCParticles.push_back(pFactory.particle(pion2TT,pion_mass,chi,ndf,pion_sigma));

          MultiTrackKinematicConstraint *  Jpsi_c = new  TwoTrackMassKinematicConstraint(Jpsi_mass);
          KinematicConstrainedVertexFitter kcvFitter;
          RefCountedKinematicTree psi2SVertexFitTree = kcvFitter.fit(vFitMCParticles, Jpsi_c);
          if (!psi2SVertexFitTree->isValid()) continue;

          psi2SVertexFitTree->movePointerToTheTop();
          RefCountedKinematicParticle psi2SCandMC = psi2SVertexFitTree->currentParticle();          
          RefCountedKinematicVertex psi2SDecayVertexMC = psi2SVertexFitTree->currentDecayVertex();
          if (!psi2SDecayVertexMC->vertexIsValid()) continue;
          if ((psi2SCandMC->currentState().mass() < 3.4) || (psi2SCandMC->currentState().mass() > 4.2) ) continue;

          double psi2S_Prob_tmp = TMath::Prob(psi2SDecayVertexMC->chiSquared(),(int)psi2SDecayVertexMC->degreesOfFreedom());
          if (psi2S_Prob_tmp<0.01) continue;

          nPsi2S++;

          for(View<pat::PackedCandidate>::const_iterator iTrack3 = thePATTrackHandle->begin(); iTrack3 != thePATTrackHandle->end(); ++iTrack3) {

            if(iTrack3->charge()==0) continue;
            if(iTrack3->pt()<0.95) continue;
            if(!(iTrack3->trackHighPurity())) continue;

            if ( IsTheSame(*iTrack3,*iMuon1) || IsTheSame(*iTrack3,*iMuon2) ) continue;
            if ( IsTheSame2(*iTrack3,*iTrack1) || IsTheSame2(*iTrack3,*iTrack2) ) continue;

            reco::TransientTrack pion3TT((*theB).build(iTrack3->pseudoTrack()));

            float chi = 0.;
            float ndf = 0.;

            TLorentzVector pion34V, Bcand4V;
            pion34V.SetXYZM(iTrack3->px(), iTrack3->py(), iTrack3->pz(), pion_mass);
            Bcand4V = psi2S4V + pion34V;
            if ( Bcand4V.M()<5.0 || Bcand4V.M()>7.0 ) continue;

	    //Now we are ready to combine!
	    // JPsi mass constraint is applied in the final Bd fit,
		     
  	    vector<RefCountedKinematicParticle> vFitMCParticles2;

            try {
	      vFitMCParticles2.push_back(pFactory.particle(muon1TT,muon_mass,chi,ndf,muon_sigma));
	      vFitMCParticles2.push_back(pFactory.particle(muon2TT,muon_mass,chi,ndf,muon_sigma));
              vFitMCParticles2.push_back(pFactory.particle(pion1TT,pion_mass,chi,ndf,pion_sigma));
              vFitMCParticles2.push_back(pFactory.particle(pion2TT,pion_mass,chi,ndf,pion_sigma));
	      vFitMCParticles2.push_back(pFactory.particle(pion3TT,pion_mass,chi,ndf,pion_sigma));
              }
            catch (...) {
              std::cout << "Exception caught ... continuing 5" << std::endl;
              continue;
            }

            MultiTrackKinematicConstraint *  Jpsi_c2 = new TwoTrackMassKinematicConstraint(Jpsi_mass);
            KinematicConstrainedVertexFitter kcvFitter2;

            RefCountedKinematicTree BVertexFitTree;
            try {
              BVertexFitTree = kcvFitter2.fit(vFitMCParticles2, Jpsi_c2);
            }
            catch (...) {
              std::cout<<" Exception caught ... continuing 6 "<<std::endl;
              continue;
            }

            if (!BVertexFitTree->isValid()) {
              //std::cout << "caught an exception in the lambdaB vertex fit" << std::endl;
              continue;
            }
	     
	    BVertexFitTree->movePointerToTheTop();		     
	     
	    RefCountedKinematicParticle BCandMC = BVertexFitTree->currentParticle();
	    RefCountedKinematicVertex BDecayVertexMC = BVertexFitTree->currentDecayVertex();
	    if (!BDecayVertexMC->vertexIsValid()){
	      //std::cout << "B MC fit vertex is not valid" << endl;
	      continue;
	    }
	     
	    if(BCandMC->currentState().mass()< 5.95 || BCandMC->currentState().mass()> 6.45) continue;
	     
	    if(BDecayVertexMC->chiSquared()<0) {
              //std::cout << " continue from negative chi2 = " << bDecayVertexMC->chiSquared() << endl;
	      continue;
	    }
		     
	    double B_Prob_tmp = TMath::Prob(BDecayVertexMC->chiSquared(),(int)BDecayVertexMC->degreesOfFreedom());
	    if(B_Prob_tmp<0.01) {
              continue;
	    }		     

            reco::Vertex bestVtx = highestptVtx;

            float primaryVertex[3] = {(float)bestVtx.x(), (float)bestVtx.y(), (float)bestVtx.z()};
            float secundaryVertex[3] = {(*BDecayVertexMC).position().x(), (*BDecayVertexMC).position().y(), (*BDecayVertexMC).position().z()};
            float flightVec[3];
            for (int i = 0; i < 3; i++) flightVec[i] = secundaryVertex[i] - primaryVertex[i];
            TVector3 flightDir(flightVec[0], flightVec[1], flightVec[2]);
            TVector3 Bmomentum(BCandMC->currentState().globalMomentum().x(), BCandMC->currentState().globalMomentum().y(), BCandMC->currentState().globalMomentum().z());

            double cos_alpha = TMath::Cos(flightDir.Angle(Bmomentum));
            if (cos_alpha < 0.95) continue;

            // get children from final B fit
	    BVertexFitTree->movePointerToTheFirstChild();
	    RefCountedKinematicParticle mu1CandMC = BVertexFitTree->currentParticle();
	    BVertexFitTree->movePointerToTheNextChild();
	    RefCountedKinematicParticle mu2CandMC = BVertexFitTree->currentParticle();
            BVertexFitTree->movePointerToTheNextChild();
            RefCountedKinematicParticle pi1CandMC = BVertexFitTree->currentParticle();
            BVertexFitTree->movePointerToTheNextChild();
            RefCountedKinematicParticle pi2CandMC = BVertexFitTree->currentParticle();
            BVertexFitTree->movePointerToTheNextChild();
            RefCountedKinematicParticle pi3CandMC = BVertexFitTree->currentParticle();
		   
	    KinematicParameters psiMu1KP = mu1CandMC->currentState().kinematicParameters();
	    KinematicParameters psiMu2KP = mu2CandMC->currentState().kinematicParameters();
	    KinematicParameters psiMupKP;
	    KinematicParameters psiMumKP;
	       
	    if ( mu1CandMC->currentState().particleCharge() > 0 ) psiMupKP = psiMu1KP;
	    if ( mu1CandMC->currentState().particleCharge() < 0 ) psiMumKP = psiMu1KP;
	    if ( mu2CandMC->currentState().particleCharge() > 0 ) psiMupKP = psiMu2KP;
            if ( mu2CandMC->currentState().particleCharge() < 0 ) psiMumKP = psiMu2KP;	 

            KinematicParameters psiPi1KP = pi1CandMC->currentState().kinematicParameters();
            KinematicParameters psiPi2KP = pi2CandMC->currentState().kinematicParameters();
            KinematicParameters psiPipKP;
            KinematicParameters psiPimKP;

            if ( pi1CandMC->currentState().particleCharge() > 0 ) psiPipKP = psiPi1KP;
            if ( pi1CandMC->currentState().particleCharge() < 0 ) psiPimKP = psiPi1KP;
            if ( pi2CandMC->currentState().particleCharge() > 0 ) psiPipKP = psiPi2KP;
            if ( pi2CandMC->currentState().particleCharge() < 0 ) psiPimKP = psiPi2KP;

            KinematicParameters pi3KP = pi3CandMC->currentState().kinematicParameters();

 	    GlobalVector Jp1vec(mu1CandMC->currentState().globalMomentum().x(),
	    			mu1CandMC->currentState().globalMomentum().y(),
 	      			mu1CandMC->currentState().globalMomentum().z());

            GlobalVector Jp2vec(mu2CandMC->currentState().globalMomentum().x(),
	      			mu2CandMC->currentState().globalMomentum().y(),
 	      			mu2CandMC->currentState().globalMomentum().z());

 	    GlobalVector p1vec( pi1CandMC->currentState().globalMomentum().x(),
	      	  		pi1CandMC->currentState().globalMomentum().y(),
 	  			pi1CandMC->currentState().globalMomentum().z());

            GlobalVector p2vec( pi2CandMC->currentState().globalMomentum().x(),
                                pi2CandMC->currentState().globalMomentum().y(),
                                pi2CandMC->currentState().globalMomentum().z());

            GlobalVector p3vec( pi3CandMC->currentState().globalMomentum().x(),
                                pi3CandMC->currentState().globalMomentum().y(),
                                pi3CandMC->currentState().globalMomentum().z());

            // fill candidate variables now
		   
            if(nB==0){		    
	      nMu  = nMu_tmp;
	      // cout<< "*Number of Muons : " << nMu_tmp << endl;
	    }		     

	    B_mass->push_back(BCandMC->currentState().mass());
	    B_px->push_back(BCandMC->currentState().globalMomentum().x());
	    B_py->push_back(BCandMC->currentState().globalMomentum().y());
	    B_pz->push_back(BCandMC->currentState().globalMomentum().z());

            piPi_mass->push_back(piPiMass);
            psiPiPi_mass->push_back( psi2SCandMC->currentState().mass());
            deltaR1->push_back(dR1_tmp);
            deltaR2->push_back(dR2_tmp);
            pointingAngle->push_back(cos_alpha);

            J_mass->push_back( J_vFit_noMC->currentState().mass() );
            J_px->push_back( J_vFit_noMC->currentState().globalMomentum().x() );
	    J_py->push_back( J_vFit_noMC->currentState().globalMomentum().y() );
	    J_pz->push_back( J_vFit_noMC->currentState().globalMomentum().z() );

            J_pt1->push_back(Jp1vec.perp());
	    J_px1->push_back(psiMu1KP.momentum().x());
	    J_py1->push_back(psiMu1KP.momentum().y());
	    J_pz1->push_back(psiMu1KP.momentum().z());
	    J_charge1->push_back(mu1CandMC->currentState().particleCharge());
 
            J_pt2->push_back(Jp2vec.perp());
	    J_px2->push_back(psiMu2KP.momentum().x());
	    J_py2->push_back(psiMu2KP.momentum().y());
	    J_pz2->push_back(psiMu2KP.momentum().z());
	    J_charge2->push_back(mu2CandMC->currentState().particleCharge());

            pi1_px->push_back(psiPi1KP.momentum().x());
            pi1_py->push_back(psiPi1KP.momentum().y());
            pi1_pz->push_back(psiPi1KP.momentum().z());
            pi1_charge->push_back(pi1CandMC->currentState().particleCharge());
            pi1_px_track->push_back(iTrack1->px());
            pi1_py_track->push_back(iTrack1->py());
            pi1_pz_track->push_back(iTrack1->pz());

            pi2_px->push_back(psiPi2KP.momentum().x());
            pi2_py->push_back(psiPi2KP.momentum().y());
            pi2_pz->push_back(psiPi2KP.momentum().z());
            pi2_charge->push_back(pi2CandMC->currentState().particleCharge());
            pi2_px_track->push_back(iTrack2->px());
            pi2_py_track->push_back(iTrack2->py());
            pi2_pz_track->push_back(iTrack2->pz());

            pi3_px->push_back(pi3KP.momentum().x());
            pi3_py->push_back(pi3KP.momentum().y());
            pi3_pz->push_back(pi3KP.momentum().z());
            pi3_charge->push_back(pi3CandMC->currentState().particleCharge());
            pi3_px_track->push_back(iTrack3->px());
            pi3_py_track->push_back(iTrack3->py());
            pi3_pz_track->push_back(iTrack3->pz());

	    J_chi2->push_back(J_vFit_vertex_noMC->chiSquared());
            psi2S_chi2->push_back( psi2SDecayVertexMC->chiSquared());
	    B_chi2->push_back(BDecayVertexMC->chiSquared());
 
	    B_Prob->push_back(B_Prob_tmp); 
            J_Prob->push_back(J_Prob_tmp);
            psi2S_Prob->push_back( psi2S_Prob_tmp);
 
	    BDecayVtxX->push_back((*BDecayVertexMC).position().x());
	    BDecayVtxY->push_back((*BDecayVertexMC).position().y());
	    BDecayVtxZ->push_back((*BDecayVertexMC).position().z());
	    BDecayVtxXE->push_back(BDecayVertexMC->error().cxx());
	    BDecayVtxYE->push_back(BDecayVertexMC->error().cyy());
	    BDecayVtxZE->push_back(BDecayVertexMC->error().czz());
	    BDecayVtxXYE->push_back(BDecayVertexMC->error().cyx());
	    BDecayVtxXZE->push_back(BDecayVertexMC->error().czx());
	    BDecayVtxYZE->push_back(BDecayVertexMC->error().czy());

 // ********************* muon-trigger-machint**************** 
/*		   
	      const pat::TriggerObjectStandAloneCollection muHLTMatches1_t1 = iMuon1->triggerObjectMatchesByFilter("hltDisplacedmumuFilterDimuon25Jpsis");
	      const pat::TriggerObjectStandAloneCollection muHLTMatches2_t1 = iMuon2->triggerObjectMatchesByFilter("hltDisplacedmumuFilterDimuon25Jpsis");
 		   
	      const pat::TriggerObjectStandAloneCollection muHLTMatches1_t2 = iMuon1->triggerObjectMatchesByFilter("hltJpsiTkVertexFilter");
	      const pat::TriggerObjectStandAloneCollection muHLTMatches2_t2 = iMuon2->triggerObjectMatchesByFilter("hltJpsiTkVertexFilter");
		   
	      const pat::TriggerObjectStandAloneCollection muHLTMatches1_t4 = iMuon1->triggerObjectMatchesByFilter("hltJpsiTkTkVertexFilterPhiKstar");
	      const pat::TriggerObjectStandAloneCollection muHLTMatches2_t4 = iMuon2->triggerObjectMatchesByFilter("hltJpsiTkTkVertexFilterPhiKstar");

	      int tri_Dim25_tmp = 0, tri_JpsiTk_tmp = 0,  tri_JpsiTkTk_tmp = 0;
		   
	      if (muHLTMatches1_t1.size() > 0 && muHLTMatches2_t1.size() > 0) tri_Dim25_tmp = 1;
	      if (muHLTMatches1_t2.size() > 0 && muHLTMatches2_t2.size() > 0) tri_JpsiTk_tmp = 1;
	      if (muHLTMatches1_t4.size() > 0 && muHLTMatches2_t4.size() > 0) tri_JpsiTkTk_tmp = 1;
		   
	      tri_Dim25->push_back( tri_Dim25_tmp );	       
	      tri_JpsiTk->push_back( tri_JpsiTk_tmp );
              tri_JpsiTkTk->push_back( tri_JpsiTkTk_tmp );
*/
    	      // ************
		  
	    mu1soft->push_back(iMuon1->isSoftMuon(bestVtx) );
	    mu2soft->push_back(iMuon2->isSoftMuon(bestVtx) );
	    mu1tight->push_back(iMuon1->isTightMuon(bestVtx) );
	    mu2tight->push_back(iMuon2->isTightMuon(bestVtx) );
	    mu1PF->push_back(iMuon1->isPFMuon());
	    mu2PF->push_back(iMuon2->isPFMuon());
	    mu1loose->push_back(muon::isLooseMuon(*iMuon1));
	    mu2loose->push_back(muon::isLooseMuon(*iMuon2));

            mumC2->push_back( glbTrackM->normalizedChi2() );
	    mumNHits->push_back( glbTrackM->numberOfValidHits() );
	    mumNPHits->push_back( glbTrackM->hitPattern().numberOfValidPixelHits() );	       
	    mupC2->push_back( glbTrackP->normalizedChi2() );
	    mupNHits->push_back( glbTrackP->numberOfValidHits() );
	    mupNPHits->push_back( glbTrackP->hitPattern().numberOfValidPixelHits() );
            mumdxy->push_back(glbTrackM->dxy(bestVtx.position()) );
	    mupdxy->push_back(glbTrackP->dxy(bestVtx.position()) );
	    mumdz->push_back(glbTrackM->dz(bestVtx.position()) );
	    mupdz->push_back(glbTrackP->dz(bestVtx.position()) );
	    muon_dca->push_back(dca);
 
	    // try refitting the primary without the tracks in the B reco candidate		   
		  
	    nB++;	       
		   
            ////////////////////////////////////
            muonParticles.clear();
	    vFitMCParticles.clear();
            vFitMCParticles2.clear();
	   
           }//pi3
         }//pi2
       }//pi1
     }//muon2
   }//muon1
 
   
   //fill the tree and clear the vectors
/*
   if (nJpsi_test > 0) {
     treeTest_->Fill();
   }
*/
   if (nB > 0 ) {
       //std::cout << "filling tree" << endl;
     tree_->Fill();
   }
   // *********

   nB = 0; nMu = 0; nJpsi = 0; nPsi2S = 0; nJpsi_test = 0;

   B_mass->clear();    B_px->clear();    B_py->clear();    B_pz->clear();

   piPi_mass->clear(); psiPiPi_mass->clear();
   deltaR1->clear(); deltaR2->clear(); pointingAngle->clear();

   J_mass->clear();  J_px->clear();  J_py->clear();  J_pz->clear();

   J_pt1->clear();  J_px1->clear();  J_py1->clear();  J_pz1->clear(), J_charge1->clear();
   J_pt2->clear();  J_px2->clear();  J_py2->clear();  J_pz2->clear(), J_charge2->clear();

   pi1_px->clear(); pi1_py->clear(); pi1_pz->clear(); pi1_charge->clear();
   pi1_px_track->clear(); pi1_py_track->clear(); pi1_pz_track->clear();
   pi2_px->clear(); pi2_py->clear(); pi2_pz->clear(); pi2_charge->clear();
   pi2_px_track->clear(); pi2_py_track->clear(); pi2_pz_track->clear();
   pi3_px->clear(); pi3_py->clear(); pi3_pz->clear(); pi3_charge->clear();
   pi3_px_track->clear(); pi3_py_track->clear(); pi3_pz_track->clear();

   J_chi2->clear(); psi2S_chi2->clear(); B_chi2->clear();
   B_Prob->clear(); J_Prob->clear(); psi2S_Prob->clear();
/*
   mu1_px_test->clear(); mu1_py_test->clear(); mu1_pz_test->clear(); mu1_charge_test->clear();
   mu2_px_test->clear(); mu2_py_test->clear(); mu2_pz_test->clear(); mu2_charge_test->clear();
   Jpsi_dca_test->clear();
   Jpsi_vx_test->clear(); Jpsi_vy_test->clear(); Jpsi_vz_test->clear();

   Jpsi_mass_test->clear(); Jpsi_prob_test->clear(); Jpsi_chi2_test->clear();
*/
   // *********

   nVtx = 0;
   priVtxX->clear(); priVtxY->clear(); priVtxZ->clear();
   priVtxXE->clear(); priVtxYE->clear(); priVtxZE->clear();
   priVtxXYE->clear(); priVtxXZE->clear(); priVtxYZE->clear();
   priVtxCL->clear();

   indexVtx->clear(); nTracksFromPV->clear();

   vRefMuP->clear(); vRefMuM->clear();
   vRefPi1->clear(); vRefPi2->clear(); vRefPi3->clear();

   flightLen->clear(); flightLenErr->clear(); flightLenSig->clear();

   BDecayVtxX->clear(); BDecayVtxY->clear(); BDecayVtxZ->clear(); 
   BDecayVtxXE->clear(); BDecayVtxYE->clear(); BDecayVtxZE->clear(); 
   BDecayVtxXYE->clear(); BDecayVtxXZE->clear(); BDecayVtxYZE->clear();  

   mumC2->clear();
   mumNHits->clear(); mumNPHits->clear();
   mupC2->clear();
   mupNHits->clear(); mupNPHits->clear();
   mumdxy->clear(); mupdxy->clear(); mumdz->clear(); mupdz->clear(); muon_dca->clear();

   tri_Dim25->clear(); tri_JpsiTk->clear(); tri_JpsiTkTk->clear();
      
   mu1soft->clear(); mu2soft->clear(); mu1tight->clear(); mu2tight->clear();
   mu1PF->clear(); mu2PF->clear(); mu1loose->clear(); mu2loose->clear(); 

}

bool Psi2Spi::IsTheSame(const pat::PackedCandidate& tk, const pat::Muon& mu){
  double DeltaEta = fabs(mu.eta()-tk.eta());
  double DeltaP   = fabs(mu.p()-tk.p());
  if (DeltaEta < 0.02 && DeltaP < 0.02) return true;
  return false;
}

bool Psi2Spi::IsTheSame2(const pat::PackedCandidate& tk1, const pat::PackedCandidate& tk2){
  double DeltaEta = fabs(tk1.eta()-tk2.eta());
  double DeltaP   = fabs(tk1.p()-tk2.p());
  if (DeltaEta < 0.02 && DeltaP < 0.02) return true;
  return false;
}


Double_t Psi2Spi::Distance(const Double_t p1[], const Double_t p2[]){
  Double_t diff[3], diff2[3];
  Double_t dist = 0;
  for (int i = 0; i < 3; i++) {
    diff[i] = p2[i]-p1[i];
    diff2[i] = diff[i]*diff[i];
    dist += diff2[i];
  }
  return TMath::Sqrt(dist);
}

Double_t Psi2Spi::DistanceError(const Double_t p1[], const Double_t err1[], const Double_t p2[], const Double_t err2[]){
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

// ------------ method called once each job just before starting event loop  ------------

void Psi2Spi::beginJob()
{

  std::cout << "Beginning analyzer job with value of doMC = " << doMC_ << std::endl;

  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>("ntuple","Bc->Psi2Spi ntuple");
/*
  treeTest_ = fs->make<TTree>("test","Jpsi->mumu");
*/
  tree_->Branch("nB",&nB,"nB/i");
  tree_->Branch("nMu",&nMu,"nMu/i");
  tree_->Branch("nJpsi",&nJpsi,"nJpsi/i");
  tree_->Branch("nPsi2S",&nPsi2S,"nPsi2S/i");

  tree_->Branch("B_mass", &B_mass);
  tree_->Branch("B_px", &B_px);
  tree_->Branch("B_py", &B_py);
  tree_->Branch("B_pz", &B_pz);

  tree_->Branch("PiPi_mass", &piPi_mass);
  tree_->Branch("PsiPiPi_mass", &psiPiPi_mass);  
  tree_->Branch("deltaR_J_p1", &deltaR1);
  tree_->Branch("deltaR_J_p2", &deltaR2);
  tree_->Branch("cosAlpha", &pointingAngle);

  tree_->Branch("J_mass", &J_mass);
  tree_->Branch("J_px", &J_px);
  tree_->Branch("J_py", &J_py);
  tree_->Branch("J_pz", &J_pz);

  tree_->Branch("J_pt1", &J_pt1);
  tree_->Branch("J_px1", &J_px1);
  tree_->Branch("J_py1", &J_py1);
  tree_->Branch("J_pz1", &J_pz1);
  tree_->Branch("J_charge1", &J_charge1);

  tree_->Branch("J_pt2", &J_pt2);
  tree_->Branch("J_px2", &J_px2);
  tree_->Branch("J_py2", &J_py2);
  tree_->Branch("J_pz2", &J_pz2);
  tree_->Branch("J_charge2", &J_charge2);

  tree_->Branch("pi1_px", &pi1_px);
  tree_->Branch("pi1_py", &pi1_py);
  tree_->Branch("pi1_pz", &pi1_pz);
  tree_->Branch("pi1_charge", &pi1_charge);
  tree_->Branch("pi1_px_track", &pi1_px_track);
  tree_->Branch("pi1_py_track", &pi1_py_track);
  tree_->Branch("pi1_pz_track", &pi1_pz_track);
  tree_->Branch("pi2_px", &pi2_px);
  tree_->Branch("pi2_py", &pi2_py);
  tree_->Branch("pi2_pz", &pi2_pz);
  tree_->Branch("pi2_charge", &pi2_charge);
  tree_->Branch("pi2_px_track", &pi2_px_track);
  tree_->Branch("pi2_py_track", &pi2_py_track);
  tree_->Branch("pi2_pz_track", &pi2_pz_track);
  tree_->Branch("pi3_px", &pi3_px);
  tree_->Branch("pi3_py", &pi3_py);
  tree_->Branch("pi3_pz", &pi3_pz);
  tree_->Branch("pi3_charge", &pi3_charge);
  tree_->Branch("pi3_px_track", &pi3_px_track);
  tree_->Branch("pi3_py_track", &pi3_py_track);
  tree_->Branch("pi3_pz_track", &pi3_pz_track);

  tree_->Branch("B_chi2", &B_chi2);
  tree_->Branch("J_chi2", &J_chi2);
  tree_->Branch("psi2S_chi2", &psi2S_chi2);

  tree_->Branch("B_Prob",    &B_Prob);
  tree_->Branch("J_Prob",  &J_Prob);
  tree_->Branch("psi2S_Prob", &psi2S_Prob);       

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
  tree_->Branch("nTracksFromPV", &nTracksFromPV);

  tree_->Branch("vRefMuP", &vRefMuP);
  tree_->Branch("vRefMuM", &vRefMuM);
  tree_->Branch("vRefPi1", &vRefPi1);
  tree_->Branch("vRefPi2", &vRefPi2);
  tree_->Branch("vRefPi3", &vRefPi3);

  tree_->Branch("nVtx",       &nVtx);
  tree_->Branch("run",        &run,       "run/I");
  tree_->Branch("event",        &event,     "event/I");
  tree_->Branch("lumiblock",&lumiblock,"lumiblock/I");

  tree_->Branch("flightLen", &flightLen);
  tree_->Branch("flightLenErr", &flightLenErr);
  tree_->Branch("flightLenSig", &flightLenSig);

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

  tree_->Branch("tri_Dim25",&tri_Dim25);
  tree_->Branch("tri_JpsiTk",&tri_JpsiTk);
  tree_->Branch("tri_JpsiTkTk",&tri_JpsiTkTk); 
 
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
/*
  treeTest_->Branch("nJpsi_test",&nJpsi_test);
  treeTest_->Branch("mu1_px_test",&mu1_px_test);
  treeTest_->Branch("mu1_py_test",&mu1_py_test);
  treeTest_->Branch("mu1_pz_test",&mu1_pz_test);
  treeTest_->Branch("mu1_charge_test",&mu1_charge_test);
  treeTest_->Branch("mu2_px_test",&mu2_px_test);
  treeTest_->Branch("mu2_py_test",&mu2_py_test);
  treeTest_->Branch("mu2_pz_test",&mu2_pz_test);
  treeTest_->Branch("mu2_charge_test",&mu2_charge_test);
  treeTest_->Branch("Jpsi_dca_test",&Jpsi_dca_test);
  treeTest_->Branch("Jpsi_vx_test",&Jpsi_vx_test);
  treeTest_->Branch("Jpsi_vy_test",&Jpsi_vy_test);
  treeTest_->Branch("Jpsi_vz_test",&Jpsi_vz_test);
  treeTest_->Branch("Jpsi_mass_test",&Jpsi_mass_test);
  treeTest_->Branch("Jpsi_prob_test",&Jpsi_prob_test);
  treeTest_->Branch("Jpsi_chi2_test",&Jpsi_chi2_test);
*/
}


// ------------ method called once each job just after ending the event loop  ------------
void Psi2Spi::endJob() {
  tree_->GetDirectory()->cd();
  tree_->Write();
//  treeTest_->GetDirectory()->cd();
//  treeTest_->Write();
}

//define this as a plug-in
DEFINE_FWK_MODULE(Psi2Spi);

