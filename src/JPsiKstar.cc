// -*- C++ -*-
//
// Package:    JPsiKstar
// Class:      JPsiKstar
// 
//=================================================
// original author:  Jhovanny Andres Mejia        |
//         created:  Fryday Sep 23                |
//         <jhovanny.andres.mejia.guisao@cern.ch> | 
//=================================================


// B0 -> JPsiKstar -> (mumu)(kpi)

// system include files
#include <memory>

// user include files
#include "myAnalyzers/JPsiKsPAT/src/JPsiKstar.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
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
JPsiKstar::JPsiKstar(const edm::ParameterSet& iConfig)
  :
  dimuon_Label(consumes<edm::View<pat::Muon>>(iConfig.getParameter<edm::InputTag>("dimuons"))),
  trakCollection_label(consumes<edm::View<pat::PackedCandidate>>(iConfig.getParameter<edm::InputTag>("Trak"))),
  primaryVertices_Label(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("primaryVertices"))),
  BSLabel_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("bslabel"))),
  builderToken_(esConsumes<TransientTrackBuilder, TransientTrackRecord>(edm::ESInputTag("", "TransientTrackBuilder"))),
  triggerResults_Label(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("TriggerResults"))),
  triggerObjects_(consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("objects"))),
//  v0PtrCollection_(consumes<reco::VertexCompositePtrCandidateCollection>(iConfig.getParameter<edm::InputTag>("secundaryVerticesPtr"))),	       

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
 
  // ************************ ****************************************************

  BDecayVtxX(0), BDecayVtxY(0), BDecayVtxZ(0), BDecayVtxXE(0), BDecayVtxYE(0), BDecayVtxZE(0),
  BDecayVtxXYE(0), BDecayVtxXZE(0), BDecayVtxYZE(0),
  
  VDecayVtxX(0), VDecayVtxY(0), VDecayVtxZ(0), VDecayVtxXE(0), VDecayVtxYE(0), VDecayVtxZE(0),
  VDecayVtxXYE(0), VDecayVtxXZE(0), VDecayVtxYZE(0), 

  // *******************************************************
  nB(0), nMu(0), nJpsi(0),
  B_mass(0), B_px(0), B_py(0), B_pz(0),
  
  kstar_mass(0), kstar_px(0), kstar_py(0), kstar_pz(0),
  kstar_pt1(0), kstar_px1(0), kstar_py1(0), kstar_pz1(0), 
  kstar_pt2(0), kstar_px2(0), kstar_py2(0), kstar_pz2(0), 

  pi1dxy(0), pi2dxy(0), pi1dz(0), pi2dz(0),
  pi1dxy_e(0), pi2dxy_e(0), pi1dz_e(0), pi2dz_e(0),
  kstar_charge1(0), kstar_charge2(0),

  J_mass(0), J_px(0), J_py(0), J_pz(0),
  J_pt1(0), J_px1(0), J_py1(0), J_pz1(0), 
  J_pt2(0), J_px2(0), J_py2(0), J_pz2(0), 
  J_charge1(0), J_charge2(0),

  Jtest_mass(0), Jtest_prob(0),

  kstar_chi2(0), J_chi2(0), B_chi2(0),
  B_Prob(0), J_Prob(0), kstar_Prob(0),

  run(0), event(0),
  lumiblock(0)

{
   //now do what ever initialization is needed
}

JPsiKstar::~JPsiKstar()
{

}

// ------------ method called to for each event  ------------
void JPsiKstar::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using std::vector;
  using namespace edm;
  using namespace reco;
  using namespace std;
  
  //*********************************
  // Get event content information
  //*********************************  

//  edm::Handle<std::vector<reco::VertexCompositePtrCandidate>> theV0PtrHandle;
//  iEvent.getByToken(v0PtrCollection_,  theV0PtrHandle);

  // Kinematic fit

  auto const &theB = iSetup.getData(builderToken_);

  edm::Handle< View<pat::PackedCandidate> > thePATTrackHandle;
  iEvent.getByToken(trakCollection_label,thePATTrackHandle);

  edm::Handle< View<pat::Muon> > thePATMuonHandle;
  iEvent.getByToken(dimuon_Label,thePATMuonHandle);

  edm::Handle<edm::TriggerResults> triggerBits;
  edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;

  iEvent.getByToken(triggerResults_Label, triggerBits);
  iEvent.getByToken(triggerObjects_, triggerObjects); 
  const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);

//  for (pat::TriggerObjectStandAlone obj : *triggerObjects) { // note: not "const &" since we want to call unpackPathNames


  for ( size_t iTrigObj = 0; iTrigObj < triggerObjects->size(); ++iTrigObj ) {
    pat::TriggerObjectStandAlone obj( triggerObjects->at( iTrigObj ) );
    obj.unpackPathNames(names);
    obj.unpackFilterLabels(iEvent, *triggerBits);
  }
  //*********************************
  //Now we get the primary vertex 
  //*********************************

  reco::Vertex bestVtx;
  reco::Vertex bestVtxBS;

  // get primary vertex
  edm::Handle<std::vector<reco::Vertex> > recVtxs;
  iEvent.getByToken(primaryVertices_Label, recVtxs);
  bestVtx = *(recVtxs->begin());
  
  priVtxX = bestVtx.x();
  priVtxY = bestVtx.y();
  priVtxZ = bestVtx.z();
  priVtxXE = bestVtx.covariance(0, 0);
  priVtxYE = bestVtx.covariance(1, 1);
  priVtxZE = bestVtx.covariance(2, 2);
  priVtxXYE = bestVtx.covariance(0, 1);
  priVtxXZE = bestVtx.covariance(0, 2);
  priVtxYZE = bestVtx.covariance(1, 2);
  priVtxCL = ChiSquaredProbability((double)(bestVtx.chi2()),(double)(bestVtx.ndof())); 

  nVtx = recVtxs->size();

  lumiblock = iEvent.id().luminosityBlock();
  run = iEvent.id().run();
  event = iEvent.id().event(); 


  //Control loop on Jpsi

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

      reco::TransientTrack muonATT((theB).build(muonTrackA));
      reco::TransientTrack muonBTT((theB).build(muonTrackB));

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

      RefCountedKinematicTree Jtest_VertexFitTree;
      try {
        Jtest_VertexFitTree = fitter.fit(muonParticles_test);
      }
      catch (...) {
        std::cout<<" Exception caught ... continuing 2 "<<std::endl;
        continue;
      }

      if (!Jtest_VertexFitTree->isValid()) {
        continue;
      }

      Jtest_VertexFitTree->movePointerToTheTop();

      RefCountedKinematicParticle Jtest_vFit_noMC = Jtest_VertexFitTree->currentParticle();
      RefCountedKinematicVertex Jtest_vFit_vertex_noMC = Jtest_VertexFitTree->currentDecayVertex();

      if(Jtest_vFit_vertex_noMC->chiSquared() < 0 ) {
        std::cout << "negative chisq from psi fit" << endl;
        continue;
      }

      double Jtest_Prob_tmp   = TMath::Prob(Jtest_vFit_vertex_noMC->chiSquared(),(int)Jtest_vFit_vertex_noMC->degreesOfFreedom());
      if(Jtest_Prob_tmp<0.01) {
        continue;
      }

      if(Jtest_vFit_noMC->currentState().mass()<3.0 || Jtest_vFit_noMC->currentState().mass()>3.2) continue;

      nJpsi++;
      Jtest_mass->push_back(Jtest_vFit_noMC->currentState().mass());
      Jtest_prob->push_back(Jtest_Prob_tmp);

      muonParticles_test.clear();

    }
  }

  //*****************************************
  //Let's begin by looking for J/psi->mu+mu-

  unsigned int nMu_tmp = thePATMuonHandle->size();
 
  for(View<pat::Muon>::const_iterator iMuon1 = thePATMuonHandle->begin(); iMuon1 != thePATMuonHandle->end(); ++iMuon1) { 
    for(View<pat::Muon>::const_iterator iMuon2 = iMuon1+1; iMuon2 != thePATMuonHandle->end(); ++iMuon2) {  
      
      if(iMuon1==iMuon2) continue;
	  
      //opposite charge 
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
      float dca = fabs( cApp.distance() );	  
      if (dca < 0. || dca > 0.5) continue;
      //cout<<" closest approach  "<<dca<<endl;

      // *****  end DCA for the 2 muons *********************

      //The mass of a muon and the insignificant mass sigma 
      //to avoid singularities in the covariance matrix.
      ParticleMass muon_mass = 0.10565837; //pdg mass
      ParticleMass psi_mass = 3.096916;
      float muon_sigma = muon_mass*1.e-6;
      //float psi_sigma = psi_mass*1.e-6;

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
	  
      RefCountedKinematicParticle J_vFit_noMC = JVertexFitTree->currentParticle();//masa del J/psi
      RefCountedKinematicVertex J_vFit_vertex_noMC = JVertexFitTree->currentDecayVertex();//vertice del J/psi
	  
      if( J_vFit_vertex_noMC->chiSquared() < 0 ) {
        std::cout << "negative chisq from psi fit" << endl;
        continue;
      }

      double J_Prob_tmp   = TMath::Prob(J_vFit_vertex_noMC->chiSquared(),(int)J_vFit_vertex_noMC->degreesOfFreedom());
      if(J_Prob_tmp<0.01) {
        continue;
      }
	  
      //some loose cuts go here

      if(J_vFit_noMC->currentState().mass()<3.0 || J_vFit_noMC->currentState().mass()>3.2) continue;
  
      for(View<pat::PackedCandidate>::const_iterator iTrack1 = thePATTrackHandle->begin(); iTrack1 != thePATTrackHandle->end(); ++iTrack1 ) {

        if(iTrack1->charge()==0) continue;
        if(fabs(iTrack1->pdgId())!=211) continue;
        if(iTrack1->pt()<0.95) continue;
        if(!(iTrack1->trackHighPurity())) continue;

        for(View<pat::PackedCandidate>::const_iterator iTrack2 = thePATTrackHandle->begin(); iTrack2 != thePATTrackHandle->end(); ++iTrack2 ) {
        //i check for ALL the combinations since the analysis is asymmetric wrt the two kstar daughters
          if(iTrack1==iTrack2) continue;
          if(iTrack2->charge()==0) continue;
          if(fabs(iTrack2->pdgId())!=211) continue;
          if(iTrack2->pt()<0.95) continue;
          if(!(iTrack2->trackHighPurity())) continue;
          if(iTrack1->charge() == iTrack2->charge()) continue;

          if ( IsTheSame(*iTrack1,*iMuon1) || IsTheSame(*iTrack1,*iMuon2) ) continue;
          if ( IsTheSame(*iTrack2,*iMuon1) || IsTheSame(*iTrack2,*iMuon2) ) continue;

	  //Now let's see if these two tracks make a vertex
	  reco::TransientTrack pion1TT((theB).build(iTrack1->pseudoTrack()));
	  reco::TransientTrack pion2TT((theB).build(iTrack2->pseudoTrack()));		     
		     
    	  ParticleMass pion_mass = 0.13957018;
          ParticleMass kaon_mass = 0.493677;
	  float pion_sigma = pion_mass*1.e-6;
          float kaon_sigma = kaon_mass*1.e-6;
		     
	  //initial chi2 and ndf before kinematic fits.
	  float chi = 0.;
	  float ndf = 0.;
	  vector<RefCountedKinematicParticle> v0dauParticles;

	  try {
	    v0dauParticles.push_back(pFactory.particle(pion1TT,kaon_mass,chi,ndf,kaon_sigma));
	    v0dauParticles.push_back(pFactory.particle(pion2TT,pion_mass,chi,ndf,pion_sigma));
	  }
	  catch(...) {
	    std::cout<<" Exception caught ... continuing 3 "<<std::endl;
	    continue;
	  }
		     
	  RefCountedKinematicTree kstarVertexFitTree;
	  try{
	    kstarVertexFitTree = fitter.fit(v0dauParticles); 
	  }
	  catch(...) {
	    std::cout<<" Exception caught ... continuing 4 "<<std::endl;                   
	    continue;
	  }
	  if (!kstarVertexFitTree->isValid()) {
            //std::cout << "invalid vertex from the Ks0 vertex fit" << std::endl;
	    continue; 
	  }

	  kstarVertexFitTree->movePointerToTheTop();	     
	  RefCountedKinematicParticle kstar_vFit_noMC = kstarVertexFitTree->currentParticle();
	  RefCountedKinematicVertex kstar_vFit_vertex_noMC = kstarVertexFitTree->currentDecayVertex();
	    
          if( kstar_vFit_vertex_noMC->chiSquared() < 0 ) { 
            //std::cout << "negative chisq from ks fit" << endl;
	    continue;
	  }

          double kstar_Prob_tmp   = TMath::Prob(kstar_vFit_vertex_noMC->chiSquared(),(int)kstar_vFit_vertex_noMC->degreesOfFreedom());
          if(kstar_Prob_tmp<0.01) {
            continue;
          }
		     
  	  //some loose cuts go here
	 
	  if(kstar_vFit_noMC->currentState().mass()<0.5 || kstar_vFit_noMC->currentState().mass()>2.5) continue;
		     
	  kstarVertexFitTree->movePointerToTheFirstChild();
	  RefCountedKinematicParticle T1CandMC = kstarVertexFitTree->currentParticle();
		     
	  kstarVertexFitTree->movePointerToTheNextChild();
	  RefCountedKinematicParticle T2CandMC = kstarVertexFitTree->currentParticle();
		       
	  //Now we are ready to combine!
		     
  	  vector<RefCountedKinematicParticle> vFitMCParticles2;

          try {
	    vFitMCParticles2.push_back(pFactory.particle(muon1TT,muon_mass,chi,ndf,muon_sigma));
	    vFitMCParticles2.push_back(pFactory.particle(muon2TT,muon_mass,chi,ndf,muon_sigma));
	    vFitMCParticles2.push_back(pFactory.particle(pion1TT,kaon_mass,chi,ndf,kaon_sigma));
            vFitMCParticles2.push_back(pFactory.particle(pion2TT,pion_mass,chi,ndf,pion_sigma));
          }
          catch (...) {
            std::cout << "Exception caught ... continuing 5" << std::endl;
            continue;
          }

          RefCountedKinematicTree BVertexFitTree;
          try {
            BVertexFitTree = fitter.fit(vFitMCParticles2);
          }
          catch (...) {
            std::cout<<" Exception caught ... continuing 6 "<<std::endl;
            continue;
          }

          if (!BVertexFitTree->isValid()) {
            //std::cout << "caught an exception in the psi vertex fit" << std::endl;
            continue;
          }
	     
	  BVertexFitTree->movePointerToTheTop();		     
	     
	  RefCountedKinematicParticle BCandMC = BVertexFitTree->currentParticle();
	  RefCountedKinematicVertex BDecayVertexMC = BVertexFitTree->currentDecayVertex();
	  if (!BDecayVertexMC->vertexIsValid()){
	    //std::cout << "B MC fit vertex is not valid" << endl;
	    continue;
	  }
	     
	  if(BCandMC->currentState().mass()< 4.0 || BCandMC->currentState().mass()> 10.0) continue;
	     
	  if(BDecayVertexMC->chiSquared()<0) {
            //std::cout << " continue from negative chi2 = " << bDecayVertexMC->chiSquared() << endl;
	  //  continue;
	  }
		     
	  double B_Prob_tmp       = TMath::Prob(BDecayVertexMC->chiSquared(),(int)BDecayVertexMC->degreesOfFreedom());
	  if(B_Prob_tmp<0.01) {
            continue;
	  }		     
		     
   	  // get children from final B fit
	  BVertexFitTree->movePointerToTheFirstChild();
	  RefCountedKinematicParticle mu1CandMC = BVertexFitTree->currentParticle();
	  BVertexFitTree->movePointerToTheNextChild();
	  RefCountedKinematicParticle mu2CandMC = BVertexFitTree->currentParticle();
	  BVertexFitTree->movePointerToTheNextChild();
	  RefCountedKinematicParticle pi1CandMC = BVertexFitTree->currentParticle();
          BVertexFitTree->movePointerToTheNextChild();
          RefCountedKinematicParticle pi2CandMC = BVertexFitTree->currentParticle();
		   
	  KinematicParameters psiMu1KP = mu1CandMC->currentState().kinematicParameters();
	  KinematicParameters psiMu2KP = mu2CandMC->currentState().kinematicParameters();
	  KinematicParameters psiMupKP;
	  KinematicParameters psiMumKP;
	       
	  if ( mu1CandMC->currentState().particleCharge() > 0 ) psiMupKP = psiMu1KP;
	  if ( mu1CandMC->currentState().particleCharge() < 0 ) psiMumKP = psiMu1KP;
	  if ( mu2CandMC->currentState().particleCharge() > 0 ) psiMupKP = psiMu2KP;
          if ( mu2CandMC->currentState().particleCharge() < 0 ) psiMumKP = psiMu2KP;	 

 	  GlobalVector Jp1vec(mu1CandMC->currentState().globalMomentum().x(),
	  mu1CandMC->currentState().globalMomentum().y(),
 	  mu1CandMC->currentState().globalMomentum().z());

          GlobalVector Jp2vec(mu2CandMC->currentState().globalMomentum().x(),
	  mu2CandMC->currentState().globalMomentum().y(),
 	  mu2CandMC->currentState().globalMomentum().z());

	  GlobalVector kstarp1vec(pi1CandMC->currentState().globalMomentum().x(),
	  pi1CandMC->currentState().globalMomentum().y(),
 	  pi1CandMC->currentState().globalMomentum().z());

          GlobalVector kstarp2vec(pi2CandMC->currentState().globalMomentum().x(),
	  pi2CandMC->currentState().globalMomentum().y(),
	  pi2CandMC->currentState().globalMomentum().z());

	  KinematicParameters kstarKaon = pi1CandMC->currentState().kinematicParameters();
	  KinematicParameters kstarPion = pi2CandMC->currentState().kinematicParameters();
	  KinematicParameters kstarDaup;
	  KinematicParameters kstarDaum;
	       
	  if ( pi1CandMC->currentState().particleCharge() > 0 ) kstarDaup = kstarKaon;
	  if ( pi1CandMC->currentState().particleCharge() < 0 ) kstarDaum = kstarKaon;
	  if ( pi2CandMC->currentState().particleCharge() > 0 ) kstarDaup = kstarPion;
	  if ( pi2CandMC->currentState().particleCharge() < 0 ) kstarDaum = kstarPion;	 

	  // fill candidate variables now
	  
	  if(nB==0){		    
	    nMu  = nMu_tmp;
	    // cout<< "*Number of Muons : " << nMu_tmp << endl;
	  } // end nB==0		     

	  B_mass->push_back(BCandMC->currentState().mass());
	  B_px->push_back(BCandMC->currentState().globalMomentum().x());
	  B_py->push_back(BCandMC->currentState().globalMomentum().y());
	  B_pz->push_back(BCandMC->currentState().globalMomentum().z());

          kstar_mass->push_back( kstar_vFit_noMC->currentState().mass() );
          kstar_px->push_back( kstar_vFit_noMC->currentState().globalMomentum().x() );
	  kstar_py->push_back( kstar_vFit_noMC->currentState().globalMomentum().y() );
	  kstar_pz->push_back( kstar_vFit_noMC->currentState().globalMomentum().z() );

          J_mass->push_back( J_vFit_noMC->currentState().mass() );
	  J_px->push_back( J_vFit_noMC->currentState().globalMomentum().x() );
	  J_py->push_back( J_vFit_noMC->currentState().globalMomentum().y() );
	  J_pz->push_back( J_vFit_noMC->currentState().globalMomentum().z() );

	  kstar_pt1->push_back(kstarp1vec.perp());
	  kstar_px1->push_back(kstarKaon.momentum().x());
	  kstar_py1->push_back(kstarKaon.momentum().y());
	  kstar_pz1->push_back(kstarKaon.momentum().z());
	  kstar_charge1->push_back(pi1CandMC->currentState().particleCharge());
 
          kstar_pt2->push_back(kstarp2vec.perp());
	  kstar_px2->push_back(kstarPion.momentum().x());
	  kstar_py2->push_back(kstarPion.momentum().y());
	  kstar_pz2->push_back(kstarPion.momentum().z());
	  kstar_charge2->push_back(pi2CandMC->currentState().particleCharge());
 
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
 
	  kstar_chi2->push_back(kstar_vFit_vertex_noMC->chiSquared());
	  J_chi2->push_back(J_vFit_vertex_noMC->chiSquared());
	  B_chi2->push_back(BDecayVertexMC->chiSquared());
 
          //double B_Prob_tmp       = TMath::Prob(bDecayVertexMC->chiSquared(),(int)bDecayVertexMC->degreesOfFreedom());
	  double J_Prob_tmp   = TMath::Prob(J_vFit_vertex_noMC->chiSquared(),(int)J_vFit_vertex_noMC->degreesOfFreedom());
	  //double kstar_Prob_tmp  = TMath::Prob(kstar_vFit_vertex_noMC->chiSquared(),(int)kstar_vFit_vertex_noMC->degreesOfFreedom());
	  B_Prob->push_back(B_Prob_tmp);
          J_Prob->push_back(J_Prob_tmp);
          kstar_Prob->push_back(kstar_Prob_tmp);
 
	  BDecayVtxX->push_back((*BDecayVertexMC).position().x());
	  BDecayVtxY->push_back((*BDecayVertexMC).position().y());
	  BDecayVtxZ->push_back((*BDecayVertexMC).position().z());
	  BDecayVtxXE->push_back(BDecayVertexMC->error().cxx());
	  BDecayVtxYE->push_back(BDecayVertexMC->error().cyy());
	  BDecayVtxZE->push_back(BDecayVertexMC->error().czz());
	  BDecayVtxXYE->push_back(BDecayVertexMC->error().cyx());
	  BDecayVtxXZE->push_back(BDecayVertexMC->error().czx());
	  BDecayVtxYZE->push_back(BDecayVertexMC->error().czy());

          VDecayVtxX->push_back( kstar_vFit_vertex_noMC->position().x() );
          VDecayVtxY->push_back( kstar_vFit_vertex_noMC->position().y() );
          VDecayVtxZ->push_back( kstar_vFit_vertex_noMC->position().z() );
          VDecayVtxXE->push_back( kstar_vFit_vertex_noMC->error().cxx() );
          VDecayVtxYE->push_back( kstar_vFit_vertex_noMC->error().cyy() );
	  VDecayVtxZE->push_back( kstar_vFit_vertex_noMC->error().czz() );
	  VDecayVtxXYE->push_back( kstar_vFit_vertex_noMC->error().cyx() );
	  VDecayVtxXZE->push_back( kstar_vFit_vertex_noMC->error().czx() );
	  VDecayVtxYZE->push_back( kstar_vFit_vertex_noMC->error().czy() );

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
 
	  pi1dxy->push_back(iTrack1->dxy());
	  pi2dxy->push_back(iTrack2->dxy());
	  pi1dz->push_back(iTrack1->dz());
	  pi2dz->push_back(iTrack2->dz());
 
	  pi1dxy_e->push_back(iTrack1->dxyError());
	  pi2dxy_e->push_back(iTrack2->dxyError());
	  pi1dz_e->push_back(iTrack1->dzError());
	  pi2dz_e->push_back(iTrack2->dzError());

	  // try refitting the primary without the tracks in the B reco candidate		   
		  
	  nB++;	       
		   
	  /////////////////////////////////////////////////
	  v0dauParticles.clear();
	  muonParticles.clear();
          vFitMCParticles2.clear();
	   
         }//for V0
       }//if V0
     }//muon2
   }//muon1

 
   if (nJpsi > 0) {
     treeTest_->Fill();
   }

   //fill the tree and clear the vectors
   if (nB > 0 ) {
     tree_->Fill();
   }
   // *********

   nB = 0; nMu = 0; nJpsi = 0;

   B_mass->clear();    B_px->clear();    B_py->clear();    B_pz->clear();
   kstar_mass->clear(); kstar_px->clear(); kstar_py->clear(); kstar_pz->clear();

   J_mass->clear();  J_px->clear();  J_py->clear();  J_pz->clear();

   kstar_pt1->clear(); kstar_px1->clear(); kstar_py1->clear(); kstar_pz1->clear(); kstar_charge1->clear(); 
   kstar_pt2->clear(); kstar_px2->clear(); kstar_py2->clear(); kstar_pz2->clear(); kstar_charge2->clear(); 

   J_pt1->clear();  J_px1->clear();  J_py1->clear();  J_pz1->clear(), J_charge1->clear();
   J_pt2->clear();  J_px2->clear();  J_py2->clear();  J_pz2->clear(), J_charge2->clear();

   kstar_chi2->clear(); J_chi2->clear(); B_chi2->clear();
   B_Prob->clear(); J_Prob->clear(); kstar_Prob->clear();

   Jtest_mass->clear(); Jtest_prob->clear();

   // *********

   nVtx = 0;
   priVtxX = 0; priVtxY = 0; priVtxZ = 0; 
   priVtxXE = 0; priVtxYE = 0; priVtxZE = 0; priVtxCL = 0;
   priVtxXYE = 0;   priVtxXZE = 0;   priVtxYZE = 0;

   BDecayVtxX->clear(); BDecayVtxY->clear(); BDecayVtxZ->clear(); 
   BDecayVtxXE->clear(); BDecayVtxYE->clear(); BDecayVtxZE->clear(); 
   BDecayVtxXYE->clear(); BDecayVtxXZE->clear(); BDecayVtxYZE->clear();  

   VDecayVtxX->clear(); VDecayVtxY->clear(); VDecayVtxZ->clear();
   VDecayVtxXE->clear(); VDecayVtxYE->clear(); VDecayVtxZE->clear();
   VDecayVtxXYE->clear(); VDecayVtxXZE->clear(); VDecayVtxYZE->clear();
  
   pi1dxy->clear(); pi2dxy->clear(); pi1dz->clear(); pi2dz->clear();
   pi1dxy_e->clear(); pi2dxy_e->clear(); pi1dz_e->clear(); pi2dz_e->clear();

   mumC2->clear();
   mumNHits->clear(); mumNPHits->clear();
   mupC2->clear();
   mupNHits->clear(); mupNPHits->clear();
   mumdxy->clear(); mupdxy->clear(); mumdz->clear(); mupdz->clear(); muon_dca->clear();

   tri_Dim25->clear(); tri_JpsiTk->clear(); tri_JpsiTkTk->clear();
      
   mu1soft->clear(); mu2soft->clear(); mu1tight->clear(); mu2tight->clear();
   mu1PF->clear(); mu2PF->clear(); mu1loose->clear(); mu2loose->clear(); 

}
/*
bool JPsiLambda::IsTheSame(const reco::Track& tk, const pat::Muon& mu){
  double DeltaEta = fabs(mu.eta()-tk.eta());
  double DeltaP   = fabs(mu.p()-tk.p());
  if (DeltaEta < 0.02 && DeltaP < 0.02) return true;
  return false;
}
*/
bool JPsiKstar::IsTheSame(const pat::PackedCandidate& tk, const pat::Muon& mu){
  double DeltaEta = fabs(mu.eta()-tk.eta());
  double DeltaP   = fabs(mu.p()-tk.p());
  if (DeltaEta < 0.02 && DeltaP < 0.02) return true;
  return false;
}

// ------------ method called once each job just before starting event loop  ------------

void JPsiKstar::beginJob()
{

  std::cout << "Beginning analyzer job with value of doMC = " << doMC_ << std::endl;

  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>("ntuple","B0->J/psi Kstar");
  treeTest_ = fs->make<TTree>("test","Jpsi->mumu");

  tree_->Branch("nB",&nB,"nB/i");
  tree_->Branch("nMu",&nMu,"nMu/i");

  tree_->Branch("B_mass", &B_mass);
  tree_->Branch("B_px", &B_px);
  tree_->Branch("B_py", &B_py);
  tree_->Branch("B_pz", &B_pz);

  tree_->Branch("kstar_mass", &kstar_mass);
  tree_->Branch("kstar_px", &kstar_px);
  tree_->Branch("kstar_py", &kstar_py);
  tree_->Branch("kstar_pz", &kstar_pz);
 
  tree_->Branch("J_mass", &J_mass);
  tree_->Branch("J_px", &J_px);
  tree_->Branch("J_py", &J_py);
  tree_->Branch("J_pz", &J_pz);

  tree_->Branch("kstar_pt1", &kstar_pt1);
  tree_->Branch("kstar_px1", &kstar_px1);
  tree_->Branch("kstar_py1", &kstar_py1);
  tree_->Branch("kstar_pz1", &kstar_pz1);
  tree_->Branch("kstar_charge1", &kstar_charge1); 
 
  tree_->Branch("kstar_pt2", &kstar_pt2);
  tree_->Branch("kstar_px2", &kstar_px2);
  tree_->Branch("kstar_py2", &kstar_py2);
  tree_->Branch("kstar_pz2", &kstar_pz2);
  tree_->Branch("kstar_charge2", &kstar_charge2);

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

  tree_->Branch("B_chi2", &B_chi2);
  tree_->Branch("kstar_chi2", &kstar_chi2);
  tree_->Branch("J_chi2", &J_chi2);

  tree_->Branch("B_Prob",    &B_Prob);
  tree_->Branch("kstar_Prob", &kstar_Prob);
  tree_->Branch("J_Prob",  &J_Prob);

  // *************************

  tree_->Branch("priVtxX",&priVtxX, "priVtxX/f");
  tree_->Branch("priVtxY",&priVtxY, "priVtxY/f");
  tree_->Branch("priVtxZ",&priVtxZ, "priVtxZ/f");
  tree_->Branch("priVtxXE",&priVtxXE, "priVtxXE/f");
  tree_->Branch("priVtxYE",&priVtxYE, "priVtxYE/f");
  tree_->Branch("priVtxZE",&priVtxZE, "priVtxZE/f");
  tree_->Branch("priVtxXYE",&priVtxXYE, "priVtxXYE/f");
  tree_->Branch("priVtxXZE",&priVtxXZE, "priVtxXZE/f");
  tree_->Branch("priVtxYZE",&priVtxYZE, "priVtxYZE/f");
  tree_->Branch("priVtxCL",&priVtxCL, "priVtxCL/f");

  tree_->Branch("nVtx",       &nVtx);
  tree_->Branch("run",        &run,       "run/I");
  tree_->Branch("event",        &event,     "event/I");
  tree_->Branch("lumiblock",&lumiblock,"lumiblock/I");

  tree_->Branch("BDecayVtxX",&BDecayVtxX);
  tree_->Branch("BDecayVtxY",&BDecayVtxY);
  tree_->Branch("BDecayVtxZ",&BDecayVtxZ);
  tree_->Branch("BDecayVtxXE",&BDecayVtxXE);
  tree_->Branch("BDecayVtxYE",&BDecayVtxYE);
  tree_->Branch("BDecayVtxZE",&BDecayVtxZE);
  tree_->Branch("BDecayVtxXYE",&BDecayVtxXYE);
  tree_->Branch("BDecayVtxXZE",&BDecayVtxXZE);
  tree_->Branch("BDecayVtxYZE",&BDecayVtxYZE);

  tree_->Branch("VDecayVtxX",&VDecayVtxX);
  tree_->Branch("VDecayVtxY",&VDecayVtxY);
  tree_->Branch("VDecayVtxZ",&VDecayVtxZ);
  tree_->Branch("VDecayVtxXE",&VDecayVtxXE);
  tree_->Branch("VDecayVtxYE",&VDecayVtxYE);
  tree_->Branch("VDecayVtxZE",&VDecayVtxZE);
  tree_->Branch("VDecayVtxXYE",&VDecayVtxXYE);
  tree_->Branch("VDecayVtxXZE",&VDecayVtxXZE);
  tree_->Branch("VDecayVtxYZE",&VDecayVtxYZE);

  tree_->Branch("pi1dxy",&pi1dxy);
  tree_->Branch("pi2dxy",&pi2dxy);
  tree_->Branch("pi1dz",&pi1dz);
  tree_->Branch("pi2dz",&pi2dz);

  tree_->Branch("pi1dxy_e",&pi1dxy_e);
  tree_->Branch("pi2dxy_e",&pi2dxy_e);
  tree_->Branch("pi1dz_e",&pi1dz_e);
  tree_->Branch("pi2dz_e",&pi2dz_e);

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

  treeTest_->Branch("Jtest_mass",&Jtest_mass);
  treeTest_->Branch("Jtest_prob",&Jtest_prob);
}


// ------------ method called once each job just after ending the event loop  ------------
void JPsiKstar::endJob() {
  tree_->GetDirectory()->cd();
  tree_->Write();
  treeTest_->GetDirectory()->cd();
  treeTest_->Write();
}

//define this as a plug-in
DEFINE_FWK_MODULE(JPsiKstar);

