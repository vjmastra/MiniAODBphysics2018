// -*- C++ -*-
//
// Package:    JPsiLambda
// Class:      JPsiLambda
// 
//=================================================
// original author:  Jhovanny Andres Mejia        |
//         created:  Fryday Sep 23                |
//         <jhovanny.andres.mejia.guisao@cern.ch> | 
//=================================================


// Lambda_b -> JPsi Lambda
// Lambda_N = B

// system include files
#include <memory>

// user include files
#include "myAnalyzers/JPsiKsPAT/src/JPsiLambda.h"

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
JPsiLambda::JPsiLambda(const edm::ParameterSet& iConfig)
  :
  dimuon_Label(consumes<edm::View<pat::Muon>>(iConfig.getParameter<edm::InputTag>("dimuons"))),
  trakCollection_label(consumes<edm::View<pat::PackedCandidate>>(iConfig.getParameter<edm::InputTag>("Trak"))),
  primaryVertices_Label(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("primaryVertices"))),
  BSLabel_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("bslabel"))),
  triggerResults_Label(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("TriggerResults"))),
  triggerObjects_(consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("objects"))),
  v0PtrCollection_(consumes<reco::VertexCompositePtrCandidateCollection>(iConfig.getParameter<edm::InputTag>("secundaryVerticesPtr"))),	       

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

  lBDecayVtxX(0), lBDecayVtxY(0), lBDecayVtxZ(0), lBDecayVtxXE(0), lBDecayVtxYE(0), lBDecayVtxZE(0),
  lBDecayVtxXYE(0), lBDecayVtxXZE(0), lBDecayVtxYZE(0),
  
  VDecayVtxX(0), VDecayVtxY(0), VDecayVtxZ(0), VDecayVtxXE(0), VDecayVtxYE(0), VDecayVtxZE(0),
  VDecayVtxXYE(0), VDecayVtxXZE(0), VDecayVtxYZE(0), 

  // *******************************************************
  nlB(0), nMu(0), nJpsi(0),
  lB_mass(0), lB_px(0), lB_py(0), lB_pz(0),
  
  lambda_mass(0), lambda_px(0), lambda_py(0), lambda_pz(0),
  lambda_pt1(0), lambda_px1(0), lambda_py1(0), lambda_pz1(0), 
  lambda_pt2(0), lambda_px2(0), lambda_py2(0), lambda_pz2(0), 

  lambda_px1_track(0), lambda_py1_track(0), lambda_pz1_track(0), 
  lambda_px2_track(0), lambda_py2_track(0), lambda_pz2_track(0), 

  pi1dxy(0), pi2dxy(0), pi1dz(0), pi2dz(0),
  pi1dxy_e(0), pi2dxy_e(0), pi1dz_e(0), pi2dz_e(0),
  lambda_charge1(0), lambda_charge2(0),

  J_mass(0), J_px(0), J_py(0), J_pz(0),
  J_pt1(0), J_px1(0), J_py1(0), J_pz1(0), 
  J_pt2(0), J_px2(0), J_py2(0), J_pz2(0), 
  J_charge1(0), J_charge2(0),

  Jtest_mass(0), Jtest_prob(0),

  lambda_chi2(0), J_chi2(0), lB_chi2(0),
  lB_Prob(0), J_Prob(0), lambda_Prob(0),

  run(0), event(0),
  lumiblock(0)

{
   //now do what ever initialization is needed
}

JPsiLambda::~JPsiLambda()
{

}

// ------------ method called to for each event  ------------
void JPsiLambda::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using std::vector;
  using namespace edm;
  using namespace reco;
  using namespace std;
  
  //*********************************
  // Get event content information
  //*********************************  

  edm::Handle<std::vector<reco::VertexCompositePtrCandidate>> theV0PtrHandle;
  iEvent.getByToken(v0PtrCollection_,  theV0PtrHandle);

// Kinematic fit
  edm::ESHandle<TransientTrackBuilder> theB; 
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB); 

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
      
      if ( theV0PtrHandle->size()>0 && thePATMuonHandle->size()>=2 ) { 
        for ( vector<VertexCompositePtrCandidate>::const_iterator iVee = theV0PtrHandle->begin();   iVee != theV0PtrHandle->end(); ++iVee ) {
         //get Lam tracks from V0 candidate
          vector<pat::PackedCandidate> v0daughters;
	  vector<Track> theDaughterTracks;
	  v0daughters.push_back( *(dynamic_cast<const pat::PackedCandidate *>(iVee->daughter(0))) );
	  v0daughters.push_back( *(dynamic_cast<const pat::PackedCandidate *>(iVee->daughter(1))) );
		     	     
          for(unsigned int j = 0; j < v0daughters.size(); ++j) {
            theDaughterTracks.push_back(v0daughters[j].pseudoTrack());
	  }
		     			
	  
	  //Now let's see if these two tracks make a vertex
	  reco::TransientTrack protonTT((*theB).build(theDaughterTracks[0]));
	  reco::TransientTrack pionTT((*theB).build(theDaughterTracks[1]));		     
		     
    	  ParticleMass pion_mass = 0.13957018;
          ParticleMass proton_mass = 0.93827208;
          ParticleMass lambda_particle_mass = 1.115683;
	  float pion_sigma = pion_mass*1.e-6;
          float proton_sigma = proton_mass*1.e-6;
	  float lambda_sigma = lambda_particle_mass*1.e-6;
		     
	  //initial chi2 and ndf before kinematic fits.
	  float chi = 0.;
	  float ndf = 0.;
	  vector<RefCountedKinematicParticle> v0dauParticles;

	  try {
	    v0dauParticles.push_back(pFactory.particle(protonTT,proton_mass,chi,ndf,proton_sigma));
	    v0dauParticles.push_back(pFactory.particle(pionTT,pion_mass,chi,ndf,pion_sigma));
	  }
	  catch(...) {
	    std::cout<<" Exception caught ... continuing 3 "<<std::endl;
	    continue;
	  }
		     
	  RefCountedKinematicTree lambdaVertexFitTree;
	  try{
	    lambdaVertexFitTree = fitter.fit(v0dauParticles); 
	  }
	  catch(...) {
	    std::cout<<" Exception caught ... continuing 4 "<<std::endl;                   
	    continue;
	  }
	  if (!lambdaVertexFitTree->isValid()) {
            //std::cout << "invalid vertex from the Ks0 vertex fit" << std::endl;
	    continue; 
	  }

	  lambdaVertexFitTree->movePointerToTheTop();	     
	  RefCountedKinematicParticle lambda_vFit_noMC = lambdaVertexFitTree->currentParticle();
	  RefCountedKinematicVertex lambda_vFit_vertex_noMC = lambdaVertexFitTree->currentDecayVertex();
	    
          if( lambda_vFit_vertex_noMC->chiSquared() < 0 ) { 
            //std::cout << "negative chisq from ks fit" << endl;
	    continue;
	  }
		     
  	  //some loose cuts go here
	    
	  //if(lambda_vFit_vertex_noMC->chiSquared()>50) continue;
	  if(lambda_vFit_noMC->currentState().mass()< 1.0 || lambda_vFit_noMC->currentState().mass()> 1.2) continue;
		     
	  lambdaVertexFitTree->movePointerToTheFirstChild();
	  RefCountedKinematicParticle T1CandMC = lambdaVertexFitTree->currentParticle();
		     
	  lambdaVertexFitTree->movePointerToTheNextChild();
	  RefCountedKinematicParticle T2CandMC = lambdaVertexFitTree->currentParticle();
		     
	  //  lambda  mass constrain
	  // do mass constrained vertex fit
	  // creating the constraint with a small sigma to put in the resulting covariance 
	  // matrix in order to avoid singularities
	  // JPsi mass constraint is applied in the final B fit
	   
	  KinematicParticleFitter csFitterlambda;
	  KinematicConstraint * lambda_c = new MassKinematicConstraint(lambda_particle_mass, lambda_sigma);
	  // add mass constraint to the lambda fit to do a constrained fit:  
		     
	  lambdaVertexFitTree = csFitterlambda.fit(lambda_c,lambdaVertexFitTree);
	  if (!lambdaVertexFitTree->isValid()){
	    //std::cout << "caught an exception in the ks mass constraint fit" << std::endl;
	    continue; 
	  }
		     
	  lambdaVertexFitTree->movePointerToTheTop();
	  RefCountedKinematicParticle lambda_vFit_withMC = lambdaVertexFitTree->currentParticle();
	     
	  //Now we are ready to combine!
		     
  	  vector<RefCountedKinematicParticle> vFitMCParticles2;

          try {
	    vFitMCParticles2.push_back(pFactory.particle(muon1TT,muon_mass,chi,ndf,muon_sigma));
	    vFitMCParticles2.push_back(pFactory.particle(muon2TT,muon_mass,chi,ndf,muon_sigma));
	    vFitMCParticles2.push_back(lambda_vFit_withMC);
          }
          catch (...) {
            std::cout << "Exception caught ... continuing 5" << std::endl;
            continue;
          }

          RefCountedKinematicTree lBVertexFitTree;
          try {
            lBVertexFitTree = fitter.fit(vFitMCParticles2);
          }
          catch (...) {
            std::cout<<" Exception caught ... continuing 6 "<<std::endl;
            continue;
          }

          if (!lBVertexFitTree->isValid()) {
            //std::cout << "caught an exception in the psi vertex fit" << std::endl;
            continue;
          }
	     
	  lBVertexFitTree->movePointerToTheTop();		     
	     
	  RefCountedKinematicParticle lBCandMC = lBVertexFitTree->currentParticle();
	  RefCountedKinematicVertex lBDecayVertexMC = lBVertexFitTree->currentDecayVertex();
	  if (!lBDecayVertexMC->vertexIsValid()){
	    //std::cout << "B MC fit vertex is not valid" << endl;
	    continue;
	  }
	     
	  if(lBCandMC->currentState().mass()< 4.0 || lBCandMC->currentState().mass()> 10.0) continue;
	     
	  if(lBDecayVertexMC->chiSquared()<0) {
            //std::cout << " continue from negative chi2 = " << bDecayVertexMC->chiSquared() << endl;
	  //  continue;
	  }
		     
	  double lB_Prob_tmp       = TMath::Prob(lBDecayVertexMC->chiSquared(),(int)lBDecayVertexMC->degreesOfFreedom());
	  if(lB_Prob_tmp<0.01) {
            continue;
	  }		     
		     
   	  // get children from final B fit
	  lBVertexFitTree->movePointerToTheFirstChild();
	  RefCountedKinematicParticle mu1CandMC = lBVertexFitTree->currentParticle();
	  lBVertexFitTree->movePointerToTheNextChild();
	  RefCountedKinematicParticle mu2CandMC = lBVertexFitTree->currentParticle();
	  lBVertexFitTree->movePointerToTheNextChild();
	  RefCountedKinematicParticle lambdaCandMC = lBVertexFitTree->currentParticle();
		   
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

	  GlobalVector lambdap1vec(T1CandMC->currentState().globalMomentum().x(),
	  T1CandMC->currentState().globalMomentum().y(),
 	  T1CandMC->currentState().globalMomentum().z());

          GlobalVector lambdap2vec(T2CandMC->currentState().globalMomentum().x(),
	  T2CandMC->currentState().globalMomentum().y(),
	  T2CandMC->currentState().globalMomentum().z());

	  KinematicParameters lambdaProtonKP = T1CandMC->currentState().kinematicParameters();
	  KinematicParameters lambdaPionKP = T2CandMC->currentState().kinematicParameters();
	  KinematicParameters lambdaDaupKP;
	  KinematicParameters lambdaDaumKP;
	       
	  if ( T1CandMC->currentState().particleCharge() > 0 ) lambdaDaupKP = lambdaProtonKP;
	  if ( T1CandMC->currentState().particleCharge() < 0 ) lambdaDaumKP = lambdaProtonKP;
	  if ( T2CandMC->currentState().particleCharge() > 0 ) lambdaDaupKP = lambdaPionKP;
	  if ( T2CandMC->currentState().particleCharge() < 0 ) lambdaDaumKP = lambdaPionKP;	 

	  // fill candidate variables now
	  
	  if(nlB==0){		    
	    nMu  = nMu_tmp;
	    // cout<< "*Number of Muons : " << nMu_tmp << endl;
	  } // end nB==0		     

	  lB_mass->push_back(lBCandMC->currentState().mass());
	  lB_px->push_back(lBCandMC->currentState().globalMomentum().x());
	  lB_py->push_back(lBCandMC->currentState().globalMomentum().y());
	  lB_pz->push_back(lBCandMC->currentState().globalMomentum().z());

          lambda_mass->push_back( lambda_vFit_noMC->currentState().mass() );
          lambda_px->push_back( lambda_vFit_noMC->currentState().globalMomentum().x() );
	  lambda_py->push_back( lambda_vFit_noMC->currentState().globalMomentum().y() );
	  lambda_pz->push_back( lambda_vFit_noMC->currentState().globalMomentum().z() );

          J_mass->push_back( J_vFit_noMC->currentState().mass() );
	  J_px->push_back( J_vFit_noMC->currentState().globalMomentum().x() );
	  J_py->push_back( J_vFit_noMC->currentState().globalMomentum().y() );
	  J_pz->push_back( J_vFit_noMC->currentState().globalMomentum().z() );

	  lambda_pt1->push_back(lambdap1vec.perp());
	  lambda_px1->push_back(lambdaProtonKP.momentum().x());
	  lambda_py1->push_back(lambdaProtonKP.momentum().y());
	  lambda_pz1->push_back(lambdaProtonKP.momentum().z());
	  lambda_px1_track->push_back(v0daughters[0].px());
	  lambda_py1_track->push_back(v0daughters[0].py());
	  lambda_pz1_track->push_back(v0daughters[0].pz());
	  lambda_charge1->push_back(T1CandMC->currentState().particleCharge());
 
          lambda_pt2->push_back(lambdap2vec.perp());
	  lambda_px2->push_back(lambdaPionKP.momentum().x());
	  lambda_py2->push_back(lambdaPionKP.momentum().y());
	  lambda_pz2->push_back(lambdaPionKP.momentum().z());
	  lambda_px2_track->push_back(v0daughters[1].px());
	  lambda_py2_track->push_back(v0daughters[1].py());
	  lambda_pz2_track->push_back(v0daughters[1].pz());
	  lambda_charge2->push_back(T2CandMC->currentState().particleCharge());
 
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
 
	  lambda_chi2->push_back(lambda_vFit_vertex_noMC->chiSquared());
	  J_chi2->push_back(J_vFit_vertex_noMC->chiSquared());
	  lB_chi2->push_back(lBDecayVertexMC->chiSquared());
 
          //double B_Prob_tmp       = TMath::Prob(bDecayVertexMC->chiSquared(),(int)bDecayVertexMC->degreesOfFreedom());
	  double J_Prob_tmp   = TMath::Prob(J_vFit_vertex_noMC->chiSquared(),(int)J_vFit_vertex_noMC->degreesOfFreedom());
	  double lambda_Prob_tmp  = TMath::Prob(lambda_vFit_vertex_noMC->chiSquared(),(int)lambda_vFit_vertex_noMC->degreesOfFreedom());
	  lB_Prob->push_back(lB_Prob_tmp);
          J_Prob->push_back(J_Prob_tmp);
          lambda_Prob->push_back(lambda_Prob_tmp);
 
	  lBDecayVtxX->push_back((*lBDecayVertexMC).position().x());
	  lBDecayVtxY->push_back((*lBDecayVertexMC).position().y());
	  lBDecayVtxZ->push_back((*lBDecayVertexMC).position().z());
	  lBDecayVtxXE->push_back(lBDecayVertexMC->error().cxx());
	  lBDecayVtxYE->push_back(lBDecayVertexMC->error().cyy());
	  lBDecayVtxZE->push_back(lBDecayVertexMC->error().czz());
	  lBDecayVtxXYE->push_back(lBDecayVertexMC->error().cyx());
	  lBDecayVtxXZE->push_back(lBDecayVertexMC->error().czx());
	  lBDecayVtxYZE->push_back(lBDecayVertexMC->error().czy());

          VDecayVtxX->push_back( lambda_vFit_vertex_noMC->position().x() );
          VDecayVtxY->push_back( lambda_vFit_vertex_noMC->position().y() );
          VDecayVtxZ->push_back( lambda_vFit_vertex_noMC->position().z() );
          VDecayVtxXE->push_back( lambda_vFit_vertex_noMC->error().cxx() );
          VDecayVtxYE->push_back( lambda_vFit_vertex_noMC->error().cyy() );
	  VDecayVtxZE->push_back( lambda_vFit_vertex_noMC->error().czz() );
	  VDecayVtxXYE->push_back( lambda_vFit_vertex_noMC->error().cyx() );
	  VDecayVtxXZE->push_back( lambda_vFit_vertex_noMC->error().czx() );
	  VDecayVtxYZE->push_back( lambda_vFit_vertex_noMC->error().czy() );

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
 
	  pi1dxy->push_back(v0daughters[0].dxy());
	  pi2dxy->push_back(v0daughters[1].dxy());
	  pi1dz->push_back(v0daughters[0].dz());
	  pi2dz->push_back(v0daughters[1].dz());
 
	  pi1dxy_e->push_back(v0daughters[0].dxyError());
	  pi2dxy_e->push_back(v0daughters[1].dxyError());
	  pi1dz_e->push_back(v0daughters[0].dzError());
	  pi2dz_e->push_back(v0daughters[1].dzError());

	  // try refitting the primary without the tracks in the B reco candidate		   
		  
	  nlB++;	       
		   
	  /////////////////////////////////////////////////
	  v0dauParticles.clear();
          v0daughters.clear();
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
   if (nlB > 0 ) {
     tree_->Fill();
   }
   // *********

   nlB = 0; nMu = 0; nJpsi = 0;

   lB_mass->clear();    lB_px->clear();    lB_py->clear();    lB_pz->clear();
   lambda_mass->clear(); lambda_px->clear(); lambda_py->clear(); lambda_pz->clear();

   J_mass->clear();  J_px->clear();  J_py->clear();  J_pz->clear();

   lambda_pt1->clear(); lambda_px1->clear(); lambda_py1->clear(); lambda_pz1->clear(); lambda_charge1->clear(); 
   lambda_pt2->clear(); lambda_px2->clear(); lambda_py2->clear(); lambda_pz2->clear(); lambda_charge2->clear(); 

   lambda_px1_track->clear(); lambda_py1_track->clear(); lambda_pz1_track->clear(); 
   lambda_px2_track->clear(); lambda_py2_track->clear(); lambda_pz2_track->clear(); 

   J_pt1->clear();  J_px1->clear();  J_py1->clear();  J_pz1->clear(), J_charge1->clear();
   J_pt2->clear();  J_px2->clear();  J_py2->clear();  J_pz2->clear(), J_charge2->clear();

   lambda_chi2->clear(); J_chi2->clear(); lB_chi2->clear();
   lB_Prob->clear(); J_Prob->clear(); lambda_Prob->clear();

   Jtest_mass->clear(); Jtest_prob->clear();

   // *********

   nVtx = 0;
   priVtxX = 0; priVtxY = 0; priVtxZ = 0; 
   priVtxXE = 0; priVtxYE = 0; priVtxZE = 0; priVtxCL = 0;
   priVtxXYE = 0;   priVtxXZE = 0;   priVtxYZE = 0;

   lBDecayVtxX->clear(); lBDecayVtxY->clear(); lBDecayVtxZ->clear(); 
   lBDecayVtxXE->clear(); lBDecayVtxYE->clear(); lBDecayVtxZE->clear(); 
   lBDecayVtxXYE->clear(); lBDecayVtxXZE->clear(); lBDecayVtxYZE->clear();  

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
bool JPsiLambda::IsTheSame(const pat::PackedCandidate& tk, const pat::Muon& mu){
  double DeltaEta = fabs(mu.eta()-tk.eta());
  double DeltaP   = fabs(mu.p()-tk.p());
  if (DeltaEta < 0.02 && DeltaP < 0.02) return true;
  return false;
}

// ------------ method called once each job just before starting event loop  ------------

void JPsiLambda::beginJob()
{

  std::cout << "Beginning analyzer job with value of doMC = " << doMC_ << std::endl;

  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>("ntuple","Lambdab->J/psi Lambda ntuple");
  treeTest_ = fs->make<TTree>("test","Jpsi->mumu");

  tree_->Branch("nlB",&nlB,"nlB/i");
  tree_->Branch("nMu",&nMu,"nMu/i");

  tree_->Branch("lB_mass", &lB_mass);
  tree_->Branch("lB_px", &lB_px);
  tree_->Branch("lB_py", &lB_py);
  tree_->Branch("lB_pz", &lB_pz);

//  tree_->Branch("PsiPiPi_mass", &psiPiPi_mass);  

  tree_->Branch("lambda_mass", &lambda_mass);
  tree_->Branch("lambda_px", &lambda_px);
  tree_->Branch("lambda_py", &lambda_py);
  tree_->Branch("lambda_pz", &lambda_pz);
 
  tree_->Branch("J_mass", &J_mass);
  tree_->Branch("J_px", &J_px);
  tree_->Branch("J_py", &J_py);
  tree_->Branch("J_pz", &J_pz);

  tree_->Branch("lambda_pt1", &lambda_pt1);
  tree_->Branch("lambda_px1", &lambda_px1);
  tree_->Branch("lambda_py1", &lambda_py1);
  tree_->Branch("lambda_pz1", &lambda_pz1);
  tree_->Branch("lambda_px1_track", &lambda_px1_track);
  tree_->Branch("lambda_py1_track", &lambda_py1_track);
  tree_->Branch("lambda_pz1_track", &lambda_pz1_track);
  tree_->Branch("lambda_charge1", &lambda_charge1); 
 
  tree_->Branch("lambda_pt2", &lambda_pt2);
  tree_->Branch("lambda_px2", &lambda_px2);
  tree_->Branch("lambda_py2", &lambda_py2);
  tree_->Branch("lambda_pz2", &lambda_pz2);
  tree_->Branch("lambda_px2_track", &lambda_px2_track);
  tree_->Branch("lambda_py2_track", &lambda_py2_track);
  tree_->Branch("lambda_pz2_track", &lambda_pz2_track);
  tree_->Branch("lambda_charge2", &lambda_charge2);

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

  tree_->Branch("lB_chi2", &lB_chi2);
  tree_->Branch("lambda_chi2", &lambda_chi2);
  tree_->Branch("J_chi2", &J_chi2);

  tree_->Branch("lB_Prob",    &lB_Prob);
  tree_->Branch("lambda_Prob", &lambda_Prob);
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

  tree_->Branch("lBDecayVtxX",&lBDecayVtxX);
  tree_->Branch("lBDecayVtxY",&lBDecayVtxY);
  tree_->Branch("lBDecayVtxZ",&lBDecayVtxZ);
  tree_->Branch("lBDecayVtxXE",&lBDecayVtxXE);
  tree_->Branch("lBDecayVtxYE",&lBDecayVtxYE);
  tree_->Branch("lBDecayVtxZE",&lBDecayVtxZE);
  tree_->Branch("lBDecayVtxXYE",&lBDecayVtxXYE);
  tree_->Branch("lBDecayVtxXZE",&lBDecayVtxXZE);
  tree_->Branch("lBDecayVtxYZE",&lBDecayVtxYZE);

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
void JPsiLambda::endJob() {
  tree_->GetDirectory()->cd();
  tree_->Write();
  treeTest_->GetDirectory()->cd();
  treeTest_->Write();
}

//define this as a plug-in
DEFINE_FWK_MODULE(JPsiLambda);

