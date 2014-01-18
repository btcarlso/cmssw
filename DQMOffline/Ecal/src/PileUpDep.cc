/*
 * \file PileUpDep.cc
 * \author Ben Carlson - CMU
 * Last Update:
 *
 */

#include "DQMOffline/Ecal/interface/PileUpDep.h"
#include "DataFormats/Scalers/interface/DcsStatus.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GtFdlWord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"

// Framework

const static int XBINS=2000;

PileUpDep::PileUpDep(const edm::ParameterSet& ps)

{
	//initialize parameters
	
	 initCaloGeometry_ = false;
	
	VertexCollection_ = consumes<reco::VertexCollection>(ps.getParameter<edm::InputTag>("VertexCollection"));

	
	basicClusterCollection_EE_ = consumes<reco::BasicClusterCollection>(ps.getParameter<edm::InputTag>("basicClusterCollection_EE"));
	basicClusterCollection_EB_ = consumes<reco::BasicClusterCollection>(ps.getParameter<edm::InputTag>("basicClusterCollection_EB"));
	superClusterCollection_EB_ = consumes<reco::SuperClusterCollection>(ps.getParameter<edm::InputTag>("superClusterCollection_EB"));
	superClusterCollection_EE_ = consumes<reco::SuperClusterCollection>(ps.getParameter<edm::InputTag>("superClusterCollection_EE"));
	
	RecHitCollection_EB_       = consumes<EcalRecHitCollection>(ps.getParameter<edm::InputTag>("RecHitCollection_EB"));
	RecHitCollection_EE_       = consumes<EcalRecHitCollection>(ps.getParameter<edm::InputTag>("RecHitCollection_EE"));

	//print_=ps.getUntrackedParameter<bool>("print"); 

}

PileUpDep::~PileUpDep(){
}

void PileUpDep::bookHistograms(DQMStore::IBooker & ibooker,
                                edm::Run const & /* iRun */,
                                edm::EventSetup const & /* iSetup */) {

  // Fetch GlobalTag information and fill the string/ME.
   // initialize
	ibooker.cd();
	ibooker.setCurrentFolder("Ecal/PileUpDep");

	std::string prof_name="bcEB_PV"; 
	std::string title="Basic clusters EB vs. PV";
	bcEB_PV=ibooker.bookProfile(prof_name,title,50,0.,50.,50,0.,350.);
	
	prof_name="bcEE_PV"; 
	title="Basic Clusters EE vs. PV"; 
	bcEE_PV=ibooker.bookProfile(prof_name,title,50,0.,50.,50,0.,350.);
	
	prof_name="scEB_PV"; 
	title="Super Clusters EB vs. PV"; 
	scEB_PV=ibooker.bookProfile(prof_name,title,50,0.,50.,50,0.,350.);
	
	prof_name="scEE_PV"; 
	title="Super Clusters EE vs. PV"; 
	scEE_PV=ibooker.bookProfile(prof_name,title,50,0.,50.,50,0.,350.);

	prof_name="scEtEB_PV"; 
	title="Super Clusters Et EB vs. PV"; 
	scEtEB_PV=ibooker.bookProfile(prof_name,title,50,0.,50.,50,0.,350.);
	
	prof_name="scEtEE_PV"; 
	title="Super Clusters Et EE vs. PV"; 
	scEtEE_PV=ibooker.bookProfile(prof_name,title,50,0.,50.,50,0.,350.);
	
	prof_name="recHitEtEB_PV"; 
	title="Reconstructed Hit Et EB vs. PV"; 
	recHitEtEB_PV=ibooker.bookProfile(prof_name,title,50,0.,50.,50,0.,350.);

	prof_name="recHitEtEE_PV"; 
	title="Reconstructed Hit Et EE vs. PV"; 
	recHitEtEE_PV=ibooker.bookProfile(prof_name,title,50,0.,50.,50,0.,350.);
	
	//construct histograms 
	
	prof_name="scHitEtEB"; 
	title="Super Cluster Hit Et EB"; 
	scHitEtEB=ibooker.book1D(prof_name,title,50,0,350);
	
	prof_name="scHitEtEE"; 
	title="Super Cluster Hit Et EE"; 
	scHitEtEE=ibooker.book1D(prof_name,title,50,0,350);
	
	prof_name="recHitEtEB"; 
	title="Reconstructed Cluster Hit Et EB"; 
	recHitEtEB=ibooker.book1D(prof_name,title,50,0,350);
	
	prof_name="recHitEtEE"; 
	title="Reconstructed Cluster Hit Et EE"; 
	recHitEtEE=ibooker.book1D(prof_name,title,50,0,350);
	
}


void PileUpDep::analyze(const edm::Event& e, const edm::EventSetup& c){

	//Analaysis Code
	
	//geometry

	if(!initCaloGeometry_){
		c.get<CaloGeometryRecord>().get(geomH);
		initCaloGeometry_=true;
	}
	
	//Vertex collection: 
	//-----------------------------------------
	edm::Handle<reco::VertexCollection> PVCollection_h;
	e.getByToken(VertexCollection_,PVCollection_h);
	if ( ! PVCollection_h.isValid() ) {
		edm::LogWarning("VertexCollection") << "VertexCollection not found"; 
	}
	
	//----------------- Basic Cluster Collection Ecal Barrel  ---------
	edm::Handle<reco::BasicClusterCollection> basicClusters_EB_h;
	e.getByToken( basicClusterCollection_EB_, basicClusters_EB_h );
//	const reco::BasicClusterCollection* theBarrelBasicClusters = basicClusters_EB_h.product () ;
	if ( ! basicClusters_EB_h.isValid() ) {
		edm::LogWarning("EBRecoSummary") << "basicClusters_EB_h not found"; 
	}
	//edm::LogInfo("PileUpOut") << "BasicClusters: " << basicClusters_EB_h->size(); 
	
	bcEB_PV->Fill(PVCollection_h->size(), basicClusters_EB_h->size()); 
	
	//----------------- Basic Cluster Collection Ecal Endcal  ---------

	edm::Handle<reco::BasicClusterCollection> basicClusters_EE_h;
	e.getByToken( basicClusterCollection_EE_, basicClusters_EE_h );
	if ( ! basicClusters_EE_h.isValid() ) {
		edm::LogWarning("EERecoSummary") << "basicClusters_EE_h not found"; 
	}

	bcEE_PV->Fill(PVCollection_h->size(), basicClusters_EE_h->size()); 

	//----------------- Super Cluster Collection Ecal Endcap  ---------
	
	edm::Handle<reco::SuperClusterCollection> superClusters_EE_h;
	e.getByToken( superClusterCollection_EE_, superClusters_EE_h );
	if ( ! superClusters_EE_h.isValid() ) {
		edm::LogWarning("EERecoSummary") << "superClusters_EE_h not found"; 
	}
	scEE_PV->Fill(PVCollection_h->size(), superClusters_EE_h->size()); 

	double scEE_Et=0; //EE Et
	
	for (reco::SuperClusterCollection::const_iterator itSC = superClusters_EE_h->begin(); 
		 itSC != superClusters_EE_h->end(); ++itSC ) {
		scEE_Et+= itSC -> energy() * sin(2.*atan( exp(- itSC->position().eta() )));
		
	}//sc-EE loop
	scEtEE_PV->Fill(PVCollection_h->size(),scEE_Et); 
	scHitEtEE->Fill(scEE_Et); 
	
	//----------------- Super Cluster Collection Ecal Barrel  ---------

	edm::Handle<reco::SuperClusterCollection> superClusters_EB_h;
	e.getByToken( superClusterCollection_EB_, superClusters_EB_h );
	if ( ! superClusters_EB_h.isValid() ) {
		edm::LogWarning("EBRecoSummary") << "superClusters_EB_h not found"; 
	}
	scEB_PV->Fill(PVCollection_h->size(), superClusters_EB_h->size()); 

	double scEB_Et=0; //EB Et 
	
	for (reco::SuperClusterCollection::const_iterator itSC = superClusters_EB_h->begin(); 
		 itSC != superClusters_EB_h->end(); ++itSC ) {
		scEB_Et+= itSC -> energy() * sin(2.*atan( exp(- itSC->position().eta() )));
		
	}//sc-EB loop
	
	scEtEB_PV->Fill(PVCollection_h->size(),scEB_Et); 
	scHitEtEB->Fill(scEB_Et); 

	
	//----------------- Reconstructed Hit Ecal barrel 
	
	edm::Handle<EcalRecHitCollection> RecHitsEB;
	e.getByToken( RecHitCollection_EB_,RecHitsEB );
	if ( ! RecHitsEB.isValid() ) {
		edm::LogWarning("EBRecoSummary") << "RecHitsEB not found"; 
	}

	//-------------------Compute scalar sum of reconstructed hit Et
	double RecHitEt_EB=0; 

	for ( EcalRecHitCollection::const_iterator itr = RecHitsEB->begin () ;
		 itr != RecHitsEB->end () ;++itr)
    {	
		//RecHitEt_EB +=itr->energy();
		
		GlobalPoint const& position  = geomH->getGeometry(itr->detid())->getPosition();
		RecHitEt_EB += itr -> energy() * sin(position.theta()) ;
	}//EB Rec Hit
	
	recHitEtEB->Fill(RecHitEt_EB); 
	recHitEtEB_PV->Fill(PVCollection_h->size(),RecHitEt_EB); 

	
	//----------------- Reconstructed Hit Ecal Endcap  

	edm::Handle<EcalRecHitCollection> RecHitsEE;
	e.getByToken( RecHitCollection_EE_, RecHitsEE );
	if ( ! RecHitsEE.isValid() ) {
		edm::LogWarning("EBRecoSummary") << "RecHitsEB not found"; 
	}
	//-------------------Compute scalar sum of reconstructed hit Et
	double RecHitEt_EE=0; 

	for ( EcalRecHitCollection::const_iterator itr = RecHitsEE->begin () ;
		 itr != RecHitsEE->end () ;++itr)
    {
		GlobalPoint const& position  = geomH->getGeometry(itr->detid())->getPosition();
		RecHitEt_EE += itr -> energy() * sin(position.theta()) ;
	}//EB Rec Hit
	
	recHitEtEE->Fill(RecHitEt_EE); 
	recHitEtEE_PV->Fill(PVCollection_h->size(),RecHitEt_EE); 
													   
													   
  return;
}

void
PileUpDep::endLuminosityBlock(const edm::LuminosityBlock& l, const edm::EventSetup& c)
{
  // int nlumi = l.id().luminosityBlock();

  // fill dcs vs lumi
  /* set those bins 0 for which bits are ON
     needed for merge off lumi histograms across files */

  return;
}

DEFINE_FWK_MODULE(PileUpDep);

