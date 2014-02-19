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
	EleTag_ = consumes<reco::GsfElectronCollection>(ps.getParameter<edm::InputTag>("EleTag")); 
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
	bcEB_PV->setAxisTitle("N_{pv}",1); 
	bcEB_PV->setAxisTitle("Basic Clusters", 2); 
	
	prof_name="bcEE_PV"; 
	title="Basic Clusters EE vs. PV"; 
	bcEE_PV=ibooker.bookProfile(prof_name,title,50,0.,50.,50,0.,350.);
	bcEE_PV->setAxisTitle("N_{pv}",1); 
	bcEE_PV->setAxisTitle("Basic Clusters", 2); 
	
	prof_name="scEB_PV"; 
	title="Super Clusters EB vs. PV"; 
	scEB_PV=ibooker.bookProfile(prof_name,title,50,0.,50.,50,0.,350.);
	scEB_PV->setAxisTitle("N_{pv}",1);
	scEB_PV->setAxisTitle("Super Clusters", 2); 
	
	prof_name="scEE_PV"; 
	title="Super Clusters EE vs. PV"; 
	scEE_PV=ibooker.bookProfile(prof_name,title,50,0.,50.,50,0.,350.);
	scEE_PV->setAxisTitle("N_{pv}",1);
	scEE_PV->setAxisTitle("Super Clusters", 2); 
	
	prof_name="scEtEB_PV"; 
	title="Super Clusters Et EB vs. PV"; 
	scEtEB_PV=ibooker.bookProfile(prof_name,title,50,0.,50.,50,0.,350.);
	scEtEB_PV->setAxisTitle("N_{pv}",1);
	scEtEB_PV->setAxisTitle("Super Cluster E_{T} [GeV]", 2); 
	
	
	prof_name="scEtEE_PV"; 
	title="Super Clusters Et EE vs. PV"; 
	scEtEE_PV=ibooker.bookProfile(prof_name,title,50,0.,50.,50,0.,350.);
	scEtEE_PV->setAxisTitle("N_{pv}",1);
	scEtEE_PV->setAxisTitle("Super Cluster E_{T} [GeV]", 2); 
	
	prof_name="recHitEtEB_PV"; 
	title="Reconstructed Hit Et EB vs. PV"; 
	recHitEtEB_PV=ibooker.bookProfile(prof_name,title,50,0.,50.,50,0.,350.);
	recHitEtEB_PV->setAxisTitle("N_{pv}",1);
	recHitEtEB_PV->setAxisTitle("Reconstructed hit E_{T} [GeV]", 2); 

	prof_name="recHitEtEE_PV"; 
	title="Reconstructed Hit Et EE vs. PV"; 
	recHitEtEE_PV=ibooker.bookProfile(prof_name,title,50,0.,50.,50,0.,350.);
	recHitEtEE_PV->setAxisTitle("N_{pv}",1);
	recHitEtEE_PV->setAxisTitle("Reconstructed hit E_{T} [GeV]", 2); 
	
	prof_name="emIso_PV"; 
	title="EM Isolation vs. PV"; 
	emIso_PV=ibooker.bookProfile(prof_name,title,50,0.,50.,50,0.,350);
	emIso_PV->setAxisTitle("N_{pv}",1);
	emIso_PV->setAxisTitle("EM_{Isolation} [GeV]", 2); 
	
	//construct histograms 
	
	prof_name="emIso"; 
	title="EM Isolation"; 
	emIso=ibooker.book1D(prof_name,title,50,0,50);
	emIso->setAxisTitle("EM_{Isolation} [GeV]",1); 
	emIso->setAxisTitle("Events",2); 
	
	prof_name="scHitEtEB"; 
	title="Super Cluster Hit Et EB"; 
	scHitEtEB=ibooker.book1D(prof_name,title,50,0,350);
	scHitEtEB->setAxisTitle("super cluster hit E_{T} [GeV]",1); 
	scHitEtEB->setAxisTitle("Events",2); 
	
	prof_name="scHitEtEE"; 
	title="Super Cluster Hit Et EE"; 
	scHitEtEE=ibooker.book1D(prof_name,title,50,0,350);
	scHitEtEE->setAxisTitle("super cluster hit E_{T} [GeV]",1); 
	scHitEtEE->setAxisTitle("Events",2); 

	
	prof_name="scHitE_EB"; 
	title="Super Cluster Hit E EB"; 
	scHitE_EB=ibooker.book1D(prof_name,title,50,0,350);
	scHitE_EB->setAxisTitle("super cluster hit E [GeV]",1);
	scHitE_EB->setAxisTitle("Events",2); 

	
	prof_name="scHitE_EE"; 
	title="Super Cluster Hit E EE"; 
	scHitE_EE=ibooker.book1D(prof_name,title,50,0,350);
	scHitE_EE->setAxisTitle("super cluster hit E [GeV]",1);
	scHitE_EE->setAxisTitle("Events",2); 
	
	//SC eta
	prof_name="scEta_EB"; 
	title="Super Cluster #eta EB"; 
	scEta_EB=ibooker.book1D(prof_name,title,50,-6,6);
	scEta_EB->setAxisTitle("#eta",1); 
	scEta_EB->setAxisTitle("Events",2); 
	
	prof_name="scEta_EE"; 
	title="Super Cluster #eta EE"; 
	scEta_EE=ibooker.book1D(prof_name,title,50,-6,6);
	scEta_EE->setAxisTitle("#eta",1); 
	scEta_EE->setAxisTitle("Events",2); 
	
	//SC phi
	prof_name="scPhi_EB"; 
	title="Super Cluster #phi EB"; 
	scPhi_EB=ibooker.book1D(prof_name,title,50,-3.14,3.14);
	scPhi_EB->setAxisTitle("super cluster #phi",1); 
	scPhi_EB->setAxisTitle("Events",2); 
	
	prof_name="scPhi_EE"; 
	title="Super Cluster #phi EE"; 
	scPhi_EE=ibooker.book1D(prof_name,title,50,-3.14,3.14);
	scPhi_EE->setAxisTitle("super cluster #phi",1); 
	scPhi_EE->setAxisTitle("Events",2); 
	
	//sc sigma eta eta / eta phi
	
	prof_name="scSigmaIetaIeta_EB"; 
	title="Super Cluster sigmaIetaIeta EB"; 
	scSigmaIetaIeta_EB=ibooker.book1D(prof_name,title,50,0,0.002);
	scSigmaIetaIeta_EB->setAxisTitle("#sigma_{i#etai#eta}",1); 
	scSigmaIetaIeta_EB->setAxisTitle("Events",2); 
	
	prof_name="scSigmaIetaIeta_EE"; 
	title="Super Cluster sigmaIetaIeta EE"; 
	scSigmaIetaIeta_EE=ibooker.book1D(prof_name,title,50,0,0.01);
	scSigmaIetaIeta_EE->setAxisTitle("#sigma_{i#etai#eta}",1); 
	scSigmaIetaIeta_EE->setAxisTitle("Events",2); 
	
	//phi
	prof_name="scSigmaIetaIphi_EB"; 
	title="Super Cluster sigmaIetaIphi EB"; 
	scSigmaIetaIphi_EB=ibooker.book1D(prof_name,title,50,0,0.002);
	scSigmaIetaIphi_EB->setAxisTitle("#sigma_{i#etai#phi}",1); 
	scSigmaIetaIphi_EB->setAxisTitle("Events",2); 
	
	prof_name="scSigmaIetaIphi_EE"; 
	title="Super Cluster sigmaIetaIphi EE"; 
	scSigmaIetaIphi_EE=ibooker.book1D(prof_name,title,50,0,0.01);
	scSigmaIetaIphi_EE->setAxisTitle("#sigma_{i#etai#phi}",1); 
	scSigmaIetaIphi_EE->setAxisTitle("Events",2); 
	
	//R9
	prof_name="r9_EB"; 
	title="r9 EB"; 
	r9_EB=ibooker.book1D(prof_name,title,50,0,1.5);
	r9_EB->setAxisTitle("R_{9}",1);
	r9_EB->setAxisTitle("Events",2); 
	
	prof_name="r9_EE"; 
	title="r9 EE"; 
	r9_EE=ibooker.book1D(prof_name,title,50,0,1.5);
	r9_EE->setAxisTitle("R_{9}",1);
	r9_EE->setAxisTitle("Events",2); 

	
	//Rec Hit
	
	prof_name="recHitEtEB"; 
	title="Reconstructed Cluster Hit Et EB"; 
	recHitEtEB=ibooker.book1D(prof_name,title,50,0,400);
	recHitEtEB->setAxisTitle("Reconstructed Hit E_{T} [GeV]",1);
	recHitEtEB->setAxisTitle("Events",2); 
	
	prof_name="recHitEtEE"; 
	title="Reconstructed Cluster Hit Et EE"; 
	recHitEtEE=ibooker.book1D(prof_name,title,50,0,400);
	recHitEtEE->setAxisTitle("Reconstructed Hit E_{T} [GeV]",1);
	recHitEtEE->setAxisTitle("Events",2); 

	
}


void PileUpDep::analyze(const edm::Event& e, const edm::EventSetup& c){

	//Analaysis Code
	
	//geometry

	if(!initCaloGeometry_){
		c.get<CaloGeometryRecord>().get(geomH);
		c.get<CaloTopologyRecord>().get(caloTop);
		
		initCaloGeometry_=true;
	}
	
	//Vertex collection: 
	//-----------------------------------------
	edm::Handle<reco::VertexCollection> PVCollection_h;
	e.getByToken(VertexCollection_,PVCollection_h);
	if ( ! PVCollection_h.isValid() ) {
		edm::LogWarning("VertexCollection") << "VertexCollection not found"; 
	}
	//-----------------gsfElectrons -------------------------
	edm::Handle<reco::GsfElectronCollection> electronCollection_h; 
	e.getByToken(EleTag_, electronCollection_h); 
	if( !electronCollection_h.isValid()){
		edm::LogWarning("EBRecoSummary") << "Electrons not found"; 
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

	//----------------- Reconstructed Hit Ecal barrel 
	
	edm::Handle<EcalRecHitCollection> RecHitsEB;
	e.getByToken( RecHitCollection_EB_,RecHitsEB );
	if ( ! RecHitsEB.isValid() ) {
		edm::LogWarning("EBRecoSummary") << "RecHitsEB not found"; 
	}
	
	//----------------- Reconstructed Hit Ecal Endcap  
	
	edm::Handle<EcalRecHitCollection> RecHitsEE;
	e.getByToken( RecHitCollection_EE_, RecHitsEE );
	if ( ! RecHitsEE.isValid() ) {
		edm::LogWarning("EBRecoSummary") << "RecHitsEB not found"; 
	}
	//----------------- Super Cluster Collection Ecal Endcap  ---------
	
	edm::Handle<reco::SuperClusterCollection> superClusters_EE_h;
	e.getByToken( superClusterCollection_EE_, superClusters_EE_h );
	if ( ! superClusters_EE_h.isValid() ) {
		edm::LogWarning("EERecoSummary") << "superClusters_EE_h not found"; 
	}
	//--------- Fill Isolation -----------------
	double IsoEcal=0; 
	
	if(electronCollection_h.isValid()){
		for (reco::GsfElectronCollection::const_iterator recoElectron =
			 electronCollection_h->begin ();
			 recoElectron != electronCollection_h->end (); recoElectron++)
		{
			//std::cout << "EM Iso: " << recoElectron->dr03EcalRecHitSumEt() << std::endl; 
			IsoEcal +=recoElectron->dr03EcalRecHitSumEt();///recoElectron->et()
		}		
		emIso_PV->Fill(PVCollection_h->size(),IsoEcal);
		emIso->Fill(IsoEcal);
	}

	//fill super clusters EE
	scEE_PV->Fill(PVCollection_h->size(), superClusters_EE_h->size()); 

	double scEE_Et=0; //EE Et (Transverse energy)
	double scEE_E=0; // EE E (energy)
	
	for (reco::SuperClusterCollection::const_iterator itSC = superClusters_EE_h->begin(); 
		 itSC != superClusters_EE_h->end(); ++itSC ) {
		scEE_Et+= itSC -> energy() * sin(2.*atan( exp(- itSC->position().eta() )));
		scEE_E+=itSC->energy();
		
		//fill super cluster endcap eta/phi
		scEta_EE->Fill(itSC->position().eta());
		scPhi_EE->Fill(itSC->position().phi());

		//get sigma eta_eta etc
		
		CaloTopology const* p_topology = caloTop.product();//get calo topology
		const EcalRecHitCollection* eeRecHits = RecHitsEE.product();
		
		reco::BasicCluster const& seedCluster(*itSC->seed());
		std::vector<float> cov = EcalClusterTools::localCovariances(seedCluster, eeRecHits, p_topology);
		float sigmaIetaIeta = cov[0];
		float sigmaIetaIphi = cov[1];

		
		float e3x3 = EcalClusterTools::e3x3(seedCluster, eeRecHits, p_topology);
		float r9 = e3x3 / itSC->energy();
		
		r9_EE->Fill(r9); 
		scSigmaIetaIeta_EE->Fill(sigmaIetaIeta);
		scSigmaIetaIphi_EE->Fill(sigmaIetaIphi);
		
		//std::cout  << " sigmaIetaIeta: " << sigmaIetaIeta << std::endl; 

		
	}//sc-EE loop
	
	scEtEE_PV->Fill(PVCollection_h->size(),scEE_Et); 
	scHitEtEE->Fill(scEE_Et); //super cluster Et historam 
	scHitE_EE->Fill(scEE_E); //super cluster energy histogram
	
	//----------------- Super Cluster Collection Ecal Barrel  ---------

	edm::Handle<reco::SuperClusterCollection> superClusters_EB_h;
	e.getByToken( superClusterCollection_EB_, superClusters_EB_h );
	if ( ! superClusters_EB_h.isValid() ) {
		edm::LogWarning("EBRecoSummary") << "superClusters_EB_h not found"; 
	}
	scEB_PV->Fill(PVCollection_h->size(), superClusters_EB_h->size()); 

	double scEB_Et=0; //EB Et 
	double scEB_E=0; //EB E 

	
	for (reco::SuperClusterCollection::const_iterator itSC = superClusters_EB_h->begin(); 
		 itSC != superClusters_EB_h->end(); ++itSC ) {
		scEB_Et+= itSC -> energy() * sin(2.*atan( exp(- itSC->position().eta() )));
		scEB_E+= itSC->energy();
		
		//fill super cluster Barrel eta/phi
		scEta_EB->Fill(itSC->position().eta());
		scPhi_EB->Fill(itSC->position().phi());
		
		//sigma ietaieta etc 
		
		CaloTopology const* p_topology = caloTop.product();//get calo topology
		const EcalRecHitCollection* ebRecHits = RecHitsEB.product();

		
		reco::BasicCluster const& seedCluster(*itSC->seed());
		std::vector<float> cov = EcalClusterTools::localCovariances(seedCluster, ebRecHits, p_topology);
		float sigmaIetaIeta = cov[0];
		float sigmaIetaIphi = cov[1];

		
		float e3x3 = EcalClusterTools::e3x3(seedCluster, ebRecHits, p_topology);
		float r9 = e3x3 / itSC->energy();
		
		r9_EB->Fill(r9);
		scSigmaIetaIeta_EB->Fill(sigmaIetaIeta);
		scSigmaIetaIphi_EB->Fill(sigmaIetaIphi);
		
	}//sc-EB loop
	
	scEtEB_PV->Fill(PVCollection_h->size(),scEB_Et); 
	scHitEtEB->Fill(scEB_Et); 
	scHitE_EB->Fill(scEB_E); 
	

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

