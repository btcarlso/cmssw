#ifndef PileUpDep_H
#define PileUpDep_H

/*
 * \file PileUpDep.h
 *
 * \author Ben Carlson - CMU
 *
*/

#include <FWCore/Framework/interface/EDAnalyzer.h>
#include <FWCore/Framework/interface/Event.h>
#include <FWCore/Framework/interface/Run.h>
#include <FWCore/Framework/interface/MakerMacros.h>
#include <FWCore/ParameterSet/interface/ParameterSet.h>
#include <FWCore/ParameterSet/interface/Registry.h>
#include <FWCore/ServiceRegistry/interface/Service.h>
#include <FWCore/MessageLogger/interface/MessageLogger.h>

#include <DQMServices/Core/interface/DQMStore.h>
#include <DQMServices/Core/interface/MonitorElement.h>
#include <DQMServices/Core/interface/DQMEDAnalyzer.h>

#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/ESDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/EgammaReco/interface/PreshowerCluster.h"
#include "DataFormats/EgammaReco/interface/PreshowerClusterFwd.h"
#include "DataFormats/EgammaReco/interface/BasicClusterFwd.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalTools.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/EcalSeverityLevelAlgo.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/EcalCleaningAlgo.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalRecHitLess.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackExtra.h"
#include "DataFormats/GsfTrackReco/interface/GsfTangent.h"
#include "DataFormats/EgammaCandidates/interface/ElectronFwd.h"
#include "DataFormats/EgammaCandidates/interface/Electron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"

#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"

#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

//DataFormats
#include "DataFormats/Scalers/interface/DcsStatus.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"


class PileUpDep: public DQMEDAnalyzer{

public:

  /// Constructor
  PileUpDep(const edm::ParameterSet& ps);

  /// Destructor
  virtual ~PileUpDep();

protected:

  /// Analyze
  void analyze(const edm::Event& e, const edm::EventSetup& c);
  void bookHistograms(DQMStore::IBooker &, edm::Run const &, edm::EventSetup const &) override;
  void endLuminosityBlock(const edm::LuminosityBlock& l, const edm::EventSetup& c);

	edm::EDGetTokenT<reco::VertexCollection> VertexCollection_; //vertex collection
	
	edm::EDGetTokenT<reco::BasicClusterCollection> basicClusterCollection_EB_; // Ecal Barrel Basic Clusters 
	edm::EDGetTokenT<reco::BasicClusterCollection> basicClusterCollection_EE_;

	edm::EDGetTokenT<reco::SuperClusterCollection> superClusterCollection_EB_;
	edm::EDGetTokenT<reco::SuperClusterCollection> superClusterCollection_EE_;
	edm::EDGetTokenT<reco::GsfElectronCollection> EleTag_;
	
	edm::EDGetTokenT<EcalRecHitCollection> RecHitCollection_EB_;
	edm::EDGetTokenT<EcalRecHitCollection> RecHitCollection_EE_;
	

private:
	
	//profiles
	MonitorElement * bcEB_PV; //basic clusters Ecal-Barrel vs Number of Primary Vertices 
	MonitorElement * bcEE_PV;
	MonitorElement * scEB_PV;
	MonitorElement * scEE_PV;
	
	MonitorElement * scEtEB_PV;//super cluster Et profiles vs Number of vertices 
	MonitorElement * scEtEE_PV;
	
	MonitorElement * recHitEtEB_PV; // reconstructed hit Et profiles vs number of vertices
	MonitorElement * recHitEtEE_PV;

	MonitorElement * emIso_PV;
	
	// histograms of reconstructed hit Et and supercluster Et
	MonitorElement * emIso; 
	MonitorElement * recHitEtEB;
	MonitorElement * recHitEtEE;
	
	MonitorElement * scHitEtEB;
	MonitorElement * scHitEtEE;
	
	MonitorElement * scHitE_EB;
	MonitorElement * scHitE_EE;
	
	//Eta
	MonitorElement * scEta_EB;
	MonitorElement * scEta_EE;
	
	//Phi
	MonitorElement * scPhi_EB;
	MonitorElement * scPhi_EE;
	
	//sc sigma eta_eta and eta phi
	
	MonitorElement *scSigmaIetaIeta_EB;
	MonitorElement *scSigmaIetaIeta_EE;

	MonitorElement *scSigmaIetaIphi_EB;
	MonitorElement *scSigmaIetaIphi_EE;
	
	//R9
	MonitorElement *r9_EB; 
	MonitorElement *r9_EE;
	
	bool initCaloGeometry_;
	edm::ESHandle<CaloGeometry> geomH;
	edm::ESHandle<CaloTopology> caloTop;


};

#endif
