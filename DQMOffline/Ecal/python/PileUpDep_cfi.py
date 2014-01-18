import FWCore.ParameterSet.Config as cms
PileUpDep = cms.EDAnalyzer('PileUpDep',
			   VertexCollection = cms.InputTag("offlinePrimaryVertices"),
                           basicClusterCollection_EE = cms.InputTag("multi5x5SuperClusters","multi5x5EndcapBasicClusters"),
			   basicClusterCollection_EB = cms.InputTag("hybridSuperClusters","hybridBarrelBasicClusters"),
                           superClusterCollection_EE = cms.InputTag("correctedMulti5x5SuperClustersWithPreshower"),
			   superClusterCollection_EB = cms.InputTag("correctedHybridSuperClusters"),
			   RecHitCollection_EB = cms.InputTag("ecalRecHit","EcalRecHitsEB"),
			   RecHitCollection_EE = cms.InputTag("ecalRecHit","EcalRecHitsEE"),

                      )

