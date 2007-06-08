#ifndef RecoMuon_MuonSeedGenerator_MuonSeedFinder_H
#define RecoMuon_MuonSeedGenerator_MuonSeedFinder_H

/** \class MuonSeedFinder
 *  
 *  Uses SteppingHelixPropagator
 *
 *  \author A. Vitelli - INFN Torino
 *  \author porting R. Bellan - INFN Torino
 *
 *  $Date: 2006/07/05 09:32:45 $
 *  $Revision: 1.7 $
 *  
 */

#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"

#include "DataFormats/TrajectorySeed/interface/TrajectorySeed.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"
#include "DataFormats/CSCRecHit/interface/CSCRecHit2D.h"

#include "RecoMuon/TransientTrackingRecHit/interface/MuonTransientTrackingRecHit.h"

#include <vector>

namespace edm {class EventSetup;}
class MagneticField;

class MuonSeedFinder {
public:
  /// Constructor
  MuonSeedFinder();

  /// Destructor
  virtual ~MuonSeedFinder(){};

  // Operations

  void add(MuonTransientTrackingRecHit::MuonRecHitPointer hit) { theRhits.push_back(hit); }
  
  std::vector<TrajectorySeed> seeds(const edm::EventSetup& eSetup) const;
  MuonTransientTrackingRecHit::ConstMuonRecHitPointer firstRecHit() const { return theRhits.front(); }
  unsigned int nrhit() const { return  theRhits.size(); }
  
private:
  //  TrackingRecHit best_cand(TrackingRecHit* rhit=0) const;
  bool createEndcapSeed(MuonTransientTrackingRecHit::MuonRecHitPointer me, 
			std::vector<TrajectorySeed>& theSeeds,
			const edm::EventSetup& eSetup) const;

  bool createEndcapSeed_OLD(MuonTransientTrackingRecHit::MuonRecHitPointer me, 
			    std::vector<TrajectorySeed>& theSeeds,
			    const edm::EventSetup& eSetup) const;
  
  float computePt(MuonTransientTrackingRecHit::ConstMuonRecHitPointer muon, const MagneticField *field) const;
  
  MuonTransientTrackingRecHit::MuonRecHitContainer theRhits;
 
  // put a parameterSet instead of
  // static SimpleConfigurable<float> theMinMomentum;
  float theMinMomentum;
};
#endif
