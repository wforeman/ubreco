///////////////////////////////////////////////////////
// BlipUtils.h
//
// Helper functions and algs. Based on RecoUtils and
// functions in AnalysisTree.
//
//////////////////////////////////////////////////////
#ifndef BLIPUTIL_H_SEEN
#define BLIPUTIL_H_SEEN

// framework
#include "canvas/Persistency/Common/Ptr.h" 
#include "canvas/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 

// LArSoft
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/WireReadout.h"

// c++
#include <vector>
#include <map>

#include "ubobj/Blip/DataTypes.h"
#include "TH1D.h"


typedef std::vector<art::Ptr<sim::SimEnergyDeposit>> SEDVec_t;

namespace BlipUtils {
 
  //###################################################
  // Functions related to blip reconstruction
  //###################################################
  //void      InitializeDetProps();
  void      FillParticleInfo(simb::MCParticle const&, blipobj::ParticleInfo&, SEDVec_t&, int plane=2);
  //void      CalcPartDrift(blipobj::ParticleInfo&, int);
  //void      CalcTotalDep(float&,int&,float&, SEDVec_t&);
  void      MakeTrueBlips(std::vector<blipobj::ParticleInfo>&, std::vector<blipobj::TrueBlip>&);
  void      GrowTrueBlip(blipobj::ParticleInfo&, blipobj::TrueBlip&);
  void      MergeTrueBlips(std::vector<blipobj::TrueBlip>&, float);
  void      GrowHitClust(blipobj::HitInfo const&, blipobj::HitClust&);
  bool      DoHitsOverlap(art::Ptr<recob::Hit> const&, art::Ptr<recob::Hit> const&);
  bool      DoHitClustsOverlap(blipobj::HitClust const&, blipobj::HitClust const&);
  bool      DoHitClustsOverlap(blipobj::HitClust const&,float,float);
  float     CalcHitClustsOverlap(blipobj::HitClust const&, blipobj::HitClust const&);
  float     CalcOverlap(float const&, float const&, float const&, float const&);
  bool      DoChannelsIntersect(int,int);
  bool      DoHitClustsMatch(blipobj::HitClust const&, blipobj::HitClust const&,float);
  blipobj::HitClust  MakeHitClust(std::vector<blipobj::HitInfo> const&);
  blipobj::Blip      MakeBlip(std::vector<blipobj::HitClust> const&);
  

  //###################################################
  // General functions 
  //###################################################
  //void    HitTruth(art::Ptr<recob::Hit> const&, int&, float&, float&, float&);
  //si_t    HitTruthIds( art::Ptr<recob::Hit> const&);
  //bool    G4IdToMCTruth( int const, art::Ptr<simb::MCTruth>&);
  double  PathLength(const simb::MCParticle&, TVector3&, TVector3&);
  double  PathLength(const simb::MCParticle&);
  bool    IsAncestorOf(int, int, bool breakAtPhots=false);
  double  DistToBoundary(const recob::Track::Point_t&);
  double  DistToLine(TVector3&, TVector3&, TVector3&);
  double  DistToLine2D(TVector2&, TVector2&, TVector2&);
  void    GetGeoBoundaries(double&,double&,double&,double&,double&,double&);
  bool    IsPointInAV(float x,float y,float z,float margin=0);
  bool    IsPointInAV(TVector3& p,float margin=0);
  bool    IsPointAtBnd(float,float,float);
  bool    IsPointAtBnd(TVector3&);
  void    NormalizeHist(TH1D*);
  float   FindMedian(std::vector<float>&);
  float   FindMean(std::vector<float>&);
  
}

#endif
