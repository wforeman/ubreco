#ifndef CLUSTER_CLUSTERMAKER_CXX
#define CLUSTER_CLUSTERMAKER_CXX

#include "ClusterMaker.h"

#include "larcore/Geometry/WireReadout.h"

namespace cluster {

  ClusterMaker::ClusterMaker() {

    _vtx_w_cm = {0.,0.,0.};
    _vtx_t_cm = {0.,0.,0.};

    // get detector specific properties
    auto const& channelMap = art::ServiceHandle<geo::WireReadout>()->Get();
    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataForJob();
    auto const detp = art::ServiceHandle<detinfo::DetectorPropertiesService>()->DataForJob(clockData);
    _wire2cm = channelMap.Plane(geo::PlaneID{0,0,0}).WirePitch();
    _time2cm = sampling_rate(clockData) / 1000.0 * detp.DriftVelocity( detp.Efield(), detp.Temperature() );

  }

  void ClusterMaker::MakeClusters(detinfo::DetectorClocksData const& clockData,
                                  const art::ValidHandle<std::vector<recob::Cluster> >& clus_h,
				  const art::FindManyP<recob::Hit>&  clus_hit_assn_v,
				  const art::ValidHandle<std::vector<recob::Vertex> >& vtx_h,
				  std::vector<::cluster::Cluster>& cluster) {


    if (loadVertex(clockData, vtx_h) == false) {
      std::cout << "NO VERTEX" << std::endl;
      return;
    }

    cluster.clear();
    cluster.reserve(clus_hit_assn_v.size());

    // loop trhough hits in each cluster
    for (size_t c=0; c < clus_hit_assn_v.size(); c++) {
      
      const std::vector<art::Ptr<recob::Hit> > hit_v = clus_hit_assn_v.at(c);
      
      // create new cluster
      ::cluster::Cluster clus;
      // and associated hits
      std::vector<::cluster::pt> pts;
      // fill hit information
      GetClusterPts(hit_v, pts);
      // assign hits to cluster
      clus.SetHits(pts);
      // and add to list of clusters
      cluster.push_back(clus);
      
    }// for all input clusters

    return;
  }

  void ClusterMaker::MakeCluster(const std::vector<art::Ptr<recob::Hit> >& hit_v,
				 ::cluster::Cluster& cluster) {

    // and associated hits
    std::vector<::cluster::pt> pts;
    // fill hit information
    GetClusterPts(hit_v, pts);
    // assign hits to cluster
    cluster.SetHits(pts);

    return;
  }

  void ClusterMaker::GetClusterPts(const std::vector<art::Ptr<recob::Hit> >& hit_v,
				   std::vector<::cluster::pt>& pt_v) {

    pt_v.clear();
    
    for (art::Ptr<recob::Hit> hit : hit_v) {
      
      auto pl = hit->WireID().Plane;

      float  q = hit->Integral();

      float hw = hit->WireID().Wire * _wire2cm;
      float ht = hit->PeakTime()    * _time2cm;
      
      float r = sqrt( (ht - _vtx_t_cm[pl]) * (ht - _vtx_t_cm[pl]) +
		      (hw - _vtx_w_cm[pl]) * (hw - _vtx_w_cm[pl]) );

      float cos_theta = (hw - _vtx_w_cm[pl]) / r;
      float sin_theta = (ht - _vtx_t_cm[pl]) / r;
      float theta_rad = atan2(sin_theta,cos_theta);
      float theta_deg = theta_rad * 180. / 3.1415;

      theta_deg += 180.;

      ::cluster::pt thispt(r,theta_deg,hw,ht,q,pl);

      thispt._idx = hit.key();
      
      pt_v.push_back(thispt);
      
    }// for all his in cluster

    return;
  }// create polar cluster


  

  bool ClusterMaker::loadVertex(detinfo::DetectorClocksData const& clockData,
                                const art::ValidHandle<std::vector<::recob::Vertex> > vtx_h) {

    if (vtx_h->size() != 1) return false;

    auto const& vtx = (*vtx_h)[0].position();


    //std::cout << "Vtx coordinates : [" << xyz[0] << ", " << xyz[1] << ", " << xyz[2] << "]" << std::endl;

    auto const& channelMap = art::ServiceHandle<geo::WireReadout>()->Get();
    auto const triggerOffset = trigger_offset(clockData);
    for (size_t pl = 0; pl < 3; pl++) {

      _vtx_w_cm[pl] = channelMap.Plane(geo::PlaneID(0,0,pl)).WireCoordinate(vtx) * _wire2cm + 0.15;
      _vtx_t_cm[pl] = vtx.X() + (triggerOffset * _time2cm) + pl*0.3;

      std::cout << "trigger offset [cm] : " << (triggerOffset * _time2cm) << std::endl;
      std::cout << "Vtx @ pl " << pl << " [" << _vtx_w_cm[pl] << ", " << _vtx_t_cm[pl] << " ]" << std::endl;

      }

    return true;
  }
  

}

#endif
