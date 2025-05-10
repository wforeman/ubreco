#ifndef SELECTIONTOOLBASE_H
#define SELECTIONTOOLBASE_H
////////////////////////////////////////////////////////////////////////
//
// Class:       IHitEfficiencyHistogramTool
// Module Type: tool
// File:        IHitEfficiencyHistogramTool.h
//
//              This provides an interface for tools which do histogramming
//              of various quantities associated to recob::Hit objects
//
// Created by David Caratelli (davidc@fnal.gov) on January 30 2019
//
////////////////////////////////////////////////////////////////////////

// art TOOLS
#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"

#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Principal/Event.h"

#include "ubreco/BlipReco/Utils/NuSelectionTypedefs.h"

#include "TTree.h"
#include <limits>

namespace selection {

  using ProxyPfpColl_t = nuselection::ProxyPfpColl_t;
  using ProxyPfpElem_t = nuselection::ProxyPfpElem_t;
  using ProxyClusColl_t = nuselection::ProxyClusColl_t;
  using ProxyClusElem_t = nuselection::ProxyClusElem_t;
  
  // a map linking the PFP Self() attribute used for hierarchy building to the PFP index in the event record
  std::map<unsigned int, unsigned int> _pfpmap;
  
  void BuildPFPMap(const ProxyPfpColl_t &pfp_pxy_col);
  void AddDaughters(const ProxyPfpElem_t &pfp_pxy,
                    const ProxyPfpColl_t &pfp_pxy_col,
                    std::vector<ProxyPfpElem_t> &slice_v);


class SelectionToolBase {

public:

    /**
     *  @brief  Virtual Destructor
     */
    virtual ~SelectionToolBase() noexcept = default;
    
    /**
     *  @brief Interface for configuring the particular algorithm tool
     *
     *  @param ParameterSet  The input set of parameters for configuration
     */
    void configure(const fhicl::ParameterSet&){};

    /**
     * @brief Selection function
     *
     * @param art::Event event record for selection
     */
    virtual bool selectEvent(art::Event const& e,
			     const std::vector<ProxyPfpElem_t>& pfp_pxy_v) = 0;

    /**
     * @brief set branches for TTree
     */
    virtual void setBranches(TTree* _tree) = 0;

    
    /**
     * @brief resetset TTree branches
     */
    virtual void resetTTree(TTree* _tree) = 0;


    /**
     * @brief set if data
     */
    void SetData(bool isdata) { fData = isdata; }


 protected:


    bool fData;


};



void BuildPFPMap(const ProxyPfpColl_t &pfp_pxy_col)
{
  _pfpmap.clear();
  unsigned int p = 0;
  for (const auto &pfp_pxy : pfp_pxy_col)
  { _pfpmap[pfp_pxy->Self()] = p; p++; }
  return;
} // BuildPFPMap

void AddDaughters(const ProxyPfpElem_t &pfp_pxy,
                  const ProxyPfpColl_t &pfp_pxy_col,
                  std::vector<ProxyPfpElem_t> &slice_v)
{
  auto daughters = pfp_pxy->Daughters();
  slice_v.push_back(pfp_pxy);
  //std::cout << "\t PFP w/ PdgCode " << pfp_pxy->PdgCode() << " has " << daughters.size() << " daughters" << std::endl;
  for (auto const &daughterid : daughters){
    if (_pfpmap.find(daughterid) == _pfpmap.end())
    {
      //std::cout << "Did not find DAUGHTERID in map! error" << std::endl;
      continue;
    }

    // const art::Ptr<recob::PFParticle> pfp_pxy(pfp_pxy_col, _pfpmap.at(daughterid) );
    auto pfp_pxy2 = pfp_pxy_col.begin();
    for (size_t j = 0; j < _pfpmap.at(daughterid); ++j)
      ++pfp_pxy2;
    // const T& pfp_pxy2 = (pfp_pxy_col.begin()+_pfpmap.at(daughterid));
    AddDaughters(*pfp_pxy2, pfp_pxy_col, slice_v);
  } // for all daughters

  return;
} // AddDaughters

} // selection namespace

#endif
