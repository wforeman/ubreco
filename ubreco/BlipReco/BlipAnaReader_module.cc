//#####################################################################
//###  BlipAnaReader analyzer module
//###
//###  Contains algorithms for reconstructing isolated, MeV-scale energy
//###  depositions in the TPC, called "blips." A TTree is made for offline
//###  analysis and plot-making. Algs will eventually be migrated into 
//###  dedicated alg/tool classes as appropriate.
//###
//###  Author: Will Foreman (wforeman_at_iit.edu)
//###  Date:   Sept 2021
//#####################################################################
#ifndef BLIPANA_H
#define BLIPANA_H

// Framework includes
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// LArSoft includes
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "cetlib/search_path.h"

// MicroBooNE-specific includes
#include "ubevt/Database/UbooneElectronLifetimeProvider.h"
#include "ubevt/Database/UbooneElectronLifetimeService.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"
#include "ubreco/BlipReco/Alg/BlipRecoAlg.h"
#include "ubreco/BlipReco/Utils/NuSelectionToolBase.h"
#include "ubreco/BlipReco/Utils/NuSelectionSCECorrections.h"
#include "ubreco/BlipReco/Utils/NuSelectionTrackShowerScoreFuncs.h"

// C++ includes
#include <cstring>
#include <utility>
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <typeinfo>
#include <cmath>

// ROOT includes
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"

// Helper templates for initializing arrays
namespace{  
  template <typename ITER, typename TYPE> 
    inline void FillWith(ITER from, ITER to, TYPE value) 
    { std::fill(from, to, value); }
  template <typename ITER, typename TYPE> 
    inline void FillWith(ITER from, size_t n, TYPE value)
    { std::fill(from, from + n, value); }
  template <typename CONT, typename V>
    inline void FillWith(CONT& data, const V& value)
    { FillWith(std::begin(data), std::end(data), value); }
}


// Set global constants and max array sizes
//const int kMaxHits    =  30000;
//const int kMaxClusts  =  10000; 
//const int kMaxTrks    =  1000;
const int kMaxBlips   = 10000;
//const int kMaxG4      = 100000;
//const int kMaxEDeps   = 10000;
//const int kMaxTrkPts  =   2000;  

class BlipAnaReader;
  
//###################################################
//  Data storage structure
//###################################################
class BlipAnaReaderTreeDataStruct 
{
  public:

  // --- TTrees
  TTree* evtTree;

  // --- Configurations and switches ---
  std::string treeName      = "anatree";

  // --- Event information ---   
  int           event;                // event number
  int           run;                  // run number
  int           subrun;               // subrun number
  unsigned int  timestamp;            // unix time of event
  float         lifetime;             // electron lifetime
  int           badchans;             // #bad chans according to wirecell
  int           longtrks;             // tracks > 5 cm
  
  // --- MCTruth neutrino info ---
  int   mctruth_nu_pdg;         // Neutrino PDG (if present, otherwise = 0)
  int   mctruth_nu_ccnc;        // CC (0) or NC (1)
  int   mctruth_nu_mode;        // interaction mode from Genie
  float mctruth_nu_vtx_x;       // vertex X 
  float mctruth_nu_vtx_y;       // vertex Y
  float mctruth_nu_vtx_z;       // vertex Z
  float mctruth_nu_KE;          // kinetic energy 
  
  // --- 3D Blip information ---
  int   nblips;                       // number of blips in event
  int   blip_id[kMaxBlips];           // blip ID / index
  int   blip_tpc[kMaxBlips];          // blip TPC
  int   blip_nplanes[kMaxBlips];      // number of planes matched (2 or 3)
  float blip_x[kMaxBlips];            // X position [cm]
  float blip_y[kMaxBlips];            // Y position [cm]
  float blip_z[kMaxBlips];            // Z position [cm]
  float blip_xSCE[kMaxBlips];         // X position SCE-corrected
  float blip_ySCE[kMaxBlips];         // Y position SCE-corrected
  float blip_zSCE[kMaxBlips];         // Z position SCE-corrected
  float blip_sigmayz[kMaxBlips];      // difference in wire intersection points
  float blip_dx[kMaxBlips];           // length along drift direction [cm]
  float blip_dw[kMaxBlips];           // length projected onto axis perpendicular to wire orientation
  float blip_size[kMaxBlips];         // rough size estimation based on time-tick extent and wire span
  int   blip_charge[kMaxBlips];       // blip charge at anode [e-]
  int   blip_chargecorr[kMaxBlips];       // blip charge at anode [e-]
  float blip_energy[kMaxBlips];       // blip reco energy [MeVee]
  float blip_energycorr[kMaxBlips];   // blip reco energy w/SCE+lifetime corrections [MeVee]
  float blip_energyTrue[kMaxBlips];   // blip truth energy [MeV]
  float blip_proxtrkdist[kMaxBlips];  // distance to nearest track
  int   blip_proxtrkid[kMaxBlips];    // index of nearest trk
  bool  blip_touchtrk[kMaxBlips];     // is blip touching track?
  int   blip_touchtrkid[kMaxBlips];   // track ID of touched track
  int   blip_clustid[kNplanes][kMaxBlips];     // cluster ID per plane
  
  // --- Reconstructed neutrino slice information (Pandora) --
  bool      nu_isNeutrino;            // neutrino slice identified by Pandora
  vfloat_t  nu_nuscore;               // neutrino score
  int       nu_pfp_pdg;               // PDG particle best matching with reco slice
  float     nu_reco_vtx_x;            // reconstructed vertex X [cm]
  float     nu_reco_vtx_y;            // reconstructed vertex Y [cm]
  float     nu_reco_vtx_z;            // reconstructed vertex Z [cm]
  vint_t    nu_trk_id;        // trackIDs for tracks in this PFP
  vfloat_t  nu_trk_score;     // track scores for tracks in this PFP
  vfloat_t  nu_shwr_score;    // shower scores in this PFP
  
  // === Function for resetting data ===
  void Clear(){ 
    event                 = -999; // --- event-wide info ---
    run                   = -999;
    subrun                = -999; 
    lifetime              = -999;
    badchans              = -99;
    longtrks              = -99;
    timestamp             = -999;
    
    mctruth_nu_pdg        = 0;
    mctruth_nu_ccnc       = -9;
    mctruth_nu_mode       = -9;
    mctruth_nu_vtx_x      = -999;
    mctruth_nu_vtx_y      = -999;
    mctruth_nu_vtx_z      = -999;
    mctruth_nu_KE         = -999;
    nu_isNeutrino         = false;
    nu_nuscore            .clear();
    nu_pfp_pdg            = -999;
    nu_reco_vtx_x         = -999;
    nu_reco_vtx_y         = -999;
    nu_reco_vtx_z         = -999;
    nu_trk_score          .clear();
    nu_trk_id             .clear();
    nu_shwr_score         .clear();  
   
    nblips                    = 0;
    FillWith(blip_id,         -9);
    FillWith(blip_tpc,        -9);
    FillWith(blip_nplanes,    -9);
    FillWith(blip_x,          -9999);
    FillWith(blip_y,          -9999);
    FillWith(blip_z,          -9999);
    FillWith(blip_sigmayz,    -9);
    FillWith(blip_dx,         -9);
    FillWith(blip_dw,        -9);
    FillWith(blip_size,       -9);
    FillWith(blip_charge,     -999);
    FillWith(blip_chargecorr,     -999);
    FillWith(blip_energy,     -999);
    FillWith(blip_energycorr,     -999);
    FillWith(blip_energyTrue, -999);
    FillWith(blip_proxtrkdist,-99);
    FillWith(blip_proxtrkid,  -9);
    FillWith(blip_touchtrk,   false);
    FillWith(blip_touchtrkid,  -9);
  }

  // === Function for initializing tree branches ===
  void MakeTree(){
    
    art::ServiceHandle<art::TFileService> tfs;
    auto vf = "std::vector<float>";
    auto vi = "std::vector<int>";
   
    evtTree = tfs->make<TTree>(treeName.c_str(),"analysis tree");
    evtTree->Branch("event",&event,"event/I");
    evtTree->Branch("run",&run,"run/I");
    evtTree->Branch("subrun",&subrun,"subrun/I");
    evtTree->Branch("timestamp",&timestamp,"timestamp/i");
    //evtTree->Branch("timestamp_hr",&timestamp_hr,"timestamp_hr/F");
    evtTree->Branch("lifetime",&lifetime,"lifetime/F");
    evtTree->Branch("badchans",&badchans,"badchans/I");
    evtTree->Branch("longtrks",&longtrks,"longtrks/I");
    evtTree->Branch("mctruth_nu_pdg",&mctruth_nu_pdg,"mctruth_nu_pdg/I");
    evtTree->Branch("mctruth_nu_ccnc",&mctruth_nu_ccnc,"mctruth_nu_ccnc/I");
    evtTree->Branch("mctruth_nu_mode",&mctruth_nu_mode,"mctruth_nu_mode/I");
    evtTree->Branch("mctruth_nu_vtx_x",&mctruth_nu_vtx_x,"mctruth_nu_vtx_x/F");
    evtTree->Branch("mctruth_nu_vtx_y",&mctruth_nu_vtx_y,"mctruth_nu_vtx_y/F");
    evtTree->Branch("mctruth_nu_vtx_z",&mctruth_nu_vtx_z,"mctruth_nu_vtx_z/F");
    evtTree->Branch("mctruth_nu_KE",&mctruth_nu_KE,"mctruth_nu_KE/F");
    evtTree->Branch("nu_isNeutrino",&nu_isNeutrino,"nu_isNeutrino/O");
    //evtTree->Branch("nu_nuscore",&nu_nuscore,"nu_nuscore/F");
    evtTree->Branch("nu_nuscore", vf, &nu_nuscore);
    evtTree->Branch("nu_pfp_pdg",&nu_pfp_pdg,"nu_pfp_pdg/I");
    evtTree->Branch("nu_reco_vtx_x",&nu_reco_vtx_x,"nu_reco_vtx_x/F");
    evtTree->Branch("nu_reco_vtx_y",&nu_reco_vtx_y,"nu_reco_vtx_y/F");
    evtTree->Branch("nu_reco_vtx_z",&nu_reco_vtx_z,"nu_reco_vtx_z/F");
    evtTree->Branch("nu_trk_id", vi, &nu_trk_id);
    evtTree->Branch("nu_trk_score", vf, &nu_trk_score);
    evtTree->Branch("nu_shwr_score", vf, &nu_shwr_score);
    evtTree->Branch("nblips",&nblips,"nblips/I");
    evtTree->Branch("blip_nplanes",blip_nplanes,"blip_nplanes[nblips]/I");
    evtTree->Branch("blip_x",blip_x,"blip_x[nblips]/F");
    evtTree->Branch("blip_y",blip_y,"blip_y[nblips]/F");
    evtTree->Branch("blip_z",blip_z,"blip_z[nblips]/F");
    evtTree->Branch("blip_xSCE",blip_x,"blip_xSCE[nblips]/F");
    evtTree->Branch("blip_ySCE",blip_y,"blip_ySCE[nblips]/F");
    evtTree->Branch("blip_zSCE",blip_z,"blip_zSCE[nblips]/F");
    //evtTree->Branch("blip_sigmayz",blip_sigmayz,"blip_sigmayz[nblips]/F");
    evtTree->Branch("blip_dx",blip_dx,"blip_dx[nblips]/F");
    evtTree->Branch("blip_dw",blip_dw,"blip_dw[nblips]/F");
    evtTree->Branch("blip_size",blip_size,"blip_size[nblips]/F");
    evtTree->Branch("blip_charge",blip_charge,"blip_charge[nblips]/I");
    evtTree->Branch("blip_charge_corr",blip_chargecorr,"blip_charge_corr[nblips]/I");
    evtTree->Branch("blip_energy",blip_energy,"blip_energy[nblips]/F");
    evtTree->Branch("blip_energy_corr",blip_energycorr,"blip_energy_corr[nblips]/F");
    evtTree->Branch("blip_energyTrue",blip_energyTrue,"blip_energyTrue[nblips]/F");
    evtTree->Branch("blip_proxtrkdist",blip_proxtrkdist,"blip_proxtrkdist[nblips]/F");
    evtTree->Branch("blip_touchtrk",blip_touchtrk,"blip_touchtrk[nblips]/O");
    for(int i=0;i<kNplanes;i++) evtTree->Branch(Form("blip_pl%i_clustid",i),blip_clustid[i],Form("blip_pl%i_clustid[nblips]/I",i));
  }
    
};//BlipAnaReaderTreeDataStruct class



//###################################################
//  BlipAnaReader class definition
//###################################################
class BlipAnaReader : public art::EDAnalyzer 
{ 
  public:
  explicit BlipAnaReader(fhicl::ParameterSet const& pset);
  virtual ~BlipAnaReader();
  
  //void beginJob();                      // called once, at start of job
  void endJob();                        // called once, at end of job
  void analyze(const art::Event& evt);  // called per event

  private:
  void    PrintBlipInfo(const blipobj::Blip&);

  // --- Data and calo objects ---
  BlipAnaReaderTreeDataStruct*  fData;

  // --- Input products
  std::string         fBlipProducer;

  // --- Counters and such ---
  bool  fIsRealData         = false;
  bool  fIsMC               = false;
  int   fNumEvents          = 0;
  int   fNum3DBlips         = 0;
  int   fNum3DBlips3Plane   = 0;
  int   fNum3DBlipsTrue     = 0;
  
  // --- Neutrino selection tools
  using ProxyPfpColl_t = selection::ProxyPfpColl_t;
  using ProxyPfpElem_t = selection::ProxyPfpElem_t;
  art::InputTag fPFPproducer;
  art::InputTag fCLSproducer; // cluster associated to PFP
  art::InputTag fSLCproducer; // slice associated to PFP
  art::InputTag fHITproducer; // hit associated to cluster
  art::InputTag fSHRproducer; // shower associated to PFP
  art::InputTag fVTXproducer; // vertex associated to PFP
  art::InputTag fPCAproducer; // PCAxis associated to PFP
  art::InputTag fMCTproducer;
  art::InputTag fTRKproducer;

  //void    BuildPFPMap(const ProxyPfpColl_t&);
  //void    AddDaughters(const ProxyPfpElem_t &pfp_pxy,
  //                  const ProxyPfpColl_t &pfp_pxy_col,
  //                  std::vector<ProxyPfpElem_t> &slice_v);

  // a map linking the PFP Self() attribute used for hierarchy building to the PFP index in the event record
  std::map<unsigned int, unsigned int> _pfpmap;

  // --- Histograms ---
  TH1D*   h_nblips;
  TH1D*   h_nblips_3plane;
  TH2D*   h_blip_zy;
  TH2D*   h_blip_zy_3plane;
  TH1D*   h_blip_nplanes;
  TH1D*   h_blip_charge;
  TH1D*   h_blip_energy;
  
  TH1D*   h_nblips_tm;
  TH2D*   h_blip_reszy;
  TH1D*   h_blip_resx;
  TH2D*   h_blip_E_tvr;
  TH2D*   h_blip_resE;

  // Initialize histograms
  void InitializeHistograms(){
    
    art::ServiceHandle<art::TFileService> tfs;
    
    float blipMax   = 500;    
    int blipBins    = 500;
    h_nblips          = tfs->make<TH1D>("nblips","Reconstructed 3D blips per event",blipBins,0,blipMax);
    h_nblips_3plane   = tfs->make<TH1D>("nblips_3plane","Reconstructed 3D blips per event (3-plane match, intersect #Delta < 1 cm)",blipBins,0,blipMax);
    h_blip_zy         = tfs->make<TH2D>("blip_zy","3D blip location;Z [cm];Y [cm]",600,-100,1100,150,-150,150);
    h_blip_zy_3plane  = tfs->make<TH2D>("blip_zy_3plane","3D blip location (3-plane match, intersect #Delta < 1 cm);Z [cm];Y [cm]",600,-100,1100,150,-150,150);
    h_blip_nplanes    = tfs->make<TH1D>("blip_nplanes","Matched planes per blip",3,1,4);
    h_blip_charge     = tfs->make<TH1D>("blip_charge","3D blips;Reconstructed charge [e-]",200,0,100e3);
    h_blip_energy     = tfs->make<TH1D>("blip_energy","3D blips;Reconstructed energy [MeV]",200,0,10);
      h_blip_zy       ->SetOption("COLZ");
      h_blip_zy_3plane->SetOption("COLZ");
    
    h_nblips_tm       = tfs->make<TH1D>("nblips_tm","Truth-matched 3D blips per event",blipBins,0,blipMax);
    h_blip_reszy      = tfs->make<TH2D>("blip_res_zy","Blip position resolution;Z_{reco} - Z_{true} [cm];Y_{reco} - Y_{true} [cm]",150,-15,15,150,-15,15);
    h_blip_resx       = tfs->make<TH1D>("blip_res_x","Blip position resolution;X_{reco} - X_{true} [cm]",150,-15,15);
    h_blip_E_tvr      = tfs->make<TH2D>("blip_E_true_vs_reco","Truth-matched 3D blips;True energy deposited [MeV];Reconstructed energy [MeVee]",200,0,10,200,0,10);
    h_blip_resE       = tfs->make<TH2D>("blip_res_energy","Energy resolution of 3D blips;Energy [MeV];#deltaE/E_{true}",200,0,5,200,-1.5,2.5);
      h_blip_E_tvr    ->SetOption("colz");
      h_blip_reszy    ->SetOption("colz");
      h_blip_resE     ->SetOption("colz");
    
 }

};//class BlipAnaReader


//###################################################
//  BlipAnaReader constructor and destructor
//###################################################
BlipAnaReader::BlipAnaReader(fhicl::ParameterSet const& pset) : 
  EDAnalyzer(pset)
  ,fData  (nullptr)
{
  // blip reconstruction algorithm class
  fBlipProducer   = pset.get<std::string>   ("BlipProducer","blipreco");

  // data tree object
  fData = new BlipAnaReaderTreeDataStruct();
  fData ->Clear();
  fData ->MakeTree();
  
  // initialize histograms
  InitializeHistograms();
  
  fPFPproducer = pset.get<art::InputTag>("PFPproducer","pandora");
  fSHRproducer = pset.get<art::InputTag>("SHRproducer","shrreco3d");
  fHITproducer = pset.get<art::InputTag>("HITproducer","pandora");
  fVTXproducer = pset.get<art::InputTag>("VTXproducer","pandora");
  fPCAproducer = pset.get<art::InputTag>("PCAproducer","pandora");
  fCLSproducer = pset.get<art::InputTag>("CLSproducer","pandora");
  fSLCproducer = pset.get<art::InputTag>("SLCproducer","pandora");
  fMCTproducer = pset.get<art::InputTag>("MCTproducer","generator");
  fTRKproducer = pset.get<art::InputTag>("TRKproducer","pandora");

    
}
BlipAnaReader::~BlipAnaReader(){}



//###################################################
//  Main event-by-event analysis
//###################################################
void BlipAnaReader::analyze(const art::Event& evt)
{ 
  //============================================
  // New event!
  //============================================
  fData            ->Clear();
  fData->event      = evt.id().event();
  fData->run        = evt.id().run();
  fData->subrun     = evt.id().subRun();
  fIsRealData       = evt.isRealData();
  fNumEvents++;

  // Get timestamp
  unsigned long long int tsval = evt.time().value();
  const unsigned long int mask32 = 0xFFFFFFFFUL;
  fData->timestamp = ( tsval >> 32 ) & mask32;

  // Retrieve lifetime
  const lariov::UBElectronLifetimeProvider& elifetime_provider 
    = art::ServiceHandle<lariov::UBElectronLifetimeService>()->GetProvider();
  float electronLifetime = elifetime_provider.Lifetime() * 1e3; // convert ms->mus
  fData->lifetime = electronLifetime;
  
  // Tell us what's going on!
  if( fNumEvents < 200 || (fNumEvents % 100) == 0 ) { 
  std::cout<<"\n"
  <<"=========== BlipAnaReader =========================\n"
  <<"Event "<<evt.id().event()<<" / run "<<evt.id().run()<<"; total: "<<fNumEvents<<"\n";
  }
  
  //===================================
  // Check neutrinos in MCTruth
  // (NuanceOffset found through simb::kNuanceOffset)
  //===================================
  art::Handle< std::vector<simb::MCTruth> > truthHandle;
  std::vector<art::Ptr<simb::MCTruth> > truthlist;
  if (evt.getByLabel("generator",truthHandle))
    art::fill_ptr_vector(truthlist, truthHandle); 
  if( truthlist.size() > 0 ) {
    auto& nu    = truthlist[0]->GetNeutrino();
    auto& part  = truthlist[0]->GetParticle(0);
    auto PDG    = fabs(part.PdgCode());
    if ( (PDG == 12) || (PDG == 14) ) {
      fData->mctruth_nu_pdg   = part.PdgCode();
      fData->mctruth_nu_ccnc  = nu.CCNC();
      fData->mctruth_nu_mode  = nu.Mode();
      fData->mctruth_nu_vtx_x = part.EndPosition()[0];
      fData->mctruth_nu_vtx_y = part.EndPosition()[1];
      fData->mctruth_nu_vtx_z = part.EndPosition()[2];
      fData->mctruth_nu_KE    =  /*GeV->MeV*/1e3 * (part.E()-part.Mass());
    }
  }
  
  
  //=======================================================
  // Check if the event contains neutrino information
  // (much of this copied from ubana/searchingfornues)
  //=======================================================
  // -- PFPs
  art::Handle< std::vector<recob::PFParticle> > pfpHandle;
  std::vector<art::Ptr<recob::PFParticle> > pfplist;
  if (evt.getByLabel(fPFPproducer,pfpHandle))
    art::fill_ptr_vector(pfplist, pfpHandle);

  // -- associated tracks/vertex
  if( pfplist.size() ) {
    // grab PFParticles in event
    ProxyPfpColl_t const &pfp_proxy = proxy::getCollection<std::vector<recob::PFParticle>>(evt, fPFPproducer,
      proxy::withAssociated<larpandoraobj::PFParticleMetadata>(fPFPproducer),
      proxy::withAssociated<recob::Cluster>(fCLSproducer),
      proxy::withAssociated<recob::Slice>(fSLCproducer),
      proxy::withAssociated<recob::Track>(fTRKproducer),
      proxy::withAssociated<recob::Vertex>(fVTXproducer),
      proxy::withAssociated<recob::PCAxis>(fPCAproducer),
      proxy::withAssociated<recob::Shower>(fSHRproducer),
      proxy::withAssociated<recob::SpacePoint>(fPFPproducer));
    selection::BuildPFPMap(pfp_proxy);
    // loopthrough PFParticles
    for (const ProxyPfpElem_t &pfp_pxy : pfp_proxy)
    {
      // get metadata for this PFP
      const auto &pfParticleMetadataList = pfp_pxy.get<larpandoraobj::PFParticleMetadata>();

      //  find neutrino candidate
      if (pfp_pxy->IsPrimary() == false) continue;
      auto PDG = fabs(pfp_pxy->PdgCode());
      fData->nu_pfp_pdg=PDG;
      if ( (PDG == 12) || (PDG == 14) )
      {
        //std::cout<<"Found a neutrino PFP\n";
        if (pfParticleMetadataList.size() != 0)
        {
          for (unsigned int j = 0; j < pfParticleMetadataList.size(); ++j)
          {
            const art::Ptr<larpandoraobj::PFParticleMetadata> &pfParticleMetadata(pfParticleMetadataList.at(j));
            auto pfParticlePropertiesMap = pfParticleMetadata->GetPropertiesMap();
            if (!pfParticlePropertiesMap.empty())
            {
              for (std::map<std::string, float>::const_iterator it = pfParticlePropertiesMap.begin(); it != pfParticlePropertiesMap.end(); ++it)
              {
                if( it->first == "IsNeutrino" ) fData->nu_isNeutrino  = it->second;
                if( it->first == "NuScore"    ) fData->nu_nuscore.push_back(it->second);
              }
            }
          }
        } // if PFP metadata exists!
  

        // Get vertex info
        double xyz[3] = {};
        auto vtx = pfp_pxy.get<recob::Vertex>();
        if (vtx.size() == 1)
        {
          // save vertex to array
          vtx.at(0)->XYZ(xyz);
          auto nuvtx = TVector3(xyz[0], xyz[1], xyz[2]);
          float _reco_nu_vtx_sce[3];
          nuselection::ApplySCECorrectionXYZ(nuvtx.X(),nuvtx.Y(),nuvtx.Z(), _reco_nu_vtx_sce);
          fData->nu_reco_vtx_x = _reco_nu_vtx_sce[0];
          fData->nu_reco_vtx_y = _reco_nu_vtx_sce[1];
          fData->nu_reco_vtx_z = _reco_nu_vtx_sce[2];
          //std::cout<<"Vertex:  "<<nuvtx.X()<<"  "<<nuvtx.Y()<<"  "<<nuvtx.Z()<<"\n";
        }
        else
        {
          std::cout << "ERROR. Found neutrino PFP w/ != 1 associated vertices..." << std::endl;
        }

        // collect PFParticle hierarchy originating from this neutrino candidate
        std::vector<ProxyPfpElem_t> slice_pfp_v;
        selection::AddDaughters(pfp_pxy, pfp_proxy, slice_pfp_v);
        //std::cout << "This slice has " << slice_pfp_v.size() << " daughter PFParticles" << std::endl;
        // create list of tracks and showers associated to this slice
        std::vector<float>  trkscore_v;
        std::vector<int>    trkid_v;
        std::vector<float>  shwrscore_v;

        for (auto pfp : slice_pfp_v)
        { 
          auto const &ass_trk_v = pfp.get<recob::Track>();
          auto const &ass_shr_v = pfp.get<recob::Shower>();
          float score = nuselection::GetTrackShowerScore(pfp);
          if (ass_trk_v.size() == 1) {
            trkid_v   .push_back(ass_trk_v.at(0)->ID());
            trkscore_v.push_back(score);
          }
          if (ass_shr_v.size() == 1) {
            shwrscore_v.push_back(score);
          }
        } // for all PFParticles in the slice
        //std::cout<<"  - "<<trkscore_v.size()<<" tracks\n";
        //std::cout<<"  - "<<shwrscore_v.size()<<" showers\n";
        for(size_t itrk = 0; itrk < trkscore_v.size(); itrk++){
          fData->nu_trk_id    .push_back(trkid_v[itrk]);
          fData->nu_trk_score .push_back(trkscore_v[itrk]);
        }
        for(size_t ishwr = 0; ishwr < shwrscore_v.size(); ishwr++){
          fData->nu_shwr_score .push_back(shwrscore_v[ishwr]);
        }
      
      }//if PDG of neutrino
    }//end loop over PFPs
  }
  
  
  //===========================================
  // Check if blip objects were saved to the event;
  //===========================================
  art::Handle< std::vector<blipobj::Blip> > blipHandle;
  std::vector<art::Ptr<blipobj::Blip> > bliplist;
  if (evt.getByLabel("blipreco",blipHandle))
    art::fill_ptr_vector(bliplist, blipHandle);
  if( bliplist.size() ) {
    
    int nblips_3p  = 0;
    int nblips_tm  = 0;
    fData->nblips = bliplist.size();
    for(size_t i=0; i<bliplist.size(); i++){
      auto& blip = bliplist[i];

      fNum3DBlips++;
      int nplanes = blip->NPlanes;
      float x = blip->Position.X();
      float y = blip->Position.Y();
      float z = blip->Position.Z();
      float E = blip->Energy;
      
      // Fill variables to be written to TTree
      fData->blip_id[i]       = blip->ID;
      fData->blip_tpc[i]      = blip->TPC;
      fData->blip_nplanes[i]  = nplanes;
      fData->blip_x[i]        = x;
      fData->blip_y[i]        = y;
      fData->blip_z[i]        = z;
      fData->blip_xSCE[i]     = blip->PositionSCE.X();
      fData->blip_ySCE[i]     = blip->PositionSCE.Y();
      fData->blip_zSCE[i]     = blip->PositionSCE.Z();
      fData->blip_sigmayz[i]  = blip->SigmaYZ;
      fData->blip_dw[i]       = blip->dYZ;
      fData->blip_charge[i]   = blip->Charge;
      fData->blip_energy[i]   = blip->Energy;
      fData->blip_chargecorr[i] = blip->ChargeCorr;
      fData->blip_energycorr[i] = blip->EnergyCorr;
      fData->blip_energyTrue[i] = blip->truth.Energy;
      fData->blip_proxtrkdist[i] = blip->ProxTrkDist;  
      fData->blip_proxtrkid[i] = blip->ProxTrkID;
      fData->blip_touchtrkid[i] = blip->TouchTrkID;
      for(size_t ipl = 0; ipl<kNplanes; ipl++){
        if( blip->clusters[ipl].NHits <= 0 ) continue;
        fData->blip_clustid[ipl][i] = blip->clusters[ipl].ID;
      }

      // Fill histograms
      h_blip_nplanes->Fill(nplanes);
      h_blip_zy->Fill(z,y);
      h_blip_charge->Fill(blip->Charge);
      h_blip_energy->Fill(blip->Energy);
      if( nplanes == 3 ) {
        fNum3DBlips3Plane++;
        nblips_3p++;
        h_blip_zy_3plane->Fill(z,y);
      }
      if( blip->truth.ID >= 0 ) {
        fNum3DBlipsTrue++;
        nblips_tm++;
        float true_x = blip->truth.Position.X();
        float true_y = blip->truth.Position.Y();
        float true_z = blip->truth.Position.Z();
        float true_E = blip->truth.Energy;
        h_blip_reszy->Fill( z - true_z, y - true_y);
        h_blip_resx->Fill( x - true_x );
        h_blip_E_tvr->Fill( true_E, E );
        h_blip_resE->Fill( true_E, (E-true_E)/true_E );
      }
    
    }// Done looping over blips loaded from event
    h_nblips        ->Fill(bliplist.size());
    h_nblips_3plane ->Fill(nblips_3p);
    h_nblips_tm     ->Fill(nblips_tm);

  }//endif event contains blips
    
    
  //====================================
  // Fill TTree
  //====================================
  fData->evtTree->Fill();

}


//###################################################
//  endJob: output useful info to screen
//###################################################
void BlipAnaReader::endJob(){
  printf("\n***********************************************\n");
  printf("BlipAnaReader Summary\n\n");
  printf("  Total events                : %i\n",        fNumEvents);
  printf("  Blips per evt, total        : %.3f\n",      (float)fNum3DBlips/fNumEvents);
  printf("                 3 planes     : %.3f\n",      (float)fNum3DBlips3Plane/fNumEvents);
  if(fIsMC){
  printf("  MC-matched blips per evt    : %.3f\n",       (float)fNum3DBlipsTrue/fNumEvents);
  }
  printf("\n***********************************************\n");
}


void BlipAnaReader::PrintBlipInfo(const blipobj::Blip& bl){
  printf("  blipID: %4i, TPC: %i, charge: %8.0i,  recoEnergy: %8.3f MeV, XYZ: %6.1f, %6.1f, %6.1f,   EdepID: %i\n",
  bl.ID,
  bl.TPC,
  (int)bl.clusters[2].Charge,
  bl.Energy,
  bl.Position.X(),bl.Position.Y(),bl.Position.Z(),
  bl.truth.ID
  );
}



DEFINE_ART_MODULE(BlipAnaReader)

#endif
