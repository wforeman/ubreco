#include "ubreco/BlipReco/Alg/BlipRecoAlg.h"
#include <fstream>

namespace blip {

  //###########################################################
  // Constructor
  //###########################################################
  BlipRecoAlg::BlipRecoAlg( fhicl::ParameterSet const& pset ) :
    fCaloAlg ( pset.get<fhicl::ParameterSet>("CaloAlg") )
  {
    this->reconfigure(pset);
   
    auto const& detProp   = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataForJob();
    auto const& clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataForJob();
    kLArDensity           = detProp.Density();
    kNominalEfield        = detProp.Efield();
    kDriftVelocity        = detProp.DriftVelocity(detProp.Efield(0),detProp.Temperature());
    kTickPeriod           = clockData.TPCClock().TickPeriod();
    kNominalRecombFactor  = ModBoxRecomb(fCalodEdx,kNominalEfield);
    kWion                 = 1000./util::kGeVToElectrons;
    
    // initialize channel list
    fBadChanMask       .resize(8256,false);
    fBadChanMaskPerEvt = fBadChanMask;
    if( fBadChanFile != "" ) {
      cet::search_path sp("FW_SEARCH_PATH");
      std::string fullname;
      sp.find_file(fBadChanFile,fullname);
      if (fullname.empty()) {
        throw cet::exception("Bad channel list not found");
      } else {
        std::ifstream inFile(fullname, std::ios::in);
        std::string line;
        while (std::getline(inFile,line)) {
          if( line.find("#") != std::string::npos ) continue;
          std::istringstream ss(line);
          int ch1, ch2;
          ss >> ch1;
          if( !(ss >> ch2) ) ch2 = ch1;
          for(int i=ch1; i<=ch2; i++) fBadChanMask[i] = true;
        }
      }
    }
    int NBadChansFromFile     = std::count(fBadChanMask.begin(),fBadChanMask.end(),true);
    
    EvtBadChanCount = 0;

    printf("******************************************\n");
    printf("Initializing BlipRecoAlg...\n");
    printf("  - Efield: %.4f kV/cm\n",detProp.Efield());
    printf("  - Drift velocity: %.4f\n",detProp.DriftVelocity(detProp.Efield(),detProp.Temperature()));
    printf("  - using dE/dx: %.2f MeV/cm\n",fCalodEdx);
    printf("  - equiv. recomb: %.4f\n",kNominalRecombFactor);
    printf("  - custom bad chans: %i\n",NBadChansFromFile);
    printf("*******************************************\n");

    // create diagnostic histograms
    art::ServiceHandle<art::TFileService> tfs;
    art::TFileDirectory hdir = tfs->mkdir("BlipRecoAlg");
   
    /*
    h_chanstatus     = hdir.make<TH1D>("chanstatus","Channel status for 'channels' list",5,0,5);
    h_chanstatus     ->GetXaxis()->SetBinLabel(1, "disconnected");
    h_chanstatus     ->GetXaxis()->SetBinLabel(2, "dead");
    h_chanstatus     ->GetXaxis()->SetBinLabel(3, "lownoise");
    h_chanstatus     ->GetXaxis()->SetBinLabel(4, "noisy");
    h_chanstatus     ->GetXaxis()->SetBinLabel(5, "good");
    
    h_hit_chanstatus     = hdir.make<TH1D>("hit_chanstatus","Channel status of hits",5,0,5);
    h_hit_chanstatus     ->GetXaxis()->SetBinLabel(1, "disconnected");
    h_hit_chanstatus     ->GetXaxis()->SetBinLabel(2, "dead");
    h_hit_chanstatus     ->GetXaxis()->SetBinLabel(3, "lownoise");
    h_hit_chanstatus     ->GetXaxis()->SetBinLabel(4, "noisy");
    h_hit_chanstatus     ->GetXaxis()->SetBinLabel(5, "good");
    */

    h_trkhits_mcfrac  = hdir.make<TH1D>("trkhits_mcfrac","Fraction of track hits matched to MC",101,0,1.01);
    h_chan_nhits      = hdir.make<TH1D>("chan_nhits","Untracked hits;TPC readout channel;Total hits",8256,0,8256);
    h_chan_nclusts    = hdir.make<TH1D>("chan_nclusts","Untracked isolated hits;TPC readout channel;Total clusts",8256,0,8256);
    h_chan_bad        = hdir.make<TH1D>("chan_bad","Channels marked as bad;TPC readout channel",8256,0,8256);
    h_recomb          = hdir.make<TH1D>("recomb","Applied recomb factor (nominal field)",150,0.40,0.70);
    h_recombSCE       = hdir.make<TH1D>("recombSCE","Applied recomb factor (SCE-corrected local field)",150,0.40,0.70);
    h_clust_nwires    = hdir.make<TH1D>("clust_nwires","Clusters (pre-cut);Wires in cluster",100,0,100);
    h_clust_timespan  = hdir.make<TH1D>("clust_timespan","Clusters (pre-cut);Time span [ticks]",300,0,300);

    int qbins = 200;
    float qmax = 100;
    //int wiresPerPlane[3]={2400,2400,3456};
    for(int i=0; i<kNplanes; i++) {
      //h_hit_maskfrac[i]       = dir_diag.make<TH1D>(Form("pl%i_hit_maskfrac",i),"",100,0,1.);
      //h_hit_maskfrac_true[i]  = dir_diag.make<TH1D>(Form("pl%i_hit_maskfrac_true",i),"",100,0,1.);
      //h_hit_mult[i]         = hdir.make<TH1D>(Form("pl%i_hit_mult",i),      Form("Plane %i;Num same-wire hits within +/- 50 ticks",i),20,0,20);
      if( i == fCaloPlane ) continue;
      //h_wire_nhits[i]       = hdir.make<TH1D>(Form("pl%i_wire_nhits",i),      Form("Plane %i untracked hits not plane-matched;TPC readout channel;Total hits",i),wiresPerPlane[i],0,wiresPerPlane[i]);
      h_clust_overlap[i]    = hdir.make<TH1D>(Form("pl%i_clust_overlap",i),   Form("Plane %i clusters;Overlap fraction",i),101,0,1.01);
      h_clust_dt[i]         = hdir.make<TH1D>(Form("pl%i_clust_dt",i),        Form("Plane %i clusters;dT [ticks]",i),200,-10,10);
      h_clust_dtfrac[i]     = hdir.make<TH1D>(Form("pl%i_clust_dtfrac",i),    Form("Plane %i clusters;Charge-weighted mean dT/RMS",i),120,-3,3);
      h_clust_q[i]     = hdir.make<TH2D>(Form("pl%i_clust_charge",i),  
        Form("Pre-cut;Plane %i cluster charge [#times10^{3} e-];Plane %i cluster charge [#times10^{3} e-]",fCaloPlane,i),
        qbins,0,qmax,qbins,0,qmax);
      h_clust_q[i]    ->SetOption("colz");
      h_clust_q_cut[i]     = hdir.make<TH2D>(Form("pl%i_clust_charge_cut",i),  
        Form("Post-cut;Plane %i cluster charge [#times10^{3} e-];Plane %i cluster charge [#times10^{3}]",fCaloPlane,i),
        qbins,0,qmax,qbins,0,qmax);
      h_clust_q_cut[i]    ->SetOption("colz");
      h_clust_score[i]    = hdir.make<TH1D>(Form("pl%i_clust_score",i),   Form("Plane %i clusters;Match score",i),101,0,1.01);
      h_clust_picky_overlap[i]   = hdir.make<TH1D>(Form("pl%i_clust_picky_overlap",i),  Form("Plane %i clusters (3 planes, intersect #Delta cut);Overlap fraction",i),101,0,1.01);
      h_clust_picky_dt[i]        = hdir.make<TH1D>(Form("pl%i_clust_picky_dt",i),       Form("Plane %i clusters (3 planes, intersect #Delta cut);dT [ticks]",i),200,-10,10);
      h_clust_picky_dtfrac[i]      = hdir.make<TH1D>(Form("pl%i_clust_picky_dtfrac",i),Form("Plane %i clusters (3 planes, intersect #Delta cut);Charge-weighted mean dT/RMS",i),120,-3,3);
      h_clust_picky_q[i]  = hdir.make<TH2D>(Form("pl%i_clust_picky_charge",i),  
        Form("3 planes, intersect #Delta < 0.5 cm;Plane %i cluster charge [#times 10^{3} e-];Plane %i cluster charge [#times 10^{3} e-]",fCaloPlane,i),
        qbins,0,qmax,qbins,0,qmax);
      h_clust_picky_q[i]     ->SetOption("colz");
      h_clust_picky_score[i]    = hdir.make<TH1D>(Form("pl%i_clust_picky_score",i),   Form("Plane %i clusters;Match score",i),101,0,1.01);
      
      h_clust_true_overlap[i]   = hdir.make<TH1D>(Form("pl%i_clust_true_overlap",i),  Form("Plane %i clusters (MC true matches);Overlap fraction",i),101,0,1.01);
      h_clust_true_dt[i]        = hdir.make<TH1D>(Form("pl%i_clust_true_dt",i),       Form("Plane %i clusters (MC true matches);dT [ticks]",i),200,-10,10);
      h_clust_true_dtfrac[i]      = hdir.make<TH1D>(Form("pl%i_clust_true_dtfrac",i),Form("Plane %i clusters (MC true matches);Charge-weighted mean dT/RMS",i),120,-3,3);
      h_clust_true_q[i]  = hdir.make<TH2D>(Form("pl%i_clust_true_charge",i),  
        Form("MC true matches;Plane %i cluster charge [#times 10^{3} e-];Plane %i cluster charge [#times 10^{3} e-]",fCaloPlane,i),
        qbins,0,qmax,qbins,0,qmax);
      h_clust_true_q[i]     ->SetOption("colz");
      h_clust_true_score[i]    = hdir.make<TH1D>(Form("pl%i_clust_true_score",i),   Form("Plane %i clusters;Match score",i),101,0,1.01);
      
      h_nmatches[i]         = hdir.make<TH1D>(Form("pl%i_nmatches",i),Form("number of plane%i matches to single collection cluster",i),20,0,20);
    }
  
    // Efficiency as a function of energy deposited on a wire
    //h_recoWireEff_denom = hdir.make<TH1D>("recoWireEff_trueCount","Collection plane;Electron energy deposited on wire [MeV];Count",40,0,2.0);
    //h_recoWireEff_num   = hdir.make<TH1D>("recoWireEff","Collection plane;Electron energy deposited on wire [MeV];Hit reco efficiency",40,0,2.0);
    //h_recoWireEffQ_denom = hdir.make<TH1D>("recoWireEffQ_trueCount","Collection plane;Charge deposited on wire [e-];Count",50,0,100000);
    //h_recoWireEffQ_num   = hdir.make<TH1D>("recoWireEffQ","Collection plane;Charge deposited on wire [e-];Hit reco efficiency",50,0,100000);

  }
  
  
  //--------------------------------------------------------------  
  //Destructor
  BlipRecoAlg::~BlipRecoAlg()
  {
  }
  
  
  //###########################################################
  // Reconfigure fcl parameters
  //###########################################################
  void BlipRecoAlg::reconfigure( fhicl::ParameterSet const& pset ){
    
    // initialize MC flags (will be determined automatically
    // later after the first event is processed)
    isMC        = false;
    isMCOverlay = false;
  
    fHitProducer        = pset.get<art::InputTag> ("HitProducer","");
    fHitProducerOG      = pset.get<art::InputTag> ("HitProducerOG","");
    fHitProducerData    = pset.get<art::InputTag> ("HitProducerData",    "gaushit::DataRecoStage1Test");
    fHitProducerOverlay = pset.get<art::InputTag> ("HitProducerOverlay", "gaushit::OverlayStage1a");
    if( fHitProducerOG=="" ) fHitProducerOG = fHitProducerData;
    if( fHitProducer=="" )   fHitProducer   = fHitProducerData;
    fHitTruthMatch      = pset.get<art::InputTag> ("HitTruthMatch",     "gaushitTruthMatch::OverlayRecoStage1b");
    fTrkProducer        = pset.get<std::string>   ("TrkProducer",       "pandoraInit");
    fGeantProducer      = pset.get<std::string>   ("GeantProducer",     "largeant");
    fSimDepProducer     = pset.get<std::string>   ("SimEDepProducer",   "ionization");
    fSimChanProducer    = pset.get<std::string>   ("SimChanProducer",   "driftWC:simpleSC");
    fSimGainFactor      = pset.get<float>         ("SimGainFactor",     0.826);
    fTrueBlipMergeDist  = pset.get<float>         ("TrueBlipMergeDist", 0.3);
    fMaxHitTrkLength    = pset.get<float>               ("MaxHitTrkLength", 5);
    fDoHitFiltering     = pset.get<bool>                ("DoHitFiltering",  false);
    fMaxHitMult         = pset.get<int>                 ("MaxHitMult",      10);
    fMaxHitAmp          = pset.get<float>               ("MaxHitAmp",       200);  
    fMinHitAmp          = pset.get<std::vector<float>>  ("MinHitAmp",       {-99e9,-99e9,-99e9});
    fMaxHitRMS          = pset.get<std::vector<float>>  ("MaxHitRMS",       { 99e9, 99e9, 99e9});
    fMinHitRMS          = pset.get<std::vector<float>>  ("MinHitRMS",       {-99e9,-99e9,-99e9});
    fMaxHitRatio        = pset.get<std::vector<float>>  ("MaxHitRatio",     { 99e9, 99e9, 99e9});
    fMinHitRatio        = pset.get<std::vector<float>>  ("MinHitRatio",     {-99e9,-99e9,-99e9});
    fMaxHitGOF          = pset.get<std::vector<float>>  ("MaxHitGOF",       { 99e9, 99e9, 99e9});
    fMinHitGOF          = pset.get<std::vector<float>>  ("MinHitGOF",       {-99e9,-99e9,-99e9});
    
    fHitClustWidthFact  = pset.get<float>         ("HitClustWidthFact", 3.0);
    fHitClustWireRange  = pset.get<int>           ("HitClustWireRange", 1);
    fMaxWiresInCluster  = pset.get<int>           ("MaxWiresInCluster", 10);
    fMaxClusterSpan     = pset.get<float>         ("MaxClusterSpan",    30);
    fMinClusterCharge   = pset.get<float>         ("MinClusterCharge",  500);
    fMaxClusterCharge   = pset.get<float>         ("MaxClusterCharge",  12e6);

    //fTimeOffsets        = pset.get<std::vector<float>>("TimeOffsets", {0.,0.,0.});
    fMatchMinOverlap    = pset.get<float>         ("ClustMatchMinOverlap",  0.5 );
    fMatchSigmaFact     = pset.get<float>         ("ClustMatchSigmaFact",   1.0);
    fMatchMaxTicks      = pset.get<float>         ("ClustMatchMaxTicks",    5.0 );
    fMatchQDiffLimit    = pset.get<float>         ("ClustMatchQDiffLimit",  15e3);
    fMatchMaxQRatio     = pset.get<float>         ("ClustMatchMaxQRatio",   4);
    fMatchMinScore      = pset.get<float>         ("ClustMatchMinScore",    -9);
    
    fMinMatchedPlanes   = pset.get<int>           ("MinMatchedPlanes",    2);
    fPickyBlips         = pset.get<bool>          ("PickyBlips",          false);
    fApplyTrkCylinderCut= pset.get<bool>          ("ApplyTrkCylinderCut", false);
    fCylinderRadius     = pset.get<float>         ("CylinderRadius",      15);
    fIgnoreDataTrks     = pset.get<bool>          ("IgnoreDataTrks",      false);

    fCaloPlane          = pset.get<int>           ("CaloPlane",           2);
    fCalodEdx           = pset.get<float>         ("CalodEdx",            2.8);
    fLifetimeCorr       = pset.get<bool>          ("LifetimeCorrection",  true);
    fSCECorr            = pset.get<bool>          ("SCECorrection",       true);
    fYZUniformityCorr   = pset.get<bool>          ("YZUniformityCorrection",true);
    fModBoxA            = pset.get<float>         ("ModBoxA",             0.93);
    fModBoxB            = pset.get<float>         ("ModBoxB",             0.212);
    
    fVetoBadChannels    = pset.get<bool>          ("VetoBadChannels",     true);
    fVetoNoisyChannels  = pset.get<bool>          ("VetoNoisyChannels",   false);
    fBadChanProducer    = pset.get<std::string>   ("BadChanProducer",     "nfspl1:badchannels");
    fBadChanFile        = pset.get<std::string>   ("BadChanFile",         "");
    fMinDeadWireGap     = pset.get<int>           ("MinDeadWireGap",      1);
    
    //fKeepAllClusts[0] = pset.get<bool>          ("KeepAllClustersInd", false);
    //fKeepAllClusts[1] = pset.get<bool>          ("KeepAllClustersInd", false);
    //fKeepAllClusts[2] = pset.get<bool>          ("KeepAllClustersCol", true);
    
    //keepAllClusts = true;
    //for(auto& config : fKeepAllClusts ) {
    //  if( config == false ) {
    //    keepAllClusts = false; 
    //    break;
    //  }
    //}
  }
  
  
  
  //###########################################################
  // Extract all the truth-level information from the event
  //###########################################################
  void BlipRecoAlg::RunBlipTruth( const art::Event& evt ) {
    
    //========================================
    // Reset things
    //=======================================
    pinfo.clear();
    trueblips.clear(); 
    
    auto const& wireReadout = art::ServiceHandle<geo::WireReadout>()->Get();
    auto const& chanFilt  = art::ServiceHandle<lariov::ChannelStatusService>()->GetProvider();

    // -- geometry
    art::ServiceHandle<geo::Geometry> geom;

    // -- G4 particles
    art::Handle< std::vector<simb::MCParticle> > pHandle;
    std::vector<art::Ptr<simb::MCParticle> > plist;
    if (evt.getByLabel(fGeantProducer,pHandle))
      art::fill_ptr_vector(plist, pHandle);
  
    // -- SimEnergyDeposits
    art::Handle<std::vector<sim::SimEnergyDeposit> > sedHandle;
    std::vector<art::Ptr<sim::SimEnergyDeposit> > sedlist;
    if (evt.getByLabel(fSimDepProducer,sedHandle)) 
      art::fill_ptr_vector(sedlist, sedHandle);
    
    // -- SimChannels (usually dropped in reco)
    art::Handle<std::vector<sim::SimChannel> > simchanHandle;
    std::vector<art::Ptr<sim::SimChannel> > simchanlist;
    if (evt.getByLabel(fSimChanProducer,simchanHandle)) 
      art::fill_ptr_vector(simchanlist, simchanHandle);
   
    if( plist.size() ) isMC = true;

    //====================================================
    // Prep the particle inventory service for MC+overlay
    //====================================================
    if( evt.isRealData() && plist.size() ) {
      isMCOverlay = true;
      art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
      pi_serv->Rebuild(evt);
      pi_serv->provider()->PrepParticleList(evt);
    }
    
    //=====================================================
    // Record PDG for every G4 Track ID
    //=====================================================
    //std::map<int,int> map_g4trkid_pdg;
    map_g4trkid_pdg.clear();
    map_g4trkid_chan.clear();
    map_g4trkid_chan_energy.clear();
    map_g4trkid_chan_charge.clear();
    for(size_t i = 0; i<plist.size(); i++) map_g4trkid_pdg[plist[i]->TrackId()] = plist[i]->PdgCode();
    //std::map<int, std::set<int>>         map_g4trkid_chan;
    //std::map<int, std::map<int,double> > map_g4trkid_chan_energy;
    //std::map<int, std::map<int,double> > map_g4trkid_chan_charge;

    //======================================================
    // Use SimChannels to make a map of the collected charge
    // for every G4 particle, instead of relying on the TDC-tick
    // matching that's done by BackTracker's other functions
    //======================================================
    std::map<int,double> map_g4trkid_charge;
    for(auto const &chan : simchanlist ) {
      int pl = (int)wireReadout.View(chan->Channel());
      for(auto const& tdcide : chan->TDCIDEMap() ) {
        for(auto const& ide : tdcide.second) {
          if( ide.trackID < 0 ) continue;
          double ne = ide.numElectrons;
          map_g4trkid_chan[ide.trackID].insert(chan->Channel());
          if( pl != fCaloPlane ) continue;
          
          // ####################################################
          // ###         behavior as of Nov 2022              ###
          // WireCell's detsim implements its gain "fudge factor" 
          // by scaling the SimChannel electrons (DocDB 31089)
          // instead of the electronics gain. So we need to correct 
          // for this effect to get accurate count of 'true' 
          // electrons collected on this channel.
          // ####################################################
          if( fSimGainFactor > 0 ) ne /= fSimGainFactor;
          map_g4trkid_charge[ide.trackID] += ne;
         
          // keep track of charge deposited per wire for efficiency plots
          // (coll plane only)
          if( chan->Channel() > 4800 ) {
            map_g4trkid_chan_charge[ide.trackID][chan->Channel()] += ne; 
            if( abs(map_g4trkid_pdg[ide.trackID]) == 11 ) 
              map_g4trkid_chan_energy[ide.trackID][chan->Channel()] += ide.energy;
          }
        
        }
      }
    
    }

    /*
    for(auto& m : map_g4trkid_chan_energy ) {
      for(auto& mm : m.second ) {
        if( mm.second > 0 ) h_recoWireEff_denom->Fill(mm.second);
      }
    }
    for(auto& m : map_g4trkid_chan_charge ) {
      for(auto& mm : m.second ) {
        if( mm.second > 0 ) h_recoWireEffQ_denom->Fill(mm.second);
      }
    }
    */
   

    //==================================================
    // Use G4 information to determine the "true" blips in this event.
    //==================================================
    if( plist.size() ) {
      pinfo.resize(plist.size());
      for(size_t i = 0; i<plist.size(); i++){
        BlipUtils::FillParticleInfo( *plist[i], pinfo[i], sedlist, fCaloPlane);
        if( map_g4trkid_charge[pinfo[i].trackId] ) pinfo[i].numElectrons = (int)map_g4trkid_charge[pinfo[i].trackId];
        pinfo[i].index = i;
      }
      BlipUtils::MakeTrueBlips(pinfo, trueblips);
      BlipUtils::MergeTrueBlips(trueblips, fTrueBlipMergeDist);
    }
    

    for(size_t i=0; i<trueblips.size(); i++){
      int g4id = trueblips[i].LeadG4ID;
      // loop over the channels and look for bad
      for(auto ch : map_g4trkid_chan[g4id] ) {
        if( chanFilt.IsBad(ch) ) {
          trueblips[i].AllChansGood = false;
          break;
        }
      }
    }

    


  }
  
  

  //###########################################################
  // Main reconstruction procedure.
  //
  // This function does EVERYTHING. The resulting collections of 
  // blipobj::HitClusts and blipobj::Blips can then be retrieved after
  // this function is run.
  //###########################################################
  void BlipRecoAlg::RunBlipReco( const art::Event& evt ) {
  
    //std::cout<<"\n"
    //<<"=========== BlipRecoAlg =========================\n"
    //<<"Event "<<evt.id().event()<<" / run "<<evt.id().run()<<"\n";
    
    //=======================================
    // Extract truth info if applicable
    //=======================================
    RunBlipTruth(evt);
    
    // Update hit collection for overlay (MCC10 changes)
    if( isMCOverlay && fHitProducer == fHitProducerData ) {
      fHitProducer    = fHitProducerOverlay;
      fHitProducerOG  = fHitProducerOverlay;
    }

    //=======================================
    // Reset things
    //=======================================
    blips.clear();
    hitclust.clear();
    hitinfo.clear();
    EvtBadChanCount = 0;
    
  
    //=======================================
    // Get data products for this event
    //========================================
    
    // --- detector properties
    auto const& detProp             = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataForJob();
    auto const& SCE_provider        = lar::providerFrom<spacecharge::SpaceChargeService>();
    auto const& lifetime_provider   = art::ServiceHandle<lariov::UBElectronLifetimeService>()->GetProvider();
    auto const& tpcCalib_provider   = art::ServiceHandle<lariov::TPCEnergyCalibService>()->GetProvider();
    auto const& chanFilt            = art::ServiceHandle<lariov::ChannelStatusService>()->GetProvider();
    
    //====================================================
    // Update map of bad channels for this event
    //====================================================
    fBadChanMaskPerEvt = fBadChanMask;
    if( fBadChanProducer != "" ) { 
      std::vector<int> badChans;
      art::Handle< std::vector<int>> badChanHandle;
      if( evt.getByLabel(fBadChanProducer, badChanHandle))
        badChans = *(badChanHandle);
      for(auto& ch : badChans ) {
        EvtBadChanCount++;
        fBadChanMaskPerEvt[ch] = true;
        h_chan_bad->Fill(ch);
      }
    }
  
    // -- geometry
    auto const& wireReadout = art::ServiceHandle<geo::WireReadout>()->Get();
    art::ServiceHandle<geo::Geometry> geom;
    
    // -- hits (from input module; can be track-masked subset of gaushit)
    art::Handle< std::vector<recob::Hit> > hitHandle;
    std::vector<art::Ptr<recob::Hit> > hitlist;
    if (evt.getByLabel(fHitProducer,hitHandle))
      art::fill_ptr_vector(hitlist, hitHandle);
    
    // -- hits (original), used for truth-matching
    art::Handle< std::vector<recob::Hit> > hitHandleOG;
    std::vector<art::Ptr<recob::Hit> > hitlistOG;
    if (evt.getByLabel(fHitProducerOG,hitHandleOG))
        art::fill_ptr_vector(hitlistOG, hitHandleOG);

    // -- tracks
    art::Handle< std::vector<recob::Track> > tracklistHandle;
    std::vector<art::Ptr<recob::Track> > tracklist;
    if (evt.getByLabel(fTrkProducer,tracklistHandle))
      art::fill_ptr_vector(tracklist, tracklistHandle);
    
    // -- associations
    art::FindManyP<recob::Track> fmtrk(hitHandle,evt,fTrkProducer);
    art::FindManyP<recob::Track> fmtrkOG(hitHandleOG,evt,fTrkProducer);
    
    //===============================================================
    // Map of each hit to its gaushit index (needed if the provided
    // hit collection is some filtered subset of gaushit, in order to
    // use gaushitTruthMatch later on)
    //===============================================================
    std::map< int, int > map_gh;
    // if input collection is already gaushit, this is trivial
    if( fHitProducer == fHitProducerOG ) {
      for(auto& h : hitlist ) map_gh[h.key()] = h.key(); 
    // ... but if not, find the matching origina hit. There's no convenient
    // hit ID, so we must loop through and compare channel/time (ugh)
    } else {
      std::map<int,std::vector<int>> map_chan_ghid;
      for(auto& gh : hitlistOG ) map_chan_ghid[gh->Channel()].push_back(gh.key());
      for(auto& h : hitlist ) {
        for(auto& igh : map_chan_ghid[h->Channel()]){
          if( hitlistOG[igh]->PeakTime() != h->PeakTime() ) continue;
          map_gh[h.key()] = igh;
          break;
        }
      }
    }
   

    //=======================================
    // Map track IDs to the index in the vector
    //=======================================
    //std::cout<<"Looping over tracks...\n";
    //std::map<size_t,size_t> map_trkid_isMC;
    map_trkid_isMC.clear();
    map_trkid_index.clear();
    std::map<size_t,size_t> map_trkid_nhits;
    std::map<size_t,size_t> map_trkid_nhitsMC;
    std::map<size_t,std::vector<int>> map_trkid_g4ids;
    for(size_t i=0; i<tracklist.size(); i++){ 
      map_trkid_index[tracklist.at(i)->ID()]    = i;
      map_trkid_isMC[tracklist.at(i)->ID()]     = false;
      map_trkid_nhits[tracklist.at(i)->ID()]    = 0;
      map_trkid_nhitsMC[tracklist.at(i)->ID()]  = 0;
      map_trkid_g4id[tracklist.at(i)->ID()]     = -9;
      map_trkid_g4ids[tracklist.at(i)->ID()]    .clear();
    }

    //=======================================
    // Fill vector of hit info
    //========================================
    hitinfo.resize(hitlist.size());
    std::map<int,std::vector<int>> planehitsMap;
    int nhits_untracked = 0;

    //std::cout<<"Looping over the hits... "<<hitlist.size()<<"\n";
    for(size_t i=0; i<hitlist.size(); i++){
      auto const& thisHit = hitlist[i];
      int   chan    = thisHit->Channel();
      int   plane   = thisHit->WireID().Plane;
      int   wire    = thisHit->WireID().Wire;
      
      hitinfo[i].hitid        = i;
      hitinfo[i].plane        = plane;
      hitinfo[i].chan         = chan;
      hitinfo[i].wire         = wire;
      hitinfo[i].tpc          = thisHit->WireID().TPC;
      hitinfo[i].amp          = thisHit->PeakAmplitude();
      hitinfo[i].rms          = thisHit->RMS();
      hitinfo[i].integralADC  = thisHit->Integral();
      hitinfo[i].sigmaintegral = thisHit->SigmaIntegral();
      hitinfo[i].sumADC       = thisHit->ROISummedADC();
      hitinfo[i].charge       = fCaloAlg.ElectronsFromADCArea(thisHit->Integral(),plane);
      hitinfo[i].peakTime     = thisHit->PeakTime();
      hitinfo[i].driftTime    = thisHit->PeakTime() - detProp.GetXTicksOffset(plane,0,0); // - fTimeOffsets[plane];
      hitinfo[i].gof          = thisHit->GoodnessOfFit() / thisHit->DegreesOfFreedom();
      if( thisHit->DegreesOfFreedom() ) hitinfo[i].gof = -9;

      // find associated track
      if( fmtrk.isValid() ) {
        if(fmtrk.at(i).size()) hitinfo[i].trkid = fmtrk.at(i)[0]->ID();
      } // if we didn't find an associated track, maybe this is a custom 
      //   masked hit collection or something, so let's try the OG hits
      if( hitinfo[i].trkid<0 && fHitProducer != fHitProducerOG ) {
        if( fmtrkOG.isValid() && map_gh.size() ) {
          if (fmtrkOG.at(map_gh[i]).size()) hitinfo[i].trkid = fmtrkOG.at(map_gh[i])[0]->ID(); 
        }
      }
      
      // add to the map
      planehitsMap[plane].push_back(i);
      if( hitinfo[i].trkid < 0  ) nhits_untracked++;
      if( hitinfo[i].trkid >= 0 ) map_trkid_nhits[hitinfo[i].trkid]++;
      //printf("  %lu   plane: %i,  wire: %i, time: %i\n",i,hitinfo[i].plane,hitinfo[i].wire,int(hitinfo[i].driftTime));

    }//endloop over hits




    //=======================================
    // Do hit truth-matching
    //========================================
    if( pinfo.size() ) {
      
      // Get backtracker associations
      art::FindMany<simb::MCParticle,anab::BackTrackerHitMatchingData> fmhh(hitHandleOG,evt,fHitTruthMatch);
    
      // Loop the hits
      for(size_t i=0; i<hitlist.size(); i++){
        
        //--------------------------------------------------
        // MicroBooNE-specific truth-matching: since SimChannels aren't
        // saved by default, the normal backtracker won't work, so instead
        // the truth-matching metadata is stored in the event
        //--------------------------------------------------
        int igh = map_gh[i];
        if( fmhh.at(igh).size() ) {
          std::vector<simb::MCParticle const*> pvec;
          std::vector<anab::BackTrackerHitMatchingData const*> btvec;
          fmhh.get(igh,pvec,btvec);
          hitinfo[i].g4energy = 0;
          hitinfo[i].g4charge = 0;
          float maxQ = -9;
          for(size_t j=0; j<pvec.size(); j++){
            hitinfo[i].g4energy += btvec.at(j)->energy;
            hitinfo[i].g4charge += btvec.at(j)->numElectrons;
            if( btvec.at(j)->numElectrons <= maxQ ) continue;
            maxQ = btvec.at(j)->numElectrons;
            hitinfo[i].g4trkid  = pvec.at(j)->TrackId();
            hitinfo[i].g4pdg    = pvec.at(j)->PdgCode();
            hitinfo[i].g4frac   = btvec.at(j)->ideNFraction;
          }
          
          // ###      uB behavior as of Nov 2022              ###
          // WireCell's detsim implements its gain "fudge factor" 
          // by scaling the SimChannel electrons instead of the gain
          // response. So we need to correct for this effect to get 
          // accurate count of 'true' electrons collected on channel.
          if( fSimGainFactor > 0 ) hitinfo[i].g4charge /= fSimGainFactor;
         
          /* 
          // Some old code to measure wire-by-wire thresholds
          if( map_g4trkid_chan_energy[hitinfo[i].g4trkid][chan] > 0 ) {
            double trueEnergyDep = map_g4trkid_chan_energy[hitinfo[i].g4trkid][chan];
            h_recoWireEff_num->Fill(trueEnergyDep);}
          if( map_g4trkid_chan_charge[hitinfo[i].g4trkid][chan] > 0 ) {
            double trueChargeDep = map_g4trkid_chan_charge[hitinfo[i].g4trkid][chan];
            h_recoWireEffQ_num->Fill(trueChargeDep);}
          */
        }
        
        // IF this hit was (a) matched to a track, and (b) matched to a truth
        // energy deposit, then keep the tally
        if( hitinfo[i].trkid >= 0 && hitinfo[i].g4energy > 0 ) {
          map_trkid_g4ids[hitinfo[i].trkid].push_back(hitinfo[i].g4trkid);
          map_trkid_nhitsMC[hitinfo[i].trkid]++;
        }
    
      }//endloop hits for truth-matching
    }//endif MC




    //================================================================
    // Mark tracks as MC if a considerable fraction of the total hits
    // were matched to truth depositions
    //===============================================================
    for(auto mtrk : map_trkid_index ){
      int nhitstrk = map_trkid_nhits[mtrk.first];
      int nhitstrkmc = map_trkid_nhitsMC[mtrk.first];

      float mcfrac = (nhitstrk>0) ? float(nhitstrkmc)/nhitstrk : 0;

      std::map<int,int> counters;
      for(auto g4id : map_trkid_g4ids[mtrk.first] ) {
        counters[g4id]++;
      }
      int bestID = -9;
      int bestIDcount = 0;
      for(auto counts : counters ) {
        if( counts.second > bestIDcount ) {
          bestID = counts.first;
          bestIDcount = counts.second;
        }
      }

      if( nhitstrk > 0 ) {
        h_trkhits_mcfrac->Fill(mcfrac);
        // Classify a track as "MC" if >50% of its hits
        // are matched to MC particle
        if( mcfrac > 0.50 ) {
          map_trkid_isMC[mtrk.first] = true;
          map_trkid_g4id[mtrk.first] = bestID;
        }
      }
    }
   
    //for(auto mtrk : map_trkid_index ){
      //std::cout<<"  Track ID "<<mtrk.first<<", L = "<<tracklist[mtrk.second]->Length()<<", isMC "<<map_trkid_isMC[mtrk.first]<<"\n";
    //}

    //=================================================================
    // Blip Reconstruction
    //================================================================
    //  
    //  Procedure
    //  [x] Look for hits that were not included in a track 
    //  [x] Filter hits based on hit width, etc
    //  [x] Merge together closely-spaced hits on same wires and adjacent wires
    //  [x] Plane-to-plane time matching
    //  [x] Wire intersection check to get XYZ
    //  [x] Create "blip" object

    // Create a series of masks that we'll update as we go along
    std::vector<bool> hitIsBad(hitlist.size(),      false);
    std::vector<bool> hitIsTracked(hitlist.size(),  false);
    std::vector<bool> hitIsClustered(hitlist.size(),false);
    
    // Basic track inclusion cut: exclude hits that were tracked
    for(size_t i=0; i<hitlist.size(); i++){
      if( hitinfo[i].trkid < 0 ) continue;
      auto it = map_trkid_index.find(hitinfo[i].trkid);
      if( it == map_trkid_index.end() ) continue;
      int trkindex = it->second;
      if( tracklist[trkindex]->Length() > fMaxHitTrkLength ) {
        hitIsTracked[i] = true;
        //hitIsGood[i] = false;
      }
    }
        

    // Filter based on hit properties. For hits that are a part of
    // multi-gaussian fits (multiplicity > 1), need to re-think this.
    if( fDoHitFiltering ) {
      for(size_t i=0; i<hitlist.size(); i++){
        if( hitIsTracked[i] ) continue;
        hitIsBad[i] = true;
        auto& hit = hitlist[i];
        int plane = hit->WireID().Plane;
        if( hitinfo[i].gof  <= fMinHitGOF[plane] && fMinHitGOF[plane] > 0) continue;
        if( hitinfo[i].gof  >= fMaxHitGOF[plane] && fMaxHitGOF[plane] > 0) continue;
        if( hit->RMS()            <= fMinHitRMS[plane] ) continue;
        if( hit->RMS()            >= fMaxHitRMS[plane] ) continue;
        if( hit->PeakAmplitude()  <= fMinHitAmp[plane] ) continue;
        if( hit->PeakAmplitude()  >= fMaxHitAmp )        continue;
        if( hit->Multiplicity()   >= fMaxHitMult )       continue;
        //float hit_ratio = hit->RMS() / hit->PeakAmplitude();
        //if( hit_ratio             < fMinHitRatio[plane] ) continue;
        //if( hit_ratio             > fMaxHitRatio[plane] ) continue;
        
        // we survived the gauntlet of cuts -- hit is good!
        hitIsBad[i] = false;
      }
    }

  
    // ---------------------------------------------------
    // Hit clustering
    // ---------------------------------------------------
    std::map<int,std::map<int,std::vector<int>>> tpc_planeclustsMap;
    for(auto const& planehits : planehitsMap){
      for(auto const& hi : planehits.second ){
        
        // select a new seed hit;
        // skip hits that are tracked, flagged as bad, or already clustered
        if( hitIsTracked[hi] || hitIsBad[hi] || hitIsClustered[hi] ) continue;

        // initialize a new cluster with this hit as seed
        std::set<int> hitIDs;
        hitIDs        .insert(hi);
        int startWire = hitinfo[hi].wire;
        int endWire   = hitinfo[hi].wire;
        hitIsClustered[hi] = true;
        bool clustIsValid = true;

        // see if we can add other hits to it; continue until 
        // no new hits can be lumped in with this clust
        int hitsAdded;
        do{
          hitsAdded = 0;  
          for(auto const& hj : planehits.second ) {
            
            // skip hits already clustered
            if( hitIsClustered[hj] ) continue;

            // skip hits outside overall cluster wire range
            int w1 = hitinfo[hj].wire - fHitClustWireRange;
            int w2 = hitinfo[hj].wire + fHitClustWireRange;
            if( w2 < startWire    || w1 > endWire ) continue;
            
            // check for proximity with every other hit added
            // to this cluster so far
            for(auto const& hii : hitIDs ) {

              if( hitinfo[hii].wire > w2 ) continue;
              if( hitinfo[hii].wire < w1 ) continue;
              float t1 = hitinfo[hj].driftTime;
              float t2 = hitinfo[hii].driftTime;
              float rms_sum = (hitinfo[hii].rms + hitinfo[hj].rms);
              if( fabs(t1-t2) > fHitClustWidthFact * rms_sum ) continue;
              
              // If a single bad hit is attempted to be added,
              // the entire cluster is tainted! Throw it out!
              if( hitIsBad[hj] ) { clustIsValid = false; break; }
      
              // if the hit we are checking is touching a track
              // take note of this so we can encode this info into
              // the cluster later on for delta-ray ID
              if( hitIsTracked[hj] ) {
                hitinfo[hii].touchTrk   = true;
                hitinfo[hii].touchTrkID = hitinfo[hj].trkid;
                continue;
              }
            
              startWire = std::min( hitinfo[hj].wire, startWire );
              endWire   = std::max( hitinfo[hj].wire, endWire );
              hitIDs.insert(hj);
              hitIsClustered[hj] = true;
              hitsAdded++;
              break;
            }
          
            if( !clustIsValid ) break;
          }
        } while ( hitsAdded!=0 && clustIsValid );
        
        if( !clustIsValid ) continue;

        std::vector<blipobj::HitInfo> hitinfoVec;
        for(auto hitID : hitIDs ) hitinfoVec.push_back(hitinfo[hitID]);

        blipobj::HitClust hc = BlipUtils::MakeHitClust(hitinfoVec);
        float span = hc.EndTime - hc.StartTime;
        h_clust_nwires->Fill(hc.NWires);
        h_clust_timespan->Fill(span);
          
        // basic cluster checks
        if( span      > fMaxClusterSpan   )   continue;
        if( hc.NWires > fMaxWiresInCluster )  continue;
        if( hc.Charge < fMinClusterCharge )   continue;
        if( hc.Charge > fMaxClusterCharge )   continue;
       
        // Exclude cluster if it is *entirely* on bad or noisy chans
        hc.NWiresNoisy  = 0;
        hc.NWiresBad    = 0;
        //for(auto const& hitID : hc.HitIDs ) {
        for(auto const& chan : hc.Chans ) {
          bool isBad  = ( chanFilt.IsBad(chan) || fBadChanMaskPerEvt[chan] );
          bool isNoisy= chanFilt.IsNoisy(chan);
          if( isBad   ) hc.NWiresBad++;
          if( isNoisy ) hc.NWiresNoisy++;
        }
        if( fVetoBadChannels    && hc.NWiresBad   == hc.NWires ) continue;
        if( fVetoNoisyChannels  && hc.NWiresNoisy == hc.NWires ) continue;

        // measure wire separation to nearest dead region
        // (0 = directly adjacent)
        for(size_t dw=1; dw<=5; dw++){
          int  w1   = hc.StartWire-dw;
          int  w2   = hc.EndWire+dw;
          bool flag = false;
          // treat edges of wireplane as "dead"
          if( w1 < 0 || w2 >= (int)wireReadout.Nwires(geo::PlaneID(0, hc.TPC, hc.Plane)) )
            flag=true;
          //otherwise, use channel filter service
          else {
            int ch1 = wireReadout.PlaneWireToChannel(geo::WireID(0, hc.TPC, hc.Plane, w1));
            int ch2 = wireReadout.PlaneWireToChannel(geo::WireID(0, hc.TPC, hc.Plane, w2));
            if( chanFilt.Status(ch1)<2 ) flag=true;
            if( chanFilt.Status(ch2)<2 ) flag=true;
          }
          if( flag ) { hc.DeadWireSep = dw-1; break; }
        }
       
        // veto this cluster if the gap between it and the
        // nearest dead wire (calculated above) isn't big enough
        if( fMinDeadWireGap > 0 && hc.DeadWireSep < fMinDeadWireGap ) continue;
        
        // **************************************
        // assign the ID, then go back and encode this 
        // cluster ID into the hit information
        // **************************************
        int idx = (int)hitclust.size();
        hc.ID = idx;
        tpc_planeclustsMap[hc.TPC][hc.Plane].push_back(idx);
        for(auto const& hitID : hc.HitIDs) hitinfo[hitID].clustid = hc.ID;
        // ... and find the associated truth-blip
        if( hc.G4IDs.size() ) {
          for(size_t j=0; j< trueblips.size(); j++){
            if( hc.G4IDs.count(trueblips[j].LeadG4ID)) {
              hc.EdepID = trueblips[j].ID; // we have a match!
              break;
            }
          }
        }
         
        // finally, add the finished cluster to the stack
        hitclust.push_back(hc);
      
      }
    }


    // =============================================================================
    // Plane matching and 3D blip formation
    // =============================================================================

    // --------------------------------------
    // Method 1A: Require match between calo plane ( typically collection) and
    //            1 or 2 induction planes. For every hitclust on the calo plane,
    //            do the following:
    //              1. Loop over hitclusts in one of the other planes (same TPC)
    //              3. Find closest-matched clust and add it to the histclust group
    //              4. Repeat for remaining plane(s)
    
    float _matchQDiffLimit= (fMatchQDiffLimit <= 0 ) ? std::numeric_limits<float>::max() : fMatchQDiffLimit;
    float _matchMaxQRatio = (fMatchMaxQRatio  <= 0 ) ? std::numeric_limits<float>::max() : fMatchMaxQRatio;
    float _matchTicks     = (fMatchMaxTicks   <= 0 ) ? std::numeric_limits<float>::max() : fMatchMaxTicks;
    float _matchSigmaFact = (fMatchSigmaFact  <= 0 ) ? std::numeric_limits<float>::max() : fMatchSigmaFact;

    for(auto& tpcMap : tpc_planeclustsMap ) { // loop on TPCs
     
      //std::cout
      //<<"Performing cluster matching in TPC "<<tpcMap.first<<", which has "<<tpcMap.second.size()<<" planes\n";
      auto& planeMap = tpcMap.second;
      if( planeMap.find(fCaloPlane) != planeMap.end() ){
        int   planeA              = fCaloPlane;
        auto&  hitclusts_planeA   = planeMap[planeA];
        //std::cout<<"using plane "<<fCaloPlane<<" as reference/calo plane ("<<planeMap[planeA].size()<<" clusts)\n";
        for(auto& i : hitclusts_planeA ) {
          auto& hcA = hitclust[i];
          
          // initiate hit-cluster group
          std::vector<blipobj::HitClust> hcGroup;
          hcGroup.push_back(hcA);

          // for each of the other planes, make a map of potential matches
          std::map<int, std::set<int>> cands;
         
          // map of cluster ID <--> match metrics
          std::map<int, float> map_clust_dtfrac;
          std::map<int, float> map_clust_dt;
          std::map<int, float> map_clust_overlap;
          std::map<int, float> map_clust_score;

          // ---------------------------------------------------
          // loop over other planes
          for(auto&  hitclusts_planeB : planeMap ) {
            int planeB = hitclusts_planeB.first;
            if( planeB == planeA ) continue;

            // Loop over all non-matched clusts on this plane
            for(auto const& j : hitclusts_planeB.second ) {
              auto& hcB = hitclust[j];
              if( hcB.isMatched ) continue;
              
              // *******************************************
              // Check that the two central wires intersect
              // *******************************************
              auto intersection = wireReadout.ChannelsIntersect(hcA.CenterChan, hcB.CenterChan);
              if( !intersection ) continue;
              // Save intersect location, so we don't have to
              // make another call to the Geometry service later
              TVector3 xloc(0,intersection->y,intersection->z);
              hcA.IntPts[hcB.ID] = xloc;
              hcB.IntPts[hcA.ID] = xloc;
                            
              // ***********************************
              // Calculate the cluster overlap
              // ***********************************
              float overlapFrac = BlipUtils::CalcHitClustsOverlap(hcA,hcB);
              
              // *******************************************
              // Calculate time difference for start/end, and
              // check that Q-weighted means are comparable
              // *******************************************
              float dt_start  = (hcB.StartTime - hcA.StartTime);
              float dt_end    = (hcB.EndTime   - hcA.EndTime);
              float dt        = ( fabs(dt_start) < fabs(dt_end) ) ? dt_start : dt_end;
              if( dt > _matchTicks * 20.) continue;

              // Charge-weighted mean:
              float sigmaT = std::sqrt(pow(hcA.RMS,2)+pow(hcB.RMS,2));
              float dtfrac = (hcB.Time - hcA.Time) / sigmaT;

              // *******************************************
              // Check relative charge between clusters
              // *******************************************
              float qdiff     = fabs(hcB.Charge-hcA.Charge);
              float ratio     = std::max(hcA.Charge,hcB.Charge)/std::min(hcA.Charge,hcB.Charge);

              // **************************************************
              // We made it through the cuts -- the match is good!
              // Combine metrics into a consolidated "score" that 
              // we can use later in the case of degenerate matches.
              // **************************************************
              float score = (1./3) * (overlapFrac + exp(-fabs(ratio-1.)) + exp(-fabs(dtfrac) ) );
              
              // If both clusters are matched to the same MC truth particle,
              // set flag to fill special diagnostic histograms...
              bool trueFlag = (hcA.EdepID > 0 && (hcA.EdepID == hcB.EdepID)) ? true : false; 
              
              h_clust_overlap[planeB]->Fill(overlapFrac);
              h_clust_dt[planeB]->Fill(dt);
              h_clust_dtfrac[planeB]->Fill(dtfrac);
              h_clust_score[planeB]->Fill(score);
              if( trueFlag ) {
                h_clust_true_overlap[planeB]->Fill(overlapFrac);
                h_clust_true_dt[planeB]->Fill(dt);
                h_clust_true_dtfrac[planeB]->Fill(dtfrac);
                h_clust_true_score[planeB]->Fill(score);
              }
              
              if( overlapFrac   < fMatchMinOverlap )    continue;
              if( fabs(dt)      > _matchTicks )         continue;
              if( fabs(dtfrac)  > _matchSigmaFact )     continue;
              
                             h_clust_q[planeB]      ->Fill(0.001*hcA.Charge,0.001*hcB.Charge);
              if( trueFlag ) h_clust_true_q[planeB] ->Fill(0.001*hcA.Charge,0.001*hcB.Charge);

              if( qdiff         > _matchQDiffLimit && ratio > _matchMaxQRatio ) continue;
              
              if( score         < fMatchMinScore )      continue;

              map_clust_dt[j]       = dt;
              map_clust_dtfrac[j]   = dtfrac;
              map_clust_overlap[j]  = overlapFrac;
              map_clust_score[j]    = score;
              cands[planeB]         .insert(j);
            
            }
              
          }//endloop over other planes
          
          // ---------------------------------------------------
          // loop over the candidates found on each plane
          // and select the one with the largest score
          if( cands.size() ) {
            for(auto& c : cands ) {
              int plane = c.first;
              h_nmatches[plane]->Fill(c.second.size());
              float bestScore   = -9;
              int   bestID      = -9;
              for(auto cid : c.second) {
                if( map_clust_score[cid] > bestScore ) {
                  bestScore = map_clust_score[cid];
                  bestID = cid;
                }
              }
              if( bestID >= 0 ) hcGroup.push_back(hitclust[bestID]);
            }
            
            // ----------------------------------------
            // make our new blip, but if it isn't valid, forget it and move on
            blipobj::Blip newBlip = BlipUtils::MakeBlip(hcGroup);
            if( !newBlip.isValid ) continue;
            if( newBlip.NPlanes < fMinMatchedPlanes ) continue;
            
            // ---------------------------------------
            // does this qualify as a "picky" blip?
            bool picky = ( newBlip.NPlanes > 2 && newBlip.SigmaYZ < 1. );

            // ----------------------------------------
            // save matching information
            for(auto& hc : hcGroup ) {
              hitclust[hc.ID].isMatched = true;
              for(auto hit : hitclust[hc.ID].HitIDs) hitinfo[hit].ismatch = true;
            
              // if this is a 3-plane blip with good intersection, fill diagnostic histos
              if( picky ) {
                int ipl = hc.Plane;
                if( ipl == fCaloPlane ) continue;
                float q1 = (float)newBlip.clusters[fCaloPlane].Charge;
                float q2 = (float)newBlip.clusters[ipl].Charge;
                h_clust_picky_overlap[ipl]->Fill(map_clust_overlap[hc.ID]);
                h_clust_picky_dtfrac[ipl] ->Fill(map_clust_dtfrac[hc.ID]);
                h_clust_picky_dt[ipl]     ->Fill(map_clust_dt[hc.ID]);
                h_clust_picky_q[ipl]      ->Fill(0.001*q1,0.001*q2);
                h_clust_picky_score[ipl]  ->Fill(map_clust_score[hc.ID]);
              }
            }

            if( fPickyBlips && !picky ) continue;

            
            // ----------------------------------------
            // apply cylinder cut 
            for(auto& trk : tracklist ){
              if( trk->Length() < fMaxHitTrkLength ) continue;
              if( fIgnoreDataTrks && !map_trkid_isMC[trk->ID()] ) continue;
              auto& a = trk->Vertex();
              auto& b = trk->End();
              TVector3 p1(a.X(), a.Y(), a.Z() );
              TVector3 p2(b.X(), b.Y(), b.Z() );
              // TO-DO: if this track starts or ends at a TPC boundary, 
              // we should extend p1 or p2 to outside the AV to avoid blind spots
              TVector3 bp = newBlip.Position;
              float dToLine = BlipUtils::DistToLine(p1,p2,bp);
              float d = dToLine;
              if( dToLine < 0 ) d = std::min( (bp-p1).Mag(), (bp-p2).Mag() );
              if( d > 0 ) {
                
                // update closest trkdist
                if( newBlip.ProxTrkDist < 0 || d < newBlip.ProxTrkDist ) {
                  newBlip.ProxTrkDist = d;
                  newBlip.ProxTrkID = trk->ID();
                }
                // need to do some math to figure out if this is in
                // the 45 degreee "cone" relative to the start/end
                if( dToLine > 0 && !newBlip.inCylinder && d < fCylinderRadius ) {
                  float angle1 = asin( d / (p1-bp).Mag() ) * 180./3.14159;
                  float angle2 = asin( d / (p2-bp).Mag() ) * 180./3.14159;
                  if( angle1 < 45. && angle2 < 45. ) newBlip.inCylinder = true;
                }
              }
            }//endloop over trks
           
            if( fApplyTrkCylinderCut && newBlip.inCylinder ) continue;
            
            // In the case that a cluster appeared to be touching a track in 
            // one of the 3 views, but after 3D evaluation is found to be positioned
            // far away from that track, we must update the "TouchTrk" status.
            if( (newBlip.ProxTrkID != newBlip.TouchTrkID )
            //    || (newBlip.ProxTrkDist > 1.5*newBlip.dX ) 
            //    || (newBlip.ProxTrkDist > 1.5*newBlip.dYZ ) ) 
            ) {
              newBlip.TouchTrkID = -9;
            }

            // ----------------------------------------
            // if we made it this far, the blip is good!
            // associate this blip with the hits and clusters within it
            newBlip.ID = blips.size();
            blips.push_back(newBlip);
            for(auto& hc : hcGroup ) {
              hitclust[hc.ID].BlipID = newBlip.ID;
              for( auto& h : hc.HitIDs ) hitinfo[h].blipid = newBlip.ID;
            }

  
          }//endif ncands > 0
        }//endloop over caloplane ("Plane A") clusters
      }//endif calo plane has clusters
    }//endloop over TPCs

    // Re-index the clusters after removing unmatched
    //if( !keepAllClusts ) {
    //  std::vector<blipobj::HitClust> hitclust_filt;
    //  for(size_t i=0; i<hitclust.size(); i++){
    //    auto& hc = hitclust[i];
    //    int blipID = hc.BlipID;
    //    if( fKeepAllClusts[hc.Plane] || blipID >= 0 ) {
    //      int idx = (int)hitclust_filt.size();
    //      hc.ID = idx;
    //      for( auto& h : hc.HitIDs ) hitinfo[h].clustid = hc.ID;
    //      if( blipID >= 0 ) blips[blipID].clusters[hc.Plane] = hc;
    //      hitclust_filt.push_back(hc);
    //    }
    //  }
    //  hitclust = hitclust_filt;
    //}


    for(size_t i=0; i<hitlist.size(); i++){
      if (hitinfo[i].trkid >= 0 ) continue;
      h_chan_nhits->Fill(wireReadout.PlaneWireToChannel(geo::WireID(0, 0, hitinfo[i].plane, hitinfo[i].wire)));
      
      int clustid = hitinfo[i].clustid;
      if( clustid >= 0 ) {
        if( hitclust[clustid].NWires > 1 ) continue;
        h_chan_nclusts->Fill(wireReadout.PlaneWireToChannel(geo::WireID(0, 0, hitinfo[i].plane, hitinfo[i].wire)));
      }
      //if( hitinfo[i].ismatch    ) continue;
      //if( hitclust[clustid].NWires > 1 ) continue;
      //h_chan_nclusts->Fill(wireReadout.PlaneWireToChannel(hitinfo[i].plane,hitinfo[i].wire));
    }

    
    //*************************************************************************
    // Loop over the vector of blips and perform calorimetry calculation.
    // Here we fill:
    //   blip.Charge
    //   blip.ChargeCorr (lifetime correction)
    //   blip.PositionSCE (SCE corrections)
    //   blip.Energy
    //   blip.EnergyCorr (lifetime + SCE corrections)
    //*************************************************************************
    for(size_t i=0; i<blips.size(); i++){
      auto& blip = blips[i];
      
      float Efield    = detProp.Efield();
      float EfieldSCE = detProp.Efield();

      // ----------------------------------------------------
      // Use designated calorimetry plane - defaults to collection. Here we also correct
      // for YZ non-uniformity across the wireplane, which is pretty standard (procedure 
      // taken from CalibrationdEdx_module).
      blip.Charge = blip.clusters[fCaloPlane].Charge;
      if( fYZUniformityCorr ) blip.Charge *= tpcCalib_provider.YZdqdxCorrection(fCaloPlane,blip.Position.Y(),blip.Position.Z());

      // ================================================================================
      // Calculate blip energy assuming T = T_beam (eventually can do more complex stuff
      // like associating blip with some nearby track/shower and using its tagged T0)
      //    Method 1: Assume a dE/dx ~ 2.8 MeV/cm for electrons, use that + local E-field
      //              to calculate recombination.
      //    Method 2: ESTAR lookup table method ala ArgoNeuT (TODO)
      // ================================================================================

      // --- Lifetime correction ---
      // Note: Without knowing real T0 of a blip, this correction is meaningless.
      //       Units of 'ms', not microseconds, hence the 1E-3 conversion factor.
      if( fLifetimeCorr && blip.Time>0 ) {
        float t = blip.Time*1e-3;
        float tau = lifetime_provider.Lifetime();
        blip.ChargeCorr = std::max(0.,(double)blip.Charge) * exp( t / tau ); 
      }
      
      // --- SCE corrections ---
      geo::Point_t point( blip.Position.X(),blip.Position.Y(),blip.Position.Z() );
      if( fSCECorr ) {
        
        // 1) Spatial correction
        //      TODO: Deal with cases where X falls outside AV (diffuse out-of-time signal)
        //            For example, maybe re-assign to center of drift volume?
        if( SCE_provider->EnableCalSpatialSCE() ) {
          geo::Vector_t loc_offset = SCE_provider->GetCalPosOffsets(point, 0);
          point.SetXYZ(point.X()-loc_offset.X(),point.Y()+loc_offset.Y(),point.Z()+loc_offset.Z());
          blip.PositionSCE.SetXYZ(point.X(),point.Y(),point.Z());
        }
      
        // 2) E-field correction
        //
        //   notes:
        //   - GetEfieldOffsets(xyz) and GetCalEfieldOffsets(xyz) return the exact
        //     same underlying E-field offset map; the only difference is the former
        //     is used in the simulation, and the latter in reconstruction (??).
        //   - The SpaceCharge service must have 'EnableCorSCE' and 'EnableCalEfieldSCE'
        //     enabled in order to use GetCalEfieldOffsets
        //   - Blips can have negative 'X' if the T0 correction isn't applied. Obviously 
        //     the SCE map will return (0,0,0) for these points. Beware!
        if( SCE_provider->EnableCalEfieldSCE() ) {
          auto const field_offset = SCE_provider->GetCalEfieldOffsets(point, 0); 
          EfieldSCE = Efield*std::hypot(1+field_offset.X(),field_offset.Y(),field_offset.Z());
        }

      }
      
      // METHOD 1 - assume a recombination
      float recomb    = ModBoxRecomb(fCalodEdx,Efield);
      float recombSCE = ModBoxRecomb(fCalodEdx,EfieldSCE);
      
      // nominal case + SCE/lifetime corrected
      blip.Energy     = blip.Charge     * (1./recomb)    * kWion;
      blip.EnergyCorr = blip.ChargeCorr * (1./recombSCE) * kWion;
      
      h_recomb        ->Fill(recomb);
      h_recombSCE     ->Fill(recombSCE);
      
      // METHOD 2 (TODO)
      //std::cout<<"Calculating ESTAR energy dep...  "<<depEl<<", "<<Efield<<"\n";
      //blips[i].EnergyESTAR = ESTAR->Interpolate(depEl, Efield); 
      
      // ================================================
      // Save the true blip into the object;
      // each cluster must match to the same energy dep
      // ================================================
      std::set<int> set_edepids;
      for(auto& hc : blip.clusters ) {
        if( !hc.isValid ) continue; 
        if( hc.EdepID < 0 ) break;
        set_edepids.insert( hc.EdepID );
      }
      if( set_edepids.size() == 1 )
        blip.truth = trueblips[*set_edepids.begin()];
      
    
    }//endloop over blip vector

  }//End main blip reco function
 
  
  
  //###########################################################
  float BlipRecoAlg::ModBoxRecomb(float dEdx, float Efield) {
    float Xi = fModBoxB * dEdx / ( Efield * kLArDensity );
    return log(fModBoxA+Xi)/Xi;
  }

  float BlipRecoAlg::dQdx_to_dEdx(float dQdx_e, float Efield){
    float beta  = fModBoxB / (kLArDensity * Efield);
    float alpha = fModBoxA;
    return ( exp( beta * kWion * dQdx_e ) - alpha ) / beta;
  }
  
  float BlipRecoAlg::Q_to_E(float Q, float Efield){
    if( Efield != kNominalEfield )  return kWion * (Q/ModBoxRecomb(fCalodEdx,Efield));
    else                            return kWion * (Q/kNominalRecombFactor);
  }
  
  //###########################################################
  

  
  //###########################################################
  void BlipRecoAlg::PrintConfig() {
  
    printf("BlipRecoAlg Configurations\n\n");
    printf("  Input hit collection      : %s\n",          fHitProducer.label().c_str());
    printf("  Input hit process         : %s\n",          fHitProducer.process().c_str());
    printf("  Input trk collection      : %s\n",          fTrkProducer.c_str());
    printf("  Max wires per cluster     : %i\n",          fMaxWiresInCluster);
    printf("  Max cluster timespan      : %.1f ticks\n",    fMaxClusterSpan);
    printf("  Min cluster overlap       : %4.1f\n",       fMatchMinOverlap);
    printf("  Clust match sigma-factor  : %4.1f\n",       fMatchSigmaFact);
    printf("  Clust match max dT        : %4.1f ticks\n", fMatchMaxTicks);
    printf("  Charge diff limit         : %.1fe3\n",      fMatchQDiffLimit/1000);
    printf("  Charge ratio maximum      : %.1f\n",        fMatchMaxQRatio);    
    printf("  Minimum match score       : %.2f\n",        fMatchMinScore);
    printf("  Ignoring data tracks?     : %i\n",          fIgnoreDataTrks);
    printf("  Track-cylinder radius     : %.1f cm\n",     fCylinderRadius);
    printf("  Applying cylinder cut?    : %i\n",          fApplyTrkCylinderCut);
    
    /*
    printf("  Min cluster overlap       : ");
    for(auto val : fMatchMinOverlap) { printf("%3.1f   ",val); } printf("\n");
    printf("  Clust match sigma factor  : ");
    for(auto val : fMatchSigmaFact) { printf("%3.1f   ",val); } printf("\n");
    printf("  Clust match max ticks     : ");
    for(auto val : fMatchMaxTicks) { printf("%3.1f   ",val); } printf("\n");
    */

    for(int i=0; i<kNplanes; i++){
    if( i == fCaloPlane ) continue;
    printf("  pl%i matches per cand      : %4.2f\n",       i,h_nmatches[i]->GetMean());
    }
    
    //printf("  Picky blip mode?          : %i\n",        fPickyBlips);
    printf("\n");
    
  }
  
}
