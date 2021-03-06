//brian doesn't suck
#include "../interface/ntupleProducer.h"

ntupleProducer::ntupleProducer(const edm::ParameterSet& iConfig):
  eventTree(0)
{
  jetTag_           = iConfig.getUntrackedParameter<edm::InputTag>("JetTag");
  jecTag_           = iConfig.getParameter<std::string>("JecTag");
  muonTag_          = iConfig.getUntrackedParameter<edm::InputTag>("MuonTag");
  electronTag_      = iConfig.getUntrackedParameter<edm::InputTag>("ElectronTag");
  photonTag_        = iConfig.getUntrackedParameter<edm::InputTag>("PhotonTag");
  genJetTag_        = iConfig.getUntrackedParameter<edm::InputTag>("GenJetTag");
  primaryVtxTag_    = iConfig.getUntrackedParameter<edm::InputTag>("PrimaryVtxTag");


  ebReducedRecHitCollection_ = iConfig.getParameter<edm::InputTag>("ebReducedRecHitCollection");
  eeReducedRecHitCollection_ = iConfig.getParameter<edm::InputTag>("eeReducedRecHitCollection");
  esReducedRecHitCollection_ = iConfig.getParameter<edm::InputTag>("esReducedRecHitCollection");

  //allMET:
  mMetRaw           = iConfig.getParameter<edm::InputTag>("srcMetRaw");
  mMetPf            = iConfig.getParameter<edm::InputTag>("srcMetPf");
  mMetType01        = iConfig.getParameter<edm::InputTag>("srcMetCorrected");
  mMetJERup         = iConfig.getParameter<edm::InputTag>("srcMetJERup");
  mMetJERdown       = iConfig.getParameter<edm::InputTag>("srcMetJERdown");
  mMetPhoup         = iConfig.getParameter<edm::InputTag>("srcMetPhoup");
  mMetPhodown       = iConfig.getParameter<edm::InputTag>("srcMetPhodown");
  mMetJetup         = iConfig.getParameter<edm::InputTag>("srcMetJetup");
  mMetJetdown       = iConfig.getParameter<edm::InputTag>("srcMetJetdown");
  mMetUncup         = iConfig.getParameter<edm::InputTag>("srcMetUncup");
  mMetUncdown       = iConfig.getParameter<edm::InputTag>("srcMetUncdown");
  mMetMVA           = iConfig.getParameter<edm::InputTag>("srcMVACorrected");
  mMVAJERup         = iConfig.getParameter<edm::InputTag>("srcMVAJERup");
  mMVAJERdown       = iConfig.getParameter<edm::InputTag>("srcMVAJERdown");
  mMVAPhoup         = iConfig.getParameter<edm::InputTag>("srcMVAPhoup");
  mMVAPhodown       = iConfig.getParameter<edm::InputTag>("srcMVAPhodown");
  mMVAJetup         = iConfig.getParameter<edm::InputTag>("srcMVAJetup");
  mMVAJetdown       = iConfig.getParameter<edm::InputTag>("srcMVAJetdown");
  mMVAUncup         = iConfig.getParameter<edm::InputTag>("srcMVAUncup");
  mMVAUncdown       = iConfig.getParameter<edm::InputTag>("srcMVAUncdown");

  mPhoUp            = iConfig.getParameter<edm::InputTag>("srcPatPhoup");
  mPhoDown          = iConfig.getParameter<edm::InputTag>("srcPatPhodown");
  mPhoMVAUp         = iConfig.getParameter<edm::InputTag>("srcPatPhoMvaup");
  mPhoMVADown       = iConfig.getParameter<edm::InputTag>("srcPatPhoMvadown");

  mNoOverlapJet     = iConfig.getParameter<edm::InputTag>("srcPatJets");
  mSelectedPatJet   = iConfig.getParameter<edm::InputTag>("srcSelectedJets");
  mSmearedPatJet    = iConfig.getParameter<edm::InputTag>("srcSmearedJets");

  rhoCorrTag_       = iConfig.getUntrackedParameter<edm::InputTag>("rhoCorrTag");
  rho25CorrTag_     = iConfig.getUntrackedParameter<edm::InputTag>("rho25CorrTag");
  rhoMuCorrTag_     = iConfig.getUntrackedParameter<edm::InputTag>("rhoMuCorrTag");
  hlTriggerResults_ = iConfig.getUntrackedParameter<string>("HLTriggerResults","TriggerResults");
  hltProcess_       = iConfig.getUntrackedParameter<string>("hltName");
  triggerPaths_     = iConfig.getUntrackedParameter<vector<string> >("triggers");

  partFlowTag_      = iConfig.getUntrackedParameter<edm::InputTag>("partFlowTag");
  skimLepton_       = iConfig.getUntrackedParameter<bool>("skimLepton");

  saveMuons_        = iConfig.getUntrackedParameter<bool>("saveMuons");
  saveJets_         = iConfig.getUntrackedParameter<bool>("saveJets");
  saveElectrons_    = iConfig.getUntrackedParameter<bool>("saveElectrons");
  saveEleCrystals_  = iConfig.getUntrackedParameter<bool>("saveEleCrystals");
  savePhotons_      = iConfig.getUntrackedParameter<bool>("savePhotons");
  savePhoCrystals_  = iConfig.getUntrackedParameter<bool>("savePhoCrystals");
  saveMET_          = iConfig.getUntrackedParameter<bool>("saveMET");
  saveMETExtra_     = iConfig.getUntrackedParameter<bool>("saveMETExtra");

  saveMoreEgammaVars_= iConfig.getUntrackedParameter<bool>("saveMoreEgammaVars");

  saveGenJets_      = iConfig.getUntrackedParameter<bool>("saveGenJets");
  saveGenParticles_ = iConfig.getUntrackedParameter<bool>("saveGenParticles");

  verboseTrigs       = iConfig.getUntrackedParameter<bool>("verboseTrigs");
  verboseMVAs        = iConfig.getUntrackedParameter<bool>("verboseMVAs");

  ecalTPFilterTag_    = iConfig.getUntrackedParameter<edm::InputTag>("ecalTPFilterTag");
  ecalBEFilterTag_    = iConfig.getUntrackedParameter<edm::InputTag>("ecalBEFilterTag");
  hcalHBHEFilterTag_  = iConfig.getUntrackedParameter<edm::InputTag>("hcalHBHEFilterTag");
  hcalLaserFilterTag_ = iConfig.getUntrackedParameter<edm::InputTag>("hcalLaserFilterTag");
  trackingFailureTag_ = iConfig.getUntrackedParameter<edm::InputTag>("trackingFailureTag");
  eeBadScFilterTag_   = iConfig.getUntrackedParameter<edm::InputTag>("eeBadScFilterTag");
  trkPOGFiltersTag1_  = iConfig.getUntrackedParameter<edm::InputTag>("trkPOGFiltersTag1");
  trkPOGFiltersTag2_  = iConfig.getUntrackedParameter<edm::InputTag>("trkPOGFiltersTag2");
  trkPOGFiltersTag3_  = iConfig.getUntrackedParameter<edm::InputTag>("trkPOGFiltersTag3");
  photonIsoCalcTag_   = iConfig.getParameter<edm::ParameterSet>("photonIsoCalcTag");
  jetPUIdAlgo_        = iConfig.getParameter<edm::ParameterSet>("jetPUIdAlgo");

  SCFPRemovalCone_     = iConfig.getUntrackedParameter<double>("isolation_cone_size_forSCremoval");

}

ntupleProducer::~ntupleProducer()
{
}

//
// member functions
//

// ------------ method called to for each event  ------------
void ntupleProducer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  // this is for CVS check test
  eventNumber  = iEvent.id().event();
  runNumber    = iEvent.id().run();
  lumiSection  = (unsigned int)iEvent.getLuminosityBlock().luminosityBlock();
  bunchCross   = (unsigned int)iEvent.bunchCrossing();
  isRealData   = iEvent.isRealData();

  edm::Handle<reco::BeamSpot> beamSpotHandle;
  iEvent.getByLabel("offlineBeamSpot", beamSpotHandle);
  reco::BeamSpot vertexBeamSpot = *beamSpotHandle;

  beamSpot->SetXYZ(vertexBeamSpot.x0(), vertexBeamSpot.y0(), vertexBeamSpot.z0());

  int vtxCount, jetCount, metCount, muCount, pfMuCount, eleCount, photonCount, pfPhotonCount, genCount, genPartCount, trigCount;
  vtxCount = jetCount = metCount = muCount = pfMuCount = eleCount = photonCount = pfPhotonCount = genCount = genPartCount = trigCount = 0;


  /////////////////////////////////////
  // Get PF candidates for later use //
  /////////////////////////////////////


  Handle<PFCandidateCollection> pfCands;
  iEvent.getByLabel(partFlowTag_,pfCands);
  const  PFCandidateCollection thePfColl = *(pfCands.product());

  Handle<PFCandidateCollection> pfCandsEleIso;
  iEvent.getByLabel("pfNoPileUp",pfCandsEleIso);
  const  PFCandidateCollection thePfCollEleIso = *(pfCandsEleIso.product());

  lazyTool.reset(new EcalClusterLazyTools(iEvent, iSetup, ebReducedRecHitCollection_, eeReducedRecHitCollection_));

  //////////////////////////
  //Get vertex information//
  //////////////////////////

  Handle<reco::VertexCollection> primaryVtcs;
  iEvent.getByLabel(primaryVtxTag_, primaryVtcs);

  for(VertexCollection::const_iterator iVtx = primaryVtcs->begin(); iVtx!= primaryVtcs->end(); ++iVtx){
    reco::Vertex myVtx = reco::Vertex(*iVtx);
    if(!myVtx.isValid() || myVtx.isFake()) continue;
    TCPrimaryVtx* vtxCon = new ((*primaryVtx)[vtxCount]) TCPrimaryVtx;
    vtxCon->SetXYZ(myVtx.x(), myVtx.y(), myVtx.z());
    vtxCon->SetNDof(myVtx.ndof());
    vtxCon->SetChi2(myVtx.chi2());
    vtxCon->SetNtracks(myVtx.nTracks());
    vtxCon->SetSumPt2Trks(sumPtSquared(myVtx));
    vtxCon->SetIsFake(myVtx.isFake());
    ++vtxCount;
  }

  unsigned ivtx = 0;
  VertexRef myVtxRef(primaryVtcs, ivtx);

  ///////////////////////
  //get jet information//
  ///////////////////////

  Handle<double> rhoCorr;
  iEvent.getByLabel(rhoCorrTag_, rhoCorr);
  rhoFactor = (float)(*rhoCorr);

  Handle<double> rho25Corr;
  iEvent.getByLabel(rho25CorrTag_, rho25Corr);
  rho25Factor = (float)(*rho25Corr);

  Handle<double> rhoMuCorr;
  iEvent.getByLabel(rhoMuCorrTag_, rhoMuCorr);
  rhoMuFactor = (float)(*rhoMuCorr);

  //cout<<" RHOS. In eta 4.4 = "<<rhoFactor<<"   in eta25 "<<rho25Factor<<"  MUs: "<<rhoMuFactor<<endl;

  if(saveJets_){

    edm::Handle<reco::JetTagCollection> bTagCollectionTCHE;
    iEvent.getByLabel("trackCountingHighEffBJetTags", bTagCollectionTCHE);
    const reco::JetTagCollection & bTagsTCHE = *(bTagCollectionTCHE.product());

    edm::Handle<reco::JetTagCollection> bTagCollectionTCHP;
    iEvent.getByLabel("trackCountingHighPurBJetTags", bTagCollectionTCHP);
    const reco::JetTagCollection & bTagsTCHP = *(bTagCollectionTCHP.product());

    edm::Handle<reco::JetTagCollection> bTagCollectionSSVHE;
    iEvent.getByLabel("simpleSecondaryVertexHighEffBJetTags", bTagCollectionSSVHE);
    const reco::JetTagCollection & bTagsSSVHE = *(bTagCollectionSSVHE.product());

    edm::Handle<reco::JetTagCollection> bTagCollectionJBP;
    iEvent.getByLabel("jetBProbabilityBJetTags", bTagCollectionJBP);
    const reco::JetTagCollection & bTagsJBP = *(bTagCollectionJBP.product());

    edm::Handle<reco::JetTagCollection> bTagCollectionCSV;
    iEvent.getByLabel("combinedSecondaryVertexBJetTags", bTagCollectionCSV);
    const reco::JetTagCollection & bTagsCSV = *(bTagCollectionCSV.product());

    edm::Handle<reco::JetTagCollection> bTagCollectionCSVMVA;
    iEvent.getByLabel("combinedSecondaryVertexBJetTags", bTagCollectionCSVMVA);
    const reco::JetTagCollection & bTagsCSVMVA = *(bTagCollectionCSVMVA.product());

    
    if(!isRealData){
      Handle<vector<pat::Jet> > patjets;
      iEvent.getByLabel(mNoOverlapJet, patjets);
      int jetCount_pat = 0;
      for (vector<pat::Jet>::const_iterator ipJet = patjets->begin(); ipJet != patjets->end(); ++ipJet) {
	TLorentzVector* jetCon_pat_d = new ((*jetCon_pat)[jetCount_pat]) TLorentzVector();
	jetCon_pat_d->SetPxPyPzE(ipJet->px(), ipJet->py(), ipJet->pz(), ipJet->energy());
	jetCount_pat++;
      }

      Handle<vector<pat::Jet> > smearedjets;
      iEvent.getByLabel(mSmearedPatJet, smearedjets);
      int jetCount_smear = 0;
      for (vector<pat::Jet>::const_iterator ipJet = smearedjets->begin(); ipJet != smearedjets->end(); ++ipJet) {
	TLorentzVector* jetCon_smear_d = new ((*jetCon_smear)[jetCount_smear]) TLorentzVector();
        jetCon_smear_d->SetPxPyPzE(ipJet->px(), ipJet->py(), ipJet->pz(), ipJet->energy());
	jetCount_smear++;
      }


    }

    //Have to have this because pu jet id requires to match to the ref_base only possible with View.
    edm::Handle<edm::View<reco::PFJet> > jets_pu;
    iEvent.getByLabel("ak5PFJets",jets_pu);

    Handle<ValueMap<float> > puJetIdMVA;
    iEvent.getByLabel(edm::InputTag("recoPuJetMva","fullDiscriminant"),puJetIdMVA);

    Handle<ValueMap<int> > puJetIdFlagMVA;
    iEvent.getByLabel(edm::InputTag("recoPuJetMva","fullId"),puJetIdFlagMVA);

    Handle<ValueMap<float> > puJetIdCut;
    iEvent.getByLabel(edm::InputTag("recoPuJetMva","cutbasedDiscriminant"),puJetIdCut);

    Handle<ValueMap<int> > puJetIdFlagCut;
    iEvent.getByLabel(edm::InputTag("recoPuJetMva","cutbasedId"),puJetIdFlagCut);

    Handle<vector<reco::PFJet> > jets;
    iEvent.getByLabel(jetTag_, jets);

    int i_jet = -1;
    for (vector<reco::PFJet>::const_iterator iJet = jets->begin(); iJet != jets->end(); ++iJet) {
      i_jet++;

      if (iJet->pt() < 10.) continue;

      TCJet* jetCon (new ((*recoJets)[jetCount]) TCJet);

      //float mva   = (*puJetIdMVA)[jets_pu->refAt(i_jet)];
      int  idflag = (*puJetIdFlagMVA)[jets_pu->refAt(i_jet)];

      //cout << "idflag: " << idflag << " mva: " << mva << endl;
      if( PileupJetIdentifier::passJetId( idflag, PileupJetIdentifier::kLoose )){
	jetCon->SetPuJetIdFlag_MVA_loose(1); }
      else jetCon->SetPuJetIdFlag_MVA_loose(0);
      
      if( PileupJetIdentifier::passJetId( idflag, PileupJetIdentifier::kMedium )) {
	jetCon->SetPuJetIdFlag_MVA_medium(1); }
      else jetCon->SetPuJetIdFlag_MVA_medium(0);

      if( PileupJetIdentifier::passJetId( idflag, PileupJetIdentifier::kTight )) {
	jetCon->SetPuJetIdFlag_MVA_tight(1); }
      else jetCon->SetPuJetIdFlag_MVA_tight(0);

      // Cut Based

      //float mva_cut   = (*puJetIdCut)[jets_pu->refAt(i_jet)];
      int  idflag_cut = (*puJetIdFlagCut)[jets_pu->refAt(i_jet)];

      //cout << "idflag: " << idflag << " mva: " << mva << endl;                                    
      if( PileupJetIdentifier::passJetId( idflag_cut, PileupJetIdentifier::kLoose )){
        jetCon->SetPuJetIdFlag_cut_loose(1); }
      else jetCon->SetPuJetIdFlag_cut_loose(0);
      
      if( PileupJetIdentifier::passJetId( idflag_cut, PileupJetIdentifier::kMedium )) {
        jetCon->SetPuJetIdFlag_cut_medium(1); }
      else jetCon->SetPuJetIdFlag_cut_medium(0);

      if( PileupJetIdentifier::passJetId( idflag_cut, PileupJetIdentifier::kTight )) {
        jetCon->SetPuJetIdFlag_cut_tight(1); }
      else jetCon->SetPuJetIdFlag_cut_tight(0);


      jetCon->SetPxPyPzE(iJet->px(), iJet->py(), iJet->pz(), iJet->energy());
      jetCon->SetVtx(0., 0., 0.);
      jetCon->SetChHadFrac(iJet->chargedHadronEnergyFraction());
      jetCon->SetNeuHadFrac(iJet->neutralHadronEnergyFraction());
      jetCon->SetChEmFrac(iJet->chargedEmEnergyFraction());
      jetCon->SetNeuEmFrac(iJet->neutralEmEnergyFraction());
      jetCon->SetNumConstit(iJet->chargedMultiplicity() + iJet->neutralMultiplicity());
      jetCon->SetNumChPart(iJet->chargedMultiplicity());

      //jetCon->SetJetFlavor(iJet->partonFlavour());

      jetCon->SetUncertaintyJES(-1);

      jetCon->SetBDiscriminatorMap("TCHE", MatchBTagsToJets(bTagsTCHE, *iJet));
      jetCon->SetBDiscriminatorMap("TCHP", MatchBTagsToJets(bTagsTCHP, *iJet));
      jetCon->SetBDiscriminatorMap("SSVHE", MatchBTagsToJets(bTagsSSVHE, *iJet));
      jetCon->SetBDiscriminatorMap("JBP", MatchBTagsToJets(bTagsJBP, *iJet));
      jetCon->SetBDiscriminatorMap("CSV", MatchBTagsToJets(bTagsCSV, *iJet));
      jetCon->SetBDiscriminatorMap("CSVMVA", MatchBTagsToJets(bTagsCSVMVA, *iJet));

      /////////////////////
      // Get Hgg Id vars //
      /////////////////////

      PileupJetIdentifier puIdentifier;
      // giving uncorrected input, must double check on this
      float jec = 1.;
      // jet corrector
      if( jecCor.get() == 0 ) {
        initJetEnergyCorrector( iSetup, iEvent.isRealData() );
      }
      jecCor->setJetPt(iJet->pt());
      jecCor->setJetEta(iJet->eta());
      jecCor->setJetA(iJet->jetArea());
      jecCor->setRho(rhoFactor);
      jec = jecCor->getCorrection();
      //cout<<"jec:\t"<<jec<<endl;
      VertexCollection::const_iterator vtx;
      const VertexCollection & vertexes = *(primaryVtcs.product());
      vtx = vertexes.begin();
      while( vtx != vertexes.end() && ( vtx->isFake() || vtx->ndof() < 4 ) ) {
        ++vtx;
      }
      if( vtx == vertexes.end() ) { vtx = vertexes.begin(); }
      puIdentifier = myPUJetID->computeIdVariables(&(*iJet), jec,  &(*vtx),  vertexes, false);
      //cout<<"betaStarClassic:\t"<<puIdentifier.betaStarClassic()<<"\t"<<"dR2Mean:\t"<<puIdentifier.dR2Mean()<<endl;
      jetCon->SetBetaStarClassic(puIdentifier.betaStarClassic());
      jetCon->SetDR2Mean(puIdentifier.dR2Mean());

      /////////////////////////
      // Associate to vertex //
      /////////////////////////

      associateJetToVertex(*iJet, primaryVtcs, jetCon);

      ++jetCount;
    }
  }

  if(saveMET_){

    edm::Handle<edm::View<pat::MET> > pfmet;
    iEvent.getByLabel(mMetPf, pfmet);
    if (pfmet->size() == 0) {    pfMET->SetMagPhi(-1,-10);    pfMET->SetSumEt(-1);    pfMET->SetSignificance(-1); }
    else { pfMET->SetMagPhi((*pfmet)[0].et(), (*pfmet)[0].phi());   pfMET->SetSumEt((*pfmet)[0].sumEt());   pfMET->SetSignificance((*pfmet)[0].significance());}
  
    edm::Handle<edm::View<pat::MET> > metraw;
    iEvent.getByLabel(mMetRaw, metraw);
    if (metraw->size() == 0) {    rawMET->SetMagPhi(-1,-10);    rawMET->SetSumEt(-1);    rawMET->SetSignificance(-1); }
    else { rawMET->SetMagPhi((*metraw)[0].et(), (*metraw)[0].phi());   rawMET->SetSumEt((*metraw)[0].sumEt());   rawMET->SetSignificance((*metraw)[0].significance());}

    edm::Handle<edm::View<pat::MET> > met01;
    iEvent.getByLabel(mMetType01, met01);
    if (met01->size() == 0) {    corrMET->SetMagPhi(-1,-10);    corrMET->SetSumEt(-1);    corrMET->SetSignificance(-1); }
    else { corrMET->SetMagPhi((*met01)[0].et(), (*met01)[0].phi());   corrMET->SetSumEt((*met01)[0].sumEt());   corrMET->SetSignificance((*met01)[0].significance());}

    edm::Handle<edm::View<pat::MET> > metMVA;
    iEvent.getByLabel(mMetMVA, metMVA);
    if (metMVA->size() == 0) {    mvaMET->SetMagPhi(-1,-10);    mvaMET->SetSumEt(-1);    mvaMET->SetSignificance(-1); }
    else { mvaMET->SetMagPhi((*metMVA)[0].et(), (*metMVA)[0].phi());   mvaMET->SetSumEt((*metMVA)[0].sumEt());   mvaMET->SetSignificance((*metMVA)[0].significance());}

  }

  if(saveMETExtra_){

    if (!isRealData){
      edm::Handle<pat::METCollection> metJERup;
      iEvent.getByLabel(mMetJERup, metJERup);
      if (metJERup->size() == 0) {    jerUpMET->SetMagPhi(-1,-10);    jerUpMET->SetSumEt(-1);    jerUpMET->SetSignificance(-1); }
      else { jerUpMET->SetMagPhi((*metJERup)[0].et(), (*metJERup)[0].phi()); jerUpMET->SetSumEt((*metJERup)[0].sumEt()); jerUpMET->SetSignificance((*metJERup)[0].significance());}   
      edm::Handle<pat::METCollection> metJERdown;
      iEvent.getByLabel(mMetJERdown, metJERdown);
      if (metJERdown->size() == 0) {    jerDownMET->SetMagPhi(-1,-10);    jerDownMET->SetSumEt(-1);    jerDownMET->SetSignificance(-1); }
      else { jerDownMET->SetMagPhi((*metJERdown)[0].et(), (*metJERdown)[0].phi()); jerDownMET->SetSumEt((*metJERdown)[0].sumEt()); jerDownMET->SetSignificance((*metJERdown)[0].significance());}
      
      edm::Handle<pat::METCollection> metMVAJERup;
      iEvent.getByLabel(mMVAJERup, metMVAJERup);
      if (metMVAJERup->size() == 0) {    jerUpMVAMET->SetMagPhi(-1,-10);    jerUpMVAMET->SetSumEt(-1);    jerUpMVAMET->SetSignificance(-1); }
      else { jerUpMVAMET->SetMagPhi((*metMVAJERup)[0].et(), (*metMVAJERup)[0].phi()); jerUpMVAMET->SetSumEt((*metMVAJERup)[0].sumEt()); jerUpMVAMET->SetSignificance((*metMVAJERup)[0].significance());}
      
      edm::Handle<pat::METCollection> metMVAJERdown;
      iEvent.getByLabel(mMVAJERdown, metMVAJERdown);
      if (metMVAJERdown->size() == 0) {    jerDownMVAMET->SetMagPhi(-1,-10);    jerDownMVAMET->SetSumEt(-1);    jerDownMVAMET->SetSignificance(-1); }
      else { jerDownMVAMET->SetMagPhi((*metMVAJERdown)[0].et(), (*metMVAJERdown)[0].phi()); jerDownMVAMET->SetSumEt((*metMVAJERdown)[0].sumEt()); jerDownMVAMET->SetSignificance((*metMVAJERdown)[0].significance());}
    }
  
    edm::Handle<pat::METCollection> metPhoup;
    iEvent.getByLabel(mMetPhoup, metPhoup);
    if (metPhoup->size() == 0) {    phoUpMET->SetMagPhi(-1,-10);    phoUpMET->SetSumEt(-1);    phoUpMET->SetSignificance(-1); }
    else { phoUpMET->SetMagPhi((*metPhoup)[0].et(), (*metPhoup)[0].phi()); phoUpMET->SetSumEt((*metPhoup)[0].sumEt()); phoUpMET->SetSignificance((*metPhoup)[0].significance());}

    edm::Handle<pat::METCollection> metPhodown;
    iEvent.getByLabel(mMetPhodown, metPhodown);
    if (metPhodown->size() == 0) {  phoDownMET->SetMagPhi(-1,-10);    phoDownMET->SetSumEt(-1);    phoDownMET->SetSignificance(-1); }
    else { phoDownMET->SetMagPhi((*metPhodown)[0].et(), (*metPhodown)[0].phi()); phoDownMET->SetSumEt((*metPhodown)[0].sumEt()); phoDownMET->SetSignificance((*metPhodown)[0].significance());}
    
    edm::Handle<pat::METCollection> metMVAPhoup;
    iEvent.getByLabel(mMVAPhoup, metMVAPhoup);
    if (metMVAPhoup->size() == 0) {    phoUpMVAMET->SetMagPhi(-1,-10);    phoUpMVAMET->SetSumEt(-1);    phoUpMVAMET->SetSignificance(-1); }
    else { phoUpMVAMET->SetMagPhi((*metMVAPhoup)[0].et(), (*metMVAPhoup)[0].phi()); phoUpMVAMET->SetSumEt((*metMVAPhoup)[0].sumEt()); phoUpMVAMET->SetSignificance((*metMVAPhoup)[0].significance());}

    edm::Handle<pat::METCollection> metMVAPhodown;
    iEvent.getByLabel(mMVAPhodown, metMVAPhodown);
    if (metMVAPhodown->size() == 0) {    phoDownMVAMET->SetMagPhi(-1,-10);    phoDownMVAMET->SetSumEt(-1);    phoDownMVAMET->SetSignificance(-1); }
    else { phoDownMVAMET->SetMagPhi((*metMVAPhodown)[0].et(), (*metMVAPhodown)[0].phi()); phoDownMVAMET->SetSumEt((*metMVAPhodown)[0].sumEt()); phoDownMVAMET->SetSignificance((*metMVAPhodown)[0].significance());}
    
    edm::Handle<pat::METCollection> metJetup;
    iEvent.getByLabel(mMetJetup, metJetup);
    if (metJetup->size() == 0) {    jetupMET->SetMagPhi(-1,-10);    jetupMET->SetSumEt(-1);    jetupMET->SetSignificance(-1); }
    else { jetupMET->SetMagPhi((*metJetup)[0].et(), (*metJetup)[0].phi()); jetupMET->SetSumEt((*metJetup)[0].sumEt()); jetupMET->SetSignificance((*metJetup)[0].significance());}
																
    edm::Handle<pat::METCollection> metJetdown;
    iEvent.getByLabel(mMetJetdown, metJetdown);
    if (metJetdown->size() == 0) {    jetdownMET->SetMagPhi(-1,-10);    jetdownMET->SetSumEt(-1);    jetdownMET->SetSignificance(-1); }
    else { jetdownMET->SetMagPhi((*metJetdown)[0].et(), (*metJetdown)[0].phi()); jetdownMET->SetSumEt((*metJetdown)[0].sumEt()); jetdownMET->SetSignificance((*metJetdown)[0].significance());}

    edm::Handle<pat::METCollection> metMVAJetup;
    iEvent.getByLabel(mMVAJetup, metMVAJetup);
    if (metMVAJetup->size() == 0) {    jetupMVAMET->SetMagPhi(-1,-10);    jetupMVAMET->SetSumEt(-1);    jetupMVAMET->SetSignificance(-1); }
    else { jetupMVAMET->SetMagPhi((*metMVAJetup)[0].et(), (*metMVAJetup)[0].phi()); jetupMVAMET->SetSumEt((*metMVAJetup)[0].sumEt()); jetupMVAMET->SetSignificance((*metMVAJetup)[0].significance());}

    edm::Handle<pat::METCollection> metMVAJetdown;
    iEvent.getByLabel(mMVAJetdown, metMVAJetdown);
    if (metMVAJetdown->size() == 0) {    jetdownMVAMET->SetMagPhi(-1,-10);    jetdownMVAMET->SetSumEt(-1);    jetdownMVAMET->SetSignificance(-1); }
    else { jetdownMVAMET->SetMagPhi((*metMVAJetdown)[0].et(), (*metMVAJetdown)[0].phi()); jetdownMVAMET->SetSumEt((*metMVAJetdown)[0].sumEt()); jetdownMVAMET->SetSignificance((*metMVAJetdown)[0].significance());}

    edm::Handle<pat::METCollection> metUncup;
    iEvent.getByLabel(mMetUncup, metUncup);
    if (metUncup->size() == 0) {    uncUpMET->SetMagPhi(-1,-10);    uncUpMET->SetSumEt(-1);    uncUpMET->SetSignificance(-1); }
    else { uncUpMET->SetMagPhi((*metUncup)[0].et(), (*metUncup)[0].phi()); uncUpMET->SetSumEt((*metUncup)[0].sumEt()); uncUpMET->SetSignificance((*metUncup)[0].significance());}

    edm::Handle<pat::METCollection> metUncdown;
    iEvent.getByLabel(mMetUncdown, metUncdown);
    if (metUncdown->size() == 0) {    uncDownMET->SetMagPhi(-1,-10);    uncDownMET->SetSumEt(-1);    uncDownMET->SetSignificance(-1); }
    else { uncDownMET->SetMagPhi((*metUncdown)[0].et(), (*metUncdown)[0].phi()); uncDownMET->SetSumEt((*metUncdown)[0].sumEt()); uncDownMET->SetSignificance((*metUncdown)[0].significance());}

    edm::Handle<pat::METCollection> metMVAUncup;
    iEvent.getByLabel(mMVAUncup, metMVAUncup);
    if (metMVAUncup->size() == 0) {    uncUpMVAMET->SetMagPhi(-1,-10);    uncUpMVAMET->SetSumEt(-1);    uncUpMVAMET->SetSignificance(-1); }
    else { uncUpMVAMET->SetMagPhi((*metMVAUncup)[0].et(), (*metMVAUncup)[0].phi()); uncUpMVAMET->SetSumEt((*metMVAUncup)[0].sumEt()); uncUpMVAMET->SetSignificance((*metMVAUncup)[0].significance());}

    edm::Handle<pat::METCollection> metMVAUncdown;
    iEvent.getByLabel(mMVAUncdown, metMVAUncdown);
    if (metMVAUncdown->size() == 0) {    uncDownMVAMET->SetMagPhi(-1,-10);    uncDownMVAMET->SetSumEt(-1);    uncDownMVAMET->SetSignificance(-1); }
    else { uncDownMVAMET->SetMagPhi((*metMVAUncdown)[0].et(), (*metMVAUncdown)[0].phi()); uncDownMVAMET->SetSumEt((*metMVAUncdown)[0].sumEt()); uncDownMVAMET->SetSignificance((*metMVAUncdown)[0].significance());}


  }


  ///////////////
  // Get muons //
  ///////////////


  if (saveMuons_) {

    Handle<vector<reco::Muon> > muons;
    iEvent.getByLabel(muonTag_, muons);

    for (vector<reco::Muon>::const_iterator iMuon = muons->begin(); iMuon != muons->end(); ++iMuon) {
      //if (!iMuon->isGlobalMuon() || iMuon->pt() < 3.) continue;
      //if (iMuon->pt() < 3.) continue;
      //if (!iMuon->isGlobalMuon()) continue;

      TCMuon* muCon = new ((*recoMuons)[muCount]) TCMuon;

      muCon->SetPxPyPzE(iMuon->px(), iMuon->py(), iMuon->pz(), iMuon->energy());
      muCon->SetCharge(iMuon->charge());

      muCon->SetIsPF(iMuon->isPFMuon());
      muCon->SetIsGLB(iMuon->isGlobalMuon());
      muCon->SetIsTRK(iMuon->isTrackerMuon());


      muCon->SetIsGood(muon::isGoodMuon(*iMuon, muon::TMOneStationTight));
      muCon->SetIsGoodLoose(muon::isGoodMuon(*iMuon, muon::TMOneStationLoose));

      if (primaryVtcs->size()>0){
        muCon->SetIsTight(muon::isTightMuon(*iMuon, *primaryVtcs->begin()));
        muCon->SetIsSoft( muon::isSoftMuon( *iMuon, *primaryVtcs->begin()));
      }
      else{
        muCon->SetIsTight(0);
        muCon->SetIsSoft(0);
      }

      muCon->SetCaloComp(iMuon->caloCompatibility());
      muCon->SetSegComp(muon::segmentCompatibility(*iMuon));
      muCon->SetNumberOfMatchedStations(iMuon->numberOfMatchedStations());
      muCon->SetNumberOfMatches(iMuon->numberOfMatches());

      if (iMuon->isGlobalMuon()){
        muCon->SetNormalizedChi2(       iMuon->globalTrack()->normalizedChi2());
        muCon->SetNumberOfValidMuonHits(iMuon->globalTrack()->hitPattern().numberOfValidMuonHits());
      }
      else{
        muCon->SetNormalizedChi2(-1);
        muCon->SetNumberOfValidMuonHits(-1);
      }

      if (iMuon->isTrackerMuon()){
        muCon->SetVtx(iMuon->track()->vx(),iMuon->track()->vy(),iMuon->track()->vz());
        muCon->SetPtError(iMuon->track()->ptError());

        muCon->SetTrackLayersWithMeasurement(iMuon->track()->hitPattern().trackerLayersWithMeasurement());
        muCon->SetPixelLayersWithMeasurement(iMuon->innerTrack()->hitPattern().pixelLayersWithMeasurement());
        muCon->SetNumberOfValidPixelHits(    iMuon->innerTrack()->hitPattern().numberOfValidPixelHits());
        muCon->SetNormalizedChi2_tracker(    iMuon->innerTrack()->normalizedChi2());
        muCon->SetNumberOfValidTrackerHits(iMuon->track()->hitPattern().numberOfValidTrackerHits());
        muCon->SetNumberOfLostPixelHits(   iMuon->track()->hitPattern().numberOfLostPixelHits());
        muCon->SetNumberOfLostTrackerHits( iMuon->track()->hitPattern().numberOfLostTrackerHits());
      }
      else{
        muCon->SetVtx(-1,-1,-1);
        muCon->SetPtError(-1);
        muCon->SetTrackLayersWithMeasurement(-1);
        muCon->SetNumberOfValidPixelHits(-1);
        muCon->SetNormalizedChi2_tracker(-1);
        muCon->SetNumberOfValidTrackerHits(-1);
        muCon->SetNumberOfLostPixelHits(-1);
        muCon->SetNumberOfLostTrackerHits(-1);
      }
      // Set isolation map values
      // Detector-based isolation
      muCon->SetIsoMap("NTracks_R03", iMuon->isolationR03().nTracks);
      muCon->SetIsoMap("EmIso_R03",   iMuon->isolationR03().emEt);
      muCon->SetIsoMap("HadIso_R03",  iMuon->isolationR03().hadEt);
      muCon->SetIsoMap("SumPt_R03",   iMuon->isolationR03().sumPt);

      muCon->SetIsoMap("NTracks_R05", iMuon->isolationR05().nTracks);
      muCon->SetIsoMap("EmIso_R05",   iMuon->isolationR05().emEt);
      muCon->SetIsoMap("HadIso_R05",  iMuon->isolationR05().hadEt);
      muCon->SetIsoMap("SumPt_R05",   iMuon->isolationR05().sumPt);

      // PF-based isolation
      muCon->SetIsoMap("pfPUPt_R03",      iMuon->pfIsolationR03().sumPUPt);
      muCon->SetIsoMap("pfPhotonEt_R03",  iMuon->pfIsolationR03().sumPhotonEt);
      muCon->SetIsoMap("pfChargedPt_R03", iMuon->pfIsolationR03().sumChargedParticlePt);
      muCon->SetIsoMap("pfChargedHadronPt_R03", iMuon->pfIsolationR03().sumChargedHadronPt);
      muCon->SetIsoMap("pfNeutralHadronEt_R03", iMuon->pfIsolationR03().sumNeutralHadronEt);

      muCon->SetIsoMap("pfPUPt_R04",      iMuon->pfIsolationR04().sumPUPt);
      muCon->SetIsoMap("pfPhotonEt_R04",  iMuon->pfIsolationR04().sumPhotonEt);
      muCon->SetIsoMap("pfChargedPt_R04", iMuon->pfIsolationR04().sumChargedParticlePt);
      muCon->SetIsoMap("pfChargedHadronPt_R04", iMuon->pfIsolationR04().sumChargedHadronPt);
      muCon->SetIsoMap("pfNeutralHadronEt_R04", iMuon->pfIsolationR04().sumNeutralHadronEt);

      muCon->SetPfIsoCharged(iMuon->pfIsolationR04().sumChargedHadronPt);
      muCon->SetPfIsoNeutral(iMuon->pfIsolationR04().sumNeutralHadronEt);
      muCon->SetPfIsoPhoton( iMuon->pfIsolationR04().sumPhotonEt);

      muCount++;
    }
  }


  ///////////////////
  // Get electrons //
  ///////////////////


  edm::Handle<reco::ConversionCollection> hConversions;
  iEvent.getByLabel("allConversions", hConversions);

  if (saveElectrons_) {

    Handle<reco::GsfElectronCollection > electrons;
    iEvent.getByLabel(electronTag_, electrons);

    Handle<reco::GsfElectronCollection > calibratedElectrons;
    iEvent.getByLabel(edm::InputTag("calibratedElectrons","calibratedGsfElectrons"), calibratedElectrons);

    edm::Handle<edm::ValueMap<float>> mvaTrigV0_handle;
    iEvent.getByLabel("mvaTrigV0", mvaTrigV0_handle);
    const edm::ValueMap<float> ele_mvaTrigV0 = (*mvaTrigV0_handle.product());

    edm::Handle<edm::ValueMap<double>> regEne_handle;
    iEvent.getByLabel(edm::InputTag("eleRegressionEnergy","eneRegForGsfEle"), regEne_handle);
    const edm::ValueMap<double> ele_regEne = (*regEne_handle.product());

    edm::Handle<edm::ValueMap<double>> regErr_handle;
    iEvent.getByLabel(edm::InputTag("eleRegressionEnergy","eneErrorRegForGsfEle"), regErr_handle);
    const edm::ValueMap<double> ele_regErr = (*regErr_handle.product());

    //This stuff is for modified isolation for close electrons,
    //following prescription here: https://twiki.cern.ch/twiki/bin/viewauth/CMS/BoostedZToEEModIso
    edm::Handle<edm::ValueMap<double> > h_modElectronIso_Tk;
    edm::Handle<edm::ValueMap<double> > h_modElectronIso_Ecal;
    edm::Handle<edm::ValueMap<double> > h_modElectronIso_HcalD1;
    iEvent.getByLabel("modElectronIso","track",      h_modElectronIso_Tk);
    iEvent.getByLabel("modElectronIso","ecal",       h_modElectronIso_Ecal);
    iEvent.getByLabel("modElectronIso","hcalDepth1", h_modElectronIso_HcalD1);
    const edm::ValueMap<double> modElectronIso_Tk     = (*h_modElectronIso_Tk.product());
    const edm::ValueMap<double> modElectronIso_Ecal   = (*h_modElectronIso_Ecal.product());
    const edm::ValueMap<double> modElectronIso_HcalD1 = (*h_modElectronIso_HcalD1.product());


    Int_t eee=0;
    for (vector<reco::GsfElectron>::const_iterator iElectron = electrons->begin(); iElectron != electrons->end(); ++iElectron) {
      eee++;
      if (iElectron->pt() < 5) continue;

      TCElectron* eleCon = new ((*recoElectrons)[eleCount]) TCElectron;

      // Basic physics object info
      eleCon->SetPtEtaPhiM(iElectron->pt(), iElectron->eta(), iElectron->phi(), iElectron->mass());
      eleCon->SetVtx(iElectron->gsfTrack()->vx(),iElectron->gsfTrack()->vy(),iElectron->gsfTrack()->vz());
      eleCon->SetCharge(iElectron->charge());

      // Fiducial variables
      eleCon->SetIsEB(iElectron->isEB());
      eleCon->SetIsEE(iElectron->isEE());
      eleCon->SetIsInGap(iElectron->isGap());


      // Electron ID variables

      //Methods that are availabel for the electrons can be found here:
      //http://cmslxr.fnal.gov/lxr/source/DataFormats/EgammaCandidates/interface/GsfElectron.h?v=CMSSW_5_3_11

      eleCon->SetR9(     iElectron->r9());
      eleCon->SetFBrem(  iElectron->fbrem());
      eleCon->SetEoP(    iElectron->eSuperClusterOverP());
      eleCon->SetEoPout( iElectron->eEleClusterOverPout());

      //eleCon->SetHadOverEm(iElectron->hadronicOverEm());
      // ***** >> Switching to officially recommended method:  <<<<<<
      eleCon->SetHadOverEm(iElectron->hcalOverEcalBc());
      // !!!!!!!!!
      // QUESTION: Does the eleIsolator below returns the recommended isolation for Hcal??
      //!!!!!!!!!!

      //cout<<"H/E compare: hcalOverEcalBc = "<<iElectron->hcalOverEcalBc()<<"   hadronicOverEm = "<<iElectron->hadronicOverEm()<<endl;

      eleCon->SetSCEta(  iElectron->superCluster()->eta());
      eleCon->SetSCPhi(  iElectron->superCluster()->phi());

      //*** --> Notice that previously some variables were defined in the IdMap, differently:
      //one has to perform a selections on analysis level to recover this behaviour:
      //(cut at over/underflow values)

      eleCon->SetSCDeltaEta(   iElectron->deltaEtaSuperClusterTrackAtVtx());
      eleCon->SetSCDeltaPhi(   iElectron->deltaPhiSuperClusterTrackAtVtx());
      eleCon->SetSigmaIEtaIEta(iElectron->sigmaIetaIeta());
      eleCon->SetSigmaIPhiIPhi(iElectron->sigmaIphiIphi());
      eleCon->SetSCEtaWidth(   iElectron->superCluster()->etaWidth());
      eleCon->SetSCPhiWidth(   iElectron->superCluster()->phiWidth());
      eleCon->SetSCEnergy(     iElectron->superCluster()->energy());
      if (iElectron->superCluster()->rawEnergy()!=0)
        eleCon->SetPreShowerOverRaw(iElectron->superCluster()->preshowerEnergy() / iElectron->superCluster()->rawEnergy());

      
      eleCon->SetE1x5(iElectron->e1x5());
      eleCon->SetE2x5(iElectron->e2x5Max());
      eleCon->SetE5x5(iElectron->e5x5());

      eleCon->SetDeltaEtaSeedCluster(iElectron->deltaEtaSeedClusterTrackAtCalo());
      eleCon->SetDeltaPhiSeedCluster(iElectron->deltaPhiSeedClusterTrackAtCalo());

      //std::vector vCov = iElectron->superCluster()->localCovariances(*(iElectron->superCluster()->seed())) ;


      eleCon->SetEoP(iElectron->eSuperClusterOverP());
      eleCon->SetPtError(iElectron->gsfTrack()->ptError());

      eleCon->SetNormalizedChi2Gsf(iElectron->gsfTrack()->normalizedChi2());

      bool validKF= false;
      reco::TrackRef myTrackRef = iElectron->closestCtfTrackRef();
      validKF = (myTrackRef.isAvailable());
      validKF = (myTrackRef.isNonnull());

      if (validKF){
        eleCon->SetTrackerLayersWithMeasurement( myTrackRef->hitPattern().trackerLayersWithMeasurement());
        eleCon->SetNormalizedChi2Kf( myTrackRef->normalizedChi2());
        eleCon->SetNumberOfValidHits(myTrackRef->numberOfValidHits());

      }
      else{
        eleCon->SetTrackerLayersWithMeasurement(-1);
        eleCon->SetNormalizedChi2Kf(-1);
        eleCon->SetNumberOfValidHits(-1);
      }


      InputTag  vertexLabel(string("offlinePrimaryVertices"));
      Handle<reco::VertexCollection> thePrimaryVertexColl;
      iEvent.getByLabel(vertexLabel,thePrimaryVertexColl);
      
      Vertex dummy;
      const Vertex *pv = &dummy;
      if (thePrimaryVertexColl->size() != 0) {
        pv = &*thePrimaryVertexColl->begin();
      } else { // create a dummy PV
        Vertex::Error e;
        e(0, 0) = 0.0015 * 0.0015;
        e(1, 1) = 0.0015 * 0.0015;
        e(2, 2) = 15. * 15.;
        Vertex::Point p(0, 0, 0);
        dummy = Vertex(p, e, 0, 0, 0);
      }

      float ip3d    = -999.0;
      float ip3derr = 1.0;
      float ip3dSig = 0.0;

      edm::ESHandle<TransientTrackBuilder> builder;
      iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", builder);
      TransientTrackBuilder thebuilder = *(builder.product());

      if (iElectron->gsfTrack().isNonnull()) {
        const double gsfsign   = ( (-iElectron->gsfTrack()->dxy(pv->position()))   >=0 ) ? 1. : -1.;

        const reco::TransientTrack &tt = thebuilder.build(iElectron->gsfTrack());
        const std::pair<bool,Measurement1D> &ip3dpv =  IPTools::absoluteImpactParameter3D(tt,*pv);
        if (ip3dpv.first) {
          ip3d = gsfsign*ip3dpv.second.value();
          ip3derr = ip3dpv.second.error();
          ip3dSig = ip3d/ip3derr;
        }
      }

      eleCon->SetIP3d(ip3d);
      eleCon->SetIP3dSig(ip3dSig);

      eleCon->SetNumberOfValidPixelHits(  iElectron->gsfTrack()->hitPattern().numberOfValidPixelHits());
      eleCon->SetNumberOfValidTrackerHits(iElectron->gsfTrack()->hitPattern().numberOfValidTrackerHits());
      eleCon->SetNumberOfLostPixelHits(   iElectron->gsfTrack()->hitPattern().numberOfLostPixelHits());
      eleCon->SetNumberOfLostTrackerHits( iElectron->gsfTrack()->hitPattern().numberOfLostTrackerHits());

      eleCon->SetInverseEnergyMomentumDiff(fabs((1/iElectron->ecalEnergy()) - (1/iElectron->trackMomentumAtVtx().R())));

      // HITINFO
      std::vector<TCElectron::HitInfo> hitmap;
      hitmap.clear();
      for( int h = 0; h < iElectron->gsfTrack()->hitPattern().numberOfHits(); h++) {
        uint32_t hit = iElectron->gsfTrack()->hitPattern().getHitPattern(h);
	TCElectron::HitInfo hinfo;
        hinfo.Layer = iElectron->gsfTrack()->hitPattern().getLayer(hit);
        hinfo.ValidFilter = iElectron->gsfTrack()->hitPattern().validHitFilter(hit);
        hinfo.PixelFilter = iElectron->gsfTrack()->hitPattern().pixelHitFilter(hit);
        hinfo.BarrelPixelFilter = iElectron->gsfTrack()->hitPattern().pixelBarrelHitFilter(hit);
	hinfo.MuonFilter = iElectron->gsfTrack()->hitPattern().muonHitFilter(hit);
        hinfo.StripFilter = iElectron->gsfTrack()->hitPattern().stripHitFilter(hit);
        hinfo.TrackerFilter = iElectron->gsfTrack()->hitPattern().trackerHitFilter(hit);

        hitmap.push_back(hinfo);

      } 
      eleCon->SetHitMap( hitmap );
      //end HITINFO        

      // Conversion information
      // See definition from here: https://twiki.cern.ch/twiki/bin/view/CMS/ConversionTools
      bool passConvVeto = !(ConversionTools::hasMatchedConversion(*iElectron,hConversions,vertexBeamSpot.position()));
      eleCon->SetPassConversionVeto(passConvVeto);
      eleCon->SetConversionMissHits(iElectron->gsfTrack()->trackerExpectedHitsInner().numberOfHits());

      eleIsolator.fGetIsolation(&(*iElectron), &thePfColl, myVtxRef, primaryVtcs);
      eleCon->SetIsoMap("pfChIso_R04", eleIsolator.getIsolationCharged());
      eleCon->SetIsoMap("pfNeuIso_R04",eleIsolator.getIsolationNeutral());
      eleCon->SetIsoMap("pfPhoIso_R04",eleIsolator.getIsolationPhoton());

      // Effective area for rho PU corrections (not sure if needed)
      float AEff03 = ElectronEffectiveArea::GetElectronEffectiveArea(ElectronEffectiveArea::kEleGammaAndNeutralHadronIso03, iElectron->eta(), ElectronEffectiveArea::kEleEAData2012);
      float AEff04 = ElectronEffectiveArea::GetElectronEffectiveArea(ElectronEffectiveArea::kEleGammaAndNeutralHadronIso04, iElectron->eta(), ElectronEffectiveArea::kEleEAData2012);
      eleCon->SetIsoMap("EffArea_R03", AEff03);
      eleCon->SetIsoMap("EffArea_R04", AEff04);
      eleCon->SetEffArea(AEff04);


      //MVA output:
      float m = ele_mvaTrigV0.get(eee-1);
      eleCon->SetMvaID(m);

      //Regression energy
      double ene = ele_regEne.get(eee-1);
      double err = ele_regErr.get(eee-1);
      eleCon->SetEnergyRegression(ene);
      eleCon->SetEnergyRegressionErr(err);

      //cout<<eee<<"  mva0 = "<<m<<endl;

      const reco::GsfElectron &iElectronTmp   ( (*calibratedElectrons)[eee-1]);

      //cout<<"ielectron , pt ="<<iElectron->pt()<<" eta="<<iElectron->eta()<<endl;
      //cout<<"  calibra , pt ="<<iElectronTmp.pt()<<" eta="<<iElectronTmp.eta()<<endl;

      TLorentzVector tmpP4;
      tmpP4.SetPtEtaPhiE(iElectronTmp.pt(), iElectronTmp.eta(), iElectronTmp.phi(), iElectronTmp.energy());
      eleCon->SetRegressionMomCombP4(tmpP4);


      eleIsolator.fGetIsolation(&(*iElectron), &thePfColl, myVtxRef, primaryVtcs);
      eleCon->SetPfIsoCharged(eleIsolator.getIsolationCharged());
      eleCon->SetPfIsoNeutral(eleIsolator.getIsolationNeutral());
      eleCon->SetPfIsoPhoton( eleIsolator.getIsolationPhoton());

      if (saveMoreEgammaVars_){
        eleCon->SetIsoMap("pfChIso_R04", eleIsolator.getIsolationCharged());
        eleCon->SetIsoMap("pfNeuIso_R04",eleIsolator.getIsolationNeutral());
        eleCon->SetIsoMap("pfPhoIso_R04",eleIsolator.getIsolationPhoton());

        eleCon->SetIsoMap("modIso_Tk",     modElectronIso_Tk.get(eee-1));
        eleCon->SetIsoMap("modIso_Ecal",   modElectronIso_Ecal.get(eee-1));
        eleCon->SetIsoMap("modIso_HcalD1", modElectronIso_HcalD1.get(eee-1));

        eleCon->SetIdMap("fabsEPDiff",fabs((1/iElectron->ecalEnergy()) - (1/iElectron->trackMomentumAtVtx().R())));


	


        // Electron Iso variables
        eleCon->SetIsoMap("EmIso_R03",  iElectron->dr03EcalRecHitSumEt());
        eleCon->SetIsoMap("HadIso_R03", iElectron->dr03HcalTowerSumEt());
        eleCon->SetIsoMap("SumPt_R03",  iElectron->dr03TkSumPt());

        eleCon->SetIsoMap("EmIso_R04",  iElectron->dr04EcalRecHitSumEt());
        eleCon->SetIsoMap("HadIso_R04", iElectron->dr04HcalTowerSumEt());
        eleCon->SetIsoMap("SumPt_R04",  iElectron->dr04TkSumPt());

        eleCon->SetIsoMap("pfPhotonEt_R03",      iElectron->pfIsolationVariables().photonIso);
        eleCon->SetIsoMap("pfChargedHadron_R03", iElectron->pfIsolationVariables().chargedHadronIso);
        eleCon->SetIsoMap("pfNeutralHadron_R03", iElectron->pfIsolationVariables().neutralHadronIso);

        eleCon->SetIsoMap("EffArea_R03", AEff03);
        eleCon->SetIsoMap("EffArea_R04", AEff04);
        // Add electron MVA ID and ISO variables
        
	electronMVA(&(*iElectron), eleCon, iEvent, iSetup, thePfCollEleIso, rhoFactor);

      }

      eleCount++;
    }
  }


  /////////////////
  // Get photons //
  /////////////////
  if (savePhotons_) {


    if(saveMETExtra_){
      edm::Handle<vector<pat::Photon> > patphoup;
      iEvent.getByLabel(mPhoUp,patphoup);
      
      int patphotonCount_up = 0;
      for (vector<pat::Photon>::const_iterator ipPhoton = patphoup->begin(); ipPhoton != patphoup->end() ; ++ipPhoton) {
	
	TLorentzVector* mypatPhoton_up = new ((*pho_up)[patphotonCount_up]) TLorentzVector();
	mypatPhoton_up->SetPxPyPzE(ipPhoton->px(), ipPhoton->py(), ipPhoton->pz(), ipPhoton->p());	
	++patphotonCount_up;
      }
      
      edm::Handle<vector<pat::Photon> > patphodown;
      iEvent.getByLabel(mPhoDown,patphodown);

      int patphotonCount_down = 0;
      for (vector<pat::Photon>::const_iterator ipPhoton = patphodown->begin(); ipPhoton != patphodown->end() ; ++ipPhoton) {
        TLorentzVector* mypatPhoton_down = new ((*pho_down)[patphotonCount_down]) TLorentzVector();
        mypatPhoton_down->SetPxPyPzE(ipPhoton->px(), ipPhoton->py(), ipPhoton->pz(), ipPhoton->p());
        ++patphotonCount_down;
      }

      edm::Handle<vector<pat::Photon> >patphoMVAdown;
      iEvent.getByLabel(mPhoMVADown,patphoMVAdown);

      int patphotonCountMVA_down = 0;
      for (vector<pat::Photon>::const_iterator ipPhoton = patphoMVAdown->begin(); ipPhoton != patphoMVAdown->end() ; ++ipPhoton) {
        TLorentzVector* mypatPhotonMVA_down = new ((*pho_mva_down)[patphotonCountMVA_down]) TLorentzVector();
        mypatPhotonMVA_down->SetPxPyPzE(ipPhoton->px(), ipPhoton->py(), ipPhoton->pz(), ipPhoton->p());
        ++patphotonCountMVA_down;
      }
      edm::Handle<vector<pat::Photon> >patphoMVAup;
      iEvent.getByLabel(mPhoMVAUp,patphoMVAup);

      int patphotonCountMVA_up = 0;
      for (vector<pat::Photon>::const_iterator ipPhoton = patphoMVAup->begin(); ipPhoton != patphoMVAup->end() ; ++ipPhoton) {
        TLorentzVector* mypatPhotonMVA_up = new ((*pho_mva_up)[patphotonCountMVA_up]) TLorentzVector();
        mypatPhotonMVA_up->SetPxPyPzE(ipPhoton->px(), ipPhoton->py(), ipPhoton->pz(), ipPhoton->p());
        ++patphotonCountMVA_up;
      }
 
    }


    //for mva id

    // ES geometry
    ESHandle<CaloGeometry> geoHandle;
    iSetup.get<CaloGeometryRecord>().get(geoHandle);
    const CaloSubdetectorGeometry *geometry = geoHandle->getSubdetectorGeometry(DetId::Ecal, EcalPreshower);
    const CaloSubdetectorGeometry *& geometry_p = geometry;

    if (geometry) topology_p.reset(new EcalPreshowerTopology(geoHandle));
    
    // make the map of rechits
    Handle<EcalRecHitCollection> ESRecHits;
    iEvent.getByLabel(esReducedRecHitCollection_,ESRecHits);
    
    rechits_map_.clear();
    if (ESRecHits.isValid()) {
      EcalRecHitCollection::const_iterator it;
      for (it = ESRecHits->begin(); it != ESRecHits->end(); ++it) {
        // remove bad ES rechits
        if (it->recoFlag()==1 || it->recoFlag()==14 || (it->recoFlag()<=10 && it->recoFlag()>=5)) continue;
        //Make the map of DetID, EcalRecHit pairs
        rechits_map_.insert(std::make_pair(it->id(), *it));
      }
    }
    
    Handle<EcalRecHitCollection> Brechit;
    iEvent.getByLabel("reducedEcalRecHitsEB",Brechit);
    //const EcalRecHitCollection* barrelRecHits= Brechit.product();

    Handle<vector<reco::Photon> > photons;
    iEvent.getByLabel(photonTag_, photons);

    edm::Handle<reco::GsfElectronCollection> hElectrons;
    iEvent.getByLabel("gsfElectrons", hElectrons);

    for (vector<reco::Photon>::const_iterator iPhoton = photons->begin(); iPhoton != photons->end() ; ++iPhoton) {
      TCPhoton* myPhoton = new ((*recoPhotons)[photonCount]) TCPhoton();

      if (savePhoCrystals_)
        {
          //Crystal Info:
          std::vector< std::pair<DetId, float> >  PhotonHit_DetIds  = iPhoton->superCluster()->hitsAndFractions();
          std::vector<TCEGamma::CrystalInfo> crystalinfo_container;
          crystalinfo_container.clear();
          TCPhoton::CrystalInfo crystal = {};
          float timing_avg =0.0;
          int ncrys   = 0;
          vector< std::pair<DetId, float> >::const_iterator detitr;

          for(detitr = PhotonHit_DetIds.begin(); detitr != PhotonHit_DetIds.end(); ++detitr)
            {

              if (((*detitr).first).det() == DetId::Ecal && ((*detitr).first).subdetId() == EcalBarrel) {
                EcalRecHitCollection::const_iterator j= Brechit->find(((*detitr).first));
                EcalRecHitCollection::const_iterator thishit;
                if ( j!= Brechit->end())  thishit = j;
                if ( j== Brechit->end()){
                  continue;
                }

                EBDetId detId  = (EBDetId)((*detitr).first);
                crystal.rawId  = thishit->id().rawId();
                crystal.energy = thishit->energy();
                crystal.time   = thishit->time();
                crystal.timeErr= thishit->timeError();
                crystal.recoFlag = thishit->recoFlag();
                crystal.ieta   = detId.ieta();
                crystal.iphi   = detId.iphi();
                if(crystal.energy > 0.1){
                  timing_avg  = timing_avg + crystal.time;
                  ncrys++;
                }
              }//end of if ((*detitr).det() == DetId::Ecal && (*detitr).subdetId() == EcalBarrel)
              crystalinfo_container.push_back(crystal);
            }//End loop over detids
          std::sort(crystalinfo_container.begin(),crystalinfo_container.end(),EnergySortCriterium);


          //Without taking into account uncertainty, this time makes no sense.
          if (ncrys !=0) timing_avg = timing_avg/(float)ncrys;
          else timing_avg = -99.;

          myPhoton->SetNCrystals(crystalinfo_container.size());

          for (unsigned int y =0; y < crystalinfo_container.size() && y < 100;y++){
            myPhoton->AddCrystal(crystalinfo_container[y]);
          }

      /*
         vector<TCPhoton::CrystalInfo> savedCrystals = myPhoton->GetCrystalVect();
         for (int y = 0; y< myPhoton->GetNCrystals();y++){
         std::cout << "savedCrystals[y].time : " << savedCrystals[y].time << std::endl; 
         std::cout << "savedCrystals[y].timeErr : " << savedCrystals[y].timeErr << std::endl;
         std::cout << "savedCrystals[y].energy : " << savedCrystals[y].energy <<std::endl;
         std::cout << "savedCrystals[y].ieta: " << savedCrystals[y].ieta << std::endl;

         std::cout << "savedCrystals[y].rawId: " << savedCrystals[y].rawId <<std::endl;
         }
         */

      //const reco::BasicCluster& seedClus = *(iPhoton->superCluster()->seed());

          //const reco::BasicCluster& seedClus = *(iPhoton->superCluster()->seed());
        }

      myPhoton->SetPxPyPzE(iPhoton->px(), iPhoton->py(), iPhoton->pz(), iPhoton->p());
      myPhoton->SetVtx(iPhoton->vx(), iPhoton->vy(), iPhoton->vz());

      myPhoton->SetMipChi2(iPhoton->mipChi2());
      myPhoton->SetMipTotEn(iPhoton->mipTotEnergy());
      myPhoton->SetMipSlope(iPhoton->mipSlope());
      myPhoton->SetMipIntercept(iPhoton->mipIntercept());
      myPhoton->SetMipNHitCone(iPhoton->mipNhitCone());
      myPhoton->SetMipIsHalo(iPhoton->mipIsHalo());

      /*
      //Mip Variables:
      cout << "mipChi2: "<< iPhoton->mipChi2() << endl;
      cout << "mipTotEnergy: " << iPhoton->mipTotEnergy() << endl;
      cout << "mipSlope: " << iPhoton->mipSlope() << endl;
      cout << "mipIntercept: " << iPhoton->mipIntercept() <<endl;
      cout << "mipNhitCone: " << iPhoton->mipNhitCone() << endl;
      cout << "mipIsHalo: " << iPhoton->mipIsHalo() << endl;
      */

      // more cluster shapes from Lazy Tools

      const EcalRecHitCollection* barrelRecHits= Brechit.product();      
      if(iPhoton->isEB()){
	vector<float> showershapes_barrel = EcalClusterTools::roundnessBarrelSuperClusters(*(iPhoton->superCluster()),*barrelRecHits,0);
	
	//std::cout << "roundness: " << (float)showershapes_barrel[0] << std::endl;
	//std::cout << "angle: " << (float)showershapes_barrel[1] << std::endl;
	
	myPhoton->SetRoundness((float)showershapes_barrel[0]);
	myPhoton->SetAngle((float)showershapes_barrel[1]);
      }
      else{myPhoton->SetRoundness(-99.); myPhoton->SetAngle(-99.);}
      
      
      vector<float> phoCov;
      const reco::CaloClusterPtr phoSeed = iPhoton->superCluster()->seed();
      phoCov = lazyTool->localCovariances(*phoSeed);

      std::pair<DetId, float> maxRH = EcalClusterTools::getMaximum( *phoSeed, &(*barrelRecHits) );
      if(maxRH.second) {
	Cluster2ndMoments moments = EcalClusterTools::cluster2ndMoments(*phoSeed, *barrelRecHits);
	//cout << "smaj: " << moments.sMaj << std::endl;
	//cout << "smin: " << moments.sMin << std::endl;
	myPhoton->SetSMin(moments.sMin);
	myPhoton->SetSMaj(moments.sMaj);
      }
      else {myPhoton->SetSMin(-99); myPhoton->SetSMaj(-99);}

      // ID variables
      //Methods that are availabel for the electrons can be found here:
      //http://cmslxr.fnal.gov/lxr/source/DataFormats/EgammaCandidates/interface/Photon.h?v=CMSSW_5_3_11
      myPhoton->SetHadOverEm(iPhoton->hadTowOverEm());
      myPhoton->SetR9(iPhoton->r9());
      myPhoton->SetTrackVeto(iPhoton->hasPixelSeed());

      myPhoton->SetSCEta(iPhoton->superCluster()->eta());
      myPhoton->SetSCPhi(iPhoton->superCluster()->phi());
      myPhoton->SetSigmaIEtaIEta(iPhoton->sigmaIetaIeta());
      myPhoton->SetSigmaIEtaIPhi(phoCov[1]);
      myPhoton->SetSigmaIPhiIPhi(phoCov[2]); 

      myPhoton->SetSCEtaWidth(  iPhoton->superCluster()->etaWidth());
      myPhoton->SetSCPhiWidth(  iPhoton->superCluster()->phiWidth());

      //these don't exist either
      //eleCon->SetSCDeltaEta( );
      //eleCon->SetSCDeltaPhi( );

      myPhoton->SetSCEnergy(iPhoton->superCluster()->energy());
      myPhoton->SetSCRawEnergy(iPhoton->superCluster()->rawEnergy());
      myPhoton->SetSCPSEnergy(iPhoton->superCluster()->preshowerEnergy());

      if (iPhoton->superCluster()->rawEnergy()!=0)
	myPhoton->SetPreShowerOverRaw(iPhoton->superCluster()->preshowerEnergy() / iPhoton->superCluster()->rawEnergy());

      myPhoton->SetE1x3(lazyTool->e3x1(*phoSeed));
      myPhoton->SetE2x2(lazyTool->e2x2(*phoSeed));
      myPhoton->SetE2x5Max(lazyTool->e2x5Max(*phoSeed));

      myPhoton->SetE1x5(iPhoton->e1x5());
      myPhoton->SetE2x5(iPhoton->e2x5());
      myPhoton->SetE5x5(iPhoton->e5x5());

      // PF Iso for photons
      phoIsolator.fGetIsolation(&(*iPhoton),&thePfColl, myVtxRef, primaryVtcs);
      myPhoton->SetPfIsoCharged(phoIsolator.getIsolationCharged());
      myPhoton->SetPfIsoNeutral(phoIsolator.getIsolationNeutral());
      myPhoton->SetPfIsoPhoton( phoIsolator.getIsolationPhoton());

      if (saveMoreEgammaVars_){
        myPhoton->SetIsoMap("chIso03",phoIsolator.getIsolationCharged());
        myPhoton->SetIsoMap("nhIso03",phoIsolator.getIsolationNeutral());
        myPhoton->SetIsoMap("phIso03",phoIsolator.getIsolationPhoton());

        // detector-based isolation
        myPhoton->SetIsoMap("EmIso_R03",  (iPhoton->ecalRecHitSumEtConeDR03()));
        myPhoton->SetIsoMap("HadIso_R03", (iPhoton->hcalTowerSumEtConeDR03()));
        myPhoton->SetIsoMap("TrkIso_R03", (iPhoton->trkSumPtHollowConeDR03()));


        myPhoton->SetIsoMap("EmIso_R04",  (iPhoton->ecalRecHitSumEtConeDR04()));
        myPhoton->SetIsoMap("HadIso_R04", (iPhoton->hcalTowerSumEtConeDR04()));
        myPhoton->SetIsoMap("TrkIso_R04", (iPhoton->trkSumPtHollowConeDR04()));
      }

      // Hcal isolation for 2012
      //myPhoton->SetIsoMap("HadIso_R03",iPhoton->hcalTowerSumEtConeDR03() +
      //        (iPhoton->hadronicOverEm() - iPhoton->hadTowOverEm())*iPhoton->superCluster()->energy()/cosh(iPhoton->superCluster()->eta()));
      //myPhoton->SetIsoMap("HadIso_R04",iPhoton->hcalTowerSumEtConeDR04() +
      //        (iPhoton->hadronicOverEm() - iPhoton->hadTowOverEm())*iPhoton->superCluster()->energy()/cosh(iPhoton->superCluster()->eta()));

      //Footprint removal
      edm::ParameterSet myConfig;
      myConfig.addUntrackedParameter("isolation_cone_size_forSCremoval",SCFPRemovalCone_);
      SuperClusterFootprintRemoval remover(iEvent,iSetup,myConfig);
      PFIsolation_struct mySCFPstruct = remover.PFIsolation(iPhoton->superCluster(),edm::Ptr<Vertex>(primaryVtcs,ivtx));
      /*
      cout<<"chargediso: "<<mySCFPstruct.chargediso<<endl;
      cout<<"chargediso_primvtx: "<<mySCFPstruct.chargediso_primvtx<<endl;
      cout<<"neutraliso: "<<mySCFPstruct.neutraliso<<endl;
      cout<<"photoniso: "<<mySCFPstruct.photoniso<<endl;
      cout<<"chargediso_rcone: "<<mySCFPstruct.chargediso_rcone<<endl;
      cout<<"chargediso_primvtx_rcone: "<<mySCFPstruct.chargediso_primvtx_rcone<<endl;
      cout<<"neutraliso_rcone: "<<mySCFPstruct.neutraliso_rcone<<endl;
      cout<<"photoniso_rcone: "<<mySCFPstruct.photoniso_rcone<<endl;
      cout<<"eta_rcone: "<<mySCFPstruct.eta_rcone<<endl;
      cout<<"phi_rcone: "<<mySCFPstruct.phi_rcone<<endl;
      cout<<"rcone_isOK: "<<mySCFPstruct.rcone_isOK<<endl;
      */
      myPhoton->SetIsoMap("SCFP_chargediso",mySCFPstruct.chargediso);
      myPhoton->SetIsoMap("SCFP_chargediso_primvtx",mySCFPstruct.chargediso_primvtx);
      myPhoton->SetIsoMap("SCFP_neutraliso",mySCFPstruct.neutraliso);
      myPhoton->SetIsoMap("SCFP_photoniso",mySCFPstruct.photoniso);
      myPhoton->SetIsoMap("SCFP_chargediso_rcone",mySCFPstruct.chargediso_rcone);
      myPhoton->SetIsoMap("SCFP_chargediso_primvtx_rcone",mySCFPstruct.chargediso_primvtx_rcone);
      myPhoton->SetIsoMap("SCFP_neutraliso_rcone",mySCFPstruct.neutraliso_rcone);
      myPhoton->SetIsoMap("SCFP_photoniso_rcone",mySCFPstruct.photoniso_rcone);
      myPhoton->SetIsoMap("SCFP_eta_rcone",mySCFPstruct.eta_rcone);
      myPhoton->SetIsoMap("SCFP_phi_rcone",mySCFPstruct.phi_rcone);
      myPhoton->SetIsoMap("SCFP_rcone_isOK",mySCFPstruct.rcone_isOK);

      //Conversion info
      bool passElectronVeto = !(ConversionTools::hasMatchedPromptElectron(iPhoton->superCluster(), hElectrons, hConversions, vertexBeamSpot.position()));
      myPhoton->SetConversionVeto(passElectronVeto);

      //Effective energy shit
      float phoESEffSigmaRR_x = 0.;
      float phoESEffSigmaRR_y = 0.;
      float phoESEffSigmaRR_z = 0.;
      
      if (ESRecHits.isValid() && (fabs(iPhoton->superCluster()->eta()) > 1.6 && fabs(iPhoton->superCluster()->eta()) < 3)) {

        vector<float> phoESHits0 = getESHits((*iPhoton).superCluster()->x(), (*iPhoton).superCluster()->y(), (*iPhoton).superCluster()->z(), rechits_map_, geometry_p, topology_p.get(), 0);

        vector<float> phoESShape = getESEffSigmaRR(phoESHits0);
        phoESEffSigmaRR_x = phoESShape[0];
        phoESEffSigmaRR_y = phoESShape[1];
        phoESEffSigmaRR_z = phoESShape[2];
      }
      myPhoton->SetESEffSigmaRR(phoESEffSigmaRR_x, phoESEffSigmaRR_y, phoESEffSigmaRR_z);


      ++photonCount;
    }
  }


  ////////////////////////
  // Get gen-level info //
  ////////////////////////


  if (!isRealData) {

    Handle<GenEventInfoProduct> GenEventInfoHandle;
    iEvent.getByLabel("generator", GenEventInfoHandle);

    evtWeight = ptHat = qScale = -1;

    if (GenEventInfoHandle.isValid()) {
      //qScale       = GenEventInfoHandle->qScale();
      //evtWeight    = GenEventInfoHandle->weight();
      ptHat        = (GenEventInfoHandle->hasBinningValues() ? GenEventInfoHandle->binningValues()[0] : 0.0);
    }


    ////////////////////
    // PU information //
    ////////////////////

    Handle<std::vector< PileupSummaryInfo > >  PUInfo;
    iEvent.getByLabel(edm::InputTag("addPileupInfo"), PUInfo);
    std::vector<PileupSummaryInfo>::const_iterator iPV;

    for(iPV = PUInfo->begin(); iPV != PUInfo->end(); ++iPV){
      if (iPV->getBunchCrossing() == 0){
        nPUVertices     = iPV->getPU_NumInteractions();
        nPUVerticesTrue = iPV->getTrueNumInteractions();
      }
    }

    //////////////////////
    // Get genParticles //
    //////////////////////

    if (saveGenParticles_) {
      Handle<GenParticleCollection> genParticleColl;
      iEvent.getByLabel("genParticles", genParticleColl);

      map<const reco::GenParticle*, TCGenParticle*> genMap;
      for (GenParticleCollection::const_iterator myParticle= genParticleColl->begin(); myParticle != genParticleColl->end(); ++myParticle) {

        ////  Leptons and photons and b's, (oh my)
        //// Z's, W's, H's, and now big juicy Gravitons
        if (
            (abs(myParticle->pdgId()) >= 11 && abs(myParticle->pdgId()) <= 16) 
            || myParticle->pdgId() == 22 
            || abs(myParticle->pdgId()) == 5 
            || abs(myParticle->pdgId()) == 23 
            || abs(myParticle->pdgId()) == 24 
            || abs(myParticle->pdgId()) == 25   //higgs
            || abs(myParticle->pdgId()) == 35   // another higgs
            || abs(myParticle->pdgId()) == 36   // more higgses
            || abs(myParticle->pdgId()) == 39   //graviton (sometimes higgs too)
            || abs(myParticle->pdgId()) == 443  //jpsi
            || abs(myParticle->pdgId()) == 553  //upsilon
           ) {
          addGenParticle(&(*myParticle), genPartCount, genMap);

        }

      }

    }


    /////////////////
    // Get genJets //
    /////////////////

    if (saveGenJets_) {

      Handle<reco::GenJetCollection> GenJets;
      iEvent.getByLabel(genJetTag_, GenJets);

      for (GenJetCollection::const_iterator iJet = GenJets->begin(); iJet!= GenJets->end(); ++iJet) {
        reco::GenJet myJet = reco::GenJet(*iJet);
        if (myJet.pt() > 10) {
          TCGenJet* jetCon = new ((*genJets)[genCount]) TCGenJet;
          jetCon->SetPxPyPzE(myJet.px(), myJet.py(), myJet.pz(), myJet.energy());
          jetCon->SetHadEnergy(myJet.hadEnergy());
          jetCon->SetEmEnergy(myJet.emEnergy());
          jetCon->SetInvEnergy(myJet.invisibleEnergy());
          jetCon->SetAuxEnergy(myJet.auxiliaryEnergy());
          jetCon->SetNumConstit(myJet.getGenConstituents().size());
          jetCon->SetJetFlavor(0);
        }
        ++genCount;
      }
    }
  }


  ///////////////////
  // Noise filters //
  ///////////////////

  //if (isRealData) {

  myNoiseFilters.isScraping = false; //isFilteredOutScraping(iEvent, iSetup, 10, 0.25);

  Handle<bool> hcalNoiseFilterHandle;
  iEvent.getByLabel(hcalHBHEFilterTag_, hcalNoiseFilterHandle);
  if (hcalNoiseFilterHandle.isValid())  myNoiseFilters.isNoiseHcalHBHE = !(Bool_t)(*hcalNoiseFilterHandle);
  else LogWarning("Filters")<<"hcal noise NOT valid  ";

  Handle<bool> hcalLaserFilterHandle;
  iEvent.getByLabel(hcalLaserFilterTag_, hcalLaserFilterHandle);
  if (hcalLaserFilterHandle.isValid())  myNoiseFilters.isNoiseHcalLaser = !(Bool_t)(*hcalLaserFilterHandle);
  else LogWarning("Filters")<<"hcal Laser NOT valid  ";

  Handle<bool> ecalTPFilterHandle;
  iEvent.getByLabel(ecalTPFilterTag_, ecalTPFilterHandle);
  if (ecalTPFilterHandle.isValid())  myNoiseFilters.isNoiseEcalTP = !(Bool_t)(*ecalTPFilterHandle);
  else LogWarning("Filters")<<"Ecal TP NOT valid  ";

  Handle<bool> ecalBEFilterHandle;
  iEvent.getByLabel(ecalBEFilterTag_, ecalBEFilterHandle);
  if (ecalBEFilterHandle.isValid())  myNoiseFilters.isNoiseEcalBE = !(Bool_t)(*ecalBEFilterHandle);
  else LogWarning("Filters")<<"Ecal BE NOT valid  ";

  Handle<bool> trackingFailureHandle;
  iEvent.getByLabel(trackingFailureTag_, trackingFailureHandle);
  if(trackingFailureHandle.isValid()) myNoiseFilters.isNoiseTracking = !(Bool_t)(*trackingFailureHandle);
  else LogWarning("Filters")<<"tracking Failure NOT valid  ";

  Handle<bool> eeBadScFilterHandle;
  iEvent.getByLabel(eeBadScFilterTag_, eeBadScFilterHandle);
  if(eeBadScFilterHandle.isValid()) myNoiseFilters.isNoiseEEBadSc = !(Bool_t)(*eeBadScFilterHandle);
  else LogWarning("Filters")<<"eeBadSc NOT valid  ";

  Handle<bool> trkPOGFiltersHandle1;
  iEvent.getByLabel(trkPOGFiltersTag1_, trkPOGFiltersHandle1);
  if(trkPOGFiltersHandle1.isValid()) myNoiseFilters.isNoisetrkPOG1 = !(Bool_t)(*trkPOGFiltersHandle1);
  else LogWarning("Filters")<<"trkPOG1 NOT valid  ";

  Handle<bool> trkPOGFiltersHandle2;
  iEvent.getByLabel(trkPOGFiltersTag2_, trkPOGFiltersHandle2);
  if(trkPOGFiltersHandle2.isValid()) myNoiseFilters.isNoisetrkPOG2 = !(Bool_t)(*trkPOGFiltersHandle2);
  else LogWarning("Filters")<<"trkPOG2 NOT valid  ";

  Handle<bool> trkPOGFiltersHandle3;
  iEvent.getByLabel(trkPOGFiltersTag3_, trkPOGFiltersHandle3);
  if(trkPOGFiltersHandle3.isValid()) myNoiseFilters.isNoisetrkPOG3 = !(Bool_t)(*trkPOGFiltersHandle3);
  else LogWarning("Filters")<<"trkPOG3 NOT valid  ";


  edm::Handle<BeamHaloSummary> TheBeamHaloSummary;
  iEvent.getByLabel("BeamHaloSummary",TheBeamHaloSummary);
  const BeamHaloSummary TheSummary = (*TheBeamHaloSummary.product() );
  if (!TheBeamHaloSummary.isValid()) LogWarning("Filters")<<"The Summary (for CSC halo) NOT valid  ";

  myNoiseFilters.isCSCTightHalo = TheSummary.CSCTightHaloId();
  myNoiseFilters.isCSCLooseHalo = TheSummary.CSCLooseHaloId();

  //LogWarning("Filters")<<"\n csc1  "<< myNoiseFilters.isCSCTightHalo<<"  csc2  "<<myNoiseFilters.isCSCLooseHalo
  //  <<" isNoiseHcal HBHE "<<myNoiseFilters.isNoiseHcalHBHE<<"  laser "<<myNoiseFilters.isNoiseHcalLaser<<"\n"
  //  <<" ecal TP  "<<myNoiseFilters.isNoiseEcalTP<<"   ecal BE  "<<myNoiseFilters.isNoiseEcalBE;

  //}

  ////////////////////////////
  // get trigger information//
  ////////////////////////////

  edm::Handle<TriggerResults> hltResults;
  triggerResultsTag_ = InputTag(hlTriggerResults_,"",hltProcess_);
  iEvent.getByLabel(triggerResultsTag_,hltResults);

  edm::Handle<trigger::TriggerEvent> hltEvent;
  triggerEventTag_ = InputTag("hltTriggerSummaryAOD","",hltProcess_);
  iEvent.getByLabel(triggerEventTag_,hltEvent);

  const TriggerNames & triggerNames = iEvent.triggerNames(*hltResults);
  hlNames = triggerNames.triggerNames();

  triggerStatus   = ULong64_t(0x0);

  for (int i=0; i < (int)hlNames.size(); ++i) {
    if (!triggerDecision(hltResults, i)) continue;

    for (int j = 0; j < (int)triggerPaths_.size(); ++j){
      if (triggerPaths_[j] == "") continue;

      if (hlNames[i].compare(0, triggerPaths_[j].length(),triggerPaths_[j]) == 0) {
        //cout << hlNames[i] << " ?= " << triggerPaths_[j] << endl;
        triggerStatus |= ULong64_t(0x01) << j;
        hltPrescale[j] = 1;

        /* if (isRealData) {
           pair<int, int> preScales;
           preScales = hltConfig_.prescaleValues(iEvent, iSetup, hlNames[i]);
           hltPrescale[j] = preScales.first*preScales.second;
           } */
      }
    }
  }

  for(unsigned int t = 1; t<hlNames.size();t++){
    analyzeTrigger(hltResults, hltEvent, hlNames[t], &trigCount);
  }

  ++nEvents;
  if (!skimLepton_){
    eventTree -> Fill();
  }
  else if(skimLepton_ && (eleCount > 0 || muCount > 0)) // possibly specify a cut in configuration
    eventTree -> Fill();

  beamSpot->Clear();
  primaryVtx    -> Clear("C");
  recoJets      -> Clear("C");
  //recoJPT       -> Clear("C");
  recoMuons     -> Clear("C");
  recoElectrons -> Clear("C");
  recoPhotons   -> Clear("C");
  triggerObjects-> Clear("C");
  genJets       -> Clear("C");
  genParticles  -> Clear("C");
}

// ------------ method called once each job just before starting event loop  ------------
void  ntupleProducer::beginJob()
{
  eventTree      = fs->make<TTree>("eventTree","eventTree");
  jobTree        = fs->make<TTree>("jobTree", "jobTree");

  primaryVtx     = new TClonesArray("TCPrimaryVtx");
  recoJets       = new TClonesArray("TCJet");
  //recoJPT        = new TClonesArray("TCJet");
  recoElectrons  = new TClonesArray("TCElectron");
  recoMuons      = new TClonesArray("TCMuon");
  recoPhotons    = new TClonesArray("TCPhoton");
  triggerObjects = new TClonesArray("TCTriggerObject");
  genJets        = new TClonesArray("TCGenJet");
  genParticles   = new TClonesArray("TCGenParticle");
  beamSpot       = new TVector3();
 
  if(saveMETExtra_){
    pho_up         = new TClonesArray("TLorentzVector");
    pho_down       = new TClonesArray("TLorentzVector");
    pho_mva_up     = new TClonesArray("TLorentzVector");
    pho_mva_down   = new TClonesArray("TLorentzVector");
  }

  jetCon_pat     = new TClonesArray("TLorentzVector");
  jetCon_smear   = new TClonesArray("TLorentzVector");
 
  pfMET.reset(  new TCMET);
  rawMET.reset(  new TCMET);  
  corrMET.reset(  new TCMET);  
  mvaMET.reset(  new TCMET);  

  if(saveMETExtra_){

    jerUpMET.reset(  new TCMET);
    jerDownMET.reset(  new TCMET);
    jerUpMVAMET.reset(  new TCMET);
    jerDownMVAMET.reset(  new TCMET);
    
    phoUpMET.reset(  new TCMET);
    phoDownMET.reset(  new TCMET);
    phoUpMVAMET.reset(  new TCMET);
    phoDownMVAMET.reset(  new TCMET);
    
    jetupMET.reset(  new TCMET);
    jetdownMET.reset(  new TCMET);
    jetupMVAMET.reset(  new TCMET);
    jetdownMVAMET.reset(  new TCMET);

    uncUpMET.reset(  new TCMET);
    uncDownMET.reset(  new TCMET); 
    uncUpMVAMET.reset(  new TCMET); 
    uncDownMVAMET.reset(  new TCMET); 
  }

  h1_numOfEvents = fs->make<TH1F>("numOfEvents", "total number of events, unskimmed", 1,0,1);

  if(saveMETExtra_){

  eventTree->Branch("pho_up", &pho_up, 6400, 0);
  eventTree->Branch("pho_down", &pho_down, 6400, 0);
  eventTree->Branch("pho_mva_up", &pho_mva_up, 6400, 0);  
  eventTree->Branch("pho_mva_down", &pho_mva_down, 6400, 0);
  }

  eventTree->Branch("jetCon_pat", &jetCon_pat, 6400,0);
  eventTree->Branch("jetCon_smear", &jetCon_smear, 6400,0);

  eventTree->Branch("recoJets",     &recoJets,       6400, 0);
  //eventTree->Branch("recoJPT",      &recoJPT,        6400, 0);
  eventTree->Branch("recoElectrons",&recoElectrons,  6400, 0);
  eventTree->Branch("recoMuons",    &recoMuons,      6400, 0);
  eventTree->Branch("recoPhotons",  &recoPhotons,    6400, 0);

  eventTree->Branch("pfMET",         pfMET.get(),     6400, 0);
  eventTree->Branch("rawMET",        rawMET.get(),     6400, 0);
  eventTree->Branch("corrMET",       corrMET.get(),     6400, 0);
  eventTree->Branch("mvaMET",        mvaMET.get(),     6400, 0);

  if(saveMETExtra_){
    eventTree->Branch("jerUpMET", jerUpMET.get(),  6400, 0);
    eventTree->Branch("jerDownMET", jerDownMET.get(),  6400, 0); 
    eventTree->Branch("jerUpMVAMET", jerUpMVAMET.get(),  6400, 0); 
    eventTree->Branch("jerDownMVAMET", jerDownMVAMET.get(),  6400, 0); 
    
    eventTree->Branch("phoUpMET", phoUpMET.get(), 6400, 0);
    eventTree->Branch("phoDownMET", phoDownMET.get(), 6400, 0);
    eventTree->Branch("phoDownMVAMET", phoUpMVAMET.get(), 6400, 0);
    eventTree->Branch("phoDownMVAMET", phoDownMVAMET.get(), 6400, 0);
    
    eventTree->Branch("jetupMET", jetupMET.get(), 6400,0);
    eventTree->Branch("jetdownMET", jetdownMET.get(), 6400,0); 
    eventTree->Branch("jetupMVAMET", jetupMVAMET.get(), 6400,0); 
    eventTree->Branch("jetdownMVAMET", jetdownMVAMET.get(), 6400,0); 
    
    eventTree->Branch("uncUpMET",uncUpMET.get(),6400,0);
    eventTree->Branch("uncDownMET",uncDownMET.get(),6400,0);
    eventTree->Branch("uncUpMVAMET",uncUpMVAMET.get(),6400,0);
    eventTree->Branch("uncDownMVAMET",uncDownMVAMET.get(),6400,0);
  }

  eventTree->Branch("genJets",      &genJets,        6400, 0);
  eventTree->Branch("genParticles", &genParticles,   6400, 0);
  eventTree->Branch("triggerObjects", &triggerObjects, 6400, 0);

  eventTree->Branch("primaryVtx",      &primaryVtx, 6400, 0);
  eventTree->Branch("beamSpot",        &beamSpot,   6400, 0);
  eventTree->Branch("nPUVertices",     &nPUVertices, "nPUVertices/I");
  eventTree->Branch("nPUVerticesTrue", &nPUVerticesTrue, "nPUVerticesTrue/F");

  eventTree->Branch("isRealData", &isRealData,  "isRealData/O");
  eventTree->Branch("runNumber",  &runNumber,   "runNumber/i");
  eventTree->Branch("eventNumber",&eventNumber, "eventNumber/l");
  eventTree->Branch("lumiSection",&lumiSection, "lumiSection/i");
  eventTree->Branch("bunchCross", &bunchCross,  "bunchCross/i");

  eventTree->Branch("ptHat",      &ptHat,       "ptHat/F");
  eventTree->Branch("qScale",     &qScale,      "qScale/F");
  eventTree->Branch("evtWeight",  &evtWeight,   "evtWeight/F");
  eventTree->Branch("rhoFactor",  &rhoFactor,   "rhoFactor/F");
  eventTree->Branch("rho25Factor",&rho25Factor, "rho25Factor/F");
  eventTree->Branch("rhoMuFactor",&rhoMuFactor, "rhoMuFactor/F");
  eventTree->Branch("triggerStatus",&triggerStatus, "triggerStatus/l");
  eventTree->Branch("hltPrescale",hltPrescale, "hltPrescale[64]/i");

  eventTree->Branch("NoiseFilters", &myNoiseFilters.isScraping, "isScraping/O:isNoiseHcalHBHE:isNoiseHcalLaser:isNoiseEcalTP:isNoiseEcalBE:isCSCTightHalo:isCSCLooseHalo:isNoiseTracking:isNoiseEEBadSc:isNoisetrkPOG1:isNoisetrkPOG2:isNoisetrkPOG3");

  jobTree->Branch("nEvents",&nEvents, "nEvents/i");
  jobTree->Branch("triggerNames", "vector<string>", &triggerPaths_);

  // Initialize HLT prescales //
  for (int i = 0; i < (int)(sizeof(hltPrescale)/sizeof(int)); ++i) hltPrescale[i] = 1;

  // Start counting number of events per job //
  nEvents = 0;

  // Photon and Electron PF Iso maker init
  phoIsolator.initializePhotonIsolation(kTRUE);
  phoIsolator.setConeSize(0.3);

  eleIsolator.initializeElectronIsolation(kTRUE);
  eleIsolator.setConeSize(0.4);

  // Initialize Jet PU ID
  myPUJetID.reset(new PileupJetIdAlgo(jetPUIdAlgo_));
}

void ntupleProducer::beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup)
{
  bool changed = true;
  hltConfig_.init(iRun, iSetup, hltProcess_, changed);
  deliveredLumi = 0;
  recordedLumi  = 0;
}

void ntupleProducer::endLuminosityBlock(const edm::LuminosityBlock& iLumi, const edm::EventSetup& iSetup)
{
  //if (isRealData) {
  if (false) {
    edm::Handle<LumiSummary> lumiSummary;
    iLumi.getByLabel("lumiProducer", lumiSummary);

    deliveredLumi  += lumiSummary->avgInsDelLumi()*93.244;
    recordedLumi   += deliveredLumi*lumiSummary->liveFrac();
  }
}


void ntupleProducer::endRun(const edm::Run& iRun, const edm::EventSetup& iSetup)
{
}


void ntupleProducer::endJob()
{
  cout<<nEvents<<endl;
  h1_numOfEvents->SetBinContent(1,nEvents);
  jobTree->Fill();
}


bool ntupleProducer::triggerDecision(edm::Handle<edm::TriggerResults> &hltResults, int iTrigger)
{
  bool triggerPassed = false;
  if(hltResults->wasrun(iTrigger) &&
      hltResults->accept(iTrigger) &&
      !hltResults->error(iTrigger) ){
    triggerPassed = true;
  }
  return triggerPassed;
}


float ntupleProducer::sumPtSquared(const Vertex & v)
{
  float sum = 0.;
  float pT;
  for (Vertex::trackRef_iterator it = v.tracks_begin(); it != v.tracks_end(); it++) {
    pT = (**it).pt();
    float epT=(**it).ptError(); pT=pT>epT ? pT-epT : 0;

    sum += pT*pT;
  }
  return sum;
}


bool ntupleProducer::isFilteredOutScraping( const edm::Event& iEvent, const edm::EventSetup& iSetup, int numtrack, double thresh)
{

  bool  accepted = false;
  float fraction = 0;
  // get GeneralTracks collection

  edm::Handle<reco::TrackCollection> tkRef;
  iEvent.getByLabel("generalTracks",tkRef);
  const reco::TrackCollection* tkColl = tkRef.product();

  int numhighpurity=0;
  reco::TrackBase::TrackQuality _trackQuality = reco::TrackBase::qualityByName("highPurity");

  if(tkColl->size()>(UInt_t)numtrack){
    reco::TrackCollection::const_iterator itk = tkColl->begin();
    reco::TrackCollection::const_iterator itk_e = tkColl->end();
    for(;itk!=itk_e;++itk){
      if(itk->quality(_trackQuality)) numhighpurity++;
    }
    fraction = (float)numhighpurity/(float)tkColl->size();
    if(fraction>thresh) accepted=true;
  } else {
    //if less than 10 Tracks accept the event anyway
    accepted= true;
  }
  return !accepted;  //if filtered out it's not accepted.
}


bool ntupleProducer::associateJetToVertex(reco::PFJet inJet, Handle<reco::VertexCollection> vtxCollection, TCJet *outJet)
{
  if(fabs(inJet.eta()) > 2.5){
    outJet->SetVtxSumPtFrac(-1);
    outJet->SetVtxSumPt(-1);
    outJet->SetVtxTrackFrac(-1);
    outJet->SetVtxNTracks(-1);
    outJet->SetVtxSumPtIndex(0);
    outJet->SetVtxCountIndex(0);

    return false;
  }

  vector<float>  associatedTrackSumPt;
  vector<float>  associatedTrackCount;
  vector<const reco::Track*> jetTracks;
  float sumTrackX, sumTrackY, sumTrackZ, sumTrackPt;
  int   nJetTracks = 0;
  int   vCount = 0;

  sumTrackX = sumTrackY = sumTrackZ  = sumTrackPt = 0;

  //const reco::TrackRefVector &tracks = inJet.associatedTracks();
  const reco::TrackRefVector &tracks = inJet.getTrackRefs();

  for (TrackRefVector::const_iterator iTrack = tracks.begin(); iTrack != tracks.end(); ++iTrack) {
    const reco::Track &jetTrack = **iTrack;

    sumTrackPt += jetTrack.pt();
    sumTrackX  += jetTrack.vx();
    sumTrackY  += jetTrack.vy();
    sumTrackZ  += jetTrack.vz();
    jetTracks.push_back(&jetTrack);
    ++nJetTracks;
  }

  if(jetTracks.size() == 0){
    outJet->SetVtxSumPtFrac(-1);
    outJet->SetVtxSumPt(0);
    outJet->SetVtxTrackFrac(-1);
    outJet->SetVtxNTracks(0);
    outJet->SetVtxSumPtIndex(0);
    outJet->SetVtxCountIndex(0);
    outJet->SetVtx(0., 0., 0.);
  } else {
    outJet->SetVtx(sumTrackX/nJetTracks, sumTrackY/nJetTracks, sumTrackZ/nJetTracks);
    
    for (VertexCollection::const_iterator iVtx = vtxCollection->begin(); iVtx!= vtxCollection->end(); ++iVtx) {
      reco::Vertex myVtx = reco::Vertex(*iVtx);
      if(!myVtx.isValid() || myVtx.isFake()) continue;
      associatedTrackSumPt.push_back(0);
      associatedTrackCount.push_back(0);

      for(Vertex::trackRef_iterator iTrackRef = myVtx.tracks_begin(); iTrackRef != myVtx.tracks_end(); ++iTrackRef){
        const edm::RefToBase<reco::Track> &myTrackRef = *iTrackRef;

        if(myTrackRef.isAvailable()){
          const reco::Track &myVertexTrack = *myTrackRef.get();

          for(vector<const reco::Track*>::const_iterator iTrack = jetTracks.begin(); iTrack != jetTracks.end(); ++iTrack){
            if (*iTrack == &myVertexTrack) {
              associatedTrackSumPt.at(vCount) += myVertexTrack.pt()/sumTrackPt;
              associatedTrackCount.at(vCount) += 1/nJetTracks;
            }
          }
        }
      }
      ++vCount;
    }

    float maxSumPtFraction = 0; float maxCountFraction = 0;
    int   vtxSumPtIndex = 0; int vtxCountIndex = 0;
    int count = 0;

    for (int i = 0; i < vCount; ++i) {
      if (associatedTrackSumPt.at(i) > maxSumPtFraction) {
        maxSumPtFraction = associatedTrackSumPt.at(i);
        vtxSumPtIndex = count + 1;
      }
      if (associatedTrackCount.at(i) > maxCountFraction) {
        maxCountFraction = associatedTrackCount.at(i);
        vtxCountIndex = count + 1;
      }
      ++count;
    }
    outJet->SetVtxSumPtFrac(maxSumPtFraction);
    outJet->SetVtxSumPt(sumTrackPt);
    outJet->SetVtxTrackFrac(maxCountFraction);
    outJet->SetVtxNTracks(nJetTracks);
    outJet->SetVtxSumPtIndex(vtxSumPtIndex);
    outJet->SetVtxCountIndex(vtxCountIndex);
  }

  return true;
}


void ntupleProducer::electronMVA(const reco::GsfElectron* iElectron, TCElectron* eleCon,
    const edm::Event& iEvent, const edm::EventSetup& iSetup, const reco::PFCandidateCollection& PFCandidates, float Rho)
{
  if (verboseMVAs) cout<<"loading up electron MVA values"<<endl;
  //**********************************************************
  //ID variables
  //**********************************************************
  bool validKF= false;
  reco::TrackRef myTrackRef = iElectron->closestCtfTrackRef();
  validKF = (myTrackRef.isAvailable());
  validKF = (myTrackRef.isNonnull());

  eleCon->SetIdMap("fbrem", (iElectron->fbrem() < -1) ? -1 : iElectron->fbrem());
  eleCon->SetIdMap("gsfChi2", (iElectron->gsfTrack()->normalizedChi2() > 200) ? 200 : iElectron->gsfTrack()->normalizedChi2());
  eleCon->SetIdMap("kfChi2", (validKF) ? myTrackRef->normalizedChi2() : 0);
  eleCon->SetIdMap("kfNLayers", (validKF) ? myTrackRef->hitPattern().trackerLayersWithMeasurement() : -1.);
  eleCon->SetIdMap("kfNLayersAll", (validKF) ? myTrackRef->numberOfValidHits() : -1.);
  eleCon->SetIdMap("dEta", (fabs(iElectron->deltaEtaSuperClusterTrackAtVtx()) > 0.06) ? 0.06 : fabs(iElectron->deltaEtaSuperClusterTrackAtVtx()));
  eleCon->SetIdMap("dPhi", iElectron->deltaPhiSuperClusterTrackAtVtx());
  eleCon->SetIdMap("dEtaAtCalo", iElectron->deltaEtaSeedClusterTrackAtCalo());
  eleCon->SetIdMap("SigmaIPhiIPhi", (isnan(iElectron->sigmaIphiIphi())) ? 0 : iElectron->sigmaIphiIphi() );
  eleCon->SetIdMap("SCEtaWidth", iElectron->superCluster()->etaWidth());
  eleCon->SetIdMap("SCPhiWidth", iElectron->superCluster()->phiWidth());
  eleCon->SetIdMap("EoP", (iElectron->eSuperClusterOverP() > 20) ? 20 : iElectron->eSuperClusterOverP());
  eleCon->SetIdMap("ome1x5oe5x5",(iElectron->e5x5()) !=0. ? 1.-(iElectron->e1x5()/iElectron->e5x5()) : -1.);
  if (eleCon->IdMap("ome1x5oe5x5") < -1) eleCon->SetIdMap("ome1x5oe5x5",-1);
  if (eleCon->IdMap("ome1x5oe5x5") > 2) eleCon->SetIdMap("ome1x5oe5x5",2);
  eleCon->SetIdMap("R9",(iElectron->r9() > 5) ? 5 : iElectron->r9());
  eleCon->SetIdMap("ooemoopV1",(1.0/iElectron->ecalEnergy()) - (1.0 / iElectron->p()));
  eleCon->SetIdMap("ooemoopV2",(1.0/iElectron->superCluster()->energy()) - (1.0 / iElectron->trackMomentumAtVtx().R()));
  eleCon->SetIdMap("eopOut",(iElectron->eEleClusterOverPout() > 20) ? 20 : iElectron->eEleClusterOverPout());
  eleCon->SetIdMap("preShowerORaw",iElectron->superCluster()->preshowerEnergy() / iElectron->superCluster()->rawEnergy());

  InputTag  vertexLabel(string("offlinePrimaryVertices"));
  Handle<reco::VertexCollection> thePrimaryVertexColl;
  iEvent.getByLabel(vertexLabel,thePrimaryVertexColl);

  Vertex dummy;
  const Vertex *pv = &dummy;
  if (thePrimaryVertexColl->size() != 0) {
    pv = &*thePrimaryVertexColl->begin();
  } else { // create a dummy PV
    Vertex::Error e;
    e(0, 0) = 0.0015 * 0.0015;
    e(1, 1) = 0.0015 * 0.0015;
    e(2, 2) = 15. * 15.;
    Vertex::Point p(0, 0, 0);
    dummy = Vertex(p, e, 0, 0, 0);
  }

  //d0
  float fMVAVar_d0 = -9999.0;
  if (iElectron->gsfTrack().isNonnull()) {
    fMVAVar_d0 = (-1.0)*iElectron->gsfTrack()->dxy(pv->position());
  } else if (iElectron->closestCtfTrackRef().isNonnull()) {
    fMVAVar_d0 = (-1.0)*iElectron->closestCtfTrackRef()->dxy(pv->position());
  } else {
    fMVAVar_d0 = -9999.0;
  }

  eleCon->SetIdMap("d0",fMVAVar_d0);

  /*
    This is added into the main part 

    //default values for IP3D
    float fMVAVar_ip3d      = -999.0;
    float fMVAVar_ip3dSig   = 0.0;


    edm::ESHandle<TransientTrackBuilder> builder;
    iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", builder);
    TransientTrackBuilder thebuilder = *(builder.product());

    if (iElectron->gsfTrack().isNonnull()) {
    const double gsfsign   = ( (-iElectron->gsfTrack()->dxy(pv->position()))   >=0 ) ? 1. : -1.;

    const reco::TransientTrack &tt = thebuilder.build(iElectron->gsfTrack());
    const std::pair<bool,Measurement1D> &ip3dpv =  IPTools::absoluteImpactParameter3D(tt,*pv);
    if (ip3dpv.first) {
    double ip3d = gsfsign*ip3dpv.second.value();
    double ip3derr = ip3dpv.second.error();
    fMVAVar_ip3d = ip3d;
    fMVAVar_ip3dSig = ip3d/ip3derr;
    }
    }

  eleCon->SetIdMap("ip3d",fMVAVar_ip3d);
  eleCon->SetIdMap("ip3dSig",fMVAVar_ip3dSig);

  ---- up to this ---------
  */

  //**********************************************************
  //Isolation variables
  //**********************************************************

  Double_t tmpChargedIso_DR0p0To0p1  = 0;
  Double_t tmpChargedIso_DR0p1To0p2  = 0;
  Double_t tmpChargedIso_DR0p2To0p3  = 0;
  Double_t tmpChargedIso_DR0p3To0p4  = 0;
  Double_t tmpChargedIso_DR0p4To0p5  = 0;
  Double_t tmpGammaIso_DR0p0To0p1  = 0;
  Double_t tmpGammaIso_DR0p1To0p2  = 0;
  Double_t tmpGammaIso_DR0p2To0p3  = 0;
  Double_t tmpGammaIso_DR0p3To0p4  = 0;
  Double_t tmpGammaIso_DR0p4To0p5  = 0;
  Double_t tmpNeutralHadronIso_DR0p0To0p1  = 0;
  Double_t tmpNeutralHadronIso_DR0p1To0p2  = 0;
  Double_t tmpNeutralHadronIso_DR0p2To0p3  = 0;
  Double_t tmpNeutralHadronIso_DR0p3To0p4  = 0;
  Double_t tmpNeutralHadronIso_DR0p4To0p5  = 0;

  //************************************************************
  //Note: Input collection is assumed to be PFNoPU collection
  //************************************************************
  for (reco::PFCandidateCollection::const_iterator iP = PFCandidates.begin();
      iP != PFCandidates.end(); ++iP) {

    double dr = sqrt(pow(iP->eta() - iElectron->eta(),2) + pow(acos(cos(iP->phi() - iElectron->phi())),2));

    Bool_t passVeto = kTRUE;
    //Charged
    if(iP->trackRef().isNonnull()) {

      //make sure charged pf candidates pass the PFNoPU condition (assumed)

      //************************************************************
      // Veto any PFmuon, or PFEle
      if (iP->particleId() == reco::PFCandidate::e || iP->particleId() == reco::PFCandidate::mu) passVeto = kFALSE;
      //************************************************************
      //************************************************************
      // Footprint Veto
      if (fabs(iElectron->superCluster()->eta()) > 1.479 && dr < 0.015) passVeto = kFALSE;
      if (iP->superClusterRef().isNonnull() &&
          iP->superClusterRef() == iElectron->superCluster()) passVeto = kFALSE;
      if (iP->gsfTrackRef().isNonnull() && iElectron->gsfTrack().isNonnull() &&
          iP->gsfTrackRef() == iElectron->gsfTrack()) passVeto = kFALSE;
      if (iP->trackRef().isNonnull() && iElectron->closestCtfTrackRef().isNonnull() &&
          iP->trackRef() == iElectron->closestCtfTrackRef()) passVeto = kFALSE;
      //************************************************************
      if (passVeto) {
        if (dr < 0.1) tmpChargedIso_DR0p0To0p1 += iP->pt();
        if (dr >= 0.1 && dr < 0.2) tmpChargedIso_DR0p1To0p2 += iP->pt();
        if (dr >= 0.2 && dr < 0.3) tmpChargedIso_DR0p2To0p3 += iP->pt();
        if (dr >= 0.3 && dr < 0.4) tmpChargedIso_DR0p3To0p4 += iP->pt();
        if (dr >= 0.4 && dr < 0.5) tmpChargedIso_DR0p4To0p5 += iP->pt();
      } //pass veto
    }
    //Gamma
    else if (iP->particleId() == reco::PFCandidate::gamma) {
      //************************************************************
      // Footprint Veto
      if (fabs(iElectron->superCluster()->eta()) > 1.479 && dr < 0.08) passVeto = kFALSE;
      if (iP->superClusterRef() == iElectron->superCluster()) passVeto = kFALSE;
      //************************************************************
      if (passVeto) {
        if (dr < 0.1) tmpGammaIso_DR0p0To0p1 += iP->pt();
        if (dr >= 0.1 && dr < 0.2) tmpGammaIso_DR0p1To0p2 += iP->pt();
        if (dr >= 0.2 && dr < 0.3) tmpGammaIso_DR0p2To0p3 += iP->pt();
        if (dr >= 0.3 && dr < 0.4) tmpGammaIso_DR0p3To0p4 += iP->pt();
        if (dr >= 0.4 && dr < 0.5) tmpGammaIso_DR0p4To0p5 += iP->pt();
      }
    }
    //NeutralHadron
    else {
      if (dr < 0.1) tmpNeutralHadronIso_DR0p0To0p1 += iP->pt();
      if (dr >= 0.1 && dr < 0.2) tmpNeutralHadronIso_DR0p1To0p2 += iP->pt();
      if (dr >= 0.2 && dr < 0.3) tmpNeutralHadronIso_DR0p2To0p3 += iP->pt();
      if (dr >= 0.3 && dr < 0.4) tmpNeutralHadronIso_DR0p3To0p4 += iP->pt();
      if (dr >= 0.4 && dr < 0.5) tmpNeutralHadronIso_DR0p4To0p5 += iP->pt();
    }
  } //loop over PF candidates

  double fMVAVar_ChargedIso_DR0p0To0p1,fMVAVar_ChargedIso_DR0p1To0p2,fMVAVar_ChargedIso_DR0p2To0p3,fMVAVar_ChargedIso_DR0p3To0p4,fMVAVar_ChargedIso_DR0p4To0p5,fMVAVar_GammaIso_DR0p0To0p1,fMVAVar_GammaIso_DR0p1To0p2,fMVAVar_GammaIso_DR0p2To0p3,
         fMVAVar_GammaIso_DR0p3To0p4,fMVAVar_GammaIso_DR0p4To0p5,fMVAVar_NeutralHadronIso_DR0p0To0p1,fMVAVar_NeutralHadronIso_DR0p1To0p2,fMVAVar_NeutralHadronIso_DR0p2To0p3,fMVAVar_NeutralHadronIso_DR0p3To0p4,fMVAVar_NeutralHadronIso_DR0p4To0p5;

  bool doPUCorrection = false;
  if (doPUCorrection) {
    fMVAVar_ChargedIso_DR0p0To0p1   = TMath::Min((tmpChargedIso_DR0p0To0p1)/iElectron->pt(), 2.5);
    fMVAVar_ChargedIso_DR0p1To0p2   = TMath::Min((tmpChargedIso_DR0p1To0p2)/iElectron->pt(), 2.5);
    fMVAVar_ChargedIso_DR0p2To0p3   = TMath::Min((tmpChargedIso_DR0p2To0p3)/iElectron->pt(), 2.5);
    fMVAVar_ChargedIso_DR0p3To0p4   = TMath::Min((tmpChargedIso_DR0p3To0p4)/iElectron->pt(), 2.5);
    fMVAVar_ChargedIso_DR0p4To0p5   = TMath::Min((tmpChargedIso_DR0p4To0p5)/iElectron->pt(), 2.5);
    fMVAVar_GammaIso_DR0p0To0p1     = TMath::Max(TMath::Min((tmpGammaIso_DR0p0To0p1 -
            Rho*ElectronEffectiveArea::GetElectronEffectiveArea(ElectronEffectiveArea::kEleGammaIsoDR0p0To0p1, iElectron->superCluster()->eta(), ElectronEffectiveArea::kEleEAData2012))/iElectron->pt(), 2.5), 0.0);
    fMVAVar_GammaIso_DR0p1To0p2     = TMath::Max(TMath::Min((tmpGammaIso_DR0p1To0p2 -
            Rho*ElectronEffectiveArea::GetElectronEffectiveArea(ElectronEffectiveArea::kEleGammaIsoDR0p1To0p2, iElectron->superCluster()->eta(), ElectronEffectiveArea::kEleEAData2012))/iElectron->pt(), 2.5), 0.0);
    fMVAVar_GammaIso_DR0p2To0p3 = TMath::Max(TMath::Min((tmpGammaIso_DR0p2To0p3 -
            Rho*ElectronEffectiveArea::GetElectronEffectiveArea(ElectronEffectiveArea::kEleGammaIsoDR0p2To0p3, iElectron->superCluster()->eta(), ElectronEffectiveArea::kEleEAData2012))/iElectron->pt(), 2.5), 0.0);
    fMVAVar_GammaIso_DR0p3To0p4 = TMath::Max(TMath::Min((tmpGammaIso_DR0p3To0p4 -
            Rho*ElectronEffectiveArea::GetElectronEffectiveArea(ElectronEffectiveArea::kEleGammaIsoDR0p3To0p4, iElectron->superCluster()->eta(), ElectronEffectiveArea::kEleEAData2012))/iElectron->pt(), 2.5), 0.0);
    fMVAVar_GammaIso_DR0p4To0p5 = TMath::Max(TMath::Min((tmpGammaIso_DR0p4To0p5 -
            Rho*ElectronEffectiveArea::GetElectronEffectiveArea(ElectronEffectiveArea::kEleGammaIsoDR0p4To0p5, iElectron->superCluster()->eta(), ElectronEffectiveArea::kEleEAData2012))/iElectron->pt(), 2.5), 0.0);
    fMVAVar_NeutralHadronIso_DR0p0To0p1 = TMath::Max(TMath::Min((tmpNeutralHadronIso_DR0p0To0p1 -
            Rho*ElectronEffectiveArea::GetElectronEffectiveArea(ElectronEffectiveArea::kEleNeutralHadronIsoDR0p0To0p1, iElectron->superCluster()->eta(), ElectronEffectiveArea::kEleEAData2012))/iElectron->pt(), 2.5), 0.0);
    fMVAVar_NeutralHadronIso_DR0p1To0p2 = TMath::Max(TMath::Min((tmpNeutralHadronIso_DR0p1To0p2 -
            Rho*ElectronEffectiveArea::GetElectronEffectiveArea(ElectronEffectiveArea::kEleNeutralHadronIsoDR0p1To0p2, iElectron->superCluster()->eta(), ElectronEffectiveArea::kEleEAData2012))/iElectron->pt(), 2.5), 0.0);
    fMVAVar_NeutralHadronIso_DR0p2To0p3 = TMath::Max(TMath::Min((tmpNeutralHadronIso_DR0p2To0p3 -
            Rho*ElectronEffectiveArea::GetElectronEffectiveArea(ElectronEffectiveArea::kEleNeutralHadronIsoDR0p2To0p3, iElectron->superCluster()->eta(), ElectronEffectiveArea::kEleEAData2012))/iElectron->pt(), 2.5), 0.0);
    fMVAVar_NeutralHadronIso_DR0p3To0p4 = TMath::Max(TMath::Min((tmpNeutralHadronIso_DR0p3To0p4 -
            Rho*ElectronEffectiveArea::GetElectronEffectiveArea(ElectronEffectiveArea::kEleNeutralHadronIsoDR0p3To0p4, iElectron->superCluster()->eta(), ElectronEffectiveArea::kEleEAData2012))/iElectron->pt(), 2.5), 0.0);
    fMVAVar_NeutralHadronIso_DR0p4To0p5 = TMath::Max(TMath::Min((tmpNeutralHadronIso_DR0p4To0p5 -
            Rho*ElectronEffectiveArea::GetElectronEffectiveArea(ElectronEffectiveArea::kEleNeutralHadronIsoDR0p4To0p5, iElectron->superCluster()->eta(), ElectronEffectiveArea::kEleEAData2012))/iElectron->pt(), 2.5), 0.0);
  } else {
    fMVAVar_ChargedIso_DR0p0To0p1   = TMath::Min((tmpChargedIso_DR0p0To0p1)/iElectron->pt(), 2.5);
    fMVAVar_ChargedIso_DR0p1To0p2   = TMath::Min((tmpChargedIso_DR0p1To0p2)/iElectron->pt(), 2.5) / 0.03;
    fMVAVar_ChargedIso_DR0p2To0p3 = TMath::Min((tmpChargedIso_DR0p2To0p3)/iElectron->pt(), 2.5) / 0.05;
    fMVAVar_ChargedIso_DR0p3To0p4 = TMath::Min((tmpChargedIso_DR0p3To0p4)/iElectron->pt(), 2.5) / 0.07;
    fMVAVar_ChargedIso_DR0p4To0p5 = TMath::Min((tmpChargedIso_DR0p4To0p5)/iElectron->pt(), 2.5) / 0.09;
    fMVAVar_GammaIso_DR0p0To0p1 = TMath::Max(TMath::Min((tmpGammaIso_DR0p0To0p1)/iElectron->pt(), 2.5), 0.0);
    fMVAVar_GammaIso_DR0p1To0p2 = TMath::Max(TMath::Min((tmpGammaIso_DR0p1To0p2)/iElectron->pt(), 2.5), 0.0) / 0.03;
    fMVAVar_GammaIso_DR0p2To0p3 = TMath::Max(TMath::Min((tmpGammaIso_DR0p2To0p3)/iElectron->pt(), 2.5), 0.0) / 0.05;
    fMVAVar_GammaIso_DR0p3To0p4 = TMath::Max(TMath::Min((tmpGammaIso_DR0p3To0p4)/iElectron->pt(), 2.5), 0.0) / 0.07;
    fMVAVar_GammaIso_DR0p4To0p5 = TMath::Max(TMath::Min((tmpGammaIso_DR0p4To0p5)/iElectron->pt(), 2.5), 0.0) / 0.09;
    fMVAVar_NeutralHadronIso_DR0p0To0p1 = TMath::Max(TMath::Min((tmpNeutralHadronIso_DR0p0To0p1)/iElectron->pt(), 2.5), 0.0);
    fMVAVar_NeutralHadronIso_DR0p1To0p2 = TMath::Max(TMath::Min((tmpNeutralHadronIso_DR0p1To0p2)/iElectron->pt(), 2.5), 0.0) / 0.03;
    fMVAVar_NeutralHadronIso_DR0p2To0p3 = TMath::Max(TMath::Min((tmpNeutralHadronIso_DR0p2To0p3)/iElectron->pt(), 2.5), 0.0) / 0.05;
    fMVAVar_NeutralHadronIso_DR0p3To0p4 = TMath::Max(TMath::Min((tmpNeutralHadronIso_DR0p3To0p4)/iElectron->pt(), 2.5), 0.0) / 0.07;
    fMVAVar_NeutralHadronIso_DR0p4To0p5 = TMath::Max(TMath::Min((tmpNeutralHadronIso_DR0p4To0p5)/iElectron->pt(), 2.5), 0.0) / 0.09;
  }

  eleCon->SetIsoMap("ChargedIso_DR0p0To0p1",fMVAVar_ChargedIso_DR0p0To0p1);
  eleCon->SetIsoMap("ChargedIso_DR0p1To0p2",fMVAVar_ChargedIso_DR0p1To0p2);
  eleCon->SetIsoMap("ChargedIso_DR0p2To0p3",fMVAVar_ChargedIso_DR0p2To0p3);
  eleCon->SetIsoMap("ChargedIso_DR0p3To0p4",fMVAVar_ChargedIso_DR0p3To0p4);
  eleCon->SetIsoMap("ChargedIso_DR0p4To0p5",fMVAVar_ChargedIso_DR0p4To0p5);
  eleCon->SetIsoMap("GammaIso_DR0p0To0p1",fMVAVar_GammaIso_DR0p0To0p1);
  eleCon->SetIsoMap("GammaIso_DR0p1To0p2",fMVAVar_GammaIso_DR0p1To0p2);
  eleCon->SetIsoMap("GammaIso_DR0p2To0p3",fMVAVar_GammaIso_DR0p2To0p3);
  eleCon->SetIsoMap("GammaIso_DR0p3To0p4",fMVAVar_GammaIso_DR0p3To0p4);
  eleCon->SetIsoMap("GammaIso_DR0p4To0p5",fMVAVar_GammaIso_DR0p4To0p5);
  eleCon->SetIsoMap("NeutralHadronIso_DR0p0To0p1",fMVAVar_NeutralHadronIso_DR0p0To0p1);
  eleCon->SetIsoMap("NeutralHadronIso_DR0p1To0p2",fMVAVar_NeutralHadronIso_DR0p1To0p2);
  eleCon->SetIsoMap("NeutralHadronIso_DR0p2To0p3",fMVAVar_NeutralHadronIso_DR0p2To0p3);
  eleCon->SetIsoMap("NeutralHadronIso_DR0p3To0p4",fMVAVar_NeutralHadronIso_DR0p3To0p4);
  eleCon->SetIsoMap("NeutralHadronIso_DR0p4To0p5",fMVAVar_NeutralHadronIso_DR0p4To0p5);

  bool preSelPassV2 = false;

  double electronTrackZ = 0;
  if (iElectron->gsfTrack().isNonnull()) {
    electronTrackZ = iElectron->gsfTrack()->dz(pv->position());
  } else if (iElectron->closestCtfTrackRef().isNonnull()) {
    electronTrackZ = iElectron->closestCtfTrackRef()->dz(pv->position());
  }

  if (fabs(iElectron->superCluster()->eta()) < 1.479) {
    if ( (
          iElectron->sigmaIetaIeta()< 0.01
          && fabs(iElectron->deltaEtaSuperClusterTrackAtVtx()) < 0.007
          && fabs(iElectron->deltaPhiSuperClusterTrackAtVtx()) < 0.15
          && iElectron->hadronicOverEm() < 0.12
          && iElectron->gsfTrack()->trackerExpectedHitsInner().numberOfLostHits() <= 1
          && fabs(electronTrackZ) < 0.1
          && ( iElectron->dr03TkSumPt()) / iElectron->pt() < 0.2
          && ( fmax(iElectron->dr03EcalRecHitSumEt() - 1.0,0.0) ) /iElectron->pt() < 0.20
          && (iElectron->dr03HcalTowerSumEt()) / iElectron->pt() < 0.20
         )
       ) {
      preSelPassV2= true;
    }
  } else { //endcap
    if ( (
          iElectron->sigmaIetaIeta()< 0.03
          && fabs(iElectron->deltaEtaSuperClusterTrackAtVtx()) < 0.009
          && fabs(iElectron->deltaPhiSuperClusterTrackAtVtx()) < 0.10
          && iElectron->hadronicOverEm() < 0.10
          && iElectron->gsfTrack()->trackerExpectedHitsInner().numberOfLostHits() <= 1
          && fabs(electronTrackZ) < 0.1
          && (iElectron->dr03TkSumPt() ) / iElectron->pt() < 0.2
          && (iElectron->dr03EcalRecHitSumEt() ) / iElectron->pt() < 0.20
          && (iElectron->dr03HcalTowerSumEt()) / iElectron->pt() < 0.20

         )
       ) {
      preSelPassV2 = true;
    }
  }
  eleCon->SetIdMap("preSelPassV2",preSelPassV2);

  bool preSelPassV1 = false;

  if (fabs(iElectron->superCluster()->eta()) < 1.479) {
    if ( (
          iElectron->sigmaIetaIeta()< 0.014
          && iElectron->hadronicOverEm() < 0.15
          && iElectron->gsfTrack()->trackerExpectedHitsInner().numberOfLostHits() == 0
          && ( iElectron->dr03TkSumPt()) / iElectron->pt() < 0.2
          && ( iElectron->dr03EcalRecHitSumEt()) /iElectron->pt() < 0.2
          && (iElectron->dr03HcalTowerSumEt()) / iElectron->pt() < 0.2
         )
       ) {
      preSelPassV1 = true;
    }
  } else { //endcap
    if ( (
          iElectron->sigmaIetaIeta()< 0.035
          && iElectron->hadronicOverEm() < 0.10
          && iElectron->gsfTrack()->trackerExpectedHitsInner().numberOfLostHits() == 0
          && (iElectron->dr03TkSumPt() ) / iElectron->pt() < 0.2
          && (iElectron->dr03EcalRecHitSumEt() ) / iElectron->pt() < 0.20
          && (iElectron->dr03HcalTowerSumEt()) / iElectron->pt() < 0.20

         )
       ) {
      preSelPassV1 = true;
    }
  }

  eleCon->SetIdMap("preSelPassV1",preSelPassV1);

  bool preSelPassV3 = false;
  if (fabs(iElectron->superCluster()->eta()) < 1.479) {
    if ( (
          iElectron->sigmaIetaIeta()< 0.014
          && iElectron->hadronicOverEm() < 0.15
          && iElectron->gsfTrack()->trackerExpectedHitsInner().numberOfLostHits() == 0
          && ( iElectron->dr03TkSumPt()) / iElectron->pt() < 0.2
          && ( iElectron->dr03EcalRecHitSumEt()) /iElectron->pt() < 0.2
          && (iElectron->dr03HcalTowerSumEt()) / iElectron->pt() < 0.2
         )
       ) {
      preSelPassV3 = true;
    }
  } else { //endcap
    if ( (
          iElectron->sigmaIetaIeta()< 0.035
          && iElectron->hadronicOverEm() < 0.10
          && iElectron->gsfTrack()->trackerExpectedHitsInner().numberOfLostHits() == 0
          && (iElectron->dr03TkSumPt() ) / iElectron->pt() < 0.2
          && (iElectron->dr03EcalRecHitSumEt()-1.0 ) / iElectron->pt() < 0.20
          && (iElectron->dr03HcalTowerSumEt()) / iElectron->pt() < 0.20

         )
       ) {
      preSelPassV3 = true;
    }
  }
  eleCon->SetIdMap("preSelPassV3",preSelPassV3);

  return;
}

void ntupleProducer::analyzeTrigger(edm::Handle<edm::TriggerResults> &hltResults, edm::Handle<trigger::TriggerEvent> &hltEvent, const std::string& triggerName, int* trigCount) {

  using namespace trigger;

  const unsigned int n(hltConfig_.size());
  const unsigned int triggerIndex(hltConfig_.triggerIndex(triggerName));

  TLorentzVector triggerLepton;

  // abort on invalid trigger name

  bool goodTrigger = false;
  for (unsigned int i =0; i< triggerPaths_.size(); i++){
    if (triggerName.find(triggerPaths_[i]) != string::npos){
      goodTrigger = true;
      break;
    }
  }

  if(!goodTrigger) return;

  //if(verboseTrigs){
  //  std::cout<<" n = "<<n<<" triggerIndex = "<<triggerIndex<<" size = "<<hltConfig_.size()<<std::endl;
  //  std::cout<<" Analyze triggerName : "<<triggerName<<std::endl;
  //}
  //std::cout<<" Analyze triggerName : "<<triggerName<<std::endl;

  if (triggerIndex>=n) {
    if(verboseTrigs){
      cout << "DimuonAna::analyzeTrigger: path "
        << triggerName << " - not found!" << endl;
    }
    return;
  }

  // modules on this trigger path
  // const unsigned int moduleIndex(hltResults->index(triggerIndex));
  const unsigned int moduleIndex(hltResults->index(triggerIndex));
  const unsigned int m(hltConfig_.size(triggerIndex));
  const vector<string>& moduleLabels(hltConfig_.moduleLabels(triggerIndex));
  if (moduleIndex != m-1) return;
  if(verboseTrigs){
    cout << "DimuonAna::analyzeTrigger: path "
      << triggerName << " [" << triggerIndex << "]" << endl;

    std::cout<<"  n = "<< n<<" triggerIndex = "<<triggerIndex<<" m = "<<m<<std::endl;
    std::cout<<" moduleLabels = "<<moduleLabels.size()<<" moduleIndex = "<<moduleIndex<<std::endl;

    // Results from TriggerResults product
    cout << " Trigger path status:"
      << " WasRun=" << hltResults->wasrun(triggerIndex)
      << " Accept=" << hltResults->accept(triggerIndex)
      << " Error =" << hltResults->error(triggerIndex)
      << endl;
    cout << " Last active module - label/type: "
      << moduleLabels[moduleIndex] << "/" << hltConfig_.moduleType(moduleLabels[moduleIndex])
      << " [" << moduleIndex << " out of 0-" << (m-1) << " on this path]"
      << endl;
  }
  assert (moduleIndex<m);

  // Results from TriggerEvent product - Attention: must look only for
  // modules actually run in this path for this event!
  std::vector < GlobalVector > passMomenta;
  for (unsigned int j=0; j<=moduleIndex; ++j) {
    const string& moduleLabel(moduleLabels[j]);
    const string  moduleType(hltConfig_.moduleType(moduleLabel));

    // check whether the module is packed up in TriggerEvent product
    //cout<<hltEvent->filterIndex(InputTag(moduleLabel,"",hltProcess_))<<endl;

    const unsigned int filterIndex(hltEvent->filterIndex(InputTag(moduleLabel,"",hltProcess_)));

    //  if ( (moduleLabel.find("Calo") == string::npos) )continue;
    //  if ( (moduleLabel.find("hltEventle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoDZ") == string::npos)
    //      && (moduleLabel.find("hltEventle17CaloId") == string::npos)
    //      && (moduleLabel.find("hltEventle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoDoubleFilter") == string::npos) ) continue;

    if(verboseTrigs){
      std::cout<<" j = "<<j<<" modLabel/moduleType = "<<moduleLabel<<"/"<<moduleType<<" filterIndex = "<<filterIndex<<" sizeF = "<<hltEvent->sizeFilters()<<std::endl;
    }
    if (filterIndex<hltEvent->sizeFilters()) {
      if(verboseTrigs){
        cout << " 'L3' (or 'L1', 'L2') filter in slot " << j << " - label/type " << moduleLabel << "/" << moduleType << endl;
      }
      const Vids& VIDS (hltEvent->filterIds(filterIndex));
      const Keys& KEYS(hltEvent->filterKeys(filterIndex));
      const size_type nI(VIDS.size());
      const size_type nK(KEYS.size());
      assert(nI==nK);
      const size_type n(max(nI,nK));
      if(verboseTrigs){
        cout << "   " << n  << " accepted 'L3' (or 'L1', 'L2') objects found: " << endl;
      }
      const TriggerObjectCollection& TOC(hltEvent->getObjects());
      for (size_type i=0; i!=n; ++i) {
        if(0==i){
          passMomenta.clear();
        }
        const TriggerObject& TO(TOC[KEYS[i]]);
        GlobalVector momentumT0(TO.px(),TO.py(),TO.pz());
        if (TO.pt() < 10) continue;
        TCTriggerObject* trigObj = new ((*triggerObjects)[*trigCount]) TCTriggerObject;

        //std::cout<<" i_KEY = "<<i<<" id = "<<TO.id()<<" typ = "<<moduleType<<std::endl;
        //if("HLTLevel1GTSeed"==moduleType){
        //allMuL1TriggerVectors.push_back(momentumT0);
        ////std::cout<<" L1 object found"<<std::endl;
        //}

        if(verboseTrigs){
          std::cout<<" i = "<<i<<" moduleLabel/moduleType : "<<moduleLabel<<"/"<<moduleType<<std::endl;
          cout << "   " << i << " " << VIDS[i] << "/" << KEYS[i] << ": "
            << TO.id() << " " << TO.pt() << " " << TO.eta() << " " << TO.phi() << " " << TO.mass()
            << endl;
        }

        trigObj->SetPtEtaPhiE(TO.pt(),TO.eta(),TO.phi(),TO.energy());
        trigObj->SetHLTName(triggerName);
        trigObj->SetModuleName(moduleLabel);
        trigObj->SetId(TO.id());

        (*trigCount)++;

      }
      //
    }
  }
  //cout<<endl;
  return;
}

float ntupleProducer::MatchBTagsToJets(const reco::JetTagCollection bTags,const reco::PFJet jet)
{
  float discrValue = -999;
  for (tag_iter iTag = bTags.begin(); iTag != bTags.end(); iTag++) {
    if (sqrt(pow(iTag->first->eta() - jet.eta(), 2) + pow(deltaPhi(iTag->first->phi(),jet.phi()), 2)) == 0.) {
      discrValue = iTag->second;
      break;
    }
  }

  return discrValue;
}

void ntupleProducer::initJetEnergyCorrector(const edm::EventSetup &iSetup, bool isData)
{
  //jet energy correction levels to apply on raw jet
  std::vector<JetCorrectorParameters> jetCorPars_;
  std::vector<std::string> jecLevels;
  jecLevels.push_back("L1FastJet");
  jecLevels.push_back("L2Relative");
  jecLevels.push_back("L3Absolute");
  if(isData) jecLevels.push_back("L2L3Residual");

  //check the corrector parameters needed according to the correction levels
  edm::ESHandle<JetCorrectorParametersCollection> parameters;
  iSetup.get<JetCorrectionsRecord>().get(jecTag_,parameters);
  for(std::vector<std::string>::const_iterator ll = jecLevels.begin(); ll != jecLevels.end(); ++ll)
  {
    const JetCorrectorParameters& ip = (*parameters)[*ll];
    jetCorPars_.push_back(ip);
  }

  //instantiate the jet corrector
  jecCor.reset(new FactorizedJetCorrector(jetCorPars_));
}


TCGenParticle* ntupleProducer::addGenParticle(const reco::GenParticle* myParticle, int& genPartCount, map<const reco::GenParticle*,TCGenParticle*>& genMap)
{
  TCGenParticle* genCon;
  map<const reco::GenParticle*,TCGenParticle*>::iterator it;
  it = genMap.find(myParticle);
  if (it == genMap.end()){
    genCon = new ((*genParticles)[genPartCount]) TCGenParticle;
    ++genPartCount;
    genMap[myParticle] = genCon;
    genCon->SetPxPyPzE(myParticle->px(), myParticle->py(), myParticle->pz(), myParticle->energy() );
    genCon->SetVtx(myParticle->vx(), myParticle->vy(), myParticle->vz());
    genCon->SetCharge(myParticle->charge());
    genCon->SetPDGId( myParticle->pdgId());
    genCon->SetStatus(myParticle->status());
    map<const reco::GenParticle*,TCGenParticle*>::iterator momIt;

    if (myParticle->numberOfMothers() == 0){
      genCon->SetMother(0);
    }else if(
        abs(myParticle->mother()->pdgId()) != 5 
        && abs(myParticle->mother()->pdgId()) != 11
        && abs(myParticle->mother()->pdgId()) != 12
        && abs(myParticle->mother()->pdgId()) != 13
        && abs(myParticle->mother()->pdgId()) != 14
        && abs(myParticle->mother()->pdgId()) != 15
        && abs(myParticle->mother()->pdgId()) != 16
        && abs(myParticle->mother()->pdgId()) != 22
        && abs(myParticle->mother()->pdgId()) != 23 
        && abs(myParticle->mother()->pdgId()) != 24 
        && abs(myParticle->mother()->pdgId()) != 25 
        && abs(myParticle->mother()->pdgId()) != 35 
        && abs(myParticle->mother()->pdgId()) != 36 
        && abs(myParticle->mother()->pdgId()) != 39
             && abs(myParticle->mother()->pdgId()) != 443  //Jpsi
             && abs(myParticle->mother()->pdgId()) != 553  //Upsilon
        )
    {
      genCon->SetMother(0);
    }else{
      momIt = genMap.find((const reco::GenParticle*)myParticle->mother());
      if (momIt == genMap.end()){
        genCon->SetMother(addGenParticle((const reco::GenParticle*)myParticle->mother(), genPartCount, genMap));
      }else{
        genCon->SetMother(momIt->second);
      }
    }
  }
  else
    genCon = it->second;

  return genCon;
}

vector<float> ntupleProducer::getESHits(double X, double Y, double Z, map<DetId, EcalRecHit> rechits_map, const CaloSubdetectorGeometry*& geometry_p, CaloSubdetectorTopology *topology_p, int row) {

  //cout<<row<<endl;

  vector<float> esHits;

  //double X = bcPtr->x();
  //double Y = bcPtr->y();
  //double Z = bcPtr->z();
  const GlobalPoint point(X,Y,Z);

  DetId esId1 = (dynamic_cast<const EcalPreshowerGeometry*>(geometry_p))->getClosestCellInPlane(point, 1);
  DetId esId2 = (dynamic_cast<const EcalPreshowerGeometry*>(geometry_p))->getClosestCellInPlane(point, 2);
  ESDetId esDetId1 = (esId1 == DetId(0)) ? ESDetId(0) : ESDetId(esId1);
  ESDetId esDetId2 = (esId2 == DetId(0)) ? ESDetId(0) : ESDetId(esId2);

  map<DetId, EcalRecHit>::iterator it;
  ESDetId next;
  ESDetId strip1;
  ESDetId strip2;

  strip1 = esDetId1;
  strip2 = esDetId2;

  EcalPreshowerNavigator theESNav1(strip1, topology_p);
  theESNav1.setHome(strip1);

  EcalPreshowerNavigator theESNav2(strip2, topology_p);
  theESNav2.setHome(strip2);

  if (row == 1) {
    if (strip1 != ESDetId(0)) strip1 = theESNav1.north();
    if (strip2 != ESDetId(0)) strip2 = theESNav2.east();
  } else if (row == -1) {
    if (strip1 != ESDetId(0)) strip1 = theESNav1.south();
    if (strip2 != ESDetId(0)) strip2 = theESNav2.west();
  }

  // Plane 2
  if (strip1 == ESDetId(0)) {
    for (unsigned int i=0; i<31; ++i) esHits.push_back(0);
  } else {

    it = rechits_map.find(strip1);
    if (it->second.energy() > 1.0e-10 && it != rechits_map.end()) esHits.push_back(it->second.energy());
    else esHits.push_back(0);
    //cout<<"center : "<<strip1<<" "<<it->second.energy()<<endl;

    // east road
    for (unsigned int i=0; i<15; ++i) {
      next = theESNav1.east();
      if (next != ESDetId(0)) {
        it = rechits_map.find(next);
        if (it->second.energy() > 1.0e-10 && it != rechits_map.end()) esHits.push_back(it->second.energy());
        else esHits.push_back(0);
        //cout<<"east "<<i<<" : "<<next<<" "<<it->second.energy()<<endl;
      } else {
        for (unsigned int j=i; j<15; ++j) esHits.push_back(0);
        break;
        //cout<<"east "<<i<<" : "<<next<<" "<<0<<endl;
      }
    }

    // west road
    theESNav1.setHome(strip1);
    theESNav1.home();
    for (unsigned int i=0; i<15; ++i) {
      next = theESNav1.west();
      if (next != ESDetId(0)) {
        it = rechits_map.find(next);
        if (it->second.energy() > 1.0e-10 && it != rechits_map.end()) esHits.push_back(it->second.energy());
        else esHits.push_back(0);
        //cout<<"west "<<i<<" : "<<next<<" "<<it->second.energy()<<endl;
      } else {
        for (unsigned int j=i; j<15; ++j) esHits.push_back(0);
        break;
        //cout<<"west "<<i<<" : "<<next<<" "<<0<<endl;
      }
    }
  }

  if (strip2 == ESDetId(0)) {
    for (unsigned int i=0; i<31; ++i) esHits.push_back(0);
  } else {

    it = rechits_map.find(strip2);
    if (it->second.energy() > 1.0e-10 && it != rechits_map.end()) esHits.push_back(it->second.energy());
    else esHits.push_back(0);
    //cout<<"center : "<<strip2<<" "<<it->second.energy()<<endl;

    // north road
    for (unsigned int i=0; i<15; ++i) {
      next = theESNav2.north();
      if (next != ESDetId(0)) {
        it = rechits_map.find(next);
        if (it->second.energy() > 1.0e-10 && it != rechits_map.end()) esHits.push_back(it->second.energy());
        else esHits.push_back(0);
        //cout<<"north "<<i<<" : "<<next<<" "<<it->second.energy()<<endl;
      } else {
        for (unsigned int j=i; j<15; ++j) esHits.push_back(0);
        break;
        //cout<<"north "<<i<<" : "<<next<<" "<<0<<endl;
      }
    }

    // south road
    theESNav2.setHome(strip2);
    theESNav2.home();
    for (unsigned int i=0; i<15; ++i) {
      next = theESNav2.south();
      if (next != ESDetId(0)) {
        it = rechits_map.find(next);
        if (it->second.energy() > 1.0e-10 && it != rechits_map.end()) esHits.push_back(it->second.energy());
        else esHits.push_back(0);
        //cout<<"south "<<i<<" : "<<next<<" "<<it->second.energy()<<endl;
      } else {
        for (unsigned int j=i; j<15; ++j) esHits.push_back(0);
        break;
        //cout<<"south "<<i<<" : "<<next<<" "<<0<<endl;
      }
    }
  }

  return esHits;
}

vector<float> ntupleProducer::getESEffSigmaRR(vector<float> ESHits0)
{
  const int nBIN = 21;
  vector<float> esShape;

  TH1F *htmpF = new TH1F("htmpF","",nBIN,0,nBIN);
  TH1F *htmpR = new TH1F("htmpR","",nBIN,0,nBIN);
  htmpF->Reset(); htmpR->Reset();

  Float_t effsigmaRR=0.;

  for(int ibin=0; ibin<((nBIN+1)/2); ++ibin) {
    if (ibin==0) {
      htmpF->SetBinContent((nBIN+1)/2,ESHits0[ibin]);
      htmpR->SetBinContent((nBIN+1)/2,ESHits0[ibin+31]);
    } else { // hits sourd the seed
      htmpF->SetBinContent((nBIN+1)/2+ibin,ESHits0[ibin]);
      htmpF->SetBinContent((nBIN+1)/2-ibin,ESHits0[ibin+15]);
      htmpR->SetBinContent((nBIN+1)/2+ibin,ESHits0[ibin+31]);
      htmpR->SetBinContent((nBIN+1)/2-ibin,ESHits0[ibin+31+15]);
    }
  }

  // ---- Effective Energy Deposit Width ---- //
  double EffWidthSigmaXX = 0.;
  double EffWidthSigmaYY = 0.;
  double totalEnergyXX = 0.;
  double totalEnergyYY = 0.;
  double EffStatsXX = 0.;
  double EffStatsYY = 0.;
  for (int id_X=1; id_X<=21; ++id_X) {
    totalEnergyXX += htmpF->GetBinContent(id_X);
    EffStatsXX += htmpF->GetBinContent(id_X)*(id_X-11)*(id_X-11);
    totalEnergyYY += htmpR->GetBinContent(id_X);
    EffStatsYY += htmpR->GetBinContent(id_X)*(id_X-11)*(id_X-11);
  }
  // If denominator == 0, effsigmaRR = 0;
  EffWidthSigmaXX = (totalEnergyXX>0.) ? sqrt(fabs(EffStatsXX / totalEnergyXX)) : 0.;
  EffWidthSigmaYY = (totalEnergyYY>0.) ? sqrt(fabs(EffStatsYY / totalEnergyYY)) : 0.;
  effsigmaRR = ((totalEnergyXX + totalEnergyYY) >0.) ? sqrt(EffWidthSigmaXX * EffWidthSigmaXX + EffWidthSigmaYY * EffWidthSigmaYY) : 0.;
  esShape.push_back(effsigmaRR);
  esShape.push_back(EffWidthSigmaXX);
  esShape.push_back(EffWidthSigmaYY);

  delete htmpF;
  delete htmpR;

  return esShape;
}

/*

TCTrack::ConversionInfo ntupleProducer::CheckForConversions(const edm::Handle<reco::ConversionCollection> &convCol,
                                                            const reco::GsfTrackRef &gsf,
                                                            const math::XYZPoint &bs, const math::XYZPoint &pv)
{
  TCTrack::ConversionInfo * convInfo = new TCTrack::ConversionInfo();
  //int iconv=-1;
  for (reco::ConversionCollection::const_iterator conv = convCol->begin(); conv!= convCol->end(); ++conv) {
    //iconv++;
    
    reco::Vertex vtx = conv->conversionVertex();
    if (vtx.isValid()) {
      if (ConversionTools::matchesConversion(gsf, *conv)) {
        
        (*convInfo).isValid = true;
        
        (*convInfo).vtxProb = TMath::Prob( vtx.chi2(), vtx.ndof() );
	math::XYZVector mom(conv->refittedPairMomentum());
        double dbsx = vtx.x() - bs.x();
        double dbsy = vtx.y() - bs.y();
        (*convInfo).lxyBS = (mom.x()*dbsx + mom.y()*dbsy)/mom.rho();
        
        double dpvx = vtx.x() - pv.x();
        double dpvy = vtx.y() - pv.y();
        (*convInfo).lxyPV = (mom.x()*dpvx + mom.y()*dpvy)/mom.rho();
        
        (*convInfo).nHitsMax=0;
        for (std::vector<uint8_t>::const_iterator it = conv->nHitsBeforeVtx().begin(); it!=conv->nHitsBeforeVtx().end(); ++it) {
          if ((*it)>(*convInfo).nHitsMax) (*convInfo).nHitsMax = (*it);
        }
        
        break;
      }
    }
  }
  return (*convInfo);
}

*/



//define this as a plug-in
DEFINE_FWK_MODULE(ntupleProducer);
