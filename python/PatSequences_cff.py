import FWCore.ParameterSet.Config as cms
from PhysicsTools.PatAlgos.tools.pfTools import *
from PhysicsTools.PatAlgos.tools.trigTools import *
from RecoJets.Configuration.RecoPFJets_cff import kt6PFJets
from CommonTools.ParticleFlow.Tools.enablePileUpCorrection import enablePileUpCorrection

##                                      ## 
## add trigger matching for the leptons ##
##                                      ##

def addTriggerMatchingForLeptons(process, postfix='') :
    # define the trigger matchers
    process.muTriggerMatchPF = cms.EDProducer( "PATTriggerMatcherDRLessByR",
                                               src     = cms.InputTag( "selectedPatMuons"+postfix ),
                                               matched = cms.InputTag( "patTrigger" ),
                                               matchedCuts = cms.string( 'type( "TriggerMuon" ) && ( path("HLT_Mu8_*") || path("HLT_Mu12_*") || path("HLT_Mu13_*") || path("HLT_DoubleMu7_*") )' ), 
                                               maxDPtRel   = cms.double( 0.5 ), # no effect here
                                               maxDeltaR   = cms.double( 0.5 ),
                                               maxDeltaEta = cms.double( 0.2 ), # no effect here
                                               # definition of matcher output
                                               resolveAmbiguities    = cms.bool( False ),
                                               resolveByMatchQuality = cms.bool( False )
                                               )
    
    process.eleTriggerMatchPF = cms.EDProducer( "PATTriggerMatcherDRLessByR",
                                                src     = cms.InputTag( "selectedPatElectrons"+postfix ),
                                                matched = cms.InputTag( "patTrigger" ),
                                                #matchedCuts = cms.string( 'type( "TriggerL1NoIsoEG" ) || type( "TriggerL1IsoEG" ) || type( "TriggerElectron" )' ),
                                                matchedCuts = cms.string( 'type( "TriggerElectron" )' ),
                                                maxDPtRel   = cms.double( 0.5 ), # no effect here
                                                maxDeltaR   = cms.double( 0.5 ),
                                                maxDeltaEta = cms.double( 0.2 ), # no effect here
                                                # definition of matcher output
                                                resolveAmbiguities    = cms.bool( False ),
                                                resolveByMatchQuality = cms.bool( False )
                                                )

    from PhysicsTools.PatAlgos.tools.coreTools import removeCleaning
    removeCleaning( process )
    setattr( process, 'muTriggerMatch' + postfix, process.muTriggerMatchPF )
    setattr( process, 'eleTriggerMatch' + postfix, process.eleTriggerMatchPF )
    switchOnTriggerMatching( process, triggerMatchers = [ 'muTriggerMatchPFlow','eleTriggerMatchPFlow' ], sequence = 'patPF2PATSequence' + postfix )
    removeCleaningFromTriggerMatching( process, sequence = 'patPF2PATSequence' + postfix )

##                   ##
## adds pat sequence ##
##                   ##

def addPatSequence(process, runOnMC, addPhotons=True) :

    #PF2PAT
    postfix = "PFlow"
    jetAlgo='AK5'
    jecSetPF = jetAlgo+'PFchs'
    jecLevels=['L1FastJet','L2Relative','L3Absolute']
    if(not runOnMC) : jecLevels.append( 'L2L3Residual' )

    #enablePileUpCorrection(process, postfix=postfix)

    #start PF2PAT
    usePF2PAT(process, runPF2PAT=True,
              jetAlgo=jetAlgo, runOnMC= runOnMC, postfix=postfix,
              jetCorrections=(jecSetPF, jecLevels), 
              typeIMetCorrections=True, pvCollection=cms.InputTag('goodOfflinePrimaryVertices')
              )

    #process.pfPileUpPFlow.checkClosestZVertex = False

    '''
    # configure muons
    getattr(process,"patMuons"+postfix).embedCaloMETMuonCorrs = False 
    getattr(process,"patMuons"+postfix).embedTcMETMuonCorrs = False
    getattr(process,"patMuons"+postfix).embedTrack = True
    getattr(process,"pfMuonsFromVertex"+postfix).dzCut = 99
    getattr(process,"pfMuonsFromVertex"+postfix).d0Cut = 99
    getattr(process,"pfSelectedMuons"+postfix).cut="pt()>3"

    # configure electrons
    useGsfElectrons(process,postfix)
    getattr(process,"patElectrons"+postfix).embedTrack = True
    getattr(process,"pfElectronsFromVertex"+postfix).dzCut = 99
    getattr(process,"pfElectronsFromVertex"+postfix).d0Cut = 99
    getattr(process,"pfSelectedElectrons"+postfix).cut="pt()>5"

    #electron ID
    process.load("ElectroWeakAnalysis.WENu.simpleEleIdSequence_cff")
    process.load("RecoEgamma.ElectronIdentification.cutsInCategoriesElectronIdentificationV06_DataTuning_cfi")
    process.load("RecoEgamma.ElectronIdentification.cutsInCategoriesElectronIdentificationV06_cfi")
    process.electronIDSequence = cms.Sequence(
        process.simpleEleIdSequence +
        process.eidVeryLoose +
        process.eidLoose +
        process.eidMedium +
        process.eidTight +
        process.eidSuperTight+
        process.eidVeryLooseMC +
        process.eidLooseMC +
        process.eidMediumMC +
        process.eidTightMC +
        process.eidSuperTightMC
        )

    applyPostfix( process, 'patElectrons', postfix ).electronIDSources = cms.PSet(
        eidVBTF95 = cms.InputTag("simpleEleId95relIso"),
        eidVBTF90 = cms.InputTag("simpleEleId90relIso"),
        eidVBTF85 = cms.InputTag("simpleEleId85relIso"),
        eidVBTF80 = cms.InputTag("simpleEleId80relIso"),
        eidVBTF70 = cms.InputTag("simpleEleId70relIso"),
        eidVBTF60 = cms.InputTag("simpleEleId60relIso"),
        eidVeryLoose = cms.InputTag("eidVeryLoose"),
        eidLoose = cms.InputTag("eidLoose"),
        eidMedium = cms.InputTag("eidMedium"),
        eidTight = cms.InputTag("eidTight"),
        eidSuperTight = cms.InputTag("eidSuperTight"),
        eidVeryLooseMC = cms.InputTag("eidVeryLooseMC"),
        eidLooseMC = cms.InputTag("eidLooseMC"),
        eidMediumMC = cms.InputTag("eidMediumMC"),
        eidTightMC = cms.InputTag("eidTightMC"),
        eidSuperTightMC = cms.InputTag("eidSuperTightMC")       
        )


    # configure jets
    enablePileUpCorrection( process, postfix=postfix)
    getattr(process,"patJetCorrFactors"+postfix).payload = jetAlgoPayLoad 
    getattr(process,"patJets"+postfix).embedPFCandidates = cms.bool(True)
    getattr(process,"patJets"+postfix).embedCaloTowers   = cms.bool(True)
    
    # use non pileup substracted rho as in the Jan2012 JEC set
    getattr(process,"patJetCorrFactors"+postfix).rho = cms.InputTag("kt6PFJets","rho")

    #configure top projections
    getattr(process,"pfNoPileUp"+postfix).enable    = True
    getattr(process,"pfNoMuon"+postfix).enable      = True
    getattr(process,"pfNoMuon"+postfix).verbose     = False
    getattr(process,"pfNoElectron"+postfix).enable  = True
    getattr(process,"pfNoTau"+postfix).enable       = False
    getattr(process,"pfNoJet"+postfix).enable       = False


    ### adding standard muons and electrons
    process.patMuons.embedTcMETMuonCorrs = False
    process.patMuons.embedCaloMETMuonCorrs = False
    process.patMuons.embedTrack = True

    process.patElectrons.pfElectronSource = 'particleFlow'
    process.patElectrons.embedTrack = True

    process.eleIsoSequence = setupPFElectronIso(process, 'gsfElectrons', 'PFIso')
    process.muIsoSequence = setupPFMuonIso(process, 'muons', 'PFIso')
    adaptPFIsoMuons( process, applyPostfix(process,"patMuons",""), 'PFIso')
    adaptPFIsoElectrons( process, applyPostfix(process,"patElectrons",""), 'PFIso')
    process.stdMuonSeq = cms.Sequence( process.pfParticleSelectionSequence +
                                       process.muIsoSequence +
                                       process.makePatMuons +
                                       process.selectedPatMuons
                                       )
    process.stdElectronSeq = cms.Sequence( process.pfParticleSelectionSequence +
                                           process.eleIsoSequence +
                                           process.makePatElectrons +
                                           process.selectedPatElectrons
                                           )
    if(runOnMC) :
        process.stdPhotonSeq = cms.Sequence( process.makePatPhotons )
    else :
        process.patPhotons.addGenMatch = cms.bool(False)
        process.patPhotons.embedGenMatch = cms.bool(False)
        process.stdPhotonSeq = cms.Sequence( process.patPhotons ) 

    #add secondary vertex mass to jets
    applyPostfix( process, 'patJets', postfix ).tagInfoSources = cms.VInputTag( cms.InputTag("secondaryVertexTagInfosAOD"+postfix) )
    applyPostfix( process, 'patJets', postfix ).userData.userFunctions = cms.vstring( "? hasTagInfo('secondaryVertex') && tagInfoSecondaryVertex('secondaryVertex').nVertices() > 0 ? tagInfoSecondaryVertex('secondaryVertex').secondaryVertex(0).p4().mass() : -999")
    applyPostfix( process, 'patJets', postfix ).userData.userFunctionLabels = cms.vstring('secvtxMass')
    '''
    
    #create the path
    process.patDefaultSequence = cms.Sequence(
       #process.electronIDSequence
        #* process.kt6PFJets25
        getattr(process,"patPF2PATSequence"+postfix)
        #* process.stdMuonSeq 
        #* process.stdElectronSeq 
        #* process.stdPhotonSeq
        #* process.patPhotons
        )
        
    print " *** PAT path has been defined"
    
