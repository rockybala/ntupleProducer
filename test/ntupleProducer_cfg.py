import os
import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
from RecoEgamma.PhotonIdentification.isolationCalculator_cfi import *

process = cms.Process("NTUPLE")

options = VarParsing.VarParsing ('analysis')
options.maxEvents = 1000
options.inputFiles = '/store/mc/Summer12_DR53X/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S10_START53_V7A-v1/0000/02CDCF05-BED2-E111-85F4-0030486740BA.root'
#options.inputFiles = '/store/data/Run2012D/SinglePhotonParked/AOD/22Jan2013-v1/30004/144D7268-4086-E211-9DC1-001E673984C1.root'

options.register("isRealData",
                 0,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 "0 if running on MC and 1 if running on Data")

options.parseArguments()

# real data or MC?
isRealData = options.isRealData

from NWU.ntupleProducer.patTemplate_cfg import *
from PhysicsTools.PatAlgos.tools.coreTools import *

if (isRealData):
    removeMCMatching(process, ['All'])

from PhysicsTools.PatAlgos.tools.metTools import *
addPfMET(process, 'PF')

from PhysicsTools.PatAlgos.tools.jetTools import *

# global tag
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')

if (isRealData):
    process.GlobalTag.globaltag = 'FT_53_V21_AN3::All'
else:
    process.GlobalTag.globaltag = 'START53_V27::All'

# Create good primary vertices for PF association
from PhysicsTools.SelectorUtils.pvSelector_cfi import pvSelector
process.goodOfflinePrimaryVertices = cms.EDFilter( "PrimaryVertexObjectFilter",
    filterParams = pvSelector.clone( minNdof = cms.double(4.0), maxZ = cms.double(24.0) ),
    src=cms.InputTag('offlinePrimaryVertices') #Standard Primary Vertex Collection
#    src=cms.InputTag('offlinePrimaryVerticesWithBS') #Primary Vertices Collection constrained by beamspot
                                                   )

# jet energy corrections
process.load('JetMETCorrections.Configuration.DefaultJEC_cff')
process.load("CondCore.DBCommon.CondDBCommon_cfi")

# Switch to PF jets
from PhysicsTools.PatAlgos.tools.jetTools import *
if isRealData == True:
    switchJetCollection(process,
                        cms.InputTag('ak5PFJets'),
                        doJTA            = True,
                        doBTagging       = False,
                        jetCorrLabel     = ('AK5PF', ['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual']),
                        doType1MET       = True,
                        genJetCollection = cms.InputTag("ak5GenJets"),
                        doJetID          = False,
                        jetIdLabel       = "ak5"
                        )
    
else:
    switchJetCollection(process,
                        cms.InputTag('ak5PFJets'),
                        doJTA            = True,
                        doBTagging       = False,
                        jetCorrLabel     = ('AK5PF', ['L1FastJet', 'L2Relative', 'L3Absolute']),
                        doType1MET       = True,
                        genJetCollection = cms.InputTag("ak5GenJets"),
                        doJetID          = False,
                        jetIdLabel       = "ak5"
                        )

from PhysicsTools.PatAlgos.tools.pfTools import *

#Select Good photons with a pt > 30. Missing Isolation and Conversion Safe Veto
process.goodPhotonsHighPtCut = cms.EDFilter("PATPhotonSelector",
                                            src = cms.InputTag("patPhotons"),
                                            cut = cms.string('pt > 30.0 && hadronicOverEm < 0.5 && r9 > 0.9 && abs(superCluster.eta) < 1.4442 && sigmaIetaIeta < 0.011 '),
                                            filter = cms.bool(False)
                                            )

if isRealData == True:
    from PhysicsTools.PatUtils.tools.runType1PFMEtUncertainties     import runType1PFMEtUncertainties
    runType1PFMEtUncertainties(process,
                               electronCollection = '',
                               photonCollection = cms.InputTag('goodPhotonsHighPtCut'),
                               muonCollection = '',
                               tauCollection = '',
                               jetCollection = cms.InputTag('patJets'),
                               jetCorrLabel = 'L2L3Residual',
                               jecUncertaintyFile = "PhysicsTools/PatUtils/data/Summer13_V1_DATA_UncertaintySources_AK5PF.txt",
                               makeType1corrPFMEt = True,
                               makeType1p2corrPFMEt = False,
                               doApplyType0corr = True,
                               doSmearJets = False,
                               addToPatDefaultSequence = False,
                               postfix = ''
                                               )
    process.patPFJetMETtype1p2Corr.jetCorrLabel = cms.string('L2L3Residual')
    process.patPFJetMETtype2Corr.jetCorrLabel = cms.string('L2L3Residual')
                                    
    from PhysicsTools.PatUtils.tools.runMVAMEtUncertainties         import runMVAMEtUncertainties
    runMVAMEtUncertainties(process,
                           electronCollection = '',
                           photonCollection = cms.InputTag('goodPhotonsHighPtCut'),
                           muonCollection = '',
                           tauCollection = '',
                           jetCollection = cms.InputTag('patJets'),
                           jecUncertaintyFile = "PhysicsTools/PatUtils/data/Summer13_V1_DATA_UncertaintySources_AK5PF.txt",
                           doSmearJets = False,
                           addToPatDefaultSequence = False,
                           postfix = ''
                           )
    process.calibratedAK5PFJetsForPFMEtMVA.correctors = cms.vstring("ak5PFL1FastL2L3Residual")
    process.pfMEtMVA.srcLeptons = cms.VInputTag('goodPhotonsHighPtCut')

else:
    from PhysicsTools.PatUtils.tools.runType1PFMEtUncertainties     import runType1PFMEtUncertainties
    runType1PFMEtUncertainties(process,
                               electronCollection = '',
                               photonCollection = cms.InputTag('goodPhotonsHighPtCut'),
                               muonCollection = '',
                               tauCollection = '',
                               jetCollection = cms.InputTag('patJets'),
                               jetCorrLabel = 'L3Absolute',
                               jecUncertaintyFile = "PhysicsTools/PatUtils/data/Summer13_V1_DATA_UncertaintySources_AK5PF.txt",
                               makeType1corrPFMEt = True,
                               makeType1p2corrPFMEt = False,
                               doApplyType0corr = True,
                               doSmearJets = True,
                               addToPatDefaultSequence = False,
                               postfix = ''
                               )
    
    from PhysicsTools.PatUtils.tools.runMVAMEtUncertainties         import runMVAMEtUncertainties
    runMVAMEtUncertainties(process,
                           electronCollection = '',
                           photonCollection = cms.InputTag('goodPhotonsHighPtCut'),
                           muonCollection = '',
                           tauCollection = '',
                           jetCollection = cms.InputTag('patJets'),
                           jecUncertaintyFile = "PhysicsTools/PatUtils/data/Summer13_V1_DATA_UncertaintySources_AK5PF.txt",
                           doSmearJets = True,
                           addToPatDefaultSequence = False,
                           postfix = ''
                           )
    process.calibratedAK5PFJetsForPFMEtMVA.correctors = cms.vstring('ak5PFL1FastL2L3')
    process.pfMEtMVA.srcLeptons = cms.VInputTag('goodPhotonsHighPtCut')
    process.pfMEtMVA.verbosity = cms.int32(0)

if (isRealData):
    #process.pfJetMETcorr.jetCorrLabel = cms.string("ak5PFL1FastL2L3Residual")

    # 53X b-jet discriminator calibration
    process.GlobalTag.toGet = cms.VPSet(
        cms.PSet(record = cms.string("BTagTrackProbability2DRcd"),
            tag = cms.string("TrackProbabilityCalibration_2D_Data53X_v2"),
            connect = cms.untracked.string("frontier://FrontierPrep/CMS_COND_BTAU")),
        cms.PSet(record = cms.string("BTagTrackProbability3DRcd"),
            tag = cms.string("TrackProbabilityCalibration_3D_Data53X_v2"),
            connect = cms.untracked.string("frontier://FrontierPrep/CMS_COND_BTAU"))
        )
else:
    # 53X b-jet discriminator calibration
    process.GlobalTag.toGet = cms.VPSet(
        cms.PSet(record = cms.string("BTagTrackProbability2DRcd"),
            tag = cms.string("TrackProbabilityCalibration_2D_Data53X_v2"),
            connect = cms.untracked.string("frontier://FrontierPrep/CMS_COND_BTAU")),
        cms.PSet(record = cms.string("BTagTrackProbability3DRcd"),
            tag = cms.string("TrackProbabilityCalibration_3D_Data53X_v2"),
            connect = cms.untracked.string("frontier://FrontierPrep/CMS_COND_BTAU"))
        )

### To get b-tags from ak5PFJets
process.load('RecoJets.JetAssociationProducers.ak5JTA_cff')
process.ak5JetTracksAssociatorAtVertex.jets = cms.InputTag("ak5PFJetsL1FastL2L3")
process.ak5JetTracksAssociatorAtCaloFace.jets = cms.InputTag("ak5PFJetsL1FastL2L3")
process.ak5JetExtender.jets = cms.InputTag("ak5PFJetsL1FastL2L3")

process.kt6PFJetsIso = process.kt6PFJets.clone()
process.kt6PFJetsIso.doRhoFastjet = True
process.kt6PFJetsIso.Rho_EtaMax = cms.double(2.5)


# for jet pileup ID variables
from RecoJets.JetProducers.PileupJetIDParams_cfi import *

process.load("RecoJets.JetProducers.PileupJetIDSequence_cff")
#from RecoJets.JetProducers.PileupJetIDSequence_cff import *

## The final MVA has to be evaluated on corrected jets, thus the JEC have to be fed to the producer. The producer can run both on corrected and uncorrected jets, provided that the parameters are properly set. The snippet below shows how to run on a collection of uncorrected jets. To run on corrected reco jets, set the applyJec flag to False and the inputIsCorrected to True. 
#
#process.recoPuJetId = puJetId.clone(
#    jets = cms.InputTag("ak5PFJets"),
    #jets = cms.InputTag("ak5PFJetsL1FastL2L3"),
#    applyJec = cms.bool(True),
#    inputIsCorrected = cms.bool(False),
#       )

##process.recoPuJetMva = puJetMva.clone(
#       jets = cms.InputTag("ak5PFJets"),
#       #jets = cms.InputTag("ak5PFJetsL1FastL2L3"),
#       jetids = cms.InputTag("recoPuJetId"),
#       applyJec = cms.bool(True),
#       inputIsCorrected = cms.bool(False),
#       )



#recoPuJetIdSqeuence = cms.Sequence(process.recoPuJetId * process.recoPuJetMva )


# global options
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.cerr.FwkReport.reportEvery = 500

'''
process.MessageLogger.categories = cms.untracked.vstring('FwkJob', 'FwkReport', 'FwkSummary', 'Root_NoDictionary', 'DataNotAvailable', 'HLTConfigData')
process.MessageLogger.destinations = cms.untracked.vstring('myOutput')
process.MessageLogger.myOutput = cms.untracked.PSet(
FwkJob              = cms.untracked.PSet(limit = cms.untracked.int32(0)),
FwkReport           = cms.untracked.PSet(limit = cms.untracked.int32(0)),
FwkSummary          = cms.untracked.PSet(limit = cms.untracked.int32(0)),
Root_NoDictionary   = cms.untracked.PSet(limit = cms.untracked.int32(0)),
DataNotAvailable    = cms.untracked.PSet(limit = cms.untracked.int32(0)),
HLTConfigData       = cms.untracked.PSet(limit = cms.untracked.int32(0))
)
'''

process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(False),
                                     SkipEvent = cms.untracked.vstring('ProductNotFound')
                                     )

# event counters
process.startCounter = cms.EDProducer("EventCountProducer")
process.endCounter = process.startCounter.clone()

##############################################
#### Met/Noise/BeamHalo filters  #############
## Following the recipe from this twiki:
## https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFilters
## Using *Tagging mode* for all filters (produces the boolean instead of filtering an event)
## Saving this boolean in the ntuples!
##############################################

## The iso-based HBHE noise filter
process.load('CommonTools/RecoAlgos/HBHENoiseFilterResultProducer_cfi')

## This might need updating for the VBF Parked - waiting for the recommendation
process.load("RecoMET.METFilters.hcalLaserEventFilter_cfi")
process.hcalLaserEventFilter.vetoByRunEventNumber=cms.untracked.bool(False)
process.hcalLaserEventFilter.vetoByHBHEOccupancy=cms.untracked.bool(True)
process.hcalLaserEventFilter.taggingMode = cms.bool(True)

## Ecal Dead Cell Filter
process.load('RecoMET.METFilters.EcalDeadCellTriggerPrimitiveFilter_cfi')
## For AOD and RECO recommendation to use recovered rechits
process.EcalDeadCellTriggerPrimitiveFilter.tpDigiCollection = cms.InputTag("ecalTPSkimNA")
process.EcalDeadCellTriggerPrimitiveFilter.taggingMode = cms.bool(True)

# The section below is for the filter on Boundary Energy. Available in AOD in CMSSW>44x
process.load('RecoMET.METFilters.EcalDeadCellBoundaryEnergyFilter_cfi')
process.EcalDeadCellBoundaryEnergyFilter.taggingMode = cms.bool(True)
process.EcalDeadCellBoundaryEnergyFilter.cutBoundEnergyDeadCellsEB=cms.untracked.double(10)
process.EcalDeadCellBoundaryEnergyFilter.cutBoundEnergyDeadCellsEE=cms.untracked.double(10)
process.EcalDeadCellBoundaryEnergyFilter.cutBoundEnergyGapEB=cms.untracked.double(100)
process.EcalDeadCellBoundaryEnergyFilter.cutBoundEnergyGapEE=cms.untracked.double(100)
process.EcalDeadCellBoundaryEnergyFilter.enableGap=cms.untracked.bool(False)
process.EcalDeadCellBoundaryEnergyFilter.limitDeadCellToChannelStatusEB = cms.vint32(12,14)
process.EcalDeadCellBoundaryEnergyFilter.limitDeadCellToChannelStatusEE = cms.vint32(12,14)
# End of Boundary Energy filter configuration

## The Good vertices collection needed by the tracking failure filter
process.goodVertices = cms.EDFilter(
      "VertexSelector",
        filter = cms.bool(False),
        src = cms.InputTag("offlinePrimaryVertices"),
#        src = cms.InputTag("offlinePrimaryVerticesWithBS"),
        cut = cms.string("!isFake && ndof > 4 && abs(z) <= 24 && position.rho < 2")
      )
## The tracking failure filter
process.load('RecoMET.METFilters.trackingFailureFilter_cfi')
process.trackingFailureFilter.taggingMode = cms.bool(True)

#Bad EE SC filter, not needed but goot to have them
process.load('RecoMET.METFilters.eeBadScFilter_cfi')
process.eeBadScFilter.tagginMode = cms.bool(True)

## The tracking POG filters
process.load('RecoMET.METFilters.trackingPOGFilters_cff')
## NOTE: to make tagging mode of the tracking POG filters (three of them), please do:
process.manystripclus53X.taggedMode = cms.untracked.bool(True)
process.manystripclus53X.forcedValue = cms.untracked.bool(False)
process.toomanystripclus53X.taggedMode = cms.untracked.bool(True)
process.toomanystripclus53X.forcedValue = cms.untracked.bool(False)
process.logErrorTooManyClusters.taggedMode = cms.untracked.bool(True)
process.logErrorTooManyClusters.forcedValue = cms.untracked.bool(False)
## Also the stored boolean for the three filters is opposite to what we usually
## have for other filters, i.e., true means rejected bad events while false means
## good events.

AllFilters = cms.Sequence(process.HBHENoiseFilterResultProducer
                          * process.hcalLaserEventFilter
                          * process.EcalDeadCellTriggerPrimitiveFilter
                          * process.EcalDeadCellBoundaryEnergyFilter
                          * process.goodVertices * process.trackingFailureFilter
                          * process.eeBadScFilter
                          * ~process.manystripclus53X #trkPOGFilter1
                          * ~process.toomanystripclus53X #trkPOGFilter2
                          * ~process.logErrorTooManyClusters #trkPOGFilter 3
                          )
##### END OF Noise Filters ############


# Electron MVA ID producer:
#process.load('EgammaAnalysis/ElectronTools/electronIdMVAProducer_cfi')

# Electron Regression (post moriond recommendation)
#process.load('EgammaAnalysis/ElectronTools/electronRegressionEnergyProducer_cfi')
#process.eleRegressionEnergy.inputElectronsTag    = cms.InputTag('gsfElectrons')
#process.eleRegressionEnergy.inputCollectionType  = cms.uint32(0)
#process.eleRegressionEnergy.useRecHitCollections = cms.bool(True)
#process.eleRegressionEnergy.produceValueMaps     = cms.bool(True)
#process.eleRegressionEnergy.energyRegressionType = cms.uint32(2)
#process.eleRegressionEnergy.regressionInputFile = cms.string("EgammaAnalysis/ElectronTools/data/eleEnergyRegWeights_WithSubClusters_VApr15.root")

# Electron Combination for calibration (post moriond recommendation)
#process.load('EgammaAnalysis/ElectronTools/calibratedElectrons_cfi')
#process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
#    calibratedElectrons = cms.PSet(
#      initialSeed = cms.untracked.uint32(1),
#      engineName = cms.untracked.string('TRandom3')
#    ),
#)
#if (isRealData):
#  process.calibratedElectrons.isMC = cms.bool(False)
#  process.calibratedElectrons.inputDataset = cms.string("22Jan2013ReReco")
#else:
#  process.calibratedElectrons.isMC = cms.bool(True)
#  process.calibratedElectrons.inputDataset = cms.string("Summer12_LegacyPaper")

#process.calibratedElectrons.updateEnergyError = cms.bool(True)
#process.calibratedElectrons.correctionsType   = cms.int32(2)
#process.calibratedElectrons.combinationType   = cms.int32(3)
#process.calibratedElectrons.lumiRatio         = cms.double(1.0)
#process.calibratedElectrons.verbose           = cms.bool(False)
#process.calibratedElectrons.synchronization   = cms.bool(False)
#process.calibratedElectrons.applyLinearityCorrection = cms.bool(True)
#process.calibratedElectrons.scaleCorrectionsInputPath = cms.string("EgammaAnalysis/ElectronTools/data/scalesMoriond.csv")
#process.calibratedElectrons.combinationRegressionInputPath = cms.string("EgammaAnalysis/ElectronTools/data/eleEnergyReg2012Weights_V1.root")


# event source
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(options.maxEvents))
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(options.inputFiles)
)

print '\n\nCommence ntuplization...\n\n'

### TFile service!
process.TFileService = cms.Service('TFileService',
    fileName = cms.string('nuTuple.root')
    #fileName = cms.string('~/EOS/V08_01_8TeV/ggHZG_M125_Pythia8_175_POWHEG_PDF7/nuTuple_9.root')
                                   )

### pfNoPU Sequence for electron MVA
process.pfPileUp = cms.EDProducer("PFPileUp",
    PFCandidates = cms.InputTag("particleFlow"),
    Enable = cms.bool(True),
    checkClosestZVertex = cms.bool(True),
    verbose = cms.untracked.bool(False),
    Vertices = cms.InputTag("offlinePrimaryVertices")
)

process.pfNoPileUp = cms.EDProducer("TPPFCandidatesOnPFCandidates",
    bottomCollection = cms.InputTag("particleFlow"),
    enable = cms.bool(True),
    topCollection = cms.InputTag("pfPileUp"),
    name = cms.untracked.string('pileUpOnPFCandidates'),
    verbose = cms.untracked.bool(False)
)

process.pfNoPUSeq = cms.Sequence(process.pfPileUp + process.pfNoPileUp)

from SHarper.HEEPAnalyzer.HEEPSelectionCuts_cfi import *
process.heepIdNoIso = cms.EDProducer("HEEPIdValueMapProducer",
                                     eleLabel = cms.InputTag("gsfElectrons"),
                                     barrelCuts = cms.PSet(heepBarrelCuts),
                                     endcapCuts = cms.PSet(heepEndcapCuts),
                                     eleIsolEffectiveAreas = cms.PSet(heepEffectiveAreas),
                                     eleRhoCorrLabel = cms.InputTag("kt6PFJets", "rho"),
                                     verticesLabel = cms.InputTag("offlinePrimaryVertices"),
                                     applyRhoCorrToEleIsol = cms.bool(True),
                                     writeIdAsInt = cms.bool(True)
                                     )
process.heepIdNoIso.barrelCuts.cuts=cms.string("et:detEta:ecalDriven:dEtaIn:dPhiIn:hadem:e2x5Over5x5:nrMissHits:dxy")
process.heepIdNoIso.endcapCuts.cuts=cms.string("et:detEta:ecalDriven:dEtaIn:dPhiIn:hadem:sigmaIEtaIEta:nrMissHits:dxy")

process.heepIdNoIsoEles = cms.EDProducer("tsw::HEEPGsfProducer", cutValueMap = cms.InputTag("heepIdNoIso"),
                                         inputGsfEles = cms.InputTag("gsfElectrons")  )

# Boosted Z ModEleIso: 1b) Calculating the modified iso. values using BstdZeeTools EDProducer

from TSWilliams.BstdZeeTools.bstdzeemodisolproducer_cff import *
process.modElectronIso = cms.EDProducer("BstdZeeModIsolProducer",
                                        bstdZeeModIsolParams, vetoGsfEles = cms.InputTag("heepIdNoIsoEles") )

### ntuple producer
process.ntupleProducer   = cms.EDAnalyzer('ntupleProducer',

                                          verboseTrigs         =    cms.untracked.bool(False),
                                          verboseMVAs          =    cms.untracked.bool(False),
                                          
                                          photonIsoCalcTag  =    cms.PSet(isolationSumsCalculator),
                                          jetPUIdAlgo       =    cms.PSet(full_5x),

                                                                                    
                                          #All MET
                                          srcPatPhoup         = cms.InputTag("shiftedGoodPhotonsHighPtCutEnUp"),
                                          srcPatPhodown       = cms.InputTag("shiftedGoodPhotonsHighPtCutEnDown"),
                                          srcPatPhoMvaup      = cms.InputTag("shiftedGoodPhotonsHighPtCutEnUpForPFMEtByMVA"),
                                          srcPatPhoMvadown    = cms.InputTag("shiftedGoodPhotonsHighPtCutEnDownForPFMEtByMVA"),
                                          
                                          srcPatJets          = cms.InputTag("patJetsNotOverlappingWithLeptonsForJetMEtUncertainty"),#in our case photons only
                                          srcSelectedJets     = cms.InputTag("selectedPatJets"),
                                          srcSmearedJets      = cms.InputTag("smearedPatJets"),
                                          srcMetRaw           = cms.InputTag("patPFMet"),                                
                                          srcMetPf            = cms.InputTag("patMETsPF"),
                                          
                                          ################################################
                                          srcMetCorrected     = cms.InputTag("patType1CorrectedPFMet"),
                                          srcMetJERup         = cms.InputTag("patType1CorrectedPFMetJetResUp"),
                                          srcMetJERdown       = cms.InputTag("patType1CorrectedPFMetJetResDown"),
                                          srcMetPhoup         = cms.InputTag("patType1CorrectedPFMetPhotonEnUp"),
                                          srcMetPhodown       = cms.InputTag("patType1CorrectedPFMetPhotonEnDown"),
                                          srcMetJetup         = cms.InputTag("patType1CorrectedPFMetJetEnUp"),
                                          srcMetJetdown       = cms.InputTag("patType1CorrectedPFMetJetEnDown"),
                                          srcMetUncup         = cms.InputTag("patType1CorrectedPFMetUnclusteredEnUp"),
                                          srcMetUncdown       = cms.InputTag("patType1CorrectedPFMetUnclusteredEnDown"),
                                          ################################################
                                          srcMVACorrected     = cms.InputTag("patPFMetMVA"),
                                          srcMVAJERup         = cms.InputTag("patPFMetMVAJetResUp"),
                                          srcMVAJERdown       = cms.InputTag("patPFMetMVAJetResDown"),
                                          srcMVAPhoup         = cms.InputTag("patPFMetMVAPhotonEnUp"),
                                          srcMVAPhodown       = cms.InputTag("patPFMetMVAPhotonEnDown"),
                                          srcMVAJetup         = cms.InputTag("patPFMetMVAJetEnUp"),
                                          srcMVAJetdown       = cms.InputTag("patPFMetMVAJetEnDown"),
                                          srcMVAUncup         = cms.InputTag("patPFMetMVAUnclusteredEnUp"),
                                          srcMVAUncdown       = cms.InputTag("patPFMetMVAUnclusteredEnDown"),
                                          

                                                                                    
                                          JetTag            =    cms.untracked.InputTag('ak5PFJetsL1FastL2L3'),
                                          JecTag            =    cms.string("AK5PF"),
                                          GenJetTag         =    cms.untracked.InputTag('ak5GenJets'),
                                          ElectronTag       =    cms.untracked.InputTag('gsfElectrons'),
                                          MuonTag           =    cms.untracked.InputTag('muons'),
                                          PhotonTag         =    cms.untracked.InputTag('photons'),
                                          PrimaryVtxTag     =    cms.untracked.InputTag('offlinePrimaryVertices'),
                                          rhoCorrTag        =    cms.untracked.InputTag('kt6PFJets',    'rho', 'RECO'),
                                          rho25CorrTag      =    cms.untracked.InputTag('kt6PFJetsIso', 'rho', 'PAT'),
                                          rhoMuCorrTag      =    cms.untracked.InputTag('kt6PFJetsCentralNeutral', 'rho','RECO'),  # specifically for muon iso
                                          
                                          partFlowTag       =  cms.untracked.InputTag("particleFlow"), #,"Cleaned"),
                                          
                                          ## MET FILTERS
                                          ecalTPFilterTag    =    cms.untracked.InputTag("EcalDeadCellTriggerPrimitiveFilter",""),
                                          ecalBEFilterTag    =    cms.untracked.InputTag("EcalDeadCellBoundaryEnergyFilter",""),
                                          hcalHBHEFilterTag  =    cms.untracked.InputTag("HBHENoiseFilterResultProducer","HBHENoiseFilterResult"),
                                          hcalLaserFilterTag =    cms.untracked.InputTag("hcalLaserEventFilter",""),
                                          trackingFailureTag =    cms.untracked.InputTag("trackingFailureFilter",""),
                                          eeBadScFilterTag   =    cms.untracked.InputTag("eeBadScFilter",""),
                                          trkPOGFiltersTag1  =    cms.untracked.InputTag("manystripclus53X",""),
                                          trkPOGFiltersTag2  =    cms.untracked.InputTag("toomanystripclus53X",""),
                                          trkPOGFiltersTag3  =    cms.untracked.InputTag("logErrorTooManyClusters",""),
                                          
                                          
                                          skimLepton        =  cms.untracked.bool(False),
                                          
                                          saveMuons         =    cms.untracked.bool(True),
                                          saveJets          =    cms.untracked.bool(True),
                                          saveElectrons     =    cms.untracked.bool(True),
                                          saveEleCrystals   =    cms.untracked.bool(False),
                                          savePhotons       =    cms.untracked.bool(True),
                                          savePhoCrystals   =    cms.untracked.bool(True),
                                          saveMoreEgammaVars=    cms.untracked.bool(True),
                                          saveMET           =    cms.untracked.bool(True),
                                          saveMETExtra      =    cms.untracked.bool(False), 
                                          saveGenJets       =    cms.untracked.bool(True),
                                          saveGenParticles  =    cms.untracked.bool(True),
                                          
                                          #for SC footprint removal
                                          
                                          isolation_cone_size_forSCremoval = cms.untracked.double(0.3),

                                          #for Ecal LazyTools and photon items
                                          ebReducedRecHitCollection = cms.InputTag("reducedEcalRecHitsEB"),
                                          eeReducedRecHitCollection = cms.InputTag("reducedEcalRecHitsEE"),
                                          esReducedRecHitCollection = cms.InputTag("reducedEcalRecHitsES"),
                                          
                                          

  hltName           =    cms.untracked.string("HLT"),
  triggers          =    cms.untracked.vstring(
                                               "HLT_Mu13_Mu8_v",
                                               "HLT_Mu17_Mu8_v",
                                               "HLT_Mu17_TkMu8_v",
                                               "HLT_Mu22_TkMu8_v",
                                               "HLT_Mu22_TkMu22_v",
                                               "HLT_IsoMu24_v",
                                               "HLT_IsoMu24_eta2p1_v",

                                               "HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v",
                                               "HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v",
                                               "HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v",
                                               "HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v",

                                               "HLT_Mu17_Ele8_CaloIdL_v",
                                               "HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v",
                                               "HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v",
                                               "HLT_Mu8_Ele17_CaloIdL_v",
                                               "HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v",
                                               "HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v",

                                               "HLT_Mu22_Photon22_CaloIdL_v",
                                               "HLT_Ele27_WP80_v",
                                               "HLT_Photon36_CaloId10_Iso50_Photon22_CaloId10_Iso50_v",
                                               "HLT_Photon36_CaloId10_Iso50_Photon22_R9Id85_v",
                                               "HLT_Photon36_R9Id85_OR_CaloId10_Iso50_Photon22_R9Id85_OR_CaloId10_Iso50_v",
                                               "HLT_Photon36_R9Id85_Photon22_CaloId10_Iso50_v",
                                               "HLT_Photon36_R9Id85_Photon22_R9Id85_v",

                                               "HLT_Photon30_R9Id90_CaloId_HE10_Iso40_EBOnly_Met25_HBHENoiseCleaned",
                                               "HLT_Photon30_R9Id90_CaloId_HE10_Iso40_EBOnly",
                                               "HLT_Photon30",
                                               "HLT_DiJet20_MJJ650_AllJets_DEta3p5_HT120_VBF",
                                               "HLT_DiJet30_MJJ700_AllJets_DEta3p5_VBF",
                                               "HLT_DiJet35_MJJ650_AllJets_DEta3p5_VBF",
                                               "HLT_DiJet35_MJJ700_AllJets_DEta3p5_VBF",
                                               "HLT_DiJet35_MJJ750_AllJets_DEta3p5_VBF"
)
)

process.ntuplePath = cms.Path(
    #MET STUFF with Syst
    process.patDefaultSequence    
    #* process.puJetIdSqeuence
    * process.goodPhotonsHighPtCut
    * process.pfType1MEtUncertaintySequence
    * process.pfMVAMEtUncertaintySequence
    * AllFilters
    #* recoPuJetIdSqeuence    
    * process.puJetIdSqeuence
    * process.goodOfflinePrimaryVertices
    * process.pfNoPUSeq
    * process.kt6PFJetsIso
    * process.ak5PFJetsL1FastL2L3
    * process.ak5JetTracksAssociatorAtVertex
    * process.btagging
#    * process.eleRegressionEnergy
#    * process.calibratedElectrons
#    * process.mvaTrigV0
    * process.heepIdNoIso
    * process.heepIdNoIsoEles
    * process.modElectronIso
    * process.ntupleProducer

)

#process.out = cms.OutputModule("PoolOutputModule",
#                               fileName = cms.untracked.string('myTuple.root'),
#                               SelectEvents   = cms.untracked.PSet( SelectEvents = cms.vstring('ntuplePath')),
#                               outputCommands = cms.untracked.vstring('keep *')
#                               )
#process.outpath = cms.EndPath(process.out)
