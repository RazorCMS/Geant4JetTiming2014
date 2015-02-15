# Auto generated configuration file
# using: 
# Revision: 1.381.2.7 
# Source: /local/reps/CMSSW/CMSSW/Configuration/PyReleaseValidation/python/ConfigBuilder.py,v 
# with command line options: Configuration/GenProduction/python/EightTeV/DiKaon_E_1to100_gun_cff.py -s GEN,SIM,DIGI,L1,DIGI2RAW,HLT --conditions auto:mc --datatier GEN-SIM-RAW --eventcontent RAWSIM -n 10 --no_exec
import FWCore.ParameterSet.Config as cms

process = cms.Process('HLT')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.GeometrySimDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Generator_cff')

#process.load('IOMC.EventVertexGenerators.VtxSmearedRealistic8TeVCollision_cfi')
#process.load('IOMC.EventVertexGenerators.VtxSmearedGauss_cfi')
process.load('IOMC.EventVertexGenerators.VtxSmearedBetafuncNominalCollision_cfi')

#process.VtxSmeared.SigmaY = 0
#process.VtxSmeared.SigmaX = 0
#process.VtxSmeared.SigmaZ = 0


process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.SimIdeal_cff')
process.load('Configuration.StandardSequences.Digi_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('Configuration.StandardSequences.DigiToRaw_cff')
process.load('HLTrigger.Configuration.HLT_GRun_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.RandomNumberGeneratorService.generator.initialSeed = 309936
process.RandomNumberGeneratorService.g4SimHits.initialSeed = 3502690

process.maxEvents = cms.untracked.PSet(
	input = cms.untracked.int32(3000)
)

# Input source
process.source = cms.Source("EmptySource")

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('\\$Revision: 1.3 $'),
    annotation = cms.untracked.string('Hgammagamma (gg->H), no tune'),
    name = cms.untracked.string('\\$Source: /local/reps/CMSSW/UserCode/yangyong/ECALG4SIM/test_xtal1.py,v $')
)

# Output definition

process.RAWSIMoutput = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    outputCommands = process.RAWSIMEventContent.outputCommands,
    fileName = cms.untracked.string('Hgammagamma_DIGI2RAW_HLT.root'),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('GEN-SIM-RAW')
    ),
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('generation_step')
    )
)

# Additional output definition

# Other statements
process.genstepfilter.triggerConditions=cms.vstring("generation_step")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:mc', '')

process.generator = cms.EDFilter("Pythia6GeneratorFilter",
				 comEnergy = cms.double(14000.0),
				 crossSection = cms.untracked.double(1.0),
				 filterEfficiency = cms.untracked.double(1.0),
				 maxEventsToPrint = cms.untracked.int32(0),
				 pythiaHepMCVerbosity = cms.untracked.bool(False),
				 pythiaPylistVerbosity = cms.untracked.int32(0),
				 PythiaParameters = cms.PSet( processParameters = cms.vstring(
	'PMAS(25,1)=125.0        !mass of Higgs',
	'MSEL=0                  ! user selection for process',
	'MSUB(102)=0             !ggH',
	'MSUB(123)=1             !ZZ fusion to H',
	'MSUB(124)=1             !WW fusion to H',
	'MDME(190,1) = 0            !W decay into dbar u',
	'MDME(191,1) = 0            !W decay into dbar c',
	'MDME(192,1) = 0            !W decay into dbar t',
	'MDME(194,1) = 0            !W decay into sbar u',
	'MDME(195,1) = 0            !W decay into sbar c',
	'MDME(196,1) = 0            !W decay into sbar t',
	'MDME(198,1) = 0            !W decay into bbar u',
	'MDME(199,1) = 0            !W decay into bbar c',
	'MDME(200,1) = 0            !W decay into bbar t',
	'MDME(206,1) = 0            !W decay into e+ nu_e',
	'MDME(207,1) = 0            !W decay into mu+ nu_mu',
	'MDME(208,1) = 0            !W decay into tau+ nu_tau',
	'MDME(210,1)=0           !Higgs decay into dd',
	'MDME(211,1)=0           !Higgs decay into uu',
	'MDME(212,1)=0           !Higgs decay into ss',
	'MDME(213,1)=0           !Higgs decay into cc',
	'MDME(214,1)=0           !Higgs decay into bb',
	'MDME(215,1)=0           !Higgs decay into tt',
	'MDME(216,1)=0           !Higgs decay into',
	'MDME(217,1)=0           !Higgs decay into Higgs decay',
	'MDME(218,1)=0           !Higgs decay into e nu e',
	'MDME(219,1)=0           !Higgs decay into mu nu mu',
	'MDME(220,1)=0           !Higgs decay into tau nu tau',
	'MDME(221,1)=0           !Higgs decay into Higgs decay',
	'MDME(222,1)=0           !Higgs decay into g g',
	'MDME(223,1)=1           !Higgs decay into gam gam',
	'MDME(224,1)=0           !Higgs decay into gam Z',
	'MDME(225,1)=0           !Higgs decay into Z Z',
	'MDME(226,1)=0          !Higgs decay into W W'),
        parameterSets = cms.vstring('processParameters')
							      )
 )


# Path and EndPath definitions
process.generation_step = cms.Path(process.pgen)
process.simulation_step = cms.Path(process.psim)
process.digitisation_step = cms.Path(process.pdigi)
process.L1simulation_step = cms.Path(process.SimL1Emulator)
process.digi2raw_step = cms.Path(process.DigiToRaw)
process.genfiltersummary_step = cms.EndPath(process.genFilterSummary)
#process.endjob_step = cms.EndPath(process.endOfProcess)
#process.RAWSIMoutput_step = cms.EndPath(process.RAWSIMoutput)


process.recoAnalyzer = cms.EDAnalyzer("RecoAnalyzer",
				      outputFile = cms.string('analysis_Hgg.root'),
				      muons = cms.untracked.InputTag("boostedMuonsStep2"),
				      genParticles = cms.untracked.InputTag("genParticles"),
				      electron = cms.untracked.InputTag("boostedElectronsStep2"),
				      photon = cms.untracked.InputTag("fsrPhotonsStep2"),
				      debugLevel = cms.int32(0)
				      )


#For Transient Track Builder
process.load('TrackingTools/TransientTrack/TransientTrackBuilder_cfi')

#For BTagging Discriminators
process.load('CMSAna.CMSNtupler.btagging_cff')

#For MC jet matching
process.load('CMSAna.CMSNtupler.jetflavorMatching_cff')
process.load("CMSAna.CMSNtupler.CMSNtupler_cfi")

process.myCMSNtupler.debugLevel = cms.int32(0)
process.myCMSNtupler.FillEGRegressionVars = cms.bool(True)
process.myCMSNtupler.FillGenOnly = cms.bool(True)
process.myCMSNtupler.FillAllGenParticles = cms.bool(True)
process.myCMSNtupler.GenJetPtMin = cms.double(15.0)
process.myCMSNtupler.JetPtMin = cms.double(20.0)
process.myCMSNtupler.outputFile = cms.string('BACONNtuple_VBFHgg.root')

#Define ntupler sequence
process.ntupler_sequence = cms.Sequence(
    process.myGenJetFlavourId*
    process.myCMSNtupler
    )
process.ntupler_step  = cms.Path(process.ntupler_sequence)



# Schedule definition
process.schedule = cms.Schedule(process.generation_step,process.genfiltersummary_step,process.simulation_step,process.digitisation_step,process.L1simulation_step,process.digi2raw_step,process.ntupler_step)

#process.out_step = cms.EndPath(process.recoAnalyzer )
#process.schedule.extend([process.out_step])


#process.schedule.extend(process.HLTSchedule)
#process.schedule.extend([process.endjob_step,process.RAWSIMoutput_step])


# filter all path with the production filter sequence
for path in process.paths:
	getattr(process,path)._seq = process.generator * getattr(process,path)._seq 

# customisation of the process.

# Automatic addition of the customisation function from HLTrigger.Configuration.customizeHLTforMC
from HLTrigger.Configuration.customizeHLTforMC import customizeHLTforMC 

#call to customisation function customizeHLTforMC imported from HLTrigger.Configuration.customizeHLTforMC
process = customizeHLTforMC(process)

# End of customisation functions
