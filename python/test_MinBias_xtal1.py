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

process.RandomNumberGeneratorService.generator.initialSeed = 589034996

process.maxEvents = cms.untracked.PSet(
	input = cms.untracked.int32(100)
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
				#fileName = cms.untracked.string('Hgammagamma_DIGI2RAW_HLT.root'),
				fileName = cms.untracked.string('MinBias_DIGI2RAW_HLT.root'),
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
				 comEnergy = cms.double(8000.0),
				 pythiaPylistVerbosity = cms.untracked.int32(0),
				 # put here the efficiency of your filter (1. if no filter) dummy
				 filterEfficiency = cms.untracked.double(1.0),
				 pythiaHepMCVerbosity = cms.untracked.bool(False),
				 # put here the cross section of your process (in pb) dummy
				 crossSection = cms.untracked.double(55000000000.0),
				 maxEventsToPrint = cms.untracked.int32(0),
				 PythiaParameters = cms.PSet(
	pythiaUESettings = cms.vstring('MSTJ(11)=3     ! Choice of the fragmentation function',
				       'MSTJ(22)=2     ! Decay those unstable particles',
				       'PARJ(71)=10 .  ! for which ctau  10 mm',
				       'MSTP(2)=1      ! which order running alphaS',
				       'MSTP(33)=0     ! no K factors in hard cross sections',
				       'MSTP(51)=10042     ! CTEQ6L1 structure function chosen',
				       'MSTP(52)=2     ! work with LHAPDF',
				       'MSTP(81)=1     ! multiple parton interactions 1 is Pythia default',
				       'MSTP(82)=4     ! Defines the multi-parton model',
				       'MSTU(21)=1     ! Check on possible errors during program execution',
				       'PARP(82)=1.8387   ! pt cutoff for multiparton interactions',
				       'PARP(89)=1960. ! sqrts for which PARP82 is set',
				       'PARP(83)=0.5   ! Multiple interactions: matter distrbn parameter',
				       'PARP(84)=0.4   ! Multiple interactions: matter distribution parameter',
				       'PARP(90)=0.16  ! Multiple interactions: rescaling power',
				       'PARP(67)=2.5    ! amount of initial-state radiation',
				       'PARP(85)=1.0  ! gluon prod. mechanism in MI',
				       'PARP(86)=1.0  ! gluon prod. mechanism in MI',
				       'PARP(62)=1.25   ! ',
				       'PARP(64)=0.2    ! ',
				       'MSTP(91)=1     !',
				       'PARP(91)=2.1   ! kt distribution',
				       'PARP(93)=15.0  ! '),
	# This is a vector of ParameterSet names to be read, in this order
	# The first two are in the include files below
	# The last one are simply my additional parameters
	parameterSets = cms.vstring('pythiaUESettings',
				    'pythiaMinBias'),
	pythiaMinBias = cms.vstring('MSEL=0         ! User defined processes',
				    'MSUB(11)=1     ! Min bias process',
				    'MSUB(12)=1     ! Min bias process',
				    'MSUB(13)=1     ! Min bias process',
				    'MSUB(28)=1     ! Min bias process',
				    'MSUB(53)=1     ! Min bias process',
				    'MSUB(68)=1     ! Min bias process',
				    'MSUB(92)=1     ! Min bias process, single diffractive',
				    'MSUB(93)=1     ! Min bias process, single diffractive',
				    'MSUB(94)=1     ! Min bias process, double diffractive',
				    'MSUB(95)=1     ! Min bias process')
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
				      outputFile = cms.string('analysis_MinBias.root'),
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
process.myCMSNtupler.outputFile = cms.string('BACONNtuple_MinBias.root')

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
