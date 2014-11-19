import FWCore.ParameterSet.Types  as CfgTypes
import FWCore.ParameterSet.Config as cms

process = cms.Process("MEAnalysisNew")

process.fwliteInput = cms.PSet(
    samples = cms.VPSet([
        cms.PSet(
                 fileName=cms.string("/Users/joosep/Documents/tth/step1/v5/tthbb.root"),
                 nickName=cms.string("tthbb_13TeV"),
                 type=cms.int32(3),
                 process=cms.int32(0)
        ),
#        cms.PSet(
#                 fileName=cms.string("/Users/joosep/Documents/tth/step1/nov18/ttjets.root"),
#                 nickName=cms.string("ttjets_13TeV"),
#                 type=cms.int32(3),
#                 process=cms.int32(1)
#        ),
        ]),
        evLimits=cms.vint32(0, -1)
)
