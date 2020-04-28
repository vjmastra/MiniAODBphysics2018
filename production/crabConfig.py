#crab submit -c crabConfig.py
#from WMCore.Configuration import Configuration
from CRABClient.UserUtilities import Configuration, getUsernameFromSiteDB
config = Configuration()

datasetnames = [
'/ParkingBPH1/Run2018A-05May2019-v1/MINIAOD',
'/ParkingBPH2/Run2018A-05May2019-v1/MINIAOD',
'/ParkingBPH3/Run2018A-05May2019-v1/MINIAOD',
'/ParkingBPH4/Run2018A-05May2019-v1/MINIAOD',
'/ParkingBPH5/Run2018A-05May2019-v1/MINIAOD',
'/ParkingBPH6/Run2018A-05May2019-v1/MINIAOD',
'/ParkingBPH1/Run2018B-05May2019-v1/MINIAOD',
'/ParkingBPH2/Run2018B-05May2019-v1/MINIAOD',
'/ParkingBPH3/Run2018B-05May2019-v1/MINIAOD',
'/ParkingBPH4/Run2018B-05May2019-v1/MINIAOD',
'/ParkingBPH5/Run2018B-05May2019-v1/MINIAOD',
'/ParkingBPH6/Run2018B-05May2019-v1/MINIAOD',
'/ParkingBPH1/Run2018C-05May2019-v1/MINIAOD',
'/ParkingBPH2/Run2018C-05May2019-v1/MINIAOD',
'/ParkingBPH3/Run2018C-05May2019-v1/MINIAOD',
'/ParkingBPH4/Run2018C-05May2019-v1/MINIAOD',
'/ParkingBPH5/Run2018C-05May2019-v1/MINIAOD',
'/ParkingBPH1/Run2018D-05May2019-v1/MINIAOD',
'/ParkingBPH2/Run2018D-05May2019-v1/MINIAOD',
'/ParkingBPH3/Run2018D-05May2019-v1/MINIAOD',
'/ParkingBPH4/Run2018D-05May2019-v1/MINIAOD',
'/ParkingBPH5/Run2018D-05May2019-v1/MINIAOD',
]

psetS = [
'../test/PsilambdaRootupler.py',
'../test/Psi2SlambdaRootupler.py'
]

#runNumber = [
#'',
#'297620,297656',
#'299420'
#]

decays = [
'PsiLambda',
'Psi2SLambda'
]

jsonfile = [
'',
'/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/ReReco/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON_MuonPhys.txt',
'/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/ReReco/Final/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON_MuonPhys.txt',
'/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/ReReco/Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON_MuonPhys.txt'
]

cond = [
'',
'ConstraintY2S',
'ConstraintChib2P'
]

#eventsPerJob = [
#10,
#20
#]

workDir = 'Psi2SLambda_BPH'
decay = decays[1]
pset = psetS[1]
#runNum = runNumber[0]
lumi = jsonfile[3] #no json: 0, 2017json: 1, 2016json: 2, 2018json: 3
#epj = eventsPerJob[0]

condM = cond[0]

datasetName = datasetnames[0]	

print "*****************"
print decay
print pset
print datasetName
print "*****************"

import datetime
timestamp = datetime.datetime.now().strftime("_%Y%m%d_%H%M%S")

dataset = filter(None, datasetName.split('/'))

config.section_('General')
config.General.transferOutputs  = True
config.General.workArea         = '%s' % (workDir)
#config.General.requestName      = dataset[0]+'_'+dataset[1]+'_'+dataset[2]+'_'+runNum+'_'+decay+'_dM1'+timestamp
config.General.requestName      = dataset[0]+'_'+dataset[1]+'_'+dataset[2]+'_'+decay+'_'+condM+timestamp
config.General.transferLogs     = False

config.section_('JobType')
config.JobType.psetName         = pset
config.JobType.pluginName       = 'Analysis'
#config.JobType.maxMemoryMB      = 2500
#config.JobType.maxJobRuntimeMin = 1000
#config.JobType.numCores         = 4
#config.JobType.outputFiles      = ['hltbits.root']
#config.JobType.priority			= 20

config.section_('Data')
config.Data.inputDataset        = datasetName
config.Data.inputDBS            = 'global'
#config.Data.totalUnits          = -1
#config.Data.unitsPerJob         = epj
#config.Data.splitting           = 'LumiBased'
config.Data.splitting           = 'Automatic'
#config.Data.runRange            = runNum
config.Data.lumiMask            = lumi
#config.Data.outLFNDirBase       = '/store/user/vmastrap/%s' % (workDir)
config.Data.outLFNDirBase       = '/store/user/%s/%s' % (getUsernameFromSiteDB(), workDir)
config.Data.publication         = False
#config.Data.ignoreLocality      = True

config.section_('Site')
config.Site.storageSite         = 'T2_IT_Bari'
#config.Site.blacklist           = ['T2_TW_NCHC', 'T2_US_Vanderbilt']
#config.Site.blacklist           = ['T1*', 'T2_BR_SPRACE', 'T2_US_Wisconsin', 'T1_RU_JINR', 'T2_RU_JINR', 'T2_EE_Estonia']
#config.Site.whitelist		= ['T2*']
