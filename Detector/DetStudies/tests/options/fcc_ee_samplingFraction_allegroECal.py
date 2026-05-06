from Gaudi.Configuration import *
from Configurables import EventDataSvc
from k4FWCore import ApplicationMgr, IOSvc

from GaudiKernel.SystemOfUnits import GeV
from math import pi
from os import environ, path

# Electron momentum in GeV
momentum = 50
# Theta min and max in degrees
thetaMin = 90.
thetaMax = 90.

# IOSvc for output
iosvc = IOSvc("IOSvc")
iosvc.Output = "fccee_samplingFraction_allegroECal.root"
iosvc.outputCommands = ["keep *"]

################## Particle gun setup

from Configurables import MomentumRangeParticleGun
pgun = MomentumRangeParticleGun("ParticleGun_Electron")
pgun.PdgCodes = [11]
pgun.MomentumMin = momentum * GeV
pgun.MomentumMax = momentum * GeV
pgun.PhiMin = 0
pgun.PhiMax = 2 * pi
pgun.ThetaMin = thetaMin * pi / 180.
pgun.ThetaMax = thetaMax * pi / 180.

from Configurables import GenAlg
genalg_pgun = GenAlg()
genalg_pgun.SignalProvider = pgun
genalg_pgun.hepmc.Path = "hepmc"

from Configurables import HepMCToEDMConverter
hepmc_converter = HepMCToEDMConverter()
hepmc_converter.hepmc.Path = "hepmc"
hepmc_converter.GenParticles.Path = "GenParticles"

# DD4hep geometry service
from Configurables import GeoSvc
detector_path = environ.get("K4GEO", "")
detectors = [
    "FCCee/ALLEGRO/compact/ALLEGRO_o1_v03/DectEmptyMaster.xml",
    "FCCee/ALLEGRO/compact/ALLEGRO_o1_v03/ECalBarrel_thetamodulemerged_calibration.xml",
]
geoservice = GeoSvc("GeoSvc",
                    detectors=[path.join(detector_path, detector) for detector in detectors],
                    OutputLevel=INFO)

# Geant4 service
from Configurables import SimG4Svc
geantservice = SimG4Svc("SimG4Svc",
                        detector="SimG4DD4hepDetector",
                        physicslist="SimG4FtfpBert",
                        actions="SimG4FullSimActions")
geantservice.g4PostInitCommands += ["/run/setCut 0.1 mm"]

# Geant4 algorithm
from Configurables import SimG4Alg, SimG4SaveCalHits
saveecaltool = SimG4SaveCalHits("saveECalBarrelHits",
                                readoutNames=["ECalBarrelModuleThetaMerged"])
saveecaltool.CaloHits.Path = "ECalBarrelHits"

from Configurables import SimG4PrimariesFromEdmTool
particle_converter = SimG4PrimariesFromEdmTool("EdmConverter")
particle_converter.GenParticles.Path = "GenParticles"

geantsim = SimG4Alg("SimG4Alg",
                    outputs=["SimG4SaveCalHits/saveECalBarrelHits"],
                    eventProvider=particle_converter,
                    OutputLevel=DEBUG)

from Configurables import SamplingFractionInLayers
hist = SamplingFractionInLayers("hists",
                                energyAxis=momentum,
                                readoutName="ECalBarrelModuleThetaMerged",
                                layerFieldName="layer",
                                activeFieldName="type",
                                activeFieldValue=0,
                                numLayers=11,
                                OutputLevel=INFO)
hist.deposits.Path = "ECalBarrelHits"

THistSvc().Output = ["rec DATAFILE='histSF_fccee_allegroECal.root' TYP='ROOT' OPT='RECREATE'"]
THistSvc().PrintAll = True
THistSvc().AutoSave = True
THistSvc().AutoFlush = False
THistSvc().OutputLevel = INFO

# CPU information
from Configurables import AuditorSvc, ChronoAuditor
chra = ChronoAuditor()
audsvc = AuditorSvc()
audsvc.Auditors = [chra]
geantsim.AuditExecute = True
hist.AuditExecute = True

# ApplicationMgr
ApplicationMgr(TopAlg=[genalg_pgun, hepmc_converter, geantsim, hist],
               EvtSel="NONE",
               EvtMax=10,
               ExtSvc=[EventDataSvc("EventDataSvc"), geoservice, geantservice, audsvc],
               OutputLevel=DEBUG)
