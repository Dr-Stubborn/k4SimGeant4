from Gaudi.Configuration import *
from Configurables import EventDataSvc
from k4FWCore import ApplicationMgr, IOSvc

from GaudiKernel.SystemOfUnits import GeV
from math import pi
from os import environ, path
import uuid

# Electron momentum in GeV
momentum = 50
# Theta and its spread in degrees
theta = 90.
thetaSpread = 10.

# Replace these constants with the values measured with
# fcc_ee_samplingFraction_allegroECal.py when available.
samplingFractions = [1.0] * 11

# IOSvc for output
iosvc = IOSvc("IOSvc")
iosvc.Output = "fccee_energyInCaloLayers_allegroECal_%ideg_%igev_%s.root" % (
    theta, momentum, uuid.uuid4().hex[0:16])
iosvc.outputCommands = ["drop *", "keep energyInLayer", "keep energyInCryo"]

################## Particle gun setup

from Configurables import MomentumRangeParticleGun
pgun = MomentumRangeParticleGun("ParticleGun_Electron")
pgun.PdgCodes = [11]
pgun.MomentumMin = momentum * GeV
pgun.MomentumMax = momentum * GeV
pgun.PhiMin = 0
pgun.PhiMax = 2 * pi
pgun.ThetaMin = (theta - thetaSpread / 2.) * pi / 180.
pgun.ThetaMax = (theta + thetaSpread / 2.) * pi / 180.

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
detector_path = environ.get("FCCDETECTORS", "")
detectors = [
    "FCCee/ALLEGRO/compact/ALLEGRO_o1_v03/DectEmptyMaster.xml",
    "FCCee/ALLEGRO/compact/ALLEGRO_o1_v03/ECalBarrel_thetamodulemerged_upstream.xml",
]
geoservice = GeoSvc("GeoSvc",
                    detectors=[path.join(detector_path, detector) for detector in detectors],
                    OutputLevel=WARNING)

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

from Configurables import CreateCaloCells
createcellsBarrel = CreateCaloCells("CreateCaloCellsBarrel",
                                    doCellCalibration=False,
                                    addPosition=True,
                                    addCellNoise=False,
                                    filterCellNoise=False)
createcellsBarrel.hits.Path = "ECalBarrelHits"
createcellsBarrel.cells.Path = "ECalBarrelCells"

from Configurables import EnergyInCaloLayers
energy_in_layers = EnergyInCaloLayers("energyInLayers",
                                      readoutName="ECalBarrelModuleThetaMerged",
                                      numLayers=11,
                                      samplingFractions=samplingFractions,
                                      OutputLevel=VERBOSE)
energy_in_layers.deposits.Path = "ECalBarrelCells"
energy_in_layers.particle.Path = "GenParticles"

# CPU information
from Configurables import AuditorSvc, ChronoAuditor
chra = ChronoAuditor()
audsvc = AuditorSvc()
audsvc.Auditors = [chra]
geantsim.AuditExecute = True
energy_in_layers.AuditExecute = True

# ApplicationMgr
ApplicationMgr(TopAlg=[genalg_pgun, hepmc_converter, geantsim, createcellsBarrel, energy_in_layers],
               EvtSel="NONE",
               EvtMax=10,
               ExtSvc=[EventDataSvc("EventDataSvc"), geoservice, geantservice, audsvc],
               OutputLevel=DEBUG)
