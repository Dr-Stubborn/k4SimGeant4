#include "EnergyInECalComponents.h"

#include "GaudiKernel/ITHistSvc.h"
#include "k4Interface/IGeoSvc.h"
#include "TTree.h"

// DD4hep
#include "DD4hep/Detector.h"
#include "DD4hep/Readout.h"

// EDM4hep
#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/SimCalorimeterHitCollection.h"

#include <algorithm>
#include <cmath>
#include <cstddef>

DECLARE_COMPONENT(EnergyInECalComponents)

EnergyInECalComponents::EnergyInECalComponents(const std::string& aName, ISvcLocator* aSvcLoc)
    : Gaudi::Algorithm(aName, aSvcLoc), m_histSvc("THistSvc", aName), m_geoSvc("GeoSvc", aName), m_tree(nullptr) {
  declareProperty("deposits", m_deposits, "Raw ECAL energy deposits (input)");
  declareProperty("particle", m_particle, "Generated single-particle event (input)");
}

EnergyInECalComponents::~EnergyInECalComponents() {}

StatusCode EnergyInECalComponents::initialize() {
  if (Gaudi::Algorithm::initialize().isFailure()) {
    return StatusCode::FAILURE;
  }

  if (!m_geoSvc) {
    error() << "Unable to locate Geometry Service! "
            << "Make sure you have GeoSvc and SimSvc in the right order in the configuration." << endmsg;
    return StatusCode::FAILURE;
  }

  if (m_geoSvc->getDetector()->readouts().find(m_readoutName) == m_geoSvc->getDetector()->readouts().end()) {
    error() << "Can't find readout <<" << m_readoutName << ">>!" << endmsg;
    return StatusCode::FAILURE;
  }

  if (m_numLayers == 0) {
    error() << "numLayers must be larger than zero." << endmsg;
    return StatusCode::FAILURE;
  }

  const auto nLayers = static_cast<std::size_t>(m_numLayers);
  m_activeLArEnergyByLayer.assign(nLayers, 0.);
  m_passiveAbsorberEnergyByLayer.assign(nLayers, 0.);
  m_readoutPcbEnergyByLayer.assign(nLayers, 0.);
  m_calorimeterEnergyByLayer.assign(nLayers, 0.);

  m_tree = new TTree("ecal_component_energy", "Event-level ECAL component energies");
  m_tree->Branch("event", &m_event, "event/I");
  m_tree->Branch("n_hits", &m_nHits, "n_hits/I");
  m_tree->Branch("n_unclassified_hits", &m_nUnclassifiedHits, "n_unclassified_hits/I");

  m_tree->Branch("particle_mass", &m_particleMass, "particle_mass/D");
  m_tree->Branch("particle_px", &m_particlePx, "particle_px/D");
  m_tree->Branch("particle_py", &m_particlePy, "particle_py/D");
  m_tree->Branch("particle_pz", &m_particlePz, "particle_pz/D");
  m_tree->Branch("particle_energy", &m_particleEnergy, "particle_energy/D");
  m_tree->Branch("particle_theta", &m_particleTheta, "particle_theta/D");

  m_tree->Branch("active_lar_energy", &m_activeLArEnergy, "active_lar_energy/D");
  m_tree->Branch("passive_absorber_energy", &m_passiveAbsorberEnergy, "passive_absorber_energy/D");
  m_tree->Branch("readout_pcb_energy", &m_readoutPcbEnergy, "readout_pcb_energy/D");
  m_tree->Branch("front_cryo_energy", &m_frontCryoEnergy, "front_cryo_energy/D");
  m_tree->Branch("back_cryo_energy", &m_backCryoEnergy, "back_cryo_energy/D");
  m_tree->Branch("side_cryo_energy", &m_sideCryoEnergy, "side_cryo_energy/D");
  m_tree->Branch("front_services_energy", &m_frontServicesEnergy, "front_services_energy/D");
  m_tree->Branch("back_services_energy", &m_backServicesEnergy, "back_services_energy/D");

  m_tree->Branch("calorimeter_total_energy", &m_calorimeterTotalEnergy, "calorimeter_total_energy/D");
  m_tree->Branch("cryo_services_total_energy", &m_cryoServicesTotalEnergy, "cryo_services_total_energy/D");
  m_tree->Branch("ecal_total_energy", &m_ecalTotalEnergy, "ecal_total_energy/D");
  m_tree->Branch("unclassified_energy", &m_unclassifiedEnergy, "unclassified_energy/D");
  m_tree->Branch("leakage_energy", &m_leakageEnergy, "leakage_energy/D");

  m_tree->Branch("active_lar_energy_by_layer", &m_activeLArEnergyByLayer);
  m_tree->Branch("passive_absorber_energy_by_layer", &m_passiveAbsorberEnergyByLayer);
  m_tree->Branch("readout_pcb_energy_by_layer", &m_readoutPcbEnergyByLayer);
  m_tree->Branch("calorimeter_energy_by_layer", &m_calorimeterEnergyByLayer);

  if (m_histSvc->regTree(m_treePath, m_tree).isFailure()) {
    error() << "Couldn't register TTree at " << m_treePath << endmsg;
    return StatusCode::FAILURE;
  }

  return StatusCode::SUCCESS;
}

StatusCode EnergyInECalComponents::execute(const EventContext& ctx) const {
  resetEvent();
  m_event = static_cast<int>(ctx.evt());

  const auto particles = m_particle.get();
  if (particles->empty()) {
    error() << "Initial particle not found!" << endmsg;
    return StatusCode::FAILURE;
  }
  if (particles->size() > 1) {
    warning() << "Found more than one initial particle; using the first entry for leakage calculation." << endmsg;
  }

  const auto particle = particles->at(0);
  m_particleMass = particle.getMass();
  m_particlePx = particle.getMomentum().x;
  m_particlePy = particle.getMomentum().y;
  m_particlePz = particle.getMomentum().z;

  const double momentum2 =
      m_particlePx * m_particlePx + m_particlePy * m_particlePy + m_particlePz * m_particlePz;
  m_particleEnergy = std::sqrt(m_particleMass * m_particleMass + momentum2);
  const double rxy = std::sqrt(m_particlePx * m_particlePx + m_particlePy * m_particlePy);
  m_particleTheta = std::atan2(rxy, m_particlePz) * 180. / std::acos(-1.);

  const auto decoder = m_geoSvc->getDetector()->readout(m_readoutName).idSpec().decoder();
  const auto deposits = m_deposits.get();
  m_nHits = static_cast<int>(deposits->size());

  for (const auto& hit : *deposits) {
    const double energy = hit.getEnergy();
    const dd4hep::DDSegmentation::CellID cellID = hit.getCellID();
    const int cryoId = decoder->get(cellID, m_cryoFieldName);
    const int typeId = decoder->get(cellID, m_typeFieldName);

    if (cryoId == 0) {
      const int layerId = decoder->get(cellID, m_layerFieldName);
      if (!addCalorimeterEnergy(typeId, layerId, energy)) {
        addUnclassifiedEnergy(energy);
      }
    } else if (cryoId == 1) {
      if (!addCryoServicesEnergy(typeId, energy)) {
        addUnclassifiedEnergy(energy);
      }
    } else {
      addUnclassifiedEnergy(energy);
    }
  }

  m_calorimeterTotalEnergy = m_activeLArEnergy + m_passiveAbsorberEnergy + m_readoutPcbEnergy;
  m_cryoServicesTotalEnergy =
      m_frontCryoEnergy + m_backCryoEnergy + m_sideCryoEnergy + m_frontServicesEnergy + m_backServicesEnergy;
  m_ecalTotalEnergy = m_calorimeterTotalEnergy + m_cryoServicesTotalEnergy + m_unclassifiedEnergy;
  m_leakageEnergy = m_particleEnergy - m_ecalTotalEnergy;

  verbose() << "ECAL component energies: active LAr = " << m_activeLArEnergy
            << " GeV, passive absorber = " << m_passiveAbsorberEnergy
            << " GeV, readout PCB = " << m_readoutPcbEnergy << " GeV, cryostat/services = "
            << m_cryoServicesTotalEnergy << " GeV, unclassified = " << m_unclassifiedEnergy
            << " GeV, leakage = " << m_leakageEnergy << " GeV" << endmsg;

  m_tree->Fill();
  return StatusCode::SUCCESS;
}

StatusCode EnergyInECalComponents::finalize() { return Gaudi::Algorithm::finalize(); }

void EnergyInECalComponents::resetEvent() const {
  m_nHits = 0;
  m_nUnclassifiedHits = 0;

  m_particleMass = 0.;
  m_particlePx = 0.;
  m_particlePy = 0.;
  m_particlePz = 0.;
  m_particleEnergy = 0.;
  m_particleTheta = 0.;

  m_activeLArEnergy = 0.;
  m_passiveAbsorberEnergy = 0.;
  m_readoutPcbEnergy = 0.;
  m_frontCryoEnergy = 0.;
  m_backCryoEnergy = 0.;
  m_sideCryoEnergy = 0.;
  m_frontServicesEnergy = 0.;
  m_backServicesEnergy = 0.;

  m_calorimeterTotalEnergy = 0.;
  m_cryoServicesTotalEnergy = 0.;
  m_ecalTotalEnergy = 0.;
  m_unclassifiedEnergy = 0.;
  m_leakageEnergy = 0.;

  std::fill(m_activeLArEnergyByLayer.begin(), m_activeLArEnergyByLayer.end(), 0.);
  std::fill(m_passiveAbsorberEnergyByLayer.begin(), m_passiveAbsorberEnergyByLayer.end(), 0.);
  std::fill(m_readoutPcbEnergyByLayer.begin(), m_readoutPcbEnergyByLayer.end(), 0.);
  std::fill(m_calorimeterEnergyByLayer.begin(), m_calorimeterEnergyByLayer.end(), 0.);
}

bool EnergyInECalComponents::addCalorimeterEnergy(int typeId, int layerId, double energy) const {
  if (layerId < 0 || layerId >= static_cast<int>(m_numLayers)) {
    return false;
  }

  const auto layerIndex = static_cast<std::size_t>(layerId);
  switch (typeId) {
  case 0:
    m_activeLArEnergy += energy;
    m_activeLArEnergyByLayer[layerIndex] += energy;
    break;
  case 1:
    m_passiveAbsorberEnergy += energy;
    m_passiveAbsorberEnergyByLayer[layerIndex] += energy;
    break;
  case 2:
    m_readoutPcbEnergy += energy;
    m_readoutPcbEnergyByLayer[layerIndex] += energy;
    break;
  default:
    return false;
  }

  m_calorimeterEnergyByLayer[layerIndex] += energy;
  return true;
}

bool EnergyInECalComponents::addCryoServicesEnergy(int typeId, double energy) const {
  switch (typeId) {
  case 1:
    m_frontCryoEnergy += energy;
    break;
  case 2:
    m_backCryoEnergy += energy;
    break;
  case 3:
    m_sideCryoEnergy += energy;
    break;
  case 4:
    m_frontServicesEnergy += energy;
    break;
  case 5:
    m_backServicesEnergy += energy;
    break;
  default:
    return false;
  }
  return true;
}

void EnergyInECalComponents::addUnclassifiedEnergy(double energy) const {
  m_unclassifiedEnergy += energy;
  ++m_nUnclassifiedHits;
}
