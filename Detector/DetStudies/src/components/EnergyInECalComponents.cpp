#include "EnergyInECalComponents.h"

#include "GaudiKernel/ITHistSvc.h"
#include "TTree.h"
#include "k4Interface/IGeoSvc.h"

// DD4hep
#include "DD4hep/Detector.h"
#include "DD4hep/Readout.h"

// EDM4hep
#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/SimCalorimeterHitCollection.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace {

constexpr std::array<double, 5> kContainmentFractions{0.50, 0.68, 0.90, 0.95, 0.99};
constexpr std::array<std::size_t, 3> kSummaryRadiusIndices{2, 3, 4};

struct Vector3 {
  double x;
  double y;
  double z;
};

struct Axis {
  Vector3 origin;
  Vector3 unitDirection;
};

struct EnergyDeposit {
  Vector3 position;
  double energy;
};

struct SampleStatistics {
  double mean;
  double median;
  double energyWeightedMean;
  double totalWeight;
};

double invalidValue() { return std::numeric_limits<double>::quiet_NaN(); }

bool isFinite(const Vector3& vector) {
  return std::isfinite(vector.x) && std::isfinite(vector.y) && std::isfinite(vector.z);
}

std::optional<Axis> makeAxis(const Vector3& origin, const Vector3& direction) {
  if (!isFinite(origin) || !isFinite(direction)) {
    return std::nullopt;
  }

  const double norm = std::hypot(direction.x, direction.y, direction.z);
  if (!std::isfinite(norm) || norm <= 0.) {
    return std::nullopt;
  }

  return Axis{origin, {direction.x / norm, direction.y / norm, direction.z / norm}};
}

double transverseDistance(const Vector3& position, const Axis& axis) {
  const Vector3 displacement{position.x - axis.origin.x, position.y - axis.origin.y, position.z - axis.origin.z};
  const Vector3 cross{displacement.y * axis.unitDirection.z - displacement.z * axis.unitDirection.y,
                      displacement.z * axis.unitDirection.x - displacement.x * axis.unitDirection.z,
                      displacement.x * axis.unitDirection.y - displacement.y * axis.unitDirection.x};
  return std::hypot(cross.x, cross.y, cross.z);
}

std::optional<std::array<double, 5>> containmentRadii(const std::vector<EnergyDeposit>& deposits, const Axis& axis,
                                                      double normalizationEnergy) {
  if (!std::isfinite(normalizationEnergy) || normalizationEnergy <= 0.) {
    return std::nullopt;
  }

  std::vector<std::pair<double, double>> radialEnergies;
  radialEnergies.reserve(deposits.size());
  for (const auto& deposit : deposits) {
    if (deposit.energy <= 0.) {
      continue;
    }
    const double radius = transverseDistance(deposit.position, axis);
    if (!std::isfinite(radius)) {
      return std::nullopt;
    }
    radialEnergies.emplace_back(radius, deposit.energy);
  }

  if (radialEnergies.empty()) {
    return std::nullopt;
  }

  std::sort(radialEnergies.begin(), radialEnergies.end(),
            [](const auto& left, const auto& right) { return left.first < right.first; });

  std::array<double, 5> radii;
  radii.fill(invalidValue());
  double cumulativeEnergy = 0.;
  std::size_t fractionIndex = 0;
  for (const auto& [radius, energy] : radialEnergies) {
    cumulativeEnergy += energy;
    while (fractionIndex < kContainmentFractions.size() &&
           cumulativeEnergy >= kContainmentFractions[fractionIndex] * normalizationEnergy) {
      radii[fractionIndex] = radius;
      ++fractionIndex;
    }
    if (fractionIndex == kContainmentFractions.size()) {
      break;
    }
  }

  if (fractionIndex != kContainmentFractions.size()) {
    return std::nullopt;
  }
  return radii;
}

bool energiesConsistent(double selectedEnergy, double normalizationEnergy) {
  const double scale = std::max({1., std::abs(selectedEnergy), std::abs(normalizationEnergy)});
  return std::abs(selectedEnergy - normalizationEnergy) <= 1.e-9 * scale;
}

SampleStatistics summarizeSamples(const std::vector<double>& radii, const std::vector<double>& energyWeights) {
  const SampleStatistics invalid{invalidValue(), invalidValue(), invalidValue(), 0.};
  if (radii.empty() || radii.size() != energyWeights.size()) {
    return invalid;
  }

  long double radiusSum = 0.;
  long double weightedRadiusSum = 0.;
  long double weightSum = 0.;
  for (std::size_t index = 0; index < radii.size(); ++index) {
    if (!std::isfinite(radii[index]) || !std::isfinite(energyWeights[index]) || energyWeights[index] <= 0.) {
      return invalid;
    }
    radiusSum += radii[index];
    weightedRadiusSum += static_cast<long double>(radii[index]) * energyWeights[index];
    weightSum += energyWeights[index];
  }

  std::vector<double> sortedRadii = radii;
  std::sort(sortedRadii.begin(), sortedRadii.end());
  const std::size_t middle = sortedRadii.size() / 2;
  const double median =
      sortedRadii.size() % 2 == 0 ? 0.5 * (sortedRadii[middle - 1] + sortedRadii[middle]) : sortedRadii[middle];

  return {static_cast<double>(radiusSum / radii.size()), median, static_cast<double>(weightedRadiusSum / weightSum),
          static_cast<double>(weightSum)};
}

} // namespace

DECLARE_COMPONENT(EnergyInECalComponents)

EnergyInECalComponents::EnergyInECalComponents(const std::string& aName, ISvcLocator* aSvcLoc)
    : Gaudi::Algorithm(aName, aSvcLoc), m_histSvc("THistSvc", aName), m_geoSvc("GeoSvc", aName), m_tree(nullptr),
      m_moliereSummaryTree(nullptr) {
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

  if (m_centroidAxisOrigin.value().size() != 3) {
    error() << "centroidAxisOrigin must contain exactly three coordinates (x, y, z)." << endmsg;
    return StatusCode::FAILURE;
  }
  const Vector3 centroidAxisOrigin{m_centroidAxisOrigin.value()[0], m_centroidAxisOrigin.value()[1],
                                   m_centroidAxisOrigin.value()[2]};
  if (!isFinite(centroidAxisOrigin)) {
    error() << "centroidAxisOrigin coordinates must be finite." << endmsg;
    return StatusCode::FAILURE;
  }
  m_centroidAxisOriginX = centroidAxisOrigin.x;
  m_centroidAxisOriginY = centroidAxisOrigin.y;
  m_centroidAxisOriginZ = centroidAxisOrigin.z;
  m_summaryUsesFrontMaterial = m_moliereIncludeFrontMaterial.value();

  const auto nLayers = static_cast<std::size_t>(m_numLayers);
  m_activeLArEnergyByLayer.assign(nLayers, 0.);
  m_passiveAbsorberEnergyByLayer.assign(nLayers, 0.);
  m_readoutPcbEnergyByLayer.assign(nLayers, 0.);
  m_calorimeterEnergyByLayer.assign(nLayers, 0.);

  m_tree = new TTree("ecal_component_energy", "Event-level ECAL component energies");
  m_tree->Branch("event", &m_event, "event/I");
  m_tree->Branch("n_hits", &m_nHits, "n_hits/I");
  m_tree->Branch("n_unclassified_hits", &m_nUnclassifiedHits, "n_unclassified_hits/I");

  m_tree->Branch("incident_particle_index", &m_incidentParticleIndex, "incident_particle_index/I");
  m_tree->Branch("incident_particle_pdg", &m_incidentParticlePdg, "incident_particle_pdg/I");
  m_tree->Branch("n_incident_particle_candidates", &m_nIncidentParticleCandidates, "n_incident_particle_candidates/I");
  m_tree->Branch("particle_mass", &m_particleMass, "particle_mass/D");
  m_tree->Branch("particle_px", &m_particlePx, "particle_px/D");
  m_tree->Branch("particle_py", &m_particlePy, "particle_py/D");
  m_tree->Branch("particle_pz", &m_particlePz, "particle_pz/D");
  m_tree->Branch("particle_energy", &m_particleEnergy, "particle_energy/D");
  m_tree->Branch("particle_theta", &m_particleTheta, "particle_theta/D");
  m_tree->Branch("particle_vertex_x_mm", &m_particleVertexX, "particle_vertex_x_mm/D");
  m_tree->Branch("particle_vertex_y_mm", &m_particleVertexY, "particle_vertex_y_mm/D");
  m_tree->Branch("particle_vertex_z_mm", &m_particleVertexZ, "particle_vertex_z_mm/D");

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

  m_tree->Branch("moliere_include_front_material", &m_moliereUsesFrontMaterial, "moliere_include_front_material/O");
  m_tree->Branch("moliere_n_hits", &m_nMoliereHits, "moliere_n_hits/I");
  m_tree->Branch("moliere_normalization_energy", &m_moliereNormalizationEnergy, "moliere_normalization_energy/D");
  m_tree->Branch("moliere_calorimeter_centroid_x_mm", &m_calorimeterCentroidX, "moliere_calorimeter_centroid_x_mm/D");
  m_tree->Branch("moliere_calorimeter_centroid_y_mm", &m_calorimeterCentroidY, "moliere_calorimeter_centroid_y_mm/D");
  m_tree->Branch("moliere_calorimeter_centroid_z_mm", &m_calorimeterCentroidZ, "moliere_calorimeter_centroid_z_mm/D");
  m_tree->Branch("moliere_truth_axis_valid", &m_truthAxisValid, "moliere_truth_axis_valid/O");
  m_tree->Branch("moliere_centroid_axis_valid", &m_centroidAxisValid, "moliere_centroid_axis_valid/O");

  m_tree->Branch("moliere_truth_r50_mm", &m_truthRadii[0], "moliere_truth_r50_mm/D");
  m_tree->Branch("moliere_truth_r68_mm", &m_truthRadii[1], "moliere_truth_r68_mm/D");
  m_tree->Branch("moliere_truth_r90_mm", &m_truthRadii[2], "moliere_truth_r90_mm/D");
  m_tree->Branch("moliere_truth_r95_mm", &m_truthRadii[3], "moliere_truth_r95_mm/D");
  m_tree->Branch("moliere_truth_r99_mm", &m_truthRadii[4], "moliere_truth_r99_mm/D");
  m_tree->Branch("moliere_centroid_r50_mm", &m_centroidRadii[0], "moliere_centroid_r50_mm/D");
  m_tree->Branch("moliere_centroid_r68_mm", &m_centroidRadii[1], "moliere_centroid_r68_mm/D");
  m_tree->Branch("moliere_centroid_r90_mm", &m_centroidRadii[2], "moliere_centroid_r90_mm/D");
  m_tree->Branch("moliere_centroid_r95_mm", &m_centroidRadii[3], "moliere_centroid_r95_mm/D");
  m_tree->Branch("moliere_centroid_r99_mm", &m_centroidRadii[4], "moliere_centroid_r99_mm/D");

  if (m_histSvc->regTree(m_treePath, m_tree).isFailure()) {
    error() << "Couldn't register TTree at " << m_treePath << endmsg;
    return StatusCode::FAILURE;
  }

  m_moliereSummaryTree = new TTree("ecal_moliere_summary", "Run-level Moliere radius summary");
  m_moliereSummaryTree->Branch("n_events_total", &m_summaryNEventsTotal, "n_events_total/L");
  m_moliereSummaryTree->Branch("moliere_include_front_material", &m_summaryUsesFrontMaterial,
                               "moliere_include_front_material/O");
  m_moliereSummaryTree->Branch("centroid_axis_origin_x_mm", &m_centroidAxisOriginX, "centroid_axis_origin_x_mm/D");
  m_moliereSummaryTree->Branch("centroid_axis_origin_y_mm", &m_centroidAxisOriginY, "centroid_axis_origin_y_mm/D");
  m_moliereSummaryTree->Branch("centroid_axis_origin_z_mm", &m_centroidAxisOriginZ, "centroid_axis_origin_z_mm/D");
  m_moliereSummaryTree->Branch("truth_axis_n_valid_events", &m_summaryTruthNValidEvents, "truth_axis_n_valid_events/L");
  m_moliereSummaryTree->Branch("centroid_axis_n_valid_events", &m_summaryCentroidNValidEvents,
                               "centroid_axis_n_valid_events/L");
  m_moliereSummaryTree->Branch("truth_axis_total_weight_energy_gev", &m_summaryTruthTotalWeightEnergy,
                               "truth_axis_total_weight_energy_gev/D");
  m_moliereSummaryTree->Branch("centroid_axis_total_weight_energy_gev", &m_summaryCentroidTotalWeightEnergy,
                               "centroid_axis_total_weight_energy_gev/D");

  const std::array<std::string, 3> summaryFractions{"90", "95", "99"};
  const auto addSummaryBranches =
      [this, &summaryFractions](const std::string& axisName, std::array<double, 3>& eventMean,
                                std::array<double, 3>& eventMedian, std::array<double, 3>& eventEnergyWeightedMean) {
        for (std::size_t index = 0; index < summaryFractions.size(); ++index) {
          const std::string prefix = axisName + "_r" + summaryFractions[index];
          const std::string meanName = prefix + "_event_mean_mm";
          const std::string medianName = prefix + "_event_median_mm";
          const std::string weightedMeanName = prefix + "_event_energy_weighted_mean_mm";
          m_moliereSummaryTree->Branch(meanName.c_str(), &eventMean[index], (meanName + "/D").c_str());
          m_moliereSummaryTree->Branch(medianName.c_str(), &eventMedian[index], (medianName + "/D").c_str());
          m_moliereSummaryTree->Branch(weightedMeanName.c_str(), &eventEnergyWeightedMean[index],
                                       (weightedMeanName + "/D").c_str());
        }
      };
  addSummaryBranches("truth_axis", m_summaryTruthEventMean, m_summaryTruthEventMedian,
                     m_summaryTruthEventEnergyWeightedMean);
  addSummaryBranches("centroid_axis", m_summaryCentroidEventMean, m_summaryCentroidEventMedian,
                     m_summaryCentroidEventEnergyWeightedMean);

  if (m_histSvc->regTree(m_moliereSummaryTreePath, m_moliereSummaryTree).isFailure()) {
    error() << "Couldn't register Moliere summary TTree at " << m_moliereSummaryTreePath << endmsg;
    return StatusCode::FAILURE;
  }

  return StatusCode::SUCCESS;
}

StatusCode EnergyInECalComponents::execute(const EventContext& ctx) const {
  const std::lock_guard<std::mutex> stateLock(m_stateMutex);
  resetEvent();
  m_event = static_cast<int>(ctx.evt());

  const auto particles = m_particle.get();
  double highestParticleEnergy = -std::numeric_limits<double>::infinity();
  bool highestEnergyTie = false;
  for (std::size_t index = 0; index < particles->size(); ++index) {
    const auto candidate = particles->at(index);
    if (candidate.getGeneratorStatus() != 1 || candidate.isCreatedInSimulation()) {
      continue;
    }

    ++m_nIncidentParticleCandidates;
    const double candidateEnergy = candidate.getEnergy();
    if (!std::isfinite(candidateEnergy)) {
      warning() << "Ignoring non-finite incident-particle candidate at collection index " << index << endmsg;
      continue;
    }

    if (m_incidentParticleIndex < 0 || candidateEnergy > highestParticleEnergy) {
      m_incidentParticleIndex = static_cast<int>(index);
      highestParticleEnergy = candidateEnergy;
      highestEnergyTie = false;
    } else if (candidateEnergy == highestParticleEnergy) {
      highestEnergyTie = true;
    }
  }

  std::optional<Axis> truthAxis;
  if (m_incidentParticleIndex >= 0) {
    if (highestEnergyTie) {
      warning() << "Multiple incident particles share the highest energy; using the lowest collection index "
                << m_incidentParticleIndex << endmsg;
    }

    const auto particle = particles->at(static_cast<std::size_t>(m_incidentParticleIndex));
    const auto momentum = particle.getMomentum();
    const auto vertex = particle.getVertex();
    m_incidentParticlePdg = particle.getPDG();
    m_particleMass = particle.getMass();
    m_particlePx = momentum.x;
    m_particlePy = momentum.y;
    m_particlePz = momentum.z;
    m_particleEnergy = highestParticleEnergy;
    m_particleVertexX = vertex.x;
    m_particleVertexY = vertex.y;
    m_particleVertexZ = vertex.z;

    const double transverseMomentum = std::hypot(m_particlePx, m_particlePy);
    const double momentumNorm = std::hypot(transverseMomentum, m_particlePz);
    if (std::isfinite(momentumNorm) && momentumNorm > 0.) {
      m_particleTheta = std::atan2(transverseMomentum, m_particlePz) * 180. / std::acos(-1.);
    }
    truthAxis =
        makeAxis({m_particleVertexX, m_particleVertexY, m_particleVertexZ}, {m_particlePx, m_particlePy, m_particlePz});
  } else {
    warning() << "No finite generator-level incident particle found; truth-axis Moliere radii are invalid for event "
              << m_event << endmsg;
  }

  const auto decoder = m_geoSvc->getDetector()->readout(m_readoutName).idSpec().decoder();
  const auto deposits = m_deposits.get();
  m_nHits = static_cast<int>(deposits->size());

  std::vector<EnergyDeposit> moliereDeposits;
  moliereDeposits.reserve(deposits->size());
  long double selectedEnergy = 0.;
  long double calorimeterEnergyForCentroid = 0.;
  long double weightedCentroidX = 0.;
  long double weightedCentroidY = 0.;
  long double weightedCentroidZ = 0.;

  for (const auto& hit : *deposits) {
    const double energy = hit.getEnergy();
    if (!std::isfinite(energy) || energy < 0.) {
      error() << "Encountered invalid ECAL deposit energy " << energy << " GeV in event " << m_event << endmsg;
      return StatusCode::FAILURE;
    }

    const dd4hep::DDSegmentation::CellID cellID = hit.getCellID();
    const int cryoId = decoder->get(cellID, m_cryoFieldName);
    const int typeId = decoder->get(cellID, m_typeFieldName);

    bool isCalorimeterDeposit = false;
    bool isFrontMaterialDeposit = false;

    if (cryoId == 0) {
      const int layerId = decoder->get(cellID, m_layerFieldName);
      isCalorimeterDeposit = addCalorimeterEnergy(typeId, layerId, energy);
      if (!isCalorimeterDeposit) {
        addUnclassifiedEnergy(energy);
      }
    } else if (cryoId == 1) {
      const bool isClassifiedCryoServices = addCryoServicesEnergy(typeId, energy);
      if (!isClassifiedCryoServices) {
        addUnclassifiedEnergy(energy);
      } else {
        isFrontMaterialDeposit = typeId == 1 || typeId == 4;
      }
    } else {
      addUnclassifiedEnergy(energy);
    }

    const bool isMoliereDeposit = isCalorimeterDeposit || (m_moliereUsesFrontMaterial && isFrontMaterialDeposit);
    if (isMoliereDeposit) {
      ++m_nMoliereHits;
    }

    if (!isMoliereDeposit) {
      continue;
    }

    const auto position = hit.getPosition();
    const Vector3 depositPosition{static_cast<double>(position.x), static_cast<double>(position.y),
                                  static_cast<double>(position.z)};
    if (!isFinite(depositPosition)) {
      error() << "Encountered non-finite selected ECAL deposit position in event " << m_event << endmsg;
      return StatusCode::FAILURE;
    }

    if (energy <= 0.) {
      continue;
    }

    if (isCalorimeterDeposit) {
      calorimeterEnergyForCentroid += energy;
      weightedCentroidX += static_cast<long double>(energy) * depositPosition.x;
      weightedCentroidY += static_cast<long double>(energy) * depositPosition.y;
      weightedCentroidZ += static_cast<long double>(energy) * depositPosition.z;
    }
    if (isMoliereDeposit) {
      selectedEnergy += energy;
      moliereDeposits.push_back({depositPosition, energy});
    }
  }

  m_calorimeterTotalEnergy = m_activeLArEnergy + m_passiveAbsorberEnergy + m_readoutPcbEnergy;
  m_cryoServicesTotalEnergy =
      m_frontCryoEnergy + m_backCryoEnergy + m_sideCryoEnergy + m_frontServicesEnergy + m_backServicesEnergy;
  m_ecalTotalEnergy = m_calorimeterTotalEnergy + m_cryoServicesTotalEnergy + m_unclassifiedEnergy;
  if (std::isfinite(m_particleEnergy)) {
    m_leakageEnergy = m_particleEnergy - m_ecalTotalEnergy;
  }

  m_moliereNormalizationEnergy = m_calorimeterTotalEnergy;
  if (m_moliereUsesFrontMaterial) {
    m_moliereNormalizationEnergy += m_frontCryoEnergy + m_frontServicesEnergy;
  }
  if (!energiesConsistent(static_cast<double>(selectedEnergy), m_moliereNormalizationEnergy)) {
    error() << "Moliere selected-hit energy " << static_cast<double>(selectedEnergy)
            << " GeV is inconsistent with normalization energy " << m_moliereNormalizationEnergy << " GeV in event "
            << m_event << endmsg;
    return StatusCode::FAILURE;
  }

  std::optional<Axis> centroidAxis;
  if (calorimeterEnergyForCentroid > 0.) {
    m_calorimeterCentroidX = static_cast<double>(weightedCentroidX / calorimeterEnergyForCentroid);
    m_calorimeterCentroidY = static_cast<double>(weightedCentroidY / calorimeterEnergyForCentroid);
    m_calorimeterCentroidZ = static_cast<double>(weightedCentroidZ / calorimeterEnergyForCentroid);
    const Vector3 centroid{m_calorimeterCentroidX, m_calorimeterCentroidY, m_calorimeterCentroidZ};
    const Vector3 centroidOrigin{m_centroidAxisOriginX, m_centroidAxisOriginY, m_centroidAxisOriginZ};
    centroidAxis = makeAxis(
        centroidOrigin, {centroid.x - centroidOrigin.x, centroid.y - centroidOrigin.y, centroid.z - centroidOrigin.z});
  }

  if (truthAxis.has_value()) {
    const auto radii = containmentRadii(moliereDeposits, *truthAxis, m_moliereNormalizationEnergy);
    if (radii.has_value()) {
      m_truthRadii = *radii;
      m_truthAxisValid = true;
      for (std::size_t index = 0; index < kSummaryRadiusIndices.size(); ++index) {
        m_truthSamples.radii[index].push_back(m_truthRadii[kSummaryRadiusIndices[index]]);
      }
      m_truthSamples.energyWeights.push_back(m_moliereNormalizationEnergy);
    }
  }

  if (centroidAxis.has_value()) {
    const auto radii = containmentRadii(moliereDeposits, *centroidAxis, m_moliereNormalizationEnergy);
    if (radii.has_value()) {
      m_centroidRadii = *radii;
      m_centroidAxisValid = true;
      for (std::size_t index = 0; index < kSummaryRadiusIndices.size(); ++index) {
        m_centroidSamples.radii[index].push_back(m_centroidRadii[kSummaryRadiusIndices[index]]);
      }
      m_centroidSamples.energyWeights.push_back(m_moliereNormalizationEnergy);
    }
  }

  verbose() << "ECAL component energies: active LAr = " << m_activeLArEnergy
            << " GeV, passive absorber = " << m_passiveAbsorberEnergy << " GeV, readout PCB = " << m_readoutPcbEnergy
            << " GeV, cryostat/services = " << m_cryoServicesTotalEnergy
            << " GeV, unclassified = " << m_unclassifiedEnergy << " GeV, leakage = " << m_leakageEnergy << " GeV"
            << endmsg;

  if (m_tree->Fill() < 0) {
    error() << "Failed to fill the event-level ECAL component energy TTree." << endmsg;
    return StatusCode::FAILURE;
  }
  ++m_summaryNEventsTotal;
  return StatusCode::SUCCESS;
}

StatusCode EnergyInECalComponents::finalize() {
  const std::lock_guard<std::mutex> stateLock(m_stateMutex);
  m_summaryTruthNValidEvents = static_cast<long long>(m_truthSamples.energyWeights.size());
  m_summaryCentroidNValidEvents = static_cast<long long>(m_centroidSamples.energyWeights.size());

  for (std::size_t index = 0; index < kSummaryRadiusIndices.size(); ++index) {
    const auto truthStatistics = summarizeSamples(m_truthSamples.radii[index], m_truthSamples.energyWeights);
    m_summaryTruthEventMean[index] = truthStatistics.mean;
    m_summaryTruthEventMedian[index] = truthStatistics.median;
    m_summaryTruthEventEnergyWeightedMean[index] = truthStatistics.energyWeightedMean;
    if (index == 0) {
      m_summaryTruthTotalWeightEnergy = truthStatistics.totalWeight;
    }

    const auto centroidStatistics = summarizeSamples(m_centroidSamples.radii[index], m_centroidSamples.energyWeights);
    m_summaryCentroidEventMean[index] = centroidStatistics.mean;
    m_summaryCentroidEventMedian[index] = centroidStatistics.median;
    m_summaryCentroidEventEnergyWeightedMean[index] = centroidStatistics.energyWeightedMean;
    if (index == 0) {
      m_summaryCentroidTotalWeightEnergy = centroidStatistics.totalWeight;
    }
  }

  if (m_moliereSummaryTree->Fill() < 0) {
    error() << "Failed to fill the Moliere summary TTree." << endmsg;
    return StatusCode::FAILURE;
  }

  info() << "Moliere summary: events = " << m_summaryNEventsTotal
         << ", truth-axis valid = " << m_summaryTruthNValidEvents
         << ", centroid-axis valid = " << m_summaryCentroidNValidEvents << endmsg;
  return Gaudi::Algorithm::finalize();
}

void EnergyInECalComponents::resetEvent() const {
  const double invalid = invalidValue();
  m_nHits = 0;
  m_nUnclassifiedHits = 0;

  m_incidentParticleIndex = -1;
  m_incidentParticlePdg = 0;
  m_nIncidentParticleCandidates = 0;

  m_particleMass = invalid;
  m_particlePx = invalid;
  m_particlePy = invalid;
  m_particlePz = invalid;
  m_particleEnergy = invalid;
  m_particleTheta = invalid;
  m_particleVertexX = invalid;
  m_particleVertexY = invalid;
  m_particleVertexZ = invalid;

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
  m_leakageEnergy = invalid;

  m_moliereUsesFrontMaterial = m_moliereIncludeFrontMaterial.value();
  m_nMoliereHits = 0;
  m_moliereNormalizationEnergy = 0.;
  m_calorimeterCentroidX = invalid;
  m_calorimeterCentroidY = invalid;
  m_calorimeterCentroidZ = invalid;
  m_truthAxisValid = false;
  m_centroidAxisValid = false;
  m_truthRadii.fill(invalid);
  m_centroidRadii.fill(invalid);

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
