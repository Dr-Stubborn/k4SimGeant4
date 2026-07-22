#ifndef DETSTUDIES_ENERGYINECALCOMPONENTS_H
#define DETSTUDIES_ENERGYINECALCOMPONENTS_H

// GAUDI
#include "Gaudi/Algorithm.h"
#include "GaudiKernel/ServiceHandle.h"

// Key4HEP
#include "k4FWCore/DataHandle.h"

#include <array>
#include <mutex>
#include <string>
#include <vector>

class IGeoSvc;
class ITHistSvc;
class TTree;

namespace edm4hep {
class MCParticleCollection;
class SimCalorimeterHitCollection;
} // namespace edm4hep

/** @class EnergyInECalComponents EnergyInECalComponents.h
 *
 *  Writes one event-level TTree entry with raw deposited energy split by ECAL
 *  material category and exact energy-containment radii around truth and
 *  calorimeter-centroid shower axes. Run-level radius statistics are written
 *  to a separate summary TTree.
 */
class EnergyInECalComponents : public Gaudi::Algorithm {
public:
  explicit EnergyInECalComponents(const std::string&, ISvcLocator*);
  virtual ~EnergyInECalComponents();

  virtual StatusCode initialize() final;
  virtual StatusCode execute(const EventContext&) const final;
  virtual StatusCode finalize() final;

private:
  struct MoliereSamples {
    std::array<std::vector<double>, 3> radii;
    std::vector<double> energyWeights;
  };

  void resetEvent() const;
  bool addCalorimeterEnergy(int typeId, int layerId, double energy) const;
  bool addCryoServicesEnergy(int typeId, double energy) const;
  void addUnclassifiedEnergy(double energy) const;

  ServiceHandle<ITHistSvc> m_histSvc;
  ServiceHandle<IGeoSvc> m_geoSvc;

  mutable k4FWCore::DataHandle<edm4hep::SimCalorimeterHitCollection> m_deposits{"ECalBarrelHits",
                                                                                Gaudi::DataHandle::Reader, this};
  mutable k4FWCore::DataHandle<edm4hep::MCParticleCollection> m_particle{"GenParticles", Gaudi::DataHandle::Reader,
                                                                         this};

  Gaudi::Property<std::string> m_readoutName{this, "readoutName", "", "Name of the detector readout"};
  Gaudi::Property<unsigned int> m_numLayers{this, "numLayers", 11, "Number of ECAL layers"};
  Gaudi::Property<std::string> m_cryoFieldName{this, "cryoFieldName", "cryo", "Readout field for cryostat flag"};
  Gaudi::Property<std::string> m_typeFieldName{this, "typeFieldName", "type", "Readout field for material type"};
  Gaudi::Property<std::string> m_layerFieldName{this, "layerFieldName", "layer", "Readout field for layer"};
  Gaudi::Property<std::string> m_treePath{this, "treePath", "/rec/ecal_component_energy",
                                          "THistSvc path for the event-level TTree"};
  Gaudi::Property<std::string> m_moliereSummaryTreePath{this, "moliereSummaryTreePath", "/rec/ecal_moliere_summary",
                                                        "THistSvc path for the Moliere summary TTree"};
  Gaudi::Property<bool> m_moliereIncludeFrontMaterial{
      this, "moliereIncludeFrontMaterial", false,
      "Include front cryostat and front services deposits in Moliere containment"};
  Gaudi::Property<std::vector<double>> m_centroidAxisOrigin{
      this, "centroidAxisOrigin", {0., 0., 0.}, "Origin (x, y, z) in mm of the calorimeter-centroid shower axis"};

  TTree* m_tree;
  TTree* m_moliereSummaryTree;
  mutable std::mutex m_stateMutex;

  mutable int m_event{0};
  mutable int m_nHits{0};
  mutable int m_nUnclassifiedHits{0};

  mutable int m_incidentParticleIndex{-1};
  mutable int m_incidentParticlePdg{0};
  mutable int m_nIncidentParticleCandidates{0};

  mutable double m_particleMass{0.};
  mutable double m_particlePx{0.};
  mutable double m_particlePy{0.};
  mutable double m_particlePz{0.};
  mutable double m_particleEnergy{0.};
  mutable double m_particleTheta{0.};
  mutable double m_particleVertexX{0.};
  mutable double m_particleVertexY{0.};
  mutable double m_particleVertexZ{0.};

  mutable double m_activeLArEnergy{0.};
  mutable double m_passiveAbsorberEnergy{0.};
  mutable double m_readoutPcbEnergy{0.};
  mutable double m_frontCryoEnergy{0.};
  mutable double m_backCryoEnergy{0.};
  mutable double m_sideCryoEnergy{0.};
  mutable double m_frontServicesEnergy{0.};
  mutable double m_backServicesEnergy{0.};

  mutable double m_calorimeterTotalEnergy{0.};
  mutable double m_cryoServicesTotalEnergy{0.};
  mutable double m_ecalTotalEnergy{0.};
  mutable double m_unclassifiedEnergy{0.};
  mutable double m_leakageEnergy{0.};

  mutable bool m_moliereUsesFrontMaterial{false};
  mutable int m_nMoliereHits{0};
  mutable double m_moliereNormalizationEnergy{0.};
  mutable double m_calorimeterCentroidX{0.};
  mutable double m_calorimeterCentroidY{0.};
  mutable double m_calorimeterCentroidZ{0.};
  mutable bool m_truthAxisValid{false};
  mutable bool m_centroidAxisValid{false};
  mutable std::array<double, 5> m_truthRadii{};
  mutable std::array<double, 5> m_centroidRadii{};

  mutable std::vector<double> m_activeLArEnergyByLayer;
  mutable std::vector<double> m_passiveAbsorberEnergyByLayer;
  mutable std::vector<double> m_readoutPcbEnergyByLayer;
  mutable std::vector<double> m_calorimeterEnergyByLayer;

  double m_centroidAxisOriginX{0.};
  double m_centroidAxisOriginY{0.};
  double m_centroidAxisOriginZ{0.};
  bool m_summaryUsesFrontMaterial{false};

  mutable long long m_summaryNEventsTotal{0};
  mutable long long m_summaryTruthNValidEvents{0};
  mutable long long m_summaryCentroidNValidEvents{0};
  mutable double m_summaryTruthTotalWeightEnergy{0.};
  mutable double m_summaryCentroidTotalWeightEnergy{0.};
  mutable std::array<double, 3> m_summaryTruthEventMean{};
  mutable std::array<double, 3> m_summaryTruthEventMedian{};
  mutable std::array<double, 3> m_summaryTruthEventEnergyWeightedMean{};
  mutable std::array<double, 3> m_summaryCentroidEventMean{};
  mutable std::array<double, 3> m_summaryCentroidEventMedian{};
  mutable std::array<double, 3> m_summaryCentroidEventEnergyWeightedMean{};
  mutable MoliereSamples m_truthSamples;
  mutable MoliereSamples m_centroidSamples;
};

#endif /* DETSTUDIES_ENERGYINECALCOMPONENTS_H */
