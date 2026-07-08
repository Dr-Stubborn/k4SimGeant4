#ifndef DETSTUDIES_ENERGYINECALCOMPONENTS_H
#define DETSTUDIES_ENERGYINECALCOMPONENTS_H

// GAUDI
#include "Gaudi/Algorithm.h"
#include "GaudiKernel/ServiceHandle.h"

// Key4HEP
#include "k4FWCore/DataHandle.h"

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
 *  material category. The category is decoded from the ECAL readout fields:
 *  cryo, type, and layer.
 */
class EnergyInECalComponents : public Gaudi::Algorithm {
public:
  explicit EnergyInECalComponents(const std::string&, ISvcLocator*);
  virtual ~EnergyInECalComponents();

  virtual StatusCode initialize() final;
  virtual StatusCode execute(const EventContext&) const final;
  virtual StatusCode finalize() final;

private:
  void resetEvent() const;
  bool addCalorimeterEnergy(int typeId, int layerId, double energy) const;
  bool addCryoServicesEnergy(int typeId, double energy) const;
  void addUnclassifiedEnergy(double energy) const;

  ServiceHandle<ITHistSvc> m_histSvc;
  ServiceHandle<IGeoSvc> m_geoSvc;

  mutable k4FWCore::DataHandle<edm4hep::SimCalorimeterHitCollection> m_deposits{
      "ECalBarrelHits", Gaudi::DataHandle::Reader, this};
  mutable k4FWCore::DataHandle<edm4hep::MCParticleCollection> m_particle{"GenParticles",
                                                                         Gaudi::DataHandle::Reader, this};

  Gaudi::Property<std::string> m_readoutName{this, "readoutName", "", "Name of the detector readout"};
  Gaudi::Property<unsigned int> m_numLayers{this, "numLayers", 11, "Number of ECAL layers"};
  Gaudi::Property<std::string> m_cryoFieldName{this, "cryoFieldName", "cryo", "Readout field for cryostat flag"};
  Gaudi::Property<std::string> m_typeFieldName{this, "typeFieldName", "type", "Readout field for material type"};
  Gaudi::Property<std::string> m_layerFieldName{this, "layerFieldName", "layer", "Readout field for layer"};
  Gaudi::Property<std::string> m_treePath{this, "treePath", "/rec/ecal_component_energy",
                                          "THistSvc path for the event-level TTree"};

  TTree* m_tree;

  mutable int m_event{0};
  mutable int m_nHits{0};
  mutable int m_nUnclassifiedHits{0};

  mutable double m_particleMass{0.};
  mutable double m_particlePx{0.};
  mutable double m_particlePy{0.};
  mutable double m_particlePz{0.};
  mutable double m_particleEnergy{0.};
  mutable double m_particleTheta{0.};

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

  mutable std::vector<double> m_activeLArEnergyByLayer;
  mutable std::vector<double> m_passiveAbsorberEnergyByLayer;
  mutable std::vector<double> m_readoutPcbEnergyByLayer;
  mutable std::vector<double> m_calorimeterEnergyByLayer;
};

#endif /* DETSTUDIES_ENERGYINECALCOMPONENTS_H */
