#include "SamplingFractionInLayers.h"

// FCCSW
#include "k4Interface/IGeoSvc.h"

// datamodel
#include "edm4hep/SimCalorimeterHitCollection.h"

#include "GaudiKernel/ITHistSvc.h"
#include "TH1F.h"
#include "TVector2.h"

#include <algorithm>
#include <cmath>

// DD4hep
#include "DD4hep/Detector.h"
#include "DD4hep/Readout.h"

DECLARE_COMPONENT(SamplingFractionInLayers)

SamplingFractionInLayers::SamplingFractionInLayers(const std::string& aName, ISvcLocator* aSvcLoc)
    : Gaudi::Algorithm(aName, aSvcLoc), m_histSvc("THistSvc", "SamplingFractionInLayers"),
      m_geoSvc("GeoSvc", "SamplingFractionInLayers"), m_totalEnergy(nullptr), m_totalActiveEnergy(nullptr),
      m_sf(nullptr) {
  declareProperty("deposits", m_deposits, "Energy deposits in sampling calorimeter (input)");
}
SamplingFractionInLayers::~SamplingFractionInLayers() {}

StatusCode SamplingFractionInLayers::initialize() {
  if (Gaudi::Algorithm::initialize().isFailure()) {
    return StatusCode::FAILURE;
  }
  // check if readouts exist
  if (m_geoSvc->getDetector()->readouts().find(m_readoutName) == m_geoSvc->getDetector()->readouts().end()) {
    error() << "Readout <<" << m_readoutName << ">> does not exist." << endmsg;
    return StatusCode::FAILURE;
  }
  info() << "Histogram energy ranges are fully data-driven with margin factor " << m_rangeMargin
         << ". Property energyAxis is ignored." << endmsg;
  if (m_energyBinWidth <= 0.) {
    error() << "Property energyBinWidth must be > 0. Current value = " << m_energyBinWidth << endmsg;
    return StatusCode::FAILURE;
  }

  // Register placeholder histograms in initialize. In finalize they are re-binned and filled from cache.
  const double initMax = static_cast<double>(m_energyBinWidth);
  for (uint i = 0; i < m_numLayers; i++) {
    m_totalEnLayers.push_back(new TH1F(("ecal_totalEnergy_layer" + std::to_string(i)).c_str(),
                                       ("Total deposited energy in layer " + std::to_string(i)).c_str(), 1, 0,
                                       initMax));
    if (m_histSvc->regHist("/rec/ecal_total_layer" + std::to_string(i), m_totalEnLayers.back()).isFailure()) {
      error() << "Couldn't register histogram" << endmsg;
      return StatusCode::FAILURE;
    }

    m_activeEnLayers.push_back(new TH1F(("ecal_activeEnergy_layer" + std::to_string(i)).c_str(),
                                        ("Deposited energy in active material, in layer " + std::to_string(i)).c_str(),
                                        1, 0, initMax));
    if (m_histSvc->regHist("/rec/ecal_active_layer" + std::to_string(i), m_activeEnLayers.back()).isFailure()) {
      error() << "Couldn't register histogram" << endmsg;
      return StatusCode::FAILURE;
    }

    m_sfLayers.push_back(new TH1F(("ecal_sf_layer" + std::to_string(i)).c_str(),
                                  ("SF for layer " + std::to_string(i)).c_str(), m_samplingFractionBins, 0, 1));
    if (m_histSvc->regHist("/rec/ecal_sf_layer" + std::to_string(i), m_sfLayers.back()).isFailure()) {
      error() << "Couldn't register histogram" << endmsg;
      return StatusCode::FAILURE;
    }
  }

  m_totalEnergy = new TH1F("ecal_totalEnergy", "Total deposited energy", 1, 0, initMax);
  if (m_histSvc->regHist("/rec/ecal_total", m_totalEnergy).isFailure()) {
    error() << "Couldn't register histogram" << endmsg;
    return StatusCode::FAILURE;
  }

  m_totalActiveEnergy = new TH1F("ecal_active", "Deposited energy in active material", 1, 0, initMax);
  if (m_histSvc->regHist("/rec/ecal_active", m_totalActiveEnergy).isFailure()) {
    error() << "Couldn't register histogram" << endmsg;
    return StatusCode::FAILURE;
  }

  m_sf = new TH1F("ecal_sf", "Sampling fraction", m_samplingFractionBins, 0, 1);
  if (m_histSvc->regHist("/rec/ecal_sf", m_sf).isFailure()) {
    error() << "Couldn't register histogram" << endmsg;
    return StatusCode::FAILURE;
  }

  m_eventTotalEnergy.clear();
  m_eventTotalActiveEnergy.clear();
  m_eventLayerTotalEnergy.clear();
  m_eventLayerActiveEnergy.clear();
  m_maxTotalEnergy = 0.;
  m_maxTotalActiveEnergy = 0.;
  m_maxLayerTotalEnergy.assign(m_numLayers, 0.);
  m_maxLayerActiveEnergy.assign(m_numLayers, 0.);

  return StatusCode::SUCCESS;
}

StatusCode SamplingFractionInLayers::execute(const EventContext&) const {
  auto decoder = m_geoSvc->getDetector()->readout(m_readoutName).idSpec().decoder();
  double sumE = 0.;
  std::vector<double> sumElayers;
  double sumEactive = 0.;
  std::vector<double> sumEactiveLayers;
  sumElayers.assign(m_numLayers, 0);
  sumEactiveLayers.assign(m_numLayers, 0);

  const auto deposits = m_deposits.get();
  for (const auto& hit : *deposits) {
    dd4hep::DDSegmentation::CellID cID = hit.getCellID();
    auto id = decoder->get(cID, m_layerFieldName);
    if (id < 0 || static_cast<uint>(id) >= m_numLayers) {
      warning() << "Skipping hit with out-of-range layer id " << id << " (configured numLayers=" << m_numLayers
                << ")" << endmsg;
      continue;
    }
    sumElayers[id] += hit.getEnergy();
    // check if energy was deposited in the calorimeter (active/passive material)
    if (id >= m_firstLayerId) {
      sumE += hit.getEnergy();
      // active material of calorimeter
      auto activeField = decoder->get(cID, m_activeFieldName);
      if (activeField == m_activeFieldValue) {
        sumEactive += hit.getEnergy();
        sumEactiveLayers[id] += hit.getEnergy();
      }
    }
  }

  // Cache per-event quantities. Histograms are filled in finalize.
  m_eventTotalEnergy.push_back(sumE);
  m_eventTotalActiveEnergy.push_back(sumEactive);
  m_eventLayerTotalEnergy.push_back(sumElayers);
  m_eventLayerActiveEnergy.push_back(sumEactiveLayers);

  m_maxTotalEnergy = std::max(m_maxTotalEnergy, sumE);
  m_maxTotalActiveEnergy = std::max(m_maxTotalActiveEnergy, sumEactive);

  for (uint i = 0; i < m_numLayers; i++) {
    m_maxLayerTotalEnergy[i] = std::max(m_maxLayerTotalEnergy[i], sumElayers[i]);
    m_maxLayerActiveEnergy[i] = std::max(m_maxLayerActiveEnergy[i], sumEactiveLayers[i]);
    if (i < m_firstLayerId) {
      debug() << "total energy deposited outside the calorimeter detector = " << sumElayers[i] << endmsg;
    } else {
      debug() << "total energy in layer " << i << " = " << sumElayers[i] << " active = " << sumEactiveLayers[i]
              << endmsg;
    }
  }
  return StatusCode::SUCCESS;
}

StatusCode SamplingFractionInLayers::finalize() {
  if (m_energyBinWidth <= 0.) {
    error() << "Property energyBinWidth must be > 0. Current value = " << m_energyBinWidth << endmsg;
    return StatusCode::FAILURE;
  }
  if (m_rangeMargin <= 0.) {
    error() << "Property rangeMargin must be > 0. Current value = " << m_rangeMargin << endmsg;
    return StatusCode::FAILURE;
  }

  const auto computeBinsAndMax = [&](double observedMax) {
    const double safeMax = observedMax * static_cast<double>(m_rangeMargin);
    const double axisMax = std::max(safeMax, static_cast<double>(m_energyBinWidth));
    const int nBins = std::max(1, static_cast<int>(std::ceil(axisMax / static_cast<double>(m_energyBinWidth))));
    const double xMax = nBins * static_cast<double>(m_energyBinWidth);
    return std::make_pair(nBins, xMax);
  };

  // Re-bin existing registered histograms with dynamic ranges and fixed bin width.
  for (uint i = 0; i < m_numLayers; i++) {
    const auto layerTotalCfg = computeBinsAndMax(m_maxLayerTotalEnergy[i]);
    m_totalEnLayers[i]->SetBins(layerTotalCfg.first, 0., layerTotalCfg.second);

    const auto layerActiveCfg = computeBinsAndMax(m_maxLayerActiveEnergy[i]);
    m_activeEnLayers[i]->SetBins(layerActiveCfg.first, 0., layerActiveCfg.second);
  }

  const auto totalCfg = computeBinsAndMax(m_maxTotalEnergy);
  m_totalEnergy->SetBins(totalCfg.first, 0., totalCfg.second);

  const auto totalActiveCfg = computeBinsAndMax(m_maxTotalActiveEnergy);
  m_totalActiveEnergy->SetBins(totalActiveCfg.first, 0., totalActiveCfg.second);

  // Fill histograms from cached event-level quantities.
  for (size_t evt = 0; evt < m_eventTotalEnergy.size(); ++evt) {
    m_totalEnergy->Fill(m_eventTotalEnergy[evt]);
    m_totalActiveEnergy->Fill(m_eventTotalActiveEnergy[evt]);
    if (m_eventTotalEnergy[evt] > 0.) {
      m_sf->Fill(m_eventTotalActiveEnergy[evt] / m_eventTotalEnergy[evt]);
    }

    for (uint i = 0; i < m_numLayers; i++) {
      m_totalEnLayers[i]->Fill(m_eventLayerTotalEnergy[evt][i]);
      m_activeEnLayers[i]->Fill(m_eventLayerActiveEnergy[evt][i]);
      if (m_eventLayerTotalEnergy[evt][i] > 0.) {
        m_sfLayers[i]->Fill(m_eventLayerActiveEnergy[evt][i] / m_eventLayerTotalEnergy[evt][i]);
      }
    }
  }

  return Gaudi::Algorithm::finalize();
}
