#ifndef SELECTIONUTILS_H
#define SELECTIONUTILS_H

/*
  Selection utilities for filtering particles and events.
  Fully RDataFrame compatible (using ROOT::VecOps::RVec).

  Includes:
    - sel_charged: selects reconstructed particles by absolute charge.
    - sel_class_filter: filters events based on their class bit.
    - sel_runs_filter: filters events by allowed run numbers.
    - get_isEl / get_isMu / get_isChargedHad / get_isNeutralHad / get_isGamma:
      classify jet constituents by particle type.

  Example usage in RDataFrame:

    df = df.Define("charged_particles",
                   "FCCAnalyses::AlephSelection::sel_charged(1)(ReconstructedParticles)")
           .Filter("FCCAnalyses::AlephSelection::sel_class_filter(16)(ClassBitset)")
           .Filter("FCCAnalyses::AlephSelection::sel_runs_filter(allowedRuns)(EventHeader)");

    df = df.Define("isMu", "FCCAnalyses::AlephSelection::get_isMu(JetConstituents)")
           .Define("n_muons_per_jet", "Sum(isMu)");
*/

#include "edm4hep/ReconstructedParticleCollection.h"
#include "edm4hep/EventHeaderCollection.h"
#include <set>
#include <bitset>
#include <cmath>
#include <ROOT/RVec.hxx>
#include "FCCAnalyses/JetConstituentsUtils.h"
#include "FCCAnalyses/ReconstructedParticle.h"
#include "FCCAnalyses/ReconstructedParticle2Track.h"
#include "FCCAnalyses/ReconstructedParticle2MC.h"
#include "edm4hep/MCParticleData.h"
#include "edm4hep/Track.h"
#include "edm4hep/TrackData.h"
#include "edm4hep/Cluster.h"
#include "edm4hep/ClusterData.h"
#include "edm4hep/CalorimeterHitData.h"
#include "edm4hep/ReconstructedParticleData.h"
#include "edm4hep/EDM4hepVersion.h"
#include "FCCAnalyses/JetClusteringUtils.h"
// #include "FCCAnalyses/ExternalRecombiner.h"
#include "fastjet/JetDefinition.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/Selector.hh"

namespace FCCAnalyses { namespace AlephSelection {

namespace rv = ROOT::VecOps;

// -----------------------------------
// Type aliases for jet constituent data
// -----------------------------------

using FCCAnalysesJetConstituents = rv::RVec<edm4hep::ReconstructedParticleData>;
using FCCAnalysesJetConstituentsData = rv::RVec<float>;

// ----------------------------
// Event and particle selectors
// ----------------------------

/// Selects charged particles based on their absolute charge.
struct sel_charged {
  const int m_charge;
  sel_charged(int arg_charge) : m_charge(arg_charge) {};

  edm4hep::ReconstructedParticleCollection
  operator()(const edm4hep::ReconstructedParticleCollection& in_coll) const {
    edm4hep::ReconstructedParticleCollection result;
    result.setSubsetCollection();

    for (const auto& i : in_coll) {
      if (std::abs(i.getCharge()) == m_charge) {
        result.push_back(i);
      }
    }
    return result;
  }
};

/// Filters events based on their class bit (RVec-compatible)
struct sel_class_filter {
  const int m_class;
  sel_class_filter(int arg_class) : m_class(arg_class) {};

  bool operator()(const ROOT::VecOps::RVec<uint32_t>& bitset_coll) const {
    if (bitset_coll.empty()) return false;
    std::bitset<32> bits(bitset_coll[0]);
    return bits[m_class - 1];
  }
};

/// Filters events by run number (RVec-compatible)
struct sel_runs_filter {
  const std::set<int>& m_runs_set;
  sel_runs_filter(const std::set<int>& arg_runs_set) : m_runs_set(arg_runs_set) {};

  bool operator()(const ROOT::VecOps::RVec<edm4hep::EventHeader>& event_header) const {
    if (event_header.empty()) return false;
    return m_runs_set.count(event_header[0].getRunNumber()) > 0;
  }
};

// --------------------------------------
// Jet constituent particle identification
// --------------------------------------

rv::RVec<FCCAnalysesJetConstituentsData>
get_isEl(const rv::RVec<FCCAnalysesJetConstituents>& jcs) {
  rv::RVec<FCCAnalysesJetConstituentsData> out;
  out.reserve(jcs.size());
  for (const auto& jet : jcs) {
    FCCAnalysesJetConstituentsData mask;
    mask.reserve(jet.size());
    for (const auto& c : jet)
      mask.push_back((std::abs(c.charge) > 0 && std::abs(c.mass - 0.000511) < 1e-5) ? 1.f : 0.f);
    out.push_back(std::move(mask));
  }
  return out;
}

rv::RVec<FCCAnalysesJetConstituentsData>
get_isMu(const rv::RVec<FCCAnalysesJetConstituents>& jcs) {
  rv::RVec<FCCAnalysesJetConstituentsData> out;
  out.reserve(jcs.size());
  for (const auto& jet : jcs) {
    FCCAnalysesJetConstituentsData mask;
    mask.reserve(jet.size());
    for (const auto& c : jet)
      mask.push_back((std::abs(c.charge) > 0 && std::abs(c.mass - 0.105658) < 1e-3) ? 1.f : 0.f);
    out.push_back(std::move(mask));
  }
  return out;
}

rv::RVec<FCCAnalysesJetConstituentsData>
get_isChargedHad(const rv::RVec<FCCAnalysesJetConstituents>& jcs) {
  rv::RVec<FCCAnalysesJetConstituentsData> out;
  out.reserve(jcs.size());
  for (const auto& jet : jcs) {
    FCCAnalysesJetConstituentsData mask;
    mask.reserve(jet.size());
    for (const auto& c : jet)
      mask.push_back((std::abs(c.charge) > 0 && std::abs(c.mass - 0.13957) < 1e-3) ? 1.f : 0.f);
    out.push_back(std::move(mask));
  }
  return out;
}

rv::RVec<FCCAnalysesJetConstituentsData>
get_isNeutralHad(const rv::RVec<FCCAnalysesJetConstituents>& jcs) {
  rv::RVec<FCCAnalysesJetConstituentsData> out;
  out.reserve(jcs.size());
  for (const auto& jet : jcs) {
    FCCAnalysesJetConstituentsData mask;
    mask.reserve(jet.size());
    for (const auto& c : jet)
#if edm4hep_VERSION > EDM4HEP_VERSION(0, 10, 5)
      mask.push_back((c.PDG == 130) ? 1.f : 0.f);
#else
      mask.push_back((c.type == 130) ? 1.f : 0.f);
#endif
    out.push_back(std::move(mask));
  }
  return out;
}

rv::RVec<FCCAnalysesJetConstituentsData>
get_isGamma(const rv::RVec<FCCAnalysesJetConstituents>& jcs) {
  rv::RVec<FCCAnalysesJetConstituentsData> out;
  out.reserve(jcs.size());
  for (const auto& jet : jcs) {
    FCCAnalysesJetConstituentsData mask;
    mask.reserve(jet.size());
    for (const auto& c : jet)
#if edm4hep_VERSION > EDM4HEP_VERSION(0, 10, 5)
      mask.push_back((c.PDG == 22) ? 1.f : 0.f);
#else
      mask.push_back((c.type == 22) ? 1.f : 0.f);
#endif
    out.push_back(std::move(mask));
  }
  return out;
}

}} // namespace FCCAnalyses::AlephSelection

#endif

