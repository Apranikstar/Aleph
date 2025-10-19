#ifndef SELECTIONUTILS_H
#define SELECTIONUTILS_H

/*
  Selection utilities for filtering particles and events.
  Fully RDataFrame compatible (using ROOT::VecOps::RVec).

  Includes:
    - sel_charged: selects reconstructed particles by absolute charge.
    - sel_class_filter: filters events based on their class bit.
    - sel_runs_filter: filters events by allowed run numbers.

  Example usage in RDataFrame:

    df = df.Define("charged_particles",
                   "FCCAnalyses::AlephSelection::sel_charged(1)(ReconstructedParticles)")
           .Filter("FCCAnalyses::AlephSelection::sel_class_filter(16)(ClassBitset)")
           .Filter("FCCAnalyses::AlephSelection::sel_runs_filter(allowedRuns)(EventHeader)")
*/

#include "edm4hep/ReconstructedParticleCollection.h"
#include "edm4hep/EventHeaderCollection.h"
#include <set>
#include <bitset>
#include <cmath>
#include <ROOT/RVec.hxx>

namespace FCCAnalyses { namespace AlephSelection {

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

}} // namespace FCCAnalyses::AlephSelection

#endif

