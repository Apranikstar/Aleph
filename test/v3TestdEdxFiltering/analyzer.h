#ifndef SELECTIONUTILS_H
#define SELECTIONUTILS_H


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
#include "FCCAnalyses/JetClusteringUtils.h"
// #include "FCCAnalyses/ExternalRecombiner.h"
#include "FCCAnalyses/MCParticle.h"

#include "edm4hep/MCParticleData.h"
#include "edm4hep/Track.h"
#include "edm4hep/TrackData.h"
#include "edm4hep/Cluster.h"
#include "edm4hep/ClusterData.h"
#include "edm4hep/CalorimeterHitData.h"
#include "edm4hep/ReconstructedParticleData.h"
#include "edm4hep/EDM4hepVersion.h"
#include "edm4hep/RecDqdx.h"

#include "fastjet/JetDefinition.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/Selector.hh"

#include <iostream>
#include <algorithm>



namespace FCCAnalyses { namespace AlephSelection {

namespace rv = ROOT::VecOps;

// -----------------------------------
// Type aliases for jet constituent data
// -----------------------------------
//////////////////////////////////////////////////////////////////////////////////////////
using FCCAnalysesJetConstituents = rv::RVec<edm4hep::ReconstructedParticleData>;
using FCCAnalysesJetConstituentsData = rv::RVec<float>;
//////////////////////////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////////////////////////
// -----------------------------------
// The following will first filter the 
// event to check if it's a qq event
// then it return the PID of qq 
// -----------------------------------

float getJetPID(const ROOT::VecOps::RVec<uint32_t>& ClassBit,
                   const ROOT::VecOps::RVec<edm4hep::MCParticleData>& particles) {
    // Check if bit 15 (16-1) is set in the first element of ClassBit
    if (ClassBit.empty() || !std::bitset<32>(ClassBit[0])[15]) {
        return -1.0f; // Not selected
    }

    // Bit 15 is true — find the first quark (|PDG| in 1..5)
    float result = -1.0f;
    for (const auto& particle : particles) {
        if (std::abs(particle.PDG) > 0 && std::abs(particle.PDG) < 6) {
            result = static_cast<float>(std::abs(particle.PDG));
            break;
        }
    }

    return result;
}
//////////////////////////////////////////////////////////////////////////////////////////




//////////////////////////////////////////////////////////////////////////////////////////
// -----------------------------------
// Code to get the type of pfcands
// inside a jet
// -----------------------------------

rv::RVec<FCCAnalysesJetConstituentsData> build_constituents_Types(
    const rv::RVec<edm4hep::ParticleIDData> &rpid,
    const std::vector<std::vector<int>> &indices)
{
    rv::RVec<FCCAnalysesJetConstituentsData> jcs;

    for (const auto &jet_index : indices)
    {
        FCCAnalysesJetConstituentsData jc;

        for (const auto &const_index : jet_index)
        {
            jc.push_back(rpid.at(const_index).type);
        }

        jcs.push_back(jc);
    }

    return jcs;
}


rv::RVec<FCCAnalysesJetConstituentsData>
get_isType(const rv::RVec<FCCAnalysesJetConstituentsData>& jcs, float type) {
    rv::RVec<FCCAnalysesJetConstituentsData> out;
    out.reserve(jcs.size());

    for (const auto& jet : jcs) {
        FCCAnalysesJetConstituentsData mask;
        mask.reserve(jet.size());

        for (const auto& c : jet) {
            if (c == type)
                mask.push_back(1);
            else
                mask.push_back(0);
        }

        out.push_back(std::move(mask));
    }

    return out;
}
//////////////////////////////////////////////////////////////////////////////////////////

// -----------------------------------
// Calculate the primary vertex
// four momentum
// -----------------------------------


TLorentzVector get_EventPrimaryVertexP4(const ROOT::VecOps::RVec<edm4hep::MCParticleData>& in, int m_genstatus = 21) {
    TLorentzVector result(-1e12, -1e12, -1e12, -1e12);
    bool found_py8 = false;

    // First, look for generatorStatus == m_genstatus (e.g., 21 for Pythia8 hard process incoming)
    for (const auto& p : in) {
        // Uncomment below to use actual generator status check
        // if (p.generatorStatus == m_genstatus) {
        if (1) { // currently always true like in your struct
            TLorentzVector res(
                p.vertex.x,
                p.vertex.y,
                p.vertex.z,
                p.time * 1.0e3 * 2.99792458e+8 // convert time to mm
            );
            result = res;
            found_py8 = true;
            break;
        }
    }

    // Fallback: look for genStatus == 2 with non-zero z vertex
    if (!found_py8) {
        for (const auto& p : in) {
            if (p.generatorStatus == 2 && std::abs(p.vertex.z) > 1.e-12) {
                TLorentzVector res(
                    p.vertex.x,
                    p.vertex.y,
                    p.vertex.z,
                    p.time * 1.0e3 * 2.99792458e+8
                );
                result = res;
                break;
            }
        }
    }

    return result;
}
//////////////////////////////////////////////////////////////////////////////////////////


// -----------------------------------
// Calculate the dEdx variables for 
// pads and wires.
// -----------------------------------



struct build_constituents_dEdx{
    rv::RVec<rv::RVec<edm4hep::RecDqdxData>>
    operator()(const rv::RVec<edm4hep::ReconstructedParticleData> &recoParticles,
             const rv::RVec<int> &_recoParticlesIndices,
             const rv::RVec<edm4hep::RecDqdxData> &dEdxCollection,
             const rv::RVec<int> &_dEdxIndicesCollection, 
             const std::vector<std::vector<int>> &jet_indices) const
    { 
        rv::RVec<rv::RVec<edm4hep::RecDqdxData>> dedx_constituents;

        // The links dEdx -> Track and RecoPart -> Track are one-directional, we need a map to store
        // Track.index -> dEdx to not have to loop everytime 
        // in addition, the object itself is stored in <Collection>
        // while the relations (=indices we need for links) are on _<Collection>
        std::unordered_map<int, edm4hep::RecDqdxData> track_index_to_dEdx;
        for (size_t i = 0; i < _dEdxIndicesCollection.size(); ++i) {
          int track_index = _dEdxIndicesCollection[i];
          edm4hep::RecDqdxData dedx_value = dEdxCollection[i];
          track_index_to_dEdx[track_index] = dedx_value;
        }

        //now, for each jet loop over the indices of the jet constituents provided by teh JetClusteringUtils
        // retrieve the associated RecoParticle
        // from there get the link to the Track from the corresponding index collection
        for (const auto &jet_const_indices : jet_indices) { //loop over jets
          rv::RVec<edm4hep::RecDqdxData> jet_dEdx;

          for (int constituent_index : jet_const_indices) { // loop over jet constituents
            const auto &recoPart = recoParticles[constituent_index];

            //loop over tracks associated to the RecoPart (should always be one in Aleph data)
            for (int track = recoPart.tracks_begin; track < recoPart.tracks_end; ++track) {
                 int track_index = _recoParticlesIndices[track]; //this should be the same index used in the link from dEdx to track

                  //find the matching dEdx in the map
                  if (track_index_to_dEdx.count(track_index)) {
                    jet_dEdx.push_back(track_index_to_dEdx[track_index]);
                  }
            }
          }
          dedx_constituents.push_back(jet_dEdx); 
        }
        return dedx_constituents;
    }
};


rv::RVec<rv::RVec<float>> get_dEdx_type(const rv::RVec<rv::RVec<edm4hep::RecDqdxData>> &dedx_vec) {
  rv::RVec<rv::RVec<float>> values;
  for (const auto &inner_vec : dedx_vec) {
    rv::RVec<float> inner_values;
    for (const auto &d : inner_vec) {
      inner_values.push_back(d.dQdx.type);
    }
    values.push_back(inner_values);
  }
  return values;
}

rv::RVec<rv::RVec<float>> get_dEdx_value(const rv::RVec<rv::RVec<edm4hep::RecDqdxData>> &dedx_vec) {
  rv::RVec<rv::RVec<float>> values;
  for (const auto &inner_vec : dedx_vec) {
    rv::RVec<float> inner_values;
    for (const auto &d : inner_vec) {
      inner_values.push_back(d.dQdx.value);
    }
    values.push_back(inner_values);
  }
  return values;
}

rv::RVec<rv::RVec<float>> get_dEdx_error(const rv::RVec<rv::RVec<edm4hep::RecDqdxData>> &dedx_vec) {
  rv::RVec<rv::RVec<float>> values;
  for (const auto &inner_vec : dedx_vec) {
    rv::RVec<float> inner_values;
    for (const auto &d : inner_vec) {
      inner_values.push_back(d.dQdx.error);
    }
    values.push_back(inner_values);
  }
  return values;
}






struct build_constituents_dEdx_filtered {
    rv::RVec<rv::RVec<edm4hep::RecDqdxData>>
    operator()(const rv::RVec<FCCAnalysesJetConstituentsData> &jet_constituents_types,
               const rv::RVec<edm4hep::ReconstructedParticleData> &recoParticles,
               const rv::RVec<int> &_recoParticlesIndices,
               const rv::RVec<edm4hep::RecDqdxData> &dEdxCollection,
               const rv::RVec<int> &_dEdxIndicesCollection,
               const std::vector<std::vector<int>> &jet_indices,
               float targetType) const {

        rv::RVec<rv::RVec<edm4hep::RecDqdxData>> dedx_filtered;

        // Build map Track.index -> dEdx object
        std::unordered_map<int, edm4hep::RecDqdxData> track_index_to_dEdx;
        for (size_t i = 0; i < _dEdxIndicesCollection.size(); ++i) {
            track_index_to_dEdx[_dEdxIndicesCollection[i]] = dEdxCollection[i];
        }

        // Loop over jets
        for (size_t j = 0; j < jet_indices.size(); ++j) {
            const auto &constituents = jet_indices[j];
            const auto &types = jet_constituents_types[j];

            rv::RVec<edm4hep::RecDqdxData> jet_dEdx_masked;

            // Loop over jet constituents
            for (size_t c = 0; c < constituents.size(); ++c) {
                if (types[c] != targetType) continue; // filter here ✅

                int rp_index = constituents[c];
                const auto &recoPart = recoParticles[rp_index];

                // Link to Track
                for (int t = recoPart.tracks_begin; t < recoPart.tracks_end; ++t) {
                    int track_index = _recoParticlesIndices[t];

                    if (track_index_to_dEdx.count(track_index)) {
                        jet_dEdx_masked.push_back(track_index_to_dEdx[track_index]);
                    }
                }
            }

            dedx_filtered.emplace_back(std::move(jet_dEdx_masked));
        }

        return dedx_filtered;
    }
};








}} // namespace FCCAnalyses::AlephSelection

#endif
