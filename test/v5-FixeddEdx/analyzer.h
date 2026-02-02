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

// -----------------------------------
// Calculate the dEdx variables for 
// pads and wires with preserved vector length
// -----------------------------------
struct build_constituents_dEdx_preserved {
    struct dEdxResult {
        rv::RVec<rv::RVec<float>> values;
        rv::RVec<rv::RVec<float>> errors;
    };
    
    dEdxResult
    operator()(const rv::RVec<FCCAnalysesJetConstituentsData> &jet_constituents_types,
               const rv::RVec<edm4hep::ReconstructedParticleData> &recoParticles,
               const rv::RVec<int> &_recoParticlesIndices,
               const rv::RVec<edm4hep::RecDqdxData> &dEdxCollection,
               const rv::RVec<int> &_dEdxIndicesCollection,
               const std::vector<std::vector<int>> &jet_indices,
               float targetType) const {
        
        rv::RVec<rv::RVec<float>> values;
        rv::RVec<rv::RVec<float>> errors;
        
        // Build map Track.index -> dEdx object
        std::unordered_map<int, edm4hep::RecDqdxData> track_index_to_dEdx;
        for (size_t i = 0; i < _dEdxIndicesCollection.size(); ++i) {
            track_index_to_dEdx[_dEdxIndicesCollection[i]] = dEdxCollection[i];
        }
        
        // Loop over jets
        for (size_t j = 0; j < jet_indices.size(); ++j) {
            const auto &constituents = jet_indices[j];
            const auto &types = jet_constituents_types[j];
            
            rv::RVec<float> jet_values;
            rv::RVec<float> jet_errors;
            jet_values.reserve(constituents.size());
            jet_errors.reserve(constituents.size());
            
            // Loop over jet constituents
            for (size_t c = 0; c < constituents.size(); ++c) {
                int rp_index = constituents[c];
                const auto &recoPart = recoParticles[rp_index];
                
                // Check if this constituent matches the target type
                if (types[c] == targetType) {
                    // Try to find dEdx for this particle
                    bool found = false;
                    for (int t = recoPart.tracks_begin; t < recoPart.tracks_end; ++t) {
                        int track_index = _recoParticlesIndices[t];
                        if (track_index_to_dEdx.count(track_index)) {
                            const auto &dEdx = track_index_to_dEdx[track_index];
                            if (dEdx.dQdx.type == 0) {
                                jet_values.push_back(dEdx.dQdx.value);
                                jet_errors.push_back(dEdx.dQdx.error);
                            } else {
                                jet_values.push_back(-9.0f);
                                jet_errors.push_back(-9.0f);
                            }
                            found = true;
                            break;
                        }
                    }
                    if (!found) {
                        // Target type but no dEdx found
                        jet_values.push_back(-9.0f);
                        jet_errors.push_back(-9.0f);
                    }
                } else {
                    // Not target type - use placeholder
                    jet_values.push_back(-9.0f);
                    jet_errors.push_back(-9.0f);
                }
            }
            
            values.emplace_back(std::move(jet_values));
            errors.emplace_back(std::move(jet_errors));
        }
        
        return {values, errors};
    }
};


struct build_constituents_dEdx_by_mass {
    struct dEdxResult {
        rv::RVec<rv::RVec<float>> pion_values;
        rv::RVec<rv::RVec<float>> pion_errors;
        rv::RVec<rv::RVec<float>> kaon_values;
        rv::RVec<rv::RVec<float>> kaon_errors;
        rv::RVec<rv::RVec<float>> proton_values;
        rv::RVec<rv::RVec<float>> proton_errors;
    };
    
    // Particle masses in GeV
    static constexpr float PION_MASS = 0.13957039f;
    static constexpr float KAON_MASS = 0.493677f;
    static constexpr float PROTON_MASS = 0.938272f;
    static constexpr float MASS_TOLERANCE = 0.010f;  // 10 MeV tolerance
    
    enum ParticleType {
        PION,
        KAON,
        PROTON,
        OTHER
    };
    
    ParticleType identifyParticle(float mass) const {
        if (std::abs(mass - PION_MASS) < MASS_TOLERANCE) return PION;
        if (std::abs(mass - KAON_MASS) < MASS_TOLERANCE) return KAON;
        if (std::abs(mass - PROTON_MASS) < MASS_TOLERANCE) return PROTON;
        return OTHER;
    }
    
    dEdxResult
    operator()(const rv::RVec<FCCAnalysesJetConstituentsData> &jet_constituents_types,
               const rv::RVec<edm4hep::ReconstructedParticleData> &recoParticles,
               const rv::RVec<int> &_recoParticlesIndices,
               const rv::RVec<edm4hep::RecDqdxData> &dEdxCollection,
               const rv::RVec<int> &_dEdxIndicesCollection,
               const std::vector<std::vector<int>> &jet_indices) const {
        
        rv::RVec<rv::RVec<float>> pion_values;
        rv::RVec<rv::RVec<float>> pion_errors;
        rv::RVec<rv::RVec<float>> kaon_values;
        rv::RVec<rv::RVec<float>> kaon_errors;
        rv::RVec<rv::RVec<float>> proton_values;
        rv::RVec<rv::RVec<float>> proton_errors;
        
        // Build map Track.index -> dEdx object
        std::unordered_map<int, edm4hep::RecDqdxData> track_index_to_dEdx;
        for (size_t i = 0; i < _dEdxIndicesCollection.size(); ++i) {
            track_index_to_dEdx[_dEdxIndicesCollection[i]] = dEdxCollection[i];
        }
        
        // Loop over jets
        for (size_t j = 0; j < jet_indices.size(); ++j) {
            const auto &constituents = jet_indices[j];
            
            rv::RVec<float> jet_pion_values;
            rv::RVec<float> jet_pion_errors;
            rv::RVec<float> jet_kaon_values;
            rv::RVec<float> jet_kaon_errors;
            rv::RVec<float> jet_proton_values;
            rv::RVec<float> jet_proton_errors;
            
            jet_pion_values.reserve(constituents.size());
            jet_pion_errors.reserve(constituents.size());
            jet_kaon_values.reserve(constituents.size());
            jet_kaon_errors.reserve(constituents.size());
            jet_proton_values.reserve(constituents.size());
            jet_proton_errors.reserve(constituents.size());
            
            // Loop over jet constituents
            for (size_t c = 0; c < constituents.size(); ++c) {
                int rp_index = constituents[c];
                const auto &recoPart = recoParticles[rp_index];
                
                // Identify particle type by mass
                ParticleType particle_type = identifyParticle(recoPart.mass);
                
                float pion_value = -9.0f, pion_error = -9.0f;
                float kaon_value = -9.0f, kaon_error = -9.0f;
                float proton_value = -9.0f, proton_error = -9.0f;
                
                if (particle_type != OTHER) {
                    // Try to find dEdx for this particle
                    float value = -9.0f;
                    float error = -9.0f;
                    
                    for (int t = recoPart.tracks_begin; t < recoPart.tracks_end; ++t) {
                        int track_index = _recoParticlesIndices[t];
                        if (track_index_to_dEdx.count(track_index)) {
                            const auto &dEdx = track_index_to_dEdx[track_index];
                            if (dEdx.dQdx.type == 0) {
                                value = dEdx.dQdx.value;
                                error = dEdx.dQdx.error;
                            }
                            break;
                        }
                    }
                    
                    // Assign to appropriate species
                    switch (particle_type) {
                        case PION:
                            pion_value = value;
                            pion_error = error;
                            break;
                        case KAON:
                            kaon_value = value;
                            kaon_error = error;
                            break;
                        case PROTON:
                            proton_value = value;
                            proton_error = error;
                            break;
                        default:
                            break;
                    }
                }
                // For OTHER particles, all remain at -9.0f
                
                // Add values for all particle types (with placeholders for non-matches)
                jet_pion_values.push_back(pion_value);
                jet_pion_errors.push_back(pion_error);
                jet_kaon_values.push_back(kaon_value);
                jet_kaon_errors.push_back(kaon_error);
                jet_proton_values.push_back(proton_value);
                jet_proton_errors.push_back(proton_error);
            }
            
            pion_values.emplace_back(std::move(jet_pion_values));
            pion_errors.emplace_back(std::move(jet_pion_errors));
            kaon_values.emplace_back(std::move(jet_kaon_values));
            kaon_errors.emplace_back(std::move(jet_kaon_errors));
            proton_values.emplace_back(std::move(jet_proton_values));
            proton_errors.emplace_back(std::move(jet_proton_errors));
        }
        
        return {pion_values, pion_errors, kaon_values, kaon_errors, proton_values, proton_errors};
    }
};

}} // namespace FCCAnalyses::AlephSelection

#endif
