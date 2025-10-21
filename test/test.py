processList = {

"QQB" : {"fraction" : 0.001},           
}


outputDir = "./output"
inputDir = "/eos/experiment/aleph/EDM4HEP/MC/1994"
nCPUS = -1
includePaths = ["analyzer.h"]


class RDFanalysis:
    def analysers(df):



        coll = {
        "GenParticles": "MCParticles",
        "PFParticles": "RecoParticles",
        "PFTracks": "EFlowTrack",
        "PFPhotons": "EFlowPhoton",
        "PFNeutralHadrons": "EFlowNeutralHadron",
        "TrackState": "_Tracks_trackStates",
        "TrackerHits": "TrackerHits",
        "CalorimeterHits": "CalorimeterHits",
        "PathLength": "EFlowTrack_L",
        "Bz": "magFieldBz",
        }


        df = df.Filter("AlephSelection::sel_class_filter(16)(ClassBitset)")


        # Define RP kinematics
        ####################################################################################################
        df = df.Define("RP_px", "ReconstructedParticle::get_px(RecoParticles)")
        df = df.Define("RP_py", "ReconstructedParticle::get_py(RecoParticles)")
        df = df.Define("RP_pz", "ReconstructedParticle::get_pz(RecoParticles)")
        df = df.Define("RP_e", "ReconstructedParticle::get_e(RecoParticles)")
        df = df.Define("RP_m", "ReconstructedParticle::get_mass(RecoParticles)")

        # Define pseudo-jets
        ####################################################################################################
        df = df.Define("pjetc", "JetClusteringUtils::set_pseudoJets(RP_px, RP_py, RP_pz, RP_e)")
        # Anti-kt clustering and jet constituents
        ####################################################################################################
        df = df.Define("_jet", "JetClustering::clustering_ee_kt(2, 2, 1, 0)(pjetc)")
        df = df.Define("jets","JetClusteringUtils::get_pseudoJets(_jet)" )
        df = df.Define("_jetc", "JetClusteringUtils::get_constituents(_jet)") 
        df = df.Define("jetc", "JetConstituentsUtils::build_constituents_cluster(RecoParticles, _jetc)")
        df = df.Define("jetConstitutentsTypes", f"AlephSelection::build_constituents_Types()(ParticleID, _jetc)")


        ############################################# Event Level Variables #######################################################
        df = df.Define("jet_p4", "JetConstituentsUtils::compute_tlv_jets(jets)" )
        df = df.Define("event_invariant_mass", "JetConstituentsUtils::InvariantMass(jet_p4[0], jet_p4[1])")
        df = df.Define("event_type", "AlephSelection::get_EventType({})".format(coll["GenParticles"]))

         ########## Picking Jet Flavors

        #df.Filter("event_type == 2") # d-quark: 1, u-quark:2, s-quark:3, c-quark:4, b-quark: 5

        # ===== VERTEX
        # MC primary vertex
        df = df.Define("pv", f'AlephSelection::get_EventPrimaryVertexP4()({coll["GenParticles"]})')


        ############################################# Particle Flow Level Variables #######################################################


        df = df.Define("pfcand_isMu",     "AlephSelection::get_isType(jetConstitutentsTypes,2)")
        df = df.Define("pfcand_isEl",     "AlephSelection::get_isType(jetConstitutentsTypes,1)")
        df = df.Define("pfcand_isGamma",  "AlephSelection::get_isType(jetConstitutentsTypes,4)")
        df = df.Define("pfcand_isChargedHad", "AlephSelection::get_isType(jetConstitutentsTypes,0)")
        df = df.Define("pfcand_isNeutralHad", "AlephSelection::get_isType(jetConstitutentsTypes,5)")


        ############################################# Kinematics and PID #######################################################

        df = df.Define("pfcand_e",        "JetConstituentsUtils::get_e(jetc)") 
        df = df.Define("pfcand_p",        "JetConstituentsUtils::get_p(jetc)") 
        df = df.Define("pfcand_theta",    "JetConstituentsUtils::get_theta(jetc)") 
        df = df.Define("pfcand_phi",      "JetConstituentsUtils::get_phi(jetc)") 
        df = df.Define("pfcand_charge",   "JetConstituentsUtils::get_charge(jetc)") 
        df = df.Define("pfcand_type",     "JetConstituentsUtils::get_type(jetc)") 
        df = df.Define("pfcand_erel",     "JetConstituentsUtils::get_erel_cluster(jets, jetc)")
        df = df.Define("pfcand_erel_log", "JetConstituentsUtils::get_erel_log_cluster(jets, jetc)")
        df = df.Define("pfcand_thetarel", "JetConstituentsUtils::get_thetarel_cluster(jets, jetc)")
        df = df.Define("pfcand_phirel",   "JetConstituentsUtils::get_phirel_cluster(jets, jetc)")

        df = df.Define("Bz", '1.5')


############################################# Track Parameters and Covariance #######################################################

        df = df.Define("pfcand_dxy",        f'JetConstituentsUtils::XPtoPar_dxy(jetc, {coll["TrackState"]}, pv, Bz)') 
        df = df.Define("pfcand_dz",         f'JetConstituentsUtils::XPtoPar_dz(jetc, {coll["TrackState"]}, pv, Bz)') 
        df = df.Define("pfcand_phi0",       f'JetConstituentsUtils::XPtoPar_phi(jetc, {coll["TrackState"]}, pv, Bz)') 
        df = df.Define("pfcand_C",          f'JetConstituentsUtils::XPtoPar_C(jetc, {coll["TrackState"]}, Bz)') 
        df = df.Define("pfcand_ct",         f'JetConstituentsUtils::XPtoPar_ct(jetc, {coll["TrackState"]}, Bz)') 
        df = df.Define("pfcand_dptdpt",     f'JetConstituentsUtils::get_omega_cov(jetc, {coll["TrackState"]})') 
        df = df.Define("pfcand_dxydxy",     f'JetConstituentsUtils::get_d0_cov(jetc, {coll["TrackState"]})') 
        df = df.Define("pfcand_dzdz",       f'JetConstituentsUtils::get_z0_cov(jetc, {coll["TrackState"]})') 
        df = df.Define("pfcand_dphidphi",   f'JetConstituentsUtils::get_phi0_cov(jetc, {coll["TrackState"]})') 
        df = df.Define("pfcand_detadeta",   f'JetConstituentsUtils::get_tanlambda_cov(jetc, {coll["TrackState"]})') 
        df = df.Define("pfcand_dxydz",      f'JetConstituentsUtils::get_d0_z0_cov(jetc, {coll["TrackState"]})') 
        df = df.Define("pfcand_dphidxy",    f'JetConstituentsUtils::get_phi0_d0_cov(jetc, {coll["TrackState"]})') 
        df = df.Define("pfcand_phidz",      f'JetConstituentsUtils::get_phi0_z0_cov(jetc, {coll["TrackState"]})') 
        df = df.Define("pfcand_phictgtheta",f'JetConstituentsUtils::get_tanlambda_phi0_cov(jetc, {coll["TrackState"]})') 
        df = df.Define("pfcand_dxyctgtheta",f'JetConstituentsUtils::get_tanlambda_d0_cov(jetc, {coll["TrackState"]})') 
        df = df.Define("pfcand_dlambdadz",  f'JetConstituentsUtils::get_tanlambda_z0_cov(jetc, {coll["TrackState"]})') 
        df = df.Define("pfcand_cctgtheta",  f'JetConstituentsUtils::get_omega_tanlambda_cov(jetc, {coll["TrackState"]})') 
        df = df.Define("pfcand_phic",       f'JetConstituentsUtils::get_omega_phi0_cov(jetc, {coll["TrackState"]})') 
        df = df.Define("pfcand_dxyc",       f'JetConstituentsUtils::get_omega_d0_cov(jetc, {coll["TrackState"]})') 
        df = df.Define("pfcand_cdz",        f'JetConstituentsUtils::get_omega_z0_cov(jetc, {coll["TrackState"]})')



############################################# Btag Variables #######################################################

        df = df.Define("pfcand_btagSip2dVal",   "JetConstituentsUtils::get_Sip2dVal_clusterV(jets, pfcand_dxy, pfcand_phi0, Bz)") 
        df = df.Define("pfcand_btagSip2dSig",   "JetConstituentsUtils::get_Sip2dSig(pfcand_btagSip2dVal, pfcand_dxydxy)") 
        df = df.Define("pfcand_btagSip3dVal",   "JetConstituentsUtils::get_Sip3dVal_clusterV(jets, pfcand_dxy, pfcand_dz, pfcand_phi0, Bz)") 
        df = df.Define("pfcand_btagSip3dSig",   "JetConstituentsUtils::get_Sip3dSig(pfcand_btagSip3dVal, pfcand_dxydxy, pfcand_dzdz)") 
        df = df.Define("pfcand_btagJetDistVal","JetConstituentsUtils::get_JetDistVal_clusterV(jets, jetc, pfcand_dxy, pfcand_dz, pfcand_phi0, Bz)") 
        df = df.Define("pfcand_btagJetDistSig","JetConstituentsUtils::get_JetDistSig(pfcand_btagJetDistVal, pfcand_dxydxy, pfcand_dzdz)")





        ############################################# Jet Level Variables #######################################################
        df=df.Define("event_njet",   "JetConstituentsUtils::count_jets(jetc)")
        df.Filter("event_njet > 1")
        ##############################################################################################################
        df = df.Define("sumTLVs", "JetConstituentsUtils::sum_tlv_constituents(jetc)")

        df = df.Define("jet_p", "ROOT::VecOps::RVec<Double_t>({sumTLVs[0].P(), sumTLVs[1].P()})")
        df = df.Define("jet_e", "ROOT::VecOps::RVec<Double_t>({sumTLVs[0].E(), sumTLVs[1].E()})")
        df = df.Define("jet_mass", "ROOT::VecOps::RVec<Double_t>({sumTLVs[0].M(), sumTLVs[1].M()})")
        df = df.Define("jet_phi", "ROOT::VecOps::RVec<Double_t>({sumTLVs[0].Phi(), sumTLVs[1].Phi()})")
        df = df.Define("jet_theta", "ROOT::VecOps::RVec<Double_t>({sumTLVs[0].Theta(), sumTLVs[1].Theta()})")
        df = df.Define("jet_pT", "ROOT::VecOps::RVec<Double_t>({sumTLVs[0].Pt(), sumTLVs[1].Pt()})")
        df = df.Define("jet_eta", "ROOT::VecOps::RVec<Double_t>({sumTLVs[0].Eta(), sumTLVs[1].Eta()})")
        # Leading jet
        df = df.Define("jet_p_leading",      "sumTLVs[0].P()")
        df = df.Define("jet_e_leading",      "sumTLVs[0].E()")
        df = df.Define("jet_mass_leading",   "sumTLVs[0].M()")
        df = df.Define("jet_phi_leading",    "sumTLVs[0].Phi()")
        df = df.Define("jet_theta_leading",  "sumTLVs[0].Theta()")
        df = df.Define("jet_pT_leading",     "sumTLVs[0].Pt()")
        df = df.Define("jet_eta_leading",    "sumTLVs[0].Eta()")
        
        # Subleading jet
        df = df.Define("jet_p_subleading",      "sumTLVs[1].P()")
        df = df.Define("jet_e_subleading",      "sumTLVs[1].E()")
        df = df.Define("jet_mass_subleading",   "sumTLVs[1].M()")
        df = df.Define("jet_phi_subleading",    "sumTLVs[1].Phi()")
        df = df.Define("jet_theta_subleading",  "sumTLVs[1].Theta()")
        df = df.Define("jet_pT_subleading",     "sumTLVs[1].Pt()")
        df = df.Define("jet_eta_subleading",    "sumTLVs[1].Eta()")


        df = df.Define("jet_nconst", "JetConstituentsUtils::count_consts(jetc)") 
        ##
        df = df.Define(f"jet_nmu",    f"JetConstituentsUtils::count_type(pfcand_isMu)") 
        df = df.Define(f"jet_nel",    f"JetConstituentsUtils::count_type(pfcand_isEl)") 
        df = df.Define(f"jet_nchad",  f"JetConstituentsUtils::count_type(pfcand_isChargedHad)") 
        df = df.Define(f"jet_ngamma", f"JetConstituentsUtils::count_type(pfcand_isGamma)") 
        df = df.Define(f"jet_nnhad",  f"JetConstituentsUtils::count_type(pfcand_isNeutralHad)")

        df = df.Define("dEdxPadsValue" , "dEdxPads.dQdx.value")
        df = df.Define("dEdxPadsError" , "dEdxPads.dQdx.error")
        df = df.Define("dEdxWiresValue" , "dEdxWires.dQdx.value")
        df = df.Define("dEdxWiresError" , "dEdxPads.dQdx.error")



        return df

    def output():

        return [
            "event_type",
            "event_invariant_mass","event_njet",  
            #"jet_mass","jet_p","jet_e", "jet_phi", "jet_theta", "jet_pT",
            "dEdxPadsValue", "dEdxPadsError", "dEdxWiresValue", "dEdxWiresError",
                "jet_p_leading",
                "jet_e_leading",
                "jet_mass_leading",
                "jet_phi_leading",
                "jet_theta_leading",
                "jet_pT_leading",
                "jet_eta_leading",
                "jet_p_subleading",
                "jet_e_subleading",
                "jet_mass_subleading",
                "jet_phi_subleading",
                "jet_theta_subleading",
                "jet_pT_subleading",
                "jet_eta_subleading",   
                
                "jet_nnhad","jet_ngamma","jet_nchad","jet_nel", "jet_nmu", "jet_nconst",
                
                "pfcand_isMu", "pfcand_isEl", "pfcand_isChargedHad", "pfcand_isGamma", "pfcand_isNeutralHad",
                "pfcand_e", "pfcand_p", "pfcand_theta", "pfcand_phi", "pfcand_charge", "pfcand_type",
                "pfcand_erel", "pfcand_erel_log", "pfcand_thetarel", "pfcand_phirel", 
 
                #"Bz",
 
                "pfcand_dxy", "pfcand_dz", "pfcand_phi0", "pfcand_C", "pfcand_ct",
                "pfcand_dptdpt", "pfcand_dxydxy", "pfcand_dzdz", "pfcand_dphidphi", "pfcand_detadeta",
                "pfcand_dxydz", "pfcand_dphidxy", "pfcand_phidz", "pfcand_phictgtheta", "pfcand_dxyctgtheta",
                "pfcand_dlambdadz", "pfcand_cctgtheta", "pfcand_phic", "pfcand_dxyc", "pfcand_cdz",
 
                "pfcand_btagSip2dVal", "pfcand_btagSip2dSig", "pfcand_btagSip3dVal", "pfcand_btagSip3dSig", 
                "pfcand_btagJetDistVal", "pfcand_btagJetDistSig",
                ]
