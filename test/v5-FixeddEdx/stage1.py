from argparse import ArgumentParser

class Analysis():

    def __init__(self, cmdline_args):
        parser = ArgumentParser(
            description='Additional analysis arguments',
            usage='Provide additional arguments after analysis script path')
        parser.add_argument('--tag', required=True, type=str,
                            help='Production tag to indicate version.')
        parser.add_argument('--doData', action='store_true',
                            help='Run on data, instead of MC (which is the default behaviour).')
        parser.add_argument('--year', default='1994',
                            help='MC/data year to run on - currently only 1994 as option.')
        parser.add_argument('--MCtype', default="zqq", type=str,
                            help='Type of MC to run on - currently only zqq as option.')
        parser.add_argument('--MCflavour', default=None, type=str,
                            help='For MC only: filter out events based on truth quark flavours. Default is none. Options: \
                            1 = dd, 2 = uu, 3 = ss, 4 = cc, 5 = bb')
        parser.add_argument('--fraction', default=1.0, type=float,
                            help='Fraction of events to run, default is 1.0 = 100%')
        # Parse additional arguments not known to the FCCAnalyses parsers
        # All command line arguments know to fccanalysis are provided in the
        # `cmdline_arg` dictionary.
        self.ana_args, _ = parser.parse_known_args(cmdline_args['remaining'])

        #Dictionary for setting output names:
        outnames_dict = {
            # proc: {flavour_id_1:{flavour_name_1}, flavour_id_2:{flavour_name_2}, ..}
            "zqq":{
                "1":"Zdd",
                "2":"Zuu",
                "3":"Zss",
                "4":"Zcc",
                "5":"Zbb",
                }
        }

        # sanity checks for the command line arguments:
        if self.ana_args.doData and self.ana_args.MCtype:
            print("----> WARNING: Incompatible input arguments: --MCtype defined with --doData, will be ignored.")

        if self.ana_args.doData and self.ana_args.MCflavour:
            print("----> WARNING: Incompatible input arguments: --MCflavour defined with --doData, will be ignored.")

        if self.ana_args.MCflavour and not self.ana_args.MCtype:
            print("----> ERROR: Requested truth flavour filter with --MCflavour without specifying --MCtype.")
            exit()
        
        if self.ana_args.MCtype and not self.ana_args.MCtype in outnames_dict:
            print("----> ERROR: Requested unknown --MCtype. Currently only zqq available.")
            exit()
        
        if not self.ana_args.doData and not self.ana_args.MCflavour:
            print(f"----> ERROR: Requested MC run but did not specify --MCflavour. Please pick one..")
            exit()
        
        if self.ana_args.MCflavour and not self.ana_args.MCflavour in outnames_dict[self.ana_args.MCtype]:
            print(f"----> ERROR: Requested unknown --MCflavour for --MCtype {self.ana_args.MCtype}. Check the dictionary.")
            exit()

        #set the input/output directories:
        if self.ana_args.doData:
            self.input_dir = "/eos/experiment/fcc/ee/analyses/case-studies/aleph/LEP1_DATA/"
            self.output_dir = f"/eos/experiment/fcc/ee/analyses/case-studies/aleph/processedData/{self.ana_args.year}/stage1/{self.ana_args.tag}"
            
            self.process_list = {
                "1994" : {"fraction" : self.ana_args.fraction},           
            }  

        else:
            self.input_dir = f"/eos/experiment/aleph/EDM4HEP/MC/{self.ana_args.year}/"
            self.output_dir = f"/eos/experiment/fcc/ee/analyses/case-studies/aleph/processedMC/{self.ana_args.year}/{self.ana_args.MCtype}/stage1/{self.ana_args.tag}"

            #set the output file name depending on resonance flavour 
            output_name = outnames_dict[self.ana_args.MCtype][self.ana_args.MCflavour]

            self.process_list = {
                "QQB" : {"fraction" : self.ana_args.fraction, "output":output_name},           
            }

        #set run options:
        self.n_threads = 64 
        self.include_paths = ["analyzer.h"]

    def analyzers(self, df):

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

        if self.ana_args.doData:
            #df = df.Filter("AlephSelection::sel_class_filter(16)(ClassBitset)   || AlephSelection::sel_class_filter(17)(ClassBitset) ")
            df = df.Filter("AlephSelection::sel_class_filter(16)(ClassBitset) ")
            df = df.Define("jetPID", "-999")
        else:
            # Using Classbit to filter out QQbar samples and then get a specific flavor of jets
            # d-quark: 1, u-quark:2, s-quark:3, c-quark:4, b-quark: 5
            df = df.Define("jetPID", f"AlephSelection::getJetPID(ClassBitset, {coll['GenParticles']})")
            df = df.Filter(f"jetPID == {self.ana_args.MCflavour}")
        
        # store the classbitset in the output
        df = df.Define("event_class", "AlephSelection::bitsetToIndices(ClassBitset)")

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

        # ===== VERTEX
        df = df.Define(
            "pv",
            "TLorentzVector(Vertices[0].position.x, Vertices[0].position.y, Vertices[0].position.z, 0.0)",
        )
        df = df.Define("VertexX", "Vertices.position.x")
        df = df.Define("VertexY", "Vertices.position.y")
        df = df.Define("VertexZ", "Vertices.position.z")

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

        df = df.Define("TrackStateFlipped",f"AlephSelection::flipD0_copy( {coll['TrackState']} )")

        df = df.Define("pfcand_dxy",        f'JetConstituentsUtils::XPtoPar_dxy(jetc, TrackStateFlipped, pv, Bz)') 
        df = df.Define("pfcand_dz",         f'JetConstituentsUtils::XPtoPar_dz(jetc, TrackStateFlipped, pv, Bz)') 
        df = df.Define("pfcand_phi0",       f'JetConstituentsUtils::XPtoPar_phi(jetc, TrackStateFlipped, pv, Bz)') 
        df = df.Define("pfcand_C",          f'JetConstituentsUtils::XPtoPar_C(jetc, TrackStateFlipped, Bz)') 
        df = df.Define("pfcand_ct",         f'JetConstituentsUtils::XPtoPar_ct(jetc, TrackStateFlipped, Bz)') 
        df = df.Define("pfcand_dptdpt",     f'JetConstituentsUtils::get_omega_cov(jetc, TrackStateFlipped)') 
        df = df.Define("pfcand_dxydxy",     f'JetConstituentsUtils::get_d0_cov(jetc, TrackStateFlipped)') 
        df = df.Define("pfcand_dzdz",       f'JetConstituentsUtils::get_z0_cov(jetc, TrackStateFlipped)') 
        df = df.Define("pfcand_dphidphi",   f'JetConstituentsUtils::get_phi0_cov(jetc, TrackStateFlipped)') 
        df = df.Define("pfcand_detadeta",   f'JetConstituentsUtils::get_tanlambda_cov(jetc, TrackStateFlipped)') 
        df = df.Define("pfcand_dxydz",      f'JetConstituentsUtils::get_d0_z0_cov(jetc, TrackStateFlipped)') 
        df = df.Define("pfcand_dphidxy",    f'JetConstituentsUtils::get_phi0_d0_cov(jetc, TrackStateFlipped)') 
        df = df.Define("pfcand_phidz",      f'JetConstituentsUtils::get_phi0_z0_cov(jetc, TrackStateFlipped)') 
        df = df.Define("pfcand_phictgtheta",f'JetConstituentsUtils::get_tanlambda_phi0_cov(jetc, TrackStateFlipped)') 
        df = df.Define("pfcand_dxyctgtheta",f'JetConstituentsUtils::get_tanlambda_d0_cov(jetc, TrackStateFlipped)') 
        df = df.Define("pfcand_dlambdadz",  f'JetConstituentsUtils::get_tanlambda_z0_cov(jetc, TrackStateFlipped)') 
        df = df.Define("pfcand_cctgtheta",  f'JetConstituentsUtils::get_omega_tanlambda_cov(jetc, TrackStateFlipped)') 
        df = df.Define("pfcand_phic",       f'JetConstituentsUtils::get_omega_phi0_cov(jetc, TrackStateFlipped)') 
        df = df.Define("pfcand_dxyc",       f'JetConstituentsUtils::get_omega_d0_cov(jetc, TrackStateFlipped)') 
        df = df.Define("pfcand_cdz",        f'JetConstituentsUtils::get_omega_z0_cov(jetc, TrackStateFlipped)')

        ############################################# Btag Variables #######################################################

        df = df.Define("pfcand_btagSip2dVal",   "JetConstituentsUtils::get_Sip2dVal_clusterV(jets, pfcand_dxy, pfcand_phi0, Bz)") 
        df = df.Define("pfcand_btagSip2dSig",   "JetConstituentsUtils::get_Sip2dSig(pfcand_btagSip2dVal, pfcand_dxydxy)") 
        df = df.Define("pfcand_btagSip3dVal",   "JetConstituentsUtils::get_Sip3dVal_clusterV(jets, pfcand_dxy, pfcand_dz, pfcand_phi0, Bz)") 
        df = df.Define("pfcand_btagSip3dSig",   "JetConstituentsUtils::get_Sip3dSig(pfcand_btagSip3dVal, pfcand_dxydxy, pfcand_dzdz)") 
        df = df.Define("pfcand_btagJetDistVal","JetConstituentsUtils::get_JetDistVal_clusterV(jets, jetc, pfcand_dxy, pfcand_dz, pfcand_phi0, Bz)") 
        df = df.Define("pfcand_btagJetDistSig","JetConstituentsUtils::get_JetDistSig(pfcand_btagJetDistVal, pfcand_dxydxy, pfcand_dzdz)")


        ############################################# Jet Level Variables and selection #######################################################
        
        df=df.Define("event_njet",   "JetConstituentsUtils::count_jets(jetc)")
        df = df.Filter("event_njet > 1")

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



        #         ## ELECTRONS — Pads
        # df = df.Define("pfcand_dEdx_pads_el_result",
        #             "AlephSelection::build_constituents_dEdx_preserved()("
        #             "jetConstitutentsTypes,"
        #             "RecoParticles, _RecoParticles_tracks.index,"
        #             "dEdxPads, _dEdxPads_track.index,"
        #             "_jetc, 1)")
        # df = df.Define("pfcand_dEdx_el_value_pads", "pfcand_dEdx_pads_el_result.values")
        # df = df.Define("pfcand_dEdx_el_error_pads", "pfcand_dEdx_pads_el_result.errors")

        # ## ELECTRONS — Wires
        # df = df.Define("pfcand_dEdx_wires_el_result",
        #             "AlephSelection::build_constituents_dEdx_preserved()("
        #             "jetConstitutentsTypes,"
        #             "RecoParticles, _RecoParticles_tracks.index,"
        #             "dEdxWires, _dEdxWires_track.index,"
        #             "_jetc, 1)")
        # df = df.Define("pfcand_dEdx_el_wires_value", "pfcand_dEdx_wires_el_result.values")
        # df = df.Define("pfcand_dEdx_el_wires_error", "pfcand_dEdx_wires_el_result.errors")

        # ## MUONS — Pads
        # df = df.Define("pfcand_dEdx_pads_mu_result",
        #             "AlephSelection::build_constituents_dEdx_preserved()("
        #             "jetConstitutentsTypes,"
        #             "RecoParticles, _RecoParticles_tracks.index,"
        #             "dEdxPads, _dEdxPads_track.index,"
        #             "_jetc, 2)")
        # df = df.Define("pfcand_dEdx_mu_value_pads", "pfcand_dEdx_pads_mu_result.values")
        # df = df.Define("pfcand_dEdx_mu_error_pads", "pfcand_dEdx_pads_mu_result.errors")

        # ## MUONS — Wires
        # df = df.Define("pfcand_dEdx_wires_mu_result",
        #             "AlephSelection::build_constituents_dEdx_preserved()("
        #             "jetConstitutentsTypes,"
        #             "RecoParticles, _RecoParticles_tracks.index,"
        #             "dEdxWires, _dEdxWires_track.index,"
        #             "_jetc, 2)")
        # df = df.Define("pfcand_dEdx_mu_wires_value", "pfcand_dEdx_wires_mu_result.values")
        # df = df.Define("pfcand_dEdx_mu_wires_error", "pfcand_dEdx_wires_mu_result.errors")

        # ## CHARGED HADRONS — Pads (π/K/p)
        # df = df.Define("pfcand_dEdx_pads_ch_result",
        #             "AlephSelection::build_constituents_dEdx_preserved()("
        #             "jetConstitutentsTypes,"
        #             "RecoParticles, _RecoParticles_tracks.index,"
        #             "dEdxPads, _dEdxPads_track.index,"
        #             "_jetc, 0)")
        # df = df.Define("pfcand_dEdx_ch_value_pads", "pfcand_dEdx_pads_ch_result.values")
        # df = df.Define("pfcand_dEdx_ch_error_pads", "pfcand_dEdx_pads_ch_result.errors")

        # ## CHARGED HADRONS — Wires
        # df = df.Define("pfcand_dEdx_wires_ch_result",
        #             "AlephSelection::build_constituents_dEdx_preserved()("
        #             "jetConstitutentsTypes,"
        #             "RecoParticles, _RecoParticles_tracks.index,"
        #             "dEdxWires, _dEdxWires_track.index,"
        #             "_jetc, 0)")
        # df = df.Define("pfcand_dEdx_ch_wires_value", "pfcand_dEdx_wires_ch_result.values")
        # df = df.Define("pfcand_dEdx_ch_wires_error", "pfcand_dEdx_wires_ch_result.errors")

        # Main result struct
        df = df.Define("pfcand_dEdx_pads_result",
                        "AlephSelection::build_constituents_dEdx_by_mass()("
                        "jetConstitutentsTypes,"
                        "RecoParticles, _RecoParticles_tracks.index,"
                        "dEdxPads, _dEdxPads_track.index,"
                        "_jetc)")
        
        # Pion variables
        df = df.Define("pfcand_dEdx_pion_value_pads", "pfcand_dEdx_pads_result.pion_values")
        df = df.Define("pfcand_dEdx_pion_error_pads", "pfcand_dEdx_pads_result.pion_errors")
        
        # Kaon variables
        df = df.Define("pfcand_dEdx_kaon_value_pads", "pfcand_dEdx_pads_result.kaon_values")
        df = df.Define("pfcand_dEdx_kaon_error_pads", "pfcand_dEdx_pads_result.kaon_errors")
        
        # Proton variables
        df = df.Define("pfcand_dEdx_proton_value_pads", "pfcand_dEdx_pads_result.proton_values")
        df = df.Define("pfcand_dEdx_proton_error_pads", "pfcand_dEdx_pads_result.proton_errors")

        df = df.Filter("jet_nconst[0] > 2")
        df = df.Filter("jet_nconst[1] > 2")

        df = df.Filter("jet_nchad[0] > 2")
        df = df.Filter("jet_nchad[1] > 2")

        df = df.Filter("jet_pT[0] > 10")
        df = df.Filter("jet_pT[1] > 10")
        
        df = df.Filter("abs(jet_theta[0]) > 0.75 && abs(jet_theta[0]) < 2.4")
        df = df.Filter("abs(jet_theta[1]) > 0.75 && abs(jet_theta[1]) < 2.4")


        return df

    def output(self):

        return [

 # Pads
    # "pfcand_dEdx_el_value_pads",
    # "pfcand_dEdx_mu_value_pads",
    # "pfcand_dEdx_ch_value_pads",

    # "pfcand_dEdx_el_error_pads",
    # "pfcand_dEdx_mu_error_pads",
    # "pfcand_dEdx_ch_error_pads",


    
    # # Wires
    # "pfcand_dEdx_el_wires_value",
    # "pfcand_dEdx_mu_wires_value",
    # "pfcand_dEdx_ch_wires_value",

    # "pfcand_dEdx_el_wires_error",
    # "pfcand_dEdx_mu_wires_error",
    # "pfcand_dEdx_ch_wires_error",
    "pfcand_dEdx_pion_value_pads",
    "pfcand_dEdx_pion_error_pads",
    "pfcand_dEdx_kaon_value_pads",
    "pfcand_dEdx_kaon_error_pads",
    "pfcand_dEdx_proton_value_pads",
    "pfcand_dEdx_proton_error_pads",









            "event_class",
            "jetPID",
            #"event_type",
            "event_invariant_mass","event_njet",  
            "jet_mass","jet_p","jet_e", "jet_phi", "jet_theta", "jet_pT",

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

                "VertexX", "VertexY", "VertexZ"
                ]
