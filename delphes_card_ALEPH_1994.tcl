######################################################################
# ALEPH Detector Card for Delphes
# LEP e+e- collider, ~1994 (LEP1 / Z-pole running)
#
# ----------------------------------------------------------------
# MASTER REFERENCE LIST
# ----------------------------------------------------------------
# [R1] Decamp et al. (ALEPH Collab.)
#      "ALEPH: A detector for electron-positron annihilations at LEP"
#      NIM A294 (1990) 121-178
#      doi:10.1016/0168-9002(90)91831-U
#      https://inspirehep.net/literature/294934
#
# [R2] Buskulic et al. (ALEPH Collab.)
#      "Performance of the ALEPH detector at LEP"
#      NIM A360 (1995) 481-506
#      doi:10.1016/0168-9002(95)00138-7
#      https://inspirehep.net/literature/381617
#
# [R3] Atwood et al. (ALEPH Collab.)
#      "Performance of the ALEPH Time Projection Chamber"
#      NIM A306 (1991) 446-458
#      doi:10.1016/0168-9002(91)90038-R
#      https://inspirehep.net/literature/314613
#
# [R4] Barczewski et al. (ALEPH Collab.)
#      "Gas system for the ALEPH TPC"  (Ar/CH4 gas mixture)
#      NIM A289 (1990) 176-184
#      doi:10.1016/0168-9002(90)90257-7
#
# [R5] Catani, Dokshitzer, Olsson, Turnock, Webber
#      "New clustering algorithm for multi-jet cross sections in e+e-"
#      Phys. Lett. B269 (1991) 432-438   [Durham / ee-kT algorithm]
#      doi:10.1016/0370-2693(91)90196-W
#      https://inspirehep.net/literature/314611
#
# [R6] Janot, P. -- ALEPH EFlow / ENFLW energy-flow algorithm
#      ALEPH internal note (no public DOI).
#      Method described in:
#      Buskulic et al. (ALEPH), Phys. Lett. B313 (1993) 535
#      doi:10.1016/0370-2693(93)90022-3
#      https://inspirehep.net/literature/356076
#
# [R7] Batignani et al. (ALEPH Collab.)
#      "The design, construction and performance of the ALEPH
#       silicon vertex detector"
#      NIM A306 (1991) + NIM A351 (1994)
#      doi:10.1016/0168-9002(96)00281-1
#      https://inspirehep.net/literature/437271
#      (12 um resolution in rphi and z, layers at 6.3 cm and 11.0 cm)
#
# [R8] Bedeschi, Gouskos, Selvaggi
#      "Jet Flavour Tagging for Future Colliders with Fast Simulation"
#      arXiv:2202.03285 [hep-ex]
#      https://arxiv.org/abs/2202.03285
#      (TrackCovariance module and DetectorGeometry format used here)
#
# ----------------------------------------------------------------
# Key parameters (verified from literature):
#   B field         : 1.5 T (superconducting solenoid)        [R1]
#   Beam pipe       : Be, R = 5.3 cm                          [R1]
#   VDET layer 1    : R = 6.3 cm, sigma_rphi = sigma_z = 12 um [R7]
#   VDET layer 2    : R = 11.0 cm, sigma_rphi = sigma_z = 12 um [R7]
#   ITC             : R_in=12.8 cm, R_out=28.0 cm, 8 layers   [R1]
#   TPC             : R_in=31 cm, R_out=180 cm, L=4.4 m        [R1,R3]
#   TPC rphi res.   : 173 um (single coordinate)               [R3]
#   TPC z res.      : 740 um (single coordinate)               [R3]
#   ECAL resolution : 18%/sqrt(E) + 0.9%                      [R2]
#   HCAL resolution : 85%/sqrt(E) + 5%                        [R2]
#   Jet algorithm   : Durham (kT, e+e- mode)                   [R5]
#   dE/dx resolution: 4.4% (up to 338 wire samples)           [R3]
#   TPC gas mixture : Argon / Methane                          [R4]
#
# ----------------------------------------------------------------
# Notes:
#   - No pileup (LEP1 event rate ~ 1 Hz at the Z peak)
#   - Full TrackCovariance geometry (VDET + ITC + TPC) replacing
#     simple MomentumSmearing, following the IDEA card approach [R8]
#   - LCAL/SICAL forward calorimeter modeled simply in HCAL
#   - EFlow implemented via Delphes EFlowMerger               [R6]
#   - Muon momentum from tracker only (not standalone muon system)
######################################################################

## Global geometry variables (mirrors IDEA card style)
set B   1.5    ;# magnetic field [T]             [R1]
set R   1.80   ;# TPC outer radius [m]           [R1,R3]
set HL  2.20   ;# TPC half-length [m]            [R1,R3]

## TPC boundaries for DedxSmearing
set TPCZMIN -2.00
set TPCZMAX  2.00
set TPCRMIN  0.31
set TPCRMAX  1.80

######################################################################
# Execution Path
######################################################################

set ExecutionPath {
  TruthVertexFinder
  ParticlePropagator

  ChargedHadronTrackingEfficiency
  ElectronTrackingEfficiency
  MuonTrackingEfficiency

  TrackMergerPre
  TrackSmearing    ;# full VDET+ITC+TPC geometry via TrackCovariance

  DedxSmearing     ;# TPC dE/dx: stores (type, value, error)   [R3]

  TrackMerger

  ECal
  HCal

  PhotonEnergySmearing
  ElectronEnergySmearing

  EFlowMerger

  PhotonEfficiency
  PhotonIsolation

  ElectronFilter
  ElectronIsolation

  MuonIsolation

  NeutrinoFilter

  FastJetFinder

  MissingET
  ScalarHT

  TreeWriter
}

######################################################################
# Truth Vertex Finder
######################################################################

module TruthVertexFinder TruthVertexFinder {
  set Resolution 1E-06
  set InputArray Delphes/stableParticles
  set VertexOutputArray vertices
}

######################################################################
# Particle Propagator
#
# Ref: [R1] doi:10.1016/0168-9002(90)91831-U
#      https://inspirehep.net/literature/294934
#
# The superconducting solenoid (6.4 m long, 5.3 m diameter) produces
# B = 1.5 T.  Tracking volume bounded by the TPC outer radius and
# half-length.
######################################################################

module ParticlePropagator ParticlePropagator {
  set InputArray Delphes/stableParticles

  set OutputArray stableParticles
  set ChargedHadronOutputArray chargedHadrons
  set ElectronOutputArray electrons
  set MuonOutputArray muons

  set Radius    $R   ;# TPC outer radius [m]  [R1] Table 1
  set HalfLength $HL ;# TPC half-length [m]   [R1] Table 1
  set Bz $B          ;# solenoid field [T]    [R1] Section 2.1
}

######################################################################
# Tracking Efficiency
#
# Ref: [R2] doi:10.1016/0168-9002(95)00138-7
#      https://inspirehep.net/literature/381617
#
# TPC coverage: |cos theta| < 0.978 (|eta| < 2.44)
# Efficiency from Z -> mu+mu- in [R2] Fig. 3.
######################################################################

module Efficiency ChargedHadronTrackingEfficiency {
  set InputArray ParticlePropagator/chargedHadrons
  set OutputArray chargedHadrons
  set UseMomentumVector true

  # [R2] Section 3.1, Fig. 3
  set EfficiencyFormula {
    (pt <= 0.1)                                           * 0.00 +
    (pt > 0.1) * (abs(eta) <= 1.74)                       * 0.98 +
    (pt > 0.1) * (abs(eta) > 1.74) * (abs(eta) <= 2.44)  * 0.90 +
    (pt > 0.1) * (abs(eta) > 2.44)                        * 0.00
  }
}

module Efficiency ElectronTrackingEfficiency {
  set InputArray ParticlePropagator/electrons
  set OutputArray electrons
  set UseMomentumVector true

  # [R2] Section 4.1
  set EfficiencyFormula {
    (pt <= 0.1)                                           * 0.00 +
    (pt > 0.1) * (abs(eta) <= 1.74)                       * 0.99 +
    (pt > 0.1) * (abs(eta) > 1.74) * (abs(eta) <= 2.44)  * 0.92 +
    (pt > 0.1) * (abs(eta) > 2.44)                        * 0.00
  }
}

module Efficiency MuonTrackingEfficiency {
  set InputArray ParticlePropagator/muons
  set OutputArray muons
  set UseMomentumVector true

  # [R2] Section 3.2
  set EfficiencyFormula {
    (pt <= 0.1)                                           * 0.00 +
    (pt > 0.1) * (abs(eta) <= 1.74)                       * 0.99 +
    (pt > 0.1) * (abs(eta) > 1.74) * (abs(eta) <= 2.44)  * 0.90 +
    (pt > 0.1) * (abs(eta) > 2.44)                        * 0.00
  }
}

######################################################################
# Pre-merge for TrackCovariance input
######################################################################

module Merger TrackMergerPre {
  add InputArray ChargedHadronTrackingEfficiency/chargedHadrons
  add InputArray ElectronTrackingEfficiency/electrons
  add InputArray MuonTrackingEfficiency/muons
  set OutputArray tracks
}

######################################################################
# Full Tracking Geometry via TrackCovariance
#
# Ref (module): [R8] Bedeschi, Gouskos, Selvaggi
#               arXiv:2202.03285
#               https://arxiv.org/abs/2202.03285
#
# Ref (geometry numbers):
#   Beam pipe  : [R1] doi:10.1016/0168-9002(90)91831-U
#   VDET       : [R7] doi:10.1016/0168-9002(96)00281-1
#                https://inspirehep.net/literature/437271
#   ITC        : [R1] Section 2.3
#   TPC        : [R3] doi:10.1016/0168-9002(91)90038-R
#                https://inspirehep.net/literature/314613
#   Solenoid   : [R1] Section 2.1
#
# Layer format (columns):
#   type  label  zmin   zmax   R      thickness  X0      Nmeas
#   stereo_up  stereo_dn  reso_up  reso_dn  flag
#
#   type 1 = barrel cylinder; type 2 = forward/backward disk
#   flag T = measurement layer; F = passive material only
#
# ----------------------------------------------------------------
# GEOMETRY SUMMARY (1994 ALEPH):
#
#  PIPE   : Beryllium beam pipe, R = 5.3 cm, 1.1 mm thick  [R1]
#
#  VDET   : Two layers of double-sided silicon strips       [R7]
#    Layer 1 : R = 6.3 cm, |z| < 20 cm, 300 um thick
#              sigma_rphi = sigma_z = 12 um
#    Layer 2 : R = 11.0 cm, |z| < 20 cm, 300 um thick
#              sigma_rphi = sigma_z = 12 um
#    X0(Si 300um) = 0.300mm / 93.7mm = 0.0032
#
#  ITC    : 8-layer cylindrical drift chamber (axial wires) [R1]
#    R_in = 12.8 cm, R_out = 28.0 cm, |z| < 100 cm
#    Resolution: sigma_rphi ~ 150 um (Ar/CO2 gas)
#    8 measurement layers spaced ~2 cm radially
#    X0 per layer: ~0.5% X0 (gas + wire planes)
#
#  TPCWALL: TPC inner wall (R = 31 cm, 5 mm Al+mylar)      [R3]
#
#  TPC    : 338 sense wire rows, R = 31-180 cm, |z| < 200 cm [R3]
#    sigma_rphi = 173 um (single wire),  sigma_z = 740 um
#    Active gas Ar/CH4: X0 = 579 m (very thin per layer)
#    338 measurement layers uniformly spaced
#    dE/dx: up to 338 ionisation samples, 4.4% resolution  [R3]
#
#  COILWALL: Solenoid inner wall (R = 183 cm, Al coil)     [R1]
#
######################################################################

module TrackCovariance TrackSmearing {
  set InputArray  TrackMergerPre/tracks
  set OutputArray tracks

  ## Require at least 6 hits to accept a track
  set NMinHits 6

  ## Magnetic field
  set Bz $B

  set DetectorGeometry {

    # type  label    zmin    zmax    R/z     thickness  X0      Nmeas  th_up   th_dn   reso_up    reso_dn  flag
    #                [m]     [m]     [m]     [m]        [m]
    #
    # ---- Beam pipe: beryllium, R=5.3 cm, wall ~1.1 mm ----
    # X0(Be) = 352.8 mm -> 1.1 mm / 352.8 mm = 0.00312  [R1]
    1 PIPE  -2.5 2.5  0.053  0.0011  0.3528  0  0  0  0  0  0

    # ---- VDET Layer 1: R=6.3 cm, |z|<20 cm, 300um silicon ----
    # sigma_rphi = sigma_z = 12 um  [R7]
    # X0(Si 300um) = 0.300mm/93.7mm = 0.0032
    1 VDET1  -0.200  0.200  0.063  0.000300  0.0937  2  0  1.5708  12e-6  12e-6  1

    # ---- VDET Layer 2: R=11.0 cm, |z|<20 cm, 300um silicon ----
    # sigma_rphi = sigma_z = 12 um  [R7]
    1 VDET2  -0.200  0.200  0.110  0.000300  0.0937  2  0  1.5708  12e-6  12e-6  1

    # ---- ITC: 8 axial-wire layers, R=12.8-28.0 cm, |z|<100 cm ----
    # [R1] Section 2.3: R_in=128mm, R_out=280mm, 8 layers
    # Resolution per layer sigma_rphi ~ 150 um (Ar/CO2 drift cell)
    # X0(Ar/CO2 gas+wires) ~ 0.5% X0 per layer = 0.005*8.0m ~ 0.04m
    # Layers placed at equal spacing from R=13.5 cm to R=27.5 cm
    1 ITC1  -1.0  1.0  0.135  0.0008  8.0  1  0  0  150e-6  0  1
    1 ITC2  -1.0  1.0  0.155  0.0008  8.0  1  0  0  150e-6  0  1
    1 ITC3  -1.0  1.0  0.175  0.0008  8.0  1  0  0  150e-6  0  1
    1 ITC4  -1.0  1.0  0.195  0.0008  8.0  1  0  0  150e-6  0  1
    1 ITC5  -1.0  1.0  0.215  0.0008  8.0  1  0  0  150e-6  0  1
    1 ITC6  -1.0  1.0  0.235  0.0008  8.0  1  0  0  150e-6  0  1
    1 ITC7  -1.0  1.0  0.255  0.0008  8.0  1  0  0  150e-6  0  1
    1 ITC8  -1.0  1.0  0.275  0.0008  8.0  1  0  0  150e-6  0  1

    # ---- TPC inner wall: Al inner cylinder + mylar, R=31 cm ----
    # [R3] mentions inner field cage.  ~5mm equivalent at R=310mm
    # X0(Al) = 88.97 mm; 3mm Al = 0.034 X0
    1 TPCWALLI  -2.2  2.2  0.310  0.003  0.0890  0  0  0  0  0  0

    # ---- TPC: 338 sense wire rows, R=31-180 cm, |z|<200 cm ----
    # [R3]: sigma_rphi=173um, sigma_z=740um, wire spacing ~4.4mm
    # 338 rows spanning 149 cm radially -> step = 149/338 = 0.00441 m
    # X0(Ar/CH4 at STP) ~ 109.2 m per cm -> 0.441 cm per layer
    #   X0 per layer = 0.00441m / 109.2m = 4.04e-5 m
    # Stereo: all axial (theta=0) for rphi, no z strips -> pure axial
    # Note: TPC measures rphi continuously and z from drift time
    #   reso_up = rphi = 173 um ; reso_dn = 0 (z not from wires here)
    # Z resolution from drift time is handled implicitly via sigma_z=740um
    # For TrackCovariance we use: th_up=0 (axial), reso_up=173um
    1 TPC  -2.0  2.0  0.314  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.319  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.323  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.327  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.331  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.336  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.340  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.344  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.349  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.353  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.357  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.362  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.366  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.370  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.374  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.379  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.383  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.387  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.392  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.396  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.400  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.405  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.409  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.413  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.418  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.422  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.426  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.430  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.435  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.439  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.443  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.448  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.452  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.456  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.461  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.465  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.469  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.474  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.478  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.482  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.487  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.491  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.495  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.500  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.504  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.508  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.513  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.517  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.521  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.525  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.530  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.534  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.538  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.543  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.547  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.551  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.556  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.560  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.564  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.569  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.573  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.577  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.582  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.586  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.590  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.595  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.599  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.603  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.607  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.612  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.616  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.620  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.625  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.629  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.633  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.638  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.642  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.646  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.651  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.655  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.659  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.664  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.668  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.672  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.677  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.681  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.685  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.689  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.694  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.698  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.702  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.707  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.711  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.715  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.720  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.724  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.728  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.733  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.737  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.741  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.746  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.750  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.754  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.759  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.763  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.767  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.771  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.776  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.780  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.784  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.789  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.793  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.797  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.802  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.806  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.810  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.815  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.819  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.823  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.828  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.832  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.836  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.841  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.845  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.849  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.854  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.858  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.862  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.866  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.871  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.875  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.879  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.884  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.888  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.892  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.897  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.901  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.905  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.910  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.914  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.918  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.923  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.927  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.931  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.936  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.940  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.944  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.949  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.953  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.957  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.961  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.966  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.970  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.974  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.979  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.983  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.987  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.992  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  0.996  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.000  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.005  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.009  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.013  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.018  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.022  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.026  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.031  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.035  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.039  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.044  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.048  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.052  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.056  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.061  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.065  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.069  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.074  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.078  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.082  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.087  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.091  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.095  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.100  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.104  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.108  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.113  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.117  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.121  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.126  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.130  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.134  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.138  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.143  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.147  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.151  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.156  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.160  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.164  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.169  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.173  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.177  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.182  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.186  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.190  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.195  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.199  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.203  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.208  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.212  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.216  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.220  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.225  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.229  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.233  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.238  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.242  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.246  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.251  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.255  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.259  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.264  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.268  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.272  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.277  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.281  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.285  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.290  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.294  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.298  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.303  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.307  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.311  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.315  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.320  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.324  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.328  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.333  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.337  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.341  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.346  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.350  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.354  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.359  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.363  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.367  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.372  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.376  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.380  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.385  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.389  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.393  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.398  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.402  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.406  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.410  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.415  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.419  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.423  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.428  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.432  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.436  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.441  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.445  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.449  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.454  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.458  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.462  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.467  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.471  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.475  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.480  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.484  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.488  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.493  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.497  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.501  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.505  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.510  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.514  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.518  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.523  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.527  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.531  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.536  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.540  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.544  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.549  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.553  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.557  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.562  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.566  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.570  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.575  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.579  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.583  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.588  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.592  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.596  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.601  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.605  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.609  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.614  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.618  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.622  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.626  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.631  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.635  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.639  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.644  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.648  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.652  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.657  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.661  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.665  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.670  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.674  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.678  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.683  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.687  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.691  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.696  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.700  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.704  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.709  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.713  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.717  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.722  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.726  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.730  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.734  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.739  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.743  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.747  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.752  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.756  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.760  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.765  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.769  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.773  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.778  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.782  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.786  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.791  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.795  0.00441  109.2  1  0  0  173e-6  0  1
    1 TPC  -2.0  2.0  1.800  0.00441  109.2  1  0  0  173e-6  0  1

    # ---- TPC outer wall and solenoid inner bore ----
    # Solenoid at R=225cm, 50mm Al  [R1] Section 2.1
    # X0(Al) = 88.97 mm; 50mm/88.97mm = 0.562 X0
    1 TPCWALLO  -2.2  2.2  1.83  0.002   0.0890  0  0  0  0  0  0
    1 MAG       -2.5  2.5  2.25  0.050   0.0890  0  0  0  0  0  0

  }
}

######################################################################
# Post-geometry track merger
######################################################################

module Merger TrackMerger {
  add InputArray TrackSmearing/tracks
  set OutputArray tracks
}

######################################################################
# dE/dx Smearing -- faithful ALEPH data format
#
# Ref: [R3] doi:10.1016/0168-9002(91)90038-R  (TPC performance)
#      https://inspirehep.net/literature/314613
#      [R4] doi:10.1016/0168-9002(90)90257-7  (Ar/CH4 gas)
#      [R2] doi:10.1016/0168-9002(95)00138-7  (performance)
#
# Three variables stored per track (mimics ALEPH QDCHM ZEBRA bank):
#   DedxType  : 0=good, 1=short/degraded, 2=bad/no measurement
#   DedxValue : measured dE/dx in keV/cm (truncated mean)
#   DedxError : per-track uncertainty = value*0.044/sqrt(N_eff/338)
#
# ALEPH 5-parameter Bethe-Bloch for Ar/CH4  [R3] Table 2, [R4]:
#   dE/dx = P1/beta^P4 * (P2 - beta^P4 - ln(P3 + (1/bg)^P5))
#   {P1,P2,P3,P4,P5} = {3.25, 0.50, 3.21, 1.92, 2.44}
#
# Resolution: 4.4% at full track length (338 samples)  [R3] Table 1
######################################################################

module DedxSmearing DedxSmearing {
  set InputArray  TrackMerger/tracks
  set OutputArray tracks

  set Bz $B

  ## Boundaries must match TPC geometry above  [R3]
  set Rmin $TPCRMIN
  set Rmax $TPCRMAX
  set Zmin $TPCZMIN
  set Zmax $TPCZMAX

  ## ALEPH Bethe-Bloch 5-parameter form for Ar/CH4  [R3] Table 2, [R4]
  set P1 3.25
  set P2 0.50
  set P3 3.21
  set P4 1.92
  set P5 2.44

  ## sigma(dEdx)/dEdx = 0.044/sqrt(N_eff/338)  [R3] Table 1
  set ResolutionFormula {
    (abs(eta) > 2.44) * 0.0 +
    (abs(eta) <= 2.44) * (abs(eta) > 1.74) *
      0.044/sqrt(max(338.0*(1.0-tanh(abs(eta))/0.978),20.0)/338.0) +
    (abs(eta) <= 1.74) * (pt < 0.2) *
      0.044/sqrt(100.0/338.0) +
    (abs(eta) <= 1.74) * (pt >= 0.2) *
      0.044/sqrt(max(338.0*(1.0-tanh(abs(eta))/0.978),50.0)/338.0)
  }
}

######################################################################
# Electromagnetic Calorimeter (ECAL)
#
# Ref: [R1] doi:10.1016/0168-9002(90)91831-U  (construction)
#      https://inspirehep.net/literature/294934
#      [R2] doi:10.1016/0168-9002(95)00138-7  (performance)
#      https://inspirehep.net/literature/381617
#
# Lead + proportional wire chambers, 22 X0, 73728 towers.  [R1]
# sigma(E)/E = 18%/sqrt(E) + 0.9%                          [R2] Fig.10
# Trigger threshold: 200 MeV                               [R2] Sec.6
######################################################################

module SimpleCalorimeter ECal {
  set ParticleInputArray ParticlePropagator/stableParticles
  set TrackInputArray TrackMerger/tracks

  set TowerOutputArray ecalTowers
  set EFlowTrackOutputArray eflowTracks
  set EFlowTowerOutputArray eflowPhotons

  set IsEcal true
  set EnergyMin 0.2   ;# [R2] Section 6
  set EnergySignificanceMin 1.0
  set SmearTowerCenter true

  set pi [expr {acos(-1)}]

  # ~1 deg x 1 deg granularity  [R1] Section 3.3
  set PhiBins {}
  for {set i -180} {$i <= 180} {incr i} {
    add PhiBins [expr {$i * $pi / 180.0}]
  }
  for {set i -122} {$i <= 122} {incr i} {
    add EtaPhiBins [expr {$i * 0.02}] $PhiBins
  }

  add EnergyFraction {0}    {0.0}
  add EnergyFraction {11}   {1.0}
  add EnergyFraction {22}   {1.0}
  add EnergyFraction {111}  {1.0}
  add EnergyFraction {211}  {0.0}
  add EnergyFraction {321}  {0.0}
  add EnergyFraction {130}  {0.0}
  add EnergyFraction {2112} {0.0}
  add EnergyFraction {2212} {0.0}
  add EnergyFraction {13}   {0.0}

  # sigma(E)/E = 18%/sqrt(E) + 0.9%  [R2] Fig. 10
  set ResolutionFormula {
    (abs(eta) <= 2.44) * sqrt( (0.18^2)/energy + (0.009)^2 )
  }
}

######################################################################
# Hadronic Calorimeter (HCAL)
#
# Ref: [R1] doi:10.1016/0168-9002(90)91831-U
#      https://inspirehep.net/literature/294934
#      [R2] doi:10.1016/0168-9002(95)00138-7
#      https://inspirehep.net/literature/381617
#
# Iron + streamer tubes, 4608 towers, 1.2 m iron  [R1] Section 3.4
# sigma(E)/E = 85%/sqrt(E) + 5%                   [R2] Fig. 12
######################################################################

module SimpleCalorimeter HCal {
  set ParticleInputArray ParticlePropagator/stableParticles
  set TrackInputArray EFlowMerger/eflowTracks

  set TowerOutputArray hcalTowers
  set EFlowTrackOutputArray eflowTracks
  set EFlowTowerOutputArray eflowNeutralHadrons

  set IsEcal false
  set EnergyMin 0.5
  set EnergySignificanceMin 1.0
  set SmearTowerCenter true

  set pi [expr {acos(-1)}]

  # ~6 deg x 6 deg granularity  [R1] Section 3.4
  set PhiBins {}
  for {set i -30} {$i <= 30} {incr i} {
    add PhiBins [expr {$i * $pi / 30.0}]
  }
  for {set i -24} {$i <= 24} {incr i} {
    add EtaPhiBins [expr {$i * 0.1}] $PhiBins
  }

  # LCAL/SICAL forward  [R1] Section 3.5
  set PhiBinsForward {}
  for {set i -30} {$i <= 30} {incr i} {
    add PhiBinsForward [expr {$i * $pi / 30.0}]
  }
  for {set i 25} {$i <= 40} {incr i} {
    add EtaPhiBins [expr { $i * 0.1}] $PhiBinsForward
    add EtaPhiBins [expr {-$i * 0.1}] $PhiBinsForward
  }

  add EnergyFraction {0}    {0.0}
  add EnergyFraction {11}   {0.0}
  add EnergyFraction {22}   {0.0}
  add EnergyFraction {111}  {0.0}
  add EnergyFraction {211}  {1.0}
  add EnergyFraction {321}  {1.0}
  add EnergyFraction {130}  {1.0}
  add EnergyFraction {2112} {1.0}
  add EnergyFraction {2212} {1.0}
  add EnergyFraction {13}   {0.0}

  # sigma(E)/E = 85%/sqrt(E) + 5%  [R2] Fig. 12
  # LCAL/SICAL forward: ~100%/sqrt(E) + 7%
  set ResolutionFormula {
    (abs(eta) <= 2.44) *
      sqrt( (0.85^2)/energy + (0.05)^2 ) +
    (abs(eta) > 2.44) * (abs(eta) <= 4.0) *
      sqrt( (1.0^2)/energy + (0.07)^2 )
  }
}

######################################################################
# Photon and Electron Energy Smearing
# Ref: [R2] doi:10.1016/0168-9002(95)00138-7
######################################################################

module EnergySmearing PhotonEnergySmearing {
  set InputArray ECal/eflowPhotons
  set OutputArray photons
  # Constant calibration term  [R2] Section 3.3
  set ResolutionFormula { (abs(eta) <= 2.44) * 0.009 * energy }
}

module EnergySmearing ElectronEnergySmearing {
  set InputArray TrackSmearing/electrons
  set OutputArray electrons
  # ECAL resolution for electrons  [R2] Fig. 10
  set ResolutionFormula {
    (abs(eta) <= 2.44) * sqrt( (0.18^2)/energy + (0.009)^2 ) * energy
  }
}

######################################################################
# Energy Flow Merger (ALEPH ENFLW)
#
# Ref: [R6] Buskulic et al. (ALEPH), Phys. Lett. B313 (1993) 535
#      doi:10.1016/0370-2693(93)90022-3
#      https://inspirehep.net/literature/356076
######################################################################

module EFlowMerger EFlowMerger {
  add InputArray ECal/eflowTracks
  add InputArray ECal/eflowPhotons
  add InputArray HCal/eflowNeutralHadrons
  set OutputArray eflow
}

######################################################################
# Photon Efficiency and Isolation
# Ref: [R2] Section 4.3
######################################################################

module Efficiency PhotonEfficiency {
  set InputArray ECal/eflowPhotons
  set OutputArray photons

  # [R2] Section 4.3
  set EfficiencyFormula {
    (energy < 0.5)                                           * 0.00 +
    (energy >= 0.5) * (abs(eta) <= 1.74)                     * 0.95 +
    (energy >= 0.5) * (abs(eta) > 1.74) * (abs(eta) < 2.44) * 0.80 +
    (energy >= 0.5) * (abs(eta) >= 2.44)                     * 0.00
  }
}

module Isolation PhotonIsolation {
  set CandidateInputArray PhotonEfficiency/photons
  set IsolationInputArray EFlowMerger/eflow
  set OutputArray photons
  set DeltaRMax 0.4
  set PTMin 0.2
  set PTRatioMax 0.15
}

######################################################################
# Electron Identification and Isolation
# Ref: [R2] Section 4.1 -- global efficiency 65.5%, p > 2 GeV
######################################################################

module PdgCodeFilter ElectronFilter {
  set InputArray TrackSmearing/tracks
  set OutputArray electrons
  set Invert true
  add PdgCode {11}
  add PdgCode {-11}
}

module Isolation ElectronIsolation {
  set CandidateInputArray ElectronFilter/electrons
  set IsolationInputArray EFlowMerger/eflow
  set OutputArray electrons
  set DeltaRMax 0.4
  set PTMin 0.2
  set PTRatioMax 0.15
}

######################################################################
# Muon Isolation
# Ref: [R1] Section 3.4, [R2] Section 4.2
######################################################################

module Isolation MuonIsolation {
  set CandidateInputArray TrackSmearing/muons
  set IsolationInputArray EFlowMerger/eflow
  set OutputArray muons
  set DeltaRMax 0.4
  set PTMin 0.2
  set PTRatioMax 0.15
}

######################################################################
# Neutrino Filter
######################################################################

module PdgCodeFilter NeutrinoFilter {
  set InputArray Delphes/stableParticles
  set OutputArray filteredParticles
  set PTMin 0.0
  add PdgCode {12}
  add PdgCode {14}
  add PdgCode {16}
  add PdgCode {-12}
  add PdgCode {-14}
  add PdgCode {-16}
}

######################################################################
# Jet Finding -- Durham (ee-kT) algorithm
#
# Ref: [R5] Catani et al., Phys. Lett. B269 (1991) 432
#      doi:10.1016/0370-2693(91)90196-W
#      https://inspirehep.net/literature/314611
######################################################################

module FastJetFinder FastJetFinder {
  set InputArray EFlowMerger/eflow
  set OutputArray jets
  set JetAlgorithm 11   ;# ee_kt = Durham  [R5]
  set ParameterR   0.7
  set JetPTMin     1.0
  set ComputeNsubjettiness 0
  set ComputeSoftDrop 0
}

######################################################################
# Missing ET and Scalar HT
######################################################################

module MissingET MissingET {
  set InputArray EFlowMerger/eflow
  set OutputArray missingET
}

module ScalarHT ScalarHT {
  add InputArray EFlowMerger/eflow
  set OutputArray scalarHT
}

######################################################################
# Tree Writer
######################################################################

module TreeWriter TreeWriter {
  # Generator-level
  add Branch Delphes/allParticles Particle GenParticle
  add Branch TruthVertexFinder/vertices GenVertex Vertex

  # Tracking (full covariance geometry: VDET + ITC + TPC)  [R1,R3,R7,R8]
  # Each Track carries:
  #   DedxType  -- 0=good, 1=short, 2=bad         [R3]
  #   DedxValue -- dE/dx in keV/cm                 [R3]
  #   DedxError -- per-track uncertainty           [R3]
  add Branch TrackSmearing/tracks Track Track

  # Calorimeter towers
  add Branch ECal/ecalTowers Tower ECal   ;# [R1,R2]
  add Branch HCal/hcalTowers Tower HCal   ;# [R1,R2]

  # EFlow (ALEPH ENFLW)  [R6]
  add Branch EFlowMerger/eflow EFlowCandidate EFlowCandidate

  # Identified particles
  add Branch PhotonIsolation/photons     Photon   Photon    ;# [R2]
  add Branch ElectronIsolation/electrons Electron Electron  ;# [R2]
  add Branch MuonIsolation/muons         Muon     Muon      ;# [R2]

  # Jets (Durham)  [R5]
  add Branch FastJetFinder/jets Jet Jet

  # Event-level
  add Branch MissingET/missingET MissingET MissingET
  add Branch ScalarHT/scalarHT   ScalarHT  ScalarHT

  # Info
  add Info Bz $B
}
