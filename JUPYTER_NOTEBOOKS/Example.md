Loading the gRASPA library
==========================


```python
import gRASPA as g
```

# Running gRASPA, the original taste with no modifications
==========================
You can set the cycles here, it will overwrite the values read from `simulation.input` file


```python
W = g.Initialize()
W.NumberOfInitializationCycles = 100
W.NumberOfProductionCycles = 0
g.RUN(W)
g.finalize(W)
```

    /home/z/anaconda3/bin/python3.12------------------GENERAL SIMULATION SETUP-------------
    Finished Checking Number of Components, There are 1 framework, 1 Adsorbates, 2 total Components
    DONE Reading Model Info from simulation.input file
    Running Cycles in the Normal Way
    -------------------------------------------------------
    device_random[0] = 2.30000 4.50000 6.70000
    ==========================================
    ====== Preparing Simulation box 0 ======
    ==========================================
    Parsing [1] Component
    -------------- READING AdsorbateComponent 0 (CO2) --------------
    ==================================================
    ACCUMULATED Probabilities:
    Translation Probability:      0.25000
    Rotation Probability:         0.50000
    Special Rotation Probability: 0.50000
    Widom Probability:            0.50000
    Reinsertion Probability:      0.75000
    Identity Swap Probability:    0.75000
    CBCF Swap Probability:        0.75000
    Swap Probability:             1.00000
    Volume Probability:           1.00000
    Gibbs Swap Probability:       1.00000
    Gibbs Volume Probability:     1.00000
    Sum of Probabilities:         1.00000
    ==================================================
    -------------- END OF READING Component 0 (CO2) --------------
    Rosen capacity: 10
    Allocated 6 double3 for reinsertion!
    HHVDW: 0.00000, HHReal: 0.00000, HGVDW: 0.00000, HGReal: 0.00000, GGVDW: 0.00000, GGReal: 0.00000, HHEwaldE: 0.00004,
     HGEwaldE: 0.00000,
     GGEwaldE: 0.00000, TailE: 0.00000, DNN_E: 0.00000
    Stored HGVDW: 0.00000, Stored HGReal: 0.00000, Stored HGEwaldE: 0.00000
    HHVDW: 0.00000, HHReal: 0.00000, HGVDW: 0.00000, HGReal: 0.00000, GGVDW: 0.00000, GGReal: 0.00000, HHEwaldE: 0.00004,
     HGEwaldE: 0.00000,
     GGEwaldE: 0.00000, TailE: 0.00000, DNN_E: 0.00000
    Stored HGVDW: 0.00000, Stored HGReal: 0.00000, Stored HGEwaldE: 0.00000
    ============================================
    == END OF PREPARATION, SIMULATION STARTS! ==
    ============================================
    ========================
    == END OF SIMULATION! ==
    ========================
    HHVDW: 0.00000, HHReal: 0.00000, HGVDW: -34420.94011, HGReal: -412.10544, GGVDW: -815.79859, GGReal: 142.47937, HHEwaldE: 0.00004,
     HGEwaldE: -8191.68084,
     GGEwaldE: 14.37537, TailE: 0.00000, DNN_E: 0.00000
    Stored HGVDW: 0.00000, Stored HGReal: 0.00000, Stored HGEwaldE: 0.00000
    
    ===========================
    ==    END OF PROGRAM!    ==
    == PRINTING MEMORY USAGE ==
    ===========================
    Total Program Size: 6561 MB
    Resident Set Size: 222 MB
    Shared Pages: 133 MB
    Text (code): 2 MB
    Data + Stack: 231 MB


    ================FUGACITY COEFFICIENT CALCULATION================
    Checking: Current Fugacity Coeff for 1 component: 1.00000
    Every Adsorbate Component has fugacity coefficient assigned, skip EOS calculation!
    ----------------- MEMORY ALLOCAION STATUS -----------------
    System allocate_sizes are: 2304, 2000
    Component allocate_sizes are: 2304, 2000
    Allocated Blocksum size: 3601, vdw_real size: 3601, fourier_size: 0
    ------------------------------------------------------------
    ------------------- SIMULATION BOX PARAMETERS -----------------
    Pressure:        0.00060
    Box Volume:      42656.20098
    Box Beta:        0.00404
    Box Temperature: 298.00000
    ---------------------------------------------------------------
    ======================== CALCULATING INITIAL STAGE ENERGY ========================
    ****** Calculating VDW + Real Energy (CPU) ******
    Host-Host   VDW: 0.00000; Real: 0.00000
    Host-Guest  VDW: 0.00000; Real: 0.00000
    Guest-Guest VDW: 0.00000; Real: 0.00000
    ********** PRINTING COMPONENT ENERGIES**********
    Compoent [0-0], VDW: 0.00000, Real: 0.00000
    Compoent [0-1], VDW: 0.00000, Real: 0.00000
    Compoent [1-1], VDW: 0.00000, Real: 0.00000
    ****** Calculating Ewald Energy (CPU) ******
    CPU Guest-Guest Fourier: 0.00000, Host-Host Fourier: 409293.46532, Framework-Guest Fourier: 0.00000
    Component: 0, SelfAtomE: 100586536.82400 (120977870.55210 kJ/mol)
    Component: 1, SelfAtomE: 0.00000 (0.00000 kJ/mol)
    Component: 0, Intra-Molecular ExclusionE: 115033008.45111 (138352992.80625 kJ/mol)
    Component: 1, Intra-Molecular ExclusionE: 0.00000 (0.00000 kJ/mol)
    HostEwald took 0.20396 sec
    Ewald Summation (total energy) on the CPU took 0.20396 secs
    Component 0, Intra Exclusion Energy: -115033008.45111 (-138352992.80625 kJ/mol)
    Component 0, Atom Self Exclusion Energy: 100586536.82400 (120977870.55210 kJ/mol)
    DEBUG: comp: 0, IntraE: -115033008.45111, SelfE: 100586536.82400
    Component 1, Intra Exclusion Energy: -13173.67573 (-15844.29972 kJ/mol)
    Component 1, Atom Self Exclusion Energy: 13215.98988 (15895.19197 kJ/mol)
    DEBUG: comp: 1, IntraE: -13173.67573, SelfE: 13215.98988
    ******   Allocating Ewald WaveVectors + StructureFactors (INITIAL STAGE ONLY)   ******
    Allocated 110592 110592 73728 space for eikxyz
    Structure Factor 0 is 0.00000 0.00000
    Structure Factor 1 is 0.00000 0.00000
    Structure Factor 2 is 0.00000 0.00000
    Structure Factor 3 is 0.00000 0.00000
    Structure Factor 4 is 0.00000 0.00000
    Structure Factor 5 is 0.00000 0.00000
    Structure Factor 6 is 0.00000 0.00000
    Structure Factor 7 is 0.00000 0.00000
    Structure Factor 8 is 0.00000 0.00000
    Structure Factor 9 is 0.00000 0.00000
    ****** DONE Allocating Ewald WaveVectors + StructureFactors(INITIAL STAGE ONLY) ******
     ****** CHECKING StructureFactors (SF) Stored on CPU vs. GPU ****** 
    CPU SF: 4140, GPU SF: 4140
    StructureFactor 0, CPU: 0.00000 0.00000, GPU: 0.00000 0.00000
    StructureFactor 1, CPU: 0.00000 0.00000, GPU: 0.00000 0.00000
    StructureFactor 2, CPU: 0.00000 0.00000, GPU: 0.00000 0.00000
    StructureFactor 3, CPU: 0.00000 0.00000, GPU: 0.00000 0.00000
    StructureFactor 4, CPU: 0.00000 0.00000, GPU: 0.00000 0.00000
    StructureFactor 5, CPU: 0.00000 0.00000, GPU: 0.00000 0.00000
    StructureFactor 6, CPU: 0.00000 0.00000, GPU: 0.00000 0.00000
    StructureFactor 7, CPU: 0.00000 0.00000, GPU: 0.00000 0.00000
    StructureFactor 8, CPU: 0.00000 0.00000, GPU: 0.00000 0.00000
    StructureFactor 9, CPU: 0.00000 0.00000, GPU: 0.00000 0.00000
     ****** CHECKING Framework StructureFactors Stored on CPU ****** 
    Framework Structure Factor 0, real: 0.00000 imag: 0.00000
    Framework Structure Factor 1, real: 0.00000 imag: 0.00000
    Framework Structure Factor 2, real: 0.00000 imag: 0.00000
    Framework Structure Factor 3, real: 0.00000 imag: 0.00000
    Framework Structure Factor 4, real: 0.00000 imag: 0.00000
    Framework Structure Factor 5, real: -0.00000 imag: 0.00000
    Framework Structure Factor 6, real: 0.00000 imag: 0.00000
    Framework Structure Factor 7, real: -0.00000 imag: -0.00000
    Framework Structure Factor 8, real: 0.00000 imag: 0.00000
    Framework Structure Factor 9, real: 0.00000 imag: -0.00000
    VDW + Real on the GPU took 0.00000 secs
    Ewald Summation (total energy) on the GPU took 0.00141 secs
    Total GPU Energy: 
    ====================== DONE CALCULATING INITIAL STAGE ENERGY ======================
    Component 1, Need to create 0 full molecule
    ======================== CALCULATING CREATE_MOLECULE STAGE ENERGY ========================
    ****** Calculating VDW + Real Energy (CPU) ******
    Host-Host   VDW: 0.00000; Real: 0.00000
    Host-Guest  VDW: 0.00000; Real: 0.00000
    Guest-Guest VDW: 0.00000; Real: 0.00000
    ********** PRINTING COMPONENT ENERGIES**********
    Compoent [0-0], VDW: 0.00000, Real: 0.00000
    Compoent [0-1], VDW: 0.00000, Real: 0.00000
    Compoent [1-1], VDW: 0.00000, Real: 0.00000
    ****** Calculating Ewald Energy (CPU) ******
    CPU Guest-Guest Fourier: 0.00000, Host-Host Fourier: 409293.46532, Framework-Guest Fourier: 0.00000
    Component: 0, SelfAtomE: 100586536.82400 (120977870.55210 kJ/mol)
    Component: 1, SelfAtomE: 0.00000 (0.00000 kJ/mol)
    Component: 0, Intra-Molecular ExclusionE: 115033008.45111 (138352992.80625 kJ/mol)
    Component: 1, Intra-Molecular ExclusionE: 0.00000 (0.00000 kJ/mol)
    HostEwald took 0.20020 sec
    Ewald Summation (total energy) on the CPU took 0.20020 secs
     ****** CHECKING StructureFactors (SF) Stored on CPU vs. GPU ****** 
    CPU SF: 4140, GPU SF: 4140
    StructureFactor 0, CPU: 0.00000 0.00000, GPU: 0.00000 0.00000
    StructureFactor 1, CPU: 0.00000 0.00000, GPU: 0.00000 0.00000
    StructureFactor 2, CPU: 0.00000 0.00000, GPU: 0.00000 0.00000
    StructureFactor 3, CPU: 0.00000 0.00000, GPU: 0.00000 0.00000
    StructureFactor 4, CPU: 0.00000 0.00000, GPU: 0.00000 0.00000
    StructureFactor 5, CPU: 0.00000 0.00000, GPU: 0.00000 0.00000
    StructureFactor 6, CPU: 0.00000 0.00000, GPU: 0.00000 0.00000
    StructureFactor 7, CPU: 0.00000 0.00000, GPU: 0.00000 0.00000
    StructureFactor 8, CPU: 0.00000 0.00000, GPU: 0.00000 0.00000
    StructureFactor 9, CPU: 0.00000 0.00000, GPU: 0.00000 0.00000
     ****** CHECKING Framework StructureFactors Stored on CPU ****** 
    Framework Structure Factor 0, real: 0.00000 imag: 0.00000
    Framework Structure Factor 1, real: 0.00000 imag: 0.00000
    Framework Structure Factor 2, real: 0.00000 imag: 0.00000
    Framework Structure Factor 3, real: 0.00000 imag: 0.00000
    Framework Structure Factor 4, real: 0.00000 imag: 0.00000
    Framework Structure Factor 5, real: -0.00000 imag: 0.00000
    Framework Structure Factor 6, real: 0.00000 imag: 0.00000
    Framework Structure Factor 7, real: -0.00000 imag: -0.00000
    Framework Structure Factor 8, real: 0.00000 imag: 0.00000
    Framework Structure Factor 9, real: 0.00000 imag: -0.00000
    VDW + Real on the GPU took 0.00000 secs
    Ewald Summation (total energy) on the GPU took 0.00131 secs
    Total GPU Energy: 
    ====================== DONE CALCULATING CREATE_MOLECULE STAGE ENERGY ======================
    Running Simulation Boxes in SERIAL, currently [0] box; pres: 10000.00000 [Pa], temp: 298.00000 [K]
    ==================================
    == RUNNING INITIALIZATION PHASE ==
    ==================================
    CBMC Uses 10 trial positions and 10 trial orientations
    Box 0, Volume: 42656.20098
    Total Volume: 42656.20098
    INITIALIZATION Cycle: 0, 0 Adsorbate Molecules, Total Energy: 0.00000  ||  Component 0 [MFI-2x2x2-P1.cif], 1 Molecules  ||  Component 1 [CO2], 0 Molecules  ||  
    ======================== MOVE STATISTICS FOR COMPONENT [1] (CO2) ========================
    =====================TRANSLATION MOVES=====================
    Translation Performed: 503
    Translation Accepted: 11
    Max Translation: 4.0044000000, 3.9798000000, 2.6766000000
    ===========================================================
    =====================SWAP MOVES=====================
    Insertion Performed:   268
    Insertion Accepted:    107
    Deletion Performed:    274
    Deletion Accepted:     89
    Reinsertion Performed: 484
    Reinsertion Accepted:  109
    ====================================================
    =====================IDENTITY SWAP MOVES=====================
    =============================================================
    ================================================================================================
    ===============================
    == INITIALIZATION PHASE ENDS ==
    ===============================
    Running Simulation Boxes in SERIAL, currently [0] box; pres: 10000.00000 [Pa], temp: 298.00000 [K]
    ==================================
    == RUNNING EQUILIBRATION PHASE ==
    ==================================
    CBMC Uses 10 trial positions and 10 trial orientations
    ===============================
    == EQUILIBRATION PHASE ENDS ==
    ===============================
    Running Simulation Boxes in SERIAL, currently [0] box; pres: 10000.00000 [Pa], temp: 298.00000 [K]
    ==================================
    ==  RUNNING PRODUCTION PHASE   ==
    ==================================
    CBMC Uses 10 trial positions and 10 trial orientations
    ===============================
    == PRODUCTION PHASE ENDS ==
    ===============================
    Work took 0.156148 seconds
    ======================================
    CHECKING FINAL ENERGY FOR SYSTEM [0]
    ======================================
    ======================== CALCULATING FINAL STAGE ENERGY ========================
    ****** Calculating VDW + Real Energy (CPU) ******
    Host-Host   VDW: 0.00000; Real: 0.00000
    Host-Guest  VDW: -34420.94011; Real: -412.10544
    Guest-Guest VDW: -815.79859; Real: 142.47937
    ********** PRINTING COMPONENT ENERGIES**********
    Compoent [0-0], VDW: 0.00000, Real: 0.00000
    Compoent [0-1], VDW: -34420.94011, Real: -412.10544
    Compoent [1-1], VDW: -815.79859, Real: 142.47937
    ****** Calculating Ewald Energy (CPU) ******
    CPU Guest-Guest Fourier: 776.03005, Host-Host Fourier: 409293.46532, Framework-Guest Fourier: -8191.68084
    Component: 0, SelfAtomE: 100586536.82400 (120977870.55210 kJ/mol)
    Component: 1, SelfAtomE: 237887.81780 (286113.45550 kJ/mol)
    Component: 0, Intra-Molecular ExclusionE: 115033008.45111 (138352992.80625 kJ/mol)
    Component: 1, Intra-Molecular ExclusionE: 237126.16311 (285197.39492 kJ/mol)
    HostEwald took 0.20168 sec
    Ewald Summation (total energy) on the CPU took 0.20168 secs
     ****** CHECKING StructureFactors (SF) Stored on CPU vs. GPU ****** 
    CPU SF: 4140, GPU SF: 4140
    StructureFactor 0, CPU: 0.00000 0.00000, GPU: 0.00000 0.00000
    StructureFactor 1, CPU: 0.00000 0.00000, GPU: 0.00000 0.00000
    StructureFactor 2, CPU: 0.00000 0.00000, GPU: 0.00000 0.00000
    StructureFactor 3, CPU: 0.00000 0.00000, GPU: 0.00000 0.00000
    StructureFactor 4, CPU: 1.12807 2.62738, GPU: 1.12807 2.62738
    StructureFactor 5, CPU: -2.13399 -0.45637, GPU: -2.13399 -0.45637
    StructureFactor 6, CPU: 0.38720 0.74671, GPU: 0.38720 0.74671
    StructureFactor 7, CPU: -1.71523 -2.01338, GPU: -1.71523 -2.01338
    StructureFactor 8, CPU: -0.04741 2.85855, GPU: -0.04741 2.85855
    StructureFactor 9, CPU: -1.00131 -1.23950, GPU: -1.00131 -1.23950
     ****** CHECKING Framework StructureFactors Stored on CPU ****** 
    Framework Structure Factor 0, real: 0.00000 imag: 0.00000
    Framework Structure Factor 1, real: 0.00000 imag: 0.00000
    Framework Structure Factor 2, real: 0.00000 imag: 0.00000
    Framework Structure Factor 3, real: 0.00000 imag: 0.00000
    Framework Structure Factor 4, real: 0.00000 imag: 0.00000
    Framework Structure Factor 5, real: -0.00000 imag: 0.00000
    Framework Structure Factor 6, real: 0.00000 imag: 0.00000
    Framework Structure Factor 7, real: -0.00000 imag: -0.00000
    Framework Structure Factor 8, real: 0.00000 imag: 0.00000
    Framework Structure Factor 9, real: 0.00000 imag: -0.00000
    VDW + Real on the GPU took 0.00052 secs
    Ewald Summation (total energy) on the GPU took 0.00149 secs
    Total GPU Energy: 
    ====================== DONE CALCULATING FINAL STAGE ENERGY ======================
    ======================================
    Random Numbers Regenerated 0 times, offset: 28202, randomsize: 333334
    DNN Feature Preparation Time: 0.00000, DNN Prediction Time: 0.00000
    DNN GPU Time: 0.00000, DNN Sort Time: 0.00000, std::sort Time: 0.00000, Featurization Time: 0.00000
    ======================== ENERGY SUMMARY (Simulation 0) =========================
     *** INITIAL STAGE *** 
    ========================================================================
    VDW [Host-Host]:            0.00000 (0.00000 [K])
    VDW [Host-Guest]:           0.00000 (0.00000 [K])
    VDW [Guest-Guest]:          0.00000 (0.00000 [K])
    Real Coulomb [Host-Host]:   0.00000 (0.00000 [K])
    Real Coulomb [Host-Guest]:  0.00000 (0.00000 [K])
    Real Coulomb [Guest-Guest]: 0.00000 (0.00000 [K])
    Ewald [Host-Host]:          0.00000 (0.00000 [K])
     --> Total Ewald [Host-Host]:
          14855765.09240 (17867389.44443 [K])
     --> Initial Ewald [Host-Host] (excluded):
          14855765.09240 (17867389.44443 [K])
    Ewald [Host-Guest]:         0.00000 (0.00000 [K])
    Ewald [Guest-Guest]:        0.00000 (0.00000 [K])
    DNN Energy:                 0.00000 (0.00000 [K])
    Tail Correction Energy:     0.00000 (0.00000 [K])
    Total Energy:               0.00000 (0.00000 [K])
    ========================================================================
     *** CREATE MOLECULE STAGE *** 
    ========================================================================
    VDW [Host-Host]:            0.00000 (0.00000 [K])
    VDW [Host-Guest]:           0.00000 (0.00000 [K])
    VDW [Guest-Guest]:          0.00000 (0.00000 [K])
    Real Coulomb [Host-Host]:   0.00000 (0.00000 [K])
    Real Coulomb [Host-Guest]:  0.00000 (0.00000 [K])
    Real Coulomb [Guest-Guest]: 0.00000 (0.00000 [K])
    Ewald [Host-Host]:          0.00000 (0.00000 [K])
     --> Total Ewald [Host-Host]:
          14855765.09240 (17867389.44443 [K])
     --> Initial Ewald [Host-Host] (excluded):
          14855765.09240 (17867389.44443 [K])
    Ewald [Host-Guest]:         0.00000 (0.00000 [K])
    Ewald [Guest-Guest]:        0.00000 (0.00000 [K])
    DNN Energy:                 0.00000 (0.00000 [K])
    Tail Correction Energy:     0.00000 (0.00000 [K])
    Total Energy:               0.00000 (0.00000 [K])
    ========================================================================
     *** RUNNING DELTA_E (CREATE MOLECULE - INITIAL) *** 
    ========================================================================
    VDW [Host-Host]:            0.00000 (0.00000 [K])
    VDW [Host-Guest]:           0.00000 (0.00000 [K])
    VDW [Guest-Guest]:          0.00000 (0.00000 [K])
    Real Coulomb [Host-Host]:   0.00000 (0.00000 [K])
    Real Coulomb [Host-Guest]:  0.00000 (0.00000 [K])
    Real Coulomb [Guest-Guest]: 0.00000 (0.00000 [K])
    Ewald [Host-Host]:          0.00000 (0.00000 [K])
    Ewald [Host-Guest]:         0.00000 (0.00000 [K])
    Ewald [Guest-Guest]:        0.00000 (0.00000 [K])
    DNN Energy:                 0.00000 (0.00000 [K])
    Tail Correction Energy:     0.00000 (0.00000 [K])
    Total Energy:               0.00000 (0.00000 [K])
    ========================================================================
     *** CHECK DELTA_E (CREATE MOLECULE - INITIAL) *** 
    ========================================================================
    VDW [Host-Host]:            0.00000 (0.00000 [K])
    VDW [Host-Guest]:           0.00000 (0.00000 [K])
    VDW [Guest-Guest]:          0.00000 (0.00000 [K])
    Real Coulomb [Host-Host]:   0.00000 (0.00000 [K])
    Real Coulomb [Host-Guest]:  0.00000 (0.00000 [K])
    Real Coulomb [Guest-Guest]: 0.00000 (0.00000 [K])
    Ewald [Host-Host]:          0.00000 (0.00000 [K])
    Ewald [Host-Guest]:         0.00000 (0.00000 [K])
    Ewald [Guest-Guest]:        0.00000 (0.00000 [K])
    DNN Energy:                 0.00000 (0.00000 [K])
    Tail Correction Energy:     0.00000 (0.00000 [K])
    Total Energy:               0.00000 (0.00000 [K])
    ========================================================================
     *** FINAL STAGE *** 
    ========================================================================
    VDW [Host-Host]:            0.00000 (0.00000 [K])
    VDW [Host-Guest]:           -34420.94011 (-41398.90057 [K])
    VDW [Guest-Guest]:          -815.79859 (-981.18077 [K])
    Real Coulomb [Host-Host]:   0.00000 (0.00000 [K])
    Real Coulomb [Host-Guest]:  -412.10544 (-495.64922 [K])
    Real Coulomb [Guest-Guest]: 142.47937 (171.36339 [K])
    Ewald [Host-Host]:          0.00000 (0.00000 [K])
     --> Total Ewald [Host-Host]:
          14855765.09240 (17867389.44443 [K])
     --> Initial Ewald [Host-Host] (excluded):
          14855765.09240 (17867389.44443 [K])
    Ewald [Host-Guest]:         -8191.68084 (-9852.33348 [K])
    Ewald [Guest-Guest]:        14.37537 (17.28960 [K])
    DNN Energy:                 0.00000 (0.00000 [K])
    Tail Correction Energy:     0.00000 (0.00000 [K])
    Total Energy:               -43683.67024 (-52539.41104 [K])
    ========================================================================
     *** RUNNING DELTA_E (FINAL - CREATE MOLECULE) *** 
    ========================================================================
    VDW [Host-Host]:            0.00000 (0.00000 [K])
    VDW [Host-Guest]:           -34420.94011 (-41398.90057 [K])
    VDW [Guest-Guest]:          -815.79859 (-981.18077 [K])
    Real Coulomb [Host-Host]:   0.00000 (0.00000 [K])
    Real Coulomb [Host-Guest]:  -412.10544 (-495.64922 [K])
    Real Coulomb [Guest-Guest]: 142.47937 (171.36339 [K])
    Ewald [Host-Host]:          0.00000 (0.00000 [K])
    Ewald [Host-Guest]:         -8191.68084 (-9852.33348 [K])
    Ewald [Guest-Guest]:        14.37537 (17.28960 [K])
    DNN Energy:                 0.00000 (0.00000 [K])
    Tail Correction Energy:     0.00000 (0.00000 [K])
    Total Energy:               -43683.67024 (-52539.41104 [K])
    ========================================================================
     *** CHECK DELTA_E (RUNNING FINAL - CREATE MOLECULE) *** 
    ========================================================================
    VDW [Host-Host]:            0.00000 (0.00000 [K])
    VDW [Host-Guest]:           -34420.94011 (-41398.90057 [K])
    VDW [Guest-Guest]:          -815.79859 (-981.18077 [K])
    Real Coulomb [Host-Host]:   0.00000 (0.00000 [K])
    Real Coulomb [Host-Guest]:  -412.10544 (-495.64922 [K])
    Real Coulomb [Guest-Guest]: 142.47937 (171.36339 [K])
    Ewald [Host-Host]:          0.00000 (0.00000 [K])
    Ewald [Host-Guest]:         -8191.68084 (-9852.33348 [K])
    Ewald [Guest-Guest]:        14.37537 (17.28960 [K])
    DNN Energy:                 0.00000 (0.00000 [K])
    Tail Correction Energy:     0.00000 (0.00000 [K])
    Total Energy:               -43683.67024 (-52539.41104 [K])
    ========================================================================
     *** ENERGY DRIFT (CPU FINAL - RUNNING FINAL) *** 
    ========================================================================
    VDW [Host-Host]:            0.00000 (0.00000 [K])
    VDW [Host-Guest]:           -0.00000 (-0.00000 [K])
    VDW [Guest-Guest]:          0.00000 (0.00000 [K])
    Real Coulomb [Host-Host]:   0.00000 (0.00000 [K])
    Real Coulomb [Host-Guest]:  0.00000 (0.00000 [K])
    Real Coulomb [Guest-Guest]: 0.00000 (0.00000 [K])
    Ewald [Host-Host]:          0.00000 (0.00000 [K])
    Ewald [Host-Guest]:         -0.00000 (-0.00000 [K])
    Ewald [Guest-Guest]:        -0.00000 (-0.00000 [K])
    DNN Energy:                 0.00000 (0.00000 [K])
    Tail Correction Energy:     0.00000 (0.00000 [K])
    Total Energy:               -0.00000 (-0.00000 [K])
    ========================================================================
     *** GPU DRIFT (GPU FINAL - CPU FINAL) *** 
    ========================================================================
    VDW [Host-Host]:            0.00000 (0.00000 [K])
    VDW [Host-Guest]:           0.00000 (0.00000 [K])
    VDW [Guest-Guest]:          -0.00000 (-0.00000 [K])
    Real Coulomb [Host-Host]:   0.00000 (0.00000 [K])
    Real Coulomb [Host-Guest]:  0.00000 (0.00000 [K])
    Real Coulomb [Guest-Guest]: -0.00000 (-0.00000 [K])
    Ewald [Host-Host]:          -0.00004 (-0.00005 [K])
    Ewald [Host-Guest]:         0.00000 (0.00000 [K])
    Ewald [Guest-Guest]:        0.00000 (0.00000 [K])
    DNN Energy:                 0.00000 (0.00000 [K])
    Tail Correction Energy:     0.00000 (0.00000 [K])
    Total Energy:               -0.00004 (-0.00005 [K])
    ========================================================================
    ================================================================================
    ======================== PRODUCTION PHASE AVERAGE ENERGIES (Simulation 0) =========================
     *** PRODUCTION PHASE AVERAGE ENERGY *** 
    ========================================================================
    VDW [Host-Host]:            0.00000 (0.00000 [K])
    VDW [Host-Guest]:           0.00000 (0.00000 [K])
    VDW [Guest-Guest]:          0.00000 (0.00000 [K])
    Real Coulomb [Host-Host]:   0.00000 (0.00000 [K])
    Real Coulomb [Host-Guest]:  0.00000 (0.00000 [K])
    Real Coulomb [Guest-Guest]: 0.00000 (0.00000 [K])
    Ewald [Host-Host]:          0.00000 (0.00000 [K])
    Ewald [Host-Guest]:         0.00000 (0.00000 [K])
    Ewald [Guest-Guest]:        0.00000 (0.00000 [K])
    DNN Energy:                 0.00000 (0.00000 [K])
    Tail Correction Energy:     0.00000 (0.00000 [K])
    Total Energy:               0.00000 (0.00000 [K])
    ========================================================================
     *** PRODUCTION PHASE AVERAGE ENERGY ERRORBAR *** 
    ========================================================================
    VDW [Host-Host]:            0.00000 (0.00000 [K])
    VDW [Host-Guest]:           0.00000 (0.00000 [K])
    VDW [Guest-Guest]:          0.00000 (0.00000 [K])
    Real Coulomb [Host-Host]:   0.00000 (0.00000 [K])
    Real Coulomb [Host-Guest]:  0.00000 (0.00000 [K])
    Real Coulomb [Guest-Guest]: 0.00000 (0.00000 [K])
    Ewald [Host-Host]:          0.00000 (0.00000 [K])
    Ewald [Host-Guest]:         0.00000 (0.00000 [K])
    Ewald [Guest-Guest]:        0.00000 (0.00000 [K])
    DNN Energy:                 0.00000 (0.00000 [K])
    Tail Correction Energy:     0.00000 (0.00000 [K])
    Total Energy:               0.00000 (0.00000 [K])
    ========================================================================
    ================================================================================
    DNN Rejection Summary:
    Translation+Rotation: 0
    Reinsertion: 0
    Insertion: 0
    Deletion: 0
    SingleSwap: 0
    DNN Drift Summary:
    Translation+Rotation: 0.00000
    Reinsertion: 0.00000
    Insertion: 0.00000
    Deletion: 0.00000
    SingleSwap: 0.00000
    Summary for simulation 0
    PseudoAtom Type: C[0], #: 0
    PseudoAtom Type: H[1], #: 0
    PseudoAtom Type: N[2], #: 0
    PseudoAtom Type: P[3], #: 0
    PseudoAtom Type: Ni[4], #: 0
    PseudoAtom Type: Ht[5], #: 0
    PseudoAtom Type: OT[6], #: 0
    PseudoAtom Type: Op[7], #: 0
    PseudoAtom Type: Ow[8], #: 0
    PseudoAtom Type: Hw[9], #: 0
    PseudoAtom Type: Lw[10], #: 0
    PseudoAtom Type: Zr[11], #: 0
    PseudoAtom Type: Cl[12], #: 0
    PseudoAtom Type: O[13], #: 1536
    PseudoAtom Type: Si[14], #: 768
    PseudoAtom Type: CH4[15], #: 0
    PseudoAtom Type: C_co2[16], #: 18
    PseudoAtom Type: O_co2[17], #: 36


# Running gRASPA, with your modifications
====================
## To add your modifications, the user needs to modify the gRASPA processes
  * To change what happens during an MC move, one needs to decompose `g.RUN()` function and add their modifications there.
  * Users are given control of these intermediate variables, such as acceptance ratio, selected component/molecule, which move to perform


```python
# run the single-box MC moves through python #
box_index = 0 # run simulation box 0 #
```


```python
# Initialize MC move, this step will initialize a stage of the simulation (Initialization/Equilibration/Production phases)
# For example, it will initialize (clean up) the vectors for storing averages if the stage is production phase
g.InitializeMC(W, box_index)
```

    ==================================
    ==  RUNNING PRODUCTION PHASE   ==
    ==================================
    CBMC Uses 10 trial positions and 10 trial orientations



```python
# Run a move, starting by randomly selecting a component in the box and a molecule in that component
g.Select_Box_Component_Molecule(W, box_index)
```


```python
# Select a move, here we tell the code to run a single-particle insertion move (no CBMC)
MoveType = g.MoveTypes.SINGLE_INSERTION
W.SystemComponents[box_index].MCMoveVariables.MoveType = int(MoveType) # write it back to the variables so that later processes will know
```


```python
W.SystemComponents[box_index].MCMoveVariables.MoveType
```




    2



## Prepare the MC move
* Then to run the new move, first we do the preparation of the mc move
* This includes: 
  * book-keeping the amount of MC move performed
  * initialize trial positions for the insertion move, etc.
* :memo: We can see from the random number offset on the GPU that the preparation of the insertion move costs 3 x double3 of random numbers
  * This corresponds to the random positions in x, y, and z directions


```python
# Prepare the MC move, single particle insertion
print(f"Initial Random offset: {W.Random.offset}\n")
g.SingleBody_Prepare(W, systemId = box_index)
print(f"After Random offset: {W.Random.offset}\n")
```

    Initial Random offset: 28202
    
    After Random offset: 28205
    


## Calculation part of the MC move
This part includes calculating the vdw+real part of the trial energy, ewald summation, tail correction, and acceptance ratio (un-biased)
a `MoveEnergy` is generated, and you can check the deltaE of this move by `DeltaE.print()` and check the singleparticle move Pacc


```python
# Move Calculation #
DeltaE = g.SingleBody_Calculation(W, systemId = box_index)
DeltaE.print()
```

    HHVDW: 0.00000, HHReal: 0.00000, HGVDW: 6448.46706, HGReal: 500.24272, GGVDW: -344.32224, GGReal: 98.20524, HHEwaldE: 0.00000,
     HGEwaldE: -365.77690,
     GGEwaldE: 17.75031, TailE: 0.00000, DNN_E: 0.00000
    Stored HGVDW: 0.00000, Stored HGReal: 0.00000, Stored HGEwaldE: 0.00000



```python
W.SystemComponents[box_index].MCMoveVariables.Overlap
W.SystemComponents[box_index].MCMoveVariables.Pacc
```




    3.922521744437901e-14



# Move Acceptance #
* Dictates whether a move is accepted by Metropolis scheme
* You can then check acceptance by looking at the temporary variables


```python
# Move Acceptance #
comp = W.SystemComponents[box_index].MCMoveVariables.component
print(f"number of molecules: {W.SystemComponents[box_index].NumberOfMolecules[comp]}")
g.SingleBody_Acceptance(W, box_index, DeltaE)
print(f"Accept?: {W.SystemComponents[box_index].MCMoveVariables.Accept}\n")
print(f"number of molecules: {W.SystemComponents[box_index].NumberOfMolecules[comp]}")
```

    number of molecules: 18
    Accept?: False
    
    number of molecules: 18


## That move was rejected and number of molecules stays the same 
* We can see Pacc is small
* DeltaE has a really high HGVDW value, indicating unfavorable host-guest interactions
## What if we force the move to be accepted?
* We can modify this by changing the Pacc to something big


```python
# The move is rejected, what if we can change that manually? #
W.SystemComponents[box_index].MCMoveVariables.Pacc = 0.95
g.SingleBody_Acceptance(W, box_index, DeltaE)
print(f"Accept?: {W.SystemComponents[box_index].MCMoveVariables.Accept}\n")
print(f"number of molecules: {W.SystemComponents[box_index].NumberOfMolecules[comp]}")
```

    Accept?: True
    
    number of molecules: 19


## Now the move is forced to be accepted!
