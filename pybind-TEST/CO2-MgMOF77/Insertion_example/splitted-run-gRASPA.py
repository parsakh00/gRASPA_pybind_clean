# gRASPA and numpy already included in utils
from utils import *

# run steps independently #
W = g.Initialize()
W.NumberOfInitializationCycles = 10 # as a test, just run 10 cycles #
# run the single-box MC moves through python #
box_index = 0 # run simulation box 0 #
g.InitializeMC(W, box_index)

# Save framework data , framework component = 0#
FrameworkData = g.GetAllAtoms(W, systemId = 0, component = 0)

start = g.omp_get_wtime()
for i in range(0, W.NumberOfInitializationCycles):
  Steps = g.Determine_Number_Of_Steps(W, box_index, i)
  for s in range(0, Steps):
    g.Select_Box_Component_Molecule(W, box_index)
    # run single MC move #
    r = W.MCMoveVariables.RandomNumber
    comp = W.MCMoveVariables.component
    mol  = W.MCMoveVariables.molecule

    # Select move #
    # for translation/rotation/deletion, if there is no molecule, skip #
    MoveType = g.MoveTypes.TRANSLATION
    if(r < W.SystemComponents[box_index].Moves[comp].TranslationProb):
      MoveType = g.MoveTypes.TRANSLATION
      if(W.SystemComponents[box_index].NumberOfMolecules[comp] == 0): continue
    elif(r < W.SystemComponents[box_index].Moves[comp].RotationProb):
      MoveType = g.MoveTypes.ROTATION
      if(W.SystemComponents[box_index].NumberOfMolecules[comp] == 0): continue
    elif(r < W.SystemComponents[box_index].Moves[comp].SwapProb):
      if(g.Get_Uniform_Random() < 0.5):
        MoveType = g.MoveTypes.SINGLE_INSERTION
      else:
        MoveType = g.MoveTypes.SINGLE_DELETION
        if(W.SystemComponents[box_index].NumberOfMolecules[comp] == 0): continue
    W.MCMoveVariables.MoveType = int(MoveType)

    # PERFORM MOVE #
    g.SingleBody_Prepare(W, systemId = box_index)
    DeltaE = g.SingleBody_Calculation(W, systemId = box_index)
    
    # Add your modifications here, change the move energy#
    TrialConfigData = g.GetTrialConfig(W, systemId = 0, component = 1, WholeConfig = True)
    # Combine Framework data and trial adsorbate data #
    CombinedData = CombineFrameworkAdsorbate(FrameworkData, TrialConfigData)
    # Save your data (FrameworkData + TrialConfigData) to POSCAR #
    WriteDataToPOSCAR(W, systemId = 0, DATA = CombinedData, filename = "TRIAL")
    # Read "TRIAL.vasp" #
    # Then get user-defined energy #
    # XXX = predict(?)
    # Modify DeltaE based on your energy #
    # You can check DeltaE via : DeltaE.print()
    #DeltaE.zero() # zero the DeltaE (vdw E, coulomb E, ...)
    #DeltaE.DNN_E = XXX

    # feed the modified DeltaE back to acceptance rules #
    g.SingleBody_Acceptance(W, box_index, DeltaE)
    # add DeltaE (move delta) to deltaE (total delta) #
    g.MoveEnergy_Add(W.SystemComponents[box_index].deltaE, DeltaE)
  # Gather statistics and averages per cycle #
  g.GatherStatisticsDuringSimulation(W, box_index, i)
g.MCEndOfPhaseSummary(W) # end of initialization phase summary #
end = g.omp_get_wtime()
print(f"simulation took {end - start} seconds\n")

g.finalize(W) # end of simulation #
