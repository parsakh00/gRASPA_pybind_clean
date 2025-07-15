# gRASPA and numpy already included in utils
from utils import *

# run steps independently #
W = g.Initialize()
W.NumberOfInitializationCycles = 100
W.NumberOfProductionCycles = 0
g.RUN(W)
g.finalize(W)

# run the single-box MC moves through python #
box_index = 0 # run simulation box 0 #
g.InitializeMC(W, box_index)

g.Select_Box_Component_Molecule(W, box_index)
# run single MC move #
r = W.MCMoveVariables.RandomNumber
comp = W.MCMoveVariables.component
mol  = W.MCMoveVariables.molecule

#NOTE: decomment each one of the three to test them#
#MoveType = g.MoveTypes.SINGLE_INSERTION
#MoveType = g.MoveTypes.SINGLE_DELETION
#MoveType = g.MoveTypes.ROTATION
MoveType = g.MoveTypes.TRANSLATION

W.MCMoveVariables.MoveType = int(MoveType)

# PERFORM MOVE #
g.SingleBody_Prepare(W, systemId = box_index)
DeltaE = g.SingleBody_Calculation(W, systemId = box_index)

# Get Adsorbate data (old)
# Get framework data, if rigid, should just run it once!#
# systemId = 0 (only one simulation box) #
# component 0 = framework (rigid), component 1 = adsorbates (CO2) for our case#
# if rigid, only need to run get FrameworkData once!#
FrameworkData = g.GetAllAtoms(W, systemId = 0, component = 0)
# Get trial positions (new), dict will get entries: "Trial_pos", "Trial_charge", ...
# If WholeConfig = False, then only the moved molecule is here #
# else, the moved molecule will be merged so that all molecules + trial molecule (moved) is copied here #
OldData = g.GetAllAtoms(W, systemId = 0, component = 1)
TrialData = g.GetTrialConfig(W, systemId = 0, component = 1, WholeConfig = False)
TrialAllData = g.GetTrialConfig(W, systemId = 0, component = 1, WholeConfig = True)
print(f"OldData: {OldData}\n")
print(f"TrialData: {TrialData}\n")
print(f"Data: {TrialAllData}\n")

Combined_OldData = CombineFrameworkAdsorbate(FrameworkData, OldData)
CombinedData = CombineFrameworkAdsorbate(FrameworkData, TrialAllData)
print(f"{CombinedData['pos'][-1]}\n")

WriteDataToPOSCAR(W, systemId = 0, DATA = FrameworkData, filename = "Framework")
WriteDataToPOSCAR(W, systemId = 0, DATA = CombinedData, filename = "TRIAL")
WriteDataToPOSCAR(W, systemId = 0, DATA = Combined_OldData, filename = "OLD")
# combine old config with trial config

g.SingleBody_Acceptance(W, box_index, DeltaE)

print(f"Accept?: {W.MCMoveVariables.Accept}\n")


