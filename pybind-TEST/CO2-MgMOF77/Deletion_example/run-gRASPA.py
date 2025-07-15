import gRASPA as g
# save simulation to this variable "V" (V stands for Variable, the dictating variable in gRASPA code that wraps everything)
# run the whole simulation #
V = g.Initialize()
V.NumberOfInitializationCycles = 10000
g.RUN(V)
g.finalize(V)

'''
# save atom data for system 0, component 1 (co2 in this case)
DataA = g.GetAllAtoms(V, 0, 1)
# Run some MC moves #
# run multiple times of rotation moves, on the 10th CO2 molecule #
g.SingleBody_Prepare(V, 0, 10, 1, 1) # prepare the MC move
DeltaE = g.SingleBody_Calculation(V, 0, 10, 1, 1) # calculate the classical energies (vdW + electrostatics)
DeltaE.print() # check the deltaE for the rotation move
g.SingleBody_Acceptance(V, 0, 10, 1, 1, DeltaE) # determine whether to accept or reject the move

# run it multiple times ...

# calculate current total energy
TotalE = g.get_total_energy(V, 0)
# print total energy
TotalE.print()
'''

