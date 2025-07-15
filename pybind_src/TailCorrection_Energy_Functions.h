//Need to consider if we use CBCF//

double TotalTailCorrection(Components& SystemComponents, size_t FFsize, double Volume)
{
  double TailE = 0.0;
  if(!SystemComponents.HasTailCorrection) return TailE;
  for(size_t i = 0; i < SystemComponents.NumberOfPseudoAtoms.size(); i++)
    for(size_t j = i; j < SystemComponents.NumberOfPseudoAtoms.size(); j++)
    {
      size_t IJ_Forward = i * FFsize + j;
      if(SystemComponents.TailCorrection[IJ_Forward].UseTail)
      {
        size_t Ni    = SystemComponents.NumberOfPseudoAtoms[i];
        size_t Nj    = SystemComponents.NumberOfPseudoAtoms[j];
        double pairE = SystemComponents.TailCorrection[IJ_Forward].Energy * static_cast<double>(Ni * Nj); if(i != j) pairE *= 2.0;
        TailE       += pairE;
        //printf("i: %zu, j: %zu, Ni: %zu, Nj: %zu, E: %.5f\n", i,j,Ni,Nj, pairE / Volume);
      }
    }
  return TailE / Volume;
}

size_t get_change_in_pseudoAtoms(Components& SystemComponents, size_t SelectedComponent, size_t Type)
{
  size_t d = 0;
  for(size_t i = 0; i < SystemComponents.NumberOfPseudoAtomsForSpecies[SelectedComponent].size(); i++)
  {
    if(Type == SystemComponents.NumberOfPseudoAtomsForSpecies[SelectedComponent][i].x)
    {
      d = SystemComponents.NumberOfPseudoAtomsForSpecies[SelectedComponent][i].y; break;
    }
  }
  return d;
}

double TailCorrectionDifference(Components& SystemComponents, size_t SelectedComponent, size_t FFsize, double Volume, int MoveType)
{
  double TailE = 0.0;
  if(!SystemComponents.HasTailCorrection) return TailE;
  int sign = 1;
  switch(MoveType)
  {
    case INSERTION:
    {
      sign = 1;
      break;
    }
    case CBCF_INSERTION:
    {
      throw std::runtime_error("TAIL CORRECTIONS NOT READY FOR CBCF MOVES!");
    }
    case CBCF_DELETION:
    {
      throw std::runtime_error("TAIL CORRECTIONS NOT READY FOR CBCF MOVES!");
    }
    case DELETION:
    {
      sign = -1;
      break;
    }
  }
  for(size_t i = 0; i < SystemComponents.NumberOfPseudoAtoms.size(); i++)
  {
    int di = get_change_in_pseudoAtoms(SystemComponents, SelectedComponent, i); di *= sign;

    for(size_t j = i; j < SystemComponents.NumberOfPseudoAtoms.size(); j++)
    {
      int dj = get_change_in_pseudoAtoms(SystemComponents, SelectedComponent, j); dj *= sign;
      size_t IJ_Forward = i * FFsize + j;
      if(SystemComponents.TailCorrection[IJ_Forward].UseTail)
      {
        int    Ni = SystemComponents.NumberOfPseudoAtoms[i];
        int    Nj = SystemComponents.NumberOfPseudoAtoms[j];
        double E  = SystemComponents.TailCorrection[IJ_Forward].Energy;
        int    dN = Ni * dj + Nj * di + di * dj; //Zhao's note: need to define a temporary variable for this first//
        double DeltaE = E * static_cast<double>(dN); if(i != j) DeltaE *= 2.0;
        TailE    +=     DeltaE;
      }
    }
  }
  return TailE / Volume;
}

//For moves that simultaneously changes more than one component//
double TailCorrectionIdentitySwap(Components& SystemComponents, size_t NEWComponent, size_t OLDComponent, size_t FFsize, double Volume)
{
  double TailE = 0.0;
  if(!SystemComponents.HasTailCorrection) return TailE;

  for(size_t i = 0; i < SystemComponents.NumberOfPseudoAtoms.size(); i++)
  {
    int di = get_change_in_pseudoAtoms(SystemComponents, NEWComponent, i); 
        di-= get_change_in_pseudoAtoms(SystemComponents, OLDComponent, i);

    for(size_t j = i; j < SystemComponents.NumberOfPseudoAtoms.size(); j++)
    {
      int dj = get_change_in_pseudoAtoms(SystemComponents, NEWComponent, j); 
          dj-= get_change_in_pseudoAtoms(SystemComponents, OLDComponent, j);

      size_t IJ_Forward = i * FFsize + j;
      if(SystemComponents.TailCorrection[IJ_Forward].UseTail)
      {
        int    Ni = SystemComponents.NumberOfPseudoAtoms[i];
        int    Nj = SystemComponents.NumberOfPseudoAtoms[j];
        double E  = SystemComponents.TailCorrection[IJ_Forward].Energy;
        int    dN = Ni * dj + Nj * di + di * dj; //Zhao's note: need to define a temporary variable for this first//
        double DeltaE = E * static_cast<double>(dN); if(i != j) DeltaE *= 2.0;
        TailE    +=     DeltaE;
      }
    }
  }
  return TailE / Volume;
}
