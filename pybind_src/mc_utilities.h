////////////////////////////////////////////////////////////
// FUNCATIONS RELATED TO ALLOCATING MORE SPACE ON THE GPU //
////////////////////////////////////////////////////////////

__global__ void AllocateMoreSpace_CopyToTemp(Atoms* d_a, Atoms temp, size_t Space, size_t SelectedComponent)
{
  for(size_t i = 0; i < Space; i++)
  {
    temp.pos[i]       = d_a[SelectedComponent].pos[i];
    temp.scale[i]     = d_a[SelectedComponent].scale[i];
    temp.charge[i]    = d_a[SelectedComponent].charge[i];
    temp.scaleCoul[i] = d_a[SelectedComponent].scaleCoul[i];
    temp.Type[i]      = d_a[SelectedComponent].Type[i];
    temp.MolID[i]     = d_a[SelectedComponent].MolID[i];
  }
}

__global__ void AllocateMoreSpace_CopyBack(Atoms* d_a, Atoms temp, size_t Space, size_t Newspace, size_t SelectedComponent)
{
  d_a[SelectedComponent].Allocate_size = Newspace;
  for(size_t i = 0; i < Space; i++)
  {
    d_a[SelectedComponent].pos[i]       = temp.pos[i];
    d_a[SelectedComponent].scale[i]     = temp.scale[i];
    d_a[SelectedComponent].charge[i]    = temp.charge[i];
    d_a[SelectedComponent].scaleCoul[i] = temp.scaleCoul[i];
    d_a[SelectedComponent].Type[i]      = temp.Type[i];
    d_a[SelectedComponent].MolID[i]     = temp.MolID[i];
  /*
    System_comp.pos[i]         = temp.pos[i];
    System_comp.scale[i]     = temp.scale[i];
    System_comp.charge[i]    = temp.charge[i];
    System_comp.scaleCoul[i] = temp.scaleCoul[i];
    System_comp.Type[i]      = temp.Type[i];
    System_comp.MolID[i]     = temp.MolID[i];
  */
  }
  //test the new allocation//
  //d_a[SelectedComponent].pos[Newspace-1].x = 0.0;
  /*if(d_a[SelectedComponent].pos[Newspace-1].x = 0.0)
  {
  }*/
}

void AllocateMoreSpace(Atoms*& d_a, size_t SelectedComponent, Components& SystemComponents)
{
  printf("Allocating more space on device\n");
  Atoms temp; // allocate a struct on the device for copying data.
  //Atoms tempSystem[SystemComponents.NComponents.x];
  size_t Copysize=SystemComponents.Allocate_size[SelectedComponent];
  size_t Morespace = 1024;
  size_t Newspace = Copysize+Morespace;
  //Allocate space on the temporary struct
  cudaMalloc(&temp.pos,       Copysize * sizeof(double3));
  cudaMalloc(&temp.scale,     Copysize * sizeof(double));
  cudaMalloc(&temp.charge,    Copysize * sizeof(double));
  cudaMalloc(&temp.scaleCoul, Copysize * sizeof(double));
  cudaMalloc(&temp.Type,      Copysize * sizeof(size_t));
  cudaMalloc(&temp.MolID,     Copysize * sizeof(size_t));
  // Copy data to temp
  AllocateMoreSpace_CopyToTemp<<<1,1>>>(d_a, temp, Copysize, SelectedComponent);
  // Allocate more space on the device pointers

  Atoms System[SystemComponents.NComponents.x]; cudaMemcpy(System, d_a, SystemComponents.NComponents.x * sizeof(Atoms), cudaMemcpyDeviceToHost);

  cudaMalloc(&System[SelectedComponent].pos,       Newspace * sizeof(double3));
  cudaMalloc(&System[SelectedComponent].scale,     Newspace * sizeof(double));
  cudaMalloc(&System[SelectedComponent].charge,    Newspace * sizeof(double));
  cudaMalloc(&System[SelectedComponent].scaleCoul, Newspace * sizeof(double));
  cudaMalloc(&System[SelectedComponent].Type,      Newspace * sizeof(size_t));
  cudaMalloc(&System[SelectedComponent].MolID,     Newspace * sizeof(size_t));
  // Copy data from temp back to the new pointers
  AllocateMoreSpace_CopyBack<<<1,1>>>(d_a, temp, Copysize, Newspace, SelectedComponent);
  //System[SelectedComponent].Allocate_size = Newspace;
}

///////////////////////////////////////////
// FUNCTIONS RELATED TO ACCEPTING A MOVE //
///////////////////////////////////////////

static inline void Update_NumberOfMolecules(Components& SystemComponents, Atoms*& d_a, size_t SelectedComponent, int MoveType)
{
  size_t Molsize = SystemComponents.Moleculesize[SelectedComponent]; //change in atom number counts
  int NumMol = -1; //default: deletion; Insertion: +1, Deletion: -1, size_t is never negative

  switch(MoveType)
  {
    case INSERTION: case SINGLE_INSERTION: case CBCF_INSERTION:
    {
      NumMol = 1;  break;
    }
    case DELETION: case SINGLE_DELETION: case CBCF_DELETION:
    {
      NumMol = -1; break;
    }
  }
  //Update Components
  SystemComponents.NumberOfMolecule_for_Component[SelectedComponent] += NumMol;
  SystemComponents.TotalNumberOfMolecules += NumMol;

  // check System size //
  size_t size = SystemComponents.Moleculesize[SelectedComponent] * SystemComponents.NumberOfMolecule_for_Component[SelectedComponent];

  SystemComponents.UpdatePseudoAtoms(MoveType,  SelectedComponent);

  if(size > SystemComponents.Allocate_size[SelectedComponent])
  {
    AllocateMoreSpace(d_a, SelectedComponent, SystemComponents);
    throw std::runtime_error("Need to allocate more space, not implemented\n");
  }
}

__global__ void Update_SINGLE_INSERTION_data(Atoms* d_a, Atoms New, size_t SelectedComponent)
{
  //Assuming single thread//
  size_t Molsize = d_a[SelectedComponent].Molsize;
  size_t UpdateLocation = d_a[SelectedComponent].size;
  for(size_t j = 0; j < Molsize; j++)
  {
    d_a[SelectedComponent].pos[UpdateLocation+j]       = New.pos[j];
    d_a[SelectedComponent].scale[UpdateLocation+j]     = New.scale[j];
    d_a[SelectedComponent].charge[UpdateLocation+j]    = New.charge[j];
    d_a[SelectedComponent].scaleCoul[UpdateLocation+j] = New.scaleCoul[j];
    d_a[SelectedComponent].Type[UpdateLocation+j]      = New.Type[j];
    d_a[SelectedComponent].MolID[UpdateLocation+j]     = New.MolID[j];
  }
  d_a[SelectedComponent].size  += Molsize;
}

__global__ void Update_deletion_data(Atoms* d_a, size_t SelectedComponent, size_t UpdateLocation, int Moleculesize, size_t LastLocation)
{
  size_t i = blockIdx.x * blockDim.x + threadIdx.x;

  //UpdateLocation should be the molecule that needs to be deleted
  //Then move the atom at the last position to the location of the deleted molecule
  //**Zhao's note** MolID of the deleted molecule should not be changed
  //**Zhao's note** if Molecule deleted is the last molecule, then nothing is copied, just change the size later.
  if(UpdateLocation != LastLocation)
  {
    for(size_t i = 0; i < Moleculesize; i++)
    {
      d_a[SelectedComponent].pos[UpdateLocation+i]       = d_a[SelectedComponent].pos[LastLocation+i];
      d_a[SelectedComponent].scale[UpdateLocation+i]     = d_a[SelectedComponent].scale[LastLocation+i];
      d_a[SelectedComponent].charge[UpdateLocation+i]    = d_a[SelectedComponent].charge[LastLocation+i];
      d_a[SelectedComponent].scaleCoul[UpdateLocation+i] = d_a[SelectedComponent].scaleCoul[LastLocation+i];
      d_a[SelectedComponent].Type[UpdateLocation+i]      = d_a[SelectedComponent].Type[LastLocation+i];
    }
  }
  //Zhao's note: the single values in d_a and System are pointing to different locations//
  //d_a is just device (cannot access on host), while System is shared (accessible on host), need to update d_a values here
  //there are two of these values: size and Allocate_size
  if(i==0)
  {
    d_a[SelectedComponent].size  -= Moleculesize; //Zhao's special note: AllData.size doesn't work... So single values are painful, need to consider pointers for single values
  }
}

__global__ void Update_deletion_data_Parallel(Atoms* d_a, size_t SelectedComponent, size_t UpdateLocation, int Moleculesize, size_t LastLocation)
{
  size_t i = blockIdx.x * blockDim.x + threadIdx.x;
  
  //UpdateLocation should be the molecule that needs to be deleted
  //Then move the atom at the last position to the location of the deleted molecule
  //**Zhao's note** MolID of the deleted molecule should not be changed
  //**Zhao's note** if Molecule deleted is the last molecule, then nothing is copied, just change the size later.
  if(UpdateLocation != LastLocation)
  {
    if(i < Moleculesize)
    {
      d_a[SelectedComponent].pos[UpdateLocation+i]       = d_a[SelectedComponent].pos[LastLocation+i];
      d_a[SelectedComponent].scale[UpdateLocation+i]     = d_a[SelectedComponent].scale[LastLocation+i];
      d_a[SelectedComponent].charge[UpdateLocation+i]    = d_a[SelectedComponent].charge[LastLocation+i];
      d_a[SelectedComponent].scaleCoul[UpdateLocation+i] = d_a[SelectedComponent].scaleCoul[LastLocation+i];
      d_a[SelectedComponent].Type[UpdateLocation+i]      = d_a[SelectedComponent].Type[LastLocation+i];
    }
  }
  //Zhao's note: the single values in d_a and System are pointing to different locations//
  //d_a is just device (cannot access on host), while System is shared (accessible on host), need to update d_a values here
  //there are two of these values: size and Allocate_size
  if(i==0)
  {
    d_a[SelectedComponent].size  -= Moleculesize; //Zhao's special note: AllData.size doesn't work... So single values are painful, need to consider pointers for single values
  }
} 

__global__ void Update_insertion_data_Parallel(Atoms* d_a, Atoms Mol, Atoms NewMol, size_t SelectedTrial, size_t SelectedComponent, size_t UpdateLocation, int Moleculesize)
{
  size_t i = blockIdx.x * blockDim.x + threadIdx.x;
  //UpdateLocation should be the last position of the dataset
  //Need to check if Allocate_size is smaller than size
  if(i == 0)
  {
    //Zhao's note: the single values in d_a and System are pointing to different locations//
    //d_a is just device (cannot access on host), while System is shared (accessible on host), need to update d_a values here
    //there are two of these values: size and Allocate_size
    d_a[SelectedComponent].size  += Moleculesize;
	  if(Moleculesize == 1) //Only first bead is inserted, first bead data is stored in NewMol
	  {
	    d_a[SelectedComponent].pos[UpdateLocation]       = NewMol.pos[SelectedTrial];
	    d_a[SelectedComponent].scale[UpdateLocation]     = NewMol.scale[SelectedTrial];
	    d_a[SelectedComponent].charge[UpdateLocation]    = NewMol.charge[SelectedTrial];
	    d_a[SelectedComponent].scaleCoul[UpdateLocation] = NewMol.scaleCoul[SelectedTrial];
	    d_a[SelectedComponent].Type[UpdateLocation]      = NewMol.Type[SelectedTrial];
	    d_a[SelectedComponent].MolID[UpdateLocation]     = NewMol.MolID[SelectedTrial];
	  }
	  else //Multiple beads: first bead + trial orientations
	  {
	    //Update the first bead, first bead data stored in position 0 of Mol //
	    d_a[SelectedComponent].pos[UpdateLocation]       = Mol.pos[0];
	    d_a[SelectedComponent].scale[UpdateLocation]     = Mol.scale[0];
	    d_a[SelectedComponent].charge[UpdateLocation]    = Mol.charge[0];
	    d_a[SelectedComponent].scaleCoul[UpdateLocation] = Mol.scaleCoul[0];
	    d_a[SelectedComponent].Type[UpdateLocation]      = Mol.Type[0];
	    d_a[SelectedComponent].MolID[UpdateLocation]     = Mol.MolID[0];
	  }
  }
  else if(i < Moleculesize)
  {
      size_t chainsize = Moleculesize - 1; size_t j = i - 1;
      size_t selectsize = SelectedTrial*chainsize+j;
      d_a[SelectedComponent].pos[UpdateLocation+j+1]       = NewMol.pos[selectsize];
      d_a[SelectedComponent].scale[UpdateLocation+j+1]     = NewMol.scale[selectsize];
      d_a[SelectedComponent].charge[UpdateLocation+j+1]    = NewMol.charge[selectsize];
      d_a[SelectedComponent].scaleCoul[UpdateLocation+j+1] = NewMol.scaleCoul[selectsize];
      d_a[SelectedComponent].Type[UpdateLocation+j+1]      = NewMol.Type[selectsize];
      d_a[SelectedComponent].MolID[UpdateLocation+j+1]     = NewMol.MolID[selectsize];
  }
}


__global__ void Update_insertion_data(Atoms* d_a, Atoms Mol, Atoms NewMol, size_t SelectedTrial, size_t SelectedComponent, size_t UpdateLocation, int Moleculesize)
{
  size_t i = blockIdx.x * blockDim.x + threadIdx.x;
  //UpdateLocation should be the last position of the dataset
  //Need to check if Allocate_size is smaller than size
  if(Moleculesize == 1) //Only first bead is inserted, first bead data is stored in NewMol
  {
    d_a[SelectedComponent].pos[UpdateLocation]       = NewMol.pos[SelectedTrial];
    d_a[SelectedComponent].scale[UpdateLocation]     = NewMol.scale[SelectedTrial];
    d_a[SelectedComponent].charge[UpdateLocation]    = NewMol.charge[SelectedTrial];
    d_a[SelectedComponent].scaleCoul[UpdateLocation] = NewMol.scaleCoul[SelectedTrial];
    d_a[SelectedComponent].Type[UpdateLocation]      = NewMol.Type[SelectedTrial];
    d_a[SelectedComponent].MolID[UpdateLocation]     = NewMol.MolID[SelectedTrial];
  }
  else //Multiple beads: first bead + trial orientations
  {
    //Update the first bead, first bead data stored in position 0 of Mol //
    d_a[SelectedComponent].pos[UpdateLocation]       = Mol.pos[0];
    d_a[SelectedComponent].scale[UpdateLocation]     = Mol.scale[0];
    d_a[SelectedComponent].charge[UpdateLocation]    = Mol.charge[0];
    d_a[SelectedComponent].scaleCoul[UpdateLocation] = Mol.scaleCoul[0];
    d_a[SelectedComponent].Type[UpdateLocation]      = Mol.Type[0];
    d_a[SelectedComponent].MolID[UpdateLocation]     = Mol.MolID[0];

    size_t chainsize = Moleculesize - 1; // For trial orientations //
    for(size_t j = 0; j < chainsize; j++) //Update the selected orientations//
    {
      size_t selectsize = SelectedTrial*chainsize+j;
      d_a[SelectedComponent].pos[UpdateLocation+j+1]       = NewMol.pos[selectsize];
      d_a[SelectedComponent].scale[UpdateLocation+j+1]     = NewMol.scale[selectsize];
      d_a[SelectedComponent].charge[UpdateLocation+j+1]    = NewMol.charge[selectsize];
      d_a[SelectedComponent].scaleCoul[UpdateLocation+j+1] = NewMol.scaleCoul[selectsize];
      d_a[SelectedComponent].Type[UpdateLocation+j+1]      = NewMol.Type[selectsize];
      d_a[SelectedComponent].MolID[UpdateLocation+j+1]     = NewMol.MolID[selectsize];
    }
  }
  //Zhao's note: the single values in d_a and System are pointing to different locations//
  //d_a is just device (cannot access on host), while System is shared (accessible on host), need to update d_a values here
  //there are two of these values: size and Allocate_size
  if(i==0)
  {
    d_a[SelectedComponent].size  += Moleculesize;
    /*
    for(size_t j = 0; j < Moleculesize; j++)
      printf("xyz: %.5f %.5f %.5f, scale/charge/scaleCoul: %.5f %.5f %.5f, Type: %lu, MolID: %lu\n", d_a[SelectedComponent].pos[UpdateLocation+j].x, d_a[SelectedComponent].pos[UpdateLocation+j].y, d_a[SelectedComponent].pos[UpdateLocation+j].z, d_a[SelectedComponent].scale[UpdateLocation+j], d_a[SelectedComponent].charge[UpdateLocation+j], d_a[SelectedComponent].scaleCoul[UpdateLocation+j], d_a[SelectedComponent].Type[UpdateLocation+j], d_a[SelectedComponent].MolID[UpdateLocation+j]);
    */
  }
}

__global__ void update_translation_position(Atoms* d_a, Atoms NewMol, size_t start_position, size_t SelectedComponent)
{
  size_t i = blockIdx.x * blockDim.x + threadIdx.x;
  d_a[SelectedComponent].pos[start_position+i] = NewMol.pos[i];
  d_a[SelectedComponent].scale[start_position+i] = NewMol.scale[i];
  d_a[SelectedComponent].charge[start_position+i] = NewMol.charge[i];
  d_a[SelectedComponent].scaleCoul[start_position+i] = NewMol.scaleCoul[i];
}

void AcceptTranslation(Variables& Vars, size_t systemId)
{
  Components& SystemComponents = Vars.SystemComponents[systemId];
  Simulations& Sims            = Vars.Sims[systemId];
  ForceField& FF               = Vars.device_FF;

  size_t& SelectedComponent = SystemComponents.TempVal.component;

  size_t Molsize = SystemComponents.Moleculesize[SelectedComponent]; //Get the size of the selected Molecule
  size_t& start_position = SystemComponents.TempVal.start_position;
  update_translation_position<<<1,Molsize>>>(Sims.d_a, Sims.New, start_position, SelectedComponent);
  if(!FF.noCharges && SystemComponents.hasPartialCharge[SelectedComponent])
  {
    Update_Vector_Ewald(Sims.Box, false, SystemComponents, SelectedComponent);
  }
}

////////////////////////////////////////////////
// GET PREFACTOR FOR INSERTION/DELETION MOVES //
////////////////////////////////////////////////

inline double GetPrefactor(Components& SystemComponents, Simulations& Sims, size_t SelectedComponent, int MoveType)
{
  double MolFraction = SystemComponents.MolFraction[SelectedComponent];
  double FugacityCoefficient = SystemComponents.FugacityCoeff[SelectedComponent];
  double NumberOfMolecules = static_cast<double>(SystemComponents.NumberOfMolecule_for_Component[SelectedComponent]);

  //If component has fractional molecule, subtract the number of molecules by 1.//
  if(SystemComponents.hasfractionalMolecule[SelectedComponent]){NumberOfMolecules-=1.0;}
  if(NumberOfMolecules < 0.0) NumberOfMolecules = 0.0;

  double preFactor = 0.0;
  //Hard to generalize the prefactors, when you consider identity swap
  //Certainly you can use insertion/deletion for identity swap, but you may lose accuracies
  //since you are multiplying then dividing by the same variables
  switch(MoveType)
  {
    case INSERTION: case SINGLE_INSERTION:
    {
      preFactor = SystemComponents.Beta * MolFraction * SystemComponents.Pressure * FugacityCoefficient * Sims.Box.Volume / (1.0+NumberOfMolecules);
      break;
    }
    case DELETION: case SINGLE_DELETION:
    {
      preFactor = (NumberOfMolecules) / (SystemComponents.Beta * MolFraction * SystemComponents.Pressure * FugacityCoefficient * Sims.Box.Volume);
      break;
    }
    case IDENTITY_SWAP:
    {
      throw std::runtime_error("Sorry, but no IDENTITY_SWAP OPTION for now. That move uses its own function\n");
    }
    case TRANSLATION: case ROTATION: case SPECIAL_ROTATION: 
    {
      preFactor = 1.0;
    }
  }
  return preFactor;
}

////////////////////////
// ACCEPTION OF MOVES //
////////////////////////

void AcceptInsertion(Variables& Vars, CBMC_Variables& InsertionVariables, size_t systemId, int MoveType)
{
  Components& SystemComponents = Vars.SystemComponents[systemId];
  Simulations& Sims            = Vars.Sims[systemId];
  ForceField& FF               = Vars.device_FF;

  //int& MoveType   = SystemComponents.TempVal.MoveType;
  size_t& SelectedComponent = SystemComponents.TempVal.component;
  bool& noCharges = FF.noCharges;
  //printf("AccInsertion, SelectedTrial: %zu, UpdateLocation: %zu\n", SelectedTrial, UpdateLocation);
  //Zhao's note: here needs more consideration: need to update after implementing polyatomic molecule
  if(MoveType == INSERTION)
  {
    //CBMC_Variables& InsertionVariables = SystemComponents.CBMC_New[0];
    size_t UpdateLocation = SystemComponents.Moleculesize[SelectedComponent] * SystemComponents.NumberOfMolecule_for_Component[SelectedComponent];
    size_t SelectedTrial = InsertionVariables.selectedTrial;
    if(SystemComponents.Moleculesize[SelectedComponent] > 1) SelectedTrial = InsertionVariables.selectedTrialOrientation;
    Update_insertion_data_Parallel<<<1,SystemComponents.Moleculesize[SelectedComponent]>>>(Sims.d_a, Sims.Old, Sims.New, SelectedTrial, SelectedComponent, UpdateLocation, (int) SystemComponents.Moleculesize[SelectedComponent]);
  }
  else if(MoveType == SINGLE_INSERTION)
  {
    Update_SINGLE_INSERTION_data<<<1,1>>>(Sims.d_a, Sims.New, SelectedComponent);
  }
  Update_NumberOfMolecules(SystemComponents, Sims.d_a, SelectedComponent, INSERTION); //true = Insertion//
  if(!noCharges && SystemComponents.hasPartialCharge[SelectedComponent])
  {
    Update_Vector_Ewald(Sims.Box, false, SystemComponents, SelectedComponent);
  }
}

inline void AcceptDeletion(Variables& Vars, size_t systemId, int MoveType)
{
  Components& SystemComponents = Vars.SystemComponents[systemId];
  Simulations& Sims            = Vars.Sims[systemId];
  ForceField& FF               = Vars.device_FF;
  bool& noCharges = FF.noCharges;

  size_t& SelectedComponent = SystemComponents.TempVal.component;

  size_t LastMolecule = SystemComponents.NumberOfMolecule_for_Component[SelectedComponent]-1;
  size_t LastLocation = LastMolecule*SystemComponents.Moleculesize[SelectedComponent];

  size_t& UpdateLocation = SystemComponents.TempVal.UpdateLocation;
  if(MoveType == SINGLE_DELETION) 
    UpdateLocation = SystemComponents.TempVal.molecule * SystemComponents.Moleculesize[SelectedComponent];

  Update_deletion_data_Parallel<<<1,SystemComponents.Moleculesize[SelectedComponent]>>>(Sims.d_a, SelectedComponent, UpdateLocation, (int) SystemComponents.Moleculesize[SelectedComponent], LastLocation);

  Update_NumberOfMolecules(SystemComponents, Sims.d_a, SelectedComponent, DELETION); //false = Deletion//
  if(!noCharges && SystemComponents.hasPartialCharge[SelectedComponent])
  {
    Update_Vector_Ewald(Sims.Box, false, SystemComponents, SelectedComponent);
  }
  //Zhao's note: the last molecule can be the fractional molecule, (fractional molecule ID is stored on the host), we need to update it as well (at least check it)//
  //The function below will only be processed if the system has a fractional molecule and the transfered molecule is NOT the fractional one //
  if((SystemComponents.hasfractionalMolecule[SelectedComponent])&&(LastMolecule == SystemComponents.Lambda[SelectedComponent].FractionalMoleculeID))
  {
    //Since the fractional molecule is moved to the place of the selected deleted molecule, update fractional molecule ID on host
    SystemComponents.Lambda[SelectedComponent].FractionalMoleculeID = SystemComponents.TempVal.molecule;
  }
}

//////////////////////////////////////////////////
// PREPARING NEW (TRIAL) LOCATIONS/ORIENTATIONS //
//////////////////////////////////////////////////

__device__ void Rotate_Quaternions(double3 &Vec, double3 RANDOM)
{
  //https://stackoverflow.com/questions/31600717/how-to-generate-a-random-quaternion-quickly
  //https://stackoverflow.com/questions/38978441/creating-uniform-random-quaternion-and-multiplication-of-two-quaternions
  //James J. Kuffner 2004//
  //Zhao's note: this one needs three random numbers//
  const double u = RANDOM.x;
  const double v = RANDOM.y;
  const double w = RANDOM.z;
  const double pi=3.14159265358979323846;//Zhao's note: consider using M_PI
  //Sadly, no double4 type available//
  const double q0 = sqrt(1-u) * std::sin(2*pi*v);
  const double q1 = sqrt(1-u) * std::cos(2*pi*v);
  const double q2 = sqrt(u)   * std::sin(2*pi*w);
  const double q3 = sqrt(u)   * std::cos(2*pi*w);

  double rot[3*3];
  const double a01=q0*q1; const double a02=q0*q2; const double a03=q0*q3;
  const double a11=q1*q1; const double a12=q1*q2; const double a13=q1*q3;
  const double a22=q2*q2; const double a23=q2*q3; const double a33=q3*q3;

  rot[0]=1.0-2.0*(a22+a33);
  rot[1]=2.0*(a12-a03);
  rot[2]=2.0*(a13+a02);
  rot[3]=2.0*(a12+a03);
  rot[4]=1.0-2.0*(a11+a33);
  rot[5]=2.0*(a23-a01);
  rot[6]=2.0*(a13-a02);
  rot[7]=2.0*(a23+a01);
  rot[8]=1.0-2.0*(a11+a22);
  const double r=Vec.x*rot[0*3+0] + Vec.y*rot[0*3+1] + Vec.z*rot[0*3+2];
  const double s=Vec.x*rot[1*3+0] + Vec.y*rot[1*3+1] + Vec.z*rot[1*3+2];
  const double c=Vec.x*rot[2*3+0] + Vec.y*rot[2*3+1] + Vec.z*rot[2*3+2];
  Vec={r, s, c};
}

__device__ void RotationAroundAxis(double3* pos, size_t i, double theta, double3 Axis)
{
  double w,s,c,rot[3*3];

  c=cos(theta);
  w=1.0-c;
  s=sin(theta);

  rot[0*3+0]=(Axis.x)*(Axis.x)*w+c;
  rot[0*3+1]=(Axis.x)*(Axis.y)*w+(Axis.z)*s;
  rot[0*3+2]=(Axis.x)*(Axis.z)*w-(Axis.y)*s;
  rot[1*3+0]=(Axis.x)*(Axis.y)*w-(Axis.z)*s;
  rot[1*3+1]=(Axis.y)*(Axis.y)*w+c;
  rot[1*3+2]=(Axis.y)*(Axis.z)*w+(Axis.x)*s;
  rot[2*3+0]=(Axis.x)*(Axis.z)*w+(Axis.y)*s;
  rot[2*3+1]=(Axis.y)*(Axis.z)*w-(Axis.x)*s;
  rot[2*3+2]=(Axis.z)*(Axis.z)*w+c;

  w=pos[i].x*rot[0*3+0]+pos[i].y*rot[0*3+1]+pos[i].z*rot[0*3+2];
  s=pos[i].x*rot[1*3+0]+pos[i].y*rot[1*3+1]+pos[i].z*rot[1*3+2];
  c=pos[i].x*rot[2*3+0]+pos[i].y*rot[2*3+1]+pos[i].z*rot[2*3+2];
  pos[i].x=w;
  pos[i].y=s;
  pos[i].z=c;
}

__global__ void get_new_position(Simulations& Sim, ForceField FF, size_t start_position, size_t SelectedComponent, double3 MaxChange, double3* RANDOM, size_t index, int MoveType)
{
  const size_t i = blockIdx.x * blockDim.x + threadIdx.x;
  const size_t real_pos = start_position + i;

  const double3 pos       = Sim.d_a[SelectedComponent].pos[real_pos];
  const double  scale     = Sim.d_a[SelectedComponent].scale[real_pos];
  const double  charge    = Sim.d_a[SelectedComponent].charge[real_pos];
  const double  scaleCoul = Sim.d_a[SelectedComponent].scaleCoul[real_pos];
  const size_t  Type      = Sim.d_a[SelectedComponent].Type[real_pos];
  const size_t  MolID     = Sim.d_a[SelectedComponent].MolID[real_pos];

  switch(MoveType)
  {
    case TRANSLATION://TRANSLATION:
    {
      Sim.New.pos[i] = pos + MaxChange * 2.0 * (RANDOM[index] - 0.5);
      Sim.New.scale[i] = scale; 
      Sim.New.charge[i] = charge; 
      Sim.New.scaleCoul[i] = scaleCoul; 
      Sim.New.Type[i] = Type; 
      Sim.New.MolID[i] = MolID;

      Sim.Old.pos[i] = pos;
      Sim.Old.scale[i] = scale; 
      Sim.Old.charge[i] = charge; 
      Sim.Old.scaleCoul[i] = scaleCoul; 
      Sim.Old.Type[i] = Type; 
      Sim.Old.MolID[i] = MolID;
      break;
    }
    case ROTATION://ROTATION:
    {
      Sim.New.pos[i] = pos - Sim.d_a[SelectedComponent].pos[start_position];
      const double3 Angle = MaxChange * 2.0 * (RANDOM[index] - 0.5);
      
      RotationAroundAxis(Sim.New.pos, i, Angle.x, {1.0, 0.0, 0.0}); //X
      RotationAroundAxis(Sim.New.pos, i, Angle.y, {0.0, 1.0, 0.0}); //Y
      RotationAroundAxis(Sim.New.pos, i, Angle.z, {0.0, 0.0, 1.0}); //Z
      
      Sim.New.pos[i] += Sim.d_a[SelectedComponent].pos[start_position];

      Sim.New.scale[i] = scale; 
      Sim.New.charge[i] = charge; 
      Sim.New.scaleCoul[i] = scaleCoul; 
      Sim.New.Type[i] = Type; 
      Sim.New.MolID[i] = MolID;

      Sim.Old.pos[i] = pos;
      Sim.Old.scale[i] = scale; 
      Sim.Old.charge[i] = charge; 
      Sim.Old.scaleCoul[i] = scaleCoul; 
      Sim.Old.Type[i] = Type; 
      Sim.Old.MolID[i] = MolID;
      break;
    }
    case SINGLE_INSERTION:
    { 
      //First ROTATION using QUATERNIONS//
      //Then TRANSLATION//
      double3 BoxLength = {Sim.Box.Cell[0], Sim.Box.Cell[4], Sim.Box.Cell[8]};
      double3 NEW_COM   = BoxLength * RANDOM[index];
      if(i == 0) Sim.New.pos[0] = NEW_COM;
      if(i > 0)
      {
        double3 Vec = pos - Sim.d_a[SelectedComponent].pos[start_position];
        Rotate_Quaternions(Vec, RANDOM[index + 1]);
        Sim.New.pos[i] = Vec + NEW_COM;
      }
      Sim.New.scale[i] = scale;
      Sim.New.charge[i] = charge;
      Sim.New.scaleCoul[i] = scaleCoul;
      Sim.New.Type[i] = Type;
      Sim.New.MolID[i] = Sim.d_a[SelectedComponent].size / Sim.d_a[SelectedComponent].Molsize;
      break;
    }
    case SINGLE_DELETION: //Just Copy the old positions//
    {
      Sim.Old.pos[i] = pos;
      Sim.Old.scale[i] = scale; 
      Sim.Old.charge[i] = charge; 
      Sim.Old.scaleCoul[i] = scaleCoul; 
      Sim.Old.Type[i] = Type; 
      Sim.Old.MolID[i] = MolID;
    }
    
    case SPECIAL_ROTATION:
    {
      // Rotation around axis (first/second atom) in the molecule //
      Sim.New.pos[i] = pos - Sim.d_a[SelectedComponent].pos[start_position];
      
      const double3 Angle = MaxChange * 2.0 * (RANDOM[index] - 0.5);
      //Take direction, normalize//
      double3 Axis;
      //if(i == 0)
      //{ 
        Axis = Sim.d_a[SelectedComponent].pos[start_position + 1] - Sim.d_a[SelectedComponent].pos[start_position];
        double norm = sqrt(dot(Axis, Axis));
        Axis *= 1.0/norm;
        //printf("Atom %lu, DistVec: %.5f %.5f %.5f\n", i, Sim.New.pos[i].x, Sim.New.pos[i].y, Sim.New.pos[i].z);
      //}
      if(i != 0 && i != 1) //Do not rotate the end points of the vector
        RotationAroundAxis(Sim.New.pos, i, 3.0 * Angle.x, Axis);
      Sim.New.pos[i] += Sim.d_a[SelectedComponent].pos[start_position];

      Sim.New.scale[i] = scale;
      Sim.New.charge[i] = charge;
      Sim.New.scaleCoul[i] = scaleCoul;
      Sim.New.Type[i] = Type;
      Sim.New.MolID[i] = MolID;

      Sim.Old.pos[i] = pos;
      Sim.Old.scale[i] = scale;
      Sim.Old.charge[i] = charge;
      Sim.Old.scaleCoul[i] = scaleCoul;
      Sim.Old.Type[i] = Type;
      Sim.Old.MolID[i] = MolID;
      //printf("Atom %lu, Angle: %.5f, Axis: %.5f %.5f %.5f, xyz: %.5f %.5f %.5f\n", i, Angle.x, Axis.x, Axis.y, Axis.z, Sim.New.pos[i].x, Sim.New.pos[i].y, Sim.New.pos[i].z);
    }
  }
  Sim.device_flag[i] = false;
}

////////////////////////////////////////////////////////////////////////////////////
// OPTIMIZING THE ACCEPTANCE (MAINTAIN AROUND 0.5) FOR TRANSLATION/ROTATION MOVES //
////////////////////////////////////////////////////////////////////////////////////

static inline void Update_Max_Translation(Components& SystemComponents, size_t Comp)
{
  if(SystemComponents.Moves[Comp].TranslationTotal == 0) return;
  SystemComponents.Moves[Comp].TranslationAccRatio = static_cast<double>(SystemComponents.Moves[Comp].TranslationAccepted)/SystemComponents.Moves[Comp].TranslationTotal;
  //printf("AccRatio is %.10f\n", MoveStats.TranslationAccRatio);
  if(SystemComponents.Moves[Comp].TranslationAccRatio > 0.5)
  {
    SystemComponents.MaxTranslation[Comp] *= 1.05;
  }
  else
  {
    SystemComponents.MaxTranslation[Comp] *= 0.95;
  }
  if(SystemComponents.MaxTranslation[Comp].x < 0.01) SystemComponents.MaxTranslation[Comp].x = 0.01;
  if(SystemComponents.MaxTranslation[Comp].y < 0.01) SystemComponents.MaxTranslation[Comp].y = 0.01;
  if(SystemComponents.MaxTranslation[Comp].z < 0.01) SystemComponents.MaxTranslation[Comp].z = 0.01;

  if(SystemComponents.MaxTranslation[Comp].x > 5.0) SystemComponents.MaxTranslation[Comp].x = 5.0;
  if(SystemComponents.MaxTranslation[Comp].y > 5.0) SystemComponents.MaxTranslation[Comp].y = 5.0;
  if(SystemComponents.MaxTranslation[Comp].z > 5.0) SystemComponents.MaxTranslation[Comp].z = 5.0;


  SystemComponents.Moves[Comp].CumTranslationAccepted += SystemComponents.Moves[Comp].TranslationAccepted;
  SystemComponents.Moves[Comp].CumTranslationTotal    += SystemComponents.Moves[Comp].TranslationTotal;
  SystemComponents.Moves[Comp].TranslationAccepted = 0;
  SystemComponents.Moves[Comp].TranslationTotal = 0;
}

static inline void Update_Max_Rotation(Components& SystemComponents, size_t Comp)
{
  if(SystemComponents.Moves[Comp].RotationTotal == 0) return;
  SystemComponents.Moves[Comp].RotationAccRatio = static_cast<double>(SystemComponents.Moves[Comp].RotationAccepted)/SystemComponents.Moves[Comp].RotationTotal;
  //printf("AccRatio is %.10f\n", MoveStats.RotationAccRatio);
  if(SystemComponents.Moves[Comp].RotationAccRatio > 0.5)
  {
    SystemComponents.MaxRotation[Comp] *= 1.05;
  }
  else
  {
    SystemComponents.MaxRotation[Comp] *= 0.95;
  }
  if(SystemComponents.MaxRotation[Comp].x < 0.01) SystemComponents.MaxRotation[Comp].x = 0.01;
  if(SystemComponents.MaxRotation[Comp].y < 0.01) SystemComponents.MaxRotation[Comp].y = 0.01;
  if(SystemComponents.MaxRotation[Comp].z < 0.01) SystemComponents.MaxRotation[Comp].z = 0.01;

  if(SystemComponents.MaxRotation[Comp].x > 3.14) SystemComponents.MaxRotation[Comp].x = 3.14;
  if(SystemComponents.MaxRotation[Comp].y > 3.14) SystemComponents.MaxRotation[Comp].y = 3.14;
  if(SystemComponents.MaxRotation[Comp].z > 3.14) SystemComponents.MaxRotation[Comp].z = 3.14;

  SystemComponents.Moves[Comp].CumRotationAccepted += SystemComponents.Moves[Comp].RotationAccepted;
  SystemComponents.Moves[Comp].CumRotationTotal    += SystemComponents.Moves[Comp].RotationTotal;
  SystemComponents.Moves[Comp].RotationAccepted = 0;
  SystemComponents.Moves[Comp].RotationTotal = 0;
}

static inline void Update_Max_SpecialRotation(Components& SystemComponents, size_t Comp)
{
  if(SystemComponents.Moves[Comp].SpecialRotationTotal == 0) return;
  SystemComponents.Moves[Comp].SpecialRotationAccRatio = static_cast<double>(SystemComponents.Moves[Comp].SpecialRotationAccepted)/SystemComponents.Moves[Comp].SpecialRotationTotal;
  //printf("AccRatio is %.10f\n", MoveStats.RotationAccRatio);
  if(SystemComponents.Moves[Comp].SpecialRotationAccRatio > 0.5)
  {
    SystemComponents.MaxSpecialRotation[Comp] *= 1.05;
  }
  else
  {
    SystemComponents.MaxSpecialRotation[Comp] *= 0.95;
  }
  if(SystemComponents.MaxSpecialRotation[Comp].x < 0.01) SystemComponents.MaxSpecialRotation[Comp].x = 0.01;
  if(SystemComponents.MaxSpecialRotation[Comp].y < 0.01) SystemComponents.MaxSpecialRotation[Comp].y = 0.01;
  if(SystemComponents.MaxSpecialRotation[Comp].z < 0.01) SystemComponents.MaxSpecialRotation[Comp].z = 0.01;

  if(SystemComponents.MaxSpecialRotation[Comp].x > 3.14) SystemComponents.MaxSpecialRotation[Comp].x = 3.14;
  if(SystemComponents.MaxSpecialRotation[Comp].y > 3.14) SystemComponents.MaxSpecialRotation[Comp].y = 3.14;
  if(SystemComponents.MaxSpecialRotation[Comp].z > 3.14) SystemComponents.MaxSpecialRotation[Comp].z = 3.14;
  SystemComponents.Moves[Comp].SpecialRotationAccepted = 0;
  SystemComponents.Moves[Comp].SpecialRotationTotal = 0;
}

static inline void Update_Max_VolumeChange(Components& SystemComponents)
{
  if(SystemComponents.VolumeMoveAttempts == 0) return;
  double AccRatio = static_cast<double>(SystemComponents.VolumeMoveAccepted) / static_cast<double>(SystemComponents.VolumeMoveAttempts);
  
  double compare_to_target_ratio = AccRatio / SystemComponents.VolumeMoveTargetAccRatio;
  if(compare_to_target_ratio>1.5) compare_to_target_ratio = 1.5;
  else if(compare_to_target_ratio<0.5) compare_to_target_ratio = 0.5;
  SystemComponents.VolumeMoveMaxChange *= compare_to_target_ratio;

  if(SystemComponents.VolumeMoveMaxChange < 0.0005)
    SystemComponents.VolumeMoveMaxChange = 0.0005;
    if(SystemComponents.VolumeMoveMaxChange > 0.5)
      SystemComponents.VolumeMoveMaxChange=0.5;

  printf("CYCLE: %zu, AccRatio: %.5f, compare_to_target_ratio: %.5f, MaxVolumeChange: %.5f\n", SystemComponents.CURRENTCYCLE, AccRatio, compare_to_target_ratio, SystemComponents.VolumeMoveMaxChange);
 
  SystemComponents.VolumeMoveTotalAccepted += SystemComponents.VolumeMoveAccepted;
  SystemComponents.VolumeMoveTotalAttempts += SystemComponents.VolumeMoveAttempts;
  SystemComponents.VolumeMoveAccepted = 0;
  SystemComponents.VolumeMoveAttempts = 0;
}
