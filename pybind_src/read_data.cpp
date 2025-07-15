#include <filesystem>
#include <stdio.h>
#include <fstream>
#include <string>
#include <vector>
#include <numeric>

#include <algorithm> //for remove_if

#include <iostream>

#include <cfloat> //for DBL_MIN

//#include <print>

//#include "data_struct.h"
#include "VDW_Coulomb.cuh"
//#include "convert_array.h"
#include "read_data.h"

//#include <torch/script.h> // One-stop header.

#define MAX2(x,y) (((x)>(y))?(x):(y))                 // the maximum of two numbers
#define MAX3(x,y,z) MAX2((x),MAX2((y),(z)))

inline std::vector<std::string> split(const std::string txt, char ch)
{
    size_t pos = txt.find(ch);
    size_t initialPos = 0;
    std::vector<std::string> strs{};

    // Decompose statement
    while (pos != std::string::npos) {

        std::string s = txt.substr(initialPos, pos - initialPos);
        if (!s.empty())
        {
            strs.push_back(s);
        }
        initialPos = pos + 1;

        pos = txt.find(ch, initialPos);
    }

    // Add the last one
    std::string s = txt.substr(initialPos, std::min(pos, txt.size()) - initialPos + 1);
    if (!s.empty())
    {
        strs.push_back(s);
    }

    return strs;
}

bool isFloat(const std::string& str) 
{
  std::istringstream iss(str);
  float f;
  iss >> f;
  return iss.eof() && !iss.fail();
}

std::string replaceTabsWithSpaces(const std::string& input, size_t spacesPerTab)
{
  std::string result;

  for (char c : input)
  {
    if (c == '\t')
    {
      // Insert spaces instead of tab
      result += std::string(spacesPerTab, ' ');
    }
    else
    {
      // Keep other characters as they are
      result += c;
    }
  }
  return result;
}

//Zhao's note: DO NOT MIX TAB AND SPACE!!!
//082823 Update, now you can mix tab and space..//
void Split_Tab_Space(std::vector<std::string>& termsScannedLined, std::string& str)
{
  if (str.find("\t", 0) != std::string::npos) //if the delimiter is tab
  {
    //termsScannedLined = split(str, '\t');
    size_t spacesPerTab = 4;
    str = replaceTabsWithSpaces(str, spacesPerTab);
  }
  termsScannedLined = split(str, ' ');
}

void FindIfInputIsThere(std::string& InputCommand, std::string& exepath)
{
  bool exist = false;
  std::vector<std::string> termsScannedLined{};
  std::string str;
  std::ifstream file(exepath + "/" + "read_data.cpp");
  while (std::getline(file, str))
  {
    if(str.find(InputCommand) != std::string::npos)
      exist = true;
  }
  if(!exist)
  {
    //printf("%s Input Command Not Found!!!!\n", InputCommand.c_str());
    throw std::runtime_error("Program Abort due to Unfound Input Command!");
  }
  file.close();
}

void Check_Inputs_In_read_data_cpp(std::string& exepath)
{
  //Check for input, if required input keywords are not there, throw runtime error//
  std::vector<std::string> RequiredInputs_Cycles = {"NumberOfInitializationCycles", "NumberOfEquilibrationCycles", "NumberOfProductionCycles", "UseMaxStep"};
  std::vector<std::string> RequiredInputs_System = {"RestartFile", "RandomSeed", "AdsorbateAllocateSpace", "NumberOfSimulations", "SingleSimulation", "Temperature", "Pressure"};
  std::vector<std::string> RequiredInputs_CBMC   = {"NumberOfTrialPositions", "NumberOfTrialOrientations"};

  std::vector<std::string> RequiredInputs_Framework = {"InputFileType", "FrameworkName", "UnitCells"};

  //Zhao's note!
  //Although there are cases where charge is not needed, for example, "ChargeMethod None", here we still require "UseChargesFromCIFFile"
  //A better solution might be to check the input keyword for "ChargeMethod", if "None", then don't check for "UseChargesFromCIFFile"
  //Need more careful design of input..
  std::vector<std::string> RequiredInputs_ForceField= {"ChargeMethod", "OverlapCriteria", "CutOffVDW", "CutOffCoulomb", "EwaldPrecision", "UseChargesFromCIFFile"};

  std::vector<std::string> RequiredInputs_Adsorbate = {"Component", "MoleculeName", "IdealGasRosenbluthWeight", "FugacityCoefficient", "CreateNumberOfMolecules"};

  std::vector<std::string> RequiredInputs;
  RequiredInputs.insert(RequiredInputs.end(), RequiredInputs_Cycles.begin(), RequiredInputs_Cycles.end());
  RequiredInputs.insert(RequiredInputs.end(), RequiredInputs_System.begin(), RequiredInputs_System.end());
  RequiredInputs.insert(RequiredInputs.end(), RequiredInputs_CBMC.begin(), RequiredInputs_CBMC.end());
  RequiredInputs.insert(RequiredInputs.end(), RequiredInputs_Framework.begin(), RequiredInputs_Framework.end());
  RequiredInputs.insert(RequiredInputs.end(), RequiredInputs_ForceField.begin(), RequiredInputs_ForceField.end());
  RequiredInputs.insert(RequiredInputs.end(), RequiredInputs_Adsorbate.begin(), RequiredInputs_Adsorbate.end());

  std::vector<int> RequiredInput_Index(RequiredInputs.size(), -1);

  //printf("Checking if all inputs are defined\n");
  std::vector<std::string> termsScannedLined{};
  //Zhao's note: Hard-coded executable name//
  termsScannedLined = split(exepath, '/');
  exepath = "/";
  for(size_t i = 0 ; i < termsScannedLined.size() - 1; i++) exepath = exepath + termsScannedLined[i] + "/";
  //printf("True path of exe is %s\n", exepath.c_str());
  std::string str;
  std::ifstream file("simulation.input");
  size_t counter = 0;
  while (std::getline(file, str))
  {
    Split_Tab_Space(termsScannedLined, str);
    if(termsScannedLined.size() == 0) continue;
    std::string InputCommand = termsScannedLined[0];
    FindIfInputIsThere(InputCommand, exepath);

    for(size_t i = 0; i < RequiredInputs.size(); i++)
    {
      if(str.find(RequiredInputs[i], 0) != std::string::npos)
      {
        RequiredInput_Index[i] = counter;
      }
    }
    counter++;
  }
  for(size_t i = 0; i < RequiredInputs.size(); i++)
  {
    if(RequiredInput_Index[i] == -1)
    {
      //printf("Cannot find *** %s *** keyword in simulation.input file!\n", RequiredInputs[i].c_str());
      throw std::runtime_error("Required Keyword Not Found! Abort!!!!");
    }
  }
  file.close();
}

bool caseInSensStringCompare(const std::string& str1, const std::string& str2)
{
    return str1.size() == str2.size() && std::equal(str1.begin(), str1.end(), str2.begin(), [](auto a, auto b) {return std::tolower(a) == std::tolower(b); });
}

void Check_Component_size(Components& SystemComponents)
{
  size_t referenceVal = SystemComponents.MoleculeName.size();
  ////printf("reference size: %zu\n", referenceVal);
  if(SystemComponents.Moleculesize.size() != referenceVal)                   throw std::runtime_error("Moleculesize does not match! reference: " + std::to_string(referenceVal) + " Vector: " + std::to_string(SystemComponents.Moleculesize.size()));
  if(SystemComponents.NumberOfMolecule_for_Component.size() != referenceVal) throw std::runtime_error("NumberOfMolecule_for_Component does not match!");
  if(SystemComponents.MolFraction.size() != referenceVal)                    throw std::runtime_error("MolFraction does not match!");
  if(SystemComponents.IdealRosenbluthWeight.size() != referenceVal)          throw std::runtime_error("IdealRosenbluthWeight does not match!");
  if(SystemComponents.FugacityCoeff.size() != referenceVal)                  throw std::runtime_error("FugacityCoeff does not match!");
  if(SystemComponents.Tc.size() != referenceVal)                             throw std::runtime_error("Tc does not match!");
  if(SystemComponents.Pc.size() != referenceVal)                             throw std::runtime_error("Pc does not match!");
  if(SystemComponents.Accentric.size() != referenceVal)                      throw std::runtime_error("Accentric does not match!");
  if(SystemComponents.rigid.size() != referenceVal)                          throw std::runtime_error("Rigidity (boolean vector) not match!");
  if(SystemComponents.hasfractionalMolecule.size() != referenceVal)          throw std::runtime_error("HasFractionalMolecule (boolean vector) not match!");
  if(SystemComponents.Lambda.size() != referenceVal)                         throw std::runtime_error("Lambda (fractional component vector) not match!");
  if(SystemComponents.Tmmc.size() != referenceVal)                         throw std::runtime_error("Tmmc (TMMC vector) not match!");
  ////printf("CreateMolecule size: %zu\n", SystemComponents.NumberOfCreateMolecules.size());
  if(SystemComponents.NumberOfCreateMolecules.size() != referenceVal)        throw std::runtime_error("Molecules need to create not match!");
}

void read_number_of_sims_from_input(size_t *NumSims, bool *SingleSim)
{
  std::vector<std::string> termsScannedLined{};
  std::string str;
  std::ifstream file("simulation.input");
  size_t tempnum = 0; bool tempsingle = false; size_t counter = 0;
  while (std::getline(file, str))
  {
    counter++;
    if (str.find("SingleSimulation", 0) != std::string::npos)
    {
      Split_Tab_Space(termsScannedLined, str);
      if(caseInSensStringCompare(termsScannedLined[1], "yes"))
      {
        tempsingle = true;
        //printf("running only one simulation\n");
      }
    }
    if (str.find("NumberOfSimulations", 0) != std::string::npos)
    {
      Split_Tab_Space(termsScannedLined, str);
      sscanf(termsScannedLined[1].c_str(), "%zu", &tempnum);
      //printf("There are %zu simulations\n", tempnum);
    }
  }
  *NumSims = tempnum; *SingleSim = tempsingle;
  file.close();
}

void read_simulation_input(Variables& Vars, bool *ReadRestart, bool *SameFrameworkEverySimulation)
{
  bool tempRestart = false;  //Whether we read restart file or not

  std::vector<std::string> termsScannedLined{};
  std::string str;
  std::ifstream file("simulation.input");
  int counter=0;

  size_t tempNComp = 0; 
  int3& NComponents = Vars.TempComponents.NComponents;

  bool tempSameFrameworkEverySimulation = true;
  bool tempSeparateFramework = false;

  while (std::getline(file, str))
  {
    counter++;
    if (str.find("UseGPUReduction", 0) != std::string::npos)
    {
      Split_Tab_Space(termsScannedLined, str);
      if(caseInSensStringCompare(termsScannedLined[1], "yes"))
      {
        Vars.TempWidom.UseGPUReduction = true;
      }
    }
    if (str.find("Useflag", 0) != std::string::npos)
    {
      Split_Tab_Space(termsScannedLined, str);
      if(caseInSensStringCompare(termsScannedLined[1], "yes"))
      {
        Vars.TempWidom.Useflag = true;
      }
    }
  
    if (str.find("RandomSeed", 0) != std::string::npos)
    {
      Split_Tab_Space(termsScannedLined, str);
      Vars.Random.RANDOMSEED = std::stoi(termsScannedLined[1]);
      //printf("Random Seed is %d\n", Vars.Random.RANDOMSEED);
    }

    if (str.find("AdsorbateAllocateSpace", 0) != std::string::npos)
    {
      Split_Tab_Space(termsScannedLined, str);
      sscanf(termsScannedLined[1].c_str(), "%zu", &Vars.Allocate_space_Adsorbate);
      //printf("Allocate space for adsorbate is %zu\n", Vars.Allocate_space_Adsorbate);
    }
    if (str.find("NumberOfInitializationCycles", 0) != std::string::npos)
    {
      Split_Tab_Space(termsScannedLined, str);
      sscanf(termsScannedLined[1].c_str(), "%d", &Vars.NumberOfInitializationCycles);
      ////printf("line is %u, there are %zu Framework Atoms\n", counter, NumberFrameworkAtom);
    }
    if (str.find("NumberOfEquilibrationCycles", 0) != std::string::npos)
    {
      Split_Tab_Space(termsScannedLined, str);
      sscanf(termsScannedLined[1].c_str(), "%d", &Vars.NumberOfEquilibrationCycles);
      ////printf("line is %u, there are %zu Framework Atoms\n", counter, NumberFrameworkAtom);
    }
    if (str.find("NumberOfProductionCycles", 0) != std::string::npos)
    {
      Split_Tab_Space(termsScannedLined, str);
      sscanf(termsScannedLined[1].c_str(), "%d", &Vars.NumberOfProductionCycles);
      ////printf("line is %u, there are %zu Framework Atoms\n", counter, NumberFrameworkAtom);
    }
    if (str.find("NumberOfTrialPositions", 0) != std::string::npos)
    {
      Split_Tab_Space(termsScannedLined, str);
      sscanf(termsScannedLined[1].c_str(), "%zu", &Vars.TempWidom.NumberWidomTrials);
      ////printf("line is %u, there are %zu Framework Atoms\n", counter, NumberFrameworkAtom);
    }
    if (str.find("NumberOfTrialOrientations", 0) != std::string::npos)
    {
      Split_Tab_Space(termsScannedLined, str);
      sscanf(termsScannedLined[1].c_str(), "%zu", &Vars.TempWidom.NumberWidomTrialsOrientations);
      ////printf("line is %u, there are %zu Framework Atoms\n", counter, NumberFrameworkAtom);
    }
    //Zhao's note: Move it somewhere else//
    if (str.find("NumberOfBlocks", 0) != std::string::npos)
    {
      Split_Tab_Space(termsScannedLined, str);
      ////printf("line is %u, there are %zu Framework Atoms\n", counter, NumberFrameworkAtom);
    }
    if (str.find("Pressure", 0) != std::string::npos)
    {
      Split_Tab_Space(termsScannedLined, str);
      Vars.TempComponents.Pressure=std::stod(termsScannedLined[1]);
    }
    if (str.find("Temperature", 0) != std::string::npos)
    {
      Split_Tab_Space(termsScannedLined, str);
      Vars.TempComponents.Temperature=std::stod(termsScannedLined[1]);
    }
    if (str.find("RestartFile", 0) != std::string::npos)
    {
      Split_Tab_Space(termsScannedLined, str);
      if(caseInSensStringCompare(termsScannedLined[1], "yes"))
      {
        tempRestart = true;
        //printf("USE CONFIGURATION FROM RESTARTINITIAL FILE\n");
      }
    }
    if (str.find("DifferentFrameworks", 0) != std::string::npos)
    {
      Split_Tab_Space(termsScannedLined, str);
      if(caseInSensStringCompare(termsScannedLined[1], "yes"))
      {
        tempSameFrameworkEverySimulation = false;
      }
    }
    //Check number of adsorbates to process, need to put them in serial order//
    //Component 0, Component 1, ...//
    if (str.find("Component " + std::to_string(tempNComp), 0) != std::string::npos)
    {
      tempNComp ++;
    }
    //Check if we need to separate framework species
    //species could be cations, linkers, attachments to the nodes, ...//
    if (str.find("SeparateFrameworkComponents", 0) != std::string::npos)
    {
      Split_Tab_Space(termsScannedLined, str);
      if(caseInSensStringCompare(termsScannedLined[1], "yes"))
      {
        tempSeparateFramework = true;
        //printf("NEED TO SEPARATE FRAMEWORK COMPONENTS\n");
      }
    }
    //If we need to separate framework species, then read the following//
    //This requires us to put the "FrameworkComponents XXX" command after "SeparateFrameworkComponents" command//
    //Zhao's note: make sure there are two spaces before the actual command//
    if(tempSeparateFramework)
    {
      if(str.find("NumberofFrameworkComponents", 0) != std::string::npos)
      {
        Split_Tab_Space(termsScannedLined, str);
        NComponents.y = std::stoi(termsScannedLined[1]);
        //printf("THERE ARE %zu SEPARATE FRAMEWORK COMPONENTS\n", NComponents.y);
      }
    }
    if(counter>200) break;
  }
  *ReadRestart   = tempRestart;
  *SameFrameworkEverySimulation = tempSameFrameworkEverySimulation;
  NComponents.z = tempNComp; //z component is the adsorbate components//
  NComponents.x = tempNComp + NComponents.y;
  ////printf("Finished Checking Number of Components, There are %zu framework, %zu Adsorbates, %zu total Components\n", NComponents.y, NComponents.z, NComponents.x);
  file.close();
}

void read_Gibbs_and_Cycle_Stats(Variables& Vars, bool& SetMaxStep, size_t& MaxStepPerCycle)
{
  size_t counter = 0;
  double temp = 0.0;
  std::vector<std::string> termsScannedLined{};
  std::string str;
  std::ifstream file("simulation.input");

  while (std::getline(file, str))
  {
    counter++;
    if (str.find("GibbsVolumeChangeProbability", 0) != std::string::npos)
    {
      Split_Tab_Space(termsScannedLined, str);
      temp=std::stod(termsScannedLined[1]);
      if(temp > 0)
      {
        Vars.GibbsStatistics.DoGibbs = true;
        Vars.GibbsStatistics.GibbsBoxProb = temp;
      }
    }
    if (str.find("NPTVolumeChangeProbability", 0) != std::string::npos)
    {
      Split_Tab_Space(termsScannedLined, str);
      temp=std::stod(termsScannedLined[1]);
      if(temp > 0)
      {
        Vars.TempComponents.PerformVolumeMove = true;
        Vars.TempComponents.VolumeMoveProbability = temp;
      }
    }
    if (str.find("UseMaxStep", 0) != std::string::npos)
    {
      Split_Tab_Space(termsScannedLined, str);
      if(caseInSensStringCompare(termsScannedLined[1], "yes"))
      {
        SetMaxStep = true;
      }
    }
    if (str.find("MaxStepPerCycle", 0) != std::string::npos)
    {
      Split_Tab_Space(termsScannedLined, str);
      sscanf(termsScannedLined[1].c_str(), "%zu", &MaxStepPerCycle);
      if(MaxStepPerCycle == 0) throw std::runtime_error("Max Steps per Cycle must be greater than ZERO!");
    }

    if (str.find("StructureFactor_Multiplier", 0) != std::string::npos)
    {
      Split_Tab_Space(termsScannedLined, str);
      sscanf(termsScannedLined[1].c_str(), "%zu", &Vars.StructureFactor_Multiplier);
    }
    
    if(counter>200) break;
  }
  if(SetMaxStep) 
  {
    printf("Setting Maximum Number of Steps for a Cycle, Max Step = %zu\n", MaxStepPerCycle);
  }
  else
  {
    printf("Running Cycles in the Normal Way\n");
  }
  Vars.GibbsStatistics.GibbsBoxStats  = {0, 0};
  Vars.GibbsStatistics.GibbsXferStats = {0, 0};
  file.close();
}

void read_FFParams_from_input(Input_Container& Input)
{
  std::vector<std::string> termsScannedLined{};
  std::string str;

  double tempOverlap = 1.0e6; double tempvdwcut = 12.0; double tempcoulcut = 12.0;
  double tempprecision = 1.0e-6;
  bool tempnoCharges = true;
  //double tempalpha = 0.26506; //Zhao's note: here we used alpha equal to the preset value from raspa3. Need to revisit RASPA-2 for the exact calculation of alpha.
  //Zhao's note: Using the heuresitic equation for converting Ewald Precision to alpha.

  std::ifstream file("simulation.input");
  while (std::getline(file, str))
  {
    if (str.find("ChargeMethod", 0) != std::string::npos)
    {
      Split_Tab_Space(termsScannedLined, str);
      if(caseInSensStringCompare(termsScannedLined[1], "Ewald"))
      {
        tempnoCharges = false;
        //printf("USE EWALD SUMMATION FOR CHARGE\n");
      }
    }
    if (str.find("OverlapCriteria", 0) != std::string::npos)
    {
      Split_Tab_Space(termsScannedLined, str);
      tempOverlap = std::stod(termsScannedLined[1]);
    }
    if (str.find("CutOffVDW", 0) != std::string::npos)
    {
      Split_Tab_Space(termsScannedLined, str);
      tempvdwcut = std::stod(termsScannedLined[1]);
    }
    if (str.find("CutOffCoulomb", 0) != std::string::npos)
    {
      Split_Tab_Space(termsScannedLined, str);
      tempcoulcut = std::stod(termsScannedLined[1]);
    }
    if (str.find("EwaldPrecision", 0) != std::string::npos)
    {
      Split_Tab_Space(termsScannedLined, str);
      tempprecision = std::stod(termsScannedLined[1]);
      //tempalpha = (1.35 - 0.15 * log(tempprecision))/tempcoulcut; // Zhao's note: heurestic equation //
    }
    if (str.find("CBMCBiasingMethod", 0) != std::string::npos)
    {
      Split_Tab_Space(termsScannedLined, str);
      if(caseInSensStringCompare(termsScannedLined[1], "LJ_And_Real_Biasing"))
      {
        Input.VDWRealBias = true; //By default, it is using LJ + Real Biasing
      }
      else if(caseInSensStringCompare(termsScannedLined[1], "LJ_Biasing"))
      {
        Input.VDWRealBias = false;
      }
    }
    if (str.find("Component ", 0) != std::string::npos) //When it reads component, skip//
      break;
  }
  //read FF array
  Input.OverlapCriteria   = tempOverlap;
  Input.CutOffVDW         = tempvdwcut*tempvdwcut;
  Input.CutOffCoul        = tempcoulcut*tempcoulcut;
  Input.EwaldPrecision    = tempprecision;
  Input.noCharges         = tempnoCharges;
  file.close();
}

void ReadFrameworkComponentMoves(Move_Statistics& MoveStats, Components& SystemComponents, size_t comp)
{
  if(SystemComponents.NComponents.y <= 1)
  { //printf("Only one Framework Component, No moves assigned\n"); 
    return;
  }
  if(comp >= SystemComponents.NComponents.y) return;
  //printf("Checking Framework Moves for Framework Component %zu\n", comp);
  std::string FrameworkComponentName = "Framework_Component_ " + std::to_string(comp); //Separate with a space, add a _ to avoid confusion with Adsorbate species//

  std::vector<std::string> termsScannedLined{};
  std::string str;

  std::ifstream file("simulation.input");

  size_t start_counter = 0; bool FOUND = false;
  std::string start_string = FrameworkComponentName;
  std::string terminate_string="END_OF_" + FrameworkComponentName;
  //first get the line number of the destinated component
  while (std::getline(file, str))
  {
    if(str.find(start_string, 0) != std::string::npos){FOUND = true; break;}
    start_counter++;
  }

  if(!FOUND)
  {
    //printf("%s not found in simulation.input\n", FrameworkComponentName.c_str()); 
    return;
  }
  
  //printf("%s starts at line number %zu\n", start_string.c_str(), start_counter);

  file.clear();
  file.seekg(0);

  size_t counter = 0;
  while (std::getline(file, str))
  {
    if(str.find(terminate_string, 0) != std::string::npos){break;}
    if(counter >= start_counter) //start reading after touching the starting line number
    if (str.find(FrameworkComponentName, 0) != std::string::npos)
    {
      Split_Tab_Space(termsScannedLined, str);
      //printf("Found Framework_Component %s in simulation.input file\n", std::to_string(comp).c_str());
    }
    if (str.find("TranslationProbability", 0) != std::string::npos)
    {
      Split_Tab_Space(termsScannedLined, str);
      MoveStats.TranslationProb=std::stod(termsScannedLined[1]);
    }
    if (str.find("RotationProbability", 0) != std::string::npos)
    {
      Split_Tab_Space(termsScannedLined, str);
      MoveStats.RotationProb=std::stod(termsScannedLined[1]);
    }
    if (str.find("RotationSpecialProbability", 0) != std::string::npos)
    {
      Split_Tab_Space(termsScannedLined, str);
      MoveStats.SpecialRotationProb=std::stod(termsScannedLined[1]);
      //printf("WARNING: Special Rotations are rotations with pre-set Rotation Axes, Rotation Axes, Angles are needed to define in def files for %s !\n", FrameworkComponentName.c_str());
    }
    //Try to add reinsertion move here//
    if (str.find("ReinsertionProbability", 0) != std::string::npos)
    {
      Split_Tab_Space(termsScannedLined, str);
      MoveStats.ReinsertionProb=std::stod(termsScannedLined[1]);
    }
    if (str.find("END_OF_Framework_Component_" + std::to_string(comp), 0) != std::string::npos)
    {
      //printf("Reach the end of %s\n", FrameworkComponentName.c_str()); 
      break;
    } 
    counter ++;
  }
  MoveStats.NormalizeProbabilities();
  MoveStats.PrintProbabilities();
  file.close();
}

void read_Ewald_Parameters_from_input(double CutOffCoul, Boxsize& Box, double precision)
{
  double tempprefactor = 138935.483496;
  Box.Prefactor = tempprefactor;
  double tempalpha = 0.26506; //Zhao's note: here we used alpha equal to the preset value from raspa3. Need to revisit RASPA-2 for the exact calculation of alpha.
  // Zhao's note: use the Ewald method in LAMMPS, need full control of alpha and kvectors //
  bool Ewald_UseLAMMPS_Setup = false; //Default using RASPA-2's setup (automatic, heurestic)
  //printf("----------------EWALD SUMMATION SETUP-----------------\n");
  //Zhao's note: add controls if the users want to use specified ewald parameters//
  std::vector<std::string> termsScannedLined{};
  std::string str;

  std::ifstream file("simulation.input");
  while (std::getline(file, str))
  {
    if (str.find("Ewald_UseLAMMPS_Setup", 0) != std::string::npos)
    {
      Split_Tab_Space(termsScannedLined, str);
      if(caseInSensStringCompare(termsScannedLined[1], "yes"))
      {
        Ewald_UseLAMMPS_Setup = true;
      }
    }
  }
  file.clear();
  file.seekg(0);
  if(Ewald_UseLAMMPS_Setup) 
  //Then read alpha and number of kvectors in xyz, calculate reciprocal cutoff//
  {
    bool AlphaFound = false;
    bool kvecFound  = false;
    while (std::getline(file, str))
    {
      if (str.find("Ewald_Alpha", 0) != std::string::npos)
      {
        Split_Tab_Space(termsScannedLined, str);
        tempalpha=std::stod(termsScannedLined[1]);
        Box.Alpha  = tempalpha;
        AlphaFound = true;
      }
      if (str.find("Ewald_kvectors", 0) != std::string::npos)
      {
        Split_Tab_Space(termsScannedLined, str);
        //printf("termsScannedLined.size(): %zu\n", termsScannedLined.size());
        if(termsScannedLined.size() != 4)
        {
          throw std::runtime_error("Need 3 values provided for the Ewald kvectors! Check your input\n");
        }
        Box.kmax.x=std::stoi(termsScannedLined[1]);
        Box.kmax.y=std::stoi(termsScannedLined[2]);
        Box.kmax.z=std::stoi(termsScannedLined[3]);
        if(Box.kmax.x == 0 || Box.kmax.y == 0 || Box.kmax.z == 0) throw std::runtime_error("Number of kvectors cannot be zero!");
        kvecFound  = true;
      }
    }
    if(AlphaFound && kvecFound)
    {
      Box.UseLAMMPSEwald = true;
      //printf("Using LAMMPS Setup for Ewald, need to specify both Alpha and number of kvectors in simulation.input file\n");
      //Calculate Reciprocal cutoff
      double lx = Box.Cell[0];
      double ly = Box.Cell[4];
      double lz = Box.Cell[8];
      double xy = 0.0; //Box.Cell[3];
      double xz = 0.0; //Box.Cell[6];
      double yz = 0.0; //Box.Cell[7];
      ////printf("Box Vals: %.5f %.5f %.5f\n %.5f %.5f %.5f\n %.5f %.5f %.5f\n", Box.Cell[0], Box.Cell[1], Box.Cell[2], Box.Cell[3], Box.Cell[4], Box.Cell[5], Box.Cell[6], Box.Cell[7], Box.Cell[8]);
      ////printf("xy: %.5f, xz: %.5f, yz: %.5f\n", xy, xz, yz);
      double ux = 2*M_PI/lx;
      double uy = 2*M_PI*(-xy)/lx/ly;
      double uz = 2*M_PI*(xy*yz - ly*xz)/lx/ly/lz;
      double vy = 2*M_PI/ly;
      double vz = 2*M_PI*(-yz)/ly/lz;
      double wz = 2*M_PI/lz;
      const double kvecx = Box.kmax.x*ux;
      const double kvecy = Box.kmax.x*uy + Box.kmax.y*vy;
      const double kvecz = Box.kmax.x*uz + Box.kmax.y*vz + Box.kmax.z*wz;
      Box.ReciprocalCutOff = MAX3(kvecx*kvecx, kvecy*kvecy, kvecz*kvecz) * 1.00001;
    }
    else throw std::runtime_error("Need to specify both **Ewald_Alpha** and **Ewald_kvectors** in your input file. Check your input!\n");
  }
  else //Doing the RASPA-2 way of automatically calculate alpha and kvecs using Precision
  {
    //printf("Using RASPA-2 way, need to specify the Precision of Ewald\n");
    double tol = sqrt(fabs(log(precision*CutOffCoul)));
    tempalpha  = sqrt(fabs(log(precision*CutOffCoul*tol)))/CutOffCoul;
    double tol1= sqrt(-log(precision*CutOffCoul*pow(2.0*tol*tempalpha, 2)));
    Box.tol1   = tol1;
    Box.Alpha  = tempalpha;
    //Zhao's note: See InitializeEwald function in RASPA-2.0 //
    Box.kmax.x = std::round(0.25 + Box.Cell[0] * tempalpha * tol1/M_PI);
    Box.kmax.y = std::round(0.25 + Box.Cell[4] * tempalpha * tol1/M_PI);
    Box.kmax.z = std::round(0.25 + Box.Cell[8] * tempalpha * tol1/M_PI);
    Box.ReciprocalCutOff = pow(1.05*static_cast<double>(MAX3(Box.kmax.x, Box.kmax.y, Box.kmax.z)), 2);
    //printf("tol: %.5f, tol1: %.5f\n", tol, tol1);
  } 
  //printf("ALpha is %.5f, Prefactor: %.5f\n", Box.Alpha, Box.Prefactor);
  //printf("kmax: %d %d %d, ReciprocalCutOff: %.5f\n", Box.kmax.x, Box.kmax.y, Box.kmax.z, Box.ReciprocalCutOff);
  //printf("------------------------------------------------------\n");
  file.close();
}

inline std::string& tolower(std::string& s)
{
    for (auto& c : s)
    {
        [[maybe_unused]] auto t = std::tolower(static_cast<unsigned char>(c));
    }

    return s;
}

inline std::string trim2(const std::string& s)
{
    auto start = s.begin();
    while (start != s.end() && std::isspace(*start)) {
        start++;
    }

    auto end = s.end();
    do {
        end--;
    } while (std::distance(start, end) > 0 && std::isspace(*end));

    return std::string(start, end + 1);
}

double Get_Shifted_Value(double epsilon, double sigma, double CutOffSquared)
{
  double scaling = 1.0;
  double arg1 = epsilon;
  double arg2 = sigma * sigma;
  double rr = CutOffSquared;
  double temp = (rr / arg2);
  double rri3 = 1.0 / ((temp * temp * temp) + 0.5 * (1.0 - scaling) * (1.0 - scaling));
  double shift = scaling * (4.0 * arg1 * (rri3 * (rri3 - 1.0)));
  return shift;
}

double Mixing_Rule_Epsilon(double ep1, double ep2)
{
  return sqrt(ep1*ep2); //Assuming Lorentz Berthelot
}

double Mixing_rule_Sigma(double sig1, double sig2)
{
  return 0.5*(sig1+sig2); //Assuming Lorentz Berthelot
}
//Zhao's note: read force field definition first, then read pseudo-atoms
//The type numbering is determined in here, then used in pseudo_atoms.def
//void ForceFieldParser(ForceField& FF, PseudoAtomDefinitions& PseudoAtom)
void ForceFieldParser(Input_Container& Input, PseudoAtomDefinitions& PseudoAtom)
{
  Input.AtomFF = std::vector<Atom_FF>();
  std::string scannedLine; std::string str;
  std::vector<std::string> termsScannedLined{};
  size_t counter = 0;
  std::ifstream FFMixfile("force_field_mixing_rules.def");
  size_t NumberOfDefinitions = 0;
  //Temporary vectors for storing the data
  // Some other temporary values
  bool shifted = false; bool tail = false;
  Atom_FF AtomFF; //placeholder for temp values to push back//
  //printf("------------------PARSING FORCE FIELD MIXING RULES----------------\n");
  // First read the pseudo atom file
  while (std::getline(FFMixfile, str))
  {
    if(counter == 1) //read shifted/truncated
    {
      Split_Tab_Space(termsScannedLined, str);
      if(termsScannedLined[0] == "shifted")
        shifted = true;
    }
    else if(counter == 3) //read tail correction
    {
      Split_Tab_Space(termsScannedLined, str);
      if(termsScannedLined[0] == "yes")
      {
        tail = true; //Zhao's note: not implemented
      }
    }
    else if(counter == 5) // read number of force field definitions
    {
      Split_Tab_Space(termsScannedLined, str);
      sscanf(termsScannedLined[0].c_str(), "%zu", &NumberOfDefinitions);
      if (NumberOfDefinitions <= 0) throw std::runtime_error("Incorrect amount of force field definitions");
      //if (NumberOfDefinitions > 100) //printf("There are A LOT OF Definitions. Do you need all of them? Okay...\n");
    }
    else if(counter >= 7) // read data for each force field definition
    {
      //printf("%s\n", str.c_str());
      Split_Tab_Space(termsScannedLined, str);

      PseudoAtom.Name.push_back(termsScannedLined[0]);

      AtomFF.Name    = termsScannedLined[0];
      AtomFF.epsilon = std::stod(termsScannedLined[2]);
      AtomFF.sigma   = std::stod(termsScannedLined[3]);
      AtomFF.shift   = shifted;
      AtomFF.tail    = tail;
      
      Input.AtomFF.push_back(AtomFF);
    }
    counter++;
    if(counter==7+NumberOfDefinitions) break; //in case there are extra empty rows, Zhao's note: I am skipping the mixing rule, assuming Lorentz-Berthelot
  }
  //printf("------------------------------------------------------------------\n");
  FFMixfile.close();
}

static inline double GetTailCorrectionValue(double epsilon, double sigma, double cutoff)
{
  //double scaling = 1.0; Zhao's note: Need more care about fractional molecules with tail corrections
  double arg1 = epsilon;
  double arg2 = sigma * sigma * sigma;
  double rr = sqrt(cutoff);
  double term1= pow(arg2, 4) / (9.0 * pow(rr, 9)); //sigma^12/(9*r^9)
  double term2= pow(arg2, 2) / (3.0 * pow(rr, 3)); //sigma^6 /(3*r^3)
  double val = 16.0 * 3.14159265358979323846 / 2.0 * arg1 * (term1 - term2);
  return val;
}

static inline void PrepareTailCorrection(size_t i, size_t j, size_t NPseudoAtoms, std::vector<Tail>& TempTail, std::vector<double>& Mix_Epsilon, std::vector<double>& Mix_Sigma, double cutoff)
{
  size_t IJ_Forward = i * NPseudoAtoms + j;
  size_t IJ_Reverse = j * NPseudoAtoms + i;
  TempTail[IJ_Forward].UseTail= true;
  TempTail[IJ_Forward].Energy = GetTailCorrectionValue(Mix_Epsilon[IJ_Forward], Mix_Sigma[IJ_Forward], cutoff);
  if(i!=j) TempTail[IJ_Reverse] = TempTail[IJ_Forward];
  //printf("TypeI: %zu, TypeJ: %zu, NPseudoAtoms: %zu, Energy: %.5f\n", i, j, NPseudoAtoms, TempTail[IJ_Forward].Energy);
}
/*
static inline size_t GetTypeForPseudoAtom(PseudoAtomDefinitions& PseudoAtom, std::string& AtomName)
{
  bool Found = false;
  size_t AtomTypeInt = 0;
  for(size_t j = 0; j < PseudoAtom.Name.size(); j++)
  {
    if(AtomName == PseudoAtom.Name[j])
    {
      AtomTypeInt = j;
      Found = true;
      break;
    }
  }
  if(!Found){throw std::runtime_error("Overwriting terms are not Found in Pseudo atoms!!! [" + AtomName + "]");}
  return AtomTypeInt;
}
*/
static inline size_t GetTypeFromFFName(Input_Container& Input, std::string& AtomName)
{
  bool Found = false;
  size_t AtomTypeInt = 0;
  for(size_t j = 0; j < Input.AtomFF.size(); j++)
  {
    if(AtomName == Input.AtomFF[j].Name)
    {
      AtomTypeInt = j;
      Found = true;
      break;
    }
  }
  if(!Found){throw std::runtime_error("Overwriting terms are not Found in Pseudo atoms!!! [" + AtomName + "]");}
  return AtomTypeInt;
}


//Function for Overwritten tail corrections
//For now, it only considers tail correction//
//Add overwritting LJ//
void OverWrite_Mixing_Rule(Input_Container& Input)
{
  std::string scannedLine; std::string str;
  std::vector<std::string> termsScannedLined{};
  std::ifstream FFMixfile("force_field_mixing_rules.def");
  bool shifted = false;
  size_t counter = 0;
  //check if vdw is shifted
  while (std::getline(FFMixfile, str))
  {
    if(counter == 1) //read shifted/truncated
    {
      Split_Tab_Space(termsScannedLined, str);
      if(termsScannedLined[0] == "shifted")
        shifted = true;
      break;
    }
    counter ++;
  }
  FFMixfile.close();

  //Check if file exist//
  std::ifstream OverWritefile("force_field.def");
  std::filesystem::path pathfile = std::filesystem::path("force_field.def");
  if (!std::filesystem::exists(pathfile))
  {
    //printf("Force Field OverWrite file not found\n");
    return;
  }
  
  size_t startline = 0; size_t Noverwrite = 0;
  counter = 0;
  //check # of mixing rules to overwrite
  while (std::getline(OverWritefile, str))
  {
    if (str.find("mixing rules to overwrite", 0) != std::string::npos) //read OverWriteSize
    {
      startline = counter;
    }
    if(startline > 0 && (counter == startline + 1)) //Read next line
    {
      Split_Tab_Space(termsScannedLined, str);
      sscanf(termsScannedLined[0].c_str(), "%zu", &Noverwrite);
      break;
    }
    counter ++;
  }
  //printf("----- OVERWRITTING VDW PARAMETERS -----\n");
  //printf("There are %zu overwritting entries, starting from line %zu\n", Noverwrite, startline);
  if(Noverwrite == 0) 
  {
    return;
  }
  OverWritefile.clear();
  OverWritefile.seekg(0);

  size_t FFsize = Input.AtomFF.size();

  counter = 0;
  while (std::getline(OverWritefile, str))
  {
    //printf("counter: %zu, %s, shifted: %s\n", counter, (counter >= (startline+3) && counter < (3 + startline + Noverwrite)) ? "true" : "false", shifted ? "true" : "false");
    if(counter >= (startline+3) && counter < (3 + startline + Noverwrite))
    {
      //printf("cutting string: %s\n", str.c_str());
      Split_Tab_Space(termsScannedLined, str);
      if(termsScannedLined.size() == 5) //5 entries = LJ mixing rule parameter length
      {
        //Ow Hw lennard-jones 1.0 1.0
        //printf("Overwritting parameters for VDW (lennard-jones): %s and %s\n", termsScannedLined[0].c_str(), termsScannedLined[1].c_str());
        size_t typeI = GetTypeFromFFName(Input, termsScannedLined[0]);
        size_t typeJ = GetTypeFromFFName(Input, termsScannedLined[1]);

        double temp_ep = std::stod(termsScannedLined[3])/1.20272430057; //K -> 10J/mol
        double temp_sig= std::stod(termsScannedLined[4]);               //Angstroem
       
        Input.Mix_Epsilon[typeI * FFsize + typeJ] = temp_ep;
        Input.Mix_Epsilon[typeJ * FFsize + typeI] = temp_ep;
        Input.Mix_Sigma[typeI * FFsize + typeJ] = temp_sig;
        Input.Mix_Sigma[typeJ * FFsize + typeI] = temp_sig;
        if(shifted)
        {
          Input.Mix_Shift[typeI * FFsize + typeJ] = Get_Shifted_Value(temp_ep, temp_sig, Input.CutOffVDW);
          Input.Mix_Shift[typeJ * FFsize + typeI] = Get_Shifted_Value(temp_ep, temp_sig, Input.CutOffVDW);
        }
      }
    }
    counter ++;
  }
  //printf("----- MIXED VDW PARAMETERS (WITH OVERWRITTEN TERMS) [10 J/mol, Angstroem]-----\n");
  for(size_t ii = 0; ii < FFsize; ii++)
    for(size_t jj = 0; jj < FFsize; jj++)
    {
      //if(ii <= jj) //printf("ii: %zu, jj: %zu, Name_i: %s, Name_j: %s, ep: %.10f, sig: %.10f, shift: %.10f\n", ii,jj,Input.AtomFF[ii].Name.c_str(), Input.AtomFF[jj].Name.c_str(), Input.Mix_Epsilon[ii * FFsize + jj], Input.Mix_Sigma[ii * FFsize + jj], Input.Mix_Shift[ii * FFsize + jj]);
    }
  //printf("----- END OF MIXED VDW PARAMETERS (WITH OVERWRITTEN TERMS) -----\n");
  OverWritefile.clear();
  OverWritefile.seekg(0);
  OverWritefile.close();
}

void OverWriteTailCorrection(Input_Container& Input)
{
  std::string scannedLine; std::string str;
  std::vector<std::string> termsScannedLined{};
  size_t counter = 0;
  size_t OverWriteSize = 0;
  size_t typeI; size_t typeJ;
  std::vector<Tail>& TempTail = Input.Mix_Tail;

  std::ifstream OverWritefile("force_field.def");
  std::filesystem::path pathfile = std::filesystem::path("force_field.def");
  if (!std::filesystem::exists(pathfile))
  {
    //printf("Force Field OverWrite file not found\n");
    return;
  }
  //printf("----------------FORCE FIELD OVERWRITTEN (TAIL CORRECTION) PARAMETERS----------------\n");
  while (std::getline(OverWritefile, str))
  { 
    if(counter == 1) //read OverWriteSize
    {
      Split_Tab_Space(termsScannedLined, str);
      sscanf(termsScannedLined[0].c_str(), "%zu", &OverWriteSize);
    }
    else if(counter >= 3 && counter < (3 + OverWriteSize)) //read Terms to OverWrite
    {
      if(str.find("#", 0) != std::string::npos) throw std::runtime_error("Found # in when processing tail correction overwrite terms. Check your number of rules to overwrite!");
      Split_Tab_Space(termsScannedLined, str);
      if(termsScannedLined.size() != 4) throw std::runtime_error("Tail Correction Overwrite need 4 terms for each entry! Check your force_field.def file!");
      if(termsScannedLined[3] == "yes") //Use Tail Correction or not//
      {
        typeI = GetTypeFromFFName(Input, termsScannedLined[0]);
        typeJ = GetTypeFromFFName(Input, termsScannedLined[1]);
        PrepareTailCorrection(typeI, typeJ, Input.AtomFF.size(), TempTail, Input.Mix_Epsilon, Input.Mix_Sigma, Input.CutOffVDW);
      }
    }
    counter ++;    if(counter==3 + OverWriteSize) break; //in case there are extra empty rows, Zhao's note: I am skipping the mixing rule, assuming Lorentz-Berthelot
  }
  //printf("------------------------------------------------------------------------------------\n");
  //Eliminate the terms that do not have tail corrections//
  OverWritefile.close();
}

//LJ, shift, and tail, General mixing rule//
void ForceField_Processing(Input_Container& Input)
{
  size_t NumberOfDefinitions = Input.AtomFF.size();
  //Re-initialize values//
  Input.Mix_Epsilon = std::vector<double>();
  Input.Mix_Sigma   = std::vector<double>();
  Input.Mix_Shift   = std::vector<double>();
  Input.Mix_Tail    = std::vector<Tail>(NumberOfDefinitions * NumberOfDefinitions);

  //Do mixing rule (assuming Lorentz-Berthelot)
  //Declare some temporary arrays
  std::vector<int>Mix_Type; //Force Field Type (Zhao's note: assuming zero, which is Lennard-Jones)
  double temp_ep = 0.0; double temp_sig = 0.0;
  for(size_t i = 0; i < NumberOfDefinitions; i++)
  {
    for(size_t j = 0; j < NumberOfDefinitions; j++)
    {
      double eps_i = Input.AtomFF[i].epsilon;
      double eps_j = Input.AtomFF[j].epsilon;
      double sig_i = Input.AtomFF[i].sigma;
      double sig_j = Input.AtomFF[j].sigma;
      temp_ep = Mixing_Rule_Epsilon(eps_i, eps_j)/1.20272430057; //Zhao's note: caveat here: need to do full energy conversion
      temp_sig= Mixing_rule_Sigma(sig_i, sig_j);
      Input.Mix_Epsilon.push_back(temp_ep);
      Input.Mix_Sigma.push_back(temp_sig);

      if(Input.AtomFF[i].shift && Input.AtomFF[j].shift)
      {
        Input.Mix_Shift.push_back(Get_Shifted_Value(temp_ep, temp_sig, Input.CutOffVDW));
      }
      else
      {
        Input.Mix_Shift.push_back(0.0);
      }
      //If atom i needs tail correction OR atom j needs tail correction OR specifically interaction IJ needs it//
      if(Input.AtomFF[i].tail && Input.AtomFF[j].tail)
      {
        PrepareTailCorrection(i, j, Input.AtomFF.size(), Input.Mix_Tail, Input.Mix_Epsilon, Input.Mix_Sigma, Input.CutOffVDW);
      }
      Input.Mix_Z.push_back(0.0);
      Input.Mix_Type.push_back(0);
    }
  }
  //For checking if mixing rule terms are correct//
  //printf("----- MIXED VDW PARAMETERS -----\n");
  //for(size_t i = 0; i < Input.Mix_Shift.size(); i++)
  //{
    //size_t ii = i/NumberOfDefinitions; size_t jj = i%NumberOfDefinitions;
    //if(ii <= jj) //printf("i: %zu, ii: %zu, jj: %zu, Name_i: %s, Name_j: %s, ep: %.10f, sig: %.10f, shift: %.10f\n", i,ii,jj,Input.AtomFF[ii].Name.c_str(), Input.AtomFF[jj].Name.c_str(), Input.Mix_Epsilon[i], Input.Mix_Sigma[i], Input.Mix_Shift[i]);
  //}
  //printf("----- END OF MIXED VDW PARAMETERS -----\n");
}

void Copy_InputLoader_Data(Variables& Vars)
{
  Vars.FF.epsilon = convert1DVectortoArray(Vars.Input.Mix_Epsilon);
  Vars.FF.sigma   = convert1DVectortoArray(Vars.Input.Mix_Sigma);
  Vars.FF.z       = convert1DVectortoArray(Vars.Input.Mix_Z);
  Vars.FF.shift   = convert1DVectortoArray(Vars.Input.Mix_Shift);
  Vars.FF.FFType  = convert1DVectortoArray(Vars.Input.Mix_Type);
  Vars.FF.size    = Vars.Input.AtomFF.size();

  Vars.FF.noCharges   = Vars.Input.noCharges;
  Vars.FF.CutOffVDW   = Vars.Input.CutOffVDW;
  Vars.FF.CutOffCoul  = Vars.Input.CutOffCoul;
  Vars.FF.VDWRealBias = Vars.Input.VDWRealBias;
  Vars.FF.OverlapCriteria = Vars.Input.OverlapCriteria;
  
  //Copy Tail Corrections//
  for(size_t j = 0; j < Vars.Input.Mix_Tail.size(); j++)
    if(Vars.Input.Mix_Tail[j].UseTail) 
      Vars.TempComponents.HasTailCorrection = true;
  Vars.TempComponents.TailCorrection    = Vars.Input.Mix_Tail;
}


void read_movies_stats_print(Components& SystemComponents, size_t sim)
{
  std::vector<std::string> termsScannedLined{};
  std::string str;
  std::ifstream file("simulation.input");

  while (std::getline(file, str))
  {
    if (str.find("MoviesEvery", 0) != std::string::npos)
    {
      Split_Tab_Space(termsScannedLined, str);
      sscanf(termsScannedLined[1].c_str(), "%zu", &SystemComponents.MoviesEvery);
    }
    if (str.find("PrintEvery", 0) != std::string::npos)
    {
      Split_Tab_Space(termsScannedLined, str);
      sscanf(termsScannedLined[1].c_str(), "%zu", &SystemComponents.PrintStatsEvery);
    }
    if (str.find("SaveOutputToFile", 0) != std::string::npos)
    {
      Split_Tab_Space(termsScannedLined, str);
      if(caseInSensStringCompare(termsScannedLined[1], "yes"))
      {
        std::filesystem::path cwd = std::filesystem::current_path();
        std::filesystem::path directoryName = cwd / "Output/";
        std::filesystem::create_directories(directoryName);
        std::string FILENAME = directoryName.string() + "System_"+\
                               std::to_string(sim) +"_"+\
                               SystemComponents.MoleculeName[0] +"_"+\
                               std::to_string(SystemComponents.NumberofUnitCells.x)+"_"+\
                               std::to_string(SystemComponents.NumberofUnitCells.y)+"_"+\
                               std::to_string(SystemComponents.NumberofUnitCells.z)+"_"+\
                               std::to_string(SystemComponents.Temperature) + "_" +\
                               std::to_string(SystemComponents.Pressure) + ".data";

        SystemComponents.OUTPUT = fopen(FILENAME.c_str(), "w");
        if(SystemComponents.OUTPUT == nullptr)
        {
          throw std::runtime_error("FAIL TO OPEN" + FILENAME + "\n");
        }
      }
    }
  }
  //printf("Writing Movies every %zu MC step(s) or cycle(s)\n", SystemComponents.MoviesEvery);
  //printf("Printing Loadings and energies every %zu MC step(s) or cycle(s)\n", SystemComponents.PrintStatsEvery);
  if(SystemComponents.OUTPUT != stderr) //printf("Saving Output to File!\n");
  file.close();
}


void PseudoAtomParser(PseudoAtomDefinitions& PseudoAtom)
{
  std::string scannedLine; std::string str;
  std::vector<std::string> termsScannedLined{};
  size_t counter = 0;
  std::ifstream PseudoAtomfile("pseudo_atoms.def");
  size_t NumberOfPseudoAtoms = 0;
  // First read the pseudo atom file
  while (std::getline(PseudoAtomfile, str))
  {
    if(counter == 1) //read number definitions
    {
      Split_Tab_Space(termsScannedLined, str);
      sscanf(termsScannedLined[0].c_str(), "%zu", &NumberOfPseudoAtoms); //printf("THERE ARE %zu PSEUDO ATOMS\n", NumberOfPseudoAtoms);
      if (NumberOfPseudoAtoms <= 0) throw std::runtime_error("Incorrect amount of pseudo-atoms");//DON'T DO TOO FEW
      //if (NumberOfPseudoAtoms > 100) //printf("You are using A LOT OF pseudo atom definitions. Do you need so many? Okay...\n");
    }
    else if(counter >= 3) // read data for each pseudo atom
    {
      Split_Tab_Space(termsScannedLined, str);
      if(termsScannedLined[0] != PseudoAtom.Name[counter-3]) throw std::runtime_error("Order of pseudo-atom and force field definition don't match!");

      PseudoAtom.Symbol.push_back(termsScannedLined[2]);
      //size_t SymbolIdx = PseudoAtom.MatchUniqueSymbolTypeFromSymbolName(termsScannedLined[2]);
      //if(SymbolIdx >= PseudoAtom.UniqueSymbol.size()) PseudoAtom.UniqueSymbol.push_back(termsScannedLined[2]);

      //PseudoAtom.SymbolIndex.push_back(SymbolIdx);
      PseudoAtom.oxidation.push_back(std::stod(termsScannedLined[4]));
      PseudoAtom.mass.push_back(std::stod(termsScannedLined[5]));
      PseudoAtom.charge.push_back(std::stod(termsScannedLined[6]));
      PseudoAtom.polar.push_back(std::stod(termsScannedLined[7]));
    }
    counter++;
    if(counter==3+NumberOfPseudoAtoms) break; //in case there are extra empty rows
  }
  //print out the values
  //printf("-------------PARSING PSEUDO ATOMS FILE-------------\n");
  for (size_t i = 0; i < NumberOfPseudoAtoms; i++)
    //printf("Name: %s, %.5f, %.5f, %.5f, %.5f\n", PseudoAtom.Name[i].c_str(), PseudoAtom.oxidation[i], PseudoAtom.mass[i], PseudoAtom.charge[i], PseudoAtom.polar[i]);
  //printf("---------------------------------------------------\n");
  PseudoAtomfile.close();
}

void PseudoAtomProcessing(Variables& Vars)
{
  PseudoAtomDefinitions& PseudoAtom = Vars.PseudoAtoms;
  size_t NumberOfPseudoAtoms = Vars.PseudoAtoms.Name.size();
  if (NumberOfPseudoAtoms != Vars.FF.size) throw std::runtime_error("Number of VDW and pseudo-atom definitions don't match!");
  for(size_t i = 0; i < NumberOfPseudoAtoms; i++)
  {
    //Match 1-to-1 list of pseudo_atom type and symbol type//
    std::string& Symbol = PseudoAtom.Symbol[i];
    size_t SymbolIdx = PseudoAtom.MatchUniqueSymbolTypeFromSymbolName(Symbol);
    if(SymbolIdx >= PseudoAtom.UniqueSymbol.size()) PseudoAtom.UniqueSymbol.push_back(Symbol);
    PseudoAtom.SymbolIndex.push_back(SymbolIdx);
  }
}

void remove_number(std::string& s)
{
  s.erase(std::remove_if(std::begin(s), std::end(s), [](auto ch) { return std::isdigit(ch); }), s.end());
}

void remove_number_at_the_end(std::string& s)
{
  auto it = std::find_if(s.rbegin(), s.rend(), [](char ch) { return !std::isdigit(ch); });
  s.erase(it.base(), s.end());
}


void DetermineFrameworkComponent(Components& SystemComponents, size_t AtomCountPerUnitcell, size_t& Atom_Comp, size_t& MolID)
{ //Default Atom_Comp = 0; MolID = 0;
  if(SystemComponents.NComponents.y <= 1) return; //Then no need to separate
  for(size_t i = 1; i < SystemComponents.FrameworkComponentDef.size(); i++) //Component//
  {
    for(size_t j = 0; j < SystemComponents.FrameworkComponentDef[i].Atom_Indices_for_Molecule.size(); j++) //Molecule//
      for(size_t k = 0; k < SystemComponents.FrameworkComponentDef[i].Atom_Indices_for_Molecule[j].size(); k++) //Atom//
        if(AtomCountPerUnitcell == SystemComponents.FrameworkComponentDef[i].Atom_Indices_for_Molecule[j][k])
        {
          Atom_Comp = i; MolID = j; return;
        }
  }
}

/*
 * Checks and corrects the ordering of atoms in system components based on atom indices.
 * The function swaps the atom properties to ensure the order from `SystemComponents`
 * matches that from `unit_AtomIndex`.
 *
 * @param SystemComponents Contains the definition of each component in the system.
 * @param unit_fpos Vector containing the fractional positions of atoms.
 * @param unit_scale Vector containing the scaling factor for atoms.
 * @param unit_charge Vector containing the charge of atoms.
 * @param unit_scaleCoul Vector containing the Coulombic scaling factor for atoms.
 * @param unit_Type Vector containing the type ID of atoms.
 * @param unit_MolID Vector containing the molecule ID for atoms.
 * @param unit_AtomIndex Vector containing the atom indices.
 */

void CheckFrameworkComponentAtomOrder(Components& SystemComponents, std::vector<std::vector<double3>>& unit_fpos, std::vector<std::vector<double>>& unit_scale, std::vector<std::vector<double>>& unit_charge, std::vector<std::vector<double>>& unit_scaleCoul, std::vector<std::vector<size_t>>& unit_Type, std::vector<std::vector<size_t>>& unit_MolID, std::vector<std::vector<size_t>>& unit_AtomIndex)
{
  if(SystemComponents.NComponents.y <= 1) return; //Then no need to Check
  for(size_t i = 1; i < SystemComponents.FrameworkComponentDef.size(); i++) //Component//
  {
    size_t AtomCounter = 0; //Record the Atom indices
    for(size_t j = 0; j < SystemComponents.FrameworkComponentDef[i].Atom_Indices_for_Molecule.size(); j++) //Molecule//
      for(size_t k = 0; k < SystemComponents.FrameworkComponentDef[i].Atom_Indices_for_Molecule[j].size(); k++) //Atom//
      {
        for(size_t l = 0; l < unit_AtomIndex[i].size(); l++)
        {
          if(SystemComponents.FrameworkComponentDef[i].Atom_Indices_for_Molecule[j][k] == unit_AtomIndex[i][l])
          { 
            size_t From = l; size_t To = AtomCounter;
            double3 temp_fpos      = unit_fpos[i][To];
            double  temp_scale     = unit_scale[i][To];
            double  temp_charge    = unit_charge[i][To];
            double  temp_scaleCoul = unit_scaleCoul[i][To];
            size_t  temp_Type      = unit_Type[i][To];
            size_t  temp_MolID     = unit_MolID[i][To];
            size_t  temp_Index     = unit_AtomIndex[i][To];
            //Copy the current one to the correct location// 
            unit_fpos[i][To]      = unit_fpos[i][From];
            unit_scale[i][To]     = unit_scale[i][From];
            unit_charge[i][To]    = unit_charge[i][From];
            unit_scaleCoul[i][To] = unit_scaleCoul[i][From];
            unit_Type[i][To]      = unit_Type[i][From];
            unit_MolID[i][To]     = unit_MolID[i][From];
            unit_AtomIndex[i][To] = unit_AtomIndex[i][From];
            //Copy the replaced data to the New location (a swap)
            unit_fpos[i][From]      = temp_fpos;
            unit_scale[i][From]     = temp_scale;
            unit_charge[i][From]    = temp_charge;
            unit_scaleCoul[i][From] = temp_scaleCoul;
            unit_Type[i][From]      = temp_Type;
            unit_MolID[i][From]     = temp_MolID;
            unit_AtomIndex[i][From] = temp_Index;
            if(j != unit_MolID[i][To])
              throw std::runtime_error("MolID from Framework Component and MolID from CIF doesn't match!!!!"); 
          }
        }
        AtomCounter ++;
      }
  }
  
} 

void CheckFrameworkCIF(Boxsize& Box, PseudoAtomDefinitions& PseudoAtom, std::string& Frameworkfile, bool UseChargeFromCIF, double3 NumberUnitCells, Components& SystemComponents)
{
  std::vector<std::string> termsScannedLined{};
  termsScannedLined = split(Frameworkfile, '.');
  std::string frameworkName = termsScannedLined[0];
  std::string CIFFile = frameworkName + ".cif";
  std::ifstream simfile(CIFFile);
  std::filesystem::path pathfile = std::filesystem::path(CIFFile);
  if (!std::filesystem::exists(pathfile))
  {
    throw std::runtime_error("CIF file [ " + CIFFile + "] not found\n");
  }
  std::string str;

  //Temp Vector For Counting Number of Pseudo Atoms in the framework//
  int2 tempint2 = {0, 0};
  std::vector<std::vector<int2>>TEMPINTTWO;
  for(size_t i = 0; i < SystemComponents.NComponents.y; i++)
  {
    std::vector<int2>INTTWO(PseudoAtom.Name.size(), tempint2);
    TEMPINTTWO.push_back(INTTWO);
  }

  size_t NMol_Framework = 0;

  //Read angles, cell lengths, and volume//
  double3 angle; //x = alpha; y = beta; z = gamma;
  double3 abc; 
  double3 axbycz; //Diagonal of the matrix//
  while (std::getline(simfile, str))
  {
    if (str.find("_cell_length_a", 0) != std::string::npos)
    {
      Split_Tab_Space(termsScannedLined, str);
      abc.x = std::stod(termsScannedLined[1]);
    }
    if (str.find("_cell_length_b", 0) != std::string::npos)
    {
      Split_Tab_Space(termsScannedLined, str);
      abc.y = std::stod(termsScannedLined[1]);
    }
    if (str.find("_cell_length_c", 0) != std::string::npos)
    {
      Split_Tab_Space(termsScannedLined, str);
      abc.z = std::stod(termsScannedLined[1]);
    }
    if (str.find("_cell_angle_alpha", 0) != std::string::npos)
    {
      Split_Tab_Space(termsScannedLined, str);
      angle.x = std::stod(termsScannedLined[1]) / (180.0/3.14159265358979323846);
    }
    if (str.find("_cell_angle_beta", 0) != std::string::npos)
    {
      Split_Tab_Space(termsScannedLined, str);
      angle.y = std::stod(termsScannedLined[1]) / (180.0/3.14159265358979323846);
    }
    if (str.find("_cell_angle_gamma", 0) != std::string::npos)
    {
      Split_Tab_Space(termsScannedLined, str);
      angle.z = std::stod(termsScannedLined[1]) / (180.0/3.14159265358979323846);
    }
  }
  //Get xy(bx), xz(cx), and yz(cy)//
  axbycz.x = abc.x;
  axbycz.y = abc.y * std::sin(angle.z);
  double tempd = (std::cos(angle.x)-std::cos(angle.z)*std::cos(angle.y))/std::sin(angle.z);
  axbycz.z = abc.z * sqrt(1 - pow(std::cos(angle.y), 2) - pow(tempd, 2));
  double bx = abc.y * std::cos(angle.z);
  double cx = abc.z * std::cos(angle.y);
  double cy = abc.z * tempd;
  Box.Cell = (double*) malloc(9 * sizeof(double));
  Box.Cell[0] = NumberUnitCells.x * axbycz.x; Box.Cell[1] = 0.0;                          Box.Cell[2] = 0.0;
  Box.Cell[3] = NumberUnitCells.y * bx;       Box.Cell[4] = NumberUnitCells.y * axbycz.y; Box.Cell[5] = 0.0;
  Box.Cell[6] = NumberUnitCells.z * cx;       Box.Cell[7] = NumberUnitCells.z * cy;       Box.Cell[8] = NumberUnitCells.z * axbycz.z;
  simfile.clear();
  simfile.seekg(0);

  // Check the atom_site keyword: get its first and last occurance//
  //.x is the start, .y is the end//
  int2 atom_site_occurance; atom_site_occurance.x = 0; atom_site_occurance.y = 0; int count=0;
  int atomsiteCount = 0;
  bool foundline = false;
  //Initialize the order of the required columns in the CIF file with -1 (meaning non-existent)
  std::vector<int>Label_x_y_z_Charge_Order(5,-1);
  while (std::getline(simfile, str))
  {
    if (str.find("_atom_site", 0) != std::string::npos)
    {
      //Zhao's note: this relies on the fact that a cif file cannot start with _atom_site as the first line//
      if(atom_site_occurance.x == 0) atom_site_occurance.x = count;
      foundline = true;
      if(str.find("_atom_site_label", 0) != std::string::npos){  Label_x_y_z_Charge_Order[0] = atomsiteCount;}
      if(str.find("_atom_site_fract_x", 0) != std::string::npos){  Label_x_y_z_Charge_Order[1] = atomsiteCount;}
      if(str.find("_atom_site_fract_y", 0) != std::string::npos){  Label_x_y_z_Charge_Order[2] = atomsiteCount;}
      if(str.find("_atom_site_fract_z", 0) != std::string::npos){  Label_x_y_z_Charge_Order[3] = atomsiteCount;}
      if(str.find("_atom_site_charge", 0) != std::string::npos){ Label_x_y_z_Charge_Order[4] = atomsiteCount;}
      atomsiteCount++;
    }
    else //if cannot find _atom_site in this line, check if the "foundline" variable (meaning that if the keyword can be found in the previous line//
    {
      if(foundline)//Zhao's note: the assumes that "atom_site" only appears in one not multiple regions in the cif file
      {
        atom_site_occurance.y = count-1;
        break;
      }
    }
    count++;
  }
  //If label/x/y/z not found, abort the program!//
  if(Label_x_y_z_Charge_Order[0] < 0 || Label_x_y_z_Charge_Order[1] < 0 || Label_x_y_z_Charge_Order[2] < 0 || Label_x_y_z_Charge_Order[3] < 0) throw std::runtime_error("Couldn't find required columns in the CIF file! Abort.");
  //printf("atom_site starts at line %d, and ends at %d\n", atom_site_occurance.x, atom_site_occurance.y);
  simfile.clear();
  simfile.seekg(0);

  //Loop Over the Atoms//
  //size_t AtomID = 0;
  int label_location  = Label_x_y_z_Charge_Order[0];
  int x_location      = Label_x_y_z_Charge_Order[1];
  int y_location      = Label_x_y_z_Charge_Order[2];
  int z_location      = Label_x_y_z_Charge_Order[3];
  int Charge_location = Label_x_y_z_Charge_Order[4];
  if(!UseChargeFromCIF) Charge_location = -1; //If we want to use Charge from the pseudo_atoms.def, sed charge_location to -1//
  //printf("label location: %d, xyz location: %d %d %d, charge: %d\n", label_location, x_location, y_location, z_location, Charge_location);
  std::vector<std::vector<double3>>super_pos;
  std::vector<std::vector<double>>super_scale;
  std::vector<std::vector<double>>super_charge;
  std::vector<std::vector<double>>super_scaleCoul;
  std::vector<std::vector<size_t>>super_Type;
  std::vector<std::vector<size_t>>super_MolID;

  //Zhao's note: temporary values for storing unit cell//
  std::vector<std::vector<double3>>unit_fpos;
  std::vector<std::vector<double>>unit_scale;
  std::vector<std::vector<double>>unit_charge;
  std::vector<std::vector<double>>unit_scaleCoul;
  std::vector<std::vector<size_t>>unit_Type;
  std::vector<std::vector<size_t>>unit_MolID;
  std::vector<std::vector<size_t>>unit_AtomIndex; //For checking the order of the atoms, used for framework components//
  for(size_t i = 0; i < SystemComponents.NComponents.y; i++)
  {
    super_pos.emplace_back();
    super_scale.emplace_back();
    super_charge.emplace_back();
    super_scaleCoul.emplace_back();
    super_Type.emplace_back();
    super_MolID.emplace_back();

    unit_fpos.emplace_back();
    unit_scale.emplace_back();
    unit_charge.emplace_back();
    unit_scaleCoul.emplace_back();
    unit_Type.emplace_back();
    unit_MolID.emplace_back();
    unit_AtomIndex.emplace_back();
    ////printf("super_pos size: %zu\n", super_pos.size());
  }
  size_t i = 0;
  ////////////////////////////////////////////////////////////////
  //Zhao's note: For making supercells, the code can only do P1!//
  ////////////////////////////////////////////////////////////////
  double3 Shift;
  Shift.x = (double)1/NumberUnitCells.x; Shift.y = (double)1/NumberUnitCells.y; Shift.z = (double)1/NumberUnitCells.z;
  size_t AtomCountPerUnitcell = 0;
  // initiate total mass as a vector of doubles with size of framework components
  std::vector<double> totalmass(SystemComponents.NComponents.y, 0.0);
  while (std::getline(simfile, str))
  {
    if(i <= atom_site_occurance.y){i++; continue;}
    //AtomID = i - atom_site_occurance.y - 1;
    Split_Tab_Space(termsScannedLined, str);
    //Check when the elements in the line is less than 4, stop when it is less than 4 (it means the atom_site region is over//
    if(termsScannedLined.size() < 4) 
    {
      ////printf("Done reading Atoms in CIF file, there are %d Atoms\n", AtomID);
      break;
    }
    ////printf("i: %zu, line: %s, splitted: %s\n", i, str.c_str(), termsScannedLined[0].c_str());
    double3 fpos;
    fpos.x = std::stod(termsScannedLined[x_location]); //fx*=Shift.x; //Divide by the number of UnitCells
    fpos.y = std::stod(termsScannedLined[y_location]); //fy*=Shift.y; 
    fpos.z = std::stod(termsScannedLined[z_location]); //fz*=Shift.z;
    double Charge = 0.0;
    if(Charge_location >= 0) Charge = std::stod(termsScannedLined[Charge_location]);

    //DETERMINE WHAT Framework COMPONENT THIS ATOM BELONGS TO//
    size_t ATOM_COMP = 0;
    size_t MoleculeID = 0;
    DetermineFrameworkComponent(SystemComponents, AtomCountPerUnitcell, ATOM_COMP, MoleculeID);

    std::string AtomName = termsScannedLined[label_location];
    //Zhao's note: try to remove numbers from the Atom labels//
    //remove_number(AtomName);
    remove_number_at_the_end(AtomName); //Just remove the ending index
    size_t AtomTypeInt = 0;
    double AtomMass = 0.0;
    //Get the type (int) for this AtomName//
    bool AtomTypeFOUND = false;
    for(size_t j = 0; j < PseudoAtom.Name.size(); j++)
    {
      if(AtomName == PseudoAtom.Name[j])
      {
        AtomTypeInt = j;
        AtomMass = PseudoAtom.mass[j];
        AtomTypeFOUND = true;
        if(!UseChargeFromCIF) Charge = PseudoAtom.charge[j];
        // add to the totalmass of the component
        totalmass[ATOM_COMP] += AtomMass;
        //Add to the number of pseudo atoms
        TEMPINTTWO[ATOM_COMP][j].x =j;
        TEMPINTTWO[ATOM_COMP][j].y ++;
        break;
      }
    }
    //printf("Atom Count %zu, Component %zu, AtomType %zu (%s), MoleculeID %zu, AtomMass %f\n", AtomCountPerUnitcell, ATOM_COMP, AtomTypeInt, AtomName.c_str(), MoleculeID, AtomMass);
    if(!AtomTypeFOUND)throw std::runtime_error("Error: Atom Label [" + AtomName + "] not defined!");
    unit_fpos[ATOM_COMP].push_back(fpos);
    unit_scale[ATOM_COMP].push_back(1.0); //For framework, use 1.0
    unit_charge[ATOM_COMP].push_back(Charge);
    unit_scaleCoul[ATOM_COMP].push_back(1.0);//For framework, use 1.0
    unit_Type[ATOM_COMP].push_back(AtomTypeInt);
    unit_MolID[ATOM_COMP].push_back(MoleculeID);
    unit_AtomIndex[ATOM_COMP].push_back(AtomCountPerUnitcell);
    AtomCountPerUnitcell ++;
    i++;
  }
  //Check total mass of framework components//
  for(size_t asd = 0; asd < totalmass.size(); asd++)
    //printf("component %zu, totalmass: %.5f\n", asd, totalmass[asd]);
  //Zhao's note: Need to sort the atom positions for the framework component to match the order in Framework_Component definition files, see Hilal's example//
  CheckFrameworkComponentAtomOrder(SystemComponents, unit_fpos, unit_scale, unit_charge, unit_scaleCoul, unit_Type, unit_MolID, unit_AtomIndex);
  
  //Zhao's note: consider separating the duplication of supercell from reading unit cell values//
  //It may cause strange bugs if the MolID for atoms are not continuous//
  for(size_t comp = 0; comp < SystemComponents.NComponents.y; comp++)
    for(size_t ix = 0; ix < NumberUnitCells.x; ix++)
      for(size_t jy = 0; jy < NumberUnitCells.y; jy++)
        for(size_t kz = 0; kz < NumberUnitCells.z; kz++)
        {
          size_t UnitCellID = (ix * NumberUnitCells.y + jy) * NumberUnitCells.z + kz;
          double3 NCell = {(double) ix, (double) jy, (double) kz};
          for(size_t Atom = 0; Atom < unit_fpos[comp].size(); Atom++)
          {
            double3 fpos = unit_fpos[comp][Atom];
            //Get supercell fx, fy, fz, and get corresponding xyz//
            double3 super_fpos = (fpos + NCell) * Shift;
            // Get real xyz from fractional xyz //
            double3 pos;
            pos.x = super_fpos.x*Box.Cell[0]+super_fpos.y*Box.Cell[3]+super_fpos.z*Box.Cell[6];
            pos.y = super_fpos.x*Box.Cell[1]+super_fpos.y*Box.Cell[4]+super_fpos.z*Box.Cell[7];
            pos.z = super_fpos.x*Box.Cell[2]+super_fpos.y*Box.Cell[5]+super_fpos.z*Box.Cell[8];
            super_pos[comp].push_back(pos);
            double scale = unit_scale[comp][Atom];
            double charge= unit_charge[comp][Atom];
            double scaleCoul= unit_scaleCoul[comp][Atom];
            size_t Type = unit_Type[comp][Atom];
            
            super_scale[comp].push_back(scale);
            super_charge[comp].push_back(charge);
            super_scaleCoul[comp].push_back(scaleCoul);
            super_Type[comp].push_back(Type);
            size_t ActualMolID = unit_MolID[comp][Atom];
            //Only duplicate MolID if the component is separated from the CIF file//
            if(SystemComponents.FrameworkComponentDef[comp].SeparatedComponent)
            {
              ActualMolID = SystemComponents.FrameworkComponentDef[comp].Number_of_Molecules_for_Framework_component * UnitCellID + ActualMolID;
            }
            super_MolID[comp].push_back(ActualMolID);
          }
        }
  //printf("Finished Reading Atoms\n");
  // initiate overall size of the framework //
  double OverallMass = 0.0;
  for(size_t i = 0; i < SystemComponents.NComponents.y; i++)
  {
    SystemComponents.HostSystem[i].pos       = (double3*) malloc(super_pos[i].size() * sizeof(double3));
    SystemComponents.HostSystem[i].scale     = (double*)  malloc(super_pos[i].size() * sizeof(double));
    SystemComponents.HostSystem[i].charge    = (double*)  malloc(super_pos[i].size() * sizeof(double));
    SystemComponents.HostSystem[i].scaleCoul = (double*)  malloc(super_pos[i].size() * sizeof(double));
    SystemComponents.HostSystem[i].Type      = (size_t*)  malloc(super_pos[i].size() * sizeof(size_t));
    SystemComponents.HostSystem[i].MolID     = (size_t*)  malloc(super_pos[i].size() * sizeof(size_t));
    
    //std::vector<size_t> MolID_i = super_MolID[i];
    //size_t NMol = std::max(MolID_i);
    std::vector<size_t>::iterator maxValueIterator = std::max_element(super_MolID[i].begin(), super_MolID[i].end());
    size_t NMol = 0;
    //Zhao's note: add protection here for reading empty cif file for an empty box//
    if(super_pos[i].size() != 0)
    {
      NMol = *maxValueIterator;
    }
    NMol ++;
    //size_t NMol = (*std::max_element(begin(super_MolID), end(super_MolID), [](size_t& a, size_t& b){ return a[i] < b[i]; }))[i];
    size_t NMol_In_Def = SystemComponents.FrameworkComponentDef[i].Number_of_Molecules_for_Framework_component;
    if(SystemComponents.FrameworkComponentDef[i].SeparatedComponent)
    {
      NMol_In_Def *= NumberUnitCells.x * NumberUnitCells.y * NumberUnitCells.z;
      // in the case there is multiple unit cells, make sure that the mass is the total mass of the box //
      //totalmass[i] *= NMol_In_Def;
    }
    double unitcells = NumberUnitCells.x * NumberUnitCells.y * NumberUnitCells.z;
    //totalmass[i] *= unitcells;

    //printf("NMol = %zu, pos_size: %zu, NMol in FrameworkDef: %zu\n", NMol, super_pos[i].size(), NMol_In_Def);

    if(NMol != NMol_In_Def)
    {
      throw std::runtime_error("In CheckFrameworkCIF function, NMol and value in FrameworkComponentDef don't match!!!!\n");
    }
    SystemComponents.HostSystem[i].size          = super_pos[i].size();
    SystemComponents.HostSystem[i].Molsize       = super_pos[i].size() / NMol;
    SystemComponents.HostSystem[i].Allocate_size = super_pos[i].size();

    SystemComponents.Moleculesize.push_back(SystemComponents.HostSystem[i].Molsize);
    SystemComponents.Allocate_size.push_back(SystemComponents.HostSystem[i].size);
    SystemComponents.NumberOfMolecule_for_Component.push_back(NMol);

    NMol_Framework += NMol;

    SystemComponents.MolecularWeight.push_back(totalmass[i]);
    // gather all the mass of the framework //
    OverallMass += totalmass[i];

    //printf("Framework Comp [%zu], size: %zu, Molsize: %zu, Allocate_size: %zu, component mass: %f\n", i, super_pos[i].size(), SystemComponents.HostSystem[i].Molsize, SystemComponents.HostSystem[i].Allocate_size, totalmass[i]);
  }
  for(size_t i = 0; i < super_pos.size(); i++)
  {
    bool FrameworkHasPartialCharge = false; double ChargeSum = 0.0;
    for(size_t j = 0; j < super_pos[i].size(); j++)
    {
      SystemComponents.HostSystem[i].pos[j] = super_pos[i][j];
      SystemComponents.HostSystem[i].scale[j] = super_scale[i][j];
      SystemComponents.HostSystem[i].charge[j] = super_charge[i][j];
      SystemComponents.HostSystem[i].scaleCoul[j] = super_scaleCoul[i][j];
      SystemComponents.HostSystem[i].Type[j] = super_Type[i][j];
      SystemComponents.HostSystem[i].MolID[j] = super_MolID[i][j];
      ChargeSum += std::abs(SystemComponents.HostSystem[i].charge[j]);
    }
    if(ChargeSum > 1e-6) FrameworkHasPartialCharge = true;
    SystemComponents.hasPartialCharge.push_back(FrameworkHasPartialCharge);
  }

  //printf("------------------CIF FILE SUMMARY------------------\n");
  //printf("CIF FILE IS: %s\n", CIFFile.c_str());
  //printf("Number of Unit Cells: %.2f %.2f %.2f\n", NumberUnitCells.x, NumberUnitCells.y, NumberUnitCells.z);
  //printf("Box size: \n%.5f %.5f %.5f\n%.5f %.5f %.5f\n%.5f %.5f %.5f\n", Box.Cell[0], Box.Cell[1], Box.Cell[2], Box.Cell[3], Box.Cell[4], Box.Cell[5], Box.Cell[6], Box.Cell[7], Box.Cell[8]);

  //Record Number of PseudoAtoms for the Framework//
  for(size_t i = 0; i < TEMPINTTWO.size(); i++)
  {
    std::vector<int2>NumberOfPseudoAtoms;
    for(size_t j = 0; j < TEMPINTTWO[i].size(); j++)
    {
      if(TEMPINTTWO[i][j].y == 0) continue;
      TEMPINTTWO[i][j].y *= NumberUnitCells.x * NumberUnitCells.y * NumberUnitCells.z;
      NumberOfPseudoAtoms.push_back(TEMPINTTWO[i][j]);
      //printf("Framework Pseudo Atom[%zu], Name: %s, #: %zu\n", TEMPINTTWO[i][j].x, PseudoAtom.Name[j].c_str(), TEMPINTTWO[i][j].y);
    }
    //printf("NumberOfPseudoAtoms size: %zu\n", NumberOfPseudoAtoms.size());
    SystemComponents.NumberOfPseudoAtomsForSpecies.push_back(NumberOfPseudoAtoms);
  }
  for(size_t i = 0; i < SystemComponents.NumberOfPseudoAtomsForSpecies.size(); i++)
    for(size_t j = 0; j < SystemComponents.NumberOfPseudoAtomsForSpecies[i].size(); j++)
    {
      //printf("Framework Component [%zu], Pseudo Atom [%zu], Name: %s, #: %zu\n", i, SystemComponents.NumberOfPseudoAtomsForSpecies[i][j].x, PseudoAtom.Name[SystemComponents.NumberOfPseudoAtomsForSpecies[i][j].x].c_str(), SystemComponents.NumberOfPseudoAtomsForSpecies[i][j].y);
    }
  //printf("Overall Mass: %f \n", OverallMass);
  //Add PseudoAtoms from the Framework to the total PseudoAtoms array//
  for(size_t i = 0; i < SystemComponents.NComponents.y; i++)
    SystemComponents.UpdatePseudoAtoms(INSERTION, i);

  SystemComponents.TotalNumberOfMolecules = NMol_Framework; //If there is a framework, framework is counted as a molecule//
  SystemComponents.NumberOfFrameworks = NMol_Framework;
  //printf("----------------------------------------------------\n");
  //throw std::runtime_error("BREAK!!!!\n");
  simfile.close();
}

void ReadFrameworkSpeciesDefinitions(Components& SystemComponents)
{
  SystemComponents.FrameworkComponentDef.resize(SystemComponents.NComponents.y);
  if(SystemComponents.NComponents.y <= 1)
  {
    SystemComponents.FrameworkComponentDef[0].Number_of_Molecules_for_Framework_component = 1;
    return; //Then no need to separate
  }
  std::string scannedLine; std::string str;
  std::vector<std::string> termsScannedLined{};
  size_t tempval = 0;
  for(size_t i = 0; i < SystemComponents.NComponents.y; i++)
  { //Zhao's note:: for framework component value = 0, the default framework component, we may want to use a file_check to see if it is separatable//
    if(i == 0) 
    {
      SystemComponents.FrameworkComponentDef[0].SeparatedComponent = false;
      SystemComponents.FrameworkComponentDef[0].Number_of_Molecules_for_Framework_component = 1;
      continue; //No need for the default framework component
    }
    std::string FileName = "Framework_Component_" + std::to_string(i) + ".def";
    //printf("Reading %s FILE\n", FileName.c_str());
    std::filesystem::path pathfile = std::filesystem::path(FileName);
    if (!std::filesystem::exists(pathfile))
    {
      throw std::runtime_error("Framework Component NOT FOUND!!!!\n");
    }
    SystemComponents.FrameworkComponentDef[i].SeparatedComponent = true;
    std::ifstream file(FileName);
    while (std::getline(file, str))
    {
      if (str.find("#", 0) != std::string::npos) continue;
      if (str.find("Framework_Component_Name", 0) != std::string::npos)
      {
        Split_Tab_Space(termsScannedLined, str);
        SystemComponents.MoleculeName.push_back(termsScannedLined[1]);
      }

      if (str.find("Number_of_Molecules_for_Framework_component", 0) != std::string::npos)
      {
        Split_Tab_Space(termsScannedLined, str);
        sscanf(termsScannedLined[1].c_str(), "%zu", &tempval);
        SystemComponents.FrameworkComponentDef[i].Number_of_Molecules_for_Framework_component = tempval;
        //SystemComponents.FrameworkComponentDef[i].Atom_Indices_for_Molecule.resize(tempval);
      }
      if (str.find("Number_of_atoms_for_each_molecule", 0) != std::string::npos)
      {
        Split_Tab_Space(termsScannedLined, str);
        sscanf(termsScannedLined[1].c_str(), "%zu", &tempval);
        SystemComponents.FrameworkComponentDef[i].Number_of_atoms_for_each_molecule = tempval;
      }
      if (str.find("Atom_Indices_for_Molecule", 0) != std::string::npos)
      {
        Split_Tab_Space(termsScannedLined, str);
        //printf("size of the indices: %zu\n", termsScannedLined.size() - 2);
        std::vector<size_t>Indices;
        for(size_t j = 2; j < termsScannedLined.size(); j++)
        {
          sscanf(termsScannedLined[j].c_str(), "%zu", &tempval);
          Indices.push_back(tempval);
        }
        SystemComponents.FrameworkComponentDef[i].Atom_Indices_for_Molecule.push_back(Indices);
        if(Indices.size() != SystemComponents.FrameworkComponentDef[i].Number_of_atoms_for_each_molecule)
          throw std::runtime_error("Number of atoms for Framework Component != Number of Atom Indices!!!!\n");
      }
    }
    file.close();
    //printf("================Framework Component [%zu] Summary================\n", i);
    //printf("Name: %s\n", SystemComponents.MoleculeName[i].c_str());
    //printf("Number of Molecules: %zu\n", SystemComponents.FrameworkComponentDef[i].Number_of_Molecules_for_Framework_component);
    //printf("Number of Atoms per Molecule: %zu\n", SystemComponents.FrameworkComponentDef[i].Number_of_atoms_for_each_molecule);
    //for(size_t j = 0; j < SystemComponents.FrameworkComponentDef[i].Number_of_Molecules_for_Framework_component; j++)
    //  for(size_t k = 0; k < SystemComponents.FrameworkComponentDef[i].Number_of_atoms_for_each_molecule; k++)
        //printf("Molecule [%zu], Index: %zu\n", j, SystemComponents.FrameworkComponentDef[i].Atom_Indices_for_Molecule[j][k]);
  }
}

void ReadFramework(Boxsize& Box, PseudoAtomDefinitions& PseudoAtom, size_t FrameworkIndex, Components& SystemComponents)
{
  //printf("------------------------PARSING FRAMEWORK DATA------------------------\n");
  bool UseChargesFromCIFFile = false;  //By default, it is using charge from pseudo_atoms.def//
  std::vector<std::string> Names = PseudoAtom.Name;
  //size_t temp = 0;
  std::string scannedLine; std::string str;
  std::vector<std::string> termsScannedLined{};
  //size_t counter = 0;
  //Determine framework name (file name)//
  std::ifstream simfile("simulation.input");
  std::string Frameworkname;
  std::string InputType="cif";

  double3 NumberUnitCells;
  bool FrameworkFound = false;

  //Check input file for "UseChargeFromCIFFile" keyword//
  while (std::getline(simfile, str))
  {
    if (str.find("UseChargesFromCIFFile", 0) != std::string::npos)
    {
      Split_Tab_Space(termsScannedLined, str);
      if(caseInSensStringCompare(termsScannedLined[1], "yes"))
      {
        UseChargesFromCIFFile = true;
        //printf("Reading Partial Charges from pseudo_atoms.def file\n");
      }
    }
  }
  simfile.clear();
  simfile.seekg(0);
    
  while (std::getline(simfile, str))
  {
    if (str.find("InputFileType", 0) != std::string::npos)
    {
      Split_Tab_Space(termsScannedLined, str);
      InputType=termsScannedLined[1];
    }
    //Zhao's note: When we are using Multiple frameworks, in simulation.input, list the frameworks one by one.//
    if (str.find("FrameworkName", 0) != std::string::npos) // get the molecule name
    {
      Split_Tab_Space(termsScannedLined, str);
      if((1+FrameworkIndex) >= termsScannedLined.size()) throw std::runtime_error("Not Enough Framework listed in input file...\n");
      Frameworkname = termsScannedLined[1 + FrameworkIndex];
      //printf("Reading Framework %zu, FrameworkName: %s\n", FrameworkIndex, Frameworkname.c_str());
    }
    if (str.find("UnitCells " + std::to_string(FrameworkIndex), 0) != std::string::npos) // get the molecule name
    {
      Split_Tab_Space(termsScannedLined, str);
      NumberUnitCells.x = std::stod(termsScannedLined[2]);
      NumberUnitCells.y = std::stod(termsScannedLined[3]);
      NumberUnitCells.z = std::stod(termsScannedLined[4]);
      //printf("Reading Framework %zu, UnitCells: %.2f %.2f %.2f\n", FrameworkIndex, NumberUnitCells.x, NumberUnitCells.y, NumberUnitCells.z);
      FrameworkFound = true;
      SystemComponents.NumberofUnitCells = {(int) NumberUnitCells.x, (int) NumberUnitCells.y, (int) NumberUnitCells.z};
    }
  }
  //If not cif or poscar, break the program!//
  if(InputType != "cif" && InputType != "poscar") throw std::runtime_error("Cannot identify framework input type [" + InputType + "]. It can only be cif or poscar!");
  if(!FrameworkFound) throw std::runtime_error("Cannot find the framework with matching index!");
  std::string FrameworkFile = Frameworkname + "." + InputType;
  std::filesystem::path pathfile = std::filesystem::path(FrameworkFile);
  if (!std::filesystem::exists(pathfile)) throw std::runtime_error("Framework File ["+ FrameworkFile + "] not found!\n");
  /////////////////////////////////////////////////////////////////////////////////////
  //Zhao's note:                                                                     //
  //If reading poscar, then you cannot have user-defined charges for every atom      //
  //the charges used for poscar are defined in pseudo_atoms.def                      //
  //To use user-defined charge for every atom, use cif format                        //
  /////////////////////////////////////////////////////////////////////////////////////
  if(InputType == "cif")
  {
    SystemComponents.MoleculeName.push_back(FrameworkFile);

    ReadFrameworkSpeciesDefinitions(SystemComponents);
    CheckFrameworkCIF(Box, PseudoAtom, FrameworkFile, UseChargesFromCIFFile, NumberUnitCells, SystemComponents);
    //printf("Reading CIF File\n");
  }
  else
  {
    throw std::runtime_error("Only supports reading CIF files\n");
  }
  //Get Volume, cubic/non-cubic of the box//
  Box.InverseCell = (double*) malloc(9 * sizeof(double));
  inverse_matrix(Box.Cell, &Box.InverseCell);
  Box.Volume = matrix_determinant(Box.Cell);
  //DETERMINE Whether Box is cubic/cuboid or not//
  Box.Cubic = true; // Start with cubic box shape, if any value (non-diagonal) is greater than 0, then set to false //
  if((fabs(Box.Cell[3]) + fabs(Box.Cell[6]) + fabs(Box.Cell[7])) > 1e-10) Box.Cubic = false;
  if(Box.Cubic)  //printf("The Simulation Box is Cubic\n");
  if(!Box.Cubic) //printf("The Simulation Box is NOT Cubic\n");
  //printf("----------------------END OF PARSING FRAMEWORK DATA----------------------\n");
  simfile.close();
}

size_t get_type_from_name(std::string Name, std::vector<std::string> PseudoAtomNames)
{
  size_t type = 0; bool match = false;
  for(size_t i = 0; i < PseudoAtomNames.size(); i++)
  {
    if(Name == PseudoAtomNames[i]){
      type = i; match = true; break;}
  }
  if(!match) throw std::runtime_error("Atom type not found in pseudo atoms definitions\n");
  return type;
}

void MoleculeDefinitionParser(Atoms& Mol, Components& SystemComponents, std::string MolName, PseudoAtomDefinitions PseudoAtom, size_t Allocate_space)
{
  //check if molecule definition file exists
  const std::string MolFileName = MolName + ".def";
  std::filesystem::path pathfile = std::filesystem::path(MolFileName);
  if (!std::filesystem::exists(pathfile))
  {
    throw std::runtime_error("Definition file not found\n");
  }
  std::string scannedLine; std::string str;
  std::vector<std::string> termsScannedLined{};
  size_t counter = 0; size_t temp_molsize = 0; size_t atomcount = 0;
  std::ifstream file(MolFileName);

  Mol.Allocate_size = Allocate_space;
  std::vector<double3> Apos(Allocate_space);
  std::vector<double>  Ascale(Allocate_space);
  std::vector<double>  Acharge(Allocate_space);
  std::vector<double>  Amass(Allocate_space);
  std::vector<double>  AscaleCoul(Allocate_space);
  std::vector<size_t>  AType(Allocate_space);
  std::vector<size_t>  AMolID(Allocate_space);
 
  double chargesum = 0.0; //a sum of charge for the atoms in the molecule, for checking charge neutrality
  double masssum = 0.0;
 
  bool temprigid = false;

  size_t PseudoAtomSize = PseudoAtom.Name.size();
  //For Tail corrections//
  int2   TempINTTWO = {0, 0};
  std::vector<int2> ANumberOfPseudoAtomsForSpecies(PseudoAtomSize, TempINTTWO);

  // skip first line
  while (std::getline(file, str))
  {
    if(counter == 1) //read Tc
    {
      Split_Tab_Space(termsScannedLined, str);
      SystemComponents.Tc.push_back(std::stod(termsScannedLined[0]));
    }
    else if(counter == 2) //read Pc
    {
      Split_Tab_Space(termsScannedLined, str);
      SystemComponents.Pc.push_back(std::stod(termsScannedLined[0]));
    }
    else if(counter == 3) //read Accentric Factor
    {
      Split_Tab_Space(termsScannedLined, str);
      SystemComponents.Accentric.push_back(std::stod(termsScannedLined[0]));
    }
    else if(counter == 5) //read molecule size
    {
      temp_molsize = 0;
      Split_Tab_Space(termsScannedLined, str);
      sscanf(termsScannedLined[0].c_str(), "%zu", &temp_molsize);
      if(temp_molsize >= Allocate_space) throw std::runtime_error("Molecule size is greater than allocated size. Break\n");
      SystemComponents.Moleculesize.push_back(temp_molsize);
      Mol.Molsize = temp_molsize; //Set Molsize for the adsorbate molecule here//
    }
    else if(counter == 9) //read if the molecule is rigid
    {
      //printf("rigid line: %s\n", str.c_str());
      Split_Tab_Space(termsScannedLined, str);
      if(caseInSensStringCompare(termsScannedLined[0], "rigid"))
      {
        temprigid = true;
        //printf("Adsorbate Component is rigid\n");
      }
      else
      {
        throw std::runtime_error("Currently Not allowing flexible molecule\n");
      }
    }
    else if(counter >= 13 && atomcount < (temp_molsize)) //read atomic positions. Zhao's note: skipped reading groups for now
    {
      //for atomic positions, read them into Mol, at the positions for the first molecule.
      //Don't set the size of Mol, it is zero, since the data loaded is only a template.
      Split_Tab_Space(termsScannedLined, str);
      if(termsScannedLined.size() == 5) //for example: 0 CH4 0.0 0.0 0.0, position provided here.
      {
        Apos[atomcount] = {std::stod(termsScannedLined[2]), std::stod(termsScannedLined[3]), std::stod(termsScannedLined[4])};
      }
      else if(termsScannedLined.size() == 2 && temp_molsize == 1) //like methane, one can do: 0 CH4, with no positions
      {
        Apos[atomcount] = {0.0, 0.0, 0.0};
      }
      else
      {
        throw std::runtime_error("Flexible molecules not implemented\n");
      }
      //Set other values for each atom
      Ascale[atomcount]  = 1.0;
      AscaleCoul[atomcount] = 1.0;
      std::string AtomName = termsScannedLined[1];
      AType[atomcount] = get_type_from_name(AtomName, PseudoAtom.Name);
      ANumberOfPseudoAtomsForSpecies[AType[atomcount]].x = AType[atomcount];
      ANumberOfPseudoAtomsForSpecies[AType[atomcount]].y ++;
      Acharge[atomcount] = PseudoAtom.charge[AType[atomcount]]; //Determine the charge after determining the type
      //printf("%s, type: %zu, Acharge = %.5f\n", AtomName.c_str(), AType[atomcount], Acharge[atomcount]);
      Amass[atomcount] = PseudoAtom.mass[AType[atomcount]];
      chargesum += PseudoAtom.charge[AType[atomcount]]; 
      masssum   += PseudoAtom.mass[AType[atomcount]];
      AMolID[atomcount] = 0;// Molecule ID = 0, since it is in the first position
      atomcount++;
    }
    else if(counter > (13 + temp_molsize +1)) //read bonding information
    {
      //printf("Bonds not implemented. Break\n"); break;
    }
    counter++; 
  }

  bool MolHasPartialCharge = false; double ChargeSum = 0.0;
  for(size_t i = 0; i < Acharge.size(); i++) 
  {
    ChargeSum += std::abs(Acharge[i]);
  }
  if(ChargeSum > 1e-6) MolHasPartialCharge = true;
  SystemComponents.hasPartialCharge.push_back(MolHasPartialCharge);

  SystemComponents.rigid.push_back(temprigid);
  if(chargesum > 1e-10)
  {
    ////printf("chargesum = %.50f, ChargeSum = %.5f\n", chargesum, ChargeSum);
    throw std::runtime_error("Molecule not neutral, bad\n");
  }

  Mol.pos       = convert1DVectortoArray(Apos);
  Mol.scale     = convert1DVectortoArray(Ascale);
  Mol.charge    = convert1DVectortoArray(Acharge);
  Mol.scaleCoul = convert1DVectortoArray(AscaleCoul);

  Mol.Type      = convert1DVectortoArray(AType);
  Mol.MolID     = convert1DVectortoArray(AMolID);

  //for(size_t i = 0; i < Mol.Molsize; i++)
    //printf("Atom [%zu]: Type [%zu], Name: %s, Mass: %f, position: %.5f %.5f %.5f\n", i, Mol.Type[i], PseudoAtom.Name[Mol.Type[i]].c_str(), PseudoAtom.mass[Mol.Type[i]], Mol.pos[i].x, Mol.pos[i].y, Mol.pos[i].z);

  //Remove Elements from ANumberOfPseudoAtomsForSpecies if the ANumberOfPseudoAtomsForSpecies.y = 0
  std::vector<int2>TEMPINTTWO;
  for(size_t i = 0; i < ANumberOfPseudoAtomsForSpecies.size(); i++)
  {
    if(ANumberOfPseudoAtomsForSpecies[i].y == 0) continue;
    TEMPINTTWO.push_back(ANumberOfPseudoAtomsForSpecies[i]);
    //printf("Adsorbate Type[%zu], Name: %s, #: %zu\n", ANumberOfPseudoAtomsForSpecies[i].x, PseudoAtom.Name[i].c_str(), ANumberOfPseudoAtomsForSpecies[i].y);
  }
  //printf("current adsorbate mass is: %f \n", masssum);
  SystemComponents.NumberOfPseudoAtomsForSpecies.push_back(TEMPINTTWO);
  SystemComponents.MolecularWeight.push_back(masssum);
  file.close();
}

double process_str_double_DBLMIN(const std::string& str)
{
  std::istringstream iss(str);
  double value;
  if (iss >> value)
  {
      if(std::abs(value) < DBL_MIN)
        return 0.0;
      else return std::stod(str);
  }
  return 0.0;
}

static inline void Read_TMMC_Initial(TMMC& tmmc, std::string& MolName)
{
  std::string scannedLine; std::string str;
  std::vector<std::string> termsScannedLined{};
  //Determine framework name (file name)//
  std::string Filename = "TMMCInitial/System_0/TMMC_Statistics.data";
  std::ifstream file(Filename);
  std::filesystem::path pathfile = std::filesystem::path(Filename);
  if (!std::filesystem::exists(pathfile))
  {
    throw std::runtime_error("TMMCInitial file not found\n");
  }
  //Zhao's note: first read the number of lines in the TMMCInitial file
  int Histogramsize = -5;
  size_t counter = 0;
  size_t header  = 5;
  //Process the header of the TMMCInitial file//
  while (std::getline(file, str))
  {
    if(counter < header)
    {
    Split_Tab_Space(termsScannedLined, str);
    if(counter == 0)
    {
      //Record the number of TMUpdateTimes
      std::string ReadMolName = termsScannedLined[2];
      if(ReadMolName != MolName) throw std::runtime_error("TMMCInitial Mol name doesn't match what is in the simulation.input file. Please check!!!");
      tmmc.TMUpdateTimes = std::stoi(termsScannedLined[5]);
    }
    else if(counter == 1)
    {
      size_t temp_MinMacrostate = static_cast<size_t>(std::stoi(termsScannedLined[3]));
      if(temp_MinMacrostate != tmmc.MinMacrostate)
        throw std::runtime_error("TMMCInitial Min Macrostate doesn't match the value in simulation.input file!!!");
    }
    else if(counter == 2)
    {
      size_t temp_MaxMacrostate = static_cast<size_t>(std::stoi(termsScannedLined[3]));
      if(temp_MaxMacrostate != tmmc.MaxMacrostate)
        throw std::runtime_error("TMMCInitial Max Macrostate doesn't match the value in simulation.input file!!!");
    }
    else if(counter == 3)
    {
      tmmc.WLFactor = std::stoi(termsScannedLined[3]);
    }
    }
    ++Histogramsize;
    counter ++;
  }
  file.clear();
  file.seekg(0);
  //printf("Lines of TMMCInitial file: %d, Histogram size: %d\n", Histogramsize+5, Histogramsize);
  if(Histogramsize<= 0) throw std::runtime_error("TMMCInitial size is zero, there is no value written in the file. Check!!!");
  //Check histogram size and bin size//
  //From 0 to 60, there are 61 bins//
  size_t temp_nbin = (Histogramsize-1) / (tmmc.MaxMacrostate - tmmc.MinMacrostate);
  if(temp_nbin != tmmc.nbinPerMacrostate)
    throw std::runtime_error("Number of bins calculated from TMMCInitial doesn't match with what is specified in the input file. Maybe check your CBCFC (fractional molecule) setup!!!");
  if((Histogramsize-1) != (tmmc.MaxMacrostate - tmmc.MinMacrostate) * temp_nbin)
     throw std::runtime_error("Number of macrostates doesn't match! Your nbin: " + std::to_string(temp_nbin) + ", Histogram size: "+ std::to_string(Histogramsize) + ", Max-Min: " + std::to_string(tmmc.MaxMacrostate - tmmc.MinMacrostate));
  counter = 0;
  while (std::getline(file, str))
  {
    Split_Tab_Space(termsScannedLined, str);
    if(counter >= header) //Skip line 5 (column names)
    {
      ////printf("str is: %s\n", str.c_str());
      //Column names: N NMol Bin CM[-1] CM[0] CM[1] WLBias ln_g TMBias lnpi Forward_lnpi Reverse_lnpi Histogram
      size_t j   = static_cast<size_t>(std::stoi(termsScannedLined[0])); //Bin index//
      //size_t N   = std::stoi(termsScannedLined[1]);
      //size_t bin = std::stoi(termsScannedLined[2]);
      //Zhao's note: some values in the TMMC data file are extremely small, add a check here, if the printed value is smaller than DBL_MIN, make it DBL_MIN
      tmmc.CMatrix[j].x    = process_str_double_DBLMIN(termsScannedLined[3]);
      tmmc.CMatrix[j].y    = process_str_double_DBLMIN(termsScannedLined[4]);
      tmmc.CMatrix[j].z    = process_str_double_DBLMIN(termsScannedLined[5]);
      tmmc.WLBias[j]       = process_str_double_DBLMIN(termsScannedLined[6]);
      tmmc.ln_g[j]         = process_str_double_DBLMIN(termsScannedLined[7]);
      tmmc.TMBias[j]       = process_str_double_DBLMIN(termsScannedLined[8]);
      tmmc.lnpi[j]         = process_str_double_DBLMIN(termsScannedLined[9]);
      tmmc.forward_lnpi[j] = process_str_double_DBLMIN(termsScannedLined[10]);
      tmmc.reverse_lnpi[j] = process_str_double_DBLMIN(termsScannedLined[11]);
      tmmc.Histogram[j]    = process_str_double_DBLMIN(termsScannedLined[12]);
    }
    counter ++;
  }
  file.close();
}

void read_component_values_from_simulation_input(Variables& Vars, Components& SystemComponents, Move_Statistics& MoveStats, size_t AdsorbateComponent, Atoms& Mol, PseudoAtomDefinitions PseudoAtom, size_t Allocate_space, size_t BoxIndex)
{
  //adsorbate component start from zero, but in the code, framework is zero-th component
  //This function also calls MoleculeDefinitionParser//
  size_t Terminatecomponent = AdsorbateComponent+1;
  std::vector<std::string> termsScannedLined{};
  std::string str;
  std::ifstream file("simulation.input");
  int counter=0; size_t start_counter = 0;
  size_t CreateMolecule = 0; double idealrosen = 0.0; double fugacoeff = 0.0; double Molfrac = 1.0; //Set Molfraction = 1.0
  bool temp_hasfracmol = false;
  int  LambdaType = SHI_MAGINN;

  TMMC temp_tmmc;

  std::string start_string = "Component " + std::to_string(AdsorbateComponent); //start when reading "Component 0" for example
  std::string terminate_string="Component " + std::to_string(Terminatecomponent);     //terminate when reading "Component 1", if we are interested in Component 0
  //first get the line number of the destinated component
  while (std::getline(file, str))
  {
    if(str.find(start_string, 0) != std::string::npos){break;}
    start_counter++;
  }
  //printf("%s starts at line number %zu\n", start_string.c_str(), start_counter); 
  file.clear();
  file.seekg(0); 
  std::string MolName;

  bool AlreadyProcessedMLPseudoAtoms = false;
  
  while (std::getline(file, str))
  {
    if(str.find(terminate_string, 0) != std::string::npos)
    {
      //printf("Found terminate string [%s]\n", terminate_string.c_str());
      break;
    }
    if(counter >= start_counter) //start reading after touching the starting line number
    {
      if (str.find(start_string, 0) != std::string::npos) // get the molecule name
      {
        Split_Tab_Space(termsScannedLined, str);
        MolName = termsScannedLined[3];
        SystemComponents.MoleculeName.push_back(MolName);
        std::cout << "-------------- READING Adsorbate" << start_string << " (" << MolName << ")" << " --------------\n";
        MoleculeDefinitionParser(Mol, SystemComponents, termsScannedLined[3], PseudoAtom, Allocate_space);
      }
      if (str.find("IdealGasRosenbluthWeight", 0) != std::string::npos)
      {
        Split_Tab_Space(termsScannedLined, str);
        idealrosen = std::stod(termsScannedLined[1]); //printf("Ideal Chain Rosenbluth Weight: %.5f\n", idealrosen);
      }
      if (str.find("TranslationProbability", 0) != std::string::npos)
      {
        Split_Tab_Space(termsScannedLined, str);
        MoveStats.TranslationProb=std::stod(termsScannedLined[1]);
        ////printf("TranslationProb: %.5f, TotalProb: %.5f\n", TranslationProb, TotalProb);
      }
      if (str.find("RotationProbability", 0) != std::string::npos)
      {
        Split_Tab_Space(termsScannedLined, str);
        MoveStats.RotationProb=std::stod(termsScannedLined[1]);
        ////printf("RotationProb: %.5f, TotalProb: %.5f\n", RotationProb, TotalProb);
      }
      if (str.find("WidomProbability", 0) != std::string::npos)
      {
        Split_Tab_Space(termsScannedLined, str);
        MoveStats.WidomProb=std::stod(termsScannedLined[1]);
        ////printf("WidomProb: %.5f, TotalProb: %.5f\n", WidomProb, TotalProb);
      }
      if (str.find("ReinsertionProbability", 0) != std::string::npos)
      {
        Split_Tab_Space(termsScannedLined, str);
        MoveStats.ReinsertionProb=std::stod(termsScannedLined[1]);
        ////printf("ReinsertionProb: %.5f, TotalProb: %.5f\n", ReinsertionProb, TotalProb);
      }
      if (str.find("IdentityChangeProbability", 0) != std::string::npos)
      {
        Split_Tab_Space(termsScannedLined, str);
        MoveStats.IdentitySwapProb=std::stod(termsScannedLined[1]);
        ////printf("IdentityChangeProb: %.5f, TotalProb: %.5f\n", IdentitySwapProb, TotalProb);
      }
      if (str.find("SwapProbability", 0) != std::string::npos)
      {
        Split_Tab_Space(termsScannedLined, str);
        MoveStats.SwapProb=std::stod(termsScannedLined[1]);
        ////printf("SwapProb: %.5f, TotalProb: %.5f\n", SwapProb, TotalProb);
      }
      if (str.find("GibbsParticleXferProbability", 0) != std::string::npos)
      {
        Split_Tab_Space(termsScannedLined, str);
        MoveStats.GibbsSwapProb=std::stod(termsScannedLined[1]);
        ////printf("GibbsParticleXferProbability: %.5f, TotalProb: %.5f\n", MoveStats.GibbsSwapProb, TotalProb);
      }

      if (str.find("CBCFProbability", 0) != std::string::npos)
      {
        Split_Tab_Space(termsScannedLined, str);
        MoveStats.CBCFProb=std::stod(termsScannedLined[1]);
        temp_hasfracmol=true;
        ////printf("CBCFProb: %.5f, TotalProb: %.5f\n", CBCFProb, TotalProb);
      }
      //Zhao's note: If using CBCF Move, choose the lambda type//
      if (MoveStats.CBCFProb > 0.0)
      {
        if (str.find("LambdaType", 0) != std::string::npos)
        {
          Split_Tab_Space(termsScannedLined, str);
          if(caseInSensStringCompare(termsScannedLined[1], "ShiMaginn"))
          {
            LambdaType = SHI_MAGINN;
          }
          else if(caseInSensStringCompare(termsScannedLined[1], "BrickCFC"))
          {
            LambdaType = BRICK_CFC;
          }
        }
      }
      if (str.find("FugacityCoefficient", 0) != std::string::npos)
      {
        Split_Tab_Space(termsScannedLined, str);
        if(caseInSensStringCompare(termsScannedLined[1], "PR-EOS"))
        {
          fugacoeff = -1.0; //Zhao's note: fugacity coefficient cannot be negative! so using negative values as flags for using Peng-Robinson EOS!
        }
        else if (isFloat(termsScannedLined[1])) //If not PR-EOS, then the user must input a value (double) //
          fugacoeff = std::stod(termsScannedLined[1]);
        else
          throw std::runtime_error("Unrecognized FugacityCoefficient Input: " + termsScannedLined[1] + ", Currently taking 'PR-EOS' or `a float`!!!");
      }
      if (str.find("MolFraction", 0) != std::string::npos)
      {
         Split_Tab_Space(termsScannedLined, str);
         Molfrac = std::stod(termsScannedLined[1]);
      }
      if (str.find("RunTMMC", 0) != std::string::npos)
      {
        Split_Tab_Space(termsScannedLined, str);
        if(caseInSensStringCompare(termsScannedLined[1], "yes"))
        {
          temp_tmmc.DoTMMC = true;
          //printf("TMMC: Running TMMC simulation\n");
        }
      }
      if (str.find("TURN_OFF_CBMC_SWAP", 0) != std::string::npos)
      {
        Split_Tab_Space(termsScannedLined, str);
        if(caseInSensStringCompare(termsScannedLined[1], "yes"))
        {
          SystemComponents.SingleSwap = true;
          //printf("SWAP WITH NO CBMC!\n");
        }
      }
      if(temp_tmmc.DoTMMC)
      {
        if (str.find("TMMCMin", 0) != std::string::npos)
        {
          Split_Tab_Space(termsScannedLined, str);
          sscanf(termsScannedLined[1].c_str(), "%zu", &temp_tmmc.MinMacrostate);
        }
        else if (str.find("TMMCMax", 0) != std::string::npos)
        {
          Split_Tab_Space(termsScannedLined, str);
          sscanf(termsScannedLined[1].c_str(), "%zu", &temp_tmmc.MaxMacrostate);
        }
        if (str.find("UseBiasOnMacrostate", 0) != std::string::npos)
        {
          Split_Tab_Space(termsScannedLined, str);
          if(caseInSensStringCompare(termsScannedLined[1], "yes"))
          {
            temp_tmmc.DoUseBias = true;
            //printf("TMMC: Biasing Insertion/Deletions\n");
          }
        }
        if (str.find("RestartTMMC", 0) != std::string::npos)
        {
          Split_Tab_Space(termsScannedLined, str);
          if(caseInSensStringCompare(termsScannedLined[1], "yes"))
          {
            temp_tmmc.TMMCRestart = true;
            //printf("Need to Read TMMC stats from TMMC_Initial Folder\n");
            //printf("WARNING: Using TMMCRestart, we thus recommend that you do NOT use Initialization steps, start by Equilibration/Production instead!!!");
          }
        }
        if (str.find("UpdateTMMCEvery", 0) != std::string::npos)
        {
          Split_Tab_Space(termsScannedLined, str);
          sscanf(termsScannedLined[1].c_str(), "%zu", &temp_tmmc.UpdateTMEvery);
          //printf("TMMC: Bias Updated from component %s every %zu MC Step(s)!\n", MolName.c_str(), temp_tmmc.UpdateTMEvery);
        }
      }

      if(SystemComponents.UseDNNforHostGuest)
      {
       //Zhao's note: Read the types of atoms that needs to be considered for the DNN Model//
        if (str.find("DNNPseudoAtoms", 0) != std::string::npos)
        {
          if(AlreadyProcessedMLPseudoAtoms)
            throw std::runtime_error("You are reading ***DNNPseudoAtoms*** twice for one component! This is not allowed. Use only one line!!!!\n");
          AlreadyProcessedMLPseudoAtoms = true;
          //This should be a general DNN setup//
          //Consider the positions of water TIP4P, 
          //you don't need to feed the position of the fictional charge to the model//
          //Since it will be re-declared later as managed mem, 
          //we just do a simple CPU mem here//
          //Initialize//
          SystemComponents.ConsiderThisAdsorbateAtom = (bool*) malloc(SystemComponents.Moleculesize[SystemComponents.NComponents.y] * sizeof(bool));
          for(size_t b = 0; b < Mol.Molsize; b++)
            SystemComponents.ConsiderThisAdsorbateAtom[b] = false;

          Split_Tab_Space(termsScannedLined, str);
          for(size_t a = 1; a < termsScannedLined.size(); a++)
          {
            std::string AtomName = termsScannedLined[a];
            size_t Type = get_type_from_name(AtomName, PseudoAtom.Name);
            //Find location of the pseudo atoms
            for(size_t b = 0; b < Mol.Molsize; b++)
            {

              if(Type == Mol.Type[b])
              {
                SystemComponents.ConsiderThisAdsorbateAtom[b] = true;
                //printf("AtomName: %s, Type: %zu, b: %zu, Consider? %s\n", AtomName.c_str(), Type, b, SystemComponents.ConsiderThisAdsorbateAtom[b] ? "true" : "false");
              }
            }
          }
        }
      }

      if (str.find("CreateNumberOfMolecules", 0) != std::string::npos) // Number of Molecules to create
      {
        Split_Tab_Space(termsScannedLined, str); 
        sscanf(termsScannedLined[1 + BoxIndex].c_str(), "%zu", &CreateMolecule);
      }
    }
    counter++;
  }
  MoveStats.GibbsVolumeMoveProb = Vars.GibbsStatistics.GibbsBoxProb;
  if(MoveStats.GibbsVolumeMoveProb > 0.0 || MoveStats.GibbsSwapProb > 0.0)
    Vars.GibbsStatistics.DoGibbs = true;
  MoveStats.VolumeMoveProb = Vars.TempComponents.VolumeMoveProbability;

  MoveStats.NormalizeProbabilities();
  MoveStats.PrintProbabilities();

  //Zhao's note: if monatomic molecule has rotation prob, break the program!//
  size_t currentCompsize = SystemComponents.Moleculesize.size();
  //printf("Current processed %zu components\n", currentCompsize);
  if(SystemComponents.Moleculesize[currentCompsize-1] == 1 && (MoveStats.RotationProb - MoveStats.TranslationProb) > 1e-10)
  {
    throw std::runtime_error("Molecule [" + SystemComponents.MoleculeName[currentCompsize-1] + "] is MONATOMIC, CANNOT DO ROTATION!\n");
  }

  if(SystemComponents.Moleculesize[currentCompsize-1] > 1 && (MoveStats.RotationProb - MoveStats.TranslationProb) <= 1e-10)
  {
    //printf("Molecule [%s] is POLYATOMIC, WE RECOMMEND USING ROTATION (Currently not using it)!\n", SystemComponents.MoleculeName[currentCompsize-1].c_str());
  }
    
  SystemComponents.NumberOfMolecule_for_Component.push_back(0); // Zhao's note: Molecules are created later in main.cpp //
  SystemComponents.Allocate_size.push_back(Allocate_space);
  if(idealrosen < 1e-150) throw std::runtime_error("Error for component {" + std::to_string(AdsorbateComponent) + "}" + "("+ MolName + "): Ideal-Rosenbluth weight not assigned (or not valid), bad. If rigid, assign 1.");
  SystemComponents.IdealRosenbluthWeight.push_back(idealrosen);
  /*
  //Zhao's note: for fugacity coefficient, if assigned in input as "PR-EOS", do Peng-Robinson
  if(fugacoeff < 1e-150)
  {
    throw std::runtime_error("Peng-rob EOS not implemented yet...");
  }
  */
   
  SystemComponents.FugacityCoeff.push_back(fugacoeff);
  //Zhao's note: for now, Molfraction = 1.0
  SystemComponents.MolFraction.push_back(Molfrac);
  SystemComponents.hasfractionalMolecule.push_back(temp_hasfracmol);
  
  LAMBDA lambda;
  //Zhao's note: for bin number = 10, there are 11 bins, the first is when lambda = 0.0, last for lambda = 1.0//
  if(temp_hasfracmol) //Prepare lambda struct if using CBCF//
  {
    lambda.newBin    = 0;
    lambda.delta     = 1.0/static_cast<double>(lambda.binsize); 
    lambda.WangLandauScalingFactor = 0.0; lambda.FractionalMoleculeID = 0;
    lambda.lambdatype = LambdaType;
    lambda.Histogram.resize(lambda.binsize + 1); lambda.biasFactor.resize(lambda.binsize + 1);
    std::fill(lambda.Histogram.begin(),  lambda.Histogram.end(),  0.0);
    std::fill(lambda.biasFactor.begin(), lambda.biasFactor.end(), 0.0);
  }

  //Zhao's note: If we are using CBCFC + TMMC, turn off Normal Swap moves//
  if(temp_hasfracmol && temp_tmmc.DoTMMC)
  {
    if((MoveStats.SwapProb - MoveStats.CBCFProb) > 1e-10)
      throw std::runtime_error("CBCFC + TMMC: Cannot use normal (non-CFC) swap moves when you have CBCFC + TMMC!!!!");
  }
 
  if(temp_tmmc.DoTMMC) //Prepare tmmc struct if using TMMC//
  { 
    if(temp_tmmc.MaxMacrostate < temp_tmmc.MinMacrostate)
    {
      throw std::runtime_error("TMMC: Bad Min/Max Macrostates for TMMC, Min has to be SMALLER THAN OR EQUAL TO Max.");
    }
    temp_tmmc.nbinPerMacrostate = 1;
    if(temp_hasfracmol) temp_tmmc.nbinPerMacrostate = lambda.Histogram.size(); 
    size_t NBin = temp_tmmc.nbinPerMacrostate * (temp_tmmc.MaxMacrostate - temp_tmmc.MinMacrostate + 1);

    temp_tmmc.CMatrix.resize(NBin);
    temp_tmmc.WLBias.resize(NBin);
    temp_tmmc.TMBias.resize(NBin);
    std::fill(temp_tmmc.TMBias.begin(), temp_tmmc.TMBias.end(), 1.0); //Initialize the TMBias//
    temp_tmmc.ln_g.resize(NBin);
    temp_tmmc.lnpi.resize(NBin);
    temp_tmmc.forward_lnpi.resize(NBin);
    temp_tmmc.reverse_lnpi.resize(NBin);
    temp_tmmc.Histogram.resize(NBin);
    //Zhao's note: one must do TMMC restart and RestartInitial together, and skip the initializations//
    if(temp_tmmc.TMMCRestart)
    {
      Read_TMMC_Initial(temp_tmmc, MolName);
      if(!SystemComponents.ReadRestart)
      throw std::runtime_error("If TMMC Restart is used, then RestartFile must be used at the same time (cannot start with no molecules)!");
    }
    
    //Zhao's note: if we set the bounds for min/max macrostate, the number of createMolecule should fall in the range//
    if(temp_tmmc.RejectOutofBound && !SystemComponents.ReadRestart)
      if(CreateMolecule < temp_tmmc.MinMacrostate || CreateMolecule > temp_tmmc.MaxMacrostate)
        throw std::runtime_error("TMMC: Number of created molecule fall out of the TMMC Macrostate range!");

    //printf("TMMC Bias for **Adsorbate** component %zu (%s) is updated every %zu step(s)!\n", AdsorbateComponent, MolName.c_str(), temp_tmmc.UpdateTMEvery);
  }

  SystemComponents.Lambda.push_back(lambda);
  SystemComponents.Tmmc.push_back(temp_tmmc);
  SystemComponents.NumberOfCreateMolecules.push_back(CreateMolecule);
  if(BoxIndex > 0) SystemComponents.NumberOfCreateMolecules[1] = CreateMolecule;
  //printf("CREATE MOL: box %zu, SystemComponents.NumberOfCreateMolecules size: %zu, creating %zu molecules\n", BoxIndex, SystemComponents.NumberOfCreateMolecules.size(), SystemComponents.NumberOfCreateMolecules[1]);

  //Finally, check if all values in SystemComponents are set properly//
  Check_Component_size(SystemComponents);
  //Initialize single values for Mol//
  Mol.size = 0;
  std::cout << "-------------- END OF READING " << start_string << " (" << MolName << ")" << " --------------\n";
  file.close();
}

//Sort vector A, and return a reorder lambda function//
auto generate_reorder_lambda(std::vector<int>& indices, std::vector<int>& A)
{
  std::sort(indices.begin(), indices.end(),
            [&](int i1, int i2) { return A[i1] < A[i2]; });

  // Helper lambda to reorder any vector according to indices
  auto reorder = [&](auto& vec) 
  {
      auto temp = vec;
      for (size_t i = 0; i < indices.size(); ++i)
          vec[i] = temp[indices[i]];
  };
  return reorder;
}

//Filter based on equal
//Has to be integer
auto generate_filter_lambda(std::vector<int>& A, int Criterion, std::vector<int>& filteredIndices)
{
  // Collect indices of elements in A that are > 0
  //std::vector<int> filteredIndices;
  for (int i = 0; i < A.size(); ++i) 
  {
    if (A[i] == Criterion)
    {
      filteredIndices.push_back(i);
    }
  }
  // Helper function to create filtered vector based on indices
  auto createFilteredVector = [&](auto& original)
  {
    typename std::remove_reference<decltype(original)>::type filtered;
    for (size_t i = 0; i < filteredIndices.size(); i++) 
    {
      size_t index = filteredIndices[i];
      //std::cout << "filtered index: " << index << '\n';
      filtered.push_back(original[index]);
    }
    return filtered;
  };
  return createFilteredVector;
}

void ReadRestartInputFileType(Components& SystemComponents)
{
  std::vector<std::string> termsScannedLined{};
  std::string str;
  std::ifstream file("simulation.input");
  while (std::getline(file, str))
  {
    if (str.find("RestartInputFileType", 0) != std::string::npos)
    {
      Split_Tab_Space(termsScannedLined, str);
      if(caseInSensStringCompare(termsScannedLined[1], "RASPA"))
      {
        SystemComponents.RestartInputFileType = RASPA_RESTART;
        //printf("Using RASPA RESTART FILE as RESTART INPUT STRUCTURE\n");
        break;
      }
      else if(caseInSensStringCompare(termsScannedLined[1], "LAMMPS"))
      {
        SystemComponents.RestartInputFileType = LAMMPS_DATA;
        //printf("Using LAMMPS DATA FILE as RESTART INPUT STRUCTURE\n");
        break;
      }
      else
      {
        throw std::runtime_error("Unrecognized input file type for restart file! gRASPA accepts 1. RASPA-2 type restart file (keyword: RASPA); 2. LAMMPS DATA file (keyword: LAMMPS). Plase Check!!!");
      }
    }
  }
  file.close();
}

static inline size_t ReadLMPDataStartComponent(Components& SystemComponents)
{
  std::vector<std::string> termsScannedLined{};
  std::string str;
  std::ifstream file("simulation.input");
  size_t StartComponent = 0;
  while (std::getline(file, str))
  {
    if (str.find("LMPData_Comp_to_Start_with", 0) != std::string::npos)
    {
      Split_Tab_Space(termsScannedLined, str);
      int TempComponent = std::stoi(termsScannedLined[1]);
      if(TempComponent >= SystemComponents.NComponents.x || TempComponent < 0)
      {
        throw std::runtime_error("Your specified Starting Component to read in LMP Data file exceeds the limit of components (0 or total number of components)! CHECK!\n");
      }
      StartComponent = static_cast<size_t>(TempComponent);
    }
    if (str.find("Read_Boxsize", 0) != std::string::npos)
    {
      Split_Tab_Space(termsScannedLined, str);
      if(caseInSensStringCompare(termsScannedLined[1], "yes"))
      {
        SystemComponents.Read_BoxsizeRestart = true;
      }
    }
  }
  file.clear();
  file.seekg(0);
  file.close();
  return StartComponent;
}

//Zhao's note: 
//Tasks: 
//1. Determine the molecule size, try to match the component //
//2. keyword about whether to read box sizes//
//3. keyword about which component to start with (sometimes we can skip the framework atoms//
void LMPDataFileParser(Boxsize& Box, Components& SystemComponents)
{
  //printf("/*----------- READING INITIAL LAMMPS DATA FILE AS INPUT!!! -----------*/\n");
  //printf("/*WARNING: Please make sure that your ATOM TYPES ARE CONSISTENT in Lammps data file and pseudo_atoms.def!!!*/\n");
  //Read the starting component to read from the LMPData file, if zero, then the framework 0 component is read, if one, then you start with one and skip zero//
  size_t StartComponent = ReadLMPDataStartComponent(SystemComponents);

  std::string scannedLine; std::string str;
  std::vector<std::string> termsScannedLined{};
  //Determine framework name (file name)//
  std::string Filename = "LMPDataInitial/System_0/init.data";
  std::ifstream file(Filename);
  std::filesystem::path pathfile = std::filesystem::path(Filename);
  if (!std::filesystem::exists(pathfile))
  {
    throw std::runtime_error("LMPData Initial file not found... You need a folder/file called: *LMPDataInitial/System_0/init.data*\n");
  }
  //Check unique components//
  std::vector<double3>pos;
  std::vector<int>type;
  std::vector<double>charge;
  std::vector<int>MolID;
  std::vector<int>CompID; //Component ID, shown by component name (after # sign) for each atom//
  std::vector<int>AtomID; //number appears as the first index in the line for each atom//
  std::vector<int>ID;     //number appears as it is in the data file//
  //Read the position/type/charge first//
  //Assume the atom types in LMP data file MATCH the pseudo-atom types in gRASPA//
  //Note that lammps index starts from 1//
  size_t counter = 0; size_t atom_count = 0;
  size_t total_atoms = 0;
  size_t start_read_line = 0;
  if(SystemComponents.Read_BoxsizeRestart)
  {
    //printf("Reading boxsizes from LMPData Initial file.\n");
    //printf("***WARNING***All box sizes are shifted to start from zero!\n");
    bool foundx = false; bool foundy = false; bool foundz = false; bool found_xy_xz_yz = false;
    bool reached_atom_types = false;
    while (std::getline(file, str))
    {
      if((str.find("xlo", 0) != std::string::npos) && (str.find("xhi", 0) != std::string::npos))
      {
        Split_Tab_Space(termsScannedLined, str);
        Box.Cell[0] = std::stod(termsScannedLined[1]) - std::stod(termsScannedLined[0]);
        foundx = true;
      }
      if((str.find("ylo", 0) != std::string::npos) && (str.find("yhi", 0) != std::string::npos))
      {
        Split_Tab_Space(termsScannedLined, str);
        Box.Cell[4] = std::stod(termsScannedLined[1]) - std::stod(termsScannedLined[0]);
        foundy = true;
      }     
      if((str.find("zlo", 0) != std::string::npos) && (str.find("zhi", 0) != std::string::npos))
      {
        Split_Tab_Space(termsScannedLined, str);
        Box.Cell[8] = std::stod(termsScannedLined[1]) - std::stod(termsScannedLined[0]);
        foundz = true;
      }
      if((str.find("xy", 0) != std::string::npos) && (str.find("xz", 0) != std::string::npos) && (str.find("yz", 0) != std::string::npos))
      {
        Split_Tab_Space(termsScannedLined, str);
        Box.Cell[3] = std::stod(termsScannedLined[0]);
        Box.Cell[6] = std::stod(termsScannedLined[1]);
        Box.Cell[7] = std::stod(termsScannedLined[2]);
        found_xy_xz_yz = true;
      }
      if((str.find("atom types", 0) != std::string::npos))
      {
        reached_atom_types = true; break;
      }
    }
    if(!(foundx && foundy && foundz) && reached_atom_types) throw std::runtime_error("Error: didn't find box size region in LMPDATA Initial file! CHECK!!!\n");
    Box.Cell[1] = 0.0; Box.Cell[2] = 0.0; Box.Cell[5] = 0.0;
    if(!found_xy_xz_yz)
    {
      Box.Cell[3] = 0.0; Box.Cell[6] = 0.0; Box.Cell[7] = 0.0;
    }
    inverse_matrix(Box.Cell, &Box.InverseCell);
    Box.Volume = matrix_determinant(Box.Cell);
    //DETERMINE Whether Box is cubic/cuboid or not//
    Box.Cubic = true; // Start with cubic box shape, if any value (non-diagonal) is greater than 0, then set to false //
    if((fabs(Box.Cell[3]) + fabs(Box.Cell[6]) + fabs(Box.Cell[7])) > 1e-10) Box.Cubic = false;
    if(Box.Cubic)  //printf("The Simulation Box is Cubic\n");
    if(!Box.Cubic) //printf("The Simulation Box is NOT Cubic\n");
    file.clear();
    file.seekg(0);
  }
  while (std::getline(file, str))
  {
    if (counter < 3 && str.find("toms", 0) != std::string::npos)
    {
      Split_Tab_Space(termsScannedLined, str);
      total_atoms = static_cast<size_t>(std::stoi(termsScannedLined[0]));
      //printf("Found Number of atoms: %zu\n", total_atoms);
    }
    if (str.find("Atoms", 0) != std::string::npos)
    {
      start_read_line = counter + 2; 
      //printf("Start from line %zu\n", start_read_line);
    }
    else if(start_read_line > 0 && counter >= start_read_line && (atom_count) < total_atoms)
    {
      Split_Tab_Space(termsScannedLined, str);
      //if(termsScannedLined.size() == 0) continue;
      //Zhao's note: you need to provide component names for each atom to be correctly read by gRASPA//
      //So We will check for length of the separated vector, if not equal, quit//
      if(termsScannedLined.size() != 10)
        throw std::runtime_error("Check the atom information at line " + std::to_string(counter + 1) + "(start from 1) ! You need 10 elements in total. \nHash symbol is counted as 1 element, after the hash, you need COMPONENT information and ATOM TYPE information!!! Abort!\n");
      AtomID.push_back(std::stoi(termsScannedLined[0]) - 1);
      MolID.push_back(std::stoi(termsScannedLined[1]) - 1);
      type.push_back(std::stoi(termsScannedLined[2]) - 1);
      charge.push_back(std::stod(termsScannedLined[3]));
      pos.push_back({std::stod(termsScannedLined[4]), std::stod(termsScannedLined[5]), std::stod(termsScannedLined[6])});
      ID.push_back(atom_count);
      CompID.push_back(SystemComponents.MatchMoleculeNameToComponentID(termsScannedLined[8]));
      atom_count ++;
    }
    counter ++;
  }
  //printf("Read %zu atoms from LAMMPS\n", pos.size());
  //Try first sort with atomid//
  auto reorder = generate_reorder_lambda(ID, AtomID);
  reorder(AtomID);
  reorder(MolID);
  reorder(type);
  reorder(charge);
  reorder(pos);
  //for(size_t i = 0; i < 3; i++)
  //  //printf("index: %zu, atomID: %d, pos.x: %.5f\n", i, AtomID[i], pos[i].x);
  //Next step: loop over each component, extract info, fill SystemComponents.HostSystem variable//
  for(size_t comp = StartComponent; comp < SystemComponents.NComponents.x; comp++)
  {
    //Declare an empty filteredIndices vector, which will be passed by ref in the next function//
    //This filteredIndices vector is required by the lambda function//
    std::vector<int>filteredIndices;
    auto filter = generate_filter_lambda(CompID, comp, filteredIndices);
    auto filtered_MolID = filter(MolID);
    auto filtered_type  = filter(type);
    auto filtered_charge= filter(charge);
    auto filtered_pos   = filter(pos);
    std::vector<double>filtered_scale(filtered_MolID.size(), 1.0);
    //for(size_t i = 0; i < 3; i++)
      //printf("component: %zu, index: %zu, atomID: %d\n", comp, i, filtered_MolID[i]);
    //Change the MolID, MolID for each component, all starts from zero!//
    size_t NAtomPerComponent = filtered_MolID.size();
    if(NAtomPerComponent % SystemComponents.Moleculesize[comp] != 0) throw std::runtime_error("Check your data file for component " + std::to_string(comp) + ", it is not divisible!!!");
    //size_t NMolPerComponent = NAtomPerComponent / SystemComponent.Molsize[comp];
    //Write to SystemComponents.HostSystem
    double3 firstAtomPos;
    for(size_t i = 0; i < filtered_pos.size(); i++)
    {
      //Component 0 is always the framework (or a part that is not moving, since not moving, no need to wrap according to PBC)
      //Starting from component 1 (this can sometimes be part of framework), it has to be subjected to wrapping of PBC
      //Consider the distance between it and the first atom//
      //wrap the first atom inside the box, then calculate the distance, if greater than 1 box, wrap//
      if(i % SystemComponents.Moleculesize[comp] == 0)
      {
        WrapInBox(filtered_pos[i], Box.Cell, Box.InverseCell, Box.Cubic);
        firstAtomPos = filtered_pos[i];
      }
      else
      {
        double3 dist = filtered_pos[i] - firstAtomPos;
        PBC(dist, Box.Cell, Box.InverseCell, Box.Cubic);
        filtered_pos[i] = dist + firstAtomPos;
      }
      SystemComponents.HostSystem[comp].pos[i] = filtered_pos[i];
      SystemComponents.HostSystem[comp].Type[i] = filtered_type[i];
      SystemComponents.HostSystem[comp].charge[i] = filtered_charge[i];
      SystemComponents.HostSystem[comp].MolID[i] = i / SystemComponents.Moleculesize[comp];
      SystemComponents.HostSystem[comp].scale[i] = filtered_scale[i];
      SystemComponents.HostSystem[comp].scaleCoul[i] = filtered_scale[i];
    }
    bool PreviouslyNotWritten = false; //Zhao's note: a flag for distinguishing whether the component is already recorded
    //For example, the framework component is already taken care for, no need to for example, UpdatePseudoAtoms for this component//
    if(SystemComponents.NumberOfMolecule_for_Component[comp] == 0) PreviouslyNotWritten = true;
    if(PreviouslyNotWritten)
    {
      SystemComponents.NumberOfMolecule_for_Component[comp] = NAtomPerComponent / SystemComponents.Moleculesize[comp];
      SystemComponents.TotalNumberOfMolecules += SystemComponents.NumberOfMolecule_for_Component[comp];
      for(size_t Nmol = 0; Nmol < SystemComponents.NumberOfMolecule_for_Component[comp]; Nmol++)
        SystemComponents.UpdatePseudoAtoms(INSERTION, comp);
      SystemComponents.HostSystem[comp].size = filtered_pos.size();
      SystemComponents.HostSystem[comp].Molsize = SystemComponents.Moleculesize[comp];
    }
  }
  //throw std::runtime_error("DONE!!!");
  file.close();
}

void RestartFileParser(Boxsize& Box, Components& SystemComponents)
{
  std::string scannedLine; std::string str;
  std::vector<std::string> termsScannedLined{};
  //Determine framework name (file name)//
  std::string Filename = "RestartInitial/System_0/restartfile";
  std::ifstream file(Filename);
  std::filesystem::path pathfile = std::filesystem::path(Filename);
  if (!std::filesystem::exists(pathfile))
  {
    throw std::runtime_error("RestartInitial file not found\n");
  }
  size_t counter = 0;
  
  //Zhao's note: MolID in our Atoms struct are relative IDs, for one component, they start with zero.
  //Yet, in Restart file, it start with an arbitrary number (equal to number of previous component molecules)
  //Need to substract it off//
  size_t PreviousCompNMol = 0;
  for(size_t i = SystemComponents.NComponents.y; i < SystemComponents.NComponents.x; i++)
  {
    size_t start = 0;
    while (std::getline(file, str))
    {
      //find range of the part we need to read// 
      if (str.find("Components " + std::to_string(i - SystemComponents.NComponents.y), 0) != std::string::npos) 
        start = counter;
      if (str.find("Maximum-rotation-change component " + std::to_string(i - SystemComponents.NComponents.y), 0) != std::string::npos)
      {  
        break; 
      }
    }
    file.clear();
    file.seekg(0);
    //printf("RASPA Restart: Now reading Component Info\n");
    //Start reading//
    while (std::getline(file, str))
    {
      if (str.find("Fractional-molecule-id component " + std::to_string(i - SystemComponents.NComponents.y), 0) != std::string::npos)
      {
        Split_Tab_Space(termsScannedLined, str);
        if(std::stoi(termsScannedLined[3]) == -1) SystemComponents.hasfractionalMolecule[i] = false;
        if(SystemComponents.hasfractionalMolecule[i])
        {
          sscanf(termsScannedLined[3].c_str(), "%zu", &SystemComponents.Lambda[i].FractionalMoleculeID);
          //printf("Fractional Molecule ID: %zu\n", SystemComponents.Lambda[i].FractionalMoleculeID);
        }
      }
      if(SystemComponents.hasfractionalMolecule[i])
      {
        if (str.find("Lambda-factors component " + std::to_string(i - SystemComponents.NComponents.y), 0) != std::string::npos)
        {
          Split_Tab_Space(termsScannedLined, str);
          SystemComponents.Lambda[i].WangLandauScalingFactor = std::stod(termsScannedLined[3]);
          //printf("WL Factor: %.5f\n", SystemComponents.Lambda[i].WangLandauScalingFactor);
        }
        if (str.find("Number-of-biasing-factors component " + std::to_string(i - SystemComponents.NComponents.y), 0) != std::string::npos)
        {
          Split_Tab_Space(termsScannedLined, str);
          sscanf(termsScannedLined[3].c_str(), "%zu", &SystemComponents.Lambda[i].binsize);
          //printf("binsize: %zu\n", SystemComponents.Lambda[i].binsize);
          if(SystemComponents.Lambda[i].binsize != SystemComponents.Lambda[i].Histogram.size()) throw std::runtime_error("CFC Bin size don't match!");
        }
        if (str.find("Biasing-factors component " + std::to_string(i - SystemComponents.NComponents.y), 0) != std::string::npos)
        {
          for(size_t j = 0; j < SystemComponents.Lambda[i].binsize; j++)
          {
            Split_Tab_Space(termsScannedLined, str);
            SystemComponents.Lambda[i].biasFactor[j] = std::stod(termsScannedLined[3 + j]); 
            //printf("Biasing Factor %zu: %.5f\n", j, SystemComponents.Lambda[i].biasFactor[j]);
          }
        }
      }
      //Read translation rotation maxes//
      if (str.find("Maximum-translation-change component " + std::to_string(i - SystemComponents.NumberOfFrameworks), 0) != std::string::npos)
      {
        Split_Tab_Space(termsScannedLined, str);
        SystemComponents.MaxTranslation[i - SystemComponents.NComponents.y].x = std::stod(termsScannedLined[3]);
        SystemComponents.MaxTranslation[i - SystemComponents.NComponents.y].y = std::stod(termsScannedLined[4]);
        SystemComponents.MaxTranslation[i - SystemComponents.NComponents.y].z = std::stod(termsScannedLined[5]);
      }
      if (str.find("Maximum-rotation-change component " + std::to_string(i - SystemComponents.NComponents.y), 0) != std::string::npos)
      {
        Split_Tab_Space(termsScannedLined, str);
        SystemComponents.MaxRotation[i - SystemComponents.NComponents.y].x = std::stod(termsScannedLined[3]);
        SystemComponents.MaxRotation[i - SystemComponents.NComponents.y].y = std::stod(termsScannedLined[4]);
        SystemComponents.MaxRotation[i - SystemComponents.NComponents.y].z = std::stod(termsScannedLined[5]);
        break;
      }
    }
    file.clear(); 
    file.seekg(0);
    //Start reading atom positions and other information//
    start = 0; counter = 0;
    while (std::getline(file, str))
    {
      if (str.find("Component: " + std::to_string(i - SystemComponents.NComponents.y), 0) != std::string::npos)
      {
        Split_Tab_Space(termsScannedLined, str);
        start = counter + 2;
        size_t temp; 
        sscanf(termsScannedLined[3].c_str(), "%zu", &temp);
        //printf("Adsorbate Component %zu, #: %zu\n", i, temp);
        SystemComponents.NumberOfMolecule_for_Component[i] = temp;
        SystemComponents.TotalNumberOfMolecules += SystemComponents.NumberOfMolecule_for_Component[i];
        ////printf("There are %zu Molecules, molsize = %zu, line is %zu\n", SystemComponents.NumberOfMolecule_for_Component[i], SystemComponents.Moleculesize[i], counter);
        break; 
      }
      counter++;
    }

    //Update The Number of pseudoAtoms//
    for(size_t Nmol = 0; Nmol < SystemComponents.NumberOfMolecule_for_Component[i]; Nmol++)
      SystemComponents.UpdatePseudoAtoms(INSERTION, i);

    if(SystemComponents.Tmmc[i].DoTMMC && SystemComponents.Tmmc[i].RejectOutofBound)
      if(SystemComponents.NumberOfMolecule_for_Component[i] < SystemComponents.Tmmc[i].MinMacrostate || SystemComponents.NumberOfMolecule_for_Component[i] > SystemComponents.Tmmc[i].MaxMacrostate)
        throw std::runtime_error("Number of molecule fall out of the TMMC Macrostate range!");

    counter  = 0;
    file.clear();
    file.seekg(0);
    if(SystemComponents.NumberOfMolecule_for_Component[i] == 0) continue;
    size_t interval = SystemComponents.NumberOfMolecule_for_Component[i]* SystemComponents.Moleculesize[i];
    double3 pos[SystemComponents.NumberOfMolecule_for_Component[i]       * SystemComponents.Moleculesize[i]];
    double  scale[SystemComponents.NumberOfMolecule_for_Component[i]     * SystemComponents.Moleculesize[i]];
    double  charge[SystemComponents.NumberOfMolecule_for_Component[i]    * SystemComponents.Moleculesize[i]];
    double  scaleCoul[SystemComponents.NumberOfMolecule_for_Component[i] * SystemComponents.Moleculesize[i]];
    size_t  Type[SystemComponents.NumberOfMolecule_for_Component[i]      * SystemComponents.Moleculesize[i]];
    size_t  MolID[SystemComponents.NumberOfMolecule_for_Component[i]     * SystemComponents.Moleculesize[i]];

    size_t atom=0;
    while (std::getline(file, str))
    {
      //Read positions, Type and MolID//
      //Position: 0, velocity: 1, force: 2, charge: 3, scaling: 4
      if((counter >= start) && (counter < start + interval))
      {
        atom=counter - start;
        if (!(str.find("Adsorbate-atom-position", 0) != std::string::npos)) throw std::runtime_error("Cannot find matching strings in the range for reading positions!");
        Split_Tab_Space(termsScannedLined, str);
        pos[atom] = {std::stod(termsScannedLined[3]), std::stod(termsScannedLined[4]), std::stod(termsScannedLined[5])};
        sscanf(termsScannedLined[1].c_str(), "%zu", &MolID[atom]);
        size_t atomid; sscanf(termsScannedLined[2].c_str(), "%zu", &atomid);
        Type[atom] = SystemComponents.HostSystem[i].Type[atomid]; //for every component, the types of atoms for the first molecule is always there, just copy it//
        ////printf("Reading Positions, atom: %zu, xyz: %.5f %.5f %.5f, Type: %zu, MolID: %zu\n", atom, pos[atom].x, pos[atom].y, pos[atom].z, Type[atom], MolID[atom]);
        //Zhao's note: adjust the MolID from absolute to relative to component//
        MolID[atom] -= PreviousCompNMol;
        //Zhao's note, 051724, Consider adding a wrap/unwrap here//
        //Calculate atom distance between atom A and the first atom in the molecule//
        //I know, many code will consider using the oxygen atom as the reference for wrapping, no, this is not general enough//
        //We will use the first atom//
        double3 first_bead_pos;
        if(atomid == 0)
        {
          first_bead_pos = pos[atom];
        }
        else if(atomid != 0)
        {
          double3 dist_vec = pos[atom] - first_bead_pos;
          PBC(dist_vec, Box.Cell, Box.InverseCell, Box.Cubic);
          pos[atom] = first_bead_pos + dist_vec;
        }
      }
      //Read charge//
      if((counter >= start + interval * 3) && (counter < start + interval * 4))
      { 
        atom = counter - (start + interval * 3);
        Split_Tab_Space(termsScannedLined, str);
        charge[atom] = std::stod(termsScannedLined[3]);
        ////printf("Reading charge, atom: %zu, charge: %.5f\n", atom, charge[atom]);
      }
      //Read scaling and scalingCoul//
      atom=0;
      if((counter >= start + interval * 4) && (counter < start + interval * 5))
      {
        atom = counter - (start + interval * 4);
        Split_Tab_Space(termsScannedLined, str);
        double  lambda     = std::stod(termsScannedLined[3]);
        double2 val; val.x = 1.0; val.y = 1.0;
        if(lambda < 1.0) val = SystemComponents.Lambda[i].SET_SCALE(lambda);
        scale[atom]     = val.x;
        scaleCoul[atom] = val.y;
        //Determine the molecule ID// 
        size_t tempMOLID = 0;
        sscanf(termsScannedLined[1].c_str(), "%zu", &tempMOLID);
        //Determine the currentBin for the fractional molecule//
        if(tempMOLID == SystemComponents.Lambda[i].FractionalMoleculeID)
        { 
          //printf("Lambda from RestartInitial is: %.100f\n", lambda);
          //floor/ceil functions
          //double smallEpsilon = 1e-5; //FOR DEBUGGING NUMERICAL ISSUES IN THE FUTURE//
          size_t currentBIN = static_cast<size_t>(lambda/SystemComponents.Lambda[i].delta);
          //Zhao's note: do a casting test (0.7 is actually 0.69999, when using static cast, when delta is 0.1, the bin is 6 (should be 7)//
          double doubBin = static_cast<double>(currentBIN) * SystemComponents.Lambda[i].delta;
          if(abs(doubBin - lambda) > 1e-3)
          {
            //printf("static cast round off too much!");
            if((doubBin - lambda) > 1e-3)  currentBIN--;
            if(-(doubBin - lambda) > 1e-3) currentBIN++;
          }
          //printf("CURRENT BIN is %zu\n", currentBIN);
          SystemComponents.Lambda[i].currentBin = currentBIN;
        }
        ////printf("Reading scaling, atom: %zu, scale: %.5f %.5f\n", atom, scale[atom], scaleCoul[atom]);
      }
      counter ++;
    }
    for(size_t j = 0; j < interval; j++)
    {
      SystemComponents.HostSystem[i].pos[j] = pos[j]; 
      SystemComponents.HostSystem[i].charge[j] = charge[j]; 
      SystemComponents.HostSystem[i].scale[j] = scale[j]; SystemComponents.HostSystem[i].scaleCoul[j] = scaleCoul[j]; 
      SystemComponents.HostSystem[i].Type[j] = Type[j]; SystemComponents.HostSystem[i].MolID[j] = MolID[j];
      ////printf("Data for %zu: %.5f %.5f %.5f %.5f %.5f %.5f %zu %zu\n", j, Host_System[i].pos[j].x, Host_System[i].pos[j].y, Host_System[i].pos[j].z, Host_System[i].charge[j], Host_System[i].scale[j], Host_System[i].scaleCoul[j], Host_System[i].Type[j], Host_System[i].MolID[j]);
    }
    SystemComponents.HostSystem[i].size = interval;
    SystemComponents.HostSystem[i].Molsize = SystemComponents.Moleculesize[i];

    PreviousCompNMol += SystemComponents.NumberOfMolecule_for_Component[i];
  }
  file.close();
}

// Read ExcessVolume and Helium Void Fraction //
void ReadVoidFraction(Variables& Vars)
{
  std::vector<std::string> termsScannedLined{};
  std::string str;
  std::ifstream file("simulation.input");
  while (std::getline(file, str))
  {
    if (str.find("HeliumVoidFraction", 0) != std::string::npos)
    {
      Split_Tab_Space(termsScannedLined, str);
      Vars.TempComponents.HeliumVoidFraction=std::stod(termsScannedLined[1]);
    }
    if (str.find("ExcessVolume", 0) != std::string::npos)
    {
      Split_Tab_Space(termsScannedLined, str);
      Vars.TempComponents.ExcessVolume=std::stod(termsScannedLined[1]);
    }
  }
  file.close();
}

void ReadDNNModelSetup(Components& SystemComponents)
{
  std::vector<std::string> termsScannedLined{};
  std::string str;
  std::ifstream file("simulation.input");
  while (std::getline(file, str))
  {
    if (str.find("UseDNNforHostGuest", 0) != std::string::npos)
    {
      Split_Tab_Space(termsScannedLined, str);
      if(caseInSensStringCompare(termsScannedLined[1], "yes"))
      {
        SystemComponents.UseDNNforHostGuest = true;
        //printf("Using DNN Model\n");
      }
      break;
    }
  }
  if(!SystemComponents.UseDNNforHostGuest) return;

  file.clear();
  file.seekg(0);

  bool foundMethod = false; bool DNNUnitFound = false;
  while (std::getline(file, str))
  {
    if (str.find("DNNMethod", 0) != std::string::npos)
    {
      Split_Tab_Space(termsScannedLined, str);
      //printf("Found DNNMethod\n");
      if(caseInSensStringCompare(termsScannedLined[1], "Allegro"))
      {
        //printf("Found Allegro\n");
        SystemComponents.UseAllegro = true; foundMethod = true;
      }
      else if(caseInSensStringCompare(termsScannedLined[1], "LCLin"))
      {
        //printf("Found LCLin\n");
        SystemComponents.UseLCLin = true; foundMethod = true;
      }
      else {throw std::runtime_error("CANNOT IDENTIFY DNN MODEL in simulation.input file");}
    }
    if (str.find("DNNEnergyUnit", 0) != std::string::npos)
    {
      Split_Tab_Space(termsScannedLined, str);
      if(caseInSensStringCompare(termsScannedLined[1], "kJ_mol"))
      {
        SystemComponents.DNNEnergyConversion = 100.0; //from kJ_mol to 10J_mol
        DNNUnitFound = true; //printf("DNN Model is using kJ/mol as the Energy Unit\n");
      }
      else if(caseInSensStringCompare(termsScannedLined[1], "eV"))
      {
        SystemComponents.DNNEnergyConversion = 9648.53074992579265; //from eV to 10J_mol
        DNNUnitFound = true; //printf("DNN Model is using eV as the Energy Unit\n");
      }
      else
      {throw std::runtime_error("Unknown Energy Unit for DNN Model");}
    }
  }
  if(SystemComponents.UseDNNforHostGuest && !DNNUnitFound)
    throw std::runtime_error("You are using DNN models but there is no ENERGY UNIT specified!!!!");

  if(SystemComponents.UseAllegro && SystemComponents.UseLCLin)
    throw std::runtime_error("CANNOT USE Li-Chiang Lin's and Allegro at the same time!!!!");
  if(!foundMethod)
    throw std::runtime_error("CANNOT FIND the DNNMethod INPUT COMMAND in simulation.input file");
  file.close();
}

//###PATCH_LCLIN_READDATA###//
//###PATCH_ALLEGRO_READDATA###//
