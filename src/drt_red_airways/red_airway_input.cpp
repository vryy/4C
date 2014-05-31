/*!----------------------------------------------------------------------
\file red_airway_input.cpp
\brief

<pre>
Maintainer: Mahmoud Ismail
            ismail@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/

#include "red_airway.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_linedefinition.H"

using namespace DRT::UTILS;


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool DRT::ELEMENTS::RedAirway::ReadElement(const std::string& eletype,
                                           const std::string& distype,
                                           DRT::INPUT::LineDefinition* linedef)
{
  const int ndim = DRT::Problem::Instance()->NDim();
  if (ndim!=3)
     dserror("Problem defined as %dd, but found Reduced dimensional AIRWAY element.",ndim);

  // read number of material model
  int material = 0;
  linedef->ExtractInt("MAT",material);
  SetMaterial(material);


  linedef->ExtractString("TYPE",elemType_);
  if (elemType_ == "Resistive" || elemType_ == "InductoResistive" || elemType_ == "ComplientResistive" || elemType_ == "RLC" || elemType_ == "ViscoElasticRLC" || elemType_ =="ConvectiveViscoElasticRLC")
  {
    linedef->ExtractString("Resistance",resistance_);
    linedef->ExtractString("ElemSolvingType",elemsolvingType_);

    double Ew, tw, A, Ts, Phis, nu,velPow;
    int generation;
    linedef->ExtractDouble("PowerOfVelocityProfile",velPow);
    linedef->ExtractDouble("WallElasticity",Ew);
    linedef->ExtractDouble("PoissonsRatio",nu);
    linedef->ExtractDouble("ViscousTs",Ts);
    linedef->ExtractDouble("ViscousPhaseShift",Phis);
    linedef->ExtractDouble("WallThickness",tw);
    linedef->ExtractDouble("Area",A);
    linedef->ExtractInt("Generation",generation);


    // Correct the velocity profile power
    // this is because the 2.0 is the minimum energy consumtive laminar profile
    if (velPow < 2.0)
      velPow = 2.0;
    elemParams_["PowerOfVelocityProfile"]= velPow;
    elemParams_["WallElasticity"]   = Ew;
    elemParams_["PoissonsRatio"]    = nu;
    elemParams_["WallThickness"]    = tw;
    elemParams_["Area"]             = A;
    elemParams_["ViscousTs"]        = Ts;
    elemParams_["ViscousPhaseShift"]= Phis;
    generation_                     = generation;
    if (linedef->HaveNamed("BranchLength"))
    {
      double l_branch = 0.0;
      linedef->ExtractDouble("BranchLength",l_branch);
      elemParams_["BranchLength"] = l_branch;
    }
    else
    {
      elemParams_["BranchLength"] = -1.0;
    }
  }
  else
  {
    dserror("Reading type of RED_AIRWAY element failed: ComplientResistive/PoiseuilleResistive/TurbulentPoiseuilleResistive/InductoResistive/RLC/ViscoElasticRLC/ConvectiveViscoElasticRLC");
    exit(1);
  }

  return true;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool DRT::ELEMENTS::RedAcinus::ReadElement(const std::string& eletype,
                                           const std::string& distype,
                                           DRT::INPUT::LineDefinition* linedef)
{
  const int ndim = DRT::Problem::Instance()->NDim();
  if (ndim!=3)
    dserror("Problem defined as %dd, but found Reduced dimensional ACINUS element.",ndim);

  // read number of material model
  int material = 0;
  linedef->ExtractInt("MAT",material);
  SetMaterial(material);


  linedef->ExtractString("TYPE",elemType_);
  if (elemType_ == "Exponential" || elemType_ == "DoubleExponential" || elemType_ == "NeoHookean"|| elemType_ == "VolumetricOgden")
  {
    double acinusVol, alveolarDuctVol, A;
    const int generation = -1;
    linedef->ExtractDouble("AcinusVolume",acinusVol);
    linedef->ExtractDouble("AlveolarDuctVolume",alveolarDuctVol);
    linedef->ExtractDouble("Area",A);

    elemParams_["AcinusVolume"]       = acinusVol;
    elemParams_["AlveolarDuctVolume"] = alveolarDuctVol;
    elemParams_["Area"]               = A;
    generation_                       = generation;

  }
  else
  {
    dserror("Reading type of Acinare element failed: Exponential/DoubleExponential/NeoHookean/VolumetricOgden");
    exit(1);
  }

  if(elemType_ == "Exponential")
  {
    double E1_0, E1_lin, E1_exp, Tau;
    linedef->ExtractDouble("E1_0" ,E1_0);
    linedef->ExtractDouble("E1_LIN",E1_lin);
    linedef->ExtractDouble("E1_EXP",E1_exp);
    linedef->ExtractDouble("TAU" , Tau);

    elemParams_["E1_0"  ] = E1_0;
    elemParams_["E1_LIN"] = E1_lin;
    elemParams_["E1_EXP"] = E1_exp;
    elemParams_["TAU"   ] = Tau;
  }
  else if(elemType_ == "DoubleExponential")
  {
    double E1_01, E1_lin1, E1_exp1, Tau1;
    linedef->ExtractDouble("E1_01" ,E1_01);
    linedef->ExtractDouble("E1_LIN1",E1_lin1);
    linedef->ExtractDouble("E1_EXP1",E1_exp1);
    linedef->ExtractDouble("TAU1" , Tau1);
    elemParams_["E1_01"  ] = E1_01;
    elemParams_["E1_LIN1"] = E1_lin1;
    elemParams_["E1_EXP1"] = E1_exp1;
    elemParams_["TAU1"   ] = Tau1;

    double E1_02, E1_lin2, E1_exp2, Tau2;
    linedef->ExtractDouble("E1_02" ,E1_02);
    linedef->ExtractDouble("E1_LIN2",E1_lin2);
    linedef->ExtractDouble("E1_EXP2",E1_exp2);
    linedef->ExtractDouble("TAU2" , Tau2);
    elemParams_["E1_02"  ] = E1_02;
    elemParams_["E1_LIN2"] = E1_lin2;
    elemParams_["E1_EXP2"] = E1_exp2;
    elemParams_["TAU2"   ] = Tau2;
  }
  else if(elemType_ == "VolumetricOgden")
  {
    double kappa, beta;
    linedef->ExtractDouble("KAPPA",kappa);
    linedef->ExtractDouble("BETA",beta);
    elemParams_["kappa"]              = kappa;
    elemParams_["beta"]               = beta;
  }
  else
  {
    // do nothing
  }

  return true;
}



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool DRT::ELEMENTS::RedInterAcinarDep::ReadElement(const std::string& eletype,
                                                   const std::string& distype,
                                                   DRT::INPUT::LineDefinition* linedef)
{
  const int ndim = DRT::Problem::Instance()->NDim();
  if (ndim!=3)
    dserror("Problem defined as %dd, but found Reduced dimensional INTER ACINAR DEPENDENCE element.",ndim);

  // read number of material model
  const int generation = -2;
  generation_                       = generation;


  return true;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool DRT::ELEMENTS::RedAirBloodScatra::ReadElement(const std::string& eletype,
                                                   const std::string& distype,
                                                   DRT::INPUT::LineDefinition* linedef)
{
  const int ndim = DRT::Problem::Instance()->NDim();
  if (ndim!=3)
    dserror("Problem defined as %dd, but found Reduced dimensional INTER ACINAR DEPENDENCE element.",ndim);

  // read number of material model
  const int generation = -2;
  generation_                       = generation;

    double diff = 0.0;
  linedef->ExtractDouble("DiffusionCoefficient",diff);
  elemParams_["DiffusionCoefficient"] = diff;

  double th = 0.0;
  linedef->ExtractDouble("WallThickness",th);
  elemParams_["WallThickness"] = th;

  double percDiffArea = 0.0;
  linedef->ExtractDouble("PercentageOfDiffusionArea",percDiffArea);
  elemParams_["PercentageOfDiffusionArea"] = percDiffArea;

  return true;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool DRT::ELEMENTS::RedAirBloodScatraLine3::ReadElement(const std::string& eletype,
                                                   const std::string& distype,
                                                   DRT::INPUT::LineDefinition* linedef)
{
  const int ndim = DRT::Problem::Instance()->NDim();
  if (ndim!=3)
    dserror("Problem defined as %dd, but found Reduced dimensional INTER ACINAR DEPENDENCE element.",ndim);

  double diff = 0.0;
  linedef->ExtractDouble("DiffusionCoefficient",diff);
  elemParams_["DiffusionCoefficient"] = diff;

  double th = 0.0;
  linedef->ExtractDouble("WallThickness",th);
  elemParams_["WallThickness"] = th;

  double percDiffArea = 0.0;
  linedef->ExtractDouble("PercentageOfDiffusionArea",percDiffArea);
  elemParams_["PercentageOfDiffusionArea"] = percDiffArea;

  return true;
}

