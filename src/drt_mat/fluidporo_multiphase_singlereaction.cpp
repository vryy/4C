/*----------------------------------------------------------------------*/
/*!
 \file fluidporo_multiphase_singlereaction.cpp

 \brief a fluid material for porous multiphase flow defining one reaction (mass sources and sinks)

   \level 3

   \maintainer  Anh-Tu Vuong
                vuong@lnm.mw.tum.de
                http://www.lnm.mw.tum.de
                089 - 289-15251
 *----------------------------------------------------------------------*/

#include "fluidporo_multiphase_singlereaction.H"

#include <vector>
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/standardtypes_cpp.H"
#include "../drt_mat/matpar_bundle.H"

/*----------------------------------------------------------------------*
 *  constructor (public)                               vuong 08/16      |
 *----------------------------------------------------------------------*/
MAT::PAR::FluidPoroSingleReaction::FluidPoroSingleReaction(Teuchos::RCP<MAT::PAR::Material> matdata) :
  Parameter(matdata),
  numscal_(matdata->GetInt("NUMSCAL")),
  numphases_(matdata->GetInt("NUMPHASES")),
  scale_(matdata->Get<std::vector<int> >("SCALE")),
  coupling_(SetCouplingType(matdata)),
  functID_(matdata->GetInt("FUNCTID")),
  isinit_(false),
  scalarnames_(numscal_),
  pressurenames_(numphases_),
  saturationnames_(numphases_),
  porosityname_("porosity")
{
}


/*----------------------------------------------------------------------*
 *  Create Material (public)                             vuong 08/16      |
*----------------------------------------------------------------------*/
void MAT::PAR::FluidPoroSingleReaction::Initialize()
{
  if(not isinit_)
  {
    // safety check
    if(Function(functID_-1).NumberComponents()!=1)
      dserror("expected only one component for single phase reaction!");

    for (int k=0;k<numscal_;k++)
    {
      // add scalar names
      {
        std::ostringstream temp;
        temp << k;
        scalarnames_[k] = "phi"+temp.str();

        if(not Function(functID_-1).IsVariable(0,scalarnames_[k]))
          Function(functID_-1).AddVariable(0,scalarnames_[k],0.0);
      }
    }

    for (int k=0;k<numphases_;k++)
    {
      // add pressure names
      {
        std::ostringstream temp;
        temp << k;
        pressurenames_[k] = "p"+temp.str();

        if(not Function(functID_-1).IsVariable(0,pressurenames_[k]))
          Function(functID_-1).AddVariable(0,pressurenames_[k],0.0);
      }

      // add saturation names
      {
        std::ostringstream temp;
        temp << k;
        saturationnames_[k] = "S"+temp.str();

        if(not Function(functID_-1).IsVariable(0,saturationnames_[k]))
          Function(functID_-1).AddVariable(0,saturationnames_[k],0.0);
      }
    }

    // add porosity
    {
      if(not Function(functID_-1).IsVariable(0,porosityname_))
        Function(functID_-1).AddVariable(0,porosityname_,0.0);
    }

    isinit_=true;
  }
  return;
}

/*----------------------------------------------------------------------*
 *  set values in function                                 vuong 08/16 |
*----------------------------------------------------------------------*/
void MAT::PAR::FluidPoroSingleReaction::EvaluateFunction(
    std::vector<double> &                    reacval,
    std::vector<std::vector<double> >&       reacderivs,
    const std::vector<double>&               pressure,
    const std::vector<double>&               saturation,
    const double&                            porosity,
    const std::vector<double>&               scalar)
{
//  dsassert(pressurenames_.size()==pressure.size(),"Invalid number of pressure values for this the fluid poro reaction material!");
//  dsassert(pressurenames_.size()==saturation.size(),"Invalid number of pressure values for this the fluid poro reaction material!");
//  dsassert(scalarnames_.size()==scalar.size(),"Invalid number of pressure values for this the fluid poro reaction material!");

  if(numphases_!=(int)pressure.size())
    dserror("Invalid number of pressure values for this the fluid poro reaction material!");
  if(numphases_!=(int)saturation.size())
    dserror("Invalid number of saturation values for this the fluid poro reaction material!");
  if(numscal_!=(int)scalar.size())
    dserror("Invalid number of scalar values for this the fluid poro reaction material!");

  std::vector<std::pair<std::string,double> > variables;
  variables.reserve(numphases_+numphases_+1);

  std::vector<std::pair<std::string,double> > constants;
  constants.reserve(numscal_);

  // set pressure values as variable
  for (int k=0;k<numphases_;k++)
    variables.push_back(std::pair<std::string,double>(pressurenames_[k],pressure[k]));

  // set saturation values as variable
  for (int k=0;k<numphases_;k++)
    variables.push_back(std::pair<std::string,double>(saturationnames_[k],saturation[k]));

  // set porosity value as variable
  variables.push_back(std::pair<std::string,double>(porosityname_,porosity));

  // set scalar values as constants
  for (int k=0;k<numscal_;k++)
    constants.push_back(std::pair<std::string,double>(scalarnames_[k],scalar[k]));

  // evaluate the reaction term
  double curval = Function(functID_-1).Evaluate(0,variables,constants);
  // evaluate derivatives
  std::vector<double> curderivs(Function(functID_-1).FctDer(0,variables,constants));
  // sum them up
  for (int k=0;k<numphases_;k++)
  {
    if((*scale_)[k]!=0)
    {
      reacval[k] += (*scale_)[k]*curval;
      for (int j=0;j<numphases_;j++)
        reacderivs[k][j] += (*scale_)[k]*curderivs[j];
    }
  }


  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
inline DRT::UTILS::VariableExprFunction& MAT::PAR::FluidPoroSingleReaction::Function(int functnum) const
{
  try
  {
    DRT::UTILS::VariableExprFunction& funct =
        dynamic_cast<DRT::UTILS::VariableExprFunction&>(DRT::Problem::Instance()->Funct(functnum));
    return funct;
  }
  catch(std::bad_cast & exp)
  {
    dserror("Cast to VarExp Function failed! For phase law definition only 'VAREXPR' functions are allowed!\n"
        "Check your input file!");
    return dynamic_cast<DRT::UTILS::VariableExprFunction&>(DRT::Problem::Instance()->Funct(functnum));
  }
}

/*----------------------------------------------------------------------*
 *  Create Material (public)                             vuong 08/16      |
*----------------------------------------------------------------------*/
Teuchos::RCP<MAT::Material> MAT::PAR::FluidPoroSingleReaction::CreateMaterial()
{
  return Teuchos::rcp(new MAT::FluidPoroSingleReaction(this));
}

/*----------------------------------------------------------------------*
 *  translate coupling type                             vuong 08/16      |
 *----------------------------------------------------------------------*/
MAT::PAR::FluidPoroSingleReaction::porofluid_reaction_coupling
MAT::PAR::FluidPoroSingleReaction::SetCouplingType( Teuchos::RCP<MAT::PAR::Material> matdata )
{
  if ( *(matdata->Get<std::string >("COUPLING")) == "scalar_by_function" )
  {
    return porofluid_reac_coup_scalarsbyfunction;
  }
  else if ( *(matdata->Get<std::string >("COUPLING")) == "no_coupling")
  {
    return porofluid_reac_coup_none;
  }
  else
  {
    return porofluid_reac_coup_none;
  }
}

/*----------------------------------------------------------------------*
  global instance of parameter class                   vuong 08/16     |
*----------------------------------------------------------------------*/
MAT::FluidPoroSingleReactionType MAT::FluidPoroSingleReactionType::instance_;

/*----------------------------------------------------------------------*
 *  Create material from given data                          vuong 08/16 |
 *----------------------------------------------------------------------*/

DRT::ParObject* MAT::FluidPoroSingleReactionType::Create(const std::vector<char> & data)
{
  MAT::FluidPoroSingleReaction* fluid_poro = new MAT::FluidPoroSingleReaction();
  fluid_poro->Unpack(data);
  return fluid_poro;
}

/*----------------------------------------------------------------------*
 *   Create empty material                                  vuong 08/16 |
 *----------------------------------------------------------------------*/
MAT::FluidPoroSingleReaction::FluidPoroSingleReaction() :
  params_(NULL)
{
}

/*----------------------------------------------------------------------*
*   Create material with parameters                         vuong 08/16 |
*----------------------------------------------------------------------*/
MAT::FluidPoroSingleReaction::FluidPoroSingleReaction(MAT::PAR::FluidPoroSingleReaction* params) :
  params_(params)
{
}

/*----------------------------------------------------------------------*
 * pack material for commuication                           vuong 08/16 |
*----------------------------------------------------------------------*/
void MAT::FluidPoroSingleReaction::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);

  // matid
  int matid = -1;
  if (params_ != NULL)
    matid = params_->Id(); // in case we are in post-process mode
  AddtoPack(data, matid);
}

/*----------------------------------------------------------------------*
* unpack material                                           vuong 08/16 |
*----------------------------------------------------------------------*/
void MAT::FluidPoroSingleReaction::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position, data, type);
  if (type != UniqueParObjectId())
    dserror("wrong instance type data");

  // matid
  int matid;
  ExtractfromPack(position,data,matid);
  params_ = NULL;
  if (DRT::Problem::Instance()->Materials() != Teuchos::null)
  if (DRT::Problem::Instance()->Materials()->Num() != 0)
  {
    const int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();
    MAT::PAR::Parameter* mat = DRT::Problem::Instance(probinst)->Materials()->ParameterById(matid);
    if (mat->Type() == MaterialType())
      params_ = static_cast<MAT::PAR::FluidPoroSingleReaction*>(mat);
    else
      dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(), MaterialType());
  }

  if (position != data.size())
  dserror("Mismatch in size of data %d <-> %d",data.size(),position);
}

/*----------------------------------------------------------------------*
 *  initialize                                              vuong 08/16 |
*----------------------------------------------------------------------*/
void MAT::FluidPoroSingleReaction::Initialize()
{
  params_->Initialize();
  return;
}

/*----------------------------------------------------------------------*
 *  set values in function                                 vuong 08/16 |
*----------------------------------------------------------------------*/
void MAT::FluidPoroSingleReaction::EvaluateReaction(
    std::vector<double> &                    reacval,
    std::vector<std::vector<double> >&       reacderivs,
    const std::vector<double>&               pressure,
    const std::vector<double>&               saturation,
    const double&                            porosity,
    const std::vector<double>&               scalar)
{
  params_->EvaluateFunction(reacval,reacderivs,pressure,saturation,porosity,scalar);

  return;
}
