/*----------------------------------------------------------------------*/
/*!
 \file fluidporo_multiphase.cpp

 \brief material for multiphase porous fluid

   \level 3

   \maintainer  Anh-Tu Vuong
                vuong@lnm.mw.tum.de
                http://www.lnm.mw.tum.de
                089 - 289-15251
 *----------------------------------------------------------------------*/




#include <vector>
#include "fluidporo_multiphase.H"
#include "fluidporo_singlephase.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_bundle.H"

/*----------------------------------------------------------------------*
 | constructor of paramter class                            vuong 08/16 |
 *----------------------------------------------------------------------*/
MAT::PAR::FluidPoroMultiPhase::FluidPoroMultiPhase(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: MatList(matdata),
  permeability_(matdata->GetDouble("PERMEABILITY"))
{
}

/*----------------------------------------------------------------------*
 | create a poro multiphase material                        vuong 08/16 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<MAT::Material> MAT::PAR::FluidPoroMultiPhase::CreateMaterial()
{
  return Teuchos::rcp(new MAT::FluidPoroMultiPhase(this));
}

/*----------------------------------------------------------------------*
 | global instance of parameter class                       vuong 08/16 |
 *----------------------------------------------------------------------*/
MAT::FluidPoroMultiPhaseType MAT::FluidPoroMultiPhaseType::instance_;

/*----------------------------------------------------------------------*
 | create material from data                                vuong 08/16 |
 *----------------------------------------------------------------------*/
DRT::ParObject* MAT::FluidPoroMultiPhaseType::Create( const std::vector<char> & data )
{
  MAT::FluidPoroMultiPhase* FluidPoroMultiPhase = new MAT::FluidPoroMultiPhase();
  FluidPoroMultiPhase->Unpack(data);
  return FluidPoroMultiPhase;
}


/*----------------------------------------------------------------------*
 | construct empty material object                          vuong 08/16 |
 *----------------------------------------------------------------------*/
MAT::FluidPoroMultiPhase::FluidPoroMultiPhase()
  : MatList(),
    paramsporo_(NULL)
{
}

/*----------------------------------------------------------------------*
 | construct the material object given material parameter   vuong 08/16 |
 *----------------------------------------------------------------------*/
MAT::FluidPoroMultiPhase::FluidPoroMultiPhase(MAT::PAR::FluidPoroMultiPhase* params)
  : MatList(params),
    paramsporo_(params)
{
}

/*----------------------------------------------------------------------*
 | reset everything                                          vuong 08/16 |
 *----------------------------------------------------------------------*/
void MAT::FluidPoroMultiPhase::Clear()
{
  paramsporo_ = NULL;
  return;
}

/*----------------------------------------------------------------------*
 | initialize                                               vuong 08/16 |
 *----------------------------------------------------------------------*/
void MAT::FluidPoroMultiPhase::Initialize()
{
  std::map<int,Teuchos::RCP<MAT::Material> >* materials;

  if (Parameter() != NULL) // params is null pointer in post-process mode
  {
    if(Parameter()->local_)
      materials = MaterialMapWrite();
    else
      materials = Parameter()->MaterialMapWrite();

    std::map<int,Teuchos::RCP<MAT::Material> >::iterator it;
    for(it=materials->begin();it!=materials->end();it++)
    {
      Teuchos::RCP<MAT::FluidPoroSinglePhase> actphase =
          Teuchos::rcp_dynamic_cast<FluidPoroSinglePhase>(it->second,true);
      actphase->Initialize();
    }
  }
  return;
}

/*----------------------------------------------------------------------*
 | Unpack data from a char vector into this class            vuong 08/16 |
 *----------------------------------------------------------------------*/
void MAT::FluidPoroMultiPhase::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm( data );
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);

  // matid
  int matid = -1;
  if (paramsporo_ != NULL) matid = paramsporo_->Id();  // in case we are in post-process mode

  AddtoPack(data,matid);

  // Pack base class material
  MAT::MatList::Pack(data);
}

/*----------------------------------------------------------------------*
 | Unpack data from a char vector into this class            vuong 08/16 |
 *----------------------------------------------------------------------*/
void MAT::FluidPoroMultiPhase::Unpack(const std::vector<char>& data)
{
  // make sure we have a pristine material
  Clear();

  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");

  // matid and recover paramsporo_
  int matid(-1);
  ExtractfromPack(position,data,matid);
  paramsporo_ = NULL;
  if (DRT::Problem::Instance()->Materials() != Teuchos::null)
    if (DRT::Problem::Instance()->Materials()->Num() != 0)
    {
      const int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();
      MAT::PAR::Parameter* mat = DRT::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
      {
        //Note: We need to do a dynamic_cast here since Chemotaxis, Reaction, and Chemo-reaction are in a diamond inheritance structure
        paramsporo_ = dynamic_cast<MAT::PAR::FluidPoroMultiPhase*>(mat);
      }
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(), MaterialType());
    }

  // extract base class material
  std::vector<char> basedata(0);
  MAT::MatList::ExtractfromPack(position,data,basedata);
  MAT::MatList::Unpack(basedata);

  // in the postprocessing mode, we do not unpack everything we have packed
  // -> position check cannot be done in this case
  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d",data.size(),position);
}

