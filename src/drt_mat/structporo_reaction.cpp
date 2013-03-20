/*----------------------------------------------------------------------*/
/*!
 \file StructPoroReaction_reaction.cpp

 \brief wrapper for structure material of porous media including reactiive reference porosity

 <pre>
   Maintainer: Anh-Tu Vuong
               vuong@lnm.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-15264
 </pre>
 *----------------------------------------------------------------------*/

#include <vector>
#include "structporo_reaction.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_bundle.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::PAR::StructPoroReaction::StructPoroReaction(Teuchos::RCP<MAT::PAR::Material> matdata) :
  StructPoro(matdata)
{
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<MAT::Material> MAT::PAR::StructPoroReaction::CreateMaterial()
{
  return Teuchos::rcp(new MAT::StructPoroReaction(this));
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::StructPoroReactionType MAT::StructPoroReactionType::instance_;

DRT::ParObject* MAT::StructPoroReactionType::Create(const std::vector<char> & data)
{
  MAT::StructPoroReaction* struct_poro = new MAT::StructPoroReaction();
  struct_poro->Unpack(data);
  return struct_poro;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::StructPoroReaction::StructPoroReaction() :
  params_(NULL),
  refporosity_(-1.0),
  dphiDphiref_(0.0),
  refporositydot_(0.0)
{
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::StructPoroReaction::StructPoroReaction(MAT::PAR::StructPoroReaction* params) :
  StructPoro(params),
  params_(params),
  refporosity_(-1.0),
  dphiDphiref_(0.0),
  refporositydot_(0.0)
{
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void MAT::StructPoroReaction::Setup(int numgp,DRT::INPUT::LineDefinition* linedef)
{
  StructPoro::Setup(numgp,linedef);
  refporosity_ = params_->initporosity_;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::StructPoroReaction::Pack(DRT::PackBuffer& data) const
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

  // refporosity_
  AddtoPack(data, refporosity_);

  // add base class material
  StructPoro::Pack(data);

}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::StructPoroReaction::Unpack(const std::vector<char>& data)
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
        params_ = static_cast<MAT::PAR::StructPoroReaction*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(), MaterialType());
    }

  // refporosity_
  ExtractfromPack(position,data,refporosity_);

  // extract base class material
  std::vector<char> basedata(0);
  ExtractfromPack(position,data,basedata);
  StructPoro::Unpack(basedata);

}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::StructPoroReaction::ComputePorosity( Teuchos::ParameterList& params,
                                       double press,
                                       double J,
                                       int gp,
                                       double& porosity,
                                       double* dphi_dp,
                                       double* dphi_dJ,
                                       double* dphi_dJdp,
                                       double* dphi_dJJ,
                                       double* dphi_dpp,
                                       bool save)
{
  //evaluate change of reference porosity due to reaction
  double cnp = params.get<double>("scalar");
  Reaction(cnp,params);

  //call base class to compute porosity
  StructPoro::ComputePorosity(
                   refporosity_,
                   press,
                   J,
                   gp,
                   porosity,
                   dphi_dp,
                   dphi_dJ,
                   dphi_dJdp,
                   dphi_dJJ,
                   dphi_dpp,
                   &dphiDphiref_,
                   save);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void MAT::StructPoroReaction::Reaction(double cnp, Teuchos::ParameterList& params)
{
  //double dt = params.get<double>("delta time",-1.0);
  double time = params.get<double>("total time",-1.0);

  if(time==-1.0)
    dserror("time step or total time not available");

 // double k = 1.0;
  double tau = 200.0*cnp;///(cnp+k);
  double limitporosity = 0.45;

  refporosity_= limitporosity - (limitporosity - params_->initporosity_)* exp(-1.0*time/tau);
  refporositydot_= (limitporosity - params_->initporosity_)/tau * exp(-1.0*time/tau);
}



