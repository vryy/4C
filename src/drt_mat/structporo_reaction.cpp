/*----------------------------------------------------------------------*/
/*!
 \file structporo_reaction.cpp

 \brief wrapper for structure material of porous media including reactive reference porosity

   \maintainer Andreas Rauch
               rauch@lnm.mw.tum.de
               http://www.lnm.mw.tum.de

\level 3
 *----------------------------------------------------------------------*/

#include <vector>
#include "structporo_reaction.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_bundle.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::PAR::StructPoroReaction::StructPoroReaction(Teuchos::RCP<MAT::PAR::Material> matdata)
    : StructPoro(matdata), dofIDReacScalar_(matdata->GetInt("DOFIDREACSCALAR"))
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

DRT::ParObject* MAT::StructPoroReactionType::Create(const std::vector<char>& data)
{
  MAT::StructPoroReaction* struct_poro = new MAT::StructPoroReaction();
  struct_poro->Unpack(data);
  return struct_poro;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::StructPoroReaction::StructPoroReaction()
    : params_(NULL), refporosity_(-1.0), dphiDphiref_(0.0), refporositydot_(0.0)
{
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::StructPoroReaction::StructPoroReaction(MAT::PAR::StructPoroReaction* params)
    : StructPoro(params),
      params_(params),
      refporosity_(-1.0),
      dphiDphiref_(0.0),
      refporositydot_(0.0)
{
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void MAT::StructPoroReaction::Setup(int numgp, DRT::INPUT::LineDefinition* linedef)
{
  StructPoro::Setup(numgp, linedef);
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
  if (params_ != NULL) matid = params_->Id();  // in case we are in post-process mode
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
  if (type != UniqueParObjectId()) dserror("wrong instance type data");

  // matid
  int matid;
  ExtractfromPack(position, data, matid);
  params_ = NULL;
  if (DRT::Problem::Instance()->Materials() != Teuchos::null)
    if (DRT::Problem::Instance()->Materials()->Num() != 0)
    {
      const int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();
      MAT::PAR::Parameter* mat =
          DRT::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        params_ = static_cast<MAT::PAR::StructPoroReaction*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }

  // refporosity_
  ExtractfromPack(position, data, refporosity_);

  // extract base class material
  std::vector<char> basedata(0);
  ExtractfromPack(position, data, basedata);
  StructPoro::Unpack(basedata);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::StructPoroReaction::ComputePorosity(Teuchos::ParameterList& params, double press,
    double J, int gp, double& porosity, double* dphi_dp, double* dphi_dJ, double* dphi_dJdp,
    double* dphi_dJJ, double* dphi_dpp, bool save)
{
  // evaluate change of reference porosity due to reaction

  // TODO: do not read from parameter list!
  if (params.isParameter("scalar"))
  {
    Teuchos::RCP<std::vector<double>> scalars =
        params.get<Teuchos::RCP<std::vector<double>>>("scalar");
    Reaction(porosity, J, scalars, params);
  }

  // call base class to compute porosity
  StructPoro::ComputePorosity(refporosity_, press, J, gp, porosity, dphi_dp, dphi_dJ, dphi_dJdp,
      dphi_dJJ, dphi_dpp, &dphiDphiref_, save);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void MAT::StructPoroReaction::ConstitutiveDerivatives(Teuchos::ParameterList& params, double press,
    double J, double porosity, double* dW_dp, double* dW_dphi, double* dW_dJ, double* dW_dphiref,
    double* W)
{
  if (porosity == 0.0)
    dserror(
        "porosity equals zero!! Wrong initial porosity? (or wrong collagen density for ecm "
        "material)");

  // evaluate change of reference porosity due to reaction

  // TODO: do not read from parameter list!
  if (params.isParameter("scalar"))
  {
    Teuchos::RCP<std::vector<double>> scalars =
        params.get<Teuchos::RCP<std::vector<double>>>("scalar");
    Reaction(porosity, J, scalars, params);
  }

  // call base class
  StructPoro::ConstitutiveDerivatives(
      params, press, J, porosity, refporosity_, dW_dp, dW_dphi, dW_dJ, dW_dphiref, W);

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void MAT::StructPoroReaction::Reaction(const double porosity, const double J,
    Teuchos::RCP<std::vector<double>> scalars, Teuchos::ParameterList& params)
{
  // double dt = params.get<double>("delta time",-1.0);
  double time = params.get<double>("total time", -1.0);

  // use scalar for this type of reaction
  double cnp = scalars->at(params_->dofIDReacScalar_);

  if (time != -1.0)
  {
    // double k = 1.0;
    double tau = 200.0 * cnp;     ///(cnp+k); 20.0/(20*cnp+1.0)
    double limitporosity = 0.45;  // 0.8;

    refporosity_ =
        limitporosity - (limitporosity - params_->initporosity_) * exp(-1.0 * time / tau);
    refporositydot_ = (limitporosity - params_->initporosity_) / tau * exp(-1.0 * time / tau);
  }
  else  //(time==-1.0) -> time not set (this happens during setup-> no reaction)
  {
    // do nothing
  }
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void MAT::StructPoroReaction::Evaluate(
    const LINALG::Matrix<3, 3>* defgrd,    ///< (i) deformation gradient
    const LINALG::Matrix<6, 1>* glstrain,  ///< (i) green lagrange strain
    Teuchos::ParameterList& params,        ///< (i) parameter list
    LINALG::Matrix<6, 1>* stress,          ///< (o) second piola kirchhoff stress
    LINALG::Matrix<6, 6>* cmat, const int eleGID)
{
  // call base class
  StructPoro::Evaluate(defgrd, glstrain, params, stress, cmat, eleGID);

  // scale stresses and cmat
  stress->Scale((1.0 - refporosity_) / (1.0 - params_->initporosity_));
  cmat->Scale((1.0 - refporosity_) / (1.0 - params_->initporosity_));
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double MAT::StructPoroReaction::RefPorosityAv() const { return refporosity_; }


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void MAT::StructPoroReaction::VisNames(std::map<std::string, int>& names)
{
  // call base class
  StructPoro::VisNames(names);
  std::string name = "reference_porosity";
  names[name] = 1;  // scalar
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool MAT::StructPoroReaction::VisData(
    const std::string& name, std::vector<double>& data, int numgp, int eleID)
{
  // call base class
  if (StructPoro::VisData(name, data, numgp, eleID)) return true;
  if (name == "reference_porosity")
  {
    if ((int)data.size() != 1) dserror("size mismatch");
    data[0] = RefPorosityAv();
    return true;
  }
  return false;
}
