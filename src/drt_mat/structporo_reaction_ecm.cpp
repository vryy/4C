/*----------------------------------------------------------------------*/
/*! \file
 \brief solid material for  porous model with special implementations
        for ECM model (dissolution of ECM)

   \maintainer Christoph Schmidt

\level 3
 *----------------------------------------------------------------------*/



#include <vector>
#include "structporo_reaction_ecm.H"
#include "poro_law.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_bundle.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::PAR::StructPoroReactionECM::StructPoroReactionECM(Teuchos::RCP<MAT::PAR::Material> matdata)
    : StructPoroReaction(matdata), densCollagen_(matdata->GetDouble("DENSCOLLAGEN"))
{
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<MAT::Material> MAT::PAR::StructPoroReactionECM::CreateMaterial()
{
  return Teuchos::rcp(new MAT::StructPoroReactionECM(this));
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::StructPoroReactionECMType MAT::StructPoroReactionECMType::instance_;

DRT::ParObject* MAT::StructPoroReactionECMType::Create(const std::vector<char>& data)
{
  MAT::StructPoroReactionECM* struct_poro = new MAT::StructPoroReactionECM();
  struct_poro->Unpack(data);
  return struct_poro;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::StructPoroReactionECM::StructPoroReactionECM()
    : refporosity_old_(-1.0), refporositydot_old_(0.0), chempot_(0), chempot_init_(0), params_(NULL)
{
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::StructPoroReactionECM::StructPoroReactionECM(MAT::PAR::StructPoroReactionECM* params)
    : StructPoroReaction(params),
      refporosity_old_(-1.0),
      refporositydot_old_(0.0),
      chempot_(0),
      chempot_init_(0),
      params_(params)
{
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void MAT::StructPoroReactionECM::Setup(int numgp, DRT::INPUT::LineDefinition* linedef)
{
  StructPoroReaction::Setup(numgp, linedef);
  refporosity_old_ = params_->initporosity_;

  double dpsidphiref = 0.0;
  Teuchos::ParameterList params;
  params_->porolaw_->ConstitutiveDerivatives(
      params, 0.0, 1.0, params_->initporosity_, refporosity_, NULL, NULL, NULL, &dpsidphiref, NULL);

  const double initphi = params_->initporosity_;
  const double deltaphi = refporosity_ - initphi;

  chempot_.resize(numgp, 0.0);
  chempot_init_.resize(numgp, 0.0);

  for (std::vector<double>::size_type i = 0; i < chempot_init_.size(); i++)
    chempot_init_[i] = -(1.0 - deltaphi / (1.0 - initphi)) / mat_->Density() * dpsidphiref;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::StructPoroReactionECM::Pack(DRT::PackBuffer& data) const
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
  AddtoPack(data, refporosity_old_);
  // refporositydot_old_
  AddtoPack(data, refporositydot_old_);
  // chempot_init_
  AddtoPack(data, chempot_init_);
  // chempot_
  AddtoPack(data, chempot_);

  // add base class material
  StructPoroReaction::Pack(data);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::StructPoroReactionECM::Unpack(const std::vector<char>& data)
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
        params_ = static_cast<MAT::PAR::StructPoroReactionECM*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }

  ExtractfromPack(position, data, refporosity_old_);
  ExtractfromPack(position, data, refporositydot_old_);
  ExtractfromPack(position, data, chempot_init_);
  ExtractfromPack(position, data, chempot_);

  // extract base class material
  std::vector<char> basedata(0);
  ExtractfromPack(position, data, basedata);
  StructPoroReaction::Unpack(basedata);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void MAT::StructPoroReactionECM::Reaction(const double porosity, const double J,
    Teuchos::RCP<std::vector<double>> scalars, Teuchos::ParameterList& params)
{
  double dt = params.get<double>("delta time", -1.0);
  // double time = params.get<double>("total time",-1.0);

  // concentration C1
  double c1 = scalars->at(params_->dofIDReacScalar_);

  if (dt < 0.0) dserror("time step not available");

  refporosity_ = 1.0 - J * c1 / params_->densCollagen_;
  refporositydot_ = (refporosity_ - refporosity_old_) / dt;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void MAT::StructPoroReactionECM::Update()
{
  refporosity_old_ = refporosity_;
  refporositydot_old_ = refporositydot_;

  StructPoroReaction::Update();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void MAT::StructPoroReactionECM::VisNames(std::map<std::string, int>& names)
{
  // call base class
  StructPoroReaction::VisNames(names);

  std::string name = "chemical_potential";
  names[name] = 1;  // scalar
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool MAT::StructPoroReactionECM::VisData(
    const std::string& name, std::vector<double>& data, int numgp, int eleID)
{
  // call base class
  if (StructPoroReaction::VisData(name, data, numgp, eleID)) return true;
  if (name == "chemical_potential")
  {
    if ((int)data.size() != 1) dserror("size mismatch");
    for (std::vector<double>::size_type i = 0; i < chempot_.size(); i++)
      data[0] += chempot_[i] / chempot_.size();
    return true;
  }
  return false;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void MAT::StructPoroReactionECM::ChemPotential(
    const LINALG::Matrix<6, 1>& glstrain,  ///< (i) green lagrange strain
    const double porosity,                 ///< (i) porosity
    const double press,                    ///< (i) pressure at gauss point
    const double J,                        ///< (i) determinant of jacobian at gauss point
    int EleID,                             ///< (i) element GID
    double& pot,                           ///< (o) chemical potential
    const int gp)
{
  dsassert(gp < (int)chempot_.size(),
      "invalid gauss point number for calculation of chemical potential!");
  dsassert(gp < (int)chempot_init_.size(),
      "invalid gauss point number for calculation of chemical potential!");

  // dummy parameter list
  Teuchos::ParameterList params;

  double psi = 0.0;
  // evaluate strain energy
  mat_->StrainEnergy(glstrain, psi, gp, EleID);

  // derivative of
  double dpsidphiref = 0.0;

  params_->porolaw_->ConstitutiveDerivatives(
      params, press, J, porosity, refporosity_, NULL, NULL, NULL, &dpsidphiref, NULL);

  pot = 1.0 / Density() * psi - 1.0 / mat_->Density() * dpsidphiref - chempot_init_[gp];
  chempot_[gp] = pot;

  return;
}
