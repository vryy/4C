/*----------------------------------------------------------------------*/
/*!
\brief
This file contains the hyperelastic toolbox. It allows summing up several summands
of several types (isotropic or anisotropic, splitted or not) to build a hyperelastic
strain energy function.

The input line should read
MAT 0   MAT_ElastHyper   NUMMAT 2 MATIDS 1 2 DENS 0

\level 1

\maintainer Fabian Braeu

*/
/*----------------------------------------------------------------------*/

#include "elasthyper.H"
#include "../drt_lib/standardtypes_cpp.H"
#include "../drt_matelast/elast_summand.H"
#include "../linalg/linalg_utils.H"
#include "../drt_lib/drt_linedefinition.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_mat/material_service.H"
#include "../drt_comm/comm_utils.H"              // for stat inverse analysis
#include "../drt_inpar/inpar_statinvanalysis.H"  // for stat inverse analysis

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::PAR::ElastHyper::ElastHyper(Teuchos::RCP<MAT::PAR::Material> matdata)
    : Parameter(matdata),
      nummat_(matdata->GetInt("NUMMAT")),
      matids_(matdata->Get<std::vector<int>>("MATIDS")),
      density_(matdata->GetDouble("DENS")),
      polyconvex_(matdata->GetInt("POLYCONVEX")),
      statiaelasthyper_(NULL)

{
  // check if sizes fit
  if (nummat_ != (int)matids_->size())
    dserror("number of materials %d does not fit to size of material vector %d", nummat_,
        matids_->size());

  // output, that polyconvexity is checked
  if (polyconvex_) std::cout << "Polyconvexity of your simulation is checked." << std::endl;

  // STAT INVERSE ANALYSIS
  // For stat inverse analysis, add all parameters to matparams_
  // set size of matparams_ here, add values in summands
  Epetra_Map dummy_map(1, 1, 0, *(DRT::Problem::Instance()->GetNPGroup()->LocalComm()));
  for (int i = MAT::ELASTIC::PAR::first; i <= MAT::ELASTIC::PAR::last; i++)
  {
    matparams_.push_back(Teuchos::rcp(new Epetra_Vector(dummy_map, true)));
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<MAT::Material> MAT::PAR::ElastHyper::CreateMaterial()
{
  return Teuchos::rcp(new MAT::ElastHyper(this));
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::PAR::ElastHyper::OptParams(std::map<std::string, int>* pnames)
{
  statiaelasthyper_->ElastOptParams(pnames);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::ElastHyperType MAT::ElastHyperType::instance_;


DRT::ParObject* MAT::ElastHyperType::Create(const std::vector<char>& data)
{
  MAT::ElastHyper* elhy = new MAT::ElastHyper();
  elhy->Unpack(data);

  return elhy;
}


/*----------------------------------------------------------------------*
 |  initialise static arrays                                 bborn 08/09|
 *----------------------------------------------------------------------*/
// 6-Voigt C-index                              0 1 2  3 4 5
const int MAT::ElastHyper::VOIGT6ROW_[6] = {0, 1, 2, 0, 1, 2};
const int MAT::ElastHyper::VOIGT6COL_[6] = {0, 1, 2, 1, 2, 0};

// tensor indices ij = 11, 12, 13, 21, 22, 23, 31, 32, 33
// C indices           00, 01, 02, 10, 11, 12, 20, 21, 22
// Access : 3*i+j
// 6-Voigt C-indices    0   3   5   3   1   4   5   4   2
const int MAT::ElastHyper::VOIGT3X3SYM_[9] = {0, 3, 5, 3, 1, 4, 5, 4, 2};


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::ElastHyper::ElastHyper()
    : isoprinc_(false),
      isomod_(false),
      anisoprinc_(false),
      anisomod_(false),
      params_(NULL),
      potsum_(0)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::ElastHyper::ElastHyper(MAT::PAR::ElastHyper* params)
    : isoprinc_(false),
      isomod_(false),
      anisoprinc_(false),
      anisomod_(false),
      params_(params),
      potsum_(0)
{
  // make sure the referenced materials in material list have quick access parameters
  std::vector<int>::const_iterator m;
  for (m = params_->matids_->begin(); m != params_->matids_->end(); ++m)
  {
    const int matid = *m;
    Teuchos::RCP<MAT::ELASTIC::Summand> sum = MAT::ELASTIC::Summand::Factory(matid);
    if (sum == Teuchos::null) dserror("Failed to allocate");
    potsum_.push_back(sum);
  }

  // For Stat Inverse Analysis
  // pointer to elasthyper
  params_->SetMaterialPtrSIA(this);

  // just in case of stat inverse analysis (so far just tested for lbfgs)
  const Teuchos::ParameterList& invp = DRT::Problem::Instance()->StatInverseAnalysisParams();
  if (DRT::INPUT::IntegralValue<INPAR::INVANA::StatInvAnalysisType>(invp, "STAT_INV_ANALYSIS") ==
      INPAR::INVANA::stat_inv_lbfgs)
  {
    // copy matparams_ to summands, to fill it with respective parameters
    // loop map of associated potential summands
    for (unsigned int p = 0; p < potsum_.size(); ++p)
    {
      potsum_[p]->CopyStatInvAnaMatParams(params_->matparams_);
      potsum_[p]->SetStatInvAnaSummandMatParams();
    }
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ElastHyper::Pack(DRT::PackBuffer& data) const
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
  AddtoPack(data, isoprinc_);
  AddtoPack(data, isomod_);
  AddtoPack(data, anisoprinc_);
  AddtoPack(data, anisomod_);

  if (params_ != NULL)  // summands are not accessible in postprocessing mode
  {
    // loop map of associated potential summands
    for (unsigned int p = 0; p < potsum_.size(); ++p)
    {
      potsum_[p]->PackSummand(data);
    }
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ElastHyper::Unpack(const std::vector<char>& data)
{
  // make sure we have a pristine material
  params_ = NULL;
  potsum_.clear();

  isoprinc_ = false;
  isomod_ = false;
  anisoprinc_ = false;
  anisomod_ = false;

  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position, data, type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");

  // matid and recover params_
  int matid;
  ExtractfromPack(position, data, matid);
  if (DRT::Problem::Instance()->Materials() != Teuchos::null)
  {
    if (DRT::Problem::Instance()->Materials()->Num() != 0)
    {
      const unsigned int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();
      MAT::PAR::Parameter* mat =
          DRT::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        params_ = static_cast<MAT::PAR::ElastHyper*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }
  }

  isoprinc_ = (bool)ExtractInt(position, data);
  isomod_ = (bool)ExtractInt(position, data);
  anisoprinc_ = (bool)ExtractInt(position, data);
  anisomod_ = (bool)ExtractInt(position, data);

  if (params_ != NULL)  // summands are not accessible in postprocessing mode
  {
    // make sure the referenced materials in material list have quick access parameters
    std::vector<int>::const_iterator m;
    for (m = params_->matids_->begin(); m != params_->matids_->end(); ++m)
    {
      const int matid = *m;
      Teuchos::RCP<MAT::ELASTIC::Summand> sum = MAT::ELASTIC::Summand::Factory(matid);
      if (sum == Teuchos::null) dserror("Failed to allocate");
      potsum_.push_back(sum);
    }

    // loop map of associated potential summands
    for (unsigned int p = 0; p < potsum_.size(); ++p)
    {
      potsum_[p]->UnpackSummand(data, position);
    }
    // in the postprocessing mode, we do not unpack everything we have packed
    // -> position check cannot be done in this case
    if (position != data.size())
      dserror("Mismatch in size of data %d <-> %d", data.size(), position);
  }

  // For Stat Inverse Analysis
  // pointer to elasthyper
  if (params_ != NULL) params_->SetMaterialPtrSIA(this);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
int MAT::ElastHyper::MatID(const unsigned index) const
{
  if ((int)index < params_->nummat_)
    return params_->matids_->at(index);
  else
  {
    dserror("Index too large");
    return -1;
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double MAT::ElastHyper::ShearMod() const
{
  // principal coefficients
  bool haveshearmod = false;
  double shearmod = 0.0;
  {
    // loop map of associated potential summands
    for (unsigned int p = 0; p < potsum_.size(); ++p)
    {
      potsum_[p]->AddShearMod(haveshearmod, shearmod);
    }
  }

  if (haveshearmod)
  {
    return shearmod;
  }
  else
  {
    dserror("Cannot provide shear modulus equivalent");
    return -1.0;
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double MAT::ElastHyper::GetYoung()
{
  double young, shear, bulk;
  young = shear = bulk = 0.;
  for (unsigned int p = 0; p < potsum_.size(); ++p) potsum_[p]->AddYoungsMod(young, shear, bulk);

  if (bulk != 0. || shear != 0.) young += 9. * bulk * shear / (3. * bulk + shear);

  return young;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ElastHyper::SetupAAA(Teuchos::ParameterList& params, const int eleGID)
{
  // loop map of associated potential summands
  for (unsigned int p = 0; p < potsum_.size(); ++p)
  {
    potsum_[p]->SetupAAA(params, eleGID);
  }

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ElastHyper::Setup(int numgp, DRT::INPUT::LineDefinition* linedef)
{
  // Setup summands
  for (unsigned int p = 0; p < potsum_.size(); ++p)
  {
    potsum_[p]->Setup(linedef);
  }

  // find out which formulations are used
  isoprinc_ = false;
  isomod_ = false;
  anisoprinc_ = false;
  anisomod_ = false;
  bool viscogeneral = false;

  for (unsigned int p = 0; p < potsum_.size(); ++p)
  {
    potsum_[p]->SpecifyFormulation(isoprinc_, isomod_, anisoprinc_, anisomod_, viscogeneral);
  }

  if (viscogeneral)
    dserror(
        "Never use viscoelastic-materials in Elasthyper-Toolbox. Use Viscoelasthyper-Toolbox "
        "instead.");

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ElastHyper::Update()
{
  // loop map of associated potential summands
  for (unsigned int p = 0; p < potsum_.size(); ++p)
  {
    potsum_[p]->Update();
  }

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ElastHyper::GetFiberVecs(std::vector<LINALG::Matrix<3, 1>>& fibervecs)
{
  if (anisoprinc_ || anisomod_)
  {
    for (unsigned int p = 0; p < potsum_.size(); ++p)
    {
      potsum_[p]->GetFiberVecs(fibervecs);
    }
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ElastHyper::EvaluateFiberVecs(
    const double newgamma, const LINALG::Matrix<3, 3>& locsys, const LINALG::Matrix<3, 3>& defgrd)
{
  if (anisoprinc_ || anisomod_)
  {
    for (unsigned int p = 0; p < potsum_.size(); ++p)
    {
      potsum_[p]->SetFiberVecs(newgamma, locsys, defgrd);
    }
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ElastHyper::InvariantsPrincipal(
    LINALG::Matrix<3, 1>& prinv, const LINALG::Matrix<6, 1>& rcg)
{
  // 1st invariant, trace
  prinv(0) = rcg(0) + rcg(1) + rcg(2);
  // 2nd invariant
  prinv(1) = 0.5 * (prinv(0) * prinv(0) - rcg(0) * rcg(0) - rcg(1) * rcg(1) - rcg(2) * rcg(2) -
                       .5 * rcg(3) * rcg(3) - .5 * rcg(4) * rcg(4) - .5 * rcg(5) * rcg(5));
  // 3rd invariant, determinant
  prinv(2) = rcg(0) * rcg(1) * rcg(2) + 0.25 * rcg(3) * rcg(4) * rcg(5) -
             0.25 * rcg(1) * rcg(5) * rcg(5) - 0.25 * rcg(2) * rcg(3) * rcg(3) -
             0.25 * rcg(0) * rcg(4) * rcg(4);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ElastHyper::InvariantsModified(LINALG::Matrix<3, 1>& modinv,  ///< modified invariants
    const LINALG::Matrix<3, 1>& prinv                                   ///< principal invariants
)
{
  // 1st invariant, trace
  modinv(0) = prinv(0) * std::pow(prinv(2), -1. / 3.);
  // 2nd invariant
  modinv(1) = prinv(1) * std::pow(prinv(2), -2. / 3.);
  // J
  modinv(2) = std::pow(prinv(2), 1. / 2.);

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ElastHyper::StretchesPrincipal(
    LINALG::Matrix<3, 1>& prstr, LINALG::Matrix<3, 3>& prdir, const LINALG::Matrix<6, 1>& rcg)
{
  // create right Cauchy-Green 2-tensor
  LINALG::Matrix<3, 3> rcgt(false);
  rcgt(0, 0) = rcg(0);
  rcgt(1, 1) = rcg(1);
  rcgt(2, 2) = rcg(2);
  rcgt(0, 1) = rcgt(1, 0) = 0.5 * rcg(3);
  rcgt(1, 2) = rcgt(2, 1) = 0.5 * rcg(4);
  rcgt(2, 0) = rcgt(0, 2) = 0.5 * rcg(5);

  // eigenvalue decomposition
  LINALG::Matrix<3, 3> prstr2;  // squared principal stretches
  LINALG::SYEV(rcgt, prstr2, prdir);

  // THE principal stretches
  for (int al = 0; al < 3; ++al) prstr(al) = std::sqrt(prstr2(al, al));

  // bye
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ElastHyper::StretchesModified(
    LINALG::Matrix<3, 1>& modstr, const LINALG::Matrix<3, 1>& prstr)
{
  // determinant of deformation gradient
  const double detdefgrad = prstr(0) * prstr(1) * prstr(2);

  // determine modified principal stretches
  modstr.Update(std::pow(detdefgrad, -1.0 / 3.0), prstr);

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool MAT::ElastHyper::HaveCoefficientsStretchesPrincipal()
{
  // set default
  bool havecoeff = false;

  // loop map of associated potential summands and see
  {
    for (unsigned int p = 0; p < potsum_.size(); ++p)
    {
      havecoeff = havecoeff or potsum_[p]->HaveCoefficientsStretchesPrincipal();
    }
  }

  // deliver
  return havecoeff;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool MAT::ElastHyper::HaveCoefficientsStretchesModified()
{
  // set default
  bool havecoeff = false;

  // loop map of associated potential summands and see
  {
    for (unsigned int p = 0; p < potsum_.size(); ++p)
    {
      havecoeff = havecoeff or potsum_[p]->HaveCoefficientsStretchesModified();
    }
  }

  // deliver
  return havecoeff;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ElastHyper::StrainEnergy(
    const LINALG::Matrix<6, 1>& glstrain, double& psi, const int eleGID)
{
  LINALG::Matrix<6, 1> id2(true);
  LINALG::Matrix<6, 1> rcg(true);
  LINALG::Matrix<6, 1> scg(true);
  LINALG::Matrix<6, 1> icg(true);
  LINALG::Matrix<6, 6> id4(true);
  LINALG::Matrix<6, 6> id4sharp(true);

  LINALG::Matrix<3, 1> prinv(true);
  LINALG::Matrix<3, 1> modinv(true);

  // evaluate kinematic quantities
  EvaluateKinQuant(glstrain, id2, scg, rcg, icg, id4, id4sharp, prinv);
  InvariantsModified(modinv, prinv);

  // loop map of associated potential summands
  for (unsigned int p = 0; p < potsum_.size(); ++p)
  {
    potsum_[p]->AddStrainEnergy(psi, prinv, modinv, glstrain, eleGID);
  }

  return;
}


/*----------------------------------------------------------------------*
 |  Evaluate for GEMM time integration                        popp 11/13|
 *----------------------------------------------------------------------*/
void MAT::ElastHyper::EvaluateGEMM(LINALG::Matrix<MAT::NUM_STRESS_3D, 1>* stress,
    LINALG::Matrix<MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D>* cmat, double* density,
    LINALG::Matrix<MAT::NUM_STRESS_3D, 1>* glstrain_m,
    LINALG::Matrix<MAT::NUM_STRESS_3D, 1>* glstrain_new,
    LINALG::Matrix<MAT::NUM_STRESS_3D, 1>* glstrain_old, LINALG::Matrix<3, 3>* rcg_new,
    LINALG::Matrix<3, 3>* rcg_old, const int eleGID)
{
#ifdef DEBUG
  if (!stress) dserror("No stress vector supplied");
  if (!cmat) dserror("No material tangent matrix supplied");
  if (!glstrain_m) dserror("No GL strains supplied");
  if (!glstrain_new) dserror("No GL strains supplied");
  if (!glstrain_old) dserror("No GL strains supplied");
#endif

  // standard material evaluate call at midpoint t_{n+1/2}
  Teuchos::ParameterList params;
  LINALG::Matrix<3, 3> defgrd(true);
  Evaluate(&defgrd, glstrain_m, params, stress, cmat, eleGID);
  *density = Density();

  //**********************************************************************
  // CHECK IF GEMM ALGORITHMIC STRESSES NEED TO BE APPLIED
  //**********************************************************************
  // increment of Cauchy-Green tensor in Voigt notation
  LINALG::Matrix<6, 1> M;
  M(0) = (*rcg_new)(0, 0) - (*rcg_old)(0, 0);
  M(1) = (*rcg_new)(1, 1) - (*rcg_old)(1, 1);
  M(2) = (*rcg_new)(2, 2) - (*rcg_old)(2, 2);
  M(3) = (*rcg_new)(0, 1) + (*rcg_new)(1, 0) - (*rcg_old)(0, 1) - (*rcg_old)(1, 0);
  M(4) = (*rcg_new)(1, 2) + (*rcg_new)(2, 1) - (*rcg_old)(1, 2) - (*rcg_old)(2, 1);
  M(5) = (*rcg_new)(0, 2) + (*rcg_new)(2, 0) - (*rcg_old)(0, 2) - (*rcg_old)(2, 0);

  // second variant of M in Voigt notation
  LINALG::Matrix<6, 1> Mtilde;
  Mtilde(0) = M(0);
  Mtilde(1) = M(1);
  Mtilde(2) = M(2);
  Mtilde(3) = 0.5 * M(3);
  Mtilde(4) = 0.5 * M(4);
  Mtilde(5) = 0.5 * M(5);

  // dot product M * Mtilde
  double Mb = M(0) * Mtilde(0) + M(1) * Mtilde(1) + M(2) * Mtilde(2) + M(3) * Mtilde(3) +
              M(4) * Mtilde(4) + M(5) * Mtilde(5);

  // second term in algorithmic stresses only exists if Mb > 0
  // see: O. Gonzalez, Exact energy and momentum conserving algorithms for
  // general models in nonlinear elasticity, CMAME, 190(2000), pp. 1763-1783
  if (Mb < 1.0e-12) return;

  //**********************************************************************
  // COMPUTE GEMM ALGORITHMIC STRESSES
  //**********************************************************************
  // some helper definitions
  LINALG::Matrix<6, 1> vecid(true);
  for (int k = 0; k < 6; ++k) vecid(k) = 1.0;
  LINALG::Matrix<6, 6> halfid(true);
  for (int k = 0; k < 3; ++k) halfid(k, k) = 1.0;
  for (int k = 3; k < 6; ++k) halfid(k, k) = 0.5;

  // strain energy function at t_{n+1} and t_{n}
  double psi = 0.0;
  double psio = 0.0;
  StrainEnergy(*glstrain_new, psi, eleGID);
  StrainEnergy(*glstrain_old, psio, eleGID);

  // derivative of strain energy function dpsi = 0.5*stress
  // double contraction dpsi : M
  double dpsiM = 0.5 * (*stress)(0) * M(0) + 0.5 * (*stress)(1) * M(1) + 0.5 * (*stress)(2) * M(2) +
                 0.5 * (*stress)(3) * M(3) + 0.5 * (*stress)(4) * M(4) + 0.5 * (*stress)(5) * M(5);

  // factor for algorithmic stresses
  double fac = 2.0 * ((psi - psio - dpsiM) / Mb);

  // algorithmic stresses
  LINALG::Matrix<6, 1> algstress(true);
  algstress.Update(fac, Mtilde, 1.0);

  //**********************************************************************
  // COMPUTE GEMM ALGORITHMIC MATERIAL TENSOR
  //**********************************************************************
  // algorithmic material tensor requires stresses at t_{n+1}
  LINALG::Matrix<6, 1> stressnew(true);
  LINALG::Matrix<6, 6> cmatnew(true);
  Evaluate(&defgrd, glstrain_new, params, &stressnew, &cmatnew, eleGID);

  // initialize algorithmic material tensor
  LINALG::Matrix<6, 6> algcmat(true);

  // part 1 (derivative of Mtilde)
  algcmat.Update(4.0 * fac, halfid, 1.0);

  // part 2a (derivative of strain energy in fac)
  LINALG::Matrix<6, 1> dfac(true);
  dfac.Update(2.0 / Mb, stressnew, 1.0);

  // part 2b (derivative of dpsiM in fac)
  LINALG::Matrix<6, 1> tmp(true);
  tmp.Multiply(*cmat, M);
  dfac.Update(-0.5 / Mb, tmp, 1.0);
  dfac.Update(-2.0 / Mb, *stress, 1.0);

  // part 2c (derivative of Mb in fac)
  tmp.Multiply(halfid, M);
  dfac.Update(-4.0 * (psi - psio - dpsiM) / (Mb * Mb), tmp, 1.0);
  dfac.Update(-4.0 * (psi - psio - dpsiM) / (Mb * Mb), Mtilde, 1.0);

  // part 2 (derivative of fac, put together parts 2a,2b and 2c)
  LINALG::Matrix<6, 6> tmpmat(true);
  tmpmat.MultiplyNT(2.0, Mtilde, dfac);
  algcmat.Update(1.0, tmpmat, 1.0);

  //**********************************************************************
  // EXTEND ORIGINAL STRESSES / CMAT WITH GEMM CONTRIBUTIONS
  //**********************************************************************
  stress->Update(1.0, algstress, 1.0);
  cmat->Update(1.0, algcmat, 1.0);

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ElastHyper::Evaluate(const LINALG::Matrix<3, 3>* defgrd,
    const LINALG::Matrix<6, 1>* glstrain, Teuchos::ParameterList& params,
    LINALG::Matrix<6, 1>* stress, LINALG::Matrix<6, 6>* cmat, const int eleGID)
{
  LINALG::Matrix<6, 1> id2(true);
  LINALG::Matrix<6, 1> rcg(true);
  LINALG::Matrix<6, 1> scg(true);
  LINALG::Matrix<6, 1> icg(true);
  LINALG::Matrix<6, 6> id4(true);
  LINALG::Matrix<6, 6> id4sharp(true);

  LINALG::Matrix<3, 1> prinv(true);
  LINALG::Matrix<3, 1> dPI(true);
  LINALG::Matrix<6, 1> ddPII(true);

  EvaluateKinQuant(*glstrain, id2, scg, rcg, icg, id4, id4sharp, prinv);
  EvaluateInvariantDerivatives(prinv, dPI, ddPII, eleGID, potsum_);

  // check if system is polyconvex (set "POLYCONVEX 1" in material input-line)
  if (params_ != NULL)
    if (params_->polyconvex_) CheckPolyconvexity(*defgrd, prinv, dPI, ddPII, params, eleGID);

  // blank resulting quantities
  // ... even if it is an implicit law that cmat is zero upon input
  stress->Clear();
  cmat->Clear();

  // build stress response and elasticity tensor
  LINALG::Matrix<NUM_STRESS_3D, 1> stressiso(true);
  LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D> cmatiso(true);

  EvaluateIsotropicStressCmat(stressiso, cmatiso, scg, id2, icg, id4sharp, prinv, dPI, ddPII);

  stress->Update(1.0, stressiso, 1.0);
  cmat->Update(1.0, cmatiso, 1.0);

  /*----------------------------------------------------------------------*/
  // coefficients in principal stretches
  const bool havecoeffstrpr = HaveCoefficientsStretchesPrincipal();
  const bool havecoeffstrmod = HaveCoefficientsStretchesModified();
  if (havecoeffstrpr or havecoeffstrmod)
  {
    ResponseStretches(*cmat, *stress, rcg, havecoeffstrpr, havecoeffstrmod, eleGID);
  }

  /*----------------------------------------------------------------------*/
  // Do all the anisotropic stuff!
  if (anisoprinc_)
  {
    LINALG::Matrix<NUM_STRESS_3D, 1> stressanisoprinc(true);
    LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D> cmatanisoprinc(true);
    EvaluateAnisotropicPrinc(stressanisoprinc, cmatanisoprinc, rcg, params, eleGID);
    stress->Update(1.0, stressanisoprinc, 1.0);
    cmat->Update(1.0, cmatanisoprinc, 1.0);
  }

  if (anisomod_)
  {
    LINALG::Matrix<NUM_STRESS_3D, 1> stressanisomod(true);
    LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D> cmatanisomod(true);
    EvaluateAnisotropicMod(stressanisomod, cmatanisomod, rcg, icg, prinv, eleGID);
    stress->Update(1.0, stressanisomod, 1.0);
    cmat->Update(1.0, cmatanisomod, 1.0);
  }

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ElastHyper::EvaluateCauchyDerivs(const LINALG::Matrix<3, 1>& prinv, const int eleGID,
    LINALG::Matrix<3, 1>& dPI, LINALG::Matrix<6, 1>& ddPII, LINALG::Matrix<10, 1>& dddPIII,
    const double* temp)
{
  for (unsigned i = 0; i < potsum_.size(); ++i)
  {
    if (isoprinc_)
    {
      potsum_[i]->AddDerivativesPrincipal(dPI, ddPII, prinv, eleGID);
      potsum_[i]->AddThirdDerivativesPrincipalIso(dddPIII, prinv, eleGID);
    }
    if (isomod_ || anisomod_ || anisoprinc_)
      dserror("not implemented for this form of strain energy function");
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ElastHyper::EvaluateCauchy(const LINALG::Matrix<3, 3>& defgrd,
    const LINALG::Matrix<3, 1>& n, const LINALG::Matrix<3, 1>& t, double& snt,
    LINALG::Matrix<3, 1>* DsntDn, LINALG::Matrix<3, 1>* DsntDt, LINALG::Matrix<9, 1>* DsntDF,
    LINALG::Matrix<9, 9>* D2sntDF2, LINALG::Matrix<9, 3>* D2sntDFDn,
    LINALG::Matrix<9, 3>* D2sntDFDt, const int eleGID, const double* temp, double* DsntDT,
    LINALG::Matrix<9, 1>* D2sntDFDT)
{
  snt = 0.0;

  const int VOIGT3X3NONSYM_[3][3] = {{0, 3, 5}, {6, 1, 4}, {8, 7, 2}};

  static LINALG::Matrix<3, 3> b(true);
  b.MultiplyNT(1.0, defgrd, defgrd, 0.0);
  static LINALG::Matrix<3, 1> bdn(true);
  bdn.Multiply(1.0, b, n, 0.0);
  static LINALG::Matrix<3, 1> bdt(true);
  bdt.Multiply(1.0, b, t, 0.0);
  const double bdndt = bdn.Dot(t);

  static LINALG::Matrix<3, 3> ib(true);
  ib.Invert(b);
  static LINALG::Matrix<3, 1> ibdn(true);
  ibdn.Multiply(1.0, ib, n, 0.0);
  static LINALG::Matrix<3, 1> ibdt(true);
  ibdt.Multiply(1.0, ib, t, 0.0);
  const double ibdndt = ibdn.Dot(t);
  const double ndt = n.Dot(t);

  static LINALG::Matrix<6, 1> bV_strain(true);
  MatrixtoStrainLikeVoigtNotation(b, bV_strain);
  static LINALG::Matrix<3, 1> prinv(true);
  InvariantsPrincipal(prinv, bV_strain);

  static LINALG::Matrix<3, 1> dPI(true);
  static LINALG::Matrix<6, 1> ddPII(true);
  static LINALG::Matrix<10, 1> dddPIII(true);
  dPI.Clear();
  ddPII.Clear();
  dddPIII.Clear();
  EvaluateCauchyDerivs(prinv, eleGID, dPI, ddPII, dddPIII, temp);

  const double prefac = 2.0 / sqrt(prinv(2));

  snt = prefac * (prinv(1) * dPI(1) * ndt + prinv(2) * dPI(2) * ndt + dPI(0) * bdndt -
                     prinv(2) * dPI(1) * ibdndt);

  if (DsntDn)
  {
    DsntDn->Update(prinv(1) * dPI(1) + prinv(2) * dPI(2), t, 0.0);  // clear DsntDn here
    DsntDn->Update(dPI(0), bdt, 1.0);
    DsntDn->Update(-prinv(2) * dPI(1), ibdt, 1.0);
    DsntDn->Scale(prefac);
  }

  if (DsntDt)
  {
    DsntDt->Update(prinv(1) * dPI(1) + prinv(2) * dPI(2), n, 0.0);  // clear DsntDt here
    DsntDt->Update(dPI(0), bdn, 1.0);
    DsntDt->Update(-prinv(2) * dPI(1), ibdn, 1.0);
    DsntDt->Scale(prefac);
  }

  // calculate stuff that is needed for evaluations of derivatives w.r.t. F
  static LINALG::Matrix<9, 1> FV(true);
  Matrix3x3to9x1(defgrd, FV);
  static LINALG::Matrix<3, 3> iF(true);
  iF.Invert(defgrd);
  static LINALG::Matrix<3, 3> iFT(true);
  iFT.UpdateT(iF);
  static LINALG::Matrix<9, 1> iFTV(true);
  Matrix3x3to9x1(iFT, iFTV);

  // calculation of dI_i/dF (derivatives of invariants of b w.r.t. deformation gradient)
  static LINALG::Matrix<3, 3> bdF(true);
  bdF.Multiply(1.0, b, defgrd, 0.0);
  static LINALG::Matrix<9, 1> bdFV(true);
  Matrix3x3to9x1(bdF, bdFV);
  static LINALG::Matrix<3, 3> ibdF(true);
  ibdF.Multiply(1.0, ib, defgrd, 0.0);
  static LINALG::Matrix<9, 1> ibdFV(true);
  Matrix3x3to9x1(ibdF, ibdFV);
  static LINALG::Matrix<9, 1> DI1DF(true);
  DI1DF.Update(2.0, FV, 0.0);
  static LINALG::Matrix<9, 1> DI2DF(true);
  DI2DF.Update(prinv(0), FV, 0.0);
  DI2DF.Update(-1.0, bdFV, 1.0);
  DI2DF.Scale(2.0);
  static LINALG::Matrix<9, 1> DI3DF(true);
  DI3DF.Update(2.0 * prinv(2), ibdFV, 0.0);

  // calculate d(b \cdot n \cdot t)/dF
  static LINALG::Matrix<3, 1> tempvec3x1(true);
  static LINALG::Matrix<1, 3> tempvec1x3(true);
  tempvec1x3.MultiplyTN(1.0, t, defgrd, 0.0);
  static LINALG::Matrix<3, 3> DbdndtDF(true);
  DbdndtDF.MultiplyNN(1.0, n, tempvec1x3, 0.0);
  tempvec1x3.MultiplyTN(1.0, n, defgrd, 0.0);
  DbdndtDF.MultiplyNN(1.0, t, tempvec1x3, 1.0);
  static LINALG::Matrix<9, 1> DbdndtDFV(true);
  Matrix3x3to9x1(DbdndtDF, DbdndtDFV);

  // calculate d(b^{-1} \cdot n \cdot t)/dF
  static LINALG::Matrix<1, 3> tdibdF(true);
  static LINALG::Matrix<1, 3> ndibdF(true);
  tdibdF.MultiplyTN(1.0, t, ibdF, 0.0);
  static LINALG::Matrix<3, 3> DibdndtDF(true);
  DibdndtDF.MultiplyNN(1.0, ibdn, tdibdF, 0.0);
  ndibdF.MultiplyTN(1.0, n, ibdF, 0.0);
  DibdndtDF.MultiplyNN(1.0, ibdt, ndibdF, 1.0);
  DibdndtDF.Scale(-1.0);
  static LINALG::Matrix<9, 1> DibdndtDFV(true);
  Matrix3x3to9x1(DibdndtDF, DibdndtDFV);

  if (temp)
    EvaluateCauchyTempDeriv(prinv, ndt, bdndt, ibdndt, temp, DsntDT, iFTV, DbdndtDFV, DibdndtDFV,
        DI1DF, DI2DF, DI3DF, D2sntDFDT);

  if (DsntDF)
  {
    // next 3 updates add partial derivative of (\sigma * n * t) w.r.t. F for constant invariants
    // 1. part is term arising from d(J^{-1})/dF
    DsntDF->Update(-prefac * (prinv(1) * dPI(1) * ndt + prinv(2) * dPI(2) * ndt + dPI(0) * bdndt -
                                 prinv(2) * dPI(1) * ibdndt),
        iFTV, 0.0);  // DsntDF is cleared here
    // 2. part is term arising from d(b * n * t)/dF
    DsntDF->Update(prefac * dPI(0), DbdndtDFV, 1.0);
    // 3. part is term arising from d(b_el^{-1} * n * t)/dF
    DsntDF->Update(-prefac * prinv(2) * dPI(1), DibdndtDFV, 1.0);
    // add d(sigma * n * t)/dI1 \otimes dI1/dF
    DsntDF->Update(prefac * (prinv(1) * ddPII(5) * ndt + prinv(2) * ddPII(4) * ndt +
                                ddPII(0) * bdndt - prinv(2) * ddPII(5) * ibdndt),
        DI1DF, 1.0);
    // add d(sigma * n * t)/dI2 \otimes dI2/dF
    DsntDF->Update(prefac * (dPI(1) * ndt + prinv(1) * ddPII(1) * ndt + prinv(2) * ddPII(3) * ndt +
                                ddPII(5) * bdndt - prinv(2) * ddPII(1) * ibdndt),
        DI2DF, 1.0);
    // add d(sigma * n * t)/dI3 \otimes dI3/dF
    DsntDF->Update(prefac * (prinv(1) * ddPII(3) * ndt + dPI(2) * ndt + prinv(2) * ddPII(2) * ndt +
                                ddPII(4) * bdndt - dPI(1) * ibdndt - prinv(2) * ddPII(3) * ibdndt),
        DI3DF, 1.0);
  }

  if (D2sntDFDn)
  {
    // next three blocks add d/dn(d(\sigma * n * t)/dF) for constant invariants
    // first part is term arising from d/dn(dJ^{-1}/dF)
    tempvec3x1.Update(prinv(1) * dPI(1) + prinv(2) * dPI(2), t, 0.0);
    tempvec3x1.Update(dPI(0), bdt, 1.0);
    tempvec3x1.Update(-prinv(2) * dPI(1), ibdt, 1.0);
    D2sntDFDn->MultiplyNT(-prefac, iFTV, tempvec3x1, 0.0);  // D2sntDFDn is cleared here

    // second part is term arising from d/dn(d(b * n * t)/dF
    const double fac = prefac * dPI(0);
    tempvec1x3.MultiplyTN(1.0, t, defgrd, 0.0);
    for (int k = 0; k < 3; ++k)
      for (int l = 0; l < 3; ++l)
        for (int z = 0; z < 3; ++z)
          (*D2sntDFDn)(VOIGT3X3NONSYM_[k][l], z) +=
              fac * (t(k, 0) * defgrd(z, l) + (k == z) * tempvec1x3(0, l));

    // third part is term arising from d/dn(d(b^{-1} * n * t)/dF
    const double fac2 = prefac * prinv(2) * dPI(1);
    for (int k = 0; k < 3; ++k)
      for (int l = 0; l < 3; ++l)
        for (int z = 0; z < 3; ++z)
          (*D2sntDFDn)(VOIGT3X3NONSYM_[k][l], z) +=
              fac2 * (ibdt(k, 0) * ibdF(z, l) + ib(z, k) * tdibdF(0, l));

    // add parts originating from d/dn(d(sigma * n * t)/dI1 \otimes dI1/dF)
    tempvec3x1.Update(prinv(1) * ddPII(5) + prinv(2) * ddPII(4), t, 0.0);
    tempvec3x1.Update(ddPII(0), bdt, 1.0);
    tempvec3x1.Update(-prinv(2) * ddPII(5), ibdt, 1.0);
    D2sntDFDn->MultiplyNT(prefac, DI1DF, tempvec3x1, 1.0);

    // add parts originating from d/dn(d(sigma * n * t)/dI2 \otimes dI2/dF)
    tempvec3x1.Update(dPI(1) + prinv(1) * ddPII(1) + prinv(2) * ddPII(3), t, 0.0);
    tempvec3x1.Update(ddPII(5), bdt, 1.0);
    tempvec3x1.Update(-prinv(2) * ddPII(1), ibdt, 1.0);
    D2sntDFDn->MultiplyNT(prefac, DI2DF, tempvec3x1, 1.0);

    // add parts originating from d/dn(d(sigma * n * t)/dI3 \otimes dI3/dF)
    tempvec3x1.Update(prinv(1) * ddPII(3) + dPI(2) + prinv(2) * ddPII(2), t, 0.0);
    tempvec3x1.Update(ddPII(4), bdt, 1.0);
    tempvec3x1.Update(-dPI(1) - prinv(2) * ddPII(3), ibdt, 1.0);
    D2sntDFDn->MultiplyNT(prefac, DI3DF, tempvec3x1, 1.0);
  }

  if (D2sntDFDt)
  {
    // next three blocks add d/dt(d(\sigma * n * t)/dF) for constant invariants
    // first part is term arising from d/dt(dJ^{-1}/dF)
    tempvec3x1.Update(prinv(1) * dPI(1) + prinv(2) * dPI(2), n, 0.0);
    tempvec3x1.Update(dPI(0), bdn, 1.0);
    tempvec3x1.Update(-prinv(2) * dPI(1), ibdn, 1.0);
    D2sntDFDt->MultiplyNT(-prefac, iFTV, tempvec3x1, 0.0);  // D2sntDFDt is cleared here

    // second part is term arising from d/dt(d(b * n * t)/dF
    const double fac = prefac * dPI(0);
    tempvec1x3.MultiplyTN(1.0, n, defgrd, 0.0);
    for (int k = 0; k < 3; ++k)
      for (int l = 0; l < 3; ++l)
        for (int z = 0; z < 3; ++z)
          (*D2sntDFDt)(VOIGT3X3NONSYM_[k][l], z) +=
              fac * (n(k, 0) * defgrd(z, l) + (k == z) * tempvec1x3(0, l));

    // third part is term arising from d/dn(d(b^{-1} * n * t)/dF
    const double fac2 = prefac * prinv(2) * dPI(1);
    for (int k = 0; k < 3; ++k)
      for (int l = 0; l < 3; ++l)
        for (int z = 0; z < 3; ++z)
          (*D2sntDFDt)(VOIGT3X3NONSYM_[k][l], z) +=
              fac2 * (ibdn(k, 0) * ibdF(z, l) + ib(z, k) * ndibdF(0, l));

    // add parts originating from d/dt(d(sigma * n * t)/dI1 \otimes dI1/dF)
    tempvec3x1.Update(prinv(1) * ddPII(5) + prinv(2) * ddPII(4), n, 0.0);
    tempvec3x1.Update(ddPII(0), bdn, 1.0);
    tempvec3x1.Update(-prinv(2) * ddPII(5), ibdn, 1.0);
    D2sntDFDt->MultiplyNT(prefac, DI1DF, tempvec3x1, 1.0);

    // add parts originating from d/dt(d(sigma * n * t)/dI2 \otimes dI2/dF)
    tempvec3x1.Update(dPI(1) + prinv(1) * ddPII(1) + prinv(2) * ddPII(3), n, 0.0);
    tempvec3x1.Update(ddPII(5), bdn, 1.0);
    tempvec3x1.Update(-prinv(2) * ddPII(1), ibdn, 1.0);
    D2sntDFDt->MultiplyNT(prefac, DI2DF, tempvec3x1, 1.0);

    // add parts originating from d/dt(d(sigma * n * t)/dI3 \otimes dI3/dF)
    tempvec3x1.Update(prinv(1) * ddPII(3) + dPI(2) + prinv(2) * ddPII(2), n, 0.0);
    tempvec3x1.Update(ddPII(4), bdn, 1.0);
    tempvec3x1.Update(-dPI(1) - prinv(2) * ddPII(3), ibdn, 1.0);
    D2sntDFDt->MultiplyNT(prefac, DI3DF, tempvec3x1, 1.0);
  }

  if (D2sntDF2)
  {
    // define and fill all tensors that can not be calculated using multiply operations first
    static LINALG::Matrix<9, 9> DiFTDF(true);
    static LINALG::Matrix<9, 9> D2bdndtDF2(true);
    static LINALG::Matrix<9, 9> D2ibdndtDF2(true);
    static LINALG::Matrix<9, 9> D2I1DF2(true);
    static LINALG::Matrix<9, 9> D2I2DF2(true);
    static LINALG::Matrix<9, 9> D2I3DF2(true);
    DiFTDF.Clear();
    D2bdndtDF2.Clear();
    D2ibdndtDF2.Clear();
    D2I1DF2.Clear();
    D2I2DF2.Clear();
    D2I3DF2.Clear();

    static LINALG::Matrix<3, 3> C(true);
    C.MultiplyTN(1.0, defgrd, defgrd, 0.0);

    for (int k = 0; k < 3; ++k)
      for (int l = 0; l < 3; ++l)
        for (int m = 0; m < 3; ++m)
          for (int a = 0; a < 3; ++a)
          {
            DiFTDF(VOIGT3X3NONSYM_[k][l], VOIGT3X3NONSYM_[m][a]) = -iF(l, m) * iF(a, k);
            D2bdndtDF2(VOIGT3X3NONSYM_[k][l], VOIGT3X3NONSYM_[m][a]) =
                (t(k, 0) * n(m, 0) + t(m, 0) * n(k, 0)) * (l == a);
            D2ibdndtDF2(VOIGT3X3NONSYM_[k][l], VOIGT3X3NONSYM_[m][a]) =
                ibdF(k, a) * (ibdt(m, 0) * ndibdF(0, l) + ibdn(m, 0) * tdibdF(0, l)) +
                ib(m, k) * (tdibdF(0, a) * ndibdF(0, l) + tdibdF(0, l) * ndibdF(0, a)) +
                ibdF(m, l) * (ibdt(k, 0) * ndibdF(0, a) + tdibdF(0, a) * ibdn(k, 0));
            D2I1DF2(VOIGT3X3NONSYM_[k][l], VOIGT3X3NONSYM_[m][a]) = 2.0 * (k == m) * (l == a);
            D2I2DF2(VOIGT3X3NONSYM_[k][l], VOIGT3X3NONSYM_[m][a]) =
                2.0 * (prinv(0) * (k == m) * (l == a) + 2.0 * defgrd(m, a) * defgrd(k, l) -
                          (k == m) * C(a, l) - defgrd(k, a) * defgrd(m, l) - b(k, m) * (l == a));
            D2I3DF2(VOIGT3X3NONSYM_[k][l], VOIGT3X3NONSYM_[m][a]) =
                2.0 * prinv(2) * (2.0 * ibdF(m, a) * ibdF(k, l) - ibdF(m, l) * ibdF(k, a));
          }

    // terms below add contributions originating from d(1st term of DsntDF)/dF
    D2sntDF2->MultiplyNT(prefac * (prinv(1) * dPI(1) * ndt + prinv(2) * dPI(2) * ndt +
                                      dPI(0) * bdndt - prinv(2) * dPI(1) * ibdndt),
        iFTV, iFTV, 0.0);  // D2sntDF2 is cleared here
    D2sntDF2->Update(-prefac * (prinv(1) * dPI(1) * ndt + prinv(2) * dPI(2) * ndt + dPI(0) * bdndt -
                                   prinv(2) * dPI(1) * ibdndt),
        DiFTDF, 1.0);
    D2sntDF2->MultiplyNT(-prefac * dPI(0), iFTV, DbdndtDFV, 1.0);
    D2sntDF2->MultiplyNT(prefac * prinv(2) * dPI(1), iFTV, DibdndtDFV, 1.0);

    D2sntDF2->MultiplyNT(-prefac * (prinv(1) * ddPII(5) * ndt + prinv(2) * ddPII(4) * ndt +
                                       ddPII(0) * bdndt - prinv(2) * ddPII(5) * ibdndt),
        iFTV, DI1DF, 1.0);
    D2sntDF2->MultiplyNT(
        -prefac * (dPI(1) * ndt + prinv(1) * ddPII(1) * ndt + prinv(2) * ddPII(3) * ndt +
                      ddPII(5) * bdndt - prinv(2) * ddPII(1) * ibdndt),
        iFTV, DI2DF, 1.0);
    D2sntDF2->MultiplyNT(
        -prefac * (prinv(1) * ddPII(3) * ndt + dPI(2) * ndt + prinv(2) * ddPII(2) * ndt +
                      ddPII(4) * bdndt - dPI(1) * ibdndt - prinv(2) * ddPII(3) * ibdndt),
        iFTV, DI3DF, 1.0);

    // terms below add contributions originating from d(2nd term of DsntDF)/dF
    D2sntDF2->MultiplyNT(-prefac * dPI(0), DbdndtDFV, iFTV, 1.0);
    D2sntDF2->Update(prefac * dPI(0), D2bdndtDF2, 1.0);
    D2sntDF2->MultiplyNT(prefac * ddPII(0), DbdndtDFV, DI1DF, 1.0);
    D2sntDF2->MultiplyNT(prefac * ddPII(5), DbdndtDFV, DI2DF, 1.0);
    D2sntDF2->MultiplyNT(prefac * ddPII(4), DbdndtDFV, DI3DF, 1.0);

    // terms below add contributions originating from d(3rd term of DsntDF)/dF
    D2sntDF2->MultiplyNT(prefac * prinv(2) * dPI(1), DibdndtDFV, iFTV, 1.0);
    D2sntDF2->Update(-prefac * prinv(2) * dPI(1), D2ibdndtDF2, 1.0);
    D2sntDF2->MultiplyNT(-prefac * prinv(2) * ddPII(5), DibdndtDFV, DI1DF, 1.0);
    D2sntDF2->MultiplyNT(-prefac * prinv(2) * ddPII(1), DibdndtDFV, DI2DF, 1.0);
    D2sntDF2->MultiplyNT(-prefac * (dPI(1) + prinv(2) * ddPII(3)), DibdndtDFV, DI3DF, 1.0);

    // terms below add contributions originating from d(4th term of DsntDF)/dF
    D2sntDF2->MultiplyNT(-prefac * (prinv(1) * ddPII(5) * ndt + prinv(2) * ddPII(4) * ndt +
                                       ddPII(0) * bdndt - prinv(2) * ddPII(5) * ibdndt),
        DI1DF, iFTV, 1.0);
    D2sntDF2->MultiplyNT(prefac * ddPII(0), DI1DF, DbdndtDFV, 1.0);
    D2sntDF2->MultiplyNT(-prefac * prinv(2) * ddPII(5), DI1DF, DibdndtDFV, 1.0);
    D2sntDF2->Update(prefac * (prinv(1) * ddPII(5) * ndt + prinv(2) * ddPII(4) * ndt +
                                  ddPII(0) * bdndt - prinv(2) * ddPII(5) * ibdndt),
        D2I1DF2, 1.0);
    D2sntDF2->MultiplyNT(prefac * (prinv(1) * dddPIII(5) * ndt + prinv(2) * dddPIII(6) * ndt +
                                      dddPIII(0) * bdndt - prinv(2) * dddPIII(5) * ibdndt),
        DI1DF, DI1DF, 1.0);
    D2sntDF2->MultiplyNT(
        prefac * (ddPII(5) * ndt + prinv(1) * dddPIII(3) * ndt + prinv(2) * dddPIII(9) * ndt +
                     dddPIII(5) * bdndt - prinv(2) * dddPIII(3) * ibdndt),
        DI1DF, DI2DF, 1.0);
    D2sntDF2->MultiplyNT(
        prefac * (prinv(1) * dddPIII(9) * ndt + ddPII(4) * ndt + prinv(2) * dddPIII(4) * ndt +
                     dddPIII(6) * bdndt - ddPII(5) * ibdndt - prinv(2) * dddPIII(9) * ibdndt),
        DI1DF, DI3DF, 1.0);

    // terms below add contributions originating from d(5th term of DsntDF)/dF
    D2sntDF2->MultiplyNT(
        -prefac * (dPI(1) * ndt + prinv(1) * ddPII(1) * ndt + prinv(2) * ddPII(3) * ndt +
                      ddPII(5) * bdndt - prinv(2) * ddPII(1) * ibdndt),
        DI2DF, iFTV, 1.0);
    D2sntDF2->MultiplyNT(prefac * ddPII(5), DI2DF, DbdndtDFV, 1.0);
    D2sntDF2->MultiplyNT(-prefac * prinv(2) * ddPII(1), DI2DF, DibdndtDFV, 1.0);
    D2sntDF2->Update(
        prefac * (dPI(1) * ndt + prinv(1) * ddPII(1) * ndt + prinv(2) * ddPII(3) * ndt +
                     ddPII(5) * bdndt - prinv(2) * ddPII(1) * ibdndt),
        D2I2DF2, 1.0);
    D2sntDF2->MultiplyNT(
        prefac * (ddPII(5) * ndt + prinv(1) * dddPIII(3) * ndt + prinv(2) * dddPIII(9) * ndt +
                     dddPIII(5) * bdndt - prinv(2) * dddPIII(3) * ibdndt),
        DI2DF, DI1DF, 1.0);
    D2sntDF2->MultiplyNT(
        prefac * (2.0 * ddPII(1) * ndt + prinv(1) * dddPIII(1) * ndt + prinv(2) * dddPIII(7) * ndt +
                     dddPIII(3) * bdndt - prinv(2) * dddPIII(1) * ibdndt),
        DI2DF, DI2DF, 1.0);
    D2sntDF2->MultiplyNT(
        prefac * (2.0 * ddPII(3) * ndt + prinv(1) * dddPIII(7) * ndt + prinv(2) * dddPIII(8) * ndt +
                     dddPIII(9) * bdndt - ddPII(1) * ibdndt - prinv(2) * dddPIII(7) * ibdndt),
        DI2DF, DI3DF, 1.0);

    // terms below add contributions originating from d(6th term of DsntDF)/dF
    D2sntDF2->MultiplyNT(
        -prefac * (prinv(1) * ddPII(3) * ndt + dPI(2) * ndt + prinv(2) * ddPII(2) * ndt +
                      ddPII(4) * bdndt - dPI(1) * ibdndt - prinv(2) * ddPII(3) * ibdndt),
        DI3DF, iFTV, 1.0);
    D2sntDF2->MultiplyNT(prefac * ddPII(4), DI3DF, DbdndtDFV, 1.0);
    D2sntDF2->MultiplyNT(-prefac * (dPI(1) + prinv(2) * ddPII(3)), DI3DF, DibdndtDFV, 1.0);
    D2sntDF2->Update(
        prefac * (prinv(1) * ddPII(3) * ndt + dPI(2) * ndt + prinv(2) * ddPII(2) * ndt +
                     ddPII(4) * bdndt - dPI(1) * ibdndt - prinv(2) * ddPII(3) * ibdndt),
        D2I3DF2, 1.0);
    D2sntDF2->MultiplyNT(
        prefac * (prinv(1) * dddPIII(9) * ndt + ddPII(4) * ndt + prinv(2) * dddPIII(4) * ndt +
                     dddPIII(6) * bdndt - ddPII(5) * ibdndt - prinv(2) * dddPIII(9) * ibdndt),
        DI3DF, DI1DF, 1.0);
    D2sntDF2->MultiplyNT(
        prefac * (2.0 * ddPII(3) * ndt + prinv(1) * dddPIII(7) * ndt + prinv(2) * dddPIII(8) * ndt +
                     dddPIII(9) * bdndt - ddPII(1) * ibdndt - prinv(2) * dddPIII(7) * ibdndt),
        DI3DF, DI2DF, 1.0);
    D2sntDF2->MultiplyNT(
        prefac * (prinv(1) * dddPIII(8) * ndt + 2.0 * ddPII(2) * ndt + prinv(2) * dddPIII(2) * ndt +
                     dddPIII(4) * bdndt - 2.0 * ddPII(3) * ibdndt - prinv(2) * dddPIII(8) * ibdndt),
        DI3DF, DI3DF, 1.0);
  }

  return;
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ElastHyper::EvaluateKinQuant(const LINALG::Matrix<6, 1>& glstrain,
    LINALG::Matrix<6, 1>& id2, LINALG::Matrix<6, 1>& scg, LINALG::Matrix<6, 1>& rcg,
    LINALG::Matrix<6, 1>& icg, LINALG::Matrix<6, 6>& id4, LINALG::Matrix<6, 6>& id4sharp,
    LINALG::Matrix<3, 1>& prinv)

{
  // build Cartesian identity 2-tensor I_{AB}
  for (int i = 0; i < 3; i++) id2(i) = 1.0;

  // right Cauchy-Green Tensor  C_{AB} = 2 * E_{AB} + I_{AB}
  // REMARK: strain-like 6-Voigt vector
  rcg.Update(2.0, glstrain, 1.0);
  rcg.Update(1.0, id2, 1.0);

  // 'contra-variant' right Cauchy-Green Tensor C^{AB}
  // REMARK: stress-like 6-Voigt vector of right CG
  scg.Update(1.0, rcg, 1.0);
  for (int i = 3; i < 6; i++) scg(i) *= 0.5;

  // principal invariants of right Cauchy-Green strain
  InvariantsPrincipal(prinv, rcg);

  // invert right Cauchy-Green tensor
  // REMARK: stress-like 6-Voigt vector
  {
    icg(0) = (rcg(1) * rcg(2) - 0.25 * rcg(4) * rcg(4)) / prinv(2);
    icg(1) = (rcg(0) * rcg(2) - 0.25 * rcg(5) * rcg(5)) / prinv(2);
    icg(2) = (rcg(0) * rcg(1) - 0.25 * rcg(3) * rcg(3)) / prinv(2);
    icg(3) = (0.25 * rcg(5) * rcg(4) - 0.5 * rcg(3) * rcg(2)) / prinv(2);
    icg(4) = (0.25 * rcg(3) * rcg(5) - 0.5 * rcg(0) * rcg(4)) / prinv(2);
    icg(5) = (0.25 * rcg(3) * rcg(4) - 0.5 * rcg(5) * rcg(1)) / prinv(2);
  }

  // set Cartesian identity 4-tensor in 6-Voigt matrix notation
  // this is fully 'contra-variant' identity tensor, ie I^{ABCD}
  // REMARK: rows are stress-like 6-Voigt
  //         columns are stress-like 6-Voigt
  for (int i = 0; i < 3; i++) id4sharp(i, i) = 1.0;
  for (int i = 3; i < 6; i++) id4sharp(i, i) = 0.5;

  // set Cartesian identity 4-tensor in 6x6-matrix notation (stress-like)
  // this is a 'mixed co- and contra-variant' identity 4-tensor, ie I^{AB}_{CD}
  // REMARK: rows are stress-like 6-Voigt
  //         columns are strain-like 6-Voigt
  for (int i = 0; i < 6; i++) id4(i, i) = 1.0;
}

/*----------------------------------------------------------------------/
 * Reads derivatives with respect to invariants and modified invariants
 * from all materials of the elasthyper-toolbox          birzle 11/2014 */
/*----------------------------------------------------------------------*/
void MAT::ElastHyper::EvaluateInvariantDerivatives(const LINALG::Matrix<3, 1>& prinv,
    LINALG::Matrix<3, 1>& dPI, LINALG::Matrix<6, 1>& ddPII, int eleGID,
    const std::vector<Teuchos::RCP<MAT::ELASTIC::Summand>>& potsum)

{
  // derivatives of principal materials
  if (isoprinc_)
  {
    // loop map of associated potential summands
    for (unsigned int p = 0; p < potsum.size(); ++p)
    {
      potsum[p]->AddDerivativesPrincipal(dPI, ddPII, prinv, eleGID);
    }
  }
  // derivatives of decoupled (volumetric or isochoric) materials
  if (isomod_)
  {
    LINALG::Matrix<3, 1> modinv;
    InvariantsModified(modinv, prinv);
    LINALG::Matrix<3, 1> dPmodI;
    LINALG::Matrix<6, 1> ddPmodII;
    // loop map of associated potential summands
    for (unsigned int p = 0; p < potsum.size(); ++p)
    {
      potsum[p]->AddDerivativesModified(dPmodI, ddPmodII, modinv, eleGID);
    }
    // convert decoupled derivatives to principal derivatives
    ConvertModToPrinc(prinv, dPmodI, ddPmodII, dPI, ddPII);
  }
}


/*----------------------------------------------------------------------/
 * Converts derivatives with respect to modified invariants in derivatives
 * with respect to principal invariants                 birzle 11/2014  */
/*----------------------------------------------------------------------*/
void MAT::ElastHyper::ConvertModToPrinc(const LINALG::Matrix<3, 1>& prinv,
    const LINALG::Matrix<3, 1>& dPmodI, const LINALG::Matrix<6, 1>& ddPmodII,
    LINALG::Matrix<3, 1>& dPI, LINALG::Matrix<6, 1>& ddPII)
{
  // Conversions to dPI
  dPI(0) += std::pow(prinv(2), -1. / 3.) * dPmodI(0);
  dPI(1) += std::pow(prinv(2), -2. / 3.) * dPmodI(1);
  dPI(2) += 0.5 * std::pow(prinv(2), -0.5) * dPmodI(2) -
            1. / 3. * prinv(0) * std::pow(prinv(2), -4. / 3.) * dPmodI(0) -
            2. / 3. * prinv(1) * std::pow(prinv(2), -5. / 3.) * dPmodI(1);

  // Conversions to ddPII
  ddPII(0) += std::pow(prinv(2), -2. / 3.) * ddPmodII(0);
  ddPII(1) += std::pow(prinv(2), -4. / 3.) * ddPmodII(1);
  ddPII(2) += (1. / 9.) * std::pow(prinv(2), -8. / 3.) * prinv(0) * prinv(0) * ddPmodII(0) +
              (4. / 9.) * prinv(0) * prinv(1) * std::pow(prinv(2), -3.) * ddPmodII(5) -
              (1. / 3.) * std::pow(prinv(2), -11. / 6.) * prinv(0) * ddPmodII(4) +
              (4. / 9.) * std::pow(prinv(2), -7. / 3.) * prinv(0) * dPmodI(0) +
              (4. / 9.) * std::pow(prinv(2), -10. / 3.) * prinv(1) * prinv(1) * ddPmodII(1) -
              (2. / 3.) * std::pow(prinv(2), -13. / 6.) * prinv(1) * ddPmodII(3) +
              (10. / 9.) * std::pow(prinv(2), -8. / 3.) * prinv(1) * dPmodI(1) +
              0.25 * std::pow(prinv(2), -1.) * ddPmodII(2) -
              0.25 * std::pow(prinv(2), -1.5) * dPmodI(2);
  ddPII(3) += -(1. / 3.) * std::pow(prinv(2), -2.) * prinv(0) * ddPmodII(5) -
              (2. / 3.) * std::pow(prinv(2), -7. / 3.) * prinv(1) * ddPmodII(1) +
              0.5 * std::pow(prinv(2), -7. / 6.) * ddPmodII(3) -
              (2. / 3.) * std::pow(prinv(2), -5. / 3.) * dPmodI(1);
  ddPII(4) += -(1. / 3.) * std::pow(prinv(2), -5. / 3.) * prinv(0) * ddPmodII(0) -
              (2. / 3.) * std::pow(prinv(2), -2.) * prinv(1) * ddPmodII(5) +
              0.5 * std::pow(prinv(2), -5. / 6.) * ddPmodII(4) -
              (1. / 3.) * std::pow(prinv(2), -4. / 3.) * dPmodI(0);
  ddPII(5) += std::pow(prinv(2), -1.) * ddPmodII(5);
}

/*----------------------------------------------------------------------*
 * Calculate the coefficients gamma and delta from the partial          *
 * derivatives w.r.t. invariants                           seitz 01/15  *
 *----------------------------------------------------------------------*/
void MAT::ElastHyper::CalculateGammaDelta(LINALG::Matrix<3, 1>& gamma, LINALG::Matrix<8, 1>& delta,
    const LINALG::Matrix<3, 1>& prinv, const LINALG::Matrix<3, 1>& dPI,
    const LINALG::Matrix<6, 1>& ddPII)
{
  // according to Holzapfel-Nonlinear Solid Mechanics p. 216
  gamma(0) = 2. * (dPI(0) + prinv(0) * dPI(1));
  gamma(1) = -2. * dPI(1);
  gamma(2) = 2. * prinv(2) * dPI(2);

  // according to Holzapfel-Nonlinear Solid Mechanics p. 261
  delta(0) = 4. * (ddPII(0) + 2. * prinv(0) * ddPII(5) + dPI(1) + prinv(0) * prinv(0) * ddPII(1));
  delta(1) = -4. * (ddPII(5) + prinv(0) * ddPII(1));
  delta(2) = 4. * (prinv(2) * ddPII(4) + prinv(0) * prinv(2) * ddPII(3));
  delta(3) = 4. * ddPII(1);
  delta(4) = -4. * prinv(2) * ddPII(3);
  delta(5) = 4. * (prinv(2) * dPI(2) + prinv(2) * prinv(2) * ddPII(2));
  delta(6) = -4. * prinv(2) * dPI(2);
  delta(7) = -4. * dPI(1);

  return;
}

/*----------------------------------------------------------------------
 * Evaluates the 2nd Piola-Kirchhoff Stress and constitutive tensor
 * with use of first and second derivatives according to invariants;
 * use principal calculation for all materials, as modified material
 * derivatives are converted before
 *                                                       birzle 11/2014 */
/*----------------------------------------------------------------------*/
void MAT::ElastHyper::EvaluateIsotropicStressCmat(LINALG::Matrix<6, 1>& stress,
    LINALG::Matrix<6, 6>& cmat, const LINALG::Matrix<6, 1>& scg, const LINALG::Matrix<6, 1>& id2,
    const LINALG::Matrix<6, 1>& icg, const LINALG::Matrix<6, 6>& id4sharp,
    const LINALG::Matrix<3, 1>& prinv, const LINALG::Matrix<3, 1>& dPI,
    const LINALG::Matrix<6, 1>& ddPII)
{
  // 2nd Piola Kirchhoff stress factors (according to Holzapfel-Nonlinear Solid Mechanics p. 216)
  LINALG::Matrix<3, 1> gamma(true);
  // constitutive tensor factors (according to Holzapfel-Nonlinear Solid Mechanics p. 261)
  LINALG::Matrix<8, 1> delta(true);

  // compose coefficients
  CalculateGammaDelta(gamma, delta, prinv, dPI, ddPII);

  // 2nd Piola Kirchhoff stress
  stress.Update(gamma(0), id2, 1.0);
  stress.Update(gamma(1), scg, 1.0);
  stress.Update(gamma(2), icg, 1.0);

  // constitutive tensor
  // contribution: Id \otimes Id
  cmat.MultiplyNT(delta(0), id2, id2, 1.0);
  // contribution: Id \otimes C + C \otimes Id
  cmat.MultiplyNT(delta(1), id2, scg, 1.0);
  cmat.MultiplyNT(delta(1), scg, id2, 1.0);
  // contribution: Id \otimes Cinv + Cinv \otimes Id
  cmat.MultiplyNT(delta(2), id2, icg, 1.0);
  cmat.MultiplyNT(delta(2), icg, id2, 1.0);
  // contribution: C \otimes C
  cmat.MultiplyNT(delta(3), scg, scg, 1.0);
  // contribution: C \otimes Cinv + Cinv \otimes C
  cmat.MultiplyNT(delta(4), scg, icg, 1.0);
  cmat.MultiplyNT(delta(4), icg, scg, 1.0);
  // contribution: Cinv \otimes Cinv
  cmat.MultiplyNT(delta(5), icg, icg, 1.0);
  // contribution: Cinv \odot Cinv
  AddtoCmatHolzapfelProduct(cmat, icg, delta(6));
  // contribution: Id4^#
  cmat.Update(delta(7), id4sharp, 1.0);


  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ElastHyper::EvaluateAnisotropicPrinc(LINALG::Matrix<6, 1>& stressanisoprinc,
    LINALG::Matrix<6, 6>& cmatanisoprinc, const LINALG::Matrix<6, 1>& rcg,
    Teuchos::ParameterList& params, const int eleGID)
{
  // loop map of associated potential summands
  for (const Teuchos::RCP<MAT::ELASTIC::Summand>& p : potsum_)
  {
    p->AddStressAnisoPrincipal(rcg, cmatanisoprinc, stressanisoprinc, params, eleGID);
  }

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ElastHyper::EvaluateAnisotropicMod(LINALG::Matrix<6, 1>& stressanisomod,
    LINALG::Matrix<6, 6>& cmatanisomod, const LINALG::Matrix<6, 1>& rcg,
    const LINALG::Matrix<6, 1>& icg, const LINALG::Matrix<3, 1>& prinv, const int eleGID)
{
  // loop map of associated potential summands
  for (const Teuchos::RCP<MAT::ELASTIC::Summand>& p : potsum_)
  {
    p->AddStressAnisoModified(rcg, icg, cmatanisomod, stressanisomod, prinv(2), eleGID);
  }

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ElastHyper::ResponseStretches(LINALG::Matrix<6, 6>& cmat, LINALG::Matrix<6, 1>& stress,
    const LINALG::Matrix<6, 1>& rcg, const bool& havecoeffstrpr, const bool& havecoeffstrmod,
    const int eleGID)
{
  // get principal stretches and directions
  LINALG::Matrix<3, 1> prstr;
  LINALG::Matrix<3, 3> prdir;
  StretchesPrincipal(prstr, prdir, rcg);
  // modified stretches
  LINALG::Matrix<3, 1> modstr;
  StretchesModified(modstr, prstr);
  // determinant of deformation gradient
  const double detdefgrad = prstr(0) * prstr(1) * prstr(2);

  // get coefficients
  LINALG::Matrix<3, 1> gamma_(true);
  LINALG::Matrix<6, 1> delta_(true);
  if (havecoeffstrpr)
  {
    // loop map of associated potential summands
    for (unsigned int p = 0; p < potsum_.size(); ++p)
    {
      potsum_[p]->AddCoefficientsStretchesPrincipal(gamma_, delta_, prstr);
    }
  }
  if (havecoeffstrmod)
  {
    // reciprocal of cubic root of determinant of deformation gradient (convenience)
    const double detdefgrad13 = std::pow(detdefgrad, -1.0 / 3.0);
    // retrieve coefficients with respect to modified principal stretches
    LINALG::Matrix<3, 1> modgamma(true);
    LINALG::Matrix<6, 1> moddelta(true);
    {
      // loop map of associated potential summands
      for (unsigned int p = 0; p < potsum_.size(); ++p)
      {
        potsum_[p]->AddCoefficientsStretchesModified(modgamma, moddelta, modstr);
      }
    }
    // convert modified coefficients to oridinary counterparts
    //
    // derivatives of modified pr. stretches WRT pr. stretches
    LINALG::Matrix<3, 3> modbypr(false);
    for (int al = 0; al < 3; ++al)
    {
      for (int be = 0; be < 3; ++be)
      {
        modbypr(al, be) = -modstr(al) / modstr(be);
      }
      modbypr(al, al) += 3.0;
    }
    modbypr.Scale(detdefgrad13 / 3.0);
    // determine unmodified coefficients gamma and add them
    gamma_.MultiplyTN(1.0, modbypr, modgamma, 1.0);
    // determine unmodified coefficients delta and add them
    //
    // rewrite mod.coeff. as 2-tensor
    LINALG::Matrix<3, 3> moddeltat(false);
    moddeltat(0, 0) = moddelta(0);
    moddeltat(1, 1) = moddelta(1);
    moddeltat(2, 2) = moddelta(2);
    moddeltat(0, 1) = moddeltat(1, 0) = moddelta(3);
    moddeltat(1, 2) = moddeltat(2, 1) = moddelta(4);
    moddeltat(2, 0) = moddeltat(0, 2) = moddelta(5);
    // Psi_{,barlam barlam} barlam_{,lam} barlam_{,lam}
    LINALG::Matrix<3, 3> aux(false);
    aux.MultiplyTN(modbypr, moddeltat);
    LINALG::Matrix<3, 3> deltat(false);
    deltat.MultiplyNN(aux, modbypr);
    // Psi_{,barlam} barlam_{,lam lam}
    for (int be = 0; be < 3; ++be)
    {
      for (int ga = 0; ga < 3; ++ga)
      {
        double deltat_bega = 0.0;
        for (int al = 0; al < 3; ++al)
        {
          deltat_bega += -modgamma(al) * modbypr(al, be) / (3.0 * prstr(ga));
          if (ga == al) deltat_bega += -modgamma(al) * detdefgrad13 / (3.0 * prstr(be));
          if (be == ga)
            deltat_bega += modgamma(al) * detdefgrad13 * prstr(al) / (3.0 * prstr(be) * prstr(be));
        }
        deltat(be, ga) += deltat_bega;
      }
    }
    // add to delta
    // Psi_{lam lam} = Psi_{,barlam barlam} barlam_{,lam} barlam_{,lam}
    //               + Psi_{,barlam} barlam_{,lam lam}
    delta_(0) += deltat(0, 0);
    delta_(1) += deltat(1, 1);
    delta_(2) += deltat(2, 2);
    delta_(3) += deltat(0, 1);
    delta_(4) += deltat(1, 2);
    delta_(5) += deltat(2, 0);
  }

  // principal 2nd Piola--Kirchhoff stress tensor, cf [1] Eq (6.47)
  LINALG::Matrix<3, 1> prsts(true);
  for (int al = 0; al < 3; ++al)
  {
    // PK2 principal stresses
    prsts(al) = gamma_(al) / prstr(al);
    // PK2 tensor in Voigt notation
    stress(0) += prsts(al) * prdir(0, al) * prdir(0, al);  // S^11
    stress(1) += prsts(al) * prdir(1, al) * prdir(1, al);  // S^22
    stress(2) += prsts(al) * prdir(2, al) * prdir(2, al);  // S^33
    stress(3) += prsts(al) * prdir(0, al) * prdir(1, al);  // S^12
    stress(4) += prsts(al) * prdir(1, al) * prdir(2, al);  // S^23
    stress(5) += prsts(al) * prdir(2, al) * prdir(0, al);  // S^31
  }

  // integration factor prfact_{al be}
  LINALG::Matrix<6, 1> prfact1(true);
  LINALG::Matrix<6, 1> prfact2(true);
  for (int albe = 0; albe < 6; ++albe)
  {
    const int al = VOIGT6ROW_[albe];
    const int be = VOIGT6COL_[albe];
    double prfact1_albe = delta_(albe) / (prstr(al) * prstr(be));
    if (albe < 3) prfact1_albe -= gamma_(al) / (prstr(be) * prstr(al) * prstr(al));
    prfact1(albe) = prfact1_albe;
    if (al != be)
    {
      if (fabs(prstr(al) - prstr(be)) < EPS6)
        prfact2(albe) = (prfact1(be) - prfact1(albe)) / 2.0;
      else
        prfact2(albe) = (prsts(be) - prsts(al)) / (prstr(be) * prstr(be) - prstr(al) * prstr(al));
    }
  }

  // add elasticity 4-tensor, cf Holzapfel [1] Eq (6.180),(6.196)
  for (int kl = 0; kl < 6; ++kl)
  {
    const int k = VOIGT6ROW_[kl];
    const int l = VOIGT6COL_[kl];
    for (int ij = 0; ij < 6; ++ij)
    {
      const int i = VOIGT6ROW_[ij];
      const int j = VOIGT6COL_[ij];
      double c_ijkl = 0.0;
      for (int albe = 0; albe < 6; ++albe)
      {
        const int al = VOIGT6ROW_[albe];
        const int be = VOIGT6COL_[albe];
        const double fact1 = prfact1(albe);
        c_ijkl += fact1 * prdir(i, al) * prdir(j, al) * prdir(k, be) * prdir(l, be);
        if (albe >= 3)
        {  // al!=be
          c_ijkl += fact1 * prdir(i, be) * prdir(j, be) * prdir(k, al) * prdir(l, al);
          const double fact2 = prfact2(albe);
          c_ijkl += fact2 * prdir(i, al) * prdir(j, be) * prdir(k, al) * prdir(l, be) +
                    fact2 * prdir(i, al) * prdir(j, be) * prdir(k, be) * prdir(l, al) +
                    fact2 * prdir(i, be) * prdir(j, al) * prdir(k, be) * prdir(l, al) +
                    fact2 * prdir(i, be) * prdir(j, al) * prdir(k, al) * prdir(l, be);
        }
      }
      cmat(ij, kl) += c_ijkl;
    }
  }
  // ready
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ElastHyper::VisNames(std::map<std::string, int>& names)
{
  if (AnisotropicPrincipal() or AnisotropicModified())
  {
    std::vector<LINALG::Matrix<3, 1>> fibervecs;
    GetFiberVecs(fibervecs);
    int vissize = fibervecs.size();
    std::string fiber;
    for (int i = 0; i < vissize; i++)
    {
      std::ostringstream s;
      s << "Fiber" << i + 1;
      fiber = s.str();
      names[fiber] = 3;  // 3-dim vector
    }
  }
  // do visualization for isotropic materials as well
  // loop map of associated potential summands
  for (unsigned int p = 0; p < potsum_.size(); ++p)
  {
    potsum_[p]->VisNames(names);
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool MAT::ElastHyper::VisData(
    const std::string& name, std::vector<double>& data, int numgp, int eleID)
{
  //
  int return_val = 0;
  if (AnisotropicPrincipal() or AnisotropicModified())
  {
    std::vector<LINALG::Matrix<3, 1>> fibervecs;
    GetFiberVecs(fibervecs);
    int vissize = fibervecs.size();
    for (int i = 0; i < vissize; i++)
    {
      std::ostringstream s;
      s << "Fiber" << i + 1;
      std::string fiber;
      fiber = s.str();
      if (name == fiber)
      {
        if ((int)data.size() != 3) dserror("size mismatch");
        data[0] = fibervecs.at(i)(0);
        data[1] = fibervecs.at(i)(1);
        data[2] = fibervecs.at(i)(2);
      }
    }
    // return true;
    return_val = 1;
  }
  // loop map of associated potential summands
  for (unsigned int p = 0; p < potsum_.size(); ++p)
  {
    return_val += potsum_[p]->VisData(name, data, numgp, eleID);
  }
  return (bool)return_val;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const MAT::ELASTIC::Summand> MAT::ElastHyper::GetPotSummandPtr(
    const INPAR::MAT::MaterialType& materialtype) const
{
  for (unsigned int p = 0; p < potsum_.size(); ++p)
  {
    if (potsum_[p]->MaterialType() == materialtype) return potsum_[p];
  }
  return Teuchos::null;
}

/*----------------------------------------------------------------------*/
/* Check polyconvexity:
 * Polyconvexity of isotropic hyperelastic material
 * (dependent on principal or modified invariants)
 * is tested with eq. (5.31) of Vera Ebbing - PHD-thesis (p. 79).
 * Partial derivatives of SEF are used for calculation of eq (5.31).
 *                                                        birzle 04/2016 */
/*----------------------------------------------------------------------*/
void MAT::ElastHyper::CheckPolyconvexity(const LINALG::Matrix<3, 3>& defgrd,
    const LINALG::Matrix<3, 1>& prinv, const LINALG::Matrix<3, 1>& dPI,
    const LINALG::Matrix<6, 1>& ddPII, Teuchos::ParameterList& params, const int eleGID)
{
  // This polyconvexity-test is just implemented for isotropic hyperelastic-materials
  // --> error if anisotropic material is tested (plastic and viscoelastic materials should not get
  // in here)
  if (anisoprinc_ || anisomod_)
    dserror(
        "This polyconvexity-check is just implemented for isotropic "
        "hyperelastic-materials (do not use for anistropic materials).");

  // principal invariants (i)
  // first strain energy derivative dPI (i)
  // second strain energy derivative ddPII (i)

  // J = sqrt(I_3) = modinv(2)
  double J = std::pow(prinv(2), 1. / 2.);

  // defgrd = F (i)
  // dfgrd = F in Voigt - Notation
  LINALG::Matrix<9, 1> dfgrd(true);
  dfgrd(0, 0) = defgrd(0, 0);
  dfgrd(1, 0) = defgrd(1, 1);
  dfgrd(2, 0) = defgrd(2, 2);
  dfgrd(3, 0) = defgrd(0, 1);
  dfgrd(4, 0) = defgrd(1, 2);
  dfgrd(5, 0) = defgrd(0, 2);
  dfgrd(6, 0) = defgrd(1, 0);
  dfgrd(7, 0) = defgrd(2, 1);
  dfgrd(8, 0) = defgrd(2, 0);

  // Cof(F) = J*F^(-T)
  LINALG::Matrix<3, 3> CoFacF(true);  // Cof(F) in Matrix-Notation
  LINALG::Matrix<9, 1> CofF(true);    // Cof(F) in Voigt-Notation
  CoFacF.Invert(defgrd);
  CoFacF.Scale(J);
  // sort in Voigt-Notation and invert!
  CofF(0, 0) = CoFacF(0, 0);
  CofF(1, 0) = CoFacF(1, 1);
  CofF(2, 0) = CoFacF(2, 2);
  CofF(6, 0) = CoFacF(0, 1);
  CofF(7, 0) = CoFacF(1, 2);
  CofF(8, 0) = CoFacF(0, 2);
  CofF(3, 0) = CoFacF(1, 0);
  CofF(4, 0) = CoFacF(2, 1);
  CofF(5, 0) = CoFacF(2, 0);

  // id4 (9x9)
  LINALG::Matrix<9, 9> ID4(true);
  for (int i = 0; i < 9; i++)
    for (int j = 0; j < 9; j++)
      if (i == j) ID4(i, j) = 1.0;

  // Frechet Derivative according to Ebbing, PhD-thesis page 79, Eq: (5.31)
  LINALG::Matrix<19, 19> FreD(true);

  // single matrices of Frechet Derivative:

  // d²P/dFdF
  // = 4 d^2\Psi/dI_1dI_1 F \otimes F + 2 \d\Psi/dI_1 *II
  LINALG::Matrix<9, 9> FreDFF(true);
  FreDFF.MultiplyNT(4 * ddPII(0), dfgrd, dfgrd, 1.0);
  FreDFF.Update(2 * dPI(0), ID4, 1.0);

  // d²P/d(cofF)d(cofF)
  // = = 4 d^2\Psi/dI_2dI_2 cof(F) \otimes cof(F) + 2 \d\Psi/dI_2 *II
  LINALG::Matrix<9, 9> FreDcFcF(true);
  FreDcFcF.MultiplyNT(4 * ddPII(1), CofF, CofF, 1.0);
  FreDcFcF.Update(2 * dPI(1), ID4, 1.0);

  // d²P/d(detF)d(detF)
  // = 2*d \Psi/dI_3 + 4*I_3*d²\Psi/dI_3dI_3
  double FreDJJ(true);
  FreDJJ += 2 * dPI(2) + 4 * prinv(2) * ddPII(2);

  // d²P/d(cofF)dF
  // = 4*d\Psi/dI_1dI_2 F /otimes CofF
  LINALG::Matrix<9, 9> FreDcFF(true);
  FreDcFF.MultiplyNT(4 * ddPII(5), dfgrd, CofF, 1.0);

  // d²P/d(detF)d(cofF)
  // = 4*J*d^2 \Psi /dI_2 dI_3 \mat{CofF}
  LINALG::Matrix<9, 1> FreDcFJ(true);
  FreDcFJ.Update(4 * J * ddPII(3), CofF, 1.0);

  // d²P/d(detF) dF = d²P/dF d(detF)
  // = 4*J*d^2 \Psi /dI_1 dI_3 \mat{F}
  LINALG::Matrix<9, 1> FreDFJ(true);
  FreDFJ.Update(4 * J * ddPII(4), dfgrd, 1.0);

  // Sort values in Frechet Derivative

  // FreD = [FreDFF   FreDcFF    FreDFJ
  //         FreDcFF  FreDcFcF   FreDcFJ
  //         FreDFJ   FreDcFJ    FreDJJ]
  for (int i = 0; i < 9; i++)
    for (int j = 0; j < 9; j++)
    {
      FreD(i, j) = FreDFF(i, j);
      FreD(i, j + 9) = FreDcFF(i, j);
      FreD(i + 9, j) = FreDcFF(i, j);
      FreD(i + 9, j + 9) = FreDcFcF(i, j);
    }

  for (int i = 0; i < 9; i++)
  {
    FreD(i + 9, 18) = FreDcFJ(i);
    FreD(18, i + 9) = FreDcFJ(i);
    FreD(i, 18) = FreDFJ(i);
    FreD(18, i) = FreDFJ(i);
  }

  FreD(18, 18) = FreDJJ;

  // EigenValues of Frechet Derivative
  LINALG::Matrix<19, 19> EWFreD(true);  // EW on diagonal
  LINALG::Matrix<19, 19> EVFreD(true);
  LINALG::SYEV(FreD, EWFreD, EVFreD);

  // Just positive EigenValues --> System is polyconvex
  for (int i = 0; i < 19; i++)
    for (int j = 0; j < 19; j++)
      if (i == j)  // values on diagonal = EigenValues
        if (EWFreD(i, i) <
            (-1.0e-10 * EWFreD.NormInf()))  // do not test < 0, but reasonable small value
        {
          std::cout << "\nWARNING: Your system is not polyconvex!" << std::endl;
          std::cout << "Polyconvexity fails at: Element-Id: " << eleGID
                    << " and Gauß-Point: " << params.get<int>("gp") << std::endl;
          std::cout << "Eigenvalues of the Frechet Derivative are: " << EWFreD << std::endl;
        }

  return;
}


/*----------------------------------------------------------------------*/
/* Fit parameters of elasthyper materials in
 * stat inverse analysis
 *                                                        birzle 05/2017 */
/*----------------------------------------------------------------------*/
void MAT::ElastHyper::ElastOptParams(std::map<std::string, int>* pnames)
{
  // loop map of associated potential summands
  for (unsigned int p = 0; p < potsum_.size(); ++p)
  {
    potsum_[p]->AddElastOptParams(pnames);
  }
  return;
}
