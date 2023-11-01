/*----------------------------------------------------------------------*/
/*! \file
\brief
This file contains the hyperelastic toolbox. It allows summing up several summands
of several types (isotropic or anisotropic, splitted or not) to build a hyperelastic
strain energy function.

The input line should read
MAT 0   MAT_ElastHyper   NUMMAT 2 MATIDS 1 2 DENS 0

\level 1


*/
/*----------------------------------------------------------------------*/

#include "baci_mat_elasthyper.H"

#include "baci_lib_globalproblem.H"
#include "baci_linalg_fixedsizematrix_voigt_notation.H"
#include "baci_mat_elasthyper_service.H"
#include "baci_mat_par_bundle.H"
#include "baci_mat_service.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::PAR::ElastHyper::ElastHyper(const Teuchos::RCP<MAT::PAR::Material>& matdata)
    : Parameter(matdata),
      nummat_(matdata->GetInt("NUMMAT")),
      matids_(matdata->Get<std::vector<int>>("MATIDS")),
      density_(matdata->GetDouble("DENS")),
      polyconvex_(matdata->GetInt("POLYCONVEX"))

{
  // check if sizes fit
  if (nummat_ != (int)matids_->size())
    dserror("number of materials %d does not fit to size of material vector %d", nummat_,
        matids_->size());

  // output, that polyconvexity is checked
  if (polyconvex_ != 0) std::cout << "Polyconvexity of your simulation is checked." << std::endl;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<MAT::Material> MAT::PAR::ElastHyper::CreateMaterial()
{
  return Teuchos::rcp(new MAT::ElastHyper(this));
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::ElastHyperType MAT::ElastHyperType::instance_;


DRT::ParObject* MAT::ElastHyperType::Create(const std::vector<char>& data)
{
  auto* elhy = new MAT::ElastHyper();
  elhy->Unpack(data);

  return elhy;
}


/*----------------------------------------------------------------------*
 |  initialise static arrays                                 bborn 08/09|
 *----------------------------------------------------------------------*/



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::ElastHyper::ElastHyper() : summandProperties_(), params_(nullptr), potsum_(0), anisotropy_() {}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::ElastHyper::ElastHyper(MAT::PAR::ElastHyper* params)
    : summandProperties_(), params_(params), potsum_(0), anisotropy_()
{
  // make sure the referenced materials in material list have quick access parameters
  std::vector<int>::const_iterator m;
  for (m = params_->matids_->begin(); m != params_->matids_->end(); ++m)
  {
    const int matid = *m;
    Teuchos::RCP<MAT::ELASTIC::Summand> sum = MAT::ELASTIC::Summand::Factory(matid);
    if (sum == Teuchos::null) dserror("Failed to allocate");
    potsum_.push_back(sum);
    sum->RegisterAnisotropyExtensions(anisotropy_);
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
  if (params_ != nullptr) matid = params_->Id();  // in case we are in post-process mode
  AddtoPack(data, matid);
  summandProperties_.Pack(data);

  anisotropy_.PackAnisotropy(data);

  if (params_ != nullptr)  // summands are not accessible in postprocessing mode
  {
    // loop map of associated potential summands
    for (const auto& p : potsum_)
    {
      p->PackSummand(data);
    }
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ElastHyper::Unpack(const std::vector<char>& data)
{
  // make sure we have a pristine material
  params_ = nullptr;
  potsum_.clear();

  std::vector<char>::size_type position = 0;

  CORE::COMM::ExtractAndAssertId(position, data, UniqueParObjectId());

  // matid and recover params_
  int matid;
  ExtractfromPack(position, data, matid);
  if (DRT::Problem::Instance()->Materials() != Teuchos::null)
  {
    if (DRT::Problem::Instance()->Materials()->Num() != 0)
    {
      const int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();
      MAT::PAR::Parameter* mat =
          DRT::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        params_ = dynamic_cast<MAT::PAR::ElastHyper*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }
  }

  summandProperties_.Unpack(position, data);

  // Pack anisotropy
  anisotropy_.UnpackAnisotropy(data, position);

  if (params_ != nullptr)  // summands are not accessible in postprocessing mode
  {
    // make sure the referenced materials in material list have quick access parameters
    std::vector<int>::const_iterator m;
    for (m = params_->matids_->begin(); m != params_->matids_->end(); ++m)
    {
      const int summand_matid = *m;
      Teuchos::RCP<MAT::ELASTIC::Summand> sum = MAT::ELASTIC::Summand::Factory(summand_matid);
      if (sum == Teuchos::null) dserror("Failed to allocate");
      potsum_.push_back(sum);
    }

    // loop map of associated potential summands
    for (auto& p : potsum_)
    {
      p->UnpackSummand(data, position);
      p->RegisterAnisotropyExtensions(anisotropy_);
    }

    // in the postprocessing mode, we do not unpack everything we have packed
    // -> position check cannot be done in this case
    if (position != data.size())
    {
      dserror("Mismatch in size of data %d <-> %d", data.size(), position);
    }
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
int MAT::ElastHyper::MatID(const unsigned index) const
{
  if ((int)index >= params_->nummat_)
  {
    dserror("Index too large");
  }

  return params_->matids_->at(index);
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
    for (const auto& p : potsum_)
    {
      p->AddShearMod(haveshearmod, shearmod);
    }
  }
  if (!haveshearmod)
  {
    dserror("Cannot provide shear modulus equivalent");
  }
  return shearmod;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double MAT::ElastHyper::GetYoung()
{
  double young;
  double shear;
  double bulk;
  young = shear = bulk = 0.;
  for (auto& p : potsum_) p->AddYoungsMod(young, shear, bulk);

  if (bulk != 0. || shear != 0.) young += 9. * bulk * shear / (3. * bulk + shear);

  return young;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ElastHyper::SetupAAA(Teuchos::ParameterList& params, const int eleGID)
{
  // loop map of associated potential summands
  for (auto& p : potsum_)
  {
    p->SetupAAA(params, eleGID);
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ElastHyper::Setup(int numgp, DRT::INPUT::LineDefinition* linedef)
{
  // Read anisotropy
  anisotropy_.SetNumberOfGaussPoints(numgp);
  anisotropy_.ReadAnisotropyFromElement(linedef);

  // Setup summands
  for (auto& p : potsum_)
  {
    p->Setup(numgp, linedef);
  }
  summandProperties_.Clear();
  ElastHyperProperties(potsum_, summandProperties_);

  if (summandProperties_.viscoGeneral)
  {
    dserror(
        "Never use viscoelastic-materials in Elasthyper-Toolbox. Use Viscoelasthyper-Toolbox "
        "instead.");
  }
}

void MAT::ElastHyper::PostSetup(Teuchos::ParameterList& params, const int eleGID)
{
  anisotropy_.ReadAnisotropyFromParameterList(params);

  // Forward PostSetup call to all summands
  for (auto& p : potsum_)
  {
    p->PostSetup(params);
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ElastHyper::Update()
{
  // loop map of associated potential summands
  for (auto& p : potsum_)
  {
    p->Update();
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ElastHyper::GetFiberVecs(std::vector<CORE::LINALG::Matrix<3, 1>>& fibervecs)
{
  if (summandProperties_.anisoprinc || summandProperties_.anisomod)
  {
    for (auto& p : potsum_)
    {
      p->GetFiberVecs(fibervecs);
    }
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ElastHyper::EvaluateFiberVecs(const double newgamma,
    const CORE::LINALG::Matrix<3, 3>& locsys, const CORE::LINALG::Matrix<3, 3>& defgrd)
{
  if (summandProperties_.anisoprinc || summandProperties_.anisomod)
  {
    for (auto& p : potsum_)
    {
      p->SetFiberVecs(newgamma, locsys, defgrd);
    }
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ElastHyper::StrainEnergy(
    const CORE::LINALG::Matrix<6, 1>& glstrain, double& psi, const int gp, const int eleGID)
{
  static CORE::LINALG::Matrix<6, 1> C_strain(true);
  C_strain.Clear();
  static CORE::LINALG::Matrix<3, 1> prinv(true);
  prinv.Clear();
  static CORE::LINALG::Matrix<3, 1> modinv(true);
  modinv.Clear();

  EvaluateRightCauchyGreenStrainLikeVoigt(glstrain, C_strain);
  CORE::LINALG::VOIGT::Strains::InvariantsPrincipal(prinv, C_strain);
  InvariantsModified(modinv, prinv);

  // loop map of associated potential summands
  for (const auto& p : potsum_)
  {
    p->AddStrainEnergy(psi, prinv, modinv, glstrain, gp, eleGID);
  }
}


/*----------------------------------------------------------------------*
 |  Evaluate for GEMM time integration                        popp 11/13|
 *----------------------------------------------------------------------*/
void MAT::ElastHyper::EvaluateGEMM(CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, 1>* stress,
    CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D>* cmat,
    const CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, 1>* glstrain_m,
    const CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, 1>* glstrain_new,
    const CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, 1>* glstrain_old,
    const CORE::LINALG::Matrix<3, 3>* rcg_new, const CORE::LINALG::Matrix<3, 3>* rcg_old,
    const int gp, const int eleGID)
{
#ifdef DEBUG
  if (stress == nullptr) dserror("No stress vector supplied");
  if (cmat == nullptr) dserror("No material tangent matrix supplied");
  if (glstrain_m == nullptr) dserror("No GL strains supplied");
  if (glstrain_new == nullptr) dserror("No GL strains supplied");
  if (glstrain_old == nullptr) dserror("No GL strains supplied");
#endif

  // standard material evaluate call at midpoint t_{n+1/2}
  Teuchos::ParameterList params;
  CORE::LINALG::Matrix<3, 3> defgrd(true);
  Evaluate(&defgrd, glstrain_m, params, stress, cmat, gp, eleGID);

  //**********************************************************************
  // CHECK IF GEMM ALGORITHMIC STRESSES NEED TO BE APPLIED
  //**********************************************************************
  // increment of Cauchy-Green tensor in Voigt notation
  CORE::LINALG::Matrix<6, 1> M;
  M(0) = (*rcg_new)(0, 0) - (*rcg_old)(0, 0);
  M(1) = (*rcg_new)(1, 1) - (*rcg_old)(1, 1);
  M(2) = (*rcg_new)(2, 2) - (*rcg_old)(2, 2);
  M(3) = (*rcg_new)(0, 1) + (*rcg_new)(1, 0) - (*rcg_old)(0, 1) - (*rcg_old)(1, 0);
  M(4) = (*rcg_new)(1, 2) + (*rcg_new)(2, 1) - (*rcg_old)(1, 2) - (*rcg_old)(2, 1);
  M(5) = (*rcg_new)(0, 2) + (*rcg_new)(2, 0) - (*rcg_old)(0, 2) - (*rcg_old)(2, 0);

  // second variant of M in Voigt notation
  CORE::LINALG::Matrix<6, 1> Mtilde;
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
  CORE::LINALG::Matrix<6, 1> vecid(true);
  for (int k = 0; k < 6; ++k) vecid(k) = 1.0;
  CORE::LINALG::Matrix<6, 6> halfid(true);
  for (int k = 0; k < 3; ++k) halfid(k, k) = 1.0;
  for (int k = 3; k < 6; ++k) halfid(k, k) = 0.5;

  // strain energy function at t_{n+1} and t_{n}
  double psi = 0.0;
  double psio = 0.0;
  StrainEnergy(*glstrain_new, psi, gp, eleGID);
  StrainEnergy(*glstrain_old, psio, gp, eleGID);

  // derivative of strain energy function dpsi = 0.5*stress
  // double contraction dpsi : M
  double dpsiM = 0.5 * (*stress)(0) * M(0) + 0.5 * (*stress)(1) * M(1) + 0.5 * (*stress)(2) * M(2) +
                 0.5 * (*stress)(3) * M(3) + 0.5 * (*stress)(4) * M(4) + 0.5 * (*stress)(5) * M(5);

  // factor for algorithmic stresses
  double fac = 2.0 * ((psi - psio - dpsiM) / Mb);

  // algorithmic stresses
  CORE::LINALG::Matrix<6, 1> algstress(true);
  algstress.Update(fac, Mtilde, 1.0);

  //**********************************************************************
  // COMPUTE GEMM ALGORITHMIC MATERIAL TENSOR
  //**********************************************************************
  // algorithmic material tensor requires stresses at t_{n+1}
  CORE::LINALG::Matrix<6, 1> stressnew(true);
  CORE::LINALG::Matrix<6, 6> cmatnew(true);
  Evaluate(&defgrd, glstrain_new, params, &stressnew, &cmatnew, gp, eleGID);

  // initialize algorithmic material tensor
  CORE::LINALG::Matrix<6, 6> algcmat(true);

  // part 1 (derivative of Mtilde)
  algcmat.Update(4.0 * fac, halfid, 1.0);

  // part 2a (derivative of strain energy in fac)
  CORE::LINALG::Matrix<6, 1> dfac(true);
  dfac.Update(2.0 / Mb, stressnew, 1.0);

  // part 2b (derivative of dpsiM in fac)
  CORE::LINALG::Matrix<6, 1> tmp(true);
  tmp.Multiply(*cmat, M);
  dfac.Update(-0.5 / Mb, tmp, 1.0);
  dfac.Update(-2.0 / Mb, *stress, 1.0);

  // part 2c (derivative of Mb in fac)
  tmp.Multiply(halfid, M);
  dfac.Update(-4.0 * (psi - psio - dpsiM) / (Mb * Mb), tmp, 1.0);
  dfac.Update(-4.0 * (psi - psio - dpsiM) / (Mb * Mb), Mtilde, 1.0);

  // part 2 (derivative of fac, put together parts 2a,2b and 2c)
  CORE::LINALG::Matrix<6, 6> tmpmat(true);
  tmpmat.MultiplyNT(2.0, Mtilde, dfac);
  algcmat.Update(1.0, tmpmat, 1.0);

  //**********************************************************************
  // EXTEND ORIGINAL STRESSES / CMAT WITH GEMM CONTRIBUTIONS
  //**********************************************************************
  stress->Update(1.0, algstress, 1.0);
  cmat->Update(1.0, algcmat, 1.0);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ElastHyper::Evaluate(const CORE::LINALG::Matrix<3, 3>* defgrd,
    const CORE::LINALG::Matrix<6, 1>* glstrain, Teuchos::ParameterList& params,
    CORE::LINALG::Matrix<6, 1>* stress, CORE::LINALG::Matrix<6, 6>* cmat, const int gp,
    const int eleGID)
{
  bool checkpolyconvexity = (params_ != nullptr and params_->polyconvex_ != 0);

  ElastHyperEvaluate(*defgrd, *glstrain, params, *stress, *cmat, gp, eleGID, potsum_,
      summandProperties_, checkpolyconvexity);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ElastHyper::EvaluateCauchyDerivs(const CORE::LINALG::Matrix<3, 1>& prinv, const int gp,
    int eleGID, CORE::LINALG::Matrix<3, 1>& dPI, CORE::LINALG::Matrix<6, 1>& ddPII,
    CORE::LINALG::Matrix<10, 1>& dddPIII, const double* temp)
{
  for (auto& i : potsum_)
  {
    if (summandProperties_.isoprinc)
    {
      i->AddDerivativesPrincipal(dPI, ddPII, prinv, gp, eleGID);
      i->AddThirdDerivativesPrincipalIso(dddPIII, prinv, gp, eleGID);
    }
    if (summandProperties_.isomod || summandProperties_.anisomod || summandProperties_.anisoprinc)
      dserror("not implemented for this form of strain energy function");
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ElastHyper::EvaluateCauchyNDirAndDerivatives(const CORE::LINALG::Matrix<3, 3>& defgrd,
    const CORE::LINALG::Matrix<3, 1>& n, const CORE::LINALG::Matrix<3, 1>& dir,
    double& cauchy_n_dir, CORE::LINALG::Matrix<3, 1>* d_cauchyndir_dn,
    CORE::LINALG::Matrix<3, 1>* d_cauchyndir_ddir, CORE::LINALG::Matrix<9, 1>* d_cauchyndir_dF,
    CORE::LINALG::Matrix<9, 9>* d2_cauchyndir_dF2, CORE::LINALG::Matrix<9, 3>* d2_cauchyndir_dF_dn,
    CORE::LINALG::Matrix<9, 3>* d2_cauchyndir_dF_ddir, int gp, int eleGID,
    const double* concentration, const double* temp, double* d_cauchyndir_dT,
    CORE::LINALG::Matrix<9, 1>* d2_cauchyndir_dF_dT)
{
  cauchy_n_dir = 0.0;

  static CORE::LINALG::Matrix<3, 3> b(true);
  b.MultiplyNT(1.0, defgrd, defgrd, 0.0);
  static CORE::LINALG::Matrix<3, 1> bdn(true);
  bdn.Multiply(1.0, b, n, 0.0);
  static CORE::LINALG::Matrix<3, 1> bddir(true);
  bddir.Multiply(1.0, b, dir, 0.0);
  const double bdnddir = bdn.Dot(dir);

  static CORE::LINALG::Matrix<3, 3> ib(true);
  ib.Invert(b);
  static CORE::LINALG::Matrix<3, 1> ibdn(true);
  ibdn.Multiply(1.0, ib, n, 0.0);
  static CORE::LINALG::Matrix<3, 1> ibddir(true);
  ibddir.Multiply(1.0, ib, dir, 0.0);
  const double ibdnddir = ibdn.Dot(dir);
  const double nddir = n.Dot(dir);

  static CORE::LINALG::Matrix<6, 1> bV_strain(true);
  CORE::LINALG::VOIGT::Strains::MatrixToVector(b, bV_strain);
  static CORE::LINALG::Matrix<3, 1> prinv(true);
  CORE::LINALG::VOIGT::Strains::InvariantsPrincipal(prinv, bV_strain);

  static CORE::LINALG::Matrix<3, 1> dPI(true);
  static CORE::LINALG::Matrix<6, 1> ddPII(true);
  static CORE::LINALG::Matrix<10, 1> dddPIII(true);
  dPI.Clear();
  ddPII.Clear();
  dddPIII.Clear();
  EvaluateCauchyDerivs(prinv, gp, eleGID, dPI, ddPII, dddPIII, temp);

  const double prefac = 2.0 / std::sqrt(prinv(2));

  cauchy_n_dir = prefac * (prinv(1) * dPI(1) * nddir + prinv(2) * dPI(2) * nddir +
                              dPI(0) * bdnddir - prinv(2) * dPI(1) * ibdnddir);

  if (d_cauchyndir_dn != nullptr)
  {
    d_cauchyndir_dn->Update(prinv(1) * dPI(1) + prinv(2) * dPI(2), dir, 0.0);
    d_cauchyndir_dn->Update(dPI(0), bddir, 1.0);
    d_cauchyndir_dn->Update(-prinv(2) * dPI(1), ibddir, 1.0);
    d_cauchyndir_dn->Scale(prefac);
  }

  if (d_cauchyndir_ddir != nullptr)
  {
    d_cauchyndir_ddir->Update(prinv(1) * dPI(1) + prinv(2) * dPI(2), n, 0.0);
    d_cauchyndir_ddir->Update(dPI(0), bdn, 1.0);
    d_cauchyndir_ddir->Update(-prinv(2) * dPI(1), ibdn, 1.0);
    d_cauchyndir_ddir->Scale(prefac);
  }

  // calculate stuff that is needed for evaluations of derivatives w.r.t. F
  static CORE::LINALG::Matrix<9, 1> FV(true);
  CORE::LINALG::VOIGT::Matrix3x3to9x1(defgrd, FV);
  static CORE::LINALG::Matrix<3, 3> iF(true);
  iF.Invert(defgrd);
  static CORE::LINALG::Matrix<3, 3> iFT(true);
  iFT.UpdateT(iF);
  static CORE::LINALG::Matrix<9, 1> iFTV(true);
  CORE::LINALG::VOIGT::Matrix3x3to9x1(iFT, iFTV);

  // calculation of dI_i/dF (derivatives of invariants of b w.r.t. deformation gradient)
  static CORE::LINALG::Matrix<3, 3> bdF(true);
  bdF.Multiply(1.0, b, defgrd, 0.0);
  static CORE::LINALG::Matrix<9, 1> bdFV(true);
  CORE::LINALG::VOIGT::Matrix3x3to9x1(bdF, bdFV);
  static CORE::LINALG::Matrix<3, 3> ibdF(true);
  ibdF.Multiply(1.0, ib, defgrd, 0.0);
  static CORE::LINALG::Matrix<9, 1> ibdFV(true);
  CORE::LINALG::VOIGT::Matrix3x3to9x1(ibdF, ibdFV);
  static CORE::LINALG::Matrix<9, 1> d_I1_dF(true);
  d_I1_dF.Update(2.0, FV, 0.0);
  static CORE::LINALG::Matrix<9, 1> d_I2_dF(true);
  d_I2_dF.Update(prinv(0), FV, 0.0);
  d_I2_dF.Update(-1.0, bdFV, 1.0);
  d_I2_dF.Scale(2.0);
  static CORE::LINALG::Matrix<9, 1> d_I3_dF(true);
  d_I3_dF.Update(2.0 * prinv(2), ibdFV, 0.0);

  // calculate d(b \cdot n \cdot t)/dF
  static CORE::LINALG::Matrix<3, 1> tempvec3x1(true);
  static CORE::LINALG::Matrix<1, 3> tempvec1x3(true);
  tempvec1x3.MultiplyTN(1.0, dir, defgrd, 0.0);
  static CORE::LINALG::Matrix<3, 3> d_bdnddir_dF(true);
  d_bdnddir_dF.MultiplyNN(1.0, n, tempvec1x3, 0.0);
  tempvec1x3.MultiplyTN(1.0, n, defgrd, 0.0);
  d_bdnddir_dF.MultiplyNN(1.0, dir, tempvec1x3, 1.0);
  static CORE::LINALG::Matrix<9, 1> d_bdnddir_dFV(true);
  CORE::LINALG::VOIGT::Matrix3x3to9x1(d_bdnddir_dF, d_bdnddir_dFV);

  // calculate d(b^{-1} \cdot n \cdot t)/dF
  static CORE::LINALG::Matrix<1, 3> dirdibdF(true);
  static CORE::LINALG::Matrix<1, 3> ndibdF(true);
  dirdibdF.MultiplyTN(1.0, dir, ibdF, 0.0);
  static CORE::LINALG::Matrix<3, 3> d_ibdnddir_dF(true);
  d_ibdnddir_dF.MultiplyNN(1.0, ibdn, dirdibdF, 0.0);
  ndibdF.MultiplyTN(1.0, n, ibdF, 0.0);
  d_ibdnddir_dF.MultiplyNN(1.0, ibddir, ndibdF, 1.0);
  d_ibdnddir_dF.Scale(-1.0);
  static CORE::LINALG::Matrix<9, 1> d_ibdnddir_dFV(true);
  CORE::LINALG::VOIGT::Matrix3x3to9x1(d_ibdnddir_dF, d_ibdnddir_dFV);

  if (temp != nullptr)
    EvaluateCauchyTempDeriv(prinv, nddir, bdnddir, ibdnddir, temp, d_cauchyndir_dT, iFTV,
        d_bdnddir_dFV, d_ibdnddir_dFV, d_I1_dF, d_I2_dF, d_I3_dF, d2_cauchyndir_dF_dT);

  if (d_cauchyndir_dF != nullptr)
  {
    // next 3 updates add partial derivative of (\sigma * n * v) w.r.t. F for constant invariants
    // 1. part is term arising from d(J^{-1})/dF
    d_cauchyndir_dF->Update(-prefac * (prinv(1) * dPI(1) * nddir + prinv(2) * dPI(2) * nddir +
                                          dPI(0) * bdnddir - prinv(2) * dPI(1) * ibdnddir),
        iFTV, 0.0);  // DsntDF is cleared here
    // 2. part is term arising from d(b * n * v)/dF
    d_cauchyndir_dF->Update(prefac * dPI(0), d_bdnddir_dFV, 1.0);
    // 3. part is term arising from d(b_el^{-1} * n * v)/dF
    d_cauchyndir_dF->Update(-prefac * prinv(2) * dPI(1), d_ibdnddir_dFV, 1.0);
    // add d(sigma * n * v)/dI1 \otimes dI1/dF
    d_cauchyndir_dF->Update(prefac * (prinv(1) * ddPII(5) * nddir + prinv(2) * ddPII(4) * nddir +
                                         ddPII(0) * bdnddir - prinv(2) * ddPII(5) * ibdnddir),
        d_I1_dF, 1.0);
    // add d(sigma * n * v)/dI2 \otimes dI2/dF
    d_cauchyndir_dF->Update(
        prefac * (dPI(1) * nddir + prinv(1) * ddPII(1) * nddir + prinv(2) * ddPII(3) * nddir +
                     ddPII(5) * bdnddir - prinv(2) * ddPII(1) * ibdnddir),
        d_I2_dF, 1.0);
    // add d(sigma * n * v)/dI3 \otimes dI3/dF
    d_cauchyndir_dF->Update(
        prefac * (prinv(1) * ddPII(3) * nddir + dPI(2) * nddir + prinv(2) * ddPII(2) * nddir +
                     ddPII(4) * bdnddir - dPI(1) * ibdnddir - prinv(2) * ddPII(3) * ibdnddir),
        d_I3_dF, 1.0);
  }

  if (d2_cauchyndir_dF_dn != nullptr)
  {
    // next three blocks add d/dn(d(\sigma * n * v)/dF) for constant invariants
    // first part is term arising from d/dn(dJ^{-1}/dF)
    tempvec3x1.Update(prinv(1) * dPI(1) + prinv(2) * dPI(2), dir, 0.0);
    tempvec3x1.Update(dPI(0), bddir, 1.0);
    tempvec3x1.Update(-prinv(2) * dPI(1), ibddir, 1.0);
    d2_cauchyndir_dF_dn->MultiplyNT(-prefac, iFTV, tempvec3x1, 0.0);

    // second part is term arising from d/dn(d(b * n * v)/dF
    const double fac = prefac * dPI(0);
    tempvec1x3.MultiplyTN(1.0, dir, defgrd, 0.0);
    for (int k = 0; k < 3; ++k)
    {
      for (int l = 0; l < 3; ++l)
      {
        for (int z = 0; z < 3; ++z)
          (*d2_cauchyndir_dF_dn)(CORE::LINALG::VOIGT::IndexMappings::NonSymToVoigt9(k, l), z) +=
              fac * (dir(k, 0) * defgrd(z, l) + static_cast<double>(k == z) * tempvec1x3(0, l));
      }
    }

    // third part is term arising from d/dn(d(b^{-1} * n * t)/dF
    const double fac2 = prefac * prinv(2) * dPI(1);
    for (int k = 0; k < 3; ++k)
    {
      for (int l = 0; l < 3; ++l)
      {
        for (int z = 0; z < 3; ++z)
          (*d2_cauchyndir_dF_dn)(CORE::LINALG::VOIGT::IndexMappings::NonSymToVoigt9(k, l), z) +=
              fac2 * (ibddir(k, 0) * ibdF(z, l) + ib(z, k) * dirdibdF(0, l));
      }
    }

    // add parts originating from d/dn(d(sigma * n * t)/dI1 \otimes dI1/dF)
    tempvec3x1.Update(prinv(1) * ddPII(5) + prinv(2) * ddPII(4), dir, 0.0);
    tempvec3x1.Update(ddPII(0), bddir, 1.0);
    tempvec3x1.Update(-prinv(2) * ddPII(5), ibddir, 1.0);
    d2_cauchyndir_dF_dn->MultiplyNT(prefac, d_I1_dF, tempvec3x1, 1.0);

    // add parts originating from d/dn(d(sigma * n * t)/dI2 \otimes dI2/dF)
    tempvec3x1.Update(dPI(1) + prinv(1) * ddPII(1) + prinv(2) * ddPII(3), dir, 0.0);
    tempvec3x1.Update(ddPII(5), bddir, 1.0);
    tempvec3x1.Update(-prinv(2) * ddPII(1), ibddir, 1.0);
    d2_cauchyndir_dF_dn->MultiplyNT(prefac, d_I2_dF, tempvec3x1, 1.0);

    // add parts originating from d/dn(d(sigma * n * t)/dI3 \otimes dI3/dF)
    tempvec3x1.Update(prinv(1) * ddPII(3) + dPI(2) + prinv(2) * ddPII(2), dir, 0.0);
    tempvec3x1.Update(ddPII(4), bddir, 1.0);
    tempvec3x1.Update(-dPI(1) - prinv(2) * ddPII(3), ibddir, 1.0);
    d2_cauchyndir_dF_dn->MultiplyNT(prefac, d_I3_dF, tempvec3x1, 1.0);
  }

  if (d2_cauchyndir_dF_ddir != nullptr)
  {
    // next three blocks add d/dt(d(\sigma * n * v)/dF) for constant invariants
    // first part is term arising from d/dt(dJ^{-1}/dF)
    tempvec3x1.Update(prinv(1) * dPI(1) + prinv(2) * dPI(2), n, 0.0);
    tempvec3x1.Update(dPI(0), bdn, 1.0);
    tempvec3x1.Update(-prinv(2) * dPI(1), ibdn, 1.0);
    d2_cauchyndir_dF_ddir->MultiplyNT(-prefac, iFTV, tempvec3x1, 0.0);

    // second part is term arising from d/dt(d(b * n * v)/dF
    const double fac = prefac * dPI(0);
    tempvec1x3.MultiplyTN(1.0, n, defgrd, 0.0);
    for (int k = 0; k < 3; ++k)
    {
      for (int l = 0; l < 3; ++l)
      {
        for (int z = 0; z < 3; ++z)
          (*d2_cauchyndir_dF_ddir)(CORE::LINALG::VOIGT::IndexMappings::NonSymToVoigt9(k, l), z) +=
              fac * (n(k, 0) * defgrd(z, l) + static_cast<double>(k == z) * tempvec1x3(0, l));
      }
    }

    // third part is term arising from d/dn(d(b^{-1} * n * v)/dF
    const double fac2 = prefac * prinv(2) * dPI(1);
    for (int k = 0; k < 3; ++k)
    {
      for (int l = 0; l < 3; ++l)
      {
        for (int z = 0; z < 3; ++z)
          (*d2_cauchyndir_dF_ddir)(CORE::LINALG::VOIGT::IndexMappings::NonSymToVoigt9(k, l), z) +=
              fac2 * (ibdn(k, 0) * ibdF(z, l) + ib(z, k) * ndibdF(0, l));
      }
    }

    // add parts originating from d/dt(d(sigma * n * v)/dI1 \otimes dI1/dF)
    tempvec3x1.Update(prinv(1) * ddPII(5) + prinv(2) * ddPII(4), n, 0.0);
    tempvec3x1.Update(ddPII(0), bdn, 1.0);
    tempvec3x1.Update(-prinv(2) * ddPII(5), ibdn, 1.0);
    d2_cauchyndir_dF_ddir->MultiplyNT(prefac, d_I1_dF, tempvec3x1, 1.0);

    // add parts originating from d/dt(d(sigma * n * v)/dI2 \otimes dI2/dF)
    tempvec3x1.Update(dPI(1) + prinv(1) * ddPII(1) + prinv(2) * ddPII(3), n, 0.0);
    tempvec3x1.Update(ddPII(5), bdn, 1.0);
    tempvec3x1.Update(-prinv(2) * ddPII(1), ibdn, 1.0);
    d2_cauchyndir_dF_ddir->MultiplyNT(prefac, d_I2_dF, tempvec3x1, 1.0);

    // add parts originating from d/dt(d(sigma * n * v)/dI3 \otimes dI3/dF)
    tempvec3x1.Update(prinv(1) * ddPII(3) + dPI(2) + prinv(2) * ddPII(2), n, 0.0);
    tempvec3x1.Update(ddPII(4), bdn, 1.0);
    tempvec3x1.Update(-dPI(1) - prinv(2) * ddPII(3), ibdn, 1.0);
    d2_cauchyndir_dF_ddir->MultiplyNT(prefac, d_I3_dF, tempvec3x1, 1.0);
  }

  if (d2_cauchyndir_dF2 != nullptr)
  {
    // define and fill all tensors that can not be calculated using multiply operations first
    static CORE::LINALG::Matrix<9, 9> d_iFT_dF(true);
    static CORE::LINALG::Matrix<9, 9> d2_bdnddir_dF2(true);
    static CORE::LINALG::Matrix<9, 9> d2_ibdnddir_dF2(true);
    static CORE::LINALG::Matrix<9, 9> d2_I1_dF2(true);
    static CORE::LINALG::Matrix<9, 9> d2_I2_dF2(true);
    static CORE::LINALG::Matrix<9, 9> d2_I3_dF2(true);
    d_iFT_dF.Clear();
    d2_bdnddir_dF2.Clear();
    d2_ibdnddir_dF2.Clear();
    d2_I1_dF2.Clear();
    d2_I2_dF2.Clear();
    d2_I3_dF2.Clear();

    static CORE::LINALG::Matrix<3, 3> C(true);
    C.MultiplyTN(1.0, defgrd, defgrd, 0.0);

    using map = CORE::LINALG::VOIGT::IndexMappings;

    for (int k = 0; k < 3; ++k)
    {
      for (int l = 0; l < 3; ++l)
      {
        for (int m = 0; m < 3; ++m)
        {
          for (int a = 0; a < 3; ++a)
          {
            d_iFT_dF(map::NonSymToVoigt9(k, l), map::NonSymToVoigt9(m, a)) = -iF(l, m) * iF(a, k);
            d2_bdnddir_dF2(map::NonSymToVoigt9(k, l), map::NonSymToVoigt9(m, a)) =
                (dir(k, 0) * n(m, 0) + dir(m, 0) * n(k, 0)) * static_cast<double>(l == a);
            d2_ibdnddir_dF2(map::NonSymToVoigt9(k, l), map::NonSymToVoigt9(m, a)) =
                ibdF(k, a) * (ibddir(m, 0) * ndibdF(0, l) + ibdn(m, 0) * dirdibdF(0, l)) +
                ib(m, k) * (dirdibdF(0, a) * ndibdF(0, l) + dirdibdF(0, l) * ndibdF(0, a)) +
                ibdF(m, l) * (ibddir(k, 0) * ndibdF(0, a) + dirdibdF(0, a) * ibdn(k, 0));
            d2_I1_dF2(map::NonSymToVoigt9(k, l), map::NonSymToVoigt9(m, a)) =
                2.0 * static_cast<double>(k == m) * static_cast<double>(l == a);
            d2_I2_dF2(map::NonSymToVoigt9(k, l), map::NonSymToVoigt9(m, a)) =
                2.0 *
                (prinv(0) * static_cast<double>(k == m) * static_cast<double>(l == a) +
                    2.0 * defgrd(m, a) * defgrd(k, l) - static_cast<double>(k == m) * C(a, l) -
                    defgrd(k, a) * defgrd(m, l) - b(k, m) * static_cast<double>(l == a));
            d2_I3_dF2(map::NonSymToVoigt9(k, l), map::NonSymToVoigt9(m, a)) =
                2.0 * prinv(2) * (2.0 * ibdF(m, a) * ibdF(k, l) - ibdF(m, l) * ibdF(k, a));
          }
        }
      }
    }

    // terms below add contributions originating from d(1st term of DsntDF)/dF
    d2_cauchyndir_dF2->MultiplyNT(prefac * (prinv(1) * dPI(1) * nddir + prinv(2) * dPI(2) * nddir +
                                               dPI(0) * bdnddir - prinv(2) * dPI(1) * ibdnddir),
        iFTV, iFTV, 0.0);  // D2sntDF2 is cleared here
    d2_cauchyndir_dF2->Update(-prefac * (prinv(1) * dPI(1) * nddir + prinv(2) * dPI(2) * nddir +
                                            dPI(0) * bdnddir - prinv(2) * dPI(1) * ibdnddir),
        d_iFT_dF, 1.0);
    d2_cauchyndir_dF2->MultiplyNT(-prefac * dPI(0), iFTV, d_bdnddir_dFV, 1.0);
    d2_cauchyndir_dF2->MultiplyNT(prefac * prinv(2) * dPI(1), iFTV, d_ibdnddir_dFV, 1.0);

    d2_cauchyndir_dF2->MultiplyNT(
        -prefac * (prinv(1) * ddPII(5) * nddir + prinv(2) * ddPII(4) * nddir + ddPII(0) * bdnddir -
                      prinv(2) * ddPII(5) * ibdnddir),
        iFTV, d_I1_dF, 1.0);
    d2_cauchyndir_dF2->MultiplyNT(
        -prefac * (dPI(1) * nddir + prinv(1) * ddPII(1) * nddir + prinv(2) * ddPII(3) * nddir +
                      ddPII(5) * bdnddir - prinv(2) * ddPII(1) * ibdnddir),
        iFTV, d_I2_dF, 1.0);
    d2_cauchyndir_dF2->MultiplyNT(
        -prefac * (prinv(1) * ddPII(3) * nddir + dPI(2) * nddir + prinv(2) * ddPII(2) * nddir +
                      ddPII(4) * bdnddir - dPI(1) * ibdnddir - prinv(2) * ddPII(3) * ibdnddir),
        iFTV, d_I3_dF, 1.0);

    // terms below add contributions originating from d(2nd term of DsntDF)/dF
    d2_cauchyndir_dF2->MultiplyNT(-prefac * dPI(0), d_bdnddir_dFV, iFTV, 1.0);
    d2_cauchyndir_dF2->Update(prefac * dPI(0), d2_bdnddir_dF2, 1.0);
    d2_cauchyndir_dF2->MultiplyNT(prefac * ddPII(0), d_bdnddir_dFV, d_I1_dF, 1.0);
    d2_cauchyndir_dF2->MultiplyNT(prefac * ddPII(5), d_bdnddir_dFV, d_I2_dF, 1.0);
    d2_cauchyndir_dF2->MultiplyNT(prefac * ddPII(4), d_bdnddir_dFV, d_I3_dF, 1.0);

    // terms below add contributions originating from d(3rd term of DsntDF)/dF
    d2_cauchyndir_dF2->MultiplyNT(prefac * prinv(2) * dPI(1), d_ibdnddir_dFV, iFTV, 1.0);
    d2_cauchyndir_dF2->Update(-prefac * prinv(2) * dPI(1), d2_ibdnddir_dF2, 1.0);
    d2_cauchyndir_dF2->MultiplyNT(-prefac * prinv(2) * ddPII(5), d_ibdnddir_dFV, d_I1_dF, 1.0);
    d2_cauchyndir_dF2->MultiplyNT(-prefac * prinv(2) * ddPII(1), d_ibdnddir_dFV, d_I2_dF, 1.0);
    d2_cauchyndir_dF2->MultiplyNT(
        -prefac * (dPI(1) + prinv(2) * ddPII(3)), d_ibdnddir_dFV, d_I3_dF, 1.0);

    // terms below add contributions originating from d(4th term of DsntDF)/dF
    d2_cauchyndir_dF2->MultiplyNT(
        -prefac * (prinv(1) * ddPII(5) * nddir + prinv(2) * ddPII(4) * nddir + ddPII(0) * bdnddir -
                      prinv(2) * ddPII(5) * ibdnddir),
        d_I1_dF, iFTV, 1.0);
    d2_cauchyndir_dF2->MultiplyNT(prefac * ddPII(0), d_I1_dF, d_bdnddir_dFV, 1.0);
    d2_cauchyndir_dF2->MultiplyNT(-prefac * prinv(2) * ddPII(5), d_I1_dF, d_ibdnddir_dFV, 1.0);
    d2_cauchyndir_dF2->Update(prefac * (prinv(1) * ddPII(5) * nddir + prinv(2) * ddPII(4) * nddir +
                                           ddPII(0) * bdnddir - prinv(2) * ddPII(5) * ibdnddir),
        d2_I1_dF2, 1.0);
    d2_cauchyndir_dF2->MultiplyNT(
        prefac * (prinv(1) * dddPIII(5) * nddir + prinv(2) * dddPIII(6) * nddir +
                     dddPIII(0) * bdnddir - prinv(2) * dddPIII(5) * ibdnddir),
        d_I1_dF, d_I1_dF, 1.0);
    d2_cauchyndir_dF2->MultiplyNT(
        prefac * (ddPII(5) * nddir + prinv(1) * dddPIII(3) * nddir + prinv(2) * dddPIII(9) * nddir +
                     dddPIII(5) * bdnddir - prinv(2) * dddPIII(3) * ibdnddir),
        d_I1_dF, d_I2_dF, 1.0);
    d2_cauchyndir_dF2->MultiplyNT(
        prefac * (prinv(1) * dddPIII(9) * nddir + ddPII(4) * nddir + prinv(2) * dddPIII(4) * nddir +
                     dddPIII(6) * bdnddir - ddPII(5) * ibdnddir - prinv(2) * dddPIII(9) * ibdnddir),
        d_I1_dF, d_I3_dF, 1.0);

    // terms below add contributions originating from d(5th term of DsntDF)/dF
    d2_cauchyndir_dF2->MultiplyNT(
        -prefac * (dPI(1) * nddir + prinv(1) * ddPII(1) * nddir + prinv(2) * ddPII(3) * nddir +
                      ddPII(5) * bdnddir - prinv(2) * ddPII(1) * ibdnddir),
        d_I2_dF, iFTV, 1.0);
    d2_cauchyndir_dF2->MultiplyNT(prefac * ddPII(5), d_I2_dF, d_bdnddir_dFV, 1.0);
    d2_cauchyndir_dF2->MultiplyNT(-prefac * prinv(2) * ddPII(1), d_I2_dF, d_ibdnddir_dFV, 1.0);
    d2_cauchyndir_dF2->Update(
        prefac * (dPI(1) * nddir + prinv(1) * ddPII(1) * nddir + prinv(2) * ddPII(3) * nddir +
                     ddPII(5) * bdnddir - prinv(2) * ddPII(1) * ibdnddir),
        d2_I2_dF2, 1.0);
    d2_cauchyndir_dF2->MultiplyNT(
        prefac * (ddPII(5) * nddir + prinv(1) * dddPIII(3) * nddir + prinv(2) * dddPIII(9) * nddir +
                     dddPIII(5) * bdnddir - prinv(2) * dddPIII(3) * ibdnddir),
        d_I2_dF, d_I1_dF, 1.0);
    d2_cauchyndir_dF2->MultiplyNT(
        prefac * (2.0 * ddPII(1) * nddir + prinv(1) * dddPIII(1) * nddir +
                     prinv(2) * dddPIII(7) * nddir + dddPIII(3) * bdnddir -
                     prinv(2) * dddPIII(1) * ibdnddir),
        d_I2_dF, d_I2_dF, 1.0);
    d2_cauchyndir_dF2->MultiplyNT(
        prefac * (2.0 * ddPII(3) * nddir + prinv(1) * dddPIII(7) * nddir +
                     prinv(2) * dddPIII(8) * nddir + dddPIII(9) * bdnddir - ddPII(1) * ibdnddir -
                     prinv(2) * dddPIII(7) * ibdnddir),
        d_I2_dF, d_I3_dF, 1.0);

    // terms below add contributions originating from d(6th term of DsntDF)/dF
    d2_cauchyndir_dF2->MultiplyNT(
        -prefac * (prinv(1) * ddPII(3) * nddir + dPI(2) * nddir + prinv(2) * ddPII(2) * nddir +
                      ddPII(4) * bdnddir - dPI(1) * ibdnddir - prinv(2) * ddPII(3) * ibdnddir),
        d_I3_dF, iFTV, 1.0);
    d2_cauchyndir_dF2->MultiplyNT(prefac * ddPII(4), d_I3_dF, d_bdnddir_dFV, 1.0);
    d2_cauchyndir_dF2->MultiplyNT(
        -prefac * (dPI(1) + prinv(2) * ddPII(3)), d_I3_dF, d_ibdnddir_dFV, 1.0);
    d2_cauchyndir_dF2->Update(
        prefac * (prinv(1) * ddPII(3) * nddir + dPI(2) * nddir + prinv(2) * ddPII(2) * nddir +
                     ddPII(4) * bdnddir - dPI(1) * ibdnddir - prinv(2) * ddPII(3) * ibdnddir),
        d2_I3_dF2, 1.0);
    d2_cauchyndir_dF2->MultiplyNT(
        prefac * (prinv(1) * dddPIII(9) * nddir + ddPII(4) * nddir + prinv(2) * dddPIII(4) * nddir +
                     dddPIII(6) * bdnddir - ddPII(5) * ibdnddir - prinv(2) * dddPIII(9) * ibdnddir),
        d_I3_dF, d_I1_dF, 1.0);
    d2_cauchyndir_dF2->MultiplyNT(
        prefac * (2.0 * ddPII(3) * nddir + prinv(1) * dddPIII(7) * nddir +
                     prinv(2) * dddPIII(8) * nddir + dddPIII(9) * bdnddir - ddPII(1) * ibdnddir -
                     prinv(2) * dddPIII(7) * ibdnddir),
        d_I3_dF, d_I2_dF, 1.0);
    d2_cauchyndir_dF2->MultiplyNT(
        prefac * (prinv(1) * dddPIII(8) * nddir + 2.0 * ddPII(2) * nddir +
                     prinv(2) * dddPIII(2) * nddir + dddPIII(4) * bdnddir -
                     2.0 * ddPII(3) * ibdnddir - prinv(2) * dddPIII(8) * ibdnddir),
        d_I3_dF, d_I3_dF, 1.0);
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ElastHyper::VisNames(std::map<std::string, int>& names)
{
  if (AnisotropicPrincipal() or AnisotropicModified())
  {
    std::vector<CORE::LINALG::Matrix<3, 1>> fibervecs;
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
  for (auto& p : potsum_)
  {
    p->VisNames(names);
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
    std::vector<CORE::LINALG::Matrix<3, 1>> fibervecs;
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
  for (auto& p : potsum_)
  {
    return_val += static_cast<int>(p->VisData(name, data, numgp, eleID));
  }
  return (bool)return_val;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const MAT::ELASTIC::Summand> MAT::ElastHyper::GetPotSummandPtr(
    const INPAR::MAT::MaterialType& materialtype) const
{
  for (const auto& p : potsum_)
  {
    if (p->MaterialType() == materialtype) return p;
  }
  return Teuchos::null;
}
