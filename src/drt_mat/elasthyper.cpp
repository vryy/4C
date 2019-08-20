/*----------------------------------------------------------------------*/
/*! \file
\brief
This file contains the hyperelastic toolbox. It allows summing up several summands
of several types (isotropic or anisotropic, splitted or not) to build a hyperelastic
strain energy function.

The input line should read
MAT 0   MAT_ElastHyper   NUMMAT 2 MATIDS 1 2 DENS 0

\level 1

\maintainer Amadeus Gebauer

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
#include "elasthyper_service.H"

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



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::ElastHyper::ElastHyper() : summandProperties_(), params_(NULL), potsum_(0) {}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::ElastHyper::ElastHyper(MAT::PAR::ElastHyper* params)
    : summandProperties_(), params_(params), potsum_(0)
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
  summandProperties_.Pack(data);

  if (params_ != NULL)  // summands are not accessible in postprocessing mode
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
  params_ = NULL;
  potsum_.clear();

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

  summandProperties_.Unpack(position, data);

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
  summandProperties_.Clear();
  ElastHyperProperties(potsum_, summandProperties_);

  if (summandProperties_.viscoGeneral)
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
  if (summandProperties_.anisoprinc || summandProperties_.anisomod)
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
  if (summandProperties_.anisoprinc || summandProperties_.anisomod)
  {
    for (unsigned int p = 0; p < potsum_.size(); ++p)
    {
      potsum_[p]->SetFiberVecs(newgamma, locsys, defgrd);
    }
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ElastHyper::StrainEnergy(
    const LINALG::Matrix<6, 1>& glstrain, double& psi, const int eleGID)
{
  static LINALG::Matrix<6, 1> C_strain(true);
  C_strain.Clear();
  static LINALG::Matrix<3, 1> prinv(true);
  prinv.Clear();
  static LINALG::Matrix<3, 1> modinv(true);
  modinv.Clear();

  EvaluateRightCauchyGreenStrainLikeVoigt(glstrain, C_strain);
  InvariantsPrincipal<VoigtNotation::strain>(prinv, C_strain);
  InvariantsModified(modinv, prinv);

  // loop map of associated potential summands
  for (const auto& p : potsum_)
  {
    p->AddStrainEnergy(psi, prinv, modinv, glstrain, eleGID);
  }
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
  bool checkpolyconvexity = (params_ != nullptr and params_->polyconvex_);

  ElastHyperEvaluate(*defgrd, *glstrain, params, *stress, *cmat, eleGID, potsum_,
      summandProperties_, checkpolyconvexity);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ElastHyper::EvaluateCauchyDerivs(const LINALG::Matrix<3, 1>& prinv, const int eleGID,
    LINALG::Matrix<3, 1>& dPI, LINALG::Matrix<6, 1>& ddPII, LINALG::Matrix<10, 1>& dddPIII,
    const double* temp)
{
  for (unsigned i = 0; i < potsum_.size(); ++i)
  {
    if (summandProperties_.isoprinc)
    {
      potsum_[i]->AddDerivativesPrincipal(dPI, ddPII, prinv, eleGID);
      potsum_[i]->AddThirdDerivativesPrincipalIso(dddPIII, prinv, eleGID);
    }
    if (summandProperties_.isomod || summandProperties_.anisomod || summandProperties_.anisoprinc)
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
  InvariantsPrincipal<MAT::VoigtNotation::strain>(prinv, bV_strain);

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
