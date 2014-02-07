/*----------------------------------------------------------------------*/
/*!
\file plasticelasthyper.cpp
\brief
This file contains the hyperelastic toolbox with application to finite
strain plasticity using a semi-smooth Newton method. It allows summing up
several summands of isotropic non-splitted type to build
a hyperelastic strain energy function.

The input line should read
MAT 0   MAT_ElastHyper   NUMMAT 0 MATIDS  DENS 0 INITYIELD 0.0 ISOHARD 0.0 EXPISOHARD 0.0 INFYIELD 0.0 KINHARD 0.0 CPL 0.0 STAB_S 0.0

<pre>
Maintainer: Alexander Seitz
            seitz@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15271
</pre>
*/

/*----------------------------------------------------------------------*/

#include "plasticelasthyper.H"
#include "../drt_lib/standardtypes_cpp.H"
#include "../drt_matelast/elast_summand.H"
#include "../linalg/linalg_utils.H"
#include "../drt_lib/drt_linedefinition.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_mat/material_service.H"
#include "Epetra_SerialDenseSolver.h"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::PAR::PlasticElastHyper::PlasticElastHyper(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: Parameter(matdata),
  nummat_(matdata->GetInt("NUMMAT")),
  matids_(matdata->Get<std::vector<int> >("MATIDS")),
  density_(matdata->GetDouble("DENS")),
  inityield_(matdata->GetDouble("INITYIELD")),
  isohard_(matdata->GetDouble("ISOHARD")),
  expisohard_(matdata->GetDouble("EXPISOHARD")),
  infyield_(matdata->GetDouble("INFYIELD")),
  kinhard_(matdata->GetDouble("KINHARD")),
  rY_11_(matdata->GetDouble("rY_11")),
  rY_22_(matdata->GetDouble("rY_22")),
  rY_33_(matdata->GetDouble("rY_33")),
  rY_12_(matdata->GetDouble("rY_12")),
  rY_23_(matdata->GetDouble("rY_23")),
  rY_13_(matdata->GetDouble("rY_13")),
  plspin_eta_(matdata->GetDouble("PL_SPIN_ETA"))
{
  // check if sizes fit
  if (nummat_ != (int)matids_->size())
    dserror("number of materials %d does not fit to size of material vector %d", nummat_, matids_->size());

  // check plastic parameter validity
  if (inityield_<=0.)
    dserror("initial yield stress must be positive");
  if (infyield_<inityield_)
    dserror("saturation yield stress must not be less than initial yield stress");
  if (expisohard_<0.)
    dserror("Nonlinear hardening exponent must be non-negative");

  if (rY_11_!=0. || rY_22_!=0. || rY_33_!=0. || rY_12_!=0. || rY_23_!=0. || rY_13_!=0.)
    if (rY_11_<=0. || rY_22_<=0. || rY_33_<=0. || rY_12_<=0. || rY_23_<=0. || rY_13_<=0.)
      dserror("Hill parameters all must be positive (incomplete set?)");


}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<MAT::Material> MAT::PAR::PlasticElastHyper::CreateMaterial()
{
  return Teuchos::rcp(new MAT::PlasticElastHyper(this));
}


MAT::PlasticElastHyperType MAT::PlasticElastHyperType::instance_;


DRT::ParObject* MAT::PlasticElastHyperType::Create( const std::vector<char> & data )
{
  MAT::PlasticElastHyper* elhy = new MAT::PlasticElastHyper();
  elhy->Unpack(data);

  return elhy;
}


/*----------------------------------------------------------------------*
 |  initialise static arrays                                 bborn 08/09|
 *----------------------------------------------------------------------*/
// 6-Voigt C-index                              0 1 2  3 4 5
const int MAT::PlasticElastHyper::VOIGT6ROW_[6] = {0,1,2, 0,1,2};
const int MAT::PlasticElastHyper::VOIGT6COL_[6] = {0,1,2, 1,2,0};

// tensor indices ij = 11, 12, 13, 21, 22, 23, 31, 32, 33
// C indices           00, 01, 02, 10, 11, 12, 20, 21, 22
// Access : 3*i+j
// 6-Voigt C-indices    0   3   5   3   1   4   5   4   2
const int MAT::PlasticElastHyper::VOIGT3X3SYM_[9] = {0,3,5, 3,1,4, 5,4,2};

const int MAT::PlasticElastHyper::VOIGT3X3_[3][3] = {{0,3,5},{6,1,4},{8,7,2}};


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::PlasticElastHyper::PlasticElastHyper()
  : params_(NULL),
    potsum_(0),
    PlAniso_(Teuchos::null),
    InvPlAniso_(Teuchos::null)

{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::PlasticElastHyper::PlasticElastHyper(MAT::PAR::PlasticElastHyper* params)
  : params_(params),
    potsum_(0),
    PlAniso_(Teuchos::null),
    InvPlAniso_(Teuchos::null)

{
  // make sure the referenced materials in material list have quick access parameters
  std::vector<int>::const_iterator m;
  for (m=params_->matids_->begin(); m!=params_->matids_->end(); ++m)
  {
    const int matid = *m;
    Teuchos::RCP<MAT::ELASTIC::Summand> sum = MAT::ELASTIC::Summand::Factory(matid);
    if (sum == Teuchos::null) dserror("Failed to allocate");
    potsum_.push_back(sum);
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::PlasticElastHyper::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm( data );
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);
  // matid
  int matid = -1;
  if (params_ != NULL) matid = params_->Id();  // in case we are in post-process mode
  AddtoPack(data,matid);
  AddtoPack(data,isoprinc_);
  AddtoPack(data,isomod_);
  AddtoPack(data,anisoprinc_);
  AddtoPack(data,anisomod_);
  AddtoPack(data,isomodvisco_);

  // hill plasticity
  bool Hill=(bool)(PlAniso_!=Teuchos::null);
  AddtoPack(data,(int)Hill);
  if (Hill)
  {
    AddtoPack(data,*PlAniso_);
    AddtoPack(data,*InvPlAniso_);
  }

  if (params_ != NULL) // summands are not accessible in postprocessing mode
  {
    // loop map of associated potential summands
    for (unsigned int p=0; p<potsum_.size(); ++p)
    {
     potsum_[p]->PackSummand(data);
    }
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::PlasticElastHyper::Unpack(const std::vector<char>& data)
{
  // make sure we have a pristine material
  params_ = NULL;
  potsum_.clear();

  isoprinc_ = false;
  isomod_ = false;
  anisoprinc_ = false;
  anisomod_ = false;
  isomodvisco_ = false;

  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");

  // matid and recover params_
  int matid;
  ExtractfromPack(position,data,matid);
  if (DRT::Problem::Instance()->Materials() != Teuchos::null)
  {
    if (DRT::Problem::Instance()->Materials()->Num() != 0)
    {
      const unsigned int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();
      MAT::PAR::Parameter* mat = DRT::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        params_ = static_cast<MAT::PAR::PlasticElastHyper*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(), MaterialType());
    }
  }

  int isoprinc;
  int isomod;
  int anisoprinc;
  int anisomod;
  int isomodvisco;

  ExtractfromPack(position,data,isoprinc);
  ExtractfromPack(position,data,isomod);
  ExtractfromPack(position,data,anisoprinc);
  ExtractfromPack(position,data,anisomod);
  ExtractfromPack(position,data,isomodvisco);

  if (isoprinc != 0) isoprinc_ = true;
  if (isomod != 0) isomod_ = true;
  if (anisoprinc != 0) anisoprinc_ = true;
  if (anisomod != 0) anisomod_ = true;
  if (isomodvisco != 0) isomodvisco_ = true;
  
  // hill plasticity information
  bool Hill=(bool)ExtractInt(position,data);
  if (Hill)
  {
    LINALG::Matrix<5,5> tmp55;
    ExtractfromPack(position,data,tmp55);
    PlAniso_=Teuchos::rcp( new LINALG::Matrix<5,5>(tmp55));
    ExtractfromPack(position,data,tmp55);
    InvPlAniso_=Teuchos::rcp( new LINALG::Matrix<5,5>(tmp55));
  }
  else
  {
    PlAniso_=Teuchos::null;
    InvPlAniso_=Teuchos::null;
  }

  if (params_ != NULL) // summands are not accessible in postprocessing mode
  {
    // make sure the referenced materials in material list have quick access parameters
    std::vector<int>::const_iterator m;
    for (m=params_->matids_->begin(); m!=params_->matids_->end(); ++m)
    {
      const int matid = *m;
      Teuchos::RCP<MAT::ELASTIC::Summand> sum = MAT::ELASTIC::Summand::Factory(matid);
      if (sum == Teuchos::null) dserror("Failed to allocate");
      potsum_.push_back(sum);
    }

    // loop map of associated potential summands
    for (unsigned int p=0; p<potsum_.size(); ++p)
    {
     potsum_[p]->UnpackSummand(data,position);
    }
    // in the postprocessing mode, we do not unpack everything we have packed
    // -> position check cannot be done in this case
    if (position != data.size())
      dserror("Mismatch in size of data %d <-> %d",data.size(),position);
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
int MAT::PlasticElastHyper::MatID(
  const unsigned index
) const
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
double MAT::PlasticElastHyper::ShearMod() const
{
  // principal coefficients
  bool haveshearmod = false;
  double shearmod = 0.0;
  {
    // loop map of associated potential summands
    for (unsigned int p=0; p<potsum_.size(); ++p)
    {
     potsum_[p]->AddShearMod(haveshearmod,shearmod);
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
void MAT::PlasticElastHyper::Setup(int numgp, DRT::INPUT::LineDefinition* linedef)
{
  // Setup summands
  for (unsigned int p=0; p<potsum_.size(); ++p)
  {
    potsum_[p]->Setup(linedef);
  }

  // find out which formulations are used

  isoprinc_ = false ;
  isomod_ = false ;
  anisoprinc_ = false ;
  anisomod_ = false;
  isomodvisco_ = false;

  for (unsigned int p=0; p<potsum_.size(); ++p)
  {
    potsum_[p]->SpecifyFormulation(isoprinc_,isomod_,anisoprinc_,anisomod_,isomodvisco_);
  }
  // in this case the mandel stress become non-symmetric and the
  // calculated derivatives have to be extended.
  if (anisomod_==true || anisoprinc_==true)
    dserror("PlasticElastHyper only for isotropic elastic material!");
  if (isomod_==true)
    dserror("PlasticElastHyper only for isotropic elastic material without volumetric split!"
        "\n(Extension without major changes possilble.)");

    // check if either zero or three fiber directions are given
    if (linedef->HaveNamed("FIBER1") || linedef->HaveNamed("FIBER2") || linedef->HaveNamed("FIBER3"))
      if (!linedef->HaveNamed("FIBER1") || !linedef->HaveNamed("FIBER2") || !linedef->HaveNamed("FIBER3"))
        dserror("so3 expects no fibers or 3 fiber directions");

    if (HaveHillPlasticity())
      SetupHillPlasticity(linedef);

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::PlasticElastHyper::SetupHillPlasticity(DRT::INPUT::LineDefinition* linedef)
{
  PlAniso_=Teuchos::rcp( new LINALG::Matrix<5,5>);
  InvPlAniso_=Teuchos::rcp( new LINALG::Matrix<5,5>);

  // anisotropy directions
  std::vector<LINALG::Matrix<3,1> > directions(3);

  // first anisotropy direction
  if (linedef->HaveNamed("FIBER1"))
  {
    std::vector<double> fiber;
     linedef->ExtractDoubleVector("FIBER1",fiber);
     double fnorm=0.;
     //normalization
     for (int i = 0; i < 3; ++i)
       fnorm += fiber[i]*fiber[i];
     fnorm = sqrt(fnorm);
     if (fnorm==0.)
       dserror("Fiber vector has norm zero");

     // fill final fiber vector
     for (int i = 0; i < 3; ++i)
       directions.at(0)(i) = fiber[i]/fnorm;
  }
  // second anisotropy direction
  if (linedef->HaveNamed("FIBER2"))
  {
    std::vector<double> fiber;
     linedef->ExtractDoubleVector("FIBER2",fiber);
     double fnorm=0.;
     //normalization
     for (int i = 0; i < 3; ++i)
       fnorm += fiber[i]*fiber[i];
     fnorm = sqrt(fnorm);
     if (fnorm==0.)
       dserror("Fiber vector has norm zero");

     // fill final fiber vector
     for (int i = 0; i < 3; ++i)
       directions.at(1)(i) = fiber[i]/fnorm;
  }
  // third anisotropy direction
  if (linedef->HaveNamed("FIBER3"))
  {
    std::vector<double> fiber;
     linedef->ExtractDoubleVector("FIBER3",fiber);
     double fnorm=0.;
     //normalization
     for (int i = 0; i < 3; ++i)
       fnorm += fiber[i]*fiber[i];
     fnorm = sqrt(fnorm);
     if (fnorm==0.)
       dserror("Fiber vector has norm zero");

     // fill final fiber vector
     for (int i = 0; i < 3; ++i)
       directions.at(2)(i) = fiber[i]/fnorm;
  }

//    std::cout << "fiber1: " << directions.at(0) << "fiber2: " << directions.at(1) << "fiber3: " << directions.at(2) << std::endl;

  // check orthogonality
  LINALG::Matrix<1,1> matrix1;
  matrix1.MultiplyTN(directions.at(0),directions.at(1));
  if (std::abs(matrix1(0,0))>1.e-16)
    dserror("fiber directions not orthogonal");
  matrix1.MultiplyTN(directions.at(0),directions.at(2));
  if (std::abs(matrix1(0,0))>1.e-16)
    dserror("fiber directions not orthogonal");
  matrix1.MultiplyTN(directions.at(2),directions.at(1));
  if (std::abs(matrix1(0,0))>1.e-16)
    dserror("fiber directions not orthogonal");

  // check right-handed trihedron
  LINALG::Matrix<3,1> A0xA1;
  A0xA1(0) = (directions.at(0)(1)*directions.at(1)(2)-directions.at(0)(2)*directions.at(1)(1));
  A0xA1(1) = (directions.at(0)(2)*directions.at(1)(0)-directions.at(0)(0)*directions.at(1)(2));
  A0xA1(2) = (directions.at(0)(0)*directions.at(1)(1)-directions.at(0)(1)*directions.at(1)(0));
  A0xA1.Update(-1.,directions.at(2),1.);
  if (A0xA1.Norm2()>1.e-8)

    dserror("fibers don't form right-handed trihedron");

  // setup anisotropy tensor
  LINALG::Matrix<6,6> PlAnisoFull(true);

  // setup structural tensor for first and second direction
  // (as the directions are orthogonal, 2 structural tensors are sufficient)
  LINALG::Matrix<3,3> M0;
  M0.MultiplyNT(directions.at(0),directions.at(0));
  LINALG::Matrix<3,3> M1;
  M1.MultiplyNT(directions.at(1),directions.at(1));
  LINALG::Matrix<3,3> M2;
  M2.MultiplyNT(directions.at(2),directions.at(2));

  double alpha1 = 2./3./params_->rY_11_/params_->rY_11_;
  double alpha2 = 2./3./params_->rY_22_/params_->rY_22_;
  double alpha3 = 2./3./params_->rY_33_/params_->rY_33_;
  double alpha4 = 1./3./params_->rY_12_/params_->rY_12_;
  double alpha5 = 1./3./params_->rY_23_/params_->rY_23_;
  double alpha6 = 1./3./params_->rY_13_/params_->rY_13_;


  ElastSymTensorMultiply(PlAnisoFull,alpha1,M0,M0,1.);
  ElastSymTensorMultiply(PlAnisoFull,alpha2,M1,M1,1.);
  ElastSymTensorMultiply(PlAnisoFull,alpha3,M2,M2,1.);
  ElastSymTensorMultiplyAddSym(PlAnisoFull,0.5*(alpha3-alpha1-alpha2),M0,M1,1.);
  ElastSymTensorMultiplyAddSym(PlAnisoFull,0.5*(alpha1-alpha2-alpha3),M1,M2,1.);
  ElastSymTensorMultiplyAddSym(PlAnisoFull,0.5*(alpha2-alpha3-alpha1),M0,M2,1.);
  AddtodMdC_gamma2(PlAnisoFull,M0,M1,alpha4*2.);
  AddtodMdC_gamma2(PlAnisoFull,M1,M2,alpha5*2.);
  AddtodMdC_gamma2(PlAnisoFull,M0,M2,alpha6*2.);

  LINALG::Matrix<6,5> red(true);
  red(0,0) = 1.;
  red(1,1) = 1.;
  red(2,0) = -1.;
  red(2,1) = -1.;
  red(3,2) = 1.;
  red(4,3) = 1.;
  red(5,4) = 1.;
  LINALG::Matrix<6,5> tmp;
  tmp.Multiply(PlAnisoFull,red);
  PlAniso_->MultiplyTN(red,tmp);
  Epetra_SerialDenseMatrix PlAniso_epetra(5,5);
  for (int i=0; i<5; i++)
    for (int j=0; j<5; j++)
      PlAniso_epetra(i,j) = (*PlAniso_)(i,j);
    Epetra_SerialDenseSolver solver;
    solver.SetMatrix(PlAniso_epetra);
    solver.Invert();
    for (int i=0; i<5; i++)
      for (int j=0; j<5; j++)
        (*InvPlAniso_)(i,j) = PlAniso_epetra(i,j);

    LINALG::Matrix<5,6> tmp56;
    tmp56.MultiplyNT(*InvPlAniso_,red);
    InvPlAniso_->Multiply(tmp56,red);
    for (int i=2; i<5; i++) (*InvPlAniso_)(i,i)*=2.;

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::PlasticElastHyper::Update()
{
  // loop map of associated potential summands
  for (unsigned int p=0; p<potsum_.size(); ++p)
  {
    potsum_[p]->Update();
  }

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::PlasticElastHyper::GetFiberVecs(std::vector<LINALG::Matrix<3,1> >& fibervecs)
{
  if (anisoprinc_ || anisomod_)
  {
    for (unsigned int p=0; p<potsum_.size(); ++p)
    {
      potsum_[p]->GetFiberVecs(fibervecs);
    }
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::PlasticElastHyper::EvaluateFiberVecs(
    const double newgamma,
    const LINALG::Matrix<3,3>& locsys,
    const LINALG::Matrix<3,3>& defgrd)
{
  if (anisoprinc_ || anisomod_)
  {
    for (unsigned int p=0; p<potsum_.size(); ++p)
    {
      potsum_[p]->SetFiberVecs(newgamma,locsys,defgrd);
    }
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::PlasticElastHyper::InvariantsPrincipal(
    LINALG::Matrix<3,1>& prinv,
    const LINALG::Matrix<6,1>& rcg)
{
  // 1st invariant, trace
  prinv(0) = rcg(0) + rcg(1) + rcg(2);
  // 2nd invariant
  prinv(1) = 0.5*( prinv(0)*prinv(0)
                   - rcg(0)*rcg(0) - rcg(1)*rcg(1) - rcg(2)*rcg(2)
                   - .5*rcg(3)*rcg(3) - .5*rcg(4)*rcg(4) - .5*rcg(5)*rcg(5) );
  // 3rd invariant, determinant
  prinv(2) = rcg(0)*rcg(1)*rcg(2)
    + 0.25 * rcg(3)*rcg(4)*rcg(5)
    - 0.25 * rcg(1)*rcg(5)*rcg(5)
    - 0.25 * rcg(2)*rcg(3)*rcg(3)
    - 0.25 * rcg(0)*rcg(4)*rcg(4);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::PlasticElastHyper::InvariantsModified(
    LINALG::Matrix<3,1>& modinv,  ///< modified invariants
    const LINALG::Matrix<3,1>& prinv  ///< principal invariants
    )
{
  // 1st invariant, trace
  modinv(0) = prinv(0)*std::pow(prinv(2),-1./3.);
  // 2nd invariant
  modinv(1) = prinv(1)*std::pow(prinv(2),-2./3.);
  // J
  modinv(2) = std::pow(prinv(2),1./2.);

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::PlasticElastHyper::StretchesPrincipal(
  LINALG::Matrix<3,1>& prstr,
  LINALG::Matrix<3,3>& prdir,
  const LINALG::Matrix<6,1>& rcg
  )
{
  // create right Cauchy-Green 2-tensor
  LINALG::Matrix<3,3> rcgt(false);
  rcgt(0,0) = rcg(0);
  rcgt(1,1) = rcg(1);
  rcgt(2,2) = rcg(2);
  rcgt(0,1) = rcgt(1,0) = 0.5*rcg(3);
  rcgt(1,2) = rcgt(2,1) = 0.5*rcg(4);
  rcgt(2,0) = rcgt(0,2) = 0.5*rcg(5);

  // eigenvalue decomposition
  LINALG::Matrix<3,3> prstr2;  // squared principal stretches
  LINALG::SYEV(rcgt,prstr2,prdir);

  // THE principal stretches
  for (int al=0; al<3; ++al) prstr(al) = std::sqrt(prstr2(al,al));

  // bye
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::PlasticElastHyper::StretchesModified(
  LINALG::Matrix<3,1>& modstr,
  const LINALG::Matrix<3,1>& prstr
  )
{
  // determinant of deformation gradient
  const double detdefgrad = prstr(0)*prstr(1)*prstr(2);

  // determine modified principal stretches
  modstr.Update(std::pow(detdefgrad,-1.0/3.0),prstr);

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool MAT::PlasticElastHyper::HaveCoefficientsStretchesPrincipal()
{
  // set default
  bool havecoeff = false;

  // loop map of associated potential summands and see
  {
    for (unsigned int p=0; p<potsum_.size(); ++p)
    {
     havecoeff = havecoeff or potsum_[p]->HaveCoefficientsStretchesPrincipal();
    }
  }
  if (havecoeff) dserror("no principal stretch materials for plasticHyperElast (yet?)");
  // deliver
  return havecoeff;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool MAT::PlasticElastHyper::HaveCoefficientsStretchesModified()
{
  // set default
  bool havecoeff = false;

  // loop map of associated potential summands and see
  {
    for (unsigned int p=0; p<potsum_.size(); ++p)
    {
      havecoeff = havecoeff or potsum_[p]->HaveCoefficientsStretchesModified();
    }
  }
  if (havecoeff) dserror("no principal stretch materials for plasticHyperElast (yet?)");

  // deliver
  return havecoeff;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::PlasticElastHyper::Evaluate(const LINALG::Matrix<3,3>* defgrd,
    const LINALG::Matrix<3,3>* invpldefgrd,
    Teuchos::ParameterList& params,
    LINALG::Matrix<6,1>* PK2stress,
    LINALG::Matrix<6,6>* cmat_ref,
    LINALG::Matrix<6,9>* dPK2dFpinv,
    LINALG::Matrix<3,3>* MandelStress,
    LINALG::Matrix<6,6>* dMdC,
    LINALG::Matrix<6,9>* dMdFpinv,
    const int eleGID
    )
{
  LINALG::Matrix<6,1> Cpi;
  LINALG::Matrix<6,1> CpiCCpi;
  LINALG::Matrix<6,1> ircg;
  LINALG::Matrix<3,1> prinv;
  LINALG::Matrix<6,1> id2V;
  LINALG::Matrix<3,3> id2;
  LINALG::Matrix<3,3> CpiC;
  LINALG::Matrix<3,3> FpiCe;
  LINALG::Matrix<9,1> CFpiCei;
  LINALG::Matrix<9,1> CFpi;
  LINALG::Matrix<3,3> FpiTC;
  LINALG::Matrix<9,1> CFpiCe;
  LINALG::Matrix<3,3> CeFpiTC;
  LINALG::Matrix<6,1> Ce;
  LINALG::Matrix<3,3> CeM;
  LINALG::Matrix<6,1> Ce2;

  LINALG::Matrix<3,1> gamma(true);
  LINALG::Matrix<8,1> delta(true);
  LINALG::Matrix<3,1> modgamma(true);
  LINALG::Matrix<5,1> moddelta(true);

  EvaluatePlastKinQuant(defgrd,invpldefgrd,Cpi,CpiCCpi,ircg,Ce,CeM,Ce2,
                        id2V,id2,CpiC,FpiCe,CFpiCei,CFpi,FpiTC,CFpiCe,CeFpiTC,prinv);
  EvaluateGammaDelta(prinv,gamma,delta);

  // blank resulting quantities
  // ... even if it is an implicit law that cmat is zero upon input
  PK2stress->Clear();
  cmat_ref->Clear();

  // build stress response and elasticity tensor
  // for potentials based on principal invariants
  if (isoprinc_)
  {
    LINALG::Matrix<NUM_STRESS_3D,1> stressisoprinc(true) ;
    LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> cmatisoprinc(true) ;
    LINALG::Matrix<NUM_STRESS_3D,9> dPK2dFpinvIsoprinc(true);
    LINALG::Matrix<3,3> MandelStressIsoprinc(true);
    LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> dMdCisoprinc(true) ;
    LINALG::Matrix<NUM_STRESS_3D,9> dMdFpinvIsoprinc(true);
    EvaluateIsotropicPrincPlast(stressisoprinc,cmatisoprinc,dPK2dFpinvIsoprinc,
                                MandelStressIsoprinc,dMdCisoprinc,dMdFpinvIsoprinc,
                                Cpi,CpiCCpi,ircg,Ce,CeM,Ce2,id2V,id2,CpiC,FpiCe,
                                *invpldefgrd,CFpiCei,CFpi,FpiTC,CFpiCe,CeFpiTC,gamma,delta);
    PK2stress->Update(1.0, stressisoprinc, 1.0);
    cmat_ref->Update(1.0,cmatisoprinc,1.0);
    dPK2dFpinv->Update(1.,dPK2dFpinvIsoprinc,1.);
    MandelStress->Update(1.,MandelStressIsoprinc,1.);
    dMdC->Update(1.,dMdCisoprinc,1.);
    dMdFpinv->Update(1.,dMdFpinvIsoprinc,1.);
  }

  if (isomod_)
  {
    dserror("isomod_ not yet implemented");
  }

  /*----------------------------------------------------------------------*/
  //Do all the anisotropic stuff!
  if (anisoprinc_)
  {
    dserror("anisoprinc_ not yet implemented");
  }

  if (anisomod_)
  {
    dserror("anisomod_ not yet implemented");
  }
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::PlasticElastHyper::EvaluatePlastKinQuant(
    const LINALG::Matrix<3,3>* defgrd,
    const LINALG::Matrix<3,3>* invpldefgrd,
    LINALG::Matrix<6,1>& Cpi,
    LINALG::Matrix<6,1>& CpiCCpi,
    LINALG::Matrix<6,1>& ircg,
    LINALG::Matrix<6,1>& Ce,
    LINALG::Matrix<3,3>& CeM,
    LINALG::Matrix<6,1>& Ce2,
    LINALG::Matrix<6,1>& id2V,
    LINALG::Matrix<3,3>& id2,
    LINALG::Matrix<3,3>& CpiC,
    LINALG::Matrix<3,3>& FpiCe,
    LINALG::Matrix<9,1>& CFpiCei,
    LINALG::Matrix<9,1>& CFpi,
    LINALG::Matrix<3,3>& FpiTC,
    LINALG::Matrix<9,1>& CFpiCe,
    LINALG::Matrix<3,3>& CeFpiTC,
    LINALG::Matrix<3,1>& prinv)
{
  id2.Clear();
  id2V.Clear();
  for (int i=0; i<3; i++)
  {
    id2V(i)  = 1.;
    id2(i,i) = 1.;
  }
  LINALG::Matrix<3,3> tmp;
  LINALG::Matrix<3,3> tmp33;
  tmp33.Multiply(*defgrd,*invpldefgrd);
  CeM.MultiplyTN(tmp33,tmp33);
  // elastic right Cauchy-Green in strain-like Voigt notation.
  LINALG::Matrix<6,1> elasticRCGv;
  for (int i=0; i<3; i++)
    elasticRCGv(i)=CeM(i,i);
  elasticRCGv(3) = (CeM(0,1)+CeM(1,0));
  elasticRCGv(4) = (CeM(2,1)+CeM(1,2));
  elasticRCGv(5) = (CeM(0,2)+CeM(2,0));
  // elastic right Cauchy-Green in stress-like Voigt notation.
  for (int i=0; i<3; i++) Ce(i) = CeM(i,i);
  Ce(3) = (CeM(0,1)+CeM(1,0))/2.;
  Ce(4) = (CeM(2,1)+CeM(1,2))/2.;
  Ce(5) = (CeM(0,2)+CeM(2,0))/2.;

  // square of elastic right Cauchy-Green in stress-like Voigt notation.
  tmp.Multiply(CeM,CeM);
  for (int i=0; i<3; i++) Ce2(i) = tmp(i,i);
  Ce2(3) = (tmp(0,1)+tmp(1,0))/2.;
  Ce2(4) = (tmp(2,1)+tmp(1,2))/2.;
  Ce2(5) = (tmp(0,2)+tmp(2,0))/2.;

  // principal invariants of elastic Cauchy-Green strain
  InvariantsPrincipal(prinv,elasticRCGv);

  // inverse plastic right Cauchy-Green
  LINALG::Matrix<3,3> CpiM;
  CpiM.MultiplyNT(*invpldefgrd,*invpldefgrd);
  // stress-like Voigt notation
  for (int i=0; i<3; i++) Cpi(i) = CpiM(i,i);
  Cpi(3) = (CpiM(0,1)+CpiM(1,0))/2.;
  Cpi(4) = (CpiM(2,1)+CpiM(1,2))/2.;
  Cpi(5) = (CpiM(0,2)+CpiM(2,0))/2.;

  // inverse RCG
  LINALG::Matrix<3,3> iRCG;
  LINALG::Matrix<3,3> RCG;
  RCG.MultiplyTN(*defgrd,*defgrd);
  iRCG.Invert(RCG);
  // stress-like Voigt notation
  for (int i=0; i<3; i++) ircg(i) = iRCG(i,i);
  ircg(3) = (iRCG(0,1)+iRCG(1,0))/2.;
  ircg(4) = (iRCG(2,1)+iRCG(1,2))/2.;
  ircg(5) = (iRCG(0,2)+iRCG(2,0))/2.;

  // C_p^-1 * C * C_p^-1
  LINALG::Matrix<3,3> CpiCCpiM;
  tmp33.Multiply(CpiM,RCG);
  CpiCCpiM.Multiply(tmp33,CpiM);
  // stress-like Voigt notation
  for (int i=0; i<3; i++) CpiCCpi(i) = CpiCCpiM(i,i);
  CpiCCpi(3) = (CpiCCpiM(0,1)+CpiCCpiM(1,0))/2.;
  CpiCCpi(4) = (CpiCCpiM(2,1)+CpiCCpiM(1,2))/2.;
  CpiCCpi(5) = (CpiCCpiM(0,2)+CpiCCpiM(2,0))/2.;

  CpiC.Multiply(CpiM,RCG);
  FpiCe.Multiply(*invpldefgrd,CeM);

  FpiTC.MultiplyTN(*invpldefgrd,RCG);
  CeFpiTC.Multiply(CeM,FpiTC);

  tmp.Multiply(RCG,*invpldefgrd);
  Matrix3x3to9x1(tmp,CFpi);
  tmp33.Multiply(tmp,CeM);
  Matrix3x3to9x1(tmp33,CFpiCe);

  tmp.Invert(CeM);
  tmp33.Multiply(*invpldefgrd,tmp);
  tmp.Multiply(RCG,tmp33);
  Matrix3x3to9x1(tmp,CFpiCei);

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::PlasticElastHyper::EvaluateGammaDelta(
    LINALG::Matrix<3,1> prinv,
    LINALG::Matrix<3,1>& gamma,
    LINALG::Matrix<8,1>& delta
    )

{
  // principal coefficients
  if (isoprinc_)
  {
    // loop map of associated potential summands
    for (unsigned int p=0; p<potsum_.size(); ++p)
    {
      potsum_[p]->AddCoefficientsPrincipal(gamma,delta,prinv);
    }
  }

  // modified coefficients
  if (isomod_)
    dserror("Plasticelasthyper only for isotropic elasticity only in coupled form (yet?)");
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::PlasticElastHyper::EvaluateIsotropicPrincPlast(
    LINALG::Matrix<6,1>& stressisoprinc,
    LINALG::Matrix<6,6>& cmatisoprinc,
    LINALG::Matrix<6,9>& dPK2dFpinvIsoprinc,
    LINALG::Matrix<3,3>& MandelStressIsoprinc,
    LINALG::Matrix<6,6>& dMdCisoprinc,
    LINALG::Matrix<6,9>& dMdFpinvIsoprinc,
    LINALG::Matrix<6,1> Cpi,
    LINALG::Matrix<6,1> CpiCCpi,
    LINALG::Matrix<6,1> ircg,
    LINALG::Matrix<6,1> Ce,
    LINALG::Matrix<3,3> CeM,
    LINALG::Matrix<6,1> Ce2,
    LINALG::Matrix<6,1> id2V,
    LINALG::Matrix<3,3> id2,
    LINALG::Matrix<3,3> CpiC,
    LINALG::Matrix<3,3> FpiCe,
    LINALG::Matrix<3,3> Fpi,
    LINALG::Matrix<9,1> CFpiCei,
    LINALG::Matrix<9,1> CFpi,
    LINALG::Matrix<3,3> FpiTC,
    LINALG::Matrix<9,1> CFpiCe,
    LINALG::Matrix<3,3> CeFpiTC,
    LINALG::Matrix<3,1> gamma,
    LINALG::Matrix<8,1> delta
    )
{
  //  // 2nd Piola Kirchhoff stresses
    stressisoprinc.Update(gamma(0), Cpi, 1.0);
    stressisoprinc.Update(gamma(1), CpiCCpi, 1.0);
    stressisoprinc.Update(gamma(2), ircg, 1.0);

    // constitutive tensor
    cmatisoprinc.MultiplyNT(delta(0),Cpi,Cpi,1.);
    cmatisoprinc.MultiplyNT(delta(1),CpiCCpi,ircg,1.);
    cmatisoprinc.MultiplyNT(delta(1),ircg,CpiCCpi,1.);
    cmatisoprinc.MultiplyNT(delta(2),Cpi,ircg,1.);
    cmatisoprinc.MultiplyNT(delta(2),ircg,Cpi,1.);
    cmatisoprinc.MultiplyNT(delta(3),CpiCCpi,CpiCCpi,1.);
    cmatisoprinc.MultiplyNT(delta(4),CpiCCpi,ircg,1.);
    cmatisoprinc.MultiplyNT(delta(4),ircg,CpiCCpi,1.);
    cmatisoprinc.MultiplyNT(delta(5),ircg,ircg,1.);
    AddtoCmatHolzapfelProduct(cmatisoprinc,ircg,delta(6));
    AddtoCmatHolzapfelProduct(cmatisoprinc,Cpi,delta(7));

    // derivative of PK2 w.r.t. inverse plastic deformation gradient
    AddtodPK2dFpinv(dPK2dFpinvIsoprinc,id2,Fpi,gamma(0));
    AddtodPK2dFpinv(dPK2dFpinvIsoprinc,CpiC,Fpi,gamma(1));
    dPK2dFpinvIsoprinc.MultiplyNT(delta(0),Cpi,CFpi,1.);
    dPK2dFpinvIsoprinc.MultiplyNT(delta(1),Cpi,CFpiCe,1.);
    dPK2dFpinvIsoprinc.MultiplyNT(delta(1),CpiCCpi,CFpi,1.);
    dPK2dFpinvIsoprinc.MultiplyNT(delta(2),Cpi,CFpiCei,1.);
    dPK2dFpinvIsoprinc.MultiplyNT(delta(2),ircg,CFpi,1.);
    dPK2dFpinvIsoprinc.MultiplyNT(delta(3),CpiCCpi,CFpiCe,1.);
    dPK2dFpinvIsoprinc.MultiplyNT(delta(4),CpiCCpi,CFpiCei,1.);
    dPK2dFpinvIsoprinc.MultiplyNT(delta(4),ircg,CFpiCe,1.);
    dPK2dFpinvIsoprinc.MultiplyNT(delta(5),ircg,CFpiCei,1.);
    AddtodPK2dFpinv(dPK2dFpinvIsoprinc,id2,FpiCe,0.5*delta(7));

    // Mandel stress
    LINALG::Matrix<6,1> Mv;
    Mv.Update(gamma(0),Ce);
    Mv.Update(gamma(1),Ce2,1.);
    Mv.Update(gamma(2),id2V,1.);
    for (int i=0; i<3; i++) MandelStressIsoprinc(i,i) += Mv(i);
    MandelStressIsoprinc(0,1) += Mv(3);
    MandelStressIsoprinc(1,0) += Mv(3);
    MandelStressIsoprinc(1,2) += Mv(4);
    MandelStressIsoprinc(2,1) += Mv(4);
    MandelStressIsoprinc(0,2) += Mv(5);
    MandelStressIsoprinc(2,0) += Mv(5);

    // derivative of Mandel stress w.r.t. GL
    AddtodMdC_gamma1(dMdCisoprinc,Fpi,gamma(0));
    AddtodMdC_gamma2(dMdCisoprinc,Fpi,FpiCe,gamma(1));
    dMdCisoprinc.MultiplyNT(delta(0),Ce,Cpi,1.);
    dMdCisoprinc.MultiplyNT(delta(1),Ce,CpiCCpi,1.);
    dMdCisoprinc.MultiplyNT(delta(1),Ce2,Cpi,1.);
    dMdCisoprinc.MultiplyNT(delta(2),Ce,ircg,1.);
    dMdCisoprinc.MultiplyNT(delta(2),id2V,Cpi,1.);
    dMdCisoprinc.MultiplyNT(delta(3),Ce2,CpiCCpi,1.);
    dMdCisoprinc.MultiplyNT(delta(4),Ce2,ircg,1.);
    dMdCisoprinc.MultiplyNT(delta(4),id2V,CpiCCpi,1.);
    dMdCisoprinc.MultiplyNT(delta(5),id2V,ircg,1.);

    // derivative of Mandel stress w.r.t. inverse plastic deformation gradient
    AddtodPK2dFpinv(dMdFpinvIsoprinc,FpiTC,id2,gamma(0));
    AddtodPK2dFpinv(dMdFpinvIsoprinc,FpiTC,CeM,gamma(1));
    AddtodPK2dFpinv(dMdFpinvIsoprinc,CeFpiTC,id2,gamma(1));
    dMdFpinvIsoprinc.MultiplyNT(delta(0),Ce,CFpi,1.);
    dMdFpinvIsoprinc.MultiplyNT(delta(1),Ce,CFpiCe,1.);
    dMdFpinvIsoprinc.MultiplyNT(delta(1),Ce2,CFpi,1.);
    dMdFpinvIsoprinc.MultiplyNT(delta(2),Ce,CFpiCei,1.);
    dMdFpinvIsoprinc.MultiplyNT(delta(2),id2V,CFpi,1.);
    dMdFpinvIsoprinc.MultiplyNT(delta(3),Ce2,CFpiCe,1.);
    dMdFpinvIsoprinc.MultiplyNT(delta(4),Ce2,CFpiCei,1.);
    dMdFpinvIsoprinc.MultiplyNT(delta(4),id2V,CFpiCe,1.);
    dMdFpinvIsoprinc.MultiplyNT(delta(5),id2V,CFpiCei,1.);

  return ;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::PlasticElastHyper::ResponseStretches(
  LINALG::Matrix<6,6>& cmat,
  LINALG::Matrix<6,1>& stress,
  const LINALG::Matrix<6,1> rcg,
  const bool& havecoeffstrpr,
  const bool& havecoeffstrmod
  )
{
  // get principal stretches and directions
  LINALG::Matrix<3,1> prstr;
  LINALG::Matrix<3,3> prdir;
  StretchesPrincipal(prstr,prdir,rcg);
  // modified stretches
  LINALG::Matrix<3,1> modstr;
  StretchesModified(modstr,prstr);
  // determinant of deformation gradient
  const double detdefgrad = prstr(0)*prstr(1)*prstr(2);

  // get coefficients
  LINALG::Matrix<3,1> gamma_(true);
  LINALG::Matrix<6,1> delta_(true);
  if (havecoeffstrpr) {
    // loop map of associated potential summands
    for (unsigned int p=0; p<potsum_.size(); ++p)
    {
      potsum_[p]->AddCoefficientsStretchesPrincipal(gamma_,delta_,prstr);
    }
  }
  if (havecoeffstrmod) {
    // reciprocal of cubic root of determinant of deformation gradient (convenience)
    const double detdefgrad13 = std::pow(detdefgrad,-1.0/3.0);
    // retrieve coefficients with respect to modified principal stretches
    LINALG::Matrix<3,1> modgamma(true);
    LINALG::Matrix<6,1> moddelta(true);
    {
      // loop map of associated potential summands
      for (unsigned int p=0; p<potsum_.size(); ++p)
      {
        potsum_[p]->AddCoefficientsStretchesModified(modgamma,moddelta,modstr);
      }
    }
    // convert modified coefficients to oridinary counterparts
    //
    // derivatives of modified pr. stretches WRT pr. stretches
    LINALG::Matrix<3,3> modbypr(false);
    for (int al=0; al<3; ++al) {
      for (int be=0; be<3; ++be) {
        modbypr(al,be) = -modstr(al)/modstr(be);
      }
      modbypr(al,al) += 3.0;
    }
    modbypr.Scale(detdefgrad13/3.0);
    // determine unmodified coefficients gamma and add them
    gamma_.MultiplyTN(1.0,modbypr,modgamma,1.0);
    // determine unmodified coefficients delta and add them
    //
    // rewrite mod.coeff. as 2-tensor
    LINALG::Matrix<3,3> moddeltat(false);
    moddeltat(0,0) = moddelta(0);
    moddeltat(1,1) = moddelta(1);
    moddeltat(2,2) = moddelta(2);
    moddeltat(0,1) = moddeltat(1,0) = moddelta(3);
    moddeltat(1,2) = moddeltat(2,1) = moddelta(4);
    moddeltat(2,0) = moddeltat(0,2) = moddelta(5);
    // Psi_{,barlam barlam} barlam_{,lam} barlam_{,lam}
    LINALG::Matrix<3,3> aux(false);
    aux.MultiplyTN(modbypr,moddeltat);
    LINALG::Matrix<3,3> deltat(false);
    deltat.MultiplyNN(aux,modbypr);
    // Psi_{,barlam} barlam_{,lam lam}
    for (int be=0; be<3; ++be) {
      for (int ga=0; ga<3; ++ga) {
        double deltat_bega = 0.0;
        for (int al=0; al<3; ++al) {
          deltat_bega += -modgamma(al)*modbypr(al,be)/(3.0*prstr(ga));
          if (ga==al)
            deltat_bega += -modgamma(al)*detdefgrad13/(3.0*prstr(be));
          if (be==ga)
            deltat_bega += modgamma(al)*detdefgrad13*prstr(al)/(3.0*prstr(be)*prstr(be));
        }
        deltat(be,ga) += deltat_bega;
      }
    }
    // add to delta
    // Psi_{lam lam} = Psi_{,barlam barlam} barlam_{,lam} barlam_{,lam}
    //               + Psi_{,barlam} barlam_{,lam lam}
    delta_(0) += deltat(0,0);
    delta_(1) += deltat(1,1);
    delta_(2) += deltat(2,2);
    delta_(3) += deltat(0,1);
    delta_(4) += deltat(1,2);
    delta_(5) += deltat(2,0);
  }

  // principal 2nd Piola--Kirchhoff stress tensor, cf [1] Eq (6.47)
  LINALG::Matrix<3,1> prsts(true);
  for (int al=0; al<3; ++al) {
    // PK2 principal stresses
    prsts(al) = gamma_(al)/prstr(al);
    // PK2 tensor in Voigt notation
    stress(0) += prsts(al)*prdir(0,al)*prdir(0,al);  // S^11
    stress(1) += prsts(al)*prdir(1,al)*prdir(1,al);  // S^22
    stress(2) += prsts(al)*prdir(2,al)*prdir(2,al);  // S^33
    stress(3) += prsts(al)*prdir(0,al)*prdir(1,al);  // S^12
    stress(4) += prsts(al)*prdir(1,al)*prdir(2,al);  // S^23
    stress(5) += prsts(al)*prdir(2,al)*prdir(0,al);  // S^31
  }

  // integration factor prfact_{al be}
  LINALG::Matrix<6,1> prfact1(true);
  LINALG::Matrix<6,1> prfact2(true);
  for (int albe=0; albe<6; ++albe) {
    const int al = VOIGT6ROW_[albe];
    const int be = VOIGT6COL_[albe];
    double prfact1_albe = delta_(albe)/(prstr(al)*prstr(be));
    if (albe<3) prfact1_albe -= gamma_(al)/(prstr(be)*prstr(al)*prstr(al));
    prfact1(albe) = prfact1_albe;
    if (al != be) {
      if (fabs(prstr(al)-prstr(be)) < EPS6)
        prfact2(albe) = (prfact1(be) - prfact1(albe))/2.0;
      else
        prfact2(albe) = (prsts(be)-prsts(al))/(prstr(be)*prstr(be)-prstr(al)*prstr(al));
    }
  }

  // add elasticity 4-tensor, cf Holzapfel [1] Eq (6.180),(6.196)
  for (int kl=0; kl<6; ++kl) {
    const int k = VOIGT6ROW_[kl];
    const int l = VOIGT6COL_[kl];
    for (int ij=0; ij<6; ++ij) {
      const int i = VOIGT6ROW_[ij];
      const int j = VOIGT6COL_[ij];
      double c_ijkl = 0.0;
      for (int albe=0; albe<6; ++albe) {
        const int al = VOIGT6ROW_[albe];
        const int be = VOIGT6COL_[albe];
        const double fact1 = prfact1(albe);
        c_ijkl += fact1*prdir(i,al)*prdir(j,al)*prdir(k,be)*prdir(l,be);
        if (albe>=3) { // al!=be
          c_ijkl += fact1*prdir(i,be)*prdir(j,be)*prdir(k,al)*prdir(l,al);
          const double fact2 = prfact2(albe);
          c_ijkl += fact2*prdir(i,al)*prdir(j,be)*prdir(k,al)*prdir(l,be)
                  + fact2*prdir(i,al)*prdir(j,be)*prdir(k,be)*prdir(l,al)
                  + fact2*prdir(i,be)*prdir(j,al)*prdir(k,be)*prdir(l,al)
                  + fact2*prdir(i,be)*prdir(j,al)*prdir(k,al)*prdir(l,be);
        }
      }
      cmat(ij,kl) += c_ijkl;
    }
  }
  // ready
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::PlasticElastHyper::VisNames(std::map<std::string,int>& names)
{
  if (AnisotropicPrincipal() or AnisotropicModified())
  {
    std::vector<LINALG::Matrix<3,1> > fibervecs;
    GetFiberVecs(fibervecs);
    int vissize = fibervecs.size();
    std::string fiber;
    for (int i = 0; i < vissize; i++)
    {
      std::ostringstream s;
      s << "Fiber" << i+1;
      fiber = s.str();
      names[fiber] = 3; // 3-dim vector
    }
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool MAT::PlasticElastHyper::VisData(const std::string& name, std::vector<double>& data, int numgp, int eleID)
{
  if (AnisotropicPrincipal() or AnisotropicModified())
  {
    std::vector<LINALG::Matrix<3,1> > fibervecs;
    GetFiberVecs(fibervecs);
    int vissize = fibervecs.size();
    for (int i = 0; i < vissize; i++)
    {
      std::ostringstream s;
      s << "Fiber" << i+1;
      std::string fiber;
      fiber = s.str();
      if (name == fiber)
      {
        if ((int)data.size()!=3)
          dserror("size mismatch");
        data[0] = fibervecs.at(i)(0);
        data[1] = fibervecs.at(i)(1);
        data[2] = fibervecs.at(i)(2);
      }
    }
    return true;
  }
  return false;
}

void MAT::PlasticElastHyper::AddtodPK2dFpinv(LINALG::Matrix<6,9>& dPK2dFpinv,
    LINALG::Matrix<3,3> A, LINALG::Matrix<3,3> B, double fac)
{
  dPK2dFpinv(0,0) += 2 * fac * A(0,0) * B(0,0);
  dPK2dFpinv(0,3) += 2 * fac * A(0,0) * B(0,1);
  dPK2dFpinv(0,5) += 2 * fac * A(0,0) * B(0,2);
  dPK2dFpinv(0,6) += 2 * fac * A(0,1) * B(0,0);
  dPK2dFpinv(0,1) += 2 * fac * A(0,1) * B(0,1);
  dPK2dFpinv(0,4) += 2 * fac * A(0,1) * B(0,2);
  dPK2dFpinv(0,8) += 2 * fac * A(0,2) * B(0,0);
  dPK2dFpinv(0,7) += 2 * fac * A(0,2) * B(0,1);
  dPK2dFpinv(0,2) += 2 * fac * A(0,2) * B(0,2);

  dPK2dFpinv(1,0) += 2 * fac * A(1,0) * B(1,0);
  dPK2dFpinv(1,3) += 2 * fac * A(1,0) * B(1,1);
  dPK2dFpinv(1,5) += 2 * fac * A(1,0) * B(1,2);
  dPK2dFpinv(1,6) += 2 * fac * A(1,1) * B(1,0);
  dPK2dFpinv(1,1) += 2 * fac * A(1,1) * B(1,1);
  dPK2dFpinv(1,4) += 2 * fac * A(1,1) * B(1,2);
  dPK2dFpinv(1,8) += 2 * fac * A(1,2) * B(1,0);
  dPK2dFpinv(1,7) += 2 * fac * A(1,2) * B(1,1);
  dPK2dFpinv(1,2) += 2 * fac * A(1,2) * B(1,2);

  dPK2dFpinv(2,0) += 2 * fac * A(2,0) * B(2,0);
  dPK2dFpinv(2,3) += 2 * fac * A(2,0) * B(2,1);
  dPK2dFpinv(2,5) += 2 * fac * A(2,0) * B(2,2);
  dPK2dFpinv(2,6) += 2 * fac * A(2,1) * B(2,0);
  dPK2dFpinv(2,1) += 2 * fac * A(2,1) * B(2,1);
  dPK2dFpinv(2,4) += 2 * fac * A(2,1) * B(2,2);
  dPK2dFpinv(2,8) += 2 * fac * A(2,2) * B(2,0);
  dPK2dFpinv(2,7) += 2 * fac * A(2,2) * B(2,1);
  dPK2dFpinv(2,2) += 2 * fac * A(2,2) * B(2,2);

  dPK2dFpinv(3,0) += fac * (A(0,0) * B(1,0) + A(1,0) * B(0,0));
  dPK2dFpinv(3,3) += fac * (A(0,0) * B(1,1) + A(1,0) * B(0,1));
  dPK2dFpinv(3,5) += fac * (A(0,0) * B(1,2) + A(1,0) * B(0,2));
  dPK2dFpinv(3,6) += fac * (A(0,1) * B(1,0) + A(1,1) * B(0,0));
  dPK2dFpinv(3,1) += fac * (A(0,1) * B(1,1) + A(1,1) * B(0,1));
  dPK2dFpinv(3,4) += fac * (A(0,1) * B(1,2) + A(1,1) * B(0,2));
  dPK2dFpinv(3,8) += fac * (A(0,2) * B(1,0) + A(1,2) * B(0,0));
  dPK2dFpinv(3,7) += fac * (A(0,2) * B(1,1) + A(1,2) * B(0,1));
  dPK2dFpinv(3,2) += fac * (A(0,2) * B(1,2) + A(1,2) * B(0,2));

  dPK2dFpinv(4,0) += fac * (A(1,0) * B(2,0) + A(2,0) * B(1,0));
  dPK2dFpinv(4,3) += fac * (A(1,0) * B(2,1) + A(2,0) * B(1,1));
  dPK2dFpinv(4,5) += fac * (A(1,0) * B(2,2) + A(2,0) * B(1,2));
  dPK2dFpinv(4,6) += fac * (A(1,1) * B(2,0) + A(2,1) * B(1,0));
  dPK2dFpinv(4,1) += fac * (A(1,1) * B(2,1) + A(2,1) * B(1,1));
  dPK2dFpinv(4,4) += fac * (A(1,1) * B(2,2) + A(2,1) * B(1,2));
  dPK2dFpinv(4,8) += fac * (A(1,2) * B(2,0) + A(2,2) * B(1,0));
  dPK2dFpinv(4,7) += fac * (A(1,2) * B(2,1) + A(2,2) * B(1,1));
  dPK2dFpinv(4,2) += fac * (A(1,2) * B(2,2) + A(2,2) * B(1,2));

  dPK2dFpinv(5,0) += fac * (A(0,0) * B(2,0) + A(2,0) * B(0,0));
  dPK2dFpinv(5,3) += fac * (A(0,0) * B(2,1) + A(2,0) * B(0,1));
  dPK2dFpinv(5,5) += fac * (A(0,0) * B(2,2) + A(2,0) * B(0,2));
  dPK2dFpinv(5,6) += fac * (A(0,1) * B(2,0) + A(2,1) * B(0,0));
  dPK2dFpinv(5,1) += fac * (A(0,1) * B(2,1) + A(2,1) * B(0,1));
  dPK2dFpinv(5,4) += fac * (A(0,1) * B(2,2) + A(2,1) * B(0,2));
  dPK2dFpinv(5,8) += fac * (A(0,2) * B(2,0) + A(2,2) * B(0,0));
  dPK2dFpinv(5,7) += fac * (A(0,2) * B(2,1) + A(2,2) * B(0,1));
  dPK2dFpinv(5,2) += fac * (A(0,2) * B(2,2) + A(2,2) * B(0,2));

  return;
}


void MAT::PlasticElastHyper::AddtodMdC_gamma1(LINALG::Matrix<6,6>& dMdC,
    LINALG::Matrix<3,3> A, double fac)
{
  dMdC(0,0) += 2. * fac * A(0,0) * A(0,0);
  dMdC(0,3) += 2. * fac * A(0,0) * A(1,0);
  dMdC(0,5) += 2. * fac * A(0,0) * A(2,0);
  dMdC(0,1) += 2. * fac * A(1,0) * A(1,0);
  dMdC(0,4) += 2. * fac * A(1,0) * A(2,0);
  dMdC(0,2) += 2. * fac * A(2,0) * A(2,0);

  dMdC(1,0) += 2. * fac * A(0,1) * A(0,1);
  dMdC(1,3) += 2. * fac * A(0,1) * A(1,1);
  dMdC(1,5) += 2. * fac * A(0,1) * A(2,1);
  dMdC(1,1) += 2. * fac * A(1,1) * A(1,1);
  dMdC(1,4) += 2. * fac * A(1,1) * A(2,1);
  dMdC(1,2) += 2. * fac * A(2,1) * A(2,1);

  dMdC(2,0) += 2. * fac * A(0,2) * A(0,2);
  dMdC(2,3) += 2. * fac * A(0,2) * A(1,2);
  dMdC(2,5) += 2. * fac * A(0,2) * A(2,2);
  dMdC(2,1) += 2. * fac * A(1,2) * A(1,2);
  dMdC(2,4) += 2. * fac * A(1,2) * A(2,2);
  dMdC(2,2) += 2. * fac * A(2,2) * A(2,2);

  dMdC(3,0) += 2. * fac * A(0,0) * A(0,1);
  dMdC(3,3) += fac * (A(0,0) * A(1,1) + A(0,1) * A(1,0));
  dMdC(3,5) += fac * (A(0,0) * A(2,1) + A(0,1) * A(2,0));
  dMdC(3,1) += 2. * fac * A(1,0) * A(1,1);
  dMdC(3,4) += fac * (A(1,0) * A(2,1) + A(1,1) * A(2,0));
  dMdC(3,2) += 2. * fac * A(2,0) * A(2,1);

  dMdC(4,0) += 2. * fac * A(0,1) * A(0,2);
  dMdC(4,3) += fac * (A(0,1) * A(1,2) + A(0,2) * A(1,1));
  dMdC(4,5) += fac * (A(0,1) * A(2,2) + A(0,2) * A(2,1));
  dMdC(4,1) += 2. * fac * A(1,1) * A(1,2);
  dMdC(4,4) += fac * (A(1,1) * A(2,2) + A(1,2) * A(2,1));
  dMdC(4,2) += 2. * fac * A(2,1) * A(2,2);

  dMdC(5,0) += 2. * fac * A(0,0) * A(0,2);
  dMdC(5,3) += fac * (A(0,0) * A(1,2) + A(0,2) * A(1,0));
  dMdC(5,5) += fac * (A(0,0) * A(2,2) + A(0,2) * A(2,0));
  dMdC(5,1) += 2. * fac * A(1,0) * A(1,2);
  dMdC(5,4) += fac * (A(1,0) * A(2,2) + A(1,2) * A(2,0));
  dMdC(5,2) += 2. * fac * A(2,0) * A(2,2);

  return;
}


void MAT::PlasticElastHyper::AddtodMdC_gamma2(LINALG::Matrix<6,6>& dMdC,
    LINALG::Matrix<3,3> A, LINALG::Matrix<3,3> B, double fac)
{
  dMdC(0,0) += 4 * fac * A(0,0) * B(0,0);
  dMdC(0,3) += fac * (2 * A(0,0) * B(1,0) + 2 * A(1,0) * B(0,0));
  dMdC(0,5) += fac * (2 * A(0,0) * B(2,0) + 2 * A(2,0) * B(0,0));
  dMdC(0,1) += 4 * fac * A(1,0) * B(1,0);
  dMdC(0,4) += fac * (2 * A(1,0) * B(2,0) + 2 * A(2,0) * B(1,0));
  dMdC(0,2) += 4 * fac * A(2,0) * B(2,0);

  dMdC(3,0) += fac * (2 * A(0,0) * B(0,1) + 2 * A(0,1) * B(0,0));
  dMdC(3,3) += fac * (A(0,0) * B(1,1) + A(1,0) * B(0,1) + A(1,1) * B(0,0) + A(0,1) * B(1,0));
  dMdC(3,5) += fac * (A(0,0) * B(2,1) + A(2,0) * B(0,1) + A(2,1) * B(0,0) + A(0,1) * B(2,0));
  dMdC(3,1) += fac * (2 * A(1,0) * B(1,1) + 2 * A(1,1) * B(1,0));
  dMdC(3,4) += fac * (A(1,0) * B(2,1) + A(2,0) * B(1,1) + A(2,1) * B(1,0) + A(1,1) * B(2,0));
  dMdC(3,2) += fac * (2 * A(2,0) * B(2,1) + 2 * A(2,1) * B(2,0));

  dMdC(5,0) += fac * (2 * A(0,0) * B(0,2) + 2 * A(0,2) * B(0,0));
  dMdC(5,3) += fac * (A(0,0) * B(1,2) + A(1,0) * B(0,2) + A(1,2) * B(0,0) + A(0,2) * B(1,0));
  dMdC(5,5) += fac * (A(0,0) * B(2,2) + A(2,0) * B(0,2) + A(2,2) * B(0,0) + A(0,2) * B(2,0));
  dMdC(5,1) += fac * (2 * A(1,0) * B(1,2) + 2 * A(1,2) * B(1,0));
  dMdC(5,4) += fac * (A(1,0) * B(2,2) + A(2,0) * B(1,2) + A(2,2) * B(1,0) + A(1,2) * B(2,0));
  dMdC(5,2) += fac * (2 * A(2,0) * B(2,2) + 2 * A(2,2) * B(2,0));

  dMdC(1,0) += 4 * fac * A(0,1) * B(0,1);
  dMdC(1,3) += fac * (2 * A(0,1) * B(1,1) + 2 * A(1,1) * B(0,1));
  dMdC(1,5) += fac * (2 * A(0,1) * B(2,1) + 2 * A(2,1) * B(0,1));
  dMdC(1,1) += 4 * fac * A(1,1) * B(1,1);
  dMdC(1,4) += fac * (2 * A(1,1) * B(2,1) + 2 * A(2,1) * B(1,1));
  dMdC(1,2) += 4 * fac * A(2,1) * B(2,1);

  dMdC(4,0) += fac * (2 * A(0,1) * B(0,2) + 2 * A(0,2) * B(0,1));
  dMdC(4,3) += fac * (A(0,1) * B(1,2) + A(1,1) * B(0,2) + A(1,2) * B(0,1) + A(0,2) * B(1,1));
  dMdC(4,5) += fac * (A(0,1) * B(2,2) + A(2,1) * B(0,2) + A(2,2) * B(0,1) + A(0,2) * B(2,1));
  dMdC(4,1) += fac * (2 * A(1,1) * B(1,2) + 2 * A(1,2) * B(1,1));
  dMdC(4,4) += fac * (A(1,1) * B(2,2) + A(2,1) * B(1,2) + A(2,2) * B(1,1) + A(1,2) * B(2,1));
  dMdC(4,2) += fac * (2 * A(2,1) * B(2,2) + 2 * A(2,2) * B(2,1));

  dMdC(2,0) += 4 * fac * A(0,2) * B(0,2);
  dMdC(2,3) += fac * (2 * A(0,2) * B(1,2) + 2 * A(1,2) * B(0,2));
  dMdC(2,5) += fac * (2 * A(0,2) * B(2,2) + 2 * A(2,2) * B(0,2));
  dMdC(2,1) += 4 * fac * A(1,2) * B(1,2);
  dMdC(2,4) += fac * (2 * A(1,2) * B(2,2) + 2 * A(2,2) * B(1,2));
  dMdC(2,2) += 4 * fac * A(2,2) * B(2,2);

  return;
}

void MAT::PlasticElastHyper::Matrix3x3to9x1(LINALG::Matrix<3,3> A, LINALG::Matrix<9,1>& Out)
{
  Out(0) = A(0,0);
  Out(1) = A(1,1);
  Out(2) = A(2,2);
  Out(3) = A(0,1);
  Out(4) = A(1,2);
  Out(5) = A(0,2);
  Out(6) = A(1,0);
  Out(7) = A(2,1);
  Out(8) = A(2,0);
  return;
}
