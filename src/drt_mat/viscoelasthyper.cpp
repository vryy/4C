/*----------------------------------------------------------------------*/
/*!
\file viscoelasthyper.cpp
\brief
This file contains the viscohyperelastic material.
This model can be applied to any hyperelastic law of the Elasthyper toolbox.
The viscos part is rate-dependent and summed up with the hyperelastic laws
to build a viscohyperelastic strain energy function.
(Description of hysteresis not added jet)

The input line should read
MAT 0   MAT_ViscoElastHyper   NUMMAT 0 MATIDS  DENS 0 GAMMA 0 INIT_MODE -1

<pre>
Maintainer: Anna Birzle
            anna.birzle@tum.de
            http://www.lnm.mw.tum.de
</pre>
*/

/*----------------------------------------------------------------------*/
#include "viscoelasthyper.H"
#include "../drt_lib/standardtypes_cpp.H"
#include "../drt_matelast/elast_summand.H"
#include "../linalg/linalg_utils.H"
#include "../drt_lib/drt_linedefinition.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_mat/material_service.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::PAR::ViscoElastHyper::ViscoElastHyper(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: Parameter(matdata),
  nummat_(matdata->GetInt("NUMMAT")),
  matids_(matdata->Get<std::vector<int> >("MATIDS")),
  density_(matdata->GetDouble("DENS"))
{
  // check if sizes fit
  if (nummat_ != (int)matids_->size())
    dserror("number of materials %d does not fit to size of material vector %d", nummat_, matids_->size());
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<MAT::Material> MAT::PAR::ViscoElastHyper::CreateMaterial()
{
  return Teuchos::rcp(new MAT::ViscoElastHyper(this));
}


MAT::ViscoElastHyperType MAT::ViscoElastHyperType::instance_;


DRT::ParObject* MAT::ViscoElastHyperType::Create( const std::vector<char> & data )
{
  MAT::ViscoElastHyper* elhy = new MAT::ViscoElastHyper();
  elhy->Unpack(data);

  return elhy;
}


/*----------------------------------------------------------------------*
 |  initialise static arrays                                 bborn 08/09|
 *----------------------------------------------------------------------*/
// 6-Voigt C-index                              0 1 2  3 4 5
const int MAT::ViscoElastHyper::VOIGT6ROW_[6] = {0,1,2, 0,1,2};
const int MAT::ViscoElastHyper::VOIGT6COL_[6] = {0,1,2, 1,2,0};

// tensor indices ij = 11, 12, 13, 21, 22, 23, 31, 32, 33
// C indices           00, 01, 02, 10, 11, 12, 20, 21, 22
// Access : 3*i+j
// 6-Voigt C-indices    0   3   5   3   1   4   5   4   2
const int MAT::ViscoElastHyper::VOIGT3X3SYM_[9] = {0,3,5, 3,1,4, 5,4,2};


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::ViscoElastHyper::ViscoElastHyper()
  : params_(NULL),
    potsum_(0)
{
	// history data 09/13
	isinitvis_=false;
	histmodrcgcurr_=Teuchos::rcp(new std::vector<LINALG::Matrix<6,1> >);
	histmodrcglast_=Teuchos::rcp(new std::vector<LINALG::Matrix<6,1> >);

}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::ViscoElastHyper::ViscoElastHyper(MAT::PAR::ViscoElastHyper* params)
  : params_(params),
    potsum_(0)

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
void MAT::ViscoElastHyper::Pack(DRT::PackBuffer& data) const
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


  if (params_ != NULL) // summands are not accessible in postprocessing mode
  {
    // loop map of associated potential summands
    for (unsigned int p=0; p<potsum_.size(); ++p)
    {
     potsum_[p]->PackSummand(data);
    }
  }

  //  pack history data 09/13
    int histsize;
    if (!Initialized())
    {
      histsize=0;
    }
    else
    {
      histsize = histmodrcglast_->size();
    }
    AddtoPack(data,histsize);  // Length of history vector(s)
    for (int var = 0; var < histsize; ++var)
      AddtoPack(data,histmodrcglast_->at(var));

}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ViscoElastHyper::Unpack(const std::vector<char>& data)
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
        params_ = static_cast<MAT::PAR::ViscoElastHyper*>(mat);
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

    // history data 09/13
    isinitvis_ = true;
    int histsize;
    ExtractfromPack(position,data,histsize);

    if (histsize == 0) isinitvis_=false;

    histmodrcgcurr_=Teuchos::rcp(new std::vector<LINALG::Matrix<6,1> >);
    histmodrcglast_=Teuchos::rcp(new std::vector<LINALG::Matrix<6,1> >);

    for (int var=0; var<histsize; var+=1)
    {
      LINALG::Matrix<6,1> tmp(true);
      histmodrcgcurr_->push_back(tmp);

      ExtractfromPack(position,data,tmp);
      histmodrcglast_->push_back(tmp);
    }

    // in the postprocessing mode, we do not unpack everything we have packed
    // -> position check cannot be done in this case
    if (position != data.size())
      dserror("Mismatch in size of data %d <-> %d",data.size(),position);
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
int MAT::ViscoElastHyper::MatID(
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
double MAT::ViscoElastHyper::ShearMod() const
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
void MAT::ViscoElastHyper::SetupAAA(Teuchos::ParameterList& params)
{
  // loop map of associated potential summands
  for (unsigned int p=0; p<potsum_.size(); ++p)
  {
    potsum_[p]->SetupAAA(params);
  }

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ViscoElastHyper::Setup(int numgp, DRT::INPUT::LineDefinition* linedef)
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
  anisomod_ = false ;
  isomodvisco_ = false ;

  for (unsigned int p=0; p<potsum_.size(); ++p)
  {
    potsum_[p]->SpecifyFormulation(isoprinc_, isomod_, anisoprinc_, anisomod_, isomodvisco_);
  }

  // Initialise/allocate history variables 09/13

  histmodrcgcurr_=Teuchos::rcp(new std::vector<LINALG::Matrix<6,1> >);
  histmodrcglast_=Teuchos::rcp(new std::vector<LINALG::Matrix<6,1> >);

  const LINALG::Matrix<6,1> emptyvec(true);
  histmodrcgcurr_->resize(numgp);
  histmodrcglast_->resize(numgp);

  for (int j=0; j<numgp; ++j)
  {
    histmodrcgcurr_->at(j) = emptyvec;
    histmodrcglast_->at(j) = emptyvec;
  }
  isinitvis_ = true;

  return;
}

/*------------------------------------------------------------------------------------------*
|  Setup internal stress variables - special for the inverse analysis (public)         09/13|
*-------------------------------------------------------------------------------------------*/
void MAT::ViscoElastHyper::ResetAll(const int numgp)
{
  // Initialise/allocate history variables 09/13

  histmodrcgcurr_=Teuchos::rcp(new std::vector<LINALG::Matrix<6,1> >);
  histmodrcglast_=Teuchos::rcp(new std::vector<LINALG::Matrix<6,1> >);

  const LINALG::Matrix<6,1> emptyvec(true);
  histmodrcgcurr_->resize(numgp);
  histmodrcglast_->resize(numgp);

  for (int j=0; j<numgp; ++j)
  {
    histmodrcgcurr_->at(j) = emptyvec;
    histmodrcglast_->at(j) = emptyvec;
  }
  isinitvis_ = true;

  return ;
}

/*-------------------------------------------------------------------------------*
|  Setup time variables - special for the inverse analysis (public)         09/13|
*-------------------------------------------------------------------------------*/
void MAT::ViscoElastHyper::SetupTimeVariables()
{
  return ;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void MAT::ViscoElastHyper::Update()
{
  // loop map of associated potential summands
  for (unsigned int p=0; p<potsum_.size(); ++p)
  {
    potsum_[p]->Update();
  }

  // Update history values 09/13
  histmodrcglast_=histmodrcgcurr_;

  const LINALG::Matrix<6,1> emptyvec(true);
  histmodrcgcurr_=Teuchos::rcp(new std::vector<LINALG::Matrix<6,1> >);

  const int numgp=histmodrcglast_->size();
  histmodrcgcurr_->resize(numgp);

  for (int j=0; j<numgp; ++j)
  {
    histmodrcgcurr_->at(j) = emptyvec;
  }

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ViscoElastHyper::GetFiberVecs(std::vector<LINALG::Matrix<3,1> >& fibervecs)
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
void MAT::ViscoElastHyper::EvaluateFiberVecs(
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
void MAT::ViscoElastHyper::InvariantsPrincipal(
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
void MAT::ViscoElastHyper::InvariantsModified(
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
void MAT::ViscoElastHyper::StretchesPrincipal(
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
void MAT::ViscoElastHyper::StretchesModified(
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
bool MAT::ViscoElastHyper::HaveCoefficientsStretchesPrincipal()
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

  // deliver
  return havecoeff;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool MAT::ViscoElastHyper::HaveCoefficientsStretchesModified()
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

  // deliver
  return havecoeff;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ViscoElastHyper::StrainEnergy(const LINALG::Matrix<6,1>& glstrain,
                                   double& psi)
{
  LINALG::Matrix<6,1> id2(true) ;
  LINALG::Matrix<6,1> rcg(true) ;
  LINALG::Matrix<6,1> scg(true) ;
  LINALG::Matrix<6,1> icg(true) ;
  LINALG::Matrix<6,6> id4(true) ;
  LINALG::Matrix<6,6> id4sharp(true) ;

  LINALG::Matrix<3,1> prinv(true);
  LINALG::Matrix<3,1> modinv(true);

  // evluate kinematic quantities
  EvaluateKinQuant(glstrain,id2,scg,rcg,icg,id4,id4sharp,prinv,modinv);

  // loop map of associated potential summands
  for (unsigned int p=0; p<potsum_.size(); ++p)
  {
    potsum_[p]->AddStrainEnergy(psi,prinv,modinv);
  }

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ViscoElastHyper::Evaluate(const LINALG::Matrix<3,3>* defgrd,
                               const LINALG::Matrix<6,1>* glstrain,
                               Teuchos::ParameterList& params,
                               LINALG::Matrix<6,1>* stress,
                               LINALG::Matrix<6,6>* cmat)
{

  LINALG::Matrix<6,1> id2(true) ;
  LINALG::Matrix<6,1> rcg(true) ;
  LINALG::Matrix<6,1> modrcg(true) ;
  LINALG::Matrix<6,1> scg(true) ;
  LINALG::Matrix<6,1> icg(true) ;
  LINALG::Matrix<6,6> id4(true) ;
  LINALG::Matrix<6,6> id4sharp(true) ;

  LINALG::Matrix<3,1> prinv(true);
  LINALG::Matrix<3,1> modinv(true);
  LINALG::Matrix<7,1> modrateinv (true);

  LINALG::Matrix<3,1> gamma(true);
  LINALG::Matrix<8,1> delta(true);
  LINALG::Matrix<3,1> modgamma(true);
  LINALG::Matrix<5,1> moddelta(true);

  LINALG::Matrix<6,1> modrcgrate(true);
  // for extension: LINALG::Matrix<6,1> modicgrate(true);
  LINALG::Matrix<8,1> modmy(true);
  LINALG::Matrix<33,1> modxi(true);

  EvaluateKinQuant(*glstrain,id2,scg,rcg,icg,id4,id4sharp,prinv,modinv);
  EvaluateGammaDelta(rcg,prinv,modinv,gamma,delta,modgamma,moddelta);
  if (isomodvisco_)
  {
	  // calculate viscos quantities
	  EvaluateKinQuantVis(rcg,modrcg,icg,params,prinv,modrcgrate,modrateinv);
	  EvaluateMyXi(modinv,modmy,modxi,modrateinv,params);
  }

  // blank resulting quantities
  // ... even if it is an implicit law that cmat is zero upon input
  stress->Clear();
  cmat->Clear();

  // build stress response and elasticity tensor
  // for potentials based on principal invariants
  // (exists not for viscos part)
  if (isoprinc_)
  {
    LINALG::Matrix<NUM_STRESS_3D,1> stressisoprinc(true) ;
    LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> cmatisoprinc(true) ;
    EvaluateIsotropicPrinc(stressisoprinc,cmatisoprinc,scg,id2,icg,id4sharp,gamma,delta);
    stress->Update(1.0, stressisoprinc, 1.0);
    cmat->Update(1.0,cmatisoprinc,1.0);
  }
  // for decoupled potentials
  if (isomod_)
  {
	  if (!isomodvisco_)
	  {
		  // just hyperelastic part
		  LINALG::Matrix<NUM_STRESS_3D,1> stressisomodiso(true) ;
		  LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> cmatisomodiso(true);
		  LINALG::Matrix<NUM_STRESS_3D,1> stressisomodvol(true) ;
		  LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> cmatisomodvol(true) ;
		  EvaluateIsotropicMod(stressisomodiso,stressisomodvol,cmatisomodiso,cmatisomodvol,rcg,id2,icg,id4,id4sharp,modinv,prinv,modgamma,moddelta);
		  stress->Update(1.0, stressisomodiso, 1.0);
		  stress->Update(1.0, stressisomodvol, 1.0);
		  cmat->Update(1.0,cmatisomodiso,1.0);
		  cmat->Update(1.0,cmatisomodvol,1.0);
	  }
	  else
	  {
		  // first hyperelastic part
		  LINALG::Matrix<NUM_STRESS_3D,1> stressisomodiso(true) ;
		  LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> cmatisomodiso(true);
		  LINALG::Matrix<NUM_STRESS_3D,1> stressisomodvol(true) ;
		  LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> cmatisomodvol(true) ;
		  EvaluateIsotropicMod(stressisomodiso,stressisomodvol,cmatisomodiso,cmatisomodvol,rcg,id2,icg,id4,id4sharp,modinv,prinv,modgamma,moddelta);
		  stress->Update(1.0, stressisomodiso, 1.0);
		  stress->Update(1.0, stressisomodvol, 1.0);
		  cmat->Update(1.0,cmatisomodiso,1.0);
		  cmat->Update(1.0,cmatisomodvol,1.0);

		  // add viscos part
		  LINALG::Matrix<NUM_STRESS_3D,1> stressisomodisovisco(true);
		  LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> cmatisomodisovisco(true);
		  LINALG::Matrix<NUM_STRESS_3D,1> stressisomodvolvisco(true) ;
		  LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> cmatisomodvolvisco(true) ;
		  EvaluateIsoModVisco(stressisomodisovisco,stressisomodvolvisco,cmatisomodisovisco,cmatisomodvolvisco,prinv,modinv,modmy,modxi,rcg,id2,icg,id4,modrcgrate);
		  stress->Update(1.0, stressisomodisovisco, 1.0);
		  stress->Update(1.0, stressisomodvolvisco, 1.0);
		  cmat->Update(1.0,cmatisomodisovisco,1.0);
		  cmat->Update(1.0,cmatisomodvolvisco,1.0);
	  }
  }


  /*----------------------------------------------------------------------*/
  // coefficients in principal stretches
  const bool havecoeffstrpr = HaveCoefficientsStretchesPrincipal();
  const bool havecoeffstrmod = HaveCoefficientsStretchesModified();
  if (havecoeffstrpr or havecoeffstrmod) {
    ResponseStretches(*cmat,*stress,rcg,havecoeffstrpr,havecoeffstrmod);
  }

  /*----------------------------------------------------------------------*/
  //Do all the anisotropic stuff!
  if (anisoprinc_)
  {
      LINALG::Matrix<NUM_STRESS_3D,1> stressanisoprinc(true) ;
      LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> cmatanisoprinc(true) ;
      EvaluateAnisotropicPrinc(stressanisoprinc,cmatanisoprinc,rcg,params);
      stress->Update(1.0, stressanisoprinc, 1.0);
      cmat->Update(1.0, cmatanisoprinc, 1.0);
  }

  if (anisomod_)
  {
      LINALG::Matrix<NUM_STRESS_3D,1> stressanisomod(true) ;
      LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> cmatanisomod(true) ;
      EvaluateAnisotropicMod(stressanisomod,cmatanisomod,rcg,icg,prinv);
      stress->Update(1.0, stressanisomod, 1.0);
      cmat->Update(1.0, cmatanisomod, 1.0);
  }

  return ;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ViscoElastHyper::EvaluateKinQuant(
    const LINALG::Matrix<6,1>& glstrain,
    LINALG::Matrix<6,1>& id2,
    LINALG::Matrix<6,1>& scg,
    LINALG::Matrix<6,1>& rcg,
    LINALG::Matrix<6,1>& icg,
    LINALG::Matrix<6,6>& id4,
    LINALG::Matrix<6,6>& id4sharp,
    LINALG::Matrix<3,1>& prinv,
    LINALG::Matrix<3,1>& modinv)

{
  // build Cartesian identity 2-tensor I_{AB}
  for (int i=0; i<3; i++) id2(i) = 1.0;

  // right Cauchy-Green Tensor  C_{AB} = 2 * E_{AB} + I_{AB}
  // REMARK: strain-like 6-Voigt vector
  rcg.Update(2.0,glstrain,1.0);
  rcg.Update(1.0, id2, 1.0);

  // 'contra-variant' right Cauchy-Green Tensor C^{AB}
  // REMARK: stress-like 6-Voigt vector of right CG
  scg.Update(1.0,rcg,1.0);
  for (int i=3; i<6; i++) scg(i) *= 0.5;

  // principal invariants of right Cauchy-Green strain
  InvariantsPrincipal(prinv,rcg);

  // invert right Cauchy-Green tensor
  // REMARK: stress-like 6-Voigt vector
  {
    icg(0) = ( rcg(1)*rcg(2) - 0.25*rcg(4)*rcg(4) ) / prinv(2);
    icg(1) = ( rcg(0)*rcg(2) - 0.25*rcg(5)*rcg(5) ) / prinv(2);
    icg(2) = ( rcg(0)*rcg(1) - 0.25*rcg(3)*rcg(3) ) / prinv(2);
    icg(3) = ( 0.25*rcg(5)*rcg(4) - 0.5*rcg(3)*rcg(2) ) / prinv(2);
    icg(4) = ( 0.25*rcg(3)*rcg(5) - 0.5*rcg(0)*rcg(4) ) / prinv(2);
    icg(5) = ( 0.25*rcg(3)*rcg(4) - 0.5*rcg(5)*rcg(1) ) / prinv(2);
  }

  // set Cartesian identity 4-tensor in 6-Voigt matrix notation
  // this is fully 'contra-variant' identity tensor, ie I^{ABCD}
  // REMARK: rows are stress-like 6-Voigt
  //         columns are stress-like 6-Voigt
  for (int i=0; i<3; i++) id4sharp(i,i) = 1.0;
  for (int i=3; i<6; i++) id4sharp(i,i) = 0.5;

  // set Cartesian identity 4-tensor in 6x6-matrix notation (stress-like)
  // this is a 'mixed co- and contra-variant' identity 4-tensor, ie I^{AB}_{CD}
  // REMARK: rows are stress-like 6-Voigt
  //         columns are strain-like 6-Voigt
  for (int i=0; i<6; i++) id4(i,i) = 1.0;

  // modified invariants of right Cauchy-Green strain
  InvariantsModified(modinv,prinv);
}

/*----------------------------------------------------------------------*/
/* Evaluate Quantities for Viscos Part							09/13	*/
/*----------------------------------------------------------------------*/
void MAT::ViscoElastHyper::EvaluateKinQuantVis(
    LINALG::Matrix<6,1> rcg,
    LINALG::Matrix<6,1>& modrcg,
    LINALG::Matrix<6,1> icg,
    Teuchos::ParameterList& params,
    LINALG::Matrix<3,1> prinv,
    LINALG::Matrix<6,1>& modrcgrate,
    LINALG::Matrix<7,1>& modrateinv
    )

{
	  // modrcg : \overline{C} = J^{-\frac{2}{3}} C
	  const double modscale = std::pow(prinv(2),-1./3.);
	  modrcg.Update(modscale,rcg);

	  // get gauss point number of this element
	  const int gp = params.get<int>("gp",-1);
	  // get time algorithmic parameters
	  double dt = params.get<double>("delta time");

      // read history
	  LINALG::Matrix<6,1> modrcglast (histmodrcglast_->at(gp));

	  // rate of Cauchy-Green Tensor \dot{C} = \frac {overline{C}^n - \overline{C}^{n-1}} {\Delta t}
	  // REMARK: strain-like 6-Voigt vector
	  modrcgrate.Update(1.0,modrcg,1.0);
	  modrcgrate.Update(-1.0,modrcglast,1.0);
	  modrcgrate.Scale(1/dt);

	  // in the first time step, set modrcgrate to zero (--> first time step is just hyperelastic, not viscos)
	  const LINALG::Matrix<6,1> emptyvec(true);
	  if(modrcglast == emptyvec)
	  {
		  modrcgrate=emptyvec;
	  }

	  // Update history of Cauchy-Green Tensor
      histmodrcgcurr_->at(gp) = modrcg;

      // Second Invariant of modrcgrate \bar{J}_2 = \frac{1}{2} \tr (\dot{\overline{C^2}}
      modrateinv(1) = 0.5*( modrcgrate(0)*modrcgrate(0) + modrcgrate(1)*modrcgrate(1) + modrcgrate(2)*modrcgrate(2)
    		                   + .5*modrcgrate(3)*modrcgrate(3) + .5*modrcgrate(4)*modrcgrate(4) + .5*modrcgrate(5)*modrcgrate(5) );

      // For further extension of material law (not necassary at the moment)
      /*
      // necassary transfer variable: LINALG::Matrix<6,1>& modicgrate
      // \overline{J}_3 = determinant of modified rate of right Cauchy-Green-Tensor
      modrateinv(2) = modrcgrate(0)*modrcgrate(1)*modrcgrate(2)
        + 0.25 * modrcgrate(3)*modrcgrate(4)*modrcgrate(5)
        - 0.25 * modrcgrate(1)*modrcgrate(5)*modrcgrate(5)
        - 0.25 * modrcgrate(2)*modrcgrate(3)*modrcgrate(3)
        - 0.25 * modrcgrate(0)*modrcgrate(4)*modrcgrate(4);

      // invert modified rate of right Cauchy-Green tensor
      // REMARK: stress-like 6-Voigt vector
      {
        modicgrate(0) = ( modrcgrate(1)*modrcgrate(2) - 0.25*modrcgrate(4)*modrcgrate(4) ) / modrateinv(2);
        modicgrate(1) = ( modrcgrate(0)*modrcgrate(2) - 0.25*modrcgrate(5)*modrcgrate(5) ) / modrateinv(2);
        modicgrate(2) = ( modrcgrate(0)*modrcgrate(1) - 0.25*modrcgrate(3)*modrcgrate(3) ) / modrateinv(2);
        modicgrate(3) = ( 0.25*modrcgrate(5)*modrcgrate(4) - 0.5*modrcgrate(3)*modrcgrate(2) ) / modrateinv(2);
        modicgrate(4) = ( 0.25*modrcgrate(3)*modrcgrate(5) - 0.5*modrcgrate(0)*modrcgrate(4) ) / modrateinv(2);
        modicgrate(5) = ( 0.25*modrcgrate(3)*modrcgrate(4) - 0.5*modrcgrate(5)*modrcgrate(1) ) / modrateinv(2);
      }*/
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ViscoElastHyper::EvaluateGammaDelta(
    const LINALG::Matrix<6,1> rcg,
    LINALG::Matrix<3,1> prinv,
    LINALG::Matrix<3,1> modinv,
    LINALG::Matrix<3,1>& gamma,
    LINALG::Matrix<8,1>& delta,
    LINALG::Matrix<3,1>& modgamma,
    LINALG::Matrix<5,1>& moddelta
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
  {

    // loop map of associated potential summands
    for (unsigned int p=0; p<potsum_.size(); ++p)
    {
      potsum_[p]->AddCoefficientsModified(modgamma,moddelta,modinv);
    }
  }
}

/*----------------------------------------------------------------------*/
/* Evaluate Factors for Viscos Quantities					09/13	*/
/*----------------------------------------------------------------------*/
void MAT::ViscoElastHyper::EvaluateMyXi(
    LINALG::Matrix<3,1> modinv,
    LINALG::Matrix<8,1>& modmy,
    LINALG::Matrix<33,1>& modxi,
    LINALG::Matrix<7,1>& modrateinv,
    Teuchos::ParameterList& params
    )
{
  // modified coefficients
  // loop map of associated potential summands
  for (unsigned int p=0; p<potsum_.size(); ++p)
  {
	  potsum_[p]->AddCoefficientsViscoModified(modinv,modmy,modxi,modrateinv,params);
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ViscoElastHyper::EvaluateIsotropicPrinc(
    LINALG::Matrix<6,1>& stressisoprinc,
    LINALG::Matrix<6,6>& cmatisoprinc,
    LINALG::Matrix<6,1> scg,
    LINALG::Matrix<6,1> id2,
    LINALG::Matrix<6,1> icg,
    LINALG::Matrix<6,6> id4sharp,
    LINALG::Matrix<3,1> gamma,
    LINALG::Matrix<8,1> delta
    )
{

  // 2nd Piola Kirchhoff stresses
  stressisoprinc.Update(gamma(0), id2, 1.0);
  stressisoprinc.Update(gamma(1), scg, 1.0);
  stressisoprinc.Update(gamma(2), icg, 1.0);

  // constitutive tensor
  // contribution: Id \otimes Id
  cmatisoprinc.MultiplyNT(delta(0), id2, id2, 1.0);
  // contribution: Id \otimes C + C \otimes Id
  cmatisoprinc.MultiplyNT(delta(1), id2, scg, 1.0);
  cmatisoprinc.MultiplyNT(delta(1), scg, id2, 1.0);
  // contribution: Id \otimes Cinv + Cinv \otimes Id
  cmatisoprinc.MultiplyNT(delta(2), id2, icg, 1.0);
  cmatisoprinc.MultiplyNT(delta(2), icg, id2, 1.0);
  // contribution: C \otimes C
  cmatisoprinc.MultiplyNT(delta(3), scg, scg, 1.0);
  // contribution: C \otimes Cinv + Cinv \otimes C
  cmatisoprinc.MultiplyNT(delta(4), scg, icg, 1.0);
  cmatisoprinc.MultiplyNT(delta(4), icg, scg, 1.0);
  // contribution: Cinv \otimes Cinv
  cmatisoprinc.MultiplyNT(delta(5), icg, icg, 1.0);
  // contribution: Cinv \odot Cinv
  AddtoCmatHolzapfelProduct(cmatisoprinc, icg, delta(6));
  // contribution: Id4^#
  cmatisoprinc.Update(delta(7), id4sharp, 1.0);

  return ;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ViscoElastHyper::EvaluateIsotropicMod(
    LINALG::Matrix<6,1>& stressisomodiso,
    LINALG::Matrix<6,1>& stressisomodvol,
    LINALG::Matrix<6,6>& cmatisomodiso,
    LINALG::Matrix<6,6>& cmatisomodvol,
    LINALG::Matrix<6,1> rcg,
    LINALG::Matrix<6,1> id2,
    LINALG::Matrix<6,1> icg,
    LINALG::Matrix<6,6> id4,
    LINALG::Matrix<6,6> id4sharp,
    LINALG::Matrix<3,1> modinv,
    LINALG::Matrix<3,1> prinv,
    LINALG::Matrix<3,1> modgamma,
    LINALG::Matrix<5,1> moddelta
    )
{
  // define necessary variables
  const double modscale = std::pow(prinv(2),-1./3.);
  // modified right Cauchy-Green
  LINALG::Matrix<6,1> modrcg(true);
  modrcg.Update(modscale,rcg);

  // 2nd Piola Kirchhoff stresses

  // isochoric contribution
  LINALG::Matrix<6,1> modstress(true);
  modstress.Update(modgamma(0), id2);
  modstress.Update(modgamma(1), modrcg, 1.0);

  // build 4-tensor for projection as 6x6 tensor
  LINALG::Matrix<6,6> Projection;
  Projection.MultiplyNT(1./3., icg, rcg);
  Projection.Update(1.0, id4, -1.0);
  // isochoric stress
  LINALG::Matrix<6,1> isostress(true);
  stressisomodiso.MultiplyNN(modscale,Projection,modstress,1.0);

  // volumetric contribution
  stressisomodvol.Update(modgamma(2)*modinv(2), icg, 1.0);

  // constitutive tensor

  //isochoric contribution
  // modified constitutive tensor
  LINALG::Matrix<6,6> modcmat;
  LINALG::Matrix<6,6> modcmat2(true);
  // contribution: Id \otimes Id
  modcmat.MultiplyNT(moddelta(0), id2, id2);
  // contribution: Id \otimes C + C \otimes Id
  modcmat.MultiplyNT(moddelta(1), id2, modrcg, 1.0);
  modcmat.MultiplyNT(moddelta(1), rcg, id2, 1.0);
  // contribution: C \otimes C
  modcmat.MultiplyNT(moddelta(2), rcg, modrcg, 1.0);
  // contribution: Id4^#
  modcmat.Update(moddelta(3), id4sharp, 1.0);
  //scaling
  modcmat.Scale(std::pow(modinv(2),-4./3.));
  //contribution: P:modC:P
  modcmat2.MultiplyNN(Projection,modcmat);
  cmatisomodiso.MultiplyNT(1.0,modcmat2,Projection,1.0);
  // contribution: 2/3*Tr(J^(-2/3)modstress) (Cinv \odot Cinv - 1/3 Cinv \otimes Cinv)
  modcmat.Clear();
  modcmat.MultiplyNT(-1.0/3.0,icg,icg);
  AddtoCmatHolzapfelProduct(modcmat, icg, 1.0);
  LINALG::Matrix<1,1> tracemat;
  tracemat.MultiplyTN(2./3.*std::pow(modinv(2),-2./3.),modstress,rcg);
  cmatisomodiso.Update(tracemat(0,0),modcmat,1.0);
  //contribution: -2/3 (Cinv \otimes S_iso + S_iso \otimes Cinv)
  cmatisomodiso.MultiplyNT(-2./3.,icg,stressisomodiso,1.0);
  cmatisomodiso.MultiplyNT(-2./3.,stressisomodiso,icg,1.0);

  //volumetric contribution
  //contribution: 2 \tilde p Cinv \otimes Cinv
  cmatisomodvol.MultiplyNT(modinv(2)* moddelta(4),icg,icg,1.0);
  //contribution: -2 J*p Cinv \odot Cinv
  AddtoCmatHolzapfelProduct(cmatisomodvol, icg, -2*modinv(2)*modgamma(2));

  return ;
}

/*----------------------------------------------------------------------*/
/* Calculates the stress and constitutive tensor for viscos part  09/13 */
/*----------------------------------------------------------------------*/
void MAT::ViscoElastHyper::EvaluateIsoModVisco(
    LINALG::Matrix<6,1>& stressisomodisovisco,
    LINALG::Matrix<6,1>& stressisomodvolvisco,
    LINALG::Matrix<6,6>& cmatisomodisovisco,
    LINALG::Matrix<6,6>& cmatisomodvolvisco,
    LINALG::Matrix<3,1> prinv,
    LINALG::Matrix<3,1> modinv,
    LINALG::Matrix<8,1> modmy,
    LINALG::Matrix<33,1> modxi,
    LINALG::Matrix<6,1> rcg,
    LINALG::Matrix<6,1> id2,
    LINALG::Matrix<6,1> icg,
    LINALG::Matrix<6,6> id4,
    LINALG::Matrix<6,1> modrcgrate
    )
{
	// define necessary variables
	const double modscale = std::pow(prinv(2),-1./3.);

	// 2nd Piola Kirchhoff stresses

	// isochoric contribution
	LINALG::Matrix<6,1> modstress(true);
	modstress.Update(modmy(1), id2);
	modstress.Update(modmy(2), modrcgrate, 1.0);
	// build 4-tensor for projection as 6x6 tensor
	LINALG::Matrix<6,6> Projection;
	Projection.MultiplyNT(1./3., icg, rcg);
	Projection.Update(1.0, id4, -1.0);
	// isochoric stress
	stressisomodisovisco.MultiplyNN(modscale,Projection,modstress,1.0);

	// volumetric contribution:
	// with visco_isoratedep: no volumetric part added --> always 0


	// Constitutive Tensor

	//isochoric contribution
	// modified constitutive tensor
	LINALG::Matrix<6,6> modcmat(true);
	LINALG::Matrix<6,6> modcmat2(true);
    // contribution:  Id \otimes \overline{\dot{C}} + \overline{\dot{C}} \otimes Id
	modcmat.MultiplyNT(modxi(1), id2, modrcgrate);
	modcmat.MultiplyNT(modxi(1), modrcgrate, id2, 1.0);
	// contribution: Id4
	modcmat.Update(modxi(2), id4, 1.0);
    //scaling
	modcmat.Scale(std::pow(modinv(2),-4./3.));
	//contribution: P:\overline{C}:P
	modcmat2.MultiplyNN(Projection,modcmat);
	cmatisomodisovisco.MultiplyNT(1.0,modcmat2,Projection,1.0);
	// contribution: 2/3*Tr(J^(-2/3)modstress) (Cinv \odot Cinv - 1/3 Cinv \otimes Cinv)
	modcmat.Clear();
	modcmat.MultiplyNT(-1.0/3.0,icg,icg);
	AddtoCmatHolzapfelProduct(modcmat, icg, 1.0);
	LINALG::Matrix<1,1> tracemat;
	tracemat.MultiplyTN(2./3.*std::pow(modinv(2),-2./3.),modstress,rcg);
	cmatisomodisovisco.Update(tracemat(0,0),modcmat,1.0);
	//contribution: -2/3 (Cinv \otimes S_iso^v + S_iso^v \otimes Cinv)
	cmatisomodisovisco.MultiplyNT(-2./3.,icg,stressisomodisovisco,1.0);
	cmatisomodisovisco.MultiplyNT(-2./3.,stressisomodisovisco,icg,1.0);

	// volumetric contribution:
	// with visco_isoratedep: no volumetric part added --> always 0

	return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ViscoElastHyper::EvaluateAnisotropicPrinc(
    LINALG::Matrix<6,1>& stressanisoprinc,
    LINALG::Matrix<6,6>& cmatanisoprinc,
    LINALG::Matrix<6,1> rcg,
    Teuchos::ParameterList& params
    )
{
  // loop map of associated potential summands
  for (unsigned int p=0; p<potsum_.size(); ++p)
  {
    potsum_[p]->AddStressAnisoPrincipal(rcg,cmatanisoprinc,stressanisoprinc,params);
  }

  return ;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ViscoElastHyper::EvaluateAnisotropicMod(
    LINALG::Matrix<6,1>& stressanisomod,
    LINALG::Matrix<6,6>& cmatanisomod,
    LINALG::Matrix<6,1> rcg,
    LINALG::Matrix<6,1> icg,
    LINALG::Matrix<3,1> prinv
    )
{

  // loop map of associated potential summands
  for (unsigned int p=0; p<potsum_.size(); ++p)
  {
    potsum_[p]->AddStressAnisoModified(rcg,icg,cmatanisomod,stressanisomod,prinv(2));
  }

  return ;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ViscoElastHyper::ResponseStretches(
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
void MAT::ViscoElastHyper::VisNames(std::map<std::string,int>& names)
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
bool MAT::ViscoElastHyper::VisData(const std::string& name, std::vector<double>& data, int numgp, int eleID)
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
