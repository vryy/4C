/*----------------------------------------------------------------------*/
/*!
\file elasthyper.cpp

<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15237
</pre>
*/

/*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include <vector>
#include "elasthyper.H"
#include "../drt_matelast/elast_summand.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::PAR::ElastHyper::ElastHyper(
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

  // make sure the referenced materials in material list have quick access parameters
  std::vector<int>::const_iterator m;
  for (m=matids_->begin(); m!=matids_->end(); ++m)
  {
    const int matid = *m;
    Teuchos::RCP<MAT::ELAST::Summand> potsum = MAT::ELAST::Summand::Factory(matid);
    if (potsum == Teuchos::null) dserror("Failed to allocate");
    pot_.insert(std::pair<int,Teuchos::RCP<MAT::ELAST::Summand> >(matid,potsum));
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::ElastHyper::ElastHyper()
  : params_(NULL)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::ElastHyper::ElastHyper(MAT::PAR::ElastHyper* params)
  : params_(params)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ElastHyper::Pack(vector<char>& data) const
{
  data.resize(0);

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);
  // matid
  int matid = -1;
  if (params_ != NULL) matid = params_->Id();  // in case we are in post-process mode
  AddtoPack(data,matid);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ElastHyper::Unpack(const std::vector<char>& data)
{
  int position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");

  // matid and recover params_
  int matid;
  ExtractfromPack(position,data,matid);
  // in post-process mode we do not have any instance of DRT::Problem
  if (DRT::Problem::NumInstances() > 0)
  {
    const int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();
    MAT::PAR::Parameter* mat = DRT::Problem::Instance(probinst)->Materials()->ParameterById(matid);
    if (mat->Type() == MaterialType())
      params_ = static_cast<MAT::PAR::ElastHyper*>(mat);
    else
      dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(), MaterialType());
  }
  else
  {
    params_ = NULL;
  }

  if (position != (int)data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ElastHyper::InvariantsPrincipal(
  LINALG::Matrix<3,1>& prinv,
  const LINALG::Matrix<6,1>& rcg
  )
{
  // 1st invariant, trace
  prinv(0) = rcg(0) + rcg(1) + rcg(2);
  // 2nd invariant
  prinv(1) = 0.5*( prinv(0)*prinv(0)
                   - rcg(0)*rcg(0) - rcg(1)*rcg(1) - rcg(2)*rcg(2)
                   - 2.0*rcg(3)*rcg(3) - 2.0*rcg(4)*rcg(4) - 2.0*rcg(5)*rcg(5) );
  // 3rd invariant, determinant
  prinv(2) = rcg(0)*rcg(1)*rcg(2)
    + 0.25 * rcg(3)*rcg(4)*rcg(5)
    - 0.25 * rcg(1)*rcg(5)*rcg(5)
    - 0.25 * rcg(2)*rcg(3)*rcg(3)
    - 0.25 * rcg(0)*rcg(4)*rcg(4);

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ElastHyper::Evaluate(
  const LINALG::Matrix<6,1>& glstrain,
  LINALG::Matrix<6,6>& cmat,
  LINALG::Matrix<6,1>& stress
  )
{
  // build Cartesian identity 2-tensor I
  LINALG::Matrix<6,1> id2(true);
  for (int i=0; i<3; i++) id2(i) = 1.0;

  // right Cauchy-Green Tensor  C = 2 * E + I
  // REMARK: strain-like 6-Voigt vector
  LINALG::Matrix<6,1> rcg(glstrain);
  rcg.Scale(2.0);
  rcg.Update(1.0, id2, 1.0);

  // principal invariants of right Cauchy-Green strain
  LINALG::Matrix<3,1> prinv;
  InvariantsPrincipal(prinv, rcg);

  // invert right Cauchy-Green tensor
  // REMARK: stress-like 6-Voigt vector
  LINALG::Matrix<6,1> icg(false);
  {
    icg(0) = ( rcg(1)*rcg(2) - 0.25*rcg(4)*rcg(4) ) / prinv(2);
    icg(1) = ( rcg(0)*rcg(2) - 0.25*rcg(5)*rcg(5) ) / prinv(2);
    icg(2) = ( rcg(0)*rcg(1) - 0.25*rcg(3)*rcg(3) ) / prinv(2);
    icg(3) = ( 0.25*rcg(5)*rcg(4) - 0.5*rcg(3)*rcg(2) ) / prinv(2);
    icg(4) = ( 0.25*rcg(3)*rcg(5) - 0.5*rcg(0)*rcg(4) ) / prinv(2);
    icg(5) = ( 0.25*rcg(3)*rcg(4) - 0.5*rcg(5)*rcg(1) ) / prinv(2);
  }

  // principal coefficients
  bool havecoeffprinc = false;
  LINALG::Matrix<3,1> gamma(true);
  LINALG::Matrix<8,1> delta(true);

  // loop map of associated potentials
  std::map<int,Teuchos::RCP<MAT::ELAST::Summand> >& pot = params_->pot_;
  std::map<int,Teuchos::RCP<MAT::ELAST::Summand> >::iterator p;
  for (p=pot.begin(); p!=pot.end(); ++p)
  {
    p->second->AddCoefficientsPrincipal(havecoeffprinc,gamma,delta,prinv);
  }

  // we make right Cauchy-Green 6-Voigt vector _stress_-like
  // needed for dyadic products underneath
  for (int i=3; i<6; i++) rcg(i) *= 0.5;

  // set Cartesian identity 4-tensor in 6x6-Voigt matrix notation
  LINALG::Matrix<6,6> id4(true);
  for (int i=0; i<3; i++) id4(i,i) = 1.0;
  for (int i=3; i<6; i++) id4(i,i) = 0.5;

  // blank resulting quantities
  // ... even if it is an implicit law that cmat is zero upon input
  stress.Clear();
  cmat.Clear();

  // build stress response and elasticity tensor
  // for potentials based on principal invariants
  if (havecoeffprinc)
  {
    // 2nd Piola Kirchhoff stresses
    stress.Update(gamma(0), id2, 1.0);
    stress.Update(gamma(1), rcg, 1.0);
    stress.Update(gamma(2), icg, 1.0);

    // constitutive tensor
    // contribution: Id \otimes Id
    cmat.MultiplyNT(delta(0), id2, id2, 1.0);
    // contribution: Id \otimes C + C \otimes Id
    cmat.MultiplyNT(delta(1), id2, rcg, 1.0);
    cmat.MultiplyNT(delta(1), rcg, id2, 1.0);
    // contribution: Id \otimes Cinv + Cinv \otimes Id
    cmat.MultiplyNT(delta(2), id2, icg, 1.0);
    cmat.MultiplyNT(delta(2), icg, id2, 1.0);
    // contribution: C \otimes C
    cmat.MultiplyNT(delta(3), rcg, rcg, 1.0);
    // contribution: C \otimes Cinv + Cinv \otimes C
    cmat.MultiplyNT(delta(4), rcg, icg, 1.0);
    cmat.MultiplyNT(delta(4), icg, rcg, 1.0);
    // contribution: Cinv \otimes Cinv
    cmat.MultiplyNT(delta(5), icg, icg, 1.0);
    // contribution: Cinv \odot Cinv
    AddtoCmatHolzapfelProduct(cmat, icg, delta(6));
    // contribution: Id4
    cmat.Update(delta(7), id4, 1.0);
  }
  
  return;
}
#endif
