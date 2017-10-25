/*----------------------------------------------------------------------*/
/*!
\file plastic_VarConstUpdate.cpp
\brief This file contains the hyperelastic toolbox with application to finite
strain plasticity using a variational constitutive update.

The input line should read
MAT 1 MAT_plasticElastHyperVCU NUMMAT 1 MATIDS 2 DENS 1.0 INITYIELD 0.45
ISOHARD 0.12924 EXPISOHARD 16.93 INFYIELD 0.715 KINHARD 0.0

\level 3

<pre>
\maintainer Alexander Seitz
            seitz@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15271
</pre>
 */

/*----------------------------------------------------------------------*/

#include "plastic_VarConstUpdate.H"
#include "../drt_lib/standardtypes_cpp.H"
#include "../drt_matelast/elast_summand.H"
#include "../linalg/linalg_utils.H"
#include "../drt_lib/drt_linedefinition.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_mat/material_service.H"
#include "Epetra_SerialDenseSolver.h"
#include "../linalg/linalg_fixedsizematrix.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::PAR::PlasticElastHyperVCU::PlasticElastHyperVCU(
    Teuchos::RCP<MAT::PAR::Material> matdata
)
: MAT::PAR::PlasticElastHyper(matdata)
{
  // polyconvexity check is just implemented for isotropic hyperlastic materials
  if (polyconvex_)
    dserror("This polyconvexity-check is just implemented for isotropic "
        "hyperelastic-materials (do not use for plastic materials).");
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<MAT::Material> MAT::PAR::PlasticElastHyperVCU::CreateMaterial()
{
  return Teuchos::rcp(new MAT::PlasticElastHyperVCU(this));
}


MAT::PlasticElastHyperVCUType MAT::PlasticElastHyperVCUType::instance_;


DRT::ParObject* MAT::PlasticElastHyperVCUType::Create( const std::vector<char> & data )
{
  MAT::PlasticElastHyperVCU* elhy = new MAT::PlasticElastHyperVCU();
  elhy->Unpack(data);

  return elhy;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::PlasticElastHyperVCU::PlasticElastHyperVCU()
: params_(NULL)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::PlasticElastHyperVCU::PlasticElastHyperVCU(MAT::PAR::PlasticElastHyperVCU* params)
: params_(params)
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
void MAT::PlasticElastHyperVCU::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm( data );
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);
  // matid
  int matid = -1;
  if (MatParams() != NULL) matid = MatParams()->Id();  // in case we are in post-process mode
  AddtoPack(data,matid);
  AddtoPack(data,(int)isoprinc_);
  AddtoPack(data,(int)isomod_);
  AddtoPack(data,(int)anisoprinc_);
  AddtoPack(data,(int)anisomod_);

  if (MatParams() != NULL) // summands are not accessible in postprocessing mode
  {
    // loop map of associated potential summands
    for (unsigned int p=0; p<potsum_.size(); ++p)
    {
     potsum_[p]->PackSummand(data);
    }
  }

  // plastic history data
  AddtoPack<3,3>(data,last_plastic_defgrd_inverse_);
  AddtoPack(data,last_alpha_isotropic_);

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::PlasticElastHyperVCU::Unpack(const std::vector<char>& data)
{
  // make sure we have a pristine material
  params_ = NULL;
  potsum_.clear();

  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");

  // matid and recover MatParams()
  int matid;
  ExtractfromPack(position,data,matid);
  if (DRT::Problem::Instance()->Materials() != Teuchos::null)
  {
    if (DRT::Problem::Instance()->Materials()->Num() != 0)
    {
      const unsigned int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();
      MAT::PAR::Parameter* mat = DRT::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        params_ = static_cast<MAT::PAR::PlasticElastHyperVCU*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(), MaterialType());
    }
  }

  isoprinc_=(bool)ExtractInt(position,data);
  isomod_=(bool)ExtractInt(position,data);
  anisoprinc_=(bool)ExtractInt(position,data);
  anisomod_=(bool)ExtractInt(position,data);

  if (MatParams() != NULL) // summands are not accessible in postprocessing mode
  {
    // make sure the referenced materials in material list have quick access parameters
    std::vector<int>::const_iterator m;
    for (m=MatParams()->matids_->begin(); m!=MatParams()->matids_->end(); ++m)
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
  }

  // plastic history data
  ExtractfromPack<3,3>(position,data,last_plastic_defgrd_inverse_);
  ExtractfromPack(position,data,last_alpha_isotropic_);

  // no need to pack this
  delta_alpha_i_ .resize(last_alpha_isotropic_.size(),0.);
  plastic_defgrd_inverse_.resize(last_plastic_defgrd_inverse_.size());

  // in the postprocessing mode, we do not unpack everything we have packed
  // -> position check cannot be done in this case
  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d",data.size(),position);

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::PlasticElastHyperVCU::Setup(int numgp, DRT::INPUT::LineDefinition* linedef)
{
  // setup the plasticelasthyper data
  PlasticElastHyper::Setup(numgp,linedef);

  // setup history
  plastic_defgrd_inverse_.resize(numgp);

  return;
}

// MAIN
void MAT::PlasticElastHyperVCU::Evaluate(const LINALG::Matrix<3,3>* defgrd,  ///< Deformation gradient
    const LINALG::Matrix<6,1>* glstrain,///< Green-Lagrange strain
    Teuchos::ParameterList& params,     ///< Container for additional information
    LINALG::Matrix<6,1>* stress,        ///< 2nd Piola-Kirchhoff stresses
    LINALG::Matrix<6,6>* cmat,          ///< Constitutive matrix
    const int eleGID)                    ///< Element GID
{
  int gp = params.get<int>("gp");

  double last_ai = last_alpha_isotropic_[gp];
  LINALG::Matrix<3,3> empty;

  // Get cetrial
  LINALG::Matrix<3,3> id2;
  for (int i=0;i<3;i++) id2(i,i)=1.0;

  LINALG::Matrix<3,3> cetrial;
  LINALG::Matrix<6,1> ee_test;
  CompElastQuant(defgrd,last_plastic_defgrd_inverse_[gp],id2,&cetrial,&ee_test);

  // get 2pk stresses
  LINALG::Matrix<6,1> etstr;
  LINALG::Matrix<6,6> etcmat;
  ElastHyper::Evaluate(NULL,&ee_test,params,&etstr,&etcmat,eleGID);

  double yf;
  double normZero =0.0;

  LINALG::Matrix<3,3> mandelStr;
  LINALG::Matrix<3,3> devMandelStr;
  YieldFunction(last_ai,normZero,id2,cetrial,etstr,&yf,devMandelStr,mandelStr);

  if (yf <= 0)
  {
    // step is elastic
    stress->Clear();
    cmat->Clear();

    LINALG::Matrix<6,1> checkStr;
    LINALG::Matrix<6,6> checkCmat;
    LINALG::Matrix<3,3> emptymat;
    PlasticElastHyper::EvaluateElast(defgrd,&emptymat,stress,cmat,gp,eleGID);
    ElastHyper::Evaluate(defgrd,&ee_test,params,&checkStr,&checkCmat,eleGID);

    // push back
    LINALG::Matrix<3,3> checkStrMat;
    for (int i=0;i<3;i++) checkStrMat(i,i) = checkStr(i);
    checkStrMat(0,1) = checkStrMat(1,0) = checkStr(3);
    checkStrMat(1,2) = checkStrMat(2,1) = checkStr(4);
    checkStrMat(0,2) = checkStrMat(2,0) = checkStr(5);

    LINALG::Matrix<3,3> tmp33;
    LINALG::Matrix<3,3> strWithPlast;
    tmp33.Multiply(last_plastic_defgrd_inverse_[gp],checkStrMat);
    strWithPlast.MultiplyNT(tmp33,last_plastic_defgrd_inverse_[gp]);

    plastic_defgrd_inverse_[gp] = last_plastic_defgrd_inverse_[gp];
    delta_alpha_i_[gp] = 0.;
  }

  else
  {
    int iterator = 0;
    const int maxiter = 200;

    LINALG::Matrix<3,3> devMandelStr_direction(devMandelStr);
    devMandelStr_direction.Scale(1./devMandelStr.Norm2());

    LINALG::Matrix<5,5> iH;

    double plastMulti = 1.e-8;
    LINALG::Matrix<3,3> dLpStart;
    dLpStart.Update(plastMulti,devMandelStr_direction);

    LINALG::Matrix<5,1> beta;
    beta(0) = dLpStart(0,0);
    beta(1) = dLpStart(1,1);
    beta(2) = dLpStart(0,1);
    beta(3) = dLpStart(1,2);
    beta(4) = dLpStart(0,2);

    LINALG::Matrix<3,3> dLp;
    bool converged = false;

    LINALG::Matrix<5,5> hess;
    LINALG::Matrix<5,1> rhs;
    LINALG::Matrix<6,6> dcedlp;
    LINALG::Matrix<9,6> dFpiDdeltaDp;
    LINALG::Matrix<6,5> P;
    P(0,0) = 1.;
    P(1,1) = 1.;
    P(2,0) = -1.;
    P(2,1) = -1.;
    P(3,2) = 1.;
    P(4,3) = 1.;
    P(5,4) = 1.;

    do
    {
      iterator++;

      dLp(0,0) = beta(0);
      dLp(1,1) = beta(1);
      dLp(2,2) = -beta(0)-beta(1);
      dLp(0,1) = dLp(1,0) = beta(2);
      dLp(2,1) = dLp(1,2) = beta(3);
      dLp(0,2) = dLp(2,0) = beta(4);

      LINALG::Matrix<5,1> rhsElast;
      LINALG::Matrix<6,1> eeOut;
      EvaluateRHS(gp,dLp,*defgrd,eeOut,rhs,rhsElast,dcedlp,dFpiDdeltaDp,params,eleGID);

      // Hessian matrix elastic component
      LINALG::Matrix<6,1> elastStress;
      LINALG::Matrix<6,6> elastCmat;
      LINALG::Matrix<6,1> elastStressDummy;
      LINALG::Matrix<6,6> elastCmatDummy;
      ElastHyper::Evaluate(NULL,&eeOut,params,&elastStress,&elastCmat,eleGID);

      LINALG::Matrix<6,6> d2ced2lpVoigt[6];
      Ce2ndDeriv(defgrd,last_plastic_defgrd_inverse_[gp],dLp,d2ced2lpVoigt);

      LINALG::Matrix<6,6> cpart_tmp;
      cpart_tmp.Multiply(elastCmat,dcedlp);
      LINALG::Matrix<6,6> cpart;
      cpart.MultiplyTN(.25,dcedlp,cpart_tmp);

      LINALG::Matrix<6,6> spart6x6;
      for (int A=0;A<6;++A)
        for (int B=0;B<6;++B)
          for (int C=0;C<6;++C)
            if (A<3)
              spart6x6(B,C) += 0.5* elastStress(A)* d2ced2lpVoigt[A](B,C);
            else
              spart6x6(B,C) += 1.*elastStress(A)* d2ced2lpVoigt[A](B,C);

      LINALG::Matrix<6,6> hessElast6x6(spart6x6);
      hessElast6x6.Update(1.,cpart,1.);

      LINALG::Matrix<6,5> tmp65;
      tmp65.Multiply(hessElast6x6,P);
      LINALG::Matrix<5,5> hessElast;
      hessElast.MultiplyTN(P,tmp65);

      LINALG::Matrix<5,1> dlp_vec;
      dlp_vec(0)=2.*beta(0)+beta(1);
      dlp_vec(1)=2.*beta(1)+beta(0);
      dlp_vec(2)=2.*beta(2);
      dlp_vec(3)=2.*beta(3);
      dlp_vec(4)=2.*beta(4);

      // dissipative component
      LINALG::Matrix<5,5> hess_a;
      for (int i=0;i<5;++i)hess_a(i,i)=2.;
      hess_a(0,1)=1.;
      hess_a(1,0)=1.;
      LINALG::Matrix<5,5> hess_aiso(hess_a);
      hess_a.Scale(1./dLp.Norm2());
      LINALG::Matrix<5,5> tmpSummandIdentity(hess_a);
      double hess_aisoScalar = Isohard();
      hess_aisoScalar *= last_alpha_isotropic_[gp] / dLp.Norm2();
      hess_aisoScalar += Isohard();
      hess_aiso.Scale(hess_aisoScalar);

      hess_a.MultiplyNT((-1.) /(dLp.Norm2()*dLp.Norm2()*dLp.Norm2()),dlp_vec,dlp_vec,1.);
      LINALG::Matrix<5,5> tmpSummandDlpVec;
      tmpSummandDlpVec.MultiplyNT((-1.) /(dLp.Norm2()*dLp.Norm2()*dLp.Norm2()),dlp_vec,dlp_vec);
      tmpSummandDlpVec.Scale(sqrt(2./3.)*Inityield());
      tmpSummandIdentity.Scale(sqrt(2./3.)*Inityield());
      hess_a.Scale(sqrt(2./3.)*Inityield());

      LINALG::Matrix<5,5> tmp55;
      tmp55.MultiplyNT(dlp_vec,dlp_vec);
      tmp55.Scale((sqrt(2./3.)*last_alpha_isotropic_[gp]*Isohard())/(dLp.Norm2()*dLp.Norm2()*dLp.Norm2()));

      hess_aiso.Update(-1.,tmp55,1.);
      hess_aiso.Scale(sqrt(2./3.));

      // Hessian matrix for nonlinear iso hardening
      LINALG::Matrix<5,5> hessIsoNL;
      for (int i=0;i<5;++i)hessIsoNL(i,i)=2.;
      hessIsoNL(0,1)=1.;
      hessIsoNL(1,0)=1.;
      double hessIsoNLscalar= Isohard();
      double new_ai = last_alpha_isotropic_[gp] + sqrt(2./3.)*dLp.Norm2();
      double k = Infyield() - Inityield();
      hessIsoNLscalar *= new_ai;
      hessIsoNLscalar += k;
      hessIsoNLscalar -= k*std::exp(-1.*Expisohard()*new_ai);
      hessIsoNL.Scale(hessIsoNLscalar/dLp.Norm2());
      LINALG::Matrix<5,5> tmpNL55;
      tmpNL55.MultiplyNT(dlp_vec,dlp_vec);
      LINALG::Matrix<5,5> tmp2(tmpNL55);
      tmpNL55.Scale(hessIsoNLscalar/(dLp.Norm2()*dLp.Norm2()*dLp.Norm2()));
      hessIsoNL.Update(-1.,tmpNL55,1.);
      double tmp2scalar=Isohard();
      tmp2scalar += Expisohard()*k*std::exp(-1.*Expisohard()*new_ai);
      tmp2.Scale((sqrt(2./3.)*tmp2scalar)/(dLp.Norm2()*dLp.Norm2()));
      hessIsoNL.Update(1.,tmp2,1.);
      hessIsoNL.Scale(sqrt(2./3.));

      // Compose the hessian matrix
      LINALG::Matrix<5,5> hess_analyt(hessElast);
      hess_analyt.Update(1.,hess_a,1.);
      hess_analyt.Update(1.,hessIsoNL,1.); // For nonlinear isotropic hardening

      if (rhs.Norm2()<1.0e-12)
        converged=true;

      iH = hess_analyt;
      LINALG::FixedSizeSerialDenseSolver<5,5,1> solver;
      solver.SetMatrix(iH);
      if (solver.Invert() != 0) dserror("Inversion failed");

      LINALG::Matrix<5,1> beta_incr;
      beta_incr.Multiply(-1.,iH,rhs,0.);

       beta.Update(1.,beta_incr,1.);

//      // Backtracking Line Search Method
//      double damping_fac=1.;
//
//      LINALG::Matrix<5,1> beta_test(beta);
//      beta_test.Update(damping_fac,beta_incr,1.);
//
//      LINALG::Matrix<3,3> dLp_test;
//      dLp_test(0,0) = beta_test(0);
//      dLp_test(1,1) = beta_test(1);
//      dLp_test(2,2) = -beta_test(0)-beta_test(1);
//      dLp_test(0,1) = dLp_test(1,0) = beta_test(2);
//      dLp_test(2,1) = dLp_test(1,2) = beta_test(3);
//      dLp_test(0,2) = dLp_test(2,0) = beta_test(4);
//
//      LINALG::Matrix<5,1> rhs_test;
//      LINALG::Matrix<5,1> rhsElast_test;
//
//      LINALG::Matrix<6,6> dummytest;
//      LINALG::Matrix<6,1> dummyeetest;
//      EvaluateRHS(gp,dLp_test,*defgrd,dummyeetest,rhs_test,rhsElast_test,dummytest,params,eleGID);
//
//      double criteria= rhs_test.Norm2() - rhs.Norm2();
//
//      if (criteria <= 0.)
//        beta = beta_test;
//      else
//      {
//        double new_criteria = 1.0;
//        int dampIter = 0;
//        int maxDampIter = 20;
//        LINALG::Matrix<5,1> beta_corr(beta);
//        while (( new_criteria > 0. )&& (dampIter<maxDampIter))
//        {
//          damping_fac *= 0.5;
//          beta_corr=beta;
//          beta_corr.Update(damping_fac,beta_incr,1.);
//
//          LINALG::Matrix<3,3> dLp_corr;
//          dLp_corr(0,0) = beta_corr(0);
//          dLp_corr(1,1) = beta_corr(1);
//          dLp_corr(2,2) = -beta_corr(0)-beta_corr(1);
//          dLp_corr(0,1) = dLp_corr(1,0) = beta_corr(2);
//          dLp_corr(2,1) = dLp_corr(1,2) = beta_corr(3);
//          dLp_corr(0,2) = dLp_corr(2,0) = beta_corr(4);
//
//          LINALG::Matrix<5,1> rhs_corr;
//          LINALG::Matrix<5,1> rhsElast_corr;
//          LINALG::Matrix<6,6> dummy;
//          LINALG::Matrix<6,1> dummy1;
//          EvaluateRHS(gp,dLp_corr,*defgrd,dummy1,rhs_corr,rhsElast_corr,dummy,params,eleGID);
//
//          new_criteria = rhs_corr.Norm2() - rhs.Norm2();
//
//          dampIter++;
//        }
//        if (dampIter==maxDampIter)
//          beta.Update(1.,beta_incr,1.);
//        else
//          beta = beta_corr;
//      }
    } while(!converged && (iterator < maxiter));

    if (!converged )
    {
      std::cout << "eleGID: " << eleGID << "gp: " << gp << std::endl;
      dserror("unconverged");
    }

    dLp(0,0) = beta(0);
    dLp(1,1) = beta(1);
    dLp(2,2) = -beta(0)-beta(1);
    dLp(0,1) = dLp(1,0) = beta(2);
    dLp(2,1) = dLp(1,2) = beta(3);
    dLp(0,2) = dLp(2,0) = beta(4);

    // Get exp, Dexp and DDexp
    LINALG::Matrix<3,3> input_dLp(dLp);
    input_dLp.Scale(-1.);
    LINALG::Matrix<3,3> expOut(input_dLp);
    LINALG::Matrix<6,6> dexpOut_mat;
    MatrixExponentialDerivativeSym3x3(input_dLp,dexpOut_mat);
    MatrixExponential3x3(expOut);

    plastic_defgrd_inverse_[gp].Multiply(last_plastic_defgrd_inverse_[gp],expOut);
    delta_alpha_i_[gp] = sqrt(2./3.)*dLp.Norm2();

    //Compute the total stresses
    LINALG::Matrix<6,6> tangent_elast;
    PlasticElastHyper::EvaluateElast(defgrd,&dLp,stress,&tangent_elast,gp,eleGID);

    LINALG::Matrix<6,9> dPK2dFpinvIsoprinc;
    Dpk2dFpi(eleGID,defgrd,&plastic_defgrd_inverse_[gp],dPK2dFpinvIsoprinc);

    LINALG::Matrix<6,6> mixedDeriv;
    mixedDeriv.Multiply(dPK2dFpinvIsoprinc,dFpiDdeltaDp);


    LINALG::Matrix<6,5> tmp2;
    tmp2.Multiply(mixedDeriv,P);

    LINALG::Matrix<6,5> mixDerivInvHess;
    mixDerivInvHess.Multiply(tmp2,iH);

    LINALG::Matrix<6,6> cmat_summand1;
    cmat_summand1.Update(1.,tangent_elast);
    LINALG::Matrix<6,6> cmat_summand2;
    cmat_summand2.MultiplyNT(-1.,mixDerivInvHess,tmp2);

    cmat->Update(cmat_summand1,cmat_summand2);
  }

  return;
}


/// update after converged time step
void MAT::PlasticElastHyperVCU::Update()
{
  // update local history data F_n <-- F_{n+1}
  for (unsigned gp=0;gp<last_plastic_defgrd_inverse_.size();++gp)
  {
    // for a real plastic update, update the plastic history here
    last_plastic_defgrd_inverse_[gp] = plastic_defgrd_inverse_[gp];
    last_alpha_isotropic_[gp] += delta_alpha_i_[gp];
  }

  return;
};




// Evaluate dCedlp
void MAT::PlasticElastHyperVCU::EvalDceDlp(
    const LINALG::Matrix<3,3> fpi,
    const LINALG::Matrix<3,3>* defgrd,
    const LINALG::Matrix<6,6> Dexp,
    const LINALG::Matrix<3,3> cetrial,
    const LINALG::Matrix<3,3> explp,
    LINALG::Matrix<6,6>& dceDdeltalp,
    LINALG::Matrix<9,6>& dFpiDdeltaDp
)
{
  // compute dcedlp this way: dcedlp = dcedfpi : dfpidlp
  // first compute dcedfpi
  LINALG::Matrix<3,3> id2;
  for (int i=0;i<3;i++) id2(i,i)+=1.0;
  LINALG::Matrix<3,3> next_fpi;
  next_fpi.Multiply(fpi,explp);
  LINALG::Matrix<3,3> rcg;
  rcg.MultiplyTN(*defgrd,*defgrd);
  LINALG::Matrix<3,3> inputB;
  inputB.Multiply(next_fpi,rcg);

  LINALG::Matrix<6,9> dcedfpi;

  LINALG::Matrix<3,3> tmp;
  tmp.MultiplyTN(next_fpi,rcg);
  AddRightNonSymmetricHolzapfelProduct(dcedfpi,tmp,id2,1.0);

  // Derivative of inverse plastic deformation gradient
  dFpiDdeltaDp.Clear();
  for (int A=0; A<3; A++)
    for (int a=0; a<3; a++)
      for (int b=0; b<3; b++)
        for (int i=0; i<6; i++)
          if (i<=2)
            dFpiDdeltaDp(VOIGT3X3NONSYM_[A][a],i) -= fpi(A,b)*Dexp(VOIGT3X3_[b][a],i);
          else
            dFpiDdeltaDp(VOIGT3X3NONSYM_[A][a],i) -= 2.*fpi(A,b)*Dexp(VOIGT3X3_[b][a],i); //fixme factor 2

  dceDdeltalp.Multiply(dcedfpi,dFpiDdeltaDp);

  for (int i=3;i<6;++i)
    for (int j=0;j<6;++j)
      dceDdeltalp(i,j)*=2.;
}


// MAP 5x1 on 9x1
void MAT::PlasticElastHyperVCU::YieldFunction(
    const double last_ai,
    const double norm_dLp,    // this is zero when yf<0 is checked
    const LINALG::Matrix<3,3> ExpEqui, // this is identity when yf < 0 is checked
    const LINALG::Matrix<3,3> cetr,
    const LINALG::Matrix<6,1> str,
    double* yieldFunc,
    LINALG::Matrix<3,3>& devMandelStr,
    LINALG::Matrix<3,3>& MandelStr)
{
  double Qi=0.0;
  double new_ai = last_ai +norm_dLp*sqrt(2./3.);  //fixme sqrt
  double k = Infyield() - Inityield();
  Qi -= Isohard()*new_ai;
  Qi -= k;
  Qi += k* std::exp(-1.*Expisohard()*new_ai);
  Qi *= sqrt(2./3.);

  double Qeq = sqrt(2./3.) * Inityield();

  LINALG::Matrix<3,3> ce;
  LINALG::Matrix<3,3> tmp;
  tmp.Multiply(cetr,ExpEqui);
  ce.Multiply(ExpEqui,tmp);

  // Compute the deviator of mandelStr, then its norm
  LINALG::Matrix<3,3> devMandelStrSumm2;

  // se strainlike
  LINALG::Matrix<6,1> se_strain(str);
  for (int i=3;i<6;i++)  se_strain(i)*=2.;

  for (int i=0;i<3;i++)
    for (int j=0;j<3;j++)
      for (int k=0;k<3;k++)
        MandelStr(i,k) += ce(i,j) * str(VOIGT3X3_[j][k]);

  // Compute the trace of the mandel stresses
  LINALG::Matrix<3,3> id2;
  for (int i=0;i<3;i++) id2(i,i)=1.0;
  double trMandelStr=0.0;
  for (int i=0;i<3;i++) trMandelStr += MandelStr(i,i);
  devMandelStrSumm2.Update(-trMandelStr/3.0,id2);
  devMandelStr.Update(MandelStr,devMandelStrSumm2);

  double NormDevMandelStr = devMandelStr.Norm2();

  // Compose the yield function
  double yf = 0.0;
  yf = NormDevMandelStr + Qi - Qeq;
  *yieldFunc = yf;

  return;
}


// Get elastic quantities Cetrial and Ee_n+1
void MAT::PlasticElastHyperVCU::CompElastQuant(
    const LINALG::Matrix<3,3>* defgrd,
    const LINALG::Matrix<3,3> fpi,
    const LINALG::Matrix<3,3> MatExp,
    LINALG::Matrix<3,3>* cetrial,
    LINALG::Matrix<6,1>* Ee
)

{
  LINALG::Matrix<3,3> fetrial;
  fetrial.Multiply(*defgrd,fpi);
  LINALG::Matrix<3,3> cetr;
  cetr.MultiplyTN(fetrial,fetrial);
  *cetrial = cetr;

  LINALG::Matrix<3,3> tmp;
  tmp.Multiply(MatExp,cetr);
  LINALG::Matrix<3,3> next_ce;
  next_ce.Multiply(tmp,MatExp);

  // Compute Ee_n+1 (first in 3x3 then map to 6x1 VOIGT)
  LINALG::Matrix<3,3> id2int;
  for (int i=0;i<3;i++) id2int(i,i)=1.0;
  LINALG::Matrix<3,3> next_ee3x3;
  next_ee3x3.Update(0.5,next_ce,-0.5,id2int);
  LINALG::Matrix<6,1> next_ee;
  for (int i=0;i<3;i++) next_ee(i) = next_ee3x3(i,i);
  next_ee(3) = 2.*next_ee3x3(0,1);
  next_ee(4) = 2.*next_ee3x3(1,2);
  next_ee(5) = 2.*next_ee3x3(0,2);

  *Ee = next_ee;

  return;
}

/*---------------------------------------------------------------------*
 | return names of visualization data (public)                         |
 *---------------------------------------------------------------------*/
void MAT::PlasticElastHyperVCU::VisNames(std::map<std::string,int>& names)
{
  std::string accumulatedstrain = "accumulatedstrain";
  names[accumulatedstrain] = 1; // scalar

}  // VisNames()


/*---------------------------------------------------------------------*
 | return visualization data (public)                                  |
 *---------------------------------------------------------------------*/
bool MAT::PlasticElastHyperVCU::VisData(
    const std::string& name,
    std::vector<double>& data,
    int numgp,
    int eleID
)
{
  if (name == "accumulatedstrain")
  {
    if ((int)data.size()!=1) dserror("size mismatch");
    double tmp=0.;
    for (unsigned gp=0; gp<last_alpha_isotropic_.size(); gp++)
      tmp+= AccumulatedStrain(gp);
    data[0] = tmp/last_alpha_isotropic_.size();
  }
  return false;

}  // VisData()


// 2nd matrix exponential derivatives with 6 parameters
void MAT::PlasticElastHyperVCU::MatrixExponentialSecondDerivativeSym3x3x6(
    const LINALG::Matrix<3,3> MatrixIn,
    LINALG::Matrix<3,3>& exp,
    LINALG::Matrix<6,6>& dexp_mat,
    LINALG::Matrix<6,6>* MatrixExp2ndDerivVoigt
)
{
  LINALG::Matrix<3,3> MatrixExp1stDeriv[6];
  LINALG::Matrix<3,3> MatrixExp2ndDeriv[6][6];

  // temporary matrices
  LINALG::Matrix<3,3> akm(true);
  LINALG::Matrix<3,3> ak(true);
  LINALG::Matrix<3,3> akmd[6];
  LINALG::Matrix<3,3> akd[6];
  LINALG::Matrix<3,3> akmdd[6][6];
  LINALG::Matrix<3,3> akdd[6][6];

  // derivatives of A w.r.t. beta's
  LINALG::Matrix<3,3> da[6];
  da[0](0,0)=1.;
  da[1](1,1)=1.;
  da[2](2,2)=1.;
  da[3](0,1) = da[3](1,0) = 0.5;
  da[4](1,2) = da[4](2,1) = 0.5;
  da[5](0,2) = da[5](2,0) = 0.5;

  // prepare
  exp.Clear();

  // start with first entry
  int k=0;
  int kmax=200;
  for (int i=0; i<3; i++)
    exp(i,i) = 1.;

  // increment
  ++k;
  akm=exp;
  ak.Multiply(1./(double)k,akm,MatrixIn);
  akd[0](0,0)=1.;
  akd[1](1,1)=1.;
  akd[2](2,2)=1.;
  akd[3](0,1) = akd[3](1,0) = 0.5;
  akd[4](1,2) = akd[4](2,1) = 0.5;
  akd[5](0,2) = akd[5](2,0) = 0.5;

  do
  {
    // add summand
    exp.Update(1.,ak,1.);

    // increment
    ++k;
    akm.Update(ak);
    ak.Multiply(1./(double)k,akm,MatrixIn);

    // 1st derivatives
    for (int i=0; i<6; i++)
    {
      MatrixExp1stDeriv[i].Update(1.,akd[i],1.);
      // increment
      akmd[i].Update(akd[i]);
      akd[i].Multiply(1./(double)k,akm,da[i]);
      akd[i].Multiply(1./(double)k,akmd[i],MatrixIn,1.);
    }

    // 2nd derivatives
    for (int i=0; i<6; i++)
      for (int j=i; j<6; j++)
        MatrixExp2ndDeriv[i][j].Update(1.,akdd[i][j],1.);

    for (int i=0; i<6; i++)
      for (int j=i; j<6; j++)
      {
        // increment
        akmdd[i][j] = akdd[i][j];
        akdd[i][j].Multiply(1./((double)k),akmd[i],da[j]);
        akdd[i][j].Multiply(1./((double)k),akmd[j],da[i],1.);
        akdd[i][j].Multiply(1./((double)k),akmdd[i][j],MatrixIn,1.);
      }

  } while (k<kmax && ak.Norm2()>1.e-16);

  if (k==kmax)
  {
    std::cout << "matrixIn: " << MatrixIn ;
    dserror("Matrix exponential unconverged with %i summands", k) ;
  }

  // Zusatz: 1. Map MatExpFirstDer from [6](3,3) to (6,6)
  for (int i=0;i<6;i++)
    for (int j=0;j<3;j++)
      for (int k=j;k<3;k++)
        dexp_mat(i,VOIGT3X3NONSYM_[j][k]) += MatrixExp1stDeriv[i](j,k);


  // 2. Map MatExp2ndDeriv from [6][6]3x3 to [6]6x6
  for (int I=0;I<6;I++)
    for (int J=I;J<6;J++)
      for (int k=0;k<3;k++)
        for (int l=k;l<3;l++)
        {
          MatrixExp2ndDerivVoigt[I](J,VOIGT3X3NONSYM_[k][l]) = MatrixExp2ndDeriv[I][J](k,l);
          MatrixExp2ndDerivVoigt[J](I,VOIGT3X3NONSYM_[k][l]) = MatrixExp2ndDeriv[I][J](k,l);
        }

  return;
}


/*---------------------------------------------------------------------------*
 | Calculate second derivative of matrix exponential via polar decomposition |
 | following Ortiz et.al. 2001                                   seitz 02/14 |
 *---------------------------------------------------------------------------*/
void MAT::PlasticElastHyperVCU::MatrixExponentialSecondDerivativeSym3x3(
    const LINALG::Matrix<3,3> MatrixIn,
    LINALG::Matrix<3,3>& exp,
    std::vector<LINALG::Matrix<3,3> >& MatrixExp1stDeriv,
    std::vector<std::vector<LINALG::Matrix<3,3> > >& MatrixExp2ndDeriv
)
{
  // temporary matrices
  const LINALG::Matrix<3,3> zeros(true);
  LINALG::Matrix<3,3> akm(true);
  LINALG::Matrix<3,3> ak(true);
  std::vector<LINALG::Matrix<3,3> > akmd(5,zeros);
  std::vector<LINALG::Matrix<3,3> > akd(5,zeros);
  std::vector<std::vector<LINALG::Matrix<3,3> > >akmdd;
  std::vector<std::vector<LINALG::Matrix<3,3> > >akdd;

  // derivatives of A w.r.t. beta's
  std::vector<LINALG::Matrix<3,3> > da(5,ak);
  da[0](0,0) = 1.; da[0](2,2) = -1.;
  da[1](1,1) = 1.; da[1](2,2) = -1.;
  da[2](0,1) = 1.; da[2](1,0) = 1.;
  da[3](1,2) = 1.; da[3](2,1) = 1.;
  da[4](0,2) = 1.; da[4](2,0) = 1.;

  // prepare
  exp.Clear();
  MatrixExp1stDeriv.resize(5,zeros);
  MatrixExp2ndDeriv.resize(5);
  akmdd.resize(5);
  akdd.resize(5);
  for (int i=0; i<5; i++)
  {
    MatrixExp2ndDeriv[i].resize(5,zeros);
    akmdd[i].resize(5,zeros);
    akdd[i].resize(5,zeros);
  }

  // start with first entry
  int k=0;
  int kmax=200;
  for (int i=0; i<3; i++)
    exp(i,i) = 1.;

  // increment
  ++k;
  akm=exp;
  ak.Multiply(1./(double)k,akm,MatrixIn);
  akd=da;

  do
  {
    // add summand
    exp.Update(1.,ak,1.);

    // increment
    ++k;
    akm.Update(ak);
    ak.Multiply(1./(double)k,akm,MatrixIn);

    // 1st derivatives
    for (int i=0; i<5; i++)
    {
      MatrixExp1stDeriv[i].Update(1.,akd[i],1.);
      // increment
      akmd[i].Update(akd[i]);
      akd[i].Multiply(1./(double)k,akm,da[i]);
      akd[i].Multiply(1./(double)k,akmd[i],MatrixIn,1.);
    }

    // 2nd derivatives
    for (int i=0; i<5; i++)
      for (int j=0; j<5; j++)
        MatrixExp2ndDeriv[i][j].Update(1.,akdd[i][j],1.);

    for (int i=0; i<5; i++)
      for (int j=0; j<5; j++)
      {
        // increment
        akmdd[i][j] = akdd[i][j];
        akdd[i][j].Multiply(1./((double)k),akmd[i],da[j]);
        akdd[i][j].Multiply(1./((double)k),akmd[j],da[i],1.);
        akdd[i][j].Multiply(1./((double)k),akmdd[i][j],MatrixIn,1.);
      }

  } while (k<kmax && ak.Norm2()>1.e-16);

  if (k==kmax)
  {
    std::cout << "matrixIn: " << MatrixIn ;
    dserror("Matrix exponential unconverged with %i summands", k) ;
  }
  return;
}


void MAT::PlasticElastHyperVCU::EvaluateRHS(
    const int gp,
    const LINALG::Matrix<3,3> dLp,
    const LINALG::Matrix<3,3> defgrd,
    LINALG::Matrix<6,1>& eeOut,
    LINALG::Matrix<5,1>& rhs,
    LINALG::Matrix<5,1>& rhsElast,
    LINALG::Matrix<6,6>& dcedlp,
    LINALG::Matrix<9,6>& dFpiDdeltaDp,
    Teuchos::ParameterList& params,
    const int eleGID
)
{
  LINALG::Matrix<3,3> zeros;
  LINALG::Matrix<6,6> zeros66;
  // set zero
  rhs.Clear();

  // Get exp, Dexp and DDexp
  LINALG::Matrix<3,3> dLpIn(dLp);
  dLpIn.Scale(-1.);
  LINALG::Matrix<3,3> expOut(dLpIn);
  MatrixExponential3x3(expOut);
  LINALG::Matrix<6,6> dexpOut_mat;
  MatrixExponentialDerivativeSym3x3(dLpIn,dexpOut_mat);

  LINALG::Matrix<3,3> fpi_incr(dLp);
  fpi_incr.Scale(-1.);
  LINALG::Matrix<6,6> derivExpMinusLP;
  MatrixExponentialDerivativeSym3x3(fpi_incr,derivExpMinusLP);
  MatrixExponential3x3(fpi_incr);

  LINALG::Matrix<3,3> fetrial;
  fetrial.Multiply(defgrd,last_plastic_defgrd_inverse_[gp]);
  LINALG::Matrix<3,3> cetrial;
  cetrial.MultiplyTN(fetrial,fetrial);

  LINALG::Matrix<3,3> fpi;
  fpi.Multiply(last_plastic_defgrd_inverse_[gp],expOut);
  LINALG::Matrix<3,3> fe;
  fe.Multiply(defgrd,fpi);
  LINALG::Matrix<3,3> ce;
  ce.MultiplyTN(fe,fe);

  for (int i=0;i<3;++i)
    eeOut(i)=0.5*(ce(i,i)-1.);
  eeOut(3)=ce(0,1);
  eeOut(4)=ce(1,2);
  eeOut(5)=ce(0,2);

  LINALG::Matrix<6,1> se;
  LINALG::Matrix<6,6> dummy;
  ElastHyper::Evaluate(NULL,&eeOut,params,&se,&dummy,eleGID);

  EvalDceDlp(last_plastic_defgrd_inverse_[gp],&defgrd,dexpOut_mat,cetrial,expOut,dcedlp,dFpiDdeltaDp);

  LINALG::Matrix<6,1> rhs6;

  LINALG::Matrix<6,1> dLp_vec;
  for (int i=0;i<3;++i) dLp_vec(i) = dLp(i,i);
  dLp_vec(3)=2.*dLp(0,1);
  dLp_vec(4)=2.*dLp(1,2);
  dLp_vec(5)=2.*dLp(0,2);


  double new_ai = last_alpha_isotropic_[gp] + sqrt(2./3.)*dLp.Norm2();
  double k = Infyield() - Inityield();
  double rhsPlastScalar = Isohard();
  rhsPlastScalar *= new_ai;
  rhsPlastScalar += k;
  rhsPlastScalar -= k*std::exp(-1.*Expisohard()*new_ai);
  rhs6.Update((rhsPlastScalar*sqrt(2./3.))/dLp.Norm2(),dLp_vec,1.);  // plastic component, nonlinear iso hardening


  rhs6.Update((Inityield()*sqrt(2./3.))/dLp.Norm2(),dLp_vec,1.);    // dissipative component

  rhs6.MultiplyTN(0.5,dcedlp,se,1.);          //elastic component
  LINALG::Matrix<6,1> rhs6Elast(rhs6);   //rhs6Elast is the elastic part of the right hand side. We need this to build the hessian numerically


  LINALG::Matrix<6,5> dAlphadBeta;
  dAlphadBeta(0,0)=1.;
  dAlphadBeta(1,1)=1.;
  dAlphadBeta(2,0)=-1.;
  dAlphadBeta(2,1)=-1.;
  dAlphadBeta(3,2)=1.;
  dAlphadBeta(4,3)=1.;
  dAlphadBeta(5,4)=1.;

  // Outputs
  rhs.MultiplyTN(1.,dAlphadBeta,rhs6,0.);
  rhsElast.MultiplyTN(1.,dAlphadBeta,rhs6Elast,0.);

  return;
}
void MAT::PlasticElastHyperVCU::EvaluatePlast(
    LINALG::Matrix<6,9>& dPK2dFpinvIsoprinc,
    const LINALG::Matrix<3,1>& gamma,
    const LINALG::Matrix<8,1>& delta,
    const LINALG::Matrix<3,3>& id2,
    const LINALG::Matrix<6,1>& Cpi,
    const LINALG::Matrix<3,3>& Fpi,
    const LINALG::Matrix<3,3>& CpiC,
    const LINALG::Matrix<9,1>& CFpi,
    const LINALG::Matrix<9,1>& CFpiCei,
    const LINALG::Matrix<6,1>& ircg,
    const LINALG::Matrix<3,3>& FpiCe,
    const LINALG::Matrix<9,1>& CFpiCe,
    const LINALG::Matrix<6,1>& CpiCCpi
)
{

  // derivative of PK2 w.r.t. inverse plastic deformation gradient
  AddRightNonSymmetricHolzapfelProduct(dPK2dFpinvIsoprinc,id2,Fpi,gamma(0));
  AddRightNonSymmetricHolzapfelProduct(dPK2dFpinvIsoprinc,CpiC,Fpi,gamma(1));
  dPK2dFpinvIsoprinc.MultiplyNT(delta(0),Cpi,CFpi,1.);
  dPK2dFpinvIsoprinc.MultiplyNT(delta(1),Cpi,CFpiCe,1.);
  dPK2dFpinvIsoprinc.MultiplyNT(delta(1),CpiCCpi,CFpi,1.);
  dPK2dFpinvIsoprinc.MultiplyNT(delta(2),Cpi,CFpiCei,1.);
  dPK2dFpinvIsoprinc.MultiplyNT(delta(2),ircg,CFpi,1.);
  dPK2dFpinvIsoprinc.MultiplyNT(delta(3),CpiCCpi,CFpiCe,1.);
  dPK2dFpinvIsoprinc.MultiplyNT(delta(4),CpiCCpi,CFpiCei,1.);
  dPK2dFpinvIsoprinc.MultiplyNT(delta(4),ircg,CFpiCe,1.);
  dPK2dFpinvIsoprinc.MultiplyNT(delta(5),ircg,CFpiCei,1.);
  AddRightNonSymmetricHolzapfelProduct(dPK2dFpinvIsoprinc,id2,FpiCe,0.5*delta(7));
}

void MAT::PlasticElastHyperVCU::EvaluateKinQuantPlast(
    const int eleGID,
    const LINALG::Matrix<3,3>* defgrd,
    const LINALG::Matrix<3,3>* fpi,
    LINALG::Matrix<3,1>& gamma,
    LINALG::Matrix<8,1>& delta,
    LINALG::Matrix<3,3>& id2,
    LINALG::Matrix<6,1>& Cpi,
    LINALG::Matrix<3,3>& CpiC,
    LINALG::Matrix<9,1>& CFpi,
    LINALG::Matrix<9,1>& CFpiCei,
    LINALG::Matrix<6,1>& ircg,
    LINALG::Matrix<3,3>& FpiCe,
    LINALG::Matrix<9,1>& CFpiCe,
    LINALG::Matrix<6,1>& CpiCCpi
)
{

  id2.Clear();
  for (int i=0;i<3;i++)id2(i,i) =1.;

  LINALG::Matrix<3,3> fe;
  fe.Multiply(*defgrd,*fpi);
  LINALG::Matrix<3,3> ce3x3;
  ce3x3.MultiplyTN(fe,fe);

  // ce here strainlike
  LINALG::Matrix<6,1> ce;
  for (int i=0;i<3;i++) ce(i) = ce3x3(i,i);
  ce(3) = 2.*ce3x3(0,1);
  ce(4) = 2.*ce3x3(1,2);
  ce(5) = 2.*ce3x3(0,2);

  LINALG::Matrix<6,1> ce_stresslike;
  for (int i=0;i<6;i++)
    if (i<3)
      ce_stresslike(i) = ce(i);
    else
      ce_stresslike(i) = 0.5*ce(i);

  LINALG::Matrix<3,1> prinv;
  InvariantsPrincipal(prinv,ce);

  LINALG::Matrix<3,1> dPI;
  LINALG::Matrix<6,1> ddPII;
  EvaluateInvariantDerivatives(prinv,dPI,ddPII,eleGID,potsum_);
  CalculateGammaDelta(gamma,delta,prinv,dPI,ddPII);

  // inverse plastic right Cauchy-Green
  LINALG::Matrix<3,3> CpiM;
  CpiM.MultiplyNT(*fpi,*fpi);
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
  LINALG::Matrix<3,3> tmp33;
  LINALG::Matrix<3,3> CpiCCpiM;
  tmp33.Multiply(CpiM,RCG);
  CpiCCpiM.Multiply(tmp33,CpiM);
  // stress-like Voigt notation
  for (int i=0; i<3; i++) CpiCCpi(i) = CpiCCpiM(i,i);
  CpiCCpi(3) = (CpiCCpiM(0,1)+CpiCCpiM(1,0))/2.;
  CpiCCpi(4) = (CpiCCpiM(2,1)+CpiCCpiM(1,2))/2.;
  CpiCCpi(5) = (CpiCCpiM(0,2)+CpiCCpiM(2,0))/2.;

  CpiC.Multiply(CpiM,RCG);
  FpiCe.Multiply(*fpi,ce3x3);

  //    FpiTC.MultiplyTN(invpldefgrd,RCG);
  //    CeFpiTC.Multiply(CeM,FpiTC);
  LINALG::Matrix<3,3> tmp;
  tmp.Multiply(RCG,*fpi);
  Matrix3x3to9x1(tmp,CFpi);
  tmp33.Multiply(tmp,ce3x3);
  Matrix3x3to9x1(tmp33,CFpiCe);



  tmp.Invert(ce3x3);
  tmp33.Multiply(*fpi,tmp);
  tmp.Multiply(RCG,tmp33);
  Matrix3x3to9x1(tmp,CFpiCei);

}
void MAT::PlasticElastHyperVCU::Dpk2dFpi(
    const int eleGID,
    const LINALG::Matrix<3,3>* defgrd,
    const LINALG::Matrix<3,3>* fpi,
    LINALG::Matrix<6,9>& dPK2dFpinvIsoprinc
)
{

  LINALG::Matrix<3,1> gamma;
  LINALG::Matrix<8,1> delta;
  LINALG::Matrix<3,3> id2;
  LINALG::Matrix<6,1> Cpi;
  LINALG::Matrix<3,3> CpiC;
  LINALG::Matrix<9,1> CFpi;
  LINALG::Matrix<9,1> CFpiCei;
  LINALG::Matrix<6,1> ircg;
  LINALG::Matrix<3,3> FpiCe;
  LINALG::Matrix<9,1> CFpiCe;
  LINALG::Matrix<6,1> CpiCCpi;
  EvaluateKinQuantPlast(eleGID,defgrd,fpi,gamma,delta,id2,Cpi,CpiC,CFpi,CFpiCei,ircg,FpiCe,CFpiCe,CpiCCpi);

  EvaluatePlast(dPK2dFpinvIsoprinc,gamma,delta,id2,Cpi,*fpi,CpiC,CFpi,CFpiCei,ircg,FpiCe,CFpiCe,CpiCCpi);
}

// Compute DpsiplastDalphaiso
void MAT::PlasticElastHyperVCU::DpsiplastDalphaiso(
    const double norm_dLp,
    const double last_alphaiso,
    const double isoHardMod,
    const double initYield,
    const double infYield,
    const double expIsoHard,
    double* dpsiplastdalphaiso
)
{
  double new_alphaiso=0.0;
  new_alphaiso = last_alphaiso + norm_dLp;
  double summand1=0.0;
  summand1 = isoHardMod*new_alphaiso;
  double summand2=0.0;
  double k= infYield - initYield;
  summand2 = k;
  double summand3=0.0;
  summand3 = k*std::exp(-expIsoHard*new_alphaiso);


  double out=0.0;
  out = summand1 + summand2 - summand3;
  *dpsiplastdalphaiso = out;
}



// Compute DDceDdLpDdLp
void MAT::PlasticElastHyperVCU::Ce2ndDeriv(
    const LINALG::Matrix<3,3>* defgrd,
    const LINALG::Matrix<3,3> fpi,
    const LINALG::Matrix<3,3> dLp,
    LINALG::Matrix<6,6>* DDceDdLpDdLpVoigt
)
{
  // Compute ce trial
  LINALG::Matrix<3,3> fetrial;
  fetrial.Multiply(*defgrd,fpi);
  LINALG::Matrix<3,3> cetrial;
  cetrial.MultiplyTN(fetrial,fetrial);

  // Get matrix exponential derivatives
  LINALG::Matrix<3,3> minus_dLp(dLp);
  minus_dLp.Scale(-1.);
  LINALG::Matrix<3,3> zeros;
  LINALG::Matrix<3,3> exp_dLp;
  LINALG::Matrix<6,6> Dexp_dLp_mat;
  LINALG::Matrix<6,6> D2exp_VOIGT[6];
  MatrixExponentialSecondDerivativeSym3x3x6(minus_dLp,exp_dLp,Dexp_dLp_mat,D2exp_VOIGT);

  LINALG::Matrix<3,3> exp_dLp_cetrial;
  exp_dLp_cetrial.Multiply(exp_dLp,cetrial);

  for (int a=0;a<3;a++)
    for (int d=a;d<3;d++)
      for (int b=0;b<3;b++)
        for (int C=0;C<6;C++)
          for (int D=0;D<6;D++)
          {
            DDceDdLpDdLpVoigt[VOIGT3X3NONSYM_[a][d]](C,D)+= (1.+(D>2))*(1.+(C>2))*(
                exp_dLp_cetrial(a,b)*D2exp_VOIGT[VOIGT3X3_[b][d]](C,D)
                +D2exp_VOIGT[VOIGT3X3_[a][b]](C,D) * exp_dLp_cetrial(d,b));
            for (int c=0;c<3;c++)
              DDceDdLpDdLpVoigt[VOIGT3X3NONSYM_[a][d]](C,D)+= (1.+(C>2))*(1.+(D>2))*(
                  Dexp_dLp_mat(VOIGT3X3_[a][b],C)*cetrial(b,c) *Dexp_dLp_mat(VOIGT3X3_[c][d],D)
                  +Dexp_dLp_mat(VOIGT3X3_[a][b],D)*cetrial(b,c)*Dexp_dLp_mat(VOIGT3X3_[c][d],C));
          }

}
