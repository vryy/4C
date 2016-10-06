/*----------------------------------------------------------------------*/
/*!
\file elast_remodelfiber.cpp

\brief Base class for fiber materials which remodel. Has a pointer to the single fiber families
the input line should read
MAT 11 ELAST_RemodelFiber NUMMAT 2 MATIDS 111 112 TDECAY 101.0 GROWTHFAC 4.951051289713897e-04 COLMASSFRAC 0.062 0.248 DEPOSITIONSTRETCH 1.062

\level 3

\maintainer Fabian Braeu

*----------------------------------------------------------------------*/
/* headers */
#include "elast_remodelfiber.H"
#include "../drt_mat/matpar_material.H"
#include "../drt_lib/standardtypes_cpp.H"
#include "../drt_lib/drt_linedefinition.H"
#include "../drt_mat/material_service.H"
#include "../linalg/linalg_utils.H"
#include "Epetra_SerialDenseSolver.h"



/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
MAT::ELASTIC::PAR::RemodelFiber::RemodelFiber(
    Teuchos::RCP<MAT::PAR::Material> matdata
)
: Parameter(matdata),
  nummat_(matdata->GetInt("NUMMAT")),
  matids_(matdata->Get<std::vector<int> >("MATIDS")),
  tdecay_(matdata->GetDouble("TDECAY")),
  k_growth_(matdata->GetDouble("GROWTHFAC")),
  init_w_col_(matdata->GetMutable<std::vector<double> >("COLMASSFRAC")),
  G_(matdata->GetDouble("DEPOSITIONSTRETCH"))
{
  // check if sizes fit
  if (nummat_ != (int)matids_->size())
    dserror("number of materials %d does not fit to size of material vector %d", nummat_, matids_->size());

  // check decay time validity
  if (tdecay_<=0.)
    dserror("decay time must be positive");
}


/*----------------------------------------------------------------------*
 |  Constructor                             (public)   fb         09/15 |
 *----------------------------------------------------------------------*/
MAT::ELASTIC::RemodelFiber::RemodelFiber(MAT::ELASTIC::PAR::RemodelFiber* params)
  : params_(params),
    potsumfiberpas_(0),
    potsumfiberact_(0)
{
  // make sure the referenced materials in material list have quick access parameters
  std::vector<int>::const_iterator m;
  for (m=params_->matids_->begin(); m!=params_->matids_->end(); ++m)
  {
    const int matid = *m;
    Teuchos::RCP<MAT::ELASTIC::Summand> sum = MAT::ELASTIC::Summand::Factory(matid);
    if (sum == Teuchos::null)
      dserror("Failed to allocate");

    if (sum->MaterialType() == INPAR::MAT::mes_coupanisoexpoactive)
      potsumfiberact_.push_back(Teuchos::rcp_static_cast<MAT::ELASTIC::CoupAnisoExpoActive>(sum));
    else
      potsumfiberpas_.push_back(sum);
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::RemodelFiber::PackSummand(DRT::PackBuffer& data) const
{
  int num_fiber = 0;
  num_fiber = last_lambda_r_.size();

  AddtoPack(data,num_fiber);

  for(int i=0;i<num_fiber;++i)
  {
    AddtoPack(data,last_lambda_r_[i]);
    AddtoPack(data,cur_lambda_r_[i]);
    AddtoPack(data,cur_rho_col_[i]);
    AddtoPack(data,last_rho_col_[i]);
    AddtoPack(data,init_rho_col_[i]);
    AddtoPack(data,stress_[i]);
    AddtoPack(data,sigmapre_[i]);
  }

  AddtoPack(data,Av_);
  AddtoPack(data,A_strain_);
  AddtoPack(data,AM_);
  AddtoPack(data,AM_orth_);
  AddtoPack(data,A9x1_);

  if (params_ != NULL) // summands are not accessible in postprocessing mode
  {
    // loop map of associated potential summands
    for (unsigned int k=0; k<potsumfiberpas_.size(); ++k)
      potsumfiberpas_[k]->PackSummand(data);

    for (unsigned int k=0; k<potsumfiberact_.size(); ++k)
      potsumfiberact_[k]->PackSummand(data);
  }

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::RemodelFiber::UnpackSummand(const std::vector<char>& data,
                                               std::vector<char>::size_type& position)
{
  int num_fiber=0;
  ExtractfromPack(position,data,num_fiber);

  last_lambda_r_.resize(num_fiber);
  cur_lambda_r_.resize(num_fiber);
  cur_rho_col_.resize(num_fiber);
  last_rho_col_.resize(num_fiber);
  init_rho_col_.resize(num_fiber);
  stress_.resize(num_fiber);
  sigmapre_.resize(num_fiber);

  for(int i=0;i<num_fiber;++i)
  {
    ExtractfromPack(position,data,last_lambda_r_[i]);
    ExtractfromPack(position,data,cur_lambda_r_[i]);
    ExtractfromPack(position,data,cur_rho_col_[i]);
    ExtractfromPack(position,data,last_rho_col_[i]);
    ExtractfromPack(position,data,init_rho_col_[i]);
    ExtractfromPack(position,data,stress_[i]);
    ExtractfromPack(position,data,sigmapre_[i]);
  }

  ExtractfromPack(position,data,Av_);
  ExtractfromPack(position,data,A_strain_);
  ExtractfromPack(position,data,AM_);
  ExtractfromPack(position,data,AM_orth_);
  ExtractfromPack(position,data,A9x1_);

  // loop map of associated potential summands
  for (unsigned int k=0; k<potsumfiberpas_.size(); ++k)
    potsumfiberpas_[k]->UnpackSummand(data,position);

  for (unsigned int k=0; k<potsumfiberact_.size(); ++k)
    potsumfiberact_[k]->UnpackSummand(data,position);

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::RemodelFiber::Setup(int numgp,double rho_tot,DRT::INPUT::LineDefinition* linedef)
{
  // setup fiber and inelastic history variable
  last_lambda_r_.resize(potsumfiberpas_.size()+potsumfiberact_.size());
  cur_lambda_r_.resize(potsumfiberpas_.size()+potsumfiberact_.size());
  cur_rho_col_.resize(potsumfiberpas_.size()+potsumfiberact_.size());
  last_rho_col_.resize(potsumfiberpas_.size()+potsumfiberact_.size());
  stress_.resize(potsumfiberpas_.size()+potsumfiberact_.size());
  init_rho_col_.resize(potsumfiberpas_.size()+potsumfiberact_.size());
  sigmapre_.resize(potsumfiberpas_.size()+potsumfiberact_.size());


  // some variables
  LINALG::Matrix<2,1> dPI(true);
  LINALG::Matrix<3,1> ddPII(true);
  LINALG::Matrix<4,1> dddPIII(true);
  LINALG::Matrix<6,1> stressactive(true);
  LINALG::Matrix<6,6> cmatactive(true);
  double stress_f_act = 0.0;

  for(unsigned k=0;k<potsumfiberpas_.size();++k)
  {
    init_rho_col_[k].resize(numgp,rho_tot * params_->init_w_col_->at(k));
    last_lambda_r_[k].resize(numgp,1.0);
    cur_lambda_r_[k].resize(numgp,1.0);
    cur_rho_col_[k].resize(numgp,init_rho_col_[k][0]);
    last_rho_col_[k].resize(numgp,init_rho_col_[k][0]);
    stress_[k].resize(numgp,1.0);

    potsumfiberpas_[k]->Setup(linedef);
  }
  for(unsigned k=potsumfiberpas_.size();k<(potsumfiberpas_.size()+potsumfiberact_.size());++k)
  {
    init_rho_col_[k].resize(numgp,rho_tot * params_->init_w_col_->at(k));
    last_lambda_r_[k].resize(numgp,1.0);
    cur_lambda_r_[k].resize(numgp,1.0);
    cur_rho_col_[k].resize(numgp,init_rho_col_[k][0]);
    last_rho_col_[k].resize(numgp,init_rho_col_[k][0]);
    stress_[k].resize(numgp,1.0);

    potsumfiberact_[k-potsumfiberpas_.size()]->Setup(linedef);
  }


  SetupStructuralTensorsGR();

  // quadratic prestretch in tensor notation
  LINALG::Matrix<3,3> GM(true);

  // identity tenosr
  LINALG::Matrix<3,3> id(true);
  for(int i=0;i<3;++i) id(i,i) = 1.0;

  for(unsigned k=0;k<potsumfiberpas_.size();++k)
  {
    GM.Update(params_->G_*params_->G_,AM_[k],0.0);

    // Cauchy prestress of new mass which is deposited during G&R (assumption: det(F_pre)=1)
    potsumfiberpas_[k]->GetDerivativesAniso(dPI,ddPII,dddPIII,GM,0);
    sigmapre_[k] = 2.0*dPI(0)*params_->G_*params_->G_;
  }
  for(unsigned k=potsumfiberpas_.size();k<(potsumfiberpas_.size()+potsumfiberact_.size());++k)
  {
    GM.Update(params_->G_*params_->G_,AM_[k],0.0);

    // Cauchy prestress of new mass which is deposited during G&R (assumption: det(F_pre)=1)
    potsumfiberact_[k-potsumfiberpas_.size()]->GetDerivativesAniso(dPI,ddPII,dddPIII,GM,0);
    sigmapre_[k] = 2.0*dPI(0)*params_->G_*params_->G_;

    // active fiber contribution
    potsumfiberact_[k-potsumfiberpas_.size()]->EvaluateActiveStressCmatAniso(id,cmatactive,stressactive,0);
    stress_f_act = stressactive.Dot(A_strain_[k]);
    sigmapre_[k] += stress_f_act;
  }

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::RemodelFiber::SetupStructuralTensorsGR()
{
  // identity tensor
  LINALG::Matrix<3,3> id(true);
  for(int i=0;i<3;++i)
    id(i,i) = 1.0;

  // fiber directions
  std::vector<LINALG::Matrix<3,1> > fibervecs;

  Av_.resize(potsumfiberpas_.size()+potsumfiberact_.size(),LINALG::Matrix<6,1>(true));
  A_strain_.resize(potsumfiberpas_.size()+potsumfiberact_.size(),LINALG::Matrix<6,1>(true));
  AM_.resize(potsumfiberpas_.size()+potsumfiberact_.size(),LINALG::Matrix<3,3>(true));
  AM_orth_.resize(potsumfiberpas_.size()+potsumfiberact_.size(),LINALG::Matrix<3,3>(true));
  A9x1_.resize(potsumfiberpas_.size()+potsumfiberact_.size(),LINALG::Matrix<9,1>(true));

  for(unsigned k=0;k<potsumfiberpas_.size();++k)
  {
    // Get fiberdirection
    potsumfiberpas_[k]->GetFiberVecs(fibervecs);

    for (int i = 0; i < 3; ++i)
      A_strain_[k](i) = Av_[k](i) = fibervecs[k](i)*fibervecs[k](i);
    Av_[k](3) = fibervecs[k](0)*fibervecs[k](1);
    Av_[k](4) = fibervecs[k](1)*fibervecs[k](2);
    Av_[k](5) = fibervecs[k](0)*fibervecs[k](2);

    A_strain_[k](3) = 2.0*Av_[k](3);
    A_strain_[k](4) = 2.0*Av_[k](4);
    A_strain_[k](5) = 2.0*Av_[k](5);

    // build structural tensor in matrix notation
    for(int i=0;i<3;++i)
      AM_[k](i,i) = Av_[k](i);
    AM_[k](0,1) = Av_[k](3);
    AM_[k](1,2) = Av_[k](4);
    AM_[k](0,2) = Av_[k](5);
    AM_[k](1,0) = Av_[k](3);
    AM_[k](2,1) = Av_[k](4);
    AM_[k](2,0) = Av_[k](5);

    // orthogonal structural tensor ( 1_{ij} - A_{ij} )
    AM_orth_[k].Update(1.0,AM_[k],0.0);
    AM_orth_[k].Update(1.0,id,-1.0);

    // build structural tensor in 9x1 vector-like notation
    for(int i=0;i<3;++i)
      A9x1_[k](i) = A_strain_[k](i);
    A9x1_[k](3)=0.5*A_strain_[k](3);
    A9x1_[k](4)=0.5*A_strain_[k](4);
    A9x1_[k](5)=0.5*A_strain_[k](5);
    A9x1_[k](6)=0.5*A_strain_[k](3);
    A9x1_[k](7)=0.5*A_strain_[k](4);
    A9x1_[k](8)=0.5*A_strain_[k](5);
  }
  for(unsigned k=potsumfiberpas_.size();k<(potsumfiberpas_.size()+potsumfiberact_.size());++k)
  {
    // Get fiberdirection
    potsumfiberact_[k-potsumfiberpas_.size()]->GetFiberVecs(fibervecs);

    for (int i = 0; i < 3; ++i)
      A_strain_[k](i) = Av_[k](i) = fibervecs[k](i)*fibervecs[k](i);
    Av_[k](3) = fibervecs[k](0)*fibervecs[k](1);
    Av_[k](4) = fibervecs[k](1)*fibervecs[k](2);
    Av_[k](5) = fibervecs[k](0)*fibervecs[k](2);

    A_strain_[k](3) = 2.0*Av_[k](3);
    A_strain_[k](4) = 2.0*Av_[k](4);
    A_strain_[k](5) = 2.0*Av_[k](5);

    // build structural tensor in matrix notation
    for(int i=0;i<3;++i)
      AM_[k](i,i) = Av_[k](i);
    AM_[k](0,1) = Av_[k](3);
    AM_[k](1,2) = Av_[k](4);
    AM_[k](0,2) = Av_[k](5);
    AM_[k](1,0) = Av_[k](3);
    AM_[k](2,1) = Av_[k](4);
    AM_[k](2,0) = Av_[k](5);

    // orthogonal structural tensor ( 1_{ij} - A_{ij} )
    AM_orth_[k].Update(1.0,AM_[k],0.0);
    AM_orth_[k].Update(1.0,id,-1.0);

    // build structural tensor in 9x1 vector-like notation
    for(int i=0;i<3;++i)
      A9x1_[k](i) = A_strain_[k](i);
    A9x1_[k](3)=0.5*A_strain_[k](3);
    A9x1_[k](4)=0.5*A_strain_[k](4);
    A9x1_[k](5)=0.5*A_strain_[k](5);
    A9x1_[k](6)=0.5*A_strain_[k](3);
    A9x1_[k](7)=0.5*A_strain_[k](4);
    A9x1_[k](8)=0.5*A_strain_[k](5);
  }

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::RemodelFiber::Update()
{
  // update history variable
  for(unsigned k=0;k<(potsumfiberpas_.size()+potsumfiberact_.size());++k)
    for(unsigned gp=0;gp<cur_rho_col_[k].size();++gp)
    {
      last_lambda_r_[k][gp] = cur_lambda_r_[k][gp];
      last_rho_col_[k][gp] = cur_rho_col_[k][gp];
    }

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::RemodelFiber::EvaluateAnisotropicStressCmat(const LINALG::Matrix<3,3> CM,
                                                               const LINALG::Matrix<3,3> iFgM,
                                                               LINALG::Matrix<6,6>& cmat,
                                                               LINALG::Matrix<6,1>& stress,
                                                               const int gp,
                                                               const int eleGID)
{
  // clear some variables
  stress.Clear();
  cmat.Clear();

  // some variables
  LINALG::Matrix<3,3> tmp(true);
  LINALG::Matrix<3,3> iCin(true);
  LINALG::Matrix<2,1> dPIe(true);
  LINALG::Matrix<3,1> ddPIIe(true);
  LINALG::Matrix<4,1> dddPIIIe(true);
  LINALG::Matrix<6,6> dCedC(true);
  LINALG::Matrix<6,6> dCedC_strain(true);
  LINALG::Matrix<6,1> pk2ev(true);
  LINALG::Matrix<6,6> cmatelastic(true);


  // right elastic Cauchy Green tensor
  LINALG::Matrix<3,3> CeM(true);
  LINALG::Matrix<6,1> Cev(true);

  // temporary variables
  std::vector<LINALG::Matrix<3,1> > fibervecs;
  LINALG::Matrix<6,6> tmp6x6(true);

  // converts stress-like to strain-like Voigt notation// build structural tensor in matrix notation
  LINALG::Matrix<6,6> strainconv(true);
  for(int i=0;i<3;++i)
    strainconv(i,i) = 1.0;
  for(int j=3;j<6;++j)
    strainconv(j,j) = 2.0;

  LINALG::Matrix<3,3> FrM(true);
  LINALG::Matrix<3,3> iFrM(true);
  LINALG::Matrix<3,3> iFin(true);
  double ilambin_sq = 0.0;

  // additional variables for active fiber contribution (smooth muscle)
  double stress_f_act = 0.0;
  LINALG::Matrix<6,1> stressactive(true);
  LINALG::Matrix<6,6> cmatactive(true);

  // there is no growth and remodeling here
  for(unsigned k=0;k<potsumfiberpas_.size();++k)
  {
    // build remodel deformation gradient (prestretch included)
    FrM.Update(cur_lambda_r_[k][gp]*(1./params_->G_),AM_[k],0.0);
    FrM.Update(1./std::sqrt(cur_lambda_r_[k][gp]*1./params_->G_),AM_orth_[k],1.0);

    // inverse remodel deformation gradient (prestretch included)
    iFrM.Invert(FrM);

    // total inelastic deformation gradient
    iFin.MultiplyNN(1.0,iFgM,iFrM,0.0);

    // elastic right Cauchy-Green in matrix notation
    tmp.MultiplyNN(1.0,CM,iFin,0.0);
    CeM.MultiplyTN(1.0,iFin,tmp,0.0);

    // get derivatives of strain energy function w.r.t. I4
    potsumfiberpas_[k]->GetDerivativesAniso(dPIe,ddPIIe,dddPIIIe,CeM,eleGID);

    // Update stress
    iCin.MultiplyNT(1.0,iFin,iFin,0.0);
    ilambin_sq = iCin.Dot(AM_[k]);
    stress.Update(2.*cur_rho_col_[k][gp]*dPIe(0)*ilambin_sq,Av_[k],1.0);

    // update elasticity tensor
    cmat.MultiplyNT(4.*cur_rho_col_[k][gp]*ddPIIe(0)*ilambin_sq*ilambin_sq,Av_[k],Av_[k],1.0);


    // fiber stress for output
    potsumfiberpas_[k]->GetFiberVecs(fibervecs);
    stress_[k][gp] = CM.Dot(AM_[k])*2.0*dPIe(0)*ilambin_sq;
  }
  for(unsigned k=potsumfiberpas_.size();k<(potsumfiberpas_.size()+potsumfiberact_.size());++k)
  {
    // build remodel deformation gradient (prestretch included)
    FrM.Update(cur_lambda_r_[k][gp]*(1./params_->G_),AM_[k],0.0);
    FrM.Update(1./std::sqrt(cur_lambda_r_[k][gp]*1./params_->G_),AM_orth_[k],1.0);

    // inverse remodel deformation gradient (prestretch included)
    iFrM.Invert(FrM);

    // total inelastic deformation gradient
    iFin.MultiplyNN(1.0,iFgM,iFrM,0.0);

    // elastic right Cauchy-Green in matrix notation
    tmp.MultiplyNN(1.0,CM,iFin,0.0);
    CeM.MultiplyTN(1.0,iFin,tmp,0.0);

    // get derivatives of strain energy function w.r.t. I4
    potsumfiberact_[k-potsumfiberpas_.size()]->GetDerivativesAniso(dPIe,ddPIIe,dddPIIIe,CeM,eleGID);

    // Update stress
    iCin.MultiplyNT(1.0,iFin,iFin,0.0);
    ilambin_sq = iCin.Dot(AM_[k]);
    stress.Update(2.*cur_rho_col_[k][gp]*dPIe(0)*ilambin_sq,Av_[k],1.0);

    // update elasticity tensor
    cmat.MultiplyNT(4.*cur_rho_col_[k][gp]*ddPIIe(0)*ilambin_sq*ilambin_sq,Av_[k],Av_[k],1.0);

    // update active stress and elasticity tensor contribution
    potsumfiberact_[k-potsumfiberpas_.size()]->EvaluateActiveStressCmatAniso(CM,cmatactive,stressactive,eleGID);

    stress.Update(cur_rho_col_[k][gp],stressactive,1.0);
    cmat.Update(cur_rho_col_[k][gp],cmatactive,1.0);


    // fiber stress for output
    potsumfiberact_[k-potsumfiberpas_.size()]->GetFiberVecs(fibervecs);
    stress_[k][gp] = CM.Dot(AM_[k])*2.0*dPIe(0)*ilambin_sq;

    // active contribution
    stress_f_act = stressactive.Dot(A_strain_[k]);
    stress_[k][gp] += CM.Dot(AM_[k])*stress_f_act;
  }

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::RemodelFiber::EvaluateDerivativesInternalNewton(
  const LINALG::Matrix<3,3>* defgrd,
  const LINALG::Matrix<3,3> id,
  const int nr_grf_proc,
  const double density,
  const int gp,
  const double dt,
  const double v,
  LINALG::Matrix<3,3> FgM,
  LINALG::Matrix<3,3> iFgM,
  LINALG::Matrix<3,3> AgM,
  LINALG::Matrix<3,3> AcirM,
  LINALG::Matrix<3,3> AradM,
  LINALG::Matrix<3,3> AaxM,
  std::vector<std::vector<double> >& dWdrho,
  std::vector<std::vector<double> >& dWdlamb,
  std::vector<double>& W,
  std::vector<std::vector<double> >& dEdrho,
  std::vector<std::vector<double> >& dEdlamb,
  std::vector<double>& E,
  const int eleGID,
  const int growthtype
  )
{
  // constant factor
  double fac = 1.0;

  // converts stress-like to strain-like Voigt notation// build structural tensor in matrix notation
  LINALG::Matrix<6,6> strainconv(true);
  for(int i=0;i<3;++i)
    strainconv(i,i) = 1.0;
  for(int j=3;j<6;++j)
    strainconv(j,j) = 2.0;

  // fiber directions
  std::vector<LINALG::Matrix<3,1> > fibervecs;

  // some variables
  double sigf = 0.0;
  double dsigfdrho = 0.0;
  double ddPIedrho = 0.0;
  double dddPIIedrho = 0.0;
  LINALG::Matrix<2,1> I_pseudo(true);
  LINALG::Matrix<1,1> dI4edlamb(true);
  LINALG::Matrix<1,1> tmp_scal(true);
  LINALG::Matrix<1,6> tmp1x6(true);
  LINALG::Matrix<3,3> tmp3x3_1(true);
  LINALG::Matrix<3,3> tmp3x3_2(true);
  LINALG::Matrix<3,3> FrM(true);
  LINALG::Matrix<3,3> iFrM(true);
  LINALG::Matrix<3,3> iFgiFrM(true);
  LINALG::Matrix<3,3> FeM(true);
  LINALG::Matrix<3,3> CeM(true);
  LINALG::Matrix<2,1> dPIe(true);
  LINALG::Matrix<3,1> ddPIIe(true);
  LINALG::Matrix<4,1> dddPIIIe(true);
  LINALG::Matrix<3,3> diFgdrhoM(true);
  LINALG::Matrix<3,3> dFgdrhoM(true);
  LINALG::Matrix<3,3> diFgdrhoiFrM(true);
  LINALG::Matrix<3,3> dCedrhoM(true);
  LINALG::Matrix<6,1> dCedrhov(true);
  LINALG::Matrix<1,1> dI4edrho(true);
  LINALG::Matrix<1,6> dsigfdCe(true);
  LINALG::Matrix<3,3> diFrdlambM(true);
  LINALG::Matrix<3,3> dFrdlambM(true);
  LINALG::Matrix<3,3> dCedlambM(true);
  LINALG::Matrix<6,1> dCedlambv(true);
  LINALG::Matrix<6,1> dCedlamb_strain(true);
  LINALG::Matrix<1,1> dsigfdlamb(true);
  LINALG::Matrix<1,1> ddPIedlamb(true);
  LINALG::Matrix<1,1> dddPIIedlamb(true);
  LINALG::Matrix<3,3> FrFgM(true);
  LINALG::Matrix<1,6> dsigfdCedlamb(true);
  LINALG::Matrix<3,3> FrnM(true);
  LINALG::Matrix<3,3> FrniFrM(true);
  LINALG::Matrix<3,3> FrdotiFrM(true);
  LINALG::Matrix<3,3> YM(true);
  LINALG::Matrix<6,1> Y_strain(true);
  LINALG::Matrix<3,3> dYdlambM(true);
  LINALG::Matrix<6,1> dYdlambv(true);
  LINALG::Matrix<3,3> dYdrhoM(true);
  LINALG::Matrix<6,1> dYdrhov(true);
  LINALG::Matrix<1,6> dsigfdCedrho(true);

  // additional variables for active fiber contribution (smooth muscle)
  double stress_f_act = 0.0;
  LINALG::Matrix<6,1> stressactive(true);
  LINALG::Matrix<6,6> cmatactive(true);
  LINALG::Matrix<3,1> Fa(true);

  // right Cauchy Green tensor
  LINALG::Matrix<3,3> CM(true);
  CM.MultiplyTN(1.0,*defgrd,*defgrd,0.0);


  // passive fiber evaluation
  for(unsigned k=0;k<potsumfiberpas_.size();++k)
  {
    // Get fiberdirection
    potsumfiberpas_[k]->GetFiberVecs(fibervecs);

    fac = dt*params_->k_growth_/sigmapre_[k];

    // build remodel deformation gradient
    FrM.Update(cur_lambda_r_[k][gp]*(1./params_->G_),AM_[k],0.0);
    FrM.Update(1./std::sqrt(cur_lambda_r_[k][gp]*1./params_->G_),AM_orth_[k],1.0);

    // inverse remodel deformation gradient
    iFrM.Invert(FrM);

    // Fg^-1 * Fr^-1
    iFgiFrM.MultiplyNN(1.0,iFgM,iFrM,0.0);

    // elastic right Cauchy Green tensor
    FeM.MultiplyNN(1.0,*defgrd,iFgiFrM,0.0);
    CeM.MultiplyTN(1.0,FeM,FeM,0.0);

    // get derivatives of strain energy function w.r.t. the fourth invariant
    potsumfiberpas_[k]->GetDerivativesAniso(dPIe,ddPIIe,dddPIIIe,CeM,eleGID);

    // Evaluate fiber Cauchy stress in updated fiber direction
    I_pseudo(0) = CeM.Dot(AM_[k]);

    sigf = 2.*I_pseudo(0)*dPIe(0);


    switch(growthtype)
    {
    case 1:
      // dFg^-1/drho
      diFgdrhoM.Update(-(1./(v*v))*(1./density),AgM,0.0);
      // dFg/drho
      dFgdrhoM.Update(1./density,AgM,0.0);
      break;
    case 0:
      // dFg^-1/drho
      diFgdrhoM.Update(-1./3.*std::pow(v,-4./3.)*(1./density),AcirM,0.0);
      diFgdrhoM.Update(-1./3.*std::pow(v,-4./3.)*(1./density),AradM,1.0);
      diFgdrhoM.Update(-1./3.*std::pow(v,-4./3.)*(1./density),AaxM,1.0);
      // dFg/drho (isotropic growth)
      dFgdrhoM.Update(1./3.*std::pow(v,-2./3.)*(1./density),AcirM,0.0);
      dFgdrhoM.Update(1./3.*std::pow(v,-2./3.)*(1./density),AradM,1.0);
      dFgdrhoM.Update(1./3.*std::pow(v,-2./3.)*(1./density),AaxM,1.0);
     break;
    default:
      dserror("growthtype has to be either 1: anisotropic growth or 0: isotropic growth");
      break;
    }

    // dCedrho (lambda_r const)
    diFgdrhoiFrM.MultiplyNN(1.0,diFgdrhoM,iFrM,0.0);
    tmp3x3_1.MultiplyTN(1.0,diFgdrhoiFrM,CM,0.0);
    dCedrhoM.MultiplyNN(1.0,tmp3x3_1,iFgiFrM,0.0);
    tmp3x3_1.MultiplyNN(1.0,CM,diFgdrhoiFrM,0.0);
    dCedrhoM.MultiplyTN(1.0,iFgiFrM,tmp3x3_1,1.0);
    MatrixtoStressVoigtNotationVector(dCedrhoM,dCedrhov);

    // d(dPIe)/drho = ddPII * dI4e/drho (lambda_r const)
    dI4edrho.MultiplyTN(1.0,A_strain_[k],dCedrhov,0.0);
    ddPIedrho = ddPIIe(0) * dI4edrho(0);

    //dsigf/drho (lambda_r const)
    dsigfdrho = 2.*dI4edrho(0)*dPIe(0);

    dsigfdrho += 2.*I_pseudo(0)*ddPIedrho;



    // dWdrho
    dWdrho[nr_grf_proc+k][nr_grf_proc+k] = 1. - fac*(sigf-sigmapre_[k]);
    dWdrho[nr_grf_proc+k][nr_grf_proc+k] += -fac*cur_rho_col_[k][gp]*dsigfdrho;

    for(unsigned m=0;m<dWdrho[nr_grf_proc+k].size();++m)
      if((int)m != nr_grf_proc+(int)k)
        dWdrho[nr_grf_proc+k][m] = -cur_rho_col_[k][gp]*fac*dsigfdrho;



    // dsigfdCe
    dsigfdCe.UpdateT(2.*dPIe(0),Av_[k],0.0);

    dsigfdCe.UpdateT(2.*I_pseudo(0)*ddPIIe(0),Av_[k],1.0);


    // dFr^-1/dlambda_r
    diFrdlambM.Update(-std::pow(cur_lambda_r_[k][gp],-2.)*params_->G_,AM_[k],0.0);
    diFrdlambM.Update(0.5*std::pow(cur_lambda_r_[k][gp],-0.5)*std::sqrt(1./params_->G_),AM_orth_[k],1.0);


    // dFr/dlambda_r
    dFrdlambM.Update(1./params_->G_,AM_[k],0.0);
    dFrdlambM.Update(-0.5*std::pow(cur_lambda_r_[k][gp],-1.5)*1./std::sqrt(1./params_->G_),AM_orth_[k],1.0);


    // dCe/dlambda_r
    tmp3x3_1.MultiplyNN(1.0,iFgM,diFrdlambM,0.0);
    tmp3x3_2.MultiplyTN(1.0,tmp3x3_1,CM,0.0);
    dCedlambM.MultiplyNN(1.0,tmp3x3_2,iFgiFrM,0.0);

    tmp3x3_2.MultiplyTN(1.0,iFgiFrM,CM,0.0);
    dCedlambM.MultiplyNN(1.0,tmp3x3_2,tmp3x3_1,1.0);
    MatrixtoStressVoigtNotationVector(dCedlambM,dCedlambv);


    // dsigf/dlambda_r
    dCedlamb_strain.Update(1.0,dCedlambv,0.0);
    for(int i=3;i<6;++i)
      dCedlamb_strain(i) = 2.*dCedlambv(i);

    dsigfdlamb.MultiplyNN(1.0,dsigfdCe,dCedlamb_strain,0.0);


    // dPI/dlambda_r
    ddPIedlamb.MultiplyTN(ddPIIe(0),A_strain_[k],dCedlambv,0.0);


    // ddPII/dlambda_r
    dddPIIedlamb.MultiplyTN(dddPIIIe(0),A_strain_[k],dCedlambv,0.0);


    // dI4e/dlambda_r
    dI4edlamb.MultiplyTN(1.0,Av_[k],dCedlamb_strain,0.0);


    // Fr * Fg
    FrFgM.MultiplyNN(1.0,FrM,FgM,0.0);

    // (dsigf/dCe)/dlambda_r
    dsigfdCedlamb.UpdateT(2.*ddPIedlamb(0),Av_[k],0.0);

    dsigfdCedlamb.UpdateT(2.*I_pseudo(0)*dddPIIedlamb(0),Av_[k],1.0);

    dsigfdCedlamb.UpdateT(2.*ddPIIe(0)*dI4edlamb(0),Av_[k],1.0);



    // dY/dlambda_r with Y = (Ce*Frdot*Fr^-1 + Fr^-T*Frdot^T*Ce)
    FrnM.Update(last_lambda_r_[k][gp]*(1./params_->G_),AM_[k],0.0);
    FrnM.Update(1./std::sqrt(last_lambda_r_[k][gp]*1./params_->G_),AM_orth_[k],1.0);
    FrniFrM.MultiplyNN(1.0,FrnM,iFrM,0.0);

    FrdotiFrM.Update(1./dt,id,0.0);
    FrdotiFrM.Update(-1./dt,FrniFrM,1.0);
    YM.MultiplyNN(1.0,CeM,FrdotiFrM,0.0);
    YM.MultiplyTN(1.0,FrdotiFrM,CeM,1.0);
    dYdlambM.MultiplyNN(1.0,dCedlambM,FrdotiFrM,0.0);
    dYdlambM.MultiplyTN(1.0,FrdotiFrM,dCedlambM,1.0);

    tmp3x3_1.MultiplyNN(-1./dt,FrnM,diFrdlambM,0.0);
    dYdlambM.MultiplyNN(1.0,CeM,tmp3x3_1,1.0);
    dYdlambM.MultiplyTN(1.0,tmp3x3_1,CeM,1.0);
    MatrixtoStressVoigtNotationVector(dYdlambM,dYdlambv);


    // dE/dlambda_r
    dEdlamb[nr_grf_proc+k][nr_grf_proc+k] = 2.*(params_->k_growth_/sigmapre_[k])*(sigf-sigmapre_[k])*dsigfdlamb(0);

    dEdlamb[nr_grf_proc+k][nr_grf_proc+k] += dsigfdlamb(0)/params_->tdecay_;

    for(int i=0;i<3;++i)
      Y_strain(i) = YM(i,i);
    Y_strain(3) = YM(0,1) + YM(1,0);
    Y_strain(4) = YM(1,2) + YM(2,1);
    Y_strain(5) = YM(0,2) + YM(2,0);
    tmp_scal.MultiplyNN(1.0,dsigfdCedlamb,Y_strain,0.0);
    dEdlamb[nr_grf_proc+k][nr_grf_proc+k] -= tmp_scal(0);

    tmp1x6.MultiplyNN(1.0,dsigfdCe,strainconv,0.0);
    tmp_scal.MultiplyNN(1.0,tmp1x6,dYdlambv,0.0);
    dEdlamb[nr_grf_proc+k][nr_grf_proc+k] -= tmp_scal(0);



    // dW/dlambda_r
    dWdlamb[nr_grf_proc+k][nr_grf_proc+k] = -cur_rho_col_[k][gp]*fac*dsigfdlamb(0);



    // dY/drho
    dYdrhoM.MultiplyNN(1.0,dCedrhoM,FrdotiFrM,0.0);
    dYdrhoM.MultiplyTN(1.0,FrdotiFrM,dCedrhoM,1.0);
    MatrixtoStressVoigtNotationVector(dYdrhoM,dYdrhov);


    // d(ddPIIe)/drho
    dddPIIedrho = dddPIIIe(0) * dI4edrho(0);


    // (dsigf/dCe)/drho
    dsigfdCedrho.UpdateT(2.*ddPIIe(0)*dI4edrho(0),Av_[k],0.0);

    dsigfdCedrho.UpdateT(2.*ddPIedrho,Av_[k],1.0);

    dsigfdCedrho.UpdateT(2.*I_pseudo(0)*dddPIIedrho,Av_[k],1.0);



    // dE/drho
    dEdrho[nr_grf_proc+k][nr_grf_proc+k] = 2.*(params_->k_growth_/sigmapre_[k])*(sigf-sigmapre_[k])*dsigfdrho;

    dEdrho[nr_grf_proc+k][nr_grf_proc+k] += dsigfdrho/params_->tdecay_;

    tmp_scal.MultiplyNN(1.0,dsigfdCedrho,Y_strain,0.0);
    dEdrho[nr_grf_proc+k][nr_grf_proc+k] -= tmp_scal(0);

    tmp1x6.MultiplyNN(1.0,dsigfdCe,strainconv,0.0);
    tmp_scal.MultiplyNN(1.0,tmp1x6,dYdrhov,0.0);
    dEdrho[nr_grf_proc+k][nr_grf_proc+k] -= tmp_scal(0);

    for(unsigned m=0;m<dEdrho[nr_grf_proc+k].size();++m)
      if((int)m != nr_grf_proc+(int)k)
        dEdrho[nr_grf_proc+k][m] = dEdrho[nr_grf_proc+k][nr_grf_proc+k];



    // single residuals
    tmp_scal.MultiplyNN(1.0,dsigfdCe,Y_strain,0.0);
    E[nr_grf_proc+k] = ((params_->k_growth_/sigmapre_[k])*(sigf-sigmapre_[k]) + 1./params_->tdecay_)*(sigf-sigmapre_[k])-tmp_scal(0);

    W[nr_grf_proc+k] = cur_rho_col_[k][gp]*(1. - fac*(sigf-sigmapre_[k])) - last_rho_col_[k][gp];
  }
  // active fiber evaluation
  for(unsigned k=potsumfiberpas_.size();k<(potsumfiberpas_.size()+potsumfiberact_.size());++k)
  {
    // Get fiberdirection
    potsumfiberact_[k-potsumfiberpas_.size()]->GetFiberVecs(fibervecs);

    fac = dt*params_->k_growth_/sigmapre_[k];

    // build remodel deformation gradient
    FrM.Update(cur_lambda_r_[k][gp]*(1./params_->G_),AM_[k],0.0);
    FrM.Update(1./std::sqrt(cur_lambda_r_[k][gp]*1./params_->G_),AM_orth_[k],1.0);

    // inverse remodel deformation gradient
    iFrM.Invert(FrM);

    // Fg^-1 * Fr^-1
    iFgiFrM.MultiplyNN(1.0,iFgM,iFrM,0.0);

    // elastic right Cauchy Green tensor
    FeM.MultiplyNN(1.0,*defgrd,iFgiFrM,0.0);
    CeM.MultiplyTN(1.0,FeM,FeM,0.0);

    // get derivatives of strain energy function w.r.t. the fourth invariant
    potsumfiberact_[k-potsumfiberpas_.size()]->GetDerivativesAniso(dPIe,ddPIIe,dddPIIIe,CeM,eleGID);

    // Evaluate fiber Cauchy stress in updated fiber direction
    I_pseudo(0) = CeM.Dot(AM_[k]);

    sigf = 2.*I_pseudo(0)*dPIe(0);

    // Evaluate fiber Cauchy stress in updated fiber direction (active contribution)
    potsumfiberact_[k-potsumfiberpas_.size()]->EvaluateActiveStressCmatAniso(*defgrd,cmatactive,stressactive,eleGID);
    stress_f_act = stressactive.Dot(A_strain_[k]);
    Fa.MultiplyNN(1.0,*defgrd,fibervecs[k],0.0);

    sigf += Fa.Norm2()*Fa.Norm2()*stress_f_act;

    switch(growthtype)
    {
    case 1:
      // dFg^-1/drho
      diFgdrhoM.Update(-(1./(v*v))*(1./density),AgM,0.0);
      // dFg/drho
      dFgdrhoM.Update(1./density,AgM,0.0);
      break;
    case 0:
      // dFg^-1/drho
      diFgdrhoM.Update(-1./3.*std::pow(v,-4./3.)*(1./density),AcirM,0.0);
      diFgdrhoM.Update(-1./3.*std::pow(v,-4./3.)*(1./density),AradM,1.0);
      diFgdrhoM.Update(-1./3.*std::pow(v,-4./3.)*(1./density),AaxM,1.0);
      // dFg/drho (isotropic growth)
      dFgdrhoM.Update(1./3.*std::pow(v,-2./3.)*(1./density),AcirM,0.0);
      dFgdrhoM.Update(1./3.*std::pow(v,-2./3.)*(1./density),AradM,1.0);
      dFgdrhoM.Update(1./3.*std::pow(v,-2./3.)*(1./density),AaxM,1.0);
     break;
    default:
      dserror("growthtype has to be either 1: anisotropic growth or 0: isotropic growth");
      break;
    }

    // dCedrho (lambda_r const)
    diFgdrhoiFrM.MultiplyNN(1.0,diFgdrhoM,iFrM,0.0);
    tmp3x3_1.MultiplyTN(1.0,diFgdrhoiFrM,CM,0.0);
    dCedrhoM.MultiplyNN(1.0,tmp3x3_1,iFgiFrM,0.0);
    tmp3x3_1.MultiplyNN(1.0,CM,diFgdrhoiFrM,0.0);
    dCedrhoM.MultiplyTN(1.0,iFgiFrM,tmp3x3_1,1.0);
    MatrixtoStressVoigtNotationVector(dCedrhoM,dCedrhov);

    // d(dPIe)/drho = ddPII * dI4e/drho (lambda_r const)
    dI4edrho.MultiplyTN(1.0,A_strain_[k],dCedrhov,0.0);
    ddPIedrho = ddPIIe(0) * dI4edrho(0);

    //dsigf/drho (lambda_r const)
    dsigfdrho = 2.*dI4edrho(0)*dPIe(0);

    dsigfdrho += 2.*I_pseudo(0)*ddPIedrho;



    // dWdrho
    dWdrho[nr_grf_proc+k][nr_grf_proc+k] = 1. - fac*(sigf-sigmapre_[k]);
    dWdrho[nr_grf_proc+k][nr_grf_proc+k] += -fac*cur_rho_col_[k][gp]*dsigfdrho;

    for(unsigned m=0;m<dWdrho[nr_grf_proc+k].size();++m)
      if((int)m != nr_grf_proc+(int)k)
        dWdrho[nr_grf_proc+k][m] = -cur_rho_col_[k][gp]*fac*dsigfdrho;



    // dsigfdCe
    dsigfdCe.UpdateT(2.*dPIe(0),Av_[k],0.0);

    dsigfdCe.UpdateT(2.*I_pseudo(0)*ddPIIe(0),Av_[k],1.0);


    // dFr^-1/dlambda_r
    diFrdlambM.Update(-std::pow(cur_lambda_r_[k][gp],-2.)*params_->G_,AM_[k],0.0);
    diFrdlambM.Update(0.5*std::pow(cur_lambda_r_[k][gp],-0.5)*std::sqrt(1./params_->G_),AM_orth_[k],1.0);


    // dFr/dlambda_r
    dFrdlambM.Update(1./params_->G_,AM_[k],0.0);
    dFrdlambM.Update(-0.5*std::pow(cur_lambda_r_[k][gp],-1.5)*1./std::sqrt(1./params_->G_),AM_orth_[k],1.0);


    // dCe/dlambda_r
    tmp3x3_1.MultiplyNN(1.0,iFgM,diFrdlambM,0.0);
    tmp3x3_2.MultiplyTN(1.0,tmp3x3_1,CM,0.0);
    dCedlambM.MultiplyNN(1.0,tmp3x3_2,iFgiFrM,0.0);

    tmp3x3_2.MultiplyTN(1.0,iFgiFrM,CM,0.0);
    dCedlambM.MultiplyNN(1.0,tmp3x3_2,tmp3x3_1,1.0);
    MatrixtoStressVoigtNotationVector(dCedlambM,dCedlambv);


    // dsigf/dlambda_r
    dCedlamb_strain.Update(1.0,dCedlambv,0.0);
    for(int i=3;i<6;++i)
      dCedlamb_strain(i) = 2.*dCedlambv(i);

    dsigfdlamb.MultiplyNN(1.0,dsigfdCe,dCedlamb_strain,0.0);


    // dPI/dlambda_r
    ddPIedlamb.MultiplyTN(ddPIIe(0),A_strain_[k],dCedlambv,0.0);


    // ddPII/dlambda_r
    dddPIIedlamb.MultiplyTN(dddPIIIe(0),A_strain_[k],dCedlambv,0.0);


    // d(Fea_norm)^2/dlambda_r
    dI4edlamb.MultiplyTN(1.0,Av_[k],dCedlamb_strain,0.0);


    // Fr * Fg
    FrFgM.MultiplyNN(1.0,FrM,FgM,0.0);


    // (dsigf/dCe)/dlambda_r
    dsigfdCedlamb.UpdateT(2.*ddPIedlamb(0),Av_[k],0.0);

    dsigfdCedlamb.UpdateT(2.*I_pseudo(0)*dddPIIedlamb(0),Av_[k],1.0);

    dsigfdCedlamb.UpdateT(2.*ddPIIe(0)*dI4edlamb(0),Av_[k],1.0);



    // dY/dlambda_r with Y = (Ce*Frdot*Fr^-1 + Fr^-T*Frdot^T*Ce)
    FrnM.Update(last_lambda_r_[k][gp]*(1./params_->G_),AM_[k],0.0);
    FrnM.Update(1./std::sqrt(last_lambda_r_[k][gp]*1./params_->G_),AM_orth_[k],1.0);
    FrniFrM.MultiplyNN(1.0,FrnM,iFrM,0.0);

    FrdotiFrM.Update(1./dt,id,0.0);
    FrdotiFrM.Update(-1./dt,FrniFrM,1.0);
    YM.MultiplyNN(1.0,CeM,FrdotiFrM,0.0);
    YM.MultiplyTN(1.0,FrdotiFrM,CeM,1.0);
    dYdlambM.MultiplyNN(1.0,dCedlambM,FrdotiFrM,0.0);
    dYdlambM.MultiplyTN(1.0,FrdotiFrM,dCedlambM,1.0);

    tmp3x3_1.MultiplyNN(-1./dt,FrnM,diFrdlambM,0.0);
    dYdlambM.MultiplyNN(1.0,CeM,tmp3x3_1,1.0);
    dYdlambM.MultiplyTN(1.0,tmp3x3_1,CeM,1.0);
    MatrixtoStressVoigtNotationVector(dYdlambM,dYdlambv);


    // dE/dlambda_r
    dEdlamb[nr_grf_proc+k][nr_grf_proc+k] = 2.*(params_->k_growth_/sigmapre_[k])*(sigf-sigmapre_[k])*dsigfdlamb(0);

    dEdlamb[nr_grf_proc+k][nr_grf_proc+k] += dsigfdlamb(0)/params_->tdecay_;

    for(int i=0;i<3;++i)
      Y_strain(i) = YM(i,i);
    Y_strain(3) = YM(0,1) + YM(1,0);
    Y_strain(4) = YM(1,2) + YM(2,1);
    Y_strain(5) = YM(0,2) + YM(2,0);
    tmp_scal.MultiplyNN(1.0,dsigfdCedlamb,Y_strain,0.0);
    dEdlamb[nr_grf_proc+k][nr_grf_proc+k] -= tmp_scal(0);

    tmp1x6.MultiplyNN(1.0,dsigfdCe,strainconv,0.0);
    tmp_scal.MultiplyNN(1.0,tmp1x6,dYdlambv,0.0);
    dEdlamb[nr_grf_proc+k][nr_grf_proc+k] -= tmp_scal(0);



    // dW/dlambda_r
    dWdlamb[nr_grf_proc+k][nr_grf_proc+k] = -cur_rho_col_[k][gp]*fac*dsigfdlamb(0);



    // dY/drho
    dYdrhoM.MultiplyNN(1.0,dCedrhoM,FrdotiFrM,0.0);
    dYdrhoM.MultiplyTN(1.0,FrdotiFrM,dCedrhoM,1.0);
    MatrixtoStressVoigtNotationVector(dYdrhoM,dYdrhov);


    // d(ddPIIe)/drho
    dddPIIedrho = dddPIIIe(0) * dI4edrho(0);


    // (dsigf/dCe)/drho
    dsigfdCedrho.UpdateT(2.*ddPIIe(0)*dI4edrho(0),Av_[k],0.0);

    dsigfdCedrho.UpdateT(2.*ddPIedrho,Av_[k],1.0);

    dsigfdCedrho.UpdateT(2.*I_pseudo(0)*dddPIIedrho,Av_[k],1.0);


    // dE/drho
    dEdrho[nr_grf_proc+k][nr_grf_proc+k] = 2.*(params_->k_growth_/sigmapre_[k])*(sigf-sigmapre_[k])*dsigfdrho;

    dEdrho[nr_grf_proc+k][nr_grf_proc+k] += dsigfdrho/params_->tdecay_;

    tmp_scal.MultiplyNN(1.0,dsigfdCedrho,Y_strain,0.0);
    dEdrho[nr_grf_proc+k][nr_grf_proc+k] -= tmp_scal(0);

    tmp1x6.MultiplyNN(1.0,dsigfdCe,strainconv,0.0);
    tmp_scal.MultiplyNN(1.0,tmp1x6,dYdrhov,0.0);
    dEdrho[nr_grf_proc+k][nr_grf_proc+k] -= tmp_scal(0);

    for(unsigned m=0;m<dEdrho[nr_grf_proc+k].size();++m)
      if((int)m != nr_grf_proc+(int)k)
        dEdrho[nr_grf_proc+k][m] = dEdrho[nr_grf_proc+k][nr_grf_proc+k];



    // single residuals
    tmp_scal.MultiplyNN(1.0,dsigfdCe,Y_strain,0.0);
    E[nr_grf_proc+k] = ((params_->k_growth_/sigmapre_[k])*(sigf-sigmapre_[k]) + 1./params_->tdecay_)*(sigf-sigmapre_[k])-tmp_scal(0);

    W[nr_grf_proc+k] = cur_rho_col_[k][gp]*(1. - fac*(sigf-sigmapre_[k])) - last_rho_col_[k][gp];
  }

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::RemodelFiber::EvaluateDerivativesCauchyGreen(
  const LINALG::Matrix<3,3>* defgrd,
  const LINALG::Matrix<3,3> id,
  const int nr_grf_proc,
  const int gp,
  const double dt,
  LINALG::Matrix<3,3> FgM,
  LINALG::Matrix<3,3> iFgM,
  std::vector<LINALG::Matrix<1,6> >& dWdC,
  std::vector<LINALG::Matrix<1,6> >& dEdC,
  const int eleGID
  )
{
  // some variables
  LINALG::Matrix<1,6> dsigfdC(true);
  LINALG::Matrix<1,6> dI4edC(true);

  LINALG::Matrix<6,6> dsigfdCedC(true);

  double sigf = 0.0;
  LINALG::Matrix<3,3> tmp3x3(true);
  LINALG::Matrix<6,6> tmp6x6(true);
  LINALG::Matrix<2,1> I_pseudo(true);
  LINALG::Matrix<3,3> FeM(true);
  LINALG::Matrix<3,3> CeM(true);
  LINALG::Matrix<2,1> dPIe(true);
  LINALG::Matrix<3,1> ddPIIe(true);
  LINALG::Matrix<4,1> dddPIIIe(true);
  std::vector<LINALG::Matrix<3,1> > fibervecs;
  LINALG::Matrix<3,3> FrM(true);
  LINALG::Matrix<3,3> iFrM(true);
  LINALG::Matrix<3,3> iFgiFrM(true);
  LINALG::Matrix<6,6> dYdC(true);
  LINALG::Matrix<1,6> dI4edCv(true);
  LINALG::Matrix<1,6> ddPIedCv(true);
  LINALG::Matrix<1,6> dddPIIedCv(true);
  LINALG::Matrix<6,6> dCedC(true);
  LINALG::Matrix<6,6> dCedC_strain(true);
  LINALG::Matrix<3,3> FrdotiFrM(true);
  LINALG::Matrix<3,3> FrnM(true);
  LINALG::Matrix<3,3> FrniFrM(true);
  LINALG::Matrix<1,6> dsigfdCe(true);
  LINALG::Matrix<3,3> YM(true);
  LINALG::Matrix<6,1> Y_strain(true);

  // additional variables for active fiber contribution (smooth muscle)
  double stress_f_act = 0.0;
  LINALG::Matrix<6,1> stressactive(true);
  LINALG::Matrix<6,6> cmatactive(true);
  LINALG::Matrix<3,1> Fa(true);

  // converts stress-like to strain-like Voigt notation// build structural tensor in matrix notation
  LINALG::Matrix<6,6> strainconv(true);
  for(int i=0;i<3;++i)
    strainconv(i,i) = 1.0;
  for(int j=3;j<6;++j)
    strainconv(j,j) = 2.0;

  // constant factor
  double fac = 1.0;


  // passive fiber evaluation
  for(unsigned k=0;k<potsumfiberpas_.size();++k)
  {
    // Get fiberdirection
    potsumfiberpas_[k]->GetFiberVecs(fibervecs);

    fac = dt*params_->k_growth_/sigmapre_[k];

    // build remodel deformation gradient
    FrM.Update(cur_lambda_r_[k][gp]*(1./params_->G_),AM_[k],0.0);
    FrM.Update(1./std::sqrt(cur_lambda_r_[k][gp]*1./params_->G_),AM_orth_[k],1.0);

    // inverse remodel deformation gradient
    iFrM.Invert(FrM);

    // Fg^-1 * Fr^-1
    iFgiFrM.MultiplyNN(1.0,iFgM,iFrM,0.0);

    // dCedC
    EvaldCedC(dCedC,iFgiFrM,0.5);
    dCedC_strain.MultiplyNN(1.0,strainconv,dCedC,0.0);

    // elastic right Cauchy Green tensor
    FeM.MultiplyNN(1.0,*defgrd,iFgiFrM,0.0);
    CeM.MultiplyTN(1.0,FeM,FeM,0.0);

    // get derivatives of strain energy function w.r.t. the fourth invariant
    potsumfiberpas_[k]->GetDerivativesAniso(dPIe,ddPIIe,dddPIIIe,CeM,eleGID);

    // Evaluate fiber Cauchy stress in updated fiber direction
    I_pseudo(0) = CeM.Dot(AM_[k]);


    // dsigf/dC (passive contribution)
    dI4edC.MultiplyTN(1.0,A_strain_[k],dCedC,0.0);
    dsigfdC.Update(2.*dPIe(0),dI4edC,0.0);

    dsigfdC.Update(2.*I_pseudo(0)*ddPIIe(0),dI4edC,1.0);


    // derivation of the growth evolution eq. w.r.t. right Cauchy Green
    dWdC[nr_grf_proc+k].Update(-fac*cur_rho_col_[k][gp],dsigfdC,0.0);



    // dY/dC
    FrnM.Update(last_lambda_r_[k][gp]*(1./params_->G_),AM_[k],0.0);
    FrnM.Update(1./std::sqrt(last_lambda_r_[k][gp]*1./params_->G_),AM_orth_[k],1.0);
    FrniFrM.MultiplyNN(1.0,FrnM,iFrM,0.0);
    FrdotiFrM.Update(1./dt,id,0.0);
    FrdotiFrM.Update(-1./dt,FrniFrM,1.0);

    tmp3x3.MultiplyNN(1.0,iFgiFrM,FrdotiFrM,0.0);
    EvaldYdC(dYdC,tmp3x3,iFgiFrM,0.5);


    // d(Norm(Fe*a))^2/dC
    dI4edCv.MultiplyTN(1.0,Av_[k],dCedC_strain,0.0);


    // d(dPIe)/dC
    ddPIedCv.Update(ddPIIe(0),dI4edC,0.0);


    // d(ddPIIe)/dC
    dddPIIedCv.Update(dddPIIIe(0),dI4edC,0.0);



    // d(sigf/dCe)/dC
    dsigfdCedC.MultiplyNN(2.*ddPIIe(0),Av_[k],dI4edCv,0.0);

    dsigfdCedC.MultiplyNN(2.,Av_[k],ddPIedCv,1.0);

    dsigfdCedC.MultiplyNN(2.*I_pseudo(0),Av_[k],dddPIIedCv,1.0);


    // dsigfdCe
    dsigfdCe.UpdateT(2.*dPIe(0),Av_[k],1.0);

    dsigfdCe.UpdateT(2.*I_pseudo(0)*ddPIIe(0),Av_[k],1.0);


    // fiber Cauchy stress (passive contribution)
    sigf = 2.*I_pseudo(0)*dPIe(0);


    // Y = (Ce*Frdot*Fr^-1 + Fr^-T*Frdot^T*Ce)
    YM.MultiplyNN(1.0,CeM,FrdotiFrM,0.0);
    YM.MultiplyTN(1.0,FrdotiFrM,CeM,1.0);

    for(int i=0;i<3;++i)
      Y_strain(i) = YM(i,i);
    Y_strain(3) = YM(0,1)+YM(1,0);
    Y_strain(4) = YM(1,2)+YM(2,1);
    Y_strain(5) = YM(0,2)+YM(2,0);



    // dE/dC
    dEdC[nr_grf_proc+k].Update(2.*(params_->k_growth_/sigmapre_[k])*(sigf-sigmapre_[k]),dsigfdC,0.0);

    dEdC[nr_grf_proc+k].Update(1./params_->tdecay_,dsigfdC,1.0);

    dEdC[nr_grf_proc+k].MultiplyTN(-1.0,Y_strain,dsigfdCedC,1.0);

    tmp6x6.MultiplyNN(1.0,strainconv,dYdC,0.0);
    dEdC[nr_grf_proc+k].MultiplyNN(-1.0,dsigfdCe,tmp6x6,1.0);

    // fiber stress for output
    stress_[k][gp] = sigf;
  }
  // active fiber evaluation
  for(unsigned k=potsumfiberpas_.size();k<(potsumfiberpas_.size()+potsumfiberact_.size());++k)
  {
    // Get fiberdirection
    potsumfiberact_[k-potsumfiberpas_.size()]->GetFiberVecs(fibervecs);

    fac = dt*params_->k_growth_/sigmapre_[k];

    // build remodel deformation gradient
    FrM.Update(cur_lambda_r_[k][gp]*(1./params_->G_),AM_[k],0.0);
    FrM.Update(1./std::sqrt(cur_lambda_r_[k][gp]*1./params_->G_),AM_orth_[k],1.0);

    // inverse remodel deformation gradient
    iFrM.Invert(FrM);

    // Fg^-1 * Fr^-1
    iFgiFrM.MultiplyNN(1.0,iFgM,iFrM,0.0);

    // dCedC
    EvaldCedC(dCedC,iFgiFrM,0.5);
    dCedC_strain.MultiplyNN(1.0,strainconv,dCedC,0.0);

    // elastic right Cauchy Green tensor
    FeM.MultiplyNN(1.0,*defgrd,iFgiFrM,0.0);
    CeM.MultiplyTN(1.0,FeM,FeM,0.0);

    // get derivatives of strain energy function w.r.t. the fourth invariant
    potsumfiberact_[k-potsumfiberpas_.size()]->GetDerivativesAniso(dPIe,ddPIIe,dddPIIIe,CeM,eleGID);

    // Evaluate fiber Cauchy stress in updated fiber direction
    I_pseudo(0) = CeM.Dot(AM_[k]);


    // dsigf/dC (passive contribution)
    dI4edC.MultiplyTN(1.0,A_strain_[k],dCedC,0.0);
    dsigfdC.Update(2.*dPIe(0),dI4edC,0.0);

    dsigfdC.Update(2.*I_pseudo(0)*ddPIIe(0),dI4edC,1.0);

    // dsigf/dC (active contribution)
    potsumfiberact_[k-potsumfiberpas_.size()]->EvaluateActiveStressCmatAniso(*defgrd,cmatactive,stressactive,eleGID);
    stress_f_act = stressactive.Dot(A_strain_[k]);
    Fa.MultiplyNN(1.0,*defgrd,fibervecs[k],0.0);

    dsigfdC.MultiplyTN(Fa.Norm2()*Fa.Norm2(),A_strain_[k],cmatactive,1.0);
    dsigfdC.UpdateT(stress_f_act,Av_[k],1.0);


    // derivation of the growth evolution eq. w.r.t. right Cauchy Green
    dWdC[nr_grf_proc+k].Update(-fac*cur_rho_col_[k][gp],dsigfdC,0.0);



    // dY/dC
    FrnM.Update(last_lambda_r_[k][gp]*(1./params_->G_),AM_[k],0.0);
    FrnM.Update(1./std::sqrt(last_lambda_r_[k][gp]*1./params_->G_),AM_orth_[k],1.0);
    FrniFrM.MultiplyNN(1.0,FrnM,iFrM,0.0);
    FrdotiFrM.Update(1./dt,id,0.0);
    FrdotiFrM.Update(-1./dt,FrniFrM,1.0);

    tmp3x3.MultiplyNN(1.0,iFgiFrM,FrdotiFrM,0.0);
    EvaldYdC(dYdC,tmp3x3,iFgiFrM,0.5);


    // d(dPIe)/dC
    ddPIedCv.Update(ddPIIe(0),dI4edC,0.0);


    // d(ddPIIe)/dC
    dddPIIedCv.Update(dddPIIIe(0),dI4edC,0.0);



    // d(sigf/dCe)/dC
    dsigfdCedC.MultiplyNN(2.*ddPIIe(0),Av_[k],dI4edCv,0.0);

    dsigfdCedC.MultiplyNN(2.,Av_[k],ddPIedCv,1.0);

    dsigfdCedC.MultiplyNN(2.*I_pseudo(0),Av_[k],dddPIIedCv,1.0);


    // dsigfdCe
    dsigfdCe.UpdateT(2.*dPIe(0),Av_[k],1.0);

    dsigfdCe.UpdateT(2.*I_pseudo(0)*ddPIIe(0),Av_[k],1.0);


    // fiber Cauchy stress (passive contribution)
    sigf = 2.*I_pseudo(0)*dPIe(0);

    // fiber Cauchy stress (active contribution)
    sigf += Fa.Norm2()*Fa.Norm2()*stress_f_act;


    // Y = (Ce*Frdot*Fr^-1 + Fr^-T*Frdot^T*Ce)
    YM.MultiplyNN(1.0,CeM,FrdotiFrM,0.0);
    YM.MultiplyTN(1.0,FrdotiFrM,CeM,1.0);

    for(int i=0;i<3;++i)
      Y_strain(i) = YM(i,i);
    Y_strain(3) = YM(0,1)+YM(1,0);
    Y_strain(4) = YM(1,2)+YM(2,1);
    Y_strain(5) = YM(0,2)+YM(2,0);



    // dE/dC
    dEdC[nr_grf_proc+k].Update(2.*(params_->k_growth_/sigmapre_[k])*(sigf-sigmapre_[k]),dsigfdC,0.0);

    dEdC[nr_grf_proc+k].Update(1./params_->tdecay_,dsigfdC,1.0);

    dEdC[nr_grf_proc+k].MultiplyTN(-1.0,Y_strain,dsigfdCedC,1.0);

    tmp6x6.MultiplyNN(1.0,strainconv,dYdC,0.0);
    dEdC[nr_grf_proc+k].MultiplyNN(-1.0,dsigfdCe,tmp6x6,1.0);

    // fiber stress for output
    stress_[k][gp] = sigf;
  }

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::RemodelFiber::AddStressCmatGrowthRemodel(
    const LINALG::Matrix<3,3>* defgrd,
    const LINALG::Matrix<3,3> id,
    const int nr_grf_tot,
    const int nr_grf_proc,
    const int gp,
    const double v,
    const double density,
    LINALG::Matrix<3,3> iFgM,
    LINALG::Matrix<3,3> AgM,
    LINALG::Matrix<3,3> AcirM,
    LINALG::Matrix<3,3> AradM,
    LINALG::Matrix<3,3> AaxM,
    std::vector<LINALG::Matrix<1,6> > drhodC,
    std::vector<LINALG::Matrix<1,6> > dlambdC,
    LINALG::Matrix<6,1>& stress,
    LINALG::Matrix<6,6>& cmat,
    const int eleGID,
    const int growthtype
)
{
  // clear variables
  stress.Clear();
  cmat.Clear();

  // some variables
  double ddPIedrho = 0.0;
  LINALG::Matrix<3,3> tmp3x3_1(true);
  LINALG::Matrix<3,3> tmp3x3_2(true);
  LINALG::Matrix<6,6> tmp6x6(true);
  LINALG::Matrix<3,3> FeM(true);
  LINALG::Matrix<3,3> CeM(true);
  LINALG::Matrix<2,1> dPIe(true);
  LINALG::Matrix<3,1> ddPIIe(true);
  LINALG::Matrix<4,1> dddPIIIe(true);
  LINALG::Matrix<3,3> FrM(true);
  LINALG::Matrix<3,3> iFrM(true);
  LINALG::Matrix<3,3> iFgiFrM(true);
  LINALG::Matrix<9,6> dCedCdrho(true);
  LINALG::Matrix<6,6> dCedC(true);
  LINALG::Matrix<6,6> dCedC_strain(true);
  LINALG::Matrix<1,6> stressuw(true);
  LINALG::Matrix<6,1> pk2ev(true);
  LINALG::Matrix<1,1> dI4edrho(true);
  LINALG::Matrix<3,3> diFgdrhoiFrM(true);
  LINALG::Matrix<3,3> dCedrhoM(true);
  LINALG::Matrix<6,1> dCedrhov(true);
  LINALG::Matrix<3,3> diFgdrhoM(true);
  LINALG::Matrix<6,6> cmatelastic(true);
  LINALG::Matrix<3,3> diFrdlambM(true);
  LINALG::Matrix<3,3> dCedlambM(true);
  LINALG::Matrix<6,1> dCedlambv(true);
  LINALG::Matrix<1,1> ddPIedlamb(true);
  LINALG::Matrix<9,6> dCedCdlamb(true);

  // additional variables for active fiber contribution (smooth muscle)
  LINALG::Matrix<6,1> stressactive(true);
  LINALG::Matrix<6,6> cmatactive(true);

  // converts stress-like to strain-like Voigt notation// build structural tensor in matrix notation
  LINALG::Matrix<6,6> strainconv(true);
  for(int i=0;i<3;++i)
    strainconv(i,i) = 1.0;
  for(int j=3;j<6;++j)
    strainconv(j,j) = 2.0;

  // right Cauchy Green tensor
  LINALG::Matrix<3,3> CM(true);
  LINALG::Matrix<3,3> iCM(true);
  CM.MultiplyTN(1.0,*defgrd,*defgrd,0.0);
  iCM.Invert(CM);

  std::vector<std::vector<LINALG::Matrix<1,6> > > dSdrho(nr_grf_tot,std::vector<LINALG::Matrix<1,6> >(nr_grf_tot,LINALG::Matrix<1,6>(true)));
  std::vector<std::vector<LINALG::Matrix<1,6> > > dSdlamb(nr_grf_tot,std::vector<LINALG::Matrix<1,6> >(nr_grf_tot,LINALG::Matrix<1,6>(true)));


  // passive fiber evaluation
  for(unsigned k=0;k<potsumfiberpas_.size();++k)
  {
    // build remodel deformation gradient
    FrM.Update(cur_lambda_r_[k][gp]*(1./params_->G_),AM_[k],0.0);
    FrM.Update(1./std::sqrt(cur_lambda_r_[k][gp]*1./params_->G_),AM_orth_[k],1.0);

    // inverse remodel deformation gradient
    iFrM.Invert(FrM);

    // Fg^-1 * Fr^-1
    iFgiFrM.MultiplyNN(1.0,iFgM,iFrM,0.0);

    // dCedC
    EvaldCedC(dCedC,iFgiFrM,0.5);
    dCedC_strain.MultiplyNN(1.0,strainconv,dCedC,0.0);

    // elastic right Cauchy Green tensor
    FeM.MultiplyNN(1.0,*defgrd,iFgiFrM,0.0);
    CeM.MultiplyTN(1.0,FeM,FeM,0.0);

    // get derivatives of strain energy function w.r.t. the fourth invariant
    potsumfiberpas_[k]->GetDerivativesAniso(dPIe,ddPIIe,dddPIIIe,CeM,eleGID);

    switch(growthtype)
    {
    case 1:
      // dFg^-1/drho
      diFgdrhoM.Update(-(1./(v*v))*(1./density),AgM,0.0);
      break;
    case 0:
      // dFg^-1/drho
      diFgdrhoM.Update(-1./3.*std::pow(v,-4./3.)*(1./density),AcirM,0.0);
      diFgdrhoM.Update(-1./3.*std::pow(v,-4./3.)*(1./density),AradM,1.0);
      diFgdrhoM.Update(-1./3.*std::pow(v,-4./3.)*(1./density),AaxM,1.0);
     break;
    default:
      dserror("growthtype has to be either 1: anisotropic growth or 0: isotropic growth");
      break;
    }

    // dCedrho (lambda_r const)
    diFgdrhoiFrM.MultiplyNN(1.0,diFgdrhoM,iFrM,0.0);
    tmp3x3_1.MultiplyTN(1.0,diFgdrhoiFrM,CM,0.0);
    dCedrhoM.MultiplyNN(1.0,tmp3x3_1,iFgiFrM,0.0);
    tmp3x3_1.MultiplyNN(1.0,CM,diFgdrhoiFrM,0.0);
    dCedrhoM.MultiplyTN(1.0,iFgiFrM,tmp3x3_1,1.0);
    MatrixtoStressVoigtNotationVector(dCedrhoM,dCedrhov);

    // d(dPIe)/drho = ddPII * dI4e/drho (lambda_r const)
    dI4edrho.MultiplyTN(1.0,A_strain_[k],dCedrhov,0.0);
    ddPIedrho = ddPIIe(0) * dI4edrho(0);


    // Update stress
    pk2ev.Update(2.*dPIe(0),Av_[k],0.0);      // evaluate elastic second Piola Kirchhoff stress
    stressuw.MultiplyTN(1.0,pk2ev,dCedC_strain,0.0);

    stress.UpdateT(cur_rho_col_[k][gp],stressuw,1.0);


    // Evaluate derivation of the second Piola Kirchhoff stress w.r.t. rho
    dSdrho[nr_grf_proc+k][nr_grf_proc+k].Update(1.0,stressuw,0.0);

    dSdrho[nr_grf_proc+k][nr_grf_proc+k].MultiplyTN(2.0*cur_rho_col_[k][gp]*ddPIedrho,A_strain_[k],dCedC,1.0);

    tmp3x3_1.UpdateT(1.0,diFgdrhoiFrM,0.0);
    tmp3x3_2.UpdateT(1.0,iFgiFrM,0.0);
    dCedCdrho.Clear();
    AddtoMatrixHalfSymProd(0.5,tmp3x3_1,tmp3x3_2,dCedCdrho);
    AddtoMatrixHalfSymProd(0.5,tmp3x3_2,tmp3x3_1,dCedCdrho);
    dSdrho[nr_grf_proc+k][nr_grf_proc+k].MultiplyTN(2.0*cur_rho_col_[k][gp]*dPIe(0),A9x1_[k],dCedCdrho,1.0);

    for(int m=0;m<nr_grf_tot;++m)
      if(m != nr_grf_proc+(int)k)
      {
        dSdrho[nr_grf_proc+k][m].Update(1.0,dSdrho[nr_grf_proc+k][nr_grf_proc+k],0.0);
        dSdrho[nr_grf_proc+k][m].Update(-1.0,stressuw,1.0);
      }



    // dFr^-1/dlambda_r
    diFrdlambM.Update(-std::pow(cur_lambda_r_[k][gp],-2.)*params_->G_,AM_[k],0.0);
    diFrdlambM.Update(0.5*std::pow(cur_lambda_r_[k][gp],-0.5)*std::sqrt(1./params_->G_),AM_orth_[k],1.0);

    // dCe/dlambda_r
    tmp3x3_1.MultiplyNN(1.0,iFgM,diFrdlambM,0.0);
    tmp3x3_2.MultiplyTN(1.0,tmp3x3_1,CM,0.0);
    dCedlambM.MultiplyNN(1.0,tmp3x3_2,iFgiFrM,0.0);

    tmp3x3_2.MultiplyTN(1.0,iFgiFrM,CM,0.0);
    dCedlambM.MultiplyNN(1.0,tmp3x3_2,tmp3x3_1,1.0);
    MatrixtoStressVoigtNotationVector(dCedlambM,dCedlambv);

    // dPI/dlambda_r
    ddPIedlamb.MultiplyTN(ddPIIe(0),A_strain_[k],dCedlambv,0.0);


    // Evaluate derivation of the second Piola Kirchhoff stress w.r.t. lambda_r
    dSdlamb[nr_grf_proc+k][nr_grf_proc+k].MultiplyTN(2.0*cur_rho_col_[k][gp]*ddPIedlamb(0),A_strain_[k],dCedC,0.0);

    tmp3x3_1.MultiplyTT(1.0,diFrdlambM,iFgM,0.0);
    tmp3x3_2.UpdateT(1.0,iFgiFrM,0.0);
    dCedCdlamb.Clear();
    AddtoMatrixHalfSymProd(0.5,tmp3x3_1,tmp3x3_2,dCedCdlamb);
    AddtoMatrixHalfSymProd(0.5,tmp3x3_2,tmp3x3_1,dCedCdlamb);
    dSdlamb[nr_grf_proc+k][nr_grf_proc+k].MultiplyTN(2.0*cur_rho_col_[k][gp]*dPIe(0),A9x1_[k],dCedCdlamb,1.0);



    // update elasticity tensor
    cmatelastic.MultiplyNT(4.*ddPIIe(0),Av_[k],Av_[k],0.0);
    tmp6x6.MultiplyTN(1.0,dCedC_strain,cmatelastic,0.0);

    cmat.MultiplyNN(cur_rho_col_[k][gp],tmp6x6,dCedC_strain,1.0);

    for(int m=0;m<nr_grf_tot;++m)
    {
      cmat.MultiplyTN(2.0,dSdrho[nr_grf_proc+k][m],drhodC[m],1.0);
      cmat.MultiplyTN(2.0,dSdlamb[nr_grf_proc+k][m],dlambdC[m],1.0);
    }
  }
  // active fiber evaluation
  for(unsigned k=potsumfiberpas_.size();k<(potsumfiberpas_.size()+potsumfiberact_.size());++k)
  {
    // build remodel deformation gradient
    FrM.Update(cur_lambda_r_[k][gp]*(1./params_->G_),AM_[k],0.0);
    FrM.Update(1./std::sqrt(cur_lambda_r_[k][gp]*1./params_->G_),AM_orth_[k],1.0);

    // inverse remodel deformation gradient
    iFrM.Invert(FrM);

    // Fg^-1 * Fr^-1
    iFgiFrM.MultiplyNN(1.0,iFgM,iFrM,0.0);

    // dCedC
    EvaldCedC(dCedC,iFgiFrM,0.5);
    dCedC_strain.MultiplyNN(1.0,strainconv,dCedC,0.0);

    // elastic right Cauchy Green tensor
    FeM.MultiplyNN(1.0,*defgrd,iFgiFrM,0.0);
    CeM.MultiplyTN(1.0,FeM,FeM,0.0);

    // get derivatives of strain energy function w.r.t. the fourth invariant
    potsumfiberact_[k-potsumfiberpas_.size()]->GetDerivativesAniso(dPIe,ddPIIe,dddPIIIe,CeM,eleGID);

    switch(growthtype)
    {
    case 1:
      // dFg^-1/drho
      diFgdrhoM.Update(-(1./(v*v))*(1./density),AgM,0.0);
      break;
    case 0:
      // dFg^-1/drho
      diFgdrhoM.Update(-1./3.*std::pow(v,-4./3.)*(1./density),AcirM,0.0);
      diFgdrhoM.Update(-1./3.*std::pow(v,-4./3.)*(1./density),AradM,1.0);
      diFgdrhoM.Update(-1./3.*std::pow(v,-4./3.)*(1./density),AaxM,1.0);
     break;
    default:
      dserror("growthtype has to be either 1: anisotropic growth or 0: isotropic growth");
      break;
    }

    // dCedrho (lambda_r const)
    diFgdrhoiFrM.MultiplyNN(1.0,diFgdrhoM,iFrM,0.0);
    tmp3x3_1.MultiplyTN(1.0,diFgdrhoiFrM,CM,0.0);
    dCedrhoM.MultiplyNN(1.0,tmp3x3_1,iFgiFrM,0.0);
    tmp3x3_1.MultiplyNN(1.0,CM,diFgdrhoiFrM,0.0);
    dCedrhoM.MultiplyTN(1.0,iFgiFrM,tmp3x3_1,1.0);
    MatrixtoStressVoigtNotationVector(dCedrhoM,dCedrhov);

    // d(dPIe)/drho = ddPII * dI4e/drho (lambda_r const)
    dI4edrho.MultiplyTN(1.0,A_strain_[k],dCedrhov,0.0);
    ddPIedrho = ddPIIe(0) * dI4edrho(0);

    // Evaluate active fiber stress (active contribution)
    potsumfiberact_[k-potsumfiberpas_.size()]->EvaluateActiveStressCmatAniso(*defgrd,cmatactive,stressactive,eleGID);


    // Update stress
    pk2ev.Update(2.*dPIe(0),Av_[k],0.0);      // evaluate elastic second Piola Kirchhoff stress
    stressuw.MultiplyTN(1.0,pk2ev,dCedC_strain,0.0);
    stressuw.UpdateT(1.0,stressactive,1.0);

    stress.UpdateT(cur_rho_col_[k][gp],stressuw,1.0);


    // Evaluate derivation of the second Piola Kirchhoff stress w.r.t. rho
    dSdrho[nr_grf_proc+k][nr_grf_proc+k].Update(1.0,stressuw,0.0);

    dSdrho[nr_grf_proc+k][nr_grf_proc+k].MultiplyTN(2.0*cur_rho_col_[k][gp]*ddPIedrho,A_strain_[k],dCedC,1.0);

    tmp3x3_1.UpdateT(1.0,diFgdrhoiFrM,0.0);
    tmp3x3_2.UpdateT(1.0,iFgiFrM,0.0);
    dCedCdrho.Clear();
    AddtoMatrixHalfSymProd(0.5,tmp3x3_1,tmp3x3_2,dCedCdrho);
    AddtoMatrixHalfSymProd(0.5,tmp3x3_2,tmp3x3_1,dCedCdrho);
    dSdrho[nr_grf_proc+k][nr_grf_proc+k].MultiplyTN(2.0*cur_rho_col_[k][gp]*dPIe(0),A9x1_[k],dCedCdrho,1.0);

    for(int m=0;m<nr_grf_tot;++m)
      if(m != nr_grf_proc+(int)k)
      {
        dSdrho[nr_grf_proc+k][m].Update(1.0,dSdrho[nr_grf_proc+k][nr_grf_proc+k],0.0);
        dSdrho[nr_grf_proc+k][m].Update(-1.0,stressuw,1.0);
      }



    // dFr^-1/dlambda_r
    diFrdlambM.Update(-std::pow(cur_lambda_r_[k][gp],-2.)*params_->G_,AM_[k],0.0);
    diFrdlambM.Update(0.5*std::pow(cur_lambda_r_[k][gp],-0.5)*std::sqrt(1./params_->G_),AM_orth_[k],1.0);

    // dCe/dlambda_r
    tmp3x3_1.MultiplyNN(1.0,iFgM,diFrdlambM,0.0);
    tmp3x3_2.MultiplyTN(1.0,tmp3x3_1,CM,0.0);
    dCedlambM.MultiplyNN(1.0,tmp3x3_2,iFgiFrM,0.0);

    tmp3x3_2.MultiplyTN(1.0,iFgiFrM,CM,0.0);
    dCedlambM.MultiplyNN(1.0,tmp3x3_2,tmp3x3_1,1.0);
    MatrixtoStressVoigtNotationVector(dCedlambM,dCedlambv);

    // dPI/dlambda_r
    ddPIedlamb.MultiplyTN(ddPIIe(0),A_strain_[k],dCedlambv,0.0);


    // Evaluate derivation of the second Piola Kirchhoff stress w.r.t. lambda_r
    dSdlamb[nr_grf_proc+k][nr_grf_proc+k].MultiplyTN(2.0*cur_rho_col_[k][gp]*ddPIedlamb(0),A_strain_[k],dCedC,0.0);

    tmp3x3_1.MultiplyTT(1.0,diFrdlambM,iFgM,0.0);
    tmp3x3_2.UpdateT(1.0,iFgiFrM,0.0);
    dCedCdlamb.Clear();
    AddtoMatrixHalfSymProd(0.5,tmp3x3_1,tmp3x3_2,dCedCdlamb);
    AddtoMatrixHalfSymProd(0.5,tmp3x3_2,tmp3x3_1,dCedCdlamb);
    dSdlamb[nr_grf_proc+k][nr_grf_proc+k].MultiplyTN(2.0*cur_rho_col_[k][gp]*dPIe(0),A9x1_[k],dCedCdlamb,1.0);



    // update elasticity tensor
    cmatelastic.MultiplyNT(4.*ddPIIe(0),Av_[k],Av_[k],0.0);
    tmp6x6.MultiplyTN(1.0,dCedC_strain,cmatelastic,0.0);

    cmat.MultiplyNN(cur_rho_col_[k][gp],tmp6x6,dCedC_strain,1.0);

    // active fiber contribution
    cmat.Update(cur_rho_col_[k][gp],cmatactive,1.0);

    for(int m=0;m<nr_grf_tot;++m)
    {
      cmat.MultiplyTN(2.0,dSdrho[nr_grf_proc+k][m],drhodC[m],1.0);
      cmat.MultiplyTN(2.0,dSdlamb[nr_grf_proc+k][m],dlambdC[m],1.0);
    }
  }

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::RemodelFiber::EvaluateGrowthAndRemodelingExpl(const LINALG::Matrix<3,3> defgrd,
                                                                 const double dt,
                                                                 LINALG::Matrix<3,3> iFgM,
                                                                 const int gp,
                                                                 const int eleGID
  )
{
  // some variables
  std::vector<LINALG::Matrix<3,1> > fibervecs;
  double sigf = 0.0;
  double stress_f_act = 0.0;
  double lamb_r_dot = 0.0;
  LINALG::Matrix<2,1> I_pseudo(true);
  LINALG::Matrix<3,3> FrM(true);
  LINALG::Matrix<3,3> iFrM(true);
  LINALG::Matrix<3,3> iFgiFrM(true);
  LINALG::Matrix<3,3> FeM(true);
  LINALG::Matrix<3,3> CeM(true);
  LINALG::Matrix<2,1> dPIe(true);
  LINALG::Matrix<3,1> ddPIIe(true);
  LINALG::Matrix<4,1> dddPIIIe(true);
  LINALG::Matrix<3,1> Fa(true);
  LINALG::Matrix<6,1> stress_act(true);
  LINALG::Matrix<6,6> cmat_act(true);

  // passive fiber evaluation
  for(unsigned k=0;k<potsumfiberpas_.size();++k)
  {
    // Get fiberdirection
    potsumfiberpas_[k]->GetFiberVecs(fibervecs);

    // build remodel deformation gradient
    FrM.Update(last_lambda_r_[k][gp]*(1./params_->G_),AM_[k],0.0);
    FrM.Update(1./std::sqrt(last_lambda_r_[k][gp]*1./params_->G_),AM_orth_[k],1.0);

    // inverse remodel deformation gradient
    iFrM.Invert(FrM);

    // Fg^-1 * Fr^-1
    iFgiFrM.MultiplyNN(1.0,iFgM,iFrM,0.0);

    // elastic right Cauchy Green tensor
    FeM.MultiplyNN(1.0,defgrd,iFgiFrM,0.0);
    CeM.MultiplyTN(1.0,FeM,FeM,0.0);

    // get derivatives of strain energy function w.r.t. the fourth invariant
    potsumfiberpas_[k]->GetDerivativesAniso(dPIe,ddPIIe,dddPIIIe,CeM,eleGID);

    // Evaluate fiber Cauchy stress in updated fiber direction
    // pseudo-invariant I4
    I_pseudo(0) = CeM.Dot(AM_[k]);

    // fiber Cauchy stress (passive contribution)
    sigf = 2.*I_pseudo(0)*dPIe(0);


    // update current reference mass density
    cur_rho_col_[k][gp] = last_rho_col_[k][gp] + last_rho_col_[k][gp]*dt*params_->k_growth_*((sigf-sigmapre_[k])/sigmapre_[k]);

    // update inelastic remodeling stretch
    lamb_r_dot = ((((cur_rho_col_[k][gp]-last_rho_col_[k][gp])/dt)/last_rho_col_[k][gp]+(1./params_->tdecay_))*(sigf-sigmapre_[k])*last_lambda_r_[k][gp])/
        (4.0*(dPIe(0)+I_pseudo(0)*ddPIIe(0))*I_pseudo(0));

    cur_lambda_r_[k][gp] = last_lambda_r_[k][gp] + dt*lamb_r_dot;
  }
  // active fiber evaluation
  for(unsigned k=potsumfiberpas_.size();k<(potsumfiberpas_.size()+potsumfiberact_.size());++k)
  {
    // Get fiberdirection
    potsumfiberact_[k-potsumfiberpas_.size()]->GetFiberVecs(fibervecs);

    // build remodel deformation gradient
    FrM.Update(last_lambda_r_[k][gp]*(1./params_->G_),AM_[k],0.0);
    FrM.Update(1./std::sqrt(last_lambda_r_[k][gp]*1./params_->G_),AM_orth_[k],1.0);

    // inverse remodel deformation gradient
    iFrM.Invert(FrM);

    // Fg^-1 * Fr^-1
    iFgiFrM.MultiplyNN(1.0,iFgM,iFrM,0.0);

    // elastic right Cauchy Green tensor
    FeM.MultiplyNN(1.0,defgrd,iFgiFrM,0.0);
    CeM.MultiplyTN(1.0,FeM,FeM,0.0);

    // get derivatives of strain energy function w.r.t. the fourth invariant
    potsumfiberact_[k-potsumfiberpas_.size()]->GetDerivativesAniso(dPIe,ddPIIe,dddPIIIe,CeM,eleGID);

    // Evaluate fiber Cauchy stress in updated fiber direction
    // pseudo-invariant I4
    I_pseudo(0) = CeM.Dot(AM_[k]);

    // fiber Cauchy stress (passive contribution)
    sigf = 2.*I_pseudo(0)*dPIe(0);

    // fiber Cauchy stress (active contribution)
    // dsigf/dC (active contribution)
    potsumfiberact_[k-potsumfiberpas_.size()]->EvaluateActiveStressCmatAniso(defgrd,cmat_act,stress_act,eleGID);
    stress_f_act = stress_act.Dot(A_strain_[k]);
    Fa.MultiplyNN(1.0,defgrd,fibervecs[k],0.0);

    sigf += Fa.Norm2()*Fa.Norm2()*stress_f_act;

    // update current reference mass density
    cur_rho_col_[k][gp] = last_rho_col_[k][gp] + last_rho_col_[k][gp]*dt*params_->k_growth_*((sigf-sigmapre_[k])/sigmapre_[k]);

    // update inelastic remodeling stretch
    lamb_r_dot = ((((cur_rho_col_[k][gp]-last_rho_col_[k][gp])/dt)/last_rho_col_[k][gp]+(1./params_->tdecay_))*(sigf-sigmapre_[k])*last_lambda_r_[k][gp])/
        (4.0*(dPIe(0)+I_pseudo(0)*ddPIIe(0))*I_pseudo(0));

    cur_lambda_r_[k][gp] = last_lambda_r_[k][gp] + dt*lamb_r_dot;
  }

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::RemodelFiber::EvaluateGrowthAndRemodelingExplMembrane(const LINALG::Matrix<3,3>& defgrd_loc,
                                                                         const LINALG::Matrix<3,3>& Q_trafo,
                                                                         const double v,
                                                                         const double dt,
                                                                         const LINALG::Matrix<3,3> iFgM,
                                                                         const int gp,
                                                                         const int eleGID)
{
  // some variables
  std::vector<LINALG::Matrix<3,1> > fibervecs;
  double sigf = 0.0;
  double stress_f_act = 0.0;
  double lamb_r_dot = 0.0;
  LINALG::Matrix<3,3> tmp(true);
  LINALG::Matrix<3,3> AM_loc(true);
  LINALG::Matrix<3,3> IAM_loc(true);
  LINALG::Matrix<3,3> CM_glob(true);
  LINALG::Matrix<3,3> C_loc(true);
  LINALG::Matrix<3,3> FrM_loc(true);
  LINALG::Matrix<3,3> iFrM_loc(true);
  LINALG::Matrix<3,3> iFg_loc(true);
  LINALG::Matrix<3,3> iFgiFrM_loc(true);
  LINALG::Matrix<3,3> FrFgM_loc(true);
  LINALG::Matrix<3,3> CeM_loc(true);
  LINALG::Matrix<3,3> CeM_glob(true);
  LINALG::Matrix<2,1> pseudo_inv(true);
  LINALG::Matrix<2,1> dPIe(true);
  LINALG::Matrix<3,1> ddPIIe(true);
  LINALG::Matrix<4,1> dddPIIIe(true);
  LINALG::Matrix<6,1> stress_act(true);
  LINALG::Matrix<6,6> cmat_act(true);

  // passive fiber evaluation
  for(unsigned k=0;k<potsumfiberpas_.size();++k)
  {
    // Get fiberdirection
    potsumfiberpas_[k]->GetFiberVecs(fibervecs);

    // fibervector in orthonormal frame on membrane surface
    LINALG::Matrix<3,1> fibervector(true);
    fibervector.MultiplyTN(1.0,Q_trafo,fibervecs[k],0.0);

    // Remark: the fibervector in local coordinates is in the membrane surface tangent plane, therefore fibervector(2)=0.0
    // structural tensor in matrix notation
    AM_loc.MultiplyNT(1.0,fibervector,fibervector,0.0);

    // tensor I-AM
    IAM_loc.Clear();
    IAM_loc(0,0) = IAM_loc(1,1) = IAM_loc(2,2) = 1.0;
    IAM_loc.Update(-1.0,AM_loc,1.0);

    // build remodel deformation gradient
    FrM_loc.Update(last_lambda_r_[k][gp]*(1./params_->G_),AM_loc,0.0);
    FrM_loc.Update(1./std::sqrt(last_lambda_r_[k][gp]*1./params_->G_),IAM_loc,1.0);

    // inverse remodel deformation gradient
    iFrM_loc.Invert(FrM_loc);

    tmp.MultiplyTN(1.0,Q_trafo,iFgM,0.0);
    iFg_loc.MultiplyNN(1.0,tmp,Q_trafo,0.0);

    // Fg^-1 * Fr^-1
    iFgiFrM_loc.MultiplyNN(1.0,iFg_loc,iFrM_loc,0.0);

    // elastic right Cauchy-Green in matrix notation
    C_loc.MultiplyTN(1.0,defgrd_loc,defgrd_loc,0.0);
    tmp.MultiplyNN(1.0,C_loc,iFgiFrM_loc,0.0);
    CeM_loc.MultiplyTN(1.0,iFgiFrM_loc,tmp,0.0);

    // impose incompressibility for component in thickness direction of elastic right Cauchy-Green
    CeM_loc(2,2) = 1.0/(CeM_loc(0,0)*CeM_loc(1,1)-CeM_loc(0,1)*CeM_loc(1,0));

    // elastic right Cauchy Green in global coordinate system
    tmp.MultiplyNN(1.0,Q_trafo,CeM_loc,0.0);
    CeM_glob.MultiplyNT(1.0,tmp,Q_trafo,0.0);
    pseudo_inv(0) = CeM_loc.Dot(AM_loc);

    // get derivatives of strain energy function w.r.t. the fourth invariant
    potsumfiberpas_[k]->GetDerivativesAniso(dPIe,ddPIIe,dddPIIIe,CeM_glob,eleGID);

    // fiber Cauchy stress (passive contribution)
    sigf = 2.*pseudo_inv(0)*dPIe(0);

    // update current reference mass density
    cur_rho_col_[k][gp] = last_rho_col_[k][gp] + last_rho_col_[k][gp]*dt*params_->k_growth_*((sigf-sigmapre_[k])/sigmapre_[k]);

    // update inelastic remodeling stretch
    lamb_r_dot = ((((cur_rho_col_[k][gp]-last_rho_col_[k][gp])/dt)/last_rho_col_[k][gp]+(1./params_->tdecay_))*(sigf-sigmapre_[k])*last_lambda_r_[k][gp])/
        (4.0*(dPIe(0)+pseudo_inv(0)*ddPIIe(0))*pseudo_inv(0));

    cur_lambda_r_[k][gp] = last_lambda_r_[k][gp] + dt*lamb_r_dot;

    stress_[k][gp] = sigf;
  }
  // active fiber evaluation
  for(unsigned k=potsumfiberpas_.size();k<(potsumfiberpas_.size()+potsumfiberact_.size());++k)
  {
    // Get fiberdirection
    potsumfiberact_[k-potsumfiberpas_.size()]->GetFiberVecs(fibervecs);

    // fibervector in orthonormal frame on membrane surface
    LINALG::Matrix<3,1> fibervector(true);
    fibervector.MultiplyTN(1.0,Q_trafo,fibervecs[k],0.0);

    // Remark: the fibervector in local coordinates is in the membrane surface tangent plane, therefore fibervector(2)=0.0
    // structural tensor in matrix notation
    AM_loc.MultiplyNT(1.0,fibervector,fibervector,0.0);

    // tensor I-AM
    IAM_loc.Clear();
    IAM_loc(0,0) = IAM_loc(1,1) = IAM_loc(2,2) = 1.0;
    IAM_loc.Update(-1.0,AM_loc,1.0);

    // build remodel deformation gradient
    FrM_loc.Update(last_lambda_r_[k][gp]*(1./params_->G_),AM_loc,0.0);
    FrM_loc.Update(1./std::sqrt(last_lambda_r_[k][gp]*1./params_->G_),IAM_loc,1.0);

    // inverse remodel deformation gradient
    iFrM_loc.Invert(FrM_loc);

    tmp.MultiplyTN(1.0,Q_trafo,iFgM,0.0);
    iFg_loc.MultiplyNN(1.0,tmp,Q_trafo,0.0);

    // Fg^-1 * Fr^-1
    iFgiFrM_loc.MultiplyNN(1.0,iFg_loc,iFrM_loc,0.0);

    // elastic right Cauchy-Green in matrix notation
    C_loc.MultiplyTN(1.0,defgrd_loc,defgrd_loc,0.0);
    tmp.MultiplyNN(1.0,C_loc,iFgiFrM_loc,0.0);
    CeM_loc.MultiplyTN(1.0,iFgiFrM_loc,tmp,0.0);

    // impose incompressibility for component in thickness direction of elastic right Cauchy-Green
    CeM_loc(2,2) = 1.0/(CeM_loc(0,0)*CeM_loc(1,1)-CeM_loc(0,1)*CeM_loc(1,0));

    // elastic right Cauchy Green in global coordinate system
    tmp.MultiplyNN(1.0,Q_trafo,CeM_loc,0.0);
    CeM_glob.MultiplyNT(1.0,tmp,Q_trafo,0.0);
    pseudo_inv(0) = CeM_loc.Dot(AM_loc);

    // get derivatives of strain energy function w.r.t. the fourth invariant
    potsumfiberact_[k-potsumfiberpas_.size()]->GetDerivativesAniso(dPIe,ddPIIe,dddPIIIe,CeM_glob,eleGID);

    // fiber Cauchy stress (passive contribution)
    sigf = 2.*pseudo_inv(0)*dPIe(0);

    // fiber Cauchy stress (active contribution)
    // right Cauchy Green in global coordinate system
    FrFgM_loc.Invert(iFgiFrM_loc);
    tmp.MultiplyTN(1.0,FrFgM_loc,CeM_loc,0.0);
    C_loc.MultiplyNN(1.0,tmp,FrFgM_loc,0.0);
    tmp.MultiplyNN(1.0,Q_trafo,C_loc,0.0);
    CM_glob.MultiplyNT(1.0,tmp,Q_trafo,0.0);

    potsumfiberact_[k-potsumfiberpas_.size()]->EvaluateActiveStressCmatAniso(CM_glob,cmat_act,stress_act,eleGID);
    stress_f_act = stress_act.Dot(A_strain_[k]);

    sigf += C_loc.Dot(AM_loc)*stress_f_act;

    // update current reference mass density
    cur_rho_col_[k][gp] = last_rho_col_[k][gp] + last_rho_col_[k][gp]*dt*params_->k_growth_*((sigf-sigmapre_[k])/sigmapre_[k]);

    // update inelastic remodeling stretch
    lamb_r_dot = ((((cur_rho_col_[k][gp]-last_rho_col_[k][gp])/dt)/last_rho_col_[k][gp]+(1./params_->tdecay_))*(sigf-sigmapre_[k])*last_lambda_r_[k][gp])/
        (4.0*(dPIe(0)+pseudo_inv(0)*ddPIIe(0))*pseudo_inv(0));

    cur_lambda_r_[k][gp] = last_lambda_r_[k][gp] + dt*lamb_r_dot;

    stress_[k][gp] = sigf;
  }

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::RemodelFiber::MatrixtoStressVoigtNotationVector(const LINALG::Matrix<3,3>& in,
                                                                   LINALG::Matrix<6,1>& out)
{
  // "stress-like" Voigt notation
  for (int i=0; i<3; i++) out(i) = in(i,i);
  out(3) = 0.5*(in(0,1) + in(1,0));
  out(4) = 0.5*(in(1,2) + in(2,1));
  out(5) = 0.5*(in(0,2) + in(2,0));

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::RemodelFiber::EvaldCedC(LINALG::Matrix<6,6>& out,
                                           const LINALG::Matrix<3,3>& in,
                                           const double fac)
{
  out(0,0) = 2. * fac * in(0,0) * in(0,0);
  out(0,3) = 2. * fac * in(0,0) * in(1,0);
  out(0,5) = 2. * fac * in(0,0) * in(2,0);
  out(0,1) = 2. * fac * in(1,0) * in(1,0);
  out(0,4) = 2. * fac * in(1,0) * in(2,0);
  out(0,2) = 2. * fac * in(2,0) * in(2,0);

  out(3,0) = 2. * fac * in(0,0) * in(0,1);
  out(3,3) = fac * (in(0,0) * in(1,1) + in(0,1) * in(1,0));
  out(3,5) = fac * (in(0,0) * in(2,1) + in(0,1) * in(2,0));
  out(3,1) = 2. * fac * in(1,0) * in(1,1);
  out(3,4) = fac * (in(1,0) * in(2,1) + in(1,1) * in(2,0));
  out(3,2) = 2. * fac * in(2,0) * in(2,1);

  out(5,0) = 2. * fac * in(0,0) * in(0,2);
  out(5,3) = fac * (in(0,0) * in(1,2) + in(0,2) * in(1,0));
  out(5,5) = fac * (in(0,0) * in(2,2) + in(0,2) * in(2,0));
  out(5,1) = 2. * fac * in(1,0) * in(1,2);
  out(5,4) = fac * (in(1,0) * in(2,2) + in(1,2) * in(2,0));
  out(5,2) = 2. * fac * in(2,0) * in(2,2);

  out(1,0) = 2. * fac * in(0,1) * in(0,1);
  out(1,3) = 2. * fac * in(0,1) * in(1,1);
  out(1,5) = 2. * fac * in(0,1) * in(2,1);
  out(1,1) = 2. * fac * in(1,1) * in(1,1);
  out(1,4) = 2. * fac * in(1,1) * in(2,1);
  out(1,2) = 2. * fac * in(2,1) * in(2,1);

  out(4,0) = 2. * fac * in(0,1) * in(0,2);
  out(4,3) = fac * (in(0,1) * in(1,2) + in(0,2) * in(1,1));
  out(4,5) = fac * (in(0,1) * in(2,2) + in(0,2) * in(2,1));
  out(4,1) = 2. * fac * in(1,1) * in(1,2);
  out(4,4) = fac * (in(1,1) * in(2,2) + in(1,2) * in(2,1));
  out(4,2) = 2. * fac * in(2,1) * in(2,2);

  out(2,0) = 2. * fac * in(0,2) * in(0,2);
  out(2,3) = 2. * fac * in(0,2) * in(1,2);
  out(2,5) = 2. * fac * in(0,2) * in(2,2);
  out(2,1) = 2. * fac * in(1,2) * in(1,2);
  out(2,4) = 2. * fac * in(1,2) * in(2,2);
  out(2,2) = 2. * fac * in(2,2) * in(2,2);

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::RemodelFiber::AddtodPK2diFg(LINALG::Matrix<6,9>& out,
                                               LINALG::Matrix<3,3> A,
                                               LINALG::Matrix<3,3> B,
                                               double fac)
{
  out(0,0) += 2 * fac * A(0,0) * B(0,0);
  out(0,3) += 2 * fac * A(0,0) * B(0,1);
  out(0,5) += 2 * fac * A(0,0) * B(0,2);
  out(0,6) += 2 * fac * A(0,1) * B(0,0);
  out(0,1) += 2 * fac * A(0,1) * B(0,1);
  out(0,4) += 2 * fac * A(0,1) * B(0,2);
  out(0,8) += 2 * fac * A(0,2) * B(0,0);
  out(0,7) += 2 * fac * A(0,2) * B(0,1);
  out(0,2) += 2 * fac * A(0,2) * B(0,2);

  out(1,0) += 2 * fac * A(1,0) * B(1,0);
  out(1,3) += 2 * fac * A(1,0) * B(1,1);
  out(1,5) += 2 * fac * A(1,0) * B(1,2);
  out(1,6) += 2 * fac * A(1,1) * B(1,0);
  out(1,1) += 2 * fac * A(1,1) * B(1,1);
  out(1,4) += 2 * fac * A(1,1) * B(1,2);
  out(1,8) += 2 * fac * A(1,2) * B(1,0);
  out(1,7) += 2 * fac * A(1,2) * B(1,1);
  out(1,2) += 2 * fac * A(1,2) * B(1,2);

  out(2,0) += 2 * fac * A(2,0) * B(2,0);
  out(2,3) += 2 * fac * A(2,0) * B(2,1);
  out(2,5) += 2 * fac * A(2,0) * B(2,2);
  out(2,6) += 2 * fac * A(2,1) * B(2,0);
  out(2,1) += 2 * fac * A(2,1) * B(2,1);
  out(2,4) += 2 * fac * A(2,1) * B(2,2);
  out(2,8) += 2 * fac * A(2,2) * B(2,0);
  out(2,7) += 2 * fac * A(2,2) * B(2,1);
  out(2,2) += 2 * fac * A(2,2) * B(2,2);

  out(3,0) += fac * (A(0,0) * B(1,0) + A(1,0) * B(0,0));
  out(3,3) += fac * (A(0,0) * B(1,1) + A(1,0) * B(0,1));
  out(3,5) += fac * (A(0,0) * B(1,2) + A(1,0) * B(0,2));
  out(3,6) += fac * (A(0,1) * B(1,0) + A(1,1) * B(0,0));
  out(3,1) += fac * (A(0,1) * B(1,1) + A(1,1) * B(0,1));
  out(3,4) += fac * (A(0,1) * B(1,2) + A(1,1) * B(0,2));
  out(3,8) += fac * (A(0,2) * B(1,0) + A(1,2) * B(0,0));
  out(3,7) += fac * (A(0,2) * B(1,1) + A(1,2) * B(0,1));
  out(3,2) += fac * (A(0,2) * B(1,2) + A(1,2) * B(0,2));

  out(4,0) += fac * (A(1,0) * B(2,0) + A(2,0) * B(1,0));
  out(4,3) += fac * (A(1,0) * B(2,1) + A(2,0) * B(1,1));
  out(4,5) += fac * (A(1,0) * B(2,2) + A(2,0) * B(1,2));
  out(4,6) += fac * (A(1,1) * B(2,0) + A(2,1) * B(1,0));
  out(4,1) += fac * (A(1,1) * B(2,1) + A(2,1) * B(1,1));
  out(4,4) += fac * (A(1,1) * B(2,2) + A(2,1) * B(1,2));
  out(4,8) += fac * (A(1,2) * B(2,0) + A(2,2) * B(1,0));
  out(4,7) += fac * (A(1,2) * B(2,1) + A(2,2) * B(1,1));
  out(4,2) += fac * (A(1,2) * B(2,2) + A(2,2) * B(1,2));

  out(5,0) += fac * (A(0,0) * B(2,0) + A(2,0) * B(0,0));
  out(5,3) += fac * (A(0,0) * B(2,1) + A(2,0) * B(0,1));
  out(5,5) += fac * (A(0,0) * B(2,2) + A(2,0) * B(0,2));
  out(5,6) += fac * (A(0,1) * B(2,0) + A(2,1) * B(0,0));
  out(5,1) += fac * (A(0,1) * B(2,1) + A(2,1) * B(0,1));
  out(5,4) += fac * (A(0,1) * B(2,2) + A(2,1) * B(0,2));
  out(5,8) += fac * (A(0,2) * B(2,0) + A(2,2) * B(0,0));
  out(5,7) += fac * (A(0,2) * B(2,1) + A(2,2) * B(0,1));
  out(5,2) += fac * (A(0,2) * B(2,2) + A(2,2) * B(0,2));

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::RemodelFiber::AddtoMatrixHalfSymProd(const double& fac,
                                                        const LINALG::Matrix<3,3>& A,
                                                        const LINALG::Matrix<3,3>& B,
                                                        LINALG::Matrix<9,6>& out)
{
  out(0,0) += 2. * fac * A(0,0) * B(0,0);
  out(0,3) += fac * (A(0,0) * B(0,1) + A(0,1) * B(0,0));
  out(0,5) += fac * (A(0,0) * B(0,2) + A(0,2) * B(0,0));
  out(0,1) += 2. * fac * A(0,1) * B(0,1);
  out(0,4) += fac * (A(0,1) * B(0,2) + A(0,2) * B(0,1));
  out(0,2) += 2. * fac * A(0,2) * B(0,2);

  out(3,0) += 2. * fac * A(0,0) * B(1,0);
  out(3,3) += fac * (A(0,0) * B(1,1) + A(0,1) * B(1,0));
  out(3,5) += fac * (A(0,0) * B(1,2) + A(0,2) * B(1,0));
  out(3,1) += 2. * fac * A(0,1) * B(1,1);
  out(3,4) += fac * (A(0,1) * B(1,2) + A(0,2) * B(1,1));
  out(3,2) += 2. * fac * A(0,2) * B(1,2);

  out(5,0) += 2. * fac * A(0,0) * B(2,0);
  out(5,3) += fac * (A(0,0) * B(2,1) + A(0,1) * B(2,0));
  out(5,5) += fac * (A(0,0) * B(2,2) + A(0,2) * B(2,0));
  out(5,1) += 2. * fac * A(0,1) * B(2,1);
  out(5,4) += fac * (A(0,1) * B(2,2) + A(0,2) * B(2,1));
  out(5,2) += 2. * fac * A(0,2) * B(2,2);

  out(6,0) += 2. * fac * A(1,0) * B(0,0);
  out(6,3) += fac * (A(1,0) * B(0,1) + A(1,1) * B(0,0));
  out(6,5) += fac * (A(1,0) * B(0,2) + A(1,2) * B(0,0));
  out(6,1) += 2. * fac * A(1,1) * B(0,1);
  out(6,4) += fac * (A(1,1) * B(0,2) + A(1,2) * B(0,1));
  out(6,2) += 2. * fac * A(1,2) * B(0,2);

  out(1,0) += 2. * fac * A(1,0) * B(1,0);
  out(1,3) += fac * (A(1,0) * B(1,1) + A(1,1) * B(1,0));
  out(1,5) += fac * (A(1,0) * B(1,2) + A(1,2) * B(1,0));
  out(1,1) += 2. * fac * A(1,1) * B(1,1);
  out(1,4) += fac * (A(1,1) * B(1,2) + A(1,2) * B(1,1));
  out(1,2) += 2. * fac * A(1,2) * B(1,2);

  out(4,0) += 2. * fac * A(1,0) * B(2,0);
  out(4,3) += fac * (A(1,0) * B(2,1) + A(1,1) * B(2,0));
  out(4,5) += fac * (A(1,0) * B(2,2) + A(1,2) * B(2,0));
  out(4,1) += 2. * fac * A(1,1) * B(2,1);
  out(4,4) += fac * (A(1,1) * B(2,2) + A(1,2) * B(2,1));
  out(4,2) += 2. * fac * A(1,2) * B(2,2);

  out(8,0) += 2. * fac * A(2,0) * B(0,0);
  out(8,3) += fac * (A(2,0) * B(0,1) + A(2,1) * B(0,0));
  out(8,5) += fac * (A(2,0) * B(0,2) + A(2,2) * B(0,0));
  out(8,1) += 2. * fac * A(2,1) * B(0,1);
  out(8,4) += fac * (A(2,1) * B(0,2) + A(2,2) * B(0,1));
  out(8,2) += 2. * fac * A(2,2) * B(0,2);

  out(7,0) += 2. * fac * A(2,0) * B(1,0);
  out(7,3) += fac * (A(2,0) * B(1,1) + A(2,1) * B(1,0));
  out(7,5) += fac * (A(2,0) * B(1,2) + A(2,2) * B(1,0));
  out(7,1) += 2. * fac * A(2,1) * B(1,1);
  out(7,4) += fac * (A(2,1) * B(1,2) + A(2,2) * B(1,1));
  out(7,2) += 2. * fac * A(2,2) * B(1,2);

  out(2,0) += 2. * fac * A(2,0) * B(2,0);
  out(2,3) += fac * (A(2,0) * B(2,1) + A(2,1) * B(2,0));
  out(2,5) += fac * (A(2,0) * B(2,2) + A(2,2) * B(2,0));
  out(2,1) += 2. * fac * A(2,1) * B(2,1);
  out(2,4) += fac * (A(2,1) * B(2,2) + A(2,2) * B(2,1));
  out(2,2) += 2. * fac * A(2,2) * B(2,2);

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::RemodelFiber::EvaldYdC(LINALG::Matrix<6,6>& out,
                                          const LINALG::Matrix<3,3>& in1,
                                          const LINALG::Matrix<3,3>& in2,
                                          const double fac)
{
  out(0,0) = 4 * fac * in1(0,0) * in2(0,0);
  out(0,3) = fac * (2 * in1(0,0) * in2(1,0) + 2 * in1(1,0) * in2(0,0));
  out(0,5) = fac * (2 * in1(0,0) * in2(2,0) + 2 * in1(2,0) * in2(0,0));
  out(0,1) = 4 * fac * in1(1,0) * in2(1,0);
  out(0,4) = fac * (2 * in1(1,0) * in2(2,0) + 2 * in1(2,0) * in2(1,0));
  out(0,2) = 4 * fac * in1(2,0) * in2(2,0);

  out(3,0) = fac * (2 * in1(0,0) * in2(0,1) + 2 * in1(0,1) * in2(0,0));
  out(3,3) = fac * (in1(0,0) * in2(1,1) + in1(0,1) * in2(1,0) + in1(1,0) * in2(0,1) + in1(1,1) * in2(0,0));
  out(3,5) = fac * (in1(0,0) * in2(2,1) + in1(0,1) * in2(2,0) + in1(2,0) * in2(0,1) + in1(2,1) * in2(0,0));
  out(3,1) = fac * (2 * in1(1,0) * in2(1,1) + 2 * in1(1,1) * in2(1,0));
  out(3,4) = fac * (in1(1,0) * in2(2,1) + in1(1,1) * in2(2,0) + in1(2,0) * in2(1,1) + in1(2,1) * in2(1,0));
  out(3,2) = fac * (2 * in1(2,0) * in2(2,1) + 2 * in1(2,1) * in2(2,0));

  out(5,0) = fac * (2 * in1(0,0) * in2(0,2) + 2 * in1(0,2) * in2(0,0));
  out(5,3) = fac * (in1(0,0) * in2(1,2) + in1(0,2) * in2(1,0) + in1(1,0) * in2(0,2) + in1(1,2) * in2(0,0));
  out(5,5) = fac * (in1(0,0) * in2(2,2) + in1(0,2) * in2(2,0) + in1(2,0) * in2(0,2) + in1(2,2) * in2(0,0));
  out(5,1) = fac * (2 * in1(1,0) * in2(1,2) + 2 * in1(1,2) * in2(1,0));
  out(5,4) = fac * (in1(1,0) * in2(2,2) + in1(1,2) * in2(2,0) + in1(2,0) * in2(1,2) + in1(2,2) * in2(1,0));
  out(5,2) = fac * (2 * in1(2,0) * in2(2,2) + 2 * in1(2,2) * in2(2,0));

  out(1,0) = 4 * fac * in1(0,1) * in2(0,1);
  out(1,3) = fac * (2 * in1(0,1) * in2(1,1) + 2 * in1(1,1) * in2(0,1));
  out(1,5) = fac * (2 * in1(0,1) * in2(2,1) + 2 * in1(2,1) * in2(0,1));
  out(1,1) = 4 * fac * in1(1,1) * in2(1,1);
  out(1,4) = fac * (2 * in1(1,1) * in2(2,1) + 2 * in1(2,1) * in2(1,1));
  out(1,2) = 4 * fac * in1(2,1) * in2(2,1);

  out(4,0) = fac * (2 * in1(0,1) * in2(0,2) + 2 * in1(0,2) * in2(0,1));
  out(4,3) = fac * (in1(0,1) * in2(1,2) + in1(0,2) * in2(1,1) + in1(1,1) * in2(0,2) + in1(1,2) * in2(0,1));
  out(4,5) = fac * (in1(0,1) * in2(2,2) + in1(0,2) * in2(2,1) + in1(2,1) * in2(0,2) + in1(2,2) * in2(0,1));
  out(4,1) = fac * (2 * in1(1,1) * in2(1,2) + 2 * in1(1,2) * in2(1,1));
  out(4,4) = fac * (in1(1,1) * in2(2,2) + in1(1,2) * in2(2,1) + in1(2,1) * in2(1,2) + in1(2,2) * in2(1,1));
  out(4,2) = fac * (2 * in1(2,1) * in2(2,2) + 2 * in1(2,2) * in2(2,1));

  out(2,0) = 4 * fac * in1(0,2) * in2(0,2);
  out(2,3) = fac * (2 * in1(0,2) * in2(1,2) + 2 * in1(1,2) * in2(0,2));
  out(2,5) = fac * (2 * in1(0,2) * in2(2,2) + 2 * in1(2,2) * in2(0,2));
  out(2,1) = 4 * fac * in1(1,2) * in2(1,2);
  out(2,4) = fac * (2 * in1(1,2) * in2(2,2) + 2 * in1(2,2) * in2(1,2));
  out(2,2) = 4 * fac * in1(2,2) * in2(2,2);

  return;
}


/*---------------------------------------------------------------------*
 | return names of visualization data (public)                         |
 *---------------------------------------------------------------------*/
void MAT::ELASTIC::RemodelFiber::VisNames(std::map<std::string,int>& names, unsigned int p)
{
  std::string inelastic_defgrd = "lambda_r";
  std::string result_inelastic_defgrad;

  for (unsigned int k=0; k<(potsumfiberpas_.size()+potsumfiberact_.size()); ++k)
  {
    std::stringstream sstm;
    sstm << inelastic_defgrd <<"_" << p <<"_" << k;
    result_inelastic_defgrad = sstm.str();

    names[result_inelastic_defgrad] = 1;
  }


  std::string fiber_cauchy_stress = "fiber_cauchy_stress";
  std::string result_fiber_cauchy_stress;

  for (unsigned int k=0; k<(potsumfiberpas_.size()+potsumfiberact_.size()); ++k)
  {
    std::stringstream sstm;
    sstm << fiber_cauchy_stress <<"_" << p <<"_" << k;
    result_fiber_cauchy_stress = sstm.str();

    names[result_fiber_cauchy_stress] = 1;
  }


  std::string cur_rho_col = "cur_rho_col";
  std::string result_cur_rho_col;

  for (unsigned int k=0; k<(potsumfiberpas_.size()+potsumfiberact_.size()); ++k)
  {
    std::stringstream sstm;
    sstm << cur_rho_col <<"_" << p <<"_" << k;
    result_cur_rho_col = sstm.str();

    names[result_cur_rho_col] = 1;
  }

  std::string result_sum_normalized_rho_col;
  if(p == 0)
  {
   result_sum_normalized_rho_col = "sum_normalized_rho_col_0_0";
   names[result_sum_normalized_rho_col] = 1;
  }

}  // VisNames()


/*---------------------------------------------------------------------*
 | return visualization data (public)                                  |
 *---------------------------------------------------------------------*/
bool MAT::ELASTIC::RemodelFiber::VisData(
  const std::string& name,
  std::vector<double>& data,
  int numgp,
  int eleID
)
{
  if ((name == "lambda_r_0_0") || (name == "lambda_r_1_0"))
  {
    if (data.size()!= 1) dserror("size mismatch");
    for(unsigned gp=0; gp<last_lambda_r_[0].size(); ++gp)
    {
      data[0] += last_lambda_r_[0][gp];
    }
    data[0] = data[0]/last_lambda_r_[0].size();

    return true;
  }
  if ((name == "lambda_r_0_1") || (name == "lambda_r_1_1"))
  {
    if (data.size()!= 1) dserror("size mismatch");
    for(unsigned gp=0; gp<last_lambda_r_[1].size(); ++gp)
    {
      data[0] += last_lambda_r_[1][gp];
    }
    data[0] = data[0]/last_lambda_r_[1].size();

    return true;
  }
  if ((name == "lambda_r_0_2") || (name == "lambda_r_1_2"))
  {
    if (data.size()!= 1) dserror("size mismatch");
    for(unsigned gp=0; gp<last_lambda_r_[2].size(); ++gp)
    {
      data[0] += last_lambda_r_[2][gp];
    }
    data[0] = data[0]/last_lambda_r_[2].size();

    return true;
  }
  if ((name == "lambda_r_0_3") || (name == "lambda_r_1_3"))
  {
    if (data.size()!= 1) dserror("size mismatch");
    for(unsigned gp=0; gp<last_lambda_r_[3].size(); ++gp)
    {
      data[0] += last_lambda_r_[3][gp];
    }
    data[0] = data[0]/last_lambda_r_[3].size();

    return true;
  }



  if ((name == "fiber_cauchy_stress_0_0") || (name == "fiber_cauchy_stress_1_0"))
  {
    if (data.size()!= 1) dserror("size mismatch");
    for(unsigned gp=0; gp<stress_[0].size(); ++gp)
    {
      data[0] += stress_[0][gp];
    }
    data[0] = data[0]/stress_[0].size();

    return true;
  }
  if ((name == "fiber_cauchy_stress_0_1") || (name == "fiber_cauchy_stress_1_1"))
  {
    if (data.size()!= 1) dserror("size mismatch");
    for(unsigned gp=0; gp<stress_[1].size(); ++gp)
    {
      data[0] += stress_[1][gp];
    }
    data[0] = data[0]/stress_[1].size();

    return true;
  }
  if ((name == "fiber_cauchy_stress_0_2") || (name == "fiber_cauchy_stress_1_2"))
  {
    if (data.size()!= 1) dserror("size mismatch");
    for(unsigned gp=0; gp<stress_[2].size(); ++gp)
    {
      data[0] += stress_[2][gp];
    }
    data[0] = data[0]/stress_[2].size();

    return true;
  }
  if ((name == "fiber_cauchy_stress_0_3") || (name == "fiber_cauchy_stress_1_3"))
  {
    if (data.size()!= 1) dserror("size mismatch");
    for(unsigned gp=0; gp<stress_[3].size(); ++gp)
    {
      data[0] += stress_[3][gp];
    }
    data[0] = data[0]/stress_[3].size();

    return true;
  }



if(name == "sum_normalized_rho_col_0_0")
{
  if (data.size()!= 1) dserror("size mismatch");
  for(unsigned i=0;i<cur_rho_col_[0].size();++i)
  {
    data[0] += (cur_rho_col_[0][i]+cur_rho_col_[1][i]+cur_rho_col_[2][i]+cur_rho_col_[3][i])/
       (init_rho_col_[0][i]+init_rho_col_[1][i]+init_rho_col_[2][i]+init_rho_col_[3][i]);
  }
  data[0] = data[0]/cur_rho_col_[0].size();

  return true;
}




if((name == "cur_rho_col_0_0") || (name == "cur_rho_col_1_0"))
{
  if (data.size()!= 1) dserror("size mismatch");
  for(unsigned i=0;i<cur_rho_col_[0].size();++i)
  {
    data[0] += cur_rho_col_[0][i];
  }
  data[0] = data[0]/cur_rho_col_[0].size();


  return true;
}
if((name == "cur_rho_col_0_1") || (name == "cur_rho_col_1_1"))
{
  if (data.size()!= 1) dserror("size mismatch");
  for(unsigned i=0;i<cur_rho_col_[0].size();++i)
  {
    data[0] += cur_rho_col_[1][i];
  }
  data[0] = data[0]/cur_rho_col_[1].size();


  return true;
}
if((name == "cur_rho_col_0_2") || (name == "cur_rho_col_1_2"))
{
  if (data.size()!= 1) dserror("size mismatch");
  for(unsigned i=0;i<cur_rho_col_[0].size();++i)
  {
    data[0] += cur_rho_col_[2][i];
  }
  data[0] = data[0]/cur_rho_col_[2].size();


  return true;
}
if((name == "cur_rho_col_0_3") || (name == "cur_rho_col_1_3"))
{
  if (data.size()!= 1) dserror("size mismatch");
  for(unsigned i=0;i<cur_rho_col_[0].size();++i)
  {
    data[0] += cur_rho_col_[3][i];
  }
  data[0] = data[0]/cur_rho_col_[3].size();


  return true;
}


dserror("The output is only implemented for four different fiber directions!!!");
return false;
}  // VisData()
