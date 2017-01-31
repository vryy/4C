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

#include "elast_utils_autodiff.H"

/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
MAT::ELASTIC::PAR::RemodelFiber::RemodelFiber(
    Teuchos::RCP<MAT::PAR::Material> matdata
)
: Parameter(matdata),
  nummat_(matdata->GetInt("NUMMAT")),
  matids_(matdata->Get<std::vector<int> >("MATIDS")),
  t_decay_(matdata->GetDouble("TDECAY")),
  k_growth_(matdata->GetDouble("GROWTHFAC")),
  init_w_col_(matdata->GetMutable<std::vector<double> >("COLMASSFRAC")),
  G_(matdata->GetDouble("DEPOSITIONSTRETCH"))
{
  // check if sizes fit
  if (nummat_ != (int)matids_->size())
    dserror("number of materials %d does not fit to size of material vector %d", nummat_, matids_->size());

  // check decay time validity
  if (t_decay_<=0.)
    dserror("decay time must be positive");
}


/*----------------------------------------------------------------------*
 |  Constructor                             (public)   fb         09/15 |
 *----------------------------------------------------------------------*/
MAT::ELASTIC::RemodelFiber::RemodelFiber(MAT::ELASTIC::PAR::RemodelFiber* params)
  : params_(params),
    potsumfiber_(0)
{
  // make sure the referenced materials in material list have quick access parameters
  std::vector<int>::const_iterator m;
  for (m=params_->matids_->begin(); m!=params_->matids_->end(); ++m)
  {
    const int matid = *m;
    Teuchos::RCP<MAT::ELASTIC::Summand> sum = MAT::ELASTIC::Summand::Factory(matid);
    if (sum == Teuchos::null) dserror("Failed to allocate");
    potsumfiber_.push_back(Teuchos::rcp(new FiberData(sum)));
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::RemodelFiber::PackSummand(DRT::PackBuffer& data) const
{
  int num_fiber = 0;
  num_fiber = potsumfiber_.size();
  int num_gp = 0;
  num_gp = potsumfiber_[0]->cur_lambr.size();
  AddtoPack(data,num_fiber);
  AddtoPack(data,num_gp);

  for(int i=0;i<num_fiber;++i) {
    AddtoPack(data,potsumfiber_[i]->cur_lambr);
    AddtoPack(data,potsumfiber_[i]->last_lambr);
    AddtoPack(data,potsumfiber_[i]->cur_rho);
    AddtoPack(data,potsumfiber_[i]->last_rho);
    AddtoPack(data,potsumfiber_[i]->AM);
    AddtoPack(data,potsumfiber_[i]->AM_orth);
    AddtoPack(data,potsumfiber_[i]->FrnM);
    AddtoPack(data,potsumfiber_[i]->diFrdlambrM);
    AddtoPack(data,potsumfiber_[i]->dFrdlambrM);
    AddtoPack(data,potsumfiber_[i]->iFrM);
    AddtoPack(data,potsumfiber_[i]->FrdotM);
    AddtoPack(data,potsumfiber_[i]->dFrdotdlambrM);
    AddtoPack(data,potsumfiber_[i]->remodel->sig_h);
    AddtoPack(data,potsumfiber_[i]->remodel->k_sig);
    AddtoPack(data,potsumfiber_[i]->remodel->t_decay);
    AddtoPack(data,potsumfiber_[i]->G);
    AddtoPack(data,cauchystress_[i]);
  }

  AddtoPack(data,init_rho_col_);

  if (params_ != NULL) // summands are not accessible in postprocessing mode
    for (unsigned int k=0; k<potsumfiber_.size(); ++k)
      potsumfiber_[k]->fiber->PackSummand(data);

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::RemodelFiber::UnpackSummand(const std::vector<char>& data,
                                               std::vector<char>::size_type& position)
{
//  // make sure we have a pristine material
//  params_ = NULL;
//  potsumfiber_.clear();

  int num_fiber=0;
  ExtractfromPack(position,data,num_fiber);
  int num_gp = 0;
  ExtractfromPack(position,data,num_gp);

  cauchystress_.resize(num_fiber);

  double sig_h = 0.0;
  double k_sig = 0.0;
  double t_decay = 0.0;
  for(int k=0;k<num_fiber;++k) {
    ExtractfromPack(position,data,potsumfiber_[k]->cur_lambr);
    ExtractfromPack(position,data,potsumfiber_[k]->last_lambr);
    ExtractfromPack(position,data,potsumfiber_[k]->cur_rho);
    ExtractfromPack(position,data,potsumfiber_[k]->last_rho);
    ExtractfromPack(position,data,potsumfiber_[k]->AM);
    ExtractfromPack(position,data,potsumfiber_[k]->AM_orth);
    ExtractfromPack(position,data,potsumfiber_[k]->FrnM);
    ExtractfromPack(position,data,potsumfiber_[k]->diFrdlambrM);
    ExtractfromPack(position,data,potsumfiber_[k]->dFrdlambrM);
    ExtractfromPack(position,data,potsumfiber_[k]->iFrM);
    ExtractfromPack(position,data,potsumfiber_[k]->FrdotM);
    ExtractfromPack(position,data,potsumfiber_[k]->dFrdotdlambrM);
    ExtractfromPack(position,data,sig_h);
    ExtractfromPack(position,data,k_sig);
    ExtractfromPack(position,data,t_decay);
    ExtractfromPack(position,data,potsumfiber_[k]->G);
    ExtractfromPack(position,data,cauchystress_[k]);

    potsumfiber_[k]->growth = Teuchos::rcp(new GrowthEvolution(k_sig,sig_h));
    potsumfiber_[k]->remodel = Teuchos::rcp(new RemodelEvolution(k_sig,sig_h,t_decay));
  }

  ExtractfromPack(position,data,init_rho_col_);

  // loop map of associated potential summands
  for (unsigned int k=0; k<potsumfiber_.size(); ++k)
    potsumfiber_[k]->fiber->UnpackSummand(data,position);

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::RemodelFiber::Setup(int numgp,double rho_tot,DRT::INPUT::LineDefinition* linedef)
{
  // setup fiber and inelastic history variable
  cauchystress_.resize(potsumfiber_.size());
  init_rho_col_.resize(potsumfiber_.size());
  for(unsigned k=0;k<potsumfiber_.size();++k) {
    init_rho_col_[k] = rho_tot * params_->init_w_col_->at(k);
    potsumfiber_[k]->FrnM.resize(numgp);
    potsumfiber_[k]->cur_lambr.resize(numgp,1.0);
    potsumfiber_[k]->cur_rho.resize(numgp,init_rho_col_[k]);
    potsumfiber_[k]->diFrdlambrM.resize(numgp);
    potsumfiber_[k]->dFrdlambrM.resize(numgp);
    potsumfiber_[k]->iFrM.resize(numgp);
    potsumfiber_[k]->FrdotM.resize(numgp);
    potsumfiber_[k]->dFrdotdlambrM.resize(numgp);
    potsumfiber_[k]->last_lambr.resize(numgp,1.0);
    potsumfiber_[k]->last_rho.resize(numgp,init_rho_col_[k]);
    potsumfiber_[k]->G = params_->G_;
    cauchystress_[k].resize(numgp,1.0);

    potsumfiber_[k]->fiber->Setup(linedef);
  }


  // some variables
  LINALG::Matrix<2,1> dPI(true);
  LINALG::Matrix<3,1> ddPII(true);
  LINALG::Matrix<4,1> dddPIII(true);
  LINALG::Matrix<6,1> stressactv(true);
  LINALG::Matrix<6,6> cmatactive(true);
  LINALG::Matrix<3,3> stressactM(true);

  SetupStructuralTensorsGR();

  // quadratic prestretch in tensor notation
  LINALG::Matrix<3,3> CpreM(true);

  // temporary pointers to check the type of the remodelfiber (active or passive)
  Teuchos::RCP<MAT::ELASTIC::CoupAnisoExpo> t1(nullptr);
  Teuchos::RCP<MAT::ELASTIC::CoupAnisoExpoActive> t2(nullptr);

  // identity matrix
  LINALG::Matrix<3,3> id(true);
  for(int i=0;i<3;++i) id(i,i) = 1.0;

  double sig_pre = 0.0;
  for(unsigned k=0;k<potsumfiber_.size();++k)
  {
    if( (t1 = Teuchos::rcp_dynamic_cast<MAT::ELASTIC::CoupAnisoExpo>(potsumfiber_[k]->fiber)).getRawPtr() ) {
      CpreM.Update(potsumfiber_[k]->G*potsumfiber_[k]->G,potsumfiber_[k]->AM,0.0);
      t1->GetDerivativesAniso(dPI,ddPII,dddPIII,CpreM,0);
      sig_pre = 2.0*dPI(0)*potsumfiber_[k]->G*potsumfiber_[k]->G;
      for(int gp=0;gp<numgp;++gp)
        cauchystress_[k][gp] = sig_pre;
      potsumfiber_[k]->growth = Teuchos::rcp(new GrowthEvolution(params_->k_growth_,sig_pre));
      potsumfiber_[k]->remodel = Teuchos::rcp(new RemodelEvolution(params_->k_growth_,sig_pre,params_->t_decay_));
    }
    else if( (t2 = Teuchos::rcp_dynamic_cast<MAT::ELASTIC::CoupAnisoExpoActive>(potsumfiber_[k]->fiber)).getRawPtr() ) {
      CpreM.Update(potsumfiber_[k]->G*potsumfiber_[k]->G,potsumfiber_[k]->AM,0.0);
      t2->GetDerivativesAniso(dPI,ddPII,dddPIII,CpreM,0);
      sig_pre = 2.0*dPI(0)*potsumfiber_[k]->G*potsumfiber_[k]->G;
      t2->EvaluateActiveStressCmatAniso(id,cmatactive,stressactv,0);
      StressVoigtNotationVectorToMatrix(stressactv,stressactM);
      sig_pre += stressactM.Dot(potsumfiber_[k]->AM);
      for(int gp=0;gp<numgp;++gp)
        cauchystress_[k][gp] = sig_pre;
      potsumfiber_[k]->growth = Teuchos::rcp(new GrowthEvolution(params_->k_growth_,sig_pre));
      potsumfiber_[k]->remodel = Teuchos::rcp(new RemodelEvolution(params_->k_growth_,sig_pre,params_->t_decay_));
    }
    else
      dserror("So far, you can only use Elast_CoupAnisoExpo and Elast_CoupAnisoExpoActive in Elast_Remodelfiber!");
  }

  // Initialize inelastic deformation gradients in FiberData (default time step size (does not have to be the real one))
  for(unsigned k=0;k<potsumfiber_.size();++k)
    for(int gp=0;gp<numgp;++gp)
      potsumfiber_[k]->UpdateNewton(gp,1.0);


  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::RemodelFiber::SetupStructuralTensorsGR()
{
  // identity tensor
  LINALG::Matrix<3,3> id(true);
  for(int i=0;i<3;++i) id(i,i) = 1.0;

  // fiber directions
  std::vector<LINALG::Matrix<3,1> > fibervecs;

  for(unsigned k=0;k<potsumfiber_.size();++k) {
    // Get fiberdirection
    potsumfiber_[k]->fiber->GetFiberVecs(fibervecs);

    // build structural tensor in matrix notation
    potsumfiber_[k]->AM.MultiplyNT(1.0,fibervecs[k],fibervecs[k],0.0);
    // orthogonal structural tensor ( 1_{ij} - A_{ij} )
    potsumfiber_[k]->AM_orth.Update(1.0,potsumfiber_[k]->AM,0.0);
    potsumfiber_[k]->AM_orth.Update(1.0,id,-1.0);
  }

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::RemodelFiber::Update()
{
  // update history variable
  for(unsigned k=0;k<potsumfiber_.size();++k)
    for(unsigned gp=0;gp<potsumfiber_[k]->cur_rho.size();++gp)
      potsumfiber_[k]->UpdateHistory(gp);

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::RemodelFiber::EvaluateAnisotropicStressCmat(LINALG::Matrix<3,3> const& CM,
                                                               LINALG::Matrix<3,3> const& iFgM,
                                                               LINALG::Matrix<6,6> & cmat,
                                                               LINALG::Matrix<6,1> & stress,
                                                               int const gp,
                                                               double const& dt,
                                                               int const eleGID)
{
  // clear some variables
  stress.Clear();
  cmat.Clear();

  for(unsigned k=0;k<potsumfiber_.size();++k) {
    potsumfiber_[k]->UpdateNewton(gp,dt);
    AddStressCmat(CM,iFgM,*(potsumfiber_[k]),gp,eleGID,stress,cmat);
  }

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::RemodelFiber::EvaluateDerivativesInternalNewton(LINALG::Matrix<3,3> const * const defgrd,
                                                                   int const nr_grf_proc,
                                                                   int const nr_grf_tot,
                                                                   int const gp,
                                                                   double const& dt,
                                                                   int const eleGID,
                                                                   LINALG::Matrix<3,3> const& iFgM,
                                                                   LINALG::Matrix<3,3> const& dFgdrhoM,
                                                                   LINALG::Matrix<3,3> const& diFgdrhoM,
                                                                   std::vector<std::vector<double> > & dWdrho,
                                                                   std::vector<std::vector<double> > & dWdlambr,
                                                                   std::vector<double> & W,
                                                                   std::vector<std::vector<double> > & dEdrho,
                                                                   std::vector<std::vector<double> > & dEdlambr,
                                                                   std::vector<double> & E)
{
  static LINALG::Matrix<3,3> CM(true);
  CM.MultiplyTN(1.0,*defgrd,*defgrd,0.0);

  for(unsigned k=0;k<potsumfiber_.size();++k) {
    potsumfiber_[k]->UpdateNewton(gp,dt);

    // Residual of evolution equations
    EvaluateEvolutionEquation(W[nr_grf_proc+k],E[nr_grf_proc+k],CM,iFgM,dt,*(potsumfiber_[k]),gp,eleGID);

    // Derivatives of evolution equations
    double dWidrhoi = 0.0;
    double dWidrhoj = 0.0;
    double dEidrho = 0.0;
    EvaluateDerivativeEvolutionEquation(dWidrhoi,dWidrhoj,dWdlambr[nr_grf_proc+k][nr_grf_proc+k],dEidrho,
        dEdlambr[nr_grf_proc+k][nr_grf_proc+k],CM,iFgM,dFgdrhoM,diFgdrhoM,dt,*(potsumfiber_[k]),gp,eleGID);

    for(int l=0;l<nr_grf_tot;++l) {
      dEdrho[nr_grf_proc+k][l] = dEidrho;
      if(l == (k+nr_grf_proc))
        dWdrho[nr_grf_proc+k][l] = dWidrhoi;
      else
        dWdrho[nr_grf_proc+k][l] = dWidrhoj;
    }
  }

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::RemodelFiber::EvaluateDerivativesCauchyGreen(LINALG::Matrix<3,3> const* const defgrd,
                                                                int const nr_grf_proc,
                                                                int const gp,
                                                                double const& dt,
                                                                LINALG::Matrix<3,3> const& iFgM,
                                                                std::vector<LINALG::Matrix<1,6> > & dWdC,
                                                                std::vector<LINALG::Matrix<1,6> > & dEdC,
                                                                int const eleGID)
{
  static LINALG::Matrix<3,3> CM(true);
  CM.MultiplyTN(1.0,*defgrd,*defgrd,0.0);

  for(unsigned k=0;k<potsumfiber_.size();++k) {
    dEdC[nr_grf_proc+k].Clear();
    dWdC[nr_grf_proc+k].Clear();
    potsumfiber_[k]->UpdateNewton(gp,dt);

    EvaluatedEvolutionEquationdC(dWdC[nr_grf_proc+k],dEdC[nr_grf_proc+k],CM,iFgM,dt,*(potsumfiber_[k]),k,gp,eleGID);
  }

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::RemodelFiber::EvaluateAdditionalGrowthRemodelCmat(LINALG::Matrix<3,3> const* const defgrd,
                                                                     int const nr_grf_proc,
                                                                     LINALG::Matrix<3,3> const& iFgM,
                                                                     LINALG::Matrix<3,3> const& diFgdrhoM,
                                                                     std::vector<LINALG::Matrix<1,6> > const& drhodC,
                                                                     std::vector<LINALG::Matrix<1,6> > const& dlambrdC,
                                                                     LINALG::Matrix<6,6> & cmat,
                                                                     int const gp,
                                                                     int const eleGID) const
{
  // clear some variables
  cmat.Clear();

  static LINALG::Matrix<3,3> CM(true);
  CM.MultiplyTN(1.0,*defgrd,*defgrd,0.0);

  static LINALG::Matrix<6,1> dSidrhoi(true);
  static LINALG::Matrix<6,1> dSidrhoj(true);
  static LINALG::Matrix<6,1> dSdlambr(true);
  for(unsigned k=0;k<potsumfiber_.size();++k) {
    EvaluateDerivatives2ndPiolaKirchhoffGrowthRemodel(dSidrhoi,dSidrhoj,dSdlambr,CM,iFgM,diFgdrhoM,*(potsumfiber_[k]),gp,eleGID);

    for(int l=0;l<drhodC.size();++l) {
      if(l == nr_grf_proc+k)
        cmat.MultiplyNN(2.0,dSidrhoi,drhodC[l],1.0);
      else
        cmat.MultiplyNN(2.0,dSidrhoj,drhodC[l],1.0);
    }
    cmat.MultiplyNN(2.0,dSdlambr,dlambrdC[nr_grf_proc+k],1.0);
  }

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::RemodelFiber::EvaluateGrowthAndRemodelingExpl(LINALG::Matrix<3,3> const& defgrd,
                                                                 double const& dt,
                                                                 LINALG::Matrix<3,3> const& iFgM,
                                                                 const int gp,
                                                                 const int eleGID)
{
  static LINALG::Matrix<3,3> CM(true);
  CM.MultiplyTN(1.0,defgrd,defgrd,0.0);
  double drhodt = 0.0;
  double dlambrdt = 0.0;

  for(unsigned k=0;k<potsumfiber_.size();++k) {
    potsumfiber_[k]->UpdateNewton(gp,dt);

    EvaluatedEvolutionEquationdt(drhodt,dlambrdt,CM,iFgM,*(potsumfiber_[k]),k,gp,eleGID);

    UpdateGrowthRemodelParameter(drhodt*dt,dlambrdt*dt,k,gp);
  }

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template<class FUNC, typename T, typename ForceAnalytical>
void MAT::ELASTIC::RemodelFiber::DerivdC(LINALG::TMatrix<T,3,3> const& CM,
                                         LINALG::TMatrix<T,3,3> const& iFinM,
                                         LINALG::TMatrix<T,3,3> const& AM,
                                         FUNC const& func,
                                         ForceAnalytical const eleGID,
                                         LINALG::TMatrix<T,3,3> & dfuncdC) const
{
  // clear some variables
  dfuncdC.Clear();

  static LINALG::FADMatrix<3,3> iFinM_fad(true);
  iFinM_fad = iFinM;

  // Setup FAD
  // first derivative
  static LINALG::FADMatrix<3,3> CM_fad(true);
  CM_fad = CM;
  CM_fad.diff(0,9);

  static LINALG::FADMatrix<3,3> CeM_fad(true);
  static LINALG::FADMatrix<3,3> tmp_fad(true);
  tmp_fad.MultiplyNN(1.0,CM_fad,iFinM_fad,0.0);
  CeM_fad.MultiplyTN(1.0,iFinM_fad,tmp_fad,0.0);

  FAD r_fad = 0.0;
  func.EvaluateFunc(r_fad,CeM_fad,eleGID);

  LINALG::Matrix<3,3> tmp(true);
  FirstDerivToMatrix(r_fad,tmp);
  dfuncdC.Update(1.0,tmp,0.0);

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template<class FUNC, typename T>
void MAT::ELASTIC::RemodelFiber::DerivdC(LINALG::TMatrix<T,3,3> const& CM,
                                         LINALG::TMatrix<T,3,3> const& iFinM,
                                         LINALG::TMatrix<T,3,3> const& AM,
                                         FUNC const& func,
                                         int const eleGID,
                                         LINALG::TMatrix<T,3,3>& dfuncdC) const
{
  // clear some variables
  dfuncdC.Clear();

  // elastic right Cauchy-Green in matrix notation
  static LINALG::TMatrix<T,3,3> tmp(true);
  static LINALG::TMatrix<T,3,3> CeM(true);
  static LINALG::TMatrix<T,3,3> FinM(true);
  static LINALG::TMatrix<T,3,3> CinM(true);
  FinM.Invert(iFinM);
  CinM.MultiplyTN(1.0,FinM,FinM,0.0);
  tmp.MultiplyNN(1.0,CM,iFinM,0.0);
  CeM.MultiplyTN(1.0,iFinM,tmp,0.0);

  // get derivatives of strain energy function w.r.t. I4
  static LINALG::TMatrix<T,2,1> dPIe(true);
  static LINALG::TMatrix<T,3,1> ddPIIe(true);
  static LINALG::TMatrix<T,4,1> dddPIIIe(true);
  func.GetDerivativesAniso(dPIe,ddPIIe,dddPIIIe,CeM,eleGID);

  dfuncdC.Update(dPIe(0)/CinM.Dot(AM),AM,0.0);

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template<class FUNC, typename T, typename ForceAnalytical>
void MAT::ELASTIC::RemodelFiber::DerivdCdC(LINALG::TMatrix<T,3,3> const& CM,
                                           LINALG::TMatrix<T,3,3> const& iFinM,
                                           LINALG::TMatrix<T,3,3> const& AM,
                                           FUNC const& func,
                                           ForceAnalytical const eleGID,
                                           LINALG::TMatrix<T,6,6>& dfuncdCdC) const
{
  // clear some variables
  dfuncdCdC.Clear();

  static LINALG::FADMatrix<3,3> iFinM_fad(true);
  iFinM_fad = iFinM;
  static LINALG::FADMatrix<3,3> AM_fad(true);
  AM_fad = AM;

  // Setup FAD
  // first derivative
  static LINALG::FADMatrix<3,3> CM_fad(true);
  CM_fad = CM;
  CM_fad.diff(0,9);

  LINALG::FADMatrix<3,3> R_fad(true);
  DerivdC(CM_fad,iFinM_fad,AM_fad,func,eleGID,R_fad);

  for(int i=0;i<3;++i)
    for(int j=0;j<6;++j)
      dfuncdCdC(i,j) = R_fad(i,i).dx(j);
  for(int j=0;j<6;++j) {
    dfuncdCdC(3,j) = R_fad(0,1).dx(j);
    dfuncdCdC(4,j) = R_fad(1,2).dx(j);
    dfuncdCdC(5,j) = R_fad(0,2).dx(j);
  }

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template<class FUNC, typename T>
void MAT::ELASTIC::RemodelFiber::DerivdCdC(LINALG::TMatrix<T,3,3> const& CM,
                                           LINALG::TMatrix<T,3,3> const& iFinM,
                                           LINALG::TMatrix<T,3,3> const& AM,
                                           FUNC const& func,
                                           int const eleGID,
                                           LINALG::TMatrix<T,6,6>& dfuncdCdC) const
{
  // clear some variables
  dfuncdCdC.Clear();

  // elastic right Cauchy-Green in matrix notation
  static LINALG::TMatrix<T,3,3> tmp(true);
  static LINALG::TMatrix<T,3,3> CeM(true);
  static LINALG::TMatrix<T,3,3> FinM(true);
  static LINALG::TMatrix<T,3,3> CinM(true);
  FinM.Invert(iFinM);
  CinM.MultiplyTN(1.0,FinM,FinM,0.0);
  tmp.MultiplyNN(1.0,CM,iFinM,0.0);
  CeM.MultiplyTN(1.0,iFinM,tmp,0.0);

  // get derivatives of strain energy function w.r.t. I4
  static LINALG::TMatrix<T,2,1> dPIe(true);
  static LINALG::TMatrix<T,3,1> ddPIIe(true);
  static LINALG::TMatrix<T,4,1> dddPIIIe(true);
  func.GetDerivativesAniso(dPIe,ddPIIe,dddPIIIe,CeM,eleGID);

  static LINALG::TMatrix<T,6,1> Av(true);
  MatrixtoStressVoigtNotationVector(AM,Av);
  dfuncdCdC.MultiplyNT(ddPIIe(0)/(CinM.Dot(AM)*CinM.Dot(AM)),Av,Av,0.0);

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::RemodelFiber::AddStressCmat(LINALG::Matrix<3,3> const& CM,
                                               LINALG::Matrix<3,3> const& iFgM,
                                               FiberData const& fiberdat,
                                               int const gp,
                                               int const eleGID,
                                               LINALG::Matrix<6,1>& stress,
                                               LINALG::Matrix<6,6>& cmat) const
{
  static LINALG::Matrix<3,3> iFinM(true);
  iFinM.MultiplyNN(1.0,iFgM,fiberdat.iFrM[gp],0.0);

  // temporary pointers to check the type of the remodelfiber (active or passive)
  Teuchos::RCP<MAT::ELASTIC::CoupAnisoExpo> t1(nullptr);
  Teuchos::RCP<MAT::ELASTIC::CoupAnisoExpoActive> t2(nullptr);

  static LINALG::Matrix<3,3> firstderivM(true);
  static LINALG::Matrix<6,1> firstderivv(true);
  static LINALG::Matrix<6,6> secderiv(true);
  static LINALG::Matrix<6,1> stressactv(true);
  static LINALG::Matrix<6,6> cmatact(true);
  if( (t1 = Teuchos::rcp_dynamic_cast<MAT::ELASTIC::CoupAnisoExpo>(fiberdat.fiber)).getRawPtr() ) {
    DerivdC(CM,iFinM,fiberdat.AM,*t1,eleGID,firstderivM);
    DerivdCdC(CM,iFinM,fiberdat.AM,*t1,eleGID,secderiv);
  }
  else if( (t2 = Teuchos::rcp_dynamic_cast<MAT::ELASTIC::CoupAnisoExpoActive>(fiberdat.fiber)).getRawPtr() ) {
    DerivdC(CM,iFinM,fiberdat.AM,*t2,eleGID,firstderivM);
    DerivdCdC(CM,iFinM,fiberdat.AM,*t2,eleGID,secderiv);
    t2->EvaluateActiveStressCmatAniso(CM,cmatact,stressactv,eleGID);
    stress.Update(fiberdat.cur_rho[gp],stressactv,1.0);
    cmat.Update(fiberdat.cur_rho[gp],cmatact,1.0);
  }

  MatrixtoStressVoigtNotationVector(firstderivM,firstderivv);
  stress.Update(2.0*fiberdat.cur_rho[gp],firstderivv,1.0);
  cmat.Update(4.0*fiberdat.cur_rho[gp],secderiv,1.0);

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template< typename T >
void MAT::ELASTIC::RemodelFiber::EvaluateLocalCauchyStress(LINALG::TMatrix<T,3,3> const& CM,
                                                           LINALG::TMatrix<T,3,3> const& iFinM,
                                                           LINALG::TMatrix<T,3,3> const& AM,
                                                           Teuchos::RCP<MAT::ELASTIC::Summand> const fiber,
                                                           int const eleGID,
                                                           T & sig) const
{
  // temporary pointers to check the type of the remodelfiber (active or passive)
  Teuchos::RCP<MAT::ELASTIC::CoupAnisoExpo> t1(nullptr);
  Teuchos::RCP<MAT::ELASTIC::CoupAnisoExpoActive> t2(nullptr);

  static LINALG::TMatrix<T,3,3> tmp(true);
  static LINALG::TMatrix<T,3,3> CeM(true);
  static LINALG::TMatrix<T,3,3> FinM(true);
  static LINALG::TMatrix<T,3,3> CinM(true);
  FinM.Invert(iFinM);
  CinM.MultiplyTN(1.0,FinM,FinM,0.0);
  tmp.MultiplyNN(1.0,CM,iFinM,0.0);
  CeM.MultiplyTN(1.0,iFinM,tmp,0.0);
  static LINALG::TMatrix<T,2,1> dPIe(true);
  static LINALG::TMatrix<T,3,1> ddPIIe(true);
  static LINALG::TMatrix<T,4,1> dddPIIIe(true);
  T dPIact = 0.0;
  if( (t1 = Teuchos::rcp_dynamic_cast<MAT::ELASTIC::CoupAnisoExpo>(fiber)).getRawPtr() ) {
    t1->GetDerivativesAniso(dPIe,ddPIIe,dddPIIIe,CeM,eleGID);
  }
  else if( (t2 = Teuchos::rcp_dynamic_cast<MAT::ELASTIC::CoupAnisoExpoActive>(fiber)).getRawPtr() ) {
    t2->GetDerivativesAniso(dPIe,ddPIIe,dddPIIIe,CeM,eleGID);
    t2->GetDerivativeAnisoActive(dPIact);
  }

  sig = 2.0*dPIe(0)*CM.Dot(AM)/CinM.Dot(AM) + dPIact;

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template< typename T >
void MAT::ELASTIC::RemodelFiber::EvaluatedsigdCe(LINALG::TMatrix<T,3,3> const& CM,
                                                 LINALG::TMatrix<T,3,3> const& iFgM,
                                                 LINALG::TMatrix<T,3,3> const& iFrM,
                                                 LINALG::TMatrix<T,3,3> const& AM,
                                                 Teuchos::RCP<MAT::ELASTIC::Summand> const fiber,
                                                 int const eleGID,
                                                 LINALG::TMatrix<T,3,3>& dsigdCe) const
{
  // clear some variables
  dsigdCe.Clear();

  // temporary pointers to check the type of the remodelfiber (active or passive)
  Teuchos::RCP<MAT::ELASTIC::CoupAnisoExpo> t1(nullptr);
  Teuchos::RCP<MAT::ELASTIC::CoupAnisoExpoActive> t2(nullptr);

  static LINALG::TMatrix<T,3,3> tmp(true);
  static LINALG::TMatrix<T,3,3> CeM(true);
  static LINALG::TMatrix<T,3,3> FinM(true);
  static LINALG::TMatrix<T,3,3> CinM(true);
  static LINALG::TMatrix<T,3,3> iFinM(true);
  static LINALG::TMatrix<T,3,3> AgrM(true);
  iFinM.MultiplyNN(1.0,iFgM,iFrM,0.0);
  FinM.Invert(iFinM);
  CinM.MultiplyTN(1.0,FinM,FinM,0.0);
  tmp.MultiplyNN(1.0,FinM,AM,0.0);
  AgrM.MultiplyNT(1.0/CinM.Dot(AM),tmp,FinM,0.0);
  tmp.MultiplyNN(1.0,CM,iFinM,0.0);
  CeM.MultiplyTN(1.0,iFinM,tmp,0.0);
  static LINALG::TMatrix<T,2,1> dPIe(true);
  static LINALG::TMatrix<T,3,1> ddPIIe(true);
  static LINALG::TMatrix<T,4,1> dddPIIIe(true);
  if( (t1 = Teuchos::rcp_dynamic_cast<MAT::ELASTIC::CoupAnisoExpo>(fiber)).getRawPtr() ) {
    t1->GetDerivativesAniso(dPIe,ddPIIe,dddPIIIe,CeM,eleGID);
  }
  else if( (t2 = Teuchos::rcp_dynamic_cast<MAT::ELASTIC::CoupAnisoExpoActive>(fiber)).getRawPtr() ) {
    t2->GetDerivativesAniso(dPIe,ddPIIe,dddPIIIe,CeM,eleGID);
  }

  dsigdCe.Update(2.0*(ddPIIe(0)*CeM.Dot(AgrM)+dPIe(0)),AgrM,0.0);

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template< typename ForceAnalytical >
void MAT::ELASTIC::RemodelFiber::EvaluatedsigdCedC(LINALG::Matrix<3,3> const& CM,
                                                   LINALG::Matrix<3,3> const& iFgM,
                                                   LINALG::Matrix<3,3> const& iFrM,
                                                   LINALG::Matrix<3,3> const& AM,
                                                   Teuchos::RCP<MAT::ELASTIC::Summand> const fiber,
                                                   ForceAnalytical const eleGID,
                                                   LINALG::Matrix<6,6>& dsigdCedC) const
{
  // clear some variables
  dsigdCedC.Clear();

  static LINALG::FADMatrix<3,3> CM_fad(true);
  CM_fad = CM;
  CM_fad.diff(0,9);
  static LINALG::FADMatrix<3,3> iFgM_fad(true);
  iFgM_fad = iFgM;
  static LINALG::FADMatrix<3,3> iFrM_fad(true);
  iFrM_fad = iFrM;
  static LINALG::FADMatrix<3,3> AM_fad(true);
  AM_fad = AM;

  static LINALG::FADMatrix<3,3> dsigdCeM_fad(true);
  EvaluatedsigdCe(CM_fad,iFgM_fad,iFrM_fad,AM_fad,fiber,eleGID,dsigdCeM_fad);

  for(int i=0;i<3;++i) {
    for(int j=0;j<3;++j)
      dsigdCedC(i,j) = dsigdCeM_fad(i,i).dx(j);
    for(int j=3;j<6;++j)
      dsigdCedC(i,j) = 0.5*(dsigdCeM_fad(i,i).dx(j) + dsigdCeM_fad(i,i).dx(j+3));
  }
  for(int j=0;j<3;++j) {
    dsigdCedC(3,j) = 0.5*(dsigdCeM_fad(0,1).dx(j) + dsigdCeM_fad(1,0).dx(j));
    dsigdCedC(4,j) = 0.5*(dsigdCeM_fad(1,2).dx(j) + dsigdCeM_fad(2,1).dx(j));
    dsigdCedC(5,j) = 0.5*(dsigdCeM_fad(0,2).dx(j) + dsigdCeM_fad(2,0).dx(j));
  }
  for(int j=3;j<6;++j) {
    dsigdCedC(3,j) = 0.25*(dsigdCeM_fad(0,1).dx(j) + dsigdCeM_fad(0,1).dx(j+3) + dsigdCeM_fad(1,0).dx(j) + dsigdCeM_fad(1,0).dx(j+3));
    dsigdCedC(4,j) = 0.25*(dsigdCeM_fad(1,2).dx(j) + dsigdCeM_fad(1,2).dx(j+3) + dsigdCeM_fad(2,1).dx(j) + dsigdCeM_fad(2,1).dx(j+3));
    dsigdCedC(5,j) = 0.25*(dsigdCeM_fad(0,2).dx(j) + dsigdCeM_fad(0,2).dx(j+3) + dsigdCeM_fad(2,0).dx(j) + dsigdCeM_fad(2,0).dx(j+3));
  }

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::RemodelFiber::EvaluatedsigdCedC(LINALG::Matrix<3,3> const& CM,
                                                   LINALG::Matrix<3,3> const& iFgM,
                                                   LINALG::Matrix<3,3> const& iFrM,
                                                   LINALG::Matrix<3,3> const& AM,
                                                   Teuchos::RCP<MAT::ELASTIC::Summand> const fiber,
                                                   int const eleGID,
                                                   LINALG::Matrix<6,6>& dsigdCedC) const
{
  // clear some variables
  dsigdCedC.Clear();

  // temporary pointers to check the type of the remodelfiber (active or passive)
  Teuchos::RCP<MAT::ELASTIC::CoupAnisoExpo> t1(nullptr);
  Teuchos::RCP<MAT::ELASTIC::CoupAnisoExpoActive> t2(nullptr);

  static LINALG::Matrix<3,3> tmp(true);
  static LINALG::Matrix<3,3> CeM(true);
  static LINALG::Matrix<3,3> FinM(true);
  static LINALG::Matrix<3,3> CinM(true);
  static LINALG::Matrix<3,3> iFinM(true);
  static LINALG::Matrix<3,3> AgrM(true);
  iFinM.MultiplyNN(1.0,iFgM,iFrM,0.0);
  FinM.Invert(iFinM);
  CinM.MultiplyTN(1.0,FinM,FinM,0.0);
  tmp.MultiplyNN(1.0,FinM,AM,0.0);
  AgrM.MultiplyNT(1.0/CinM.Dot(AM),tmp,FinM,0.0);
  tmp.MultiplyNN(1.0,CM,iFinM,0.0);
  CeM.MultiplyTN(1.0,iFinM,tmp,0.0);
  static LINALG::Matrix<2,1> dPIe(true);
  static LINALG::Matrix<3,1> ddPIIe(true);
  static LINALG::Matrix<4,1> dddPIIIe(true);
  if( (t1 = Teuchos::rcp_dynamic_cast<MAT::ELASTIC::CoupAnisoExpo>(fiber)).getRawPtr() ) {
    t1->GetDerivativesAniso(dPIe,ddPIIe,dddPIIIe,CeM,eleGID);
  }
  else if( (t2 = Teuchos::rcp_dynamic_cast<MAT::ELASTIC::CoupAnisoExpoActive>(fiber)).getRawPtr() ) {
    t2->GetDerivativesAniso(dPIe,ddPIIe,dddPIIIe,CeM,eleGID);
  }

  static LINALG::Matrix<6,1> Agrv(true);
  static LINALG::Matrix<6,1> Av(true);
  MatrixtoStressVoigtNotationVector(AgrM,Agrv);
  MatrixtoStressVoigtNotationVector(AM,Av);
  dsigdCedC.MultiplyNT(2.0/CinM.Dot(AM)*(dddPIIIe(0)*CM.Dot(AM)/CinM.Dot(AM) + 2.0*ddPIIe(0)),Agrv,Av,0.0);

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template< typename ForceAnalytical >
void MAT::ELASTIC::RemodelFiber::EvaluateDerivativesCauchyGrowth(LINALG::Matrix<3,3> const& CM,
                                                                 LINALG::Matrix<3,3> const& iFgM,
                                                                 LINALG::Matrix<3,3> const& dFgdrhoM,
                                                                 LINALG::Matrix<3,3> const& diFgdrhoM,
                                                                 FiberData const& fiberdat,
                                                                 int const gp,
                                                                 ForceAnalytical const eleGID,
                                                                 double & dsigdrho,
                                                                 LINALG::Matrix<3,3> & dsigdCedrhoM) const
{
  // clear some variables
  dsigdrho = 0.0;
  dsigdCedrhoM.Clear();

  static LINALG::FADMatrix<3,3> iFinM_fad(true);
  static LINALG::FADMatrix<3,3> iFgM_fad(true);
  static LINALG::FADMatrix<3,3> iFrM_fad(true);
  static LINALG::FADMatrix<3,3> CM_fad(true);
  static LINALG::FADMatrix<3,3> AM_fad(true);
  iFgM_fad = iFgM;
  iFgM_fad.diff(0,9);
  iFrM_fad = fiberdat.iFrM[gp];
  CM_fad = CM;
  AM_fad = fiberdat.AM;
  iFinM_fad.MultiplyNN(1.0,iFgM_fad,iFrM_fad,0.0);

  FAD sig_fad = 0.0;
  EvaluateLocalCauchyStress(CM_fad,iFinM_fad,AM_fad,fiberdat.fiber,eleGID,sig_fad);

  LINALG::Matrix<3,3> dsigdiFgM(true);
  FirstDerivToMatrix(sig_fad,dsigdiFgM);
  dsigdrho = dsigdiFgM.Dot(diFgdrhoM);


  static LINALG::FADMatrix<3,3> dsigdCeM_fad(true);
  static LINALG::Matrix<6,9> dsigdCediFg(true);
  EvaluatedsigdCe(CM_fad,iFgM_fad,iFrM_fad,AM_fad,fiberdat.fiber,eleGID,dsigdCeM_fad);

  for(int i=0;i<3;++i)
    for(int j=0;j<9;++j)
      dsigdCediFg(i,j) = dsigdCeM_fad(i,i).dx(j);
  for(int j=0;j<9;++j)
    dsigdCediFg(3,j) = 0.5*(dsigdCeM_fad(0,1).dx(j) + dsigdCeM_fad(1,0).dx(j));
  for(int j=0;j<9;++j)
    dsigdCediFg(4,j) = 0.5*(dsigdCeM_fad(1,2).dx(j) + dsigdCeM_fad(2,1).dx(j));
  for(int j=0;j<9;++j)
    dsigdCediFg(5,j) = 0.5*(dsigdCeM_fad(0,2).dx(j) + dsigdCeM_fad(2,0).dx(j));

  static LINALG::Matrix<9,1> diFgdrho9x1(true);
  static LINALG::Matrix<6,1> tmp6x1(true);
  Matrix3x3toVector9x1(diFgdrhoM,diFgdrho9x1);
  tmp6x1.MultiplyNN(1.0,dsigdCediFg,diFgdrho9x1,0.0);
  StressVoigtNotationVectorToMatrix(tmp6x1,dsigdCedrhoM);

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::RemodelFiber::EvaluateDerivativesCauchyGrowth(LINALG::Matrix<3,3> const& CM,
                                                                 LINALG::Matrix<3,3> const& iFgM,
                                                                 LINALG::Matrix<3,3> const& dFgdrhoM,
                                                                 LINALG::Matrix<3,3> const& diFgdrhoM,
                                                                 FiberData const& fiberdat,
                                                                 int const gp,
                                                                 int const eleGID,
                                                                 double & dsigdrho,
                                                                 LINALG::Matrix<3,3> & dsigdCedrhoM) const
{
  // clear some variables
  dsigdrho = 0.0;
  dsigdCedrhoM.Clear();

  static LINALG::Matrix<3,3> CAM(true);
  CAM.MultiplyNN(1.0,CM,fiberdat.AM,0.0);
  static LINALG::Matrix<3,3> tmp(true);
  static LINALG::Matrix<3,3> CeM(true);
  static LINALG::Matrix<3,3> iFinM(true);
  static LINALG::Matrix<3,3> FgM(true);
  static LINALG::Matrix<3,3> FinM(true);
  static LINALG::Matrix<3,3> CinM(true);
  static LINALG::Matrix<3,3> CinAM(true);
  FgM.Invert(iFgM);
  iFinM.MultiplyNN(1.0,iFgM,fiberdat.iFrM[gp],0.0);
  FinM.Invert(iFinM);
  CinM.MultiplyTN(1.0,FinM,FinM,0.0);
  CinAM.MultiplyTN(1.0,CinM,fiberdat.AM,0.0);
  tmp.MultiplyNN(1.0,CM,iFinM,0.0);
  CeM.MultiplyTN(1.0,iFinM,tmp,0.0);

  // temporary pointers to check the type of the remodelfiber (active or passive)
  Teuchos::RCP<MAT::ELASTIC::CoupAnisoExpo> t1(nullptr);
  Teuchos::RCP<MAT::ELASTIC::CoupAnisoExpoActive> t2(nullptr);

  static LINALG::Matrix<2,1> dPIe(true);
  static LINALG::Matrix<3,1> ddPIIe(true);
  static LINALG::Matrix<4,1> dddPIIIe(true);
  if( (t1 = Teuchos::rcp_dynamic_cast<MAT::ELASTIC::CoupAnisoExpo>(fiberdat.fiber)).getRawPtr() )
    t1->GetDerivativesAniso(dPIe,ddPIIe,dddPIIIe,CeM,eleGID);
  else if( (t2 = Teuchos::rcp_dynamic_cast<MAT::ELASTIC::CoupAnisoExpoActive>(fiberdat.fiber)).getRawPtr() )
    t2->GetDerivativesAniso(dPIe,ddPIIe,dddPIIIe,CeM,eleGID);

  static LINALG::Matrix<3,3> dsigdiFgM(true);
  dsigdiFgM.MultiplyNT(4.0*ddPIIe(0)*CM.Dot(fiberdat.AM)/(CinM.Dot(fiberdat.AM)*CinM.Dot(fiberdat.AM)),CAM,FgM,0.0);
  dsigdiFgM.MultiplyNT(4.0*dPIe(0)*CM.Dot(fiberdat.AM)/(CinM.Dot(fiberdat.AM)*CinM.Dot(fiberdat.AM)),CinAM,FgM,1.0);
  dsigdrho = dsigdiFgM.Dot(diFgdrhoM);


  static LINALG::Matrix<3,3> AgrM(true);
  static LINALG::Matrix<3,3> CAFgTM(true);
  static LINALG::Matrix<3,3> CinAFgTM(true);
  tmp.MultiplyNN(1.0,FinM,fiberdat.AM,0.0);
  AgrM.MultiplyNT(1.0/CinM.Dot(fiberdat.AM),tmp,FinM,0.0);
  CAFgTM.MultiplyNT(1.0,CAM,FgM,0.0);
  CinAFgTM.MultiplyNT(1.0,CinAM,FgM,0.0);

  dsigdCedrhoM.Update(4.0/CinM.Dot(fiberdat.AM)*(dddPIIIe(0)*CAFgTM.Dot(diFgdrhoM)*CeM.Dot(AgrM) +
                                                 ddPIIe(0)*CM.Dot(fiberdat.AM)*CinAFgTM.Dot(diFgdrhoM)/CinM.Dot(fiberdat.AM) +
                                                 ddPIIe(0)*CAFgTM.Dot(diFgdrhoM)),AgrM,0.0);

  static LINALG::Matrix<3,3> FrM(true);
  FrM.Invert(fiberdat.iFrM[gp]);
  static LINALG::Matrix<3,3> FrdFgdrhoM(true);
  static LINALG::Matrix<3,3> dAgrdrhoM(true);
  FrdFgdrhoM.MultiplyNN(1.0,FrM,dFgdrhoM,0.0);
  tmp.MultiplyNN(1.0,FrdFgdrhoM,fiberdat.AM,0.0);
  dAgrdrhoM.MultiplyNT(1.0/CinM.Dot(fiberdat.AM),tmp,FinM,0.0);
  dAgrdrhoM.MultiplyNT(1.0/CinM.Dot(fiberdat.AM),FinM,tmp,1.0);
  dAgrdrhoM.Update(2.0*CinAFgTM.Dot(diFgdrhoM)/CinM.Dot(fiberdat.AM),AgrM,1.0);

  dsigdCedrhoM.Update(2.0*(ddPIIe(0)*CeM.Dot(AgrM)+dPIe(0)),dAgrdrhoM,1.0);

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template< typename ForceAnalytical >
void MAT::ELASTIC::RemodelFiber::EvaluateDerivativesCauchyRemodel(LINALG::Matrix<3,3> const& CM,
                                                                  LINALG::Matrix<3,3> const& iFgM,
                                                                  FiberData const& fiberdat,
                                                                  int const gp,
                                                                  ForceAnalytical const eleGID,
                                                                  double & dsigdlambr,
                                                                  LINALG::Matrix<3,3> & dsigdCedlambrM) const
{
  // clear some variables
  dsigdlambr = 0.0;
  dsigdCedlambrM.Clear();

  static LINALG::FADMatrix<3,3> iFinM_fad(true);
  static LINALG::FADMatrix<3,3> iFgM_fad(true);
  static LINALG::FADMatrix<3,3> iFrM_fad(true);
  static LINALG::FADMatrix<3,3> CM_fad(true);
  static LINALG::FADMatrix<3,3> AM_fad(true);
  iFgM_fad = iFgM;
  iFrM_fad = fiberdat.iFrM[gp];
  iFrM_fad.diff(0,9);
  CM_fad = CM;
  AM_fad = fiberdat.AM;
  iFinM_fad.MultiplyNN(1.0,iFgM_fad,iFrM_fad,0.0);

  FAD sig_fad = 0.0;
  EvaluateLocalCauchyStress(CM_fad,iFinM_fad,AM_fad,fiberdat.fiber,eleGID,sig_fad);

  LINALG::Matrix<3,3> dsigdiFrM(true);
  FirstDerivToMatrix(sig_fad,dsigdiFrM);
  dsigdlambr = dsigdiFrM.Dot(fiberdat.diFrdlambrM[gp]);


  static LINALG::FADMatrix<3,3> dsigdCeM_fad(true);
  static LINALG::Matrix<6,9> dsigdCediFr(true);
  EvaluatedsigdCe(CM_fad,iFgM_fad,iFrM_fad,AM_fad,fiberdat.fiber,eleGID,dsigdCeM_fad);

  for(int i=0;i<3;++i)
    for(int j=0;j<9;++j)
      dsigdCediFr(i,j) = dsigdCeM_fad(i,i).dx(j);
  for(int j=0;j<9;++j)
    dsigdCediFr(3,j) = 0.5*(dsigdCeM_fad(0,1).dx(j) + dsigdCeM_fad(1,0).dx(j));
  for(int j=0;j<9;++j)
    dsigdCediFr(4,j) = 0.5*(dsigdCeM_fad(1,2).dx(j) + dsigdCeM_fad(2,1).dx(j));
  for(int j=0;j<9;++j)
    dsigdCediFr(5,j) = 0.5*(dsigdCeM_fad(0,2).dx(j) + dsigdCeM_fad(2,0).dx(j));

  static LINALG::Matrix<9,1> diFrdlambr9x1(true);
  static LINALG::Matrix<6,1> tmp6x1(true);
  Matrix3x3toVector9x1(fiberdat.diFrdlambrM[gp],diFrdlambr9x1);
  tmp6x1.MultiplyNN(1.0,dsigdCediFr,diFrdlambr9x1,0.0);
  StressVoigtNotationVectorToMatrix(tmp6x1,dsigdCedlambrM);

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::RemodelFiber::EvaluateDerivativesCauchyRemodel(LINALG::Matrix<3,3> const& CM,
                                                                  LINALG::Matrix<3,3> const& iFgM,
                                                                  FiberData const& fiberdat,
                                                                  int const gp,
                                                                  int const eleGID,
                                                                  double & dsigdlambr,
                                                                  LINALG::Matrix<3,3> & dsigdCedlambrM) const
{
  // clear some variables
  dsigdlambr = 0.0;
  dsigdCedlambrM.Clear();

  static LINALG::Matrix<3,3> CAFinTM(true);
  static LINALG::Matrix<3,3> tmp1(true);
  static LINALG::Matrix<3,3> CeM(true);
  static LINALG::Matrix<3,3> iFinM(true);
  static LINALG::Matrix<3,3> FgM(true);
  static LINALG::Matrix<3,3> FrM(true);
  static LINALG::Matrix<3,3> FinM(true);
  static LINALG::Matrix<3,3> AFinTM(true);
  static LINALG::Matrix<3,3> FinAFinTM(true);
  static LINALG::Matrix<3,3> CinM(true);
  static LINALG::Matrix<3,3> CinAM(true);
  FrM.Invert(fiberdat.iFrM[gp]);
  FgM.Invert(iFgM);
  iFinM.MultiplyNN(1.0,iFgM,fiberdat.iFrM[gp],0.0);
  FinM.Invert(iFinM);
  AFinTM.MultiplyNT(1.0,fiberdat.AM,FinM,0.0);
  FinAFinTM.MultiplyNN(1.0,FinM,AFinTM,0.0);
  CAFinTM.MultiplyNN(1.0,CM,AFinTM,0.0);
  CinM.MultiplyTN(1.0,FinM,FinM,0.0);
  CinAM.MultiplyTN(1.0,CinM,fiberdat.AM,0.0);
  tmp1.MultiplyNN(1.0,CM,iFinM,0.0);
  CeM.MultiplyTN(1.0,iFinM,tmp1,0.0);

  // temporary pointers to check the type of the remodelfiber (active or passive)
  Teuchos::RCP<MAT::ELASTIC::CoupAnisoExpo> t1(nullptr);
  Teuchos::RCP<MAT::ELASTIC::CoupAnisoExpoActive> t2(nullptr);
  static LINALG::Matrix<2,1> dPIe(true);
  static LINALG::Matrix<3,1> ddPIIe(true);
  static LINALG::Matrix<4,1> dddPIIIe(true);
  if( (t1 = Teuchos::rcp_dynamic_cast<MAT::ELASTIC::CoupAnisoExpo>(fiberdat.fiber)).getRawPtr() )
    t1->GetDerivativesAniso(dPIe,ddPIIe,dddPIIIe,CeM,eleGID);
  else if( (t2 = Teuchos::rcp_dynamic_cast<MAT::ELASTIC::CoupAnisoExpoActive>(fiberdat.fiber)).getRawPtr() )
    t2->GetDerivativesAniso(dPIe,ddPIIe,dddPIIIe,CeM,eleGID);

  static LINALG::Matrix<3,3> dsigdiFrM(true);
  dsigdiFrM.MultiplyTN(4.0*ddPIIe(0)*CM.Dot(fiberdat.AM)/(CinM.Dot(fiberdat.AM)*CinM.Dot(fiberdat.AM)),iFgM,CAFinTM,0.0);
  dsigdiFrM.MultiplyTN(4.0*dPIe(0)*CM.Dot(fiberdat.AM)/(CinM.Dot(fiberdat.AM)*CinM.Dot(fiberdat.AM)),FrM,FinAFinTM,1.0);
  dsigdlambr = dsigdiFrM.Dot(fiberdat.diFrdlambrM[gp]);


  static LINALG::Matrix<3,3> AgrM(true);
  static LINALG::Matrix<3,3> tmp(true);
  static LINALG::Matrix<3,3> FrTFinAFinTM(true);
  static LINALG::Matrix<3,3> iFgTCAFinTM(true);
  tmp.MultiplyNN(1.0,FinM,fiberdat.AM,0.0);
  AgrM.MultiplyNT(1.0/CinM.Dot(fiberdat.AM),tmp,FinM,0.0);
  FrTFinAFinTM.MultiplyTN(1.0,FrM,FinAFinTM,0.0);
  iFgTCAFinTM.MultiplyTN(1.0,iFgM,CAFinTM,0.0);

  dsigdCedlambrM.Update(4.0/CinM.Dot(fiberdat.AM)*(dddPIIIe(0)*iFgTCAFinTM.Dot(fiberdat.diFrdlambrM[gp])*CeM.Dot(AgrM) +
                                     ddPIIe(0)*CM.Dot(fiberdat.AM)*FrTFinAFinTM.Dot(fiberdat.diFrdlambrM[gp])/CinM.Dot(fiberdat.AM) +
                                     ddPIIe(0)*iFgTCAFinTM.Dot(fiberdat.diFrdlambrM[gp])),AgrM,0.0);

  static LINALG::Matrix<3,3> dFrdlambrFgM(true);
  static LINALG::Matrix<3,3> dAgrdlambrM(true);
  dFrdlambrFgM.MultiplyNN(1.0,fiberdat.dFrdlambrM[gp],FgM,0.0);
  tmp.MultiplyNN(1.0,dFrdlambrFgM,fiberdat.AM,0.0);
  dAgrdlambrM.MultiplyNT(1.0/CinM.Dot(fiberdat.AM),tmp,FinM,0.0);
  dAgrdlambrM.MultiplyNT(1.0/CinM.Dot(fiberdat.AM),FinM,tmp,1.0);
  dAgrdlambrM.Update(2.0*FrTFinAFinTM.Dot(fiberdat.diFrdlambrM[gp])/CinM.Dot(fiberdat.AM),AgrM,1.0);

  dsigdCedlambrM.Update(2.0*(ddPIIe(0)*CeM.Dot(AgrM)+dPIe(0)),dAgrdlambrM,1.0);

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template< typename T, typename ForceAnalytical >
void MAT::ELASTIC::RemodelFiber::EvaluatedsigdC(LINALG::TMatrix<T,3,3> const& CM,
                                                LINALG::TMatrix<T,3,3> const& iFinM,
                                                LINALG::TMatrix<T,3,3> const& AM,
                                                Teuchos::RCP<MAT::ELASTIC::Summand> const fiber,
                                                ForceAnalytical const eleGID,
                                                LINALG::TMatrix<T,3,3>& dsigdC) const
{
  // clear some variables
  dsigdC.Clear();

  static LINALG::FADMatrix<3,3> iFinM_fad(true);
  iFinM_fad = iFinM;
  static LINALG::FADMatrix<3,3> AM_fad(true);
  AM_fad = AM;

  static LINALG::FADMatrix<3,3> CM_fad(true);
  CM_fad = CM;
  CM_fad.diff(0,9);

  FAD sig_fad = 0.0;
  EvaluateLocalCauchyStress(CM_fad,iFinM_fad,AM_fad,fiber,eleGID,sig_fad);

  LINALG::TMatrix<T,3,3> tmp(true);
  FirstDerivToMatrix(sig_fad,dsigdC);
std::cout <<dsigdC <<std::endl;
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template< typename T >
void MAT::ELASTIC::RemodelFiber::EvaluatedsigdC(LINALG::TMatrix<T,3,3> const& CM,
                                                LINALG::TMatrix<T,3,3> const& iFinM,
                                                LINALG::TMatrix<T,3,3> const& AM,
                                                Teuchos::RCP<MAT::ELASTIC::Summand> const fiber,
                                                int const eleGID,
                                                LINALG::TMatrix<T,3,3>& dsigdC) const
{
  // clear some variables
  dsigdC.Clear();

  // temporary pointers to check the type of the remodelfiber (active or passive)
  Teuchos::RCP<MAT::ELASTIC::CoupAnisoExpo> t1(nullptr);
  Teuchos::RCP<MAT::ELASTIC::CoupAnisoExpoActive> t2(nullptr);

  static LINALG::TMatrix<T,3,3> FinM(true);
  FinM.Invert(iFinM);
  static LINALG::TMatrix<T,3,3> CinM(true);
  CinM.MultiplyTN(1.0,FinM,FinM,0.0);
  static LINALG::TMatrix<T,3,3> CeM(true);
  static LINALG::TMatrix<T,3,3> tmp(true);
  tmp.MultiplyNN(1.0,CM,iFinM,0.0);
  CeM.MultiplyTN(1.0,iFinM,tmp,0.0);

  static LINALG::Matrix<2,1> dPIe(true);
  static LINALG::Matrix<3,1> ddPIIe(true);
  static LINALG::Matrix<4,1> dddPIIIe(true);
  if( (t1 = Teuchos::rcp_dynamic_cast<MAT::ELASTIC::CoupAnisoExpo>(fiber)).getRawPtr() )
    t1->GetDerivativesAniso(dPIe,ddPIIe,dddPIIIe,CeM,eleGID);
  else if( (t2 = Teuchos::rcp_dynamic_cast<MAT::ELASTIC::CoupAnisoExpoActive>(fiber)).getRawPtr() )
    t2->GetDerivativesAniso(dPIe,ddPIIe,dddPIIIe,CeM,eleGID);

  dsigdC.Update(2.0/CinM.Dot(AM)*(ddPIIe(0)*CM.Dot(AM)/CinM.Dot(AM) + dPIe(0)),AM,0.0);

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::RemodelFiber::EvaluateEvolutionEquation(double & rg,
                                                           double & rr,
                                                           LINALG::Matrix<3,3> const& CM,
                                                           LINALG::Matrix<3,3> const& iFgM,
                                                           double const& dt,
                                                           FiberData const& fiberdat,
                                                           int const gp,
                                                           int const eleGID) const
{
  static LINALG::Matrix<3,3> iFinM(true);
  iFinM.MultiplyNN(1.0,iFgM,fiberdat.iFrM[gp],0.0);
  double sig = 0.0;
  EvaluateLocalCauchyStress(CM,iFinM,fiberdat.AM,fiberdat.fiber,eleGID,sig);

  // Growth evolution equation
  fiberdat.growth->EvaluateFunc(rg,sig,fiberdat.cur_rho[gp],fiberdat.last_rho[gp],dt,eleGID);

  static LINALG::Matrix<3,3> dsigdCe(true);
  EvaluatedsigdCe(CM,iFgM,fiberdat.iFrM[gp],fiberdat.AM,fiberdat.fiber,eleGID,dsigdCe);

  static LINALG::Matrix<3,3> YM(true);
  static LINALG::Matrix<3,3> FrdotiFrM(true);
  static LINALG::Matrix<3,3> tmp(true);
  static LINALG::Matrix<3,3> CeM(true);
  tmp.MultiplyNN(1.0,CM,iFinM,0.0);
  CeM.MultiplyTN(1.0,iFinM,tmp,0.0);
  FrdotiFrM.MultiplyNN(1.0,fiberdat.FrdotM[gp],fiberdat.iFrM[gp],0.0);
  YM.MultiplyNN(1.0,CeM,FrdotiFrM,0.0);

  // Remodel evolution equation
  fiberdat.remodel->EvaluateFunc(rr,sig,YM,dsigdCe,eleGID);

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::RemodelFiber::EvaluateDerivativeEvolutionEquation(double & dWidrhoi,
                                                                     double & dWidrhoj,
                                                                     double & dWdlambr,
                                                                     double & dEdrho,
                                                                     double & dEdlambr,
                                                                     LINALG::Matrix<3,3> const& CM,
                                                                     LINALG::Matrix<3,3> const& iFgM,
                                                                     LINALG::Matrix<3,3> const& dFgdrhoM,
                                                                     LINALG::Matrix<3,3> const& diFgdrhoM,
                                                                     double const& dt,
                                                                     FiberData const& fiberdat,
                                                                     int const gp,
                                                                     int const eleGID) const
{
  // Derivative of growth evolution eq. w.r.t. the mass density
  static LINALG::Matrix<3,3> iFinM(true);
  iFinM.MultiplyNN(1.0,iFgM,fiberdat.iFrM[gp],0.0);
  double sig = 0.0;
  EvaluateLocalCauchyStress(CM,iFinM,fiberdat.AM,fiberdat.fiber,eleGID,sig);

  double dsigdrho = 0.0;
  static LINALG::Matrix<3,3> dsigdCedrhoM(true);
  EvaluateDerivativesCauchyGrowth(CM,iFgM,dFgdrhoM,diFgdrhoM,fiberdat,gp,eleGID,dsigdrho,dsigdCedrhoM);

  fiberdat.growth->EvaluatedFuncidrhoi(dWidrhoi,sig,dsigdrho,fiberdat.cur_rho[gp],dt,eleGID);
  fiberdat.growth->EvaluatedFuncidrhoj(dWidrhoj,sig,dsigdrho,fiberdat.cur_rho[gp],dt,eleGID);

  // Derivative of growth evolution eq. w.r.t. the inelastic remodel stretch
  static LINALG::Matrix<3,3> dsigdCedlambrM(true);
  double dsigdlambr = 0.0;
  EvaluateDerivativesCauchyRemodel(CM,iFgM,fiberdat,gp,eleGID,dsigdlambr,dsigdCedlambrM);

  fiberdat.growth->EvaluatedFuncidlambr(dWdlambr,dsigdlambr,fiberdat.cur_rho[gp],dt,eleGID);

  // Derivative of remodel evolution eq. w.r.t. the inelastic remodel stretch
  static LINALG::Matrix<3,3> tmp1(true);
  static LINALG::Matrix<3,3> tmp2(true);
  static LINALG::Matrix<3,3> FrdotiFrM(true);
  static LINALG::Matrix<3,3> CeM(true);
  static LINALG::Matrix<3,3> YM(true);
  tmp1.MultiplyNN(1.0,CM,iFinM,0.0);
  CeM.MultiplyTN(1.0,iFinM,tmp1,0.0);
  FrdotiFrM.MultiplyNN(1.0,fiberdat.FrdotM[gp],fiberdat.iFrM[gp],0.0);
  YM.MultiplyNN(1.0,CeM,FrdotiFrM,0.0);

  static LINALG::Matrix<3,3> dYdlambrM(true);
  static LINALG::Matrix<3,3> dFrdotdlambriFrM(true);
  static LINALG::Matrix<3,3> FrdotdiFrdlambrM(true);
  tmp1.MultiplyTN(1.0,iFinM,CM,0.0);
  tmp2.MultiplyNN(1.0,tmp1,iFgM,0.0);
  tmp1.MultiplyNN(1.0,tmp2,fiberdat.diFrdlambrM[gp],0.0);
  dYdlambrM.MultiplyNN(1.0,tmp1,FrdotiFrM,0.0);
  dYdlambrM.MultiplyTN(1.0,tmp1,FrdotiFrM,1.0);
  dFrdotdlambriFrM.MultiplyNN(1.0,fiberdat.dFrdotdlambrM[gp],fiberdat.iFrM[gp],0.0);
  dYdlambrM.MultiplyNN(1.0,CeM,dFrdotdlambriFrM,1.0);
  FrdotdiFrdlambrM.MultiplyNN(1.0,fiberdat.FrdotM[gp],fiberdat.diFrdlambrM[gp],0.0);
  dYdlambrM.MultiplyNN(1.0,CeM,FrdotdiFrdlambrM,1.0);

  static LINALG::Matrix<3,3> dsigdCeM(true);
  EvaluatedsigdCe(CM,iFgM,fiberdat.iFrM[gp],fiberdat.AM,fiberdat.fiber,eleGID,dsigdCeM);

  fiberdat.remodel->EvaluatedFuncidlambri(dEdlambr,sig,dsigdlambr,YM,dYdlambrM,dsigdCeM,dsigdCedlambrM,eleGID);

  // Derivative of remodel evolution eq. w.r.t. the mass density
  static LINALG::Matrix<3,3> dYdrhoM(true);
  tmp1.MultiplyTT(1.0,fiberdat.iFrM[gp],diFgdrhoM,0.0);
  tmp2.MultiplyNN(1.0,tmp1,CM,0.0);
  tmp1.MultiplyNN(1.0,tmp2,iFinM,0.0);
  dYdrhoM.MultiplyNN(1.0,tmp1,FrdotiFrM,0.0);
  dYdrhoM.MultiplyTN(1.0,tmp1,FrdotiFrM,1.0);

  fiberdat.remodel->EvaluatedFuncidrho(dEdrho,sig,dsigdrho,YM,dYdrhoM,dsigdCeM,dsigdCedrhoM,eleGID);


  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::RemodelFiber::EvaluatedEvolutionEquationdC(LINALG::Matrix<1,6> & dWdC,
                                                              LINALG::Matrix<1,6> & dEdC,
                                                              LINALG::Matrix<3,3> const& CM,
                                                              LINALG::Matrix<3,3> const& iFgM,
                                                              double const& dt,
                                                              FiberData const& fiberdat,
                                                              int const k,
                                                              int const gp,
                                                              int const eleGID)
{
  // Growth law
  static LINALG::Matrix<3,3> iFinM(true);
  iFinM.MultiplyNN(1.0,iFgM,fiberdat.iFrM[gp],0.0);

  static LINALG::Matrix<3,3> dsigdC(true);
  EvaluatedsigdC(CM,iFinM,fiberdat.AM,fiberdat.fiber,eleGID,dsigdC);
  static LINALG::Matrix<6,1> dsigdCv(true);
  MatrixtoStressVoigtNotationVector(dsigdC,dsigdCv);

  fiberdat.growth->EvaluatedFuncidC(dWdC,dsigdCv,fiberdat.cur_rho[gp],dt,eleGID);


  // Remodel law
  double sig = 0.0;
  EvaluateLocalCauchyStress(CM,iFinM,fiberdat.AM,fiberdat.fiber,eleGID,sig);
  cauchystress_[k][gp] = sig;
  static LINALG::Matrix<6,6> dsigdCedC(true);
  EvaluatedsigdCedC(CM,iFgM,fiberdat.iFrM[gp],fiberdat.AM,fiberdat.fiber,eleGID,dsigdCedC);
  static LINALG::Matrix<3,3> dsigdCeM(true);
  static LINALG::Matrix<9,1> dsigdCe9x1(true);
  EvaluatedsigdCe(CM,iFgM,fiberdat.iFrM[gp],fiberdat.AM,fiberdat.fiber,eleGID,dsigdCeM);
  Matrix3x3toVector9x1(dsigdCeM,dsigdCe9x1);

  static LINALG::Matrix<3,3> tmp(true);
  static LINALG::Matrix<3,3> FrdotiFrM(true);
  static LINALG::Matrix<3,3> CeM(true);
  static LINALG::Matrix<3,3> YM(true);
  static LINALG::Matrix<6,1> Y_strain(true);
  tmp.MultiplyNN(1.0,CM,iFinM,0.0);
  CeM.MultiplyTN(1.0,iFinM,tmp,0.0);
  FrdotiFrM.MultiplyNN(1.0,fiberdat.FrdotM[gp],fiberdat.iFrM[gp],0.0);
  YM.MultiplyNN(1.0,CeM,FrdotiFrM,0.0);
  MatrixtoStrainVoigtNotationVector(YM,Y_strain);

  LINALG::Matrix<9,6> dYdC(true);
  static LINALG::Matrix<3,3> iFinTM(true);
  static LINALG::Matrix<3,3> iFrTFrdotTiFinTM(true);
  iFinTM.UpdateT(1.0,iFinM,0.0);
  iFrTFrdotTiFinTM.MultiplyTT(1.0,FrdotiFrM,iFinM,0.0);
  MAT::AddLeftNonSymmetricHolzapfelProduct(dYdC,iFinTM,iFrTFrdotTiFinTM,0.5);

  fiberdat.remodel->EvaluatedFuncidC(dEdC,sig,dsigdCv,Y_strain,dYdC,dsigdCe9x1,dsigdCedC,eleGID);

  return;
}





/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::RemodelFiber::EvaluatedEvolutionEquationdt(double & drhodt,
                                                              double & dlambrdt,
                                                              LINALG::Matrix<3,3> const& CM,
                                                              LINALG::Matrix<3,3> const& iFgM,
                                                              FiberData const& fiberdat,
                                                              int const k,
                                                              int const gp,
                                                              int const eleGID)
{
  // Time derivative of the mass density
  static LINALG::Matrix<3,3> iFinM(true);
  iFinM.MultiplyNN(1.0,iFgM,fiberdat.iFrM[gp],0.0);
  double sig = 0.0;
  EvaluateLocalCauchyStress(CM,iFinM,fiberdat.AM,fiberdat.fiber,eleGID,sig);
  cauchystress_[k][gp] = sig;

  fiberdat.growth->Evaluatedrhodt(drhodt,sig,fiberdat.cur_rho[gp],eleGID);

  // Time derivative of the remodel stretch
  static LINALG::Matrix<3,3> tmp(true);
  static LINALG::Matrix<3,3> FrdotredM(true);
  static LINALG::Matrix<3,3> CeM(true);
  static LINALG::Matrix<3,3> YredM(true);
  tmp.MultiplyNN(1.0,CM,iFinM,0.0);
  CeM.MultiplyTN(1.0,iFinM,tmp,0.0);
  static LINALG::Matrix<3,3> dsigdCeM(true);
  EvaluatedsigdCe(CM,iFgM,fiberdat.iFrM[gp],fiberdat.AM,fiberdat.fiber,eleGID,dsigdCeM);

  FrdotredM.Update(1.0/fiberdat.G,fiberdat.AM,0.0);
  FrdotredM.Update(-0.5*std::pow(fiberdat.cur_lambr[gp],-3./2.)*std::pow(fiberdat.G,0.5),fiberdat.AM_orth,1.0);
  tmp.MultiplyNN(1.0,CeM,FrdotredM,0.0);
  YredM.MultiplyNN(1.0,tmp,fiberdat.iFrM[gp],0.0);

  fiberdat.remodel->Evaluatedlambrdt(dlambrdt,sig,YredM,dsigdCeM,eleGID);

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template < typename ForceAnalytical >
void MAT::ELASTIC::RemodelFiber::EvaluateDerivatives2ndPiolaKirchhoffGrowthRemodel(LINALG::Matrix<6,1> & dSidrhoi,
                                                                                   LINALG::Matrix<6,1> & dSidrhoj,
                                                                                   LINALG::Matrix<6,1> & dSdlambr,
                                                                                   LINALG::Matrix<3,3> const& CM,
                                                                                   LINALG::Matrix<3,3> const& iFgM,
                                                                                   LINALG::Matrix<3,3> const& diFgdrhoM,
                                                                                   FiberData const& fiberdat,
                                                                                   int const gp,
                                                                                   ForceAnalytical const eleGID) const
{
  // Derivative w.r.t. the mass density
  static LINALG::FADMatrix<3,3> CM_fad(true);
  CM_fad = CM;
  static LINALG::FADMatrix<3,3> iFgM_fad(true);
  iFgM_fad = iFgM;
  iFgM_fad.diff(0,18);
  static LINALG::FADMatrix<3,3> iFrM_fad(true);
  iFrM_fad = fiberdat.iFrM[gp];
  iFrM_fad.diff(9,18);
  static LINALG::FADMatrix<3,3> iFinM_fad(true);
  iFinM_fad.MultiplyNN(1.0,iFgM_fad,iFrM_fad,0.0);
  static LINALG::FADMatrix<3,3> AM_fad(true);
  AM_fad = fiberdat.AM;

  // temporary pointers to check the type of the remodelfiber (active or passive)
  Teuchos::RCP<MAT::ELASTIC::CoupAnisoExpo> t1(nullptr);
  Teuchos::RCP<MAT::ELASTIC::CoupAnisoExpoActive> t2(nullptr);

  static LINALG::FADMatrix<3,3> firstderivM_fad(true);
  static LINALG::Matrix<6,1> Sactv(true);
  static LINALG::Matrix<6,6> cmatact(true);
  static LINALG::FADMatrix<3,3> S_fad(true);
  S_fad.Clear();
  if( (t1 = Teuchos::rcp_dynamic_cast<MAT::ELASTIC::CoupAnisoExpo>(fiberdat.fiber)).getRawPtr() ) {
    DerivdC(CM_fad,iFinM_fad,AM_fad,*t1,eleGID,firstderivM_fad);
  }
  else if( (t2 = Teuchos::rcp_dynamic_cast<MAT::ELASTIC::CoupAnisoExpoActive>(fiberdat.fiber)).getRawPtr() ) {
    DerivdC(CM_fad,iFinM_fad,AM_fad,*t2,eleGID,firstderivM_fad);
    t2->EvaluateActiveStressCmatAniso(CM,cmatact,Sactv,eleGID);
  }

  S_fad.Update(2.0*fiberdat.cur_rho[gp],firstderivM_fad,1.0);

  static LINALG::Matrix<6,9> dSdiFg(true);
  for(int i=0;i<3;++i)
    for(int j=0;j<9;++j)
      dSdiFg(i,j) = S_fad(i,i).dx(j);
  for(int j=0;j<9;++j)
    dSdiFg(3,j) = 0.5*(S_fad(0,1).dx(j) + S_fad(1,0).dx(j));
  for(int j=0;j<9;++j)
    dSdiFg(4,j) = 0.5*(S_fad(1,2).dx(j) + S_fad(2,1).dx(j));
  for(int j=0;j<9;++j)
    dSdiFg(5,j) = 0.5*(S_fad(0,2).dx(j) + S_fad(2,0).dx(j));

  static LINALG::Matrix<9,1> diFgdrho9x1(true);
  Matrix3x3toVector9x1(diFgdrhoM,diFgdrho9x1);
  static LINALG::Matrix<3,3> firstderivM(true);
  static LINALG::Matrix<6,1> firstderivv(true);
  firstderivM = firstderivM_fad.ConverttoDouble();
  MatrixtoStressVoigtNotationVector(firstderivM,firstderivv);
  dSidrhoj.MultiplyNN(1.0,dSdiFg,diFgdrho9x1,0.0);
  dSidrhoi.Update(1.0,dSidrhoj,0.0);
  dSidrhoi.Update(2.0,firstderivv,1.0);
  dSidrhoi.Update(1.0,Sactv,1.0);


  // Derivative w.r.t. the inelastic remodel stretch
  static LINALG::Matrix<6,9> dSdiFr(true);
  for(int i=0;i<3;++i)
    for(int j=0;j<9;++j)
      dSdiFr(i,j) = S_fad(i,i).dx(j+9);
  for(int j=0;j<9;++j)
    dSdiFr(3,j) = 0.5*(S_fad(0,1).dx(j+9) + S_fad(1,0).dx(j+9));
  for(int j=0;j<9;++j)
    dSdiFr(4,j) = 0.5*(S_fad(1,2).dx(j+9) + S_fad(2,1).dx(j+9));
  for(int j=0;j<9;++j)
    dSdiFr(5,j) = 0.5*(S_fad(0,2).dx(j+9) + S_fad(2,0).dx(j+9));

  static LINALG::Matrix<9,1> diFrdlambr9x1(true);
  Matrix3x3toVector9x1(fiberdat.diFrdlambrM[gp],diFrdlambr9x1);
  dSdlambr.MultiplyNN(1.0,dSdiFr,diFrdlambr9x1,0.0);

  return;
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::RemodelFiber::EvaluateDerivatives2ndPiolaKirchhoffGrowthRemodel(LINALG::Matrix<6,1> & dSidrhoi,
                                                                                   LINALG::Matrix<6,1> & dSidrhoj,
                                                                                   LINALG::Matrix<6,1> & dSdlambr,
                                                                                   LINALG::Matrix<3,3> const& CM,
                                                                                   LINALG::Matrix<3,3> const& iFgM,
                                                                                   LINALG::Matrix<3,3> const& diFgdrhoM,
                                                                                   FiberData const& fiberdat,
                                                                                   int const gp,
                                                                                   int const eleGID) const
{
  // Clear some variables
  dSidrhoi.Clear();
  dSidrhoj.Clear();
  dSdlambr.Clear();

  // Derivative w.r.t. the mass density
  static LINALG::Matrix<3,3> tmp(true);
  static LINALG::Matrix<3,3> CeM(true);
  static LINALG::Matrix<3,3> iFinM(true);
  static LINALG::Matrix<3,3> FgM(true);
  static LINALG::Matrix<3,3> FinM(true);
  static LINALG::Matrix<3,3> CinM(true);
  static LINALG::Matrix<3,3> CAFgTM(true);
  static LINALG::Matrix<3,3> CinAFgTM(true);
  FgM.Invert(iFgM);
  iFinM.MultiplyNN(1.0,iFgM,fiberdat.iFrM[gp],0.0);
  FinM.Invert(iFinM);
  CinM.MultiplyTN(1.0,FinM,FinM,0.0);
  tmp.MultiplyTN(1.0,CinM,fiberdat.AM,0.0);
  CinAFgTM.MultiplyNT(1.0,tmp,FgM,0.0);
  tmp.MultiplyTN(1.0,CM,fiberdat.AM,0.0);
  CAFgTM.MultiplyNT(1.0,tmp,FgM,0.0);
  tmp.MultiplyNN(1.0,CM,iFinM,0.0);
  CeM.MultiplyTN(1.0,iFinM,tmp,0.0);

  // temporary pointers to check the type of the remodelfiber (active or passive)
  Teuchos::RCP<MAT::ELASTIC::CoupAnisoExpo> t1(nullptr);
  Teuchos::RCP<MAT::ELASTIC::CoupAnisoExpoActive> t2(nullptr);

  static LINALG::Matrix<2,1> dPIe(true);
  static LINALG::Matrix<3,1> ddPIIe(true);
  static LINALG::Matrix<4,1> dddPIIIe(true);
  static LINALG::Matrix<6,1> stressactv(true);
  static LINALG::Matrix<6,6> cmatact(true);
  if( (t1 = Teuchos::rcp_dynamic_cast<MAT::ELASTIC::CoupAnisoExpo>(fiberdat.fiber)).getRawPtr() )
    t1->GetDerivativesAniso(dPIe,ddPIIe,dddPIIIe,CeM,eleGID);
  else if( (t2 = Teuchos::rcp_dynamic_cast<MAT::ELASTIC::CoupAnisoExpoActive>(fiberdat.fiber)).getRawPtr() ) {
    t2->GetDerivativesAniso(dPIe,ddPIIe,dddPIIIe,CeM,eleGID);
    t2->EvaluateActiveStressCmatAniso(CM,cmatact,stressactv,eleGID);
    dSidrhoi.Update(1.0,stressactv,0.0);
  }

  static LINALG::Matrix<6,1> Av(true);
  MatrixtoStressVoigtNotationVector(fiberdat.AM,Av);
  dSidrhoj.Update(4.0*fiberdat.cur_rho[gp]/(CinM.Dot(fiberdat.AM)*CinM.Dot(fiberdat.AM))*
                  (ddPIIe(0)*CAFgTM.Dot(diFgdrhoM) + dPIe(0)*CinAFgTM.Dot(diFgdrhoM)),Av,0.0);
  dSidrhoi.Update(1.0,dSidrhoj,1.0);
  dSidrhoi.Update(2.0*dPIe(0)/CinM.Dot(fiberdat.AM),Av,1.0);


  // Derivative w.r.t. the inelastic remodel stretch
  static LINALG::Matrix<3,3> iFgTCAFinTM(true);
  static LINALG::Matrix<3,3> FrM(true);
  static LINALG::Matrix<3,3> AFinTM(true);
  static LINALG::Matrix<3,3> FrTFinAFinTM(true);
  FrM.Invert(fiberdat.iFrM[gp]);
  AFinTM.MultiplyNT(1.0,fiberdat.AM,FinM,0.0);
  tmp.MultiplyNN(1.0,FinM,AFinTM,0.0);
  FrTFinAFinTM.MultiplyTN(1.0,FrM,tmp,0.0);
  tmp.MultiplyNN(1.0,CM,AFinTM,0.0);
  iFgTCAFinTM.MultiplyTN(1.0,iFgM,tmp,0.0);

  dSdlambr.Update(4.0*fiberdat.cur_rho[gp]/(CinM.Dot(fiberdat.AM)*CinM.Dot(fiberdat.AM))*
                  (ddPIIe(0)*iFgTCAFinTM.Dot(fiberdat.diFrdlambrM[gp]) + dPIe(0)*FrTFinAFinTM.Dot(fiberdat.diFrdlambrM[gp])),Av,0.0);

  return;
}


/*---------------------------------------------------------------------*
 | return names of visualization data (public)                         |
 *---------------------------------------------------------------------*/
void MAT::ELASTIC::RemodelFiber::VisNames(std::map<std::string,int>& names, unsigned int p)
{
  std::string inelastic_defgrd = "lambda_r";
  std::string result_inelastic_defgrad;
  for (unsigned int k=0; k<potsumfiber_.size(); ++k) {
    std::stringstream sstm;
    sstm << inelastic_defgrd <<"_" << p <<"_" << k;
    result_inelastic_defgrad = sstm.str();

    names[result_inelastic_defgrad] = 1;
  }


  std::string fiber_cauchy_stress = "fiber_cauchy_stress";
  std::string result_fiber_cauchy_stress;
  for (unsigned int k=0; k<potsumfiber_.size(); ++k) {
    std::stringstream sstm;
    sstm << fiber_cauchy_stress <<"_" << p <<"_" << k;
    result_fiber_cauchy_stress = sstm.str();

    names[result_fiber_cauchy_stress] = 1;
  }


  std::string cur_rho_col = "cur_rho_col";
  std::string result_cur_rho_col;
  for (unsigned int k=0; k<potsumfiber_.size(); ++k) {
    std::stringstream sstm;
    sstm << cur_rho_col <<"_" << p <<"_" << k;
    result_cur_rho_col = sstm.str();

    names[result_cur_rho_col] = 1;
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
  if ((name == "lambda_r_0_0") || (name == "lambda_r_1_0")) {
    if (data.size()!= 1) dserror("size mismatch");
    for(unsigned gp=0; gp<potsumfiber_[0]->last_lambr.size(); ++gp) {
      data[0] += potsumfiber_[0]->last_lambr[gp];
    }
    data[0] = data[0]/potsumfiber_[0]->last_lambr.size();

    return true;
  }
  if ((name == "lambda_r_0_1") || (name == "lambda_r_1_1")) {
    if (data.size()!= 1) dserror("size mismatch");
    for(unsigned gp=0; gp<potsumfiber_[1]->last_lambr.size(); ++gp) {
      data[0] += potsumfiber_[1]->last_lambr[gp];
    }
    data[0] = data[0]/potsumfiber_[1]->last_lambr.size();

    return true;
  }
  if ((name == "lambda_r_0_2") || (name == "lambda_r_1_2")) {
    if (data.size()!= 1) dserror("size mismatch");
    for(unsigned gp=0; gp<potsumfiber_[2]->last_lambr.size(); ++gp) {
      data[0] += potsumfiber_[2]->last_lambr[gp];
    }
    data[0] = data[0]/potsumfiber_[2]->last_lambr.size();

    return true;
  }
  if ((name == "lambda_r_0_3") || (name == "lambda_r_1_3")) {
    if (data.size()!= 1) dserror("size mismatch");
    for(unsigned gp=0; gp<potsumfiber_[3]->last_lambr.size(); ++gp) {
      data[0] += potsumfiber_[3]->last_lambr[gp];
    }
    data[0] = data[0]/potsumfiber_[3]->last_lambr.size();

    return true;
  }


  if ((name == "fiber_cauchy_stress_0_0") || (name == "fiber_cauchy_stress_1_0")) {
    if (data.size()!= 1) dserror("size mismatch");
    for(unsigned gp=0; gp<cauchystress_[0].size(); ++gp) {
      data[0] += cauchystress_[0][gp];
    }
    data[0] = data[0]/cauchystress_[0].size();

    return true;
  }
  if ((name == "fiber_cauchy_stress_0_1") || (name == "fiber_cauchy_stress_1_1")) {
    if (data.size()!= 1) dserror("size mismatch");
    for(unsigned gp=0; gp<cauchystress_[1].size(); ++gp) {
      data[0] += cauchystress_[1][gp];
    }
    data[0] = data[0]/cauchystress_[1].size();

    return true;
  }
  if ((name == "fiber_cauchy_stress_0_2") || (name == "fiber_cauchy_stress_1_2")) {
    if (data.size()!= 1) dserror("size mismatch");
    for(unsigned gp=0; gp<cauchystress_[2].size(); ++gp) {
      data[0] += cauchystress_[2][gp];
    }
    data[0] = data[0]/cauchystress_[2].size();

    return true;
  }
  if ((name == "fiber_cauchy_stress_0_3") || (name == "fiber_cauchy_stress_1_3")) {
    if (data.size()!= 1) dserror("size mismatch");
    for(unsigned gp=0; gp<cauchystress_[3].size(); ++gp) {
      data[0] += cauchystress_[3][gp];
    }
    data[0] = data[0]/cauchystress_[3].size();

    return true;
  }


if((name == "cur_rho_col_0_0") || (name == "cur_rho_col_1_0")) {
  if (data.size()!= 1) dserror("size mismatch");
  for(unsigned i=0;i<potsumfiber_[0]->cur_rho.size();++i) {
    data[0] += potsumfiber_[0]->cur_rho[i];
  }
  data[0] = data[0]/potsumfiber_[0]->cur_rho.size();


  return true;
}
if((name == "cur_rho_col_0_1") || (name == "cur_rho_col_1_1")) {
  if (data.size()!= 1) dserror("size mismatch");
  for(unsigned i=0;i<potsumfiber_[1]->cur_rho.size();++i) {
    data[0] += potsumfiber_[1]->cur_rho[i];
  }
  data[0] = data[0]/potsumfiber_[1]->cur_rho.size();


  return true;
}
if((name == "cur_rho_col_0_2") || (name == "cur_rho_col_1_2")) {
  if (data.size()!= 1) dserror("size mismatch");
  for(unsigned i=0;i<potsumfiber_[2]->cur_rho.size();++i) {
    data[0] += potsumfiber_[2]->cur_rho[i];
  }
  data[0] = data[0]/potsumfiber_[2]->cur_rho.size();


  return true;
}
if((name == "cur_rho_col_0_3") || (name == "cur_rho_col_1_3")) {
  if (data.size()!= 1) dserror("size mismatch");
  for(unsigned i=0;i<potsumfiber_[3]->cur_rho.size();++i) {
    data[0] += potsumfiber_[3]->cur_rho[i];
  }
  data[0] = data[0]/potsumfiber_[3]->cur_rho.size();


  return true;
}

dserror("The output is only implemented for four different fiber directions!!!");
return false;
}  // VisData()
