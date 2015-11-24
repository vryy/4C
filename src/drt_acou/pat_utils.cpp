/*!----------------------------------------------------------------------
\file pat_utils.cpp

<pre>
Maintainer: Svenja Schoeder
            schoeder@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/staff/svenja-schoeder/
            089 - 289-15271
</pre>
*----------------------------------------------------------------------*/

#include "pat_utils.H"
#include "pat_imagereconstruction.H"

#include "../drt_inv_analysis/invana_utils.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_timestepping/timintmstep.H"


/*----------------------------------------------------------------------*/
ACOU::PATSearchDirection::PATSearchDirection(INPAR::ACOU::OptimizationType opti)
{
  opti_ = opti;
}

/*----------------------------------------------------------------------*/
void ACOU::PATSearchDirection::Setup(Teuchos::RCP<Epetra_Map> map, Teuchos::RCP<Epetra_Map> uniquemap, int numvecs)
{
  if(opti_==INPAR::ACOU::inv_lbfgs)
  {
    const Teuchos::ParameterList& invp = DRT::Problem::Instance()->StatInverseAnalysisParams();

    actsize_ = 0;
    numvecs_ = numvecs;
    ssize_ = invp.get<int>("SIZESTORAGE")*numvecs;

    sstore_ = Teuchos::rcp(new TIMINT::TimIntMStep<Epetra_Vector>(0, 0, map.get(), true));
    ystore_ = Teuchos::rcp(new TIMINT::TimIntMStep<Epetra_Vector>(0, 0, map.get(), true));

    map_ = map;
    uniquemap_ = uniquemap;

    oldparams_ = Teuchos::rcp(new Epetra_MultiVector(*map_,numvecs_,false));
    oldgradient_ = Teuchos::rcp(new Epetra_MultiVector(*map_,numvecs_,false));
  }
  // nothing to do for steepest descent

  return;
}

/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_MultiVector> ACOU::PATSearchDirection::ComputeDirection(Teuchos::RCP<Epetra_MultiVector> gradient, Teuchos::RCP<Epetra_MultiVector> params, int iter)
{
  Teuchos::RCP<Epetra_MultiVector> direction = Teuchos::rcp(new Epetra_MultiVector(gradient->Map(),gradient->NumVectors(),false));
  if(opti_==INPAR::ACOU::inv_lbfgs)
  {
    if(iter==0)
    {
      oldgradient_->Update(1.0,*gradient,0.0);
      oldparams_->Update(1.0,*params,0.0);
      direction->Update(-1.0,*gradient,0.0);
    }
    else
    {
      // store vectors
      if(iter*numvecs_<=ssize_)
      {
        actsize_+=numvecs_;

        sstore_->Resize(-actsize_+1,0,map_.get(),false);
        ystore_->Resize(-actsize_+1,0,map_.get(),false);

        Epetra_MultiVector s(*map_, numvecs_,false);
        s.Update(1.0,*params,-1.0,*oldparams_,0.0);
        for(int i=0; i<numvecs_; ++i)
          sstore_->UpdateSteps(*s(i));

        s.Update(1.0,*gradient,-1.0,*oldgradient_,false);
        for(int i=0; i<numvecs_; ++i)
          ystore_->UpdateSteps(*s(i));
      }

      // compute direction
      direction->Update(1.0,*gradient,0.0);
      std::vector<double> alpha;

      // loop steps
      for (int i=0; i>-actsize_; i-=numvecs_)
      {
        double a=0.0;
        double b=0.0;
        double aa=0.0;
        double bb=0.0;
        int ind=0;
        for (int j=numvecs_; j>0; j--)
        {
          INVANA::MVDotProduct(*(*ystore_)(i-j+1),*(*sstore_)(i-j+1),*uniquemap_,&a);
          INVANA::MVDotProduct(*(*sstore_)(i-j+1),*(*direction)(ind),*uniquemap_,&b);
          ind++;
          aa += a;
          bb += b;
        }

        alpha.push_back(1/aa*bb);
        ind=0;
        for (int j=numvecs_; j>0; j--)
        {
          (*direction)(ind)->Update(-1.0*alpha.back(), *(*ystore_)(i-j+1),1.0 );
          ind++;
        }
      }

      for (int i=-actsize_+1; i<=0; i+=numvecs_)
      {
        double a=0.0;
        double b=0.0;
        double aa=0.0;
        double bb=0.0;
        double beta=0.0;

        for (int j=0; j<numvecs_; j++)
        {
          INVANA::MVDotProduct(*(*ystore_)(i+j),*(*sstore_)(i+j),*uniquemap_,&a);
          INVANA::MVDotProduct(*(*ystore_)(i+j),*(*direction)(j),*uniquemap_,&b);
          aa += a;
          bb += b;
        }
        beta=1/aa*bb;
        double alphac=alpha.back();
        alpha.pop_back();

        for (int j=0; j<numvecs_; j++)
          (*direction)(j)->Update(alphac-beta, *(*sstore_)(i+j),1.0 );
      }

      // minimization not maximization
      direction->Scale(-1.0);
    }
  }
  else
  {
    if(0)
      direction->Update(-1000.0,*gradient,0.0);
    else
    {
      double maxval = 0.0;
      gradient->MaxValue(&maxval);
      double minval = 0.0;
      gradient->MinValue(&minval);
      if(abs(minval)>maxval)
        direction->Update(0.1/minval,*gradient,0.0);
      else
        direction->Update(-0.1/maxval,*gradient,0.0);
    }
  }
  return direction;
}

/*----------------------------------------------------------------------*/
ACOU::PATLineSearch::PATLineSearch(Teuchos::RCP<PatImageReconstruction> imagereconstruction)
{
  itermax_ = imagereconstruction->acouparams_->sublist("PA IMAGE RECONSTRUCTION").get<int>("INV_LS_MAX_RUN");
  alpha_max_ = 15.0;

  c1_ =  1.0e-12;
  c2_ = 0.9;

  imagereconstruction_ =  imagereconstruction;
  myrank_ = imagereconstruction_->myrank_;
}

/*----------------------------------------------------------------------*/
void ACOU::PATLineSearch::Init(double J0, Teuchos::RCP<Epetra_MultiVector> gradient, Teuchos::RCP<Epetra_MultiVector> direction, Teuchos::RCP<Epetra_MultiVector> state, Teuchos::RCP<Epetra_Map> uniquemap)
{
  // initialize objective function value
  J_0_ = J0;
  J_i_ = J0;
  J_im1_ = J0;

  // initialize vectorial quantities
  dir_ = direction;
  step_ = Teuchos::rcp(new Epetra_MultiVector(dir_->Map(),dir_->NumVectors()));
  state_ = Teuchos::rcp(new Epetra_MultiVector(*state));

  // initialize unique map
  uniquemap_ = uniquemap;

  // initialize norm of derivative of line search function
  imagereconstruction_->CalculateGradDirNorm(*dir_,*uniquemap_,&normgradphi_0_);
  normgradphi_i_ = 0.0;

  // set step lengths
  alpha_i_ = 1.0;
  alpha_im1_ = 0.0;
  alpha_x_ = 0.0;

  return;
}

/*----------------------------------------------------------------------*/
bool ACOU::PATLineSearch::Run()
{
  if(!myrank_)
    std::cout<<"*************** RUN LINE SEARCH: J_0 "<<J_0_<<" normgradphi_0 "<<normgradphi_0_<<std::endl;

  for(int i=0; i<itermax_; ++i)
  {
    if(!myrank_)
      std::cout<<"*************** line search iteration "<<i<<" of maximal "<<itermax_<<" line search iterations, alpha_i_ "<<alpha_i_<<", J_0_ "<<J_0_<<", J_i_ "<< J_i_<<", normgradphi_0_ "<<normgradphi_0_<<", normgradphi_i_ "<<normgradphi_i_<<std::endl;

    // update parameters
    step_->Update(1.0,*state_,0.0);
    step_->Update(alpha_i_,*dir_,1.0);
    imagereconstruction_->ReplaceParams(step_);

    // solve forward problem
    if( imagereconstruction_->GetSequenze()==0 || imagereconstruction_->GetSequenze()==-1 )
    {
      imagereconstruction_->SolveStandardScatra();
      imagereconstruction_->SolveStandardAcou();
    }

    // evaluate objective function
    J_i_ = imagereconstruction_->EvalulateObjectiveFunction();

    // check first condition
    if( J_i_ > J_0_ + c1_ * alpha_i_ * normgradphi_0_ || (J_i_ >= J_im1_ && i>0) )
    {
      if(!myrank_)
        std::cout<<"*************** line search condition 1 met, J_i "<<J_i_<<", J_0 "<<J_0_<<std::endl;
      alpha_x_ = Zoom(alpha_im1_,alpha_i_,J_im1_,false);
      break;
    }
    else if(!myrank_)
      std::cout<<"*************** line search condition 1 NOT met, J_i "<<J_i_<<", J_0 "<<J_0_<<std::endl;

    // solve adjoint problem
    if( imagereconstruction_->GetSequenze()==0 || imagereconstruction_->GetSequenze()==-1 )
    {
      imagereconstruction_->SolveAdjointAcou();
      imagereconstruction_->SolveAdjointScatra();
    }

    // calculate gradient
    imagereconstruction_->EvaluateGradient();
    imagereconstruction_->CalculateGradDirNorm(*dir_,*uniquemap_,&normgradphi_i_);

    // check second condition
    if(std::abs(normgradphi_i_) <= -c2_*normgradphi_0_)
    {
      alpha_x_ = alpha_i_;
      if(!myrank_)
        std::cout<<"*************** line search condition 2 met, |\\/phi|_i "<<normgradphi_i_<<", |\\/phi|_0 "<<normgradphi_0_<<std::endl;
      break;
    }
    else if(!myrank_)
      std::cout<<"*************** line search condition 2 NOT met, |\\/phi|_i "<<normgradphi_i_<<", |\\/phi|_0 "<<normgradphi_0_<<std::endl;

    // check third condition
    if(normgradphi_i_>= 0)
    {
      if(!myrank_)
        std::cout<<"*************** line search condition 3 met"<<std::endl;
      alpha_x_ = Zoom(alpha_i_,alpha_im1_,J_im1_,true);
      break;
    }
    else if(!myrank_)
      std::cout<<"*************** line search condition 3 not met"<<std::endl;

    // update alphas
    alpha_im1_ = alpha_i_;
    alpha_i_ = PredictStepLength();
  }

  if(alpha_x_ != 0.0)
  {
    if(!myrank_)
      std::cout<<"*************** line search succeeded"<<std::endl;
    return true;
  }
  else
  {
    if(!myrank_)
      std::cout<<"*************** line search failed"<<std::endl;
    return false;
  }
}

/*----------------------------------------------------------------------*/
double ACOU::PATLineSearch::PredictStepLength()
{
  // step length must become longer when we are here!
  double alpha = 0.0;

  // calculate the coefficients of the quartic polynomial
  // double d = J_0_;
  double c = normgradphi_0_;
  double b = -1.0/alpha_i_/alpha_i_*(alpha_i_*normgradphi_i_+2.0*normgradphi_0_*alpha_i_-3.0*J_i_+3*J_0_);
  double a = 1.0/alpha_i_/alpha_i_/alpha_i_*(J_i_-J_0_-normgradphi_0_*alpha_i_-b*alpha_i_*alpha_i_);

  // calculate the minima of the quartic polynomial
  double radi = b*b/9.0/a/a-c/3.0/a;
  if(radi>0.0)
  {
    double alpha1 = -b/3.0/a+sqrt(radi);
    double alpha2 = -b/3.0/a-sqrt(radi);
    // check if the results suit us
    if(alpha1>alpha_i_&&alpha1<10.0*alpha_i_)
      alpha = alpha1;
    else if(alpha2>alpha_i_&&alpha2<10.0*alpha_i_)
      alpha = alpha2;
    else if(alpha1>10.0*alpha_i_&&alpha2>10.0*alpha_i_)
      alpha = 10.0*alpha_i_;
    else if(alpha1<alpha_i_&&alpha2<alpha_i_)
      alpha = 2.0*alpha_i_;
    else
      alpha = 2.0*alpha_i_;
  }
  else
  {
    // quadratic interpolation
    alpha = -normgradphi_0_/2.0/(J_i_-J_0_-normgradphi_0_*alpha_i_)*alpha_i_;
    if(alpha<alpha_i_)
      alpha = 2.0*alpha_i_;
    else if(alpha>10.0*alpha_i_)
      alpha = 10.0*alpha_i_;
  }
  return alpha;
}

/*----------------------------------------------------------------------*/
double ACOU::PATLineSearch::Zoom(double alpha_lo, double alpha_hi, double J_alpha_lo, bool turn)
{
  double alpha_j = 0.0;
  double J_j = 0.0;
  double normgradphi_j = 0.0;

  // zoom iterations
  for(int j=0; j<itermax_; ++j)
  {
    // update alpha
    if(turn)
      alpha_j = alpha_hi+(alpha_lo-alpha_hi)/3.0;
    else
      alpha_j = (alpha_lo + alpha_hi) / 2.0;

    // output for user
    if(!myrank_)
      std::cout<<"*************** zoom iteration "<<j<<": alpha_lo "<<alpha_lo<<" alpha_hi "<<alpha_hi<<" alpha_j "<<alpha_j<<std::endl;

    // update parameters
    step_->Update(1.0,*state_,0.0);
    step_->Update(alpha_j,*dir_,1.0);
    imagereconstruction_->ReplaceParams(step_);

    // solve forward problem
    imagereconstruction_->SolveStandardScatra();
    imagereconstruction_->SolveStandardAcou();

    // evaluate objective function
    J_j = imagereconstruction_->EvalulateObjectiveFunction();

    // output for user
    if(!myrank_)
      std::cout<<"J_j "<<J_j<<" J_0_ "<<J_0_<<" J_0_+c1_... "<<J_0_ + c1_ * alpha_j * normgradphi_0_<<" J_alpha_lo "<<J_alpha_lo<<std::endl;

    if( J_j > J_0_ + c1_ * alpha_j * normgradphi_0_ || J_j >= J_alpha_lo )
    {
      if(!myrank_)
        std::cout<<"*************** zoom condition 1 met"<<", J_j "<<J_j<<", J_0 "<<J_0_<<std::endl;
      alpha_hi = alpha_j;
    }
    else
    {
      // solve adjoint problem
      if(imagereconstruction_->GetSequenze()==0 || imagereconstruction_->GetSequenze()==-1 )
      {
        imagereconstruction_->SolveAdjointAcou();
        imagereconstruction_->SolveAdjointScatra();
      }

      // calculate gradient
      imagereconstruction_->EvaluateGradient();
      imagereconstruction_->CalculateGradDirNorm(*dir_,*uniquemap_,&normgradphi_j);

      // check second condition
      if(std::abs(normgradphi_j) <= -c2_*normgradphi_0_)
      {
        if(!myrank_)
          std::cout<<"*************** zoom condition 2 met, |\\/phi|_j "<<normgradphi_j<<", |\\/phi|_0 "<<normgradphi_0_<<std::endl;
        return alpha_j;
      }

      // check third condition
      if( normgradphi_j * (alpha_hi - alpha_lo) >= 0 )
      {
        if(!myrank_)
          std::cout<<"*************** zoom condition 3 met"<<std::endl;
        alpha_hi = alpha_lo;
      }

      alpha_lo = alpha_j;
    }
  }

  return 0.0;
}
