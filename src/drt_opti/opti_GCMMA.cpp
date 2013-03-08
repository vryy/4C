/*!------------------------------------------------------------------------------------------------*
\file opti_GCMMA.cpp

\brief Optimization algorithm GCMMA

<pre>
Maintainer: Martin Winklmaier
            winklmaier@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
 *------------------------------------------------------------------------------------------------*/


#include "opti_GCMMA.H"

#include "../drt_inpar/inpar_topopt.H"
#include "../drt_inpar/inpar_parameterlist_utils.H"
#include "../drt_io/io.H"
#include "../drt_io/io_gmsh.H"
#include "../drt_lib/drt_discret.H"
#include "../headers/definitions.h"
#include <Epetra_SerialDenseSolver.h>
#include <stdio.h>
#include <stdlib.h>

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
OPTI::GCMMA::GCMMA(
    Teuchos::RCP<DRT::Discretization> discret,
    const Teuchos::ParameterList& params,
    Teuchos::RCP<Epetra_Vector> x,
    int numConstraints,
    Teuchos::RCP<Epetra_Vector> x_min,
    Teuchos::RCP<Epetra_Vector> x_max,
    Teuchos::RCP<IO::DiscretizationWriter>& output
) :
discret_(discret),
params_(params),
total_iter_(0),
outer_iter_(0),
inner_iter_(0),
max_total_iter_(params.get<int>("MAX_ITER")),
max_inner_iter_(params.get<int>("MAX_INNER_ITER")),
max_outer_iter_(params.get<int>("MAX_GRAD_ITER")),
m_(numConstraints),
n_loc_(x->MyLength()),
n_(x->GlobalLength()),
x_(Teuchos::rcp(new Epetra_Vector(*x))),
x_old_(Teuchos::rcp(new Epetra_Vector(*x))),
x_old2_(Teuchos::rcp(new Epetra_Vector(*x))),
x_mma_(Teuchos::rcp(new Epetra_Vector(*x))),
x_diff_min_(params.get<double>("X_DIFF_MIN")),
obj_(0.0),
obj_deriv_(Teuchos::rcp(new Epetra_Vector(x->Map()))),
obj_appr_(0.0),
constr_(Teuchos::rcp(new Epetra_SerialDenseVector(m_))),
constr_deriv_(Teuchos::rcp(new Epetra_MultiVector(x->Map(),m_))),
constr_appr_(Teuchos::rcp(new Epetra_SerialDenseVector(m_))),
p0_(Teuchos::rcp(new Epetra_Vector(x->Map()))),
q0_(Teuchos::rcp(new Epetra_Vector(x->Map()))),
r0_(0.0),
P_(Teuchos::rcp(new Epetra_MultiVector(x->Map(),m_))),
Q_(Teuchos::rcp(new Epetra_MultiVector(x->Map(),m_))),
b_(Teuchos::rcp(new Epetra_SerialDenseVector(m_))),
rho0_(0.01),
rho_(Teuchos::rcp(new Epetra_SerialDenseVector(m_))),
rho0min_(params.get<double>("RHOMIN")),
rhomin_(Teuchos::rcp(new Epetra_SerialDenseVector(m_))),
y_mma_(Teuchos::rcp(new Epetra_SerialDenseVector(m_))),
z_mma_(0.0),
xsi_(Teuchos::rcp(new Epetra_Vector(x->Map()))),
eta_(Teuchos::rcp(new Epetra_Vector(x->Map()))),
lam_(Teuchos::rcp(new Epetra_SerialDenseVector(m_))),
mu_(Teuchos::rcp(new Epetra_SerialDenseVector(m_))),
zet_(1.0),
a0_(1.0),
a_(Teuchos::rcp(new Epetra_SerialDenseVector(m_))),
c_(Teuchos::rcp(new Epetra_SerialDenseVector(m_))),
d_(Teuchos::rcp(new Epetra_SerialDenseVector(m_))),
tol_sub_(params.get<double>("TOL_SUB")),
tol_kkt_(params.get<double>("TOL_KKT")),
facmin_(params.get<double>("FACMIN")),
s_(Teuchos::rcp(new Epetra_SerialDenseVector(m_))),
upres_(params_.get<int>("UPRES")),
output_(output)
{
  if (m_>100)
    dserror("current implementation inefficient for large number of constraints due to array structure and used solver");

  if (x_min==Teuchos::null)
  {
    x_min_ = Teuchos::rcp(new Epetra_Vector(x_->Map(),true));

    if (discret_->Comm().MyPID()==0)
      printf("WARNING: Initialized lower boundary for optimization variables with zeros\n");
  }
  else
    x_min_ = Teuchos::rcp(new Epetra_Vector(*x_min));

  if (x_max==Teuchos::null)
  {
    x_max_ = Teuchos::rcp(new Epetra_Vector(x_->Map(),false));
    x_max_->PutScalar(1.0);

    if (discret_->Comm().MyPID()==0)
      printf("WARNING: Initialized upper boundary for optimization variables with ones\n");
  }
  else
    x_max_ = Teuchos::rcp(new Epetra_Vector(*x_max));

  // test case modification
  if ((DRT::INPUT::IntegralValue<INPAR::TOPOPT::OptiTestCases>(params_,"TESTCASE") == INPAR::TOPOPT::optitest_snake_one_constr) or
      (DRT::INPUT::IntegralValue<INPAR::TOPOPT::OptiTestCases>(params_,"TESTCASE") == INPAR::TOPOPT::optitest_snake_multiple_constr))
  {
    x_min_->PutScalar(-2.0);
    x_max_->PutScalar(2.0);
  }

  if ((not x_min_->Map().SameAs(x_max_->Map())) or
      (not x_min_->Map().SameAs(x_->Map())))
    dserror("non-matching maps for boundaries of optimization variables");

  x_diff_ = Teuchos::rcp(new Epetra_Vector(*x_max_));
  x_diff_->Update(-1.0,*x_min_,1.0);

  asymp_min_ = Teuchos::rcp(new Epetra_Vector(*x_min_));
  asymp_max_ = Teuchos::rcp(new Epetra_Vector(*x_max_));

  alpha_ = Teuchos::rcp(new Epetra_Vector(*x_min_));
  beta_ = Teuchos::rcp(new Epetra_Vector(*x_max_));

  double* p1 = a_->Values();
  double* p2 = c_->Values();
  double* p3 = d_->Values();
  double* p4 = rho_->Values();
  double* p5 = rhomin_->Values();

  for (int i=0;i<m_;i++)
  {
    *p1 = 0.0;
    *p2 = 1000*a0_;
    *p3 = a0_;
    *p4 = rho0_;
    *p5 = rho0min_;
    p1++;
    p2++;
    p3++;
    p4++;
    p5++;
  }
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> OPTI::GCMMA::Iterate(
    double& objective,
    Teuchos::RCP<Epetra_Vector> objectivegrad,
    Teuchos::RCP<Epetra_SerialDenseVector> constraints,
    Teuchos::RCP<Epetra_MultiVector> constraintsgrad
)
{
  InitIter(objective,objectivegrad,constraints,constraintsgrad);

  InitSubSolve();

  SubSolve();

  Update();

  return x_mma_; // current solution
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void OPTI::GCMMA::InitIter(
    double& objective,
    Teuchos::RCP<Epetra_Vector> objectivegrad,
    Teuchos::RCP<Epetra_SerialDenseVector> constraints,
    Teuchos::RCP<Epetra_MultiVector> constraintsgrad
)
{
  if ((inner_iter_==0) and // new outer_iter
      (objectivegrad == Teuchos::null or constraintsgrad == Teuchos::null))
    dserror("gradients must be given in new outer iteration");

  // initialization
  if (outer_iter_ == 0)
  {
    // reset of optimization variable for test cases
    if ((DRT::INPUT::IntegralValue<INPAR::TOPOPT::OptiTestCases>(params_,"TESTCASE") == INPAR::TOPOPT::optitest_snake_one_constr) or
        (DRT::INPUT::IntegralValue<INPAR::TOPOPT::OptiTestCases>(params_,"TESTCASE") == INPAR::TOPOPT::optitest_snake_multiple_constr))
    {
      if (x_->Map().NumGlobalElements()%3!=0)
        dserror("cannot be tested");

      int l = x_->Map().NumGlobalElements()/3;

      for(int i=0;i<l;i++)
      {
        double alphai = (3*(i+1)-2*l)*M_PI/(6*l);

        if (x_->Map().MyGID(i))     (*x_)[x_->Map().LID(i)] = cos(alphai + M_PI/12);
        if (x_->Map().MyGID(i+l))   (*x_)[x_->Map().LID(i+l)] = sin(alphai + M_PI/12);
        if (x_->Map().MyGID(i+2*l)) (*x_)[x_->Map().LID(i+2*l)] = sin(2*alphai + M_PI/6);
      }
    }

    // initialization of old values
    *x_old_ = *x_;
    *x_old2_ = *x_;
    *x_mma_ = *x_;
  }


  // reset of objective function, constraints and their derivatives for test cases
  if ((DRT::INPUT::IntegralValue<INPAR::TOPOPT::OptiTestCases>(params_,"TESTCASE") == INPAR::TOPOPT::optitest_snake_one_constr) or
      (DRT::INPUT::IntegralValue<INPAR::TOPOPT::OptiTestCases>(params_,"TESTCASE") == INPAR::TOPOPT::optitest_snake_multiple_constr))
  {
    double locobj = 0.0;
    Epetra_SerialDenseVector locconstr(m_);

    int l = x_mma_->Map().NumGlobalElements()/3;

    double delta = 0.1;
    double alpha = 0.0;
    double locg = 0.0;double locx1 = 0.0;double locx2 = 0.0; double locx3 = 0.0;
    for (int i=0;i<l;i++)
    {
      alpha = M_PI*(3*(i+1)-2*l)/(6*l);

      if (x_->Map().MyGID(i))
      {
        int lid = x_mma_->Map().LID(i);
        double x = (*x_mma_)[lid];
        locobj += cos(alpha)*x;
        locconstr(0) += x*x;
        if (inner_iter_==0) (*objectivegrad)[lid] = cos(alpha);
        if (i==0) {locg+=x*x;locx1=x;}
      }

      if (x_->Map().MyGID(i+l))
      {
        int lid = x_mma_->Map().LID(i+l);
        double x = (*x_mma_)[lid];
        locobj += sin(alpha)*x;
        locconstr(0) += x*x;
        if (inner_iter_==0) (*objectivegrad)[lid] = sin(alpha);
        if (i==0) {locg+=x*x;locx2=x;}
      }

      if (x_->Map().MyGID(i+2*l))
      {
        int lid = x_mma_->Map().LID(i+2*l);
        double x = (*x_mma_)[lid];
        locobj += -0.1*x;
        locconstr(0) += x*x;
        if (inner_iter_==0) (*objectivegrad)[lid] = -0.1;
        if (i==0) locx3=x;
      }
    }
    discret_->Comm().SumAll(&locobj,&objective,1);
    discret_->Comm().SumAll(locconstr.Values(),constraints->Values(),m_);


    double g = 0.0;
    discret_->Comm().SumAll(&locg,&g,1);
    g -= 1;
    g /= delta;

    double h = 0.0;double x1 = 0.0;double x2 = 0.0; double x3 = 0.0;
    discret_->Comm().SumAll(&locx1,&x1,1);
    discret_->Comm().SumAll(&locx2,&x2,1);
    discret_->Comm().SumAll(&locx3,&x3,1);
    h = (x3-2*x1*x2)/delta;

    (*constraints)(0) -= l + 1.0e-5;

    if (DRT::INPUT::IntegralValue<INPAR::TOPOPT::OptiTestCases>(params_,"TESTCASE") == INPAR::TOPOPT::optitest_snake_multiple_constr)
    {
      (*constraints)(1) = +g + pow(g,7) - 2 - 1.0e-5;
      (*constraints)(2) = -g - pow(g,7) - 2 - 1.0e-5;
      (*constraints)(3) = +h + pow(h,7) - 2 - 1.0e-5;
      (*constraints)(4) = -h - pow(h,7) - 2 - 1.0e-5;
    }

    if (inner_iter_==0)
    {
      (*constraintsgrad)(0)->Update(2.0,*x_mma_,0.0);

      if (DRT::INPUT::IntegralValue<INPAR::TOPOPT::OptiTestCases>(params_,"TESTCASE") == INPAR::TOPOPT::optitest_snake_multiple_constr)
      {
        if (x_->Map().MyGID(0))
        {
          int lid = x_mma_->Map().LID(0);
          (*(*constraintsgrad)(1))[lid] = +2*x1/delta * (1+7*pow(g,6));
          (*(*constraintsgrad)(2))[lid] = -2*x1/delta * (1+7*pow(g,6));
          (*(*constraintsgrad)(3))[lid] = -2*x2/delta * (1+7*pow(h,6));
          (*(*constraintsgrad)(4))[lid] = +2*x2/delta * (1+7*pow(h,6));
        }

        if (x_->Map().MyGID(l))
        {
          int lid = x_mma_->Map().LID(l);
          (*(*constraintsgrad)(1))[lid] = +2*x2/delta * (1+7*pow(g,6));
          (*(*constraintsgrad)(2))[lid] = -2*x2/delta * (1+7*pow(g,6));
          (*(*constraintsgrad)(3))[lid] = -2*x1/delta * (1+7*pow(h,6));
          (*(*constraintsgrad)(4))[lid] = +2*x1/delta * (1+7*pow(h,6));
        }

        if (x_->Map().MyGID(2*l))
        {
          int lid = x_mma_->Map().LID(2*l);
          (*(*constraintsgrad)(3))[lid] = +1.0/delta * (1+7*pow(h,6));
          (*(*constraintsgrad)(4))[lid] = -1.0/delta * (1+7*pow(h,6));
        }
      }
    }
  }


  // update counters
  inner_iter_++;
  total_iter_++;

  // new outer iter -> values and derivatives, update asymptotes and re-initialize rho
  if (inner_iter_==1)
  {
    // update counters
    outer_iter_++;

    // reset constraints and objective values
    *constr_ = *constraints;
    *constr_deriv_ = *constraintsgrad;
    obj_ = objective;
    *obj_deriv_ = *objectivegrad;
//    cout << "new obj is " << obj_ << endl;
//    cout << "new obj deriv is " << *obj_deriv_ << endl;
//    cout << "new constr are " << *constr_ << endl;
//    cout << "new constr deriv are " << *constr_deriv_ << endl;

    // reset optimization variables
    *x_old2_ = *x_old_;
    *x_old_ = *x_;
    *x_ = *x_mma_;

    Asymptotes();

    InitRho();
  }
  else // new inner iter -> new values, update rho
  {
//cout << "new obj is " << objective << endl;
//cout << "new constr are " << *constraints << endl;
    UpdateRho(objective,constraints);
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void OPTI::GCMMA::Asymptotes()
{
  // simplified asymptotes in first two iterations
  if (outer_iter_<3)
  {
    *asymp_min_ = *x_;
    *asymp_max_ = *x_;

    asymp_min_->Update(-0.5,*x_diff_,1.0);
    asymp_max_->Update(+0.5,*x_diff_,1.0);
  }
  else
  {
    /*
     * componentwise:
     *
     * asy_min = x - fac*(x_old-asy_min)
     * asy_min = std::max(asy_min, x-10*x_diff)
     * asy_min = min(asy_min, x-0.01*x_diff)
     *
     * asy_max = x + fac*(asy_max-x_old)
     * asy_max = min(asy_max, x+10*x_diff)
     * asy_max = std::max(asy_max, x+0.01*x_diff)
     *
     * fac = 1.0 / 0.7 / 1.2 depending on sign of (x-x_old)(x_old-x_old2)
     */
    double* xval = x_->Values();
    double* xold = x_old_->Values();
    double* xold2 = x_old2_->Values();
    double* asy_min = asymp_min_->Values();
    double* asy_max = asymp_max_->Values();
    double* xdiff = x_diff_->Values();

    for (int i=0;i<n_loc_;i++)
    {
      double val = (*xval-*xold)*(*xold-*xold2);

      double fac = 1.0;
      if (val<0)
        fac = 0.7;
      else if (val>0)
        fac = 1.2;

      *asy_min = *xval - fac*(*xold-*asy_min);
      *asy_max = *xval + fac*(*asy_max-*xold);

      *asy_min = std::max(*asy_min,*xval-10**xdiff);
      *asy_min = std::min(*asy_min,*xval-0.01**xdiff);

      *asy_max = std::min(*asy_max,*xval+10**xdiff);
      *asy_max = std::max(*asy_max,*xval+0.01**xdiff);

      xval++;
      xold++;
      xold2++;
      asy_min++;
      asy_max++;
      xdiff++;
    }
  }
//  cout << "low is " << *asymp_min_ << endl;
//  cout << "upp is " << *asymp_max_ << endl;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void OPTI::GCMMA::InitRho()
{
  /*
   * in matrix form:
   *
   * rho0 = std::max(rho0min,(abs(dJ/dx)^T * xdiff)/(10*size(x)))
   * J: Objective function
   *
   * rho = std::max(rhomin,(abs(dF/dx)^T * xdiff)/(10*size(x)))
   * F: Constraint function(s)
   *
   */

  // set rho of objective function
  double rho0loc = 0.0;
  double* obj_deriv = obj_deriv_->Values();
  double* xdiff = x_diff_->Values();

  for (int i=0;i<n_loc_;i++)
  {
    rho0loc += abs(*obj_deriv)**xdiff;

    obj_deriv++;
    xdiff++;
  }

  // communicate values
  discret_->Comm().SumAll(&rho0loc,&rho0_,1);

  // apply minimal value
  rho0_ = std::max(rho0min_,rho0_/(10*n_));


  // set rho of constraints
  Epetra_SerialDenseVector rholoc(m_);
  double* rho = rholoc.Values();

  for (int j=0;j<m_;j++)
  {
    *rho = 0.0;

    Epetra_Vector constr_der(View,*constr_deriv_,j);
    double* constr_deriv = constr_der.Values();
    xdiff = x_diff_->Values();

    for (int i=0;i<n_loc_;i++)
    {
      *rho += abs(*constr_deriv)**xdiff;

      constr_deriv++;
      xdiff++;
    }
    rho++;
  }

  // communicate values
  discret_->Comm().SumAll(rholoc.Values(),rho_->Values(),m_);

  // apply minimal value
  rho = rho_->Values();
  double* rhomin = rhomin_->Values();

  for (int j=0;j<m_;j++)
  {
    *rho = std::max(*rhomin,*rho/(10*n_));

    rho++;
    rhomin++;
  }
//  cout << "rho0 is " << rho0_ << endl;
//  cout << "rho is " << *rho_ << endl;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void OPTI::GCMMA::UpdateRho(
    double& objective,
    Teuchos::RCP<Epetra_SerialDenseVector> constraints
)
{
  /*
   * fac = (x_mma-x)^2*(asy_max-asy_min)/ ((asy_max-x_mma)*(x_mma-asy_min)*max(x_diff_min,xdiff))
   *
   * rho0 = std::min(1.1*(rho0+(J-J_app)/fac, 10*rho0)
   * rho = std::min(1.1*(rho+(F-F_app)/fac, 10*rho)
   *
   * J/F_app Approximation of J/F
   */

  double fac = 0.0;
  double facloc = 0.0;

  double* x = x_->Values();
  double* xmma = x_mma_->Values();
  double* xdiff = x_diff_->Values();
  double* asy_min = asymp_min_->Values();
  double* asy_max = asymp_max_->Values();

  for (int i=0;i<n_loc_;i++)
  {
    facloc += (*xmma-*x)/(*asy_max-*xmma) * (*xmma-*x)/(*xmma-*asy_min) * (*asy_max-*asy_min)/std::max(x_diff_min_,*xdiff);

    x++;
    xmma++;
    asy_max++;
    asy_min++;
    xdiff++;
  }

  // communicate values
  discret_->Comm().SumAll(&facloc,&fac,1);

  fac = std::max(fac,facmin_);

  if (objective > obj_appr_ + 0.5*tol_sub_)
  {
    rho0_ = std::min(1.1*(rho0_ + (objective-obj_appr_)/fac),10*rho0_);
  }

  double* constr = constraints->Values();
  double* constr_appr = constr_appr_->Values();
  double* rho = rho_->Values();

  for (int i=0;i<m_;i++)
  {
    if (*constr > *constr_appr + 0.5*tol_sub_)
    {
      *rho = std::min(1.1*((*rho) + (*constr-*constr_appr)/fac),10*(*rho));
    }

    constr++;
    constr_appr++;
    rho++;
  }
//  cout << "rho0 is " << rho0_ << endl;
//  cout << "rho is " << *rho_ << endl;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool OPTI::GCMMA::Converged(
    double& objective,
    Teuchos::RCP<Epetra_Vector> objectivegrad,
    Teuchos::RCP<Epetra_SerialDenseVector> constraints,
    Teuchos::RCP<Epetra_MultiVector> constraintsgrad
)
{
  // reset of objective function, constraints and their derivatives for test cases
  if ((DRT::INPUT::IntegralValue<INPAR::TOPOPT::OptiTestCases>(params_,"TESTCASE") == INPAR::TOPOPT::optitest_snake_one_constr) or
      (DRT::INPUT::IntegralValue<INPAR::TOPOPT::OptiTestCases>(params_,"TESTCASE") == INPAR::TOPOPT::optitest_snake_multiple_constr))
  {
    double locobj = 0.0;
    Epetra_SerialDenseVector locconstr(m_);

    int l = x_mma_->Map().NumGlobalElements()/3;

    double delta = 0.1;
    double alpha = 0.0;
    double locg = 0.0;double locx1 = 0.0;double locx2 = 0.0; double locx3 = 0.0;
    for (int i=0;i<l;i++)
    {
      alpha = M_PI*(3*(i+1)-2*l)/(6*l);

      if (x_->Map().MyGID(i))
      {
        int lid = x_mma_->Map().LID(i);
        double x = (*x_mma_)[lid];
        locobj += cos(alpha)*x;
        locconstr(0) += x*x;
        if (i==0) {locg+=x*x;locx1=x;}
      }

      if (x_->Map().MyGID(i+l))
      {
        int lid = x_mma_->Map().LID(i+l);
        double x = (*x_mma_)[lid];
        locobj += sin(alpha)*x;
        locconstr(0) += x*x;
        if (i==0) {locg+=x*x;locx2=x;}
      }

      if (x_->Map().MyGID(i+2*l))
      {
        int lid = x_mma_->Map().LID(i+2*l);
        double x = (*x_mma_)[lid];
        locobj += -0.1*x;
        locconstr(0) += x*x;
        if (i==0) locx3=x;
      }
    }
    discret_->Comm().SumAll(&locobj,&objective,1);
    discret_->Comm().SumAll(locconstr.Values(),constraints->Values(),m_);

    double g = 0.0;
    discret_->Comm().SumAll(&locg,&g,1);
    g -= 1;
    g /= delta;

    double h = 0.0;double x1 = 0.0;double x2 = 0.0; double x3 = 0.0;
    discret_->Comm().SumAll(&locx1,&x1,1);
    discret_->Comm().SumAll(&locx2,&x2,1);
    discret_->Comm().SumAll(&locx3,&x3,1);
    h = (x3-2*x1*x2)/delta;

    (*constraints)(0) -= l + 1.0e-5;

    if (DRT::INPUT::IntegralValue<INPAR::TOPOPT::OptiTestCases>(params_,"TESTCASE") == INPAR::TOPOPT::optitest_snake_multiple_constr)
    {
      (*constraints)(1) = +g + pow(g,7) - 2 - 1.0e-5;
      (*constraints)(2) = -g - pow(g,7) - 2 - 1.0e-5;
      (*constraints)(3) = +h + pow(h,7) - 2 - 1.0e-5;
      (*constraints)(4) = -h - pow(h,7) - 2 - 1.0e-5;
    }


    int dummy=0;
    if ((InnerConvergence(objective,constraints,dummy)==true) and (outer_iter_>0))
    {
      for (int i=0;i<l;i++)
      {
        alpha = M_PI*(3*(i+1)-2*l)/(6*l);

        if (x_->Map().MyGID(i))
        {
          int lid = x_mma_->Map().LID(i);
          (*objectivegrad)[lid] = cos(alpha);
        }

        if (x_->Map().MyGID(i+l))
        {
          int lid = x_mma_->Map().LID(i+l);
          (*objectivegrad)[lid] = sin(alpha);
        }

        if (x_->Map().MyGID(i+2*l))
        {
          int lid = x_mma_->Map().LID(i+2*l);
          (*objectivegrad)[lid] = -0.1;
        }

        (*constraintsgrad)(0)->Update(2.0,*x_mma_,0.0);
      }

      if (DRT::INPUT::IntegralValue<INPAR::TOPOPT::OptiTestCases>(params_,"TESTCASE") == INPAR::TOPOPT::optitest_snake_multiple_constr)
      {
        if (x_->Map().MyGID(0))
        {
          int lid = x_mma_->Map().LID(0);
          (*(*constraintsgrad)(1))[lid] = +2*x1/delta * (1+7*pow(g,6));
          (*(*constraintsgrad)(2))[lid] = -2*x1/delta * (1+7*pow(g,6));
          (*(*constraintsgrad)(3))[lid] = -2*x2/delta * (1+7*pow(h,6));
          (*(*constraintsgrad)(4))[lid] = +2*x2/delta * (1+7*pow(h,6));
        }

        if (x_->Map().MyGID(l))
        {
          int lid = x_mma_->Map().LID(l);
          (*(*constraintsgrad)(1))[lid] = +2*x2/delta * (1+7*pow(g,6));
          (*(*constraintsgrad)(2))[lid] = -2*x2/delta * (1+7*pow(g,6));
          (*(*constraintsgrad)(3))[lid] = -2*x1/delta * (1+7*pow(h,6));
          (*(*constraintsgrad)(4))[lid] = +2*x1/delta * (1+7*pow(h,6));
        }

        if (x_->Map().MyGID(2*l))
        {
          int lid = x_mma_->Map().LID(2*l);
          (*(*constraintsgrad)(3))[lid] = +1.0/delta * (1+7*pow(h,6));
          (*(*constraintsgrad)(4))[lid] = -1.0/delta * (1+7*pow(h,6));
        }
      }
    }
  }


  bool converged = false;
  int num = -1; // number of failed conditions in inner iteration


  if (discret_->Comm().MyPID()==0)
  {
    printf("Checking convergence of optimization algorithm GCMMA\n");
    printf("+-----------------------+-----------+-----------+\n");
    printf("|     Condition         |   Value   | Max-value |\n");
  }

  if ((InnerConvergence(objective,constraints,num)==true) and (outer_iter_>0)) // new outer iteration
  {
    inner_iter_ = 0;

    // check if outer iteration converged
    if (KKTCond(objective,objectivegrad,constraints,constraintsgrad))
      converged = true;

    Epetra_Vector inc(*x_mma_);
    inc.Update(-1.0,*x_,1.0);

    double inc2norm = 0.0;
    inc.Norm2(&inc2norm);
    double incinfnorm = 0.0;
    inc.NormInf(&incinfnorm);

    if (discret_->Comm().MyPID()==0)
    {
      printf("| Increment [L2-norm]   |%10.3E |%10.3E |\n",inc2norm,tol_kkt_);
      printf("| Increment [LInf-norm] |%10.3E |%10.3E |\n",incinfnorm,tol_kkt_);
    }
  }
  else // new inner iteration -> no global convergence
    converged = false;

  if (discret_->Comm().MyPID()==0)
  {
    printf("+-----------------------+-----------+-----------+\n");
    printf("| Total iteration       |      %4d |      %4d |\n",total_iter_+1,max_total_iter_);

    if (inner_iter_==0)
      printf("| Outer iteration       |      %4d |      %4d |\n",outer_iter_+1,max_outer_iter_);
    else
      printf("| Outer iteration       |      %4d |      %4d |\n",outer_iter_,max_outer_iter_);

    printf("| Inner iteration       |      %4d |      %4d |\n",inner_iter_,max_inner_iter_);

    if (num>0)
      printf("| Failing conditions    |      %4d |      %4d |\n",num,m_+1);

    printf("+-----------------------+-----------+-----------+\n");
  }

  // stop if total iteration counter or outer iteration counter reaches their
  // respective maximum number of iterations
  if ((total_iter_==max_total_iter_) or (outer_iter_ == max_outer_iter_))
  {
    if (discret_->Comm().MyPID()==0)
      printf("WARNING: GCMMA optimization algorithm did not converge\n");

    converged = true;
  }



  return converged;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void OPTI::GCMMA::FinishIteration(
    double& objective,
    Teuchos::RCP<Epetra_SerialDenseVector> constraints,
    bool& doGradient
)
{
  // reset of objective function and constraints for test cases
  if ((DRT::INPUT::IntegralValue<INPAR::TOPOPT::OptiTestCases>(params_,"TESTCASE") == INPAR::TOPOPT::optitest_snake_one_constr) or
      (DRT::INPUT::IntegralValue<INPAR::TOPOPT::OptiTestCases>(params_,"TESTCASE") == INPAR::TOPOPT::optitest_snake_multiple_constr))
  {
    double locobj = 0.0;
    Epetra_SerialDenseVector locconstr(m_);

    int l = x_mma_->Map().NumGlobalElements()/3;

    double delta = 0.1;
    double alpha = 0.0;
    double locg = 0.0;double locx1 = 0.0;double locx2 = 0.0; double locx3 = 0.0;
    for (int i=0;i<l;i++)
    {
      alpha = M_PI*(3*(i+1)-2*l)/(6*l);

      if (x_->Map().MyGID(i))
      {
        int lid = x_mma_->Map().LID(i);
        double x = (*x_mma_)[lid];
        locobj += cos(alpha)*x;
        locconstr(0) += x*x;
        if (i==0) {locg+=x*x;locx1=x;}
      }

      if (x_->Map().MyGID(i+l))
      {
        int lid = x_mma_->Map().LID(i+l);
        double x = (*x_mma_)[lid];
        locobj += sin(alpha)*x;
        locconstr(0) += x*x;
        if (i==0) {locg+=x*x;locx2=x;}
      }

      if (x_->Map().MyGID(i+2*l))
      {
        int lid = x_mma_->Map().LID(i+2*l);
        double x = (*x_mma_)[lid];
        locobj += -0.1*x;
        locconstr(0) += x*x;
        if (i==0) locx3=x;
      }
    }
    discret_->Comm().SumAll(&locobj,&objective,1);
    discret_->Comm().SumAll(locconstr.Values(),constraints->Values(),m_);

    double g = 0.0;
    discret_->Comm().SumAll(&locg,&g,1);
    g -= 1;
    g /= delta;

    double h = 0.0;double x1 = 0.0;double x2 = 0.0; double x3 = 0.0;
    discret_->Comm().SumAll(&locx1,&x1,1);
    discret_->Comm().SumAll(&locx2,&x2,1);
    discret_->Comm().SumAll(&locx3,&x3,1);
    h = (x3-2*x1*x2)/delta;

    (*constraints)(0) -= l + 1.0e-5;

    if (DRT::INPUT::IntegralValue<INPAR::TOPOPT::OptiTestCases>(params_,"TESTCASE") == INPAR::TOPOPT::optitest_snake_multiple_constr)
    {
      (*constraints)(1) = +g + pow(g,7) - 2 - 1.0e-5;
      (*constraints)(2) = -g - pow(g,7) - 2 - 1.0e-5;
      (*constraints)(3) = +h + pow(h,7) - 2 - 1.0e-5;
      (*constraints)(4) = -h - pow(h,7) - 2 - 1.0e-5;
    }
  }


  // gradient required for new outer iterations
  if (outer_iter_==0)
    doGradient = true;

  int dummy=0;
  if (InnerConvergence(objective,constraints,dummy))
    doGradient = true;
  else
    doGradient = false;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool OPTI::GCMMA::KKTCond(
    double& objective,
    Teuchos::RCP<Epetra_Vector> objectivegrad,
    Teuchos::RCP<Epetra_SerialDenseVector> constraints,
    Teuchos::RCP<Epetra_MultiVector> constraintsgrad
)
{
  /*
   * Global residual of primal dual interior point algorithm consists of:
   *
   * resx = eta - xsi + dF/dx^T * lam
   * resy = c + d.*y - mu - lam
   * resz = a0 - zet - a^T * lam
   * reslam = F - z_mma_*a - y + s
   * resmu = mu.*y
   * res = lam.*s
   * resxsi = xsi.*(x-x_min)
   * reseta = eta.*(x_max-x)
   * reszet = z*z_mma
   *
   */
  Epetra_Vector resX(*objectivegrad);
  resX.Update(-1.0,*xsi_,1.0);
  resX.Update(+1.0,*eta_,1.0);

  Epetra_SerialDenseVector resy(m_);
  double* resyptr = resy.Values();

  double resz = a0_ - zet_;

  Epetra_SerialDenseVector reslam(m_);
  double* reslamptr = reslam.Values();

  Epetra_SerialDenseVector resmu(m_);
  double* resmuptr = resmu.Values();

  Epetra_SerialDenseVector res(m_);
  double* resptr = res.Values();

  double* lam = lam_->Values();
  double* c = c_->Values();
  double* d = d_->Values();
  double* y = y_mma_->Values();
  double* mu = mu_->Values();
  double* a = a_->Values();
  double* constr = constraints->Values();
  double* s = s_->Values();

  for (int i=0;i<m_;i++)
  {
    Epetra_Vector constr_deriv(View,*constraintsgrad,i);
    resX.Update(*lam,constr_deriv,1.0);

    *resyptr = *c + *d**y-*mu-*lam;
    resz -= *a**lam;
    *reslamptr = *constr - z_mma_**a - *y + *s;
    *resmuptr = *mu**y;
    *resptr = *lam**s;

    resyptr++;
    reslamptr++;
    resmuptr++;
    resptr++;
    c++;
    d++;
    y++;
    mu++;
    lam++;
    a++;
    constr++;
    s++;
  }

  Epetra_Vector resXsi(x_->Map());
  double* resXsiptr = resXsi.Values();

  Epetra_Vector resEta(x_->Map());
  double* resEtaptr = resEta.Values();

  double* xsi = xsi_->Values();
  double* eta = eta_->Values();
  double* x = x_mma_->Values();
  double* xmin = x_min_->Values();
  double* xmax = x_max_->Values();

  for (int i=0;i<n_loc_;i++)
  {
    *resXsiptr = *xsi*(*x-*xmin);
    *resEtaptr = *eta*(*xmax-*x);

    resXsiptr++;
    resEtaptr++;
    xsi++;
    eta++;
    x++;
    xmax++;
    xmin++;
  }

  double reszet = zet_*z_mma_;

  double resnorm = 0.0;
  resnorm = Res2Norm(&resX,&resXsi,&resEta,&resy,&resmu,&reslam,&res,&resz,&reszet);

  double resinf = 0.0;
  resinf = ResInfNorm(&resX,&resXsi,&resEta,&resy,&resmu,&reslam,&res,&resz,&reszet);

  if (discret_->Comm().MyPID()==0)
  {
    printf("+-----------------------+-----------+-----------+\n");
    printf("| Res/KKT [L2-norm]     |%10.3E |%10.3E |\n",resnorm,tol_kkt_);
    printf("| Res/KKT [LInf-norm]   |%10.3E |%10.3E |\n",resinf,tol_kkt_);
  }

  if (resnorm<tol_kkt_)
    return true;

  return false;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool OPTI::GCMMA::InnerConvergence(
    double& objective,
    Teuchos::RCP<Epetra_SerialDenseVector> constraints,
    int& numNotFinished
)
{
  // initially new outer iter
  if (outer_iter_==0)
    return true;

  if ((inner_iter_==max_inner_iter_) and (discret_->Comm().MyPID()==0))
  {
    printf("WARNING: inner GCMMA optimization loop did not converge\n");

    return true;
  }

  // if all approximations (constraints and objective) are smaller than the
  // values itself (modulo a small tolerance), then inner iteration is finished
  bool finished = true;
  numNotFinished = 0;

  if (obj_appr_ + tol_sub_ < objective)
  {cout << "objective not finished" << endl;
    finished = false;
    numNotFinished++;
  }

  double* constr_appr = constr_appr_->Values();
  double* constr = constraints->Values();

  for (int i=0;i<m_;i++)
  {
    if (*constr_appr + tol_sub_ < *constr)
    {cout << "constraint " << i << " not finished" << endl;
      finished = false;
      numNotFinished++;
    }

    constr_appr++;
    constr++;
  }

  return finished;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void OPTI::GCMMA::InitSubSolve()
{
  /*
   * alpha = (1-albeta)*asymp_min + albeta*x
   * beta = (1-albeta)*asymp_max + albeta*x
   */
  double albeta = 0.1;

  alpha_->Update(1.0 - albeta,*asymp_min_,0.0);
  alpha_->Update(albeta,*x_,1.0);

  beta_->Update(1.0 - albeta,*asymp_max_,0.0);
  beta_->Update(albeta,*x_,1.0);

  double* alpha = alpha_->Values();
  double* beta = beta_->Values();
  double* xmin = x_min_->Values();
  double* xmax = x_max_->Values();

  for (int i=0;i<n_loc_;i++)
  {
    *alpha = std::max(*alpha,*xmin);
    *beta = std::min(*beta,*xmax);

    alpha++;
    beta++;
    xmin++;
    xmax++;
  }


  /*
   * ux1 = asymp_max - x
   * ux2 = (asymp_max - x).^2
   * xl1 = x - asymp_min
   * xl2 = (x - asymp_min).^2
   */
  double* x = x_->Values();
  double* asymp_max = asymp_max_->Values();
  double* asymp_min = asymp_min_->Values();

  Epetra_Vector ux1(x_->Map(),false);
  Epetra_Vector ux2(x_->Map(),false);
  Epetra_Vector xl1(x_->Map(),false);
  Epetra_Vector xl2(x_->Map(),false);

  double* ux1ptr = ux1.Values();
  double* ux2ptr = ux2.Values();
  double* xl1ptr = xl1.Values();
  double* xl2ptr = xl2.Values();

  for (int i=0;i<n_loc_;i++)
  {
    *ux1ptr = *asymp_max-*x;
    *ux2ptr = *ux1ptr**ux1ptr;
    *xl1ptr = *x-*asymp_min;
    *xl2ptr = *xl1ptr**xl1ptr;

    ux1ptr++;
    ux2ptr++;
    xl1ptr++;
    xl2ptr++;
    x++;
    asymp_max++;
    asymp_min++;
  }


  /* compute
   * p0 = (max(dJ/dx,0) + fac*p0q0 + rho0/max(x_diff,x_diff_min)).*ux2
   * q0 = (max(-dJ/dx,0) + fac*p0q0 + rho0/max(x_diff,x_diff_min)).*xl2
   * r0 = J - p0./ux1 - q0./xl1
   *
   * with
   * p0q0 = (max(dJ/dx,0) + (max(-dJ/dx,0)
   * fac << 1
   */
  double* p0 = p0_->Values();
  double* q0 = q0_->Values();
  double* obj_deriv = obj_deriv_->Values();

  double p0q0; // tmp variable
  double xdiffinv;

  ux1ptr = ux1.Values();
  ux2ptr = ux2.Values();
  xl1ptr = xl1.Values();
  xl2ptr = xl2.Values();

  double* xdiff = x_diff_->Values();

  double fac = 0.001;

  double r0loc = 0.0;

  for (int i=0;i<n_loc_;i++)
  {
    *p0 = std::max(*obj_deriv,0.0);
    *q0 = std::max(-*obj_deriv,0.0);

    xdiffinv = 1.0/std::max(*xdiff,x_diff_min_);

    p0q0 = *p0 + *q0;
    *p0 = (*p0 + fac*p0q0 + rho0_*xdiffinv)**ux2ptr;
    *q0 = (*q0 + fac*p0q0 + rho0_*xdiffinv)**xl2ptr;

    r0loc -= *p0/(*ux1ptr) + *q0/(*xl1ptr);

    p0++;
    q0++;
    obj_deriv++;
    xdiff++;
    ux1ptr++;
    ux2ptr++;
    xl1ptr++;
    xl2ptr++;
  }
  discret_->Comm().SumAll(&r0loc,&r0_,1);

  r0_ += obj_;


  /* compute
   * P = (max(dF/dx,0) + fac*pq + rho/max(x_diff,x_diff_min)).*ux2
   * Q = (max(-dF/dx,0) + fac*pq + rho/max(x_diff,x_diff_min)).*xl2
   * b = F - P^T/ux1 - Q^T/xl1
   *
   * with
   * pq = (max(dF/dx,0) + (max(-dF/dx,0)
   * fac << 1
   */
  Epetra_SerialDenseVector bloc(m_);
  double* b = bloc.Values();

  double* rho = rho_->Values();

  for (int i=0;i<m_;i++)
  {
    Epetra_Vector p_i(View,*P_,i);
    Epetra_Vector q_i(View,*Q_,i);
    Epetra_Vector constr_deriv(View,*constr_deriv_,i);

    double* p = p_i.Values();
    double* q = q_i.Values();
    double* constr_der = constr_deriv.Values();

    double pq = 0.0; // tmp variable

    xdiff = x_diff_->Values();

    ux1ptr = ux1.Values();
    ux2ptr = ux2.Values();
    xl1ptr = xl1.Values();
    xl2ptr = xl2.Values();

    for (int j=0;j<n_loc_;j++)
    {
      *p = std::max(*constr_der,0.0);
      *q = std::max(-*constr_der,0.0);

      xdiffinv = 1.0/std::max(*xdiff,x_diff_min_);

      pq = *p + *q;
      *p = (*p + fac*pq + *rho*xdiffinv)**ux2ptr;
      *q = (*q + fac*pq + *rho*xdiffinv)**xl2ptr;

      p++;
      q++;
      constr_der++;
      xdiff++;
      ux1ptr++;
      ux2ptr++;
      xl1ptr++;
      xl2ptr++;
    }

    rho++;
  }
  for (int i=0;i<m_;i++)
  {
    Epetra_Vector p_i(View,*P_,i);
    Epetra_Vector q_i(View,*Q_,i);
    double* p = p_i.Values();
    double* q = q_i.Values();
    ux1ptr = ux1.Values();
    xl1ptr = xl1.Values();
    for (int j=0;j<n_loc_;j++)
    {
      *b -= *p/(*ux1ptr) + *q/(*xl1ptr);
      p++;
      q++;
      ux1ptr++;
      xl1ptr++;
    }
    b++;
  }
  discret_->Comm().SumAll(bloc.Values(),b_->Values(),m_);
  *b_ += *constr_;
  b_->Scale(-1.0);

//  cout << "init subsolve: " << endl;
//  cout << "m is " << m_ << ", n is " << n_ << ", epsimin is " << tol_sub_ << endl;
//  cout << "low asy is " << *asymp_min_ << "upp asy is " << *asymp_max_ << endl;
//  cout << "alpha is " << *alpha_ << "beta is " << *beta_ << endl;
//  cout << "p0 is " << *p0_ << "q0 is " << *q0_ << endl;
//  cout << "p is " << *P_ << "q is " << *Q_ << endl;
//  cout << "a0 is " << a0_ << ", a is " << *a_ << ", b is " << *b_ << ", c is " << *c_ << ", d is " << *d_ << endl;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void OPTI::GCMMA::SubSolve()
{
  /*
   * Initialisation
   *
   * x_mma = 0.5*(alpha+beta)
   * y = ones
   * lam = ones
   * mu = std::max(ones,0.5*c)
   * s = ones
   * z_mma = 1.0
   * xsi = std::max(ones,ones/(x-alpha))
   * eta = std::max(ones,ones/(beta-x))
   */
  x_mma_->Update(0.5,*alpha_,0.0);
  x_mma_->Update(0.5,*beta_,1.0);

  double* y = y_mma_->Values();
  double* lam = lam_->Values();
  double* mu = mu_->Values();
  double* s = s_->Values();
  double* c = c_->Values();

  for (int i=0;i<m_;i++)
  {
    *y = 1.0;
    *lam = 1.0;
    *mu = std::max(1.0, *c/2.0);
    *s = 1.0;

    y++;
    lam++;
    mu++;
    s++;
    c++;
  }


  z_mma_ = 1.0;


  double* xsi = xsi_->Values();
  double* eta = eta_->Values();
  double* x = x_mma_->Values();
  double* alpha = alpha_->Values();
  double* beta = beta_->Values();
  for (int i=0;i<n_loc_;i++)
  {
    *xsi = std::max(1.0, 1.0/(*x-*alpha));
    *eta = std::max(1.0, 1.0/(*beta-*x));

    xsi++;
    eta++;
    x++;
    alpha++;
    beta++;
  }




  int total_iter = 0;

  double tol_sub = 1.0;
  double tol_fac = 0.1; // factor from one tolerance to next (smaller) one
  if (tol_fac>0.99)
    dserror("factor for tolerance adaption shall be significant smaller than one!");

  bool tol_reached = false;
  /*
   * Initialisation finished
   */

  while (tol_reached == false)
  {
    // first step: compute norms of all composed residuals
    double resnorm = 0.0;
    double resinf = 0.0;
    ResApp(resnorm,resinf,tol_sub);

    int inner_iter = 0;
    int max_sub_inner_iter_ = 200;

    // iteration loop of inner algorithm
    while (inner_iter<max_sub_inner_iter_)
    {
      inner_iter++;
      total_iter++;

      /*
       * helper
       *
       * plam = P^T * lam
       * qlam = Q^T * lam
       */
      Epetra_Vector plam(*p0_);
      Epetra_Vector qlam(*q0_);
      double* lam = lam_->Values();
      for (int j=0;j<m_;j++)
      {
        Epetra_Vector P(View,*P_,j);
        Epetra_Vector Q(View,*Q_,j);

        plam.Update(*lam,P,1.0);
        qlam.Update(*lam,Q,1.0);

        lam++;
      }


      /*
       * helper data
       *
       * delx = plam./ux2 - qlam./xl2 - tol_sub/(x-alpha) + tol_sub/(beta-x)
       * diagx = 2*(plam/(ux1^3) + qlamptr/(xl1^3)) + xsi/(x-alpha) + eta/(beta-x)
       */
      Epetra_Vector delx(x_->Map(),false);
      Epetra_Vector diagx(x_->Map(),false);

      Epetra_Vector uxinv1(x_->Map(),false);
      Epetra_Vector xlinv1(x_->Map(),false);

      double ux1 = 0.0;
      double xl1 = 0.0;

      double ux2 = 0.0;
      double xl2 = 0.0;

      double* uxinv1ptr = uxinv1.Values();
      double* xlinv1ptr = xlinv1.Values();

      double* x = x_mma_->Values();
      double* asymp_max = asymp_max_->Values();
      double* asymp_min = asymp_min_->Values();

      double* delxptr = delx.Values();
      double* diagxptr = diagx.Values();

      xsi = xsi_->Values();
      eta = eta_->Values();

      alpha = alpha_->Values();
      beta = beta_->Values();


      double* plamptr = plam.Values();
      double* qlamptr = qlam.Values();

      for (int i=0;i<n_loc_;i++)
      {
        ux1 = *asymp_max-*x;
        xl1 = *x-*asymp_min;

        ux2 = ux1*ux1;
        xl2 = xl1*xl1;

        *uxinv1ptr = 1.0/ux1;
        *xlinv1ptr = 1.0/xl1;

        double dpsidx = *plamptr/ux2 - *qlamptr/xl2;

        *delxptr = dpsidx - tol_sub/(*x-*alpha) + tol_sub/(*beta-*x);
        *diagxptr = 2*(*plamptr/(ux2*ux1) + *qlamptr/(xl2*xl1)) + *xsi/(*x-*alpha) + *eta/(*beta-*x);

        uxinv1ptr++;
        xlinv1ptr++;
        delxptr++;
        diagxptr++;
        asymp_max++;
        asymp_min++;
        x++;
        plamptr++;
        qlamptr++;
        alpha++;
        beta++;
        xsi++;
        eta++;
      }


      /*
       * helper data
       *
       * GG = P./(ux1^2) - Q./(xl1^2)
       */
      Epetra_MultiVector GG(x_->Map(),m_);

      for (int i=0;i<m_;i++)
      {
        Epetra_Vector gg(View,GG,i);

        Epetra_Vector P(View,*P_,i);
        Epetra_Vector Q(View,*Q_,i);

        double* ggptr = gg.Values();
        double* pptr = P.Values();
        double* qptr = Q.Values();
        uxinv1ptr = uxinv1.Values();
        xlinv1ptr = xlinv1.Values();

        for (int j=0;j<n_loc_;j++)
        {
          *ggptr = *pptr*(*uxinv1ptr**uxinv1ptr) - *qptr*(*xlinv1ptr**xlinv1ptr);

          ggptr++;
          pptr++;
          qptr++;
          uxinv1ptr++;
          xlinv1ptr++;
        }
      }


      /*
       * dely = c + d.*y - lam - tol_sub/y;
       * delz = a0 - tol_sub/z_mma - a^T*lam
       * dellam = P^T/ux + Q^T/xl - z_mma*a - y - b + tol_sub/lam
       * diagy = d + mu/y
       * diaglamyi = s./lam + 1.0/diagy;
       */
      Epetra_SerialDenseVector dely(m_);
      double* delyptr = dely.Values();

      double delz = a0_ - tol_sub/z_mma_;

      Epetra_SerialDenseVector dellam(m_);
      double* dellamptr = dellam.Values();

      Epetra_SerialDenseVector diagy(m_);
      double* diagyptr = diagy.Values();

      Epetra_SerialDenseVector diaglamyi(m_);
      double* diaglamyiptr = diaglamyi.Values();

      c = c_->Values();
      double* d = d_->Values();
      y = y_mma_->Values();
      lam = lam_->Values();
      double* a = a_->Values();
      double gvec1 = 0.0;
      double gvec2 = 0.0;
      double* b = b_->Values();
      mu = mu_->Values();
      s = s_->Values();

      for (int i=0;i<m_;i++)
      {
        *delyptr = *c + *d**y - *lam - tol_sub/(*y);

        delz -= *a**lam;

        uxinv1.Dot(*(*P_)(i),&gvec1);
        xlinv1.Dot(*(*Q_)(i),&gvec2);
        *dellamptr = gvec1 + gvec2 - *a*z_mma_ - *y - *b + tol_sub/(*lam);

        *diagyptr = *d + *mu/(*y);

        *diaglamyiptr = *s/(*lam) + 1.0/(*diagyptr);

        delyptr++;
        dellamptr++;
        diagyptr++;
        diaglamyiptr++;
        c++;
        d++;
        y++;
        lam++;
        a++;
        b++;
        mu++;
        s++;
      }


      Epetra_SerialDenseVector dlam(m_);
      double dz = 0.0;
      Epetra_Vector dx(x_->Map());

      if (m_ < x_->GlobalLength())
      {
        /*
         *
         * Build final system matrix and solve (size m+1)
         *
         * rhs = bb = [dellam + dely./diagy - GG^T*delx./diagx
         *             delz]
         * sysmat = AA = [GG^T * diagx(diagmat) * GG + diaglamyi (diagmat) || a
         *                a^T                                              || -zet/z_mma]
         * solut = AA^-1 * bb
         */
        Epetra_SerialDenseVector bb(m_+1);
        double* bbval = bb.A();

        dellamptr = dellam.Values();
        delyptr = dely.Values();
        diagyptr = diagy.Values();

        for (int i=0;i<m_;i++)
        {
          double locval = 0.0;
          double val = 0.0;

          Epetra_Vector gg(View,GG,i);
          double* ggptr = gg.Values();
          delxptr = delx.Values();
          diagxptr = diagx.Values();

          for (int j=0;j<n_loc_;j++)
          {
            locval += *ggptr**delxptr/(*diagxptr);

            ggptr++;
            delxptr++;
            diagxptr++;
          }
          discret_->Comm().SumAll(&locval,&val,1);

          *bbval = *dellamptr + *delyptr/(*diagyptr) - val;

          bbval++;
          dellamptr++;
          delyptr++;
          diagyptr++;
        }
        *bbval = delz;


        Epetra_SerialDenseMatrix AA(m_+1,m_+1);
        double* aa = AA.A();
        a = a_->Values();
        diaglamyiptr = diaglamyi.Values();

        for (int icol=0;icol<m_;icol++)
        {
//          Alam = spdiags(diaglamyi,0,m,m) + GG*spdiags(diagxinv,0,n,n)*GG';
          Epetra_Vector colGG(View,GG,icol);
          double* colgg;

          for (int irow=0;irow<m_;irow++)
          {
            colgg = colGG.Values(); // reset here every loop

            Epetra_Vector rowGG(View,GG,irow);
            double* rowgg = rowGG.Values();

            double locval = 0.0;
            diagxptr = diagx.Values();

            for (int j=0;j<n_loc_;j++)
            {
              locval += *rowgg /(*diagxptr) * *colgg;

              rowgg++;
              diagxptr++;
              colgg++;
            }

            discret_->Comm().SumAll(&locval,aa,1);

            if (irow==icol)
              *aa += *diaglamyiptr;

            aa++;
          }

          *aa = *a;

          a++;
          aa++;
          diaglamyiptr++;
        }

        a = a_->Values();
        for (int i=0;i<m_;i++)
        {
          *aa = *a;

          a++;
          aa++;
        }

        *aa = -zet_/z_mma_;


        Epetra_SerialDenseVector solut(m_+1);

        Epetra_SerialDenseSolver solver;
        solver.SetMatrix(AA);
        solver.SetVectors(solut,bb);
        solver.Solve();


        /*
         *
         * dlam = solut(1:m) // all except for the last entry / first m entries
         * dz = solut(end) // last entry
         * dx = (-GG^T*dlam - delx)./diagx
         *
         */
        double* dlamptr = dlam.Values();
        double* sol = solut.A();

        for (int i=0;i<m_;i++)
        {
          *dlamptr = *sol;

          dlamptr++;
          sol++;
        }


        dz = *sol;


        dlamptr = dlam.Values();
        for (int i=0;i<m_;i++)
        {
          Epetra_Vector gg(View,GG,i);
          dx.Update(-*dlamptr,gg,1.0);

          dlamptr++;
        }

        double* dxptr = dx.Values();
        delxptr = delx.Values();
        diagxptr = diagx.Values();

        for (int i=0;i<n_loc_;i++)
        {
          *dxptr -= *delxptr;
          *dxptr /= *diagxptr;

          dxptr++;
          delxptr++;
          diagxptr++;
        }
      }
      else
      {
        dserror("Not implemented! Efficient data structures and solver required here");
//        see matlab code:
//        diaglamyiinv = eem./diaglamyi;
//        dellamyi = dellam + dely./diagy;
//        Axx = spdiags(diagx,0,n,n) + GG'*spdiags(diaglamyiinv,0,m,m)*GG;
//        azz = zet/z + a'*(a./diaglamyi);
//        axz = -GG'*(a./diaglamyi);
//        bx = delx + GG'*(dellamyi./diaglamyi);
//        bz  = delz - a'*(dellamyi./diaglamyi);
//        AA = [Axx   axz
//              axz'  azz ];
//        bb = [-bx' -bz]';
//        solut = AA\bb;
//        dx  = solut(1:n);
//        dz = solut(n+1);
//        dlam = (GG*dx)./diaglamyi - dz*(a./diaglamyi) + dellamyi./diaglamyi;
      }


      /*
       *
       * dxsi = -xsi + (tol_sub - xsi*dx)/(x-alpha)
       * deta = -eta + (tol_sub + eta*dx)/(beta-x)
       *
       */
      Epetra_Vector dxsi(x_->Map());
      Epetra_Vector deta(x_->Map());

      double* dxsiptr = dxsi.Values();
      double* detaptr = deta.Values();
      xsi = xsi_->Values();
      eta = eta_->Values();
      double* dxptr = dx.Values();
      x = x_mma_->Values();
      alpha = alpha_->Values();
      beta = beta_->Values();

      for (int i=0;i<n_loc_;i++)
      {
        *dxsiptr = -*xsi + (tol_sub - *xsi**dxptr)/(*x-*alpha);
        *detaptr = -*eta + (tol_sub + *eta**dxptr)/(*beta-*x);

        dxsiptr++;
        detaptr++;
        xsi++;
        eta++;
        dxptr++;
        x++;
        alpha++;
        beta++;
      }


      /*
       *
       * dy = - dely./diagy + dlam./diagy
       * dmu = -mu + tol_sub/y - mu.*dy./y
       * ds = -s + tol_sub/lam - s.*dlam./lam
       *
       */
      Epetra_SerialDenseVector dy(m_);
      Epetra_SerialDenseVector dmu(m_);
      Epetra_SerialDenseVector ds(m_);

      double* dyptr = dy.Values();
      double* dmuptr = dmu.Values();
      double* dsptr = ds.Values();
      delyptr = dely.Values();
      diagyptr = diagy.Values();
      double* dlamptr = dlam.Values();
      mu = mu_->Values();
      y = y_mma_->Values();
      s = s_->Values();
      lam = lam_->Values();

      for (int i=0;i<m_;i++)
      {
        *dyptr = -*delyptr/(*diagyptr) + *dlamptr/(*diagyptr);
        *dmuptr = -*mu + tol_sub/(*y) - (*mu**dyptr)/(*y);
        *dsptr = -*s + tol_sub/(*lam) - (*s**dlamptr)/(*lam);

        dyptr++;
        delyptr++;
        diagyptr++;
        dlamptr++;
        dmuptr++;
        mu++;
        y++;
        dsptr++;
        s++;
        lam++;
      }


      double dzet = -zet_ + tol_sub/z_mma_ - zet_*dz/z_mma_;


      double val = std::min(dz/z_mma_, dzet/zet_);

      y = y_mma_->Values();
      dyptr = dy.Values();
      lam = lam_->Values();
      dlamptr = dlam.Values();
      mu = mu_->Values();
      dmuptr = dmu.Values();
      s = s_->Values();
      dsptr = ds.Values();

      for (int i=0;i<m_;i++)
      {
        val = std::min(val, *dyptr/(*y));
        val = std::min(val, *dlamptr/(*lam));
        val = std::min(val, *dmuptr/(*mu));
        val = std::min(val, *dsptr/(*s));

        y++;
        dyptr++;
        lam++;
        dlamptr++;
        mu++;
        dmuptr++;
        s++;
        dsptr++;
      }

      xsi = xsi_->Values();
      dxsiptr = dxsi.Values();
      eta = eta_->Values();
      detaptr = deta.Values();
      dxptr = dx.Values();
      x = x_mma_->Values();
      alpha = alpha_->Values();
      beta = beta_->Values();

      for (int i=0;i<n_loc_;i++)
      {
        val = std::min(val, *dxsiptr/(*xsi));
        val = std::min(val, *detaptr/(*eta));
        val = std::min(val, *dxptr/(*x-*alpha));
        val = std::min(val, *dxptr/(*x-*beta));

        xsi++;
        dxsiptr++;
        eta++;
        detaptr++;
        dxptr++;
        x++;
        alpha++;
        beta++;
      }

      double fac = -1.01;
      if (fac>-1.0) dserror("unsensible factor");

      // min becomes max since fac<0
      val = std::max(1.0,fac*val);
      double steg = 0.0;
      discret_->Comm().MaxAll(&val,&steg,1);
      steg = 1.0/steg;
      /*
       * steg = 1.0/max(1.0
       *                fac*min(dz/z_mma
       *                        dzet/zet
       *                        dy/y
       *                        dlam/lam
       *                        dmu/mu
       *                        ds/s
       *                        dxsi/xsi
       *                        deta/eta
       *                        dx/(x-alpha)
       *                        dx/(-beta)
       *                       )
       *               )
       */



      int it = 0;
      int max_it = 50;
      double resnew = 2*resnorm;


      // save old values
      Epetra_Vector x_mma_old(*x_mma_);
      Epetra_Vector xsi_old(*xsi_);
      Epetra_Vector eta_old(*eta_);

      Epetra_SerialDenseVector y_old(*y_mma_);
      Epetra_SerialDenseVector lam_old(*lam_);
      Epetra_SerialDenseVector mu_old(*mu_);
      Epetra_SerialDenseVector s_old(*s_);
      double z_mma_old = z_mma_;
      double zet_old = zet_;


      while ((resnorm<resnew) and (it<max_it))
      {
        it++;

        x_mma_->Update(steg,dx,1.0,x_mma_old,0.0);
        xsi_->Update(steg,dxsi,1.0,xsi_old,0.0);
        eta_->Update(steg,deta,1.0,eta_old,0.0);

        y = y_mma_->Values();
        double* yptr = y_old.Values();
        dyptr = dy.Values();
        lam = lam_->Values();
        double* lamptr = lam_old.Values();
        dlamptr = dlam.Values();
        mu = mu_->Values();
        double* muptr = mu_old.Values();
        dmuptr = dmu.Values();
        s = s_->Values();
        double* sptr = s_old.Values();
        dsptr = ds.Values();

        for (int i=0;i<m_;i++)
        {
          *y = *yptr + steg**dyptr;
          *lam = *lamptr + steg**dlamptr;
          *mu = *muptr + steg**dmuptr;
          *s = *sptr + steg**dsptr;

          y++;
          yptr++;
          dyptr++;
          lam++;
          lamptr++;
          dlamptr++;
          mu++;
          muptr++;
          dmuptr++;
          s++;
          sptr++;
          dsptr++;
        }

        z_mma_ = z_mma_old + steg*dz;
        zet_ = zet_old + steg*dzet;


        // compute residuals
        ResApp(resnew,resinf,tol_sub);

        steg = steg/2.0;

        if ((it==max_it) and (discret_->Comm().MyPID()==0))
          printf("Reached maximal number of iterations in most inner loop of primal dual interior point optimization algorithm\n");
      }

      resnorm = resnew;

      // it would be sufficient to compute resinf only once here and not in
      // every iteration as done above. but the effort is neglegible
      if (resinf < 0.9*tol_sub)
        break;
    }

    if ((inner_iter==max_sub_inner_iter_) and (discret_->Comm().MyPID()==0))
      printf("Reached maximal number of iterations in inner loop of primal dual interior point optimization algorithm\n");

    if (tol_sub > 1.001*tol_sub_)
      tol_sub *= tol_fac;
    else
      tol_reached = true;
  }
//  cout << "after subsolv:" << endl;
//  cout << "x is " << *x_mma_ << endl;
//  cout << "y is " << *y_mma_ << ", z is " << z_mma_ << ", lam is " << *lam_ << endl;
//  cout << "xsi is " << *xsi_ << endl;
//  cout << "eta is " << *eta_ << endl;
//  cout << "mu is " << *mu_ << ", zet is " << zet_ << ", s is " << *s_ << endl;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void OPTI::GCMMA::ResApp(
    double& resnorm,
    double& resinf,
    const double& tol_sub
)
{
  /*
   * helper
   *
   * plam = P^T*lam
   * qlam = Q^T*lam
   */
  Epetra_Vector plam(*p0_);
  Epetra_Vector qlam(*q0_);

  double* lam = lam_->Values();

  for (int j=0;j<m_;j++)
  {
    Epetra_Vector P(View,*P_,j);
    Epetra_Vector Q(View,*Q_,j);

    plam.Update(*lam,P,1.0);
    qlam.Update(*lam,Q,1.0);

    lam++;
  }


  /*
   * resx = plam/(ux1*ux1) - qlam/(xl1*xl1) - xsi + eta;
   */
  Epetra_Vector resX(x_->Map());
  double* resxptr = resX.Values();

  Epetra_Vector uxinv1(x_->Map(),false);
  Epetra_Vector xlinv1(x_->Map(),false);
  double* uxinv1ptr = uxinv1.Values();
  double* xlinv1ptr = xlinv1.Values();

  double ux1 = 0.0; // one entry of Epetra_Vector, just used locally
  double xl1 = 0.0; // one entry of Epetra_Vector, just used locally
  double* x = x_mma_->Values();
  double* asymp_max = asymp_max_->Values();
  double* asymp_min = asymp_min_->Values();
  double* plamptr = plam.Values();
  double* qlamptr = qlam.Values();
  double* xsi = xsi_->Values();
  double* eta = eta_->Values();

  for (int i=0;i<n_loc_;i++)
  {
    ux1 = *asymp_max-*x;
    xl1 = *x-*asymp_min;

    *uxinv1ptr = 1.0/ux1;
    *xlinv1ptr = 1.0/xl1;

    double dpsidx = *plamptr/(ux1*ux1) - *qlamptr/(xl1*xl1);
    *resxptr = dpsidx - *xsi + *eta;

    uxinv1ptr++;
    xlinv1ptr++;
    resxptr++;
    x++;
    asymp_max++;
    asymp_min++;
    plamptr++;
    qlamptr++;
    xsi++;
    eta++;
  }


  /*
   * resy = c + d.*y - mu - lam
   * resz = a0 - zet - a^T*lam
   * reslam = P^T/ux1 + Q^T/xl - a*z_mma_ - y + s - b
   * resmu = mu.*y - tol_sub
   * res = lam.*s - tol_sub
   */
  Epetra_SerialDenseVector resy(m_);
  double* resyptr = resy.Values();

  double resz = a0_ - zet_;

  Epetra_SerialDenseVector reslam(m_);
  double* reslamptr = reslam.Values();

  Epetra_SerialDenseVector resmu(m_);
  double* resmuptr = resmu.Values();

  Epetra_SerialDenseVector res(m_);
  double* resptr = res.Values();

  double* c = c_->Values();
  double* d = d_->Values();
  double* y = y_mma_->Values();
  double* mu = mu_->Values();
  lam = lam_->Values();
  double* a = a_->Values();
  double* s = s_->Values();
  double* b = b_->Values();
  double gvec1 = 0.0;
  double gvec2 = 0.0;

  for (int j=0;j<m_;j++)
  {
    *resyptr = *c + *d**y - *mu - *lam;

    resz -= *a**lam;

    uxinv1.Dot(*(*P_)(j),&gvec1);
    xlinv1.Dot(*(*Q_)(j),&gvec2);
    *reslamptr = gvec1 + gvec2 - *a*z_mma_ - *y + *s - *b;

    *resmuptr = *mu**y - tol_sub;

    *resptr = *lam**s - tol_sub;

    resyptr++;
    reslamptr++;
    resmuptr++;
    resptr++;
    c++;
    d++;
    y++;
    mu++;
    lam++;
    a++;
    s++;
    b++;
  }


  /*
   * resxsi = xsi.*(x-alpha) - tol_sub
   * reseta = eta.*(beta-x) - tol_sub
   * reszet = zet*z_mma - tol_sub
   */
  Epetra_Vector resXsi(x_->Map());
  Epetra_Vector resEta(x_->Map());

  double* resxsiptr = resXsi.Values();
  double* resetaptr = resEta.Values();

  xsi = xsi_->Values();
  eta = eta_->Values();
  x = x_mma_->Values();
  double* alpha = alpha_->Values();
  double* beta = beta_->Values();

  for (int i=0;i<n_loc_;i++)
  {
    *resxsiptr = *xsi*(*x-*alpha) - tol_sub;
    *resetaptr = *eta*(*beta-*x) - tol_sub;

    resxsiptr++;
    resetaptr++;
    xsi++;
    eta++;
    x++;
    alpha++;
    beta++;
  }


  double reszet = zet_*z_mma_ - tol_sub;


  resnorm = Res2Norm(&resX,&resXsi,&resEta,&resy,&resmu,&reslam,&res,&resz,&reszet);
  resinf = ResInfNorm(&resX,&resXsi,&resEta,&resy,&resmu,&reslam,&res,&resz,&reszet);
}


double OPTI::GCMMA::Res2Norm(
    Epetra_Vector* res1,
    Epetra_Vector* res2,
    Epetra_Vector* res3,
    Epetra_SerialDenseVector* res4,
    Epetra_SerialDenseVector* res5,
    Epetra_SerialDenseVector* res6,
    Epetra_SerialDenseVector* res7,
    double* res8,
    double* res9
)
{
  double resnorm = *res8**res8 + *res9**res9;
  resnorm += std::pow(res4->Norm2(),2) + pow(res5->Norm2(),2) + pow(res6->Norm2(),2) + pow(res7->Norm2(),2);

  double locnorm = 0.0;
  res1->Norm2(&locnorm);
  resnorm += locnorm*locnorm;
  res2->Norm2(&locnorm);
  resnorm += locnorm*locnorm;
  res3->Norm2(&locnorm);
  resnorm += locnorm*locnorm;

  return sqrt(resnorm);
}


double OPTI::GCMMA::ResInfNorm(
    Epetra_Vector* res1,
    Epetra_Vector* res2,
    Epetra_Vector* res3,
    Epetra_SerialDenseVector* res4,
    Epetra_SerialDenseVector* res5,
    Epetra_SerialDenseVector* res6,
    Epetra_SerialDenseVector* res7,
    double* res8,
    double* res9
)
{
  double locnorm = 0.0;
  double resinf = std::max(abs(*res8),abs(*res9));

  resinf = std::max(resinf,res4->NormInf());
  resinf = std::max(resinf,res5->NormInf());
  resinf = std::max(resinf,res6->NormInf());
  resinf = std::max(resinf,res7->NormInf());

  res1->NormInf(&locnorm);
  resinf = std::max(locnorm,resinf);
  res2->NormInf(&locnorm);
  resinf = std::max(locnorm,resinf);
  res3->NormInf(&locnorm);
  resinf = std::max(locnorm,resinf);

  return resinf;
}


void OPTI::GCMMA::Update(
)
{
  // helper vectors
  Epetra_Vector uxinv1(x_->Map());
  Epetra_Vector xlinv1(x_->Map());

  double* asymp_max = asymp_max_->Values();
  double* asymp_min = asymp_min_->Values();
  double* x = x_mma_->Values();

  double* uxinv1ptr = uxinv1.Values();
  double* xlinv1ptr = xlinv1.Values();


  for (int i=0;i<n_loc_;i++)
  {
    *uxinv1ptr = 1.0/(*asymp_max-*x);
    *xlinv1ptr = 1.0/(*x-*asymp_min);

    asymp_max++;
    asymp_min++;
    x++;
    uxinv1ptr++;
    xlinv1ptr++;
  }

  // set new approximation of objective value
  double value = 0.0;

  uxinv1.Dot(*p0_,&value);
  obj_appr_ = r0_ + value;

  xlinv1.Dot(*q0_,&value);
  obj_appr_ += value;


  /*
   * set new approximation of constraints
   *
   * constr_app = b + P^T/ux + Q^T/xl
   */
  *constr_appr_ = *b_;
  constr_appr_->Scale(-1.0);

  double* constr = constr_appr_->Values();
  for (int i=0;i<m_;i++)
  {
    uxinv1.Dot(*(*P_)(i),&value);
    *constr += value;

    xlinv1.Dot(*(*Q_)(i),&value);
    *constr += value;

    constr++;
  }

//  cout << "obj_appr is " << obj_appr_ << endl;
//  cout << "constr_appr is " << *constr_appr_ << endl;
}



void OPTI::GCMMA::Output()
{
  // step number and time
  output_->NewStep(total_iter_,(double)total_iter_);
  // velocity/pressure vector
  output_->WriteVector("x_mma",x_mma_);

  if ((DRT::INPUT::IntegralValue<bool>(params_,"GMSH_OUTPUT")==true) and (total_iter_%upres_ == 0))
    OutputToGmsh(); // TODO look at this: total_iter_,false);

  // write domain decomposition for visualization (only once!)
  if (total_iter_==upres_)
    output_->WriteElementData(true);

  output_->WriteVector("x", x_);
  output_->WriteVector("x_old", x_old_);
  output_->WriteVector("x_old2",x_old2_);

  output_->WriteDouble("obj",obj_);
  output_->WriteVector("obj_deriv",obj_deriv_);

  double* constr = constr_->Values();
  for (int i=0;i<m_;i++)
  {
    std::ostringstream s; s << i;
    std::string name = "constr" + s.str();
    output_->WriteDouble(name,*constr);
    constr++;
  }
  output_->WriteVector("constr_deriv",constr_deriv_);

  return;
}


void OPTI::GCMMA::OutputToGmsh()
{
  // turn on/off screen output for writing process of Gmsh postprocessing file
  const bool screen_out = false;

  // create Gmsh postprocessing file
  const std::string filename = IO::GMSH::GetNewFileNameAndDeleteOldFiles("optimization_field", total_iter_, 500, screen_out, discret_->Comm().MyPID());
  std::ofstream gmshfilecontent(filename.c_str());

  {
    // add 'View' to Gmsh postprocessing file
    gmshfilecontent << "View \" " << "X \" {" << endl;
    // draw scalar field 'Phinp' for every element
    IO::GMSH::ScalarFieldToGmsh(discret_,x_,gmshfilecontent);
    gmshfilecontent << "};" << endl;
  }

  gmshfilecontent.close();
  if (screen_out) std::cout << " done" << endl;
}




