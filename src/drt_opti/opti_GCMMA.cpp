/*---------------------------------------------------------------------*/
/*! \file

\brief Optimization algorithm Globally Convergent Method of Moving Asyptotes (GCMMA)


\level 3

*/
/*---------------------------------------------------------------------*/


#include "opti_GCMMA.H"

#include "topopt_utils.H"
#include "../drt_inpar/inpar_topopt.H"
#include "../drt_inpar/inpar_parameterlist_utils.H"
#include "../drt_io/io.H"
#include "../drt_io/io_control.H"
#include "../drt_io/io_gmsh.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../linalg/linalg_utils_sparse_algebra_math.H"
#include "../headers/definitions.h"
#include <Epetra_SerialDenseSolver.h>
#include <stdio.h>
#include <stdlib.h>



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
OPTI::GCMMA::GCMMA(Teuchos::RCP<DRT::Discretization> discret, const Teuchos::ParameterList& params,
    Teuchos::RCP<Epetra_Vector> x, int numConstraints, Teuchos::RCP<Epetra_Vector> x_min,
    Teuchos::RCP<Epetra_Vector> x_max, Teuchos::RCP<IO::DiscretizationWriter>& output)
    : discret_(discret),
      params_(params.sublist("TOPOLOGY OPTIMIZER")),
      total_iter_(0),
      outer_iter_(0),
      inner_iter_(0),
      max_total_iter_(params_.get<int>("MAX_ITER")),
      max_inner_iter_(params_.get<int>("MAX_INNER_ITER")),
      max_outer_iter_(params_.get<int>("MAX_GRAD_ITER")),
      max_sub_iter_(params_.get<int>("MAX_SUB_ITER")),
      max_inner_sub_iter_(params_.get<int>("MAX_INNER_SUB_ITER")),
      m_(numConstraints),
      n_loc_(x->MyLength()),
      n_(x->GlobalLength()),
      dens_type_(DRT::INPUT::IntegralValue<INPAR::TOPOPT::DensityField>(params, "DENS_TYPE")),
      solver_sub_(DRT::INPUT::IntegralValue<INPAR::TOPOPT::SolverType>(params, "GCMMA_SOLVER")),
      x_(Teuchos::rcp(new Epetra_Vector(*x))),
      x_old_(Teuchos::rcp(new Epetra_Vector(*x))),
      x_old2_(Teuchos::rcp(new Epetra_Vector(*x))),
      x_mma_(Teuchos::rcp(new Epetra_Vector(*x))),
      obj_(0.0),
      obj_deriv_(Teuchos::rcp(new Epetra_Vector(x->Map()))),
      obj_appr_(0.0),
      constr_(Teuchos::rcp(new Epetra_SerialDenseVector(m_))),
      constr_deriv_(Teuchos::rcp(new Epetra_MultiVector(x->Map(), m_))),
      constr_appr_(Teuchos::rcp(new Epetra_SerialDenseVector(m_))),
      p0_(Teuchos::rcp(new Epetra_Vector(x->Map()))),
      q0_(Teuchos::rcp(new Epetra_Vector(x->Map()))),
      r0_(0.0),
      P_(Teuchos::rcp(new Epetra_MultiVector(x->Map(), m_))),
      Q_(Teuchos::rcp(new Epetra_MultiVector(x->Map(), m_))),
      b_(Teuchos::rcp(new Epetra_SerialDenseVector(m_))),
      rho0_(0.01),  // initial value never used
      rho_(Teuchos::rcp(new Epetra_SerialDenseVector(m_))),
      rho0min_(params_.get<double>("RHOMIN")),
      rhomin_(Teuchos::rcp(new Epetra_SerialDenseVector(m_))),
      rho_fac1_(params_.get<double>("RHO_FAC1")),
      rho_fac2_(params_.get<double>("RHO_FAC2")),
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
      tol_sub_(params_.get<double>("TOL_SUB")),
      tol_sub_adaptiv_(DRT::INPUT::IntegralValue<int>(params_, "TOL_SUB_ADAPTIV")),
      tol_sub_quot_fac_(params_.get<double>("TOL_SUB_QUOT_FAC")),
      tol_sub_min_(params_.get<double>("TOL_SUB_MIN")),
      tol_kkt_(params_.get<double>("TOL_KKT")),
      fac_stepsize_(params_.get<double>("fac_stepsize")),
      asy_fac1_(params_.get<double>("asymptotes_fac1")),
      asy_fac2_(params_.get<double>("asymptotes_fac2")),
      fac_x_bd_(params_.get<double>("fac_x_boundaries")),
      fac_sub_reg_(params_.get<double>("fac_sub_reg")),
      inc2_tol_(params_.get<double>("inc2_tol")),
      gamma_up_(params_.get<double>("GAMMA_UP")),
      gamma_down_(params_.get<double>("GAMMA_DOWN")),
      facmin_(params_.get<double>("FACMIN")),
      s_(Teuchos::rcp(new Epetra_SerialDenseVector(m_))),
      upres_(params_.get<int>("RESULTSEVRY")),
      uprestart_(params.get<int>("RESTARTEVRY")),
      output_(output)
{
  if (m_ > 100)
    dserror(
        "current implementation inefficient for large number of constraints due to array structure "
        "and used solver");

  if ((tol_sub_quot_fac_ < 1e-6) || (tol_sub_quot_fac_ > 1e-1))
    dserror("unsensible quotient between KKT and subproblem tolerances");
  if ((tol_sub_min_ < 1e-12) || (tol_sub_min_ > 1e-3))
    dserror("unsensible minimal tolerance for subproblem");
  if ((inc2_tol_ < 1e-25) || (inc2_tol_ > 1e-5)) dserror("unsensible zero-off-boundary");
  if (fac_stepsize_ > -1.0 || fac_stepsize_ < -1.1) dserror("unsensible step size factor");
  if (asy_fac1_ < 2.0 || asy_fac1_ > 100000.0) dserror("unsensible factor for updating asymptotes");
  if (asy_fac2_ < 1.0e-3 || asy_fac2_ > 1.0e-1)
    dserror("unsensible factor for updating asymptotes");
  if (rho_fac1_ < 1.01 || rho_fac1_ > 100) dserror("unsensible factor for updating rho");
  if (rho_fac2_ < 5 || rho_fac2_ > 5000) dserror("unsensible factor for updating rho");
  if (fac_x_bd_ < 0.01 || fac_x_bd_ > 0.5)
    dserror("unsensible factor for computation of boundaries for optimization variable");
  if (fac_sub_reg_ < 0 || fac_sub_reg_ > 1.0e-1)
    dserror("unsensible regularization factor in setup of subproblem");

  if (x_min == Teuchos::null)
  {
    x_min_ = Teuchos::rcp(new Epetra_Vector(x_->Map(), true));

    if (discret_->Comm().MyPID() == 0)
      printf("WARNING: Initialized lower boundary for optimization variables with zeros\n");
  }
  else
    x_min_ = Teuchos::rcp(new Epetra_Vector(*x_min));

  if (x_max == Teuchos::null)
  {
    x_max_ = Teuchos::rcp(new Epetra_Vector(x_->Map(), false));
    x_max_->PutScalar(1.0);

    if (discret_->Comm().MyPID() == 0)
      printf("WARNING: Initialized upper boundary for optimization variables with ones\n");
  }
  else
    x_max_ = Teuchos::rcp(new Epetra_Vector(*x_max));

  // test case modification
  if ((DRT::INPUT::IntegralValue<INPAR::TOPOPT::OptiCase>(params_, "TESTCASE") ==
          INPAR::TOPOPT::optitest_snake_one_constr) or
      (DRT::INPUT::IntegralValue<INPAR::TOPOPT::OptiCase>(params_, "TESTCASE") ==
          INPAR::TOPOPT::optitest_snake_multiple_constr))
  {
    x_min_->PutScalar(-2.0);
    x_max_->PutScalar(2.0);
  }

  if ((not x_min_->Map().SameAs(x_max_->Map())) or (not x_min_->Map().SameAs(x_->Map())))
    dserror("non-matching maps for boundaries of optimization variables");

  x_diff_ = Teuchos::rcp(new Epetra_Vector(*x_max_));
  x_diff_->Update(-1.0, *x_min_, 1.0);

  // minimal distance of boundaries for x-values
  const double xdiffmin = params_.get<double>("X_DIFF_MIN");
  double* xdiff = x_diff_->Values();
  for (int j = 0; j < n_loc_; j++)
  {
    *xdiff = std::max(*xdiff, xdiffmin);

    xdiff++;
  }

  asymp_min_ = Teuchos::rcp(new Epetra_Vector(*x_min_));
  asymp_max_ = Teuchos::rcp(new Epetra_Vector(*x_max_));

  alpha_ = Teuchos::rcp(new Epetra_Vector(*x_min_));
  beta_ = Teuchos::rcp(new Epetra_Vector(*x_max_));

  double* p1 = a_->Values();
  double* p2 = c_->Values();
  double* p3 = d_->Values();
  double* p4 = rho_->Values();
  double* p5 = rhomin_->Values();

  for (int i = 0; i < m_; i++)
  {
    *p1 = params_.get<double>("a_init");
    *p2 = params_.get<double>("c_init") * a0_;
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
Teuchos::RCP<Epetra_Vector> OPTI::GCMMA::Iterate(double& objective,
    Teuchos::RCP<Epetra_Vector> objectivegrad, Teuchos::RCP<Epetra_SerialDenseVector> constraints,
    Teuchos::RCP<Epetra_MultiVector> constraintsgrad)
{
  InitIter(objective, objectivegrad, constraints, constraintsgrad);

  InitSubSolve();

  SubSolve();

  Update();

  return x_mma_;  // current solution
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void OPTI::GCMMA::InitIter(double& objective, Teuchos::RCP<Epetra_Vector> objectivegrad,
    Teuchos::RCP<Epetra_SerialDenseVector> constraints,
    Teuchos::RCP<Epetra_MultiVector> constraintsgrad)
{
  if ((inner_iter_ == 0) and  // new outer_iter
      (objectivegrad == Teuchos::null or constraintsgrad == Teuchos::null))
    dserror("gradients must be given in new outer iteration");

  // initialization
  if (outer_iter_ == 0)
  {
    // reset of optimization variable for test cases
    if ((DRT::INPUT::IntegralValue<INPAR::TOPOPT::OptiCase>(params_, "TESTCASE") ==
            INPAR::TOPOPT::optitest_snake_one_constr) or
        (DRT::INPUT::IntegralValue<INPAR::TOPOPT::OptiCase>(params_, "TESTCASE") ==
            INPAR::TOPOPT::optitest_snake_multiple_constr))
    {
      if (x_->Map().NumGlobalElements() % 3 != 0) dserror("cannot be tested");

      int l = x_->Map().NumGlobalElements() / 3;

      for (int i = 0; i < l; i++)
      {
        double alphai = (3 * (i + 1) - 2 * l) * M_PI / (6 * l);

        if (x_->Map().MyGID(i)) (*x_)[x_->Map().LID(i)] = cos(alphai + M_PI / 12);
        if (x_->Map().MyGID(i + l)) (*x_)[x_->Map().LID(i + l)] = sin(alphai + M_PI / 12);
        if (x_->Map().MyGID(i + 2 * l))
          (*x_)[x_->Map().LID(i + 2 * l)] = sin(2 * alphai + M_PI / 6);
      }
    }

    // initialization of old values
    *x_old_ = *x_;
    *x_old2_ = *x_;
    *x_mma_ = *x_;
  }

  // reset of objective function, constraints and their derivatives for test cases
  if ((DRT::INPUT::IntegralValue<INPAR::TOPOPT::OptiCase>(params_, "TESTCASE") ==
          INPAR::TOPOPT::optitest_snake_one_constr) or
      (DRT::INPUT::IntegralValue<INPAR::TOPOPT::OptiCase>(params_, "TESTCASE") ==
          INPAR::TOPOPT::optitest_snake_multiple_constr))
  {
    double locobj = 0.0;
    Epetra_SerialDenseVector locconstr(m_);

    int l = x_mma_->Map().NumGlobalElements() / 3;

    double delta = 0.1;
    double alpha = 0.0;
    double locg = 0.0;
    double locx1 = 0.0;
    double locx2 = 0.0;
    double locx3 = 0.0;
    for (int i = 0; i < l; i++)
    {
      alpha = M_PI * (3 * (i + 1) - 2 * l) / (6 * l);

      if (x_->Map().MyGID(i))
      {
        int lid = x_mma_->Map().LID(i);
        double x = (*x_mma_)[lid];
        locobj += cos(alpha) * x;
        locconstr(0) += x * x;
        if (inner_iter_ == 0) (*objectivegrad)[lid] = cos(alpha);
        if (i == 0)
        {
          locg += x * x;
          locx1 = x;
        }
      }

      if (x_->Map().MyGID(i + l))
      {
        int lid = x_mma_->Map().LID(i + l);
        double x = (*x_mma_)[lid];
        locobj += sin(alpha) * x;
        locconstr(0) += x * x;
        if (inner_iter_ == 0) (*objectivegrad)[lid] = sin(alpha);
        if (i == 0)
        {
          locg += x * x;
          locx2 = x;
        }
      }

      if (x_->Map().MyGID(i + 2 * l))
      {
        int lid = x_mma_->Map().LID(i + 2 * l);
        double x = (*x_mma_)[lid];
        locobj += -0.1 * x;
        locconstr(0) += x * x;
        if (inner_iter_ == 0) (*objectivegrad)[lid] = -0.1;
        if (i == 0) locx3 = x;
      }
    }
    discret_->Comm().SumAll(&locobj, &objective, 1);
    discret_->Comm().SumAll(locconstr.Values(), constraints->Values(), m_);


    double g = 0.0;
    discret_->Comm().SumAll(&locg, &g, 1);
    g -= 1;
    g /= delta;

    double h = 0.0;
    double x1 = 0.0;
    double x2 = 0.0;
    double x3 = 0.0;
    discret_->Comm().SumAll(&locx1, &x1, 1);
    discret_->Comm().SumAll(&locx2, &x2, 1);
    discret_->Comm().SumAll(&locx3, &x3, 1);
    h = (x3 - 2 * x1 * x2) / delta;

    (*constraints)(0) -= l + 1.0e-5;

    if (DRT::INPUT::IntegralValue<INPAR::TOPOPT::OptiCase>(params_, "TESTCASE") ==
        INPAR::TOPOPT::optitest_snake_multiple_constr)
    {
      (*constraints)(1) = +g + pow(g, 7) - 2 - 1.0e-5;
      (*constraints)(2) = -g - pow(g, 7) - 2 - 1.0e-5;
      (*constraints)(3) = +h + pow(h, 7) - 2 - 1.0e-5;
      (*constraints)(4) = -h - pow(h, 7) - 2 - 1.0e-5;
    }

    if (inner_iter_ == 0)
    {
      (*constraintsgrad)(0)->Update(2.0, *x_mma_, 0.0);

      if (DRT::INPUT::IntegralValue<INPAR::TOPOPT::OptiCase>(params_, "TESTCASE") ==
          INPAR::TOPOPT::optitest_snake_multiple_constr)
      {
        if (x_->Map().MyGID(0))
        {
          int lid = x_mma_->Map().LID(0);
          (*(*constraintsgrad)(1))[lid] = +2 * x1 / delta * (1 + 7 * pow(g, 6));
          (*(*constraintsgrad)(2))[lid] = -2 * x1 / delta * (1 + 7 * pow(g, 6));
          (*(*constraintsgrad)(3))[lid] = -2 * x2 / delta * (1 + 7 * pow(h, 6));
          (*(*constraintsgrad)(4))[lid] = +2 * x2 / delta * (1 + 7 * pow(h, 6));
        }

        if (x_->Map().MyGID(l))
        {
          int lid = x_mma_->Map().LID(l);
          (*(*constraintsgrad)(1))[lid] = +2 * x2 / delta * (1 + 7 * pow(g, 6));
          (*(*constraintsgrad)(2))[lid] = -2 * x2 / delta * (1 + 7 * pow(g, 6));
          (*(*constraintsgrad)(3))[lid] = -2 * x1 / delta * (1 + 7 * pow(h, 6));
          (*(*constraintsgrad)(4))[lid] = +2 * x1 / delta * (1 + 7 * pow(h, 6));
        }

        if (x_->Map().MyGID(2 * l))
        {
          int lid = x_mma_->Map().LID(2 * l);
          (*(*constraintsgrad)(3))[lid] = +1.0 / delta * (1 + 7 * pow(h, 6));
          (*(*constraintsgrad)(4))[lid] = -1.0 / delta * (1 + 7 * pow(h, 6));
        }
      }
    }
  }

  // update counters
  inner_iter_++;
  total_iter_++;

  // new outer iter -> values and derivatives, update asymptotes and re-initialize rho
  if (inner_iter_ == 1)
  {
    // update counters
    outer_iter_++;

    // reset constraints and objective values
    *constr_ = *constraints;
    *constr_deriv_ = *constraintsgrad;
    obj_ = objective;
    *obj_deriv_ = *objectivegrad;
    //    std::cout.precision(16);
    //    std::cout << "new obj is " << obj_ << std::endl;
    //    std::cout << "new obj deriv is " << *obj_deriv_ << std::endl;
    //    std::cout << "new constr are " << *constr_ << std::endl;
    //    std::cout << "new constr deriv are " << *constr_deriv_ << std::endl;
    //
    //    std::ostringstream filename1;
    //    std::ostringstream filename2;
    //    std::ostringstream filename3;
    //    std::ostringstream filename4;
    //    filename1 << "/home/winklmaier/workspace/output/1in2out/failing_constraint/obj_" <<
    //    total_iter_ << ".mtl"; filename2 <<
    //    "/home/winklmaier/workspace/output/1in2out/failing_constraint/obj_der_" << total_iter_ <<
    //    ".mtl"; filename3 <<
    //    "/home/winklmaier/workspace/output/1in2out/failing_constraint/constr_" << total_iter_ <<
    //    ".mtl"; filename4 <<
    //    "/home/winklmaier/workspace/output/1in2out/failing_constraint/constr_der_" << total_iter_
    //    << ".mtl"; LINALG::PrintVectorInMatlabFormat(filename2.str(),*obj_deriv_);
    //    LINALG::PrintVectorInMatlabFormat(filename4.str(),*(*constr_deriv_)(0));
    //    LINALG::PrintSerialDenseMatrixInMatlabFormat(filename3.str(),*constr_);
    //
    //    std::ofstream f1;
    //    f1.precision(16);
    //    f1.open(filename1.str().c_str(),std::fstream::trunc);
    //    f1 << obj_;
    //    f1.close();

    // reset optimization variables
    *x_old2_ = *x_old_;
    *x_old_ = *x_;
    *x_ = *x_mma_;

    Asymptotes();

    InitRho();
  }
  else  // new inner iter -> new values, update rho
  {
    //    std::cout.precision(16);
    //    std::cout << "new obj is " << objective << std::endl;
    //    std::cout << "new constr are " << *constraints << std::endl;
    //
    //    std::ostringstream filename1;
    //    std::ostringstream filename2;
    //    filename1 << "/home/winklmaier/workspace/output/1in2out/failing_constraint/obj_" <<
    //    total_iter_ << ".mtl"; filename2 <<
    //    "/home/winklmaier/workspace/output/1in2out/failing_constraint/constr_" << total_iter_ <<
    //    ".mtl"; LINALG::PrintSerialDenseMatrixInMatlabFormat(filename2.str(),*constraints);
    //
    //    std::ofstream f1;
    //    f1.open(filename1.str().c_str(),std::fstream::ate | std::fstream::app);
    //    f1 << objective;
    //    f1.close();

    UpdateRho(objective, constraints);
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void OPTI::GCMMA::Asymptotes()
{
  // simplified asymptotes in first two iterations
  if (outer_iter_ < 3)
  {
    *asymp_min_ = *x_;
    *asymp_max_ = *x_;

    asymp_min_->Update(-0.5, *x_diff_, 1.0);
    asymp_max_->Update(+0.5, *x_diff_, 1.0);
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
     * fac = 1.0 / 0.7 / 2.3 depending on sign of (x-x_old)(x_old-x_old2)
     */
    double* xval = x_->Values();
    double* xold = x_old_->Values();
    double* xold2 = x_old2_->Values();
    double* asy_min = asymp_min_->Values();
    double* asy_max = asymp_max_->Values();
    double* xdiff = x_diff_->Values();

    for (int i = 0; i < n_loc_; i++)
    {
      double val = (*xval - *xold) * (*xold - *xold2);

      double fac = 1.0;
      if (val < -inc2_tol_)
        fac = gamma_down_;
      else if (val > inc2_tol_)
        fac = gamma_up_;

      *asy_min = *xval - fac * (*xold - *asy_min);
      *asy_max = *xval + fac * (*asy_max - *xold);

      *asy_min = std::max(*asy_min, *xval - asy_fac1_ * *xdiff);
      *asy_min = std::min(*asy_min, *xval - asy_fac2_ * *xdiff);

      *asy_max = std::min(*asy_max, *xval + asy_fac1_ * *xdiff);
      *asy_max = std::max(*asy_max, *xval + asy_fac2_ * *xdiff);

      xval++;
      xold++;
      xold2++;
      asy_min++;
      asy_max++;
      xdiff++;
    }
  }
  //  std::cout << "low is " << *asymp_min_ << std::endl;
  //  std::cout << "upp is " << *asymp_max_ << std::endl;
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

  for (int i = 0; i < n_loc_; i++)
  {
    rho0loc += abs(*obj_deriv) * *xdiff;

    obj_deriv++;
    xdiff++;
  }

  // communicate values
  discret_->Comm().SumAll(&rho0loc, &rho0_, 1);

  // apply minimal value
  rho0_ = std::max(rho0min_, rho0_ / (10 * n_));


  // set rho of constraints
  Epetra_SerialDenseVector rholoc(m_);
  double* rho = rholoc.Values();

  for (int j = 0; j < m_; j++)
  {
    *rho = 0.0;

    Epetra_Vector constr_der(View, *constr_deriv_, j);
    double* constr_deriv = constr_der.Values();
    xdiff = x_diff_->Values();

    for (int i = 0; i < n_loc_; i++)
    {
      *rho += abs(*constr_deriv) * *xdiff;

      constr_deriv++;
      xdiff++;
    }
    rho++;
  }

  // communicate values
  discret_->Comm().SumAll(rholoc.Values(), rho_->Values(), m_);

  // apply minimal value
  rho = rho_->Values();
  double* rhomin = rhomin_->Values();

  for (int j = 0; j < m_; j++)
  {
    *rho = std::max(*rhomin, *rho / (10 * n_));

    rho++;
    rhomin++;
  }
  //  std::cout << "rho0 is " << rho0_ << std::endl;
  //  std::cout << "rho is " << *rho_ << std::endl;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void OPTI::GCMMA::UpdateRho(double& objective, Teuchos::RCP<Epetra_SerialDenseVector> constraints)
{
  /*
   * fac = (x_mma-x)^2*(asy_max-asy_min)/ ((asy_max-x_mma)*(x_mma-asy_min)/xdiff)
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

  for (int i = 0; i < n_loc_; i++)
  {
    facloc += (*xmma - *x) / (*asy_max - *xmma) * (*xmma - *x) / (*xmma - *asy_min) *
              (*asy_max - *asy_min) / (*xdiff);

    x++;
    xmma++;
    asy_max++;
    asy_min++;
    xdiff++;
  }

  // communicate values
  discret_->Comm().SumAll(&facloc, &fac, 1);

  fac = std::max(fac, facmin_);

  if (objective > obj_appr_ + 0.5 * tol_sub_)
  {
    rho0_ = std::min(rho_fac1_ * (rho0_ + (objective - obj_appr_) / fac), rho_fac2_ * rho0_);
  }

  double* constr = constraints->Values();
  double* constr_appr = constr_appr_->Values();
  double* rho = rho_->Values();

  for (int i = 0; i < m_; i++)
  {
    if (*constr > *constr_appr + 0.5 * tol_sub_)
    {
      *rho = std::min(rho_fac1_ * ((*rho) + (*constr - *constr_appr) / fac), rho_fac2_ * (*rho));
    }

    constr++;
    constr_appr++;
    rho++;
  }
  //  std::cout << "rho0 is " << rho0_ << std::endl;
  //  std::cout << "rho is " << *rho_ << std::endl;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool OPTI::GCMMA::Converged(double& objective, Teuchos::RCP<Epetra_Vector> objectivegrad,
    Teuchos::RCP<Epetra_SerialDenseVector> constraints,
    Teuchos::RCP<Epetra_MultiVector> constraintsgrad)
{
  // reset of objective function, constraints and their derivatives for test cases
  if ((DRT::INPUT::IntegralValue<INPAR::TOPOPT::OptiCase>(params_, "TESTCASE") ==
          INPAR::TOPOPT::optitest_snake_one_constr) or
      (DRT::INPUT::IntegralValue<INPAR::TOPOPT::OptiCase>(params_, "TESTCASE") ==
          INPAR::TOPOPT::optitest_snake_multiple_constr))
  {
    double locobj = 0.0;
    Epetra_SerialDenseVector locconstr(m_);

    int l = x_mma_->Map().NumGlobalElements() / 3;

    double delta = 0.1;
    double alpha = 0.0;
    double locg = 0.0;
    double locx1 = 0.0;
    double locx2 = 0.0;
    double locx3 = 0.0;
    for (int i = 0; i < l; i++)
    {
      alpha = M_PI * (3 * (i + 1) - 2 * l) / (6 * l);

      if (x_->Map().MyGID(i))
      {
        int lid = x_mma_->Map().LID(i);
        double x = (*x_mma_)[lid];
        locobj += cos(alpha) * x;
        locconstr(0) += x * x;
        if (i == 0)
        {
          locg += x * x;
          locx1 = x;
        }
      }

      if (x_->Map().MyGID(i + l))
      {
        int lid = x_mma_->Map().LID(i + l);
        double x = (*x_mma_)[lid];
        locobj += sin(alpha) * x;
        locconstr(0) += x * x;
        if (i == 0)
        {
          locg += x * x;
          locx2 = x;
        }
      }

      if (x_->Map().MyGID(i + 2 * l))
      {
        int lid = x_mma_->Map().LID(i + 2 * l);
        double x = (*x_mma_)[lid];
        locobj += -0.1 * x;
        locconstr(0) += x * x;
        if (i == 0) locx3 = x;
      }
    }
    discret_->Comm().SumAll(&locobj, &objective, 1);
    discret_->Comm().SumAll(locconstr.Values(), constraints->Values(), m_);

    double g = 0.0;
    discret_->Comm().SumAll(&locg, &g, 1);
    g -= 1;
    g /= delta;

    double h = 0.0;
    double x1 = 0.0;
    double x2 = 0.0;
    double x3 = 0.0;
    discret_->Comm().SumAll(&locx1, &x1, 1);
    discret_->Comm().SumAll(&locx2, &x2, 1);
    discret_->Comm().SumAll(&locx3, &x3, 1);
    h = (x3 - 2 * x1 * x2) / delta;

    (*constraints)(0) -= l + 1.0e-5;

    if (DRT::INPUT::IntegralValue<INPAR::TOPOPT::OptiCase>(params_, "TESTCASE") ==
        INPAR::TOPOPT::optitest_snake_multiple_constr)
    {
      (*constraints)(1) = +g + pow(g, 7) - 2 - 1.0e-5;
      (*constraints)(2) = -g - pow(g, 7) - 2 - 1.0e-5;
      (*constraints)(3) = +h + pow(h, 7) - 2 - 1.0e-5;
      (*constraints)(4) = -h - pow(h, 7) - 2 - 1.0e-5;
    }


    int dummy = 0;
    if ((InnerConvergence(objective, constraints, dummy) == true) and (outer_iter_ > 0))
    {
      for (int i = 0; i < l; i++)
      {
        alpha = M_PI * (3 * (i + 1) - 2 * l) / (6 * l);

        if (x_->Map().MyGID(i))
        {
          int lid = x_mma_->Map().LID(i);
          (*objectivegrad)[lid] = cos(alpha);
        }

        if (x_->Map().MyGID(i + l))
        {
          int lid = x_mma_->Map().LID(i + l);
          (*objectivegrad)[lid] = sin(alpha);
        }

        if (x_->Map().MyGID(i + 2 * l))
        {
          int lid = x_mma_->Map().LID(i + 2 * l);
          (*objectivegrad)[lid] = -0.1;
        }

        (*constraintsgrad)(0)->Update(2.0, *x_mma_, 0.0);
      }

      if (DRT::INPUT::IntegralValue<INPAR::TOPOPT::OptiCase>(params_, "TESTCASE") ==
          INPAR::TOPOPT::optitest_snake_multiple_constr)
      {
        if (x_->Map().MyGID(0))
        {
          int lid = x_mma_->Map().LID(0);
          (*(*constraintsgrad)(1))[lid] = +2 * x1 / delta * (1 + 7 * pow(g, 6));
          (*(*constraintsgrad)(2))[lid] = -2 * x1 / delta * (1 + 7 * pow(g, 6));
          (*(*constraintsgrad)(3))[lid] = -2 * x2 / delta * (1 + 7 * pow(h, 6));
          (*(*constraintsgrad)(4))[lid] = +2 * x2 / delta * (1 + 7 * pow(h, 6));
        }

        if (x_->Map().MyGID(l))
        {
          int lid = x_mma_->Map().LID(l);
          (*(*constraintsgrad)(1))[lid] = +2 * x2 / delta * (1 + 7 * pow(g, 6));
          (*(*constraintsgrad)(2))[lid] = -2 * x2 / delta * (1 + 7 * pow(g, 6));
          (*(*constraintsgrad)(3))[lid] = -2 * x1 / delta * (1 + 7 * pow(h, 6));
          (*(*constraintsgrad)(4))[lid] = +2 * x1 / delta * (1 + 7 * pow(h, 6));
        }

        if (x_->Map().MyGID(2 * l))
        {
          int lid = x_mma_->Map().LID(2 * l);
          (*(*constraintsgrad)(3))[lid] = +1.0 / delta * (1 + 7 * pow(h, 6));
          (*(*constraintsgrad)(4))[lid] = -1.0 / delta * (1 + 7 * pow(h, 6));
        }
      }
    }
  }


  bool converged = false;
  int num = -1;  // number of failed conditions in inner iteration


  if (discret_->Comm().MyPID() == 0)
  {
    printf("Checking convergence of optimization algorithm GCMMA\n");
    printf("+-----------------------+-----------+-----------+\n");
    printf("|     Condition         |   Value   | Max-value |\n");
  }

  if ((InnerConvergence(objective, constraints, num) == true) and
      (outer_iter_ > 0))  // new outer iteration
  {
    inner_iter_ = 0;

    Epetra_Vector inc(*x_mma_);
    inc.Update(-1.0, *x_, 1.0);

    double inc2norm = 0.0;
    inc.Norm2(&inc2norm);
    double incinfnorm = 0.0;
    inc.NormInf(&incinfnorm);

    // check if outer iteration converged
    if ((KKTCond(objective, objectivegrad, constraints, constraintsgrad)) and (inc2norm < tol_kkt_))
      converged = true;

    if (discret_->Comm().MyPID() == 0)
    {
      printf("| Increment [L2-norm]   |%10.3E |%10.3E |\n", inc2norm, tol_kkt_);
      printf("| Increment [LInf-norm] |%10.3E |%10.3E |\n", incinfnorm, tol_kkt_);
    }
  }
  else  // new inner iteration -> no global convergence
    converged = false;

  if (discret_->Comm().MyPID() == 0)
  {
    printf("+-----------------------+-----------+-----------+\n");
    printf("| Total iteration       |      %4d |      %4d |\n", total_iter_ + 1, max_total_iter_);

    if (inner_iter_ == 0)
      printf("| Outer iteration       |      %4d |      %4d |\n", outer_iter_ + 1, max_outer_iter_);
    else
      printf("| Outer iteration       |      %4d |      %4d |\n", outer_iter_, max_outer_iter_);

    printf("| Inner iteration       |      %4d |      %4d |\n", inner_iter_, max_inner_iter_);

    if (num > 0) printf("| Failing conditions    |      %4d |      %4d |\n", num, m_ + 1);

    printf("+-----------------------+-----------+-----------+\n");
  }

  // stop if total iteration counter or outer iteration counter reaches their
  // respective maximum number of iterations
  if ((total_iter_ == max_total_iter_) or (outer_iter_ == max_outer_iter_))
  {
    if (discret_->Comm().MyPID() == 0)
      printf("WARNING: GCMMA optimization algorithm did not converge\n");

    converged = true;
  }



  return converged;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void OPTI::GCMMA::FinishIteration(
    double& objective, Teuchos::RCP<Epetra_SerialDenseVector> constraints, bool& doGradient)
{
  // reset of objective function and constraints for test cases
  if ((DRT::INPUT::IntegralValue<INPAR::TOPOPT::OptiCase>(params_, "TESTCASE") ==
          INPAR::TOPOPT::optitest_snake_one_constr) or
      (DRT::INPUT::IntegralValue<INPAR::TOPOPT::OptiCase>(params_, "TESTCASE") ==
          INPAR::TOPOPT::optitest_snake_multiple_constr))
  {
    double locobj = 0.0;
    Epetra_SerialDenseVector locconstr(m_);

    int l = x_mma_->Map().NumGlobalElements() / 3;

    double delta = 0.1;
    double alpha = 0.0;
    double locg = 0.0;
    double locx1 = 0.0;
    double locx2 = 0.0;
    double locx3 = 0.0;
    for (int i = 0; i < l; i++)
    {
      alpha = M_PI * (3 * (i + 1) - 2 * l) / (6 * l);

      if (x_->Map().MyGID(i))
      {
        int lid = x_mma_->Map().LID(i);
        double x = (*x_mma_)[lid];
        locobj += cos(alpha) * x;
        locconstr(0) += x * x;
        if (i == 0)
        {
          locg += x * x;
          locx1 = x;
        }
      }

      if (x_->Map().MyGID(i + l))
      {
        int lid = x_mma_->Map().LID(i + l);
        double x = (*x_mma_)[lid];
        locobj += sin(alpha) * x;
        locconstr(0) += x * x;
        if (i == 0)
        {
          locg += x * x;
          locx2 = x;
        }
      }

      if (x_->Map().MyGID(i + 2 * l))
      {
        int lid = x_mma_->Map().LID(i + 2 * l);
        double x = (*x_mma_)[lid];
        locobj += -0.1 * x;
        locconstr(0) += x * x;
        if (i == 0) locx3 = x;
      }
    }
    discret_->Comm().SumAll(&locobj, &objective, 1);
    discret_->Comm().SumAll(locconstr.Values(), constraints->Values(), m_);

    double g = 0.0;
    discret_->Comm().SumAll(&locg, &g, 1);
    g -= 1;
    g /= delta;

    double h = 0.0;
    double x1 = 0.0;
    double x2 = 0.0;
    double x3 = 0.0;
    discret_->Comm().SumAll(&locx1, &x1, 1);
    discret_->Comm().SumAll(&locx2, &x2, 1);
    discret_->Comm().SumAll(&locx3, &x3, 1);
    h = (x3 - 2 * x1 * x2) / delta;

    (*constraints)(0) -= l + 1.0e-5;

    if (DRT::INPUT::IntegralValue<INPAR::TOPOPT::OptiCase>(params_, "TESTCASE") ==
        INPAR::TOPOPT::optitest_snake_multiple_constr)
    {
      (*constraints)(1) = +g + pow(g, 7) - 2 - 1.0e-5;
      (*constraints)(2) = -g - pow(g, 7) - 2 - 1.0e-5;
      (*constraints)(3) = +h + pow(h, 7) - 2 - 1.0e-5;
      (*constraints)(4) = -h - pow(h, 7) - 2 - 1.0e-5;
    }
  }


  // gradient required for new outer iterations
  if (outer_iter_ == 0) doGradient = true;

  int dummy = 0;
  if (InnerConvergence(objective, constraints, dummy))
    doGradient = true;
  else
    doGradient = false;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool OPTI::GCMMA::KKTCond(double& objective, Teuchos::RCP<Epetra_Vector> objectivegrad,
    Teuchos::RCP<Epetra_SerialDenseVector> constraints,
    Teuchos::RCP<Epetra_MultiVector> constraintsgrad)
{
  /*
   * Global residual of primal dual interior point algorithm consists of:
   *
   * resx = eta - xsi + dF/dx^T * lam + df_0/dx
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
  resX.Update(-1.0, *xsi_, 1.0);
  resX.Update(+1.0, *eta_, 1.0);

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

  for (int i = 0; i < m_; i++)
  {
    Epetra_Vector constr_deriv(View, *constraintsgrad, i);
    resX.Update(*lam, constr_deriv, 1.0);

    *resyptr = *c + *d * *y - *mu - *lam;
    resz -= *a * *lam;
    *reslamptr = *constr - z_mma_ * *a - *y + *s;
    *resmuptr = *mu * *y;
    *resptr = *lam * *s;

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

  for (int i = 0; i < n_loc_; i++)
  {
    *resXsiptr = *xsi * (*x - *xmin);
    *resEtaptr = *eta * (*xmax - *x);

    resXsiptr++;
    resEtaptr++;
    xsi++;
    eta++;
    x++;
    xmax++;
    xmin++;
  }

  double reszet = zet_ * z_mma_;

  double resnorm = 0.0;
  resnorm = Res2Norm(&resX, &resXsi, &resEta, &resy, &resmu, &reslam, &res, &resz, &reszet);

  double resinf = 0.0;
  resinf = ResInfNorm(&resX, &resXsi, &resEta, &resy, &resmu, &reslam, &res, &resz, &reszet);

  if (discret_->Comm().MyPID() == 0)
  {
    printf("+-----------------------+-----------+-----------+\n");
    printf("| Res/KKT [L2-norm]     |%10.3E |%10.3E |\n", resnorm, tol_kkt_);
    printf("| Res/KKT [LInf-norm]   |%10.3E |%10.3E |\n", resinf, tol_kkt_);
  }

  if (resnorm < tol_kkt_) return true;

  // adopt tolerance of subproblem to improve convergence
  // Comparison of tol_sub and kkt-residual:
  // - if subproblem tolerance is too small, convergence is slow
  // - if it is too large (= near kkt residuum), no more convergence
  // Moreover, it must not become too small. Otherwise subproblem will not converge below tol
  if ((tol_sub_adaptiv_ == true) && (tol_sub_ > resnorm * tol_sub_quot_fac_))
    tol_sub_ = std::max(0.1 * tol_sub_, tol_sub_min_);

  return false;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool OPTI::GCMMA::InnerConvergence(
    double& objective, Teuchos::RCP<Epetra_SerialDenseVector> constraints, int& numNotFinished)
{
  // initially new outer iter
  if (outer_iter_ == 0) return true;

  if (inner_iter_ == max_inner_iter_)
  {
    if (discret_->Comm().MyPID() == 0)
      printf("WARNING: inner GCMMA optimization loop did not converge\n");

    return true;
  }

  // if all approximations (constraints and objective) are smaller than the
  // values itself (modulo a small tolerance), then inner iteration is finished
  bool finished = true;
  numNotFinished = 0;

  //  std::cout << "obj: diff is " << objective-obj_appr_-tol_sub_ << std::endl;
  //  std::cout << "objective is " << objective << std::endl;
  //  std::cout << "approximated objective is " << obj_appr_ << std::endl;

  if (obj_appr_ + tol_sub_ < objective)
  {
    finished = false;
    numNotFinished++;
  }

  double* constr_appr = constr_appr_->Values();
  double* constr = constraints->Values();

  for (int i = 0; i < m_; i++)
  {
    //    std::cout << "constr: diff is " << *constr-*constr_appr-tol_sub_ << std::endl;
    //    std::cout << "constraint is " << (*constraints)(0) << std::endl;
    //    std::cout << "approximated constraint is " << *constr_appr << std::endl;

    if (*constr_appr + tol_sub_ < *constr)
    {
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
  alpha_->Update(1.0 - fac_x_bd_, *asymp_min_, 0.0);
  alpha_->Update(fac_x_bd_, *x_, 1.0);

  beta_->Update(1.0 - fac_x_bd_, *asymp_max_, 0.0);
  beta_->Update(fac_x_bd_, *x_, 1.0);

  double* alpha = alpha_->Values();
  double* beta = beta_->Values();
  double* xmin = x_min_->Values();
  double* xmax = x_max_->Values();

  for (int i = 0; i < n_loc_; i++)
  {
    *alpha = std::max(*alpha, *xmin);
    *beta = std::min(*beta, *xmax);

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

  Epetra_Vector ux1(x_->Map(), false);
  Epetra_Vector ux2(x_->Map(), false);
  Epetra_Vector xl1(x_->Map(), false);
  Epetra_Vector xl2(x_->Map(), false);

  double* ux1ptr = ux1.Values();
  double* ux2ptr = ux2.Values();
  double* xl1ptr = xl1.Values();
  double* xl2ptr = xl2.Values();

  for (int i = 0; i < n_loc_; i++)
  {
    *ux1ptr = *asymp_max - *x;
    *ux2ptr = *ux1ptr * *ux1ptr;
    *xl1ptr = *x - *asymp_min;
    *xl2ptr = *xl1ptr * *xl1ptr;

    ux1ptr++;
    ux2ptr++;
    xl1ptr++;
    xl2ptr++;
    x++;
    asymp_max++;
    asymp_min++;
  }


  /* compute
   * p0 = (max(dJ/dx,0) + fac*p0q0 + rho0/x_diff).*ux2
   * q0 = (max(-dJ/dx,0) + fac*p0q0 + rho0/x_diff).*xl2
   * r0 = J - p0./ux1 - q0./xl1
   *
   * with
   * p0q0 = (max(dJ/dx,0) + (max(-dJ/dx,0)
   * fac << 1
   */
  double* p0 = p0_->Values();
  double* q0 = q0_->Values();
  double* obj_deriv = obj_deriv_->Values();

  double p0q0 = 0.0;  // tmp variable

  ux1ptr = ux1.Values();
  ux2ptr = ux2.Values();
  xl1ptr = xl1.Values();
  xl2ptr = xl2.Values();

  double* xdiff = x_diff_->Values();

  double r0loc = 0.0;

  for (int i = 0; i < n_loc_; i++)
  {
    *p0 = std::max(*obj_deriv, 0.0);
    *q0 = std::max(-*obj_deriv, 0.0);

    p0q0 = *p0 + *q0;
    *p0 = (*p0 + fac_sub_reg_ * p0q0 + rho0_ / (*xdiff)) * *ux2ptr;
    *q0 = (*q0 + fac_sub_reg_ * p0q0 + rho0_ / (*xdiff)) * *xl2ptr;

    r0loc -= *p0 / (*ux1ptr) + *q0 / (*xl1ptr);

    p0++;
    q0++;
    obj_deriv++;
    xdiff++;
    ux1ptr++;
    ux2ptr++;
    xl1ptr++;
    xl2ptr++;
  }
  discret_->Comm().SumAll(&r0loc, &r0_, 1);

  r0_ += obj_;


  /* compute
   * P = (max(dF/dx,0) + fac*pq + rho/x_diff).*ux2
   * Q = (max(-dF/dx,0) + fac*pq + rho/x_diff).*xl2
   * b = F - P^T/ux1 - Q^T/xl1
   *
   * with
   * pq = (max(dF/dx,0) + (max(-dF/dx,0)
   * fac << 1
   */
  Epetra_SerialDenseVector bloc(m_);
  double* b = bloc.Values();

  double* rho = rho_->Values();

  for (int i = 0; i < m_; i++)
  {
    Epetra_Vector p_i(View, *P_, i);
    Epetra_Vector q_i(View, *Q_, i);
    Epetra_Vector constr_deriv(View, *constr_deriv_, i);

    double* p = p_i.Values();
    double* q = q_i.Values();
    double* constr_der = constr_deriv.Values();

    double pq = 0.0;  // tmp variable

    xdiff = x_diff_->Values();

    ux1ptr = ux1.Values();
    ux2ptr = ux2.Values();
    xl1ptr = xl1.Values();
    xl2ptr = xl2.Values();

    for (int j = 0; j < n_loc_; j++)
    {
      *p = std::max(*constr_der, 0.0);
      *q = std::max(-*constr_der, 0.0);

      pq = *p + *q;
      *p = (*p + fac_sub_reg_ * pq + *rho / (*xdiff)) * *ux2ptr;
      *q = (*q + fac_sub_reg_ * pq + *rho / (*xdiff)) * *xl2ptr;

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
  for (int i = 0; i < m_; i++)
  {
    Epetra_Vector p_i(View, *P_, i);
    Epetra_Vector q_i(View, *Q_, i);
    double* p = p_i.Values();
    double* q = q_i.Values();
    ux1ptr = ux1.Values();
    xl1ptr = xl1.Values();
    for (int j = 0; j < n_loc_; j++)
    {
      *b -= *p / (*ux1ptr) + *q / (*xl1ptr);
      p++;
      q++;
      ux1ptr++;
      xl1ptr++;
    }
    b++;
  }
  discret_->Comm().SumAll(bloc.Values(), b_->Values(), m_);
  *b_ += *constr_;
  b_->Scale(-1.0);

  //  std::cout << "init subsolve: " << std::endl;
  //  std::cout << "m is " << m_ << ", n is " << n_ << ", epsimin is " << tol_sub_ << std::endl;
  //  std::cout << "low asy is " << *asymp_min_ << "upp asy is " << *asymp_max_ << std::endl;
  //  std::cout << "alpha is " << *alpha_ << "beta is " << *beta_ << std::endl;
  //  std::cout << "p0 is " << *p0_ << "q0 is " << *q0_ << std::endl;
  //  std::cout << "p is " << *P_ << "q is " << *Q_ << std::endl;
  //  std::cout << "a0 is " << a0_ << ", a is " << *a_ << ", b is " << *b_ << ", c is " << *c_ << ",
  //  d is " << *d_ << std::endl;

  //  const std::string outname(DRT::Problem::Instance()->OutputControlFile()->FileName());
  //
  //  std::ostringstream filename1;
  //  std::ostringstream filename2;
  //  std::ostringstream filename3;
  //  std::ostringstream filename4;
  //  std::ostringstream filename5;
  //  std::ostringstream filename6;
  //  std::ostringstream filename7;
  //  std::ostringstream filename8;
  //  std::ostringstream filename9;
  //  std::ostringstream filename10;
  //
  //  filename1 << outname << "_" << total_iter_ << "_lower.mtl";
  //  filename2 << outname << "_" << total_iter_ << "_upper.mtl";
  //  filename3 << outname << "_" << total_iter_ << "_alpha.mtl";
  //  filename4 << outname << "_" << total_iter_ << "_beta.mtl";
  //  filename5 << outname << "_" << total_iter_ << "_P.mtl";
  //  filename6 << outname << "_" << total_iter_ << "_Q.mtl";
  //  filename7 << outname << "_" << total_iter_ << "_p0.mtl";
  //  filename8 << outname << "_" << total_iter_ << "_q0.mtl";
  //  filename9 << outname << "_" << total_iter_ << "_a.mtl";
  //  filename10 << outname << "_" << total_iter_ << "_d.mtl";
  //
  //  LINALG::PrintVectorInMatlabFormat(filename1.str(),*asymp_min_);
  //  LINALG::PrintVectorInMatlabFormat(filename2.str(),*asymp_max_);
  //  LINALG::PrintVectorInMatlabFormat(filename3.str(),*alpha_);
  //  LINALG::PrintVectorInMatlabFormat(filename4.str(),*beta_);
  //  LINALG::PrintVectorInMatlabFormat(filename5.str(),*(*P_)(0));
  //  LINALG::PrintVectorInMatlabFormat(filename6.str(),*(*Q_)(0));
  //  LINALG::PrintVectorInMatlabFormat(filename7.str(),*(*p0_)(0));
  //  LINALG::PrintVectorInMatlabFormat(filename8.str(),*(*q0_)(0));
  //  LINALG::PrintSerialDenseMatrixInMatlabFormat(filename9.str(),*a_);
  //  LINALG::PrintSerialDenseMatrixInMatlabFormat(filename10.str(),*d_);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void OPTI::GCMMA::SubSolve()
{
  if (solver_sub_ == INPAR::TOPOPT::gcmma_solver_orig)  // old solver
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
    x_mma_->Update(0.5, *alpha_, 0.0);
    x_mma_->Update(0.5, *beta_, 1.0);

    double* y = y_mma_->Values();
    double* lam = lam_->Values();
    double* mu = mu_->Values();
    double* s = s_->Values();
    double* c = c_->Values();

    for (int i = 0; i < m_; i++)
    {
      *y = 1.0;
      *lam = 1.0;
      *mu = std::max(1.0, *c / 2.0);
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
    for (int i = 0; i < n_loc_; i++)
    {
      *xsi = std::max(1.0, 1.0 / (*x - *alpha));
      *eta = std::max(1.0, 1.0 / (*beta - *x));

      xsi++;
      eta++;
      x++;
      alpha++;
      beta++;
    }


    int total_iter = 0;
    double tol_sub = 1.0;
    bool tol_reached = false;

    while (tol_reached == false)  // this loop has to finish, no max_iter required
    {
      // first step: compute norms of all composed residuals
      double resnorm = 0.0;
      double resinf = 0.0;
      ResApp(resnorm, resinf, tol_sub);

      int inner_iter = 0;

      // iteration loop of inner algorithm
      while (inner_iter < max_sub_iter_)
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
        for (int j = 0; j < m_; j++)
        {
          Epetra_Vector P(View, *P_, j);
          Epetra_Vector Q(View, *Q_, j);

          plam.Update(*lam, P, 1.0);
          qlam.Update(*lam, Q, 1.0);

          lam++;
        }


        /*
         * helper data
         *
         * delx = plam./ux2 - qlam./xl2 - tol_sub/(x-alpha) + tol_sub/(beta-x)
         * diagx = 2*(plam/(ux1^3) + qlamptr/(xl1^3)) + xsi/(x-alpha) + eta/(beta-x)
         */
        Epetra_Vector delx(x_->Map(), false);
        Epetra_Vector diagx(x_->Map(), false);

        Epetra_Vector uxinv1(x_->Map(), false);
        Epetra_Vector xlinv1(x_->Map(), false);

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

        for (int i = 0; i < n_loc_; i++)
        {
          ux1 = *asymp_max - *x;
          xl1 = *x - *asymp_min;

          ux2 = ux1 * ux1;
          xl2 = xl1 * xl1;

          *uxinv1ptr = 1.0 / ux1;
          *xlinv1ptr = 1.0 / xl1;

          double dpsidx = *plamptr / ux2 - *qlamptr / xl2;

          *delxptr = dpsidx - tol_sub / (*x - *alpha) + tol_sub / (*beta - *x);
          *diagxptr = 2 * (*plamptr / (ux2 * ux1) + *qlamptr / (xl2 * xl1)) + *xsi / (*x - *alpha) +
                      *eta / (*beta - *x);

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
        Epetra_MultiVector GG(x_->Map(), m_);

        for (int i = 0; i < m_; i++)
        {
          Epetra_Vector gg(View, GG, i);

          Epetra_Vector P(View, *P_, i);
          Epetra_Vector Q(View, *Q_, i);

          double* ggptr = gg.Values();
          double* pptr = P.Values();
          double* qptr = Q.Values();
          uxinv1ptr = uxinv1.Values();
          xlinv1ptr = xlinv1.Values();

          for (int j = 0; j < n_loc_; j++)
          {
            *ggptr = *pptr * (*uxinv1ptr * *uxinv1ptr) - *qptr * (*xlinv1ptr * *xlinv1ptr);

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

        double delz = a0_ - tol_sub / z_mma_;

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

        for (int i = 0; i < m_; i++)
        {
          *delyptr = *c + *d * *y - *lam - tol_sub / (*y);

          delz -= *a * *lam;

          uxinv1.Dot(*(*P_)(i), &gvec1);
          xlinv1.Dot(*(*Q_)(i), &gvec2);
          *dellamptr = gvec1 + gvec2 - *a * z_mma_ - *y - *b + tol_sub / (*lam);

          *diagyptr = *d + *mu / (*y);

          *diaglamyiptr = *s / (*lam) + 1.0 / (*diagyptr);

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
           * Build final system matrix and solve (size m+1)
           *
           * rhs = bb = [dellam + dely./diagy - GG^T*delx./diagx
           *             delz]
           * sysmat = AA = [GG^T * diagx(diagmat) * GG + diaglamyi (diagmat) || a
           *                a^T                                              || -zet/z_mma]
           * solut = AA^-1 * bb
           */
          Epetra_SerialDenseVector bb(m_ + 1);
          double* bbval = bb.A();

          dellamptr = dellam.Values();
          delyptr = dely.Values();
          diagyptr = diagy.Values();

          for (int i = 0; i < m_; i++)
          {
            double locval = 0.0;
            double val = 0.0;

            Epetra_Vector gg(View, GG, i);
            double* ggptr = gg.Values();
            delxptr = delx.Values();
            diagxptr = diagx.Values();

            for (int j = 0; j < n_loc_; j++)
            {
              locval += *ggptr * *delxptr / (*diagxptr);

              ggptr++;
              delxptr++;
              diagxptr++;
            }
            discret_->Comm().SumAll(&locval, &val, 1);

            *bbval = *dellamptr + *delyptr / (*diagyptr) - val;

            bbval++;
            dellamptr++;
            delyptr++;
            diagyptr++;
          }
          *bbval = delz;


          Epetra_SerialDenseMatrix AA(m_ + 1, m_ + 1);
          double* aa = AA.A();
          a = a_->Values();
          diaglamyiptr = diaglamyi.Values();

          for (int icol = 0; icol < m_; icol++)  // Column loop for m x m Matrix
          {
            // Alam = spdiags(diaglamyi,0,m,m) + GG*spdiags(diagxinv,0,n,n)*GG';
            Epetra_Vector colGG(View, GG, icol);
            double* colgg;

            for (int irow = 0; irow < m_; irow++)  // row loop for m x m Matrix
            {
              colgg = colGG.Values();  // reset here every loop

              Epetra_Vector rowGG(View, GG, irow);
              double* rowgg = rowGG.Values();

              double locval = 0.0;
              diagxptr = diagx.Values();

              // Remark: interpret x'*D*x as weighted dotproduct (x,x) with weights=diag(D)
              for (int j = 0; j < n_loc_; j++)  // dot Product loop for GG*diagpsiinv*GG'.
              {
                locval += *rowgg / (*diagxptr) * *colgg;

                rowgg++;
                diagxptr++;
                colgg++;
              }

              discret_->Comm().SumAll(&locval, aa, 1);

              if (irow == icol) *aa += *diaglamyiptr;

              aa++;
            }

            *aa = *a;

            a++;
            aa++;
            diaglamyiptr++;
          }

          a = a_->Values();
          for (int i = 0; i < m_; i++)
          {
            *aa = *a;

            a++;
            aa++;
          }

          *aa = -zet_ / z_mma_;


          Epetra_SerialDenseVector solut(m_ + 1);

          Epetra_SerialDenseSolver solver;
          solver.SetMatrix(AA);
          solver.SetVectors(solut, bb);
          solver.Solve();


          /*
           * dlam = solut(1:m) // all except for the last entry / first m entries
           * dz = solut(end) // last entry
           * dx = (-GG^T*dlam - delx)./diagx
           */
          double* dlamptr = dlam.Values();
          double* sol = solut.A();

          for (int i = 0; i < m_; i++)
          {
            *dlamptr = *sol;

            dlamptr++;
            sol++;
          }


          dz = *sol;


          dlamptr = dlam.Values();
          for (int i = 0; i < m_; i++)
          {
            Epetra_Vector gg(View, GG, i);
            dx.Update(-*dlamptr, gg, 1.0);

            dlamptr++;
          }

          double* dxptr = dx.Values();
          delxptr = delx.Values();
          diagxptr = diagx.Values();

          for (int i = 0; i < n_loc_; i++)
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
          // see matlab code:
          // diaglamyiinv = eem./diaglamyi;
          // dellamyi = dellam + dely./diagy;
          // Axx = spdiags(diagx,0,n,n) + GG'*spdiags(diaglamyiinv,0,m,m)*GG;
          // azz = zet/z + a'*(a./diaglamyi);
          // axz = -GG'*(a./diaglamyi);
          // bx = delx + GG'*(dellamyi./diaglamyi);
          // bz  = delz - a'*(dellamyi./diaglamyi);
          // AA = [Axx   axz
          //       axz'  azz ];
          // bb = [-bx' -bz]';
          // solut = AA\bb;
          // dx  = solut(1:n);
          // dz = solut(n+1);
          // dlam = (GG*dx)./diaglamyi - dz*(a./diaglamyi) + dellamyi./diaglamyi;
        }


        /*
         * dxsi = -xsi + (tol_sub - xsi*dx)/(x-alpha)
         * deta = -eta + (tol_sub + eta*dx)/(beta-x)
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

        for (int i = 0; i < n_loc_; i++)
        {
          *dxsiptr = -*xsi + (tol_sub - *xsi * *dxptr) / (*x - *alpha);
          *detaptr = -*eta + (tol_sub + *eta * *dxptr) / (*beta - *x);

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
         * dy = - dely./diagy + dlam./diagy
         * dmu = -mu + tol_sub/y - mu.*dy./y
         * ds = -s + tol_sub/lam - s.*dlam./lam
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

        for (int i = 0; i < m_; i++)
        {
          *dyptr = -*delyptr / (*diagyptr) + *dlamptr / (*diagyptr);
          *dmuptr = -*mu + tol_sub / (*y) - (*mu * *dyptr) / (*y);
          *dsptr = -*s + tol_sub / (*lam) - (*s * *dlamptr) / (*lam);

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


        double dzet = -zet_ + tol_sub / z_mma_ - zet_ * dz / z_mma_;


        double val = std::min(dz / z_mma_, dzet / zet_);

        y = y_mma_->Values();
        dyptr = dy.Values();
        lam = lam_->Values();
        dlamptr = dlam.Values();
        mu = mu_->Values();
        dmuptr = dmu.Values();
        s = s_->Values();
        dsptr = ds.Values();

        for (int i = 0; i < m_; i++)
        {
          val = std::min(val, *dyptr / (*y));
          val = std::min(val, *dlamptr / (*lam));
          val = std::min(val, *dmuptr / (*mu));
          val = std::min(val, *dsptr / (*s));

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

        for (int i = 0; i < n_loc_; i++)
        {
          val = std::min(val, *dxsiptr / (*xsi));
          val = std::min(val, *detaptr / (*eta));
          val = std::min(val, *dxptr / (*x - *alpha));
          val = std::min(val, *dxptr / (*x - *beta));

          xsi++;
          dxsiptr++;
          eta++;
          detaptr++;
          dxptr++;
          x++;
          alpha++;
          beta++;
        }

        // min becomes max since fac<0
        val = std::max(1.0, fac_stepsize_ * val);
        double stepsize = 0.0;
        discret_->Comm().MaxAll(&val, &stepsize, 1);
        stepsize = 1.0 / stepsize;
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
        double resnew = 2 * resnorm;


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


        while ((resnorm < resnew) and (it < max_inner_sub_iter_))
        {
          it++;

          x_mma_->Update(stepsize, dx, 1.0, x_mma_old, 0.0);
          xsi_->Update(stepsize, dxsi, 1.0, xsi_old, 0.0);
          eta_->Update(stepsize, deta, 1.0, eta_old, 0.0);

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

          for (int i = 0; i < m_; i++)
          {
            *y = *yptr + stepsize * *dyptr;
            *lam = *lamptr + stepsize * *dlamptr;
            *mu = *muptr + stepsize * *dmuptr;
            *s = *sptr + stepsize * *dsptr;

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

          z_mma_ = z_mma_old + stepsize * dz;
          zet_ = zet_old + stepsize * dzet;


          // compute residuals
          ResApp(resnew, resinf, tol_sub);

          stepsize = stepsize / 2.0;

          if ((it == max_inner_sub_iter_) and (discret_->Comm().MyPID() == 0))
            printf(
                "Reached maximal number of iterations in most inner loop of primal dual interior "
                "point optimization algorithm\n");
        }

        resnorm = resnew;

        // it would be sufficient to compute resinf only once here and not in
        // every iteration as done above. but the effort is neglegible
        if (resinf < 0.9 * tol_sub) break;
      }

      if ((inner_iter == max_sub_iter_) and (discret_->Comm().MyPID() == 0))
        printf(
            "Reached maximal number of iterations in inner loop of primal dual interior point "
            "optimization algorithm\n");

      if (tol_sub > 1.001 * tol_sub_)
        tol_sub *= 0.1;
      else
        tol_reached = true;
    }
  }
  else if (solver_sub_ == INPAR::TOPOPT::gcmma_solver_gauss)  // new solver
  {
    // AB HIER NEUER SOLVER tag newsolver
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
    // set outputs

    if (m_ > 1) dserror("Warning: Algorithm has bug for m>1 until now!");

    x_mma_->Update(0.5, *alpha_, 0.0);
    x_mma_->Update(0.5, *beta_, 1.0);

    double* y = y_mma_->Values();
    double* lam = lam_->Values();
    double* mu = mu_->Values();
    double* s = s_->Values();
    double* c = c_->Values();

    // std::cout <<"c is " << *c << std::endl;

    for (int i = 0; i < m_; i++)
    {
      *y = 1.0;
      *lam = 1.0;
      *mu = std::max(1.0, *c / 2.0);
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
    for (int i = 0; i < n_loc_; i++)
    {
      *xsi = std::max(1.0, 1.0 / (*x - *alpha));
      *eta = std::max(1.0, 1.0 / (*beta - *x));

      xsi++;
      eta++;
      x++;
      alpha++;
      beta++;
    }

    int total_iter = 0;
    double tol_sub = 1.0;
    bool tol_reached = false;

    while (tol_reached == false)  // this loop has to finish, no max_iter required
    {
      // first step: compute norms of all composed residuals
      double resnorm = 0.0;
      double resinf = 0.0;
      ResApp(resnorm, resinf, tol_sub);

      int inner_iter = 0;

      // iteration loop of inner algorithm
      // std::cout << " eine neue innere schleife startet hier " << std::endl;
      while (inner_iter < max_sub_iter_)
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
        for (int j = 0; j < m_; j++)
        {
          Epetra_Vector P(View, *P_, j);
          Epetra_Vector Q(View, *Q_, j);

          plam.Update(*lam, P, 1.0);
          qlam.Update(*lam, Q, 1.0);

          lam++;
        }

        // AB HIER BACKUP ALTES VERFARHEN MANCHE GRÖßEN WERDEN NOCH BENÖTIGT, deshalb ins neue
        // Verfahren gecopypastet
        /*
         * helper data
         *
         * delx = plam./ux2 - qlam./xl2 - tol_sub/(x-alpha) + tol_sub/(beta-x)
         * diagx = 2*(plam/(ux1^3) + qlamptr/(xl1^3)) + xsi/(x-alpha) + eta/(beta-x)
         * xmina = xmma - alfa
         * bminx = beta - xmma
         */
        Epetra_Vector diagx(x_->Map(), false);
        Epetra_Vector uxinv1(x_->Map(), false);
        Epetra_Vector xlinv1(x_->Map(), false);
        Epetra_Vector dpsidx(x_->Map(), false);
        Epetra_Vector psivec(x_->Map(), false);
        Epetra_Vector diagpsi(x_->Map(), false);

        double ux1 = 0.0;
        double xl1 = 0.0;
        double ux2 = 0.0;
        double xl2 = 0.0;

        double* uxinv1ptr = uxinv1.Values();
        double* xlinv1ptr = xlinv1.Values();

        double* x = x_mma_->Values();
        double* asymp_max = asymp_max_->Values();
        double* asymp_min = asymp_min_->Values();
        double* diagxptr = diagx.Values();
        double* dpsidxptr = dpsidx.Values();
        double* psivecptr = psivec.Values();

        xsi = xsi_->Values();
        eta = eta_->Values();
        alpha = alpha_->Values();
        beta = beta_->Values();

        double* plamptr = plam.Values();
        double* qlamptr = qlam.Values();
        double* diagpsiptr = diagpsi.Values();

        // ux1, ux2, xl1, xl2, diagx, dpsidx
        // dpsidx = plam./ux2 - qlam./xl2
        for (int i = 0; i < n_loc_; i++)
        {
          ux1 = *asymp_max - *x;
          xl1 = *x - *asymp_min;
          ux2 = ux1 * ux1;
          xl2 = xl1 * xl1;

          *uxinv1ptr = 1.0 / ux1;
          *xlinv1ptr = 1.0 / xl1;
          *dpsidxptr = *plamptr / ux2 - *qlamptr / xl2;
          *diagxptr = 2 * (*plamptr / (ux2 * ux1) + *qlamptr / (xl2 * xl1));
          *psivecptr = *diagxptr;
          *diagpsiptr = *psivecptr;

          uxinv1ptr++;
          xlinv1ptr++;
          diagxptr++;
          diagpsiptr++;
          psivecptr++;
          dpsidxptr++;
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

        // xmina bminx
        Epetra_Vector xmina(x_->Map(), false);
        Epetra_Vector bminx(x_->Map(), false);

        x = x_mma_->Values();
        double* alpha = alpha_->Values();
        double* beta = beta_->Values();

        double* xminaptr = xmina.Values();
        double* bminxptr = bminx.Values();

        for (int i = 0; i < n_loc_; i++)
        {
          *xminaptr = *x - *alpha;
          *bminxptr = *beta - *x;

          x++;
          alpha++;
          beta++;
          xminaptr++;
          bminxptr++;
        }

        /*
         * helper data
         *
         * GG = P./(ux1^2) - Q./(xl1^2)
         */
        Epetra_MultiVector GG(x_->Map(), m_);

        for (int i = 0; i < m_; i++)
        {
          Epetra_Vector gg(View, GG, i);

          Epetra_Vector P(View, *P_, i);
          Epetra_Vector Q(View, *Q_, i);

          double* ggptr = gg.Values();
          double* pptr = P.Values();
          double* qptr = Q.Values();
          uxinv1ptr = uxinv1.Values();
          xlinv1ptr = xlinv1.Values();

          for (int j = 0; j < n_loc_; j++)
          {
            *ggptr = *pptr * (*uxinv1ptr * *uxinv1ptr) - *qptr * (*xlinv1ptr * *xlinv1ptr);

            ggptr++;
            pptr++;
            qptr++;
            uxinv1ptr++;
            xlinv1ptr++;
          }
        }

        // HIER BEGINN NEUES VERFAHREN

        // diagxsi = spdiags(xsi,0,n,n)
        /* "diagxsi = xsi" ist nur zur Umbenennung gedacht, da man mit dem einen eher den
         * Vektor, mit dem anderen aber eher die Diagonalmatrix assoziiert. Das Ganze kann
         * also ruhig wieder umgebaut werden, falls man es aus Effizienzgründen streichen will.
         */
        Epetra_Vector diagxsi(x_->Map(), false);
        Epetra_Vector diageta(x_->Map(), false);
        Epetra_Vector diagpsiinv(x_->Map(), false);
        Epetra_Vector diagxmina(x_->Map(), false);
        Epetra_Vector diagbminx(x_->Map(), false);

        double* diagxsiptr = diagxsi.Values();
        double* diagetaptr = diageta.Values();
        double* diagpsiinvptr = diagpsiinv.Values();
        double* xsiptr = xsi_->Values();
        double* etaptr = eta_->Values();
        psivecptr = psivec.Values();
        double* diagxminaptr = diagxmina.Values();
        double* diagbminxptr = diagbminx.Values();
        xminaptr = xmina.Values();
        bminxptr = bminx.Values();

        for (int i = 0; i < n_loc_; i++)
        {
          *diagxsiptr = *xsiptr;
          *diagetaptr = *etaptr;
          *diagpsiinvptr = 1.0 / (*psivecptr);
          *diagxminaptr = *xminaptr;
          *diagbminxptr = *bminxptr;

          diagxsiptr++;
          diagetaptr++;
          xsiptr++;
          etaptr++;
          diagpsiinvptr++;
          psivecptr++;
          diagxminaptr++;
          diagbminxptr++;
          xminaptr++;
          bminxptr++;
        }

        Epetra_SerialDenseVector diagd(*d_);
        Epetra_SerialDenseVector diagdinv(*d_);
        Epetra_SerialDenseVector diagmu(*mu_);
        Epetra_SerialDenseVector diags(*s_);
        Epetra_SerialDenseVector diagymma(*y_mma_);
        Epetra_SerialDenseVector diaglam(*lam_);

        double* diagdptr = diagd.Values();
        double* dptr = d_->Values();
        double* diagdinvptr = diagdinv.Values();
        double* diagmuptr = diagmu.Values();
        double* muptr = mu_->Values();
        double* diagsptr = diags.Values();
        double* sptr = s_->Values();
        double* diagymmaptr = diagymma.Values();
        double* ymmaptr = y_mma_->Values();
        double* diaglamptr = diaglam.Values();
        double* lamptr = lam_->Values();

        for (int i = 0; i < m_; i++)
        {
          *diagdptr = *dptr;
          *diagdinvptr = 1.0 / (*dptr);
          *diagmuptr = *muptr;
          *diagsptr = *sptr;
          *diagymmaptr = *ymmaptr;
          *diaglamptr = *lamptr;

          diagdptr++;
          dptr++;
          diagdinvptr++;
          diagmuptr++;
          muptr++;
          diagsptr++;
          sptr++;
          diagymmaptr++;
          ymmaptr++;
          diaglamptr++;
          lamptr++;
        }

        /* Ab hier Hilfsgrößen für rechte Seite
         * und Aufstellen der rechten Seite
         */
        Epetra_Vector delx(x_->Map(), false);
        Epetra_Vector delxsi(x_->Map(), false);
        Epetra_Vector deleta(x_->Map(), false);

        double* delxptr = delx.Values();
        dpsidxptr = dpsidx.Values();
        xsiptr = xsi_->Values();
        etaptr = eta_->Values();

        double* delxsiptr = delxsi.Values();
        double* alphaptr = alpha_->Values();
        double* deletaptr = deleta.Values();
        double* betaptr = beta_->Values();

        double* xmmaptr = x_mma_->Values();

        for (int i = 0; i < n_loc_; i++)
        {
          *delxptr = *dpsidxptr - *xsiptr + *etaptr;
          *delxsiptr = *xsiptr * (*xmmaptr - *alphaptr) - tol_sub;
          *deletaptr = *etaptr * (*betaptr - *xmmaptr) - tol_sub;

          delxptr++;
          dpsidxptr++;
          xsiptr++;
          etaptr++;
          delxsiptr++;
          alphaptr++;
          deletaptr++;
          betaptr++;
          xmmaptr++;
        }

        Epetra_SerialDenseVector dely(*y_mma_);
        Epetra_SerialDenseVector dellam(*y_mma_);
        Epetra_SerialDenseVector delmu(*y_mma_);
        Epetra_SerialDenseVector dels(*y_mma_);

        double* delyptr = dely.Values();
        double* cptr = c_->Values();
        dptr = d_->Values();
        ymmaptr = y_mma_->Values();
        lamptr = lam_->Values();
        muptr = mu_->Values();
        sptr = s_->Values();

        double* dellamptr = dellam.Values();
        double* aptr = a_->Values();
        double* bptr = b_->Values();
        double* delmuptr = delmu.Values();
        double* delsptr = dels.Values();



        for (int i = 0; i < m_; i++)
        {
          double gvec1 = 0.0;
          double gvec2 = 0.0;
          uxinv1.Dot(*(*P_)(i), &gvec1);
          xlinv1.Dot(*(*Q_)(i), &gvec2);

          *delyptr = *cptr + *dptr * (*ymmaptr) - *lamptr - *muptr;
          *dellamptr = gvec1 + gvec2 - *aptr * (z_mma_) - *ymmaptr + *sptr - *bptr;
          *delmuptr = *muptr * (*ymmaptr) - tol_sub;
          *delsptr = *lamptr * (*sptr) - tol_sub;

          delyptr++;
          cptr++;
          dptr++;
          ymmaptr++;
          lamptr++;
          muptr++;
          dellamptr++;
          aptr++;
          bptr++;
          delmuptr++;
          delsptr++;
          sptr++;
        }

        double delz = 0.0;
        double delzet = 0.0;
        double dot_a_lam = 0.0;
        dot_a_lam = a_->Dot(*lam_);
        delz = a0_ - zet_ - dot_a_lam;
        delzet = zet_ * z_mma_ - tol_sub;


        /* Aufstellen der Blöcke der Matrix:
         * Zuerst ein paar Hilfsgrößen, dann die Matrixblöcke
         */
        Epetra_Vector temp1(x_->Map(), false);

        double* temp1ptr = temp1.Values();
        diagxminaptr = diagxmina.Values();
        diagpsiinvptr = diagpsiinv.Values();
        diagxsiptr = diagxsi.Values();

        for (int i = 0; i < n_loc_; i++)
        {
          *temp1ptr = 1.0 / (*diagxminaptr + *diagpsiinvptr * (*diagxsiptr));

          temp1ptr++;
          diagxminaptr++;
          diagpsiinvptr++;
          diagxsiptr++;
        }

        Epetra_Vector help1(x_->Map(), false);
        Epetra_Vector help1inv(x_->Map(), false);

        double* help1ptr = help1.Values();
        double* help1invptr = help1inv.Values();
        diagbminxptr = diagbminx.Values();
        diagpsiinvptr = diagpsiinv.Values();
        diagetaptr = diageta.Values();
        diagxsiptr = diagxsi.Values();
        temp1ptr = temp1.Values();

        for (int i = 0; i < n_loc_; i++)
        {
          *help1ptr =
              *diagbminxptr + *diagpsiinvptr * (*diagetaptr) -
              *diagpsiinvptr * (*diagxsiptr) * (*temp1ptr) * (*diagpsiinvptr) * (*diagetaptr);
          *help1invptr = 1.0 / (*help1ptr);

          help1ptr++;
          diagbminxptr++;
          diagpsiinvptr++;
          diagetaptr++;
          diagxsiptr++;
          temp1ptr++;
          help1invptr++;
        }

        /* help2 = diagd*(diagymma + diagdinv*diagmu)*GG*diagpsiinv; % ist eine m x n
         * Matrix(Multivektor). Interpretiere diesen Term wie folgt: Zunächst (durch die
         * Multiplikation von links) werden die ZEILEN von GG mit den entsprechenden Zahlen aus der
         * linken Matrix skaliert. Danach werden die entsprechenden Zeilen SPALTEN so skaliert. In
         * Ahnlehnung an das Konzept der Multivektoren folgt also die folgende
         * Programmierung.
         * Definiere leftmat die Matrix links von GG, nur hilfsgröße
         */

        Epetra_MultiVector help2(x_->Map(), m_);

        diagdptr = diagd.Values();
        diagymmaptr = diagymma.Values();
        diagdinvptr = diagdinv.Values();
        diagmuptr = diagmu.Values();

        for (int i = 0; i < m_; i++)
        {
          Epetra_Vector gg(View, GG, i);
          Epetra_Vector help2temp(View, help2, i);

          double leftmat = 0.0;
          double* help2tempptr = help2temp.Values();
          double* ggptr = gg.Values();
          diagpsiinvptr = diagpsiinv.Values();

          leftmat = *diagdptr * (*diagymmaptr + *diagdinvptr * (*diagmuptr));

          for (int j = 0; j < n_loc_; j++)
          {
            *help2tempptr = leftmat * (*ggptr) * (*diagpsiinvptr);

            help2tempptr++;
            ggptr++;
            diagpsiinvptr++;
          }

          diagdptr++;
          diagymmaptr++;
          diagdinvptr++;
          diagmuptr++;
        }


        Epetra_MultiVector help3(x_->Map(), m_);

        for (int i = 0; i < m_; i++)
        {
          Epetra_Vector help2temp(View, help2, i);
          Epetra_Vector help3temp(View, help3, i);

          double* help3tempptr = help3temp.Values();
          double* help2tempptr = help2temp.Values();

          for (int j = 0; j < n_loc_; j++)
          {
            *help3tempptr = -*help2tempptr;

            help2tempptr++;
            help3tempptr++;
          }
        }

        /* Gehe für diesen Term analog zu "help2" vor: Also zuerst
         * Zusammenfassen aller Diagonalmatrizen, Dann Multiplikation der
         * Diagonalmatrizen mit GG'; also skalierung der SPALTEN von GG mit dem
         * jeweiligen Eintrag der Diagonalmatrix.
         * help4 = diagpsiinv*diageta*GG' - temp1*diagpsiinv*diageta*diagpsiinv*diagxsi*GG';
         */
        Epetra_MultiVector help4(x_->Map(), m_);

        for (int i = 0; i < m_; i++)
        {
          Epetra_Vector help4temp(View, help4, i);
          Epetra_Vector gg(View, GG, i);

          double* help4tempptr = help4temp.Values();
          double* diagpsiinvptr = diagpsiinv.Values();
          double* diagetaptr = diageta.Values();
          double* ggptr = gg.Values();
          double* temp1ptr = temp1.Values();
          double* diagxsiptr = diagxsi.Values();

          for (int j = 0; j < n_loc_; j++)
          {
            *help4tempptr = (*diagpsiinvptr) * (*diagetaptr) * (*ggptr) -
                            (*temp1ptr) * (*diagpsiinvptr) * (*diagetaptr) * (*diagpsiinvptr) *
                                (*diagxsiptr) * (*ggptr);

            help4tempptr++;
            diagpsiinvptr++;
            diagetaptr++;
            ggptr++;
            temp1ptr++;
            diagxsiptr++;
          }
        }

        /* Vorgehen wieder analog zu oben. Ergebnis ist eine vollbesetzte m x m Matrix
         * help5 = diagdinv*diagmu + diagd*(diagymma + diagdinv*diagmu)*(-GG*diagpsiinv*GG' -
         * diagdinv); %test ausgabe help5 = zeros(m,m); % initialisieren m x m Matrizen sind klein
         * und werden deshalb werden diese hier gelegentlich seriell ausgeführt
         */

        Epetra_SerialDenseMatrix help5(m_, m_);
        double* help5ptr = help5.A();

        diagdinvptr = diagdinv.Values();
        diagdptr = diagd.Values();
        diagymmaptr = diagymma.Values();
        diagmuptr = diagmu.Values();

        // Schleife über alle Spalten
        for (int icol = 0; icol < m_; icol++)
        {
          Epetra_Vector Gcol(View, GG, icol);

          // Schleife über alle Zeilen
          for (int irow = 0; irow < m_; irow++)
          {
            double* Gcolptr = Gcol.Values();

            Epetra_Vector Grow(View, GG, irow);
            double* Growptr = Grow.Values();

            double* diagpsiinvptr = diagpsiinv.Values();

            double tempvar = 0.0;

            // Skalarproduktschleife
            for (int i = 0; i < n_loc_; i++)
            {
              tempvar -= *Gcolptr * *Growptr * *diagpsiinvptr;

              Growptr++;
              Gcolptr++;
              diagpsiinvptr++;
            }
            discret_->Comm().SumAll(&tempvar, help5ptr, 1);

            // Fallunterscheidung zum draufaddieren von Diagonalmatrix
            if (irow == icol) *help5ptr -= *diagdinvptr;

            // Multiplikation mit Diagonalmatrix
            *help5ptr *= *diagdptr * (*diagymmaptr + *diagdinvptr * *diagmuptr);

            // Fallunterscheidung zu draufaddieren von weiterer Diagonalmatrix
            if (irow == icol) *help5ptr += *diagdinvptr * *diagmuptr;

            help5ptr++;
          }
          diagdinvptr++;
          diagdptr++;
          diagymmaptr++;
          diagmuptr++;
        }

        // Ergibt einen m x 1 Vektor
        // help6 = -diagd*(diagymma + diagdinv*diagmu)*a;
        Epetra_SerialDenseVector help6(*y_mma_);

        double* help6ptr = help6.Values();
        diagdptr = diagd.Values();
        diagymmaptr = diagymma.Values();
        diagdinvptr = diagdinv.Values();
        diagmuptr = diagmu.Values();
        aptr = a_->Values();

        for (int i = 0; i < m_; i++)
        {
          *help6ptr = -(*diagdptr) * (*diagymmaptr + (*diagdinvptr) * (*diagmuptr)) * (*aptr);

          help6ptr++;
          diagdptr++;
          diagymmaptr++;
          diagdinvptr++;
          diagmuptr++;
          aptr++;
        }

        // ergibt mxn Multivektor. Help2,3 sind mxn Multis, alles andere ist
        // diag(nxn)
        // help7 = help3 + help2*temp1*diagpsiinv*diagxsi;

        Epetra_MultiVector help7(x_->Map(), m_);
        for (int i = 0; i < m_; i++)
        {
          Epetra_Vector help2temp(View, help2, i);
          Epetra_Vector help3temp(View, help3, i);
          Epetra_Vector help7temp(View, help7, i);

          double* help2tempptr = help2temp.Values();
          double* help3tempptr = help3temp.Values();
          double* help7tempptr = help7temp.Values();

          double* temp1ptr = temp1.Values();
          double* diagpsiinvptr = diagpsiinv.Values();
          double* diagxsiptr = diagxsi.Values();

          for (int j = 0; j < n_loc_; j++)
          {
            *help7tempptr =
                *help3tempptr + (*help2tempptr) * (*temp1ptr) * (*diagpsiinvptr) * (*diagxsiptr);

            help2tempptr++;
            help3tempptr++;
            help7tempptr++;
            temp1ptr++;
            diagpsiinvptr++;
            diagxsiptr++;
          }
        }

        // ergibt mxm Matrix(Dense). help5=mxm; help2=mxn; temp1=diag(nxn);
        // help7=mxn; help4=nxm; help1inv=diag(nxn)
        // help8 = help5 + help2*temp1*diagpsiinv*diagxsi*GG' - help7*help1inv*help4;
        Epetra_SerialDenseMatrix help8(m_, m_);
        double* help8ptr = help8.A();
        help5ptr = help5.A();

        // Schleife über alle Spalten
        for (int icol = 0; icol < m_; icol++)
        {
          Epetra_Vector colvec1(View, GG, icol);
          Epetra_Vector colvec2(View, help4, icol);

          // Schleife über alle Zeilen
          for (int irow = 0; irow < m_; irow++)
          {
            double* colvec1ptr = colvec1.Values();
            double* colvec2ptr = colvec2.Values();

            Epetra_Vector rowvec1(View, help2, irow);
            double* rowvec1ptr = rowvec1.Values();

            Epetra_Vector rowvec2(View, help7, irow);
            double* rowvec2ptr = rowvec2.Values();

            double* temp1ptr = temp1.Values();
            double* diagpsiinvptr = diagpsiinv.Values();
            double* help1invptr = help1inv.Values();
            double* diagxsiptr = diagxsi.Values();

            double tempvar = 0.0;

            // Skalarproduktschleife
            for (int i = 0; i < n_loc_; i++)
            {
              tempvar += (*rowvec1ptr * *temp1ptr * *diagpsiinvptr * *diagxsiptr * *colvec1ptr) -
                         (*rowvec2ptr * *help1invptr * *colvec2ptr);

              rowvec1ptr++;
              temp1ptr++;
              diagpsiinvptr++;
              colvec1ptr++;
              rowvec2ptr++;
              help1invptr++;
              colvec2ptr++;
              diagxsiptr++;
            }
            discret_->Comm().SumAll(&tempvar, help8ptr, 1);

            // Draufaddieren von anderer mxm Matrix
            *help8ptr += *help5ptr;

            help8ptr++;
            help5ptr++;
          }
        }


        //  ergibt diag(mxm)
        // help9 = diagd*(diagymma + diagdinv*diagmu);
        Epetra_SerialDenseVector help9(*y_mma_);

        double* help9ptr = help9.Values();
        diagdptr = diagd.Values();
        diagymmaptr = diagymma.Values();
        diagdinvptr = diagdinv.Values();
        diagmuptr = diagmu.Values();

        for (int i = 0; i < m_; i++)
        {
          *help9ptr = *diagdptr * (*diagymmaptr + *diagdinvptr * *diagmuptr);
        }

        // A11 = diagpsi
        Epetra_Vector A11(x_->Map(), false);
        double* A11ptr = A11.Values();
        diagpsiptr = diagpsi.Values();

        for (int i = 0; i < n_loc_; i++)
        {
          *A11ptr = *diagpsiptr;

          A11ptr++;
          diagpsiptr++;
        }

        // A15 = -eye(n); % einfach: for... A15(i)=1; also als Vektor erstellen
        Epetra_Vector A15(x_->Map(), false);
        double* A15ptr = A15.Values();

        for (int i = 0; i < n_loc_; i++)
        {
          *A15ptr = -1.0;
          A15ptr++;
        }

        // A16 = eye(n); % - " -
        Epetra_Vector A16(x_->Map(), false);
        double* A16ptr = A16.Values();

        for (int i = 0; i < n_loc_; i++)
        {
          *A16ptr = 1.0;
          A16ptr++;
        }

        // A18 = GG'; % einfach übernehmen
        Epetra_MultiVector A18(x_->Map(), m_);

        for (int i = 0; i < m_; i++)
        {
          Epetra_Vector A18temp(View, A18, i);
          Epetra_Vector GGtemp(View, GG, i);

          double* A18tempptr = A18temp.Values();
          double* GGtempptr = GGtemp.Values();

          for (int j = 0; j < n_loc_; j++)
          {
            *A18tempptr = *GGtempptr;

            A18tempptr++;
            GGtempptr++;
          }
        }

        // A22 = diagd;
        Epetra_SerialDenseVector A22(*y_mma_);

        diagdptr = diagd.Values();
        double* A22ptr = A22.Values();

        for (int i = 0; i < m_; i++)
        {
          *A22ptr = *diagdptr;

          A22ptr++;
          diagdptr++;
        }

        // A23 = -ones(m,1) //bzw -eye(m)
        // A28 = -ones(m,1) //bzw -eye(m)
        Epetra_SerialDenseVector A23(*y_mma_);
        double* A23ptr = A23.Values();
        Epetra_SerialDenseVector A28(*y_mma_);
        double* A28ptr = A28.Values();

        for (int i = 0; i < m_; i++)
        {
          *A23ptr = -1.0;
          *A28ptr = -1.0;
          A23ptr++;
          A28ptr++;
        }

        // A33 = -diagdinv;
        Epetra_SerialDenseVector A33(*y_mma_);
        double* A33ptr = A33.Values();
        diagdinvptr = diagdinv.Values();

        for (int i = 0; i < m_; i++)
        {
          *A33ptr = -*diagdinvptr;

          A33ptr++;
          diagdinvptr++;
        }

        // A35 =  GG*diagpsiinv;
        Epetra_MultiVector A35(x_->Map(), m_);

        for (int i = 0; i < m_; i++)
        {
          Epetra_Vector ggtemp(View, GG, i);
          Epetra_Vector A35temp(View, A35, i);

          double* ggtempptr = ggtemp.Values();
          double* A35tempptr = A35temp.Values();
          double* diagpsiinvptr = diagpsiinv.Values();

          for (int j = 0; j < n_loc_; j++)
          {
            *A35tempptr = *ggtempptr * *diagpsiinvptr;

            A35tempptr++;
            ggtempptr++;
            diagpsiinvptr++;
          }
        }

        // A36 =  -GG*diagpsiinv;
        Epetra_MultiVector A36(x_->Map(), m_);

        for (int i = 0; i < m_; i++)
        {
          Epetra_Vector ggtemp(View, GG, i);
          Epetra_Vector A36temp(View, A36, i);

          double* ggtempptr = ggtemp.Values();
          double* A36tempptr = A36temp.Values();
          double* diagpsiinvptr = diagpsiinv.Values();

          for (int j = 0; j < n_loc_; j++)
          {
            *A36tempptr = -*ggtempptr * *diagpsiinvptr;

            A36tempptr++;
            ggtempptr++;
            diagpsiinvptr++;
          }
        }

        // A37 = -a;
        Epetra_SerialDenseVector A37(*y_mma_);
        double* A37ptr = A37.Values();
        aptr = a_->Values();

        for (int i = 0; i < m_; i++)
        {
          *A37ptr = -*aptr;

          A37ptr++;
          aptr++;
        }


        // A38 = -GG*diagpsiinv*GG' - diagdinv; % muss Umgeschrieben werden, code kann teilweise von
        // help8 übernommen werden, ergibt mxm dense
        Epetra_SerialDenseMatrix A38(m_, m_);
        double* A38ptr = A38.A();
        diagdinvptr = diagd.Values();

        // Schleife über alle Spalten
        for (int icol = 0; icol < m_; icol++)
        {
          // Schleife über alle Zeilen
          for (int irow = 0; irow < m_; irow++)
          {
            Epetra_Vector ggcol(View, GG, icol);
            double* ggcolptr = ggcol.Values();

            Epetra_Vector ggrow(View, GG, irow);
            double* ggrowptr = ggrow.Values();

            double* diagpsiinvptr = diagpsiinv.Values();

            double tempvar = 0.0;

            // Skalarproduktschleife
            for (int i = 0; i < n_loc_; i++)
            {
              tempvar -= *ggrowptr * *diagpsiinvptr * *ggcolptr;

              ggrowptr++;
              diagpsiinvptr++;
              ggcolptr++;
            }
            discret_->Comm().SumAll(&tempvar, A38ptr, 1);

            // Fallunterscheidung Draufaddieren von anderer mxm Matrix
            if (icol == irow) *A38ptr -= *diagdinvptr;

            A38ptr++;
          }
          diagdinvptr++;
        }


        // A39 = ones(m,1);
        Epetra_SerialDenseVector A39(*y_mma_);
        double* A39ptr = A39.Values();

        for (int i = 0; i < m_; i++)
        {
          *A39ptr = 1.0;
          A39ptr++;
        }

        // A44 = -1;
        double A44 = -1.0;

        // A48 = -a';
        Epetra_SerialDenseVector A48(*y_mma_);

        double* A48ptr = A48.Values();
        aptr = a_->Values();

        for (int i = 0; i < m_; i++)
        {
          *A48ptr = -*aptr;

          A48ptr++;
          aptr++;
        }


        // A55 = diagxmina + diagpsiinv*diagxsi; % als Vektor übernehmen, minimal umschreiben
        Epetra_Vector A55(x_->Map(), false);

        double* A55ptr = A55.Values();
        diagxminaptr = diagxmina.Values();
        diagpsiinvptr = diagpsiinv.Values();
        diagxsiptr = diagxsi.Values();

        for (int i = 0; i < n_loc_; i++)
        {
          *A55ptr = *diagxminaptr + *diagpsiinvptr * *diagxsiptr;

          A55ptr++;
          diagxminaptr++;
          diagpsiinvptr++;
          diagxsiptr++;
        }


        // A56 = -diagpsiinv*diagxsi;
        Epetra_Vector A56(x_->Map(), false);

        double* A56ptr = A56.Values();
        diagpsiinvptr = diagpsiinv.Values();
        diagxsiptr = diagxsi.Values();

        for (int i = 0; i < n_loc_; i++)
        {
          *A56ptr = -*diagpsiinvptr * *diagxsiptr;

          A56ptr++;
          diagpsiinvptr++;
          diagxsiptr++;
        }


        // A58 = -diagpsiinv*diagxsi*GG';
        Epetra_MultiVector A58(x_->Map(), m_);

        for (int i = 0; i < m_; i++)
        {
          Epetra_Vector A58temp(View, A58, i);
          Epetra_Vector ggtemp(View, GG, i);

          double* A58tempptr = A58temp.Values();
          double* diagpsiinvptr = diagpsiinv.Values();
          double* diagxsiptr = diagxsi.Values();
          double* ggtempptr = ggtemp.Values();

          for (int j = 0; j < n_loc_; j++)
          {
            *A58tempptr = -*diagpsiinvptr * *diagxsiptr * *ggtempptr;

            A58tempptr++;
            diagpsiinvptr++;
            diagxsiptr++;
            ggtempptr++;
          }
        }


        // A66 = help1; nxn(diag)
        Epetra_Vector A66(x_->Map(), false);

        help1ptr = help1.Values();
        double* A66ptr = A66.Values();

        for (int i = 0; i < n_loc_; i++)
        {
          *A66ptr = *help1ptr;

          A66ptr++;
          help1ptr++;
        }


        // A68 = help4; nxm

        Epetra_MultiVector A68(x_->Map(), m_);

        for (int i = 0; i < m_; i++)
        {
          Epetra_Vector A68temp(View, A68, i);
          Epetra_Vector help4temp(View, help4, i);

          double* help4tempptr = help4temp.Values();
          double* A68tempptr = A68temp.Values();

          for (int j = 0; j < n_loc_; j++)
          {
            *A68tempptr = *help4tempptr;

            A68tempptr++;
            help4tempptr++;
          }
        }

        // A77 = zet;
        double A77 = zet_;


        // A78 = -a'*zmma;
        Epetra_SerialDenseVector A78(*y_mma_);

        aptr = a_->Values();
        double* A78ptr = A78.Values();

        for (int i = 0; i < m_; i++)
        {
          *A78ptr = -*aptr * z_mma_;

          A78ptr++;
          aptr++;
        }


        // A79 = zeros(1,m); % Braucht man für kleines LGS
        Epetra_SerialDenseVector A79(*y_mma_);
        double* A79ptr = A79.Values();

        for (int i = 0; i < m_; i++)
        {
          *A79ptr = 0.0;
          A79ptr++;
        }


        // A87 = help6;
        Epetra_SerialDenseVector A87(*y_mma_);

        double* A87ptr = A87.Values();
        help6ptr = help6.Values();

        for (int i = 0; i < m_; i++)
        {
          *A87ptr = *help6ptr;

          A87ptr++;
          help6ptr++;
        }


        // A88 = help8;
        Epetra_SerialDenseMatrix A88(m_, m_);

        double* A88ptr = A88.A();
        help8ptr = help8.A();

        for (int icol = 0; icol < m_; icol++)
        {
          for (int irow = 0; irow < m_; irow++)
          {
            *A88ptr = *help8ptr;

            A88ptr++;
            help8ptr++;
          }
        }


        // A89 = help9;
        Epetra_SerialDenseVector A89(*y_mma_);

        double* A89ptr = A89.Values();
        help9ptr = help9.Values();

        for (int i = 0; i < m_; i++)
        {
          *A89ptr = *help9ptr;

          A89ptr++;
          help9ptr++;
        }


        // A97 = zeros(m,1); %braucht man für kleines Gleichungssystem
        Epetra_SerialDenseVector A97(*y_mma_);
        double* A97ptr = A97.Values();

        for (int i = 0; i < m_; i++)
        {
          *A97ptr = 0.0;

          A97ptr++;
        }


        // A98 = diags;
        Epetra_SerialDenseVector A98(*y_mma_);

        double* A98ptr = A98.Values();
        diagsptr = diags.Values();

        for (int i = 0; i < m_; i++)
        {
          *A98ptr = *diagsptr;

          A98ptr++;
          diagsptr++;
        }


        // A99 = diaglam;
        Epetra_SerialDenseVector A99(*y_mma_);

        double* A99ptr = A99.Values();
        diaglamptr = diaglam.Values();

        for (int i = 0; i < m_; i++)
        {
          *A99ptr = *diaglamptr;

          A99ptr++;
          diaglamptr++;
        }

        // noch ein paar inverse Matrizen...
        Epetra_Vector A11inv(x_->Map(), false);
        Epetra_Vector A55inv(x_->Map(), false);
        Epetra_Vector A66inv(x_->Map(), false);

        double* A11invptr = A11inv.Values();
        double* A55invptr = A55inv.Values();
        double* A66invptr = A66inv.Values();
        A11ptr = A11.Values();
        A55ptr = A55.Values();
        A66ptr = A66.Values();

        for (int i = 0; i < n_loc_; i++)
        {
          *A11invptr = 1.0 / (*A11ptr);
          *A55invptr = 1.0 / (*A55ptr);
          *A66invptr = 1.0 / (*A66ptr);

          A11invptr++;
          A55invptr++;
          A66invptr++;
          A11ptr++;
          A55ptr++;
          A66ptr++;
        }


        Epetra_SerialDenseVector A22inv(*y_mma_);
        Epetra_SerialDenseVector A33inv(*y_mma_);

        double* A22invptr = A22inv.Values();
        double* A33invptr = A33inv.Values();
        A22ptr = A22.Values();
        A33ptr = A33.Values();

        for (int i = 0; i < m_; i++)
        {
          *A22invptr = 1.0 / (*A22ptr);
          *A33invptr = 1.0 / (*A33ptr);

          A22invptr++;
          A33invptr++;
          A22ptr++;
          A33ptr++;
        }

        double A44inv = 1.0 / A44;


        // Hilfsgrößen zum Aufstellen der rechten Seite///////////////////////////

        // rhs1 = temp1*diagpsiinv*diageta*(-delxsi + diagpsiinv*diagxsi*delx); % muss geändert
        // werden. ergibt nx1
        Epetra_Vector rhs1(x_->Map());

        double* rhs1ptr = rhs1.Values();
        temp1ptr = temp1.Values();
        diagpsiinvptr = diagpsiinv.Values();
        diagetaptr = diageta.Values();
        delxsiptr = delxsi.Values();
        diagxsiptr = diagxsi.Values();
        delxptr = delx.Values();

        for (int i = 0; i < n_loc_; i++)
        {
          *rhs1ptr = (*temp1ptr) * (*diagpsiinvptr) * (*diagetaptr) *
                     (-*delxsiptr + (*diagpsiinvptr) * (*diagxsiptr) * (*delxptr));

          rhs1ptr++;
          temp1ptr++;
          diagpsiinvptr++;
          diagetaptr++;
          delxsiptr++;
          diagxsiptr++;
          delxptr++;
        }


        // rhs2 = diagd*(diagymma + diagdinv*diagmu)*(-dellam + GG*diagpsiinv*delx - diagdinv*dely);
        // % muss geändert werden, ergibt mx1
        Epetra_SerialDenseVector rhs2(*y_mma_);

        double* rhs2ptr = rhs2.Values();
        diagdptr = diagd.Values();
        diagymmaptr = diagymma.Values();
        diagdinvptr = diagdinv.Values();
        diagmuptr = diagmu.Values();
        dellamptr = dellam.Values();
        delyptr = dely.Values();

        for (int i = 0; i < m_; i++)
        {
          Epetra_Vector ggtemp(View, GG, i);
          double tempvar = 0.0;

          double* ggtempptr = ggtemp.Values();
          double* diagpsiinvptr = diagpsiinv.Values();
          double* delxptr = delx.Values();

          for (int j = 0; j < n_loc_; j++)
          {
            tempvar += *ggtempptr * *diagpsiinvptr * *delxptr;

            ggtempptr++;
            diagpsiinvptr++;
            delxptr++;
          }
          discret_->Comm().SumAll(&tempvar, rhs2ptr, 1);

          *rhs2ptr = *diagdptr * (*diagymmaptr + *diagdinvptr * *diagmuptr) *
                     (-*dellamptr + *rhs2ptr - *diagdinvptr * *delyptr);

          rhs2ptr++;
          diagdptr++;
          diagymmaptr++;
          diagdinvptr++;
          diagmuptr++;
          dellamptr++;
          delyptr++;
        }



        // b1 = -delx;
        Epetra_Vector b1(x_->Map(), false);

        double* b1ptr = b1.Values();
        delxptr = delx.Values();

        for (int i = 0; i < n_loc_; i++)
        {
          *b1ptr = -*delxptr;

          b1ptr++;
          delxptr++;
        }


        // b2 = -dely;
        Epetra_SerialDenseVector b2(*y_mma_);

        double* b2ptr = b2.Values();
        delyptr = dely.Values();

        for (int i = 0; i < m_; i++)
        {
          *b2ptr = -*delyptr;

          b2ptr++;
          delyptr++;
        }


        // b3 = -dellam+GG*diagpsiinv*delx - diagdinv*dely;
        Epetra_SerialDenseVector b3(*y_mma_);

        double* b3ptr = b3.Values();
        dellamptr = dellam.Values();
        diagdinvptr = diagd.Values();
        delyptr = dely.Values();

        for (int i = 0; i < m_; i++)
        {
          Epetra_Vector ggtemp(View, GG, i);
          double* ggtempptr = ggtemp.Values();
          double* diagpsiinvptr = diagpsiinv.Values();
          double* delxptr = delx.Values();

          double tempvar = 0.0;

          for (int j = 0; j < n_loc_; j++)
          {
            tempvar += *ggtempptr * *diagpsiinvptr * *delxptr;

            ggtempptr++;
            diagpsiinvptr++;
            delxptr++;
          }
          discret_->Comm().SumAll(&tempvar, b3ptr, 1);

          *b3ptr = -*dellamptr + *b3ptr - *diagdinvptr * *delyptr;

          b3ptr++;
          dellamptr++;
          diagdinvptr++;
          delyptr++;
        }


        // b4 = -delz;
        double b4 = -delz;


        // b5 = -delxsi + diagpsiinv*diagxsi*delx;
        Epetra_Vector b5(x_->Map(), false);

        double* b5ptr = b5.Values();
        delxsiptr = delxsi.Values();
        diagpsiinvptr = diagpsiinv.Values();
        diagxsiptr = diagxsi.Values();
        delxptr = delx.Values();

        for (int i = 0; i < n_loc_; i++)
        {
          *b5ptr = -(*delxsiptr) + (*diagpsiinvptr) * (*diagxsiptr) * (*delxptr);

          b5ptr++;
          delxsiptr++;
          diagpsiinvptr++;
          diagxsiptr++;
          delxptr++;
        }


        // b6 = -deleta - diagpsiinv*diageta*delx + rhs1;
        Epetra_Vector b6(x_->Map(), false);

        double* b6ptr = b6.Values();
        deletaptr = deleta.Values();
        diagpsiinvptr = diagpsiinv.Values();
        diagetaptr = diageta.Values();
        delxptr = delx.Values();
        rhs1ptr = rhs1.Values();

        for (int i = 0; i < n_loc_; i++)
        {
          *b6ptr = -(*deletaptr) - (*diagpsiinvptr) * (*diagetaptr) * (*delxptr) + *rhs1ptr;

          b6ptr++;
          deletaptr++;
          diagpsiinvptr++;
          diagetaptr++;
          delxptr++;
          rhs1ptr++;
        }


        // b7 = -delzet - zmma*delz;
        double b7 = -delzet - (z_mma_ * delz);


        // b8 = -delmu + diagdinv*diagmu*dely + rhs2 - (help7*help1inv*b6 + help2*temp1*b5); ergibt
        // m x 1
        Epetra_SerialDenseVector b8(*y_mma_);

        double* b8ptr = b8.Values();
        delmuptr = delmu.Values();
        diagdinvptr = diagdinv.Values();
        diagmuptr = diagmu.Values();
        delyptr = dely.Values();
        rhs2ptr = rhs2.Values();

        for (int i = 0; i < m_; i++)
        {
          Epetra_Vector help2temp(View, help2, i);
          Epetra_Vector help7temp(View, help7, i);

          double* help2tempptr = help2temp.Values();
          double* help7tempptr = help7temp.Values();
          double* help1invptr = help1inv.Values();
          double* b6ptr = b6.Values();
          double* temp1ptr = temp1.Values();
          double* b5ptr = b5.Values();


          double tempvar = 0.0;

          for (int j = 0; j < n_loc_; j++)
          {
            tempvar += (*help7tempptr) * (*help1invptr) * (*b6ptr) +
                       (*help2tempptr) * (*temp1ptr) * (*b5ptr);

            help2tempptr++;
            help7tempptr++;
            help1invptr++;
            b6ptr++;
            temp1ptr++;
            b5ptr++;
          }
          discret_->Comm().SumAll(&tempvar, b8ptr, 1);

          *b8ptr = -(*delmuptr) + (*diagdinvptr) * (*diagmuptr) * (*delyptr) + *rhs2ptr - *b8ptr;

          b8ptr++;
          delmuptr++;
          diagdinvptr++;
          diagmuptr++;
          delyptr++;
          rhs2ptr++;
        }


        // b9 = -dels;
        Epetra_SerialDenseVector b9(*y_mma_);

        double* b9ptr = b9.Values();
        delsptr = dels.Values();

        for (int i = 0; i < m_; i++)
        {
          *b9ptr = -(*delsptr);

          b9ptr++;
          delsptr++;
        }


        // AUFSTELLEN DER KLEINEN MATRIX UND RECHTER SEITE, DANN SOLVER.
        // HINWEIS: DIE VORKONDITIONIERUNG IST ESSENTIELL, ZUMINDEST IN MATLAB
        Epetra_SerialDenseMatrix AA(2 * m_ + 1, 2 * m_ + 1);
        double* AAptr = AA.A();

        for (int irow = 0; irow < (2 * m_ + 1); irow++)
        {
          for (int icol = 0; icol < (2 * m_ + 1); icol++)
          {
            *AAptr = 0.0;

            AAptr++;
          }
        }


        // A77 (1 x 1) double
        AA(0, 0) = A77;

        // A78 (1 x m) SerialDenseVector
        A78ptr = A78.Values();

        for (int i = 1; i < m_ + 1; i++)
        {
          AA(0, i) = *A78ptr;

          A78ptr++;
        }

        // A79 (1 x m) SerialDenseVector
        A79ptr = A79.Values();

        for (int i = m_ + 1; i < 2 * m_ + 1; i++)
        {
          AA(0, i) = *A79ptr;

          A79ptr++;
        }

        // A87 (m x 1) SerialDenseVector)
        A87ptr = A87.Values();

        for (int i = 1; i < m_ + 1; i++)
        {
          AA(i, 0) = *A87ptr;

          A87ptr++;
        }

        // A88 (m x m) SerialDenseMatrix
        for (int irow = 0; irow < m_; irow++)
        {
          for (int icol = 0; icol < m_; icol++)
          {
            AA((1 + irow), (1 + icol)) = A88(irow, icol);
          }
        }

        // A89 (m x m) SerialDenseVector, weil diagonal
        A89ptr = A89.Values();

        for (int i = 0; i < m_; i++)
        {
          AA(1 + i, m_ + 1 + i) = *A89ptr;

          A89ptr++;
        }

        // A97 (m x 1) SerialDenseVector
        A97ptr = A97.Values();

        for (int i = m_ + 1; i < 2 * m_ + 1; i++)
        {
          AA(i, 0) = *A97ptr;

          A97ptr++;
        }

        // A98 (m x m) SerialDenseVector, weil diagonal
        A98ptr = A98.Values();

        for (int i = 0; i < m_; i++)
        {
          AA(m_ + 1 + i, 1 + i) = *A98ptr;

          A98ptr++;
        }

        // A98 (m x m) SerialDenseVector, weil diagonal
        A99ptr = A99.Values();

        for (int i = 0; i < m_; i++)
        {
          AA(m_ + 1 + i, m_ + 1 + i) = *A99ptr;

          A99ptr++;
        }

        // Aufstellen rechte Seite
        Epetra_SerialDenseVector bb(2 * m_ + 1);

        // b7 (1 x 1)
        bb(0) = b7;

        // b8 (m x 1)
        b8ptr = b8.Values();

        for (int i = 0; i < m_; i++)
        {
          bb(i + 1) = *b8ptr;

          b8ptr++;
        }

        // b9 (m x 1)
        b9ptr = b9.Values();

        for (int i = 0; i < m_; i++)
        {
          bb(i + m_ + 1) = *b9ptr;

          b9ptr++;
        }


        // VORKONDITIONIERUNG (SEHR WICHTIG !!!) (es werden einfach alle Zeilen mit dem
        // Betragsmaximalen Element skaliert)
        for (int irow = 0; irow < 2 * m_ + 1; irow++)
        {
          // Berechne Betragsmaximales Element
          double maxeintrag = 0.0;
          for (int icol = 0; icol < 2 * m_ + 1; icol++)
          {
            if (std::abs(AA(irow, icol)) > maxeintrag) maxeintrag = AA(irow, icol);
          }

          // Skaliere ganze Zeile mit "maxeintrag", WICHTIG: AUCH RECHTE SEITE WIRD SKALIERT
          for (int icol = 0; icol < 2 * m_ + 1; icol++)
            AA(irow, icol) = AA(irow, icol) / maxeintrag;
          bb(irow) = bb(irow) / maxeintrag;
        }

        Epetra_SerialDenseVector solut(2 * m_ + 1);

        Epetra_SerialDenseSolver solver;
        solver.SetMatrix(AA);
        solver.SetVectors(solut, bb);
        solver.Solve();

        // Neudefinition der Lösungskomponenten
        double sol7 = solut(0);

        Epetra_SerialDenseVector sol8(*y_mma_);
        Epetra_SerialDenseVector sol9(*y_mma_);

        for (int i = 0; i < m_; i++)
        {
          sol8(i) = solut(i + 1);
          sol9(i) = solut(m_ + 1 + i);
        }


        // AB HIER RÜCKWÄRTSSUBSTITUTION

        // sol6 = A66inv*(b6 - A68*sol8);

        // Zwischenergebnis: sol68*sol8
        Epetra_Vector A68sol8(x_->Map(), true);

        double* sol8ptr = sol8.Values();

        for (int i = 0; i < m_; i++)
        {
          Epetra_Vector A68temp(View, A68, i);

          double* A68tempptr = A68temp.Values();
          double* A68sol8ptr = A68sol8.Values();

          for (int j = 0; j < n_loc_; j++)
          {
            *A68sol8ptr += *A68tempptr * *sol8ptr;

            A68sol8ptr++;
            A68tempptr++;
          }
          sol8ptr++;
        }

        Epetra_Vector sol6(x_->Map(), false);

        double* sol6ptr = sol6.Values();
        double* A68sol8ptr = A68sol8.Values();
        b6ptr = b6.Values();
        A66invptr = A66inv.Values();

        for (int i = 0; i < n_loc_; i++)
        {
          *sol6ptr = (*A66invptr) * (*b6ptr - *A68sol8ptr);

          sol6ptr++;
          A68sol8ptr++;
          b6ptr++;
          A66invptr++;
        }


        // sol5 = A55inv*(b5 - A56*sol6 - A58*sol8);

        // Zwischenergebnis: sol58*sol8
        Epetra_Vector A58sol8(x_->Map(), true);
        sol8ptr = sol8.Values();

        for (int i = 0; i < m_; i++)
        {
          Epetra_Vector A58temp(View, A58, i);

          double* A58tempptr = A58temp.Values();
          double* A58sol8ptr = A58sol8.Values();

          for (int j = 0; j < n_loc_; j++)
          {
            *A58sol8ptr += *A58tempptr * *sol8ptr;

            A58sol8ptr++;
            A58tempptr++;
          }
          sol8ptr++;
        }

        Epetra_Vector sol5(x_->Map(), false);

        double* sol5ptr = sol5.Values();
        double* A58sol8ptr = A58sol8.Values();
        b5ptr = b5.Values();
        A55invptr = A55inv.Values();
        A56ptr = A56.Values();
        sol6ptr = sol6.Values();

        for (int i = 0; i < n_loc_; i++)
        {
          *sol5ptr = (*A55invptr) * (*b5ptr - (*A56ptr) * (*sol6ptr) - *A58sol8ptr);

          sol5ptr++;
          A58sol8ptr++;
          b5ptr++;
          A55invptr++;
          A56ptr++;
          sol6ptr++;
        }

        // sol4 = A44inv*(b4 - A48*sol8);

        // Hilfsgröße: A48*sol8;
        A48ptr = A48.Values();
        sol8ptr = sol8.Values();
        double A48sol8 = 0.0;

        for (int i = 0; i < m_; i++)
        {
          A48sol8 += (*A48ptr) * (*sol8ptr);

          A48ptr++;
          sol8ptr++;
        }

        double sol4 = A44inv * (b4 - A48sol8);


        // sol3 = A33inv*(b3 - A35*sol5 - A36*sol6 - A37*sol7 - A38*sol8 - A39*sol9);

        // Zwischenergebnis: A35*sol5, A36*sol6
        Epetra_Vector A35sol5(x_->Map(), true);
        Epetra_Vector A36sol6(x_->Map(), true);

        double* A35sol5ptr = A35sol5.Values();
        double* A36sol6ptr = A36sol6.Values();

        for (int i = 0; i < m_; i++)
        {
          Epetra_Vector A35temp(View, A35, i);
          Epetra_Vector A36temp(View, A36, i);

          double* A35tempptr = A35temp.Values();
          double* A36tempptr = A36temp.Values();
          double* sol5ptr = sol5.Values();
          double* sol6ptr = sol6.Values();

          double tempvar35 = 0.0;
          double tempvar36 = 0.0;

          for (int j = 0; j < n_loc_; j++)
          {
            tempvar35 += *A35tempptr * *sol5ptr;
            tempvar36 += *A36tempptr * *sol6ptr;

            A35tempptr++;
            A36tempptr++;
            sol5ptr++;
            sol6ptr++;
          }
          discret_->Comm().SumAll(&tempvar35, A35sol5ptr, 1);
          discret_->Comm().SumAll(&tempvar36, A36sol6ptr, 1);

          A35sol5ptr++;
          A36sol6ptr++;
        }

        // Zwischenergebnis: A38*sol8;
        Epetra_SerialDenseVector A38sol8(m_);

        for (int irow = 0; irow < m_; irow++)
        {
          for (int icol = 0; icol < m_; icol++)
          {
            A38sol8(irow) += A38(irow, icol) * sol8(icol);
          }
        }

        Epetra_Vector sol3(x_->Map(), false);

        double* sol3ptr = sol3.Values();
        A33invptr = A33inv.Values();
        b3ptr = b3.Values();
        A35sol5ptr = A35sol5.Values();
        A36sol6ptr = A36sol6.Values();
        A37ptr = A37.Values();
        double* A38sol8ptr = A38sol8.Values();
        A39ptr = A39.Values();
        double* sol9ptr = sol9.Values();

        for (int i = 0; i < m_; i++)
        {
          *sol3ptr = (*A33invptr) * (*b3ptr - *A35sol5ptr - *A36sol6ptr - (*A37ptr) * (sol7) -
                                        *A38sol8ptr - (*A39ptr) * (*sol9ptr));

          sol3ptr++;
          A33invptr++;
          b3ptr++;
          A35sol5ptr++;
          A36sol6ptr++;
          A37ptr++;
          A38sol8ptr++;
          A39ptr++;
          sol9ptr++;
        }


        // sol2 = A22inv*(b2 - A23*sol3 - A28*sol8);
        Epetra_SerialDenseVector sol2(m_);

        double* sol2ptr = sol2.Values();
        A22invptr = A22inv.Values();
        b2ptr = b2.Values();
        A23ptr = A23.Values();
        sol3ptr = sol3.Values();
        A28ptr = A28.Values();
        sol8ptr = sol8.Values();

        for (int i = 0; i < m_; i++)
        {
          *sol2ptr = (*A22invptr) * (*b2ptr - (*A23ptr) * (*sol3ptr) - (*A28ptr) * (*sol8ptr));

          sol2ptr++;
          A22ptr++;
          b2ptr++;
          A23ptr++;
          sol3ptr++;
          A28ptr++;
          sol8ptr++;
        }


        // sol1 = A11inv*(b1 - A15*sol5 - A16*sol6 - A18*sol8);

        // Zwischenergebnis: sol18*sol8
        Epetra_Vector A18sol8(x_->Map(), true);

        sol8ptr = sol8.Values();

        for (int i = 0; i < m_; i++)
        {
          Epetra_Vector A18temp(View, A18, i);

          double* A18tempptr = A18temp.Values();
          double* A18sol8ptr = A18sol8.Values();

          for (int j = 0; j < n_loc_; j++)
          {
            *A18sol8ptr += *A18tempptr * *sol8ptr;

            A18sol8ptr++;
            A18tempptr++;
          }
          sol8ptr++;
        }

        Epetra_Vector sol1(x_->Map(), false);

        double* sol1ptr = sol1.Values();
        A11invptr = A11inv.Values();
        b1ptr = b1.Values();
        A15ptr = A15.Values();
        sol5ptr = sol5.Values();
        A16ptr = A16.Values();
        sol6ptr = sol6.Values();
        double* A18sol8ptr = A18sol8.Values();

        for (int i = 0; i < n_loc_; i++)
        {
          *sol1ptr = (*A11invptr) *
                     (*b1ptr - (*A15ptr) * (*sol5ptr) - (*A16ptr) * (*sol6ptr) - *A18sol8ptr);

          sol1ptr++;
          A11invptr++;
          b1ptr++;
          A15ptr++;
          sol5ptr++;
          A16ptr++;
          sol6ptr++;
          A18sol8ptr++;
        }


        // Lösungsvektoren noch umbenennen

        Epetra_Vector dx(x_->Map(), false);
        Epetra_Vector dxsi(x_->Map(), false);
        Epetra_Vector deta(x_->Map(), false);

        double* dxptr = dx.Values();
        double* dxsiptr = dxsi.Values();
        double* detaptr = deta.Values();

        sol1ptr = sol1.Values();
        sol5ptr = sol5.Values();
        sol6ptr = sol6.Values();

        for (int i = 0; i < n_loc_; i++)
        {
          *dxptr = *sol1ptr;
          *dxsiptr = *sol5ptr;
          *detaptr = *sol6ptr;

          dxptr++;
          dxsiptr++;
          detaptr++;
          sol1ptr++;
          sol5ptr++;
          sol6ptr++;
        }

        Epetra_SerialDenseVector dy(*y_mma_);
        Epetra_SerialDenseVector dlam(*y_mma_);
        Epetra_SerialDenseVector dmu(*y_mma_);
        Epetra_SerialDenseVector ds(*y_mma_);

        double* dyptr = dy.Values();
        double* dlamptr = dlam.Values();
        double* dmuptr = dmu.Values();
        double* dsptr = ds.Values();

        sol2ptr = sol2.Values();
        sol8ptr = sol8.Values();
        sol3ptr = sol3.Values();
        sol9ptr = sol9.Values();

        for (int i = 0; i < m_; i++)
        {
          *dyptr = *sol2ptr;
          *dlamptr = *sol8ptr;
          *dmuptr = *sol3ptr;
          *dsptr = *sol9ptr;

          dyptr++;
          dlamptr++;
          dmuptr++;
          dsptr++;
          sol2ptr++;
          sol8ptr++;
          sol3ptr++;
          sol9ptr++;
        }

        double dz = sol7;
        double dzet = sol4;

        // HIER ENDE NEUES VERFAHREN

        double val = std::min(dz / z_mma_, dzet / zet_);

        y = y_mma_->Values();
        dyptr = dy.Values();
        lam = lam_->Values();
        dlamptr = dlam.Values();
        mu = mu_->Values();
        dmuptr = dmu.Values();
        s = s_->Values();
        dsptr = ds.Values();


        for (int i = 0; i < m_; i++)
        {
          val = std::min(val, *dyptr / (*y));
          val = std::min(val, *dlamptr / (*lam));
          val = std::min(val, *dmuptr / (*mu));
          val = std::min(val, *dsptr / (*s));

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

        for (int i = 0; i < n_loc_; i++)
        {
          val = std::min(val, *dxsiptr / (*xsi));
          val = std::min(val, *detaptr / (*eta));
          val = std::min(val, *dxptr / (*x - *alpha));
          val = std::min(val, *dxptr / (*x - *beta));

          xsi++;
          dxsiptr++;
          eta++;
          detaptr++;
          dxptr++;
          x++;
          alpha++;
          beta++;
        }

        // min becomes max since fac<0
        val = std::max(1.0, fac_stepsize_ * val);
        double stepsize = 0.0;
        discret_->Comm().MaxAll(&val, &stepsize, 1);
        stepsize = 1.0 / stepsize;
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
        double resnew = 2 * resnorm;


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


        while ((resnorm < resnew) and (it < max_inner_sub_iter_))
        {
          it++;

          x_mma_->Update(stepsize, dx, 1.0, x_mma_old, 0.0);
          xsi_->Update(stepsize, dxsi, 1.0, xsi_old, 0.0);
          eta_->Update(stepsize, deta, 1.0, eta_old, 0.0);
          // std::cout << "eta is " << *eta_ << std::endl;

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

          for (int i = 0; i < m_; i++)
          {
            *y = *yptr + stepsize * *dyptr;
            *lam = *lamptr + stepsize * *dlamptr;
            *mu = *muptr + stepsize * *dmuptr;
            *s = *sptr + stepsize * *dsptr;

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

          z_mma_ = z_mma_old + stepsize * dz;
          zet_ = zet_old + stepsize * dzet;


          // compute residuals
          ResApp(resnew, resinf, tol_sub);

          stepsize = stepsize / 2.0;

          if ((it == max_inner_sub_iter_) and (discret_->Comm().MyPID() == 0))
            printf(
                "Reached maximal number of iterations in most inner loop of primal dual interior "
                "point optimization algorithm\n");
        }

        resnorm = resnew;

        // it would be sufficient to compute resinf only once here and not in
        // every iteration as done above. but the effort is neglegible
        if (resinf < 0.9 * tol_sub) break;
      }

      if ((inner_iter == max_sub_iter_) and (discret_->Comm().MyPID() == 0))
        printf(
            "Reached maximal number of iterations in inner loop of primal dual interior point "
            "optimization algorithm\n");

      if (tol_sub > 1.001 * tol_sub_)
        tol_sub *= 0.1;
      else
        tol_reached = true;
    }
  }
  //  std::cout << "after subsolv:" << std::endl;
  //  std::cout << "x is " << *x_mma_ << std::endl;
  //  std::cout << "y is " << *y_mma_ << ", z is " << z_mma_ << ", lam is " << *lam_ << std::endl;
  //  std::cout << "xsi is " << *xsi_ << std::endl;
  //  std::cout << "eta is " << *eta_ << std::endl;
  //  std::cout << "mu is " << *mu_ << ", zet is " << zet_ << ", s is " << *s_ << std::endl;

  //  const std::string outname(DRT::Problem::Instance()->OutputControlFile()->FileName());
  //
  //  std::ostringstream filename1;
  //  std::ostringstream filename2;
  //  std::ostringstream filename3;
  //  std::ostringstream filename4;
  //  std::ostringstream filename5;
  //  std::ostringstream filename6;
  //  std::ostringstream filename7;
  //  std::ostringstream filename8;
  //  std::ostringstream filename9;
  //
  //  filename1 << outname << "_" << total_iter_ << "_x.mtl";
  //  filename2 << outname << "_" << total_iter_ << "_xi.mtl";
  //  filename3 << outname << "_" << total_iter_ << "_eta.mtl";
  //  filename4 << outname << "_" << total_iter_ << "_y.mtl";
  //  filename5 << outname << "_" << total_iter_ << "_lam.mtl";
  //  filename6 << outname << "_" << total_iter_ << "_mu.mtl";
  //  filename7 << outname << "_" << total_iter_ << "_z.mtl";
  //  filename8 << outname << "_" << total_iter_ << "_zeta.mtl";
  //  filename9 << outname << "_" << total_iter_ << "_s.mtl";
  //
  //  LINALG::PrintVectorInMatlabFormat(filename1.str(),*x_mma_);
  //  LINALG::PrintVectorInMatlabFormat(filename2.str(),*xsi_);
  //  LINALG::PrintVectorInMatlabFormat(filename3.str(),*eta_);
  //  LINALG::PrintSerialDenseMatrixInMatlabFormat(filename4.str(),*y_mma_);
  //  LINALG::PrintSerialDenseMatrixInMatlabFormat(filename5.str(),*lam_);
  //  LINALG::PrintSerialDenseMatrixInMatlabFormat(filename6.str(),*mu_);
  //  LINALG::PrintSerialDenseMatrixInMatlabFormat(filename7.str(),*s_);
  //
  //  Epetra_SerialDenseVector tmp_z(Copy, &z_mma_, 1);
  //  Epetra_SerialDenseVector tmp_zeta(Copy, &zet_, 1);
  //  LINALG::PrintSerialDenseMatrixInMatlabFormat(filename8.str(),tmp_z);
  //  LINALG::PrintSerialDenseMatrixInMatlabFormat(filename9.str(),tmp_zeta);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void OPTI::GCMMA::ResApp(double& resnorm, double& resinf, const double& tol_sub)
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

  for (int j = 0; j < m_; j++)
  {
    Epetra_Vector P(View, *P_, j);
    Epetra_Vector Q(View, *Q_, j);

    plam.Update(*lam, P, 1.0);
    qlam.Update(*lam, Q, 1.0);

    lam++;
  }


  /*
   * resx = plam/(ux1*ux1) - qlam/(xl1*xl1) - xsi + eta;
   */
  Epetra_Vector resX(x_->Map());
  double* resxptr = resX.Values();

  Epetra_Vector uxinv1(x_->Map(), false);
  Epetra_Vector xlinv1(x_->Map(), false);
  double* uxinv1ptr = uxinv1.Values();
  double* xlinv1ptr = xlinv1.Values();

  double ux1 = 0.0;  // one entry of Epetra_Vector, just used locally
  double xl1 = 0.0;  // one entry of Epetra_Vector, just used locally
  double* x = x_mma_->Values();
  double* asymp_max = asymp_max_->Values();
  double* asymp_min = asymp_min_->Values();
  double* plamptr = plam.Values();
  double* qlamptr = qlam.Values();
  double* xsi = xsi_->Values();
  double* eta = eta_->Values();

  for (int i = 0; i < n_loc_; i++)
  {
    ux1 = *asymp_max - *x;
    xl1 = *x - *asymp_min;

    *uxinv1ptr = 1.0 / ux1;
    *xlinv1ptr = 1.0 / xl1;

    double dpsidx = *plamptr / (ux1 * ux1) - *qlamptr / (xl1 * xl1);
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

  for (int j = 0; j < m_; j++)
  {
    *resyptr = *c + *d * *y - *mu - *lam;

    resz -= *a * *lam;

    uxinv1.Dot(*(*P_)(j), &gvec1);
    xlinv1.Dot(*(*Q_)(j), &gvec2);
    *reslamptr = gvec1 + gvec2 - *a * z_mma_ - *y + *s - *b;

    *resmuptr = *mu * *y - tol_sub;

    *resptr = *lam * *s - tol_sub;

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

  for (int i = 0; i < n_loc_; i++)
  {
    *resxsiptr = *xsi * (*x - *alpha) - tol_sub;
    *resetaptr = *eta * (*beta - *x) - tol_sub;

    resxsiptr++;
    resetaptr++;
    xsi++;
    eta++;
    x++;
    alpha++;
    beta++;
  }


  double reszet = zet_ * z_mma_ - tol_sub;


  resnorm = Res2Norm(&resX, &resXsi, &resEta, &resy, &resmu, &reslam, &res, &resz, &reszet);
  resinf = ResInfNorm(&resX, &resXsi, &resEta, &resy, &resmu, &reslam, &res, &resz, &reszet);
}


double OPTI::GCMMA::Res2Norm(Epetra_Vector* res1, Epetra_Vector* res2, Epetra_Vector* res3,
    Epetra_SerialDenseVector* res4, Epetra_SerialDenseVector* res5, Epetra_SerialDenseVector* res6,
    Epetra_SerialDenseVector* res7, double* res8, double* res9)
{
  double resnorm = *res8 * *res8 + *res9 * *res9;
  resnorm += std::pow(res4->Norm2(), 2) + pow(res5->Norm2(), 2) + pow(res6->Norm2(), 2) +
             pow(res7->Norm2(), 2);

  double locnorm = 0.0;
  res1->Norm2(&locnorm);
  resnorm += locnorm * locnorm;
  res2->Norm2(&locnorm);
  resnorm += locnorm * locnorm;
  res3->Norm2(&locnorm);
  resnorm += locnorm * locnorm;

  return sqrt(resnorm);
}


double OPTI::GCMMA::ResInfNorm(Epetra_Vector* res1, Epetra_Vector* res2, Epetra_Vector* res3,
    Epetra_SerialDenseVector* res4, Epetra_SerialDenseVector* res5, Epetra_SerialDenseVector* res6,
    Epetra_SerialDenseVector* res7, double* res8, double* res9)
{
  double locnorm = 0.0;
  double resinf = std::max(abs(*res8), abs(*res9));

  resinf = std::max(resinf, res4->NormInf());
  resinf = std::max(resinf, res5->NormInf());
  resinf = std::max(resinf, res6->NormInf());
  resinf = std::max(resinf, res7->NormInf());

  res1->NormInf(&locnorm);
  resinf = std::max(locnorm, resinf);
  res2->NormInf(&locnorm);
  resinf = std::max(locnorm, resinf);
  res3->NormInf(&locnorm);
  resinf = std::max(locnorm, resinf);

  return resinf;
}


void OPTI::GCMMA::Update()
{
  // helper vectors
  Epetra_Vector uxinv1(x_->Map());
  Epetra_Vector xlinv1(x_->Map());

  double* asymp_max = asymp_max_->Values();
  double* asymp_min = asymp_min_->Values();
  double* x = x_mma_->Values();

  double* uxinv1ptr = uxinv1.Values();
  double* xlinv1ptr = xlinv1.Values();


  for (int i = 0; i < n_loc_; i++)
  {
    *uxinv1ptr = 1.0 / (*asymp_max - *x);
    *xlinv1ptr = 1.0 / (*x - *asymp_min);

    asymp_max++;
    asymp_min++;
    x++;
    uxinv1ptr++;
    xlinv1ptr++;
  }

  // set new approximation of objective value
  double value = 0.0;

  uxinv1.Dot(*p0_, &value);
  obj_appr_ = r0_ + value;

  xlinv1.Dot(*q0_, &value);
  obj_appr_ += value;


  /*
   * set new approximation of constraints
   *
   * constr_app = b + P^T/ux + Q^T/xl
   */
  *constr_appr_ = *b_;
  constr_appr_->Scale(-1.0);

  double* constr = constr_appr_->Values();
  for (int i = 0; i < m_; i++)
  {
    uxinv1.Dot(*(*P_)(i), &value);
    *constr += value;

    xlinv1.Dot(*(*Q_)(i), &value);
    *constr += value;

    constr++;
  }

  //  std::cout << "obj_appr is " << obj_appr_ << std::endl;
  //  std::cout << "constr_appr is " << *constr_appr_ << std::endl;
}



void OPTI::GCMMA::Output()
{
  // step number and time
  output_->NewStep(total_iter_, (double)total_iter_);

  // for result test distinction between element and node based vector required
  // for restart and standard output not, so only current objective variable distincted
  if (dens_type_ == INPAR::TOPOPT::dens_node_based)
    output_->WriteVector("x_mma_n", x_mma_);
  else if (dens_type_ == INPAR::TOPOPT::dens_ele_based)
    output_->WriteVector("x_mma_e", x_mma_);
  else
    dserror("undefined type of optimization field");

  if ((DRT::INPUT::IntegralValue<bool>(DRT::Problem::Instance()->IOParams(), "OUTPUT_GMSH") ==
          true) and
      (DRT::INPUT::IntegralValue<bool>(params_, "GMSH_OUTPUT") == true) and
      (total_iter_ % upres_ == 0))
    OutputToGmsh();  // TODO look at this: total_iter_,false);

  // write domain decomposition for visualization (only once!)
  if (total_iter_ == upres_) output_->WriteElementData(true);

  output_->WriteVector("x", x_);
  output_->WriteVector("x_old", x_old_);
  output_->WriteVector("x_old2", x_old2_);

  output_->WriteDouble("obj", obj_);
  output_->WriteVector("obj_deriv", obj_deriv_);

  // optimization algorithm data - part 2 (Epetra_SerialDenseVector)
  Teuchos::RCP<std::vector<double>> constrvec = Teuchos::rcp(new std::vector<double>(m_));
  double* constr = constr_->Values();
  std::vector<double>::iterator constrvecit = constrvec->begin();
  for (int i = 0; i < m_; i++)
  {
    *constrvecit = *constr;
    constr++;
    constrvecit++;
  }
  output_->WriteRedundantDoubleVector("constr", constrvec);
  output_->WriteVector("constr_deriv", constr_deriv_);

  if (uprestart_ != 0 && total_iter_ > 0 && total_iter_ % uprestart_ == 0)  // add restart data
  {
    // iteration counter
    output_->WriteInt("out_it", outer_iter_);
    output_->WriteInt("in_it", inner_iter_);

    // optimization algorithm data - part 1 (Epetra_Vector)
    output_->WriteVector("xsi", xsi_);
    output_->WriteVector("eta", eta_);
    output_->WriteVector("asymp_max", asymp_max_);
    output_->WriteVector("asymp_min", asymp_min_);

    // optimization algorithm data - part 2 (Epetra_SerialDenseVector)
    Teuchos::RCP<std::vector<double>> yvec = Teuchos::rcp(new std::vector<double>(m_));
    Teuchos::RCP<std::vector<double>> lamvec = Teuchos::rcp(new std::vector<double>(m_));
    Teuchos::RCP<std::vector<double>> muvec = Teuchos::rcp(new std::vector<double>(m_));
    Teuchos::RCP<std::vector<double>> svec = Teuchos::rcp(new std::vector<double>(m_));
    Teuchos::RCP<std::vector<double>> rhovec = Teuchos::rcp(new std::vector<double>(m_));
    Teuchos::RCP<std::vector<double>> constrappvec = Teuchos::rcp(new std::vector<double>(m_));

    double* yptr = y_mma_->Values();
    double* lamptr = lam_->Values();
    double* muptr = mu_->Values();
    double* sptr = s_->Values();
    double* rhoptr = rho_->Values();
    double* constrappptr = constr_appr_->Values();

    std::vector<double>::iterator yvecit = yvec->begin();
    std::vector<double>::iterator lamvecit = lamvec->begin();
    std::vector<double>::iterator muvecit = muvec->begin();
    std::vector<double>::iterator svecit = svec->begin();
    std::vector<double>::iterator rhovecit = rhovec->begin();
    std::vector<double>::iterator constrappvecit = constrappvec->begin();

    for (int i = 0; i < m_; i++)
    {
      *yvecit = *yptr;
      *lamvecit = *lamptr;
      *muvecit = *muptr;
      *svecit = *sptr;
      *rhovecit = *rhoptr;
      *constrappvecit = *constrappptr;

      yptr++;
      yvecit++;
      lamptr++;
      lamvecit++;
      muptr++;
      muvecit++;
      sptr++;
      svecit++;
      rhoptr++;
      rhovecit++;
      constrappptr++;
      constrappvecit++;
    }

    output_->WriteRedundantDoubleVector("y_mma", yvec);
    output_->WriteRedundantDoubleVector("lam", lamvec);
    output_->WriteRedundantDoubleVector("mu", muvec);
    output_->WriteRedundantDoubleVector("s", svec);
    output_->WriteRedundantDoubleVector("rho", rhovec);
    output_->WriteRedundantDoubleVector("constr_app", constrappvec);


    // optimization algorithm data - part 3 (doubles)
    output_->WriteDouble("zet", zet_);
    output_->WriteDouble("z_mma", z_mma_);
    output_->WriteDouble("rho0", rho0_);
    output_->WriteDouble("obj_app", obj_appr_);
  }

  return;
}


void OPTI::GCMMA::OutputToGmsh()
{
  // turn on/off screen output for writing process of Gmsh postprocessing file
  const bool screen_out = false;

  // create Gmsh postprocessing file
  const std::string filename = IO::GMSH::GetNewFileNameAndDeleteOldFiles(
      "optimization_field", total_iter_, 500, screen_out, discret_->Comm().MyPID());
  std::ofstream gmshfilecontent(filename.c_str());

  {
    // add 'View' to Gmsh postprocessing file
    gmshfilecontent << "View \" "
                    << "X \" {" << std::endl;

    // draw scalar field 'x' for every element for element or node-based field
    if (dens_type_ == INPAR::TOPOPT::dens_node_based)
      IO::GMSH::ScalarFieldToGmsh(discret_, x_, gmshfilecontent);
    else if (dens_type_ == INPAR::TOPOPT::dens_ele_based)
      IO::GMSH::ScalarElementFieldToGmsh(discret_, x_, gmshfilecontent);
    else
      dserror("not implemented");

    gmshfilecontent << "};" << std::endl;
  }

  gmshfilecontent.close();
  if (screen_out) std::cout << " done" << std::endl;
}


Teuchos::RCP<IO::DiscretizationReader> OPTI::GCMMA::ReadRestart(int step)
{
  // discretization reader with opti filename (restart number reduced by one here for correct file
  // names)
  Teuchos::RCP<IO::DiscretizationReader> reader =
      Teuchos::rcp(new IO::DiscretizationReader(discret_,
          Teuchos::rcp(new IO::InputControl(
              TOPOPT::modifyFilename(output_->Output()->FileName(), "", (bool)step))),
          step));

  // iteration counter
  total_iter_ = reader->ReadInt("step");
  outer_iter_ = reader->ReadInt("out_it");
  inner_iter_ = reader->ReadInt("in_it");

  // optimization algorithm data - part 1 (Epetra_Vector)
  if (dens_type_ == INPAR::TOPOPT::dens_node_based)
    reader->ReadVector(x_mma_, "x_mma_n");
  else if (dens_type_ == INPAR::TOPOPT::dens_ele_based)
    reader->ReadVector(x_mma_, "x_mma_e");
  else
    dserror("undefined type of optimization field");
  reader->ReadVector(x_, "x");
  reader->ReadVector(x_old_, "x_old");
  reader->ReadVector(x_old2_, "x_old2");
  reader->ReadVector(xsi_, "xsi");
  reader->ReadVector(eta_, "eta");
  reader->ReadVector(asymp_max_, "asymp_max");
  reader->ReadVector(asymp_min_, "asymp_min");
  reader->ReadVector(obj_deriv_, "obj_deriv");
  reader->ReadMultiVector(constr_deriv_, "constr_deriv");

  // optimization algorithm data - part 2 (Epetra_SerialDenseVector)
  Teuchos::RCP<std::vector<double>> yvec = Teuchos::rcp(new std::vector<double>(m_));
  Teuchos::RCP<std::vector<double>> lamvec = Teuchos::rcp(new std::vector<double>(m_));
  Teuchos::RCP<std::vector<double>> muvec = Teuchos::rcp(new std::vector<double>(m_));
  Teuchos::RCP<std::vector<double>> svec = Teuchos::rcp(new std::vector<double>(m_));
  Teuchos::RCP<std::vector<double>> rhovec = Teuchos::rcp(new std::vector<double>(m_));
  Teuchos::RCP<std::vector<double>> constrvec = Teuchos::rcp(new std::vector<double>(m_));
  Teuchos::RCP<std::vector<double>> constrappvec = Teuchos::rcp(new std::vector<double>(m_));

  reader->ReadRedundantDoubleVector(yvec, "y_mma");
  reader->ReadRedundantDoubleVector(lamvec, "lam");
  reader->ReadRedundantDoubleVector(muvec, "mu");
  reader->ReadRedundantDoubleVector(svec, "s");
  reader->ReadRedundantDoubleVector(rhovec, "rho");
  reader->ReadRedundantDoubleVector(constrvec, "constr");
  reader->ReadRedundantDoubleVector(constrappvec, "constr_app");

  double* yptr = y_mma_->Values();
  double* lamptr = lam_->Values();
  double* muptr = mu_->Values();
  double* sptr = s_->Values();
  double* rhoptr = rho_->Values();
  double* constrptr = constr_->Values();
  double* constrappptr = constr_appr_->Values();

  std::vector<double>::iterator yvecit = yvec->begin();
  std::vector<double>::iterator lamvecit = lamvec->begin();
  std::vector<double>::iterator muvecit = muvec->begin();
  std::vector<double>::iterator svecit = svec->begin();
  std::vector<double>::iterator rhovecit = rhovec->begin();
  std::vector<double>::iterator constrvecit = constrvec->begin();
  std::vector<double>::iterator constrappvecit = constrappvec->begin();

  for (int i = 0; i < m_; i++)
  {
    *yptr = *yvecit;
    *lamptr = *lamvecit;
    *muptr = *muvecit;
    *sptr = *svecit;
    *rhoptr = *rhovecit;
    *constrptr = *constrvecit;
    *constrappptr = *constrappvecit;

    yptr++;
    yvecit++;
    lamptr++;
    lamvecit++;
    muptr++;
    muvecit++;
    sptr++;
    svecit++;
    rhoptr++;
    rhovecit++;
    constrptr++;
    constrvecit++;
    constrappptr++;
    constrappvecit++;
  }

  // optimization algorithm data - part 3 (doubles)
  zet_ = reader->ReadDouble("zet");
  z_mma_ = reader->ReadDouble("z_mma");
  rho0_ = reader->ReadDouble("rho0");
  obj_ = reader->ReadDouble("obj");
  obj_appr_ = reader->ReadDouble("obj_app");

  return reader;
}
