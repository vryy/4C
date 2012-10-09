/*!------------------------------------------------------------------------------------------------*
\file opti_GCMMA.cpp

\brief 

<pre>
Maintainer: Martin Winklmaier
            winklmaier@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
 *------------------------------------------------------------------------------------------------*/


#include "opti_GCMMA.H"

#include "../drt_lib/drt_discret.H"
#include "../headers/definitions.h"
#include <Epetra_SerialDenseSolver.h>


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
OPTI::GCMMA::GCMMA(
    Teuchos::RCP<const DRT::Discretization> discret,
    const Teuchos::ParameterList& params,
    Teuchos::RCP<Epetra_Vector> x,
    int numConstraints,
    Teuchos::RCP<Epetra_Vector> x_min,
    Teuchos::RCP<Epetra_Vector> x_max
) :
discret_(discret),
total_iter_(0),
outer_iter_(0),
inner_iter_(0),
max_total_iter_(params.get<int>("MAX_ITER")),
max_inner_iter_(min(max_total_iter_,100)),
m_(numConstraints),
n_(x->MyLength()),
x_(rcp(new Epetra_Vector(*x))),
x_old_(rcp(new Epetra_Vector(*x))),
x_old2_(rcp(new Epetra_Vector(*x))),
x_mma_(rcp(new Epetra_Vector(*x))),
x_diff_min_(1.0e-5),
obj_(0.0),
obj_deriv_(rcp(new Epetra_Vector(x->Map()))),
obj_appr_(0.0),
constr_(rcp(new Epetra_SerialDenseVector(m_))),
constr_deriv_(rcp(new Epetra_MultiVector(x->Map(),m_))),
constr_appr_(rcp(new Epetra_SerialDenseVector(m_))),
p0_(rcp(new Epetra_Vector(x->Map()))),
q0_(rcp(new Epetra_Vector(x->Map()))),
r0_(0.0),
P_(rcp(new Epetra_MultiVector(x->Map(),m_))),
Q_(rcp(new Epetra_MultiVector(x->Map(),m_))),
b_(rcp(new Epetra_SerialDenseVector(m_))),
rho0_(0.01),
rho_(rcp(new Epetra_SerialDenseVector(m_))),
rho0min_(1.0e-6),
rhomin_(rcp(new Epetra_SerialDenseVector(m_))),
y_mma_(rcp(new Epetra_SerialDenseVector(m_))),
z_mma_(0.0),
xsi_(rcp(new Epetra_Vector(x->Map()))),
eta_(rcp(new Epetra_Vector(x->Map()))),
lam_(rcp(new Epetra_SerialDenseVector(m_))),
mu_(rcp(new Epetra_SerialDenseVector(m_))),
zet_(1.0),
a0_(1.0),
a_(rcp(new Epetra_SerialDenseVector(m_))),
c_(rcp(new Epetra_SerialDenseVector(m_))),
d_(rcp(new Epetra_SerialDenseVector(m_))),
tol_sub_min_(1.0e-09),
tol_(1.0e-5),
res_norm_(1.0e-6),
inc_norm_(1.0e-6),
facmin_(1.0e-10),
s_(rcp(new Epetra_SerialDenseVector(m_)))
{
  if (m_>100)
    dserror("current implementation inefficient for large number of constraints due to array structure and used solver");

  if (x_min==Teuchos::null)
  {
    x_min_ = Teuchos::rcp(new Epetra_Vector(x_->Map(),true));
    cout << endl << "WARNING: Initialized lower boundary for optimization variables with zeros" << endl << endl;
  }
  else
    x_min_ = rcp(new Epetra_Vector(*x_min));

  if (x_max==Teuchos::null)
  {
    x_max_ = Teuchos::rcp(new Epetra_Vector(x_->Map(),false));
    cout << endl << "WARNING: Initialized upper boundary for optimization variables with ones" << endl << endl;
    x_max_->PutScalar(1.0);
  }
  else
    x_max_ = rcp(new Epetra_Vector(*x_max));

  // TODO remove these two lines!!!
  x_min_->PutScalar(-2.0);
  x_max_->PutScalar(2.0);

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
  // currently gradients are always computed
  if ((inner_iter_==0) and // new outer_iter
      (objectivegrad == Teuchos::null or constraintsgrad == Teuchos::null))
    dserror("gradients must be given in new outer iteration");

  // initialization of old values
  if (outer_iter_ == 0)
  {
    // TODO test case
    if (x_->Map().NumMyElements()%3!=0)
      dserror("cannot be tested");

    int l = x_->Map().NumMyElements()/3;
    double val[3];

    for(int i=0;i<l;i++)
    {
      double alphai = (3*(i+1)-2*l)*M_PI/(6*l);

      val[0] = cos(alphai + M_PI/12);
      val[1] = sin(alphai + M_PI/12);
      val[2] = sin(2*alphai + M_PI/6);

      int indices[] = {i,i+l,i+2*l};
      x_->ReplaceMyValues(3,(double*)val,(int*)indices); // lnodeid = ldofid
    }



    *x_old_ = *x_;
    *x_old2_ = *x_;
    *x_mma_ = *x_;
  }



  // TODO testing example -> remove
  objective = 0.0;
  (*constraints)[0] = 0.0;

  int l = x_mma_->Map().NumMyElements()/3;

  double alpha = 0.0;
  for (int i=0;i<l;i++)
  {
    alpha = M_PI*(3*(i+1)-2*l)/(6*l);

    objective += cos(alpha)*(*x_mma_)[i] + sin(alpha)*(*x_mma_)[i+l] - 0.1*(*x_mma_)[i+2*l];

    if (inner_iter_==0)
    {
      (*objectivegrad)[i] = cos(alpha);
      (*objectivegrad)[i+l] = sin(alpha);
      (*objectivegrad)[i+2*l] = -0.1;
    }

    (*constraints)[0] += (*x_mma_)[i]*(*x_mma_)[i] + (*x_mma_)[i+l]*(*x_mma_)[i+l] + (*x_mma_)[i+2*l]*(*x_mma_)[i+2*l];

    if (inner_iter_==0)
    {
      (*constraintsgrad)[0][i] = 2*(*x_mma_)[i];
      (*constraintsgrad)[0][i+l] = 2*(*x_mma_)[i+l];
      (*constraintsgrad)[0][i+2*l] = 2*(*x_mma_)[i+2*l];
    }
  }

  (*constraints)[0] -= l + 1.0e-5;



  // update counters
  inner_iter_++;
  total_iter_++;

  if (inner_iter_==1) // new outer iter
  {
    // update counters
    outer_iter_++;

    // reset constraints and objective values
    double* constraint = constraints->Values();
    double* constr = constr_->Values();
    for (int i=0;i<m_;i++)
    {
      *constr = *constraint;

      constr++;
      constraint++;
    }
    *constr_deriv_ = *constraintsgrad;

    obj_ = objective;
    *obj_deriv_ = *objectivegrad;
    cout << "new obj is " << obj_ << endl;
    cout << "new obj deriv is " << *obj_deriv_ << endl;
    cout << "new constr are " << (*constr_)[0] << endl;
    cout << "new constr deriv are " << *constr_deriv_ << endl;

    // reset optimization variables
    *x_old2_ = *x_old_;
    *x_old_ = *x_;
    *x_ = *x_mma_;

    Asymptotes();

    InitRho();
  }
  else
  {
cout << "new obj is " << objective << endl;
cout << "new constr are " << (*constraints)[0] << endl;
    UpdateRho(objective,constraints);
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void OPTI::GCMMA::Asymptotes()
{
  // compute new asymptotes
  if (outer_iter_<3)
  {
    *asymp_min_ = *x_;
    *asymp_max_ = *x_;

    asymp_min_->Update(-0.5,*x_diff_,1.0);
    asymp_max_->Update(+0.5,*x_diff_,1.0);
  }
  else
  {
    double* xval = x_->Values();
    double* xold = x_old_->Values();
    double* xold2 = x_old2_->Values();
    double* asy_min = asymp_min_->Values();
    double* asy_max = asymp_max_->Values();
    double* xdiff = x_diff_->Values();

    for (int i=0;i<n_;i++)
    {
      double val = (*xval-*xold)*(*xold-*xold2);

      double fac = 1.0;
      if (val<0)
        fac = 0.7;
      else if (val>0)
        fac = 1.2;

      *asy_min = *xval - fac*(*xold-*asy_min);
      *asy_max = *xval + fac*(*asy_max-*xold);

      *asy_min = max(*asy_min,*xval-10**xdiff);
      *asy_min = min(*asy_min,*xval-0.01**xdiff);

      *asy_max = min(*asy_max,*xval+10**xdiff);
      *asy_max = max(*asy_max,*xval+0.01**xdiff);

      xval++;
      xold++;
      xold2++;
      asy_min++;
      asy_max++;
      xdiff++;
    }
  }
  cout << "low is " << *asymp_min_ << endl;
  cout << "upp is " << *asymp_max_ << endl;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void OPTI::GCMMA::InitRho()
{
  // set rho of objective function
  double rho0loc = 0.0;
  double* obj_deriv = obj_deriv_->Values();
  double* xdiff = x_diff_->Values();

  for (int i=0;i<n_;i++)
  {
    rho0loc += abs(*obj_deriv)**xdiff;

    obj_deriv++;
    xdiff++;
  }

  // communicate values
  discret_->Comm().SumAll(&rho0loc,&rho0_,1);

  // apply minimal value
  rho0_ = max(rho0min_,rho0_/(10*n_));


  // set rho of constraints
  Epetra_SerialDenseVector rholoc(m_);
  double* rho = rholoc.Values();

  for (int j=0;j<m_;j++)
  {
    *rho = 0.0;

    Epetra_Vector constr_der(View,*constr_deriv_,j);
    double* constr_deriv = constr_der.Values();
    xdiff = x_diff_->Values();

    for (int i=0;i<n_;i++)
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
    *rho = max(*rhomin,*rho/(10*n_));

    rho++;
    rhomin++;
  }
  cout << "rho0 is " << rho0_ << endl;
  cout << "rho is " << (*rho_)[0] << endl;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void OPTI::GCMMA::UpdateRho(
    double& objective,
    Teuchos::RCP<Epetra_SerialDenseVector> constraints
)
{
  double fac = 0.0;
  double facloc = 0.0;

  double* x = x_->Values();
  double* xmma = x_mma_->Values();
  double* xdiff = x_diff_->Values();
  double* asy_min = asymp_min_->Values();
  double* asy_max = asymp_max_->Values();

  for (int i=0;i<n_;i++)
  {
    facloc += (*xmma-*x)/(*asy_max-*xmma) * (*xmma-*x)/(*xmma-*asy_min) * (*asy_max-*asy_min)/max(x_diff_min_,*xdiff);

    x++;
    xmma++;
    asy_max++;
    asy_min++;
    xdiff++;
  }

  // communicate values
  discret_->Comm().SumAll(&facloc,&fac,1);

  fac = max(fac,facmin_);

  if (objective > obj_appr_ + 0.5*tol_sub_min_)
  {
    rho0_ = min(1.1*(rho0_ + (objective-obj_appr_)/fac),10*rho0_);
  }

  double* constr = constraints->Values();
  double* constr_appr = constr_appr_->Values();
  double* rho = rho_->Values();

  for (int i=0;i<m_;i++)
  {
    if (*constr > *constr_appr + 0.5*tol_sub_min_)
    {
      *rho = min(1.1*((*rho) + (*constr-*constr_appr)/fac),10*(*rho));
    }

    constr++;
    constr_appr++;
    rho++;
  }
  cout << "rho0 is " << rho0_ << endl;
  cout << "rho is " << (*rho_)[0] << endl;
}


bool OPTI::GCMMA::Converged(
    double& objective,
    Teuchos::RCP<Epetra_Vector> objectivegrad,
    Teuchos::RCP<Epetra_SerialDenseVector> constraints,
    Teuchos::RCP<Epetra_MultiVector> constraintsgrad
)
{
  // TODO testing example -> remove
  objective = 0.0;
  (*constraints)[0] = 0.0;

  int l = x_mma_->Map().NumMyElements()/3;

  double alpha;
  for (int i=0;i<l;i++)
  {
    alpha = M_PI*(3*(i+1)-2*l)/(6*l);

    objective += cos(alpha)*(*x_mma_)[i] + sin(alpha)*(*x_mma_)[i+l] - 0.1*(*x_mma_)[i+2*l];
    (*constraints)[0] += (*x_mma_)[i]*(*x_mma_)[i] + (*x_mma_)[i+l]*(*x_mma_)[i+l] + (*x_mma_)[i+2*l]*(*x_mma_)[i+2*l];
  }

  (*constraints)[0] -= l + 1.0e-5;

  if ((InnerConvergence(objective,constraints)==true) and (outer_iter_>0))
  {
    for (int i=0;i<l;i++)
    {
      alpha = M_PI*(3*(i+1)-2*l)/(6*l);

      (*objectivegrad)[i] = cos(alpha);
      (*objectivegrad)[i+l] = sin(alpha);
      (*objectivegrad)[i+2*l] = -0.1;

      (*constraintsgrad)[0][i] = 2*(*x_mma_)[i];
      (*constraintsgrad)[0][i+l] = 2*(*x_mma_)[i+l];
      (*constraintsgrad)[0][i+2*l] = 2*(*x_mma_)[i+2*l];
    }
  }




//  cout << "test if iter counters works well: total iter always increases by 1," << endl;
//  cout << "either inner or outer iter increase by 1, in second case inner iter shall be 1" << endl;
//  cout << "total iter is " << total_iter_ << endl;
//  cout << "inner iter is " << inner_iter_ << endl;
//  cout << "outer iter is " << outer_iter_ << endl;

  if (total_iter_==max_total_iter_)
  {
    cout << "WARNING: GCMMA optimization algorithm did not converge" << endl;
    return true;
  }


  bool converged = false;

  if ((InnerConvergence(objective,constraints)==true) and (outer_iter_>0))
  {
    inner_iter_ = 0;

    if (KKTCond(objective,objectivegrad,constraints,constraintsgrad)==false)
      converged = false;
    else
      converged = true;
  }
  else
    converged = false;

  // TODO check maybe also increment

  return converged;
}


void OPTI::GCMMA::FinishIteration(
    double& objective,
    Teuchos::RCP<Epetra_SerialDenseVector> constraints,
    bool& doGradient
)
{
  // TODO testing example -> remove
  objective = 0.0;
  (*constraints)[0] = 0.0;

  int l = x_mma_->Map().NumMyElements()/3;

  double alpha;
  for (int i=0;i<l;i++)
  {
    alpha = M_PI*(3*(i+1)-2*l)/(6*l);

    objective += cos(alpha)*(*x_mma_)[i] + sin(alpha)*(*x_mma_)[i+l] - 0.1*(*x_mma_)[i+2*l];
    (*constraints)[0] += (*x_mma_)[i]*(*x_mma_)[i] + (*x_mma_)[i+l]*(*x_mma_)[i+l] + (*x_mma_)[i+2*l]*(*x_mma_)[i+2*l];
  }

  (*constraints)[0] -= l + 1.0e-5;



  if (outer_iter_==0)
    doGradient = true;

  if (InnerConvergence(objective,constraints))
    doGradient = true;
  else
    doGradient = false;
}


bool OPTI::GCMMA::KKTCond(
    double& objective,
    Teuchos::RCP<Epetra_Vector> objectivegrad,
    Teuchos::RCP<Epetra_SerialDenseVector> constraints,
    Teuchos::RCP<Epetra_MultiVector> constraintsgrad
)
{
  /*
   * Global residual of primal dual interior point algorithm consists of:
   * res x/y/z
   * res lam/mu/std
   * res xsi/eta/zet
   */
  double resnorm = 0.0;
  double value = 0.0;



  Teuchos::RCP<Epetra_Vector> resX = rcp(new Epetra_Vector(*objectivegrad));
  resX->Update(-1.0,*xsi_,1.0);
  resX->Update(+1.0,*eta_,1.0);

  double resy = 0.0; // whole residual is of size m_
  double resz = a0_ - zet_;
  double reslam = 0.0; // whole residual is of size m_
  double resmu = 0.0; // whole residual is of size m_
  double res = 0.0; // whole residual is of size m_


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
    resX->Update(*lam,constr_deriv,1.0);

    resy = *c + *d**y-*mu-*lam;
    resz -= *a**lam;
    reslam = *constr - z_mma_**a - *y + *s;
    resmu = *mu**y;
    res = *lam**s;

    resnorm += resy*resy + reslam*reslam + resmu*resmu + res*res;

    c++;
    d++;
    y++;
    mu++;
    lam++;
    a++;
    constr++;
    s++;
 }

  resX->Norm2(&value);
  resnorm += value*value;

  resnorm += resz*resz;



  Epetra_Vector resXsi(x_->Map());
  Epetra_Vector resEta(x_->Map());

  double* resxsi = resXsi.Values();
  double* reseta = resEta.Values();
  double* xsi = xsi_->Values();
  double* eta = eta_->Values();
  double* x = x_mma_->Values();
  double* xmin = x_min_->Values();
  double* xmax = x_max_->Values();


  for (int i=0;i<n_;i++)
  {
    *resxsi = *xsi*(*x-*xmin);
    *reseta = *eta*(*xmax-*x);

    resxsi++;
    reseta++;
    xsi++;
    eta++;
    x++;
    xmax++;
    xmin++;
  }

  resXsi.Norm2(&value);
  resnorm += value*value;

  resEta.Norm2(&value);
  resnorm += value*value;



  double reszet = zet_*z_mma_;
  resnorm += reszet*reszet;




  resnorm = sqrt(resnorm);
  cout << "resnorm is " << resnorm << endl;
  if (resnorm<tol_)
    return true;

  return false;
}


bool OPTI::GCMMA::InnerConvergence(
    double& objective,
    Teuchos::RCP<Epetra_SerialDenseVector> constraints
)
{
  if (outer_iter_==0)
    return true;

  if (inner_iter_==max_inner_iter_)
  {
    cout << "WARNING: inner GCMMA optimization loop did not converge" << endl;
    return true;
  }

  if (obj_appr_ + tol_sub_min_ < objective)
    return false;

  double* constr_appr = constr_appr_->Values();
  double* constr = constraints->Values();

  for (int i=0;i<m_;i++)
  {
    if (*constr_appr + tol_sub_min_ < *constr)
      return false;

    constr_appr++;
    constr++;
  }

  return true;
}


void OPTI::GCMMA::InitSubSolve()
{
  double albeta = 0.1;

  alpha_->Update(1.0 - albeta,*asymp_min_,0.0);
  alpha_->Update(albeta,*x_,1.0);

  beta_->Update(1.0 - albeta,*asymp_max_,0.0);
  beta_->Update(albeta,*x_,1.0);

  double* alpha = alpha_->Values();
  double* beta = beta_->Values();
  double* xmin = x_min_->Values();
  double* xmax = x_max_->Values();

  for (int i=0;i<n_;i++)
  {
    *alpha = max(*alpha,*xmin);
    *beta = min(*beta,*xmax);

    alpha++;
    beta++;
    xmin++;
    xmax++;
  }


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

  for (int i=0;i<n_;i++)
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

  r0_ = obj_;

  for (int i=0;i<n_;i++)
  {
    *p0 = max(*obj_deriv,0.0);
    *q0 = max(-*obj_deriv,0.0);

    xdiffinv = 1.0/max(*xdiff,x_diff_min_);

    p0q0 = *p0 + *q0;
    *p0 = (*p0 + fac*p0q0 + rho0_*xdiffinv)**ux2ptr;
    *q0 = (*q0 + fac*p0q0 + rho0_*xdiffinv)**xl2ptr;

    r0_ -= *p0/(*ux1ptr) + *q0/(*xl1ptr);

    p0++;
    q0++;
    obj_deriv++;
    xdiff++;
    ux1ptr++;
    ux2ptr++;
    xl1ptr++;
    xl2ptr++;
  }


  double* b = b_->Values();
  double* constr = constr_->Values();

  double* rho = rho_->Values();
  for (int i=0;i<m_;i++)
  {
    Epetra_Vector p_i(View,*P_,i);
    Epetra_Vector q_i(View,*Q_,i);
    Epetra_Vector constr_deriv(View,*constr_deriv_,i);

    double* p = p_i.Values();
    double* q = q_i.Values();
    double* constr_der = constr_deriv.Values();

    double pq; // tmp variable

    *b = *constr;

    xdiff = x_diff_->Values();

    ux1ptr = ux1.Values();
    ux2ptr = ux2.Values();
    xl1ptr = xl1.Values();
    xl2ptr = xl2.Values();

    for (int j=0;j<n_;j++)
    {
      *p = max(*constr_der,0.0);
      *q = max(-*constr_der,0.0);

      xdiffinv = 1.0/max(*xdiff,x_diff_min_);

      pq = *p + *q;
      *p = (*p + fac*pq + *rho*xdiffinv)**ux2ptr;
      *q = (*q + fac*pq + *rho*xdiffinv)**xl2ptr;

      *b -= *p/(*ux1ptr) + *q/(*xl1ptr);

      p++;
      q++;
      constr_der++;
      xdiff++;
      ux1ptr++;
      ux2ptr++;
      xl1ptr++;
      xl2ptr++;
    }

    *b = -*b;

    b++;
    constr++;
    rho++;
  }

//  cout << "init subsolve: " << endl;
//  cout << "m is " << m_ << ", n is " << n_ << ", epsimin is " << tol_sub_min_ << endl;
//  cout << "low asy is " << *asymp_min_ << "upp asy is " << *asymp_max_ << endl;
//  cout << "alpha is " << *alpha_ << "beta is " << *beta_ << endl;
//  cout << "p0 is " << *p0_ << "q0 is " << *q0_ << endl;
//  cout << "p is " << *P_ << "q is " << *Q_ << endl;
//  cout << "a0 is " << a0_ << ", a is " << a_[0] << ", b is " << b_[0] << ", c is " << c_[0] << ", d is " << d_[0] << endl;
}


void OPTI::GCMMA::SubSolve()
{
  // Initialization
  x_mma_->Update(0.5,*alpha_,0.0);
  x_mma_->Update(0.5,*beta_,1.0);

  double* y = y_mma_->Values();
  double* lam = lam_->Values();
  for (int i=0;i<m_;i++)
  {
    *y = 1.0;
    *lam = 1.0;

    y++;
    lam++;
  }


  z_mma_ = 1.0;


  double* xsi = xsi_->Values();
  double* eta = eta_->Values();
  double* x = x_mma_->Values();
  double* alpha = alpha_->Values();
  double* beta = beta_->Values();
  for (int i=0;i<n_;i++)
  {
    *xsi = max(1.0, 1.0/(*x-*alpha));
    *eta = max(1.0, 1.0/(*beta-*x));

    xsi++;
    eta++;
    x++;
    alpha++;
    beta++;
  }


  double* mu = mu_->Values();
  double* s = s_->Values();
  double* c = c_->Values();
  for (int i=0;i<m_;i++)
  {
    *mu = max(1.0, *c/2.0);
    *s = 1.0;

    mu++;
    s++;
    c++;
  }

  int total_iter = 0;

  double tol_sub = 1.0;
  double tol_fac = 0.1; // factor from one tolerance to next (smaller) one
  if (tol_fac>0.99)
    dserror("factor for tolerance adaption shall be significant smaller than one!");

  bool tol_reached = false;

  while (tol_reached == false)
  {
    Epetra_Vector ux1(x_->Map(),false);
    Epetra_Vector ux2(x_->Map(),false);
    Epetra_Vector xl1(x_->Map(),false);
    Epetra_Vector xl2(x_->Map(),false);

    Epetra_Vector uxinv1(x_->Map(),false);
    Epetra_Vector xlinv1(x_->Map(),false);

    double* ux1ptr = ux1.Values();
    double* ux2ptr = ux2.Values();
    double* xl1ptr = xl1.Values();
    double* xl2ptr = xl2.Values();

    double* uxinv1ptr = uxinv1.Values();
    double* xlinv1ptr = xlinv1.Values();

    double* x = x_mma_->Values();
    double* asymp_max = asymp_max_->Values();
    double* asymp_min = asymp_min_->Values();

    Teuchos::RCP<Epetra_Vector> resX = rcp(new Epetra_Vector(x_->Map()));
    double* resxptr = resX->Values();

    xsi = xsi_->Values();
    eta = eta_->Values();

    for (int i=0;i<n_;i++)
    {
      *ux1ptr = *asymp_max-*x;
      *ux2ptr = *ux1ptr**ux1ptr;
      *xl1ptr = *x-*asymp_min;
      *xl2ptr = *xl1ptr**xl1ptr;

      *uxinv1ptr = 1.0/(*ux1ptr);
      *xlinv1ptr = 1.0/(*xl1ptr);


      double plam = (*p0_)[i];
      double qlam = (*q0_)[i];

      for (int j=0;j<m_;j++)
      {
        plam += (*P_)[j][i]*(*lam_)[j];
        qlam += (*Q_)[j][i]*(*lam_)[j];
      }

      double dpsidx = plam/(*ux2ptr) - qlam/(*xl2ptr);
      *resxptr = dpsidx - *xsi + *eta;

      ux1ptr++;
      ux2ptr++;
      xl1ptr++;
      xl2ptr++;
      uxinv1ptr++;
      xlinv1ptr++;
      resxptr++;
      x++;
      asymp_max++;
      asymp_min++;
      xsi++;
      eta++;
    }


    Epetra_SerialDenseVector resy(m_);
    double* c = c_->Values();
    double* d = d_->Values();
    double* y = y_mma_->Values();
    double* mu = mu_->Values();
    double* lam = lam_->Values();

    for (int j=0;j<m_;j++)
    {
      resy[j] = *c + *d**y - *mu - *lam;

      c++;
      d++;
      y++;
      mu++;
      lam++;
    }

    double resz = a0_ - zet_;
    double* a = a_->Values();
    lam = lam_->Values();

    for (int j=0;j<m_;j++)
    {
      resz -= *a**lam;

      a++;
      lam++;
    }


    Epetra_SerialDenseVector reslam(m_);

    Epetra_SerialDenseVector gvec1(m_);
    Epetra_SerialDenseVector gvec2(m_);

    uxinv1.Dot(*P_,gvec1.Values());
    xlinv1.Dot(*Q_,gvec2.Values());

    double* reslamptr = reslam.Values();
    double* gvec1ptr = gvec1.Values();
    double* gvec2ptr = gvec2.Values();
    a = a_->Values();
    y = y_mma_->Values();
    double* s = s_->Values();
    double* b = b_->Values();

    for (int j=0;j<m_;j++)
    {
      *reslamptr = *gvec1ptr + *gvec2ptr - *a*z_mma_ - *y + *s - *b;

      reslamptr++;
      gvec1ptr++;
      gvec2ptr++;
      a++;
      y++;
      s++;
      b++;
    }


    Teuchos::RCP<Epetra_Vector> resXsi = rcp(new Epetra_Vector(x_->Map()));
    Teuchos::RCP<Epetra_Vector> resEta = rcp(new Epetra_Vector(x_->Map()));

    double* resxsiptr = resXsi->Values();
    double* resetaptr = resEta->Values();

    xsi = xsi_->Values();
    eta = eta_->Values();
    x = x_mma_->Values();
    alpha = alpha_->Values();
    beta = beta_->Values();

    for (int i=0;i<n_;i++)
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


    Epetra_SerialDenseVector resmu(m_);
    Epetra_SerialDenseVector res(m_);

    double* resmuptr = resmu.Values();
    double* resptr = res.Values();
    mu = mu_->Values();
    y = y_mma_->Values();
    lam = lam_->Values();
    s = s_->Values();

    for (int j=0;j<m_;j++)
    {
      *resmuptr = *mu**y - tol_sub;
      *resptr = *lam**s - tol_sub;

      resmuptr++;
      mu++;
      y++;
      lam++;
      s++;
    }


    double reszet = zet_*z_mma_ - tol_sub;


    // sum up all residuals
    double locnorm = 0.0;
    double resnorm = 0.0;
    double resinf = abs(resz);

    resX->Norm2(&locnorm);
    resnorm += locnorm*locnorm;
    resX->NormInf(&locnorm);
    resinf = max(locnorm,resinf);

    double* resyptr = resy.Values();
    for (int j=0;j<m_;j++)
    {
      resnorm += *resyptr**resyptr;
      resinf = max(resinf,abs(*resyptr));
      resyptr++;
    }

    resnorm += resz*resz;

    reslamptr = reslam.Values();
    for (int j=0;j<m_;j++)
    {
      resnorm += *reslamptr**reslamptr;
      resinf = max(resinf,abs(*reslamptr));
      reslamptr++;
    }

    resXsi->Norm2(&locnorm);
    resnorm += locnorm*locnorm;
    resXsi->NormInf(&locnorm);
    resinf = max(locnorm,resinf);

    resEta->Norm2(&locnorm);
    resnorm += locnorm*locnorm;
    resEta->NormInf(&locnorm);
    resinf = max(locnorm,resinf);

    resmuptr = resmu.Values();
    for (int j=0;j<m_;j++)
    {
      resnorm += *resmuptr**resmuptr;
      resinf = max(resinf,abs(*resmuptr));
      resmuptr++;
    }

    resnorm += reszet*reszet;
    resinf = max(resinf,abs(reszet));

    resptr = res.Values();
    for (int j=0;j<m_;j++)
    {
      resnorm += *resptr**resptr;
      resinf = max(resinf,abs(*resptr));
      resptr++;
    }


    resnorm = sqrt(resnorm);

    int inner_iter = 0;
    int max_sub_inner_iter_ = 200;

    while (inner_iter<max_sub_inner_iter_)
    {
      inner_iter++;
      total_iter++;

      Epetra_Vector delx(x_->Map(),false);
      Epetra_Vector diagx(x_->Map(),false);
      Epetra_Vector ux3(x_->Map(),false);
      Epetra_Vector xl3(x_->Map(),false);

      Epetra_Vector uxinv2(x_->Map(),false);
      Epetra_Vector xlinv2(x_->Map(),false);

      double* ux1ptr = ux1.Values();
      double* ux2ptr = ux2.Values();
      double* ux3ptr = ux3.Values();
      double* xl1ptr = xl1.Values();
      double* xl2ptr = xl2.Values();
      double* xl3ptr = xl3.Values();

      double* uxinv1ptr = uxinv1.Values();
      double* uxinv2ptr = uxinv2.Values();
      double* xlinv1ptr = xlinv1.Values();
      double* xlinv2ptr = xlinv2.Values();

      double* x = x_mma_->Values();
      double* asymp_max = asymp_max_->Values();
      double* asymp_min = asymp_min_->Values();

      double* delxptr = delx.Values();
      double* diagxptr = diagx.Values();

      xsi = xsi_->Values();
      eta = eta_->Values();

      alpha = alpha_->Values();
      beta = beta_->Values();

      for (int i=0;i<n_;i++)
      {
        *ux1ptr = *asymp_max-*x;
        *ux2ptr = *ux1ptr**ux1ptr;
        *ux3ptr = *ux2ptr**ux1ptr;
        *xl1ptr = *x-*asymp_min;
        *xl2ptr = *xl1ptr**xl1ptr;
        *xl3ptr = *xl2ptr**xl1ptr;

        *uxinv1ptr = 1.0/(*ux1ptr);
        *uxinv2ptr = 1.0/(*ux1ptr**ux1ptr);
        *xlinv1ptr = 1.0/(*xl1ptr);
        *xlinv2ptr = 1.0/(*xl1ptr**xl1ptr);


        double plam = (*p0_)[i];
        double qlam = (*q0_)[i];

        for (int j=0;j<m_;j++)
        {
          plam += (*P_)[j][i]*(*lam_)[j];
          qlam += (*Q_)[j][i]*(*lam_)[j];
        }

        double dpsidx = plam/(*ux2ptr) - qlam/(*xl2ptr);

        *delxptr = dpsidx - tol_sub/(*x-*alpha) + tol_sub/(*beta-*x);
        *diagxptr = 2*(plam/(*ux3ptr) + qlam/(*xl3ptr)) + *xsi/(*x-*alpha) + *eta/(*beta-*x);

        ux1ptr++;
        ux2ptr++;
        ux3ptr++;
        xl1ptr++;
        xl2ptr++;
        xl3ptr++;
        uxinv1ptr++;
        uxinv2ptr++;
        xlinv1ptr++;
        xlinv2ptr++;
        delxptr++;
        diagxptr++;
        asymp_max++;
        asymp_min++;
        x++;
        alpha++;
        beta++;
        xsi++;
        eta++;
      }


      uxinv1.Dot(*P_,gvec1.Values());
      xlinv1.Dot(*Q_,gvec2.Values());

      Epetra_MultiVector GG(x_->Map(),m_);

      for (int i=0;i<m_;i++)
      {
        Epetra_Vector gg(View,GG,i);

        Epetra_Vector P(View,*P_,i);
        Epetra_Vector Q(View,*Q_,i);

        double* ggptr = gg.Values();
        double* pptr = P.Values();
        double* qptr = Q.Values();
        uxinv2ptr = uxinv2.Values();
        xlinv2ptr = xlinv2.Values();

        for (int j=0;j<n_;j++)
        {
          *ggptr += *pptr**uxinv2ptr - *qptr**xlinv2ptr;

          ggptr++;
          pptr++;
          qptr++;
          uxinv2ptr++;
          xlinv2ptr++;
        }
      }


      Epetra_SerialDenseVector dely(m_);

      double* delyptr = dely.Values();
      c = c_->Values();
      d = d_->Values();
      y = y_mma_->Values();
      lam = lam_->Values();

      for (int i=0;i<m_;i++)
      {
        *delyptr = *c + *d**y - *lam - tol_sub/(*y);

        delyptr++;
        c++;
        d++;
        y++;
        lam++;
      }


      double delz = a0_ - tol_sub/z_mma_;

      a = a_->Values();
      lam = lam_->Values();

      for (int i=0;i<m_;i++)
      {
        delz -= *a**lam;

        a++;
        lam++;
      }


      Epetra_SerialDenseVector dellam(m_);

      double* dellamptr = dellam.Values();
      gvec1ptr = gvec1.Values();
      gvec2ptr = gvec2.Values();
      a = a_->Values();
      y = y_mma_->Values();
      b = b_->Values();
      lam = lam_->Values();

      for (int i=0;i<m_;i++)
      {
        *dellamptr = *gvec1ptr + *gvec2ptr - *a*z_mma_ - *y - *b + tol_sub/(*lam);

        dellamptr++;
        gvec1ptr++;
        gvec2ptr++;
        a++;
        y++;
        b++;
        lam++;
      }


      Epetra_SerialDenseVector diagy(m_);

      double* diagyptr = diagy.Values();
      d = d_->Values();
      mu = mu_->Values();
      y = y_mma_->Values();

      for (int i=0;i<m_;i++)
      {
        *diagyptr = *d + *mu/(*y);

        diagyptr++;
        d++;
        mu++;
        y++;
      }


      Epetra_SerialDenseVector diaglamyi(m_);

      double* diaglamyiptr = diaglamyi.Values();
      s = s_->Values();
      lam = lam_->Values();
      diagyptr = diagy.Values();

      for (int i=0;i<m_;i++)
      {
        *diaglamyiptr = *s/(*lam) + 1.0/(*diagyptr);

        diaglamyiptr++;
        s++;
        lam++;
        diagyptr++;
      }


      Epetra_SerialDenseVector dlam(m_);
      double dz = 0.0;
      Epetra_Vector dx(x_->Map());

      if (m_ < x_->GlobalLength())
      {
        Epetra_SerialDenseVector blam(m_);

        double* blamptr = blam.Values();
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

          for (int j=0;j<n_;j++)
          {
            locval += *ggptr**delxptr/(*diagxptr);

            ggptr++;
            delxptr++;
            diagxptr++;
          }
          discret_->Comm().SumAll(&locval,&val,1);

          *blamptr = *dellamptr + *delyptr/(*diagyptr) - val;

          blamptr++;
          dellamptr++;
          delyptr++;
          diagyptr++;
        }


        Epetra_SerialDenseVector bb(m_+1);
        double* bbval = bb.A();
        blamptr = blam.Values();

        for (int i=0;i<m_;i++)
        {
          *bbval = *blamptr;
          bbval++;
          blamptr++;
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

            for (int j=0;j<n_;j++)
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


        double* dlamptr = dlam.Values();
        double* sol = solut.A();

        for (int i=0;i<m_;i++) // first m entries of sol
        {
          *dlamptr = *sol;

          dlamptr++;
          sol++;
        }


        dz = *sol; // last entry (=m+1) of sol


        dlamptr = dlam.Values();
        for (int i=0;i<m_;i++)
        {
          Epetra_Vector gg(View,GG,i);
          dx.Update(-*dlamptr,gg,1.0); // dx = -GG'*dlam

          dlamptr++;
        }

        double* dxptr = dx.Values();
        delxptr = delx.Values();
        diagxptr = diagx.Values();

        for (int i=0;i<n_;i++)
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

      for (int i=0;i<n_;i++)
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
        dlamptr++;
      }


      double dzet = -zet_ + tol_sub/z_mma_ - zet_*dz/z_mma_;


      double val = min(dz/z_mma_, dzet/zet_);

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
        val = min(val, *dyptr/(*y));
        val = min(val, *dlamptr/(*lam));
        val = min(val, *dmuptr/(*mu));
        val = min(val, *dsptr/(*s));

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

      for (int i=0;i<n_;i++)
      {
        val = min(val, *dxsiptr/(*xsi));
        val = min(val, *detaptr/(*eta));
        val = min(val, *dxptr/(*x-*alpha));
        val = min(val, *dxptr/(*x-*beta));

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
      val = max(1.0,fac*val);
      double steg = 1.0/val;


      int it = 0;
      int max_it = 50;
      double resnew = 2*resnorm;


      // save old values
      Epetra_Vector x_mma_old(*x_mma_);
      Epetra_Vector xsi_old(*xsi_);
      Epetra_Vector eta_old(*eta_);

      Epetra_SerialDenseVector y_old(m_);
      Epetra_SerialDenseVector lam_old(m_);
      Epetra_SerialDenseVector mu_old(m_);
      Epetra_SerialDenseVector s_old(m_);
      double z_mma_old = z_mma_;
      double zet_old = zet_;

      y = y_mma_->Values();
      lam = lam_->Values();
      mu = mu_->Values();
      s = s_->Values();
      double* yptr = y_old.Values();
      double* lamptr = lam_old.Values();
      double* muptr = mu_old.Values();
      double* sptr = s_old.Values();

      for (int j=0;j<m_;j++)
      {
        *yptr = *y;
        *lamptr = *lam;
        *muptr = *mu;
        *sptr = *s;

        yptr++;
        lamptr++;
        muptr++;
        sptr++;
        y++;
        lam++;
        mu++;
        s++;
      }

      while ((resnorm<resnew) and (it<max_it))
      {
        it++;

        x_mma_->Update(steg,dx,1.0,x_mma_old,0.0);
        xsi_->Update(steg,dxsi,1.0,xsi_old,0.0);
        eta_->Update(steg,deta,1.0,eta_old,0.0);

        y = y_mma_->Values();
        yptr = y_old.Values();
        dyptr = dy.Values();
        lam = lam_->Values();
        lamptr = lam_old.Values();
        dlamptr = dlam.Values();
        mu = mu_->Values();
        muptr = mu_old.Values();
        dmuptr = dmu.Values();
        s = s_->Values();
        sptr = s_old.Values();
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


        double* ux1ptr = ux1.Values();
        double* ux2ptr = ux2.Values();
        double* xl1ptr = xl1.Values();
        double* xl2ptr = xl2.Values();

        double* uxinv1ptr = uxinv1.Values();
        double* xlinv1ptr = xlinv1.Values();

        double* x = x_mma_->Values();
        double* asymp_max = asymp_max_->Values();
        double* asymp_min = asymp_min_->Values();

        Teuchos::RCP<Epetra_Vector> resX = rcp(new Epetra_Vector(x_->Map()));
        double* resxptr = resX->Values();

        xsi = xsi_->Values();
        eta = eta_->Values();

        for (int i=0;i<n_;i++)
        {
          *ux1ptr = *asymp_max-*x;
          *ux2ptr = *ux1ptr**ux1ptr;
          *xl1ptr = *x-*asymp_min;
          *xl2ptr = *xl1ptr**xl1ptr;

          *uxinv1ptr = 1.0/(*ux1ptr);
          *xlinv1ptr = 1.0/(*xl1ptr);


          double plam = (*p0_)[i];
          double qlam = (*q0_)[i];

          for (int j=0;j<m_;j++)
          {
            plam += (*P_)[j][i]*(*lam_)[j];
            qlam += (*Q_)[j][i]*(*lam_)[j];
          }

          double dpsidx = plam/(*ux2ptr) - qlam/(*xl2ptr);

          *resxptr = dpsidx - *xsi + *eta;

          ux1ptr++;
          ux2ptr++;
          xl1ptr++;
          xl2ptr++;
          uxinv1ptr++;
          xlinv1ptr++;
          resxptr++;
          x++;
          asymp_max++;
          asymp_min++;
          xsi++;
          eta++;
        }


        Epetra_SerialDenseVector resy(m_);
        double* c = c_->Values();
        double* d = d_->Values();
        double* y = y_mma_->Values();
        double* mu = mu_->Values();
        double* lam = lam_->Values();

        for (int j=0;j<m_;j++)
        {
          resy[j] = *c + *d**y - *mu - *lam;

          c++;
          d++;
          y++;
          mu++;
          lam++;
        }


        double resz = a0_ - zet_;

        double* a = a_->Values();
        lam = lam_->Values();

        for (int j=0;j<m_;j++)
        {
          resz -= *a**lam;

          a++;
          lam++;
        }


        Epetra_SerialDenseVector reslam(m_);

        Epetra_SerialDenseVector gvec1(m_);
        Epetra_SerialDenseVector gvec2(m_);

        uxinv1.Dot(*P_,gvec1.Values());
        xlinv1.Dot(*Q_,gvec2.Values());

        double* reslamptr = reslam.Values();
        double* gvec1ptr = gvec1.Values();
        double* gvec2ptr = gvec2.Values();
        a = a_->Values();
        y = y_mma_->Values();
        double* s = s_->Values();
        double* b = b_->Values();

        for (int j=0;j<m_;j++)
        {
          *reslamptr = *gvec1ptr + *gvec2ptr - *a*z_mma_ - *y + *s - *b;

          reslamptr++;
          gvec1ptr++;
          gvec2ptr++;
          a++;
          y++;
          s++;
          b++;
        }


        Teuchos::RCP<Epetra_Vector> resXsi = rcp(new Epetra_Vector(x_->Map()));
        Teuchos::RCP<Epetra_Vector> resEta = rcp(new Epetra_Vector(x_->Map()));

        double* resxsiptr = resXsi->Values();
        double* resetaptr = resEta->Values();

        xsi = xsi_->Values();
        eta = eta_->Values();
        x = x_mma_->Values();
        alpha = alpha_->Values();
        beta = beta_->Values();

        for (int i=0;i<n_;i++)
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


        Epetra_SerialDenseVector resmu(m_);
        Epetra_SerialDenseVector res(m_);

        double* resmuptr = resmu.Values();
        double* resptr = res.Values();
        mu = mu_->Values();
        y = y_mma_->Values();
        lam = lam_->Values();
        s = s_->Values();

        for (int j=0;j<m_;j++)
        {
          *resmuptr = (*mu)*(*y) - tol_sub;
          *resptr = *lam**s - tol_sub;

          resmuptr++;
          mu++;
          y++;
          lam++;
          s++;
        }


        double reszet = zet_*z_mma_ - tol_sub;


        // sum up all residuals
        resnew = 0.0;
        resinf = abs(resz);

        resX->Norm2(&locnorm);
        resnew += locnorm*locnorm;
        resX->NormInf(&locnorm);
        resinf = max(locnorm,resinf);

        double* resyptr = resy.Values();
        for (int j=0;j<m_;j++)
        {
          resnew += *resyptr**resyptr;
          resinf = max(resinf,abs(*resyptr));
          resyptr++;
        }

        resnew += resz*resz;

        reslamptr = reslam.Values();
        for (int j=0;j<m_;j++)
        {
          resnew += *reslamptr**reslamptr;
          resinf = max(resinf,abs(*reslamptr));
          reslamptr++;
        }

        resXsi->Norm2(&locnorm);
        resnew += locnorm*locnorm;
        resXsi->NormInf(&locnorm);
        resinf = max(locnorm,resinf);

        resEta->Norm2(&locnorm);
        resnew += locnorm*locnorm;
        resEta->NormInf(&locnorm);
        resinf = max(locnorm,resinf);

        resmuptr = resmu.Values();
        for (int j=0;j<m_;j++)
        {
          resnew += *resmuptr**resmuptr;
          resinf = max(resinf,abs(*resmuptr));
          resmuptr++;
        }

        resnew += reszet*reszet;
        resinf = max(resinf,abs(reszet));

        resptr = res.Values();
        for (int j=0;j<m_;j++)
        {
          resnew += *resptr**resptr;
          resinf = max(resinf,abs(*resptr));
          resptr++;
        }


        resnew = sqrt(resnew);
        steg = steg/2.0;

        if (it==max_it)
          cout << "reached maximal number of iterations in most inner loop of primal dual interior point optimization algorithm" << endl;
      }

      resnorm = resnew;
      // es würde genügen resinf hier (und damit nur einmal) zu berechnen!!!

      if (resinf < 0.9*tol_sub)
        break;
    }

    if (inner_iter==max_sub_inner_iter_)
      cout << "reached maximal number of iterations in inner loop of primal dual interior point optimization algorithm" << endl;

    if (tol_sub > 1.001*tol_sub_min_)
      tol_sub *= tol_fac;
    else
      tol_reached = true;
  }
  cout << "after subsolv:" << endl;
  cout << "x is " << *x_mma_ << endl;
//  cout << "y is " << y_mma_[0] << ", z is " << z_mma_ << ", lam is " << lam_[0] << endl;
//  cout << "xsi is " << *xsi_ << endl;
//  cout << "eta is " << *eta_ << endl;
//  cout << "mu is " << mu_[0] << ", zet is " << zet_ << ", s is " << s_[0] << endl;
}


void OPTI::GCMMA::Update(
)
{
  // Update data

  // helper vectors
  Epetra_Vector uxinv(x_->Map());
  Epetra_Vector xlinv(x_->Map());

  double* asymp_max = asymp_max_->Values();
  double* asymp_min = asymp_min_->Values();
  double* x = x_mma_->Values();

  double* uxinvptr = uxinv.Values();
  double* xlinvptr = xlinv.Values();


  for (int i=0;i<n_;i++)
  {
    *uxinvptr = 1.0/(*asymp_max-*x);
    *xlinvptr = 1.0/(*x-*asymp_min);

    asymp_max++;
    asymp_min++;
    x++;
    uxinvptr++;
    xlinvptr++;
  }

  // set new approximation of objective value
  double value = 0.0;

  uxinv.Dot(*p0_,&value);
  obj_appr_ = r0_ + value;

  xlinv.Dot(*q0_,&value);
  obj_appr_ += value;

  // set new approximation of constraints
  Epetra_SerialDenseVector values(m_);

  double* valptr = values.Values();
  uxinv.Dot(*P_,valptr);

  double* constrptr = constr_appr_->Values();
  double* b = b_->Values();
  valptr = values.Values();

  for (int i=0;i<m_;i++)
  {
    *constrptr = -*b + *valptr;

    constrptr++;
    b++;
    valptr++;
  }

  valptr = values.Values();
  xlinv.Dot(*Q_,valptr);

  constrptr = constr_appr_->Values();
  valptr = values.Values();

  for (int i=0;i<m_;i++)
  {
    *constrptr += *valptr;

    constrptr++;
    valptr++;
  }
}



