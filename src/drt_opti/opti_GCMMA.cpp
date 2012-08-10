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
params_(params),
outer_iter_(0),
inner_iter_(0),
max_outer_iter_(100),
max_inner_iter_(100),
m_(numConstraints),
n_(x->MyLength()),
x_(rcp(new Epetra_Vector(*x))),
x_old_(rcp(new Epetra_Vector(*x))),
x_old2_(rcp(new Epetra_Vector(*x))),
x_mma_(rcp(new Epetra_Vector(*x))),
x_diff_min_(1.0e-5),
obj_(0.0),
obj_mma_(0.0),
obj_deriv_(rcp(new Epetra_Vector(x->Map()))),
obj_appr_(0.0),
constr_(new double[m_]),
constr_mma_(new double[m_]),
constr_deriv_(rcp(new Epetra_MultiVector(x->Map(),m_))),
constr_appr_(new double[m_]),
p0_(rcp(new Epetra_Vector(x->Map()))),
q0_(rcp(new Epetra_Vector(x->Map()))),
r0_(0.0),
P_(rcp(new Epetra_MultiVector(x->Map(),m_))),
Q_(rcp(new Epetra_MultiVector(x->Map(),m_))),
r_(new double[m_]),
rho0_(0.01),
rho_(new double[m_]),
rho0min_(1.0e-6),
rhomin_(new double[m_]),
y_mma_(new double[m_]),
z_mma_(0.0),
xsi_(rcp(new Epetra_Vector(x->Map()))),
eta_(rcp(new Epetra_Vector(x->Map()))),
lam_(new double[m_]),
mu_(new double[m_]),
zet_(0.0),
a0_(1.0),
a_(new double[m_]),
c_(new double[m_]),
d_(new double[m_]),
tol_sub_min_(1.0e-09),
tol_(1.0e-5),
res_norm_(1.0e-6),
inc_norm_(1.0e-6),
facmin_(1.0e-10),
s_(new double[m_])
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

  double* p1 = a_;
  double* p2 = c_;
  double* p3 = d_;
  double* p4 = rho_;
  double* p5 = rhomin_;

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
void OPTI::GCMMA::Iterate(
    double& objective,
    Teuchos::RCP<Epetra_Vector> objectivegrad,
    double* constraints,
    Teuchos::RCP<Epetra_MultiVector> constraintsgrad,
    bool& innerIterConverged
)
{
  if (inner_iter_ == 0)
  {
    if (constraints == NULL or constraintsgrad == Teuchos::null)
      dserror("constraints and its gradients must be given in new outer iteration");

    constr_ = constraints;
    *constr_deriv_ = *constraintsgrad;

    obj_ = objective;
    *obj_deriv_ = *objectivegrad;

  }

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



    double* constr_mma = constr_mma_;
    double* constr = constr_;

    for (int i=0;i<m_;i++)
    {
      *constr_mma = *constr;

      constr++;
      constr_mma++;
    }

    obj_mma_ = obj_;

    x_old_ = x_;
    x_old2_ = x_;
    x_mma_ = x_;
  }



  // TODO testing example -> remove
  obj_ = 0.0;
  constr_[0] = 0.0;

  int l = x_->Map().NumMyElements()/3;
  double alpha;
  for (int i=0;i<l;i++)
  {
    alpha = M_PI*(3*(i+1)-2*l)/(6*l);

    obj_ += cos(alpha)*(*x_)[i] + sin(alpha)*(*x_)[i+l] - 0.1*(*x_)[i+2*l];

    (*obj_deriv_)[i] = cos(alpha);
    (*obj_deriv_)[i+l] = sin(alpha);
    (*obj_deriv_)[i+2*l] = -0.1;

    constr_[0] += (*x_)[i]*(*x_)[i] + (*x_)[i+l]*(*x_)[i+l] + (*x_)[i+2*l]*(*x_)[i+2*l];

    (*constr_deriv_)[0][i] = 2*(*x_)[i];
    (*constr_deriv_)[0][i+l] = 2*(*x_)[i+l];
    (*constr_deriv_)[0][i+2*l] = 2*(*x_)[i+2*l];
  }

  constr_[0] -= l + 1.0e-5;



  if (inner_iter_==0) // new outer iteration
  {
    outer_iter_++;
    inner_iter_++;

    PrepareOuterIter();
  }
  else
  {
    inner_iter_++;

    PrepareInnerIter();
  }

  InitSubSolve();

  SubSolve();
//  %
//  % Solving the subproblem by a primal-dual Newton method
//  [xmma,ymma,zmma,lam,xsi,eta,mu,zet,s] = ...
//  subsolv(m,n,epsimin,low,upp,alfa,beta,p0,q0,P,Q,a0,a,b,c,d);
//  %
//  % Calculations of f0app and fapp.
//  ux1 = upp-xmma;
//  xl1 = xmma-low;
//  uxinv = eeen./ux1;
//  xlinv = eeen./xl1;
//  f0app = r0 + p0'*uxinv + q0'*xlinv;
//  fapp  =  r +   P*uxinv +   Q*xlinv;

  Update(innerIterConverged);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void OPTI::GCMMA::PrepareOuterIter()
{
  Asymptotes();

  UpdateRho();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void OPTI::GCMMA::PrepareInnerIter()
{
  UpdateRho();
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

      *asy_max = min(*asy_max,*xval+0.01**xdiff);
      *asy_max = max(*asy_max,*xval+10**xdiff);

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
void OPTI::GCMMA::UpdateRho()
{
  if (inner_iter_==1)
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
    double* rholoc = new double[m_];
    double* rho = rholoc;

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
    discret_->Comm().SumAll(rholoc,rho_,m_);

    // apply minimal value
    rho = rho_;
    double* rhomin = rhomin_;

    for (int j=0;j<m_;j++)
    {
      *rho = max(*rhomin,*rho/(10*n_));

      rho++;
      rhomin++;
    }
  }
  else
  {
    double fac = 0.0;

    double* x = x_->Values();
    double* xmma = x_mma_->Values();
    double* xmax = x_max_->Values();
    double* xmin = x_min_->Values();
    double* xdiff = x_diff_->Values();
    double* asy_min = asymp_min_->Values();
    double* asy_max = asymp_max_->Values();
    for (int i=0;i<n_;i++)
    {
      fac += (*xmma-*x)/(*asy_max-*xmma) * (*xmma-*x)/(*xmma-*asy_min) * (*xmax-*xmin)/max(x_diff_min_,*xdiff);

      x++;
      xmax++;
      xmin++;
    }
    fac = max(fac,facmin_);

    if (obj_ > obj_appr_ + 0.5*tol_sub_min_)
    {
      rho0_ = min(1.1*(rho0_ + (obj_-obj_appr_)/fac),10*rho0_);
    }

    double* constr = constr_;
    double* constr_appr = constr_appr_;
    double* rho = rho_;

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
  }
}


bool OPTI::GCMMA::Converged()
{
//  /*
//   * Usually it notfinished is initially set true and it is checked if
//   * convergence is reached
//   *
//   * Since here in optimization usually only the objective function is checked
//   * and the density increment usually not, it is done here the other way:
//   * setting notfinished false and checking whether convergence is not reached
//   *
//   * winklmaier 12/11
//   */
//  iter_++;
//  if (max_iter_<1) dserror("maximal number of optimization steps smaller than 1");
//
//  bool finished = true;
//
//  if (iter_==1)
//  {
//    finished = false;
//    // TODO initial output
//  }
//  else
//  {
//    if (iter_>max_iter_)
//    {
//      finished = true;
//      // TODO output
//    }
//    else // TODO check convergence also if maxiter reached?
//    {
//      if ((conv_check_type_==INPAR::TOPOPT::inc) or
//          (conv_check_type_==INPAR::TOPOPT::inc_and_res))
//      {
//        Epetra_Vector inc(*optimizer_->RowMap(),false);
////        inc.Update(1.0,*optimizer_->DensityIp(),-1.0,*optimizer_->DensityI(),0.0);
//
//        double incvelnorm;
//        inc.Norm2(&incvelnorm);
//
//        if (incvelnorm>inc_tol_)
//          finished = false;
//      }
//
//      if ((conv_check_type_==INPAR::TOPOPT::res) or
//          (conv_check_type_==INPAR::TOPOPT::inc_and_res))
//      {
//        if (fabs(objective_ip_-objective_i_)>res_tol_)
//          finished = false;
//      }
//      // TODO output
//    }
//  }
  // no check in first iteration
  if (outer_iter_==0)
    return false;

  // check only new outer iterations
  if (inner_iter_!=0)
    return false;

  // compute residual for original variables x
  Teuchos::RCP<Epetra_Vector> resX = rcp(new Epetra_Vector(*obj_deriv_));

  resX->Update(-1.0,*xsi_,1.0);
  resX->Update(+1.0,*eta_,1.0);

  double* lam = lam_;

  for (int i=0;i<m_;i++)
  {
    Epetra_Vector constr_deriv(View,*constr_deriv_,i);

    resX->Update(*lam,constr_deriv,1.0);

    lam++;
  }



  double* resy = new double[m_];
  double* resyptr = resy;
  double* c = c_;
  double* d = d_;
  double* y = y_mma_;
  double* mu = mu_;
  lam = lam_;

  for (int i=0;i<m_;i++)
  {
    *resyptr = *c + *d**y-*mu-*lam;

    resyptr++;
    c++;
    d++;
    y++;
    mu++;
    lam++;
  }

  double resz = a0_ - zet_;
  double* a = a_;
  lam = lam_;

  for (int i=0;i<m_;i++)
  {
    resz -= *a**lam;

    a++;
    lam++;
  }


  double* relam = new double[m_];
  double* constr = constr_;
  a = a_;
  y = y_mma_;
  double* s = s_;

  for (int i=0;i<m_;i++)
  {
    *relam = *constr - z_mma_**a - *y + *s;

    relam++;
    constr++;
    a++;
    y++;
    s++;
  }


  double* rexsi;
  double* reeta;
  double* xsi = xsi_->Values();
  double* eta = eta_->Values();
  double* x = x_->Values();
  double* xmin = x_min_->Values();
  double* xmax = x_max_->Values();


  for (int i=0;i<n_;i++)
  {
    *rexsi = *xsi*(*x-*xmin);
    *reeta = *eta*(*xmax-*x);

    rexsi++;
    reeta++;
    xsi++;
    x++;
    xmax++;
    xmin++;
  }


  double* remu;
  mu = mu_;
  y = y_mma_;

  for (int i=0;i<m_;i++)
  {
    *remu = *mu**y;

    remu++;
    mu++;
    y++;
    s++;
  }


  double reszet = zet_*z_mma_;


  double* res;
  s = s_;
  lam = lam_;

  for (int i=0;i<m_;i++)
  {
    *res = *lam**s;

    res++;
    s++;
    lam++;
  }

  double resnorm;
  resX->Norm2(&resnorm);
  resnorm += resnorm*resnorm + resz*resz;

  resyptr = resy;
  for (int i=0;i<m_;i++)
  {
    resnorm += *resyptr**resyptr;
    resyptr++;
  }

  resnorm = sqrt(resnorm);

  if (resnorm<tol_)
    return true;

  return false;
}


bool OPTI::GCMMA::InnerConvergence()
{
  if (obj_appr_ + tol_sub_min_ < obj_)
    return false;

  double* constr_appr = constr_appr_;
  double* constr = constr_;

  for (int i=0;i<m_;i++)
  {
    if (*constr_appr + tol_sub_min_ < *constr)
      return false;
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


  double* r = r_;
  double* constr = constr_;

  double* rho = rho_;
  for (int i=0;i<m_;i++)
  {
    Epetra_Vector p_i(View,*P_,i);
    Epetra_Vector q_i(View,*Q_,i);
    Epetra_Vector constr_deriv(View,*constr_deriv_,i);

    double* p = p_i.Values();
    double* q = q_i.Values();
    double* constr_der = constr_deriv.Values();

    double pq; // tmp variable

    *r = *constr;

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

      *r -= *p/(*ux1ptr) + *q/(*xl1ptr);

      p++;
      q++;
      constr_der++;
      xdiff++;
      ux1ptr++;
      ux2ptr++;
      xl1ptr++;
      xl2ptr++;
    }

    *r = -*r;

    r++;
    constr++;
    rho++;
  }

//  cout << "m is " << m_ << ", n is " << n_ << ", epsimin is " << tol_sub_ << endl;
//  cout << "low asy is " << *asymp_min_ << "upp asy is " << *asymp_max_ << endl;
//  cout << "alpha is " << *alpha_ << "beta is " << *beta_ << endl;
//  cout << "p0 is " << p0_ << "q0 is " << q0_ << endl;
//  cout << "p is " << P_ << "q is " << Q_ << endl;
//  cout << "a0 is " << a0_ << ", a is " << a_[0] << ", b is " << r[0] << ", c is " << c_[0] << ", d is " << d_[0] << endl;
}


void OPTI::GCMMA::SubSolve()
{
  // Initialization
  x_mma_->Update(0.5,*alpha_,0.0);
  x_mma_->Update(0.5,*beta_,1.0);

  double* y = y_mma_;
  double* lam = lam_;
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


  double* mu = mu_;
  double* s = s_;
  double* c = c_;
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

    double* x = x_->Values();
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


      double plam = (*p0_)[n_];
      double qlam = (*q0_)[n_];

      for (int j=0;j<m_;j++)
      {
        plam += (*P_)[j][i]*lam_[j];
        qlam += (*Q_)[j][i]*lam_[j];
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
    }


    double* resy = new double[m_];
    double* c = c_;
    double* d = d_;
    double* y = y_mma_;
    double* mu = mu_;
    double* lam = lam_;

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

    double* a = a_;
    lam = lam_;

    for (int j=0;j<m_;j++)
    {
      resz -= *a**lam;

      a++;
      lam++;
    }


    double* reslam = new double[m_];

    double* gvec1 = new double[m_];
    double* gvec2 = new double[m_];

    uxinv1.Dot(*P_,gvec1);
    xlinv1.Dot(*Q_,gvec2);

    double* reslamptr = reslam;
    double* gvec1ptr = gvec1;
    double* gvec2ptr = gvec2;
    a = a_;
    y = y_mma_;
    double* s = s_;
    double* r = r_;

    for (int j=0;j<m_;j++)
    {
      *reslamptr = *gvec1ptr + *gvec2ptr - *a*z_mma_ - *y + *s - *r;

      reslamptr++;
      gvec1ptr++;
      gvec2ptr++;
      a++;
      y++;
      s++;
      r++;
    }


    Teuchos::RCP<Epetra_Vector> resXsi = rcp(new Epetra_Vector(x_->Map()));
    Teuchos::RCP<Epetra_Vector> resEta = rcp(new Epetra_Vector(x_->Map()));

    double* resxsiptr = resXsi->Values();
    double* resetaptr = resEta->Values();

    xsi = xsi_->Values();
    eta = eta_->Values();
    x = x_->Values();
    alpha = alpha_->Values();
    beta = beta_->Values();

    for (int i=0;i<n_;i++)
    {
      *resxsiptr = *xsi*(*x-*alpha) - tol_sub;
      *resetaptr = *eta*(*beta-*x) - tol_sub;

      resxsiptr++;
      xsi++;
      eta++;
      x++;
      alpha++;
      beta++;
    }


    double* resmu = new double[m_];
    double* res = new double[m_];

    double* resmuptr = resmu;
    double* resptr = res;
    mu = mu_;
    y = y_mma_;
    lam = lam_;
    s = s_;

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

    double* resyptr = resy;
    for (int j=0;j<m_;j++)
    {
      resnorm += *resyptr**resyptr;
      resinf = max(resinf,abs(*resyptr));
      resyptr++;
    }

    resnorm += resz*resz;

    reslamptr = reslam;
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

    resmuptr = resmu;
    for (int j=0;j<m_;j++)
    {
      resnorm += *resmuptr**resmuptr;
      resinf = max(resinf,abs(*resmuptr));
      resmuptr++;
    }

    resnorm += reszet*reszet;
    resinf = max(resinf,abs(reszet));

    resptr = res;
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

      double* x = x_->Values();
      double* asymp_max = asymp_max_->Values();
      double* asymp_min = asymp_min_->Values();

      double* delxptr = delx.Values();
      double* diagxptr = diagx.Values();

      xsi = xsi_->Values();
      eta = eta_->Values();

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


        double plam = (*p0_)[n_];
        double qlam = (*q0_)[n_];

        for (int j=0;j<m_;j++)
        {
          plam += (*P_)[j][i]*lam_[j];
          qlam += (*Q_)[j][i]*lam_[j];
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
      }


      uxinv1.Dot(*P_,gvec1);
      xlinv1.Dot(*Q_,gvec2);

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


      double* dely = new double[m_];

      double* delyptr = dely;
      c = c_;
      d = d_;
      y = y_mma_;
      lam = lam_;

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

      a = a_;
      lam = lam_;

      for (int i=0;i<m_;i++)
      {
        delz -= *a**lam;

        a++;
        lam++;
      }


      double* dellam = new double[m_];

      double* dellamptr = dellam;
      gvec1ptr = gvec1;
      gvec2ptr = gvec2;
      a = a_;
      y = y_mma_;
      r = r_;
      lam = lam_;

      for (int i=0;i<m_;i++)
      {
        *dellamptr = *gvec1ptr + *gvec2ptr - *a*z_mma_ - *y - *r + tol_sub/(*lam);

        dellamptr++;
        gvec1ptr++;
        gvec2ptr++;
        a++;
        y++;
        r++;
        lam++;
      }


      double* diagy = new double[m_];

      double* diagyptr = diagy;
      d = d_;
      mu = mu_;
      y = y_mma_;

      for (int i=0;i<m_;i++)
      {
        *diagyptr = *d + *mu/(*y);

        diagyptr++;
        d++;
        mu++;
        y++;
      }


      double* diaglamyi = new double[m_];

      double* diaglamyiptr = diaglamyi;
      s = s_;
      lam = lam_;
      diagyptr = diagy;

      for (int i=0;i<m_;i++)
      {
        *diaglamyiptr = *s/(*lam) + 1.0/(*diagyptr);

        diaglamyiptr++;
        s++;
        lam++;
        diagyptr++;
      }



      double* dlam = new double[m_];
      double dz = 0.0;
      Epetra_Vector dx(x_->Map());

      if (m_ < x_->GlobalLength())
      {
        double* blam = new double[m_];

        double* blamptr = blam;
        dellamptr = dellam;
        delyptr = dely;
        diagyptr = diagy;

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
        blamptr = blam;

        for (int i=0;i<m_;i++)
        {
          *bbval = *blamptr;
          bbval++;
          blamptr++;
        }
        *bbval = delz;

        Epetra_SerialDenseMatrix AA(m_+1,m_+1);
        double* aa = AA.A();
        a = a_;
        for (int icol=0;icol<m_;icol++)
        {
          Epetra_Vector colGG(View,GG,icol);
          double* colgg;

          for (int irow=0;irow<m_;irow++)
          {
            colgg = colGG.Values();

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
            diaglamyiptr++;
          }
          *aa = *a;

          a++;
          aa++;
        }

        a = a_;
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


        double* dlamptr = dlam;
        double* sol = solut.A();

        for (int i=0;i<m_;i++)
        {
          *dlamptr = *sol;

          dlamptr++;
          sol++;
        }

        dz = *dlamptr;


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
      x = x_->Values();
      alpha = alpha_->Values();
      beta = beta_->Values();

      for (int i=0;i<m_;i++)
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


      double* dy = new double[m_];
      double* dmu = new double[m_];
      double* ds = new double[m_];

      double* dyptr = dy;
      double* dmuptr = dmu;
      double* dsptr = ds;
      delyptr = dely;
      diagyptr = diagy;
      double* dlamptr = dlam;
      mu = mu_;
      y = y_mma_;
      s = s_;
      lam = lam_;
      dlamptr = dlam;

      for (int i=0;i<m_;i++)
      {
        *dyptr = -*delyptr/(*diagyptr) + *dlamptr/(*diagy);
        *dmuptr = -*mu + tol_sub/(*y) - (*mu**dyptr)/(*y);
        *dsptr = -*s + tol_sub/(*lam) - (*s**dlam)/(*lam);

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
        dlam++;
      }


      double dzet = -zet_ + tol_sub/z_mma_ - zet_*dz/z_mma_;



      double val = min(z_mma_/dz, zet_/dzet);

      y = y_mma_;
      dyptr = dy;
      lam = lam_;
      dlamptr = dlam;
      mu = mu_;
      dmuptr = dmu;
      s = s_;
      dsptr = ds;

      for (int i=0;i<m_;i++)
      {
        val = min(val, *y/(*dyptr));
        val = min(val, *lam/(*dlamptr));
        val = min(val, *mu/(*dmuptr));
        val = min(val, *s/(*dsptr));

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
      x = x_->Values();
      alpha = alpha_->Values();
      beta = beta_->Values();

      for (int i=0;i<n_;i++)
      {
        val = min(val, *xsi/(*dxsiptr));
        val = min(val, *eta/(*detaptr));
        val = min(val, *dxptr/(*x-*alpha));
        val = min(val, *dxptr/(*beta-*x));

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

      // min becomes max since fac<0
      val = max(1.0,fac*val);
      double steg = 1.0/val;


      int it = 0;
      int max_it = 50;
      double resnew = 2*resnorm;

      while ((resnorm<resnew) and (it<max_it))
      {
        it++;

        x_->Update(steg,dx,1.0);
        xsi_->Update(steg,dxsi,1.0);
        eta_->Update(steg,deta,1.0);

        y = y_mma_;
        dyptr = dy;
        lam = lam_;
        dlamptr = dlam;
        mu = mu_;
        dmuptr = dmu;
        s = s_;
        dsptr = ds;
        for (int i=0;i<m_;i++)
        {
          *y = *y + steg**dyptr;
          *lam = *lam + steg**dlamptr;
          *mu = *mu + steg**dmuptr;
          *s = *s + steg**dsptr;
        }

        z_mma_ = z_mma_ + steg*dz;
        zet_ = zet_ + steg*dzet;



        double* ux1ptr = ux1.Values();
        double* ux2ptr = ux2.Values();
        double* xl1ptr = xl1.Values();
        double* xl2ptr = xl2.Values();

        double* uxinv1ptr = uxinv1.Values();
        double* xlinv1ptr = xlinv1.Values();

        double* x = x_->Values();
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
            plam += (*P_)[j][i]*lam_[j];
            qlam += (*Q_)[j][i]*lam_[j];
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
        }


        double* resy = new double[m_];
        double* c = c_;
        double* d = d_;
        double* y = y_mma_;
        double* mu = mu_;
        double* lam = lam_;

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

        double* a = a_;
        lam = lam_;

        for (int j=0;j<m_;j++)
        {
          resz -= *a**lam;

          a++;
          lam++;
        }


        double* reslam = new double[m_];

        double* gvec1 = new double[m_];
        double* gvec2 = new double[m_];

        uxinv1.Dot(*P_,gvec1);
        xlinv1.Dot(*Q_,gvec2);

        double* reslamptr = reslam;
        double* gvec1ptr = gvec1;
        double* gvec2ptr = gvec2;
        a = a_;
        y = y_mma_;
        double* s = s_;
        double* r = r_;

        for (int j=0;j<m_;j++)
        {
          *reslamptr = *gvec1ptr + *gvec2ptr - *a*z_mma_ - *y + *s - *r;

          reslamptr++;
          gvec1ptr++;
          gvec2ptr++;
          a++;
          y++;
          s++;
          r++;
        }


        Teuchos::RCP<Epetra_Vector> resXsi = rcp(new Epetra_Vector(x_->Map()));
        Teuchos::RCP<Epetra_Vector> resEta = rcp(new Epetra_Vector(x_->Map()));

        double* resxsiptr = resXsi->Values();
        double* resetaptr = resEta->Values();

        xsi = xsi_->Values();
        eta = eta_->Values();
        x = x_->Values();
        alpha = alpha_->Values();
        beta = beta_->Values();

        for (int i=0;i<n_;i++)
        {
          *resxsiptr = *xsi*(*x-*alpha) - tol_sub;
          *resetaptr = *eta*(*beta-*x) - tol_sub;

          resxsiptr++;
          xsi++;
          eta++;
          x++;
          alpha++;
          beta++;
        }


        double* resmu = new double[m_];
        double* res = new double[m_];

        double* resmuptr = resmu;
        double* resptr = res;
        mu = mu_;
        y = y_mma_;
        lam = lam_;
        s = s_;

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

        double* resyptr = resy;
        for (int j=0;j<m_;j++)
        {
          resnew += *resyptr**resyptr;
          resinf = max(resinf,abs(*resyptr));
          resyptr++;
        }

        resnew += resz*resz;

        reslamptr = reslam;
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

        resmuptr = resmu;
        for (int j=0;j<m_;j++)
        {
          resnew += *resmuptr**resmuptr;
          resinf = max(resinf,abs(*resmuptr));
          resmuptr++;
        }

        resnew += reszet*reszet;
        resinf = max(resinf,abs(reszet));

        resptr = res;
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
}


void OPTI::GCMMA::Update(
    bool& innerIterConverged
)
{
  innerIterConverged = InnerConvergence();

  if (innerIterConverged)
  {
    inner_iter_ = 0;
  }

}



