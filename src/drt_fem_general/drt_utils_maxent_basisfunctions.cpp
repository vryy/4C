/*!----------------------------------------------------------------------
\file drt_utils_maxent_basisfunctions..cpp

<pre>
Maintainer: Keijo Nissen
            nissen@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15253
</pre>

*----------------------------------------------------------------------*/

#include <iomanip>
#include "drt_utils_maxent_basisfunctions.H"
#include "../linalg/linalg_serialdensematrix.H"
#include "../linalg/linalg_serialdensevector.H"

/*============================================================================*\
 * Maximum-Entropy problem class                                              *
\*============================================================================*/

/*--------------------------------------------------------------------------*
 | ctor MaxentProblem                                             nis Mar12 |
 *--------------------------------------------------------------------------*/
DRT::MESHFREE::MaxEntProblem::MaxEntProblem(Teuchos::ParameterList const & params)
{
  // determine type of prior
  prior_ = DRT::INPUT::IntegralValue<INPAR::MESHFREE::priortype>(params,"T_PRIOR");
  // determine skewness of prior
  skew_  = DRT::INPUT::IntegralValue<INPAR::MESHFREE::priorskew>(params,"T_SKEW");
  // dilation parameter given as variance of prior function
  var_   = params.get<double>("T_VARIANCE");
  // tolerance at which Newton is considered to be converged
  NewtonTol_ = params.get<double>("NEWTON_TOL");
  // maximum number of Newton steps
  NewtonMax_ = params.get<double>("NEWTON_MAX");
  // Tolerance at which prior considered numerically zero - not "const" for non-convex domains
  double rangeTol = params.get<double>("T_RANGE_TOL");
  SetRange(rangeTol);
  // determine type of compliance condition
  cmpl_  = DRT::INPUT::IntegralValue<INPAR::MESHFREE::compltype>(params,"T_COMPLIANCE");
  // determine parameters for compliance condition
//  cmpl_param_ = params.get<Teuchos::ParameterList*>("compltype_param");
  // determine if partition-of-unity is applied
  const bool pu         = DRT::INPUT::IntegralValue<int>(params,"PARTITION_OF_UNITY");
  // negativity weight of basis function; 0=non-negative
  const double neg      = params.get<double>("NEGATIVITY");

  dual1d_ = Teuchos::rcp(new DualProblem<1>(pu,neg,this));
  dual2d_ = Teuchos::rcp(new DualProblem<2>(pu,neg,this));
  dual3d_ = Teuchos::rcp(new DualProblem<3>(pu,neg,this));
};

/*--------------------------------------------------------------------------*
 | set range_ of basis function according to prior types          nis Mar12 |
 *--------------------------------------------------------------------------*/
void DRT::MESHFREE::MaxEntProblem::SetRange(double rangeTol)
{
  // for gaussian priors
  if (prior_==INPAR::MESHFREE::p_gauss) {
    // for symmetric gaussian priors
    if (skew_==INPAR::MESHFREE::p_sym) {
      range_ = sqrt(- var_ * log(rangeTol));
    } // end if prior symmetry
    // no assymetric priors implemented
    else
      dserror("No other than symmetric max-ent prior implemented, yet.");
  } // end if prior type
  // no other prior types implemented
  else
    dserror("No other than Gaussian max-ent prior implemented, yet.");

  //cout << "Range of basis functions: " << range_ << endl;

  return;
}

/*============================================================================*\
 * Dual problem member class                                                  *
\*============================================================================*/

/*--------------------------------------------------------------------------*
 | ctor DualProblem                                               nis Mar12 |
 *--------------------------------------------------------------------------*/
template<int dim>
DRT::MESHFREE::MaxEntProblem::DualProblem<dim>::DualProblem(
  bool const & pu, double const & neg, MaxEntProblem* that
)
{
  // determine type of gap function/projection if necessary
  // TODO: gap function
  // if ((!neg) and (dim>1))
  //  enum gaptype gap    = params->get<enum gaptype>("gaptype");
  // assing correct functions to function pointers according to

  // make this-pointer of MaxEntProblem (that) known to DualProblem
  this_ = that;

  // set funtion pointers once
  if (pu)
    if (!neg) {
      UpdateParams = &DualProblem::DualParamsNonNegPU;
      GetDerivs = &DualProblem::DerivsNonNegPU;
    }
    else
      dserror("No other than non-negative basis functions implemented.");
  else
    dserror("No other than reduced problem with partition-of-unity for basis functions implemented.");
  return;
}

/*--------------------------------------------------------------------------*
 | evaluating maxent basis functions and derivatives              nis Jan12 |
 *--------------------------------------------------------------------------*/

template<int dim>
int DRT::MESHFREE::MaxEntProblem::DualProblem<dim>::maxent_basisfunction(
  LINALG::SerialDenseVector &           funct_ , // basis functions values
  LINALG::SerialDenseMatrix &           deriv_ , // spatial derivatives of basis functions
  LINALG::SerialDenseMatrix const &     diffx    // distance vector between node and integration point
  )
{
  // get number of nodes
  const int  na  = diffx.N();
  const bool der = deriv_.M();

  // get prior functions
  LINALG::SerialDenseVector q(na);
  LINALG::SerialDenseMatrix dxq;
  if (der)
    dxq.Reshape(dim,na);
  SetPriorFunctDeriv(q,dxq,diffx);

  // get consistency conditions
  LINALG::SerialDenseMatrix c(dim,na);
  std::vector<LINALG::SerialDenseMatrix> dxc;
  if (der) {
    dxc.resize(dim);
    for (size_t i=0; i<dim; i++)
      dxc[i].Reshape(dim,na);
  }
  SetComplCondFunctDeriv(c,dxc,diffx);

  // perform Newton-Raphson optimisation
  LINALG::Matrix<dim,1>   lam(true);   // initialisation needed
  LINALG::Matrix<dim,1>   r(false);    // no initialisation needed
  LINALG::Matrix<dim,dim> Jinv(false); // no initialisation needed
  for(int i=0;i<this_->NewtonMax_;i++) {
    // get parameters of dual problem
    (this->*UpdateParams)(funct_,r,Jinv,na,q,c,lam);
    // perform Newton step
    lam.MultiplyNN(1.0,Jinv,r,1.0);
    if (r.Norm2()<this_->NewtonTol_)
      break;
  }
  (this->*UpdateParams)(funct_,r,Jinv,na,q,c,lam);

  if (isnan(funct_.Norm2()) or isnan(r.Norm2()))
    dserror("NaN in MaxEnt-Optimization.");

  // if Newton-Raphson optimisation not converged
  int err = 0;
  if (r.Norm2()>this_->NewtonTol_)
    dserror("No convergence in MaxEnt-Optimization.");
//    err = 1;

  // compute basis function derivatives if required
  if (der) {
    (this->*GetDerivs)(funct_,deriv_,lam,Jinv,q,dxq,c,dxc);
  }

  // return error flag for Newton-Raphson optimisation
  return err;
}

/*--------------------------------------------------------------------------*
 | function to compute the prior function and its derivatives.    nis Jan12 |
 *--------------------------------------------------------------------------*/
template<int dim>
void DRT::MESHFREE::MaxEntProblem::DualProblem<dim>::SetPriorFunctDeriv(
  LINALG::SerialDenseVector & q    ,
  LINALG::SerialDenseMatrix & dxq  , // scaled by 1/q
  LINALG::SerialDenseMatrix const & diffx
  ) const
{
  // for gaussian priors
  if (this_->prior_==INPAR::MESHFREE::p_gauss) {
    // for symmetric gaussian priors
    if (this_->skew_==INPAR::MESHFREE::p_sym) {
      for (int i=0; i<diffx.N(); i++){
        double temp = 0;
        for (int j=0; j<dim; j++){
          temp += diffx(j,i)*diffx(j,i);
        }
        q(i) = exp(-(1/this_->var_) * temp);
      }
      if (dxq.M())
        dxq.Update(2/(this_->var_),diffx,0.);
    } // end if prior symmetry
    // no assymetric priors implemented
    else
      dserror("No other than symmetric max-ent prior implemented, yet.");
  } // end if prior type
  // no other prior types implemented
  else
    dserror("No other than Gaussian max-ent prior implemented, yet.");

//  cout << "diffx = ";
//  for (int i=0; i<diffx.M(); i++){
//    for (int j=0; j<diffx.N(); j++)
//      cout << diffx(i,j) << " ";
//    cout << endl;
//  }
//  cout << endl;
//
//  cout << "q = ";
//  for (int i=0; i<diffx.N(); i++)
//    cout << q(i) << " ";
//  cout << endl;
//  cout << endl;
//
//
//  cout << "dxq = ";
//  for (int i=0; i<diffx.M(); i++){
//    for (int j=0; j<diffx.N(); j++)
//      cout << dxq(i,j) << " ";
//    cout << endl;
//  }
//  cout << endl;

  return;
}

/*--------------------------------------------------------------------------*
 | function to compute the compliance condition and derivatives.  nis Jan12 |
 *--------------------------------------------------------------------------*/
template<int dim>
void DRT::MESHFREE::MaxEntProblem::DualProblem<dim>::SetComplCondFunctDeriv(
  LINALG::SerialDenseMatrix              & c    , // (dim x numnodes)-matrix
  std::vector<LINALG::SerialDenseMatrix> & dxc  , // dim-vector of (dim x numnodes)-matrix
  const LINALG::SerialDenseMatrix        & diffx  // (dim x numnodes)-matrix
  ) const
{
  // linear consistency
  if (this_->cmpl_==INPAR::MESHFREE::c_linear) {
    c = diffx;
    if (dxc.size()) {
      for (size_t i=0;i<dim;i++)
        for (int j=0;j<dxc[0].N();j++)
          dxc[i](i,j) = -1;
    }
  }
  // no other comliance conditions implemented
  else
    dserror("No other than linear consistency condition implemented, yet.");

//  cout << "c = ";
//  for (int i=0; i<diffx.M(); i++){
//    for (int j=0; j<diffx.N(); j++)
//      cout << c(i,j) << " ";
//    cout << endl;
//  }
//  cout << endl;
//
//  for (int h=0; h<diffx.M(); h++){
//    cout << "dxc[" << h <<"] = ";
//    for (int i=0; i<diffx.M(); i++){
//      for (int j=0; j<diffx.N(); j++)
//        cout << dxc[h](i,j) << " ";
//      cout << endl;
//    }
//    cout << endl;
//  }

  return;
}

/*--------------------------------------------------------------------------*
 | dual problem parameters for nonnegative basis functions with             |
 | partition-of-unity constraint                         (public) nis Jan12 |
 *--------------------------------------------------------------------------*/
template<int dim>
void DRT::MESHFREE::MaxEntProblem::DualProblem<dim>::DualParamsNonNegPU(
  LINALG::SerialDenseVector       & funct_,
  LINALG::Matrix<dim,1>           & r     ,
  LINALG::Matrix<dim,dim>         & J     ,
  const int                       & na    ,
  const LINALG::SerialDenseVector & q     ,
  LINALG::SerialDenseMatrix & c     , // THIS SHOULD BE CONST!!! (Impossible
                                      // because of LINALG::Matrix-view on its
                                      // columns for ca)
  const LINALG::Matrix<dim,1>     & lam
  )
{
  // temporary double for sum of basis functions.
  double fsum=0;
  // initialization of r and J
  r.Clear(); // since value is used in sum
  J.Clear(); // since value is used in sum/multiplication
  for (int i=0; i<na; i++) {
    // temporary vector for compliance of node a
    const LINALG::Matrix<dim,1> ca(c[i],true);
    // funct_i = q_i exp(lam_j c_ji) [/ (q_k exp(-lam_l c_lk))]
    funct_[i] = q[i] * exp(-lam.Dot(ca));
    // r_i  = funct_j c_ij
    r.Update(funct_[i], ca, 1.);
    // J_ij = funct_k c_ik c_jk [- r_i r_j]
    J.MultiplyNT(funct_[i],ca,ca,1);
    // fsum = q_k exp(-lam_l c_lk)
    fsum += funct_[i];
  }
  funct_.Scale(1/fsum); // [/ (q_k exp(-lam_l c_lk))]
  r.Scale(1/fsum);      // [/ (q_k exp(-lam_l c_lk))]
  J.MultiplyNT(-1,r,r,1/fsum); // [/ (q_k exp(-lam_l c_lk)) - r_i r_j]
  J.Invert();
};

/*==========================================================================*/
//! struct of functions to compute spatial derivatives of basis function
/*==========================================================================*/

/*--------------------------------------------------------------------------*
 | spatial derivatives for nonnegative basis functions with                 |
 | partition-of-unity constraint                         (public) nis Jan12 |
 *--------------------------------------------------------------------------*/
template<int dim>
void DRT::MESHFREE::MaxEntProblem::DualProblem<dim>::DerivsNonNegPU(
  const LINALG::SerialDenseVector              & funct_, // basis functions values
        LINALG::SerialDenseMatrix              & deriv_, // spatial derivatives of basis functions
  const LINALG::Matrix<dim,1>                  & lam   , // argmax of dual problem
  const LINALG::Matrix<dim,dim>                & Jinv  , // inverted Hessian of dual function at maximum
  const LINALG::SerialDenseVector              & q     , // prior functions
  const LINALG::SerialDenseMatrix              & dxq   , // spatial derivatives of prior functions
  const LINALG::SerialDenseMatrix              & c     , // constraints
  const std::vector<LINALG::SerialDenseMatrix> & dxc     // spatial derivatives of constraints
  )
{
  /*
   *  for non-negative basis functions with partition-of-unity property:
   *  ==================================================================
   *
   *  d(b_a)/dx = ð(b_a)/ðx + ð(b_a)/dł · d(ł)/dx
   *            = ð(b_a)/ðx - ð(b_a)/dł · (J)^-1 · d(r)/dx
   *
   *  ð(b_a)/ðx = db_dx = b_a (f_{a,i} - (b_b f_{b,i}))
   *  ð(b_a)/dł = db_dl = b_a c_a^k
   *  d(r)/dx   = dr_dx = c_b^j (ð(b_b)/ðx)^i + b_b c_{b,i}^j
   *                   (= c_b^j  b_b f_{b,i}  + b_b c_{b,i}^j)
   *
   *                                                                          //          I
   *                                                                          //  ,-------^-------.
   *  b_{a,i} =   b_a (f_{a,i} - (b_b f_{b,i}))                               //  b_a * (f_ai - Ia)
   *                                                                          //                    IIc
   *                                                                          //               ,-----^-----.
   *            - (b_a c_a^k) (J_inv^kj) (b_b c_b^j f_{b,i} + b_b c_{b,i}^j)  // IIa * J_inv * (IIc1 + IIc2)
   *                                                                          //       `---------v---------´
   *                                                                          //                IIb
   */

  size_t na=q.Length();

  LINALG::SerialDenseMatrix f_ai(dxq);               // to be: q_{a,i}/q_a + lam^l*c_{a,i}^l
  LINALG::SerialDenseMatrix temp_IIa(c);             // to be: b_a*c_a^k
  LINALG::SerialDenseMatrix temp_IIb(dim,dim);       // to be: Jinv * IIc
  LINALG::SerialDenseMatrix temp_IIc(dim,dim);       // to be: IIc1 + IIc2
  LINALG::SerialDenseVector temp(dim);               // to be either: b_b * f_{a,i}
                                                     //           or: b_b * c_{a,(i)}^j [j as leading dimension]
  LINALG::Matrix<dim,dim> tempf_IIb(temp_IIb,true);  // fixedsize view on temp_IIb
  LINALG::Matrix<dim,dim> tempf_IIc(temp_IIc,true);  // fixedsize view on temp_IIc
  LINALG::Matrix<dim,1>   tempf(temp,true);          // fixedsize view on temp

  //TODO: improve enrolling at compile time
  for(size_t i=0; i<dim; i++) {
    f_ai.Update(lam(i,0), dxc[i], 1.);               // f_ai = q_{a,i}/q_a + lam_j c_{a,i}^j
    temp.Multiply('N','N',1.,dxc[i],funct_,0.);      // temp = IIc2_j = c_{a,i}^j * b_a
    for(size_t j=0; j<dim; j++)
      tempf_IIc(i,j) = tempf(j,0);                   // IIc = IIc2 [ + IIc1 still missing]
  }

  temp.Multiply('N','N',1.0,f_ai,funct_,0.0);        // temp = Ia_i = f_{a,i} * b_a

  for(size_t i=0; i<dim; i++) {
    for(size_t j=0; j<na; j++) {
      temp_IIa(i,j) *= funct_[j];                    // temp_IIa   = b_a*c_a^k
      deriv_(i,j) = funct_[j]*(f_ai(i,j)-temp[i]);   // deriv_ = I = (f_ai - Ia) [ + ... still missing]
    }
  }

  temp_IIc.Multiply('N','T',1.0,c,deriv_,1.0);       // IIc = IIc1 + IIc2

  tempf_IIb.MultiplyNN(-1.,Jinv,tempf_IIc);          // IIb = Jinv * IIc

//  cout << "lam = ";
//  for (int i=0; i<dim; i++){
//    cout << lam(i,0) << " ";
//  }
//  cout << endl;
//  cout << endl;
//
//  cout << "Jinv = ";
//  for (int i=0; i<dim; i++){
//    for (int j=0; j<dim; j++)
//      cout << Jinv(i,j) << " ";
//    cout << endl;
//  }
//  cout << endl;
//
//  cout << "f_ai = ";
//  for (int i=0; i<dim; i++){
//    for (int j=0; j<na; j++)
//      cout << f_ai(i,j) << " ";
//    cout << endl;
//  }
//  cout << endl;
//
//  cout << "dp_dx = ";
//  for (int i=0; i<dim; i++){
//    for (int j=0; j<na; j++)
//      cout << deriv_(i,j) << " ";
//    cout << endl;
//  }
//  cout << endl;
//
//  cout << "dp_dl = ";
//  for (int i=0; i<dim; i++){
//    for (int j=0; j<na; j++)
//      cout << temp_IIa(i,j) << " ";
//    cout << endl;
//  }
//  cout << endl;
//
//  cout << "dr_dx = ";
//  for (int i=0; i<dim; i++){
//    for (int j=0; j<dim; j++)
//      cout << temp_IIc(i,j) << " ";
//    cout << endl;
//  }
//  cout << endl;
//
//  cout << "dl_dx = ";
//  for (int i=0; i<dim; i++){
//    for (int j=0; j<dim; j++)
//      cout << temp_IIb(i,j) << " ";
//    cout << endl;
//  }
//  cout << endl;

  deriv_.Multiply('N','N',1.,temp_IIb,temp_IIa,1.); // deriv = I + IIa * IIb

}
