/*!---------------------------------------------------------------------------
\file drt_utils_maxent_basisfunctions..cpp

<pre>
Maintainer: Keijo Nissen
            nissen@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15253
</pre>

*---------------------------------------------------------------------------*/

#include <iomanip>
#include "drt_utils_maxent_basisfunctions.H"
#include "../linalg/linalg_serialdensematrix.H"
#include "../linalg/linalg_serialdensevector.H"

/*==========================================================================*\
 * Maximum-Entropy problem class                                            *
\*==========================================================================*/

/*--------------------------------------------------------------------------*
 | ctor MaxentProblem                                             nis Mar12 |
 *--------------------------------------------------------------------------*/
DRT::MESHFREE::MaxEntApprox::MaxEntApprox(Teuchos::ParameterList const & params)
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

  // Tolerance at which prior considered numerically zero - not "const" for
  // non-convex domains
  double rangeTol = params.get<double>("T_RANGE_TOL"); //
  SetRange(rangeTol);

  // determine type of compliance condition
  cmpl_  = DRT::INPUT::IntegralValue<INPAR::MESHFREE::compltype>(params,"T_COMPLIANCE");

  // determine parameters for compliance condition
//  cmpl_param_ = params.get<Teuchos::ParameterList*>("compltype_param");

  // determine if partition-of-unity is applied
  const bool pu         = DRT::INPUT::IntegralValue<int>(params,"PARTITION_OF_UNITY");

  // negativity weight of basis function; 0=non-negative
  const double neg      = params.get<double>("NEGATIVITY");

  // create dual problems of this type
  //   note: lower dimensional problems needed for evaluation on boundary
  dual1d_ = Teuchos::rcp(new DualProblem<1>(pu,neg,this));
  dual2d_ = Teuchos::rcp(new DualProblem<2>(pu,neg,this));
  dual3d_ = Teuchos::rcp(new DualProblem<3>(pu,neg,this));
};

/*--------------------------------------------------------------------------*
 | set range_ of basis function according to prior types          nis Mar12 |
 *--------------------------------------------------------------------------*/
void DRT::MESHFREE::MaxEntApprox::SetRange(double rangeTol)
{
  // for gaussian priors
  if (prior_==INPAR::MESHFREE::p_gauss) {
    // for symmetric gaussian priors
    if (skew_==INPAR::MESHFREE::p_sym) {
      range_ = sqrt(- var_ * log(rangeTol));
    } // end if prior symmetry
    // no assymetric priors implemented
    else {
      dserror("No other than symmetric max-ent prior implemented, yet.");
    }
  } // end if prior type
  // no other prior types implemented
  else {
    dserror("No other than Gaussian max-ent prior implemented, yet.");
  }

  //cout << "Range of basis functions: " << range_ << endl;

  return;
}

/*==========================================================================*\
 * Dual problem member class                                                *
\*==========================================================================*/

/*--------------------------------------------------------------------------*
 | ctor DualProblem                                               nis Mar12 |
 *--------------------------------------------------------------------------*/
template<int dim>
DRT::MESHFREE::MaxEntApprox::DualProblem<dim>::DualProblem(
  bool const    pu,
  double const  neg,
  MaxEntApprox* that
)
{
  // determine type of gap function/projection if necessary
  // TODO: gap function
  // if ((!neg) and (dim>1))
  //  enum gaptype gap    = params->get<enum gaptype>("gaptype");
  // assing correct functions to function pointers according to

  // make this-pointer of MaxEntApprox (that) known to DualProblem
  this_ = that;

  // set funtion pointers once
  if (pu) {
    if (!neg) {
      dualprob_ = Teuchos::rcp(new DualStandard);
    }
    else {
      dserror("No other than non-negative basis functions implemented.");
    }
  }
  else {
    dserror("No other than reduced problem with partition-of-unity for basis functions implemented.");
  }

  return;
}

/*--------------------------------------------------------------------------*
 | evaluating maxent basis functions and derivatives              nis Jan12 |
 *--------------------------------------------------------------------------*/

template<int dim>
int DRT::MESHFREE::MaxEntApprox::DualProblem<dim>::maxent_basisfunction(
  LINALG::SerialDenseVector &           funct , // basis functions values
  LINALG::SerialDenseMatrix &           deriv , // spatial derivatives of basis functions
  LINALG::SerialDenseMatrix const &     diffx   // distance vector between node and integration point
  )
{
  // get number of nodes
  const int  na  = diffx.N();
  const bool der = (deriv.M()!=0);

  // check some dimensions
  if (diffx.M()!=dim) dserror("Number of rows of diffx is %i but must be  %i (dim)!",diffx.M(),dim);
  if (funct.Length()!=na) funct.LightSize(na);
  if (der and (deriv.M()!=dim or deriv.N()!=na)) deriv.LightShape(dim,na);

  // initilize vector/matrix for prior functions
  LINALG::SerialDenseVector q(na);
  LINALG::SerialDenseMatrix dxq(0,0);
  if (der) {
    dxq.Reshape(dim,na);
  }
  // get prior functions and spatial derivatives
  SetPriorFunctDeriv(q,dxq,diffx);

  // initilize matrices for consistency/compliance conditions
  LINALG::SerialDenseMatrix c(dim,na);
  std::vector<LINALG::SerialDenseMatrix> dxc(0);
  if (der) {
    dxc.resize(dim);
    for (size_t i=0; i<dim; i++) {
      dxc[i].Reshape(dim,na);
    }
  }
  // get consistency/compliance conditions and spatial derivatives
  SetComplCondFunctDeriv(c,dxc,diffx);

  // perform Newton-Raphson optimisation
  LINALG::Matrix<dim,1>   lam(true);   // initialisation needed
  LINALG::Matrix<dim,1>   r(false);    // no initialisation needed
  LINALG::Matrix<dim,dim> Jinv(false); // no initialisation needed

  for(int i=0;i<this_->NewtonMax_;i++) {
    // get parameters of dual problem
    dualprob_->UpdateParams(funct,r,Jinv,na,q,c,lam);
    // perform Newton step
    lam.MultiplyNN(1.0,Jinv,r,1.0);

    // check if converged
    if (r.Norm2()<this_->NewtonTol_) {
      break;
    }
  }
  // update dual problem parameters at new lambda
  dualprob_->UpdateParams(funct,r,Jinv,na,q,c,lam);

  // NaN-check
  if (isnan(funct.Norm2()) or isnan(r.Norm2())) {
    dserror("NaN in MaxEnt-Optimization.");
  }

  // check if Newton-Raphson optimisation really converged
  int err = 0;
  if (r.Norm2()>this_->NewtonTol_) {
    dserror("No convergence in MaxEnt-Optimization.");
  }
//    err = 1;

  // compute spatial derivatives of basis function if required
  if (der) {
    dualprob_->GetDerivs(funct,deriv,lam,Jinv,q,dxq,c,dxc);
  }

  // return error flag for Newton-Raphson optimisation - bit useless though...
  return err;
}

/*--------------------------------------------------------------------------*
 | function to compute the prior function and its derivatives.    nis Jan12 |
 *--------------------------------------------------------------------------*/
template<int dim>
void DRT::MESHFREE::MaxEntApprox::DualProblem<dim>::SetPriorFunctDeriv(
  LINALG::SerialDenseVector       & q    ,
  LINALG::SerialDenseMatrix       & dxq  , // scaled by 1/q
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
      if (dxq.M()) {
        dxq.Update(2/(this_->var_),diffx,0.);
      }
    } // end if prior symmetry
    // no assymetric priors implemented
    else {
      dserror("No other than symmetric max-ent prior implemented, yet.");
    }
  } // end if prior type
  // no other prior types implemented
  else {
    dserror("No other than Gaussian max-ent prior implemented, yet.");
  }

  return;
}

/*--------------------------------------------------------------------------*
 | function to compute the compliance condition and derivatives.  nis Jan12 |
 *--------------------------------------------------------------------------*/
template<int dim>
void DRT::MESHFREE::MaxEntApprox::DualProblem<dim>::SetComplCondFunctDeriv(
  LINALG::SerialDenseMatrix              & c    , // (dim x numnodes)-matrix
  std::vector<LINALG::SerialDenseMatrix> & dxc  , // dim-vector of (dim x numnodes)-matrix
  const LINALG::SerialDenseMatrix        & diffx  // (dim x numnodes)-matrix
  ) const
{
  // linear consistency
  if (this_->cmpl_==INPAR::MESHFREE::c_linear) {
    c = diffx;
    if (dxc.size()) {
      for (size_t i=0;i<dim;i++) {
        for (int j=0;j<dxc[0].N();j++) {
          dxc[i](i,j) = -1;
        }
      }
    }
  }
  // no other comliance conditions implemented
  else {
    dserror("No other than linear consistency condition implemented, yet.");
  }

  return;
}

/*==========================================================================*/
//! functions of specific DualProblemTypes: standard dual problem (Arroyo2006)
/*==========================================================================*/

/*--------------------------------------------------------------------------*
 | dual problem parameters for nonnegative basis functions with             |
 | partition-of-unity constraint                         (public) nis Jan12 |
 *--------------------------------------------------------------------------*/
template<int dim>
void DRT::MESHFREE::MaxEntApprox::DualProblem<dim>::DualStandard::UpdateParams(
  LINALG::SerialDenseVector       & funct,
  LINALG::Matrix<dim,1>           & r    ,
  LINALG::Matrix<dim,dim>         & J    ,
  const int                         na   , // number of nodes
  const LINALG::SerialDenseVector & q    ,
  const LINALG::SerialDenseMatrix & c    ,
  const LINALG::Matrix<dim,1>     & lam    // unknown Lagrange multiplier
  )
{
  /*------------------------------------------------------------------------*
   * initialization
   *------------------------------------------------------------------------*/

  // temporary double for sum of basis functions.
  double fsum=0;
  // initialization of r and J
  r.Clear(); // since value is used in sum
  J.Clear(); // since value is used in sum/multiplication

  /*------------------------------------------------------------------------*
   * for-loop over all nodes
   *------------------------------------------------------------------------*/
  for (int i=0; i<na; i++) {
    // temporary vector for compliance of node a
    const LINALG::Matrix<dim,1> ca(const_cast<double*>(c[i]),true);
    // functi = q_i exp(lam_j c_ji) [/ (q_k exp(-lam_l c_lk))]
    funct[i] = q[i] * exp(-lam.Dot(ca));
    // r_i  = funct_j c_ij
    r.Update(funct[i], ca, 1.);
    // J_ij = funct_k c_ik c_jk [- r_i r_j]
    J.MultiplyNT(funct[i],ca,ca,1.);
    // fsum = q_k exp(-lam_l c_lk)
    fsum += funct[i];
  }

  /*------------------------------------------------------------------------*
   * final scaling
   *------------------------------------------------------------------------*/
  funct.Scale(1/fsum);           // ["/ (q_k exp(-lam_l c_lk))"]
  r.Scale(1/fsum);               // ["/ (q_k exp(-lam_l c_lk))"]
  J.MultiplyNT(-1.0,r,r,1/fsum); // ["/ (q_k exp(-lam_l c_lk)) - r_i r_j"]

  J.Invert();
};

/*--------------------------------------------------------------------------*
 | spatial derivatives for nonnegative basis functions with                 |
 | partition-of-unity constraint                         (public) nis Jan12 |
 *--------------------------------------------------------------------------*/
template<int dim>
void DRT::MESHFREE::MaxEntApprox::DualProblem<dim>::DualStandard::GetDerivs(
  const LINALG::SerialDenseVector              & funct, // basis functions values
        LINALG::SerialDenseMatrix              & deriv, // spatial derivatives of basis functions
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
   *  d(b_a)/dx = ð(b_a)/ðx + ð(b_a)/ðl · d(l)/dx
   *            = ð(b_a)/ðx - ð(b_a)/ðl · (J)^-1 · d(r)/dx
   *
   *  with functions:
   *  ---------------
   *  ð(b_a)/ðx^i = b_a (f_{a,x^i} - (b_b f_{b,x^i}))
   *  ð(b_a)/ðl^k = b_a c_a^k
   *  d(r^j)/dx^i = c_b^j ð(b_b)/ðx^i    + b_b c_{b,x^i}^j
   *             (= c_b^j b_b f_{b,x^i}  + b_b c_{b,x^i}^j)
   *
   *  expanded and with classified terms of implementation:                   //          I
   *  ----------------------------------------------------                    //  ,-------^-------.
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
    temp.Multiply('N','N',1.,dxc[i],funct,0.);       // temp = IIc2_j = c_{a,i}^j * b_a
    for(size_t j=0; j<dim; j++) {
      tempf_IIc(i,j) = tempf(j,0);                   // IIc = IIc2 [ + IIc1 still missing]
    }
  }

  temp.Multiply('N','N',1.0,f_ai,funct,0.0);         // temp = Ia_i = f_{a,i} * b_a

  for(size_t j=0; j<na; j++) {
    for(size_t i=0; i<dim; i++) {
      temp_IIa(i,j) *= funct[j];                    // temp_IIa   = b_a*c_a^k
      deriv(i,j) = funct[j]*(f_ai(i,j)-temp[i]);    // deriv = I = (f_ai - Ia) [ + ... still missing]
    }
  }

  temp_IIc.Multiply('N','T',1.0,c,deriv,1.0);       // IIc = IIc1 + IIc2

  tempf_IIb.MultiplyNN(-1.,Jinv,tempf_IIc);         // IIb = Jinv * IIc

  deriv.Multiply('N','N',1.,temp_IIb,temp_IIa,1.);  // deriv = I + IIa * IIb

}
