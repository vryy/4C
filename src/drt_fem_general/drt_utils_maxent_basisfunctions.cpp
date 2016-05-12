/*!---------------------------------------------------------------------------
\file drt_utils_maxent_basisfunctions.cpp

\brief Implementation of maxent basis functions

<pre>
\level 2
\maintainer Keijo Nissen
            http://www.lnm.mw.tum.de
            089 - 289-15253
</pre>

*---------------------------------------------------------------------------*/

#include <iomanip>
#include <string>
#include "drt_utils_maxent_basisfunctions.H"
#include "../drt_meshfree_discret/drt_meshfree_utils.H"
#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_serialdensematrix.H"
#include "../linalg/linalg_serialdensevector.H"
#include "../drt_inpar/inpar_meshfree.H"

/*==========================================================================*\
 * Maximum-Entropy problem class                                            *
\*==========================================================================*/

/*--------------------------------------------------------------------------*
 | ctor MaxentProblem                                             nis Mar12 |
 *--------------------------------------------------------------------------*/
DRT::MESHFREE::MaxEntApprox::MaxEntApprox(Teuchos::ParameterList const & params, int type)
  : type_(type)
{
  //----------------------------------------------------------------------
  // initialize members
  //----------------------------------------------------------------------

  // dilation parameter given as variance of prior function
  var_   = params.get<double>("VARIANCE");

  // determine type of prior
  prior_ = DRT::INPUT::IntegralValue<INPAR::MESHFREE::priortype>(params,"PRIOR");

  // determine type of compliance condition
//  for (int i=0; i<2; ++i)
//    cmpl_[i] =
//
//    DRT::INPUT::IntegralValue<INPAR::MESHFREE::compltype>(params,(std::to_string(i)+"_COMPLIANCE"));
  if (type_==0)
    cmpl_ = DRT::INPUT::IntegralValue<INPAR::MESHFREE::compltype>(params,("S_COMPLIANCE"));
  else if(type_==1)
    cmpl_ = DRT::INPUT::IntegralValue<INPAR::MESHFREE::compltype>(params,("W_COMPLIANCE"));
  else
    dserror("The basis function type has to be either (0) for solution or (1) for weighting basis fucntions!");

  // tolerance at which Newton is considered to be converged
  newtontol_ = params.get<double>("NEWTON_TOL");

  // maximum number of Newton steps
  newtonmaxiter_ = params.get<double>("NEWTON_MAXITER");

  //----------------------------------------------------------------------
  // initialize base class members (need members var_, prior_, and skew_)
  //----------------------------------------------------------------------

  // set range at which prior considered numerically zero
  SetRange(params.get<double>("RANGE_TOL"));

  //----------------------------------------------------------------------
  // get further input variables
  //----------------------------------------------------------------------

  // determine if partition-of-unity is applied
  const bool pu         = DRT::INPUT::IntegralValue<int>(params,"PARTITION_OF_UNITY");

  // negativity weight of basis function; 0=non-negative
  const double neg      = params.get<double>("NEGATIVITY");

  //----------------------------------------------------------------------
  // create dual problems of this type
  //   note: lower dimensional problems needed for evaluation on boundary
  //----------------------------------------------------------------------

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
    range_ = sqrt(- var_ * log(rangeTol));
  } // end if prior type
  // no other prior types implemented
  else {
    dserror("No other than Gaussian max-ent prior implemented, yet.");
  }

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
  MaxEntApprox* maxent
)
{
  // determine type of gap function/projection if necessary
  // TODO: gap function
  // if ((!neg) and (dim>1))
  //  enum gaptype gap    = params->get<enum gaptype>("gaptype");
  // assing correct functions to function pointers according to

  maxent_ = maxent;

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
  Teuchos::RCP<LINALG::SerialDenseMatrix>  diffx, // distance vector between node and integration point
  Teuchos::RCP<LINALG::SerialDenseVector>& funct, // basis functions values
  Teuchos::RCP<LINALG::SerialDenseMatrix>& deriv, // spatial derivatives of basis functions
  Teuchos::RCP<const LINALG::SerialDenseVector> upsilon, //!< eigen value at cell points for inf-flux basis weighting functions
  Teuchos::RCP<LINALG::SerialDenseVector> sfunct, // basis solution functions values in case of Bubnov Galerkin
  Teuchos::RCP<LINALG::SerialDenseMatrix> sderiv // spatial derivatives of basis solution functions in case of Bubnov Galerkin
  )
{
  // get number of nodes
  const int  na  = diffx->N();

  // check and if necessary adapt dimension of diffx
  if (diffx->M()!=dim)
  {
    int num_reduced = ReduceCloudDimensionOnFaces(diffx);
    if (diffx->M()!=dim)
      dserror("Covariance matrix has %i zero eigenvalues but needs to have %i for dimension reduction to %i.\nSo far, boundaries have to be plane/straight",num_reduced,num_reduced+(diffx->M()-dim),dim);
  }

  // check if we actually need to compute basis functions
  if ((maxent_->type_==0) and (maxent_->cmpl_!=INPAR::MESHFREE::c_linear))
    dserror("Solution basis function should have linear consistency!\nIf you REALLY want to change this, remove implicite Petrov-Galerkin check in 'else' case.");
  else if ((maxent_->type_==1 and (maxent_->cmpl_==INPAR::MESHFREE::c_linear))
           and (sfunct!=Teuchos::null))
  {
    funct = sfunct;
    deriv = sderiv;
    return 0;
  }

  // check and if necessary adapt sizes of funct and deriv
  bool der = (deriv!=Teuchos::null);
  if (funct->Length()!=na) funct->LightSize(na);
  if (der and (deriv->M()!=dim or deriv->N()!=na)) deriv->LightShape(dim,na);

  // initilize vector/matrix for prior functions
  LINALG::SerialDenseVector q(na);
  LINALG::SerialDenseMatrix dxq(0,0);
  if (der)
    dxq.Reshape(dim,na);

  // get prior functions and spatial derivatives
  SetPriorFunctDeriv(q,dxq,*diffx);

  // initilize matrices for consistency/compliance conditions
  LINALG::SerialDenseMatrix c(dim,na);
  std::vector<LINALG::SerialDenseMatrix> dxc(0);
  if (der)
  {
    dxc.resize(dim);
    for (size_t i=0; i<dim; i++)
      dxc[i].Reshape(dim,na);
  }

  // get consistency/compliance conditions and spatial derivatives
  SetComplCondFunctDeriv(c,dxc,*diffx,maxent_->type_,upsilon);

  // perform Newton-Raphson optimisation
  LINALG::Matrix<dim,1>   lam(true);   // initialisation needed
  LINALG::Matrix<dim,1>   r(false);    // no initialisation needed
  LINALG::Matrix<dim,dim> Jinv(false); // no initialisation needed

  // loop index initialized outside to be usable in error message if necessary
  int i;
  for(i=0;i<maxent_->newtonmaxiter_;i++)
  {
    // get parameters of dual problem
    int err = dualprob_->UpdateParams(*funct,r,Jinv,na,q,c,lam);
    if (err>0)
      return err;

    // perform Newton step
    lam.MultiplyNN(-1.0,Jinv,r,1.0);

    // check if converged
    if (r.Norm2()<maxent_->newtontol_)
      break;
  }

  // update dual problem parameters at new lambda
  dualprob_->UpdateParams(*funct,r,Jinv,na,q,c,lam);

  // NaN-check
  if (std::isnan(funct->Norm2()) or std::isnan(r.Norm2())) {
    dserror("NaN in MaxEnt-Optimization.");
  }

  // check if Newton-Raphson optimisation really converged
  if (r.Norm2()>maxent_->newtontol_) {
    dserror("No convergence in MaxEnt-Optimization after %i of %i iterations. |r|= %.2e > %.2e",i,maxent_->newtonmaxiter_,r.Norm2(),maxent_->newtontol_);
  }

  // compute spatial derivatives of basis function if required
  if (der)
  {
    int err = dualprob_->GetDerivs(*funct,*deriv,lam,Jinv,q,dxq,c,dxc);
    if (err>0)
      return err;
  }

  // everything went fine
  return 0;
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
  if (maxent_->prior_==INPAR::MESHFREE::p_gauss) {
    for (int i=0; i<diffx.N(); i++){
      double temp = 0.0;
      for (int j=0; j<dim; j++){
        temp += diffx(j,i)*diffx(j,i);
      }
      q(i) = exp(-(1.0/maxent_->var_) * temp);
    }
    if (dxq.M()) {
      dxq.Update(2.0/(maxent_->var_),diffx,0.0);
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
  LINALG::SerialDenseMatrix                   & c    , // (dim x numnodes)-matrix
  std::vector<LINALG::SerialDenseMatrix>      & dxc  , // dim-vector of (dim x numnodes)-matrix
  const LINALG::SerialDenseMatrix             & diffx, // (dim x numnodes)-matrix
  const int                                     type,
  Teuchos::RCP<const LINALG::SerialDenseVector> ups    // eigenvector at evaluation point for inf-flux weighting basis functions
  ) const
{
  // auxiliary variable: number of nodes
  int na = c.N();
  // auxiliary variable: 2-norm of upsilon
  double upsnorm = 0.0;
  if (ups!=Teuchos::null)
    upsnorm = ups->Norm2();

  // in case of linear consistency:
  // ==============================
  if (maxent_->cmpl_==INPAR::MESHFREE::c_linear or upsnorm<1e-14)
  {
    // linear consistency terms correspond to distance matrix
    c = diffx;

    // compute spatial derivatives of linear consistency terms if necessary
    if (dxc.size())
    {
      for (size_t i=0;i<dim;i++)
        for (int j=0;j<na;j++)
          dxc[i](i,j) = -1.0;
    }
  }
  // in case of streamline based convective consistency:
  // ===================================================
  else if (maxent_->cmpl_==INPAR::MESHFREE::c_stream)
  {
    // check whether eigenvector at evaluation point has already been evaluated
    if (ups==Teuchos::null) dserror("First compute eigenvectors at evaluation point by means of solution basis function before evaluation inf-flux weighting basis functions!");

    // clear consistency matrix
    c.Zero();

    // compute streamline convective consistency term
    LINALG::SerialDenseVector temp(na,false);
    temp.Multiply('T','N', 1.0, diffx, *ups, 0.0);
    for (int j=0; j<na; ++j)
      c(0,j) = (exp(temp(j))-1); // (exp(ups.Norm2()*10.0/2.0 - 1));

    // compute spatial derivatives of streamline convective consistency term if necessary
    if (dxc.size())
    {
      LINALG::SerialDenseMatrix& dxc0 = dxc[0];
      for (int j=0; j<na; ++j)
        for (int k=0; k<dim; ++k)
          dxc0(k,j) = - exp(temp(j)) * (*ups)(k);///(exp(ups->Norm2()*10.0/2.0 - 1));
    }

    // switch over dimension for linear consistency perpendicular to convection
    switch (dim)
    {
    case 2:
    {
      // get direction perpendicular to upsilon
      LINALG::Matrix<dim,1> perp(false);
      perp(0,0) =  (*ups)(1)/upsnorm;
      perp(1,0) = -(*ups)(0)/upsnorm;

      // compute linear consistency terms perpendicular to upsilon
      for (int j=0; j<na; ++j)
        c(1,j) = perp(0,0)*diffx(0,j) + perp(1,0)*diffx(1,j);

      // compute spatial derivatives of linear consistency terms perpendicular to upsilon if necessary
      if (dxc.size())
      {
        LINALG::SerialDenseMatrix& dxc1 = dxc[1];
        for (int j=0; j<na; ++j)
          for (int k=0; k<dim; ++k)
            dxc1(k,j) = - perp(k,0);
      }

      break;
    }
    case 3:
    {
      // initialization
      LINALG::Matrix<dim,2> perp(false);

      // note:
      // ======
      // To enforce linear consistency perpendicular to upsilon, we need to
      // impose linear consistency in two linearly independant directions in
      // the plane normal to upsilon. It does not matter which two. However,
      // orthogonality is assumed to improve robustness. Thus, our strategy is
      // to find any direction perpendicular to upsilon and then get the
      // second by cross product of upsilon and the first.

      // get first direction perpendicular to upsilon
      if (std::abs((*ups)(2)/upsnorm-1)>1e-3)
      {
        // for upsilon not too z-dominant (1e-3 is arbitrary value):
        const double scale = std::sqrt((*ups)(0)*(*ups)(0) + (*ups)(1)*(*ups)(1));
        perp(0,0) =  (*ups)(1)/scale;
        perp(1,0) = -(*ups)(0)/scale;
        perp(2,0) = 0.0;
      }
      else
      {
        // for z-dominant upsilon:
        const double scale = std::sqrt((*ups)(1)*(*ups)(1) + (*ups)(2)*(*ups)(2));
        perp(0,0) = 0.0;
        perp(1,0) =  (*ups)(2)/scale;
        perp(2,0) = -(*ups)(1)/scale;
      }

      // get second direction by cross-product
      perp(0,1) = (*ups)(1)/upsnorm*perp(2,0) - (*ups)(2)/upsnorm*perp(1,0);
      perp(1,1) = (*ups)(2)/upsnorm*perp(0,0) - (*ups)(0)/upsnorm*perp(2,0);
      perp(2,1) = (*ups)(0)/upsnorm*perp(1,0) - (*ups)(1)/upsnorm*perp(0,0);

      // compute linear consistency terms perpendicular to upsilon
      for (int i=1; i<dim; ++i)
        for (int j=0; j<na; ++j)
          c(i,j) = perp(0,i-1)*diffx(0,j) + perp(1,i-1)*diffx(1,j) + perp(2,i-1)*diffx(2,j);

      // compute spatial derivatives of linear consistency terms perpendicular to upsilon if necessary
      if (dxc.size())
      {
        for (int i=1; i<dim; ++i)
        {
          LINALG::SerialDenseMatrix& dxci = dxc[i];
          for (int j=0; j<na; ++j)
            for (int k=0; k<dim; ++k)
              dxci(k,j) = - perp(k,i-1);
        }
      }

      break;
    }
    default:
      dserror("Unknown dimension for streamline convective consistency.");
    }

  }
  // no other compliance conditions implemented
  else {
    dserror("Unknown consistency condition for MaxEnt dual problem.");
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
int DRT::MESHFREE::MaxEntApprox::DualProblem<dim>::DualStandard::UpdateParams(
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
  double fsum = 0.0;
  // initialization of r and J
  r.Clear(); // since value is used in sum
  J.Clear(); // since value is used in sum/multiplication

  /*------------------------------------------------------------------------*
   * for-loop over all nodes
   *------------------------------------------------------------------------*/
  for (int i=0; i<na; i++)
  {
    // temporary vector for compliance of node a
    const LINALG::Matrix<dim,1> ca(const_cast<double*>(c[i]),true);
    // functi = q_i exp(lam_j c_ji) ["/ (q_k exp(lam_l c_lk))" will follow]
    funct[i] = q[i] * exp(lam.Dot(ca));
    // r_i  = funct_j c_ij
    r.Update(funct[i], ca, 1.0);
    // J_ij = funct_k c_ik c_jk ["- r_i r_j" will follow]
    J.MultiplyNT(funct[i],ca,ca,1.0);
    // fsum = q_k exp(-lam_l c_lk)
    fsum += funct[i];
  }

  /*------------------------------------------------------------------------*
   * final scaling
   *------------------------------------------------------------------------*/

  funct.Scale(1.0/fsum);           // ["/ (q_k exp(lam_l c_lk))"]
  r.Scale(1.0/fsum);               // ["/ (q_k exp(lam_l c_lk))"]
  J.MultiplyNT(-1.0,r,r,1.0/fsum); // ["/ (q_k exp(lam_l c_lk)) - r_i r_j"]

  // return with error code 1 if determinant of J is too close to zero
  if (std::abs(J.Determinant())<1e-25)
    return 1;

  J.Invert();

  // everything went fine
  return 0;
};

/*--------------------------------------------------------------------------*
 | spatial derivatives for nonnegative basis functions with                 |
 | partition-of-unity constraint                         (public) nis Jan12 |
 *--------------------------------------------------------------------------*/
template<int dim>
int DRT::MESHFREE::MaxEntApprox::DualProblem<dim>::DualStandard::GetDerivs(
  const LINALG::SerialDenseVector              & funct, // basis functions values
        LINALG::SerialDenseMatrix              & deriv, // spatial derivatives of basis functions
  const LINALG::Matrix<dim,1>                  & lam ,  // argmax of dual problem
  const LINALG::Matrix<dim,dim>                & Jinv,  // inverted Hessian of dual function at maximum
  const LINALG::SerialDenseVector              & q   ,  // prior functions
  const LINALG::SerialDenseMatrix              & dxq ,  // spatial derivatives of prior functions
  const LINALG::SerialDenseMatrix              & c   ,  // constraints
  const std::vector<LINALG::SerialDenseMatrix> & dxc    // spatial derivatives of constraints
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
    f_ai.Update(lam(i,0), dxc[i], 1.0);              // f_ai = q_{a,i}/q_a + lam_j c_{a,i}^j
    temp.Multiply('N','N',1.0,dxc[i],funct,0.0);     // temp = IIc2_j = c_{a,i}^j * b_a
    for(size_t j=0; j<dim; j++) {
      tempf_IIc(i,j) = tempf(j,0);                   // IIc = IIc2 [ + IIc1 still missing]
    }
  }

  temp.Multiply('N','N',1.0,f_ai,funct,0.0);         // temp = Ia_i = f_{a,i} * b_a

  for(size_t j=0; j<na; j++) {
    for(size_t i=0; i<dim; i++) {
      temp_IIa(i,j) *= funct[j];                     // temp_IIa   = b_a*c_a^k
      deriv(i,j) = funct[j]*(f_ai(i,j)-temp[i]);     // deriv = I = (f_ai - Ia) [ + ... still missing]
    }
  }

  temp_IIc.Multiply('N','T',1.0,c,deriv,1.0);        // IIc = IIc1 + IIc2

  tempf_IIb.MultiplyNN(-1.0,Jinv,tempf_IIc);         // IIb = - Jinv * IIc

  deriv.Multiply('T','N',1.0,temp_IIb,temp_IIa,1.0); // (deriv = I + IIa * IIb)^T

  return 0;
}
