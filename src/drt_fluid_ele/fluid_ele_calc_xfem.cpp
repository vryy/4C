/*----------------------------------------------------------------------*/
/*!
\file fluid_ele_calc_xfem.cpp

\brief Internal implementation of XFluid element interface coupling

\level 2

<pre>
\maintainer Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15240
</pre>
*/
/*----------------------------------------------------------------------*/
#include <Teuchos_TimeMonitor.hpp>

#include <fstream>

#include "../drt_cut/cut_boundarycell.H"
#include "../drt_cut/cut_position.H"
#include "../drt_cut/cut_volumecell.H"

#include "../drt_geometry/position_array.H"

#include "../linalg/linalg_utils.H"

#include "fluid_ele.H"
#include "fluid_ele_parameter_xfem.H"
#include "fluid_ele_calc_xfem.H"

#include "../drt_inpar/inpar_xfem.H"
#include "../drt_inpar/inpar_cut.H"

#include "../drt_mat/newtonianfluid.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_function.H"

#include "../drt_lib/drt_condition_utils.H"

#include "../drt_xfem/xfem_condition_manager.H"

#include <Teuchos_TimeMonitor.hpp>

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::FluidEleCalcXFEM<distype> * DRT::ELEMENTS::FluidEleCalcXFEM<distype>::Instance( bool create )
{
  static FluidEleCalcXFEM<distype> * instance;
  if ( create )
  {
    if ( instance==NULL )
    {
      instance = new FluidEleCalcXFEM<distype>();
    }
  }
  else
  {
    if ( instance!=NULL )
      delete instance;
    instance = NULL;
  }
  return instance;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleCalcXFEM<distype>::Done()
{
  // delete this pointer! Afterwards we have to go! But since this is a
  // cleanup call, we can do it this way.
  Instance( false );
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::FluidEleCalcXFEM<distype>::FluidEleCalcXFEM()
  : DRT::ELEMENTS::FluidEleCalc<distype>::FluidEleCalc(),
    densaf_master_(0.0),
    densaf_slave_(0.0),
    viscaf_master_(0.0),
    viscaf_slave_(0.0),
    gamma_m_(0.0),
    gamma_s_(0.0),
    evelaf_(true),
    epreaf_(true),
    eveln_(true),
    epren_(true),
    ivelint_jump_(true),
    itraction_jump_(true),
    proj_tangential_(true),
    LB_proj_matrix_(true),
    ivelintn_jump_(true),
    itractionn_jump_(true),
    velint_s_(true),
    velintn_s_(true),
    rst_(true),
    normal_(true),
    x_side_(true),
    x_gp_lin_(true)
{
  // we use the standard parameter list here, since there are not any additional
  // xfem-specific parameters required in this derived class
  my::fldpara_=DRT::ELEMENTS::FluidEleParameterXFEM::Instance();
  fldparaxfem_ = static_cast<DRT::ELEMENTS::FluidEleParameterXFEM*>(my::fldpara_);
}


namespace DRT
{
namespace ELEMENTS
{

/*-------------------------------------------------------------------------------*
          Evaluate routine for cut elements of XFEM  (public)
*-------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int FluidEleCalcXFEM<distype>::EvaluateXFEM(DRT::ELEMENTS::Fluid*                             ele,
                                            DRT::Discretization &                             discretization,
                                            const std::vector<int> &                          lm,
                                            Teuchos::ParameterList&                           params,
                                            Teuchos::RCP<MAT::Material> &                     mat,
                                            Epetra_SerialDenseMatrix&                         elemat1_epetra,
                                            Epetra_SerialDenseMatrix&                         elemat2_epetra,
                                            Epetra_SerialDenseVector&                         elevec1_epetra,
                                            Epetra_SerialDenseVector&                         elevec2_epetra,
                                            Epetra_SerialDenseVector&                         elevec3_epetra,
                                            const std::vector<DRT::UTILS::GaussIntegration> & intpoints,
                                            const GEO::CUT::plain_volumecell_set &            cells,
                                            bool                                              offdiag )
{
  int err=0;

  for( std::vector<DRT::UTILS::GaussIntegration>::const_iterator i=intpoints.begin();i!=intpoints.end();++i )
  {
    const DRT::UTILS::GaussIntegration intpoints_cell = *i;
    err = my::Evaluate( ele, discretization, lm, params, mat,
                   elemat1_epetra, elemat2_epetra,
                   elevec1_epetra, elevec2_epetra, elevec3_epetra,
                   intpoints_cell, offdiag);
    if(err)
      return err;
  }
  return err;
}

/*-------------------------------------------------------------------------------*
          Evaluate routine for cut elements of XFEM  (public)
*-------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int FluidEleCalcXFEM<distype>::IntegrateShapeFunctionXFEM(
    DRT::ELEMENTS::Fluid*                             ele,
    DRT::Discretization &                             discretization,
    const std::vector<int> &                          lm,
    Epetra_SerialDenseVector&                         elevec1_epetra,
    const std::vector<DRT::UTILS::GaussIntegration> & intpoints,
    const GEO::CUT::plain_volumecell_set &            cells
)
{
  int err=0;

  for( std::vector<DRT::UTILS::GaussIntegration>::const_iterator i=intpoints.begin();i!=intpoints.end();++i )
  {
    const DRT::UTILS::GaussIntegration gint = *i;
    err = my::IntegrateShapeFunction( ele, discretization, lm,
                   elevec1_epetra,
                   gint);
    if(err)
      return err;
  }

  return err;
}


/// error computation
template <DRT::Element::DiscretizationType distype>
int FluidEleCalcXFEM<distype>::ComputeError(
                         DRT::ELEMENTS::Fluid*         ele,
                         Teuchos::ParameterList&       params,
                         Teuchos::RCP<MAT::Material>&  mat,
                         DRT::Discretization&          discretization,
                         std::vector<int>&             lm,
                         Epetra_SerialDenseVector&     ele_dom_norms)
{
  // integrations points and weights
  // more GP than usual due to (possible) cos/exp fcts in analytical solutions
  // degree 5
  const DRT::UTILS::GaussIntegration intpoints(distype, 8);
  return ComputeError( ele, params, mat,
                       discretization, lm,
                       ele_dom_norms, intpoints);
}

template <DRT::Element::DiscretizationType distype>
int FluidEleCalcXFEM<distype>::ComputeError(
                         DRT::ELEMENTS::Fluid*                ele,
                         Teuchos::ParameterList&              params,
                         Teuchos::RCP<MAT::Material>&         mat,
                         DRT::Discretization&                 discretization,
                         std::vector<int>&                    lm,
                         Epetra_SerialDenseVector&            ele_dom_norms, // squared element domain norms
                         const DRT::UTILS::GaussIntegration & intpoints)
{
  // analytical solution
  LINALG::Matrix<my::nsd_,1>         u_analyt(true);
  LINALG::Matrix<my::nsd_,my::nsd_>  grad_u_analyt(true);
  double p_analyt = 0.0;

  // error
  LINALG::Matrix<my::nsd_,1>        u_err(true);
  LINALG::Matrix<my::nsd_,my::nsd_> grad_u_err(true);
  double p_err = 0.0;

  const int calcerr = DRT::INPUT::get<INPAR::FLUID::CalcError>(params,"calculate error");
  const int calcerrfunctno = DRT::INPUT::get<int>(params,"error function number");

  const double t = my::fldparatimint_->Time();

  //set element id
  my::eid_ = ele->Id();


  //----------------------------------------------------------------------------
  //   Extract velocity/pressure from global vectors
  //----------------------------------------------------------------------------

  // fill the local element vector/matrix with the global values
  LINALG::Matrix<my::nsd_,my::nen_> evelaf(true);
  LINALG::Matrix<my::nen_,1>        epreaf(true);
  this->ExtractValuesFromGlobalVector(discretization,lm, *my::rotsymmpbc_, &evelaf, &epreaf, "u and p at time n+1 (converged)" );

  //----------------------------------------------------------------------------
  //                         ELEMENT GEOMETRY
  //----------------------------------------------------------------------------

  // get node coordinates
  GEO::fillInitialPositionArray<distype,my::nsd_, LINALG::Matrix<my::nsd_,my::nen_> >(ele,my::xyze_);

  //----------------------------------------------------------------
  // Now do the nurbs specific stuff (for isogeometric elements)
  //----------------------------------------------------------------
  if(my::isNurbs_)
  {
//    // access knots and weights for this element
//    bool zero_size = DRT::NURBS::GetMyNurbsKnotsAndWeights(discretization,ele,my::myknots_,my::weights_);
//
//    // if we have a zero sized element due to a interpolated point -> exit here
//    if(zero_size)
//      return(0);
    dserror("compute error not implemented for nurbs");
  } // Nurbs specific stuff

  if (ele->IsAle())
  {
    LINALG::Matrix<my::nsd_,my::nen_>       edispnp(true);
    this->ExtractValuesFromGlobalVector(discretization,lm, *my::rotsymmpbc_, &edispnp, NULL,"dispnp");

    // get new node positions for isale
    my::xyze_ += edispnp;
  }

//------------------------------------------------------------------
//                       INTEGRATION LOOP
//------------------------------------------------------------------

  for ( DRT::UTILS::GaussIntegration::iterator iquad=intpoints.begin(); iquad!=intpoints.end(); ++iquad )
  {
    // evaluate shape functions and derivatives at integration point
    my::EvalShapeFuncAndDerivsAtIntPoint(iquad.Point(),iquad.Weight());

    // get velocity at integration point
    // (values at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
    my::velint_.Multiply(evelaf,my::funct_);

    // get velocity derivatives at integration point
    // (values at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
    my::vderxy_.MultiplyNT(evelaf,my::derxy_);

    // get pressure at integration point
    // (value at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
    double preint = my::funct_.Dot(epreaf);

    // get coordinates at integration point
    LINALG::Matrix<my::nsd_,1> xyzint(true);
    xyzint.Multiply(my::xyze_,my::funct_);

    // get viscosity
    if (mat->MaterialType() == INPAR::MAT::m_fluid)
    {
      const MAT::NewtonianFluid* actmat = static_cast<const MAT::NewtonianFluid*>(mat.get());

      // get constant kinematic viscosity
      my::visc_ = actmat->Viscosity()/actmat->Density();
    }
    else dserror("Material is not Newtonian Fluid");

    AnalyticalReference(
             calcerr,          ///< which reference solution
             calcerrfunctno,   ///< error function number
             u_analyt,         ///< exact velocity
             grad_u_analyt,    ///< exact velocity gradient
             p_analyt,         ///< exact pressure
             xyzint,           ///< xyz position of gaussian point
             t,
             mat
   );

    // compute difference between analytical solution and numerical solution
    p_err = preint - p_analyt;
    u_err.Update(1.0, my::velint_, -1.0, u_analyt, 0.0);
    grad_u_err.Update(1.0, my::vderxy_, -1.0, grad_u_analyt, 0.0);

    // error on pre-defined functional
    // here G=sin(x)( u,x - u,x exact )
    double funcerr = ( sin(xyzint(0,0)) * ( grad_u_analyt(0,0)-my::vderxy_(0,0) ) )*my::fac_;

    // standard domain errors
    // 1.   || u - u_h ||_L2(Omega)              =   standard L2-norm for velocity
    // 2.   || grad( u - u_h ) ||_L2(Omega)      =   standard H1-seminorm for velocity
    // 3.   || u - u_h ||_H1(Omega)              =   standard H1-norm for velocity
    //                                           =   sqrt( || u - u_h ||^2_L2(Omega) + || grad( u - u_h ) ||^2_L2(Omega) )
    // 4.   || p - p_h ||_L2(Omega)              =   standard L2-norm for pressure
    //
    // viscosity-scaled domain errors
    // 5.   || nu^(+1/2) grad( u - u_h ) ||_L2(Omega)      =   visc-scaled H1-seminorm for velocity
    //                                                     =   nu^(+1/2) * || grad( u - u_h ) ||_L2(Omega) (for homogeneous visc)
    // 6.   || nu^(-1/2) (p - p_h) ||_L2(Omega)            =   visc-scaled L2-norm for pressure
    //                                                     =   nu^(-1/2) * || p - p_h ||_L2(Omega) (for homogeneous visc)
    // 7.   || sigma^(+1/2) ( u - u_h ) ||_L2(Omega)       =   sigma-scaled L2-norm for velocity
    //                                                     =   sigma^(+1/2) * || u - u_h ||_L2(Omega) (for homogeneous sigma)
    // 8.   || Phi^(+1/2) (p - p_h) ||_L2(Omega)           =   Phi-scaled L2-norm for pressure
    //                                                     =   Phi^(+1/2) * || p - p_h ||_L2(Omega) (for homogeneous Phi)
    // with Phi^{-1} = sigma*CP^2 + |beta|*CP + nu + (|beta|*CP/sqrt(sigma*CP^2 + nu))^2, see Massing,Schott,Wall Oseen paper
    //
    // 9. functional G=sin(x)( u,x - u,x exact ) (Sudhakar)


    double u_err_squared      = 0.0;
    double grad_u_err_squared = 0.0;
    double p_err_squared      = 0.0;

    // evaluate squared errors at gaussian point
    for (int isd=0;isd<my::nsd_;isd++)
    {
      u_err_squared += u_err(isd)*u_err(isd)*my::fac_;

      for(int jsd=0; jsd<my::nsd_;jsd++)
      {
        grad_u_err_squared += grad_u_err(isd,jsd)*grad_u_err(isd,jsd)*my::fac_;
      }

    }

    p_err_squared = p_err*p_err*my::fac_;

    double Poincare_const = 1.0; // scales as upper bound for mesh size
    double beta_maximum   = 1.0; // maximal advective velocity in domain
    double sigma = 1./my::fldparatimint_->TimeFac(); // sigma scaling in Oseen results from time discretization

    double Phi_tmp = beta_maximum * Poincare_const/sqrt(sigma * Poincare_const*Poincare_const + my::visc_);
    double Phi_inverse_squared = sigma * Poincare_const*Poincare_const + beta_maximum * Poincare_const + my::visc_ + Phi_tmp*Phi_tmp;

    // standard domain errors
    ele_dom_norms[0] += u_err_squared;
    ele_dom_norms[1] += grad_u_err_squared;
    ele_dom_norms[2] += u_err_squared + grad_u_err_squared;
    ele_dom_norms[3] += p_err_squared;

    // viscosity-scaled domain errors
    ele_dom_norms[4] += my::visc_ * grad_u_err_squared;
    ele_dom_norms[5] += 1.0/my::visc_ * p_err_squared;
    ele_dom_norms[6] += sigma * u_err_squared;
    ele_dom_norms[7] += 1.0/Phi_inverse_squared * p_err_squared;

    // error for predefined functional
    ele_dom_norms[8] += funcerr;


  } // loop gaussian points

  return 0;
}

template <DRT::Element::DiscretizationType distype>
void FluidEleCalcXFEM<distype>::AnalyticalReference(
    const int                               calcerr,     ///< which reference solution
    const int                               calcerrfunctno, ///< error function number
    LINALG::Matrix<my::nsd_,1> &            u,           ///< exact jump vector (coupled)
    LINALG::Matrix<my::nsd_,my::nsd_> &     grad_u,      ///< exact velocity gradient
    double &                                p,           ///< exact pressure
    LINALG::Matrix<my::nsd_,1> &            xyzint,      ///< xyz position of gaussian point
    const double &                          t,           ///< time
    Teuchos::RCP<MAT::Material>             mat
)
{
  // Compute analytical solution
  switch(calcerr)
  {
  case INPAR::FLUID::beltrami_stat_stokes:
  case INPAR::FLUID::beltrami_stat_navier_stokes:
  case INPAR::FLUID::beltrami_instat_stokes:
  case INPAR::FLUID::beltrami_instat_navier_stokes:
  {

    // function evaluation requires a 3D position vector!!
    double position[3];

    if(my::nsd_ == 3)
    {
      position[0] = xyzint(0);
      position[1] = xyzint(1);
      position[2] = xyzint(2);
    }
    else dserror("invalid nsd %d", my::nsd_);

    // evaluate velocity and pressure
    Teuchos::RCP<DRT::UTILS::Function> function = Teuchos::null;

    // evaluate the velocity gradient
    Teuchos::RCP<DRT::UTILS::Function> function_grad = Teuchos::null;

    // evaluate velocity and pressure
    // evaluate the velocity gradient

    function      = Teuchos::rcp(new DRT::UTILS::BeltramiUP( mat));
    function_grad = Teuchos::rcp(new DRT::UTILS::BeltramiGradU( mat));

    if(my::nsd_==3)
    {
      u(0) = function->Evaluate(0,position,t,NULL);
      u(1) = function->Evaluate(1,position,t,NULL);
      u(2) = function->Evaluate(2,position,t,NULL);
      p = function->Evaluate(3,position,t,NULL);
    }
    else dserror("case 'kimmoin_stat' is a 3D specific case");


    if(my::nsd_==3)
    {
      grad_u(0,0) = function_grad->Evaluate(0,position,t,NULL); // u,x
      grad_u(0,1) = function_grad->Evaluate(1,position,t,NULL); // u,y
      grad_u(0,2) = function_grad->Evaluate(2,position,t,NULL); // u,z

      grad_u(1,0) = function_grad->Evaluate(3,position,t,NULL); // v,x
      grad_u(1,1) = function_grad->Evaluate(4,position,t,NULL); // v,y
      grad_u(1,2) = function_grad->Evaluate(5,position,t,NULL); // v,z

      grad_u(2,0) = function_grad->Evaluate(6,position,t,NULL); // w,x
      grad_u(2,1) = function_grad->Evaluate(7,position,t,NULL); // w,y
      grad_u(2,2) = function_grad->Evaluate(8,position,t,NULL); // w,z
    }
    else dserror("case 'kimmoin_stat' is a 3D specific case");
  }
  break;

  case INPAR::FLUID::beltrami_flow:
  {
    if (my::nsd_ == 3)
    {
      const double a      = M_PI/4.0;
      const double d      = M_PI/2.0;

      double x = xyzint(0);
      double y = xyzint(1);
      double z = xyzint(2);

      double visc = my::visc_;

      // calculation of velocities and pressure
      u(0) = -a * ( exp(a*x) * sin(a*y + d*z) + exp(a*z) * cos(a*x + d*y) ) * exp(-visc*d*d*t);
      u(1) = -a * ( exp(a*y) * sin(a*z + d*x) + exp(a*x) * cos(a*y + d*z) ) * exp(-visc*d*d*t);
      u(2) = -a * ( exp(a*z) * sin(a*x + d*y) + exp(a*y) * cos(a*z + d*x) ) * exp(-visc*d*d*t);
      p    = -a*a/2.0 *
              ( exp(2.0*a*x) + exp(2.0*a*y) + exp(2.0*a*z)
                + 2.0 * sin(a*x + d*y) * cos(a*z + d*x) * exp(a*(y+z))
                + 2.0 * sin(a*y + d*z) * cos(a*x + d*y) * exp(a*(z+x))
                + 2.0 * sin(a*z + d*x) * cos(a*y + d*z) * exp(a*(x+y))
              )* exp(-2.0*visc*d*d*t);

      // velocity gradients
      grad_u(0,0) = -a * ( a * exp(a*x) * sin(a*y + d*z) - a * exp(a*z) * sin(a*x + d*y) ) * exp(-visc*d*d*t); //u,x
      grad_u(0,1) = -a * ( a * exp(a*x) * cos(a*y + d*z) - d * exp(a*z) * sin(a*x + d*y) ) * exp(-visc*d*d*t); //u,y
      grad_u(0,2) = -a * ( d * exp(a*x) * cos(a*y + d*z) + a * exp(a*z) * cos(a*x + d*y) ) * exp(-visc*d*d*t); //u,z

      grad_u(1,0) = -a * ( d * exp(a*y) * cos(a*z + d*x) + a * exp(a*x) * cos(a*y + d*z) ) * exp(-visc*d*d*t); //v,x
      grad_u(1,1) = -a * ( a * exp(a*y) * sin(a*z + d*x) - a * exp(a*x) * sin(a*y + d*z) ) * exp(-visc*d*d*t); //v,y
      grad_u(1,2) = -a * ( a * exp(a*y) * cos(a*z + d*x) - d * exp(a*x) * sin(a*y + d*z) ) * exp(-visc*d*d*t); //v,z

      grad_u(2,0) = -a * ( a * exp(a*z) * cos(a*x + d*y) - d * exp(a*y) * sin(a*z + d*x) ) * exp(-visc*d*d*t); //w,x
      grad_u(2,1) = -a * ( d * exp(a*z) * cos(a*x + d*y) + a * exp(a*y) * cos(a*z + d*x) ) * exp(-visc*d*d*t); //w,y
      grad_u(2,2) = -a * ( a * exp(a*z) * sin(a*x + d*y) - a * exp(a*y) * sin(a*z + d*x) ) * exp(-visc*d*d*t); //w,z
    }
    else dserror("action 'calc_fluid_beltrami_error' is a 3D specific action");
  }
  break;

  case INPAR::FLUID::kimmoin_stat_stokes:
  case INPAR::FLUID::kimmoin_stat_navier_stokes:
  case INPAR::FLUID::kimmoin_instat_stokes:
  case INPAR::FLUID::kimmoin_instat_navier_stokes:
  {

    // function evaluation requires a 3D position vector!!
    double position[3];

    if(my::nsd_ == 3)
    {
      position[0] = xyzint(0);
      position[1] = xyzint(1);
      position[2] = xyzint(2);
    }
    else dserror("invalid nsd %d", my::nsd_);

    // evaluate velocity and pressure
    Teuchos::RCP<DRT::UTILS::Function> function = Teuchos::null;

    // evaluate the velocity gradient
    Teuchos::RCP<DRT::UTILS::Function> function_grad = Teuchos::null;

    bool is_stationary = false;

    // evaluate velocity and pressure
    // evaluate the velocity gradient
    if(calcerr == INPAR::FLUID::kimmoin_stat_stokes or
       calcerr == INPAR::FLUID::kimmoin_stat_navier_stokes)
    {
      is_stationary = true;
    }
    else if(calcerr == INPAR::FLUID::kimmoin_instat_stokes or
            calcerr == INPAR::FLUID::kimmoin_instat_navier_stokes)
    {
      is_stationary = false;
    }

    function      = Teuchos::rcp(new DRT::UTILS::KimMoinUP( mat, is_stationary));
    function_grad = Teuchos::rcp(new DRT::UTILS::KimMoinGradU( mat, is_stationary));

    if(my::nsd_==3)
    {
      u(0) = function->Evaluate(0,position,t,NULL);
      u(1) = function->Evaluate(1,position,t,NULL);
      u(2) = function->Evaluate(2,position,t,NULL);
      p = function->Evaluate(3,position,t,NULL);
    }
    else dserror("case 'kimmoin_stat' is a 3D specific case");


    if(my::nsd_==3)
    {
      grad_u(0,0) = function_grad->Evaluate(0,position,t,NULL); // u,x
      grad_u(0,1) = function_grad->Evaluate(1,position,t,NULL); // u,y
      grad_u(0,2) = function_grad->Evaluate(2,position,t,NULL); // u,z

      grad_u(1,0) = function_grad->Evaluate(3,position,t,NULL); // v,x
      grad_u(1,1) = function_grad->Evaluate(4,position,t,NULL); // v,y
      grad_u(1,2) = function_grad->Evaluate(5,position,t,NULL); // v,z

      grad_u(2,0) = function_grad->Evaluate(6,position,t,NULL); // w,x
      grad_u(2,1) = function_grad->Evaluate(7,position,t,NULL); // w,y
      grad_u(2,2) = function_grad->Evaluate(8,position,t,NULL); // w,z
    }
    else dserror("case 'kimmoin_stat' is a 3D specific case");
  }
  break;

  case INPAR::FLUID::shear_flow:
  {
    const double maxvel = 1.0;
    const double hight = 1.0;

    // y=0 is located in the middle of the domain
    if (my::nsd_ == 2)
    {
      p = 1.0;
      u(0) = xyzint(1)*maxvel + hight/2*maxvel;
      u(1) = 0.0;
    }
    if (my::nsd_ == 3)
    {
      p = 0.0;
      u(0) = xyzint(1)*maxvel + hight/2*maxvel;
      u(1) = 0.0;
      u(2) = 0.0;
    }

  }
  break;

  case INPAR::FLUID::gravitation:
  {
    const double gravity = 10.0;
    const double hight = 1.0;

    // 2D: rectangle 1.0x1.0
    // 3D: cube 1.0x1.0x1.0
    // y=0 is located in the middle of the domain
    if (my::nsd_ == 2)
    {
      p = -xyzint(1)*gravity + hight/2*gravity;
      u(0) = 0.0;
      u(1) = 0.0;
    }
    if (my::nsd_ == 3)
    {
      p = -xyzint(1)*gravity + hight/2*gravity;
      u(0) = 0.0;
      u(1) = 0.0;
      u(2) = 0.0;
    }
  }
  break;

  case INPAR::FLUID::channel2D:
  {
    const double maxvel=1.25;
    const double hight = 1.0;
    const double visc = 1.0;
    const double pressure_gradient = 10.0;

    // u_max = 1.25
    // y=0 is located in the middle of the channel
    if (my::nsd_ == 2)
    {
      p = 1.0;
      //p = -10*xyzint(0)+20;
      u(0) = maxvel -((hight*hight)/(2.0*visc)*pressure_gradient*(xyzint(1)/hight)*(xyzint(1)/hight));
      u(1) = 0.0;
    }
    else
      dserror("3D analytical solution is not implemented yet");
  }
  break;

  case INPAR::FLUID::jeffery_hamel_flow:
  {
    //LINALG::Matrix<3,1> physpos(true);
    //GEO::elementToCurrentCoordinates(distype, xyzint, xsi_, physpos);

    // function evaluation requires a 3D position vector!!
    double position[3];
    position[0] = xyzint(0);
    position[1] = xyzint(1);
    position[2] = 0.0;

    if (1.0 < position[0] and position[0] < 2.0 and 0.0 < position[1] and position[1] < position[0])
    {
      const double u_exact_x = DRT::Problem::Instance()->Funct(0).Evaluate(0,position,t,NULL);
      const double u_exact_y = DRT::Problem::Instance()->Funct(0).Evaluate(1,position,t,NULL);
      u(0) = u_exact_x;
      u(1) = u_exact_y;
    }
  }
  break;

  case INPAR::FLUID::byfunct:
  {

    // function evaluation requires a 3D position vector!!
    double position[3];

    if (my::nsd_ == 2)
    {

      position[0] = xyzint(0);
      position[1] = xyzint(1);
      position[2] = 0.0;
    }
    else if(my::nsd_ == 3)
    {
      position[0] = xyzint(0);
      position[1] = xyzint(1);
      position[2] = xyzint(2);
    }
    else dserror("invalid nsd %d", my::nsd_);

    if(my::nsd_ == 2)
    {
      const double u_exact_x = DRT::Problem::Instance()->Funct(calcerrfunctno-1).Evaluate(0,position,t,NULL);
      const double u_exact_y = DRT::Problem::Instance()->Funct(calcerrfunctno-1).Evaluate(1,position,t,NULL);
      const double p_exact   = DRT::Problem::Instance()->Funct(calcerrfunctno-1).Evaluate(2,position,t,NULL);

      u(0) = u_exact_x;
      u(1) = u_exact_y;
      p    = p_exact;


      std::vector<double> uder_exact_x = DRT::Problem::Instance()->Funct(calcerrfunctno-1).FctDer(0,position,t,NULL);
      std::vector<double> uder_exact_y = DRT::Problem::Instance()->Funct(calcerrfunctno-1).FctDer(1,position,t,NULL);
      //std::vector<double> pder_exact   = DRT::Problem::Instance()->Funct(func_no-1).FctDer(2,position,t,1,NULL);

      if(uder_exact_x.size())
      {
        grad_u(0,0)=uder_exact_x[0];
        grad_u(0,1)=uder_exact_x[1];
      }

      if(uder_exact_y.size())
      {
        grad_u(1,0)=uder_exact_y[0];
        grad_u(1,1)=uder_exact_y[1];
      }

    }
    else if(my::nsd_==3)
    {
      const double u_exact_x = DRT::Problem::Instance()->Funct(calcerrfunctno-1).Evaluate(0,position,t,NULL);
      const double u_exact_y = DRT::Problem::Instance()->Funct(calcerrfunctno-1).Evaluate(1,position,t,NULL);
      const double u_exact_z = DRT::Problem::Instance()->Funct(calcerrfunctno-1).Evaluate(2,position,t,NULL);
      const double p_exact   = DRT::Problem::Instance()->Funct(calcerrfunctno-1).Evaluate(3,position,t,NULL);

      u(0) = u_exact_x;
      u(1) = u_exact_y;
      u(2) = u_exact_z;
      p    = p_exact;

      std::vector<double> uder_exact_x = DRT::Problem::Instance()->Funct(calcerrfunctno-1).FctDer(0,position,t,NULL);
      std::vector<double> uder_exact_y = DRT::Problem::Instance()->Funct(calcerrfunctno-1).FctDer(1,position,t,NULL);
      std::vector<double> uder_exact_z = DRT::Problem::Instance()->Funct(calcerrfunctno-1).FctDer(2,position,t,NULL);

      if(uder_exact_x.size())
      {
        grad_u(0,0)=uder_exact_x[0];
        grad_u(0,1)=uder_exact_x[1];
        grad_u(0,2)=uder_exact_x[2];
      }

      if(uder_exact_y.size())
      {
        grad_u(1,0)=uder_exact_y[0];
        grad_u(1,1)=uder_exact_y[1];
        grad_u(1,2)=uder_exact_y[2];
      }

      if(uder_exact_z.size())
      {
        grad_u(2,0)=uder_exact_z[0];
        grad_u(2,1)=uder_exact_z[1];
        grad_u(2,2)=uder_exact_z[2];
      }

//      u(0) = 5.0+30.0*position[1];
//      u(1) = 0.0;
//      u(2) = 0.0;
//
//      p = 4.0;
//
//      grad_u(0,0) = 0.0;       grad_u(0,1) = 30.0;       grad_u(0,2) = 0.0;
//      grad_u(1,0) = 0.0;       grad_u(1,1) =  0.0;       grad_u(1,2) = 0.0;
//      grad_u(2,0) = 0.0;       grad_u(2,1) =  0.0;       grad_u(2,2) = 0.0;

    }
    else dserror("invalid dimension");
  }
  break;

  default:
    dserror("analytical solution is not defined");
    break;
  }
}


/*--------------------------------------------------------------------------------
 * compute interface error norms
 *--------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int FluidEleCalcXFEM<distype>::ComputeErrorInterface(
    DRT::ELEMENTS::Fluid *                                              ele,               ///< fluid element
    DRT::Discretization &                                               dis,               ///< background discretization
    const std::vector<int> &                                            lm,                ///< element local map
    const Teuchos::RCP<XFEM::ConditionManager> &                        cond_manager,      ///< XFEM condition manager
    Teuchos::RCP<MAT::Material>&                                        mat,               ///< material
    Epetra_SerialDenseVector&                                           ele_interf_norms,  /// squared element interface norms
    const std::map<int, std::vector<GEO::CUT::BoundaryCell*> > &        bcells,            ///< boundary cells
    const std::map<int, std::vector<DRT::UTILS::GaussIntegration> > &   bintpoints,        ///< boundary integration points
    const GEO::CUT::plain_volumecell_set &                              vcSet,             ///< set of plain volume cells
    Teuchos::ParameterList&                                             params             ///< parameter list
)
{
#ifdef DEBUG
  if(cond_manager == Teuchos::null) dserror("set the condition manager!");
#endif

  const int calcerr = DRT::INPUT::get<INPAR::FLUID::CalcError>(params,"calculate error");
  const int calcerrfunctno = DRT::INPUT::get<int>(params,"error function number");

  const double t = my::fldparatimint_->Time();


  //----------------------------------------------------------------------------
  //                         ELEMENT GEOMETRY
  //----------------------------------------------------------------------------

  // ---------------------------------------------------------------------
  // get initial node coordinates for element
  // ---------------------------------------------------------------------
  // get node coordinates
  GEO::fillInitialPositionArray< distype, my::nsd_, LINALG::Matrix<my::nsd_,my::nen_> >( ele, my::xyze_ );

  // ---------------------------------------------------------------------
  // get additional state vectors for ALE case: grid displacement and vel.
  // ---------------------------------------------------------------------

  LINALG::Matrix<my::nsd_, my::nen_> edispnp(true);
  LINALG::Matrix<my::nsd_, my::nen_> egridv(true);

  if (ele->IsAle()) my::GetGridDispVelALE(dis, lm, edispnp, egridv);
  // add displacement when fluid nodes move in the ALE case
  if (ele->IsAle()) my::xyze_ += edispnp;


  // ---------------------------------------------------------------------

  /// element coordinates in EpetraMatrix
  Epetra_SerialDenseMatrix ele_xyze(my::nsd_,my::nen_);
  for ( int i=0; i<my::nen_; ++i )
  {
    for(int j=0; j<my::nsd_; j++)
      ele_xyze(j,i) = my::xyze_( j, i );
  }

  // ---------------------------------------------------------------------
  // get velocity state vectors
  // ---------------------------------------------------------------------

  // get element-wise velocity/pressure field
  LINALG::Matrix<my::nsd_,my::nen_> evelaf(true);
  LINALG::Matrix<my::nen_,1> epreaf(true);
  my::ExtractValuesFromGlobalVector(dis, lm, *my::rotsymmpbc_, &evelaf, &epreaf, "u and p at time n+1 (converged)");


  // ---------------------------------------------------------------------
  // set element advective field for Oseen problems
  // ---------------------------------------------------------------------
  if (my::fldpara_->PhysicalType()==INPAR::FLUID::oseen) my::SetAdvectiveVelOseen(ele);


  // ---------------------------------------------------------------------
  // get the element volume
  // ---------------------------------------------------------------------

  // evaluate shape functions and derivatives at element center
  my::EvalShapeFuncAndDerivsAtEleCenter();
  // set element area or volume
  const double vol = my::fac_;

  //-----------------------------------------------------------------------------------
  //         evaluate element length, stabilization factors and average weights
  //-----------------------------------------------------------------------------------

  // element length
  double h_k = 0.0;
  double inv_hk = 0.0;

  // take a volume based element length
  h_k = ComputeVolEqDiameter(vol);
  inv_hk = 1.0/h_k;


  // evaluate shape function derivatives
  bool eval_deriv = true;


  //-----------------------------------------------------------------------------------
  //         initialize analytical solution vectors and error variables
  //-----------------------------------------------------------------------------------

  // analytical solution
  LINALG::Matrix<my::nsd_,1>         u_analyt(true);
  LINALG::Matrix<my::nsd_,my::nsd_>  grad_u_analyt(true);
  double p_analyt = 0.0;

  // error
  LINALG::Matrix<my::nsd_,1>        u_err(true);
  LINALG::Matrix<my::nsd_,my::nsd_> grad_u_err(true);
  double p_err = 0.0;

  LINALG::Matrix<my::nsd_,1> flux_u_err(true);
  LINALG::Matrix<my::nsd_,1> flux_p_err(true);


  //--------------------------------------------
  // loop intersecting sides
  //--------------------------------------------
  // map of side-element id and Gauss points
  for ( std::map<int, std::vector<DRT::UTILS::GaussIntegration> >::const_iterator i=bintpoints.begin();
      i!=bintpoints.end();
      ++i )
  {
    //-----------------------------------------------------------------------------------

    // interface normal vector, pointing from background domain into the interface
    LINALG::Matrix<3,1> normal(true);
    // gauss-point coordinates
    LINALG::Matrix<3,1> x_side(true);

    // we need an interface to the boundary element (for projection)
    Teuchos::RCP<DRT::ELEMENTS::XFLUID::SlaveElementInterface<distype> > si;

    // location array of boundary element
    DRT::Element::LocationArray cutla( 1 );

    // pointer to boundary element
    DRT::Element * side = NULL;

    // location array of element to couple with (only used for embedded fluid problems)
    DRT::Element::LocationArray coupl_la( 1 );

    // coordinates of boundary element
    Epetra_SerialDenseMatrix side_xyze;

    //-----------------------------------------------------------------------------------
    // only used for couplings:

    // coupling object between background element and each coupling element (side for xfluid-sided coupling, element for other couplings)
    Teuchos::RCP<DRT::ELEMENTS::XFLUID::SlaveElementInterface<distype> > ci;

    // pointer to coupling element
    DRT::Element * coupl_ele = NULL;

    // coupling element coordinates
    Epetra_SerialDenseMatrix coupl_xyze;

    //-----------------------------------------------------------------------------------

    //-----------------------------------------------------------------------------------

    int coup_sid = i->first; // global coupling side id

    // get the coupling strategy for coupling of two fields
    const XFEM::EleCoupCond & coupcond = cond_manager->GetCouplingCondition(coup_sid, my::eid_);
    const INPAR::XFEM::EleCouplingCondType & cond_type = coupcond.first;

    const std::vector<DRT::UTILS::GaussIntegration> & cutintpoints = i->second;

    std::map<int, std::vector<GEO::CUT::BoundaryCell*> >::const_iterator j = bcells.find( coup_sid );
    if ( j==bcells.end() )
      dserror( "missing boundary cell" );

    const std::vector<GEO::CUT::BoundaryCell*> & bcs = j->second;
    if ( bcs.size()!=cutintpoints.size() )
      dserror( "boundary cell integration rules mismatch" );


    //---------------------------------------------------------------------------------
    // set flags used for coupling with given levelset/mesh coupling side
    bool is_ls_coupling_side   = cond_manager->IsLevelSetCoupling(coup_sid);
    bool is_mesh_coupling_side = cond_manager->IsMeshCoupling(coup_sid);

    Teuchos::RCP<DRT::Discretization> cutter_dis = cond_manager->GetCutterDis(coup_sid);

#ifdef DEBUG
    if( is_ls_coupling_side and  is_mesh_coupling_side) dserror("side cannot be a levelset-coupling side and a mesh coupling side at once: side %i", coup_sid);
    if(!is_ls_coupling_side and !is_mesh_coupling_side) dserror("side is neither a levelset-coupling side nor a mesh coupling side: side %i", coup_sid);
#endif
    //-----------------------------------------------------------------------------------

    //---------------------------------------------------------------------------------
    // prepare the coupling objects
    if(is_mesh_coupling_side)
    {
      // get the side element and its coordinates for projection of Gaussian points
      side = cond_manager->GetSide( coup_sid );
      GEO::InitialPositionArray(side_xyze,side);

      // create auxiliary coupling object for the boundary element, in order to perform projection
      si = DRT::ELEMENTS::XFLUID::SlaveElementInterface<distype>::CreateSlaveElementRepresentation( side, side_xyze );

      // set displacement of side
      side->LocationVector(*cutter_dis,cutla,false);
      si->AddSlaveEleDisp(*cutter_dis,cutla[0].lm_);

      if( cond_type == INPAR::XFEM::CouplingCond_SURF_WEAK_DIRICHLET or
          cond_type == INPAR::XFEM::CouplingCond_SURF_FSI_PART or
          cond_type == INPAR::XFEM::CouplingCond_SURF_CRACK_FSI_PART or
          cond_type == INPAR::XFEM::CouplingCond_SURF_NAVIER_SLIP)
      {
        si->SetInterfaceJumpStatenp(*cutter_dis, "ivelnp", cutla[0].lm_);
      }

      if (cond_type == INPAR::XFEM::CouplingCond_SURF_FLUIDFLUID)
      {
        // force to get the embedded element, even if background-sided coupling is active
        coupl_ele = cond_manager->GetCondElement(coup_sid);

        GEO::InitialPositionArray(coupl_xyze,coupl_ele);

        ci = DRT::ELEMENTS::XFLUID::SlaveElementInterface<distype>::CreateSlaveElementRepresentation(coupl_ele,coupl_xyze);

        // set velocity (and pressure) of coupling/slave element at current time step
        const int coup_idx = cond_manager->GetCouplingIndex(coup_sid,my::eid_);
        coupl_ele->LocationVector(*cond_manager->GetCouplingByIdx(coup_idx)->GetCondDis(),coupl_la,false);
        ci->SetSlaveState(*cond_manager->GetCouplingByIdx(coup_idx)->GetCondDis(),coupl_la[0].lm_);
      }
    }

    if (cond_manager->HasAveragingStrategy(INPAR::XFEM::Xfluid_Sided))
    {
      h_k = ComputeCharEleLength(ele, ele_xyze, cond_manager, vcSet, bcells, bintpoints);
      inv_hk = 1.0 / h_k;
    }

    //--------------------------------------------
    // loop boundary cells w.r.t current cut side
    //--------------------------------------------
    for ( std::vector<DRT::UTILS::GaussIntegration>::const_iterator i=cutintpoints.begin();
        i!=cutintpoints.end();
        ++i )
    {
      const DRT::UTILS::GaussIntegration & gi = *i;
      GEO::CUT::BoundaryCell * bc = bcs[i - cutintpoints.begin()]; // get the corresponding boundary cell

      //--------------------------------------------
      // loop gausspoints w.r.t current boundary cell
      //--------------------------------------------
      for ( DRT::UTILS::GaussIntegration::iterator iquad=gi.begin(); iquad!=gi.end(); ++iquad )
      {
        double drs = 0.0; // transformation factor between reference cell and linearized boundary cell

        const LINALG::Matrix<2,1> eta( iquad.Point() ); // xi-coordinates with respect to side

        LINALG::Matrix<3,1> rst(true); // local coordinates w.r.t background element

        LINALG::Matrix<3,1> x_gp_lin(true); // gp in xyz-system on linearized interface

        // compute transformation factor, normal vector and global Gauss point coordiantes
        if(bc->Shape() != DRT::Element::dis_none) // Tessellation approach
        {
          ComputeSurfaceTransformation(drs, x_gp_lin, normal, bc, eta);
        }
        else // MomentFitting approach
        {
          drs = 1.0;
          normal = bc->GetNormalVector();
          const double* gpcord = iquad.Point();
          for (int idim=0;idim<3;++idim)
          {
            x_gp_lin(idim,0) = gpcord[idim];
          }
        }

        // find element local position of gauss point
        GEO::CUT::Position<distype> pos( my::xyze_, x_gp_lin );
        pos.Compute();
        rst = pos.LocalCoordinates();

        //        if (!levelset_cut)
        //        {
        //          // project gaussian point from linearized interface to warped side (get/set local side coordinates in SideImpl)
        //          LINALG::Matrix<2,1> xi_side(true);
        //
        //          //          side_impl[sid]->ProjectOnSide(x_gp_lin, x_side, xi_side);
        //          si->ProjectOnSide(x_gp_lin, x_side, xi_side);
        //        }
        //        else x_side = x_gp_lin;


        //TODO unify ProjectOnSide and Evaluate for different spatial dimensions of boundary slave element and volumetric slave element
        if (is_mesh_coupling_side)
        {
          // project gaussian point from linearized interface to warped side (get/set local side coordinates in SideImpl)
          LINALG::Matrix<3,1> xi_side(true);
          // project on boundary element
          si->ProjectOnSide(x_gp_lin, x_side, xi_side);

          if (cond_type == INPAR::XFEM::CouplingCond_SURF_FLUIDFLUID)
            ci->Evaluate( x_side ); // evaluate embedded element's shape functions at gauss-point coordinates
        }
        else if (is_ls_coupling_side)
        {
          // TODO: do we need this here?
          //          if(cond_manager->IsCoupling( coup_sid, my::eid_ ))
          //            ci->Evaluate( x_gp_lin ); // evaluate embedded element's shape functions at gauss-point coordinates
        }


        const double surf_fac = drs*iquad.Weight();

        //--------------------------------------------
        // evaluate shape functions (and derivatives)

        if(eval_deriv)
        {
          EvalFuncAndDeriv( rst );
        }
        else
        {
          DRT::UTILS::shape_function<distype>( rst, my::funct_ );
        }


        // get velocity at integration point
        // (values at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
        my::velint_.Multiply(evelaf,my::funct_);

        // get velocity derivatives at integration point
        // (values at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
        my::vderxy_.MultiplyNT(evelaf,my::derxy_);

        // get pressure at integration point
        // (value at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
        double press = my::funct_.Dot(epreaf);

        //----------------------------------------------
        // get convective velocity at integration point
        my::SetConvectiveVelint(ele->IsAle());

#if(0)
        // FOR INTERPOLATION ERRORS
        // I_h(u_analyt) - u_analyt
        // grad(I_h(u_analyt) - u_analyt)
        // I_h(p_analyt) - p_analyt
        //---------------------------------------------------------------
        LINALG::Matrix<my::nsd_,1> u_analyt_interpolated(true);
        LINALG::Matrix<my::nsd_,my::nsd_> grad_u_analyt_interpolated(true);
        double p_analyt_interpolated = 0.0;

        LINALG::Matrix<my::nsd_,my::nen_> u_analyt_element(true);
        LINALG::Matrix<my::nen_,1> p_analyt_element(true);


        //---------------------------------------------------------------
        for(int i=0; i<my::nen_; ++i)
        {
          LINALG::Matrix<my::nsd_,1> u_analyt_node(true);
          LINALG::Matrix<my::nsd_,my::nsd_> grad_u_analyt_node(true);
          double p_analyt_node = 0.0;


          LINALG::Matrix<my::nsd_,1> x_node(true);
          x_node(0) = my::xyze_(0,i);
          x_node(1) = my::xyze_(1,i);
          x_node(2) = my::xyze_(2,i);

          AnalyticalReference(
               calcerr,          ///< which reference solution
               u_analyt_node,         ///< exact velocity (onesided), exact jump vector (coupled)
               grad_u_analyt_node,    ///< exact velocity gradient
               p_analyt_node,         ///< exact pressure
               x_node,           ///< xyz position of node
               t,                ///< time t
               mat
           );

          u_analyt_element(0,i) = u_analyt_node(0);
          u_analyt_element(1,i) = u_analyt_node(1);
          u_analyt_element(2,i) = u_analyt_node(2);

          p_analyt_element(i)=p_analyt_node;
        }

        u_analyt_interpolated.Multiply(u_analyt_element, my::funct_);
        grad_u_analyt_interpolated.MultiplyNT(u_analyt_element,my::derxy_);

        p_analyt_interpolated = my::funct_.Dot(p_analyt_element);
#endif

        //--------------------------------------------
        // compute errors

        LINALG::Matrix<my::nsd_,1>        u_analyt(true);      // boundary condition to enforce (xfsi), interfacial jump to enforce (fluidfluid)
        LINALG::Matrix<my::nsd_,my::nsd_> grad_u_analyt(true);
        p_analyt = 0.0;

        AnalyticalReference(
            calcerr,          ///< which reference solution
            calcerrfunctno,   ///< error function number
            u_analyt,         ///< exact velocity (onesided), exact jump vector (coupled)
            grad_u_analyt,    ///< exact velocity gradient
            p_analyt,         ///< exact pressure
            x_gp_lin,         ///< xyz position of gaussian point
            t,                ///< time t
            mat
        );

        LINALG::Matrix<my::nsd_,1> velint_s;
        if (is_mesh_coupling_side)
        {
          si->GetInterfaceVelnp(velint_s);
        }

        if (cond_type == INPAR::XFEM::CouplingCond_SURF_FLUIDFLUID)
        {
          u_err.Update(1.0, my::velint_, -1.0, velint_s, 0.0);

          LINALG::Matrix<my::nsd_,my::nsd_> grad_u_side(true);
          ci->GetInterfaceVelGradnp(grad_u_side);

          grad_u_err.Update(1.0, my::vderxy_, -1.0, grad_u_side, 0.0);

          double press_coupl = 0.0;
          ci->GetInterfacePresnp(press_coupl);
          //p_err = p_background - p_emb;
          p_err = press - press_coupl;
        }
        else
        {
          u_err.Update(1.0, my::velint_, -1.0, u_analyt, 0.0);
          grad_u_err.Update(1.0, my::vderxy_, -1.0, grad_u_analyt, 0.0);
          p_err = press - p_analyt;

#if(0)
          // FOR interpolation error
          u_err.Update(1.0, u_analyt_interpolated, -1.0, u_analyt, 0.0);
          grad_u_err.Update(1.0, grad_u_analyt_interpolated, -1.0, grad_u_analyt, 0.0);
          p_err = p_analyt_interpolated - p_analyt;
#endif
        }

        flux_u_err.Multiply(grad_u_err,normal);
        flux_p_err.Update(p_err,normal,0.0);

        /*
         * Scaling of interface error norms:
         *
         *                       (1)           (2)          (3)
         *                    /  \mu    \rho             h * \rho         \
         *  NIT :=  \gamma * |    --  +  -- * |u|_inf  + ----------------- |
         *                    \   h      6               12 * \theta * dt /
         *
         *                             interface errors
         *  1.   || nu/h * (u_b - u_e - u_jump) ||_L_2(Gamma)        =  broken H1/2 Sobolev norm for boundary/coupling condition
         *  2.   || nu^(+1/2) grad( u_b - u_e )*n ||_H-1/2(Gamma)    =  standard H-1/2 Sobolev norm for normal flux (velocity part)
         *  3.   || (p_b - p_e)*n ||_H-1/2(Gamma)                    =  standard H-1/2 Sobolev norm for normal flux (pressure part)
         *  4.   || (u*n)_inflow (u - u*) ||_L2(Gamma)               =  L^2 Sobolev norm for inflow boundary/coupling condition
         *  5.   || (sigma*h+|u|+nu/h)^(+1/2) (u - u*)*n ||_L2(Gamma) = L^2 Sobolev norm for mass conservation coupling condition
         */
        double u_err_squared        = 0.0;
        double u_err_squared_normal = 0.0;
        double flux_u_err_squared   = 0.0;
        double flux_p_err_squared   = 0.0;

        // evaluate squared errors at gaussian point
        for (int isd=0;isd<my::nsd_;isd++)
        {
          u_err_squared        += u_err(isd)*u_err(isd)*surf_fac;
          u_err_squared_normal += u_err(isd)*normal(isd)*normal(isd)*u_err(isd)*surf_fac;
          flux_u_err_squared   += flux_u_err(isd)*flux_u_err(isd)*surf_fac;
          flux_p_err_squared   += flux_p_err(isd)*flux_p_err(isd)*surf_fac;
        }

        // interface errors
        double nit_stabfac = 0.0;

        //Needs coupling condition to get kappas!
        const double kappa_m = 1.0;
        const double kappa_s = 0.0;
        double visc_stab_fac = 0.0;
        cond_manager->Get_ViscPenalty_Stabfac(coup_sid, ele,kappa_m,kappa_s, inv_hk,fldparaxfem_,visc_stab_fac);
        NIT_Compute_FullPenalty_Stabfac(nit_stabfac,normal,h_k,kappa_m,kappa_s,my::convvelint_,velint_s,visc_stab_fac,true);

        const double veln_normal = my::convvelint_.Dot(normal); // TODO: shift this to routine
        double NIT_inflow_stab = std::max(0.0,-veln_normal);

        ele_interf_norms[0] += visc_stab_fac * u_err_squared;
        ele_interf_norms[1] += h_k * my::visc_ * flux_u_err_squared;
        ele_interf_norms[2] += h_k * flux_p_err_squared;
        ele_interf_norms[3] += NIT_inflow_stab * u_err_squared;
        ele_interf_norms[4] += nit_stabfac * u_err_squared_normal;

      } // end loop gauss points of boundary cell
    } // end loop boundary cells of side

  } // end loop cut sides

  return 0;
}


/*--------------------------------------------------------------------------------
 * add mixed/hybrid stress-based LM interface condition to element matrix and rhs
 *--------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void FluidEleCalcXFEM<distype>::ElementXfemInterfaceHybridLM(
    DRT::ELEMENTS::Fluid *                                            ele,                      ///< fluid element
    DRT::Discretization &                                             dis,                      ///< background discretization
    const std::vector<int> &                                          lm,                       ///< element local map
    const Teuchos::RCP<XFEM::ConditionManager> &                      cond_manager,             ///< XFEM condition manager
    const std::vector<DRT::UTILS::GaussIntegration> &                 intpoints,                ///< element gauss points
    const std::map<int, std::vector<GEO::CUT::BoundaryCell*> >&       bcells,                   ///< boundary cells
    const std::map<int, std::vector<DRT::UTILS::GaussIntegration> > & bintpoints,               ///< boundary integration points
    const std::map<int, std::vector<int> > &                          patchcouplm,              ///< lm vectors for coupling elements, key= global coupling side-Id
    std::map<int, std::vector<Epetra_SerialDenseMatrix> > &           side_coupling,            ///< side coupling matrices
    Teuchos::ParameterList&                                           params,                   ///< parameter list
    Teuchos::RCP<MAT::Material>&                                      mat,                      ///< material
    Epetra_SerialDenseMatrix&                                         elemat1_epetra,           ///< local system matrix of intersected element
    Epetra_SerialDenseVector&                                         elevec1_epetra,           ///< local element vector of intersected element
    Epetra_SerialDenseMatrix&                                         Cuiui,                    ///< coupling matrix of a side with itself
    const GEO::CUT::plain_volumecell_set &                            vcSet                    ///< set of plain volume cells
)
{
#ifdef DEBUG
  if(cond_manager == Teuchos::null) dserror("set the condition manager!");
#endif

  //--------------------------------------------------------
  // determine, whether this is a Cauchy stress-based (MHCS)
  // or viscous stress-based LM approach (MHVS)
  //--------------------------------------------------------
  INPAR::XFEM::CouplingMethod coupling_method = fldparaxfem_->GetCouplingMethod();

  // check for valid boundary integration type
  switch (coupling_method)
  {
  case INPAR::XFEM::Hybrid_LM_viscous_stress:
  case INPAR::XFEM::Hybrid_LM_Cauchy_stress:
    break;
  case INPAR::XFEM::Nitsche:
    dserror("Wrong evaluation routine for Nitsche coupling. Try ElementXfemInterfaceNIT/NIT2 instead.");
    break;
  default:
    dserror("Landed in evaluation routine for stress-based LM, given an unknown or unsupported coupling method.");
    break;
  }

  const bool is_MHVS = (coupling_method == INPAR::XFEM::Hybrid_LM_viscous_stress);


  // REMARK: to avoid confusion -
  // 'side' = boundary element, part of cutdis (can be warped)
  // 'boundary cell' = belongs to volume-cell

  // do we need convective stabilization?
  bool add_conv_stab( fldparaxfem_->XffConvStabScaling() != INPAR::XFEM::XFF_ConvStabScaling_none ||
      fldparaxfem_->ConvStabScaling() != INPAR::XFEM::ConvStabScaling_none);

  //----------------------------------------------------------------------------
  //                         ELEMENT GEOMETRY
  //----------------------------------------------------------------------------

  // ---------------------------------------------------------------------
  // get initial node coordinates for element
  // ---------------------------------------------------------------------
  // get node coordinates
  GEO::fillInitialPositionArray< distype, my::nsd_, LINALG::Matrix<my::nsd_,my::nen_> >( ele, my::xyze_ );

  // ---------------------------------------------------------------------
  // get additional state vectors for ALE case: grid displacement and vel.
  // ---------------------------------------------------------------------

  LINALG::Matrix<my::nsd_, my::nen_> edispnp(true);
  LINALG::Matrix<my::nsd_, my::nen_> egridv(true);

  if (ele->IsAle()) my::GetGridDispVelALE(dis, lm, edispnp, egridv);


  // ---------------------------------------------------------------------

  /// element coordinates in EpetraMatrix
  Epetra_SerialDenseMatrix ele_xyze(my::nsd_,my::nen_);
  for ( int i=0; i<my::nen_; ++i )
  {
    for(int j=0; j<my::nsd_; j++)
      ele_xyze(j,i) = my::xyze_( j, i );
  }

  // ---------------------------------------------------------------------
  // get velocity state vectors
  // ---------------------------------------------------------------------

  // get element-wise velocity/pressure field for current time step
  LINALG::Matrix<my::nsd_,my::nen_> evelaf(true);
  LINALG::Matrix<my::nen_,1> epreaf(true);
  my::ExtractValuesFromGlobalVector(dis, lm, *my::rotsymmpbc_, &evelaf, &epreaf, "velaf");

  // get element-wise velocity/pressure field for previous time step
  LINALG::Matrix<my::nsd_,my::nen_> eveln(true);
  LINALG::Matrix<my::nen_,1> epren(true);
  if (my::fldparatimint_->IsNewOSTImplementation())
    my::ExtractValuesFromGlobalVector(dis, lm, *my::rotsymmpbc_, &eveln, &epren, "veln");

  // ---------------------------------------------------------------------
  // set element advective field for Oseen problems
  // ---------------------------------------------------------------------
  if (my::fldpara_->PhysicalType()==INPAR::FLUID::oseen) my::SetAdvectiveVelOseen(ele);

  // compute characteristic element length based on the background element
  const double h_k = ComputeCharEleLength(ele,ele_xyze,cond_manager,vcSet,bcells,bintpoints);

  //--------------------------------------------------------
  // declaration of matrices & rhs
  //--------------------------------------------------------

  // sub-blocks of matrix K_{\sigma\sigma} (-->K_ss)
  LINALG::Matrix<my::nen_,my::nen_>  bK_ss(true);        // N * N^T
  LINALG::Matrix<my::nen_,my::nen_>  invbK_ss(true);     // inverse of bK_ss, (N * N^T)^-1
  LINALG::Matrix<my::nen_,my::nen_>  halfInvbK_ss(true); // inverse scaled by 1/2

  // The block matrices K_... result from the volume integrals on the cut element.
  // In case of a viscous stress-based approach (MHVS), there is no term like K_sp,
  // that couples the Lagrange multiplier stresses with the pressure.
  // Instead, we have contributions G_up and G_pu from surface coupling terms, that don't play a role
  // in the condensation of the stress-based LM (in contrast to velocity coupling terms s-u and u-s).
  // For compatibility reasons, we include column- & row-blocks for the pressure in K_su and K_us.
  // In the case of a MHVS-approach, these remain empty, as we have G_up and G_pu.

  LINALG::BlockMatrix<LINALG::Matrix<my::nen_,my::nen_>,numstressdof_,my::numdofpernode_> K_su; // s-u, s-p
  LINALG::BlockMatrix<LINALG::Matrix<my::nen_,my::nen_>,my::numdofpernode_,numstressdof_> K_us; // u-s, p-s

  LINALG::BlockMatrix<LINALG::Matrix<my::nen_,my::nen_>,numstressdof_,numstressdof_>      invK_ss; // K_ss^-1
  LINALG::BlockMatrix<LINALG::Matrix<my::nen_,1>,numstressdof_,1>                         rhs_s;

  // Only for MHVS:

  LINALG::BlockMatrix<LINALG::Matrix<my::nen_,my::nen_>,my::nsd_,my::nsd_>      K_uu;
  LINALG::BlockMatrix<LINALG::Matrix<my::nen_,1>,my::nsd_,1>                    rhs_uu;

  // surface-based pressure terms, analogous to Nitsche
  LINALG::BlockMatrix<LINALG::Matrix<my::nen_,my::nen_>,my::nsd_,1>             G_up;
  LINALG::BlockMatrix<LINALG::Matrix<my::nen_,my::nen_>,1,my::nsd_>             G_pu;

  // rhs-contributions from interface integration
  LINALG::BlockMatrix<LINALG::Matrix<my::nen_,1>,my::nsd_,1>                    rhs_up;
  LINALG::Matrix<my::nen_,1>                                                    rhs_pu(true);

  //--------------------------------------------------------
  // build matrices K (based on volume terms)
  //--------------------------------------------------------

  // in case of MHVS we have a stabilizing parameter!
  double mhvs_param = 1.0;
  if (is_MHVS)
  {
    // MHVS-stabilization parameter n
    // REMARK:
    // NIT_visc_stab_fac =   gamma * mu * C^2
    // (C^2 includes characteristic length scale);
    // the analogous MHVS-factor = 2 * n * mu * meas(\Gamma)/meas(\Omega_K)
    // gamma <--> 2 * n
    double mhvs_param = fldparaxfem_->NITStabScaling() / 2.0;
    if ( fabs(mhvs_param) < 1.e-8 )
      dserror("MHVS stabilizing parameter n appears in denominator. Please avoid choosing 0.");
  }

  // build volumetric coupling matrices
  HybridLM_Build_VolBased(
      intpoints, vcSet,
      evelaf, epreaf,
      bK_ss, invbK_ss, K_su, rhs_s, K_us, K_uu, rhs_uu,
      is_MHVS, mhvs_param);

  /*--------------------------------------------------------
    build surface coupling terms
  --------------------------------------------------------*/

  //-----------------------------------------------------------------------------------

  // side coupling implementation between background element and each cutting side OR
  // embedded element
  std::map<int,Teuchos::RCP<DRT::ELEMENTS::XFLUID::HybridLMInterface<distype> > > ci;

  //-----------------------------------------------------------------------------------
  //                     application-specific flags & parameters
  //-----------------------------------------------------------------------------------

  // map of boundary element gids and coupling matrices, [0]: Gsui, [1]: Guis
  std::map<int, std::vector<Epetra_SerialDenseMatrix> > Cuiui_coupling;

  std::vector<int> patchelementslm;


  // auxiliary coupling implementation for terms
  // find all the intersecting elements of actele
  std::set<int> begids;
  for (std::map<int,  std::vector<GEO::CUT::BoundaryCell*> >::const_iterator bc=bcells.begin();
       bc!=bcells.end(); ++bc )
  {
    const int coup_sid = bc->first;

    // get the coupling strategy for coupling of two fields
    const INPAR::XFEM::AveragingStrategy averaging_strategy = cond_manager->GetAveragingStrategy(coup_sid, my::eid_);

    begids.insert(coup_sid);

    if(!cond_manager->IsCoupling( coup_sid, my::eid_ )) continue; // no coupling with current side

    if(cond_manager->IsLevelSetCoupling(coup_sid)) dserror("PatchLocationVector for level-set coupling not supported for hybrid-lm methods yet");

    // get coupling matrices for the current side (boundary element)
    std::vector<Epetra_SerialDenseMatrix> & Cuiui_matrices = Cuiui_coupling[coup_sid]; // create new vector of Coupling matrices

    std::map<int, std::vector<int> >::const_iterator j = patchcouplm.find( coup_sid );
    if ( j==patchcouplm.end() )
      dserror( "missing side" );

    const std::vector<int> & patchlm = j->second;

    // get number of dofs for coupling side/element
    const size_t ndof_i = j->second.size();

    patchelementslm.reserve( patchelementslm.size() + ndof_i);
    patchelementslm.insert(patchelementslm.end(), patchlm.begin(), patchlm.end());

    if(averaging_strategy != INPAR::XFEM::Xfluid_Sided)
      dserror("Embedded-sided or Mean or Harmonic coupling for stress-based hybrid LM approach is not yet available!");

    Cuiui_matrices.resize(2);
    Cuiui_matrices[0].Shape(my::nen_*numstressdof_,ndof_i); //Gsui (coupling between background elements sigma and current side!)
    Cuiui_matrices[1].Shape(ndof_i,my::nen_*numstressdof_); //Guis

  }


  // map of Nitsche-interfaces for building convective stabilization terms
  std::map<int,Teuchos::RCP<DRT::ELEMENTS::XFLUID::NitscheInterface<distype> > > si_nit;

  // map of boundary element gids and coupling contributions from convective stabilization terms
  std::map<int, std::vector<Epetra_SerialDenseMatrix> > side_coupling_extra;

  // reshape coupling matrices for convective stabilization terms
  if ( add_conv_stab || my::fldparatimint_->IsNewOSTImplementation() )
  {
    HybridLM_CreateSpecialContributionMatrices(cond_manager,begids,side_coupling_extra);
  }


  //--------------------------------------------
  // loop intersecting sides
  //--------------------------------------------
  // map of side-element id and Gauss points
  for ( std::map<int, std::vector<DRT::UTILS::GaussIntegration> >::const_iterator i=bintpoints.begin();
        i!=bintpoints.end();
        ++i )
  {
    //-----------------------------------------------------------------------------------

    // interface normal vector, pointing from background domain into the interface
    LINALG::Matrix<3,1> normal(true);
    // gauss-point coordinates
    LINALG::Matrix<3,1> x_side(true);

    // we need an interface to the boundary element (for projection)
    Teuchos::RCP<DRT::ELEMENTS::XFLUID::SlaveElementInterface<distype> > si;

    // location array of boundary element
    DRT::Element::LocationArray cutla( 1 );

    // pointer to boundary element
    DRT::Element * side = NULL;

    // coordinates of boundary element
    Epetra_SerialDenseMatrix side_xyze;

    //-----------------------------------------------------------------------------------
    // only used for couplings:

    // pointer to coupling element
    DRT::Element * coupl_ele = NULL;

    // coupling element coordinates
    Epetra_SerialDenseMatrix coupl_xyze;

    //-----------------------------------------------------------------------------------


    int coup_sid = i->first;

    // get the coupling strategy for coupling of two fields
    const INPAR::XFEM::AveragingStrategy averaging_strategy = cond_manager->GetAveragingStrategy(coup_sid, my::eid_);
    const XFEM::EleCoupCond & coupcond = cond_manager->GetCouplingCondition(coup_sid, my::eid_);
    const INPAR::XFEM::EleCouplingCondType & cond_type = coupcond.first;

    const int coup_idx = cond_manager->GetCouplingIndex(coup_sid, my::eid_);
    Teuchos::RCP<XFEM::CouplingBase> coupling = cond_manager->GetCouplingByIdx(coup_idx);

    const std::vector<DRT::UTILS::GaussIntegration> & cutintpoints = i->second;

    // get side's boundary cells
    std::map<int, std::vector<GEO::CUT::BoundaryCell*> >::const_iterator j = bcells.find( coup_sid );
    if ( j==bcells.end() )
      dserror( "missing boundary cell" );

    const std::vector<GEO::CUT::BoundaryCell*> & bcs = j->second;
    if ( bcs.size()!=cutintpoints.size() )
      dserror( "boundary cell integration rules mismatch" );

    //---------------------------------------------------------------------------------
    // set flags used for coupling with given levelset/mesh coupling side
    bool is_ls_coupling_side   = cond_manager->IsLevelSetCoupling(coup_sid);
    bool is_mesh_coupling_side = cond_manager->IsMeshCoupling(coup_sid);

    Teuchos::RCP<DRT::Discretization> cutter_dis = cond_manager->GetCutterDis(coup_sid);

#ifdef DEBUG
    if( is_ls_coupling_side and  is_mesh_coupling_side) dserror("side cannot be a levelset-coupling side and a mesh coupling side at once: side %i", coup_sid);
    if(!is_ls_coupling_side and !is_mesh_coupling_side) dserror("side is neither a levelset-coupling side nor a mesh coupling side: side %i", coup_sid);
#endif

    //-----------------------------------------------------------------------------------
    Teuchos::RCP<XFEM::MeshCouplingFSI> mc_fsi = Teuchos::null;
    bool assemble_iforce = false;

    //---------------------------------------------------------------------------------
    // prepare the coupling objects
    //---------------------------------------------------------------------------------
    // prepare the coupling objects
    if(is_mesh_coupling_side)
    {
      mc_fsi = Teuchos::rcp_dynamic_cast<XFEM::MeshCouplingFSI>(coupling);
      if(mc_fsi != Teuchos::null) assemble_iforce = true;

      // get the side element and its coordinates for projection of Gaussian points
      side = cond_manager->GetSide( coup_sid );
      GEO::InitialPositionArray(side_xyze,side);

      // create auxiliary coupling object for the boundary element, in order to perform projection
      si = DRT::ELEMENTS::XFLUID::SlaveElementInterface<distype>::CreateSlaveElementRepresentation( side, side_xyze );

      // set displacement of side
      side->LocationVector(*cutter_dis,cutla,false);
      si->AddSlaveEleDisp(*cutter_dis,cutla[0].lm_);

      if(cond_type == INPAR::XFEM::CouplingCond_SURF_WEAK_DIRICHLET or
         cond_type == INPAR::XFEM::CouplingCond_SURF_FSI_PART or
         cond_type == INPAR::XFEM::CouplingCond_SURF_CRACK_FSI_PART or
         cond_type == INPAR::XFEM::CouplingCond_SURF_NAVIER_SLIP)
      {
        si->SetInterfaceJumpStatenp(*cutter_dis, "ivelnp", cutla[0].lm_);
        if (my::fldparatimint_->IsNewOSTImplementation())
          si->SetInterfaceJumpStaten(*cutter_dis, "iveln", cutla[0].lm_);
      }
    }

    Teuchos::RCP<DRT::Discretization> coupl_dis_ = cond_manager->GetCouplingDis(coup_sid);

    if(!(is_ls_coupling_side and !cond_manager->IsCoupling( coup_sid, my::eid_ ))) // not level-set-WDBC case
    {
      if (averaging_strategy == INPAR::XFEM::Embedded_Sided or
          averaging_strategy == INPAR::XFEM::Mean) // for coupling-sided and two-sided coupling
        dserror("embedded or two-sided coupling not supported");
      else
      {
        // TODO get the coupling element / the coupling side!!!
        coupl_ele = cond_manager->GetSide(coup_sid);
      }

      GEO::InitialPositionArray(coupl_xyze,coupl_ele);
    }

    if(!cond_manager->IsCoupling( coup_sid, my::eid_ ))
    {
      if(is_ls_coupling_side) //... for problems with cut interface defined by level-set field, currently only one-sided
      {
        ci[coup_sid] = DRT::ELEMENTS::XFLUID::HybridLMInterface<distype>::CreateHybridLMCoupling_XFluidWDBC(
            fldparaxfem_->IsViscousAdjointSymmetric());
      }
      else if(is_mesh_coupling_side)
      {
        ci[coup_sid] = DRT::ELEMENTS::XFLUID::HybridLMInterface<distype>::CreateHybridLMCoupling_XFluidWDBC(
            coupl_ele,coupl_xyze,fldparaxfem_->IsViscousAdjointSymmetric());
      }

    }
    else //coupling
    {
      std::map<int,std::vector<Epetra_SerialDenseMatrix> >::iterator c = side_coupling.find( coup_sid );

      std::vector<Epetra_SerialDenseMatrix> & side_matrices = c->second;

      if (side_matrices.size() != 3)
        dserror("Obtained only %d side coupling matrices. 3 required.", side_matrices.size());

      // coupling matrices between background element and one! side
      Epetra_SerialDenseMatrix & C_uiu  = side_matrices[0];
      Epetra_SerialDenseMatrix & C_uui  = side_matrices[1];
      Epetra_SerialDenseMatrix & rhC_ui = side_matrices[2];

      // coupling matrices between one side and itself via the element Kss
      std::map<int,std::vector<Epetra_SerialDenseMatrix> >::iterator c2 = Cuiui_coupling.find( coup_sid );
      std::vector<Epetra_SerialDenseMatrix> & Cuiui_matrices = c2->second;
      Epetra_SerialDenseMatrix & eleGsui = Cuiui_matrices[0];
      Epetra_SerialDenseMatrix & eleGuis = Cuiui_matrices[1];

      if (averaging_strategy == INPAR::XFEM::Embedded_Sided or
          averaging_strategy == INPAR::XFEM::Mean) // for coupling-sided and two-sided coupling
        dserror("embedded or two-sided coupling not supported");

      ci[coup_sid] = DRT::ELEMENTS::XFLUID::HybridLMInterface<distype>::CreateHybridLMCoupling_XFluidSided(
          coupl_ele,coupl_xyze,C_uiu,C_uui,rhC_ui,eleGsui,eleGuis,fldparaxfem_->IsViscousAdjointSymmetric());
    }

    if(cond_manager->IsCoupling( coup_sid, my::eid_ ))
    {
      std::map<int, std::vector<int> >::const_iterator k = patchcouplm.find( coup_sid );
      const std::vector<int> & coupl_lm = k->second;

      // set velocity (and pressure) of coupling/slave element at current time step
      ci[coup_sid]->SetSlaveState(*coupl_dis_,coupl_lm);

      // note: old state is handled by nitsche coupling si_nit
    }


    if(!(is_ls_coupling_side and !cond_manager->IsCoupling( coup_sid, my::eid_ ))) // not level-set-WDBC case
    {
      std::map<int, std::vector<int> >::const_iterator k = patchcouplm.find( coup_sid );
      const std::vector<int> & coupl_lm = k->second;

      // add displacement of coupling element at current time step
      ci[coup_sid]->AddSlaveEleDisp(*coupl_dis_,coupl_lm);
    }


    // define interface force vector w.r.t side
    Epetra_SerialDenseVector iforce;
    iforce.Size(cutla[0].lm_.size());

    // we need an instance of Nitsche-evaluation class for evaluation of
    // inflow terms and for evaluation of terms for the previous time step
    // (new OST)
    if (add_conv_stab || my::fldparatimint_->IsNewOSTImplementation())
    {
      if(is_mesh_coupling_side)
      {
        // create si_nit based on side and side_xyze and not w.r.t coup_ele, as no derivatives are used in these coupling terms
        // see also HybridLM_CreateSpecialContributionMatrices()
        if (cond_manager->IsCoupling( coup_sid, my::eid_ )) //... for two-sided problems
        {
          // coupling matrices between background element and one! side
          std::map<int,std::vector<Epetra_SerialDenseMatrix> >::iterator c = side_coupling_extra.find( coup_sid );
          std::vector<Epetra_SerialDenseMatrix> & side_matrices_extra = c->second;
          Epetra_SerialDenseMatrix & C_uiu  = side_matrices_extra[0];
          Epetra_SerialDenseMatrix & C_uui  = side_matrices_extra[1];
          Epetra_SerialDenseMatrix & rhC_ui = side_matrices_extra[2];
          Epetra_SerialDenseMatrix & C_uiui = side_matrices_extra[3];

          si_nit[coup_sid] = DRT::ELEMENTS::XFLUID::NitscheInterface<distype>::CreateNitscheCoupling_XFluidSided(
              side, side_xyze, elemat1_epetra, C_uiu, C_uui, C_uiui, elevec1_epetra, rhC_ui, *fldparaxfem_);
        }
        else
        {
          si_nit[coup_sid] = DRT::ELEMENTS::XFLUID::NitscheInterface<distype>::CreateNitscheCoupling_XFluidWDBC(
              side, side_xyze, elemat1_epetra, elevec1_epetra, *fldparaxfem_);
        }

        // set velocity for current time step
        si_nit[coup_sid]->SetSlaveState(*cutter_dis,cutla[0].lm_);

        // set displacement of side for current time step
        si_nit[coup_sid]->AddSlaveEleDisp(*cutter_dis,cutla[0].lm_);

      }
      else if(is_ls_coupling_side)
      {
        if (cond_manager->IsCoupling( coup_sid, my::eid_ )) //... for two-sided problems
        {
          dserror("convective terms for hybrid lm coupling not implemented yet for level-set cuts");
        }
        else
        {
          si_nit[coup_sid] = DRT::ELEMENTS::XFLUID::NitscheInterface<distype>::CreateNitscheCoupling_XFluidWDBC(
            elemat1_epetra,elevec1_epetra,*fldparaxfem_);
        }
      }
      else dserror("no mesh-/level-set coupling object for coupling sid %i", coup_sid);
    }

    // Set State for current and previous time
    if(cond_manager->IsCoupling( coup_sid, my::eid_ ))
    {
      if (my::fldparatimint_->IsNewOSTImplementation())
      {
        // set velocity for previous time step
        si_nit[coup_sid]->SetSlaveStaten(*cutter_dis,cutla[0].lm_);
      }
    }


    //--------------------------------------------
    // loop boundary cells w.r.t current cut side
    //--------------------------------------------
    for ( std::vector<DRT::UTILS::GaussIntegration>::const_iterator i=cutintpoints.begin();
          i!=cutintpoints.end();
          ++i )
    {
      const DRT::UTILS::GaussIntegration & gi = *i;
      GEO::CUT::BoundaryCell * bc = bcs[i - cutintpoints.begin()]; // get the corresponding boundary cell

      //--------------------------------------------
      // loop gausspoints w.r.t current boundary cell
      //--------------------------------------------
      for ( DRT::UTILS::GaussIntegration::iterator iquad=gi.begin(); iquad!=gi.end(); ++iquad )
      {
        double drs = 0.0; // transformation factor between reference cell and linearized boundary cell

        const LINALG::Matrix<2,1> eta( iquad.Point() ); // xi-coordinates with respect to side

        LINALG::Matrix<3,1> rst(true); // local coordinates w.r.t background element

        LINALG::Matrix<3,1> x_gp_lin(true); // gp in xyz-system on linearized interface

        // compute transformation factor, normal vector and global Gauss point coordiantes
        if(bc->Shape() != DRT::Element::dis_none) // Tessellation approach
        {
          ComputeSurfaceTransformation(drs, x_gp_lin, normal, bc, eta);
        }
        else // MomentFitting approach
        {
          drs = 1.0;
          normal = bc->GetNormalVector();
          const double* gpcord = iquad.Point();
          for (int idim=0;idim<3;++idim)
          {
            x_gp_lin(idim,0) = gpcord[idim];
          }
        }

        // find element local position of gauss point
        GEO::CUT::Position<distype> pos( my::xyze_, x_gp_lin );
        pos.Compute();
        rst = pos.LocalCoordinates();

        //TODO unify ProjectOnSide and Evaluate for different spatial dimensions of boundary slave element and volumetric slave element
        if (is_mesh_coupling_side)
        {
          // project gaussian point from linearized interface to warped side (get/set local side coordinates in SideImpl)
          LINALG::Matrix<3,1> xi_side(true);
          // project on boundary element
          si->ProjectOnSide(x_gp_lin, x_side, xi_side);

          if (averaging_strategy == INPAR::XFEM::Embedded_Sided or
              averaging_strategy == INPAR::XFEM::Mean)
            dserror("embedded or two-sided weighting not supported"); // evaluate embedded element's shape functions at gauss-point coordinates
          else
          {
            ci.at(coup_sid)->Evaluate( xi_side ); // evaluate side's shape functions at gauss-point coordinates
            if(add_conv_stab || my::fldparatimint_->IsNewOSTImplementation())
              si_nit.at(coup_sid)->Evaluate( xi_side ); // evaluate side's shape functions at gauss-point coordinates
          }
        }
        else if (is_ls_coupling_side)
        {
          if(cond_manager->IsCoupling( coup_sid, my::eid_ ))
            dserror("coupling for level-sets not supported here"); // evaluate embedded element's shape functions at gauss-point coordinates
        }


        // integration factors
        const double surf_fac = drs*iquad.Weight();
        const double timefacfac = surf_fac * my::fldparatimint_->TimeFac();

        //--------------------------------------------

        // evaluate shape functions (and derivatives)

        if (assemble_iforce) // evaluate derivatives to assemble iforce vector
        {
          EvalFuncAndDeriv( rst );
        }
        else
        {
          DRT::UTILS::shape_function<distype>( rst, my::funct_ );
        }


        // get velocity at integration point
        // (values at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
        my::velint_.Multiply(evelaf,my::funct_);


        //-----------------------------------------------------------------------------
        // define the prescribed interface jump vectors for velocity and traction

        LINALG::Matrix<my::nsd_,1> ivelint_jump (true);
        LINALG::Matrix<my::nsd_,1> itraction_jump (true);
        LINALG::Matrix<my::nsd_,my::nsd_> proj_tangential (true);
        LINALG::Matrix<my::nsd_,my::nsd_> LB_proj_matrix (true);

        double kappa_m=0.0;
        double kappa_s=0.0;
        double visc_m=0.0;
#ifdef DEBUG
        // Only Navier Slip used kappa_m,kappa_s and visc_m defined just before!
        // To use Navier Slip specify them correct!
        if(cond_type == INPAR::XFEM::CouplingCond_SURF_NAVIER_SLIP || cond_type == INPAR::XFEM::CouplingCond_LEVELSET_NAVIER_SLIP)
          dserror("ElementXfemInterfaceHybridLM with Navier Slip, what to do with kappa_m/kappa_s for the dyn_visc in the traction_jump?");
#endif

        GetInterfaceJumpVectors(
            coupcond,
            coupling,
            ivelint_jump,
            itraction_jump,
            proj_tangential,
            LB_proj_matrix,
            x_gp_lin,
            normal,
            si,
            rst,
            kappa_m,
            kappa_s,
            visc_m
        );

        if(cond_type == INPAR::XFEM::CouplingCond_LEVELSET_NEUMANN or
            cond_type == INPAR::XFEM::CouplingCond_SURF_NEUMANN)
        {
          //-----------------------------------------------------------------------------
          // evaluate the Neumann boundary condition term
          EvaluateNeumann(
              timefacfac,             ///< theta*dt
              my::funct_,             ///< coupling master shape functions
              itraction_jump,         ///< prescribed interface traction, jump height for coupled problems
              elevec1_epetra          ///< element rhs vector
          );

          if (my::fldparatimint_->IsNewOSTImplementation())
          {
            dserror("how to deal with Neumann boundary condition and new OSTImplementation");
          }
        }
        else // standard Hybrid lm terms
        {
          //--------------------------------------------

          bK_ss.MultiplyNT( my::funct_, my::funct_ );

          /*                      \
            - |  (virt tau) * n^f , Du  |
               \                      */

          HybridLM_Evaluate_SurfBased(
              ci[coup_sid],
              bK_ss,
              K_su,
              K_us,
              rhs_s,
              epreaf,
              K_uu,
              rhs_uu,
              G_up,
              G_pu,
              rhs_up,
              rhs_pu,
              normal,
              timefacfac,
              ivelint_jump,
              itraction_jump,
              cond_manager->IsCoupling( coup_sid, my::eid_ ),
              is_MHVS);

          //--------------------------------------------
          // evaluate additional inflow/convective stabilization terms

          if (add_conv_stab)
          {
            double NIT_full_stab_fac = 0.0;
            const double NIT_visc_stab_fac = 0.0;

            my::SetConvectiveVelint(ele->IsAle());


            //-----------------------------------------------------------------------------

            // define the coupling between two not matching grids
            // for fluidfluidcoupling
            // domain Omega^m := Coupling master (XFluid)
            // domain Omega^s := Alefluid( or monolithic: structure) ( not available for non-coupling (Dirichlet) )

            // [| v |] := vm - vs


            if (cond_manager->IsCoupling( coup_sid, my::eid_ ))
            {
              LINALG::Matrix<my::nsd_,1> velint_s;
              ci[coup_sid]->GetInterfaceVelnp(velint_s);


              //Get Material parameters for the master side!
              GetMaterialParametersVolumeCell(mat,densaf_master_,viscaf_master_,gamma_m_);

              bool non_xfluid_coupling;
              double kappa_m;
              double kappa_s;
              cond_manager->GetAverageWeights(coup_sid, ele, kappa_m, kappa_s, non_xfluid_coupling);

              NIT_Compute_FullPenalty_Stabfac(
                NIT_full_stab_fac,  ///< to be filled: full Nitsche's penalty term scaling (viscous+convective part)
                normal,
                h_k,
                kappa_m, //weights (only existing for Nitsche currently!!)
                kappa_s, //weights (only existing for Nitsche currently!!)
                my::convvelint_,
                velint_s,
                NIT_visc_stab_fac   ///< Nitsche's viscous scaling part of penalty term
              );

              si_nit.at(coup_sid)->ApplyConvStabTerms(
                ci[coup_sid],
                my::funct_,
                my::velint_,
                normal,
                my::densaf_,  //Look into this term when changing HybridLM
                NIT_full_stab_fac,
                timefacfac,
                ivelint_jump,
                cond_type
              );
            }
            else // non-coupling
            {
              LINALG::Matrix<my::nsd_,1> velint_s;
              ci[coup_sid]->GetInterfaceVelnp(velint_s);

              //Get Material parameters for the master side!
              GetMaterialParametersVolumeCell(mat,densaf_master_,viscaf_master_,gamma_m_);

              bool non_xfluid_coupling;
              double kappa_m;
              double kappa_s;
              cond_manager->GetAverageWeights(coup_sid, ele, kappa_m, kappa_s, non_xfluid_coupling);

              NIT_Compute_FullPenalty_Stabfac(
                NIT_full_stab_fac,  ///< to be filled: full Nitsche's penalty term scaling (viscous+convective part)
                normal,
                h_k,
                kappa_m, //weights (only existing for Nitsche currently!!)
                kappa_s, //weights (only existing for Nitsche currently!!)
                my::convvelint_,
                velint_s,
                NIT_visc_stab_fac   ///< Nitsche's viscous scaling part of penalty term
              );

              si_nit.at(coup_sid)->ApplyConvStabTerms(
                ci[coup_sid],
                my::funct_,
                my::velint_,
                normal,
                my::densaf_,
                NIT_full_stab_fac,
                timefacfac,
                ivelint_jump,
                cond_type
              );
            } // if coupling
          } // if add conv stab

          if (my::fldparatimint_->IsNewOSTImplementation())
          {
            dserror("New OST for HybridLM not implemented - check out the code below this dserror!");
//            // get velocity at integration point
//            // (values at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
//            my::velintn_.Multiply(eveln,my::funct_);
//
//            // get velocity derivatives at integration point
//            // (values at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
//            my::vderxyn_.MultiplyNT(eveln,my::derxy_);
//
//            //-----------------------------------------------------------------------------
//            // evaluate the coupling terms for coupling with current side
//            // (or embedded element through current side)
//            // time step n
//            const double kappa_m = 1.0;
//            const double kappa_s = 0.0;
//            GetMaterialParametersVolumeCell(mat,densaf_master_,viscaf_master_,gamma_m_);
//
//            // REMARK: do not add adjoint and penalty terms at t_n for hybrid LM approach!
//            // (these are Nitsche-terms! find best settings for Nitsche's method first!)
//            LINALG::Matrix<my::nsd_,1> ivelintn_jump (true);
//            LINALG::Matrix<my::nsd_,1> itractionn_jump(true);
//
//            //Get Configuration Map (finally we should modify the configuration map here in a way that it fits hybrid LM approach)
//            std::map<INPAR::XFEM::CoupTerm, std::pair<bool,double> >& hlm_configmap_n = coupling->GetConfigurationmap();
//
//            si_nit.at(coup_sid)->NIT_evaluateCouplingOldState(
//              normal,
//              surf_fac * (my::fldparatimint_->Dt()-my::fldparatimint_->TimeFac()), // scaling of rhs depending on time discretization scheme
//              false,
//              viscaf_master_,              // dynvisc viscosity in background fluid
//              viscaf_slave_,               // dynvisc viscosity in embedded fluid
//              kappa_m,                     // mortaring weighting
//              kappa_s,                     // mortaring weighting
//              my::densn_,                  // fluid density
//              my::funct_,                  // bg shape functions
//              my::derxy_,                  // bg shape function gradient
//              my::vderxyn_,                // bg grad u^n
//              my::funct_.Dot(epren),       // bg p^n
//              my::velintn_,                // bg u^n
//              ivelintn_jump,
//              itractionn_jump,
//              hlm_configmap_n,
//              INPAR::XFEM::PreviousState_only_consistency
//            );
          }
        } // hybrid lm method

        if (!assemble_iforce)
          continue;

        //--------------------------------------------
        // calculate interface forces for XFSI

        // get velocity derivatives at integration point
        // (values at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
        my::vderxy_.MultiplyNT(evelaf,my::derxy_);

        // get pressure at integration point
        // (value at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
        double press = my::funct_.Dot(epreaf);

        //-------------------------------
        // traction vector w.r.t fluid domain, resulting stresses acting on the fluid surface
        // t= (-p*I + 2mu*eps(u))*n^f
        LINALG::Matrix<my::nsd_,1> traction(true);

        BuildTractionVector( traction, press, normal );

        ci[coup_sid]->ComputeInterfaceForce(iforce, traction, surf_fac );

      } // end loop gauss points of boundary cell
    } // end loop boundary cells of side

    if(assemble_iforce)
      AssembleInterfaceForce(mc_fsi->IForcecol(), *cutter_dis, cutla[0].lm_, iforce);

  } // end loop cut sides

  /*--------------------------------------------------------
    build final element and coupling matrices
  --------------------------------------------------------*/

  // compute inverse K_ss^-1
  LINALG::FixedSizeSerialDenseSolver<my::nen_, my::nen_> invsolver;
  invsolver.SetMatrix(invbK_ss);
  invsolver.Invert();

  // the non-diagonal entries (shear stresses) lead to the factor 2 in the
  // K_ss matrix; inversion leads to 1/2 in the matrix blocks of the shear stress
  halfInvbK_ss.Update(0.5, invbK_ss, 0.0);

  invK_ss.AddView(Sigmaxx, Sigmaxx, invbK_ss);
  invK_ss.AddView(Sigmaxy, Sigmaxy, halfInvbK_ss);
  invK_ss.AddView(Sigmaxz, Sigmaxz, halfInvbK_ss);
  invK_ss.AddView(Sigmayy, Sigmayy, invbK_ss);
  invK_ss.AddView(Sigmayz, Sigmayz, halfInvbK_ss);
  invK_ss.AddView(Sigmazz, Sigmazz, invbK_ss);

  // create views
  LINALG::Matrix<my::numdofpernode_*my::nen_,my::numdofpernode_*my::nen_> elemat(elemat1_epetra, true);
  LINALG::Matrix<my::numdofpernode_*my::nen_,1>                           elevec(elevec1_epetra, true);

  // now the matrix products involving the inverse matrix will be computed!

  // REMARK: at this step, the K matrices already include contributions from surface terms G_us, G_su
  LINALG::BlockMatrix<LINALG::Matrix<my::nen_,my::nen_>, my::numdofpernode_,numstressdof_> KusinvKss;

  // (K_us + G_us) K_ss^-1 (MHVS) or G_us K_ss^-1 (MHCS)
  KusinvKss.Multiply(K_us,invK_ss);

  // (K_us + G_us) K_ss^-1 (K_su + G_su) (MHVS) or G_us  K_ss^-1 (K_su + G_su + K_sp) (MHCS)
  LINALG::BlockMatrix<LINALG::Matrix<my::nen_, my::nen_>, my::numdofpernode_, my::numdofpernode_> KusinvKssKsu;
  KusinvKssKsu.Multiply(KusinvKss, K_su);

  // (K_us + G_us) K_ss^-1 rhs_s (MHVS) or G_us K_ss^-1 rhs_s (MHCS)
  LINALG::BlockMatrix<LINALG::Matrix<my::nen_,1>,my::numdofpernode_,1> KusinvKssrhs_s;

  KusinvKssrhs_s.Multiply(KusinvKss,rhs_s);

  //REMARK: The following term goes into the interface element's rhs-vector rhC_ui!
  /*
   *  G_uis  *(K_ss^-1 rhs_s)
   *  |___________________|
   *        |
   *        |
   *    -->side-impl
   */


  /*--------------------------------------------------------
    setup element matrix
  --------------------------------------------------------*/
  /*
   * Matrix of intersected fluid element for VISCOUS stress-based LM
   *  _                                                                          _
   * |                                                                            |
   *             _
   *      K_uu + K_uu - (K_us + G_us) K_ss^-1 (K_su + G_su) | K_up + G_up
   *
   *  ----------------------------------------------------- |------------
   *
   *      K_pu + G_pu                                       | K_pp
   * |_                                                                          _|
   *
   * Matrix of intersected fluid element for CAUCHY stress-based LM
   *
   *  _                                                                          _
   * |                                                                            |
   *
   *      K_uu  - G_us  K_ss^-1 (K_su + G_su) | K_up - G_us K_ss^-1 K_sp
   *
   *  ----------------------------------------|------------
   *
   *      K_pu                                | K_pp
   * |_                                                                          _|
   *
   */

  // -KusinvKssKsu --> -(K_us + G_us) K_ss^-1 (K_su + G_su) (MHVS) or -G_us  K_ss^-1 (K_su + G_su) -G_us K_ss^-1 K_sp (MHCS)

  // complete the background element matrix with the calculated entries:

  // number of row blocks
  const unsigned numbrow = my::nsd_; // no (p,u)-block for MHVS & MHCS
  // number of column blocks
  const unsigned numdofpernode = my::numdofpernode_; // avoid possible linker error on some compilers
  const unsigned numbcol = (is_MHVS ? numbrow : numdofpernode); // we have (u,p)-block for MHCS (-G_us K_ss^-1 K_sp)

  // loop over row blocks (only velocities)
  for (unsigned ibr=0; ibr<numbrow; ++ibr)
  {
    // loop over column blocks (only velocities for MHVS!)
    for (unsigned ibc=0; ibc<numbcol; ++ibc)
    {
      // add -KusinvKssKsu
      if ( KusinvKssKsu.IsUsed(ibr,ibc) )
      {
        LINALG::Matrix<my::nen_, my::nen_> & bKusinvKssKsu = * KusinvKssKsu(ibr,ibc);

        for (int ir=0; ir<my::nen_; ++ir)
        {
          // row position in the final element matrix
          unsigned row = ibr + ir * my::numdofpernode_;

          for (int ic = 0; ic<my::nen_; ++ic)
          {
            // column position in the final element matrix
            unsigned col = ibc + ic * my::numdofpernode_;

            // - (K_us + G_us) K_ss^-1 (K_su + G_su ) (MHVS) or
            // - G_us  K_ss^-1 (K_su + G_su + K_sp) (MHCS)
            elemat(row, col) -= bKusinvKssKsu(ir,ic);
          }
        }
      } // -KusinvKssKsu

      if (!is_MHVS) continue;

      /*
       * ONLY MHVS:
       * add contribution from viscous term scaled by (-1/n) (n: MHVS-parameter):
       * _
       * K_uu = (-1/n) * K_uu^{viscous}
       *
       * */

      if (K_uu.IsUsed(ibr,ibc))
      {
        LINALG::Matrix<my::nen_, my::nen_> & bK_uu = * K_uu(ibr,ibc);

        for (int ir=0; ir<my::nen_; ++ir)
        {
          // row position in the final element matrix
          unsigned velrow = ibr + ir * my::numdofpernode_;

          for (int ic=0; ic<my::nen_; ++ic)
          {
            // column position in the final element matrix
            unsigned velcol = ibc + ic * my::numdofpernode_;

            // + K_uu
            elemat(velrow, velcol) += bK_uu(ir,ic);
          }
        }
      } // K_uu
    }//end column block loop

    if (!is_MHVS) continue;

    // ONLY MHVS: velocity-pressure coupling entries G_up
    // loop over row blocks
    if ( G_up.IsUsed(ibr,0) )
    {
      LINALG::Matrix<my::nen_, my::nen_> & bGup = * G_up(ibr, 0);

      for (int ic=0; ic<my::nen_; ++ic)
      {
        // column position in final element matrix
        unsigned prescol = my::nsd_ + ic * my::numdofpernode_;

        // loop over rows of velocity-pressure submatrix
        for (int ir=0; ir<my::nen_; ++ir)
        {
          // row position in the final element matrix
          unsigned velrow = ibr + ir * my::numdofpernode_;

          // + G_up
          elemat(velrow, prescol) += bGup(ir, ic);
        }
      }
    } // G_up
  }// end row block loop

  if (is_MHVS)
  {
    // ONLY MHVS: pressure-velocity coupling entries G_pu
    // loop over column blocks
    for (unsigned ibc = 0; ibc < numbcol; ++ibc)
    {
      if ( G_pu.IsUsed(0,ibc) )
      {
        LINALG::Matrix<my::nen_, my::nen_> & bGpu = * G_pu(0, ibc);

        // pressure-velocity entries
        for (int ir=0; ir<my::nen_; ++ir)
        {
          // row position in final element matrix
          unsigned presrow = my::nsd_ + ir * my::numdofpernode_;

          for (int ic=0; ic<my::nen_; ++ic)
          {
            // column position in final element matrix
            unsigned velcol = ibc + ic * my::numdofpernode_;

            // + Gpu
            elemat(presrow,velcol) += bGpu(ir,ic);
          }
        }
      } // G_pu
    } // end column block loop
  }
  //element matrix complete!

  /*--------------------------------------------------------
    setup rhs-vector
  --------------------------------------------------------*/
  /*
   *  for MHVS:
   *  - (K_us + G_us) K_ss^-1 * rhs_s + rhs_uu + rhs_up
   *          +
   *  rhs_pu + rhs_pui
   *  |_______________|
   *        |
   *       united in rhs_pu!
   *
   *  for MHCS:
   *  - G_us K_ss^-1 * rhs_s
   *
   */
  // loop over row blocks
  for (unsigned ibr = 0; ibr<numbrow; ++ibr)
  {
    if ( KusinvKssrhs_s.IsUsed(ibr,0) )
    {
      LINALG::Matrix<my::nen_,1> & bKusinvKssrhs_s = * KusinvKssrhs_s(ibr,0);

      for (int ir=0; ir<my::nen_; ++ir)
      {
        // row position in final element matrix
        unsigned int velrow = ibr + ir * my::numdofpernode_;

        // - (K_us + G_us) K_ss^-1 * rhs_s (MHVS) or
        // - G_us K_ss^-1 * rhs_s (MHCS)
        elevec(velrow,0) -= bKusinvKssrhs_s(ir,0);
      }
    } // -KusinvKssrhs_s

    // ONLY MHVS!
    if (! is_MHVS) continue;

    if ( rhs_uu.IsUsed(ibr,0) )
    {
      LINALG::Matrix<my::nen_,1> & brhs_uu = * rhs_uu(ibr,0);

      for (int ir=0; ir<my::nen_; ++ir)
      {
        // row position in final element matrix
        unsigned velrow = ibr + ir * my::numdofpernode_;

        // + rhs_uu
        elevec(velrow,0) += brhs_uu(ir,0);
      }
    } // rhs_uu

    if ( rhs_up.IsUsed(ibr,0) )
    {
      LINALG::Matrix<my::nen_, 1> & brhs_up = * rhs_up(ibr,0);

      for (int ir=0; ir<my::nen_; ++ir)
      {
        // row position in final element matrix
        unsigned int velrow = ibr + ir * my::numdofpernode_;
        // + rhs_up
        elevec(velrow,0) += brhs_up(ir,0);
      }
    } // rhs_up
  } // end row block loop

  if (is_MHVS)
  {
    // add rhs_pu
    for (int ir=0; ir<my::nen_; ++ir)
    {
      // row position in final element matrix
      unsigned int presrow = my::nsd_ + ir * my::numdofpernode_;
      // + rhs_pu + rhs_pui
      elevec(presrow,0) += rhs_pu(ir,0);
    } // rhs_pu
  }

  // coupling contributions are added matrix & rhs-vector for background element!

  // in case that no coupling objects are available, we are done here
  if (side_coupling.empty())
    return;

  //-------------------------------------------------
  // finalize creation of side coupling terms
  // Cuiu, Cuui, rhCui & Gsui, Guis (-->Cuiui)
  //-------------------------------------------------

  // build interface coupling matrices - therefore iterate through the interface elements
  for (typename std::map<int,Teuchos::RCP<DRT::ELEMENTS::XFLUID::HybridLMInterface<distype> > >::iterator sit=ci.begin();
      sit!=ci.end(); ++sit)
  {
    Teuchos::RCP<DRT::ELEMENTS::XFLUID::HybridLMInterface<distype> > si = sit->second;
    const int coup_sid = sit->first;

    // creation of Cuiu,Cuui,rhCui,Guis and Gsui:
    si->HybridLM_buildFinalCouplingMatrices(invK_ss,KusinvKss,K_su,rhs_s);

    // add contributions from convective stabilization, if active
    if (add_conv_stab)
    {
      std::map<int,std::vector<Epetra_SerialDenseMatrix> >::iterator c = side_coupling.find( coup_sid );
      std::vector<Epetra_SerialDenseMatrix> & side_matrices = c->second;
      std::map<int,std::vector<Epetra_SerialDenseMatrix> >::iterator cc = side_coupling_extra.find( coup_sid );
      std::vector<Epetra_SerialDenseMatrix> & side_matrices_extra = cc->second;
#ifdef DEBUG
      if (side_matrices.size() != 3)
        dserror("Obtained only %d side coupling matrices. 3 required.", side_matrices.size());
      if (side_matrices_extra.size() != 4)
        dserror("Obtained only %d conv. side coupling matrices. 4 required.", side_matrices_extra.size());
#endif

      for (int i=0; i<3; ++i)
      {
#ifdef DEBUG
        if (side_matrices[i].M() != side_matrices_extra[i].M() ||
            side_matrices[i].N() != side_matrices_extra[i].N())
          dserror("Mismatch in matrix dimensions of convective stabilization matrix and MHCS/MHVS coupling matrix");
#endif
        side_matrices[i] += side_matrices_extra[i];
      }

      // in case of new OST and with active convective stab. terms, we have already
      //
      continue;
    }

    // add contributions from old time step to RHS
    if (my::fldparatimint_->IsNewOSTImplementation())
    {
      std::map<int,std::vector<Epetra_SerialDenseMatrix> >::iterator c = side_coupling.find( coup_sid );
      std::vector<Epetra_SerialDenseMatrix> & side_matrices = c->second;
      std::map<int,std::vector<Epetra_SerialDenseMatrix> >::iterator cc = side_coupling_extra.find( coup_sid );
      std::vector<Epetra_SerialDenseMatrix> & side_matrices_extra = cc->second;
#ifdef DEBUG
      if (side_matrices.size() != 3)
        dserror("Obtained only %d side coupling matrices. 3 required.", side_matrices.size());
      if (side_matrices_extra.size() != 4)
        dserror("Obtained only %d conv. side coupling matrices. 4 required.", side_matrices_extra.size());
      if (side_matrices[2].M() != side_matrices_extra[2].M() ||
          side_matrices[2].N() != side_matrices_extra[2].N())
        dserror("Mismatch in matrix dimensions of convective stabilization matrix and MHCS/MHVS coupling matrix");
#endif

      // we only add the RHS contribution
      side_matrices[2] += side_matrices_extra[2];
    }
  } // end loop over map of coupling matrices

  //-------------------------------------------------
  // build C_uiui
  //-------------------------------------------------

  /*
   * "patchelementslm" includes u,v,w,p- DOFs
   * Therefore, the final matrices G_sui and G_uis have pressure columns/rows.
   * The pressure entries are set to zero.
   * The task is now, to build the final G_uis and G_sui matrices
   * out of the element submatrices collected in vector 'Cuiui_matrices'.
   * In there, the G_sui & G_uis contributions from the sides are collected!
   */
  Epetra_SerialDenseMatrix G_sui(numstressdof_*my::nen_,patchelementslm.size());
  Epetra_SerialDenseMatrix G_uis(patchelementslm.size(),numstressdof_*my::nen_);
  Epetra_SerialDenseMatrix Cuiui_conv(patchelementslm.size(),patchelementslm.size());

  // transform the block matrix invK_ss to an EpetraSerialDenseMatrix,
  // to be later multiplied with G_sui & G_uis!
  Epetra_SerialDenseMatrix InvKss(my::nen_*numstressdof_, my::nen_*numstressdof_);

  //--------------------------------------------
  // Build InvKss ( K_ss^(-1) )
  //--------------------------------------------
  for (int ibc=0; ibc<numstressdof_; ++ibc)
  {
    for (int ibr=0; ibr<numstressdof_; ++ibr)
    {
      if (invK_ss.IsUsed(ibr,ibc))
      {
        LINALG::Matrix<my::nen_,my::nen_> & binvK_ss = *invK_ss(ibr,ibc);
        for ( int ic=0; ic<my::nen_; ++ic )
        {
          unsigned col = ibc + numstressdof_ * ic;
          for ( int ir=0; ir<my::nen_; ++ir )
          {
            unsigned row = ibr + numstressdof_ * ir;
            InvKss(row,col) -= binvK_ss(ir,ic);
          }
        }
      }
    }
  }

  //--------------------------------------------
  // assemble G_sui & G_uis
  //--------------------------------------------
  int ipatchsizesbefore = 0;

  for (std::map<int, std::vector<Epetra_SerialDenseMatrix> >::const_iterator m=Cuiui_coupling.begin();
       m!=Cuiui_coupling.end(); ++m)
  {
    const std::vector<Epetra_SerialDenseMatrix> & Cuiui_matrices = m->second;

    // assemble Gsui
    for ( int ibc=0; ibc<Cuiui_matrices[0].N(); ++ibc )
    {
      for ( int ibr=0; ibr<numstressdof_*my::nen_; ++ibr )
      {
        G_sui(ibr,ibc+ipatchsizesbefore) = Cuiui_matrices[0](ibr,ibc);
      }
    }

    // assemble Guis
    for ( int ibc=0; ibc <numstressdof_*my::nen_; ++ibc )
    {
      for ( int ibr=0; ibr<Cuiui_matrices[1].M(); ++ibr )
      {
        G_uis(ibr+ipatchsizesbefore, ibc) = Cuiui_matrices[1](ibr,ibc);
      }
    }

    if (add_conv_stab)
    {
      const int coup_sid = m->first;
      std::map<int,std::vector<Epetra_SerialDenseMatrix> >::const_iterator c = side_coupling_extra.find( coup_sid );
      const std::vector<Epetra_SerialDenseMatrix>  & Cuiui_conv_matrices = c->second;

      for ( int ibc=0; ibc<Cuiui_conv_matrices[3].N(); ++ibc )
      {
        for ( int ibr=0; ibr<Cuiui_conv_matrices[3].M(); ++ibr )
        {
          Cuiui_conv(ibr+ipatchsizesbefore,ibc+ipatchsizesbefore) = Cuiui_conv_matrices[3](ibr,ibc);
        }
      }
    }

    ipatchsizesbefore += Cuiui_matrices[0].N();

  }

  Epetra_SerialDenseMatrix GuisInvKss(patchelementslm.size(),numstressdof_*my::nen_);

  // G_uis * K_ss^-1
  GuisInvKss.Multiply('N', 'N', 1.0, G_uis, InvKss, 1.0);

  // Cuiui <--> (-)G_uis * K_ss^-1 * G_sui
  Cuiui.Multiply('N', 'N', 1.0, GuisInvKss, G_sui, 1.0);

  if (add_conv_stab )
    Cuiui += Cuiui_conv;
}

/*--------------------------------------------------------------------------------
 * setup volume-based terms
 * (mixed/hybrid viscous stress-based LM approach)
 *-------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void FluidEleCalcXFEM<distype>::HybridLM_Build_VolBased(
    const std::vector<DRT::UTILS::GaussIntegration> &                                           intpoints,
    const GEO::CUT::plain_volumecell_set&                                                       cells,
    const LINALG::Matrix<my::nsd_, my::nen_> &                                                  evelaf,            ///< element velocity
    const LINALG::Matrix<my::nen_,1>&                                                           epreaf,            ///< element pressure
    LINALG::Matrix<my::nen_,my::nen_>&                                                          bK_ss,             ///< block K_ss matrix
    LINALG::Matrix<my::nen_,my::nen_>&                                                          invbK_ss,          ///< inverse of block K_ss matrix
    LINALG::BlockMatrix<LINALG::Matrix<my::nen_,my::nen_>,numstressdof_,my::numdofpernode_> &   K_su,              ///< K_su matrix
    LINALG::BlockMatrix<LINALG::Matrix<my::nen_,1>,numstressdof_,1> &                           rhs_s,             ///< rhs_s vector
    LINALG::BlockMatrix<LINALG::Matrix<my::nen_,my::nen_>,my::numdofpernode_,numstressdof_> &   K_us,              ///< K_us matrix
    LINALG::BlockMatrix<LINALG::Matrix<my::nen_,my::nen_>,my::nsd_,my::nsd_> &                  K_uu,              ///< K_uu matrix
    LINALG::BlockMatrix<LINALG::Matrix<my::nen_,1>,my::nsd_, 1> &                               rhs_uu,            ///< rhs_u(u) vector
    const bool                                                                                  is_MHVS,           ///< viscous (true) or Cauchy (false) stress-based LM
    const double                                                                                mhvs_param         ///< stabilizing parameter for viscous stress-based LM
)
{
  // full L2-projection means integration over the full background element,
  // not only the physical part
  if (fldparaxfem_->HybridLM_L2Proj() == INPAR::XFEM::Hybrid_LM_L2_Proj_full)
  {
    // get the standard set of gauss-points from the intersected element
    for (DRT::UTILS::GaussIntegration::const_iterator iquad = my::intpoints_.begin();
        iquad != my::intpoints_.end(); ++iquad)
    {
      my::EvalShapeFuncAndDerivsAtIntPoint(iquad.Point(), iquad.Weight());

      if (is_MHVS)
      {
        MHVS_Evaluate_VolBased(evelaf, bK_ss, invbK_ss, K_su, rhs_s, K_us, K_uu, rhs_uu, mhvs_param);
      }
      else
      {
        MHCS_Evaluate_VolBased(evelaf, epreaf, bK_ss, invbK_ss, K_su, rhs_s);
      }
    }
  }
  else
  {
    for( std::vector<DRT::UTILS::GaussIntegration>::const_iterator i=intpoints.begin();i!=intpoints.end();++i )
    {
      const DRT::UTILS::GaussIntegration intcell = *i;
      for ( DRT::UTILS::GaussIntegration::iterator iquad=intcell.begin(); iquad!=intcell.end(); ++iquad )
      {
        // evaluate shape functions and derivatives at integration point
        my::EvalShapeFuncAndDerivsAtIntPoint(iquad.Point(),iquad.Weight());
        if (is_MHVS)
        {
          MHVS_Evaluate_VolBased(evelaf, bK_ss, invbK_ss, K_su, rhs_s, K_us, K_uu, rhs_uu, mhvs_param);
        }
        else
        {
          MHCS_Evaluate_VolBased(evelaf, epreaf, bK_ss, invbK_ss, K_su, rhs_s);
        }
      }
    }
  }

  return;
}

/*--------------------------------------------------------------------------------
 * evaluate volume-based terms for current gauss point
 * (mixed/hybrid Cauchy stress-based LM approach)
 *-------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void FluidEleCalcXFEM<distype>::MHCS_Evaluate_VolBased(
    const LINALG::Matrix<my::nsd_,my::nen_>&                                                  evelaf,   ///< element velocity
    const LINALG::Matrix<my::nen_,1>&                                                         epreaf,   ///< element pressure
    LINALG::Matrix<my::nen_,my::nen_>&                                                        bK_ss,    ///< block K_ss matrix
    LINALG::Matrix<my::nen_,my::nen_>&                                                        invbK_ss, ///< inverse of block K_ss matrix
    LINALG::BlockMatrix<LINALG::Matrix<my::nen_,my::nen_>,numstressdof_,my::numdofpernode_>&  K_su,     ///< K_su matrix
    LINALG::BlockMatrix<LINALG::Matrix<my::nen_,1>,numstressdof_,1> &                         rhs_s     ///< rhs_s vector
)
{
  LINALG::Matrix<my::nen_,1> dx;
  LINALG::Matrix<my::nen_,1> dy;
  LINALG::Matrix<my::nen_,1> dz;

  LINALG::Matrix<my::nen_,my::nen_> conv_x;
  LINALG::Matrix<my::nen_,my::nen_> conv_y;
  LINALG::Matrix<my::nen_,my::nen_> conv_z;

  //----------------------------------------------------------------------
  // set time-integration factors for left- and right-hand side
  // (two right-hand-side factors: general and for residuals)
  //----------------------------------------------------------------------

  const double viscfac = 1.0/(2.0*my::visceff_);

  // get velocity at integration point
  // (values at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
  my::velint_.Multiply(evelaf,my::funct_);

  // get velocity derivatives at integration point
  // (values at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
  my::vderxy_.MultiplyNT(evelaf,my::derxy_);

  // get pressure at integration point
  // (value at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
  double press = my::funct_.Dot(epreaf);

  // time integration factor & spatial integration factor
  const double timefacfac = my::fldparatimint_->TimeFac() * my::fac_;


  //--------------------------------------------

  for ( int i=0; i<my::nen_; ++i )
  {
    dx( i ) = my::derxy_( 0, i );
    dy( i ) = my::derxy_( 1, i );
    dz( i ) = my::derxy_( 2, i );
  }

  // block - K_ss
  bK_ss.MultiplyNT( my::funct_, my::funct_ );


  conv_x.MultiplyNT( my::funct_, dx );
  conv_y.MultiplyNT( my::funct_, dy );
  conv_z.MultiplyNT( my::funct_, dz );

    /*                     \
  - |  virt tau , eps(Dtau)  |
    \                     */

  invbK_ss.Update( -viscfac*timefacfac, bK_ss, 1.0 );

    /*                 \
   | virt tau , eps(Du) |
    \                 */


  // K_su

  K_su( Sigmaxx, Velx )->Update( timefacfac, conv_x, 1.0 );
  K_su( Sigmaxy, Velx )->Update( timefacfac, conv_y, 1.0 );
  K_su( Sigmayx, Vely )->Update( timefacfac, conv_x, 1.0 );
  K_su( Sigmaxz, Velx )->Update( timefacfac, conv_z, 1.0 );
  K_su( Sigmazx, Velz )->Update( timefacfac, conv_x, 1.0 );
  K_su( Sigmayy, Vely )->Update( timefacfac, conv_y, 1.0 );
  K_su( Sigmayz, Vely )->Update( timefacfac, conv_z, 1.0 );
  K_su( Sigmazy, Velz )->Update( timefacfac, conv_y, 1.0 );
  K_su( Sigmazz, Velz )->Update( timefacfac, conv_z, 1.0 );

  // r_su

  rhs_s( Sigmaxx, 0 )->Update( - timefacfac* my::vderxy_(0, 0)                     , my::funct_, 1.0 );
  rhs_s( Sigmaxy, 0 )->Update( - timefacfac*(my::vderxy_(0, 1) + my::vderxy_(1, 0)), my::funct_, 1.0 );
  rhs_s( Sigmaxz, 0 )->Update( - timefacfac*(my::vderxy_(0, 2) + my::vderxy_(2, 0)), my::funct_, 1.0 );
  rhs_s( Sigmayy, 0 )->Update( - timefacfac* my::vderxy_(1, 1)                     , my::funct_, 1.0 );
  rhs_s( Sigmayz, 0 )->Update( - timefacfac*(my::vderxy_(1, 2) + my::vderxy_(2, 1)), my::funct_, 1.0 );
  rhs_s( Sigmazz, 0 )->Update( - timefacfac* my::vderxy_(2, 2)                     , my::funct_, 1.0 );

  // stressbar-pressure coupling
  /*
                  /                    \
                 |                      |
               - | tr(virt tau^e) , p I |
                 |                      |
                  \                    /
  */

  // K_sp
  K_su( Sigmaxx, Pres )->Update( -viscfac*timefacfac, bK_ss, 1.0 );
  K_su( Sigmayy, Pres )->Update( -viscfac*timefacfac, bK_ss, 1.0 );
  K_su( Sigmazz, Pres )->Update( -viscfac*timefacfac, bK_ss, 1.0 );

  // r_sp
  rhs_s( Sigmaxx, 0 )->Update( viscfac*timefacfac*press, my::funct_, 1.0 );
  rhs_s( Sigmayy, 0 )->Update( viscfac*timefacfac*press, my::funct_, 1.0 );
  rhs_s( Sigmazz, 0 )->Update( viscfac*timefacfac*press, my::funct_, 1.0 );

  return;

} //EvaluateMatricesMSH

/*--------------------------------------------------------------------------------
 * evaluate volume-based matrices K for current gauss point
 * (mixed/hybrid viscous stress-based LM approach)
 *-------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void FluidEleCalcXFEM<distype>::MHVS_Evaluate_VolBased(
    const LINALG::Matrix<my::nsd_,my::nen_> &                                                 evelaf,
    LINALG::Matrix<my::nen_,my::nen_> &                                                       bK_ss,
    LINALG::Matrix<my::nen_,my::nen_> &                                                       invbK_ss,
    LINALG::BlockMatrix<LINALG::Matrix<my::nen_,my::nen_>,numstressdof_,my::numdofpernode_> & K_su,
    LINALG::BlockMatrix<LINALG::Matrix<my::nen_,1>,numstressdof_,1> &                         rhs_s,
    LINALG::BlockMatrix<LINALG::Matrix<my::nen_,my::nen_>,my::numdofpernode_,numstressdof_> & K_us,
    LINALG::BlockMatrix<LINALG::Matrix<my::nen_,my::nen_>,my::nsd_,my::nsd_> &                K_uu,
    LINALG::BlockMatrix<LINALG::Matrix<my::nen_,1>,my::nsd_, 1> &                             rhs_uu,
    const double &                                                                            mhvs_param
)
{
  // velocities at current gauss point
  my::velint_.Multiply(evelaf, my::funct_);

  // velocity gradient at current gauss point
  my::vderxy_.MultiplyNT(evelaf, my::derxy_);

  // compute shape function derivatives:
  // get derivatives of nodal shape function vector w. r. t. x,y,z
  LINALG::Matrix<my::nen_,1> dx;
  LINALG::Matrix<my::nen_,1> dy;
  LINALG::Matrix<my::nen_,1> dz;

  // derxy(coordinate,node)
  for (int i = 0; i < my::nen_; ++i)
  {
    dx(i) = my::derxy_(0,i);
    dy(i) = my::derxy_(1,i);
    dz(i) = my::derxy_(2,i);
  }

  // fill the (nen_ x nen_) matrix block of K_ss
  bK_ss.MultiplyNT(my::funct_, my::funct_); //N * N^T

  // scaling with inverse dynamic effective viscosity
  const double viscfac = -1.0 / (2.0 * my::visceff_);

  // time integration factor & spatial integration factor, scaled with mhvs-parameter 1/n
  const double timefacfac = my::fldparatimint_->TimeFac() * my::fac_ / mhvs_param;


  /* K_ss
   *
   *       1        /                \
   * (-)  ---   *  |  \tau, \sigma    |
   *   (2n \mu)     \                /
   */
  invbK_ss.Update(viscfac * timefacfac, bK_ss, 1.0);

  /*
   *  K_su
   *  from:
   *
   *  1     /                    \
   *  - *  |  \tau, \epsilon(u)   |
   *  n     \                    /
   *
   */

  // create blocks
  LINALG::Matrix<my::nen_,my::nen_> NdNdxT;
  LINALG::Matrix<my::nen_,my::nen_> NdNdyT;
  LINALG::Matrix<my::nen_,my::nen_> NdNdzT;

  NdNdxT.MultiplyNT(my::funct_, dx);
  NdNdyT.MultiplyNT(my::funct_, dy);
  NdNdzT.MultiplyNT(my::funct_, dz);

  // add main diagonal blocks
  K_su(Sigmaxx, Velx)->Update(timefacfac, NdNdxT, 1.0);
  K_su(Sigmayy, Vely)->Update(timefacfac, NdNdyT, 1.0);
  K_su(Sigmazz, Velz)->Update(timefacfac, NdNdzT, 1.0);

  // add off-diagonal blocks
  K_su(Sigmaxy, Velx)->Update(timefacfac, NdNdyT, 1.0);
  K_su(Sigmaxz, Velx)->Update(timefacfac, NdNdzT, 1.0);
  K_su(Sigmayx, Vely)->Update(timefacfac, NdNdxT, 1.0);
  K_su(Sigmayz, Vely)->Update(timefacfac, NdNdzT, 1.0);
  K_su(Sigmazx, Velz)->Update(timefacfac, NdNdxT, 1.0);
  K_su(Sigmazy, Velz)->Update(timefacfac, NdNdyT, 1.0);

  /*
   * rhs_s contribution
   * from:
   *
   *  1     /                   \
   *  -  * |  \tau, \epsilon(u)  |
   *  n     \                   /
   *
   */

  rhs_s(Sigmaxx, 0)->Update(-timefacfac * my::vderxy_(0,0),my::funct_, 1.0);
  rhs_s(Sigmaxy, 0)->Update(-timefacfac * (my::vderxy_(1,0) + my::vderxy_(0,1)),my::funct_, 1.0);
  rhs_s(Sigmaxz, 0)->Update(-timefacfac * (my::vderxy_(0,2) + my::vderxy_(2,0)),my::funct_, 1.0);
  rhs_s(Sigmayy, 0)->Update(-timefacfac * my::vderxy_(1,1),my::funct_,1.0);
  rhs_s(Sigmayz, 0)->Update(-timefacfac * (my::vderxy_(1,2) + my::vderxy_(2,1)),my::funct_, 1.0);
  rhs_s(Sigmazz, 0)->Update(-timefacfac * my::vderxy_(2,2),my::funct_, 1.0);

  // nothing like K_sp here, as this is a purely viscous stress-based approach

  /*
   * coupling matrix K_us which results from testing the LM-stress field with
   * test strains computed from the test velocities
   *
   *     1     /                         \
   *     -  * |  \epsilon(v), \sigma      | * \alpha
   *     n     \                         /
   *
   */

  // as the viscous stresses are condensed, there is no contribution of this term to
  // the RHS vector
  LINALG::Matrix<my::nen_,my::nen_> dNdxNT;
  LINALG::Matrix<my::nen_,my::nen_> dNdyNT;
  LINALG::Matrix<my::nen_,my::nen_> dNdzNT;

  dNdxNT.UpdateT(NdNdxT);
  dNdyNT.UpdateT(NdNdyT);
  dNdzNT.UpdateT(NdNdzT);

  // leads to terms, that are analogous to a symmetric/non-symmetric Nitsche-formulation
  // REMARK: behaves unstable for betau=-1.0 in fluid-fluid problems, so keep that in mind!
  // betau (-)1 <--> symmetric Nitsche
  // betau (+)1 <--> non-symmetric Nitsche
  // the stabilizing parameter n has been applied to timefacfac

  const double alpha = fldparaxfem_->IsViscousAdjointSymmetric() ? 1.0 : -1.0;

  // add main diagonal submatrices
  K_us(Velx, Sigmaxx)->Update(alpha * timefacfac, dNdxNT, 1.0);
  K_us(Vely, Sigmayy)->Update(alpha * timefacfac, dNdyNT, 1.0);
  K_us(Velz, Sigmazz)->Update(alpha * timefacfac, dNdzNT, 1.0);

  // add off-diagonal blocks
  K_us(Velx, Sigmaxy)->Update(alpha * timefacfac, dNdyNT, 1.0);
  K_us(Velx, Sigmaxz)->Update(alpha * timefacfac, dNdzNT, 1.0);
  K_us(Vely, Sigmayx)->Update(alpha * timefacfac, dNdxNT, 1.0);
  K_us(Vely, Sigmayz)->Update(alpha * timefacfac, dNdzNT, 1.0);
  K_us(Velz, Sigmazx)->Update(alpha * timefacfac, dNdxNT, 1.0);
  K_us(Velz, Sigmazy)->Update(alpha * timefacfac, dNdyNT, 1.0);

  // computation of additional stress term, scaled with inverse MHVS-parameter

  // build K_uu coupling matrix
  /*
   *   /                           \     1
   * -|  \epsilon(v), \epsilon(u)   | *  - * 2 \mu * \alpha
   *   \                           /     n
   *
   */

  // factor 2 from above is cancelled out
  const double visc_timefac_mhvs = -alpha * my::visceff_ * timefacfac;

  std::vector<const LINALG::Matrix<my::nen_,1> *> dN;
  dN.push_back(&dx);
  dN.push_back(&dy);
  dN.push_back(&dz);

  LINALG::Matrix<my::nen_,my::nen_> dNidxj(true);
  LINALG::Matrix<my::nen_,my::nen_> dNjdxj(true);

  for (int idim = 0; idim < my::nsd_; ++idim)
  {
    for (int jdim = 0; jdim < my::nsd_; ++jdim)
    {
      dNidxj.MultiplyNT(*dN[jdim],*dN[idim]);
      K_uu(idim, jdim)->Update(visc_timefac_mhvs, dNidxj, 1.0);
      dNjdxj.MultiplyNT(*dN[jdim],*dN[jdim]);
      K_uu(idim, idim)->Update(visc_timefac_mhvs, dNjdxj, 1.0);
      rhs_uu(idim, 0)->Update(-visc_timefac_mhvs * (my::vderxy_(idim, jdim) + my::vderxy_(jdim, idim)), *dN[jdim], 1.0);
    }
  }
}

/*--------------------------------------------------------------------------------
 * build surface-based terms for current gauss point
 * (mixed/hybrid Cauchy or viscous stress-based LM approach)
 *-------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void FluidEleCalcXFEM<distype>::HybridLM_Evaluate_SurfBased(
    Teuchos::RCP<DRT::ELEMENTS::XFLUID::HybridLMInterface<distype> > &                        si,
    const LINALG::Matrix<my::nen_,my::nen_> &                                                 bK_ss,
    LINALG::BlockMatrix<LINALG::Matrix<my::nen_,my::nen_>,numstressdof_,my::numdofpernode_> & K_su,
    LINALG::BlockMatrix<LINALG::Matrix<my::nen_,my::nen_>,my::numdofpernode_,numstressdof_> & K_us,
    LINALG::BlockMatrix<LINALG::Matrix<my::nen_,1>,numstressdof_,1> &                         rhs_s,
    const LINALG::Matrix<my::nen_,1> &                                                        epreaf,
    LINALG::BlockMatrix<LINALG::Matrix<my::nen_,my::nen_>,my::nsd_,my::nsd_> &                K_uu,
    LINALG::BlockMatrix<LINALG::Matrix<my::nen_,1>,my::nsd_, 1> &                             rhs_uu,
    LINALG::BlockMatrix<LINALG::Matrix<my::nen_,my::nen_>,my::nsd_,1> &                       G_up,
    LINALG::BlockMatrix<LINALG::Matrix<my::nen_,my::nen_>,1,my::nsd_> &                       G_pu,
    LINALG::BlockMatrix<LINALG::Matrix<my::nen_,1>,my::nsd_,1> &                              rhs_up,
    LINALG::Matrix<my::nen_,1> &                                                              rhs_pu,
    const LINALG::Matrix<my::nsd_,1> &                                                        normal,
    const double &                                                                            timesurffac,
    const LINALG::Matrix<my::nsd_,1> &                                                        ivelint_jump,
    const LINALG::Matrix<my::nsd_,1> &                                                        itraction_jump,
    const bool                                                                                eval_side_coupling,
    const bool                                                                                is_MHVS
)
{
  K_us(Velx, Sigmaxx)->Update(-timesurffac * normal(Velx), bK_ss, 1.0);
  K_us(Velx, Sigmaxy)->Update(-timesurffac * normal(Vely), bK_ss, 1.0);
  K_us(Velx, Sigmaxz)->Update(-timesurffac * normal(Velz), bK_ss, 1.0);
  K_us(Vely, Sigmayx)->Update(-timesurffac * normal(Velx), bK_ss, 1.0);
  K_us(Vely, Sigmayy)->Update(-timesurffac * normal(Vely), bK_ss, 1.0);
  K_us(Vely, Sigmayz)->Update(-timesurffac * normal(Velz), bK_ss, 1.0);
  K_us(Velz, Sigmazx)->Update(-timesurffac * normal(Velx), bK_ss, 1.0);
  K_us(Velz, Sigmazy)->Update(-timesurffac * normal(Vely), bK_ss, 1.0);
  K_us(Velz, Sigmazz)->Update(-timesurffac * normal(Velz), bK_ss, 1.0);


  // K_su - add the blocks from the surface contribution
  // viscous stress-tested interface continuity term
  /*
  *
  *     /                             \
  *   - |       \tau_{ij} * n_j, u^i   |
  *     \                             /
  */
  K_su(Sigmaxx, Velx)->Update(-timesurffac * normal(Velx), bK_ss, 1.0);
  K_su(Sigmaxy, Velx)->Update(-timesurffac * normal(Vely), bK_ss, 1.0);
  K_su(Sigmaxz, Velx)->Update(-timesurffac * normal(Velz), bK_ss, 1.0);
  K_su(Sigmayx, Vely)->Update(-timesurffac * normal(Velx), bK_ss, 1.0);
  K_su(Sigmayy, Vely)->Update(-timesurffac * normal(Vely), bK_ss, 1.0);
  K_su(Sigmayz, Vely)->Update(-timesurffac * normal(Velz), bK_ss, 1.0);
  K_su(Sigmazx, Velz)->Update(-timesurffac * normal(Velx), bK_ss, 1.0);
  K_su(Sigmazy, Velz)->Update(-timesurffac * normal(Vely), bK_ss, 1.0);
  K_su(Sigmazz, Velz)->Update(-timesurffac * normal(Velz), bK_ss, 1.0);

  //Add surface integral contribution to rhs_s

  // from diagonal terms
  rhs_s(Sigmaxx, 0)->Update(timesurffac * normal(Velx) * my::velint_(Velx), my::funct_, 1.0);
  rhs_s(Sigmayy, 0)->Update(timesurffac * normal(Vely) * my::velint_(Vely), my::funct_, 1.0);
  rhs_s(Sigmazz, 0)->Update(timesurffac * normal(Velz) * my::velint_(Velz), my::funct_, 1.0);

  // from off-diagonal terms
  rhs_s(Sigmaxy, 0)->Update(timesurffac * (normal(Vely) * my::velint_(Velx) + normal(Velx) * my::velint_(Vely)), my::funct_, 1.0);
  rhs_s(Sigmaxz, 0)->Update(timesurffac * (normal(Velz) * my::velint_(Velx) + normal(Velx) * my::velint_(Velz)), my::funct_, 1.0);
  rhs_s(Sigmayz, 0)->Update(timesurffac * (normal(Velz) * my::velint_(Vely) + normal(Vely) * my::velint_(Velz)), my::funct_, 1.0);

  // get pressure at current integration point
  double press = my::funct_.Dot(epreaf);

  // MHVS terms
  if (is_MHVS)
  {
    // interface pressure term
    /*
    *
    *     /          \
    *    |   v, p n   |
    *     \          /
    */
    G_up(Velx, 0)->Update(timesurffac * normal(0), bK_ss, 1.0);
    G_up(Vely, 0)->Update(timesurffac * normal(1), bK_ss, 1.0);
    G_up(Velz, 0)->Update(timesurffac * normal(2), bK_ss, 1.0);

    // velocity residual rhs_up
    rhs_up(Velx,0)->Update(-timesurffac * normal(Velx) * press, my::funct_, 1.0);
    rhs_up(Vely,0)->Update(-timesurffac * normal(Vely) * press, my::funct_, 1.0);
    rhs_up(Velz,0)->Update(-timesurffac * normal(Velz) * press, my::funct_, 1.0);

    // pressure-tested interface continuity term
    /*
    *
    *      /         \
    *    -|  q n, u   |
    *      \         /
    */
    G_pu(0, Velx)->Update(-timesurffac * normal(0), bK_ss, 1.0);
    G_pu(0, Vely)->Update(-timesurffac * normal(1), bK_ss, 1.0);
    G_pu(0, Velz)->Update(-timesurffac * normal(2), bK_ss, 1.0);

    // pressure residual rhs_pu
    // this results from -(q, u_i * n_i)_{\Gamma} (pressure-tested kinematic continuity)
    const double normalvel = my::velint_.Dot(normal);
    rhs_pu.Update(timesurffac*normalvel, my::funct_, 1.0);
  }

  // the terms involving side-DOF are treated by the side implementation class,
  // as the sides have their own shape functions!
  if (eval_side_coupling)
  {
    if (is_MHVS)
      si->MHVS_buildCouplingMatrices(normal, timesurffac, my::funct_, rhs_s, press, rhs_pu, ivelint_jump, itraction_jump);
    else
      si->MHCS_buildCouplingMatrices(normal, timesurffac, my::funct_, rhs_s, ivelint_jump, itraction_jump);
  }
  else
  {
    // from velocity jump tested with viscous test stresses
    /*
    *
    *     /                    _     \
    *   +|  \tau_{ij} * n_j,   u^i    |
    *     \                          /
    */

    // add surface integral contribution to rhs_s

    // from diagonal terms
    rhs_s(Sigmaxx, 0)->Update(-timesurffac * normal(Velx) * ivelint_jump(Velx), my::funct_, 1.0);
    rhs_s(Sigmayy, 0)->Update(-timesurffac * normal(Vely) * ivelint_jump(Vely), my::funct_, 1.0);
    rhs_s(Sigmazz, 0)->Update(-timesurffac * normal(Velz) * ivelint_jump(Velz), my::funct_, 1.0);

    // from off-diagonal terms
    rhs_s(Sigmaxy, 0)->Update(-timesurffac * (normal(Vely) * ivelint_jump(Velx) + normal(Velx) * ivelint_jump(Vely)), my::funct_, 1.0);
    rhs_s(Sigmaxz, 0)->Update(-timesurffac * (normal(Velz) * ivelint_jump(Velx) + normal(Velx) * ivelint_jump(Velz)), my::funct_, 1.0);
    rhs_s(Sigmayz, 0)->Update(-timesurffac * (normal(Velz) * ivelint_jump(Vely) + normal(Vely) * ivelint_jump(Velz)), my::funct_, 1.0);

    if (! is_MHVS) return;

    // ONLY MHVS:
    // pressure-tested kinematic continuity term
    /*
     *
     *      /      _  \
     *    +|  q n, u   |
     *      \         /
     */
    const double normalvel = ivelint_jump.Dot(normal);

    rhs_pu.Update(-timesurffac * normalvel, my::funct_, 1.0);
  }
}

/*--------------------------------------------------------------------------------
 *--------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void FluidEleCalcXFEM<distype>::ElementXfemInterfaceNIT(
    DRT::ELEMENTS::Fluid *                                              ele,
    DRT::Discretization &                                               dis,
    const std::vector<int> &                                            lm,
    const Teuchos::RCP<XFEM::ConditionManager> &                        cond_manager,      ///< XFEM condition manager
    const std::map<int, std::vector<GEO::CUT::BoundaryCell*> > &        bcells,
    const std::map<int, std::vector<DRT::UTILS::GaussIntegration> > &   bintpoints,
    const std::map<int, std::vector<int> > &                            patchcouplm,       ///< lm vectors for coupling elements, key= global coupling side-Id
    Teuchos::ParameterList&                                             params,
    Teuchos::RCP<MAT::Material>&                                        mat_master,       ///< material for the background
    Teuchos::RCP<MAT::Material>&                                        mat_slave,        ///< material for the coupled side
    Epetra_SerialDenseMatrix&                                           elemat1_epetra,
    Epetra_SerialDenseVector&                                           elevec1_epetra,
    const GEO::CUT::plain_volumecell_set &                              vcSet,
    std::map<int, std::vector<Epetra_SerialDenseMatrix> > &             side_coupling,
    Epetra_SerialDenseMatrix&                                           Cuiui
  )
{
#ifdef DEBUG
  if(cond_manager == Teuchos::null) dserror("set the condition manager!");
#endif


  //----------------------------------------------------------------------------
  //                         ELEMENT GEOMETRY
  //----------------------------------------------------------------------------


  // ---------------------------------------------------------------------
  // get initial node coordinates for element
  // ---------------------------------------------------------------------
  // get node coordinates
  GEO::fillInitialPositionArray< distype, my::nsd_, LINALG::Matrix<my::nsd_,my::nen_> >( ele, my::xyze_ );

  // ---------------------------------------------------------------------
  // get additional state vectors for ALE case: grid displacement and vel.
  // ---------------------------------------------------------------------

  my::edispnp_.Clear();
  my::egridv_.Clear();

  if (ele->IsAle()) my::GetGridDispVelALE(dis, lm, my::edispnp_, my::egridv_);


  // ---------------------------------------------------------------------

  /// element coordinates in EpetraMatrix
  Epetra_SerialDenseMatrix ele_xyze(my::nsd_,my::nen_);
  for ( int i=0; i<my::nen_; ++i )
  {
    for(int j=0; j<my::nsd_; ++j)
      ele_xyze(j,i) = my::xyze_( j, i );
  }

  // ---------------------------------------------------------------------
  // get velocity state vectors
  // ---------------------------------------------------------------------

  // get element-wise velocity/pressure field for current time step
  evelaf_.Clear();
  epreaf_.Clear();
  my::ExtractValuesFromGlobalVector(dis, lm, *my::rotsymmpbc_, &evelaf_, &epreaf_, "velaf");

  // get element-wise velocity/pressure field for previous time step
  eveln_.Clear();
  epren_.Clear();
  if (my::fldparatimint_->IsNewOSTImplementation())
    my::ExtractValuesFromGlobalVector(dis, lm, *my::rotsymmpbc_, &eveln_, &epren_, "veln");

  // ---------------------------------------------------------------------
  // set element advective field for Oseen problems
  // ---------------------------------------------------------------------
  if (my::fldpara_->PhysicalType()==INPAR::FLUID::oseen) my::SetAdvectiveVelOseen(ele);


  //-----------------------------------------------------------------------------------
  //                     application-specific flags & parameters
  //-----------------------------------------------------------------------------------

  // map of boundary element gids, to coupling matrices Cuiui
  std::map<int, std::vector<Epetra_SerialDenseMatrix> > Cuiui_coupling;

  //-----------------------------------------------------------------------------------
  //            preparation of Cuiui-coupling matrices for each side
  //-----------------------------------------------------------------------------------

  // create Cuiui coupling matrices for each coupling side (we don't have Cuiui for standard Dirichlet problems...)

  // loop all the intersecting sides of actele
  for (std::map<int,  std::vector<GEO::CUT::BoundaryCell*> >::const_iterator bc=bcells.begin();
      bc!=bcells.end(); ++bc )
  {
    const int coup_sid = bc->first;

    if(!cond_manager->IsCoupling( coup_sid, my::eid_ )) continue; // no couplings to be evaluated for current side

    // get coupling matrices for the current side (boundary element)
    std::vector<Epetra_SerialDenseMatrix> & Cuiui_matrices = Cuiui_coupling[coup_sid]; // create new vector of Coupling matrices

    std::map<int, std::vector<int> >::const_iterator j = patchcouplm.find( coup_sid );
    if ( j==patchcouplm.end() )
      dserror( "missing side" );

    // get number of dofs for coupling side/element
    const size_t ndof_i = j->second.size();

    Cuiui_matrices.resize(1);
    Cuiui_matrices[0].Shape(ndof_i,ndof_i); //Cuiui
  }

  //-----------------------------------------------------------------------------------
  //         evaluate element length, stabilization factors and average weights
  //-----------------------------------------------------------------------------------

  // initialize the characteristic element length and it's inverse
  double h_k = 0.0;
  double inv_hk = 1.0;

  //-----------------------------------------------------------------------------------
  // compute characteristic element length for background element in case of background-sided coupling

  if (cond_manager->HasAveragingStrategy(INPAR::XFEM::Xfluid_Sided))
  {
    h_k = ComputeCharEleLength(ele, ele_xyze, cond_manager, vcSet, bcells, bintpoints);
    inv_hk = 1.0 / h_k;
  }

  //Get materials for both master and slave side
  // If non-constant density or viscosity wants to be calculated on the different sides. A review over how the scalar variables
  // are set at the surface should be made.
  GetMaterialParametersVolumeCell(mat_master,densaf_master_,viscaf_master_,gamma_m_);
  if(mat_slave != Teuchos::null)
  {
    GetMaterialParametersVolumeCell(mat_slave,densaf_slave_,viscaf_slave_,gamma_s_);
    //Security check:
    if(gamma_s_!=gamma_m_)
    {
      std::cout << "Surface tension for master side: "<< gamma_m_ << ", is not equal to surface tension on slave side:" << gamma_s_ << std::endl;
      dserror("Non-matching surface tension provided for Master and Slave side.");
    }
  }

  //-----------------------------------------------------------------------------------
  //      surface integral --- loop sides
  //-----------------------------------------------------------------------------------
  // map of side-element id and Gauss points
  for ( std::map<int, std::vector<DRT::UTILS::GaussIntegration> >::const_iterator i=bintpoints.begin();
      i!=bintpoints.end();
      ++i )
  {
    TEUCHOS_FUNC_TIME_MONITOR("FluidEleCalcXFEM::GaussIntegrationloop");

    //-----------------------------------------------------------------------------------

    // interface normal vector, pointing from background domain into the interface
    normal_.Clear();
    // gauss-point coordinates
    x_side_.Clear();

    // we need an interface to the boundary element (for projection)
    Teuchos::RCP<DRT::ELEMENTS::XFLUID::SlaveElementInterface<distype> > si;

    // location array of boundary element
    DRT::Element::LocationArray cutla( 1 );

    // pointer to boundary element
    DRT::Element * side = NULL;

    // coordinates of boundary element
    Epetra_SerialDenseMatrix side_xyze;

    //-----------------------------------------------------------------------------------
    // only used for couplings:

    // coupling object between background element and each coupling element (side for xfluid-sided coupling, element for other couplings)
    Teuchos::RCP<DRT::ELEMENTS::XFLUID::NitscheInterface<distype> > ci;

    // pointer to coupling element
    DRT::Element * coupl_ele = NULL;

    // coupling element coordinates
    Epetra_SerialDenseMatrix coupl_xyze;

    //-----------------------------------------------------------------------------------

    int coup_sid = i->first; // global coupling side id

    // get the coupling strategy for coupling of two fields
    const XFEM::EleCoupCond & coupcond = cond_manager->GetCouplingCondition(coup_sid, my::eid_);
    const INPAR::XFEM::EleCouplingCondType & cond_type = coupcond.first;

    const int coup_idx = cond_manager->GetCouplingIndex(coup_sid, my::eid_);
    Teuchos::RCP<XFEM::CouplingBase> coupling = cond_manager->GetCouplingByIdx(coup_idx);


    const std::vector<DRT::UTILS::GaussIntegration> & cutintpoints = i->second;

    std::map<int, std::vector<GEO::CUT::BoundaryCell*> >::const_iterator j = bcells.find( coup_sid );
    if ( j==bcells.end() )
      dserror( "missing boundary cell" );

    const std::vector<GEO::CUT::BoundaryCell*> & bcs = j->second;
    if ( bcs.size()!=cutintpoints.size() )
      dserror( "boundary cell integration rules mismatch" );

    //-----------------------------------------------------------------------------------
    // define average weights

    bool non_xfluid_coupling;
    double kappa_m;
    double kappa_s;

    cond_manager->GetAverageWeights(coup_sid, ele, kappa_m, kappa_s, non_xfluid_coupling);

    //---------------------------------------------------------------------------------
    // set flags used for coupling with given levelset/mesh coupling side
    bool is_ls_coupling_side   = cond_manager->IsLevelSetCoupling(coup_sid);
    bool is_mesh_coupling_side = cond_manager->IsMeshCoupling(coup_sid);

    Teuchos::RCP<DRT::Discretization> cutter_dis = cond_manager->GetCutterDis(coup_sid);

#ifdef DEBUG
    if( is_ls_coupling_side and  is_mesh_coupling_side) dserror("side cannot be a levelset-coupling side and a mesh coupling side at once: side %i", coup_sid);
    if(!is_ls_coupling_side and !is_mesh_coupling_side) dserror("side is neither a levelset-coupling side nor a mesh coupling side: side %i", coup_sid);
#endif

    //-----------------------------------------------------------------------------------
    Teuchos::RCP<XFEM::MeshCouplingFSI> mc_fsi = Teuchos::null;
    bool assemble_iforce = false;

    //---------------------------------------------------------------------------------
    // prepare the coupling objects
    if(is_mesh_coupling_side)
    {
      mc_fsi = Teuchos::rcp_dynamic_cast<XFEM::MeshCouplingFSI>(coupling);
      if(mc_fsi != Teuchos::null) assemble_iforce = true;

      // get the side element and its coordinates for projection of Gaussian points
      side = cond_manager->GetSide( coup_sid );
      GEO::InitialPositionArray(side_xyze,side);

      // create auxiliary coupling object for the boundary element, in order to perform projection
      si = DRT::ELEMENTS::XFLUID::SlaveElementInterface<distype>::CreateSlaveElementRepresentation( side, side_xyze );

      // set displacement of side
      side->LocationVector(*cutter_dis,cutla,false);
      si->AddSlaveEleDisp(*cutter_dis,cutla[0].lm_);

      if(cond_type == INPAR::XFEM::CouplingCond_SURF_WEAK_DIRICHLET or
         cond_type == INPAR::XFEM::CouplingCond_SURF_FSI_PART or
         cond_type == INPAR::XFEM::CouplingCond_SURF_CRACK_FSI_PART or
         cond_type == INPAR::XFEM::CouplingCond_SURF_NAVIER_SLIP)
      {
        si->SetInterfaceJumpStatenp(*cutter_dis, "ivelnp", cutla[0].lm_);
        if (my::fldparatimint_->IsNewOSTImplementation())
          si->SetInterfaceJumpStaten(*cutter_dis, "iveln", cutla[0].lm_);
      }
    }

    Teuchos::RCP<DRT::Discretization> coupl_dis_ = cond_manager->GetCouplingDis(coup_sid);

    if(!(is_ls_coupling_side and !cond_manager->IsCoupling( coup_sid, my::eid_ ))) // not level-set-WDBC case
    {
      coupl_ele = cond_manager->GetCouplingElement(coup_sid, ele);
      if (coupl_ele == NULL)
        dserror("Failed to obtain coupling element for global coup_sid %d", coup_sid);
      GEO::InitialPositionArray(coupl_xyze,coupl_ele);
    }

    if(!cond_manager->IsCoupling( coup_sid, my::eid_ ))
    {
      if(is_ls_coupling_side) //... for problems with cut interface defined by level-set field, currently only one-sided
      {
        ci = DRT::ELEMENTS::XFLUID::NitscheInterface<distype>::CreateNitscheCoupling_XFluidWDBC(
          elemat1_epetra,elevec1_epetra,*fldparaxfem_);
      }
      else if(is_mesh_coupling_side)
      {
        ci = DRT::ELEMENTS::XFLUID::NitscheInterface<distype>::CreateNitscheCoupling_XFluidWDBC(
            coupl_ele,coupl_xyze,elemat1_epetra,elevec1_epetra,*fldparaxfem_);
      }
    }
    else // coupling
    {
      std::map<int,std::vector<Epetra_SerialDenseMatrix> >::iterator c = side_coupling.find( coup_sid );

      std::vector<Epetra_SerialDenseMatrix> & side_matrices = c->second;

      // coupling matrices between background element and one! side
      Epetra_SerialDenseMatrix & C_uiu  = side_matrices[0];
      Epetra_SerialDenseMatrix & C_uui  = side_matrices[1];
      Epetra_SerialDenseMatrix & rhC_ui = side_matrices[2];

      // coupling matrices between one side and itself via the element Kss
      std::map<int,std::vector<Epetra_SerialDenseMatrix> >::iterator c2 = Cuiui_coupling.find( coup_sid );
      std::vector<Epetra_SerialDenseMatrix> & Cuiui_matrices = c2->second;
      Epetra_SerialDenseMatrix & eleCuiui = Cuiui_matrices[0];

      if (non_xfluid_coupling)
      {
        // create interface for the embedded element and the associated side
        ci = DRT::ELEMENTS::XFLUID::NitscheInterface<distype>::CreateNitscheCoupling_TwoSided(
            coupl_ele,coupl_xyze,elemat1_epetra,C_uiu,C_uui,eleCuiui,elevec1_epetra,rhC_ui,*fldparaxfem_);
      }
      else // ... for xfluid-sided coupling
      {
        ci = DRT::ELEMENTS::XFLUID::NitscheInterface<distype>::CreateNitscheCoupling_XFluidSided(
            coupl_ele,coupl_xyze,elemat1_epetra,C_uiu,C_uui,eleCuiui,elevec1_epetra,rhC_ui,*fldparaxfem_);
      }
    }

    if(cond_manager->IsCoupling( coup_sid, my::eid_ ))
    {
      std::map<int, std::vector<int> >::const_iterator k = patchcouplm.find( coup_sid );
      const std::vector<int> & coupl_lm = k->second;

      // set velocity (and pressure) of coupling/slave element at current time step
      ci->SetSlaveState(*coupl_dis_,coupl_lm);

      // set velocity (and pressure) of coupling element at old time step
      if (my::fldparatimint_->IsNewOSTImplementation())
        ci->SetSlaveStaten(*coupl_dis_,coupl_lm);
    }


    if(!(is_ls_coupling_side and !cond_manager->IsCoupling( coup_sid, my::eid_ ))) // not level-set-WDBC case
    {
      std::map<int, std::vector<int> >::const_iterator k = patchcouplm.find( coup_sid );
      const std::vector<int> & coupl_lm = k->second;

      // add displacement of coupling element at current time step
      ci->AddSlaveEleDisp(*coupl_dis_,coupl_lm);
    }



    if(cond_manager->IsCoupling( coup_sid, my::eid_ ) and non_xfluid_coupling)
    {
      //---------------------------------------------------------------------------------
      // compute characteristic element length for the case of embedded-sided coupling

      // char. length defined by local eigenvalue problem
      if (fldparaxfem_->ViscStabTracEstimate() == INPAR::XFEM::ViscStab_TraceEstimate_eigenvalue)
      {
        inv_hk = cond_manager->Get_TraceEstimate_MaxEigenvalue(coup_sid);
#if(0)
        std::cout.precision(15);
        std::cout << "C_T/hk (formula): "
            << NIT_getTraceEstimateConstant(ele_distype)/ComputeCharEleLength(coupl_ele, coupl_xyze, cond_manager, vcSet, bcells, bintpoints, emb, side);
        << " max_eigenvalue ~ C_T/hk: "
            <<  my::fldpara_->Get_TraceEstimate_MaxEigenvalue(coup_sid)
            << " max_eigenvalue*h_k = C_T: "<< my::fldpara_->Get_TraceEstimate_MaxEigenvalue(coup_sid)*h_k << " vs: C_T (formula) " << NIT_getTraceEstimateConstant(ele_distype) << std::endl;
#endif
        h_k = 1.0 / inv_hk;
      }
      else // ... char. length defined otherwise
      {
        // compute characteristic element length based on the embedded element
        h_k = ComputeCharEleLength(coupl_ele, coupl_xyze, cond_manager, vcSet, bcells, bintpoints, ci, side);
        inv_hk = 1.0 / h_k;
      }
    }


    //---------------------------------------------------------------------------------
    // compute viscous part of Nitsche's penalty term scaling for Nitsche's method
    // based on the inverse characteristic element length
    //---------------------------------------------------------------------------------

    double NIT_visc_stab_fac = 0.0;
    cond_manager->Get_ViscPenalty_Stabfac(coup_sid, ele,kappa_m,kappa_s, inv_hk,fldparaxfem_,NIT_visc_stab_fac);

    // define interface force vector w.r.t side (for XFSI)
    Epetra_SerialDenseVector iforce;
    iforce.Size(cutla[0].lm_.size());

    //---------------------------------------------------------------------------------
    // loop boundary cells w.r.t current cut side
    //---------------------------------------------------------------------------------
    for ( std::vector<DRT::UTILS::GaussIntegration>::const_iterator i=cutintpoints.begin();
          i!=cutintpoints.end();
          ++i )
    {
      const DRT::UTILS::GaussIntegration & gi = *i;
      GEO::CUT::BoundaryCell * bc = bcs[i - cutintpoints.begin()]; // get the corresponding boundary cell

      //-------------------------------------------------------------------------------
      // loop gausspoints w.r.t current boundary cell
      //-------------------------------------------------------------------------------
      //-------------------------------------------------------------------------------
      for ( DRT::UTILS::GaussIntegration::iterator iquad=gi.begin(); iquad!=gi.end(); ++iquad )
      {
        double drs = 0.0; // transformation factor between reference cell and linearized boundary cell

        const LINALG::Matrix<2,1> eta( iquad.Point() ); // xi-coordinates with respect to side



        // compute transformation factor, normal vector and global Gauss point coordinates
        if (bc->Shape() != DRT::Element::dis_none) // Tessellation approach
        {
          TEUCHOS_FUNC_TIME_MONITOR( "FluidEleCalcXFEM::ComputeSurfaceTransformation" );

          ComputeSurfaceTransformation(drs, x_gp_lin_, normal_, bc, eta);
        }
        else // MomentFitting approach
        {
          drs = 1.0;
          normal_ = bc->GetNormalVector();
          const double* gpcord = iquad.Point();
          for (int idim=0;idim<3;++idim)
          {
            x_gp_lin_(idim,0) = gpcord[idim];
          }
        }

        {
          TEUCHOS_FUNC_TIME_MONITOR("FluidEleCalcXFEM::GEO::CUT::Position");

          // find element local position of gauss point
          GEO::CUT::Position<distype> pos( my::xyze_, x_gp_lin_ );
          pos.Compute();
          rst_ = pos.LocalCoordinates();
        }

        if (is_mesh_coupling_side)
        {
          // project gaussian point from linearized interface to warped side (get/set local side coordinates in SideImpl)
          LINALG::Matrix<3,1> xi_side;

          // project on boundary element
          si->ProjectOnSide(x_gp_lin_, x_side_, xi_side);

          if (non_xfluid_coupling)
            ci->Evaluate( x_side_ ); // evaluate embedded element's shape functions at gauss-point coordinates
          else
            ci->Evaluate( xi_side ); // evaluate side's shape functions at gauss-point coordinates
        }
        else if (is_ls_coupling_side)
        {
          if(cond_manager->IsCoupling( coup_sid, my::eid_ ))
            ci->Evaluate( x_gp_lin_ ); // evaluate embedded element's shape functions at gauss-point coordinates
        }

        // integration factors
        const double surf_fac = drs*iquad.Weight();
        const double timefacfac = surf_fac * my::fldparatimint_->TimeFac();

        // evaluate background element shape functions
        EvalFuncAndDeriv( rst_ );

        // get velocity at integration point
        // (values at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
        my::velint_.Multiply(evelaf_,my::funct_);

        // get velocity derivatives at integration point
        // (values at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
        my::vderxy_.MultiplyNT(evelaf_,my::derxy_);

        // get pressure at integration point
        // (value at n+1)
        double press = my::funct_.Dot(epreaf_);

        //----------------------------------------------
        // get convective velocity at integration point
        my::SetConvectiveVelint(ele->IsAle());

        //-----------------------------------------------------------------------------
        // compute stabilization factors

        double NIT_full_stab_fac = 0.0;

        //Extract slave velocity at Gausspoint
        ci->GetInterfaceVelnp(velint_s_);

        NIT_Compute_FullPenalty_Stabfac(
            NIT_full_stab_fac,             ///< to be filled: full Nitsche's penalty term scaling (viscous+convective part)
            normal_,
            h_k,
            kappa_m,
            kappa_s,
            my::convvelint_,
            velint_s_,
            NIT_visc_stab_fac             ///< Nitsche's viscous scaling part of penalty term
            );


        //-----------------------------------------------------------------------------
        // define the prescribed interface jump vectors for velocity and traction

        ivelint_jump_.Clear();
        itraction_jump_.Clear();
        proj_tangential_.Clear();
        LB_proj_matrix_.Clear();

        GetInterfaceJumpVectors(
            coupcond,
            coupling,
            ivelint_jump_,
            itraction_jump_,
            proj_tangential_,
            LB_proj_matrix_,
            x_gp_lin_,
            normal_,
            si,
            rst_,
            kappa_m,
            viscaf_master_,
            viscaf_slave_
        );

        if(cond_type == INPAR::XFEM::CouplingCond_LEVELSET_NEUMANN or
           cond_type == INPAR::XFEM::CouplingCond_SURF_NEUMANN)
        {
          //-----------------------------------------------------------------------------
          // evaluate the Neumann boundary condition term
          EvaluateNeumann(
              timefacfac,             ///< theta*dt
              my::funct_,             ///< coupling master shape functions
              itraction_jump_,        ///< prescribed interface traction, jump height for coupled problems
              elevec1_epetra          ///< element rhs vector
          );

          if (my::fldparatimint_->IsNewOSTImplementation())
          {
            dserror("how to deal with Neumann boundary condition and new OSTImplementation");
          }
        }
        else
        {
          TEUCHOS_FUNC_TIME_MONITOR( "FluidEleCalcXFEM::NIT_evaluateCoupling" );

          //Get Configuration Map
          std::map<INPAR::XFEM::CoupTerm, std::pair<bool,double> >& configmap =
              coupling->GetConfigurationmap(kappa_m,viscaf_master_,viscaf_slave_,NIT_visc_stab_fac, NIT_full_stab_fac,x_gp_lin_,coupcond.second);

          //-----------------------------------------------------------------------------
          // evaluate the coupling terms for coupling with current side
          // (or embedded element through current side)
          // time step n+1

          const bool isImplPressureNewOst(my::fldparatimint_->IsFullImplPressureAndCont() && my::fldparatimint_->IsNewOSTImplementation());

          const double pres_timefacfac(isImplPressureNewOst ? my::fldparatimint_->Dt() * surf_fac : timefacfac);

          ci->NIT_evaluateCoupling(
            normal_,                     // normal vector
            timefacfac,                  // theta*dt*fac
            pres_timefacfac,             // impl. pressure with new OST: dt * fac, else theta*dt*fac
            viscaf_master_,              // dynvisc viscosity in background fluid
            viscaf_slave_,               // dynvisc viscosity in embedded fluid
            my::densaf_,                 // fluid density
            my::funct_,                  // bg shape functions
            my::derxy_,                  // bg shape function gradient
            my::vderxy_,                 // bg grad u^n+1
            press,                       // bg p^n+1
            my::velint_,                 // bg u^n+1
            ivelint_jump_,               // prescribed interface velocity, Dirichlet values or jump height for coupled problems
            itraction_jump_,             // traction jump at interface (i.e. [| -pI + \mu*[\nabla u + (\nabla u)^T]  |] \cdot n)
            proj_tangential_,            // tangential projection matrix
            LB_proj_matrix_,             // prescribed projection matrix for laplace-beltrami problems
            configmap                    // Configuration Map
          );

          if (my::fldparatimint_->IsNewOSTImplementation())
          {
            // get velocity at integration point
            // (values at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
            my::velintn_.Multiply(eveln_,my::funct_);

            // get velocity derivatives at integration point
            // (values at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
            my::vderxyn_.MultiplyNT(eveln_,my::derxy_);

            ivelintn_jump_.Clear();
            itractionn_jump_.Clear();

            // Safety check
            if(cond_type == INPAR::XFEM::CouplingCond_LEVELSET_NAVIER_SLIP or
               cond_type == INPAR::XFEM::CouplingCond_SURF_NAVIER_SLIP       )
            {

              if (my::fldparatimint_->IsNewOSTImplementation())
              {
                dserror("How to deal with NavierSlip boundary condition and new OSTImplementation?");
              }

            }

            if(cond_type != INPAR::XFEM::CouplingCond_LEVELSET_TWOPHASE and
                cond_type != INPAR::XFEM::CouplingCond_LEVELSET_COMBUSTION)
            {
              GetInterfaceJumpVectorsOldState(
                  coupcond,
                  coupling,
                  ivelintn_jump_,
                  itractionn_jump_,
                  x_gp_lin_,
                  normal_,
                  si,
                  my::funct_.Dot(epren_),       // bg p^n
                  rst_
              );
            }
            else
            {
              GetInterfaceJumpVectorsOldState(
                  coupcond,
                  coupling,
                  ivelintn_jump_,
                  itractionn_jump_,
                  x_gp_lin_,
                  normal_,
                  ci,
                  my::funct_.Dot(epren_),       // bg p^n
                  rst_
              );
            }

            //-----------------------------------------------------------------------------
            // evaluate the coupling terms for coupling with current side
            // (or embedded element through current side)
            // time step n

            // REMARK: evaluation of all Nitsche terms for the previous solution
            // only makes sense for stationary interfaces!

            double NIT_full_stab_fac_n = 0.0;
            if (fldparaxfem_->InterfaceTermsPreviousState() == INPAR::XFEM::PreviousState_full)
            {
              velintn_s_.Clear();
              ci->GetInterfaceVeln(velintn_s_);

              NIT_Compute_FullPenalty_Stabfac(
                NIT_full_stab_fac_n,  ///< to be filled: full Nitsche's penalty term scaling (viscous+convective part)
                normal_,
                h_k,
                kappa_m, //weights (only existing for Nitsche currently!!)
                kappa_s, //weights (only existing for Nitsche currently!!)
                my::convvelintn_,
                velintn_s_,
                NIT_visc_stab_fac   ///< Nitsche's viscous scaling part of penalty term
              );
            }

            //Get Configuration Map
            std::map<INPAR::XFEM::CoupTerm, std::pair<bool,double> >& configmap_n =
                coupling->GetConfigurationmap(kappa_m,viscaf_master_,viscaf_slave_,NIT_visc_stab_fac, NIT_full_stab_fac,x_gp_lin_,coupcond.second);

            const double timefacfacn = surf_fac * (my::fldparatimint_->Dt()-my::fldparatimint_->TimeFac());
            ci->NIT_evaluateCouplingOldState(
              normal_,
              timefacfacn,
              isImplPressureNewOst,
              viscaf_master_,              // dynvisc viscosity in background fluid
              viscaf_slave_,               // dynvisc viscosity in embedded fluid
              my::densn_,                  // fluid density
              my::funct_,                  // bg shape functions
              my::derxy_,                  // bg shape function gradient
              my::vderxyn_,                // bg grad u^n
              my::funct_.Dot(epren_),       // bg p^n
              my::velintn_,                 // bg u^n
              ivelintn_jump_,               // velocity jump at interface (i.e. [| u |])
              proj_tangential_,            // tangential projection matrix
              itractionn_jump_,             // traction jump at interface (i.e. [| -pI + \mu*[\nabla u + (\nabla u)^T]  |] \cdot n)
              configmap_n
            );
          }
        }

        if (!assemble_iforce)
          continue;

        //-----------------------------------------------------------------------------
        // calculate interface forces for XFSI

        //-------------------------------
        // traction vector w.r.t fluid domain, resulting stresses acting on the fluid surface
        // t= (-p*I + 2mu*eps(u))*n^f
        LINALG::Matrix<my::nsd_,1> traction;

        BuildTractionVector( traction, press, normal_ );
        ci->ComputeInterfaceForce(iforce, traction, surf_fac );

      } // end loop gauss points of boundary cell
    } // end loop boundary cells of side

    if(assemble_iforce)
      AssembleInterfaceForce(mc_fsi->IForcecol(), *cutter_dis, cutla[0].lm_, iforce);
  } // end loop cut sides


  //-----------------------------------------------------------------------------------
  // build Cuiui coupling matrix (includes patch of Cuiui matrices for all sides)
  //-----------------------------------------------------------------------------------

  NIT_BuildPatchCuiui(Cuiui, Cuiui_coupling);

  return;
}



/*--------------------------------------------------------------------------------
 * get the interface jump vectors for velocity and traction at the Gaussian point
 *--------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void FluidEleCalcXFEM<distype>::GetInterfaceJumpVectors(
    const XFEM::EleCoupCond & coupcond,                                      ///< coupling condition for given interface side
    Teuchos::RCP<XFEM::CouplingBase> coupling,                               ///< coupling object
    LINALG::Matrix<my::nsd_,1>& ivelint_jump,                                ///< prescribed interface jump vector for velocity
    LINALG::Matrix<my::nsd_,1>& itraction_jump,                              ///< prescribed interface jump vector for traction
    LINALG::Matrix<my::nsd_,my::nsd_>&  proj_tangential,                     ///< tangential projection matrix
    LINALG::Matrix<my::nsd_,my::nsd_>&  LB_proj_matrix,                      ///< prescribed projection matrix for laplace-beltrami problems
    const LINALG::Matrix<my::nsd_,1>& x,                                     ///< global coordinates of Gaussian point
    const LINALG::Matrix<my::nsd_,1>& normal,                                ///< normal vector at Gaussian point
    Teuchos::RCP<DRT::ELEMENTS::XFLUID::SlaveElementInterface<distype> > si, ///< side implementation for cutter element
    LINALG::Matrix<3,1>& rst,                                                ///< local coordinates of GP for bg element
    double& kappa_m,                                                         ///< fluid sided weighting
    double& visc_m,                                                          ///< fluid sided weighting
    double& visc_s                                                           ///< slave sided dynamic viscosity
)
{
  TEUCHOS_FUNC_TIME_MONITOR("FluidEleCalcXFEM::GetInterfaceJumpVectors");

  // [| v |] := vm - vs

  const INPAR::XFEM::EleCouplingCondType & cond_type = coupcond.first; ///< condition type for given interface side
  const DRT::Condition* cond = coupcond.second;                        ///< condition to be evaluated

  switch (cond_type)
  {
  case INPAR::XFEM::CouplingCond_SURF_WEAK_DIRICHLET:
  {
    const std::string* evaltype = cond->Get<std::string>("evaltype");

    if(*evaltype == "funct_gausspoint")
    {
      // evaluate function at Gaussian point at current time
      coupling->EvaluateCouplingConditions(ivelint_jump,itraction_jump,x,cond); //itraction_jump.Clear() called here...
    }
    else
    {
      // evaluate function at nodes at current time
      si->GetInterfaceJumpVelnp(ivelint_jump);
    }

    break;
  }
  case INPAR::XFEM::CouplingCond_SURF_NEUMANN:
  case INPAR::XFEM::CouplingCond_LEVELSET_WEAK_DIRICHLET:
  case INPAR::XFEM::CouplingCond_LEVELSET_NEUMANN:
  {
    // evaluate condition function at Gaussian point
    coupling->EvaluateCouplingConditions(ivelint_jump,itraction_jump,x,cond);
    break;
  }
  case INPAR::XFEM::CouplingCond_SURF_FSI_PART:
  case INPAR::XFEM::CouplingCond_SURF_CRACK_FSI_PART:
  {
    // evaluate function at nodes at current time
    si->GetInterfaceJumpVelnp(ivelint_jump);
    break;
  }
  case INPAR::XFEM::CouplingCond_SURF_FLUIDFLUID:
  case INPAR::XFEM::CouplingCond_SURF_FSI_MONO:
  {
    // nothing to evaluate as continuity coupling conditions have to be evaluated
    break;
  }
  case INPAR::XFEM::CouplingCond_SURF_FPI_MONO:
  {
    Teuchos::rcp_dynamic_cast<XFEM::MeshCouplingFPI>(coupling)->EvaluateCouplingConditions<distype>(proj_tangential,normal);
    break;
  }
  case INPAR::XFEM::CouplingCond_SURF_NAVIER_SLIP:
  {


    bool eval_dirich_at_gp = (*(cond->Get<std::string>("evaltype")) == "funct_gausspoint");

    // The velocity is evaluated twice in this framework...
    Teuchos::rcp_dynamic_cast<XFEM::MeshCouplingNavierSlip>(coupling)->EvaluateCouplingConditions(ivelint_jump,itraction_jump,proj_tangential,x,normal,cond,eval_dirich_at_gp,kappa_m,visc_m,visc_s);

    if(!eval_dirich_at_gp)
    {

      si->GetInterfaceJumpVelnp(ivelint_jump);

    }

    break;
  }
  case INPAR::XFEM::CouplingCond_LEVELSET_NAVIER_SLIP:
  {

    Teuchos::rcp_dynamic_cast<XFEM::LevelSetCouplingNavierSlip>(coupling)->EvaluateCouplingConditions<distype>(ivelint_jump,itraction_jump,x,cond,proj_tangential,my::eid_,my::funct_,my::derxy_,normal,kappa_m,visc_m,visc_s);

    break;
  }
  case INPAR::XFEM::CouplingCond_LEVELSET_TWOPHASE:
  {
    // where [*] = (*)^m - (*)^s = (*)^+ - (*)^-
    // n = n^m = n^+
    // [sigma*n] = gamma * curv * n   with curv = div(grad(phi)/||grad(phi)||)

    double surf_coeff    = gamma_m_;

    if(gamma_m_ != 0.0)
    {
      Teuchos::rcp_dynamic_cast<XFEM::LevelSetCouplingTwoPhase>(coupling)->EvaluateTractionDiscontinuity<distype>(itraction_jump,LB_proj_matrix,my::eid_,my::funct_,my::derxy_,normal,surf_coeff);
    }

    break;
  }
  case INPAR::XFEM::CouplingCond_LEVELSET_COMBUSTION:
  {
    // TODO: evaluate the ivelint_jump and the itraction_jump
    break;
  }
  // Neumann boundary conditions for Mesh and Levelset
  default: dserror("invalid type of condition %i, which prescribed interface vectors have to be set?", cond_type); break;
  }

  // Create a projection matrix.
  //  If it is a Navier-Slip coupling, the matrix is provided from the Evaluation.
  //   Furthermore, if it is a Laplace-Beltrami way of calculating the surface tension,
  //   do not fill the matrix as it contains the "projection matrix" for LB implementation.
  if(cond_type != INPAR::XFEM::CouplingCond_LEVELSET_NAVIER_SLIP
      and cond_type != INPAR::XFEM::CouplingCond_SURF_NAVIER_SLIP)
  {
    //Create normal projection matrix.
    LINALG::Matrix<my::nsd_,my::nsd_> eye(true);
    for(int i =0; i<my::nsd_; ++i)
      eye(i,i)=1;
    for(int i =0; i<my::nsd_; ++i)
    {
      for(int j =0; j<my::nsd_; ++j)
      {
        proj_tangential(i,j)          = eye(i,j) - normal(i,0) * normal(j,0);
      }
    }
  }

  return;
}

/*--------------------------------------------------------------------------------
 * get the interface jump vectors for velocity and traction at the Gaussian point
 * for previous time step
 *--------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void FluidEleCalcXFEM<distype>::GetInterfaceJumpVectorsOldState(
    const XFEM::EleCoupCond & coupcond,                                      ///< coupling condition for given interface side
    Teuchos::RCP<XFEM::CouplingBase> coupling,                               ///< coupling object
    LINALG::Matrix<my::nsd_,1>& ivelintn_jump,                                ///< prescribed interface jump vector for velocity
    LINALG::Matrix<my::nsd_,1>& itractionn_jump,                              ///< prescribed interface jump vector for traction
    const LINALG::Matrix<my::nsd_,1>& x,                                     ///< global coordinates of Gaussian point
    const LINALG::Matrix<my::nsd_,1>& normal,                                ///< normal vector at Gaussian point
    Teuchos::RCP<DRT::ELEMENTS::XFLUID::SlaveElementInterface<distype> > si, ///< side implementation for cutter element
    const double &                           presn_m,                        ///< coupling master pressure
    LINALG::Matrix<3,1>& rst                                                 ///< local coordinates of GP for bg element
)
{

  // [| v |] := vm - vs

  const INPAR::XFEM::EleCouplingCondType & cond_type = coupcond.first; ///< condition type for given interface side
  const DRT::Condition* cond = coupcond.second;                        ///< condition to be evaluated

  switch (cond_type)
  {
  case INPAR::XFEM::CouplingCond_SURF_WEAK_DIRICHLET:
  {
    const std::string* evaltype = cond->Get<std::string>("evaltype");

    if(*evaltype == "funct_gausspoint")
    {
      // evaluate function at Gaussian point at current time
      coupling->EvaluateCouplingConditionsOldState(ivelintn_jump,itractionn_jump,x,cond);
    }
    else
    {
      // evaluate function at nodes for previous time
      si->GetInterfaceJumpVeln(ivelintn_jump);
    }

    break;
  }
  case INPAR::XFEM::CouplingCond_SURF_NEUMANN:
  case INPAR::XFEM::CouplingCond_LEVELSET_WEAK_DIRICHLET:
  case INPAR::XFEM::CouplingCond_LEVELSET_NEUMANN:
  {
    // evaluate condition function at Gaussian point
    coupling->EvaluateCouplingConditionsOldState(ivelintn_jump,itractionn_jump,x,cond);
    break;
  }
  case INPAR::XFEM::CouplingCond_SURF_NAVIER_SLIP:
  case INPAR::XFEM::CouplingCond_LEVELSET_NAVIER_SLIP:
  {
    dserror("Navier Slip Condition not implemented for NEWOst yet!");
    //here you would need the dyn_visc for summing up vel_jump and traction_jump...
    break;
  }
  case INPAR::XFEM::CouplingCond_SURF_FSI_PART:
  case INPAR::XFEM::CouplingCond_SURF_CRACK_FSI_PART:
  {
    // evaluate function at nodes at current time
    si->GetInterfaceJumpVeln(ivelintn_jump);
    break;
  }
  case INPAR::XFEM::CouplingCond_SURF_FLUIDFLUID:
  case INPAR::XFEM::CouplingCond_SURF_FSI_MONO:
  {
    // nothing to evaluate as continuity coupling conditions have to be evaluated
    break;
  }
  case INPAR::XFEM::CouplingCond_SURF_FPI_MONO:
  {
    dserror("Fluid Poro Structure Interaction not implemented for NEWOst yet!");
    break;
  }
  case INPAR::XFEM::CouplingCond_LEVELSET_TWOPHASE:
  {

    // Spatial velocity gradient for slave side
    LINALG::Matrix<my::nsd_,my::nsd_> vderxyn_s(true);
    si->GetInterfaceVelGradn(vderxyn_s);


    //Calculate the old jump using the reconstructed values of p and u for the new interface position.
    // i.e.

    // itractionn_jump = [|  -pI + \mu*[\nabla u + (\nabla u)^T]  |] \cdot n

    // Pressure part
    double presn_s = 0.0;
    si->GetInterfacePresn(presn_s);

    itractionn_jump.Update(-(presn_m - presn_s),normal,0.0);

    // Shear tensor part
    //===================
    LINALG::Matrix<my::nsd_,my::nsd_> tmp_matrix(true);
    tmp_matrix.Update(viscaf_master_,my::vderxyn_,-viscaf_slave_,vderxyn_s);

    //Initialize dummy variable
    LINALG::Matrix<my::nsd_,1> tmp_vector(true);

    //Normal
    tmp_vector.Multiply(tmp_matrix,normal);
    itractionn_jump.Update(1.0,tmp_vector,1.0);

    //Transposed
    tmp_vector.MultiplyTN(tmp_matrix,normal);
    itractionn_jump.Update(1.0,tmp_vector,1.0);
    //===================

    break;
  }
  case INPAR::XFEM::CouplingCond_LEVELSET_COMBUSTION:
  {
    // TODO: evaluate the ivelint_jump and the itraction_jump
    break;
  }
  // Neumann boundary conditions for Mesh and Levelset
  default: dserror("invalid type of condition %i, which prescribed interface vectors have to be set?", cond_type); break;
  }

  return;
}




/*--------------------------------------------------------------------------------
 * build the patch coupling matrix Cuiui containing Cuiui for all cutting sides
 *--------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void FluidEleCalcXFEM<distype>::NIT_BuildPatchCuiui(
    Epetra_SerialDenseMatrix &                              Cuiui,            ///< ui-ui patch coupling matrix containing Cuiui for all cutting sides
    std::map<int, std::vector<Epetra_SerialDenseMatrix> >&  Cuiui_coupling    ///< Cuiui matrices for all cutting sides
)
{
  //TEUCHOS_FUNC_TIME_MONITOR("FluidEleCalcXFEM::NIT_BuildPatchCuiui");

  // build patch-Cuiui matrix
  int ipatchsizesbefore = 0;
  for (std::map<int, std::vector<Epetra_SerialDenseMatrix> >::const_iterator m=Cuiui_coupling.begin();
      m!=Cuiui_coupling.end(); ++m)
  {

    int coup_sid = m->first;
    std::vector<Epetra_SerialDenseMatrix> & Cuiui_mats = Cuiui_coupling[coup_sid];

    // Cuiui matrices in Cuiui_mats[0]

    // assemble Cuiui
    for ( int ic=0; ic<Cuiui_mats[0].N(); ++ic) // Cuiui includes only ui,ui coupling, not (ui,p) ...
    {
      for ( int ir=0; ir<Cuiui_mats[0].M(); ++ir )
      {
        Cuiui(ir+ipatchsizesbefore,ic+ipatchsizesbefore) = Cuiui_mats[0](ir,ic);
      }
    }

    ipatchsizesbefore += Cuiui_mats[0].N();

  }

  return;
}

/*--------------------------------------------------------------------------------
 *    compute stabilization factor for the Nitsche's penalty term
 *--------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void FluidEleCalcXFEM<distype>::NIT_Compute_FullPenalty_Stabfac(
    double &                           NIT_full_stab_fac,             ///< to be filled: full Nitsche's penalty term scaling (viscous+convective part)
    const LINALG::Matrix<my::nsd_,1>&  normal,                        ///< interface-normal vector
    const double                       h_k,                           ///< characteristic element length
    const double                       kappa_m,                        ///< Weight parameter (parameter +/master side)
    const double                       kappa_s,                        ///< Weight parameter (parameter -/slave  side)
    const LINALG::Matrix<my::nsd_,1>&  velint_m,                      ///< Master side velocity at gauss-point
    const LINALG::Matrix<my::nsd_,1>&  velint_s,                      ///< Slave side velocity at gauss-point
    const double                       NIT_visc_stab_fac,             ///< Nitsche's viscous scaling part of penalty term
    bool                               error_calc                     ///< when called in error calculation, don't add the inflow terms
)
{
  //TEUCHOS_FUNC_TIME_MONITOR("FluidEleCalcXFEM::NIT_Compute_FullPenalty_Stabfac");

  //------------------------------------------------------------------------------
  // compute the full Nitsche parameter

  /*
   * Depending on the flow regime, the factor alpha of Nitsches penalty term
   * (\alpha * [v],[u]) can take various forms.
   * Based on INPAR::XFEM::MassConservationCombination, we choose:
   *
   *                       (1)           (2)          (3)
   *                    /  \mu    \rho             h * \rho         \
   *  NIT :=  \gamma * |    --  +  -- * |u|_inf  + ----------------- |
   *                    \   h      6               12 * \theta * dt /
   *
   *       OR:
   *                    /  \mu    \rho             h * \rho         \
   *  NIT :=  max      |    --  ;  -- * |u|_inf  ; ----------------- | *\gamma
   *                    \   h      6               12 * \theta * dt /
   *
   *          (1) NIT_visc_stab_fac = \gamma * \mu/h
   *
   *          (2) convective contribution
   *
   *          (3) transient contribution
   *
   * see Schott and Rasthofer, 'A face-oriented stabilized Nitsche-type extended variational
   * multiscale method for incompressible two-phase flow', Int. J. Numer. Meth. Engng, 2014
   *
   *
   * If INPAR::XFEM::MassConservationScaling_only_visc is set, we choose only (1),
   * no matter if the combination option max or sum is active!
   *
   */


  // (1)
  NIT_full_stab_fac = NIT_visc_stab_fac;

  if (fldparaxfem_->MassConservationScaling() == INPAR::XFEM::MassConservationScaling_full)
  {
    //TODO: Raffaela: which velocity has to be evaluated for these terms in ALE? the convective velocity or the velint?
    double velnorminf_m = velint_m.NormInf(); // relative convective velocity
    double velnorminf_s = velint_s.NormInf();

    // take the maximum of viscous & convective contribution or the sum?
    if (fldparaxfem_->MassConservationCombination() == INPAR::XFEM::MassConservationCombination_max)
    {
      NIT_full_stab_fac = std::max(NIT_full_stab_fac,fldparaxfem_->NITStabScaling() * (kappa_m*densaf_master_*fabs(velnorminf_m)+kappa_s*densaf_slave_*fabs(velnorminf_s)) / 6.0);
      if (! my::fldparatimint_->IsStationary())
        NIT_full_stab_fac= std::max( fldparaxfem_->NITStabScaling() * h_k * ( kappa_m*densaf_master_ + kappa_s*densaf_slave_ ) / (12.0 * my::fldparatimint_->TimeFac()),NIT_full_stab_fac);
    }
    else // the sum
    {
      // (2)
      NIT_full_stab_fac += fldparaxfem_->NITStabScaling() * (kappa_m*densaf_master_*fabs(velnorminf_m)+kappa_s*densaf_slave_*fabs(velnorminf_s)) / 6.0; //THIS ONE NEEDS CHANGING!

      // (3)
      if (! my::fldparatimint_->IsStationary())
        NIT_full_stab_fac += fldparaxfem_->NITStabScaling() * h_k * ( kappa_m*densaf_master_ + kappa_s*densaf_slave_ ) / (12.0 * my::fldparatimint_->TimeFac());
    }
  }
  else if (fldparaxfem_->MassConservationScaling() != INPAR::XFEM::MassConservationScaling_only_visc)
    dserror("Unknown scaling choice in calculation of Nitsche's penalty parameter");

  if (my::fldpara_->IsConservative() and (fldparaxfem_->XffConvStabScaling() != INPAR::XFEM::XFF_ConvStabScaling_none
                                       or fldparaxfem_->ConvStabScaling()    != INPAR::XFEM::ConvStabScaling_none ) )
  {
    dserror("convective stabilization is not available for conservative form of Navier-Stokes, but possible to implement!");
  }

  //----------------------------------------------------------------------------------------------
  // add inflow terms to ensure coercivity at inflow boundaries in the convective limit

  if ((fldparaxfem_->ConvStabScaling() == INPAR::XFEM::ConvStabScaling_none &&
       fldparaxfem_->XffConvStabScaling() == INPAR::XFEM::XFF_ConvStabScaling_none) ||
       fldparaxfem_->XffConvStabScaling() == INPAR::XFEM::XFF_ConvStabScaling_only_averaged || error_calc)
    return;

  const double veln_normal = velint_m.Dot(normal);

  double NIT_inflow_stab = 0.0;

  if (fldparaxfem_->XffConvStabScaling() == INPAR::XFEM::XFF_ConvStabScaling_upwinding)
  {
    NIT_inflow_stab = fabs(veln_normal)*0.5;
  }
  else
  {
    INPAR::XFEM::ConvStabScaling conv_stab_scaling = fldparaxfem_->ConvStabScaling();
    if (conv_stab_scaling == INPAR::XFEM::ConvStabScaling_abs_inflow)
    {
      //      | u*n |
      NIT_inflow_stab = fabs(veln_normal);
    }
    else if (conv_stab_scaling == INPAR::XFEM::ConvStabScaling_inflow)
    {
      //      ( -u*n ) if (u*n)<0 (inflow), conv_stabfac >= 0
      NIT_inflow_stab = std::max(0.0,-veln_normal);
    }
    else
      dserror("No valid INPAR::XFEM::ConvStabScaling for xfluid/xfsi problems");
  }

  NIT_inflow_stab *= densaf_master_; //my::densaf_;

  // Todo (kruse): it is planned to add the inflow contributions independent from the max. option!
  // This version is only kept to shift the adaption of test results to a single commit.
  if (fldparaxfem_->MassConservationCombination() == INPAR::XFEM::MassConservationCombination_max)
  {
    NIT_full_stab_fac = std::max(NIT_full_stab_fac,NIT_inflow_stab);
  }
  else if (fldparaxfem_->MassConservationCombination() == INPAR::XFEM::MassConservationCombination_sum)
  {
    NIT_full_stab_fac += NIT_inflow_stab;
  }
  else
    dserror("Unknown combination type in calculation of Nitsche's penalty parameter.");

  return;
}

/*--------------------------------------------------------------------------------
 * prepare coupling matrices, that include contributions from convective stabilization
 * and contributions from previous time steps (rhs)
 *--------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void FluidEleCalcXFEM<distype>::HybridLM_CreateSpecialContributionMatrices(
    const Teuchos::RCP<XFEM::ConditionManager> &            cond_manager,            ///< XFEM condition manager
    std::set<int> &                                         begids,                  ///< ids of intersecting boundary elements
    std::map<int, std::vector<Epetra_SerialDenseMatrix> > & side_coupling_extra       ///< contributions to coupling matrices from convective stabilizations
)
{
  if (fldparaxfem_->GetCouplingMethod() != INPAR::XFEM::Hybrid_LM_Cauchy_stress &&
      fldparaxfem_->GetCouplingMethod() != INPAR::XFEM::Hybrid_LM_viscous_stress)
    dserror("Do not call this method with a non-Lagrange multiplier based approach!");

  for (std::set<int>::const_iterator bgid=begids.begin(); bgid!=begids.end(); ++bgid)
  {
    const int coup_sid = *bgid;

    if(!cond_manager->IsCoupling( coup_sid, my::eid_ )) continue; // no coupling with current side

    if(cond_manager->IsLevelSetCoupling(coup_sid)) dserror("HybridLM_CreateSpecialContributionMatrices for level-set coupling not supported yet");

    Teuchos::RCP<DRT::Discretization> cutter_dis = Teuchos::null;
    if(cond_manager->IsMeshCoupling(coup_sid)) cutter_dis = cond_manager->GetCutterDis(coup_sid);

    DRT::Element * side = cond_manager->GetSide(coup_sid); // for each boundary element there is one corresponding side

    std::vector<int> patchlm;
    std::vector<int> patchlmowner;
    std::vector<int> patchlmstride;
    side->LocationVector(*cutter_dis, patchlm, patchlmowner, patchlmstride);

    // get coupling matrices for the current side (boundary element)
    std::vector<Epetra_SerialDenseMatrix> & side_matrices_extra = side_coupling_extra[coup_sid];

    side_matrices_extra.resize(4);
    side_matrices_extra[0].Shape(patchlm.size(), my::nen_*my::numdofpernode_); //Cuiu
    side_matrices_extra[1].Shape(my::nen_*my::numdofpernode_, patchlm.size()), //Cuui
    side_matrices_extra[2].Shape(patchlm.size(),1);                            //rhs_Cui
    side_matrices_extra[3].Shape(patchlm.size(), patchlm.size());              //Cuiui
  }
}


/*--------------------------------------------------------------------------------
 * compute transformation factor for surface integration, normal, local and global gp coordinates
 *--------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void FluidEleCalcXFEM<distype>::ComputeSurfaceTransformation(
    double &                    drs,         ///< surface transformation factor
    LINALG::Matrix<3,1> &       x_gp_lin,    ///< global coordiantes of gaussian point
    LINALG::Matrix<3,1> &       normal,      ///< normal vector on boundary cell
    GEO::CUT::BoundaryCell *    bc,          ///< boundary cell
    const LINALG::Matrix<2,1> & eta          ///< local coordinates of gaussian point w.r.t boundarycell
)
{

  normal.Clear();

  // get normal vector on linearized boundary cell, x-coordinates of gaussian point and surface transformation factor
  switch ( bc->Shape() )
  {
  case DRT::Element::tri3:
  {
    bc->Transform<DRT::Element::tri3>(eta, x_gp_lin, normal, drs);
    break;
  }
  case DRT::Element::quad4:
  {
    bc->Transform<DRT::Element::quad4>(eta, x_gp_lin, normal, drs);
    break;
  }
  default:
    throw std::runtime_error( "unsupported integration cell type" );
  }

  return;
}



/*----------------------------------------------------------------------*
 | evaluate shape functions and derivatives at given local coordinates  |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void FluidEleCalcXFEM<distype>::EvalFuncAndDeriv( LINALG::Matrix<3,1> &  rst )
{
  // evaluate shape functions
  DRT::UTILS::shape_function<distype>( rst, my::funct_ );

  // evaluate the derivatives of shape functions
  DRT::UTILS::shape_function_deriv1<distype>(rst,my::deriv_);
  my::xjm_.MultiplyNT(my::deriv_,my::xyze_);
  my::det_ = my::xji_.Invert(my::xjm_);

  // compute global first derivates
  my::derxy_.Multiply(my::xji_,my::deriv_);

  return;
}


/*----------------------------------------------------------------------*
 | build traction vector w.r.t fluid domain                             |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void FluidEleCalcXFEM<distype>::BuildTractionVector(
    LINALG::Matrix<my::nsd_,1> &  traction,   ///< traction vector
    double &                      press,      ///< pressure at gaussian point
    LINALG::Matrix<my::nsd_,1> &  normal      ///< normal vector
)
{
  // compute the stresses at the current Gaussian point for computing the interface force
  LINALG::Matrix<my::nsd_,my::nsd_> two_eps;
  for(int i=0; i<my::nsd_; ++i)
  {
    for(int j=0; j<my::nsd_; ++j)
    {
      two_eps(j,i) = my::vderxy_(i,j) + my::vderxy_(j,i);
    }
  }
  //-------------------------------

  // t = ( -pI + 2mu eps(u) )*n^f
  traction.Multiply(two_eps, normal);

  // add the pressure part and scale the viscous part with the viscosity
  traction.Update( -press, normal, viscaf_master_);

  return;
}

/*----------------------------------------------------------------------*
 | assemble side's interface force                                      |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void FluidEleCalcXFEM<distype>::AssembleInterfaceForce(
    Teuchos::RCP<Epetra_Vector>            iforcecol, ///< interface force column vector
    DRT::Discretization &                  cutdis,    ///< cut discretization
    std::vector<int> &                     lm,        ///< local dof map
    Epetra_SerialDenseVector &             iforce     ///< interface force vector
)
{
  //TEUCHOS_FUNC_TIME_MONITOR("FluidEleCalcXFEM::AssembleInterfaceForce");

  const Epetra_Map* dofcolmap = cutdis.DofColMap();

  for (int idof = 0; idof < (int)(lm.size()); ++idof)
  {
    int gdof = lm[idof];

    // f^i = ( N^i, t ) = ( N^i, (-pI+2mu*eps(u))*n )
    (*iforcecol)[dofcolmap->LID(gdof)] += iforce[idof];
  }

  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void FluidEleCalcXFEM<distype>::EvaluateNeumann(
  const double &                        timefacfac,             ///< theta*dt
  const LINALG::Matrix<my::nen_,1> &    funct_m,                ///< coupling master shape functions
  const LINALG::Matrix<my::nsd_,1> &    itraction_jump,         ///< prescribed interface traction, jump height for coupled problems
  Epetra_SerialDenseMatrix&             elevec1_epetra          ///< element vector
)
{
  const int master_numdof = my::nsd_+1;
  LINALG::Matrix<master_numdof*my::nen_,1> rhC_um(elevec1_epetra.A(), true);

  // funct_m * timefac * fac
  LINALG::Matrix<my::nen_,1> funct_m_timefacfac(funct_m);
  funct_m_timefacfac.Scale(timefacfac);

  //-----------------------------------------------------------------
  // standard consistency Neumann term

       /*            \
    - |    v  ,   t   |   with t = [sigma * n]
       \             /     */

  // loop over velocity components
  for (int ivel = 0; ivel < my::nsd_; ++ ivel)
  {
    //-----------------------------------------------
    //    - (vm, t)
    //-----------------------------------------------
    for (int ir = 0; ir<my::nen_; ++ir)
    {
      const unsigned row = ir*(my::nsd_+1) + ivel;
      rhC_um(row,0) += funct_m_timefacfac(ir)*itraction_jump(ivel);
    }
  } // end loop over velocity components
}


/*--------------------------------------------------------------------------------
 *--------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void FluidEleCalcXFEM<distype>::CalculateContinuityXFEM(
    DRT::ELEMENTS::Fluid *               ele,            ///< fluid element
    DRT::Discretization &                dis,            ///< discretization
    const std::vector<int> &             lm,             ///< local map
    Epetra_SerialDenseVector&            elevec1_epetra, ///< element vector
    const DRT::UTILS::GaussIntegration & intpoints       ///< integration points
  )
{
  LINALG::Matrix<my::numdofpernode_*my::nen_,1> elevec1(elevec1_epetra,true);
  LINALG::Matrix<my::numdofpernode_,my::nen_>   tmpvel;
  my::eid_ = ele->Id();

  // get node coordinates and number of elements per node
  GEO::fillInitialPositionArray<distype,my::nsd_,LINALG::Matrix<my::nsd_,my::nen_> >(ele,my::xyze_);

  //------------------------------------------------------------------------
  //  start loop over integration points
  //------------------------------------------------------------------------
  for ( DRT::UTILS::GaussIntegration::const_iterator iquad=intpoints.begin(); iquad!=intpoints.end(); ++iquad )
  {
    // evaluate shape functions and derivatives at integration point
    my::EvalShapeFuncAndDerivsAtIntPoint(iquad.Point(),iquad.Weight());

     // Summe über alle Knoten
     for(int ui=0; ui<my::nen_; ++ui)
     {
       for (int idim = 0; idim <my::nsd_; ++idim)
       {
         // Bloecke ueber Knoten
         const int fui = my::numdofpernode_*ui;
         /* continuity term */
         /*
              /           \
             |             |
             | nabla o Du  |
             |             |
              \           /
         */

         elevec1(fui+idim) += my::fac_*my::derxy_(idim,ui);
       }
     }
  }

  return;
}

/*--------------------------------------------------------------------------------
 *--------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void FluidEleCalcXFEM<distype>::CalculateContinuityXFEM(
    DRT::ELEMENTS::Fluid *      ele,                ///< fluid element
    DRT::Discretization &       dis,                ///< discretization
    const std::vector<int> &    lm,                 ///< local map
    Epetra_SerialDenseVector&   elevec1_epetra      ///< element vector
  )
{
  CalculateContinuityXFEM(ele,
                          dis,
                          lm,
                          elevec1_epetra,
                          my::intpoints_);
}

/*--------------------------------------------------------------------------------
 * compute characteristic element length
 *--------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
double FluidEleCalcXFEM<distype>::ComputeCharEleLength(
    DRT::Element *                                                        ele,                   ///< fluid element
    Epetra_SerialDenseMatrix &                                            ele_xyze,              ///< element coordinates
    const Teuchos::RCP<XFEM::ConditionManager> &                          cond_manager,          ///< XFEM condition manager
    const GEO::CUT::plain_volumecell_set &                                vcSet,                 ///< volumecell sets for volume integration
    const std::map<int, std::vector<GEO::CUT::BoundaryCell*> > &          bcells,                ///< bcells for boundary cell integration
    const std::map<int, std::vector<DRT::UTILS::GaussIntegration> > &     bintpoints,            ///< integration points for boundary cell integration
    Teuchos::RCP<DRT::ELEMENTS::XFLUID::SlaveElementInterface<distype> >  emb,                   ///< pointer to the embedded coupling implementation
    DRT::Element *                                                        face                   ///< side element in 3D
)
{
  //TEUCHOS_FUNC_TIME_MONITOR("FluidEleCalcXFEM::ComputeCharEleLength");

  const INPAR::XFEM::ViscStab_hk visc_stab_hk = fldparaxfem_->ViscStabHK();

  const int coup_sid = bintpoints.begin()->first;
  const INPAR::XFEM::AveragingStrategy averaging_strategy = cond_manager->GetAveragingStrategy(coup_sid,ele->Id());
  if (emb == Teuchos::null and averaging_strategy == INPAR::XFEM::Embedded_Sided)
    dserror("no coupling interface available, however Embedded_Sided coupling is activated!");

  // characteristic element length to be computed
  double h_k = 0.0;

  // measure of the face (surface in 3D, line in 2D) or measure of the cut-face
  double meas_surf = 0.0;

  // measure of the element volume or measure of the physical cut part of the element volume
  double meas_vol = 0.0;

  switch (visc_stab_hk)
  {
  //---------------------------------------------------
  // volume-equivalent diameter
  //---------------------------------------------------
  case INPAR::XFEM::ViscStab_hk_vol_equivalent:
  {
    // evaluate shape functions and derivatives at element center
    if(averaging_strategy == INPAR::XFEM::Embedded_Sided)
    {
      // evaluate shape functions and derivatives at element center w.r.t embedded element
      meas_vol = emb->EvalShapeFuncAndDerivsAtEleCenter();
    }
    else
    {
      // evaluate shape functions and derivatives at element center w.r.t background element
      my::EvalShapeFuncAndDerivsAtEleCenter();
      meas_vol = my::fac_;
    }

    // compute h_k as volume-equivalent diameter and directly return the value
    return h_k = ComputeVolEqDiameter(meas_vol);

    break;
  }
  //---------------------------------------------------
  // compute h_k as physical/cut volume divided by physical partial/cut surface measure
  // ( used to estimate the cut-dependent inverse estimate on cut elements, not useful for sliver and/or dotted cut situations)
  //---------------------------------------------------
  case INPAR::XFEM::ViscStab_hk_cut_vol_div_by_cut_surf:
  {
    if(averaging_strategy == INPAR::XFEM::Embedded_Sided)
      dserror("ViscStab_hk_cut_vol_div_by_cut_surf not reasonable for Embedded_Sided_Coupling!");

    // compute the cut surface measure
    meas_surf = ComputeMeasCutSurf(bintpoints, bcells);

    if (fabs(meas_surf) < 1.e-8)  dserror("Element contribution to interface has zero size.");

    // compute the cut volume measure
    for( GEO::CUT::plain_volumecell_set::const_iterator i=vcSet.begin();i!=vcSet.end();++i )
    {
      GEO::CUT::VolumeCell* vc = *i;
      meas_vol += vc->Volume();
    }

    if(meas_vol < 0.0) dserror(" measure of cut partial volume is smaller than 0.0: %f Attention with increasing Nitsche-Parameter!!!", meas_vol);

    break;
  }
  //---------------------------------------------------
  // full element volume divided by physical partial/cut surface measure ( used to estimate the cut-dependent inverse estimate on cut elements, however, avoids problems with sliver cuts, not useful for dotted cuts)
  //---------------------------------------------------
  case INPAR::XFEM::ViscStab_hk_ele_vol_div_by_cut_surf:
  {
    if(averaging_strategy == INPAR::XFEM::Embedded_Sided)
      dserror("ViscStab_hk_ele_vol_div_by_cut_surf not reasonable for Embedded_Sided_Coupling!");

    // compute the cut surface measure
    meas_surf = ComputeMeasCutSurf(bintpoints, bcells);

    // evaluate shape functions and derivatives at element center
    // compute the full element volume measure
    my::EvalShapeFuncAndDerivsAtEleCenter();
    meas_vol = my::fac_;

    break;
  }
  //---------------------------------------------------
  // full element volume divided by surface measure ( used for uncut situations, standard weak Dirichlet boundary/coupling conditions)
  //---------------------------------------------------
  case INPAR::XFEM::ViscStab_hk_ele_vol_div_by_ele_surf:
  {
    if(averaging_strategy != INPAR::XFEM::Embedded_Sided)
      dserror("ViscStab_hk_ele_vol_div_by_ele_surf just reasonable for Embedded_Sided_Coupling!");

    //---------------------------------------------------
    // find the corresponding local id of the face of the element
    // REMARK: this is quite slow, however at the moment the easiest way to get the local id
    //         here is space for improvement

    const int coup_idx = cond_manager->GetCouplingIndex(coup_sid, my::eid_);
    Teuchos::RCP<XFEM::MeshVolCoupling> mc_vol = Teuchos::rcp_dynamic_cast<XFEM::MeshVolCoupling>(cond_manager->GetCouplingByIdx(coup_idx));
    if (mc_vol == Teuchos::null) dserror("ComputeCharEleLength-ViscStab_hk_ele_vol_div_by_ele_surf: Cast to MeshVolCoupling failed!");

    const int lid = mc_vol->GetFaceLidOfEmbeddedElement(face->Id());
    //---------------------------------------------------

    // compute the uncut element's surface measure
    meas_surf = ComputeMeasFace(ele, ele_xyze, lid);

    // evaluate shape functions and derivatives at element center w.r.t embedded element
    // compute the full element volume measure
    meas_vol = emb->EvalShapeFuncAndDerivsAtEleCenter();

    break;
  }
  //---------------------------------------------------
  case INPAR::XFEM::ViscStab_hk_ele_vol_div_by_max_ele_surf:
  //---------------------------------------------------
  {
    if(averaging_strategy == INPAR::XFEM::Embedded_Sided)
      dserror("ViscStab_hk_ele_vol_div_by_max_ele_surf not reasonable for Embedded_Sided_Coupling!");

    // compute the uncut element's surface measure
    const int numfaces = DRT::UTILS::getNumberOfElementFaces(ele->Shape());

    // loop all surfaces
    for(int lid=0; lid< numfaces; ++lid)
    {
      meas_surf = std::max(meas_surf, ComputeMeasFace(ele, ele_xyze, lid));
    }

    // evaluate shape functions and derivatives at element center
    // compute the full element volume measure
    my::EvalShapeFuncAndDerivsAtEleCenter();
    meas_vol = my::fac_;

    break;
  }
  default:
    dserror("unknown type of characteristic element length");
    break;
  }

  //--------------------------------------
  // compute the final element length if fraction-based computation and not returned yet
  h_k = meas_vol / meas_surf;
  //--------------------------------------

  // check plausibility
  if(h_k < 1e-14) dserror("the characteristic element length is zero or smaller, it has not been set properly!");

  return h_k;
}

/*--------------------------------------------------------------------------------
 * compute the measure of the elements surface with given local id
 *--------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
double FluidEleCalcXFEM<distype>::ComputeMeasFace(
    DRT::Element *            ele,              ///< fluid element
    Epetra_SerialDenseMatrix& ele_xyze,         ///< element coordinates
    const int                 local_face_id     ///< the local id of the face w.r.t the fluid element
)
{
  // get the shape of the face
  DRT::Element::DiscretizationType face_shape = DRT::UTILS::getEleFaceShapeType(ele->Shape(), local_face_id);

  // get the current node coordinates, extract them from the element's node coordinates
  const int numnode_face = DRT::UTILS::getNumberOfElementNodes(face_shape);
  Epetra_SerialDenseMatrix xyze_face( my::nsd_, numnode_face);

  // map for numbering of nodes of the surfaces
  std::vector< std::vector<int> > map = DRT::UTILS::getEleNodeNumberingFaces(ele->Shape());

  // extract the surface's node coordinates from the element's nodes coordinates
  for(int n=0; n<numnode_face; ++n)
  {
    const int node_lid = map[local_face_id][n];
    for(int idim = 0; idim < my::nsd_; ++idim)
      xyze_face(idim,n) = ele_xyze(idim,node_lid);
  }

  // the metric tensor and the area of an infintesimal surface element
  Epetra_SerialDenseMatrix  metrictensor(my::nsd_-1,my::nsd_-1);
  double                    drs = 0.0;

  if(my::nsd_ != 3) dserror("don't call this function for non-3D examples, adapt the following for 2D!");

  DRT::UTILS::GaussRule2D gaussrule = DRT::UTILS::intrule2D_undefined;
  switch(face_shape)
  {
  case DRT::Element::quad4:
  case DRT::Element::quad8:
  case DRT::Element::quad9:
    gaussrule = DRT::UTILS::intrule_quad_1point;
    break;
  case DRT::Element::tri3:
  case DRT::Element::tri6:
    gaussrule = DRT::UTILS::intrule_tri_1point;
    break;
  default:
    dserror("shape type unknown!\n"); break;
  }

  double meas_face = 0.0;

  /*----------------------------------------------------------------------*
    |               start loop over integration points                     |
   *----------------------------------------------------------------------*/
  const DRT::UTILS::IntegrationPoints2D  intpoints(gaussrule);
  for (int gpid=0; gpid<intpoints.nquad; ++gpid)
  {
    const double e0 = intpoints.qxg[gpid][0];
    const double e1 = intpoints.qxg[gpid][1];

    Epetra_SerialDenseMatrix deriv( my::nsd_-1, numnode_face);

    // get shape functions and derivatives in the plane of the element
    DRT::UTILS::shape_function_2D_deriv1(deriv, e0, e1, face_shape);

    // compute measure tensor for surface element and the infinitesimal
    // area element drs for the integration
    DRT::UTILS::ComputeMetricTensorForSurface(xyze_face,deriv,metrictensor,&drs);

    meas_face += intpoints.qwgt[gpid] * drs;

  }

  return meas_face;
}

/*--------------------------------------------------------------------------------
 * pre-compute the measure of all side's surface cutting the element
 *--------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
double FluidEleCalcXFEM<distype>::ComputeMeasCutSurf(
    const std::map<int, std::vector<DRT::UTILS::GaussIntegration> > &   bintpoints,        ///< boundary cell integration points
    const std::map<int, std::vector<GEO::CUT::BoundaryCell*> > &        bcells             ///< boundary cells
    )
{
  double surf = 0.0;

  //--------------------------------------------
  // loop intersecting sides
  // map of side-element id and Gauss points
  for ( std::map<int, std::vector<DRT::UTILS::GaussIntegration> >::const_iterator i=bintpoints.begin();
        i!=bintpoints.end();
        ++i )
  {
    int sid = i->first;
    const std::vector<DRT::UTILS::GaussIntegration> & cutintpoints = i->second;

    // get side's boundary cells
    std::map<int, std::vector<GEO::CUT::BoundaryCell*> >::const_iterator j = bcells.find( sid );
    if ( j==bcells.end() )
      dserror( "missing boundary cell" );

    const std::vector<GEO::CUT::BoundaryCell*> & bcs = j->second;
    if ( bcs.size()!=cutintpoints.size() )
      dserror( "boundary cell integration rules mismatch" );

    //--------------------------------------------
    // loop boundary cells w.r.t current cut side
    //--------------------------------------------
    for ( std::vector<DRT::UTILS::GaussIntegration>::const_iterator i=cutintpoints.begin();
          i!=cutintpoints.end();
          ++i )
    {
      const DRT::UTILS::GaussIntegration & gi = *i;
      GEO::CUT::BoundaryCell * bc = bcs[i - cutintpoints.begin()]; // get the corresponding boundary cell

      //--------------------------------------------
      // loop gausspoints w.r.t current boundary cell
      //--------------------------------------------
      for ( DRT::UTILS::GaussIntegration::iterator iquad=gi.begin(); iquad!=gi.end(); ++iquad )
      {
        double drs = 0.0; // transformation factor between reference cell and linearized boundary cell

        const LINALG::Matrix<2,1> eta( iquad.Point() ); // xi-coordinates with respect to side

        LINALG::Matrix<3,1> normal(true);

        LINALG::Matrix<3,1> x_gp_lin(true); // gp in xyz-system on linearized interface

        // compute transformation factor, normal vector and global Gauss point coordiantes
        if(bc->Shape() != DRT::Element::dis_none) // Tessellation approach
        {
          ComputeSurfaceTransformation(drs, x_gp_lin, normal, bc, eta);
        }
        else // MomentFitting approach
        {
          drs = 1.0;
          normal = bc->GetNormalVector();
          const double* gpcord = iquad.Point();
          for (int idim=0;idim<3;++idim)
          {
            x_gp_lin(idim,0) = gpcord[idim];
          }
        }

        const double surf_fac = drs*iquad.Weight();

        surf += surf_fac;

      } //loop gausspoints w.r.t current boundary cell
    } // loop boundary cells
  } // loop intersecting sides

  return surf;
}

template <DRT::Element::DiscretizationType distype>
void FluidEleCalcXFEM<distype>::GetMaterialParametersVolumeCell( Teuchos::RCP<const MAT::Material>  material,
    double &                           densaf, //done
    double &                           viscaf,   //done
    double &                           gamma   //done
    )
{

  //TEUCHOS_FUNC_TIME_MONITOR("FluidEleCalcXFEM::GetMaterialParametersVolumeCell");

  //Initiate dummy variables:
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // get element-wise velocity/pressure field for current time step
  my::evelaf_.Clear();
  //Scatra field
  my::escabofoaf_.Clear();
  my::escaaf_.Clear();
  my::escaam_.Clear();
  // set thermodynamic pressure at n+1/n+alpha_F and n+alpha_M/n and
  // its time derivative at n+alpha_M/n+1 (LOMA specific!!)
  const double thermpressaf   = 1.0;
  const double thermpressam   = 1.0;
  const double thermpressdtaf = 0.0;
  const double thermpressdtam = 0.0;
  const double vol            = 0.0;
  //Values of material parameters at not needed time steps.
  double densn  = 0.0;
  double densam = 0.0;
  double viscn  = 0.0;
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  //Get material from FluidEleCalc routine.
  // If non-constant density or viscosity wants to be calculated on the different sides. A review over how the scalar variables
  // are set at the surface should be made.
  my::GetMaterialParams(material,my::evelaf_,my::escaaf_,my::escaam_,my::escabofoaf_,thermpressaf,thermpressam,thermpressdtaf,thermpressdtam,vol,densam,densaf,densn,viscaf,viscn,gamma);

  return;

}

  } // end namespace ELEMENTS
} // end namespace DRT



// Ursula is responsible for this comment!
template class DRT::ELEMENTS::FluidEleCalcXFEM<DRT::Element::hex8>;
template class DRT::ELEMENTS::FluidEleCalcXFEM<DRT::Element::hex20>;
template class DRT::ELEMENTS::FluidEleCalcXFEM<DRT::Element::hex27>;
template class DRT::ELEMENTS::FluidEleCalcXFEM<DRT::Element::tet4>;
template class DRT::ELEMENTS::FluidEleCalcXFEM<DRT::Element::tet10>;
template class DRT::ELEMENTS::FluidEleCalcXFEM<DRT::Element::wedge6>;
//template class DRT::ELEMENTS::FluidEleCalcXFEM<DRT::Element::pyramid5>;


