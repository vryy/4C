/*----------------------------------------------------------------------*/
/*!
\file fluid_ele_calc_xfem.cpp

\brief Internal implementation of XFluid element interface coupling

<pre>
Maintainer: Raffaela Kruse /Benedikt Schott
            kruse@lnm.mw.tum.de
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
  : DRT::ELEMENTS::FluidEleCalc<distype>::FluidEleCalc()
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

  const INPAR::CUT::VCellGaussPts vcellgausspts = fldparaxfem_->VolumeCellGaussPoints();

  if( vcellgausspts == INPAR::CUT::VCellGaussPts_Tessellation ) // standard "Tessellation"
  {
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
  }
  else if( vcellgausspts==INPAR::CUT::VCellGaussPts_MomentFitting ) // standard "MomentFitting" method
  {
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
  }
  else if( vcellgausspts==INPAR::CUT::VCellGaussPts_DirectDivergence)  // DirectDivergence approach
  {

    LINALG::Matrix<my::numdofpernode_*my::nen_,my::numdofpernode_*my::nen_> elemat1(elemat1_epetra,true);
    //LINALG::Matrix<(nsd_+1)*nen_,(nsd_+1)*nen_> elemat2(elemat2_epetra,true);
    LINALG::Matrix<my::numdofpernode_*my::nen_,            1> elevec1(elevec1_epetra,true);

    for( std::vector<DRT::UTILS::GaussIntegration>::const_iterator i=intpoints.begin();i!=intpoints.end();++i )
    {
      const DRT::UTILS::GaussIntegration intcell = *i;
      GEO::CUT::VolumeCell * vc = cells[i-intpoints.begin()];

      // This happens when the volume of the cell is too small that the integration method
      // predicts negative volume for this
      if( vc->IsNegligiblySmall() )
        continue;

      //----------------------------------------------------------------------
      //integration over the main gauss points to get the required integral
      //----------------------------------------------------------------------
      int mainPtno = 0;
      for ( DRT::UTILS::GaussIntegration::iterator iquad=intcell.begin(); iquad!=intcell.end(); ++iquad )
      {
        Epetra_SerialDenseMatrix elematTemp1(my::numdofpernode_*my::nen_,my::numdofpernode_*my::nen_); //Sudhakar : Is this efficient?
        Epetra_SerialDenseMatrix elematTemp2(my::numdofpernode_*my::nen_,my::numdofpernode_*my::nen_);
        Epetra_SerialDenseVector elevecTemp1(my::numdofpernode_*my::nen_);

        // get internal Gaussian rule for every main Gauss point
        DRT::UTILS::GaussIntegration gint = vc->GetInternalRule( mainPtno );
        mainPtno++;

        //----------------------------------------------------------------------
        //integration over the internal gauss points - to get modified integrand
        //----------------------------------------------------------------------

        err = my::Evaluate( ele, discretization, lm, params, mat,
                        elematTemp1, elematTemp2,
                        elevecTemp1, elevec2_epetra, elevec3_epetra,
                        gint, offdiag );


        if(err)
          return err;

        LINALG::Matrix<my::numdofpernode_*my::nen_,my::numdofpernode_*my::nen_> elem1(elematTemp1,true);
        LINALG::Matrix<my::numdofpernode_*my::nen_,my::numdofpernode_*my::nen_> elem2(elematTemp2,true);
        LINALG::Matrix<my::numdofpernode_*my::nen_,1> elev1(elevecTemp1,true);

        elemat1.Update(iquad.Weight(), elem1, 1.0);
        //elemat2.Update(1.0, elem2, 1.0);
        elevec1.Update(iquad.Weight(), elev1, 1.0);
      }
    }
  }
  else dserror("unsupported type of VCellGaussPts");

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

  const INPAR::CUT::VCellGaussPts vcellgausspts = fldparaxfem_->VolumeCellGaussPoints();

  if( vcellgausspts==INPAR::CUT::VCellGaussPts_Tessellation ) // standard "Tessellation" or "MomentFitting" method
  {
    for( std::vector<DRT::UTILS::GaussIntegration>::const_iterator i=intpoints.begin();i!=intpoints.end();++i )
    {
      const DRT::UTILS::GaussIntegration gint = *i;
      err = my::IntegrateShapeFunction( ele, discretization, lm,
                     elevec1_epetra,
                     gint);
      if(err)
        return err;
    }
  }
  else if( vcellgausspts==INPAR::CUT::VCellGaussPts_MomentFitting ) // standard "Tessellation" or "MomentFitting" method
  {
    for( std::vector<DRT::UTILS::GaussIntegration>::const_iterator i=intpoints.begin();i!=intpoints.end();++i )
    {
      const DRT::UTILS::GaussIntegration gint = *i;
      err = my::IntegrateShapeFunction( ele, discretization, lm,
                     elevec1_epetra,
                     gint);
      if(err)
        return err;
    }
  }
  else if( vcellgausspts==INPAR::CUT::VCellGaussPts_DirectDivergence )  // DirectDivergence approach
  {

    LINALG::Matrix<my::numdofpernode_*my::nen_,            1> elevec1(elevec1_epetra,true);

    for( std::vector<DRT::UTILS::GaussIntegration>::const_iterator i=intpoints.begin();i!=intpoints.end();++i )
    {
      const DRT::UTILS::GaussIntegration intcell = *i;
      GEO::CUT::VolumeCell * vc = cells[i-intpoints.begin()];

      //----------------------------------------------------------------------
      //integration over the main gauss points to get the required integral
      //----------------------------------------------------------------------
      int mainPtno = 0;
      for ( DRT::UTILS::GaussIntegration::iterator iquad=intcell.begin(); iquad!=intcell.end(); ++iquad )
      {
        Epetra_SerialDenseVector elevecTemp1(my::numdofpernode_*my::nen_);

        // get internal Gaussian rule for every main Gauss point
        DRT::UTILS::GaussIntegration gint = vc->GetInternalRule( mainPtno );
        mainPtno++;

        //----------------------------------------------------------------------
        //integration over the internal gauss points - to get modified integrand
        //----------------------------------------------------------------------

        err = my::IntegrateShapeFunction( ele, discretization, lm, elevecTemp1, gint);


        if(err)
          return err;

        LINALG::Matrix<my::numdofpernode_*my::nen_,            1> elev1(elevecTemp1,true);

        elevec1.Update(iquad.Weight(), elev1, 1.0);
      }
    }
  }
  else dserror("unsupported type of VCellGaussPts");

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
    // 4.   || p - p_h ||_L2(Omega)              =   standard L2-norm for for pressure
    //
    // viscosity-scaled domain errors
    // 5.   || nu^(+1/2) grad( u - u_h ) ||_L2(Omega)      =   visc-scaled H1-seminorm for velocity
    //                                                     =   nu^(+1/2) * || grad( u - u_h ) ||_L2(Omega) (for homogeneous visc)
    // 6.   || nu^(-1/2) (p - p_h) ||_L2(Omega)            =   visc-scaled L2-norm for for pressure
    //                                                     =   nu^(-1/2) * || p - p_h ||_L2(Omega) (for homogeneous visc)

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

    // standard domain errors
    ele_dom_norms[0] += u_err_squared;
    ele_dom_norms[1] += grad_u_err_squared;
    ele_dom_norms[2] += u_err_squared + grad_u_err_squared;
    ele_dom_norms[3] += p_err_squared;

    // viscosity-scaled domain errors
    ele_dom_norms[4] += my::visc_ * grad_u_err_squared;
    ele_dom_norms[5] += 1.0/my::visc_ * p_err_squared;

    // error for predefined functional
    ele_dom_norms[6] += funcerr;


  } // loop gaussian points

  return 0;
}

template <DRT::Element::DiscretizationType distype>
void FluidEleCalcXFEM<distype>::AnalyticalReference(
    const int                               calcerr,     ///< which reference solution
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

    bool is_stationary = false;

    // evaluate velocity and pressure
    // evaluate the velocity gradient
    if(calcerr == INPAR::FLUID::beltrami_stat_stokes or
       calcerr == INPAR::FLUID::beltrami_stat_navier_stokes)
    {
      is_stationary = true;
    }
    else if(calcerr == INPAR::FLUID::beltrami_instat_stokes or
            calcerr == INPAR::FLUID::beltrami_instat_navier_stokes)
    {
      is_stationary = false;
    }

    function      = Teuchos::rcp(new DRT::UTILS::BeltramiUP( mat, is_stationary));
    function_grad = Teuchos::rcp(new DRT::UTILS::BeltramiGradU( mat, is_stationary));

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

  case INPAR::FLUID::byfunct1:
  {
    const int func_no = 1;


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
      const double u_exact_x = DRT::Problem::Instance()->Funct(func_no-1).Evaluate(0,position,t,NULL);
      const double u_exact_y = DRT::Problem::Instance()->Funct(func_no-1).Evaluate(1,position,t,NULL);
      const double p_exact   = DRT::Problem::Instance()->Funct(func_no-1).Evaluate(2,position,t,NULL);

      u(0) = u_exact_x;
      u(1) = u_exact_y;
      p    = p_exact;
    }
    else if(my::nsd_==3)
    {
      const double u_exact_x = DRT::Problem::Instance()->Funct(func_no-1).Evaluate(0,position,t,NULL);
      const double u_exact_y = DRT::Problem::Instance()->Funct(func_no-1).Evaluate(1,position,t,NULL);
      const double u_exact_z = DRT::Problem::Instance()->Funct(func_no-1).Evaluate(2,position,t,NULL);
      const double p_exact   = DRT::Problem::Instance()->Funct(func_no-1).Evaluate(3,position,t,NULL);

      u(0) = u_exact_x;
      u(1) = u_exact_y;
      u(2) = u_exact_z;
      p    = p_exact;

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
    Teuchos::RCP<MAT::Material>&                                        mat,               ///< material
    Epetra_SerialDenseVector&                                           ele_interf_norms,  /// squared element interface norms
    DRT::Discretization &                                               cutdis,            ///< cut discretization
    const std::map<int, std::vector<GEO::CUT::BoundaryCell*> > &        bcells,            ///< boundary cells
    const std::map<int, std::vector<DRT::UTILS::GaussIntegration> > &   bintpoints,        ///< boundary integration points
    Teuchos::ParameterList&                                             params,            ///< parameter list
    const GEO::CUT::plain_volumecell_set&                               vcSet              ///< volumecell sets in this element
)
{
  const int calcerr = DRT::INPUT::get<INPAR::FLUID::CalcError>(params,"calculate error");

  const double t = my::fldparatimint_->Time();


  //----------------------------------------------------------------------------
  //                         ELEMENT GEOMETRY
  //----------------------------------------------------------------------------

  // get node coordinates
  GEO::fillInitialPositionArray< distype, my::nsd_, LINALG::Matrix<my::nsd_,my::nen_> >( ele, my::xyze_ );
  // evaluate shape functions and derivatives at element center
  my::EvalShapeFuncAndDerivsAtEleCenter();
  // set element area or volume
  const double vol = my::fac_;

  // get element-wise velocity/pressure field
  LINALG::Matrix<my::nsd_,my::nen_> evelaf(true);
  LINALG::Matrix<my::nen_,1> epreaf(true);
  my::ExtractValuesFromGlobalVector(dis, lm, *my::rotsymmpbc_, &evelaf, &epreaf, "u and p at time n+1 (converged)");


  //----------------------------------------------------------------------------
  //      surface integral --- build Cuiui, Cuui, Cuiu and Cuu matrix and rhs
  //----------------------------------------------------------------------------

  DRT::Element::LocationArray cutla( 1 );

  LINALG::Matrix<3,1> normal;
  LINALG::Matrix<3,1> x_side;

  // coupling implementation between background element and each cut side (std::map<sid, side_impl)
  std::map<int, Teuchos::RCP<DRT::ELEMENTS::XFLUID::SlaveElementInterface<distype> > > side_impl;
  Teuchos::RCP<DRT::ELEMENTS::XFLUID::SlaveElementInterface<distype> > si;

  // find all the intersecting elements of actele
  std::set<int> begids;
  for (std::map<int,  std::vector<GEO::CUT::BoundaryCell*> >::const_iterator bc=bcells.begin();
       bc!=bcells.end(); ++bc )
  {
    int sid = bc->first;
    begids.insert(sid);
  }

  // map of boundary element gids and coupling matrices, [0]: Cuiui matrix
  std::map<int, std::vector<Epetra_SerialDenseMatrix> > Cuiui_coupling;

  // lm vector of all intersecting boundary elements that intersect the current background element
  std::vector<int> patchelementslmv;
  std::vector<int> patchelementslmowner;

  // create location vectors for intersecting boundary elements and reshape coupling matrices
  PatchLocationVector(begids,cutdis,patchelementslmv,patchelementslmowner, Cuiui_coupling);


  //-----------------------------------------------------------------------------------
  //         evaluate element length, stabilization factors and average weights
  //-----------------------------------------------------------------------------------

  // element length
  double h_k = 0.0;

  // take a volume based element length
  h_k = ComputeVolEqDiameter(vol);


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
    int sid = i->first;
    const std::vector<DRT::UTILS::GaussIntegration> & cutintpoints = i->second;

    // get side's boundary cells
    std::map<int, std::vector<GEO::CUT::BoundaryCell*> >::const_iterator j = bcells.find( sid );
    if ( j==bcells.end() )
      dserror( "missing boundary cell" );

    const std::vector<GEO::CUT::BoundaryCell*> & bcs = j->second;
    if ( bcs.size()!=cutintpoints.size() )
      dserror( "boundary cell integration rules mismatch" );

    // side and location vector
    DRT::Element * side = cutdis.gElement( sid );
    side->LocationVector(cutdis,cutla,false);

    // side geometry
    const int numnodes = side->NumNode();
    DRT::Node ** nodes = side->Nodes();
    Epetra_SerialDenseMatrix side_xyze( 3, numnodes );
    for ( int i=0; i<numnodes; ++i )
    {
      const double * x = nodes[i]->X();
      std::copy( x, x+3, &side_xyze( 0, i ) );
    }

    si = DRT::ELEMENTS::XFLUID::NitscheInterface<distype>::CreateSlaveElementRepresentation(side,side_xyze);

    side_impl[sid] = si;

    // get velocity at integration point of boundary dis
    si->SetSlaveVel(cutdis,"ivelnp",cutla[0].lm_);

    // set displacement of side
    si->AddSlaveEleDisp(cutdis,"idispnp",cutla[0].lm_);


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

#ifdef BOUNDARYCELL_TRANSFORMATION_OLD

        si->Evaluate(eta,x_side,normal,drs);

        // find element local position of gauss point at interface
        GEO::CUT::Position<distype> pos( my::xyze_, x_side );
        pos.Compute();
        rst = pos.LocalCoordinates();

#else
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
          for (int idim=0;idim<3;idim++)
          {
            x_gp_lin(idim,0) = gpcord[idim];
          }
        }

        // find element local position of gauss point
        GEO::CUT::Position<distype> pos( my::xyze_, x_gp_lin );
        pos.Compute();
        rst = pos.LocalCoordinates();

        // project gaussian point from linearized interface to warped side (get/set local side coordinates in SideImpl)
        LINALG::Matrix<2,1> xi_side(true);
        si->ProjectOnSide(x_gp_lin, x_side, xi_side);
 #endif

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

        //--------------------------------------------
        // compute errors

        LINALG::Matrix<my::nsd_,1>        u_analyt(true);      // boundary condition to enforce (xfsi), interfacial jump to enforce (fluidfluid)
        LINALG::Matrix<my::nsd_,my::nsd_> grad_u_analyt(true);
        p_analyt = 0.0;

        AnalyticalReference(
                     calcerr,          ///< which reference solution
                     u_analyt,         ///< exact velocity (onesided), exact jump vector (coupled)
                     grad_u_analyt,    ///< exact velocity gradient
                     p_analyt,         ///< exact pressure
                     x_side,           ///< xyz position of gaussian point which lies on the real side, projected from linearized interface
                     t,                ///< time t
                     mat
                     );

        //--------------------------------------------

        // prescribed velocity vector at weak Dirichlet boundary
        LINALG::Matrix<my::nsd_,1> ivelint_WDBC_JUMP(true);
        si->GetInterfaceVel(ivelint_WDBC_JUMP);

        LINALG::Matrix<my::nsd_,1> check_diff(true);
        check_diff.Update(1.0, u_analyt, -1.0, ivelint_WDBC_JUMP);


        // This check is obsolete, since there is a interpolation difference
        // between analytical solution and the enforced boundary condition
        if(check_diff.Norm2() > 1e-12)
        {
//            cout.precision(12);
//            cout << "u_analyt" << u_analyt << endl;
//            cout << "ivelint_WDBC_JUMP" <<  ivelint_WDBC_JUMP << endl;
          //dserror("do you want to enforce other boundary conditions than the analytical domain solution?");
        }

        u_err.Update(1.0, my::velint_, -1.0, u_analyt, 0.0);


        grad_u_err.Update(1.0, my::vderxy_, -1.0, grad_u_analyt, 0.0);
        p_err = press - p_analyt;

        flux_u_err.Multiply(grad_u_err,normal);
        flux_p_err.Update(p_err,normal,0.0);


        // interface errors
        // 1.   || nu^(+1/2) (u - u*) ||_H1/2(Gamma)             =  broken H1/2 Sobolev norm for boundary/coupling condition
        // 2.   || nu^(+1/2) grad( u - u_h )*n ||_H-1/2(Gamma)   =  standard H-1/2 Sobolev norm for normal flux (velocity part)
        // 3.   || nu^(-1/2) (p - p_h)*n ||_H-1/2(Gamma)         =  standard H-1/2 Sobolev norm for normal flux (pressure part)

        double u_err_squared      = 0.0;
        double flux_u_err_squared = 0.0;
        double flux_p_err_squared = 0.0;

        // evaluate squared errors at gaussian point
        for (int isd=0;isd<my::nsd_;isd++)
        {
          u_err_squared += u_err(isd)*u_err(isd)*surf_fac;
          flux_u_err_squared += flux_u_err(isd)*flux_u_err(isd)*surf_fac;
          flux_p_err_squared += flux_p_err(isd)*flux_p_err(isd)*surf_fac;
        }

        // interface errors
        ele_interf_norms[0] += 1.0/h_k * my::visc_ * u_err_squared;
        ele_interf_norms[1] += h_k     * my::visc_ * flux_u_err_squared;
        ele_interf_norms[2] += h_k     / my::visc_ * flux_p_err_squared;

      } // end loop gauss points of boundary cell
    } // end loop boundary cells of side

  } // end loop cut sides

  return 0;
}
/*--------------------------------------------------------------------------------
 * compute interface error norms for XFF-problems
 *--------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int FluidEleCalcXFEM<distype>::ComputeErrorInterfaceXFluidFluid(
    DRT::ELEMENTS::Fluid *                                              ele,               ///< fluid element
    DRT::Discretization &                                               dis,               ///< background discretization
    const std::vector<int> &                                            lm,                ///< element local map
    Teuchos::RCP<MAT::Material>&                                        mat,               ///< material
    Epetra_SerialDenseVector&                                           ele_interf_norms,  /// squared element interface norms
    DRT::Discretization &                                               cutdis,            ///< cut discretization
    DRT::Discretization &                                               embdis,            ///< embedded discretization
    const std::map<int, std::vector<GEO::CUT::BoundaryCell*> > &        bcells,            ///< boundary cells
    const std::map<int, std::vector<DRT::UTILS::GaussIntegration> > &   bintpoints,        ///< boundary integration points
    Teuchos::ParameterList&                                             params,            ///< parameter list
    const GEO::CUT::plain_volumecell_set&                               vcSet,              ///< volumecell sets in this element,
    std::map<int,int> &                                                 boundary_emb_gid_map
)
{
  const int calcerr = DRT::INPUT::get<INPAR::FLUID::CalcError>(params,"calculate error");

  const double t = my::fldparatimint_->Time();

  //----------------------------------------------------------------------------
  //                         ELEMENT GEOMETRY
  //----------------------------------------------------------------------------

  // get node coordinates
  GEO::fillInitialPositionArray< distype, my::nsd_, LINALG::Matrix<my::nsd_,my::nen_> >( ele, my::xyze_ );
  // evaluate shape functions and derivatives at element center
  my::EvalShapeFuncAndDerivsAtEleCenter();
  // set element area or volume
  const double vol = my::fac_;

  // get element-wise velocity/pressure field
  LINALG::Matrix<my::nsd_,my::nen_> evelaf(true);
  LINALG::Matrix<my::nen_,1> epreaf(true);
  my::ExtractValuesFromGlobalVector(dis, lm, *my::rotsymmpbc_, &evelaf, &epreaf, "u and p at time n+1 (converged)");

  //----------------------------------------------------------------------------
  //      surface integral --- build Cuiui, Cuui, Cuiu and Cuu matrix and rhs
  //----------------------------------------------------------------------------

  DRT::Element::LocationArray alela( 1 );
  DRT::Element::LocationArray cutla( 1 );

  LINALG::Matrix<3,1> normal;
  LINALG::Matrix<3,1> x_side;

  // interface for coupling between background element and cutting embedded element
  Teuchos::RCP<DRT::ELEMENTS::XFLUID::SlaveElementInterface<distype> > emb;
  // auxiliary 2D coupling interface (side)
  Teuchos::RCP<DRT::ELEMENTS::XFLUID::SlaveElementInterface<distype> > si;

  // find all the intersecting elements of actele
  std::set<int> begids;
  for (std::map<int,  std::vector<GEO::CUT::BoundaryCell*> >::const_iterator bc=bcells.begin();
       bc!=bcells.end(); ++bc )
  {
    int sid = bc->first;
    begids.insert(sid);
  }

  //-----------------------------------------------------------------------------------
  //         evaluate element length, stabilization factors and average weights
  //-----------------------------------------------------------------------------------

  // element length
  double h_k = 0.0;

  // take a volume based element length
  h_k = ComputeVolEqDiameter(vol);

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
    int sid = i->first;
    const std::vector<DRT::UTILS::GaussIntegration> & cutintpoints = i->second;

    // get side's boundary cells
    std::map<int, std::vector<GEO::CUT::BoundaryCell*> >::const_iterator j = bcells.find( sid );
    if ( j==bcells.end() )
      dserror( "missing boundary cell" );

    const std::vector<GEO::CUT::BoundaryCell*> & bcs = j->second;
    if ( bcs.size()!=cutintpoints.size() )
      dserror( "boundary cell integration rules mismatch" );

    // side and location vector
    DRT::Element * side = cutdis.gElement( sid );
    side->LocationVector(cutdis,cutla,false);
    // the corresponding embedded element
    DRT::Element * emb_ele = embdis.gElement( boundary_emb_gid_map.find(sid)->second );
    emb_ele->LocationVector(embdis,alela,false);

    // embedded geometry
    const int emb_numnodes = emb_ele->NumNode();
    DRT::Node ** emb_nodes = emb_ele->Nodes();
    Epetra_SerialDenseMatrix emb_xyze( 3, emb_numnodes );
    for ( int i=0; i<emb_numnodes; ++i )
    {
      const double * x = emb_nodes[i]->X();
      std::copy( x, x+3, &emb_xyze( 0, i ) );
    }

    // side geometry
    const int numnodes = side->NumNode();
    DRT::Node ** nodes = side->Nodes();
    Epetra_SerialDenseMatrix side_xyze( 3, numnodes );
    for ( int i=0; i<numnodes; ++i )
    {
      const double * x = nodes[i]->X();
      std::copy( x, x+3, &side_xyze( 0, i ) );
    }

    // create coupling interface class
    emb = DRT::ELEMENTS::XFLUID::SlaveElementInterface<distype>::CreateSlaveElementRepresentation(emb_ele,emb_xyze);
    si = DRT::ELEMENTS::XFLUID::SlaveElementInterface<distype>::CreateSlaveElementRepresentation(side,side_xyze);

    // set side element velocity at integration point of boundary dis
    si->SetSlaveVel(cutdis,"ivelnp",cutla[0].lm_);
    // set embedded element velocity at integration point
    emb->SetSlaveVel(embdis,"velaf",alela[0].lm_);

    // set displacement of side
    si->AddSlaveEleDisp(cutdis,"idispnp",cutla[0].lm_);

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

#ifdef BOUNDARYCELL_TRANSFORMATION_OLD

        si->Evaluate(eta,x_side,normal,drs);

        // find element local position of gauss point at interface
        GEO::CUT::Position<distype> pos( my::xyze_, x_side );
        pos.Compute();
        rst = pos.LocalCoordinates();

#else
        LINALG::Matrix<3,1> x_gp_lin(true); // gp in xyz-system on linearized interface

        // compute transformation factor, normal vector and global Gauss point coordiantes
        if (bc->Shape() != DRT::Element::dis_none) // Tessellation approach
        {
          ComputeSurfaceTransformation(drs, x_gp_lin, normal, bc, eta);
        }
        else // MomentFitting approach
        {
          drs = 1.0;
          normal = bc->GetNormalVector();
          const double* gpcord = iquad.Point();
          for (int idim=0;idim<3;idim++)
          {
            x_gp_lin(idim,0) = gpcord[idim];
          }
        }

        // find element local position of gauss point
        GEO::CUT::Position<distype> pos( my::xyze_, x_gp_lin );
        pos.Compute();
        rst = pos.LocalCoordinates();

        // project gaussian point from linearized interface to warped side (get/set local side coordinates in SideImpl)
        LINALG::Matrix<2,1> xi_side(true);
        si->ProjectOnSide(x_gp_lin, x_side, xi_side);
 #endif

        const double surf_fac = drs*iquad.Weight();

        // -----------------------------------------
        // evaluate embedded element shape functions
        emb->Evaluate( x_side );

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

        //--------------------------------------------
        // compute errors

        LINALG::Matrix<my::nsd_,1>        u_analyt(true);      // boundary condition to enforce (xfsi), interfacial jump to enforce (fluidfluid)
        LINALG::Matrix<my::nsd_,my::nsd_> grad_u_analyt(true);
        p_analyt = 0.0;


        AnalyticalReference(
                     calcerr,          ///< which reference solution
                     u_analyt,         ///< exact velocity (onesided), exact jump vector (coupled)
                     grad_u_analyt,    ///< exact velocity gradient
                     p_analyt,         ///< exact pressure
                     x_side,           ///< xyz position of gaussian point which lies on the real side, projected from linearized interface
                     t,                ///< time t
                     mat);

        //--------------------------------------------
        // zero velocity jump for fluidfluidcoupling
        LINALG::Matrix<my::nsd_,1> ivelint_WDBC_JUMP(true);

        // from side element (not embedded!)
        LINALG::Matrix<my::nsd_,1> ivelint(true);
        si->GetInterfaceVel(ivelint);

        u_err.Update(1.0, my::velint_, -1.0, ivelint, 0.0); // u_backgr - u_emb
        u_err.Update(-1.0, ivelint_WDBC_JUMP, 1.0);         // u_backgr - u_emb - u_jump

        LINALG::Matrix<my::nsd_,my::nsd_> grad_u_side(true);
        emb->GetInterfaceVelGrad(grad_u_side);

        grad_u_err.Update(1.0, my::vderxy_, -1.0, grad_u_side, 0.0);

        double press_emb = 0.0;
        emb->GetInterfacePres(press_emb);
        //p_err = p_background - p_emb;
        p_err = press - press_emb;

        flux_u_err.Multiply(grad_u_err,normal);
        flux_p_err.Update(p_err,normal,0.0);

        // interface errors
        // 1.   || nu^(+1/2) (u_b - u_e - u_jump) ||_H1/2(Gamma)             =  broken H1/2 Sobolev norm for boundary/coupling condition
        // 2.   || nu^(+1/2) grad( u_b - u_e )*n ||_H-1/2(Gamma)   =  standard H-1/2 Sobolev norm for normal flux (velocity part)
        // 3.   || nu^(-1/2) (p_b - p_e)*n ||_H-1/2(Gamma)         =  standard H-1/2 Sobolev norm for normal flux (pressure part)

        double u_err_squared      = 0.0;
        double flux_u_err_squared = 0.0;
        double flux_p_err_squared = 0.0;

        // evaluate squared errors at gaussian point
        for (int isd=0;isd<my::nsd_;isd++)
        {
          u_err_squared += u_err(isd)*u_err(isd)*surf_fac;
          flux_u_err_squared += flux_u_err(isd)*flux_u_err(isd)*surf_fac;
          flux_p_err_squared += flux_p_err(isd)*flux_p_err(isd)*surf_fac;
        }

        // interface errors
        ele_interf_norms[0] += 1.0/h_k * my::visc_ * u_err_squared;
        ele_interf_norms[1] += h_k     * my::visc_ * flux_u_err_squared;
        ele_interf_norms[2] += h_k     / my::visc_ * flux_p_err_squared;
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
    const std::vector<DRT::UTILS::GaussIntegration> &                 intpoints,                ///< element gauss points
    DRT::Discretization &                                             cutdis,                   ///< boundary discretization
    const std::map<int, std::vector<GEO::CUT::BoundaryCell*> >&       bcells,                   ///< boundary cells
    const std::map<int, std::vector<DRT::UTILS::GaussIntegration> > & bintpoints,               ///< boundary integration points
    std::map<int, std::vector<Epetra_SerialDenseMatrix> > &           side_coupling,            ///< side coupling matrices
    Teuchos::ParameterList&                                           params,                   ///< parameter list
    Epetra_SerialDenseMatrix&                                         elemat1_epetra,           ///< local system matrix of intersected element
    Epetra_SerialDenseVector&                                         elevec1_epetra,           ///< local element vector of intersected element
    Epetra_SerialDenseMatrix&                                         Cuiui,                    ///< coupling matrix of a side with itself
    const GEO::CUT::plain_volumecell_set &                            vcSet,                    ///< set of plain volume cells
    bool                                                              fluidfluidcoupling        ///< indicates fluid-fluid coupling context
)
{
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

  // flag, that indicates whether we have side coupling terms Cuiu, Cuui,...
  // generally, this is the case for fluid-fluid coupling and monolithic XFSI
  // we distinguish by examination of side_coupling
  const bool eval_side_coupling = !side_coupling.empty();

  // plausibility check
  if (!eval_side_coupling && fluidfluidcoupling)
    dserror("You have to evaluate side coupling terms in fluid-fluid coupling");

  // REMARK: to avoid confusion -
  // 'side' = boundary element, part of cutdis (can be warped)
  // 'boundary cell' = belongs to volume-cell

  // do we need convective stabilization?
  bool add_conv_stab = false;

  if (fluidfluidcoupling)
  {
    add_conv_stab = ( fldparaxfem_->XffConvStabScaling() != INPAR::XFEM::XFF_ConvStabScaling_none );
  }
  else
    add_conv_stab = ( fldparaxfem_->ConvStabScaling() != INPAR::XFEM::ConvStabScaling_none );

  //----------------------------------------------------------------------------
  //                         ELEMENT GEOMETRY
  //----------------------------------------------------------------------------

  // get nodal coordinates
  GEO::fillInitialPositionArray<distype, my::nsd_, LINALG::Matrix<my::nsd_, my::nen_> > (ele,my::xyze_);

  ///< element coordinates in EpetraMatrix
  Epetra_SerialDenseMatrix ele_xyze(my::nsd_,my::nen_);
  for ( int i=0; i<my::nen_; ++i )
  {
    for(int j=0; j<my::nsd_; j++)
      ele_xyze(j,i) = my::xyze_( j, i );
  }

  // extract the current velocity & pressure values for current element
  LINALG::Matrix<my::nsd_,my::nen_>  evelaf(true);
  LINALG::Matrix<my::nen_,1>         epreaf(true);
  my::ExtractValuesFromGlobalVector(dis, lm, *my::rotsymmpbc_, &evelaf, &epreaf, "velaf");

  // reconstruct interface force vector
  const Teuchos::RCP<Epetra_Vector> iforcecol = params.get<Teuchos::RCP<Epetra_Vector> >("iforcenp", Teuchos::null);
  bool assemble_iforce = false;
  if (iforcecol != Teuchos::null) assemble_iforce = true;

  // compute characteristic element length based on the background element
  const double h_k = ComputeCharEleLength(ele,ele_xyze,vcSet,bcells,bintpoints);

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

  // L2-projection between stress fields on whole or on physical element volume?
  INPAR::XFEM::Hybrid_LM_L2_Proj hybrid_lm_l2_proj = fldparaxfem_->HybridLM_L2Proj();

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
    double mhvs_param = fldparaxfem_->ViscStabGamma() / 2.0;
    if ( fabs(mhvs_param) < 1.e-8 )
      dserror("MHVS stabilizing parameter n appears in denominator. Please avoid choosing 0.");
  }

  // build volumetric coupling matrices
  HybridLM_Build_VolBased(
      hybrid_lm_l2_proj,
      intpoints, vcSet,
      evelaf, epreaf,
      bK_ss, invbK_ss, K_su, rhs_s, K_us, K_uu, rhs_uu,
      is_MHVS, mhvs_param);

  /*--------------------------------------------------------
    build surface coupling terms
  --------------------------------------------------------*/

  DRT::Element::LocationArray cutla( 1 );

  LINALG::Matrix<3,1> normal(true);
  LINALG::Matrix<3,1> x_side(true);

  // side coupling implementation between background element and each cut side (std::map<sid, side_impl>)
  // for the coupling terms obtained by hybrid LM approach
  std::map<int, Teuchos::RCP<DRT::ELEMENTS::XFLUID::HybridLMInterface<distype> > > side_impl;
  Teuchos::RCP<DRT::ELEMENTS::XFLUID::HybridLMInterface<distype> > si;

  // auxiliary coupling implementation for terms
  // find all the intersecting elements of actele
  std::set<int> begids;
  for (std::map<int,  std::vector<GEO::CUT::BoundaryCell*> >::const_iterator bc=bcells.begin();
       bc!=bcells.end(); ++bc )
  {
    int sid = bc->first;
    begids.insert(sid);
  }

  // map of boundary element gids and coupling matrices, [0]: Gsui, [1]: Guis
  std::map<int, std::vector<Epetra_SerialDenseMatrix> > Cuiui_coupling;

  // map of Nitsche-interfaces for building convective stabilization terms
  std::map<int,Teuchos::RCP<DRT::ELEMENTS::XFLUID::NitscheInterface<distype> > > si_nit;

  // map of boundary element gids and coupling contributions from convective stabilization terms
  std::map<int, std::vector<Epetra_SerialDenseMatrix> > conv_side_coupling;

  // lm vector of all intersecting boundary elements that intersect the current background element
  std::vector<int> patchelementslm;
  std::vector<int> patchelementslmowner;

  // create location vectors for intersecting boundary elements and reshape coupling matrices
  PatchLocationVector(begids,cutdis,patchelementslm,patchelementslmowner, Cuiui_coupling);

  // reshape coupling matrices for convective stabilization terms
  if (add_conv_stab && eval_side_coupling)
  {
    HybridLM_CreateConvStabMatrices(begids,cutdis,conv_side_coupling);
  }
  // evaluate shape function derivatives
  bool eval_deriv = false;

  if (!fluidfluidcoupling) eval_deriv = true; // evaluate derivatives to evaluate traction vector

  //--------------------------------------------
  // loop intersecting sides
  //--------------------------------------------
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

    // side and location vector
    DRT::Element * side = cutdis.gElement( sid );
    side->LocationVector(cutdis,cutla,false);

    // side geometry
    const int numnodes = side->NumNode();
    DRT::Node ** nodes = side->Nodes();
    Epetra_SerialDenseMatrix side_xyze( 3, numnodes );
    for ( int i=0; i<numnodes; ++i )
    {
      const double * x = nodes[i]->X();
      std::copy( x, x+3, &side_xyze( 0, i ) );
    }

    // create side impl
    if (eval_side_coupling)
    {
      std::map<int,std::vector<Epetra_SerialDenseMatrix> >::iterator c = side_coupling.find( sid );
      std::vector<Epetra_SerialDenseMatrix> & side_matrices = c->second;

      if (side_matrices.size() != 3)
        dserror("Obtained only %d side coupling matrices. 3 required.", side_matrices.size());

      // coupling matrices between background element and one! side
      Epetra_SerialDenseMatrix & C_uiu  = side_matrices[0];
      Epetra_SerialDenseMatrix & C_uui  = side_matrices[1];
      Epetra_SerialDenseMatrix & rhC_ui = side_matrices[2];

      // coupling matrices between one side and itself via the element Kss
      std::map<int,std::vector<Epetra_SerialDenseMatrix> >::iterator c2 = Cuiui_coupling.find( sid );
      std::vector<Epetra_SerialDenseMatrix> & Cuiui_matrices = c2->second;
      Epetra_SerialDenseMatrix & eleGsui = Cuiui_matrices[0];
      Epetra_SerialDenseMatrix & eleGuis = Cuiui_matrices[1];

      si = DRT::ELEMENTS::XFLUID::HybridLMInterface<distype>::CreateHybridLMCoupling_XFluidSided(side,side_xyze,C_uiu,C_uui,rhC_ui,eleGsui,eleGuis,fldparaxfem_->IsViscousAdjointSymmetric());
    }
    else
    {
      si = DRT::ELEMENTS::XFLUID::HybridLMInterface<distype>::CreateHybridLMCoupling_XFluidWDBC(side,side_xyze);
    }

    // get velocity at integration point of boundary dis
    si->SetSlaveVel(cutdis,"ivelnp",cutla[0].lm_);
    // set displacement of side
    si->AddSlaveEleDisp(cutdis,"idispnp",cutla[0].lm_);

    // store
    side_impl[sid] = si;

    // define interface force vector w.r.t side
    Epetra_SerialDenseVector iforce;
    iforce.Size(cutla[0].lm_.size());

    if (add_conv_stab)
    {
      if (fluidfluidcoupling)
      {
        // coupling matrices between background element and one! side
        std::map<int,std::vector<Epetra_SerialDenseMatrix> >::iterator c = conv_side_coupling.find( sid );
        std::vector<Epetra_SerialDenseMatrix> & conv_side_matrices = c->second;
        Epetra_SerialDenseMatrix & C_uiu  = conv_side_matrices[0];
        Epetra_SerialDenseMatrix & C_uui  = conv_side_matrices[1];
        Epetra_SerialDenseMatrix & rhC_ui = conv_side_matrices[2];
        Epetra_SerialDenseMatrix & C_uiui = conv_side_matrices[3];

         si_nit[sid] = DRT::ELEMENTS::XFLUID::NitscheInterface<distype>::CreateNitscheCoupling_XFluidSided(
            side, side_xyze, C_uiu, C_uui, rhC_ui, C_uiui, true);
      }
      else if (! eval_side_coupling)
      {
        si_nit[sid] = DRT::ELEMENTS::XFLUID::NitscheInterface<distype>::CreateNitscheCoupling_XFluidWDBC(
                      side, side_xyze);
      }
      else
        dserror("Convective stabilization not yet available for monolithic XFSI!");
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


#ifdef BOUNDARYCELL_TRANSFORMATION_OLD

        si->Evaluate(eta,x_side,normal,drs);

        // find element local position of gauss point at interface
        GEO::CUT::Position<distype> pos( my::xyze_, x_side );
        pos.Compute();
        rst = pos.LocalCoordinates();

#else
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
          for (int idim=0;idim<3;idim++)
          {
            x_gp_lin(idim,0) = gpcord[idim];
          }
        }

        // find element local position of gauss point
        GEO::CUT::Position<distype> pos( my::xyze_, x_gp_lin );
        pos.Compute();
        rst = pos.LocalCoordinates();

        // project gaussian point from linearized interface to warped side (get/set local side coordinates in SideImpl)
        LINALG::Matrix<2,1> xi_side(true);
        si->ProjectOnSide(x_gp_lin, x_side, xi_side);
 #endif

        const double surf_fac = drs*iquad.Weight();

        const double fac = surf_fac * my::fldparatimint_->TimeFac();

        //--------------------------------------------

        // evaluate shape functions (and derivatives)

        if (eval_deriv)
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


        //--------------------------------------------

        bK_ss.MultiplyNT( my::funct_, my::funct_ );

               /*                      \
            - |  (virt tau) * n^f , Du  |
               \                      */

        HybridLM_Evaluate_SurfBased(
          si,
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
          fac,
          eval_side_coupling,
          is_MHVS);

        //--------------------------------------------
        // evaluate additional inflow/convective stabilization terms
        if (add_conv_stab)
        {
          double NIT_full_stab_fac = 0.0;
          double avg_conv_stab_fac = 0.0;
          double velgrad_interface_fac = 0.0;
          double gamma_ghost_penalty = 0.0;
          double presscoupling_interface_fac = 0.0;
          double gamma_press_coupling = 0.0;
          const double NIT_visc_stab_fac = 0.0;

          INPAR::XFEM::XFF_ConvStabScaling xff_conv_stab_scaling = fldparaxfem_->XffConvStabScaling();

          if (fluidfluidcoupling)
          {
            NIT_ComputeFullStabfacFluidFluid(
                NIT_full_stab_fac,             ///< to be filled: full Nitsche's penalty term scaling (viscous+convective part)
                avg_conv_stab_fac,             ///< scaling of the convective average coupling term for fluidfluid problems
                velgrad_interface_fac,
                presscoupling_interface_fac,
                xff_conv_stab_scaling,
                my::velint_.Dot(normal),
                my::velint_,
                h_k,
                NIT_visc_stab_fac, ///< Nitsche's viscous scaling part of penalty term
                gamma_ghost_penalty,
                gamma_press_coupling
                );

             si_nit.at(sid)->ApplyConvStabTerms(
                 si,
                 elemat1_epetra,
                 elevec1_epetra,
                 my::funct_,
                 my::velint_,
                 NIT_full_stab_fac,
                 avg_conv_stab_fac,
                 fac,
                 xff_conv_stab_scaling);
          }
          else if (! eval_side_coupling)
          {
            INPAR::XFEM::ConvStabScaling conv_stab_scaling = fldparaxfem_->ConvStabScaling();
            NIT_ComputeFullStabfac(
                NIT_full_stab_fac, ///< to be filled: full Nitsche's penalty term scaling (viscous+convective part)
                NIT_visc_stab_fac, ///< Nitsche's viscous scaling part of penalty term
                conv_stab_scaling, ///< type of convective stabilization for xfluid
                my::velint_.Dot(normal)
            );

            si_nit.at(sid)->ApplyConvStabTerms(
                si,
                elemat1_epetra,
                elevec1_epetra,
                my::funct_,
                my::velint_,
                NIT_full_stab_fac,
                avg_conv_stab_fac,
                fac);
          }
        }

        //--------------------------------------------
        // calculate interface forces for XFSI
        if (assemble_iforce)
        {
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

          si->ComputeInterfaceForce(iforce, traction, surf_fac );
        } // buildInterfaceForce


      } // end loop gauss points of boundary cell
    } // end loop boundary cells of side

    if(assemble_iforce) AssembleInterfaceForce(iforcecol, cutdis, cutla[0].lm_, iforce);

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

  // in case of WDBC, we are done here
  if (! eval_side_coupling)
    return;

  //-------------------------------------------------
  // finalize creation of side coupling terms
  // Cuiu, Cuui, rhCui & Gsui, Guis (-->Cuiui)
  //-------------------------------------------------

  // build interface coupling matrices - therefore iterate through the interface elements
  for (typename std::map<int,Teuchos::RCP<DRT::ELEMENTS::XFLUID::HybridLMInterface<distype> > >::iterator sit=side_impl.begin();
      sit!=side_impl.end(); ++sit)
  {
    Teuchos::RCP<DRT::ELEMENTS::XFLUID::HybridLMInterface<distype> > si = sit->second;

    // creation of Cuiu,Cuui,rhCui,Guis and Gsui:
    si->HybridLM_buildFinalCouplingMatrices(invK_ss,KusinvKss,K_su,rhs_s);

    // add contributions from convective stabilization, if active
    if (add_conv_stab)
    {
      const int sid = sit->first;
      std::map<int,std::vector<Epetra_SerialDenseMatrix> >::iterator c = side_coupling.find( sid );
      std::vector<Epetra_SerialDenseMatrix> & side_matrices = c->second;
      if (side_matrices.size() != 3)
        dserror("Obtained only %d side coupling matrices. 3 required.", side_matrices.size());

      std::map<int,std::vector<Epetra_SerialDenseMatrix> >::iterator cc = conv_side_coupling.find( sid );
      std::vector<Epetra_SerialDenseMatrix> & conv_side_matrices = cc->second;
      if (conv_side_matrices.size() != 4)
        dserror("Obtained only %d conv. side coupling matrices. 4 required.", conv_side_matrices.size());

      for (int i=0; i<3; ++i)
      {
        if (side_matrices[i].M() != conv_side_matrices[i].M() ||
            side_matrices[i].N() != conv_side_matrices[i].N())
          dserror("Mismatch in matrix dimensions of convective stabilization matrix and MHCS/MHVS coupling matrix");
        side_matrices[i] += conv_side_matrices[i];
      }
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
      const int sid = m->first;
      std::map<int,std::vector<Epetra_SerialDenseMatrix> >::const_iterator c = conv_side_coupling.find( sid );
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
    const INPAR::XFEM::Hybrid_LM_L2_Proj                                                        hybrid_lm_l2_proj,
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

  INPAR::CUT::VCellGaussPts vcellgausspts = fldparaxfem_->VolumeCellGaussPoints();

  // full L2-projection means integration over the full background element,
  // not only the physical part
  if (hybrid_lm_l2_proj == INPAR::XFEM::Hybrid_LM_L2_Proj_full)
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
  else if (vcellgausspts == INPAR::CUT::VCellGaussPts_Tessellation or vcellgausspts == INPAR::CUT::VCellGaussPts_MomentFitting)
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
  else if (vcellgausspts == INPAR::CUT::VCellGaussPts_DirectDivergence )  // DirectDivergence method
  {
    for( std::vector<DRT::UTILS::GaussIntegration>::const_iterator i=intpoints.begin();i!=intpoints.end();++i )
    {
      const DRT::UTILS::GaussIntegration intcell = *i;
      GEO::CUT::VolumeCell * vc = cells[i-intpoints.begin()];

      //----------------------------------------------------------------------
      //integration over the main gauss points to get the required integral
      //----------------------------------------------------------------------
      int mainPtno = 0;
      for ( DRT::UTILS::GaussIntegration::iterator iquad=intcell.begin(); iquad!=intcell.end(); ++iquad )
      {
        // to store values temporarily
        LINALG::BlockMatrix<LINALG::Matrix<my::nen_,my::nen_>,numstressdof_,my::numdofpernode_>     K_suTemp;
        LINALG::BlockMatrix<LINALG::Matrix<my::nen_,my::nen_>,my::numdofpernode_,numstressdof_>     K_usTemp;
        LINALG::BlockMatrix<LINALG::Matrix<my::nen_,1>,numstressdof_,1>                             rhs_sTemp;
        LINALG::BlockMatrix<LINALG::Matrix<my::nen_,my::nen_>,my::nsd_,my::nsd_>                    K_uuTemp;
        LINALG::BlockMatrix<LINALG::Matrix<my::nen_,1>,my::nsd_, 1>                                 rhs_uuTemp;
        LINALG::Matrix<my::nen_,my::nen_>                                                           invbK_ssTemp;

        // get internal Gaussian rule for every main Gauss point
        DRT::UTILS::GaussIntegration gint = vc->GetInternalRule(mainPtno);
        mainPtno++;

        //----------------------------------------------------------------------
        //integration over the internal gauss points - to get modified integrand
        //----------------------------------------------------------------------
        for ( DRT::UTILS::GaussIntegration::iterator quadint=gint.begin(); quadint!=gint.end(); ++quadint )
        {
          my::EvalShapeFuncAndDerivsAtIntPoint( quadint.Point(), quadint.Weight() );
          if (is_MHVS)
          {
            MHVS_Evaluate_VolBased(
                evelaf, bK_ss, invbK_ssTemp, K_suTemp, rhs_sTemp, K_usTemp, K_uuTemp, rhs_uuTemp, mhvs_param);
          }
          else
          {
            MHCS_Evaluate_VolBased(
                evelaf, epreaf, bK_ss, invbK_ssTemp, K_suTemp, rhs_sTemp);
          }
        }

        // Update
        invbK_ss.Update( iquad.Weight(), invbK_ssTemp, 1.0);

        for (int ivel = 0; ivel < my::nsd_; ++ ivel)
        {
          for (int jvel = 0; jvel < my::nsd_; ++jvel)
          {
            K_su(stressIndex(ivel,jvel),ivel)->Update( iquad.Weight(), *K_suTemp(stressIndex(ivel,jvel),ivel), 1.0);
          }
        }

        for (int isigma = 0; isigma < my::numderiv2_; ++ isigma)
          rhs_s(isigma,0)->Update( iquad.Weight(), *rhs_sTemp(isigma,0), 1.0);

        // in case of a Cauchy stress-based approach, fill the pressure column block of K_su
        // and we're done
        if (! is_MHVS)
        {
          K_su( Sigmaxx, Pres )->Update( iquad.Weight(), *K_suTemp( Sigmaxx, Pres ), 1.0 );
          K_su( Sigmayy, Pres )->Update( iquad.Weight(), *K_suTemp( Sigmayy, Pres ), 1.0 );
          K_su( Sigmazz, Pres )->Update( iquad.Weight(), *K_suTemp( Sigmazz, Pres ), 1.0 );

          continue;
        }

        for (int ivel = 0; ivel < my::nsd_; ++ ivel)
        {
          for (int jvel = 0; jvel < my::nsd_; ++jvel)
          {
            K_us(ivel, stressIndex(ivel,jvel))->Update( iquad.Weight(), *K_usTemp(ivel, stressIndex(ivel,jvel)), 1.0);

            K_uu(ivel, jvel)->Update( iquad.Weight(), *K_uuTemp(ivel,jvel), 1.0);
          }
          rhs_uu(ivel,0)->Update( iquad.Weight(), *rhs_uuTemp(ivel,0),1.0);
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
    const LINALG::Matrix<3,1> &                                                               normal,
    const double &                                                                            timesurffac,
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
      si->MHVS_buildCouplingMatrices(normal, timesurffac, my::funct_, rhs_s, press, rhs_pu);
    else
      si->MHCS_buildCouplingMatrices(normal, timesurffac, my::funct_, rhs_s);
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
    LINALG::Matrix<my::nsd_,1> velint_WDBC(true);
    si->GetInterfaceVel(velint_WDBC);

    // add surface integral contribution to rhs_s

    // from diagonal terms
    rhs_s(Sigmaxx, 0)->Update(-timesurffac * normal(Velx) * velint_WDBC(Velx), my::funct_, 1.0);
    rhs_s(Sigmayy, 0)->Update(-timesurffac * normal(Vely) * velint_WDBC(Vely), my::funct_, 1.0);
    rhs_s(Sigmazz, 0)->Update(-timesurffac * normal(Velz) * velint_WDBC(Velz), my::funct_, 1.0);

    // from off-diagonal terms
    rhs_s(Sigmaxy, 0)->Update(-timesurffac * (normal(Vely) * velint_WDBC(Velx) + normal(Velx) * velint_WDBC(Vely)), my::funct_, 1.0);
    rhs_s(Sigmaxz, 0)->Update(-timesurffac * (normal(Velz) * velint_WDBC(Velx) + normal(Velx) * velint_WDBC(Velz)), my::funct_, 1.0);
    rhs_s(Sigmayz, 0)->Update(-timesurffac * (normal(Velz) * velint_WDBC(Vely) + normal(Vely) * velint_WDBC(Velz)), my::funct_, 1.0);

    if (! is_MHVS) return;

    // ONLY MHVS:
    // pressure-tested kinematic continuity term
    /*
     *
     *      /      _  \
     *    +|  q n, u   |
     *      \         /
     */
    const double normalvel = velint_WDBC.Dot(normal);

    rhs_pu.Update(-timesurffac * normalvel, my::funct_, 1.0);
  }
}

/*--------------------------------------------------------------------------------
 * add Nitsche (NIT) interface condition to element matrix and rhs
 *--------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void FluidEleCalcXFEM<distype>::ElementXfemInterfaceNIT(
    DRT::ELEMENTS::Fluid *                                              ele,               ///< fluid element
    DRT::Discretization &                                               dis,               ///< background discretization
    const std::vector<int> &                                            lm,                ///< element local map
    DRT::Discretization &                                               cutdis,            ///< cut discretization
    const std::map<int, std::vector<GEO::CUT::BoundaryCell*> > &        bcells,            ///< boundary cells
    const std::map<int, std::vector<DRT::UTILS::GaussIntegration> > &   bintpoints,        ///< boundary integration points
    std::map<int, std::vector<Epetra_SerialDenseMatrix> > &             side_coupling,     ///< side coupling matrices
    Teuchos::ParameterList&                                             params,            ///< parameter list
    Epetra_SerialDenseMatrix&                                           elemat1_epetra,    ///< element matrix
    Epetra_SerialDenseVector&                                           elevec1_epetra,    ///< element vector
    Epetra_SerialDenseMatrix&                                           Cuiui,             ///< ui-ui coupling matrix
    const GEO::CUT::plain_volumecell_set&                               vcSet,             ///< volumecell sets in this element
    bool                                                                fluidfluidcoupling ///< Is this xfluidfluid problem?
  )
{
  const Teuchos::RCP<Epetra_Vector> iforcecol = params.get<Teuchos::RCP<Epetra_Vector> >("iforcenp", Teuchos::null);

  bool assemble_iforce = false;
  if(iforcecol != Teuchos::null) assemble_iforce = true;

  //----------------------------------------------------------------------------
  //                         ELEMENT GEOMETRY
  //----------------------------------------------------------------------------

  // get node coordinates
  GEO::fillInitialPositionArray< distype, my::nsd_, LINALG::Matrix<my::nsd_,my::nen_> >( ele, my::xyze_ );

  ///< element coordinates in EpetraMatrix
  Epetra_SerialDenseMatrix ele_xyze(my::nsd_,my::nen_);
  for ( int i=0; i<my::nen_; ++i )
  {
    for (int j=0; j<my::nsd_; j++)
      ele_xyze(j,i) = my::xyze_( j, i );
  }

  // get element-wise velocity/pressure field
  LINALG::Matrix<my::nsd_,my::nen_> evelaf(true);
  LINALG::Matrix<my::nen_,1> epreaf(true);
  my::ExtractValuesFromGlobalVector(dis, lm, *my::rotsymmpbc_, &evelaf, &epreaf, "velaf");

  // flag, that indicates whether we have side coupling terms Cuiu, Cuui,...
  // generally, this is the case for fluid-fluid coupling and monolithic XFSI
  // we distinguish by examination of side_coupling
  bool eval_side_coupling = !side_coupling.empty();

  //----------------------------------------------------------------------------
  //      surface integral --- build Cuiui, Cuui, Cuiu and Cuu matrix and rhs
  //----------------------------------------------------------------------------

  DRT::Element::LocationArray cutla( 1 );

  LINALG::Matrix<3,1> normal(true); // normal vector w.r.t the fluid domain, points from the fluid into the structure
  LINALG::Matrix<3,1> x_side(true);

  // side coupling implementation between background element and each cut side
  std::map<int, Teuchos::RCP<DRT::ELEMENTS::XFLUID::NitscheInterface<distype> > > side_impl;
  Teuchos::RCP<DRT::ELEMENTS::XFLUID::NitscheInterface<distype> > si;

  // map of boundary element gids and coupling matrices, [0]: Cuiui matrix
  std::map<int, std::vector<Epetra_SerialDenseMatrix> > Cuiui_coupling;
  // we don't have Cuiui for standard Dirichlet problems...
  if (eval_side_coupling)
  {
    // find all the intersecting elements of actele
    std::set<int> begids;
    for (std::map<int,  std::vector<GEO::CUT::BoundaryCell*> >::const_iterator bc=bcells.begin();
         bc!=bcells.end(); ++bc )
    {
      int sid = bc->first;
      begids.insert(sid);
    }

    // lm vector of all intersecting boundary elements that intersect the current background element
    std::vector<int> patchelementslmv;
    std::vector<int> patchelementslmowner;

    // create location vectors for intersecting boundary elements and reshape coupling matrices
    PatchLocationVector(begids,cutdis,patchelementslmv,patchelementslmowner, Cuiui_coupling);
  }

  //-----------------------------------------------------------------------------------
  //         evaluate element length, stabilization factors and average weights
  //-----------------------------------------------------------------------------------

  //--------------------------------------------
  // element length based on the background element
  //--------------------------------------------
  const double h_k = ComputeCharEleLength(ele, ele_xyze, vcSet, bcells, bintpoints);

  //--------------------------------------------
  // compute viscous part of Nitsche's penalty term scaling for Nitsche's method (dimensionless)
  //--------------------------------------------
  const double NIT_visc_stab_fac = NIT_ComputeNitscheStabfac( ele->Shape(), h_k, -1 );

  //--------------------------------------------
  // loop intersecting sides
  //--------------------------------------------
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

    // side and location vector
    DRT::Element * side = cutdis.gElement( sid );
    side->LocationVector(cutdis,cutla,false);

    // side geometry
    const int numnodes = side->NumNode();
    DRT::Node ** nodes = side->Nodes();
    Epetra_SerialDenseMatrix side_xyze( 3, numnodes );
    for ( int i=0; i<numnodes; ++i )
    {
      const double * x = nodes[i]->X();
      std::copy( x, x+3, &side_xyze( 0, i ) );
    }

    // create side impl

    if (eval_side_coupling)
    {
      std::map<int,std::vector<Epetra_SerialDenseMatrix> >::iterator c = side_coupling.find( sid );
      std::vector<Epetra_SerialDenseMatrix> & side_matrices = c->second;

      // do consistency check
      if (side_matrices.size() == 3)
      {
        // coupling matrices between background element and one side
        Epetra_SerialDenseMatrix & C_uiu  = side_matrices[0];
        Epetra_SerialDenseMatrix & C_uui  = side_matrices[1];
        Epetra_SerialDenseMatrix & rhC_ui = side_matrices[2];

        // coupling matrices between one side and itself
        std::map<int,std::vector<Epetra_SerialDenseMatrix> >::iterator c2 = Cuiui_coupling.find( sid );
        std::vector<Epetra_SerialDenseMatrix> & Cuiui_matrices = c2->second;
        Epetra_SerialDenseMatrix & eleCuiui = Cuiui_matrices[0];

        si = DRT::ELEMENTS::XFLUID::NitscheInterface<distype>::CreateNitscheCoupling_XFluidSided(
            side,side_xyze,C_uiu,C_uui,rhC_ui,eleCuiui,fldparaxfem_->IsViscousAdjointSymmetric());
      }
      else
      {
        dserror("Got an invalid number of %d side coupling matrices for Nitsche-coupling.", side_matrices.size());
      }
    }
    else
    {
      si = DRT::ELEMENTS::XFLUID::NitscheInterface<distype>::CreateNitscheCoupling_XFluidWDBC(side,side_xyze);
    }

    side_impl[sid] = si;

    // get velocity at integration point of boundary dis
    si->SetSlaveVel(cutdis,"ivelnp",cutla[0].lm_);

    // set displacement of side
    si->AddSlaveEleDisp(cutdis,"idispnp",cutla[0].lm_);

    // define interface force vector w.r.t side
    Epetra_SerialDenseVector iforce;
    iforce.Size(cutla[0].lm_.size());

    //--------------------------------------------
    // loop boundary cells w.r.t current cut side
    //--------------------------------------------
    for ( std::vector<DRT::UTILS::GaussIntegration>::const_iterator i=cutintpoints.begin();
          i!=cutintpoints.end();
          ++i )
    {
      const DRT::UTILS::GaussIntegration & gi = *i;
      GEO::CUT::BoundaryCell * bc = bcs[i - cutintpoints.begin()]; // get the corresponding boundary cell, bc-orientation is outward-pointing from fluid to structure

      //--------------------------------------------
      // loop gausspoints w.r.t current boundary cell
      //--------------------------------------------
      for ( DRT::UTILS::GaussIntegration::iterator iquad=gi.begin(); iquad!=gi.end(); ++iquad )
      {
        double drs = 0.0; // transformation factor between reference cell and linearized boundary cell

        const LINALG::Matrix<2,1> eta( iquad.Point() ); // xi-coordinates with respect to side

        LINALG::Matrix<3,1> rst(true); // local coordinates w.r.t background element

#ifdef BOUNDARYCELL_TRANSFORMATION_OLD

        si->Evaluate(eta,x_side,normal,drs);

        // find element local position of gauss point at interface
        GEO::CUT::Position<distype> pos( my::xyze_, x_side );
        pos.Compute();
        rst = pos.LocalCoordinates();

#else
        LINALG::Matrix<3,1> x_gp_lin(true); // gp in xyz-system on linearized interface

        // compute transformation factor, normal vector and global Gauss point coordiantes
        if (bc->Shape() != DRT::Element::dis_none) // Tessellation approach
        {
          ComputeSurfaceTransformation(drs, x_gp_lin, normal, bc, eta);
        }
        else // MomentFitting approach
        {
          drs = 1.0;
          normal = bc->GetNormalVector();
          const double* gpcord = iquad.Point();
          for (int idim=0;idim<3;idim++)
          {
            x_gp_lin(idim,0) = gpcord[idim];
          }
        }

        // find element local position of gauss point
        GEO::CUT::Position<distype> pos( my::xyze_, x_gp_lin );
        pos.Compute();
        rst = pos.LocalCoordinates();

        // project gaussian point from linearized interface to warped side (get/set local side coordinates in SideImpl)
        LINALG::Matrix<2,1> xi_side(true);
        si->ProjectOnSide(x_gp_lin, x_side, xi_side);
 #endif

        const double surf_fac = drs*iquad.Weight();
        const double timefacfac = surf_fac * my::fldparatimint_->TimeFac();

        //--------------------------------------------

        // evaluate shape functions (and derivatives)
        EvalFuncAndDeriv( rst );

        // get velocity at integration point
        // (values at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
        my::velint_.Multiply(evelaf,my::funct_);

        // get velocity derivatives at integration point
        // (values at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
        my::vderxy_.MultiplyNT(evelaf,my::derxy_);

        // get pressure at integration point
        // (value at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
        double press = my::funct_.Dot(epreaf);

        //--------------------------------------------
        // compute stabilization factors

        double NIT_full_stab_fac = 0.0;
        double avg_conv_stab_fac = 0.0;

        if (fluidfluidcoupling)
        {
          // convective stabilization type for xfluidfluid
          INPAR::XFEM::XFF_ConvStabScaling xff_conv_stab_scaling = fldparaxfem_->XffConvStabScaling();

          double velgrad_interface_fac = 0.0;
          double gamma_ghost_penalty = 0.0;
          double presscoupling_interface_fac = 0.0;
          double gamma_press_coupling = 0.0;

          NIT_ComputeFullStabfacFluidFluid(
              NIT_full_stab_fac,             ///< to be filled: full Nitsche's penalty term scaling (viscous+convective part)
              avg_conv_stab_fac,             ///< to be filled: scaling of the convective average coupling term for fluidfluid problems
              velgrad_interface_fac,         ///< to be filled: stabilization parameter for velocity-gradients penalty term
              presscoupling_interface_fac,   ///< to be filled: stabilization parameter for pressure-coupling penalty term
              xff_conv_stab_scaling,
              my::velint_.Dot(normal),
              my::velint_,
              h_k,
              NIT_visc_stab_fac,
              gamma_ghost_penalty,
              gamma_press_coupling
              );
        }
        else
        {
          // convective stabilization type for xfluid
          INPAR::XFEM::ConvStabScaling conv_stab_scaling = fldparaxfem_->ConvStabScaling();

          NIT_ComputeFullStabfac(
              NIT_full_stab_fac, ///< to be filled: full Nitsche's penalty term scaling (viscous+convective part)
              NIT_visc_stab_fac, ///< Nitsche's viscous scaling part of penalty term
              conv_stab_scaling, ///< type of convective stabilization for xfluid
              my::velint_.Dot(normal)
          );
        }

        //--------------------------------------------
        // define average weights -
        // this is the routine for xfluid-sided
        // Nitsche, so kappa_s is set to 0
        //--------------------------------------------

        const double kappa1 = 1.0;
        const double kappa2 = 0.0;

        si->NIT_buildCouplingMatrices(
          elemat1_epetra,          // standard bg-bg-matrix
          elevec1_epetra,          // standard bg-rhs
          normal,                  // normal vector
          timefacfac,              // theta*dt
          my::visceff_,            // dynamic viscosity in background fluid
          my::visceff_,            // dynamic viscosity in embedded fluid
          kappa1,                  // mortaring weighting
          kappa2,                  // mortaring weighting
          NIT_full_stab_fac,       // full Nitsche's penalty term scaling (viscous+convective part)
          avg_conv_stab_fac,       // scaling of the convective average coupling term for fluidfluid problems
          my::funct_,              // bg shape functions
          my::derxy_,              // bg deriv
          my::vderxy_,             // bg deriv^n
          press,                   // bg p^n
          my::velint_              // bg u^n
        );


        //--------------------------------------------
        // calculate interface forces for XFSI
        if (assemble_iforce)
        {

          // get pressure at integration point
          // (value at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
          double press = my::funct_.Dot(epreaf);

          //-------------------------------
          // traction vector w.r.t fluid domain, resulting stresses acting on the fluid surface
          // t= (-p*I + 2mu*eps(u))*n^f
          LINALG::Matrix<my::nsd_,1> traction(true);

          BuildTractionVector( traction, press, normal );

          si->ComputeInterfaceForce(iforce, traction, surf_fac );

        } // buildInterfaceForce
      } // end loop gauss points of boundary cell
    } // end loop boundary cells of side

    if(assemble_iforce) AssembleInterfaceForce(iforcecol, cutdis, cutla[0].lm_, iforce);

  } // end loop cut sides


  //----------------------------------------------------------------------------
  // build Cuiui coupling matrix (includes patch of Cuiui matrices for all sides)
  //----------------------------------------------------------------------------

  if ( eval_side_coupling )
    NIT_BuildPatchCuiui(Cuiui, Cuiui_coupling);

  return;
}



/*--------------------------------------------------------------------------------
 *--------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void FluidEleCalcXFEM<distype>::ElementXfemInterfaceNIT2(
    DRT::ELEMENTS::Fluid *                                              ele,
    DRT::Discretization &                                               dis,
    const std::vector<int> &                                            lm,
    DRT::Discretization &                                               cutdis,
    const std::map<int, std::vector<GEO::CUT::BoundaryCell*> > &        bcells,
    const std::map<int, std::vector<DRT::UTILS::GaussIntegration> > &   bintpoints,
    std::map<int, std::vector<Epetra_SerialDenseMatrix> > &             side_coupling,
    Teuchos::ParameterList&                                             params,
    DRT::Discretization &                                               alediscret,
    std::map<int,int> &                                                 boundary_emb_gid_map,
    Epetra_SerialDenseMatrix&                                           elemat1_epetra,
    Epetra_SerialDenseVector&                                           elevec1_epetra,
    Epetra_SerialDenseMatrix&                                           Cuiui,
    const GEO::CUT::plain_volumecell_set &                              vcSet
  )
{
  //----------------------------------------------------------------------------
  //                         ELEMENT GEOMETRY
  //----------------------------------------------------------------------------

  // get node coordinates
  GEO::fillInitialPositionArray< distype, my::nsd_, LINALG::Matrix<my::nsd_,my::nen_> >( ele, my::xyze_ );

  ///< element coordinates in EpetraMatrix
  Epetra_SerialDenseMatrix ele_xyze(my::nsd_,my::nen_);
  for ( int i=0; i<my::nen_; ++i )
  {
    for(int j=0; j<my::nsd_; j++)
      ele_xyze(j,i) = my::xyze_( j, i );
  }

  // get element-wise velocity/pressure field
  LINALG::Matrix<my::nsd_,my::nen_> evelaf(true);
  LINALG::Matrix<my::nen_,1> epreaf(true);
  my::ExtractValuesFromGlobalVector(dis, lm, *my::rotsymmpbc_, &evelaf, &epreaf, "velaf");

  //----------------------------------------------------------------------------
  //      surface integral --- build Cuiui, Cuui, Cuiu and Cuu matrix and rhs
  //----------------------------------------------------------------------------

  DRT::Element::LocationArray alela( 1 );
  DRT::Element::LocationArray cutla( 1 );

  LINALG::Matrix<3,1> normal(true);
  LINALG::Matrix<3,1> x_side(true);

  bool fluidfluidcoupling = false;

  // embedded element coupling implementation between background element and each cutting embedded element (std::map<sid, emb_impl)
  // side coupling implementation between background element and each cut side (std::map<sid, side_impl)
  std::map<int, Teuchos::RCP<DRT::ELEMENTS::XFLUID::NitscheInterface<distype> > > emb_impl;
  Teuchos::RCP<DRT::ELEMENTS::XFLUID::NitscheInterface<distype> > emb;
  std::map<int, Teuchos::RCP<DRT::ELEMENTS::XFLUID::SlaveElementInterface<distype> > > side_impl;
  Teuchos::RCP<DRT::ELEMENTS::XFLUID::SlaveElementInterface<distype> > si;

  // find all the intersecting elements of actele
  std::set<int> begids;
  for (std::map<int,  std::vector<GEO::CUT::BoundaryCell*> >::const_iterator bc=bcells.begin();
      bc!=bcells.end(); ++bc )
  {
    int sid = bc->first;
    begids.insert(sid);
  }

  // map of boundary element gids, to coupling matrices Cuiui
  std::map<int, std::vector<Epetra_SerialDenseMatrix> > Cuiui_coupling;

  // lm vector of all intersecting elements (boundary elements which intersect the current element)
  std::vector<int> patchelementslmv;
  std::vector<int> patchelementslmowner;

  // create location vectors for intersecting embedded elements and reshape coupling matrices
  PatchLocationVector(begids,alediscret,patchelementslmv,patchelementslmowner,Cuiui_coupling,boundary_emb_gid_map);

  //-----------------------------------------------------------------------------------
  //         evaluate element length, stabilization factors and average weights
  //-----------------------------------------------------------------------------------

  INPAR::XFEM::CouplingStrategy coupling_strategy = fldparaxfem_->GetCouplingStrategy();

  // initialize the characteristic element length
  double h_k = 0.0;

  // compute characteristic element length based on the background element
  if(coupling_strategy != INPAR::XFEM::Embedded_Sided_Coupling)
    h_k = ComputeCharEleLength(ele, ele_xyze, vcSet, bcells, bintpoints);

  //------------------------------
  // define average weights

  double kappa1 = 0.0;

  if (coupling_strategy == INPAR::XFEM::Embedded_Sided_Coupling)
  {
    kappa1 = 0.0;
  }
  else if (coupling_strategy == INPAR::XFEM::Two_Sided_Coupling)
  {
    if( h_k <= 0.0 ) dserror("element length is <= 0.0");

    kappa1 = 0.5;
  }
  else if (coupling_strategy == INPAR::XFEM::Xfluid_Sided_Coupling)
  {
    dserror("Use ElementXfemInterfaceNIT for purely xfluid-sided coupling.");
  }
  else dserror("coupling strategy not known");

  if ( kappa1 > 1.0 || kappa1 < 0.0) dserror("Nitsche weights for inverse estimate kappa1 lies not in [0,1]: %d", kappa1);

  double kappa2 = 1.0-kappa1;


  //----------------------------------------------------------------------------
  //      surface integral --- loop sides
  //----------------------------------------------------------------------------
  // map of side-element id and Gauss points
  for ( std::map<int, std::vector<DRT::UTILS::GaussIntegration> >::const_iterator i=bintpoints.begin();
      i!=bintpoints.end();
      ++i )
  {
    int sid = i->first;
    const std::vector<DRT::UTILS::GaussIntegration> & cutintpoints = i->second;

    std::map<int, std::vector<GEO::CUT::BoundaryCell*> >::const_iterator j = bcells.find( sid );
    if ( j==bcells.end() )
      dserror( "missing boundary cell" );

    const std::vector<GEO::CUT::BoundaryCell*> & bcs = j->second;
    if ( bcs.size()!=cutintpoints.size() )
      dserror( "boundary cell integration rules mismatch" );

    DRT::Element * side = cutdis.gElement( sid );
    side->LocationVector(cutdis,cutla,false);
    DRT::Element * emb_ele = alediscret.gElement( boundary_emb_gid_map.find(sid)->second );
    emb_ele->LocationVector(alediscret,alela,false);

    const int emb_numnodes = emb_ele->NumNode();
    DRT::Node ** emb_nodes = emb_ele->Nodes();
    Epetra_SerialDenseMatrix emb_xyze( 3, emb_numnodes );
    for ( int i=0; i<emb_numnodes; ++i )
    {
      const double * x = emb_nodes[i]->X();
      std::copy( x, x+3, &emb_xyze( 0, i ) );
    }

    const int numnodes = side->NumNode();
    DRT::Node ** nodes = side->Nodes();
    Epetra_SerialDenseMatrix side_xyze( 3, numnodes );
    for ( int i=0; i<numnodes; ++i )
    {
      const double * x = nodes[i]->X();
      std::copy( x, x+3, &side_xyze( 0, i ) );
    }

    std::map<int,std::vector<Epetra_SerialDenseMatrix> >::iterator c = side_coupling.find( sid );

    std::vector<Epetra_SerialDenseMatrix> & side_matrices = c->second;

    // Todo: the number of coupling matrices between element and intersecting sides
    // is not a valid criterion to distinguish between fluidfluid-coupling and monolithic XFSI!
    if ( side_matrices.size()==3 )
      fluidfluidcoupling = true;


    if (fluidfluidcoupling)
    {
      // coupling matrices between background element and one! side
      Epetra_SerialDenseMatrix & C_uiu  = side_matrices[0];
      Epetra_SerialDenseMatrix & C_uui  = side_matrices[1];
      Epetra_SerialDenseMatrix & rhC_ui = side_matrices[2];

      // coupling matrices between one side and itself via the element Kss
      std::map<int,std::vector<Epetra_SerialDenseMatrix> >::iterator c2 = Cuiui_coupling.find( sid );
      std::vector<Epetra_SerialDenseMatrix> & Cuiui_matrices = c2->second;
      Epetra_SerialDenseMatrix & eleCuiui = Cuiui_matrices[0];

      // create interface for the embedded element and the associated side
      emb = DRT::ELEMENTS::XFLUID::NitscheInterface<distype>::CreateNitscheCoupling_TwoSided(
          emb_ele,emb_xyze,C_uiu,C_uui,rhC_ui,eleCuiui,fldparaxfem_->IsViscousAdjointSymmetric());
      si = DRT::ELEMENTS::XFLUID::SlaveElementInterface<distype>::CreateSlaveElementRepresentation(
          side,side_xyze);

    }
    else
    {
      // Todo: It is planned to use this method also for two-phase flow...
      dserror("InterfaceNitscheTwoSided should not be called for non-fluidfluidcoupling!");
    }

    emb_impl[sid] = emb;
    side_impl[sid] = si;

    // get velocity at integration point of boundary dis
    emb->SetSlaveVel(alediscret,"velaf",alela[0].lm_);

    // set displacement of embedded element & side
    emb->AddSlaveEleDisp(alediscret,"dispnp",alela[0].lm_);
    si->AddSlaveEleDisp(cutdis,"idispnp",cutla[0].lm_);

    //--------------------------------------------
    // compute viscous part of Nitsche's penalty term scaling for Nitsche's method (dimensionless)
    //--------------------------------------------

    double NIT_visc_stab_fac = 0.0;

    // set the embedded element length dependent on side in case of Emb Embedded_Sided_Coupling
    if (coupling_strategy == INPAR::XFEM::Embedded_Sided_Coupling)
    {
      // compute characteristic element length based on the embedded element
      h_k = ComputeCharEleLength(emb_ele, emb_xyze, vcSet, bcells, bintpoints, emb, side);

      // compute Nitsche stabilization factor based on the embedded element
      NIT_visc_stab_fac = NIT_ComputeNitscheStabfac( emb_ele->Shape(), h_k, sid );
    }
    else
    {
      // compute Nitsche stabilization factor based on the background element
      NIT_visc_stab_fac = NIT_ComputeNitscheStabfac( ele->Shape(), h_k, -1 );
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

#ifdef BOUNDARYCELL_TRANSFORMATION_OLD

        si->Evaluate(eta,x_side,normal,drs);

        // find element local position of gauss point at interface
        GEO::CUT::Position<distype> pos( my::xyze_, x_side );
        pos.Compute();
        rst = pos.LocalCoordinates();

#else
        LINALG::Matrix<3,1> x_gp_lin(true); // gp in xyz-system on linearized interface

        // compute transformation factor, normal vector and global Gauss point coordiantes
        if (bc->Shape() != DRT::Element::dis_none) // Tessellation approach
        {
          ComputeSurfaceTransformation(drs, x_gp_lin, normal, bc, eta);
        }
        else // MomentFitting approach
        {
          drs = 1.0;
          normal = bc->GetNormalVector();
          const double* gpcord = iquad.Point();
          for (int idim=0;idim<3;idim++)
          {
            x_gp_lin(idim,0) = gpcord[idim];
          }
        }

        // find element local position of gauss point
        GEO::CUT::Position<distype> pos( my::xyze_, x_gp_lin );
        pos.Compute();
        rst = pos.LocalCoordinates();

        // project gaussian point from linearized interface to warped side (get/set local side coordinates in SideImpl)
        LINALG::Matrix<2,1> xi_side(true);
        si->ProjectOnSide(x_gp_lin, x_side, xi_side);
 #endif

        const double surf_fac = drs*iquad.Weight();

        const double timefacfac = surf_fac * my::fldparatimint_->TimeFac();

        //--------------------------------------------
        // evaluate shape functions (and derivatives)

        // evaluate embedded element shape functions
        emb->Evaluate( x_side );

        // evaluate background element shape functions
        EvalFuncAndDeriv( rst );

        // get velocity at integration point
        // (values at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
        my::velint_.Multiply(evelaf,my::funct_);

        // get velocity derivatives at integration point
        // (values at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
        my::vderxy_.MultiplyNT(evelaf,my::derxy_);

        // get pressure at integration point
        // (value at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
        double press = my::funct_.Dot(epreaf);

        //--------------------------------------------
        // compute stabilization factors

        double NIT_full_stab_fac = 0.0;
        double avg_conv_stab_fac = 0.0;

        double gamma_velgrad = 0.0;
        bool velgrad_interface_stab = fldparaxfem_->IsVelGradInterfaceStab();
        if (velgrad_interface_stab)
          gamma_velgrad = fldparaxfem_->GetVelGradInterfaceStabFac();

        double gamma_press_coupling = 0.0;
        bool presscoupling_interface_stab = fldparaxfem_->IsPressureGradientInterfaceStab();
        if (presscoupling_interface_stab)
          gamma_press_coupling = fldparaxfem_->GetPressureGradientInterfaceStabFac();

        if (fluidfluidcoupling)
        {
          // convective stabilization type for xfluidfluid
          INPAR::XFEM::XFF_ConvStabScaling xff_conv_stab_scaling = fldparaxfem_->XffConvStabScaling();

          NIT_ComputeFullStabfacFluidFluid(
              NIT_full_stab_fac,             ///< to be filled: full Nitsche's penalty term scaling (viscous+convective part)
              avg_conv_stab_fac,             ///< to be filled: scaling of the convective average coupling term for fluidfluid problems
              gamma_velgrad,                 ///< to be filled: stabilization parameter for velocity-gradients penalty term
              gamma_press_coupling,          ///< to be filled: stabilization parameter for pressure-coupling penalty term
              xff_conv_stab_scaling,
              my::velint_.Dot(normal),
              my::velint_,
              h_k,
              NIT_visc_stab_fac,             ///< Nitsche's viscous scaling part of penalty term
              gamma_velgrad,
              gamma_press_coupling
              );

          //zero velocity jump for fluidfluidcoupling
          emb->NIT_buildCouplingMatrices(
              elemat1_epetra,              // standard bg-bg-matrix
              elevec1_epetra,              // standard bg-rhs
              normal,                      // normal vector
              timefacfac,                  // theta*dt
              my::visceff_,                // dynvisc viscosity in background fluid
              my::visceff_,                // dynvisc viscosity in embedded fluid
              kappa1,                      // mortaring weighting
              kappa2,                      // mortaring weighting
              NIT_full_stab_fac,           // full Nitsche's penalty term scaling (viscous+convective part)
              avg_conv_stab_fac,           // scaling of the convective average coupling term for fluidfluid problems
              my::funct_,                  // bg shape functions
              my::derxy_,                  // bg deriv
              my::vderxy_,                 // bg deriv^n
              press,                       // bg p^n
              my::velint_,                 // bg u^n
              xff_conv_stab_scaling        // Inflow term strategies for xfluidfluid
            );


          // if desired, evaluate additional interface terms, that penalize velocity gradient jumps
          if (velgrad_interface_stab || presscoupling_interface_stab)
            emb->NIT_applyAdditionalInterfacePenaltyTerms(
                elemat1_epetra,              // standard bg-bg-matrix
                elevec1_epetra,              // standard bg-rhs
                normal,                      // normal vector
                timefacfac,                  // theta*dt
                my::funct_,                  // bg shape functions
                my::derxy_,                  // bg deriv
                my::vderxy_,                 // bg deriv^n
                press,                       // bg p^n
                velgrad_interface_stab,      // penalty term for velocity gradients at the interface
                gamma_velgrad,               // stabilization fac for velocity gradients at the interface
                presscoupling_interface_stab,// penalty term for pressure coupling at the interface // tODo: use gamma here!
                gamma_press_coupling         // stabilization fac for pressure coupling at the interface
            );


        }
        else  dserror(" no two sided mortaring for non-fluidfluidcoupling");

      } // end loop gauss points of boundary cell
    } // end loop boundary cells of side
  } // end loop cut sides


  //----------------------------------------------------------------------------
  // build Cuiui coupling matrix (includes patch of Cuiui matrices for all sides)
  //----------------------------------------------------------------------------

  if ( fluidfluidcoupling )
    NIT_BuildPatchCuiui(Cuiui, Cuiui_coupling);

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
  // build patch-Cuiui matrix
  int ipatchsizesbefore = 0;
  for (std::map<int, std::vector<Epetra_SerialDenseMatrix> >::const_iterator m=Cuiui_coupling.begin();
      m!=Cuiui_coupling.end(); ++m)
  {

    int bid = m->first;
    std::vector<Epetra_SerialDenseMatrix> & Cuiui_mats = Cuiui_coupling[bid];

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
 * compute viscous part of Nitsche's penalty term scaling for Nitsche's method
 *--------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
double FluidEleCalcXFEM<distype>::NIT_ComputeNitscheStabfac(
    const DRT::Element::DiscretizationType ele_distype, ///< the discretization type of the element w.r.t which the stabilization factor is computed
    const double h_k,                                   ///< the characteristic element length
    const int sid                                       ///< Id of the side element for embedded sided fluidfluid formulations
)
{


  //------------------------------
  // scaling factors for Nitsche's standard stabilization term
  /*
        viscous_Nitsche-part
      /                                        \        /                           i        \
     |  NIT_visc_stab_fac *  [ v ] , [ Du ]     | =  - |   NIT_visc_stab_fac * [ v ], [ u ]   |
      \                                        /        \                                    /


      with

            NIT_visc_stab_fac =   gamma * mu * C^2

      where
           - gamma is a dimensionless user-defined parameter
           - mu is the viscosity
           - C is an estimate of the following trace inequality with  C^2 ~ 1/h
             and C depends on element type, shape and polynomial degree


                    trace inequality to be satisfied for Nitsche's method

                    || grad v_h ||         <=   C  *  || grad v_h ||
                                  Gamma_K                           K

           - C^2 \approx lambda can be estimated as the maximal eigenvalue lambda of the generalized eigenvalue problem (A*x = lambda * B*x)

                    < grad w_h,  grad v_h >          =  lambda *  ( grad w_h, grad v_h)
                                           Gamma_K                                     K

           - alternatively we can estimate C^2 as C^2 \approx \rho_safety * CT / h_K (see generalized hp-framework by Erik Burman 2006)

                      with CT depending on the dimension d and the polynomial degree p

                               CT = p(p+1)/2 * (2+1/p)^d ( for tensor-product elements (hexahedral,quadrilateral...)
                               CT = (p+1)*(p+d)/d        ( for simplices (lines, triangles, tetrahedral...)

                      and  1/h_K \approx meas(Gamma_K)/meas(K)
                      and rho_safety a safety factor depending on the element-type which is useful when the surface cuts the element in an inclined situation or when the surface
                          contains small acute angles or corners
  */


  // get the type how to estimate the scaling from the trace inequality
  INPAR::XFEM::ViscStab_TraceEstimate visc_stab_trace_estimate = fldparaxfem_->ViscStabTracEstimate();

  //------------------------------
  // to be filled viscous part of Nitsche's penalty term scaling for Nitsche's method
  double NIT_visc_stab_fac = 0.0;

  // compute the final viscous scaling part of Nitsche's penalty term
  if(visc_stab_trace_estimate == INPAR::XFEM::ViscStab_TraceEstimate_CT_div_by_hk)
  {
    if(fabs(h_k) < 1e-14) dserror("the characteristic element length is zero, it has not been set properly!");

    // get an estimate of the hp-depending constant C_T satisfying the trace inequality w.r.t the corresponding element (compare different weightings)
    double C_T = NIT_getTraceEstimateConstant(ele_distype);

    // get a safety factor dependent on the the element shape for cut elements, since the surface can cut the
    // background element in a bad way when the surface includes sharp corners or small acute angles
    //TODO: still to be set
    double rho_safety = 1.0;

    // build the final viscous scaling
    NIT_visc_stab_fac = rho_safety * C_T / h_k;
  }
  else if(visc_stab_trace_estimate == INPAR::XFEM::ViscStab_TraceEstimate_eigenvalue)
  {

    // get an estimate of the hp-depending constant C_T satisfying the trace inequality from a previously solved eigenvalue problem
    NIT_visc_stab_fac = my::fldpara_->Get_TraceEstimate_MaxEigenvalue(sid);

#if(0)
    std::cout.precision(15);
    std::cout << "C_T/hk (formula): "
              << NIT_getTraceEstimateConstant(ele_distype)/h_k
              << " max_eigenvalue ~ C_T/hk: "
              <<  my::fldpara_->Get_TraceEstimate_MaxEigenvalue(sid)
              << " max_eigenvalue*h_k = C_T: "<< my::fldpara_->Get_TraceEstimate_MaxEigenvalue(sid)*h_k << " vs: C_T (formula) " << NIT_getTraceEstimateConstant(ele_distype) << std::endl;
#endif

  }
  else dserror("unknown trace-inequality-estimate type for viscous part of Nitsche's penalty term");

  // get final viscous scaling part of Nitsche's penalty term from a previously solved eigenvalue problem

  //--------------------------------------
  // scale the viscous part of the penalty scaling with the maximal viscosity of both sides
  NIT_visc_stab_fac *= my::visceff_;

  //--------------------------------------
  // scale the viscous part of the penalty scaling with the dimensionless user defined Nitsche-parameter gamma
  NIT_visc_stab_fac *= fldparaxfem_->ViscStabGamma();


  return NIT_visc_stab_fac;
}


/*--------------------------------------------------------------------------------
 * get the constant which satisfies the trace inequality depending on the spatial dimension and polynomial order of the element
 *--------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
double FluidEleCalcXFEM<distype>::NIT_getTraceEstimateConstant(
    const DRT::Element::DiscretizationType ele_distype
)
{
  /*
  => return
       CT = p(p+1)/2 * (2+1/p)^d  OR an estimate obtained form solving an EVP ( for tensor-product elements (hexahedral,quadrilateral...)
       CT = (p+1)*(p+d)/d                                                     ( for simplices (lines, triangles, tetrahedral...)
     constant depending on element-type/polynomial order, see below

  - C is an estimate of the following trace inequality with  C^2 ~ 1/h
     and C depends on element type, shape and polynomial degree


            trace inequality to be satisfied for Nitsche's method

            || grad v_h ||         <=   C  *  || grad v_h ||
                          Gamma_K                           K

   - we can estimate C^2 as C^2 \approx \rho_safety * CT / h_K (see generalized hp-framework by Erik Burman 2006)

              with CT depending on the dimension d and the polynomial degree p of grad(v_h)

              - For tensor-product elements (hexahedral,quadrilateral,...)
                   CT = p(p+1)/2 * (2+1/p)^d
                      -> REMARK: the theoretical estimate of the constant CT for hexahedral and quadrilateral elements seems to overestimate the actual constant
                      -> therefore we take the approximate values from solving a local eigenvalue problem...)
              - for simplices (lines, triangles, tetrahedral...)
                   CT = (p+1)*(p+d)/d        ( this estimate is sharp and leads to the same constants as solving an EVP)

              and  1/h_K \approx meas(Gamma_K)/meas(K)
              and rho_safety a safety factor depending on the element-type which is useful when the surface cuts the element in an inclined situation or when the surface
                  contains small acute angles or corners

   - Remark regarding the dimension of the problem (pseudo-2D simulations):
     -> important for CT is the degree of the element-wise polynomial which is clear for pure 2D and 3D problems
     -> for pseudo-2D with one element-layer of 3D element in z-direction, but strong Dirichlet conditions which eliminate the z-contributions
        in the polynomial, we have to adapt the constant from d=3 to d=2 although a 3D element is used
        (use the "is_pseudo_2D"-flag)
        => this ensures same results of a 2D simulation and a pseudo-2D simulation using the 3D implementation in combination with
           strong Dirichlet conditions in z-direction
     -> an important property is that the elongation of the element in z-direction does not play a role when surface-element-measure based
        hk-definitions are used since then the relation meas(Gamma_K)/meas(K) remains constant
     -> using non-strong Dirichlet conditions for the enforcement of u_z = 0, CT for 3D is required which then leads
        to differences between a 2D implementation and pseudo-2D simulation using the 3D implementation
     -> IMPORTANT: to forget to set the "is_pseudo_2D"-flag to "true" for pseudo-2D examples leads to an overestimate of CT,
        and so does not result in stability problems, it might lead to a slightly worse convergence of the linear solver,
        however, usually it has not a to worse effect.
*/

  // TODO: introduce 2D-flag
  bool is_pseudo_2D = fldparaxfem_->IsPseudo2D();

  DRT::Element::DiscretizationType trace_inequality_distype;

  if(!is_pseudo_2D)
  {
    // the standard case for real 2D or 3D simulations
    trace_inequality_distype = ele_distype;
  }
  else
  {
    // modification for pseudo 2D simulations
    switch (ele_distype)
    {
    case DRT::Element::hex8:
      trace_inequality_distype = DRT::Element::quad4; break;    // hex8 -> quad4 reduction
    case DRT::Element::hex20:
      trace_inequality_distype = DRT::Element::quad8; break;    // hex20-> quad8 reduction
    case DRT::Element::hex27:
      trace_inequality_distype = DRT::Element::quad9; break;    // hex27-> quad9 reduction
    case DRT::Element::wedge15:
      trace_inequality_distype = DRT::Element::tri6; break;     // wedge15 -> tri6 reduction (tri6 elements in 2D plane)
    case DRT::Element::wedge6:
      trace_inequality_distype = DRT::Element::tri3; break;     // wedge6 -> tri3 reduction (tri3 elements in 2D plane)
    default:
    {
      dserror("not a valid pseudo 2D element-type - what to do?");
      trace_inequality_distype = DRT::Element::dis_none;
      break;
    }
    };
  }

//  return 2.0;

  // switch over the distype which determines the right polynomial degree and dimension
  switch (trace_inequality_distype)
  {
  // for triangluar/tetradhedral elements:
  //    -> the grad-operator reduces the polynomial degree p(grad(v_h)) = p(v_h)-1
  // for quadrilateral/hexahedral and wedge/pyramid elements:
  //    -> the grad-operator does not reduce the polynomial degree due to the mixed polynomials p(grad(v_h)) = p(v_h)

  /*
  // CT = p(p+1)/2 * (2+1/p)^d
  // 3D hexahedral elements
  case DRT::Element::hex8:      return 27.0;                 break;  /// d=3, p(v_h) = 1 -> p(grad(v_h)) = 1   =>  CT=27.0
  case DRT::Element::hex20:     return 46.875;               break;  /// d=3, p(v_h) = 2 -> p(grad(v_h)) = 2   =>  CT=46.875 (375/8)
  case DRT::Element::hex27:     return 46.875;               break;  /// d=3, p(v_h) = 2 -> p(grad(v_h)) = 2   =>  CT=46.875 (375/8)
  // 2D quadrilateral elements
  case DRT::Element::quad4:     return 9.0;                  break;  /// d=2, p(v_h) = 1 -> p(grad(v_h)) = 1   =>  CT=9.0
  case DRT::Element::quad8:     return 9.0;                  break;  /// d=2, p(v_h) = 1 -> p(grad(v_h)) = 1   =>  CT=9.0
  case DRT::Element::quad9:     return 18.75;                break;  /// d=2, p(v_h) = 2 -> p(grad(v_h)) = 2   =>  CT=18.75  (75/4)
  */
  // REMARK: the theroretical estimate of the constant CT for hexahedral and quadrilateral elements seems to overestimate the actual constant
  // -> therefore we take the approximate values from solving a local eigenvalue problem
  // estimates from solving the eigenvalue problem on regular elements
  case DRT::Element::hex8:      return 1.59307;              break;  /// d=3, p(v_h) = 1 -> p(grad(v_h)) = 1   =>  CT= from eigenvalue problem
  case DRT::Element::hex20:     return 4.10462;              break;  /// d=3, p(v_h) = 2 -> p(grad(v_h)) = 2   =>  CT= from eigenvalue problem
  case DRT::Element::hex27:     return 4.27784;              break;  /// d=3, p(v_h) = 2 -> p(grad(v_h)) = 2   =>  CT= from eigenvalue problem
  // 2D quadrilateral elements
  case DRT::Element::quad4:     return 1.43426;              break;  /// d=2, p(v_h) = 1 -> p(grad(v_h)) = 1   =>  CT= from eigenvalue problem
  case DRT::Element::quad8:     return 4.06462;              break;  /// d=2, p(v_h) = 1 -> p(grad(v_h)) = 1   =>  CT= from eigenvalue problem
  case DRT::Element::quad9:     return 4.19708;              break;  /// d=2, p(v_h) = 2 -> p(grad(v_h)) = 2   =>  CT= from eigenvalue problem
  //----------------------------------------------------------------
  // CT = (p+1)*(p+d)/d (this estimate leads to the same results as solving a local eigenvalue problem)
  // 3D tetrahedral elements
  case DRT::Element::tet4:      return 1.0;                  break;  /// d=3, p(v_h) = 1 -> p(grad(v_h)) = 0   =>  CT=1.0
  case DRT::Element::tet10:     return 2.6666666666666666;   break;  /// d=3, p(v_h) = 2 -> p(grad(v_h)) = 1   =>  CT=2.6666666666666666 (8/3)
  // 2D triangular elements
  case DRT::Element::tri3:      return 1.0;                  break;  /// d=2, p(v_h) = 1 -> p(grad(v_h)) = 0   =>  CT=1.0
  case DRT::Element::tri6:      return 3.0;                  break;  /// d=2, p(v_h) = 2 -> p(grad(v_h)) = 1   =>  CT=3.0
  // 1D line elements
  case DRT::Element::line2:     return 1.0;                  break;  /// d=1, p(v_h) = 1 -> p(grad(v_h)) = 0   =>  CT=1.0
  case DRT::Element::line3:     return 4.0;                  break;  /// d=1, p(v_h) = 1 -> p(grad(v_h)) = 0   =>  CT=4.0
  //----------------------------------------------------------------
  // 3D wedge/pyramid elements, the current estimates are taken from the maximum of hex and tet elements, the correct value has to switch between the faces!
  case DRT::Element::pyramid5:  std::cout << "WARNING: calibrate this value!" << std::endl; return 1.59307;              break;  /// d=3, p(v_h) = 1 -> p(grad(v_h)) = 1   =>  CT taken from hex8
  case DRT::Element::wedge6:    std::cout << "WARNING: calibrate this value!" << std::endl; return 1.59307;              break;  /// d=3, p(v_h) = 1 -> p(grad(v_h)) = 1   =>  CT taken from hex8
  case DRT::Element::wedge15:   std::cout << "WARNING: calibrate this value!" << std::endl; return 4.10462;              break;  /// d=3, p(v_h) = 2 -> p(grad(v_h)) = 2   =>  CT taken from hex20
  default:
    dserror("constant for trace inequality not specified for this element type yet: % i", trace_inequality_distype); break;
  };

  return 0.0;
}


/*--------------------------------------------------------------------------------
 *      compute stabilization factor for the Nitsche's penalty term (xfluid/xfsi)
 *--------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void FluidEleCalcXFEM<distype>::NIT_ComputeFullStabfac(
    double &                         NIT_full_stab_fac,     ///< to be filled: full Nitsche's penalty term scaling (viscous+convective part)
    const double                     NIT_visc_stab_fac,     ///< Nitsche's viscous scaling part of penalty term
    INPAR::XFEM::ConvStabScaling     conv_stab_scaling,     ///< type of convective stabilization for xfluid
    const double                     veln_normal            ///< interface-normal velocity contribution
)
{
  if(my::fldpara_->IsConservative() and conv_stab_scaling != INPAR::XFEM::ConvStabScaling_none)
  {
    dserror("convective stabilization is not available for conservative form of Navier-Stokes, but possible to implement!");
  }

  // additional stabilization for convective stabilization
  double conv_stab_fac = 0.0;

  if (conv_stab_scaling == INPAR::XFEM::ConvStabScaling_abs_normal_vel)
  {
    //      | u*n |
    conv_stab_fac =  fabs(veln_normal);
  }
  else if (conv_stab_scaling == INPAR::XFEM::ConvStabScaling_max_abs_normal_vel)
  {
    //    max(1.0,| u*n |)
    conv_stab_fac =  std::max(1.0, fabs(veln_normal));
  }
  else if (conv_stab_scaling == INPAR::XFEM::ConvStabScaling_inflow)
  {
    //      ( -u*n ) if (u*n)<0 (inflow), conv_stabfac >= 0
    conv_stab_fac = std::max(0.0,-veln_normal);
  }
  else if (conv_stab_scaling == INPAR::XFEM::ConvStabScaling_none)
  {
    conv_stab_fac = 0.0;
  }
  else
  dserror("No valid INPAR::XFEM::ConvStabScaling for xfluid/xfsi problems");

  // Nitsche penalty scaling, Combined convective and viscous scaling

  //=================================================================================
  // definition in Burman 2007
  // Interior penalty variational multiscale method for the incompressible Navier-Stokes equation:
  // Monitoring artificial dissipation
  /*
  //      viscous_Nitsche-part, convective inflow part
  //
  //                    mu                                               /       _      \
  //  max( gamma_Nit * ----  , gamma_conv * | u_h * n | )       *       |  u_h - u, v_h  |
  //                    h_k                                              \              /
  */

  // final stabilization factors (viscous + convective + transient parts)
  NIT_full_stab_fac = std::max(NIT_visc_stab_fac, conv_stab_fac*my::densaf_);
  //=================================================================================

  return;
}

/*--------------------------------------------------------------------------------
 *    compute stabilization factor for Nitsche's method
 *--------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void FluidEleCalcXFEM<distype>::NIT_ComputeFullStabfacFluidFluid(
    double &                         NIT_full_stab_fac,             ///< to be filled: full Nitsche's penalty term scaling (viscous+convective part)
    double &                         avg_conv_stab_fac,             ///< to be filled: scaling of the convective average coupling term for fluidfluid problems
    double &                         velgrad_interface_fac,         ///< to be filled: stabilization parameter for velocity-gradients penalty term
    double &                         press_coupling_fac,            ///< to be filled: stabilization parameter for pressure-coupling penalty term
    INPAR::XFEM::XFF_ConvStabScaling xff_conv_stab_scaling,         ///< type of convective stabilization for fluid-fluid problem
    const double                     veln_normal,                   ///< interface-normal velocity contribution
    const LINALG::Matrix<my::nsd_,1> velint,                        ///< interface velocity
    const double                     h_k,                           ///< characteristic element length
    const double                     NIT_visc_stab_fac,             ///< Nitsche's viscous scaling part of penalty term
    const double                     gamma_ghost_penalty,           ///< ghost-penalty parameter
    const double                     gamma_press_coupling           ///< factor of stabilization parameter for pressure-coupling penalty term
)
{
  if(my::fldpara_->IsConservative() and xff_conv_stab_scaling != INPAR::XFEM::XFF_ConvStabScaling_none)
  {
    dserror("convective stabilization is not available for conservative form of Navier-Stokes, but possible to implement!");
  }

  //----------------------------------------------------------
  // additional average term stabilization for convective regime

  /*
  //    /                                           \        /                        i                   \
  //   |  gamma_conv * ( u * n^e ) *  { v } , [ Du ] | =   - |  gamma_conv*  ( u * n^e ) * { v }, [ u ]   |
  //    \                                           /        \                                            /
  //
  */

  if(xff_conv_stab_scaling == INPAR::XFEM::XFF_ConvStabScaling_onesidedinflow or
     xff_conv_stab_scaling == INPAR::XFEM::XFF_ConvStabScaling_averaged or
     xff_conv_stab_scaling == INPAR::XFEM::XFF_ConvStabScaling_onesidedinflow_max_penalty or
     xff_conv_stab_scaling == INPAR::XFEM::XFF_ConvStabScaling_averaged_max_penalty)
  {
    // (beta .normal_e)({v},||u||)
    avg_conv_stab_fac = my::densaf_*veln_normal;
  }
  else if (xff_conv_stab_scaling == INPAR::XFEM::XFF_ConvStabScaling_none)
  {
    avg_conv_stab_fac = 0.0;
  }
  else
    dserror("No valid INPAR::XFEM::XFF_ConvStabScaling for xfluidfluid/xffsi problems");


  //----------------------------------------------------------------

  // set the viscous part for the Nitsche penalty term
  NIT_full_stab_fac = NIT_visc_stab_fac;

  // if additional penalty factor is on..
  // stabfac_visc = max(my::densaf_*|u*n|_inf,stabfac_visc);
  if (xff_conv_stab_scaling == INPAR::XFEM::XFF_ConvStabScaling_averaged_max_penalty or
      xff_conv_stab_scaling == INPAR::XFEM::XFF_ConvStabScaling_onesidedinflow_max_penalty)
  {
    NIT_full_stab_fac = std::max(my::densaf_*fabs(velint.NormInf()), NIT_full_stab_fac);
  }

  // if ConvStabScaling_inflow the Nitsche penalty term obtains also a convective contribution
  // (\rho 0.5 |beta.n^e| ||v||,||u||)
  if (xff_conv_stab_scaling == INPAR::XFEM::XFF_ConvStabScaling_onesidedinflow or
      xff_conv_stab_scaling == INPAR::XFEM::XFF_ConvStabScaling_onesidedinflow_max_penalty)
  {
    NIT_full_stab_fac += fabs(veln_normal)*my::densaf_*0.5;
  }


  //----------------------------------------------------------------


  // stabilization parameter for velocity-gradients penalty term
  velgrad_interface_fac = gamma_ghost_penalty*my::visceff_*h_k;


  // stabilization parameter for pressure-coupling penalty term
  // no scaling with density here. The scaling with the density follows...


  // parameter \gamma_p of pressure-stabilizing term
  /*
   *   /                                  \
   *   | \gamma_p * \rho^-1 * [ Dq] [ Dp ] |
   *   \                                  /
   */
  press_coupling_fac = gamma_press_coupling;
  double kinvisc = my::visceff_/my::densaf_;

  //  nu-weighted definition
  // \gamma_p = \alpha_p * h_k^2 / ( 1+ \nu/h_k )
  press_coupling_fac /= (1.0 + (kinvisc/h_k));

  // to have a consistent formulation we need only one density factor for
  // pressure coupling. That is because we have two times pressure (test
  // functions and the shape function) in the formuation. If we do not
  // cross out one density, we would multiply the term two times with
  // density, which is not correct.
  press_coupling_fac /= my::densaf_;

  return;
}

/*--------------------------------------------------------------------------------
 * prepare coupling matrices, that include contributions from convective stabilization
 *--------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void FluidEleCalcXFEM<distype>::HybridLM_CreateConvStabMatrices(
    std::set<int> &                                         begids,                  ///< ids of intersecting boundary elements
    DRT::Discretization &                                   cutdis,                  ///< cut discretization
    std::map<int, std::vector<Epetra_SerialDenseMatrix> > & conv_side_coupling       ///< contributions to coupling matrices from convective stabilizations
)
{
  if (fldparaxfem_->GetCouplingMethod() != INPAR::XFEM::Hybrid_LM_Cauchy_stress &&
      fldparaxfem_->GetCouplingMethod() != INPAR::XFEM::Hybrid_LM_viscous_stress)
    dserror("Do not call this method with a non-Lagrange multiplier based approach!");

  for (std::set<int>::const_iterator bgid=begids.begin(); bgid!=begids.end(); ++bgid)
  {
    DRT::Element * side = cutdis.gElement(*bgid); // for each boundary element there is one corresponding side
    std::vector<int> patchlm;
    std::vector<int> patchlmowner;
    std::vector<int> patchlmstride;
    side->LocationVector(cutdis, patchlm, patchlmowner, patchlmstride);

    // get coupling matrices for the current side (boundary element)
    std::vector<Epetra_SerialDenseMatrix> & conv_side_matrices = conv_side_coupling[*bgid];

    conv_side_matrices.resize(4);
    conv_side_matrices[0].Shape(patchlm.size(), my::nen_*my::numdofpernode_); //Cuiu
    conv_side_matrices[1].Shape(my::nen_*my::numdofpernode_, patchlm.size()), //Cuui
    conv_side_matrices[2].Shape(patchlm.size(),1);                            //rhs_Cui
    conv_side_matrices[3].Shape(patchlm.size(), patchlm.size());              //Cuiui
  }
}

/*--------------------------------------------------------------------------------
 * create location vector w.r.t patch of intersecting boundary elements and reshape coupling matrices
 *--------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void FluidEleCalcXFEM<distype>::PatchLocationVector(
    std::set<int> &                                         begids,                  ///< ids of intersecting boundary elements
    DRT::Discretization &                                   cutdis,                  ///< cut discretization
    std::vector<int> &                                      patchelementslmv,        ///< lm vector for patch of boundary elements
    std::vector<int> &                                      patchelementslmowner,    ///< lmowner vector for patch of boundary elements
    std::map<int, std::vector<Epetra_SerialDenseMatrix> > & Cuiui_coupling           ///< coupling matrices
)
{
  for (std::set<int>::const_iterator bgid=begids.begin(); bgid!=begids.end(); ++bgid)
  {
    DRT::Element * side = cutdis.gElement(*bgid); // for each boundary element there is one corresponding side
    std::vector<int> patchlm;
    std::vector<int> patchlmowner;
    std::vector<int> patchlmstride;
    side->LocationVector(cutdis, patchlm, patchlmowner, patchlmstride);

    patchelementslmv.reserve( patchelementslmv.size() + patchlm.size());
    patchelementslmv.insert(patchelementslmv.end(), patchlm.begin(), patchlm.end());

    patchelementslmowner.reserve( patchelementslmowner.size() + patchlmowner.size());
    patchelementslmowner.insert( patchelementslmowner.end(), patchlmowner.begin(), patchlmowner.end());

    // get coupling matrices for the current side (boundary element)
    std::vector<Epetra_SerialDenseMatrix> & Cuiui_matrices = Cuiui_coupling[*bgid];

    INPAR::XFEM::CouplingMethod coupl_meth = fldparaxfem_->GetCouplingMethod();

    if (coupl_meth == INPAR::XFEM::Hybrid_LM_viscous_stress ||
        coupl_meth == INPAR::XFEM::Hybrid_LM_Cauchy_stress)
    {
      Cuiui_matrices.resize(2);
      Cuiui_matrices[0].Shape(my::nen_*numstressdof_,patchlm.size()); //Gsui (coupling between background elements sigma and current side!)
      Cuiui_matrices[1].Shape(patchlm.size(),my::nen_*numstressdof_); //Guis
    }
    else if (coupl_meth == INPAR::XFEM::Nitsche)
    {
      Cuiui_matrices.resize(1);
      Cuiui_matrices[0].Shape(patchlm.size(),patchlm.size()); //Cuiui
    }
    else dserror("Type of coupling method unknown.");
  }

  return;
}

/*--------------------------------------------------------------------------------
 * create location vector w.r.t patch of intersecting boundary elements and
 * reshape coupling matrices (for embedded coupling)
 *--------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void FluidEleCalcXFEM<distype>::PatchLocationVector(
    std::set<int> &                                         begids,                  ///< ids of intersecting boundary elements
    DRT::Discretization &                                   alediscret,              ///< cut discretization
    std::vector<int> &                                      patchelementslmv,        ///< lm vector for patch of boundary elements
    std::vector<int> &                                      patchelementslmowner,    ///< lmowner vector for patch of boundary elements
    std::map<int, std::vector<Epetra_SerialDenseMatrix> > & Cuiui_coupling,          ///< coupling matrices
    std::map<int,int> &                                     boundary_emb_gid_map     ///< map between boundary sid and corresponding embedded element id
)
{
  for (std::set<int>::const_iterator bgid=begids.begin(); bgid!=begids.end(); ++bgid)
  {
    DRT::Element * emb_ele = alediscret.gElement(boundary_emb_gid_map.find(*bgid)->second);

    std::vector<int> patchlm;
    std::vector<int> patchlmowner;
    std::vector<int> patchlmstride;
    emb_ele->LocationVector(alediscret, patchlm, patchlmowner, patchlmstride);

    patchelementslmv.reserve( patchelementslmv.size() + patchlm.size());
    patchelementslmv.insert(patchelementslmv.end(), patchlm.begin(), patchlm.end());

    patchelementslmowner.reserve( patchelementslmowner.size() + patchlmowner.size());
    patchelementslmowner.insert( patchelementslmowner.end(), patchlmowner.begin(), patchlmowner.end());

    // get coupling matrices for the current side (boundary element)
    std::vector<Epetra_SerialDenseMatrix> & Cuiui_matrices = Cuiui_coupling[*bgid];

    INPAR::XFEM::CouplingMethod coupl_meth = fldparaxfem_->GetCouplingMethod();
    if (coupl_meth == INPAR::XFEM::Nitsche)
    {
      Cuiui_matrices.resize(1);
      Cuiui_matrices[0].Shape(patchlm.size(),patchlm.size()); //Cuiui
    }
    else if (coupl_meth == INPAR::XFEM::Hybrid_LM_viscous_stress ||
             coupl_meth == INPAR::XFEM::Hybrid_LM_Cauchy_stress)
    {
      dserror("Embedded-sided coupling for stress-based hybrid LM approach is not yet available!");
    }
    else dserror("Not supported coupling method.");
  }

  return;
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
  LINALG::Matrix<my::nsd_,my::nsd_> eps(true);
  for(int i=0; i<my::nsd_; i++)
  {
    for(int j=0; j<my::nsd_; j++)
    {
      eps(i,j) = 0.5 * (my::vderxy_(i,j) + my::vderxy_(j,i));
    }
  }
  //-------------------------------

  // t = ( -pI + 2mu eps(u) )*n^f
  traction.Clear();
  traction.Multiply(eps, normal);

  // add the pressure part and scale the viscous part with the viscosity
  traction.Update( -press, normal, 2.0*my::visceff_);

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
  const Epetra_Map* dofcolmap = cutdis.DofColMap();

  for (int idof = 0; idof < (int)(lm.size()); ++idof)
  {
      int gdof = lm[idof];

      // f^i = ( N^i, t ) = ( N^i, (-pI+2mu*eps(u))*n )
      (*iforcecol)[dofcolmap->LID(gdof)] += iforce[idof];
  }

  return;
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
    const GEO::CUT::plain_volumecell_set &                                vcSet,                 ///< volumecell sets for volume integration
    const std::map<int, std::vector<GEO::CUT::BoundaryCell*> > &          bcells,                ///< bcells for boundary cell integration
    const std::map<int, std::vector<DRT::UTILS::GaussIntegration> > &     bintpoints,            ///< integration points for boundary cell integration
    Teuchos::RCP<DRT::ELEMENTS::XFLUID::SlaveElementInterface<distype> >  emb,                   ///< pointer to the embedded coupling implementation
    DRT::Element *                                                        face                   ///< side element in 3D
)
{
  const INPAR::XFEM::ViscStab_hk visc_stab_hk = fldparaxfem_->ViscStabHK();

  const INPAR::XFEM::CouplingStrategy coupling_strategy = fldparaxfem_->GetCouplingStrategy();

  if (emb == Teuchos::null and coupling_strategy == INPAR::XFEM::Embedded_Sided_Coupling)
    dserror("no EmbCoupling available, however Embedded_Sided_Coupling is activated!");

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
    if(coupling_strategy == INPAR::XFEM::Embedded_Sided_Coupling)
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
    if(coupling_strategy == INPAR::XFEM::Embedded_Sided_Coupling)
      dserror("ViscStab_hk_cut_vol_div_by_cut_surf not reasonable for Embedded_Sided_Coupling!");

    // compute the cut surface measure
    meas_surf = ComputeMeasCutSurf(bintpoints, bcells);

    if (fabs(meas_surf) < 1.e-8)  dserror("Element contribution to interface has zero size.");

    // compute the cut volume measure
    for( GEO::CUT::plain_volumecell_set::const_iterator i=vcSet.begin();i!=vcSet.end();i++ )
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
    if(coupling_strategy == INPAR::XFEM::Embedded_Sided_Coupling)
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
    if(coupling_strategy != INPAR::XFEM::Embedded_Sided_Coupling)
      dserror("ViscStab_hk_ele_vol_div_by_ele_surf just reasonable for Embedded_Sided_Coupling!");

    //---------------------------------------------------
    // find the corresponding local id of the face of the element
    // REMARK: this is quite slow, however at the moment the easiest way to get the local id
    //         here is space for improvement

    const int lid = fldparaxfem_->Get_face_lid_of_embedded_ele(face->Id());
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
    if(coupling_strategy == INPAR::XFEM::Embedded_Sided_Coupling)
      dserror("ViscStab_hk_ele_vol_div_by_max_ele_surf not reasonable for Embedded_Sided_Coupling!");

    // compute the uncut element's surface measure
    const int numfaces = DRT::UTILS::getNumberOfElementFaces(ele->Shape());

    // loop all surfaces
    for(int lid=0; lid< numfaces; lid++)
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
  if( h_k <= 0.0 ) dserror("element length is <= 0.0");

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
  for(int n=0; n<numnode_face; n++)
  {
    const int node_lid = map[local_face_id][n];
    for(int idim = 0; idim < my::nsd_; idim++)
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
  for (int gpid=0; gpid<intpoints.nquad; gpid++)
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

#ifdef BOUNDARYCELL_TRANSFORMATION_OLD
        dserror("at the moment not available -> fix it");
//        si->Evaluate(eta,x_side,normal,drs);

#else
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
          for (int idim=0;idim<3;idim++)
          {
            x_gp_lin(idim,0) = gpcord[idim];
          }
        }

 #endif
        const double surf_fac = drs*iquad.Weight();

        surf += surf_fac;

      } //loop gausspoints w.r.t current boundary cell
    } // loop boundary cells
  } // loop intersecting sides

  return surf;
}

  } // end namespace ELEMENTS
} // end namespace DRT



// Ursula is responsible for this comment!
template class DRT::ELEMENTS::FluidEleCalcXFEM<DRT::Element::hex8>;
template class DRT::ELEMENTS::FluidEleCalcXFEM<DRT::Element::hex20>;
template class DRT::ELEMENTS::FluidEleCalcXFEM<DRT::Element::hex27>;
template class DRT::ELEMENTS::FluidEleCalcXFEM<DRT::Element::tet4>;
template class DRT::ELEMENTS::FluidEleCalcXFEM<DRT::Element::tet10>;
//template class DRT::ELEMENTS::FluidEleCalcXFEM<DRT::Element::wedge6>;
//template class DRT::ELEMENTS::FluidEleCalcXFEM<DRT::Element::pyramid5>;


