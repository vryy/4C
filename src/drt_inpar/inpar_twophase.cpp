/*----------------------------------------------------------------------*/
/*!
\file inpar_twophase.cpp

\brief Input parameters for combustion

\level 2

<pre>
\maintainer Magnus Winter
            winter@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
</pre>
*/

/*----------------------------------------------------------------------*/



#include "drt_validparameters.H"
#include "inpar_twophase.H"



void INPAR::TWOPHASE::SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  using namespace DRT::INPUT;
  using Teuchos::tuple;
  using Teuchos::setStringToIntegralParameter;

  Teuchos::ParameterList& twophasedyn = list->sublist("TWO PHASE FLOW",false,"");
  IntParameter("NUMSTEP",10,"Number of Time Steps",&twophasedyn);
  DoubleParameter("TIMESTEP",0.01,"Time increment dt",&twophasedyn);
  DoubleParameter("MAXTIME",0.0,"Total simulation time",&twophasedyn);
  DoubleParameter("CONVTOL",1E-6,"Tolerance for convergence check",&twophasedyn);
  IntParameter("RESULTSEVRY",1,"Increment for writing solution",&twophasedyn);
  IntParameter("RESTARTEVRY",1,"Increment for writing restart",&twophasedyn);
  IntParameter("ITEMAX",1,"Maximum number of iterations in levelset-fluid loop",&twophasedyn);
  BoolParameter("WRITE_CENTER_OF_MASS","No","Write center of mass to file",&twophasedyn);
  //TODO: Is this parameter even used?
  BoolParameter("RESTART_SCATRA_INPUT","No","Use ScaTra field from .dat-file instead",&twophasedyn);

  Teuchos::ParameterList& twophasedyn_smeared = twophasedyn.sublist("SMEARED",false,"");
  DoubleParameter("INTERFACE_THICKNESS",0.0,"Thickness of interface for multiphase flow",&twophasedyn_smeared);
  BoolParameter("ENHANCED_GAUSSRULE","No","Set higher order gaussrule within the interface layer.",&twophasedyn_smeared);


  Teuchos::ParameterList& twophase_surftens = twophasedyn.sublist("SURFACE TENSION",false,"");
  setStringToIntegralParameter<int>("SURFTENSAPPROX","surface_tension_approx_none","Type of surface tension approximation",
      tuple<std::string>(
          "surface_tension_approx_none",                         //none
          "surface_tension_approx_fixed_curvature",              //prescribed curvature, gamma = curv*gamma
          //"surface_tension_approx_divgrad",                      //calcs curvature at GP using the smoothed grad_phi and smoothed grad_phi for normal  //Do Not Migrate.
          "surface_tension_approx_divgrad_normal",               //calcs curvature at GP using the smoothed grad_phi and normal on Boundary Cell for normal
          "surface_tension_approx_nodal_curvature",              //calcs curvature at nodes and normal on Boundary Cell for normal
          "surface_tension_approx_laplacebeltrami"),              //standard Laplace-Beltrami (see e.g. Fries 2009)
          //"surface_tension_approx_laplacebeltrami_smoothed"),    //         Laplace-Beltrami, includes additional projection based on the smoothed normal (see e.g. Gross, Reusken 2009)
          tuple<int>(
              surface_tension_approx_none,
              surface_tension_approx_fixed_curvature,
              //surface_tension_approx_divgrad,
              surface_tension_approx_divgrad_normal,
              surface_tension_approx_nodal_curvature,
              surface_tension_approx_laplacebeltrami),
//               surface_tension_approx_laplacebeltrami_smoothed),
              &twophase_surftens);

  setStringToIntegralParameter<int>("SMOOTHGRADPHI","smooth_grad_phi_l2_projection","Type of smoothing for grad(phi)",
      tuple<std::string>(
          //"smooth_grad_phi_none",
          "smooth_grad_phi_meanvalue",
          "smooth_grad_phi_leastsquares_3D",
          "smooth_grad_phi_leastsquares_2Dx",
          "smooth_grad_phi_leastsquares_2Dy",
          "smooth_grad_phi_leastsquares_2Dz",
          "smooth_grad_phi_l2_projection"),
          tuple<int>(
              //smooth_grad_phi_none,
              smooth_grad_phi_meanvalue,
              smooth_grad_phi_leastsquares_3D,
              smooth_grad_phi_leastsquares_2Dx,
              smooth_grad_phi_leastsquares_2Dy,
              smooth_grad_phi_leastsquares_2Dz,
              smooth_grad_phi_l2_projection),
              &twophase_surftens);
  BoolParameter("L2_PROJECTION_SECOND_DERIVATIVES","No","L2 Projection Second Derivatives of Level Set",&twophase_surftens);

  setStringToIntegralParameter<int>("NODAL_CURVATURE","l2_projected","Type of calculation of nodal curvature value",
      tuple<std::string>(
          "l2_projected",
          "averaged"),
          tuple<int>(
              l2_projected,
              averaged),
              &twophase_surftens);

  setStringToIntegralParameter<int>("LAPLACE_BELTRAMI","matrix_mixed_smoothed","Type of calculation of Laplace-Beltrami projection matrix",
      tuple<std::string>(
          "matrix_non_smoothed",
          "matrix_smoothed",
          "matrix_mixed_smoothed"),
          tuple<int>(
              matrix_non_smoothed,
              matrix_smoothed,
              matrix_mixed_smoothed),
              &twophase_surftens);

  DoubleParameter("SMOOTHING_PARAMETER",0.0,"Diffusion Coefficient for Smoothing",&twophase_surftens); //Added diffusion for L2_projection

}
