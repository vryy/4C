/*!
\file turbulence_statistic_manager.H

\brief Manage the computation of averages for several 
canonical flows like channel flow, flow around a square 
cylinder, flow in a lid driven cavity etc.

The manager is intended to remove as much of the averaging
overhead as possible from the time integration method.

<pre>

Maintainer: Peter Gamnitzer
            gamnitzer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15235

</pre>
*/
#ifdef CCADISCRET
#ifndef TURBULENCE_STATISTIC_MANAGER_H
#define TURBULENCE_STATISTIC_MANAGER_H

#include "../drt_fluid/fluid_genalpha_integration.H"
#include "../drt_fluid/fluidimplicitintegration.H"
#include "../drt_fluid/fluid_utils.H"

#include "../drt_fluid/turbulence_statistics_mean_general.H"
#include "../drt_fluid/turbulence_statistics_ccy.H"
#include "../drt_fluid/turbulence_statistics_cha.H"
#include "../drt_fluid/turbulence_statistics_ldc.H"
#include "../drt_fluid/turbulence_statistics_bfs.H"
#include "../drt_fluid/turbulence_statistics_sqc.H"


using namespace std;
using namespace Teuchos;

namespace FLD
{

class FluidGenAlphaIntegration;
class FluidImplicitTimeInt;

class TurbulenceStatisticManager
{
public:

  /*!
  \brief Genalpha Constructor (public)

  */
  TurbulenceStatisticManager(FluidGenAlphaIntegration& timeint);

  /*!
  \brief One-Step-Theta Constructor (public)

  */
  TurbulenceStatisticManager(FluidImplicitTimeInt& timeint);

  /*!
  \brief Destructor

  */
  virtual ~TurbulenceStatisticManager();


  //! @name managing methods 

  /*!
  \brief Store values computed during the element call

  (including Cs, visceff etc for a dynamic 
   Smagorinsky model)

  */
  void StoreElementValues(int step);

  /*!
  \brief Include current quantities in the time 
  averaging procedure

  */
  void DoTimeSample(int step, double time, const double eosfac);

  /*!
  \brief Write (dump) the statistics to a file

  */
  void DoOutput(IO::DiscretizationWriter& output, int step, const double eosfac);

  /*!
  \brief Clear all statistics collected up to now

  */
  void Reset();

  /*!
  \brief Restart collection of statistics

  */
  void Restart(IO::DiscretizationReader& reader,int step);

  //@}

private:
  /*!
  \brief Time integration independent setup called by Constructor

  */
  void Setup();


  //! time step size
  double                               dt_;
  //! parameter for generalised alpha --- alpha "momentum"
  double                               alphaM_;
  //! parameter for generalised alpha --- alpha "force"
  double                               alphaF_;
  //! third parameter for generalised alpha
  double                               gamma_;

  //! reference to material density to rescale pressure
  double &                             density_;

  //! parameters for sampling/dumping period
  //! start of sampling at step samstart
  int                                  samstart_;
  //! stop sampling at step samstop
  int                                  samstop_;
  //! incremental dump every dumperiod steps or standalone 
  //! records (0)
  int                                  dumperiod_  ;

  //! the discretization
  RCP<DRT::Discretization>             discret_    ;

  //! parameterlist of the discretization including time params, 
  //! stabilization params and turbulence sublist
  ParameterList&                       params_     ;

  //! parameterlist specially designed for the evaluation of 
  //! gausspoint statistics
  ParameterList                        eleparams_  ;

  //! decides whether we use an Eulerian or an ALE formulation
  bool                                 alefluid_   ;

  //! all my solution vectors needed for element evaluation
  RCP<Epetra_Vector>                   myaccnp_    ;
  RCP<Epetra_Vector>                   myaccn_     ;
  RCP<Epetra_Vector>                   myaccam_    ;

  RCP<Epetra_Vector>                   myvelnp_    ;
  RCP<Epetra_Vector>                   myveln_     ;
  RCP<Epetra_Vector>                   myvelaf_    ;

  RCP<Epetra_Vector>                   myvedeam_   ;
  RCP<Epetra_Vector>                   myvedenp_   ;
  RCP<Epetra_Vector>                   myveden_    ;

  RCP<Epetra_Vector>                   mydispnp_   ;
  RCP<Epetra_Vector>                   mydispn_    ;

  RCP<Epetra_Vector>                   mygridveln_ ;
  RCP<Epetra_Vector>                   mygridvelaf_;

  RCP<Epetra_Vector>                   mysubgrvisc_;

  RCP<Epetra_Vector>                   myforce_    ;


  //! specifies the special flow
  enum SpecialFlow
  { no_special_flow                 ,
    channel_flow_of_height_2        ,
    loma_channel_flow_of_height_2   ,
    lid_driven_cavity               ,
    loma_lid_driven_cavity          ,
    backward_facing_step            ,
    loma_backward_facing_step       ,
    square_cylinder                 ,
    square_cylinder_nurbs           ,
    rotating_circular_cylinder_nurbs}    flow_;

  //! toggle evaluation of dynamic Smagorinsky/Smagorinsky with 
  //! wall damping quantities
  bool                                 smagorinsky_;

  //! toggle evaluation of subgrid quantities, dissipation rates etc
  //! this is only possible for the genalpha implementation since 
  //! we need a corresponding element implementation
  bool                                 subgrid_dissipation_;

  //! mean values of velocity and pressure, independent of special flow
  //! averaging takes place in time, if hom. directions have been specified
  //! additionally along these lines
  RCP<TurbulenceStatisticsGeneralMean> statistics_general_mean_;

  //! turbulence statistics for turbulent channel flow
  RCP<TurbulenceStatisticsCha>         statistics_channel_;

  //! turbulence statistics for a rotating circular cylinder
  RCP<TurbulenceStatisticsCcy>         statistics_ccy_;

  //! turbulence statistics for lid-driven cavity
  RCP<TurbulenceStatisticsLdc>         statistics_ldc_;

  //! turbulence statistics for backward-facing step
  RCP<TurbulenceStatisticsBfs>         statistics_bfs_;

  //! turbulence statistics for square cylinder
  RCP<TurbulenceStatisticsSqc>         statistics_sqc_;

}; // end class TurbulenceStatisticManager

} // end namespace FLD

#endif  // #ifndef TURBULENCE_STATISTIC_MANAGER_H
#endif  // #ifdef CCADISCRET
