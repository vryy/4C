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

#include "Epetra_Vector.h"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

using namespace std;
using namespace Teuchos;

//forward declarations
namespace DRT
{
  class Discretization;
  class DofSet;
}

namespace IO
{
  class DiscretizationReader;
  class DiscretizationWriter;
}

namespace COMBUST
{
  class TurbulenceStatisticsORACLES;
  class TurbulenceStatisticsBcf;
}

namespace FLD
{
// forward declarations
class FluidGenAlphaIntegration;
class FluidImplicitTimeInt;
class CombustFluidImplicitTimeInt;
class TurbulenceStatisticsGeneralMean;
class TurbulenceStatisticsCha;
//class TurbulenceStatisticsBcf;
class TurbulenceStatisticsCcy;
class TurbulenceStatisticsLdc;
class TurbulenceStatisticsBfs;
class TurbulenceStatisticsSqc;


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
  \brief Combustion Constructor (public)

  */
  TurbulenceStatisticManager(CombustFluidImplicitTimeInt& timeint);

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
  \brief Store values computed during the element call

  (including subfilter stresses for a scale
  similarity model)

  */
  void StoreNodalValues(int step, const RCP<Epetra_Vector>   stress12);

  /*!
  \brief Include current quantities in the time
  averaging procedure

  Include current quantities in the time
  averaging procedure. The velnp parameters with
  the default of Teuchos::null can be used in
  case:
  1. discretmatchingvelnp: the discretization has a varying
                           number of dofs.
  2. velnp:                instead of updating the original
                           velnp to which myvelnp_ points a
                           new velnp is created.

  \param step (in)       : current time step
  \param time (in)       : current time
  \param eosfac (in)     : equation of state factor
  \param velnp (in)      : velnp built with dofrowmap from standarddofset
  \param force (in)      : trueresidual vector built with standarddofset
  \param discretmatchingvelnp (in) : velnp built with dofrowmap from discretization

  */
  void DoTimeSample(int step                                                              ,
                    double time                                                           ,
                    const double eosfac                                                   );
//                    Teuchos::RCP<const Epetra_Vector> velnp                = Teuchos::null,
//                    Teuchos::RCP<const Epetra_Vector> force                = Teuchos::null,
//                    Teuchos::RCP<const DRT::DofSet>   stddofset            = Teuchos::null,
//                    Teuchos::RCP<const Epetra_Vector> discretmatchingvelnp = Teuchos::null,
//                    Teuchos::RCP<const Epetra_Vector> phinp                = Teuchos::null);

  /*!
  \brief function called by the combustion time integration

  */
  void DoTimeSample(int                             step,
                    double                          time,
                    Teuchos::RCP<Epetra_Vector>     velnp,
                    Teuchos::RCP<Epetra_Vector>     force,
                    Teuchos::RCP<Epetra_Vector>     phi,
                    Teuchos::RCP<const DRT::DofSet> stddofset
                    );

  /*!
  \brief Write (dump) the statistics to a file

  */
  void DoOutput(
       IO::DiscretizationWriter& output,
       int step,
       const double inflow = false);

  /*!
  \brief Restart collection of statistics

  */
  void Restart(IO::DiscretizationReader& reader,int step);

  /*!
  \brief Restart scatra-specific collection of statistics

  */
  void RestartScaTra(IO::DiscretizationReader& scatrareader,int step);

  /*!
  \brief Add results from scalar transport fields to statistics

  */
  void AddScaTraResults(
      RCP<DRT::Discretization> scatradis,
      RCP<Epetra_Vector> myphinp);

  /*!
  \brief   Write (dump) the scatra-specific mean fields to the result file

  */
  void DoOutputForScaTra(IO::DiscretizationWriter& output, int step);

  //@}

  /*!
  \brief remote access method to general mean statistics manager

  */
  const RCP<TurbulenceStatisticsGeneralMean> GetTurbulenceStatisticsGeneralMean() const
    { return statistics_general_mean_; }

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
  double                               density_;

  //! parameters for sampling/dumping period
  //! start of sampling at step samstart
  int                                  samstart_;
  //! stop sampling at step samstop
  int                                  samstop_;
  //! incremental dump every dumperiod steps or standalone
  //! records (0)
  int                                  dumperiod_  ;

  //! the fluid discretization
  RCP<DRT::Discretization>             discret_    ;
  //! the scatra discretization
  RCP<DRT::Discretization>             scatradis_  ;

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

  RCP<Epetra_Vector>                   myscanp_    ;
  //! scatra result vector (defined on the scatra dofrowmap!)
  RCP<Epetra_Vector>                   myfullphinp_;

  RCP<Epetra_Vector>                   mydispnp_   ;
  RCP<Epetra_Vector>                   mydispn_    ;

  RCP<Epetra_Vector>                   mygridveln_ ;
  RCP<Epetra_Vector>                   mygridvelaf_;

  RCP<Epetra_Vector>                   myforce_    ;

  RCP<Epetra_MultiVector>              myfilteredvel_;
  RCP<Epetra_MultiVector>              myfilteredreystr_;
  RCP<Epetra_Vector>                   myfsvelaf_    ;


  //! specifies the special flow
  enum SpecialFlow
  { no_special_flow                 ,
    channel_flow_of_height_2        ,
    loma_channel_flow_of_height_2   ,
    lid_driven_cavity               ,
    loma_lid_driven_cavity          ,
    backward_facing_step            ,
    loma_backward_facing_step       ,
    combust_oracles                 ,
    square_cylinder                 ,
    square_cylinder_nurbs           ,
    bubbly_channel_flow             ,
    rotating_circular_cylinder_nurbs,
    rotating_circular_cylinder_nurbs_scatra,
    time_averaging                           }    flow_;

  //! toggle evaluation of dynamic Smagorinsky/Smagorinsky with
  //! wall damping quantities
  bool                                 smagorinsky_;

  //! toggle evaluation of scale similarity quantities
  bool                                 scalesimilarity_;

  //! toggle evaluation of scale multifractal quantities
  bool                                 multifractal_;

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

  //! turbulence statistics for turbulent multi phase channel flow
  RCP<COMBUST::TurbulenceStatisticsBcf> statistics_channel_multiphase_;

  //! turbulence statistics for a rotating circular cylinder
  RCP<TurbulenceStatisticsCcy>         statistics_ccy_;

  //! turbulence statistics for lid-driven cavity
  RCP<TurbulenceStatisticsLdc>         statistics_ldc_;

  //! turbulence statistics for backward-facing step
  RCP<TurbulenceStatisticsBfs>         statistics_bfs_;

  //! turbulence statistics for ORACLES
  RCP<COMBUST::TurbulenceStatisticsORACLES> statistics_oracles_;

  //! turbulence statistics for square cylinder
  RCP<TurbulenceStatisticsSqc>         statistics_sqc_;

}; // end class TurbulenceStatisticManager

} // end namespace FLD

#endif  // #ifndef TURBULENCE_STATISTIC_MANAGER_H
#endif  // #ifdef CCADISCRET
