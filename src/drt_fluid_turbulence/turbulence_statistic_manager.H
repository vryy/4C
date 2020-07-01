/*----------------------------------------------------------------------*/
/*! \file

\brief Manage the computation of averages for several
canonical flows like channel flow, flow around a square
cylinder, flow in a lid driven cavity etc.

The manager is intended to remove as much of the averaging
overhead as possible from the time integration method.


\level 2

*/
/*----------------------------------------------------------------------*/

#ifndef TURBULENCE_STATISTIC_MANAGER_H
#define TURBULENCE_STATISTIC_MANAGER_H

#include "Epetra_Vector.h"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "../drt_inpar/inpar_fluid.H"
#include "../drt_inpar/inpar_scatra.H"


// forward declarations
namespace DRT
{
  class Discretization;
  class DofSet;
}  // namespace DRT

namespace IO
{
  class DiscretizationReader;
  class DiscretizationWriter;
}  // namespace IO

namespace SCATRA
{
  class ScaTraTimIntImpl;
}

namespace FLD
{
  // forward declarations
  class FluidImplicitTimeInt;
  class XWall;
  class TurbulenceStatisticsCha;
  class TurbulenceStatisticsCcy;
  class TurbulenceStatisticsLdc;
  class TurbulenceStatisticsBfs;
  class TurbulenceStatisticsSqc;
  class TurbulenceStatisticsHit;
  class TurbulenceStatisticsTgv;
  class TurbulenceStatisticsGeneralMean;
  class TurbulenceStatisticsPh;
  class TurbulenceStatisticsBfda;
  namespace UTILS
  {
    class StressManager;
  }


  class TurbulenceStatisticManager
  {
   public:
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
    \brief Get current velnp for meshtying, since it may have changed

    */
    void GetCurrentVelnp(Teuchos::RCP<Epetra_Vector>);

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
    void DoTimeSample(int step, const double eosfac, const double thermpressaf = 0.0,
        const double thermpressam = 0.0, const double thermpressdtaf = 0.0,
        const double thermpressdtam = 0.0);

    /*!
    \brief ???

    */
    void DoTimeSample(int step, Teuchos::RCP<Epetra_Vector> velnp,
        Teuchos::RCP<Epetra_Vector> force, Teuchos::RCP<Epetra_Vector> phi,
        Teuchos::RCP<const DRT::DofSet> stddofset);

    /*!
    \brief Write (dump) the statistics to a file

    */
    void DoOutput(IO::DiscretizationWriter& output, int step, const bool inflow = false);

    /*!
    \brief Restart collection of statistics

    */
    void ReadRestart(IO::DiscretizationReader& reader, int step);

    /*!
    \brief Restart scatra-specific collection of statistics

    */
    void ReadRestartScaTra(IO::DiscretizationReader& scatrareader, int step);

    /*!
    \brief Provide access to scalar transport field

    */
    void AddScaTraField(Teuchos::RCP<SCATRA::ScaTraTimIntImpl> scatra_timeint);

    /*!
    \brief   Write (dump) the scatra-specific mean fields to the result file

    */
    void DoOutputForScaTra(IO::DiscretizationWriter& output, int step);

    //@}

    /*!
    \brief remote access method to general mean statistics manager

    */
    const Teuchos::RCP<TurbulenceStatisticsGeneralMean> GetTurbulenceStatisticsGeneralMean() const
    {
      return statistics_general_mean_;
    }

    bool WithScaTra() { return withscatra_; }

   private:
    /*!
    \brief Time integration independent setup called by Constructor

    */
    void Setup();

    //! time step size
    double dt_;

    //! parameters for sampling/dumping period
    //! start of sampling at step samstart
    int samstart_;
    //! stop sampling at step samstop
    int samstop_;
    //! incremental dump every dumperiod steps or standalone
    //! records (0)
    int dumperiod_;

    //! the fluid discretization
    Teuchos::RCP<DRT::Discretization> discret_;
    //! the scatra discretization
    Teuchos::RCP<DRT::Discretization> scatradis_;

    //! parameterlist of the discretization including time params,
    //! stabilization params and turbulence sublist
    Teuchos::RCP<Teuchos::ParameterList> params_;
    //! parameterlist of the scatra discretization including time params,
    //! stabilization params and turbulence sublist
    Teuchos::RCP<Teuchos::ParameterList> scatraparams_;
    //! since model parameters (Prt, Csgs, etc.) are recomputed during
    //! simulation and change their values, we need a pointer to this
    //! list and cannot use a copy stored in scatraparams_
    Teuchos::RCP<Teuchos::ParameterList> scatraextraparams_;
    Teuchos::RCP<Teuchos::ParameterList> scatratimeparams_;

    //! parameterlist specially designed for the evaluation of
    //! gausspoint statistics
    Teuchos::ParameterList eleparams_;

    //! decides whether we use an Eulerian or an ALE formulation
    bool alefluid_;

    //! all my solution vectors needed for element evaluation
    Teuchos::RCP<Epetra_Vector> myaccnp_;
    Teuchos::RCP<Epetra_Vector> myaccn_;
    Teuchos::RCP<Epetra_Vector> myaccam_;

    Teuchos::RCP<Epetra_Vector> myvelnp_;
    Teuchos::RCP<Epetra_Vector> myveln_;
    Teuchos::RCP<Epetra_Vector> myvelaf_;

    Teuchos::RCP<Epetra_Vector> myhist_;

    Teuchos::RCP<Epetra_Vector> myscaaf_;
    Teuchos::RCP<Epetra_Vector> myscaam_;

    Teuchos::RCP<Epetra_Vector> mydispnp_;
    Teuchos::RCP<Epetra_Vector> mydispn_;

    Teuchos::RCP<Epetra_Vector> mygridvelaf_;

    Teuchos::RCP<Epetra_Vector> myforce_;

    Teuchos::RCP<Epetra_Vector> myfsvelaf_;
    Teuchos::RCP<Epetra_Vector> myfsscaaf_;
    // xwall object is required for evaluating inner element planes of channel
    Teuchos::RCP<FLD::XWall> myxwall_;

    Teuchos::RCP<FLD::UTILS::StressManager> mystressmanager_;

    //! scatra result vector (defined on the scatra dofrowmap!)
    Teuchos::RCP<Epetra_Vector> myphinp_;
    Teuchos::RCP<Epetra_Vector> myphiaf_;
    Teuchos::RCP<Epetra_Vector> myphiam_;
    Teuchos::RCP<Epetra_Vector> myfsphi_;
    Teuchos::RCP<Epetra_Vector> myscatrahist_;
    Teuchos::RCP<Epetra_Vector> myphidtam_;

    //! number of dofset associated with velocity-related dofs on scatra discretization
    int myscatrandsvel_;

    //! specifies the special flow
    enum SpecialFlow
    {
      no_special_flow,
      channel_flow_of_height_2,
      loma_channel_flow_of_height_2,
      scatra_channel_flow_of_height_2,
      lid_driven_cavity,
      loma_lid_driven_cavity,
      backward_facing_step,
      loma_backward_facing_step,
      backward_facing_step2,
      square_cylinder,
      square_cylinder_nurbs,
      bubbly_channel_flow,
      rotating_circular_cylinder_nurbs,
      rotating_circular_cylinder_nurbs_scatra,
      decaying_homogeneous_isotropic_turbulence,
      forced_homogeneous_isotropic_turbulence,
      scatra_forced_homogeneous_isotropic_turbulence,
      taylor_green_vortex,
      periodic_hill,
      blood_fda_flow,
      time_averaging
    } flow_;

    //! flag to include scatra field
    bool withscatra_;

    //! toggle additional evaluations for turbulence models
    enum INPAR::FLUID::TurbModelAction turbmodel_;

    //! toggle evaluation of subgrid quantities, dissipation rates etc
    //! this is only possible for the genalpha implementation since
    //! we need a corresponding element implementation
    bool subgrid_dissipation_;

    //! toggle statistics output of additional inflow
    bool inflow_;

    //! mean values of velocity and pressure, independent of special flow
    //! averaging takes place in time, if hom. directions have been specified
    //! additionally along these lines
    bool out_mean_;

    //! name of statistics output file, despite the ending
    const std::string statistics_outfilename_;

    Teuchos::RCP<TurbulenceStatisticsGeneralMean> statistics_general_mean_;

    //! turbulence statistics for turbulent channel flow
    Teuchos::RCP<TurbulenceStatisticsCha> statistics_channel_;

    //! turbulence statistics for a rotating circular cylinder
    Teuchos::RCP<TurbulenceStatisticsCcy> statistics_ccy_;

    //! turbulence statistics for lid-driven cavity
    Teuchos::RCP<TurbulenceStatisticsLdc> statistics_ldc_;

    //! turbulence statistics for backward-facing step
    Teuchos::RCP<TurbulenceStatisticsBfs> statistics_bfs_;

    //! turbulence statistics for periodic hill
    Teuchos::RCP<TurbulenceStatisticsPh> statistics_ph_;

    //! turbulence statistics for blood fda flow
    Teuchos::RCP<TurbulenceStatisticsBfda> statistics_bfda_;

    //! turbulence statistics for square cylinder
    Teuchos::RCP<TurbulenceStatisticsSqc> statistics_sqc_;

    //! turbulence statistics for homogeneous isotropic turbulence
    Teuchos::RCP<TurbulenceStatisticsHit> statistics_hit_;

    //! turbulence statistics for Taylor-Green Vortex
    Teuchos::RCP<TurbulenceStatisticsTgv> statistics_tgv_;

  };  // end class TurbulenceStatisticManager

}  // end namespace FLD

#endif  // #ifndef TURBULENCE_STATISTIC_MANAGER_H
