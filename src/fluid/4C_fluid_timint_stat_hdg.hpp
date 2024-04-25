/*-----------------------------------------------------------*/
/*! \file

\brief Stationary fluid problem with HDG discretization


\level 2

*/
/*-----------------------------------------------------------*/

#ifndef FOUR_C_FLUID_TIMINT_STAT_HDG_HPP
#define FOUR_C_FLUID_TIMINT_STAT_HDG_HPP
// TODO als fix
// fluid_timint_stat_hdg
// because it is not working


#include "4C_config.hpp"

#include "4C_fluid_timint_stat.hpp"

FOUR_C_NAMESPACE_OPEN


namespace FLD
{
  class TimIntStationaryHDG : public virtual TimIntStationary
  {
   public:
    /// Standard Constructor
    TimIntStationaryHDG(const Teuchos::RCP<DRT::Discretization>& actdis,
        const Teuchos::RCP<CORE::LINALG::Solver>& solver,
        const Teuchos::RCP<Teuchos::ParameterList>& params,
        const Teuchos::RCP<IO::DiscretizationWriter>& output, bool alefluid = false);

    /*!
    \brief initialization

    */
    void Init() override;

    /*!
    \brief Set the part of the righthandside belonging to the last
           timestep

       Stationary:

                     mom: hist_ = 0.0
                    (con: hist_ = 0.0)


    */
    /*!
      \brief Reset state vectors
       */
    void Reset(bool completeReset = false, int numsteps = 1, int iter = -1) override;

    void SetOldPartOfRighthandside() override;

    /*!
      \brief Set custom parameters in the respective time integration class (Loma, RedModels...)

    */
    void SetCustomEleParamsAssembleMatAndRHS(Teuchos::ParameterList& eleparams) override;

    /*!
    \brief Set states in the time integration schemes: differs between GenAlpha and the others

    */
    void SetStateTimInt() override;

    /*!
    \brief Call discret_->ClearState() after assembly (HDG needs to read from state vectors...)

    */
    void ClearStateAssembleMatAndRHS() override;

    /*!
    \brief set initial flow field for analytical test problems

    */
    void SetInitialFlowField(
        const INPAR::FLUID::InitialField initfield, const int startfuncno) override;

    /*!

    \brief parameter (fix over a time step) are set in this method.
           Therefore, these parameter are accessible in the fluid element
           and in the fluid boundary element

    */
    void SetElementTimeParameter() override;

    //@}

   protected:
    //! @name velocity gradient, velocity and pressure at time n+1, n, n-1
    //!  and n+alpha_F for element interior in HDG
    Teuchos::RCP<Epetra_Vector> intvelnp_;

   private:
    ///< Keep track of whether we do the first assembly because we reconstruct the local HDG
    ///< solution as part of assembly
    bool first_assembly_;


  };  // class TimIntStationaryHDG

}  // namespace FLD

FOUR_C_NAMESPACE_CLOSE

#endif
