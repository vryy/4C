/*-----------------------------------------------------------*/
/*! \file

\brief Stationary fluid problem with HDG discretization


\level 2

*/
/*-----------------------------------------------------------*/

#ifndef FLD_TIMINT_STAT_HDG_H
#define FLD_TIMINT_STAT_HDG_H  // TODO als fix fluid_timint_stat_hdg because it is not working


#include "fluid_timint_stat.H"


namespace FLD
{
  class TimIntStationaryHDG : public virtual TimIntStationary
  {
   public:
    /// Standard Constructor
    TimIntStationaryHDG(const Teuchos::RCP<DRT::Discretization>& actdis,
        const Teuchos::RCP<LINALG::Solver>& solver,
        const Teuchos::RCP<Teuchos::ParameterList>& params,
        const Teuchos::RCP<IO::DiscretizationWriter>& output, bool alefluid = false);

    /// Destructor
    virtual ~TimIntStationaryHDG() {}

    /*!
    \brief initialization

    */
    virtual void Init();

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
    virtual void Reset(bool completeReset = false, int numsteps = 1, int iter = -1);

    void SetOldPartOfRighthandside();

    /*!
      \brief Set custom parameters in the respective time integration class (Loma, RedModels...)

    */
    virtual void SetCustomEleParamsAssembleMatAndRHS(Teuchos::ParameterList& eleparams);

    /*!
    \brief Set states in the time integration schemes: differs between GenAlpha and the others

    */
    void SetStateTimInt();

    /*!
    \brief Call discret_->ClearState() after assembly (HDG needs to read from state vectors...)

    */
    virtual void ClearStateAssembleMatAndRHS();

    /*!
    \brief set initial flow field for analytical test problems

    */
    virtual void SetInitialFlowField(
        const INPAR::FLUID::InitialField initfield, const int startfuncno);

    /*!

    \brief parameter (fix over a time step) are set in this method.
           Therefore, these parameter are accessible in the fluid element
           and in the fluid boundary element

    */
    void SetElementTimeParameter();

    //@}

   protected:
    //! @name velocity gradient, velocity and pressure at time n+1, n, n-1
    //!  and n+alpha_F for element interior in HDG
    Teuchos::RCP<Epetra_Vector> intvelnp_;

   private:
    ///< Keep track of whether we do the first assembly because we reconstruct the local HDG
    ///< solution as part of assembly
    bool firstAssembly_;


  };  // class TimIntStationaryHDG

}  // namespace FLD

#endif  // #ifndef FLD_TIMINT_STAT_HDG_H
