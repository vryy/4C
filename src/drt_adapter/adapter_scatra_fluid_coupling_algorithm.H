/*----------------------------------------------------------------------*/
/*!
\file adapter_scatra_fluid_coupling_algorithm.H

\brief Basis of all algorithms that perform a coupling between Navier-Stokes
       and (active or passive) scalar transport equations

\level 1
<pre>
\maintainer Rui Fang
            fang@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089 - 289-15251
</pre>
*/
/*----------------------------------------------------------------------*/


#ifndef ADAPTER_SCATRA_FLUID_COUPLING_ALGORITHM_H
#define ADAPTER_SCATRA_FLUID_COUPLING_ALGORITHM_H

#include "adapter_algorithmbase.H"
#include "../drt_adapter/ad_fld_base_algorithm.H"
#include "../drt_adapter/adapter_scatra_base_algorithm.H"
#include "../drt_inpar/inpar_scatra.H"

namespace ADAPTER
{
  class MortarVolCoupl;

  /// basis coupling algorithm for scalar transport with fluid velocity field
  /*!

    Base class for scalar transport problems coupled to Navier-Stokes velocity field.
    Derives from FluidBaseAlgorithm and ScaTraBaseAlgorithm.
    There are different subclasses that implement different coupling schemes
    (one-way coupling for the transport of passive scalars and
     fully coupled schemes for other multiphysics applications like electrochemistry).

    \author gjb
    \date 07/08
   */
  class ScaTraFluidCouplingAlgorithm : public AlgorithmBase,
                                       public ADAPTER::FluidBaseAlgorithm,
                                       public ADAPTER::ScaTraBaseAlgorithm
  {
   public:
    /// constructor using a Epetra_Comm
    ScaTraFluidCouplingAlgorithm(const Epetra_Comm& comm,  ///< communicator
        const Teuchos::ParameterList& prbdyn,              ///< problem-specific parameters
        bool isale,                        ///< do we need an ALE formulation of the fields?
        const std::string scatra_disname,  ///< scatra discretization name
        const Teuchos::ParameterList& solverparams);

    /// virtual destructor to support polymorph destruction
    virtual ~ScaTraFluidCouplingAlgorithm(){};

    /// setup this class
    virtual void Setup();

    /// init this class
    virtual void Init(const Teuchos::ParameterList& prbdyn,  ///< parameter list for global problem
        const Teuchos::ParameterList&
            scatradyn,  ///< parameter list for scalar transport subproblem
        const Teuchos::ParameterList& solverparams,  ///< parameter list for scalar transport solver
        const std::string& disname = "scatra",       ///< name of scalar transport discretization
        const bool isale = false                     ///< ALE flag
    );

    /// outer level time loop (to be implemented by deriving classes)
    virtual void TimeLoop() = 0;

    /// read restart data
    virtual void ReadRestart(int step  ///< step number where the calculation is continued
    );

   protected:
    /// perform algorithm specific initialization stuff
    virtual void DoAlgorithmSpecificInit(){};

    /// provide access to algorithm parameters
    virtual const Teuchos::ParameterList& AlgoParameters() { return params_; }

    /// interpolate fluid quantity to a scatra one (e.g. via volmortar)
    const Teuchos::RCP<const Epetra_Vector> FluidToScatra(
        const Teuchos::RCP<const Epetra_Vector> fluidvector) const;

    /// interpolate scatra quantity to a fluid one (e.g. via volmortar)
    const Teuchos::RCP<const Epetra_Vector> ScatraToFluid(
        const Teuchos::RCP<const Epetra_Vector> scatravector) const;

   private:
    /// setup adapters for transport on boundary if necessary
    void SetupFieldCoupling(const std::string fluid_disname, const std::string scatra_disname);

    /// flag for type of field coupling (i.e. matching or volmortar)
    INPAR::SCATRA::FieldCoupling fieldcoupling_;

    //! volume coupling (using mortar) adapter
    Teuchos::RCP<::ADAPTER::MortarVolCoupl> volcoupl_fluidscatra_;

    /// problem-specific parameter list
    const Teuchos::ParameterList& params_;

    /// name of scatra discretization
    const std::string scatra_disname_;

   private:
    //! flag indicating if class is setup
    bool issetup_;

    //! flag indicating if class is initialized
    bool isinit_;

   protected:
    //! returns true if Setup() was called and is still valid
    bool IsSetup() { return issetup_; };

    //! returns true if Init(..) was called and is still valid
    bool IsInit() { return isinit_; };

    //! check if \ref Setup() was called
    void CheckIsSetup()
    {
      if (not IsSetup()) dserror("Setup() was not called.");
    };

    //! check if \ref Init() was called
    void CheckIsInit()
    {
      if (not IsInit()) dserror("Init(...) was not called.");
    };

   public:
    //! set flag true after setup or false if setup became invalid
    void SetIsSetup(bool trueorfalse) { issetup_ = trueorfalse; };

    //! set flag true after init or false if init became invalid
    void SetIsInit(bool trueorfalse) { isinit_ = trueorfalse; };
  };

}  // namespace ADAPTER

#endif  // ADAPTER_SCATRA_FLUID_COUPLING_ALGORITHM_H
