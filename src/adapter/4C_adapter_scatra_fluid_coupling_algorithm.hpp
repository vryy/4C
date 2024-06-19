/*----------------------------------------------------------------------*/
/*! \file

\brief Basis of all algorithms that perform a coupling between Navier-Stokes
       and (active or passive) scalar transport equations

\level 1


*/
/*----------------------------------------------------------------------*/


#ifndef FOUR_C_ADAPTER_SCATRA_FLUID_COUPLING_ALGORITHM_HPP
#define FOUR_C_ADAPTER_SCATRA_FLUID_COUPLING_ALGORITHM_HPP

#include "4C_config.hpp"

#include "4C_adapter_algorithmbase.hpp"
#include "4C_adapter_fld_base_algorithm.hpp"
#include "4C_adapter_scatra_base_algorithm.hpp"
#include "4C_coupling_adapter_volmortar.hpp"
#include "4C_inpar_scatra.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Adapter
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
                                       public Adapter::FluidBaseAlgorithm,
                                       public Adapter::ScaTraBaseAlgorithm
  {
   public:
    /// constructor using a Epetra_Comm
    ScaTraFluidCouplingAlgorithm(const Epetra_Comm& comm,  ///< communicator
        const Teuchos::ParameterList& prbdyn,              ///< problem-specific parameters
        bool isale,                        ///< do we need an ALE formulation of the fields?
        const std::string scatra_disname,  ///< scatra discretization name
        const Teuchos::ParameterList& solverparams);

    /// setup this class
    void setup() override;

    /// init this class
    void init() override;

    /// outer level time loop (to be implemented by deriving classes)
    virtual void TimeLoop() = 0;

    /// read restart data
    void read_restart(int step  ///< step number where the calculation is continued
        ) override;

   protected:
    /// perform algorithm specific initialization stuff
    virtual void do_algorithm_specific_init(){};

    /// provide access to algorithm parameters
    virtual const Teuchos::ParameterList& algo_parameters() { return params_; }

    /// interpolate fluid quantity to a scatra one (e.g. via volmortar)
    Teuchos::RCP<const Epetra_Vector> fluid_to_scatra(
        const Teuchos::RCP<const Epetra_Vector> fluidvector) const;

    /// interpolate scatra quantity to a fluid one (e.g. via volmortar)
    Teuchos::RCP<const Epetra_Vector> scatra_to_fluid(
        const Teuchos::RCP<const Epetra_Vector> scatravector) const;

   private:
    /// setup adapters for transport on boundary if necessary
    void setup_field_coupling(const std::string fluid_disname, const std::string scatra_disname);

    /// flag for type of field coupling (i.e. matching or volmortar)
    Inpar::ScaTra::FieldCoupling fieldcoupling_;

    //! volume coupling (using mortar) adapter
    Teuchos::RCP<Core::Adapter::MortarVolCoupl> volcoupl_fluidscatra_;

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
    //! returns true if setup() was called and is still valid
    bool is_setup() { return issetup_; };

    //! returns true if init(..) was called and is still valid
    bool is_init() { return isinit_; };

    //! check if \ref setup() was called
    void check_is_setup()
    {
      if (not is_setup()) FOUR_C_THROW("setup() was not called.");
    };

    //! check if \ref init() was called
    void check_is_init()
    {
      if (not is_init()) FOUR_C_THROW("init(...) was not called.");
    };

   public:
    //! set flag true after setup or false if setup became invalid
    void set_is_setup(bool trueorfalse) { issetup_ = trueorfalse; };

    //! set flag true after init or false if init became invalid
    void set_is_init(bool trueorfalse) { isinit_ = trueorfalse; };
  };

}  // namespace Adapter

FOUR_C_NAMESPACE_CLOSE

#endif
