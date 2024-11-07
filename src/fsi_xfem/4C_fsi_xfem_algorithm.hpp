// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_FSI_XFEM_ALGORITHM_HPP
#define FOUR_C_FSI_XFEM_ALGORITHM_HPP


#include "4C_config.hpp"

#include "4C_adapter_algorithmbase.hpp"
#include "4C_adapter_field_wrapper.hpp"
#include "4C_linalg_vector.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE

namespace Adapter
{
  class StructurePoroWrapper;
  class AleFpsiWrapper;
}  // namespace Adapter

namespace FLD
{
  class XFluid;
}

/*----------------------------------------------------------------------*
 |                                                         schott 08/14 |
 *----------------------------------------------------------------------*/
//! FSI: Fluid-Structure Interaction
namespace FSI
{
  //! XFSI algorithm base
  //!
  //!  Base class of XFSI algorithms. Derives from AlgorithmBase.
  //!
  //!  \author schott
  //!  \date 08/14
  class AlgorithmXFEM : public Adapter::AlgorithmBase
  {
   public:
    //! create using a Epetra_Comm
    AlgorithmXFEM(const Epetra_Comm& comm, const Teuchos::ParameterList& timeparams,
        const Adapter::FieldWrapper::Fieldtype type);


    //! setup
    virtual void setup();

    //! outer level time loop (to be implemented by deriving classes)
    virtual void timeloop() = 0;

    /// initialise XFSI system
    virtual void setup_system() = 0;

    //! read restart data
    void read_restart(int step  //!< step number where the calculation is continued
        ) override = 0;

    //--------------------------------------------------------------------------//
    //! @name Access to single fields

    //! access to structural & poro field
    const std::shared_ptr<Adapter::StructurePoroWrapper>& structure_poro()
    {
      return structureporo_;
    }

    //! access to fluid field
    const std::shared_ptr<FLD::XFluid>& fluid_field() { return fluid_; }

    //! access to ale field
    const std::shared_ptr<Adapter::AleFpsiWrapper>& ale_field() { return ale_; }

    //! is an monolithic ale computations
    bool have_ale() { return (ale_field() != nullptr); }

    //! number of physical fields to solve involved
    int num_fields() { return num_fields_; }

    //@}

   protected:
    //--------------------------------------------------------------------------//
    //! @name Time loop building blocks

    //! start a new time step
    void prepare_time_step() override = 0;

    //! calculate stresses, strains, energies
    virtual void prepare_output(bool force_prepare);

    //! take current results for converged and save for next time step
    void update() override;

    //! write output
    void output() override = 0;

    //@}

    //--------------------------------------------------------------------------//
    //! @name Underlying fields

    //! underlying structure / poro of the FSI/FPSI problem
    std::shared_ptr<Adapter::StructurePoroWrapper> structureporo_;

    //! underlying fluid of the FSI problem
    std::shared_ptr<FLD::XFluid> fluid_;

    // underlying ale of the FSI problem
    std::shared_ptr<Adapter::AleFpsiWrapper> ale_;

    //@}

    //--------------------------------------------------------------------------//
    //! @name block ids of the monolithic system
    int num_fields_;
    int structp_block_;
    int fluid_block_;
    int fluidp_block_;
    int ale_i_block_;

    //@}

   private:
  };  // Algorithm
}  // namespace FSI


/*----------------------------------------------------------------------*/
FOUR_C_NAMESPACE_CLOSE

#endif
