/*----------------------------------------------------------------------*/
/*! \file

\brief Basis of monolithic XFSI algorithm that performs a coupling between the
       structural field equation and XFEM fluid field equations

\level 2

*/


#ifndef FSI_ALGORITHMXFEM_H
#define FSI_ALGORITHMXFEM_H


#include "../drt_adapter/adapter_algorithmbase.H"
#include <Epetra_Vector.h>
#include "../drt_adapter/ad_field_wrapper.H"

namespace DRT
{
  class Discretization;
}

namespace ADAPTER
{
  class StructurePoroWrapper;
  class AleFpsiWrapper;
}  // namespace ADAPTER

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
  class AlgorithmXFEM : public ADAPTER::AlgorithmBase
  {
   public:
    //! create using a Epetra_Comm
    AlgorithmXFEM(const Epetra_Comm& comm, const Teuchos::ParameterList& timeparams,
        const ADAPTER::FieldWrapper::Fieldtype type);

    //! virtual destructor to support polymorph destruction
    virtual ~AlgorithmXFEM();

    //! setup
    virtual void Setup();

    //! outer level time loop (to be implemented by deriving classes)
    virtual void Timeloop() = 0;

    /// initialise XFSI system
    virtual void SetupSystem() = 0;

    //! read restart data
    virtual void ReadRestart(int step  //!< step number where the calculation is continued
        ) = 0;

    //--------------------------------------------------------------------------//
    //! @name Access to single fields

    //! access to structural & poro field
    const Teuchos::RCP<ADAPTER::StructurePoroWrapper>& StructurePoro() { return structureporo_; }

    //! access to fluid field
    const Teuchos::RCP<::FLD::XFluid>& FluidField() { return fluid_; }

    //! access to ale field
    const Teuchos::RCP<::ADAPTER::AleFpsiWrapper>& AleField() { return ale_; }

    //! is an monolithic ale computations
    bool HaveAle() { return (AleField() != Teuchos::null); }

    //! number of physical fields to solve involved
    int NumFields() { return num_fields_; }

    //@}

   protected:
    //--------------------------------------------------------------------------//
    //! @name Time loop building blocks

    //! start a new time step
    virtual void PrepareTimeStep() = 0;

    //! calculate stresses, strains, energies
    virtual void PrepareOutput();

    //! take current results for converged and save for next time step
    virtual void Update();

    //! write output
    virtual void Output() = 0;

    //@}

    //--------------------------------------------------------------------------//
    //! @name Underlying fields

    //! underlying structure / poro of the FSI/FPSI problem
    Teuchos::RCP<::ADAPTER::StructurePoroWrapper> structureporo_;

    //! underlying fluid of the FSI problem
    Teuchos::RCP<::FLD::XFluid> fluid_;

    // underlying ale of the FSI problem
    Teuchos::RCP<ADAPTER::AleFpsiWrapper> ale_;

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
#endif  // FSI_ALGORITHMXFEM_H
