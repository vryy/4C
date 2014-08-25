/*----------------------------------------------------------------------*/
/*!
\file fsi_xfem_algorithm.H

\brief Basis of monolithic XFSI algorithm that performs a coupling between the
       structural field equation and XFEM fluid field equations

<pre>
Maintainer: Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
*/


#ifndef FSI_ALGORITHMXFEM_H
#define FSI_ALGORITHMXFEM_H


#include "../drt_adapter/adapter_algorithmbase.H"
#include <Epetra_Vector.h>



namespace DRT
{
  class Discretization;
}

namespace ADAPTER
{
  class Structure;
  class FSIStructureWrapper;
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
    explicit AlgorithmXFEM(const Epetra_Comm& comm, const Teuchos::ParameterList& timeparams);

    //! virtual destructor to support polymorph destruction
    virtual ~AlgorithmXFEM();

    //! outer level time loop (to be implemented by deriving classes)
    virtual void Timeloop() = 0;

    /// initialise XFSI system
    virtual void SetupSystem() = 0;

    //! read restart data
    virtual void ReadRestart(
      int step  //!< step number where the calculation is continued
      ) = 0;

    //--------------------------------------------------------------------------//
    //! @name Access to single fields

    //! access to structural field
    const Teuchos::RCP<ADAPTER::FSIStructureWrapper>& StructureField(){return structure_;}

    //! access to fluid field
    const Teuchos::RCP< ::ADAPTER::Fluid>& FluidField(){return fluid_;}

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

    //! underlying structure of the FSI problem
    Teuchos::RCP< ::ADAPTER::FSIStructureWrapper> structure_;

    //! underlying fluid of the FSI problem
    Teuchos::RCP< ::ADAPTER::Fluid> fluid_;

    //@}

  private:

  };  // Algorithm
}  // namespace FSI


/*----------------------------------------------------------------------*/
#endif  // FSI_ALGORITHMXFEM_H
