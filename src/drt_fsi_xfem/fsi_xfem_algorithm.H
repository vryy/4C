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


/*----------------------------------------------------------------------*
 | definitions                                             schott 08/13 |
 *----------------------------------------------------------------------*/
#ifndef FSI_ALGORITHMXFEM_H
#define FSI_ALGORITHMXFEM_H


/*----------------------------------------------------------------------*
 | headers                                                 schott 08/13 |
 *----------------------------------------------------------------------*/
#include "../drt_adapter/ad_str_structure.H"
#include "../drt_adapter/adapter_algorithmbase.H"
#include "../drt_adapter/ad_fld_base_algorithm.H"


#include <Epetra_Vector.h>


/*----------------------------------------------------------------------*
 | forward declarations                                    schott 08/13 |
 *----------------------------------------------------------------------*/
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
 |                                                         schott 08/13 |
 *----------------------------------------------------------------------*/
//! TSI: Thermo-Structure Interaction
namespace FSI
{
  //! TSI algorithm base
  //!
  //!
  //!  Base class of TSI algorithms. Derives from StructureBaseAlgorithm and
  //!  ThermoBaseAlgorithm with temperature field.
  //!  There can (and will) be different subclasses that implement different
  //!  coupling schemes.
  //!
  //!  \warning The order of calling the two BaseAlgorithm-constructors (that
  //!  is the order in which we list the base classes) is important here! In the
  //!  constructors control file entries are written. And these entries define
  //!  the order in which the filters handle the Discretizations, which in turn
  //!  defines the dof number ordering of the Discretizations... Don't get
  //!  confused. Just always list structure, thermo. In that order.
  //!
  //!  \author u.kue
  //!  \date 02/08
  class AlgorithmXFEM : public ADAPTER::AlgorithmBase
  {

  public:

    //! create using a Epetra_Comm
    explicit AlgorithmXFEM(const Epetra_Comm& comm, const Teuchos::ParameterList& timeparams);

    //! virtual destructor to support polymorph destruction
    virtual ~AlgorithmXFEM();

    //! outer level time loop (to be implemented by deriving classes)
    virtual void Timeloop() = 0;

    /// initialise TSI system
    virtual void SetupSystem() = 0;

    //! read restart data
    virtual void ReadRestart(
      int step  //!< step number where the calculation is continued
      ) = 0;

    //! access to structural field
    const Teuchos::RCP<ADAPTER::FSIStructureWrapper>& StructureField(){return structure_;}

    //! access to fluid field
    const Teuchos::RCP< ::ADAPTER::Fluid>& FluidField(){return fluid_;}

  protected:

    //! @name Time loop building blocks

    //! start a new time step
    virtual void PrepareTimeStep() = 0;

    //! calculate stresses, strains, energies
//    void PrepareOutput()
//    {
//      StructureField()->PrepareOutput();
//    }

    //! take current results for converged and save for next time step
    virtual void Update();

    //! write output
    virtual void Output();

    //! communicate displacement vector to thermal field to enable their
    //! visualisation on the deformed body
    void OutputDeformationInThr(
      Teuchos::RCP<const Epetra_Vector> dispnp,
      Teuchos::RCP<DRT::Discretization> structdis
      );

    //@}


    //! @name Access methods

    //! velocity calculation given the displacements (like in FSI)
    Teuchos::RCP<Epetra_Vector> CalcVelocity(
      Teuchos::RCP<const Epetra_Vector> dispnp
      );

    //! displacements at time n+1 for thermal output
    Teuchos::RCP<Epetra_MultiVector> dispnp_;

    //@}


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
