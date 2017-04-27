/*----------------------------------------------------------------------*/
/*!
\file fsi_algorithm.H

\brief Basis of all FSI algorithms

\maintainer Andreas Rauch

\level 1
*/
/*----------------------------------------------------------------------*/



#ifndef FSI_ALGORITHM_H
#define FSI_ALGORITHM_H

#include <Epetra_Comm.h>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

#include "../drt_adapter/ad_fld_fluid_ale.H"
#include "../drt_adapter/adapter_algorithmbase.H"


/// forward declarations
namespace ADAPTER
{
  class FSIStructureWrapper;
  class StructureBaseAlgorithmNew;
}


namespace FSI
{

  /// FSI algorithm base
  /*!

    Base class of FSI algorithms with generalized fluid field. There can (and
    will) be different subclasses that implement different coupling schemes.

    \note The generalized fluid field hides any ale or xfem handling of the
    variable fluid domain. This is the appropriate base class if direct access
    to ale or xfem is not required. If a coupling algorithm needs an ale
    field, MonolithicBase is the better choice for a base class.

    \warning The order of calling the three BaseAlgorithm-constructors (that
    is the order in which we list the base classes) is important here! In the
    constructors control file entries are written. And these entries define
    the order in which the filters handle the Discretizations, which in turn
    defines the dof number ordering of the Discretizations... Don't get
    confused. Just always list structure, fluid, ale. In that order.

    \author u.kue
    \date 02/08
   */
  class Algorithm : public ADAPTER::AlgorithmBase
  {
  public:

    /// create using a Epetra_Comm
    explicit Algorithm(const Epetra_Comm& comm);

    /// virtual destructor to support polymorph destruction
    virtual ~Algorithm();

    /// setup this object
    virtual void Setup();

    /// access to structure field
    const Teuchos::RCP<ADAPTER::FSIStructureWrapper>& StructureField(){return structure_;}

    /// access to fluid field
    const Teuchos::RCP<ADAPTER::FluidMovingBoundary>& MBFluidField(){return fluid_;}

    /// read restart data
    virtual void ReadRestart(int step);

  protected:

    //! @name Time loop building blocks

    /// start a new time step
    virtual void PrepareTimeStep();

    /// take current results for converged and save for next time step
    virtual void Update();

    /// calculate stresses, strains, energies
    virtual void PrepareOutput();

    /// write output
    virtual void Output();

    //@}

    //! @name Transfer helpers
    virtual Teuchos::RCP<Epetra_Vector> StructToFluid(Teuchos::RCP<Epetra_Vector> iv) const;
    virtual Teuchos::RCP<Epetra_Vector> FluidToStruct(Teuchos::RCP<Epetra_Vector> iv) const;
    virtual Teuchos::RCP<Epetra_Vector> StructToFluid(Teuchos::RCP<const Epetra_Vector> iv) const;
    virtual Teuchos::RCP<Epetra_Vector> FluidToStruct(Teuchos::RCP<const Epetra_Vector> iv) const;
    //@}

    /// return the structure fluid coupling object
    ADAPTER::Coupling& StructureFluidCoupling();

    /// return const version of structure fluid coupling object
    const ADAPTER::Coupling& StructureFluidCoupling() const;

  protected:

    /// underlying structure of the FSI problem
    Teuchos::RCP< ::ADAPTER::FSIStructureWrapper>    structure_;

    /// underlying fluid of the FSI problem
    Teuchos::RCP< ::ADAPTER::FluidMovingBoundary>    fluid_;

    /// RCP pointer to the base algorithm of the structure
    Teuchos::RCP<ADAPTER::StructureBaseAlgorithmNew> adapterbase_ptr_;

    /// use deprecated old structural time integration. todo Has to be removed !
    bool use_old_structure_;

  private:

    /// coupling of structure and fluid at the interface
    Teuchos::RCP<ADAPTER::Coupling>               coupsf_;
  };
}

#endif
