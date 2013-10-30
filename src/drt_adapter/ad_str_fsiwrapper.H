/*----------------------------------------------------------------------*/
/*!
\file ad_str_fsiwrapper.H

\brief Structural adapter for FSI problems containing the interface
       and methods dependent on the interface

<pre>
Maintainer: Georg Hammerl
            hammerl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15237
</pre>
*/

#ifndef FSI_STRUCTURE_WRAPPER_H
#define FSI_STRUCTURE_WRAPPER_H

#include "ad_str_wrapper.H"

// forward declarations
namespace STR
{
  namespace AUX
  {
    class MapExtractor;
  }
}


namespace ADAPTER
{

  class FSIStructureWrapper : public StructureWrapper
  {
  public:

    // constructor
    explicit FSIStructureWrapper(Teuchos::RCP<Structure> structure);


    /// communication object at the interface
    virtual Teuchos::RCP<const STR::AUX::MapExtractor> Interface() const
    {
      return interface_;
    }

    //! switch structure field to block matrix in fsi simulations
    virtual void UseBlockMatrix();

    /// linear structure solve with just a interface load
    ///
    /// The very special solve done in steepest descent relaxation
    /// calculation (and matrix free Newton Krylov).
    ///
    /// \note Can only be called after a valid structural solve.
    virtual Teuchos::RCP<Epetra_Vector> RelaxationSolve(Teuchos::RCP<Epetra_Vector> iforce);

    /// @name Extract interface values

    /// extract interface displacements at \f$t_{n}\f$
    virtual Teuchos::RCP<Epetra_Vector> ExtractInterfaceDispn();

    /// extract interface displacements at \f$t_{n+1}\f$
    virtual Teuchos::RCP<Epetra_Vector> ExtractInterfaceDispnp();

    /// Predictor for interface displacements
    virtual Teuchos::RCP<Epetra_Vector> PredictInterfaceDispnp();

    /// @name Apply interface forces

    /// apply interface forces to structural solver
    ///
    /// This prepares a new solve of the structural field within one time
    /// step. The middle values are newly created.
    ///
    /// \note This is not yet the most efficient implementation.
    virtual void ApplyInterfaceForces(Teuchos::RCP<Epetra_Vector> iforce);

		//! rebuild FSI interface from structure side
    virtual void RebuildInterface();


  private:

    /// the interface map setup for interface <-> full translation
    Teuchos::RCP<STR::AUX::MapExtractor> interface_;

    /// predictor type
    int predictor_;

  };
}

#endif
