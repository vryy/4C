/*----------------------------------------------------------------------*/
/*!
 * \file tsi_partitioned.H
\brief Partitioned TSI base

<pre>
Maintainer: Caroline Danowski
            danowski@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15253
</pre>
*/

/*----------------------------------------------------------------------*
 | definitions                                               dano 11/09 |
 *----------------------------------------------------------------------*/
#ifdef CCADISCRET

#ifndef TSI_PARTITIONED_H
#define TSI_PARTITIONED_H


/*----------------------------------------------------------------------*
 | headers                                                   dano 11/09 |
 *----------------------------------------------------------------------*/
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

#include <AztecOO.h>
#include <Epetra_CrsGraph.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_Import.h>
#include <Epetra_LinearProblem.h>
#include <Epetra_Map.h>
#include <Epetra_RowMatrix.h>
#include <Epetra_Vector.h>
#include <NOX.H>
#include <NOX_Epetra.H>
#include <NOX_Epetra_Interface_Required.H>
#include <NOX_Epetra_LinearSystem_AztecOO.H>

#include <vector>

#include "tsi_algorithm.H"


/*----------------------------------------------------------------------*
 | TSI: Thermo Structure Interaction                         dano 08/09 |
 *----------------------------------------------------------------------*/
//!  This is the place where we bring temperature and structure together.
namespace TSI
{

namespace UTILS
{
  class DebugWriter;
}

//! Partitioned TSI base
//!
//!  This is the framework for partitioned TSI. The actual work is done by
//!  subclasses.
//!
//!  \author u.kue
//!  \date 04/08
  class Partitioned : public Algorithm,
                      // 23.11.09 not for TSI!!!
                      public NOX::Epetra::Interface::Required
  {
  public:

    //! setup TSI algorithm
    explicit Partitioned(Epetra_Comm& comm);

    //! outer level TSI time loop
    void Timeloop(const Teuchos::RCP<NOX::Epetra::Interface::Required>& interface);

    //! compute TSI interface residual S^{-1}(F(d)) - d
    bool computeF(const Epetra_Vector &x, Epetra_Vector &F, const FillType fillFlag);

  protected:

    //! @name Transfer helpers

    Teuchos::RCP<Epetra_Vector> StructToThermo(Teuchos::RCP<Epetra_Vector> iv);
    Teuchos::RCP<Epetra_Vector> ThermoToStruct(Teuchos::RCP<Epetra_Vector> iv);

    //@}

    //! @name Operators implemented by subclasses

    //! composed TSI operator
    virtual void TSIOp(
      const Epetra_Vector &x,
      Epetra_Vector &F,
      const FillType fillFlag
      );

// 23.11.09 no interface!!!!!!!!!!
    //! interface thermo operator
    virtual Teuchos::RCP<Epetra_Vector> ThermoOp(
      Teuchos::RCP<Epetra_Vector> idisp,
      const FillType fillFlag
      );

// 23.11.09 no interface!!!!!!!!!!
    //! interface structural operator
    virtual Teuchos::RCP<Epetra_Vector> StructOp(
      Teuchos::RCP<Epetra_Vector> iforce,
      const FillType fillFlag
      );

    //@}

    //! @name Encapsulation of interface unknown
    //! default is displacement, but subclasses might change that
    virtual Teuchos::RCP<Epetra_Vector> InitialGuess();

    //@}


// 23.11.09 adapt this part to node based coupling!!!!
//  TSI do not has an interface --> therefore changing...

    //! @name Access methods for subclasses

  //  in TSI the velocity could only be needed in conjunction with material law
  //  //! interface velocity calculation given the interface displacements
  //  Teuchos::RCP<Epetra_Vector> InterfaceVelocity(
  //    Teuchos::RCP<const Epetra_Vector> idispnp
  //    ) const;

    //! current interface displacements
    //!
    //!  Extract structural displacement at t(n+1)
    Teuchos::RCP<Epetra_Vector> InterfaceDisp();

    //! current interface forces
    //!
    //!  Extract fluid force at t(n+1)
    Teuchos::RCP<Epetra_Vector> InterfaceForce();

    //@}

    //! create convergance tests
    virtual void CreateStatusTest(
      Teuchos::ParameterList& nlParams,
      Teuchos::RCP<NOX::Epetra::Group> grp,
      Teuchos::RCP<NOX::StatusTest::Combo> converged
      );

    Teuchos::RCP<UTILS::DebugWriter> MyDebugWriter() const { return debugwriter_; }

  private:

    //! create linear solver framework
    Teuchos::RCP<NOX::Epetra::LinearSystem> CreateLinearSystem(
      Teuchos::ParameterList& nlParams,
      const Teuchos::RCP<NOX::Epetra::Interface::Required>& interface,
      NOX::Epetra::Vector& noxSoln,
      Teuchos::RCP<NOX::Utils> utils
      );

    //! create convergance tests including testing framework
    Teuchos::RCP<NOX::StatusTest::Combo> CreateStatusTest(
      Teuchos::ParameterList& nlParams,
      Teuchos::RCP<NOX::Epetra::Group> grp
      );

    //! setup list with default parameters
    void SetDefaultParameters(
      const Teuchos::ParameterList& tsidyn,
      Teuchos::ParameterList& list
      );

    //! interface displacement from time step begin
    Teuchos::RCP<Epetra_Vector> idispn_;

    //! interface temperature from time step begin
    Teuchos::RCP<Epetra_Vector> itempn_;
  // ??????????
  //  //! connection of interface dofs for finite differences
  //  Teuchos::RCP<Epetra_CrsGraph> rawGraph_;


    //! counters on how many times the residuum was called in a time step
    //!
    //!  NOX knows different types of residuum calls depending on
    //!  circumstances (normal, finite difference, matrix free
    //!  jacobi). It is possible to do approximations depending on the
    //!  type.
    std::vector<int> counter_;

    //! number of residuum calculations per nonlinear solve in one time step
    std::vector<int> linsolvcount_;

    //! print parameters and stuff
    //!
    //!  \warning This variable is only valid during when the time loop runs.
    Teuchos::RCP<NOX::Utils> utils_;

    //! parameters handed in to NOX
    Teuchos::ParameterList noxparameterlist_;

  protected:
// 23.11.09 needed in TSI???
    //! matrix free itemax
    int mfresitemax_;
  private:
// 23.11.09 in TSI all nodes are matching nodes!!!
    //! nodes at the fluid-structure interface match
    bool matchingnodes_;

    //! special debugging output
    Teuchos::RCP<UTILS::DebugWriter> debugwriter_;

  };  // class Partitioned

}  // namespace::TSI

/*----------------------------------------------------------------------*/
#endif // TSI_PARTITIONED_H
#endif // CCADISCRET
