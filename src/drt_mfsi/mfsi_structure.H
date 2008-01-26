/*----------------------------------------------------------------------*/
/*!
\file mfsi_structure.H

\brief Structure field adapter

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>
*/
/*----------------------------------------------------------------------*/

#ifdef CCADISCRET

#ifndef MFSI_STRUCTURE_H
#define MFSI_STRUCTURE_H

#include <Teuchos_RefCountPtr.hpp>
#include <Teuchos_ParameterList.hpp>

#include <Epetra_Vector.h>
#include <Epetra_Map.h>
#include <Epetra_Operator.h>

#include "../drt_structure/strugenalpha.H"

#include "../drt_fsi/fsi_utils_mapextractor.H"


namespace MFSI {

  /// general structural field interface for monolithic FSI
  /*!

  The point is to keep monolithic FSI as far apart from our field
  solvers as possible. Each structure field solver we want to use should
  get its own subclass of this. The monolithic FSI algorithm should be
  able to extract all the information from the structure field it needs
  using this interface.

  \warning A cleanup is needed badly.

  \author u.kue
  \date 11/07
  */
  class Structure
  {
  public:

    virtual ~Structure();

    virtual Teuchos::RCP<const Epetra_Vector> InitialGuess() = 0;
    virtual Teuchos::RCP<Epetra_Vector> RHS() = 0;
    virtual Teuchos::RCP<Epetra_Vector> Dispnp() = 0;
    virtual Teuchos::RCP<Epetra_Vector> Disp() = 0;
    virtual Teuchos::RCP<Epetra_Vector> Dispm() = 0;
    virtual Teuchos::RCP<const Epetra_Map> DofRowMap() = 0;
    virtual Teuchos::RCP<Epetra_Operator> SysMat() = 0;
    virtual Teuchos::RCP<DRT::Discretization> Discretization() = 0;
    virtual double DispIncrFactor() = 0;

    //virtual Teuchos::RCP<Epetra_Vector> InterfaceDispn() = 0;

    virtual Teuchos::RCP<Epetra_Vector> FluidCondRHS() = 0;
    virtual Teuchos::RCP<Epetra_Vector> MeshCondRHS() = 0;

    virtual void PrepareTimeStep() = 0;
    virtual void Evaluate(Teuchos::RCP<const Epetra_Vector> disp) = 0;
    virtual void Update() = 0;
    virtual void Output() = 0;

    virtual void SetInterfaceMap(Teuchos::RCP<Epetra_Map> im) = 0;
    virtual Teuchos::RCP<Epetra_Map> InterfaceMap() = 0;
    virtual Teuchos::RCP<Epetra_Map> InnerDisplacementRowMap() = 0;

    /// communication object at the interface
    virtual const FSI::UTILS::MapExtractor& Interface() const = 0;

    virtual const Epetra_Map& DomainMap() = 0;
  };


  /// adapter to generalized alpha algorithm
  class StructureAdapter : public Structure
  {
  public:

    StructureAdapter(Teuchos::RCP<Teuchos::ParameterList> params,
                     Teuchos::RCP<DRT::Discretization> dis,
                     Teuchos::RCP<LINALG::Solver> solver,
                     Teuchos::RCP<IO::DiscretizationWriter> output);

    Teuchos::RCP<const Epetra_Vector> InitialGuess();
    Teuchos::RCP<Epetra_Vector> RHS();
    Teuchos::RCP<Epetra_Vector> Dispnp();
    Teuchos::RCP<Epetra_Vector> Disp();
    Teuchos::RCP<Epetra_Vector> Dispm();
    Teuchos::RCP<const Epetra_Map> DofRowMap();
    Teuchos::RCP<Epetra_Operator> SysMat();
    Teuchos::RCP<DRT::Discretization> Discretization();
    double DispIncrFactor();

    virtual Teuchos::RCP<Epetra_Vector> FluidCondRHS();
    virtual Teuchos::RCP<Epetra_Vector> MeshCondRHS();

    //Teuchos::RCP<Epetra_Vector> InterfaceDispn();

    void PrepareTimeStep();
    void Evaluate(Teuchos::RCP<const Epetra_Vector> disp);
    void Update();
    void Output();

    void SetInterfaceMap(Teuchos::RCP<Epetra_Map> im);
    Teuchos::RCP<Epetra_Map> InterfaceMap() { return interface_.CondDofMap(); }
    Teuchos::RCP<Epetra_Map> InnerDisplacementRowMap() { return interface_.OtherDofMap(); }

    /// communication object at the interface
    virtual const FSI::UTILS::MapExtractor& Interface() const { return interface_; }

    virtual const Epetra_Map& DomainMap();

  private:

    /// the actual structure algorithm
    StruGenAlpha structure_;

    /// the interface map setup for interface <-> full translation
    FSI::UTILS::MapExtractor interface_;

    //! @name local copies of input parameters

    Teuchos::RCP<DRT::Discretization> dis_;
    Teuchos::RCP<Teuchos::ParameterList> params_;
    Teuchos::RCP<LINALG::Solver> solver_;
    Teuchos::RCP<IO::DiscretizationWriter> output_;

    //@}
  };
}

#endif
#endif
