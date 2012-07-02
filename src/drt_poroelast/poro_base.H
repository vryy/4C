/*----------------------------------------------------------------------*/
/*!
 \file poro_base.H

 \brief  Basis of all porous media algorithms

 <pre>
   Maintainer: Anh-Tu Vuong
               vuong@lnm.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-15264
 </pre>
 *-----------------------------------------------------------------------*/
#ifndef PORO_BASE_H_
#define PORO_BASE_H_

/*----------------------------------------------------------------------*
 | headers                                                              |
 *----------------------------------------------------------------------*/
#include <Epetra_Time.h>
#include <Epetra_Vector.h>

#include "../drt_adapter/adapter_algorithmbase.H"
#include "../linalg/linalg_mapextractor.H"

/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
namespace DRT
{
  class Condition;
}

namespace LINALG
{
  class SparseMatrix;
  class SparseOperator;
  class MapExtractor;

  class BlockSparseMatrixBase;
  class Solver;
}

namespace ADAPTER
{
  class Fluid;
  class Coupling;
  class FSIStructureWrapper;
}

namespace POROELAST
{
  /// monolithic Poroelasticity algorithm base

  class PoroBase: public ADAPTER::AlgorithmBase
  {
    public:

      //! create using a Epetra_Comm
      explicit PoroBase(const Epetra_Comm& comm, const Teuchos::ParameterList& timeparams);

      //! virtual destructor to support polymorph destruction
      virtual ~PoroBase();

      /// read restart data
      virtual void ReadRestart(int restart);

      //! outer level time loop (to be implemented by deriving classes)
      virtual void TimeLoop() = 0;

      /// initialise system
      virtual void SetupSystem() = 0;

      //! @name Time loop building blocks

      //! start a new time step
      virtual void PrepareTimeStep();

      //! take current results for converged and save for next time step
      virtual void Update();

      //! calculate stresses, strains, energies
      virtual void PrepareOutput();

      //! perform result tests
      virtual void TestResults(const Epetra_Comm& comm);

      //! access to structural field
      const Teuchos::RCP<ADAPTER::FSIStructureWrapper>& StructureField(){return structure_;}

      //! access to fluid field
      ::ADAPTER::Fluid& FluidField(){return *fluid_;}

    protected:

      //! @name Underlying fields

      //! underlying structure of the FSI problem
      Teuchos::RCP< ::ADAPTER::FSIStructureWrapper>    structure_;

      //! underlying fluid of the FSI problem
      Teuchos::RCP< ::ADAPTER::Fluid>        fluid_;

      //@}

      //! @name Transfer helpers

      //! field transform
      virtual Teuchos::RCP<Epetra_Vector> StructureToFluidField(Teuchos::RCP<
          const Epetra_Vector> iv) const;

      //! field transform
      virtual Teuchos::RCP<Epetra_Vector> FluidToStructureField(Teuchos::RCP<
          const Epetra_Vector> iv) const;

      virtual Teuchos::RCP<Epetra_Vector> StructureToFluidAtInterface(Teuchos::RCP<
          const Epetra_Vector> iv) const;

      virtual Teuchos::RCP<Epetra_Vector> FluidToStructureAtInterface(Teuchos::RCP<
          const Epetra_Vector> iv) const;

      void BuidNoPenetrationMap();

      ADAPTER::Coupling& FluidStructureCoupling() { return *coupfs_; }

      //! Extractor used for constraint structure
      LINALG::MapExtractor consplitter_;

      //! map of gids of fluid  DOFs and respective structure DOFs
      Teuchos::RCP<std::map<int, int> > nopenetrationmap_;

      //! Extractor used for no penetration condition
      LINALG::MapExtractor nopenetration_;

      //! vector containing no penetration - conditions
      std::vector<DRT::Condition*> nopencond_;

      //! vector containing globla IDs of dofs with no penetration condition
      Teuchos::RCP<std::set<int> >condIDs_;

      /// coupling of fluid and structure (interface only)
      Teuchos::RCP<ADAPTER::Coupling> icoupfs_;

    private:
      //! coupling of fluid and structure (whole field)
      Teuchos::RCP<ADAPTER::Coupling> coupfs_;

      //@}

  }; // MonolithicBase
}

#endif /* PORO_BASE_H_ */
