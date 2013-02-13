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

#include "../drt_lib/drt_dserror.H"
#include "poroelast_utils.H"

/*----------------------------------------------------------------------*
 | forward declarations                                                 |
 *----------------------------------------------------------------------*/
namespace DRT
{
  class Condition;
  class Discretization;
}

namespace LINALG
{
  class MapExtractor;
  class SparseMatrix;
}

namespace ADAPTER
{
  class FluidPoro;
  class Coupling;
  class FSIStructureWrapper;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
/// POROELAST: namespace for all porous media problems
namespace POROELAST
{
  //forward declaration
  class NoPenetrationConditionHandle;

  /// base class for porous media algorithm
  class PoroBase: public ADAPTER::AlgorithmBase
  {
    public:

      //! create using a Epetra_Comm
      explicit PoroBase(const Epetra_Comm& comm, const Teuchos::ParameterList& timeparams);

      //! virtual destructor to support polymorph destruction
      virtual ~PoroBase();

      /// read restart data
      virtual void ReadRestart(int restart);

      //! outer level time loop
      virtual void TimeLoop();

      /// initialize system
      virtual void SetupSystem() = 0;

      //! perform result tests
      virtual void TestResults(const Epetra_Comm& comm);

      //! access to structural field
      const Teuchos::RCP<ADAPTER::FSIStructureWrapper>& StructureField(){return structure_;}

      //! access to fluid field
      ::ADAPTER::FluidPoro& FluidField(){return *fluid_;}

      //! solve time step (depending on algorithm)
      virtual void Solve() = 0;

      /// setup solver (for monolithic only)
      virtual bool SetupSolver(){return false;};

      //! @name Time loop building blocks

      //! start a new time step
      virtual void PrepareTimeStep();

      //! take current results for converged and save for next time step
      virtual void Update();

      //! calculate stresses, strains, energies
      virtual void PrepareOutput();

      //! output
      virtual void Output();

      //@}

      //! Evaluate porosity at surface
      void CalculateSurfPoro(const string& condstring);

    protected:

      //! @name Transfer helpers

      //! field transform
      Teuchos::RCP<Epetra_Vector> StructureToFluidField(Teuchos::RCP<
          const Epetra_Vector> iv) const;

      //! field transform
      Teuchos::RCP<Epetra_Vector> FluidToStructureField(Teuchos::RCP<
          const Epetra_Vector> iv) const;

      //@}

      //! @name Transfer methods

      //! set fluid solution on structure
      void SetFluidSolution();

      //! set structure solution on fluid
      void SetStructSolution();

      //@}


      //! return coupling object
      ADAPTER::Coupling& FluidStructureCoupling() { return *coupfs_; }

      //! Extractor used for constraint structure
      Teuchos::RCP< LINALG::MapExtractor > consplitter_;

      //! Extractor used for poro structure interaction
      Teuchos::RCP< LINALG::MapExtractor > psiextractor_;

      //! helper class for no penetration condition
      Teuchos::RCP<NoPenetrationConditionHandle> noPenHandle_;

      //! flag for partial integration condition of porous fluid continuity equation
      bool partincond_;

      //! flag for pressure coupling condition
      bool presintcond_;

    private:

      //! setup of everything, that is needed for the volumetric coupling
      void SetupProxiesAndCoupling();

      //! check for special poro conditions and set bools
      void CheckForPoroConditions();

      //!flag indicating not fully overlapping fluid and structure discretization
      bool submeshes_;

      //! coupling of fluid and structure (whole field)
      Teuchos::RCP<ADAPTER::Coupling> coupfs_;

      //! @name Underlying fields

      //! underlying structure of the poro problem
      Teuchos::RCP< ::ADAPTER::FSIStructureWrapper>    structure_;

      //! underlying fluid of the poro problem
      Teuchos::RCP< ::ADAPTER::FluidPoro>        fluid_;

      //@}

  }; // PoroBase

  //! helper class for no penetration condition
  class NoPenetrationConditionHandle
  {
  public:

    //! create using a Epetra_Comm
    explicit NoPenetrationConditionHandle( std::vector<DRT::Condition*> nopencond):
    condIDs_(Teuchos::null),
    condVector_(Teuchos::null),
    condRHS_(Teuchos::null),
    nopencond_(nopencond),
    nopenetration_(Teuchos::null),
    hascond_(false),
    fluidfluidConstraintMatrix_(Teuchos::null),
    fluidstructureConstraintMatrix_(Teuchos::null),
    structVelConstraintMatrix_(Teuchos::null)
    {
      if(nopencond_.size())
      {
        hascond_ = true;
        condIDs_ = Teuchos::rcp(new std::set<int>());
      }
    };

    ~NoPenetrationConditionHandle(){};

    //! build map containing dofs with no penetration condition (fluid)
    void BuidNoPenetrationMap(const Epetra_Comm& comm, Teuchos::RCP<const Epetra_Map> dofRowMap);

    //! apply rhs terms of no penetration condition to global rhs vector
    void ApplyCondRHS(Teuchos::RCP<Epetra_Vector> iterinc,
                      Teuchos::RCP<Epetra_Vector> rhs);

    //! return no penetration map extractor
    Teuchos::RCP< const LINALG::MapExtractor > Extractor(){return nopenetration_;};

    //! return vector containing global IDs of dofs with no penetration condition
    Teuchos::RCP<std::set<int> > CondIDs(){return condIDs_;};

    //! return vector containing global IDs of dofs with no penetration condition
    Teuchos::RCP<Epetra_Vector> CondVector(){return condVector_;};

    //! check if a no penetration condition exists
    bool HasCond(){return hascond_;};

    //! return condrhs
    Teuchos::RCP<Epetra_Vector> RHS(){return condRHS_;};

    //! clear everything that is needed for coupling
    void Clear(POROELAST::coupltype coupltype=POROELAST::undefined);

    //! setup coupling matrixes and vecors
    void Setup(Teuchos::RCP<const Epetra_Map> dofRowMap, const Epetra_Map* dofRowMapFluid);

    //! return constraint matrix, that fits to coupling type
    Teuchos::RCP<LINALG::SparseMatrix> ConstraintMatrix(POROELAST::coupltype coupltype);

    //! return constraint matrix for structure velocity coupling
    Teuchos::RCP<LINALG::SparseMatrix> StructVelConstraintMatrix(POROELAST::coupltype coupltype);

  private:

    //! set containing global IDs of dofs with no penetration condition
    Teuchos::RCP<std::set<int> > condIDs_;

    //! vector marking dofs with no penetration condition
    Teuchos::RCP<Epetra_Vector> condVector_;

    //! vector containing rhs terms from no penetration condition
    Teuchos::RCP<Epetra_Vector> condRHS_;

    //! vector containing no penetration - conditions
    std::vector<DRT::Condition*> nopencond_;

    //! Extractor used for no penetration condition
    Teuchos::RCP< LINALG::MapExtractor > nopenetration_;

    //! flag indicating if a no penetration condition exists
    bool hascond_;

    //! @name coupling matrices
    Teuchos::RCP<LINALG::SparseMatrix> fluidfluidConstraintMatrix_;
    Teuchos::RCP<LINALG::SparseMatrix> fluidstructureConstraintMatrix_;
    Teuchos::RCP<LINALG::SparseMatrix> structVelConstraintMatrix_;
    //@}
  };
}

#endif /* PORO_BASE_H_ */
