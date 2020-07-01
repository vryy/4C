/*----------------------------------------------------------------------*/
/*! \file

 \brief  Basis of all porous media algorithms

 \level 2

 *-----------------------------------------------------------------------*/
#ifndef PORO_BASE_H_
#define PORO_BASE_H_

/*----------------------------------------------------------------------*
 | headers                                                              |
 *----------------------------------------------------------------------*/
#include <Epetra_Time.h>

#include "../drt_adapter/adapter_algorithmbase.H"
#include "../drt_adapter/ad_field.H"

#include "../drt_lib/drt_dserror.H"
#include "poroelast_utils.H"

#include "../linalg/linalg_mapextractor.H"

/*----------------------------------------------------------------------*
 | forward declarations                                                 |
 *----------------------------------------------------------------------*/
namespace DRT
{
  class Condition;
  class Discretization;
}  // namespace DRT

namespace LINALG
{
  class MapExtractor;
  class SparseMatrix;
  class BlockSparseMatrixBase;
}  // namespace LINALG

namespace ADAPTER
{
  class FluidPoro;
  class Coupling;
  class FPSIStructureWrapper;
  class MortarVolCoupl;
}  // namespace ADAPTER

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
namespace POROELAST
{
  // forward declaration
  class NoPenetrationConditionHandle;

  /// base class for porous media algorithm
  class PoroBase : public ADAPTER::AlgorithmBase, public ADAPTER::Field
  {
   public:
    //! create using a Epetra_Comm
    explicit PoroBase(const Epetra_Comm& comm, const Teuchos::ParameterList& timeparams);

    //! virtual destructor to support polymorph destruction
    virtual ~PoroBase();

    /// read restart data
    virtual void ReadRestart(const int step);

    //! outer level time loop
    virtual void TimeLoop();

    /// initialize system
    virtual void SetupSystem() = 0;

    //! perform result tests
    virtual void TestResults(const Epetra_Comm& comm);

    //! build combined dirichlet map for the monolithic problem
    virtual void BuildCombinedDBCMap()
    {
      dserror(
          "BuildCombinedDBCMap() not implemented in base class. must be implemented in sub "
          "classes.");
    };

    //! @name access methods

    //! access to structural field
    const Teuchos::RCP<ADAPTER::FPSIStructureWrapper>& StructureField() { return structure_; }

    //! access to fluid field
    const Teuchos::RCP<ADAPTER::FluidPoro>& FluidField() { return fluid_; }

    //! composed system matrix
    virtual Teuchos::RCP<LINALG::SparseMatrix> SystemMatrix()
    {
      dserror("SystemMatrix() only available for monolithic schemes!");
      return Teuchos::null;
    };

    //! block system matrix
    virtual Teuchos::RCP<LINALG::BlockSparseMatrixBase> BlockSystemMatrix()
    {
      dserror("BlockSystemMatrix() only available for monolithic schemes!");
      return Teuchos::null;
    };

    //! full monolithic dof row map
    virtual Teuchos::RCP<const Epetra_Map> DofRowMap()
    {
      dserror("DofRowMap() only available for monolithic schemes!");
      return Teuchos::null;
    };

    //! dof row map of Structure field
    virtual Teuchos::RCP<const Epetra_Map> DofRowMapStructure() = 0;

    //! dof row map of Fluid field
    virtual Teuchos::RCP<const Epetra_Map> DofRowMapFluid() = 0;

    //! extractor to communicate between full monolithic map and block maps
    virtual Teuchos::RCP<const LINALG::MultiMapExtractor> Extractor() const
    {
      dserror("ExtractorPointer only available for monolithic schemes!");
      return Teuchos::null;
    };

    //! unique map of all dofs that should be constrained with DBC
    virtual Teuchos::RCP<const Epetra_Map> CombinedDBCMap() const
    {
      dserror("CombinedDBCMap() only available for monolithic schemes!");
      return Teuchos::null;
    };

    //! return rhs of poro problem
    virtual Teuchos::RCP<const Epetra_Vector> RHS()
    {
      dserror("RHS() only available for monolithic schemes!");
      return Teuchos::null;
    };

    //@}

    //! update all fields at x^n+1_i+1 with x^n+1_i+1 = x_n+1_i + iterinc
    void UpdateStateIncrementally(
        Teuchos::RCP<const Epetra_Vector> iterinc  //!< increment between iteration i and i+1
        ) override
    {
      dserror("UpdateStateIncrementally() only available for monolithic schemes!");
    };

    //! evaluate all fields at x^n+1_i+1 with x^n+1_i+1 = x_n+1_i + iterinc
    virtual void Evaluate(
        Teuchos::RCP<const Epetra_Vector> iterinc  //!< increment between iteration i and i+1
    )
    {
      dserror("Evaluate() only available for monolithic schemes!");
    };

    //! evaluate all fields at x^n+1_i+1 with x^n+1_i+1 = x_n+1_i + iterinc
    virtual void Evaluate(
        Teuchos::RCP<const Epetra_Vector> sx, Teuchos::RCP<const Epetra_Vector> fx)
    {
      dserror("Evaluate(sx,fx) only available for monolithic schemes!");
    };

    //! evaluate all fields at x^n+1_i+1 with x^n+1_i+1 = x_n+1_i + iterinc
    virtual void Evaluate(
        Teuchos::RCP<const Epetra_Vector> iterinc,  //!< increment between iteration i and i+1
        bool firstiter)
    {
      dserror("Evaluate() only available for monolithic schemes!");
    };

    //! evaluate all fields at x^n+1_i+1 with x^n+1_i+1 = x_n+1_i + iterinc
    virtual void Evaluate(
        Teuchos::RCP<const Epetra_Vector> sx, Teuchos::RCP<const Epetra_Vector> fx, bool firstiter)
    {
      dserror("Evaluate(sx,fx) only available for monolithic schemes!");
    };

    //! solve time step (depending on algorithm)
    virtual void Solve() = 0;

    //! perform one time step (setup + solve + output)
    virtual void DoTimeStep() = 0;

    //! setup solver (for monolithic only)
    virtual bool SetupSolver() { return false; };

    virtual void SetupRHS(bool firstcall = false)
    {
      dserror("SetupRHS() only available for monolithic schemes!");
    };

    //! @name Time loop building blocks

    //! start a new time step
    virtual void PrepareTimeStep();

    //! take current results for converged and save for next time step
    virtual void Update();

    //! calculate stresses, strains, energies
    virtual void PrepareOutput();

    //! output
    virtual void Output(bool forced_writerestart = false);

    //@}

    //! return whether the poro discretization contains submeshes (i.e. it is coupled with a pure
    //! solid)
    bool HasSubmeshes() { return submeshes_; };

    //! Evaluate porosity at surface
    void CalculateSurfPoro(const std::string& condstring);

    //! return coupling object
    ADAPTER::Coupling& FluidStructureCoupling() { return *coupfs_; }

   protected:
    //! @name Transfer helpers

    //! field transform
    Teuchos::RCP<Epetra_Vector> StructureToFluidField(Teuchos::RCP<const Epetra_Vector> iv);

    //! field transform
    Teuchos::RCP<Epetra_Vector> FluidToStructureField(Teuchos::RCP<const Epetra_Vector> iv);

    //@}

    //! @name Transfer methods

    //! set fluid solution on structure
    void SetFluidSolution();

    //! set structure solution on fluid
    void SetStructSolution();

    //@}

    //! Extractor used for constraint structure
    Teuchos::RCP<LINALG::MapExtractor> consplitter_;

    //! Extractor used for poro structure interaction
    Teuchos::RCP<LINALG::MapExtractor> psiextractor_;

    //! helper class for no penetration condition
    Teuchos::RCP<NoPenetrationConditionHandle> noPenHandle_;

    //! flag for partial integration condition of porous fluid continuity equation
    bool partincond_;

    //! flag for pressure coupling condition
    bool presintcond_;

    //! flag for additional porosity degree of freedom
    bool porositydof_;

    Teuchos::RCP<LINALG::MapExtractor> porositysplitter_;

    //! true if the poroelast problem is only part of a larger problem (like e.g. in FPSI)
    bool PartOfMultifieldProblem_;

    //! @name Volume Mortar stuff

    //! flag for matchinggrid
    const bool matchinggrid_;

    //! volume coupling (using mortar) adapter
    Teuchos::RCP<::ADAPTER::MortarVolCoupl> volcoupl_;
    //@}

    //! flag for old time integration
    // TODO: delete this once unused!!
    const bool oldstructimint_;

    // flag indicating if nitsche contact is active
    bool nit_contact_;

   private:
    //! setup of everything, that is needed for the volumetric coupling
    // void SetupProxies();
    void SetupCoupling();

    //! add dof set of structure/fluid discretization to fluid/structure discretization
    void ReplaceDofSets();

    //! check for special poro conditions and set bools
    void CheckForPoroConditions();

    //! flag indicating not fully overlapping fluid and structure discretization
    bool submeshes_;

    //! coupling of fluid and structure (whole field)
    Teuchos::RCP<ADAPTER::Coupling> coupfs_;

    //! @name Underlying fields

    //! underlying structure of the poro problem
    Teuchos::RCP<::ADAPTER::FPSIStructureWrapper> structure_;

    //! underlying fluid of the poro problem
    Teuchos::RCP<::ADAPTER::FluidPoro> fluid_;

    //@}

  };  // PoroBase

  //! helper class for no penetration condition
  class NoPenetrationConditionHandle
  {
   public:
    //! create using a Epetra_Comm
    explicit NoPenetrationConditionHandle(std::vector<DRT::Condition*> nopencond)
        : condIDs_(Teuchos::null),
          condVector_(Teuchos::null),
          condRHS_(Teuchos::null),
          nopencond_(nopencond),
          nopenetration_(Teuchos::null),
          hascond_(false),
          fluidfluidConstraintMatrix_(Teuchos::null),
          fluidstructureConstraintMatrix_(Teuchos::null),
          structVelConstraintMatrix_(Teuchos::null)
    {
      if (nopencond_.size())
      {
        hascond_ = true;
        condIDs_ = Teuchos::rcp(new std::set<int>());
      }
    };

    ~NoPenetrationConditionHandle(){};

    //! build map containing dofs with no penetration condition (fluid)
    void BuidNoPenetrationMap(const Epetra_Comm& comm, Teuchos::RCP<const Epetra_Map> dofRowMap);

    //! apply rhs terms of no penetration condition to global rhs vector
    void ApplyCondRHS(Teuchos::RCP<Epetra_Vector> iterinc, Teuchos::RCP<Epetra_Vector> rhs);

    //! return no penetration map extractor
    Teuchos::RCP<const LINALG::MapExtractor> Extractor() { return nopenetration_; };

    //! return vector containing global IDs of dofs with no penetration condition
    Teuchos::RCP<std::set<int>> CondIDs() { return condIDs_; };

    //! return vector containing global IDs of dofs with no penetration condition
    Teuchos::RCP<Epetra_Vector> CondVector() { return condVector_; };

    //! check if a no penetration condition exists
    bool HasCond() { return hascond_; };

    //! return condrhs
    Teuchos::RCP<Epetra_Vector> RHS() { return condRHS_; };

    //! clear everything that is needed for coupling
    void Clear(POROELAST::coupltype coupltype = POROELAST::undefined);

    //! setup coupling matrixes and vecors
    void Setup(Teuchos::RCP<const Epetra_Map> dofRowMap, const Epetra_Map* dofRowMapFluid);

    //! return constraint matrix, that fits to coupling type
    Teuchos::RCP<LINALG::SparseMatrix> ConstraintMatrix(POROELAST::coupltype coupltype);

    //! return constraint matrix for structure velocity coupling
    Teuchos::RCP<LINALG::SparseMatrix> StructVelConstraintMatrix(POROELAST::coupltype coupltype);

   private:
    //! set containing global IDs of dofs with no penetration condition
    Teuchos::RCP<std::set<int>> condIDs_;

    //! vector marking dofs with no penetration condition
    Teuchos::RCP<Epetra_Vector> condVector_;

    //! vector containing rhs terms from no penetration condition
    Teuchos::RCP<Epetra_Vector> condRHS_;

    //! vector containing no penetration - conditions
    std::vector<DRT::Condition*> nopencond_;

    //! Extractor used for no penetration condition
    Teuchos::RCP<LINALG::MapExtractor> nopenetration_;

    //! flag indicating if a no penetration condition exists
    bool hascond_;

    //! @name coupling matrices
    Teuchos::RCP<LINALG::SparseMatrix> fluidfluidConstraintMatrix_;
    Teuchos::RCP<LINALG::SparseMatrix> fluidstructureConstraintMatrix_;
    Teuchos::RCP<LINALG::SparseMatrix> structVelConstraintMatrix_;
    //@}
  };
}  // namespace POROELAST

#endif /* PORO_BASE_H_ */
