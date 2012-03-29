/*----------------------------------------------------------------------*/
/*!
 \file poroelast_monolithic.H

 \brief  Basis of all monolithic poroelasticity algorithms

 <pre>
   Maintainer: Anh-Tu Vuong
               vuong@lnm.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-15264
 </pre>
 *------------------------------------------------------------------------------------------------*/

#ifdef CCADISCRET
#ifndef POROELAST_MONOLITHIC_H_
#define POROELAST_MONOLITHIC_H_

/*----------------------------------------------------------------------*
 | headers                                                              |
 *----------------------------------------------------------------------*/
#include <Epetra_Time.h>

#include "../drt_adapter/adapter_structure.H"
#include "../drt_adapter/adapter_fluid_base_algorithm.H"
#include "../drt_adapter/adapter_algorithmbase.H"


/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
namespace LINALG
{
  class SparseMatrix;
  class MapExtractor;

  class BlockSparseMatrixBase;
  class Solver;
}

namespace ADAPTER
{
  class Coupling;
}

namespace POROELAST
{
  /// monolithic Poroelasticity algorithm base

  class MonolithicBase: public ADAPTER::AlgorithmBase,
      public ADAPTER::StructureBaseAlgorithm,
      public ADAPTER::FluidBaseAlgorithm
  {
    public:

      //! create using a Epetra_Comm
      explicit MonolithicBase(const Epetra_Comm& comm);

      //! virtual destructor to support polymorph destruction
      virtual ~MonolithicBase();

      /// read restart data
      virtual void ReadRestart(int step ///< step number where the calculation is continued
          );

      //! @name Time loop building blocks

      //! start a new time step
      virtual void PrepareTimeStep();

      //! take current results for converged and save for next time step
      virtual void Update();

      //! calculate stresses, strains, energies
      virtual void PrepareOutput()
      {
        StructureField().PrepareOutput();
      }

    protected:
      //! @name Transfer helpers

      //! field transform
      virtual Teuchos::RCP<Epetra_Vector> StructureToFluidField(Teuchos::RCP<
          const Epetra_Vector> iv) const;

      //! field transform
      virtual Teuchos::RCP<Epetra_Vector> FluidToStructureField(Teuchos::RCP<
          const Epetra_Vector> iv) const;

      void BuidNoPenetrationMap();

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

    private:
      //! coupling of fluid and structure (whole field)
    Teuchos::RCP<ADAPTER::Coupling> coupfa_;

      //@}

  }; // MonolithicBase


  /// base class of all monolithic Poroelasticity algorithms
  class Monolithic: public MonolithicBase
  {

    public:
      //! create using a Epetra_Comm
      explicit Monolithic(const Epetra_Comm& comm,
          const Teuchos::ParameterList& sdynparams);

      /*! do the setup for the monolithic system

       <pre>

       1.) setup coupling
       2.) get maps for all blocks in the system (and for the whole system as well)
       create combined map
       3.) create system matrix

       </pre>

       \note We want to do this setup after reading the restart information, not
       directly in the constructor. This is necessary since during restart (if
       ReadMesh is called), the dofmaps for the blocks might get invalid.
       */
      // Setup the monolithic Poroelasticity system
      virtual void SetupSystem();

      /// outer level Poroelasticity time loop
      void TimeLoop(const Teuchos::ParameterList& sdynparams //!< structural dynamic input parameters
          );

      //! setup composed right hand side from field solvers
      // implemented in FSI in the "split level"
      void SetupRHS();

      //! setup composed system matrix from field solvers
      void SetupSystemMatrix(const Teuchos::ParameterList& sdynparams);

      //! composed system matrix
      Teuchos::RCP<LINALG::BlockSparseMatrixBase> SystemMatrix() const
      {
        return systemmatrix_;
      }

      //! solve linear Poroelasticity system
      void LinearSolve();

      //! @name Apply current field state to system

      //! Evaluate mechanical-fluid system matrix
      void ApplyStrCouplMatrix(Teuchos::RCP<LINALG::SparseMatrix> k_sf, //!< mechanical-fluid stiffness matrix
          const Teuchos::ParameterList& sdynparams);

      //! Evaluate fluid-mechanical system matrix
      void ApplyFluidCouplMatrix(Teuchos::RCP<LINALG::SparseMatrix> k_fs, //!< fluid-mechanical tangent matrix
          const Teuchos::ParameterList& sdynparams);

      //@}

      //! evaluate all fields at x^n+1 with x^n+1 = x_n + stepinc
      virtual void Evaluate(Teuchos::RCP<const Epetra_Vector> stepinc //!< increment between time step n and n+1
          );

      //! extract initial guess from fields
      //! returns \f$\Delta x_{n+1}^{<k>}\f$
      virtual void InitialGuess(Teuchos::RCP<Epetra_Vector> ig);

      //! is convergence reached of iterative solution technique?
      //! keep your fingers crossed...
      //! \author lw  \date 12/07
      bool Converged();

      //! outer iteration loop
      void NewtonFull(const Teuchos::ParameterList& sdyn);

      //! @name Output

      //! print to screen information about residual forces and displacements
      //! \author lw (originally) \date 12/07
      void PrintNewtonIter();

      //! contains text to PrintNewtonIter
      //! \author lw (originally) \date 12/07
      void PrintNewtonIterText(FILE* ofile //!< output file handle
          );

      //! contains header to PrintNewtonIter
      //! \author lw (originally) \date 12/07
      void PrintNewtonIterHeader(FILE* ofile //!< output file handle
          );

      //! print statistics of converged Newton-Raphson iteration
      void PrintNewtonConv();

      //! write output
      void Output();

      //@}

      //!finite difference check of stiffness matrix
      void PoroFDCheck();

      /// Evaluate no penetration condition
      void EvaluateCondition(Teuchos::RCP<LINALG::SparseMatrix> Sysmat,
                             Teuchos::RCP<Epetra_Vector> Cond_RHS,
                             int coupltype=0);

    protected:

      //! @name Access methods

      //! velocity calculation given the displacements (like in FSI)
      Teuchos::RCP<Epetra_Vector> CalcVelocity(
          Teuchos::RCP<const Epetra_Vector> dispnp);

      //@}

      //! convergence check for Newton solver
      bool ConvergenceCheck(int itnum, int itmax, double ittol);

      //! extract the three field vectors from a given composed vector
      /*!
       x is the sum of all increments up to this point.
       \param x  (i) composed vector that contains all field vectors
       \param sx (o) structural vector (e.g. displacements)
       \param tx (o) thermal vector (e.g. temperatures)
       */
      virtual void ExtractFieldVectors(Teuchos::RCP<const Epetra_Vector> x,
          Teuchos::RCP<const Epetra_Vector>& sx, Teuchos::RCP<
              const Epetra_Vector>& fx);

      //! @name Access methods for subclasses

      //! full monolithic dof row map
      Teuchos::RCP<const Epetra_Map> DofRowMap() const
      {
        return blockrowdofmap_.FullMap();
      }

      //! set full monolithic dof row map
      /*!
       A subclass calls this method (from its constructor) and thereby
       defines the number of blocks, their maps and the block order. The block
       maps must be row maps by themselves and must not contain identical GIDs.
       */
      void SetDofRowMaps(
          const std::vector<Teuchos::RCP<const Epetra_Map> >& maps);

      //! combined DBC map
      //! unique map of all dofs that should be constrained with DBC
      Teuchos::RCP<Epetra_Map> CombinedDBCMap();

      //! extractor to communicate between full monolithic map and block maps
      const LINALG::MultiMapExtractor& Extractor() const
      {
        return blockrowdofmap_;
      }

      //@}

      //! @name General purpose algorithm members
      //@{

      bool solveradapttol_; //!< adapt solver tolerance
      double solveradaptolbetter_; //!< tolerance to which is adpated ????
      Teuchos::RCP<LINALG::Solver> solver_; //!< linear algebraic solver

      //@}

      //! @name Printing and output
      //@{

      int printscreen_; //!< print infos to standard out every printscreen_ steps
      bool printiter_; //!< print intermediate iterations during solution
      bool printerrfile_; //!< print infos to error file #errfile_
      FILE* errfile_; //!< error file handle

      //@}

      //! @name Global vectors
      Teuchos::RCP<Epetra_Vector> zeros_; //!< a zero vector of full length
      //@}

      enum INPAR::STR::DynamicType strmethodname_; //!< enum for STR time integartion

    private:

      //! dof row map splitted in (field) blocks
      LINALG::MultiMapExtractor blockrowdofmap_;

      //! build block vector from field vectors, e.g. rhs, increment vector
      void SetupVector(Epetra_Vector &f, //!< vector of length of all dofs
          Teuchos::RCP<const Epetra_Vector> sv, //!< vector containing only structural dofs
          Teuchos::RCP<const Epetra_Vector> fv //!< vector containing only fluid dofs
          );

      //! block systemmatrix
      Teuchos::RCP<LINALG::BlockSparseMatrixBase> systemmatrix_;

      //! @name Iterative solution technique

      enum INPAR::POROELAST::ConvNorm normtypeinc_; //!< convergence check for residual temperatures
      enum INPAR::POROELAST::ConvNorm normtypefres_; //!< convergence check for residual forces
      enum INPAR::POROELAST::BinaryOp combincfres_; //!< binary operator to combine temperatures and forces

      double tolinc_; //!< tolerance residual temperatures
      double tolfres_; //!< tolerance force residual
      int iter_; //!< iteration step
      int itermax_; //!< maximally permitted iterations
      int itermin_; //!< minimally requested iteration
      double normrhs_; //!< norm of residual forces
      double norminc_; //!< norm of residual unknowns

      double normrhsfluid_; //!< norm of residual forces (fluid)
      double normincfluid_; //!< norm of residual unknowns (fluid)

      double normrhsstruct_; //!< norm of residual forces (structure)
      double normincstruct_; //!< norm of residual unknowns (structure)

      Epetra_Time timer_; //!<timer for solution technique

      //@}

      //! @name Various global forces

      Teuchos::RCP<Epetra_Vector> rhs_; //!< rhs of Poroelasticity system
      // Teuchos::RCP<Epetra_Vector> rhs_fs_;//!< rhs of fluid-structure coupling
      // Teuchos::RCP<Epetra_Vector> rhs_sf_;//!< rhs of structure-fluid coupling
      Teuchos::RCP<Epetra_Vector> cond_rhs_; //!< additional rhs-Terms for no penetration condition

      Teuchos::RCP<Epetra_Vector> freact_; //!< reaction force
      Teuchos::RCP<Epetra_Vector> iterinc_; //!< increment between Newton steps k and k+1
      //!< \f$\Delta{x}^{<k>}_{n+1}\f$
      Teuchos::RCP<Epetra_Vector> veln_; //!< global structure velocities
      //!< \f${V}_{n+1}\f$
      //!< at \f$t_{n+1}\f$
      Teuchos::RCP<Epetra_Vector> dispn_; //!< global structure displacements
      //!< \f${d}_{n+1}\f$
      //!< at \f$t_{n+1}\f$

      //@}

  }; // Monolithic

} // namespace POROELAST


/*----------------------------------------------------------------------*/
#endif /* POROELAST_MONOLITHIC_H_ */
#endif  // CCADISCRET
