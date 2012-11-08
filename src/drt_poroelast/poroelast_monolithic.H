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

#ifndef POROELAST_MONOLITHIC_H_
#define POROELAST_MONOLITHIC_H_

/*----------------------------------------------------------------------*
 | headers                                                              |
 *----------------------------------------------------------------------*/
#include "../drt_inpar/inpar_poroelast.H"
#include "../drt_inpar/inpar_structure.H"
#include "poro_base.H"

#include "../linalg/linalg_mapextractor.H"

/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
namespace LINALG
{
  class SparseMatrix;
  class SparseOperator;

  class BlockSparseMatrixBase;
  class Solver;
}

namespace POROELAST
{
  /// base class of all monolithic Poroelasticity algorithms
  class Monolithic: public PoroBase
  {

    public:
      //! create using a Epetra_Comm
      explicit Monolithic(const Epetra_Comm& comm,
          const Teuchos::ParameterList& timeparams
          );

      virtual ~Monolithic(){};

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
      //void TimeLoop();

      //! setup composed right hand side from field solvers
      virtual void SetupRHS(bool firstcall=false) ;

      //! setup composed system matrix from field solvers
      //void SetupSystemMatrix();

      /// setup composed system matrix from field solvers
      virtual void SetupSystemMatrix() { SetupSystemMatrix(*systemmatrix_); }

      /// setup composed system matrix from field solvers
      virtual void SetupSystemMatrix(LINALG::BlockSparseMatrixBase& mat) ;

      //! composed system matrix
      Teuchos::RCP<LINALG::BlockSparseMatrixBase> SystemMatrix() const
      {
        return systemmatrix_;
      }

      //! solve linear system
      void LinearSolve();

      //! create linear solver (setup of parameter lists, etc...)
      void CreateLinearSolver();

      //! @name Apply current field state to system

      //! Evaluate mechanical-fluid system matrix
      void ApplyStrCouplMatrix(Teuchos::RCP<LINALG::SparseOperator> k_sf //!< mechanical-fluid stiffness matrix
          );

      //! Evaluate fluid-mechanical system matrix
      void ApplyFluidCouplMatrix(Teuchos::RCP<LINALG::SparseOperator> k_fs //!< fluid-mechanical tangent matrix
            );

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

      //! newton iteration
      void NewtonFull();

      //! solve one time step
      void Solve();

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

      //@}

      //!finite difference check of stiffness matrix
      //void PoroFDCheck();

      /// Evaluate no penetration condition
      void EvaluateCondition(Teuchos::RCP<LINALG::SparseOperator> Sysmat,
                             Teuchos::RCP<Epetra_Vector> Cond_RHS,
                             int coupltype=0);

      /// recover Lagrange multiplier \f$\lambda_\Gamma\f$ at the interface at the end of each time step
      /// (i.e. condensed forces onto the structure) needed for rhs in next time step
      virtual void RecoverLagrangeMultiplier(){return;};

    protected:

      //! @name Access methods

      //! velocity calculation given the displacements (like in FSI)
      Teuchos::RCP<Epetra_Vector> CalcVelocity(
          Teuchos::RCP<const Epetra_Vector> dispnp);

      //@}

      //! setup newton solver
      void SetupNewton();

      //! convergence check for Newton solver
      virtual void BuildCovergenceNorms();

      //! convergence check for Newton solver
      bool ConvergenceCheck(int itnum, int itmax, double ittol);

      //! extract the field vectors from a given composed vector
      /*!
       x is the sum of all increments up to this point.
       \param x  (i) composed vector that contains all field vectors
       \param sx (o) structural vector (e.g. displacements)
       \param fx (o) fluid vector (e.g. velocities and pressure)
       */
      virtual void ExtractFieldVectors(Teuchos::RCP<const Epetra_Vector> x,
          Teuchos::RCP<const Epetra_Vector>& sx, Teuchos::RCP<
              const Epetra_Vector>& fx,
              bool firstcall=false);

      //! @name Access methods for subclasses

      //! full monolithic dof row map
      Teuchos::RCP<const Epetra_Map> DofRowMap() const {return blockrowdofmap_.FullMap();}

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
      virtual Teuchos::RCP<Epetra_Map> CombinedDBCMap();

      //! extractor to communicate between full monolithic map and block maps
      const LINALG::MultiMapExtractor& Extractor() const
      {
        return blockrowdofmap_;
      }

      //@}

      //! @name General purpose algorithm members
      //@{

      bool solveradapttol_;                   //!< adapt solver tolerance
      double solveradaptolbetter_;            //!< tolerance to which is adpated ????
      Teuchos::RCP<LINALG::Solver> solver_;   //!< linear algebraic solver

      //@}

      //! @name Printing and output
      //@{

      int printscreen_;     //!< print infos to standard out every printscreen_ steps
      bool printiter_;      //!< print intermediate iterations during solution
      bool printerrfile_;   //!< print infos to error file #errfile_
      FILE* errfile_;       //!< error file handle

      //@}

      //! @name Global vectors
      Teuchos::RCP<Epetra_Vector> zeros_; //!< a zero vector of full length
      //@}

      enum INPAR::STR::DynamicType strmethodname_; //!< enum for STR time integration

      //! block systemmatrix
      Teuchos::RCP<LINALG::BlockSparseMatrixBase> systemmatrix_;

      Teuchos::RCP<Epetra_Vector> rhs_; //!< rhs of Poroelasticity system

      //! dof row map (not splitted)
      Teuchos::RCP<Epetra_Map> fullmap_;

      //! dof row map splitted in (field) blocks
      LINALG::MultiMapExtractor blockrowdofmap_;

      Teuchos::RCP<LINALG::SparseMatrix> StructFluidCouplingMatrix();

      Teuchos::RCP<LINALG::SparseMatrix> FluidStructCouplingMatrix();

      Teuchos::RCP<LINALG::BlockSparseMatrixBase> StructFluidCouplingBlockMatrix();

      Teuchos::RCP<LINALG::BlockSparseMatrixBase> FluidStructCouplingBlockMatrix();


      //void UseBlockCoulpingMatrix();

      //coupling matrixes
      Teuchos::RCP<LINALG::SparseOperator> k_sf_;
      Teuchos::RCP<LINALG::SparseOperator> k_fs_;

    private:

      //! build block vector from field vectors, e.g. rhs, increment vector
      void SetupVector(Epetra_Vector &f, //!< vector of length of all dofs
          Teuchos::RCP<const Epetra_Vector> sv, //!< vector containing only structural dofs
          Teuchos::RCP<const Epetra_Vector> fv //!< vector containing only fluid dofs
          );

      //! @name Iterative solution technique

      enum INPAR::POROELAST::ConvNorm normtypeinc_;   //!< convergence check for residual temperatures
      enum INPAR::POROELAST::ConvNorm normtypefres_;  //!< convergence check for residual forces
      enum INPAR::POROELAST::BinaryOp combincfres_;   //!< binary operator to combine temperatures and forces

      double tolinc_;     //!< tolerance residual increment
      double tolfres_;    //!< tolerance force residual

      int itermax_;     //!< maximally permitted iterations
      int itermin_;     //!< minimally requested iteration
      double normrhs_;  //!< norm of residual forces
      double norminc_;  //!< norm of residual unknowns

      double normrhsfluidvel_;  //!< norm of residual forces (fluid velocity)
      double normincfluidvel_;  //!< norm of residual unknowns (fluid velocity)
      double normrhsfluidpres_; //!< norm of residual forces (fluid velocity)
      double normincfluidpres_; //!< norm of residual unknowns (fluid velocity)
      double normrhsfluid_;     //!< norm of residual forces (fluid )
      double normincfluid_;     //!< norm of residual unknowns (fluid )

      double normrhsstruct_;    //!< norm of residual forces (structure)
      double normincstruct_;    //!< norm of residual unknowns (structure)

      Epetra_Time timer_;       //!<timer for solution technique

      int iter_;                //!< iteration step

      //@}

      //! @name Various global forces

      Teuchos::RCP<Epetra_Vector> cond_rhs_;  //!< additional rhs-Terms for no penetration condition

      Teuchos::RCP<Epetra_Vector> freact_;    //!< reaction force
      Teuchos::RCP<Epetra_Vector> iterinc_;   //!< increment between Newton steps k and k+1
      //!< \f$\Delta{x}^{<k>}_{n+1}\f$

      //@}

  }; // Monolithic

} // namespace POROELAST


/*----------------------------------------------------------------------*/
#endif /* POROELAST_MONOLITHIC_H_ */
