/*!----------------------------------------------------------------------
\file cardiovascular0d_manager.H

\brief Monolithic coupling of 3D structure 0D cardiovascular flow models

\level 2

<pre>
\maintainer Marc Hirschvogel
            hirschvogel@mhpc.mw.tum.de
            http://www.mhpc.mw.tum.de
            089 - 289-10363
</pre>

*----------------------------------------------------------------------*/

#ifndef CARDIOVASCULAR0D_MANAGER_H
#define CARDIOVASCULAR0D_MANAGER_H

#include <Teuchos_RCP.hpp>
#include <Epetra_Vector.h>
#include <Epetra_Operator.h>
#include <Epetra_RowMatrix.h>
#include <Teuchos_ParameterList.hpp>

#include "../drt_inpar/inpar_cardiovascular0d.H"

namespace IO
{
  class DiscretizationReader;
}

namespace DRT
{
  class Condition;
}

namespace DRT
{
  class Discretization;
}

namespace LINALG
{
  class SparseMatrix;
  class SparseOperator;
  class MapExtractor;
  class MultiMapExtractor;
  class Solver;
}

namespace UTILS
{

  //forward declarations
  class Cardiovascular0D;
  class Cardiovascular0DDofSet;

  class Cardiovascular0DManager
  {
  public:

    /*!
      \brief Constructor of cardiovascular0d manager
    */
    Cardiovascular0DManager
    (
      Teuchos::RCP<DRT::Discretization> disc,  ///< standard discretization
      Teuchos::RCP<const Epetra_Vector> disp,  ///< current displacement
      Teuchos::ParameterList strparams,  ///<  parameterlist from structural time integration algorithm
      Teuchos::ParameterList cv0dparams,  ///<  parameterlist from cardiovascular0d
      LINALG::Solver& solver  ///< Solver to solve linear subproblem in iteration
    );
    /*!
        \brief Destructor

     */
    ~Cardiovascular0DManager()
    {
        return;
    };

    /*!
      \brief Assemble cardiovascular0d stiffness and rhs contributions to full coupled problem
    */
    void EvaluateForceStiff
    (
      const double time,  ///< time at end of time step
      Teuchos::RCP<const Epetra_Vector> disp,  ///< displacement at end of time step
      Teuchos::ParameterList scalelist
    );

    /*!
     \brief Return cardiovascular0d rhs norm at generalized midpoint $t_{n+\theta}$
    */
    double GetCardiovascular0DRHSNorm() const
    {
      double foo;
      cardiovascular0drhsm_->Norm2(&foo);
      return foo;
    };

    /*!
     \brief Return cardiovascular0d rhs norm at generalized midpoint $t_{n+\theta}$
    */
    double GetCardiovascular0DDofIncrNorm() const
    {
      double foo;
      cv0ddofincrement_->Norm2(&foo);
      return foo;
    };

    /*!
     \brief Return cardiovascular0d rhs norm at generalized midpoint $t_{n+\theta}$
    */
    int GetCardiovascular0DLinSolveError() const
    {
      return linsolveerror_;
    };

    /*!
         \brief Update cardiovascular0d dofs
    */
    void UpdateTimeStep();

    /*!
         \brief Update cardiovascular0d dofs
    */
    void ResetStep();

    /// Add a vector as residual increment to the cardiovascular0d dof vector
    void UpdateCv0DDof
    (
      Teuchos::RCP<Epetra_Vector> cv0ddofincrement  ///< vector to add
    );

    ///
    void EvaluateNeumannCardiovascular0DCoupling
    (
      Teuchos::RCP<Epetra_Vector> actpres
    );

    /*!
         \brief Return cardiovascular0d rhs at generalized midpoint $t_{n+\theta}$
    */
    Teuchos::RCP<Epetra_Vector> GetCardiovascular0DRHS() const
    {
      return cardiovascular0drhsm_;
    }

    /*!
     \brief Return EpetraMap that determined distribution of Cardiovascular0D functions and pressures over processors
    */
    Teuchos::RCP<Epetra_Map> GetCardiovascular0DMap() const
    {
      return cardiovascular0dmap_;
    };

    //! Return the additional rectangular matrix, constructed for pressure evaluation
    Teuchos::RCP<LINALG::SparseMatrix> GetMatDcardvasc0dDd() //const
    {
      return mat_dcardvasc0d_dd_;
    };

//    //! Return the transpose of the additional rectangular matrix, constructed for pressure evaluation
//    Teuchos::RCP<LINALG::SparseMatrix> GetMatDcardvasc0dDdTrans() //const
//    {
//      return mat_dcardvasc0d_dd_->Transpose();
//    };

    //! Return the additional rectangular matrix, constructed for pressure evaluation
    Teuchos::RCP<LINALG::SparseMatrix> GetMatDstructDcv0ddof() //const
    {
      return mat_dstruct_dcv0ddof_;
    };

    //! Return the additional rectangular matrix, constructed for pressure evaluation
    Teuchos::RCP<LINALG::SparseMatrix> GetCardiovascular0DStiffness() //const
    {
      return cardiovascular0dstiffness_;
    };

    /*!
      \brief Return dof vector
    */
    Teuchos::RCP<Epetra_Vector> Get0DDofVector() const
    {
      return cv0ddofn_;
    };

    /*!
      \brief Return dof vector of last converged step
    */
    Teuchos::RCP<Epetra_Vector> Get0DDofVectorOld() const
    {
      return cv0ddof_;
    };

    /*!
     \brief Return if there are Cardiovascular0Ds
    */
    bool HaveCardiovascular0D() const
    {
      return havecardiovascular0d_;
    };

    /*!
     \brief Read restart information
    */
    void ReadRestart(
      IO::DiscretizationReader& reader,
      const double& time);

    /*!
     \brief Return structural input parameter list
    */
    Teuchos::ParameterList& StrParams() { return strparams_; }

    /*!
     \brief Return cardiovascular0d input parameter list
    */
    Teuchos::ParameterList& Cardvasc0DParams() { return cv0dparams_; }

    Teuchos::RCP<LINALG::Solver>& GetSolver() { return solver_; }

    /// Reset reference base values for restart computations
    void SetRefVolValue
    (
      Teuchos::RCP<Epetra_Vector> newrefval  ///< new reference base values
    );
    void SetRefFluxValue
    (
      Teuchos::RCP<Epetra_Vector> newrefval  ///< new reference base values
    );
    void SetRefDFluxValue
    (
      Teuchos::RCP<Epetra_Vector> newrefval  ///< new reference base values
    );
    void SetRefDDFluxValue
    (
      Teuchos::RCP<Epetra_Vector> newrefval  ///< new reference base values
    );

    /// Reset dofs
    void SetCv0DDofVector
    (
      Teuchos::RCP<Epetra_Vector> newdof  ///< new Cardiovascular0D dofs
    )
    {
      cv0ddofn_->Update(1.0,*newdof,0.0);
      cv0ddof_->Update(1.0,*newdof,0.0);
      return;
    }

    void PrintPresFlux(bool init) const;

    /// Return Reference base values to write restart
    Teuchos::RCP<Epetra_Vector> GetRefVolValue() const
    {
      return v_;
    }
    Teuchos::RCP<Epetra_Vector> GetRefFluxValue() const
    {
      return Q_;
    }
    Teuchos::RCP<Epetra_Vector> GetRefDFluxValue() const
    {
      return dQ_;
    }
    Teuchos::RCP<Epetra_Vector> GetRefDDFluxValue() const
    {
      return ddQ_;
    }

    //! switch Cardiovascular0D matrix to block matrix
    void UseBlockMatrix(Teuchos::RCP<const LINALG::MultiMapExtractor> domainmaps,
                        Teuchos::RCP<const LINALG::MultiMapExtractor> rangemaps);


    void SolverSetup
    (
      LINALG::Solver& solver,
      Teuchos::ParameterList params
    );


    int Solve
    (
      Teuchos::RCP<LINALG::SparseMatrix> stiff,  ///< stiffness matrix
      Teuchos::RCP<Epetra_Vector> dispinc,      ///< displacement increment to compute
      const Teuchos::RCP<Epetra_Vector> rhsstandard  ///< standard right hand side
    );


  private:

    // don't want = operator, cctor and destructor
    Cardiovascular0DManager operator = (const Cardiovascular0DManager& old);
    Cardiovascular0DManager(const Cardiovascular0DManager& old);


    Teuchos::RCP<DRT::Discretization> actdisc_; ///< discretization where elements of cardiovascular0d boundary live in
    Teuchos::RCP<Cardiovascular0DDofSet> cardiovascular0ddofset_; ///< degrees of freedom of pressures
    Teuchos::RCP<Epetra_Map> cardiovascular0dmap_;  ///< unique map of Cardiovascular0D values
    Teuchos::RCP<Epetra_Map> redcardiovascular0dmap_;   ///< fully redundant map of Cardiovascular0D values
    Teuchos::RCP<Epetra_Export> cardvasc0dimpo_;  ///< importer for fully redundant Cardiovascular0D vector into distributed one
    Teuchos::RCP<Epetra_Vector> actpres_;
    Teuchos::RCP<Epetra_Vector> cv0ddof_; ///< cvdof vector at t_{n}
    Teuchos::RCP<Epetra_Vector> cv0ddofn_; ///< cvdof vector at t_{n+1}
    Teuchos::RCP<Epetra_Vector> cv0ddofm_; ///< cvdof vector at mid-point
    Teuchos::RCP<Epetra_Vector> dcv0ddof_; ///< cvdof rate vector at t_{n}
    Teuchos::RCP<Epetra_Vector> dcv0ddofn_; ///< cvdof rate vector at t_{n+1}
    Teuchos::RCP<Epetra_Vector> dcv0ddofm_; ///< cvdof rate vector at mid-point
    Teuchos::RCP<Epetra_Vector> cv0ddofincrement_; ///< increment of cvdof
    Teuchos::RCP<Epetra_Vector> v_; ///< vol vector at t_{n}
    Teuchos::RCP<Epetra_Vector> vn_; ///< vol vector at t_{n+1}
    Teuchos::RCP<Epetra_Vector> vm_; ///< vol vector at mid-point
    Teuchos::RCP<Epetra_Vector> Q_; ///< flux vector at t_{n}
    Teuchos::RCP<Epetra_Vector> Qn_; ///< flux vector at t_{n+1}
    Teuchos::RCP<Epetra_Vector> Qm_; ///< flux vector at mid-point
    Teuchos::RCP<Epetra_Vector> dQ_; ///< flux rate vector at t_{n}
    Teuchos::RCP<Epetra_Vector> dQn_; ///< flux rate vector at t_{n+1}
    Teuchos::RCP<Epetra_Vector> dQm_; ///< flux rate vector at mid-point
    Teuchos::RCP<Epetra_Vector> ddQ_; ///< time deriv of flux rate vector at t_{n}
    Teuchos::RCP<Epetra_Vector> ddQn_; ///< time deriv of flux rate vector at t_{n+1}
    Teuchos::RCP<Epetra_Vector> ddQm_; ///< time deriv of flux rate vector at mid-point
    Teuchos::RCP<Epetra_Vector> cardiovascular0drhsm_; ///< Cardiovascular0D rhs vector at mid-point
    Teuchos::RCP<Epetra_Vector> cardvasc0d_rhs_cv0ddof_; ///< Cardiovascular0D rhs associated with p
    Teuchos::RCP<Epetra_Vector> cardvasc0d_rhs_dcv0ddof_; ///< Cardiovascular0D rhs associated with pressure rate dp/dt
    Teuchos::RCP<Epetra_Vector> cardvasc0d_rhs_Q_; ///< Cardiovascular0D rhs associated with q
    Teuchos::RCP<Epetra_Vector> cardvasc0d_rhs_dQ_; ///< Cardiovascular0D rhs associated with s
    Teuchos::RCP<Epetra_Vector> cardvasc0d_rhs_ddQ_; ///< Cardiovascular0D rhs associated with ds/dt
    Teuchos::RCP<Epetra_Vector> cardvasc0d_rhs_1_; ///< Cardiovascular0D rhs associated with 1
    Teuchos::RCP<Epetra_Vector> compvolm_; ///< compartment volume vector (V_at, V_ar, V_ven) at mid-point
    Teuchos::RCP<Epetra_Vector> compvol_; ///< compartment volume vector (V_at, V_ar, V_ven) at t_n
    int Cardiovascular0DID_;  ///< smallest Cardiovascular0D bc id
    int offsetID_;  ///< smallest Cardiovascular0D bc id
    int numCardiovascular0DID_;  ///< number of Cardiovascular0D bcs
    std::vector<int> currentID;
    Teuchos::RCP<LINALG::SparseMatrix> mat_dcardvasc0d_dd_;  ///< additional rectangular matrix
    Teuchos::RCP<LINALG::SparseMatrix> mat_dstruct_dcv0ddof_;  ///< additional rectangular matrix
    Teuchos::RCP<LINALG::SparseMatrix> cardiovascular0dstiffness_;  ///< additional rectangular matrix
    bool havecardiovascular0d_;  ///< are there Cardiovascular0D bcs at all?
    Teuchos::RCP<Cardiovascular0D> cardvasc0d_windkesselonly_;
    Teuchos::RCP<Cardiovascular0D> cardvasc0d_arterialproxdist_;
    Teuchos::RCP<Cardiovascular0D> cardvasc0d_arterialvenoussyspulcoupled_;
    int myrank_;
    Teuchos::RCP<LINALG::Solver> solver_; ///< solver for linear standard linear system
    int counter_;  ///< counts how often #Solve is called
    bool isadapttol_;  ///< adaptive tolerance for solver?
    double adaptolbetter_; ///< adaptive tolerance for solver useful?
    double tolres_struct_; ///< tolerace for structural residual
    double tolres_cardvasc0d_; ///< tolerace for cardiovascular0d residual
    INPAR::CARDIOVASCULAR0D::Cardvasc0DSolveAlgo algochoice_;
    Teuchos::RCP<Epetra_Vector> dirichtoggle_;  ///< \b only for compatability: dirichlet toggle -- monitor its target change!
    Teuchos::RCP<LINALG::MapExtractor> dbcmaps_;  ///< map for Dirichlet DOFs
    Teuchos::RCP<Epetra_Vector> zeros_;  //!< a zero vector of full length
    double theta_;
    bool enhanced_output_;
    double pstime_;
    double totaltime_;
    int linsolveerror_;
    Teuchos::ParameterList strparams_; ///< structure input parameters
    Teuchos::ParameterList cv0dparams_; ///< 0D cardiovascular input parameters
  }; //class
}
#endif /*0DCARDIOVASCULAR_MANAGER_H*/
