/*!----------------------------------------------------------------------
\file windkessel_manager.H

\brief Monolithic coupling of 3D structure 0D Windkessel models

\level 2

<pre>
\maintainer Marc Hirschvogel
            hirschvogel@mhpc.mw.tum.de
            http://www.mhpc.mw.tum.de
            089 - 289-10363
</pre>

*----------------------------------------------------------------------*/

#ifndef WINDKESSEL_MANAGER_H
#define WINDKESSEL_MANAGER_H

#include <Teuchos_RCP.hpp>
#include <Epetra_Vector.h>
#include <Epetra_Operator.h>
#include <Epetra_RowMatrix.h>
#include <Teuchos_ParameterList.hpp>

#include "../drt_inpar/inpar_windkessel.H"

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
  class Windkessel;
  class WindkesselDofSet;

  class WindkesselManager
  {
  public:

    /*!
      \brief Constructor of windkessel manager
    */
    WindkesselManager
    (
      Teuchos::RCP<DRT::Discretization> disc,  ///< standard discretization
      Teuchos::RCP<Epetra_Vector> disp,  ///< current displacement
      Teuchos::ParameterList strparams,  ///<  parameterlist from structural time integration algorithm
      Teuchos::ParameterList wkparams,  ///<  parameterlist from windkessel
      LINALG::Solver& solver,  ///< Solver to solve linear subproblem in iteration
      Teuchos::RCP<LINALG::MapExtractor> dbcmaps  ///< Map extractor for Dirichlet DOFs
    );
    /*!
        \brief Destructor

     */
    ~WindkesselManager()
    {
        return;
    };

    /*!
      \brief Assemble windkessel stiffness and rhs contributions to full coupled problem
    */
    void StiffnessAndInternalForces
    (
      const double time,  ///< time at end of time step
      Teuchos::RCP<Epetra_Vector> displast,  ///< displacement at beginning of time step
      Teuchos::RCP<Epetra_Vector> disp,  ///< displacement at end of time step
      Teuchos::ParameterList scalelist
    );

    /*!
     \brief Return windkessel rhs norm at generalized midpoint $t_{n+\theta}$
    */
    double GetWindkesselRHSNorm() const
    {
      double foo;
      windkesselrhsm_->Norm2(&foo);
      return foo;
    };

    /*!
     \brief Return windkessel rhs norm at generalized midpoint $t_{n+\theta}$
    */
    double GetWindkesselDofIncrNorm() const
    {
      double foo;
      wkdofincrement_->Norm2(&foo);
      return foo;
    };

    /*!
     \brief Return windkessel rhs norm at generalized midpoint $t_{n+\theta}$
    */
    int GetWindkesselLinSolveError() const
    {
      return linsolveerror_;
    };

    /*!
         \brief Update windkessel dofs
    */
    void UpdateTimeStep();

    /*!
         \brief Update windkessel dofs
    */
    void ResetStep();

    /// Add a vector as residual increment to the windkessel dof vector
    void UpdateWkDof
    (
      Teuchos::RCP<Epetra_Vector> wkdofincrement  ///< vector to add
    );

    ///
    void EvaluateNeumannWindkesselCoupling
    (
      Teuchos::RCP<Epetra_Vector> actpres
    );

    /*!
         \brief Return windkessel rhs at generalized midpoint $t_{n+\theta}$
    */
    Teuchos::RCP<Epetra_Vector> GetWindkesselRHS() const
    {
      return windkesselrhsm_;
    }

    /*!
     \brief Return EpetraMap that determined distribution of Windkessel functions and pressures over processors
    */
    Teuchos::RCP<Epetra_Map> GetWindkesselMap() const
    {
      return windkesselmap_;
    };

    //! Return the additional rectangular matrix, constructed for pressure evaluation
    Teuchos::RCP<LINALG::SparseOperator> GetMatDwindkDd() //const
    {
      return mat_dwindk_dd_;
    };

    //! Return the additional rectangular matrix, constructed for pressure evaluation
    Teuchos::RCP<LINALG::SparseOperator> GetMatDstructDwkdof() //const
    {
      return mat_dstruct_dwkdof_;
    };

    //! Return the additional rectangular matrix, constructed for pressure evaluation
    Teuchos::RCP<LINALG::SparseMatrix> GetWindkesselStiffness() //const
    {
      return windkesselstiffness_;
    };

    /*!
      \brief Return dof vector
    */
    Teuchos::RCP<Epetra_Vector> GetWkDofVector() const
    {
      return wkdofn_;
    };

    /*!
      \brief Return dof vector of last converged step
    */
    Teuchos::RCP<Epetra_Vector> GetWkDofVectorOld() const
    {
      return wkdof_;
    };

    /*!
     \brief Return if there are Windkessels
    */
    bool HaveWindkessel() const
    {
      return havewindkessel_;
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
     \brief Return windkessel input parameter list
    */
    Teuchos::ParameterList& WkParams() { return wkparams_; }

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
    void SetWkDofVector
    (
      Teuchos::RCP<Epetra_Vector> newdof  ///< new Windkessel dofs
    )
    {
      wkdofn_->Update(1.0,*newdof,0.0);
      wkdof_->Update(1.0,*newdof,0.0);
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

    //! switch Windkessel matrix to block matrix
    void UseBlockMatrix(Teuchos::RCP<const LINALG::MultiMapExtractor> domainmaps,
                        Teuchos::RCP<const LINALG::MultiMapExtractor> rangemaps);


    void SolverSetup
    (
      Teuchos::RCP<DRT::Discretization> discr,
      LINALG::Solver& solver,
      Teuchos::RCP<LINALG::MapExtractor> dbcmaps,
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
    WindkesselManager operator = (const WindkesselManager& old);
    WindkesselManager(const WindkesselManager& old);


    Teuchos::RCP<DRT::Discretization> actdisc_; ///< discretization where elements of windkessel boundary live in
    Teuchos::RCP<WindkesselDofSet> windkesseldofset_; ///< degrees of freedom of pressures
    Teuchos::RCP<Epetra_Map> windkesselmap_;  ///< unique map of Windkessel values
    Teuchos::RCP<Epetra_Map> redwindkesselmap_;   ///< fully redundant map of Windkessel values
    Teuchos::RCP<Epetra_Export> windkimpo_;  ///< importer for fully redundant Windkessel vector into distributed one
    Teuchos::RCP<Epetra_Vector> actpres_;
    Teuchos::RCP<Epetra_Vector> wkdof_; ///< wkdof vector at t_{n}
    Teuchos::RCP<Epetra_Vector> wkdofn_; ///< wkdof vector at t_{n+1}
    Teuchos::RCP<Epetra_Vector> wkdofm_; ///< wkdof vector at mid-point
    Teuchos::RCP<Epetra_Vector> dwkdof_; ///< wkdof rate vector at t_{n}
    Teuchos::RCP<Epetra_Vector> dwkdofn_; ///< wkdof rate vector at t_{n+1}
    Teuchos::RCP<Epetra_Vector> dwkdofm_; ///< wkdof rate vector at mid-point
    Teuchos::RCP<Epetra_Vector> wkdofincrement_; ///< increment of wkdof
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
    Teuchos::RCP<Epetra_Vector> windkesselrhsm_; ///< Windkessel rhs vector at mid-point
    Teuchos::RCP<Epetra_Vector> windk_rhs_wkdof_; ///< Windkessel rhs associated with p
    Teuchos::RCP<Epetra_Vector> windk_rhs_dwkdof_; ///< Windkessel rhs associated with pressure rate dp/dt
    Teuchos::RCP<Epetra_Vector> windk_rhs_Q_; ///< Windkessel rhs associated with q
    Teuchos::RCP<Epetra_Vector> windk_rhs_dQ_; ///< Windkessel rhs associated with s
    Teuchos::RCP<Epetra_Vector> windk_rhs_ddQ_; ///< Windkessel rhs associated with ds/dt
    Teuchos::RCP<Epetra_Vector> windk_rhs_1_; ///< Windkessel rhs associated with 1
    Teuchos::RCP<Epetra_Vector> compvolm_; ///< compartment volume vector (V_at, V_ar, V_ven) at mid-point
    Teuchos::RCP<Epetra_Vector> compvol_; ///< compartment volume vector (V_at, V_ar, V_ven) at t_n
    int WindkesselID_;  ///< smallest Windkessel bc id
    int offsetID_;  ///< smallest Windkessel bc id
    int numWindkesselID_;  ///< number of Windkessel bcs
    std::vector<int> currentID;
    Teuchos::RCP<LINALG::SparseOperator> mat_dwindk_dd_;  ///< additional rectangular matrix
    Teuchos::RCP<LINALG::SparseOperator> mat_dstruct_dwkdof_;  ///< additional rectangular matrix
    Teuchos::RCP<LINALG::SparseMatrix> windkesselstiffness_;  ///< additional rectangular matrix
    bool havewindkessel_;  ///< are there Windkessel bcs at all?
    Teuchos::RCP<Windkessel> wk_std_;
    Teuchos::RCP<Windkessel> wk_heartvalvearterial_proxdist_;
    Teuchos::RCP<Windkessel> wk_heartvalvecardiovascular_full_;
    int myrank_;
    Teuchos::RCP<LINALG::Solver> solver_; ///< solver for linear standard linear system
    int counter_;  ///< counts how often #Solve is called
    bool isadapttol_;  ///< adaptive tolerance for solver?
    double adaptolbetter_; ///< adaptive tolerance for solver useful?
    double tolres_struct_; ///< tolerace for structural residual
    double tolres_windk_; ///< tolerace for windkessel residual
    INPAR::WINDKESSEL::WindkSolveAlgo algochoice_;
    Teuchos::RCP<Epetra_Vector> dirichtoggle_;  ///< \b only for compatability: dirichlet toggle -- monitor its target change!
    Teuchos::RCP<LINALG::MapExtractor> dbcmaps_;  ///< map for Dirichlet DOFs
    Teuchos::RCP<Epetra_Vector> zeros_;  //!< a zero vector of full length
    double theta_;
    bool enhanced_output_;
    INPAR::STR::PreStress pstype_;
    double pstime_;
    double totaltime_;
    int linsolveerror_;
    Teuchos::ParameterList strparams_; ///< structure input parameters
    Teuchos::ParameterList wkparams_; ///< windkessel input parameters
  }; //class
}
#endif /*WINDKESSELMANAGER_H*/
