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
#include "../drt_inpar/inpar_structure.H"

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
  class Structure;
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
      Teuchos::RCP<Epetra_Vector> fint,  ///< vector of internal structural forces
      Teuchos::RCP<LINALG::SparseOperator> stiff,  ///< structural stiffness matrix
      Teuchos::ParameterList scalelist
    );

    /*!
     \brief Return cardiovascular0d rhs norm at generalized midpoint $t_{n+\theta}$
    */
    double GetCardiovascular0DRHSNorm() const
    {
      double foo;
      cardvasc0d_res_m_->Norm2(&foo);
      return foo;
    };

    /*!
     \brief Return cardiovascular0d rhs norm at generalized midpoint $t_{n+\theta}$
    */
    double GetCardiovascular0DRHSInfNorm() const
    {
      double foo;
      cardvasc0d_res_m_->NormInf(&foo);
      return foo;
    };

    /*!
     \brief Return cardiovascular0d dof inbcr norm
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
         \brief check periodic state of cardiovascular model
    */
    void CheckPeriodic();


    bool IsRealtiveEqualTo(const double A, const double B, const double Ref);

    bool ModuloIsRealtiveZero(const double value, const double modulo, const double Ref);

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
      Teuchos::ParameterList params,
      const Teuchos::RCP<Epetra_Vector> actpres,
      Teuchos::RCP<Epetra_Vector> systemvector, ///< structural rhs
      Teuchos::RCP<LINALG::SparseOperator> systemmatrix ///< structural stiffness matrix
    );

    /*!
         \brief Return cardiovascular0d rhs at generalized midpoint $t_{n+\theta}$
    */
    Teuchos::RCP<Epetra_Vector> GetCardiovascular0DRHS() const
    {
      return cardvasc0d_res_m_;
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
    Teuchos::RCP<Epetra_Vector> Get0D_dof_np() const
    {
      return cv0ddof_np_;
    };

    /*!
      \brief Return dof vector of last converged step
    */
    Teuchos::RCP<Epetra_Vector> Get0D_dof_n() const
    {
      return cv0ddof_n_;
    };

    /*!
      \brief Return dof vector of last converged step
    */
    Teuchos::RCP<Epetra_Vector> Get0D_dof_m() const
    {
      return cv0ddof_m_;
    };

    /*!
      \brief Return vol vector
    */
    Teuchos::RCP<Epetra_Vector> Get0D_vol_np() const
    {
      return v_np_;
    };


    Teuchos::RCP<Epetra_Vector> Get0D_df_np() const
    {
      return cardvasc0d_df_np_;
    };

    /*!
      \brief Return dof vector of last converged step
    */
    Teuchos::RCP<Epetra_Vector> Get0D_df_n() const
    {
      return cardvasc0d_df_n_;
    };

    Teuchos::RCP<Epetra_Vector> Get0D_f_np() const
    {
      return cardvasc0d_f_np_;
    };

    /*!
      \brief Return dof vector of last converged step
    */
    Teuchos::RCP<Epetra_Vector> Get0D_f_n() const
    {
      return cardvasc0d_f_n_;
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
    void Set0D_v_n
    (
      Teuchos::RCP<Epetra_Vector> newval  ///< new reference base values
    )
    {
      v_n_->Update(1.0,*newval,0.0);
    }

    /// set df_n, f_n
    void Set0D_df_n
    (
      Teuchos::RCP<Epetra_Vector> newval  ///< new Cardiovascular0D dofs
    )
    {
      cardvasc0d_df_n_->Update(1.0,*newval,0.0);
      return;
    }
    void Set0D_f_n
    (
      Teuchos::RCP<Epetra_Vector> newval  ///< new Cardiovascular0D dofs
    )
    {
      cardvasc0d_f_n_->Update(1.0,*newval,0.0);
      return;
    }

    /// Reset dofs
    void Set0D_dof_n
    (
      Teuchos::RCP<Epetra_Vector> newdof  ///< new Cardiovascular0D dofs
    )
    {
      cv0ddof_np_->Update(1.0,*newdof,0.0);
      cv0ddof_n_->Update(1.0,*newdof,0.0);
      return;
    }

    void PrintPresFlux(bool init) const;



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
      const Teuchos::RCP<Epetra_Vector> rhsstandard,  ///< standard right hand side
      const double k_ptc ///< for 3D-0D PTC
    );

    Teuchos::RCP<Cardiovascular0D> GetCardvasc0D4ElementWindkessel()
    {
      return cardvasc0d_4elementwindkessel_;
    }

    Teuchos::RCP<Cardiovascular0D> GetCardvasc0DArterialProxDist()
    {
      return cardvasc0d_arterialproxdist_;
    }

    Teuchos::RCP<Cardiovascular0D> GetCardvasc0DSysPulCirculation()
    {
      return cardvasc0d_syspulcirculation_;
    }

    Teuchos::RCP<Cardiovascular0D> GetCardvascRespir0DSysPulPeriphCirculation()
    {
      return cardvascrespir0d_syspulperiphcirculation_;
    }

    bool GetIsPeriodic() const
    {
      return is_periodic_;
    };

    double Get_k_ptc() const
    {
      return k_ptc_;
    };

    void Modify_k_ptc(const double sum, const double fac)
    {
      // increase PTC factor
      if(k_ptc_==0.0) k_ptc_+=sum;
      else k_ptc_*=fac;
    };

    void Reset_k_ptc()
    {
      // reset PTC factor - for adaptivity!
      k_ptc_=0.0;
    };


  private:

    // don't want = operator, cctor and destructor
    Cardiovascular0DManager operator = (const Cardiovascular0DManager& old);
    Cardiovascular0DManager(const Cardiovascular0DManager& old);


    Teuchos::RCP<DRT::Discretization> actdisc_; ///< discretization where elements of cardiovascular0d boundary live in
    int myrank_;
    Teuchos::RCP<LINALG::MapExtractor> dbcmaps_;  ///< map for Dirichlet DOFs
    Teuchos::RCP<Cardiovascular0DDofSet> cardiovascular0ddofset_; ///< degrees of freedom of pressures
    Teuchos::RCP<Epetra_Map> cardiovascular0dmap_;  ///< unique map of Cardiovascular0D values
    Teuchos::RCP<Epetra_Map> redcardiovascular0dmap_;   ///< fully redundant map of Cardiovascular0D values
    Teuchos::RCP<Epetra_Export> cardvasc0dimpo_;  ///< importer for fully redundant Cardiovascular0D vector into distributed one
    Teuchos::RCP<Epetra_Vector> cv0ddofincrement_; ///< increment of cvdof
    Teuchos::RCP<Epetra_Vector> cv0ddof_n_; ///< cvdof vector at t_{n}
    Teuchos::RCP<Epetra_Vector> cv0ddof_np_; ///< cvdof vector at t_{n+1}
    Teuchos::RCP<Epetra_Vector> cv0ddof_m_; ///< cvdof vector at mid-point
    Teuchos::RCP<Epetra_Vector> dcv0ddof_m_; ///< cvdof rate vector at mid-point
    Teuchos::RCP<Epetra_Vector> v_n_; ///< vol vector at t_{n}
    Teuchos::RCP<Epetra_Vector> v_np_; ///< vol vector at t_{n+1}
    Teuchos::RCP<Epetra_Vector> v_m_; ///< vol vector at mid-point
    Teuchos::RCP<Epetra_Vector> cv0ddof_T_N_; ///< cvdof vector at periodic time T_{N}
    Teuchos::RCP<Epetra_Vector> cv0ddof_T_NP_; ///< cvdof vector at periodic time T_{N+1}
    Teuchos::RCP<Epetra_Vector> cardvasc0d_res_m_; ///< Cardiovascular0D full rhs vector, at t_{n+theta}
    Teuchos::RCP<Epetra_Vector> cardvasc0d_df_n_; ///< Cardiovascular0D rhs part associated with time derivaties, at t_{n}
    Teuchos::RCP<Epetra_Vector> cardvasc0d_df_np_; ///< Cardiovascular0D rhs part associated with time derivaties, at t_{n+1}
    Teuchos::RCP<Epetra_Vector> cardvasc0d_df_m_; ///< Cardiovascular0D rhs part associated with time derivaties, at t_{n+theta}
    Teuchos::RCP<Epetra_Vector> cardvasc0d_f_n_; ///< Cardiovascular0D rhs part associated with non-derivatives, at t_{n}
    Teuchos::RCP<Epetra_Vector> cardvasc0d_f_np_; ///< Cardiovascular0D rhs part associated with non-derivatives, at t_{n+1}
    Teuchos::RCP<Epetra_Vector> cardvasc0d_f_m_; ///< Cardiovascular0D rhs part associated with non-derivatives, at t_{n+theta}
    const double T_period_; ///< periodic time
    const double eps_periodic_; ///< tolerance for periodic state
    bool is_periodic_; ///< true, if periodic state is reached, false otherwise
    double cycle_error_; ///< perdiodicity error
    int numCardiovascular0DID_;  ///< number of Cardiovascular0D bcs
    int Cardiovascular0DID_;  ///< smallest Cardiovascular0D bc id
    int offsetID_;  ///< smallest Cardiovascular0D bc id
    std::vector<int> currentID_;
    bool havecardiovascular0d_;  ///< are there Cardiovascular0D bcs at all?
    Teuchos::RCP<Cardiovascular0D> cardvasc0d_model_;
    Teuchos::RCP<Cardiovascular0D> cardvasc0d_4elementwindkessel_;
    Teuchos::RCP<Cardiovascular0D> cardvasc0d_arterialproxdist_;
    Teuchos::RCP<Cardiovascular0D> cardvasc0d_syspulcirculation_;
    Teuchos::RCP<Cardiovascular0D> cardvascrespir0d_syspulperiphcirculation_;
    Teuchos::RCP<LINALG::Solver> solver_; ///< solver for linear standard linear system
    Teuchos::RCP<LINALG::SparseMatrix> cardiovascular0dstiffness_;  ///< additional rectangular matrix
    Teuchos::RCP<LINALG::SparseMatrix> mat_dcardvasc0d_dd_;  ///< additional rectangular matrix
    Teuchos::RCP<LINALG::SparseMatrix> mat_dstruct_dcv0ddof_;  ///< additional rectangular matrix
    int counter_;  ///< counts how often #Solve is called
    bool isadapttol_;  ///< adaptive tolerance for solver?
    double adaptolbetter_; ///< adaptive tolerance for solver useful?
    double tolres_struct_; ///< tolerace for structural residual
    double tolres_cardvasc0d_; ///< tolerace for cardiovascular0d residual
    INPAR::CARDIOVASCULAR0D::Cardvasc0DSolveAlgo algochoice_;
    Teuchos::RCP<Epetra_Vector> dirichtoggle_;  ///< \b only for compatability: dirichlet toggle -- monitor its target change!
    Teuchos::RCP<Epetra_Vector> zeros_;  //!< a zero vector of full length
    const double theta_;
    const bool enhanced_output_;
    const bool ptc_3d0d_;
    double k_ptc_;
    double totaltime_;
    int linsolveerror_;
    Teuchos::ParameterList strparams_; ///< structure input parameters
    Teuchos::ParameterList cv0dparams_; ///< 0D cardiovascular input parameters
    INPAR::STR::IntegrationStrategy intstrat_; ///< structural time-integration strategy (old vs. standard)

  }; //class
}
#endif /*CARDIOVASCULAR0D_MANAGER_H*/
