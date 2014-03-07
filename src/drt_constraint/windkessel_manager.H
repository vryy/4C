/*!----------------------------------------------------------------------
\file windkessel_manager.cpp

\brief Class controlling Windkessel boundary conditions and containing the necessary data

**************************************************************************************************************************
Monolithic coupling of structure and a three-element Windkessel governed by
c dp/dt - c r2 dQ/dt + (p-p_0)/r1 - (1 + r2/r1) Q(d) = 0
[c: compliance, r1: first resistance, r2: second resistance, Q = -dV/dt: flux, p: pressure variable]
                            ____
                       ____|R_1 |___
                      |    |____|   |
----Q >----|_R_2_|----|             |--->----| p_0
                      |______|C|____|
                             | |

There are three different versions:
a) the standard model,

b) a special trimodular one where every parameter is pressure-dependent in terms of
PARAM(p) = (PARAM_b - PARAM_c)*0.5*(1.0 - tanh[(p-p_bc)/k_p] ) + PARAM_c + (PARAM_a - PARAM_b)*(1.0 - tanh[(p-p_ab)/k_p] )

c) a heart-specific model governing the arterial pressure with a three-element Windkessel with an additional valve law
(resistive Windkessel) infront of it, opening when arterial pressure is reached in the ventricels and closing when flux is
beginning to be reversed (cf. Sainte-Marie et. al. "Modeling and estimation of the cardiac electromechanical activity",
Comp. & Struct. 84 (2006) 1743-1759)

The Windkessel is momolithically coupled with the standard structural dynamics governing equation

M a + C v + f_int(d) - f_ext(d,p) = 0,

with Q being a function of the displacement vector d and f_ext additionally being a function of the Windkessel pressure p.
**************************************************************************************************************************

<pre>
Maintainer: Marc Hirschvogel
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
      \brief Constructor of Windkessel manager, allocating the Windkessels
    */
    WindkesselManager
    (
      Teuchos::RCP<DRT::Discretization> disc,  ///< standard discretization
      Teuchos::RCP<Epetra_Vector> disp,  ///< current displacement
      Teuchos::ParameterList param,  ///<  parameterlist from time integration algorithm
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
      \brief Change stiffness matrix and force vector according to the Windkessel functions.
      Values of pressure are taken from intern variable.
      Difference between current and prescribed values is calculated and stored as well.
    */
    void StiffnessAndInternalForces
    (
      const double time,  ///< time at end of time step
      Teuchos::RCP<Epetra_Vector> displast,  ///< displacement at beginning of time step
      Teuchos::RCP<Epetra_Vector> disp,  ///< displacement at end of time step
      Teuchos::ParameterList scalelist
    );

    /*!
     \brief Return norm of difference between actual and Windkessel values
    */
    double GetWindkesselRHSNorm() const
    {
      double foo;
      windkesselrhsm_->Norm2(&foo);
      return foo;
    };

    /*!
         \brief Return number of Windkessel functions
    */
    int GetNumberOfWindkesselFunctions() const
    {
        return numWindkesselID_;
    };

    /*!
     \brief Scale all pressures by a double d
    */
    void ScalePres
    (
      double d   ///< scale factor
    )
    {
        pn_->Scale(d);
        return;
    };

    /*!
         \brief Update Windkessel variables
    */
    void UpdateTimeStep();

    /// Add a vector as residual increment to the vector of pressures
    void UpdatePres
    (
      Teuchos::RCP<Epetra_Vector> presincrement  ///< vector to add
    );

    ///
    void EvaluateNeumannWindkesselCoupling
    (
      Teuchos::RCP<Epetra_Vector> actpres
    );

    void SetPres
    (
      Teuchos::RCP<Epetra_Vector> actpres
    );

    /*!
         \brief Return differences between prescribed and actual value of Windkessel number i
    */
    double GetWindkesselRHS
    (
      int i ///< ID of Windkessel of interest
    ) const
    {
        return (*windkesselrhsm_)[i];
    }

    /// return vector of differences between prescribed and actual values
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
    Teuchos::RCP<LINALG::SparseOperator> GetMatDstructDp() //const
    {
      return mat_dstruct_dp_;
    };

    //! Return the additional rectangular matrix, constructed for pressure evaluation
    Teuchos::RCP<LINALG::SparseMatrix> GetWindkesselStiffness() //const
    {
      return windkesselstiffness_;
    };

    /*!
      \brief Return pressure for Windkessel i
    */
    double GetPres
    (
      int i ///< ID of Windkessel of interest
    ) const
    {
      return (*pn_)[i];
    };

    /*!
      \brief Return pressure vector
    */
    Teuchos::RCP<Epetra_Vector> GetPresVector() const
    {
      return pn_;
    };

    /*!
      \brief Return pressure of last converged step
    */
    Teuchos::RCP<Epetra_Vector> GetPresVectorOld() const
    {
      return p_;
    };

    /*!
     \brief Return if there are Windkessels
    */
    bool HaveWindkessel() const
    {
      return havewindkessel_;
    };

    /*!
     \brief Return current value
    */
    double GetActVolValue
    (
      int i ///< ID of Windkessel of interest
    ) const
    {
        return (*vn_)[i];
    };

    /// Reset reference base values for restart computations
    void SetRefBaseValues
    (
      Teuchos::RCP<Epetra_Vector> newrefvals,  ///< new reference base values
      const double& time   ///< current time
    );

    /// Reset pressures
    void SetPresVector
    (
      Teuchos::RCP<Epetra_Vector> newpres  ///< new pressures
    )
    {
      pn_->Update(1.0,*newpres,0.0);
      p_->Update(1.0,*newpres,0.0);
      return;
    }

    void PrintPresFlux() const;

    /// Return Reference base values to write restart
    Teuchos::RCP<Epetra_Vector> GetRefVolValue() const
    {
      return v_;
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


    void Solve
    (
      Teuchos::RCP<LINALG::SparseMatrix> stiff,  ///< stiffness matrix
      Teuchos::RCP<Epetra_Vector> dispinc,      ///< displacement increment to compute
      const Teuchos::RCP<Epetra_Vector> rhsstandard  ///< standard right hand side
    );


  private:

    // don't want = operator, cctor and destructor
    WindkesselManager operator = (const WindkesselManager& old);
    WindkesselManager(const WindkesselManager& old);


    Teuchos::RCP<DRT::Discretization> actdisc_; ///< discretization where elements of Windkessel boundary live in
    Teuchos::RCP<WindkesselDofSet> windkesseldofset_; ///< degrees of freedom of pressures
    Teuchos::RCP<Epetra_Map> windkesselmap_;  ///< unique map of Windkessel values
    Teuchos::RCP<Epetra_Map> redwindkesselmap_;   ///< fully redundant map of Windkessel values
    Teuchos::RCP<Epetra_Export> windkimpo_;  ///< importer for fully redundant Windkessel vector into distributed one
    Teuchos::RCP<Epetra_Vector> actpres_;
    Teuchos::RCP<Epetra_Vector> p_; ///< pres vector at t_{n}
    Teuchos::RCP<Epetra_Vector> pn_; ///< pres vector at t_{n+1}
    Teuchos::RCP<Epetra_Vector> pm_; ///< pres vector at mid-point
    Teuchos::RCP<Epetra_Vector> dp_; ///< pres rate vector at t_{n}
    Teuchos::RCP<Epetra_Vector> dpn_; ///< pres rate vector at t_{n+1}
    Teuchos::RCP<Epetra_Vector> dpm_; ///< pres rate vector at mid-point
    Teuchos::RCP<Epetra_Vector> v_; ///< vol vector at t_{n}
    Teuchos::RCP<Epetra_Vector> vn_; ///< vol vector at t_{n+1}
    Teuchos::RCP<Epetra_Vector> vm_; ///< vol vector at mid-point
    Teuchos::RCP<Epetra_Vector> dv_; ///< time deriv of vol vector at t_{n}
    Teuchos::RCP<Epetra_Vector> dvn_; ///< time deriv of vol vector at t_{n+1}
    Teuchos::RCP<Epetra_Vector> dvm_; ///< time deriv of vol vector at mid-point
    Teuchos::RCP<Epetra_Vector> q_; ///< flux vector at t_{n}
    Teuchos::RCP<Epetra_Vector> qn_; ///< flux vector at t_{n+1}
    Teuchos::RCP<Epetra_Vector> qm_; ///< flux vector at mid-point
    Teuchos::RCP<Epetra_Vector> dq_; ///< time deriv of flux vector at t_{n}
    Teuchos::RCP<Epetra_Vector> dqn_; ///< time deriv of flux vector at t_{n+1}
    Teuchos::RCP<Epetra_Vector> dqm_; ///< time deriv of flux vector at mid-point
    Teuchos::RCP<Epetra_Vector> s_; ///< flux rate vector at t_{n}
    Teuchos::RCP<Epetra_Vector> sn_; ///< flux rate vector at t_{n+1}
    Teuchos::RCP<Epetra_Vector> sm_; ///< flux rate vector at mid-point
    Teuchos::RCP<Epetra_Vector> ds_; ///< time deriv of flux rate vector at t_{n}
    Teuchos::RCP<Epetra_Vector> dsn_; ///< time deriv of flux rate vector at t_{n+1}
    Teuchos::RCP<Epetra_Vector> dsm_; ///< time deriv of flux rate vector at mid-point
    Teuchos::RCP<Epetra_Vector> windkesselrhsm_; ///< Windkessel rhs vector at mid-point
    Teuchos::RCP<Epetra_Vector> windk_rhs_p_; ///< Windkessel rhs associated with p
    Teuchos::RCP<Epetra_Vector> windk_rhs_dp_; ///< Windkessel rhs associated with pressure rate dp/dt
    Teuchos::RCP<Epetra_Vector> windk_rhs_q_; ///< Windkessel rhs associated with q
    Teuchos::RCP<Epetra_Vector> windk_rhs_dq_; ///< Windkessel rhs associated with dq/dt
    Teuchos::RCP<Epetra_Vector> windk_rhs_s_; ///< Windkessel rhs associated with s
    Teuchos::RCP<Epetra_Vector> windk_rhs_ds_; ///< Windkessel rhs associated with ds/dt
    Teuchos::RCP<Epetra_Vector> windk_rhs_1_; ///< Windkessel rhs associated with 1
    int WindkesselID_;  ///< smallest Windkessel bc id
    int offsetID_;  ///< smallest Windkessel bc id
    int numWindkesselID_;  ///< number of Windkessel bcs
    std::vector<int> currentID;
    Teuchos::RCP<LINALG::SparseOperator> mat_dwindk_dd_;  ///< additional rectangular matrix
    Teuchos::RCP<LINALG::SparseOperator> mat_dstruct_dp_;  ///< additional rectangular matrix
    Teuchos::RCP<LINALG::SparseMatrix> windkesselstiffness_;  ///< additional rectangular matrix
    bool havewindkessel_;  ///< are there Windkessel bcs at all?
    Teuchos::RCP<Windkessel> wk_std_;
    Teuchos::RCP<Windkessel> wk_trimodular_;
    Teuchos::RCP<Windkessel> wk_heartvalvearterial_;
    int myrank_;
    Teuchos::RCP<LINALG::Solver> solver_; ///< solver for linear standard linear system
    int counter_;  ///< counts how often #Solve is called
    Teuchos::RCP<Epetra_Vector> dirichtoggle_;  ///< \b only for compatability: dirichlet toggle -- monitor its target change!
    Teuchos::RCP<LINALG::MapExtractor> dbcmaps_;  ///< map for Dirichlet DOFs
    Teuchos::RCP<Epetra_Vector> zeros_;  //!< a zero vector of full length
    double theta_;

  }; //class
}
#endif /*WINDKESSELMANAGER_H*/
