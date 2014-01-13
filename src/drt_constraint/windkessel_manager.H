/*!----------------------------------------------------------------------
\file windkessel_manager.cpp

\brief Class controlling Windkessel boundary conditions and containing the necessary data

**************************************************************************************************************************
Monolithic coupling of a three-element Windkessel governed by

either the standard linear version in p
(1) c dp/dt - c r2 dq/dt + p/r1 - (1 + r2/r1) q(d) = 0

or a special nonlinear heart version to mimic opened and closed valves
(2) c(p) dp/dt - c(p) r2(p) dq/dt + p/r1(p) - (1 + r2(p)/r1(p)) q(d) = 0

with
c(p) = (c_closed - c_opened)*0.5*(1.0 - tanh[(p-p_open)/k_p] ) + c_opened
r1(p) = (r1_closed - r1_opened)*0.5*(1.0 - tanh[(p-p_open)/k_p] ) + r1_opened
r2(p) = (r2_closed - r2_opened)*0.5*(1.0 - tanh[(p-p_open)/k_p] ) + r2_opened

[c: compliance, r1: first resistance, r2: second resistance, q = -dV/dt: flux, p: pressure variable]

and the standard structural dynamics governing equation

M a + C v + f_int(d) - f_ext(d,p) = 0,

with q being a function of the displacement vector d and f_ext additionally being a function of the Windkessel pressure p.
**************************************************************************************************************************

<pre>
Maintainer: Marc Hirschvogel
            hirschvogel@lnm.mw.tum.de
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
        presn_->Scale(d);
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
    Teuchos::RCP<LINALG::SparseOperator> GetCoupOffdiagVolD() //const
    {
      return coupoffdiag_vol_d_;
    };

    //! Return the additional rectangular matrix, constructed for pressure evaluation
    Teuchos::RCP<LINALG::SparseOperator> GetCoupOffdiagFextP() //const
    {
      return coupoffdiag_fext_p_;
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
      return (*presn_)[i];
    };

    /*!
      \brief Return pressure vector
    */
    Teuchos::RCP<Epetra_Vector> GetPresVector() const
    {
      return presn_;
    };

    /*!
      \brief Return pressure of last converged step
    */
    Teuchos::RCP<Epetra_Vector> GetPresVectorOld() const
    {
      return pres_;
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
        return (*voln_)[i];
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
      presn_->Update(1.0,*newpres,0.0);
      pres_->Update(1.0,*newpres,0.0);
      return;
    }

    void PrintPresFlux() const;

    /// Return Reference base values to write restart
    Teuchos::RCP<Epetra_Vector> GetRefVolValue() const
    {
      return vol_;
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
    Teuchos::RCP<Epetra_Vector> pres_; ///< pres vector at t_{n}
    Teuchos::RCP<Epetra_Vector> presn_; ///< pres vector at t_{n+1}
    Teuchos::RCP<Epetra_Vector> presm_; ///< pres vector at mid-point
    Teuchos::RCP<Epetra_Vector> presrate_; ///< pres rate vector at t_{n}
    Teuchos::RCP<Epetra_Vector> presraten_; ///< pres rate vector at t_{n+1}
    Teuchos::RCP<Epetra_Vector> presratem_; ///< pres rate vector at mid-point
    Teuchos::RCP<Epetra_Vector> flux_; ///< flux vector at t_{n}
    Teuchos::RCP<Epetra_Vector> fluxn_; ///< flux vector at t_{n+1}
    Teuchos::RCP<Epetra_Vector> fluxm_; ///< flux vector at mid-point
    Teuchos::RCP<Epetra_Vector> fluxrate_; ///< flux vector at t_{n}
    Teuchos::RCP<Epetra_Vector> fluxraten_; ///< flux vector at t_{n+1}
    Teuchos::RCP<Epetra_Vector> fluxratem_; ///< flux vector at mid-point
    Teuchos::RCP<Epetra_Vector> vol_; ///< flux vector at t_{n}
    Teuchos::RCP<Epetra_Vector> voln_; ///< flux vector at t_{n+1}
    Teuchos::RCP<Epetra_Vector> volm_; ///< flux vector at mid-point
    Teuchos::RCP<Epetra_Vector> windkesselrhsm_; ///< Windkessel rhs vector at mid-point
    Teuchos::RCP<Epetra_Vector> windk_rhs_p_; ///< Windkessel rhs associated with p
    Teuchos::RCP<Epetra_Vector> windk_rhs_dpdt_; ///< Windkessel rhs associated with dp/dt
    Teuchos::RCP<Epetra_Vector> windk_rhs_q_; ///< Windkessel rhs associated with q
    Teuchos::RCP<Epetra_Vector> windk_rhs_dqdt_; ///< Windkessel rhs associated with dq/dt
    Teuchos::RCP<Epetra_Vector> presnprint_;
    Teuchos::RCP<Epetra_Vector> fluxnprint_;
    int WindkesselID_;  ///< smallest Windkessel bc id
    int offsetID_;  ///< smallest Windkessel bc id
    int numWindkesselID_;  ///< number of Windkessel bcs
    std::vector<int> currentID;
    Teuchos::RCP<LINALG::SparseOperator> coupoffdiag_vol_d_;  ///< additional rectangular matrix
    Teuchos::RCP<LINALG::SparseOperator> coupoffdiag_fext_p_;  ///< additional rectangular matrix
    Teuchos::RCP<LINALG::SparseMatrix> windkesselstiffness_;  ///< additional rectangular matrix
    bool havewindkessel_;  ///< are there Windkessel bcs at all?
    RCP<Windkessel> rcr_;
    RCP<Windkessel> rcr_nlnheart_;
    int myrank_;
    Teuchos::RCP<LINALG::Solver> solver_; ///< solver for linear standard linear system
    int counter_;  ///< counts how often #Solve is called
    Teuchos::RCP<Epetra_Vector> dirichtoggle_;  ///< \b only for compatability: dirichlet toggle -- monitor its target change!
    Teuchos::RCP<LINALG::MapExtractor> dbcmaps_;  ///< map for Dirichlet DOFs
    Teuchos::RCP<Epetra_Vector> zeros_;  //!< a zero vector of full length

  }; //class
}
#endif /*WINDKESSELMANAGER_H*/
