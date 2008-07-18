/*!----------------------------------------------------------------------
\file constraint_manager.H

\brief Class controlling constraints and containing the necessary data

<pre>
Maintainer: Thomas Kloeppel
            kloeppel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>

*----------------------------------------------------------------------*/

#ifndef CONSTRAINT_MANAGER_H
#define CONSTRAINT_MANAGER_H

#include <Teuchos_RCP.hpp>
#include <Epetra_Vector.h>
#include <Epetra_Operator.h>
#include <Epetra_RowMatrix.h>


#include <Teuchos_ParameterList.hpp>
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/linalg_solver.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/linalg_systemmatrix.H"
#include "constraint.H"
#include "multipointconstraint.H"

using namespace std;
using namespace Teuchos;

namespace UTILS
{
  /*!
  \brief Class controlling constraints and containing the necessary data
  
  \author Kloeppel (kloeppel@lnm.mw.tum.de)
  */
  class ConstrManager
  {
  public:
  
    //! @name Constructors and destructors and related methods
  
    /*!
    \brief Standard Constructor
  
    */
  
    ConstrManager
    (
      RCP<DRT::Discretization> disc,  ///< standard discretization
      RCP<Epetra_Vector> disp,  ///< current displacement
      ParameterList param  ///<  parameterlist from timeintegration algorithm
    );
    /*!
        \brief Destructor
  
     */
    ~ConstrManager()
    {
        return;
    };
  
    //@}
  
    //! @name Public Member Functions
    /*!
      \brief Change stiffness matrix and force vector according to the constraints.
      Values of lagrange multiplier are taken from intern variable.
      Difference between current and prescribed values is calculated and stored as well.
    */
    void StiffnessAndInternalForces
    (
      const double time,  ///< time at end of time step
      RCP<Epetra_Vector> disp,  ///< displacement at end of time step
      RCP<Epetra_Vector> fint,  ///< vector of internal forces
      RCP<LINALG::SparseMatrix> stiff  ///< stiffness matrix
    );
  
    /*!
     \brief Return norm of difference between actual and constraint values
    */
    double GetErrorNorm()
    {
      double foo;
      constrainterr_->Norm2(&foo);
      return foo;
    };
  
    /*!
         \brief Return number of constraints
    */
    int GetNumberOfConstraints()
    {
        return numConstrID_;
    };
  
    /*!
             \brief Scale lagrange multiplier by a double d
    */
    void ScaleLagrMult(double d)
    {
        lagrMultVec_->Scale(d);
        return;
    };
  
    /*!
         \brief Update lagrange multiplier lambda_{n+1}=lambda_{n}+factor*(volerr)
    */
    void UpdateLagrMult(double factor);
  
    void UpdateLagrMult( RCP<Epetra_Vector> vect);
  
    /*!
         \brief Compute difference between current and prescribed values at a given time and a given displacement
    */
    void ComputeError(double time,RCP<Epetra_Vector> disp);
  
    /*!
         \brief Return differences between prescribed and actual value of constraint number i
    */
    double GetError(int i)
    {
        return (*constrainterr_)[i];
    }
  
    /// return vector of differences between prescribed and actual values
    RCP<Epetra_Vector> GetError()
    {
      return constrainterr_;
    }
    /*!
         \brief Add Increment to lagrange multiplier
    */
    void UpdateLagrMult();
  
    /*!
     \brief Return EpetraMap that determined distributation of constraints and lagrange
     multiplier over processors
    */
    RCP<Epetra_Map> GetConstraintMap()
    {
      return constrmap_;
    };
  
    RCP<LINALG::SparseMatrix> GetConstrMatrix()
    {
      return constrMatrix_;
    };
  
    /*!
      \brief Return lagrange multiplier for constraint i
    */
    double GetLagrMult(int i)
    {
      return (*lagrMultVec_)[i];
    };
  
    /*!
      \brief Return lagrange multiplier vector
    */
    void GetLagrMultVector(RCP<Epetra_Vector> tmpmult)
    {
      LINALG::Export(*lagrMultVec_,*tmpmult);
      return ;
    };
  
    /*!
     \brief Return if there are constraints
    */
    bool HaveConstraint()
    {
        return haveconstraint_;
    };
    /*!
         \brief Return if there are monitors
    */
    bool HaveMonitor()
    {
        return havemonitor_;
    };
  
    /*!
     \brief Return current value of volume i
    */
    double GetCurrValue(int i)
    {
        return (*actvalues_)[i];
    };
  
    /*
         \brief Print out the values of current monitor values
     */
    void PrintMonitorValues();
  
    /*
       \brief Compute values described by a monitor boundary condition
    */
    void ComputeMonitorValues(RCP<Epetra_Vector> disp);
    //@}
  private:
  
      //! name@ Private Member Functions
  
    /*! \brief Sum values of a completely redundant vector into a distributed vector of the same size
    */
    void ImportResults
    (
      RCP<Epetra_Vector>& vect_dist, ///< distributed to store values in
      RCP<Epetra_Vector>& vect_redu  ///< fully redundant and assembled vector containing element reults
    );
  
    /*!
      \brief synchronize processors after evaluating constraints by minimizing values returned woith a parameterlist
      over all processors
    */
  
    void SynchronizeMinConstraint(ParameterList& params,
            RCP<Epetra_Vector>& vect,const char* resultstr);
  
    /*!
        \brief Setup number of constraint ID's by evaluating the parameterlist returned by elements
    */
    void ManageIDs(ParameterList& params,
         int& minID,
         int& maxID,
         int& numID);
  
    /*!
           \brief Setup number of constraint ID's by evaluating the parameterlist returned by elements
    */
    void ManageIDs(ParameterList& params,
        int& minID,
        int& maxID,
        int& numID,
        vector<int>& MPCids);
  
  
    //@}
  
  private:
  
      // don't want = operator, cctor and destructor
  
    ConstrManager operator = (const ConstrManager& old);
    ConstrManager(const ConstrManager& old);
  
  
    //! name@ Private Attributes
    RCP<DRT::Discretization> actdisc_;
    RCP<Epetra_Map> constrmap_;
    RCP<Epetra_Map> redconstrmap_;
    RCP<Epetra_Map> monitormap_;
    RCP<Epetra_Map> redmonmap_;
    RCP<Epetra_Vector> referencevalues_;
    RCP<Epetra_Vector> initialvalues_;
    RCP<Epetra_Vector> actvalues_;
    RCP<Epetra_Vector> constrainterr_;
    RCP<Epetra_Vector> monitorvalues_;
    RCP<Epetra_Vector> initialmonvalues_;
    int minConstrID_;
    int maxConstrID_;
    int numConstrID_;
    int numMonitorID_;
    int maxMonitorID_;
    int minMonitorID_;
    RCP<Epetra_Vector> fact_;
    RCP<Epetra_Vector> lagrMultVec_;
    RCP<LINALG::SparseMatrix> constrMatrix_;
    bool haveconstraint_;
    bool havemonitor_;
    double uzawaparam_;
  
    RCP<Constraint> volconstr3d_;
    RCP<Constraint> areaconstr3d_;
    RCP<Constraint> areaconstr2d_;
    RCP<MPConstraint> mpconplane3d_;
    RCP<MPConstraint> mpconline2d_;
  
    RCP<Constraint> volmonitor3d_;
    RCP<Constraint> areamonitor3d_;
    RCP<Constraint> areamonitor2d_;
    //@}
  
  
  }; //class
}
#endif /*CONSTRAINTMANAGER_H*/
