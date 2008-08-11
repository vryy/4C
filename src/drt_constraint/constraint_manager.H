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
  
    /*!
      \brief Constructor of constraint manager, allocating the constraints
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
  
    /*!
      \brief Change stiffness matrix and force vector according to the constraints.
      Values of lagrange multiplier are taken from intern variable.
      Difference between current and prescribed values is calculated and stored as well.
    */
    void StiffnessAndInternalForces
    (
      const double time,  ///< time at end of time step
      RCP<Epetra_Vector> displast,  ///< displacement at beginning of time step
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
     \brief Scale all lagrange multipliers by a double d
    */
    void ScaleLagrMult
    (
      double d   ///< scale factor
    )
    {
        lagrMultVec_->Scale(d);
        return;
    };
  
    /*!
         \brief Update lagrange multiplier \f$\lambda_{n+1}=\lambda_{n}+factor*\f$(volerr)
    */
    void UpdateLagrMult
    (
      double factor
    );
    
    /// Add a vector to the vector of lagrnage multipliers 
    void UpdateLagrMult
    (
      RCP<Epetra_Vector> vect  ///< vector to add
    );
  
    /*!
         \brief Compute difference between current and prescribed values at a given time and a given displacement
    */
    void ComputeError
    (
      double time,   /// time, at which the error is to compute at
      RCP<Epetra_Vector> disp  /// displacement vector at the given time  
    );
  
    /*!
         \brief Return differences between prescribed and actual value of constraint number i
    */
    double GetError
    (
      int i ///< ID of constraint of interest
    )
    {
        return (*constrainterr_)[i];
    }
  
    /// return vector of differences between prescribed and actual values
    RCP<Epetra_Vector> GetError()
    {
      return constrainterr_;
    }
  
    /*!
     \brief Return EpetraMap that determined distributation of constraints and lagrange
     multiplier over processors
    */
    RCP<Epetra_Map> GetConstraintMap()
    {
      return constrmap_;
    };
  
    //! Return the additional rectangular matrix, construted for lagrange multiplier evaluation
    RCP<LINALG::SparseMatrix> GetConstrMatrix()
    {
      return constrMatrix_;
    };
  
    /*!
      \brief Return lagrange multiplier for constraint i
    */
    double GetLagrMult
    (
      int i ///< ID of constraint of interest
    )
    {
      return (*lagrMultVec_)[i];
    };
  
    /*!
      \brief Return lagrange multiplier vector
    */
    void GetLagrMultVector
    (
      RCP<Epetra_Vector> tmpmult ///< vector to store lagrange multiplier in
    )
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
     \brief Return current value 
    */
    double GetCurrValue
    (
      int i ///< ID of constraint of interest
    )
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
    void ComputeMonitorValues
    (
      RCP<Epetra_Vector> disp ///< current diplacement
    );
    
    /// Reset reference base values for restart computations
    void SetRefBaseValues
    (
      RCP<Epetra_Vector> newrefvals,  ///< new reference base values
      const double& time   ///< current time
    );
    
    RCP<Epetra_Vector> GetRefBaseValues()
    {
      return refbasevalues_;
    }
    
  private:
    
    /*! \brief Sum values of a completely redundant vector into a distributed vector of the same size
    */
    void ImportResults
    (
      RCP<Epetra_Vector>& vect_dist, ///< distributed to store values in
      RCP<Epetra_Vector>& vect_redu,  ///< fully redundant and assembled vector containing element reults
      bool zero=true ///< is the distributed vector to be set to zero before adding the redundant vector?
    );
  
    /*!
      \brief synchronize processors after evaluating constraints by minimizing values returned woith a parameterlist
      over all processors
    */
    void SynchronizeMinConstraint
    (
      ParameterList& params,   ///< parameter list, values are stored in
      RCP<Epetra_Vector>& vect, ///< vector to store results in
      const char* resultstr  ///< name of value to take minimum of
    );
  
  
  private:
  
      // don't want = operator, cctor and destructor
  
    ConstrManager operator = (const ConstrManager& old);
    ConstrManager(const ConstrManager& old);
  
  
    RCP<DRT::Discretization> actdisc_; ///< discretization, elements to constraint live in 
    RCP<Epetra_Map> constrmap_;  ///< unique map of constraint values
    RCP<Epetra_Map> redconstrmap_;   ///< fully redundant map of constraint values
    RCP<Epetra_Map> monitormap_;  ///< unique map of monitor values
    RCP<Epetra_Map> redmonmap_;  ///< fully redundant map of monitor values
    RCP<Epetra_Vector> referencevalues_;  ///< reference at current time stepto constrain values to
    RCP<Epetra_Vector> refbasevalues_;   ///< reference base values at activation time of constrained structures 
    RCP<Epetra_Vector> actvalues_;   ///< current values of constrained structures
    RCP<Epetra_Vector> constrainterr_;   ///< vector with deflection between reference and current values
    RCP<Epetra_Vector> monitorvalues_;  ///< current values of monitored structures
    RCP<Epetra_Vector> initialmonvalues_;   ///< initial values of monitored structures
    int minConstrID_;  ///< smallest constraint boundary condition ID
    int numConstrID_;  ///< number of constraint boundary conditions
    int numMonitorID_;  ///< smallest monitor boundary condition ID
    int minMonitorID_;  ///< number monitor boundary condition ID
    RCP<Epetra_Vector> fact_;  ///< vector with current time curve values
    RCP<Epetra_Vector> lagrMultVec_;  ///< lagrange multipliers
    RCP<LINALG::SparseMatrix> constrMatrix_;  ///< additional rectangular matrix 
    bool haveconstraint_;  ///< are there constraints at all?
    bool havemonitor_;  ///< are there monitor conditions?
    double uzawaparam_;  ///< parameter of Uzawa algorithm (only for the case the linear uzawa is not used) 
  
    RCP<Constraint> volconstr3d_;  ///< 3d volume constraints defined on surfaces 
    RCP<Constraint> areaconstr3d_;  ///< 3d area constraints defined on surfaces
    RCP<Constraint> areaconstr2d_;  ///< 2d area constraints defined on lines
    RCP<MPConstraint> mpconplane3d_;  ///< 3d multipoint constraint prescribing the motion of a node relatively to a plane
    RCP<MPConstraint> mpconline2d_; ///< 2d multipoint constraint prescribing the motion of a node relatively to a straight line
  
    RCP<Constraint> volmonitor3d_; ///< 3d volume monitors defined on surfaces
    RCP<Constraint> areamonitor3d_;  ///< 3d area monitors defined on surfaces
    RCP<Constraint> areamonitor2d_;  ///< 2d area monitors defined on lines
    
  
  
  }; //class
}
#endif /*CONSTRAINTMANAGER_H*/
