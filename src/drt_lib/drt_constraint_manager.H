/*!----------------------------------------------------------------------
\file drt_constraint_manager.H

\brief Class controlling constraints and containing the necessary data

<pre>
Maintainer: Thomas Kloeppel
            kloeppel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>

*----------------------------------------------------------------------*/

#ifndef DRT_CONSTRAINT_MANAGER_H
#define DRT_CONSTRAINT_MANAGER_H

#include <Teuchos_RCP.hpp>
#include <Epetra_Vector.h>
#include <Epetra_Operator.h>
#include <Epetra_RowMatrix.h>


#include <Teuchos_ParameterList.hpp>
#include "./drt_discret.H"
#include "./linalg_solver.H"
#include "./linalg_utils.H"
#include "linalg_systemmatrix.H"


using namespace std;
using namespace Teuchos;



/*!
\brief Class controlling constraints and containing the necessary data
*/
class ConstrManager
{
public:

    class StruGenAlpha;
    //! @name Constructors and destructors and related methods

    /*!
    \brief Standard Constructor

    \param disc (in): Discretization
    \param disp (in): current displacement to compute initial constraints
    */

    ConstrManager(DRT::Discretization& disc,
        RCP<Epetra_Vector> disp,
        ParameterList param);

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

      \param time (in):      time at end of time step
      \param disp (in):      Displacement at end of time step
      \param fint (in/out):     Vector of internal forces
      \param stiff (in/out):   Stiffness matrix
    */
    void StiffnessAndInternalForces(const double time,
          RCP<Epetra_Vector> disp,
          RCP<Epetra_Vector> fint,
          RCP<LINALG::SparseMatrix> stiff);

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
         \brief Scale increment of lagrange multiplier increment by a double d
    */
    void ScaleLagrIncr(double d)
    {
        lagrMultInc_->Scale(d);
        return;
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
         \brief Update increment of lagrange multiplier Delta lambda_{n+1}=Delta lambda_{n}+factor*(vect+volerr)
    */
    void UpdateLagrIncr(double factor, Epetra_Vector vect);

    /*!
         \brief Update lagrange multiplier lambda_{n+1}=lambda_{n}+factor*(volerr)
    */
    void UpdateLagrMult(double factor);
    /*!
         \brief Return lagrange multiplier increment
    */
    double GetLagrIncr(int i)
    {
        return (*lagrMultInc_)[i];
    };

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
      \brief Set all lagrange multiplier for restart, assuming that Epetra_Maps fit together
    */
    void SetLagrMultVector(RCP<Epetra_Vector> restartvec)
    {
      LINALG::Export(*restartvec,*lagrMultVec_);
      return;
    };

    /*!
      \brief Set uzawa parameter (used for dynamic adaptation and restart)
    */
    void SetUzawaParameter(double restartval)
    {
      uzawaparam_=restartval;
      return;
    }    
    
    /*!
       \brief Return the current value of the uzawa parameter
    */
    double GetUzawaParameter()
    {
      return uzawaparam_;
    }
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

    /*!
     \brief Computes inner product of constraint multivector times lagrange multiplier increment
     vector as required for the uzawa algorithm to compute new displacement increment

     \param dotprod (out):      resulting value of B^T * Delta lambda
     */
    void ComputeConstrTimesLagrIncr(RCP<Epetra_Vector> dotprod);

    /*!
      \brief Return inner product of constraint multivector times displacement increment

      \param disi (in)      :      displacement increment Delta u
      \param dotprod (out)     :      resulting value of B * Delta u
      */
    void ComputeConstrTimesDisi(Epetra_Vector disi,RCP<Epetra_Vector> dotprod);

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
    //methods used by constructor and other public functions
    /*!
      \brief synchronize processors after evaluating constraints by summing up values returned by elements
      over all processors
    */

    void SynchronizeSumConstraint(ParameterList& params,
            RCP<Epetra_Vector>& vect,const char* resultstr, const int numID, const int minID);

    /*!
      \brief synchronize processors after evaluating constraints by minimizing values returned woith a parameterlist
      over all processors
    */

    void SynchronizeMinConstraint(ParameterList& params,
            RCP<Epetra_Vector>& vect,const char* resultstr);
    
    /*!
        \brief Setup number of constraint ID's by evaluating the parameterlist returned by elements 
    */
    
    void ManageIDs(ParameterList& params,int& minID, int& maxID, int& numID);
    //@}

private:

    // don't want = operator, cctor and destructor

    ConstrManager operator = (const ConstrManager& old);
    ConstrManager(const ConstrManager& old);


    //! name@ Private Attributes
    DRT::Discretization& actdisc_;
    RCP<Epetra_Map> constrmap_;
    RCP<Epetra_Map> monitormap_;
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
    RCP<Epetra_Vector> lagrMultInc_;
    //RCP<Epetra_MultiVector> constrVec_;
    RCP<LINALG::SparseMatrix> constrMatrix_;
    bool haveconstraint_;
    bool havevolconstr_;
    bool haveareaconstr3D_;
    bool haveareaconstr2D_;
    bool havemonitor_;
    bool haveareamonitor3D_;
    bool haveareamonitor2D_;
    bool havevolmonitor_;
    double uzawaparam_;
    //@}


}; //class

#endif /*CONSTRAINTMANAGER_H*/
