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
#include "../drt_lib/drt_utils.H"


using namespace std;
using namespace Teuchos;



/*!
\brief Class controlling constraints and containing the necessary data

\author Kloeppel (kloeppel@lnm.mw.tum.de)
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

    ConstrManager(RCP<DRT::Discretization> disc,
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
    //methods used by public functions
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
    
    //! creating a new discretization based on conditions containing constraint elements
    RCP<DRT::Discretization> CreateDiscretizationFromCondition(
            vector< DRT::Condition* >      constrcond, ///< conditions as discretization basis  
            const string&             discret_name, ///< name of new discretization
            const string&             element_name ///< name of element type to create
            );
    
    //! Reorder MPC nodes based on condition input 
    void ReorderConstraintNodes(
        vector<int>& nodeids,  ///< reordered node ids
        const DRT::Condition*      condname); ///< condition to deal with
    
    ///replace numdofs in elements of constraint discretization
    void ReplaceNumDof(
            const RCP<DRT::Discretization> sourcedis,       ///< basis discretization
            const RCP<DRT::Discretization> constraintdis    ///< discretization with elements to modify
            ) const;
    
    /// find col node map so that we can evaluate the constraint elements
    RCP<Epetra_Map> ComputeNodeColMap(
            const RCP<DRT::Discretization> sourcedis,  ///< standard discretization we want to redistribute
            const RCP<DRT::Discretization> constraintdis ///< constraint discretization prescribing ghosting
            ) const;
    
    //! fill vectors with MPC Ids and amplitudes and update parameter list
    void SetupMPC(
        vector<double>& amplit,   ///< vector for MPC amplitudes (as given in the input file)
        vector<int>& IDs, ///< vector for MPC IDs
        const vector<DRT::Condition*>& constrcond, ///< conditions to deal with
        ParameterList& param ///< parameter list given to the elements
        );
    
    //! Evaluate a discretization and assemble to a rectangular matrix
    void Evaluate( 
        RCP<DRT::Discretization> disc,  ///< discretization to evaluate
        Teuchos::ParameterList&        params,  ///< parameter list to communicate between elements and discretization
        Teuchos::RCP<LINALG::SparseOperator> systemmatrix1, ///< sparse matrix that may be filled by assembly of element contributions
        Teuchos::RCP<LINALG::SparseOperator> systemmatrix2, ///< sparse (rectangular) matrix that may be filled by assembly of element contributions
        Teuchos::RCP<Epetra_Vector>    systemvector1, //< distributed vector that may be filled by aasembly of element contributions
        Teuchos::RCP<Epetra_Vector>    systemvector2, //< distributed vector that may be filled by aasembly of element contributions
        Teuchos::RCP<Epetra_Vector>    systemvector3, //< distributed vector that may be filled by aasembly of element contributions
        vector<DRT::Condition*>& constrcond);  //< conditions determing where to assemble vector into rectangular matrix
    
    //! Evaluate conditions and assemble to a rectangular matrix
    void EvaluateCondition( 
        RCP<DRT::Discretization> disc,  ///< discretization to evaluate
        Teuchos::ParameterList&        params,  ///< parameter list to communicate between elements and discretization
        Teuchos::RCP<LINALG::SparseOperator> systemmatrix1, ///< sparse matrix that may be filled by assembly of element contributions
        Teuchos::RCP<LINALG::SparseOperator> systemmatrix2, ///< sparse (rectangular) matrix that may be filled by assembly of element contributions
        Teuchos::RCP<Epetra_Vector>    systemvector1, //< distributed vector that may be filled by aasembly of element contributions
        Teuchos::RCP<Epetra_Vector>    systemvector2, //< distributed vector that may be filled by aasembly of element contributions
        Teuchos::RCP<Epetra_Vector>    systemvector3, //< distributed vector that may be filled by aasembly of element contributions
        vector<DRT::Condition*>& constrcond);  //< conditions determing where to assemble vector into rectangular matrix
    //@}

private:

    // don't want = operator, cctor and destructor

    ConstrManager operator = (const ConstrManager& old);
    ConstrManager(const ConstrManager& old);


    //! name@ Private Attributes
    RCP<DRT::Discretization> actdisc_;
    RCP<DRT::Discretization> constraintdis_;
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
    RCP<LINALG::SparseMatrix> constrMatrix_;
    bool haveconstraint_;
    bool havevolconstr_;
    bool haveareaconstr3D_;
    bool haveareaconstr2D_;
    bool havenodeconstraint_;
    bool havemonitor_;
    bool haveareamonitor3D_;
    bool haveareamonitor2D_;
    bool havevolmonitor_;
    double uzawaparam_;
    //@}


}; //class

#endif /*CONSTRAINTMANAGER_H*/
