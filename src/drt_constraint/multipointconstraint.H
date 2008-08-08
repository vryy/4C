/*!----------------------------------------------------------------------
\file multipointconstraint.H

\brief Basic constraint class, dealing with multi point constraints 

<pre>
Maintainer: Thomas Kloeppel
            kloeppel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>

*----------------------------------------------------------------------*/

#ifndef MPCONSTRAINT_H
#define MPCONSTRAINT_H

#include "constraint.H"

using namespace std;
using namespace Teuchos;




namespace UTILS
{
  /*!
  \brief This class can handle multi point constraints. It is derived from the basic constraint class 
  and reimplements the evaluate routine.
  
  \author Kloeppel (kloeppel@lnm.mw.tum.de)
  */
  class MPConstraint : UTILS::Constraint
  {
  public:
  
    //! @name Constructors and destructors and related methods
  
    /*!
    \brief Standard Constructor
    */
  
    MPConstraint
    (
      RCP<DRT::Discretization> discr,  ///< Discretization constraint lives on
      const string& CondName, ///< Name of condition to creat constraint from
      int& minID,  ///< minimum constraint or monitor ID so far
      int& maxID  ///< maximum constraint or monitor ID so far   
    );  
  
    /*!
        \brief Destructor
  
     */
    ~MPConstraint()
    {
        return;
    };
  
    /// Set state of the underlying constraint discretization
    void SetConstrState
    (
      const string& state,  ///< name of state to set
      RCP<Epetra_Vector> V  ///< values to set
    )
    {
      if (constrtype_!=none)
      {
        constraintdis_->ClearState();
        constraintdis_->SetState(state,V);
      }
    };
      
    //@}
  
    //! @name Public Member Functions
  
    //! Initialize routine, to compute the right initial values, all constraints are based on
    void Initialize
    ( 
      Teuchos::ParameterList&  params,  ///< parameter list to communicate between elements and discretization
      Teuchos::RCP<Epetra_Vector> systemvector  ///< distributed vector that may be filled by assembly of element contributions
    ); 
    
    //! Evaluate routine to call from outside. In here the right action is determined and the #EvaluateConstraint routine is called 
    void Evaluate
    ( 
      Teuchos::ParameterList&  params,  ///< parameter list to communicate between elements and discretization
      Teuchos::RCP<LINALG::SparseOperator> systemmatrix1, ///< sparse matrix that may be filled by assembly of element contributions
      Teuchos::RCP<LINALG::SparseOperator> systemmatrix2, ///< sparse (rectangular) matrix that may be filled by assembly of element contributions
      Teuchos::RCP<Epetra_Vector> systemvector1, ///< distributed vector that may be filled by assembly of element contributions
      Teuchos::RCP<Epetra_Vector> systemvector2, ///< distributed vector that may be filled by assembly of element contributions
      Teuchos::RCP<Epetra_Vector> systemvector3  ///< distributed vector that may be filled by assembly of element contributions
    ); 
    
    bool HaveConstraint()
    {
        return constrtype_!=none;
    };
  
    //@}
  
  private:
      
    // don't want = operator, cctor and destructor
  
    MPConstraint operator = (const MPConstraint& old);
    MPConstraint(const MPConstraint& old);
  
  
    //! additional discretization consisting of constraint elements
    RCP<DRT::Discretization> constraintdis_;
  
    
    //! @name Private Member Functions
      
    //! Return the ConstrType based on the condition name
    ConstrType GetConstrType(const string& Name);///<condition name
    
    //! Evaluate constraint discretization and assemble the results
    void EvaluateConstraint( 
        RCP<DRT::Discretization> disc,  ///< discretization to evaluate
        Teuchos::ParameterList&        params,  ///< parameter list to communicate between elements and discretization
        Teuchos::RCP<LINALG::SparseOperator> systemmatrix1, ///< sparse matrix that may be filled by assembly of element contributions
        Teuchos::RCP<LINALG::SparseOperator> systemmatrix2, ///< sparse (rectangular) matrix that may be filled by assembly of element contributions
        Teuchos::RCP<Epetra_Vector>    systemvector1, ///< distributed vector that may be filled by aasembly of element contributions
        Teuchos::RCP<Epetra_Vector>    systemvector2, ///< distributed vector that may be filled by aasembly of element contributions
        Teuchos::RCP<Epetra_Vector>    systemvector3); ///< distributed vector that may be filled by aasembly of element contributions
    
    
    //! creating a new discretization based on conditions containing constraint elements
    RCP<DRT::Discretization> CreateDiscretizationFromCondition
    (
      RCP<DRT::Discretization> actdisc,
      vector< DRT::Condition* >      constrcond, ///< conditions as discretization basis  
      const string&             discret_name, ///< name of new discretization
      const string&             element_name ///< name of element type to create
    );
     
    //! Reorder MPC nodes based on condition input 
    void ReorderConstraintNodes(
         vector<int>& nodeids,  ///< reordered node ids
         const DRT::Condition* condname); ///< condition to deal with
     
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
     
     
  
    //@}
  
  }; //class
}
#endif /*MPCONSTRAINT_H*/
