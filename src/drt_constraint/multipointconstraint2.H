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

#ifndef MPCONSTRAINT2_H
#define MPCONSTRAINT2_H

#include "multipointconstraint.H"

using namespace std;
using namespace Teuchos;


namespace UTILS
{
  /*!
  \brief This class can handle twodimensional multi point constraints. It is derived from the basic constraint class 
  and reimplements the evaluate routine.
  
  \author Kloeppel (kloeppel@lnm.mw.tum.de)
  */
  class MPConstraint2 : public MPConstraint
  {
  public:
  
    //! @name Constructors and destructors and related methods
  
    /*!
    \brief Standard Constructor
    */
  
    MPConstraint2
    (
      RCP<DRT::Discretization> discr,  ///< Discretization constraint lives on
      const string& CondName, ///< Name of condition to creat constraint from
      int& minID,  ///< minimum constraint or monitor ID so far
      int& maxID  ///< maximum constraint or monitor ID so far   
    );  
  
    /*!
        \brief Destructor
  
     */
      
    //@}
  
    //! @name Public Member Functions
  
    /// initialization routine called by the manager ctor to get correct reference base values and activating the right conditions at the beginning
    void Initialize
    (
      ParameterList&        params,  ///< parameter list to communicate between elements and discretization
      RCP<Epetra_Vector>    systemvector  ///< distributed vector that may be filled by assembly of element contributions
    );
    
    /// initialization routine called at restart to activate the right conditions
    void Initialize
    (
      const double& time  ///< current time
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
  
    MPConstraint2 operator = (const MPConstraint2& old);
    MPConstraint2(const MPConstraint2& old);
  
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
    
  
    //@}
  
  }; //class
}
#endif /*MPCONSTRAINT_H*/
