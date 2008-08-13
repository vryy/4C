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
  \brief This pure virtual class is the commen interface for multi point constraints. 
  It is derived from the basic constraint class.
  
  \author Kloeppel (kloeppel@lnm.mw.tum.de)
  */
  class MPConstraint : public Constraint
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
        \brief Alternative Constructor
    */
    MPConstraint
    (
      RCP<DRT::Discretization> discr,  ///< Discretization constraint lives on
      const string& CondName ///< Name of condition to creat constraint from
    );  
  
    /*!
        \brief Destructor
    */
    virtual ~MPConstraint()
    {
      ;
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
        map<int,RCP<DRT::Discretization> >::iterator discrit;
        for(discrit=constraintdis_.begin();discrit!=constraintdis_.end();++discrit)
        {
          (discrit->second)->ClearState();
          (discrit->second)->SetState(state,V);
        }
      }
    };
      
    //@}
  
    //! @name Public Member Functions
  
    /// initialization routine called by the manager ctor to get correct reference base values and activating the right conditions at the beginning
    virtual void Initialize
    (
      ParameterList&        params,  ///< parameter list to communicate between elements and discretization
      RCP<Epetra_Vector>    systemvector3  ///< distributed vector that may be filled by assembly of element contributions
    )=0;
    
    /// initialization routine called at restart to activate the right conditions
    virtual void Initialize
    (
      const double& time  ///< current time
    )=0;
    
    //! Evaluate routine to call from outside. In here the right action is determined and the #EvaluateConstraint routine is called 
    virtual void Evaluate
    ( 
      Teuchos::ParameterList&  params,  ///< parameter list to communicate between elements and discretization
      Teuchos::RCP<LINALG::SparseOperator> systemmatrix1, ///< sparse matrix that may be filled by assembly of element contributions
      Teuchos::RCP<LINALG::SparseOperator> systemmatrix2, ///< sparse (rectangular) matrix that may be filled by assembly of element contributions
      Teuchos::RCP<Epetra_Vector> systemvector1, ///< distributed vector that may be filled by assembly of element contributions
      Teuchos::RCP<Epetra_Vector> systemvector2, ///< distributed vector that may be filled by assembly of element contributions
      Teuchos::RCP<Epetra_Vector> systemvector3  ///< distributed vector that may be filled by assembly of element contributions
    )=0; 
    
    bool HaveConstraint()
    {
        return constrtype_!=none;
    };
  
    //@}
  
  protected:
      
    // don't want = operator, cctor and destructor
  
//    MPConstraint operator = (const MPConstraint& old);
    MPConstraint(const MPConstraint& old);
  
  
    //! additional discretization consisting of constraint elements
    map<int,RCP<DRT::Discretization> > constraintdis_;
  
    
    //! @name Private Member Functions
      
    //! Evaluate constraint discretization and assemble the results
    virtual void EvaluateConstraint
    ( 
        RCP<DRT::Discretization> disc,  ///< discretization to evaluate
        Teuchos::ParameterList&        params,  ///< parameter list to communicate between elements and discretization
        Teuchos::RCP<LINALG::SparseOperator> systemmatrix1, ///< sparse matrix that may be filled by assembly of element contributions
        Teuchos::RCP<LINALG::SparseOperator> systemmatrix2, ///< sparse (rectangular) matrix that may be filled by assembly of element contributions
        Teuchos::RCP<Epetra_Vector>    systemvector1, ///< distributed vector that may be filled by aasembly of element contributions
        Teuchos::RCP<Epetra_Vector>    systemvector2, ///< distributed vector that may be filled by aasembly of element contributions
        Teuchos::RCP<Epetra_Vector>    systemvector3 ///< distributed vector that may be filled by aasembly of element contributions
    )=0;
    
    //! creating a new discretization based on conditions containing constraint elements
    virtual map<int,RCP<DRT::Discretization> > CreateDiscretizationFromCondition
    (
      RCP<DRT::Discretization> actdisc,
      vector< DRT::Condition* >      constrcond, ///< conditions as discretization basis  
      const string&             discret_name, ///< name of new discretization
      const string&             element_name, ///< name of element type to create
      int& startID
    )=0;
        
    ///replace numdofs in elements of constraint discretization
    void ReplaceNumDof
    (
      const RCP<DRT::Discretization> sourcedis,       ///< basis discretization
      const RCP<DRT::Discretization> constraintdis    ///< discretization with elements to modify
    ) const;
     
    /// find col node map so that we can evaluate the constraint elements
    RCP<Epetra_Map> ComputeNodeColMap(
             const RCP<DRT::Discretization> sourcedis,  ///< standard discretization we want to redistribute
             const RCP<DRT::Discretization> constraintdis ///< constraint discretization prescribing ghosting
             ) const;
  
    void PrintVec(ostream& os, const vector<int> actvec)
    {
      vector<int>::const_iterator i;
      for (i=actvec.begin(); i<actvec.end(); ++i) {
        os << *i << ",";
      }
      os << endl;
    }

    void PrintVec(ostream& os, const vector<double> actvec)
    {
      vector<double>::const_iterator i;
      for (i=actvec.begin(); i<actvec.end(); ++i) {
        os << *i << ",";
      }
      os << endl;
    }

    void PrintSet(ostream& os, const set<int> actset)
    {
      set<int>::iterator i;
      for (i=actset.begin(); i != actset.end(); ++i) {
        os << *i << ",";
      }
      os << endl;
      
    }
    //@}
  
  }; //class
}
#endif /*MPCONSTRAINT_H*/
