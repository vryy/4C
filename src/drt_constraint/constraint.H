/*!----------------------------------------------------------------------
\file constraint.H

\brief Basic constraint class, dealing with constraints living on boundaries

<pre>
Maintainer: Thomas Kloeppel
            kloeppel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>

*----------------------------------------------------------------------*/

#ifndef CONSTRAINT_H
#define CONSTRAINT_H

#include <Teuchos_RCP.hpp>
#include <Epetra_Vector.h>
#include <Teuchos_ParameterList.hpp>

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/linalg_solver.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/linalg_systemmatrix.H"


using namespace std;
using namespace Teuchos;

/*!
\brief UTILS: namespace constaining general utils for BACI, as constraints, ...

*/
namespace UTILS
{
  /*!
  \brief Basic constraint class, dealing with constraint and monitor boundary conditions. 
  This class cannot handle multi point constraints, they will be dealt with by a derived class. 
  
  \author Kloeppel (kloeppel@lnm.mw.tum.de)
  */
  class Constraint
  
  {
  public:
  
    //! Constraint types
    enum ConstrType
    {
      none,
      volconstr3d,
      areaconstr3d,
      areaconstr2d,
      volmonitor3d,
      areamonitor3d,
      areamonitor2d,
      mpcnodeonplane3d,
      mpcnodeonline2d
    };
    
    /*!
    \brief Constructor of a constraint based on a conditions with a given name. It also
    takes care of the constraint IDs.
    */
  
    Constraint
    (
      RCP<DRT::Discretization> discr,  ///< Discretization constraint lives on
      const string& CondName, ///< Name of condition to creat constraint from
      int& minID,  ///< minimum constraint or monitor ID so far
      int& maxID  ///< maximum constraint or monitor ID so far   
    );  
  
    /*!
        \brief Destructor
  
     */
    ~Constraint()
    {
        return;
    };
  
    /*!
     \brief Return if there are constraints
    */
    bool HaveConstraint()
    {
        return constrtype_!=none;
    };
  
    /// Set state of the underlying discretization
    void SetState
    (
      const string& state,  ///< name of state to set
      RCP<Epetra_Vector> V  ///< values to set
    )
    {
      actdisc_->SetState(state,V);
    };
    
    void Initialize
    (
      ParameterList&        params,  ///< parameter list to communicate between elements and discretization
      RCP<Epetra_Vector>    systemvector3  ///< distributed vector that may be filled by assembly of element contributions
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
    
    /*
       \brief Compute values described by a monitor boundary condition
    */
    void ComputeMonitorValues
    (
      RCP<Epetra_Vector> disp ///< current displacement
    );
    
    ConstrType Type()
    {
      return constrtype_;
    }
    

  protected:
    
    RCP<DRT::Discretization> actdisc_; /// standard discretization
    vector<DRT::Condition*> constrcond_; /// conditions, that define the constraint (all of the same kind)
    ConstrType constrtype_; /// constraint type
    map<int,double>  inittimes_; /// map with times at which constraint is supposed to become active
    map<int,bool> activecons_; /// map with indicator if constraints are active
      
  private:
    
    // don't want = operator, cctor and destructor
  
    Constraint operator = (const Constraint& old);
    Constraint(const Constraint& old);
  
    //! Return the ConstrType based on the condition name
    ConstrType GetConstrType
    (
      const string& Name///<condition name
    );
    
    //! Evaluate constraint conditions and assemble the results
    void EvaluateConstraint
    ( 
      Teuchos::ParameterList&        params,  ///< parameter list to communicate between elements and discretization
      Teuchos::RCP<LINALG::SparseOperator> systemmatrix1, ///< sparse matrix that may be filled by assembly of element contributions
      Teuchos::RCP<LINALG::SparseOperator> systemmatrix2, ///< sparse (rectangular) matrix that may be filled by assembly of element contributions
      Teuchos::RCP<Epetra_Vector>    systemvector1, ///< distributed vector that may be filled by aasembly of element contributions
      Teuchos::RCP<Epetra_Vector>    systemvector2, ///< distributed vector that may be filled by aasembly of element contributions
      Teuchos::RCP<Epetra_Vector>    systemvector3  ///< distributed vector that may be filled by aasembly of element contributions
      );
    
    //! Compute and assemble initial constraint values (depending on user specific activation times)
    void InitializeConstraint
    (
      Teuchos::ParameterList&        params,  ///< parameter list to communicate between elements and discretization
      Teuchos::RCP<Epetra_Vector>    systemvector ///< distributed vector that may be filled by aasembly of element contributions
    ); 
  }; //class
} //namespace

#endif /*CONSTRAINT_H*/
