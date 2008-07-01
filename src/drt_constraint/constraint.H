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
\brief Basic constraint class, dealing with constraint and monitor boundary conditions. 
This class cannot handle multi point constraints, they will be dealt with by a derived class. 

\author Kloeppel (kloeppel@lnm.mw.tum.de)
*/

class Constraint
{
public:

  //! @name Constructors and destructors and related methods

  /*!
  \brief Standard Constructor
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

  //@}

  //! @name Public Member Functions

  /*!
   \brief Return if there are constraints
  */
  bool HaveConstraint()
  {
      return constrtype_!=none;
  };

  //! Evaluate conditions and assemble to a rectangular matrix
  void Evaluate
  ( 
    Teuchos::ParameterList&  params,  ///< parameter list to communicate between elements and discretization
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix1, ///< sparse matrix that may be filled by assembly of element contributions
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix2, ///< sparse (rectangular) matrix that may be filled by assembly of element contributions
    Teuchos::RCP<Epetra_Vector> systemvector1, //< distributed vector that may be filled by aasembly of element contributions
    Teuchos::RCP<Epetra_Vector> systemvector2, //< distributed vector that may be filled by aasembly of element contributions
    Teuchos::RCP<Epetra_Vector> systemvector3  //< distributed vector that may be filled by aasembly of element contributions
  ); 
  
  /*
     \brief Compute values described by a monitor boundary condition
  */
  void ComputeMonitorValues
  (
    RCP<Epetra_Vector> disp ///< current displacement
  );
  //@}

private:
    
  enum ConstrType
  {
    none,
    volconstr3d,
    areaconstr3d,
    areaconstr2d,
    volmonitor3d,
    areamonitor3d,
    areamonitor2d
  };
  
  // don't want = operator, cctor and destructor

  Constraint operator = (const Constraint& old);
  Constraint(const Constraint& old);


  //! name@ Private Attributes
  RCP<DRT::Discretization> actdisc_;
  vector<DRT::Condition*> constrcond_;
  ConstrType constrtype_;
  //@}
  
  //! name@ Private Member Functions
    
  //! Return the ConstrType based on the condition name
  ConstrType GetConstrType(const string& Name);///<condition name
  
  //! Evaluate conditions and assemble to a rectangular matrix
  void EvaluateConstraint( 
      RCP<DRT::Discretization> disc,  ///< discretization to evaluate
      Teuchos::ParameterList&        params,  ///< parameter list to communicate between elements and discretization
      Teuchos::RCP<LINALG::SparseOperator> systemmatrix1, ///< sparse matrix that may be filled by assembly of element contributions
      Teuchos::RCP<LINALG::SparseOperator> systemmatrix2, ///< sparse (rectangular) matrix that may be filled by assembly of element contributions
      Teuchos::RCP<Epetra_Vector>    systemvector1, //< distributed vector that may be filled by aasembly of element contributions
      Teuchos::RCP<Epetra_Vector>    systemvector2, //< distributed vector that may be filled by aasembly of element contributions
      Teuchos::RCP<Epetra_Vector>    systemvector3); //< distributed vector that may be filled by aasembly of element contributions
  //@}




}; //class

#endif /*CONSTRAINTMANAGER_H*/
