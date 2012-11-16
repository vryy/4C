/*!----------------------------------------------------------------------
\file multipointconstraint3penalty.H

\brief Basic constraint class, dealing with multi point constraints

<pre>
Maintainer: Thomas Kloeppel
            kloeppel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>

*----------------------------------------------------------------------*/

#ifndef MPCONSTRAINT3PENALTY_H
#define MPCONSTRAINT3PENALTY_H

#include "multipointconstraint.H"

using namespace std;
using namespace Teuchos;



namespace UTILS
{
  /*!
  \brief This class can handle multi point constraints in 3D.
  It is derived from the basic multipointconstraint class.

  \author Kloeppel (kloeppel@lnm.mw.tum.de)
  */
  class MPConstraint3Penalty : public UTILS::MPConstraint
  {
  public:

    /*!
    \brief Standard Constructor
    */
    MPConstraint3Penalty
    (
        RCP<DRT::Discretization> discr,  ///< Discretization constraint lives on
        const string& CondName ///< Name of condition to create constraint from
    );
    
    /// unused
    void Initialize
    (
      Teuchos::ParameterList&        params,  ///< parameter list to communicate between elements and discretization
      RCP<Epetra_Vector>    systemvector3  ///< distributed vector that may be filled by assembly of element contributions
    );
    
    /// initialization routine called by the manager ctor
    void Initialize
    (
      Teuchos::ParameterList&        params  ///< parameter list to communicate between elements and discretization
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
      RCP<LINALG::SparseOperator> systemmatrix1, ///< sparse matrix that may be filled by assembly of element contributions
      RCP<LINALG::SparseOperator> systemmatrix2, ///< sparse (rectangular) matrix that may be filled by assembly of element contributions
      RCP<Epetra_Vector> systemvector1, ///< distributed vector that may be filled by assembly of element contributions
      RCP<Epetra_Vector> systemvector2, ///< distributed vector that may be filled by assembly of element contributions
      RCP<Epetra_Vector> systemvector3  ///< distributed vector that may be filled by assembly of element contributions
    );

  private:

    // don't want = operator, cctor
    MPConstraint3Penalty operator = (const MPConstraint3Penalty& old);
    MPConstraint3Penalty(const MPConstraint3Penalty& old);

    //! Return the ConstrType based on the condition name
    ConstrType GetConstrType(const string& Name);///<condition name

    //! Evaluate constraint discretization and assemble the results
    void EvaluateConstraint(
        RCP<DRT::Discretization> disc,  ///< discretization to evaluate
        Teuchos::ParameterList&        params,  ///< parameter list to communicate between elements and discretization
        RCP<LINALG::SparseOperator> systemmatrix1, ///< sparse matrix that may be filled by assembly of element contributions
        RCP<LINALG::SparseOperator> systemmatrix2, ///< sparse (rectangular) matrix that may be filled by assembly of element contributions
        RCP<Epetra_Vector>    systemvector1, ///< distributed vector that may be filled by assembly of element contributions
        RCP<Epetra_Vector>    systemvector2, ///< distributed vector that may be filled by assembly of element contributions
        RCP<Epetra_Vector>    systemvector3); ///< distributed vector that may be filled by assembly of element contributions

    //! Initialize constraint discretization and assemble the results to the refbasevector
    void EvaluateError
    (
      RCP<DRT::Discretization> disc,  ///< discretization to evaluate
      Teuchos::ParameterList&        params,  ///< parameter list to communicate between elements and discretization
      RCP<Epetra_Vector>    systemvector3, ///< distributed vector that may be filled by aasembly of element contributions
      bool init =  false
    );

    //! creating a new discretization based on conditions containing constraint elements
    map<int,RCP<DRT::Discretization> > CreateDiscretizationFromCondition
    (
      RCP<DRT::Discretization> actdisc,
      vector< DRT::Condition* >      constrcond, ///< conditions as discretization basis
      const string&             discret_name, ///< name of new discretization
      const string&             element_name, ///< name of element type to create
      int& startID
    );

    // projected attributes
    map<int,bool> absconstraint_;///<maps condition ID to indicator if absolute values are to use for controlling
    map<int,int> eletocondID_;///<maps element ID to condition ID, to allow use of other maps
    map<int,int> eletocondvecindex_;///<maps element ID to condition index in vector #constrcond_
    map<int,double> penalties_; ///< maps condition ID to penalty factor
    RCP<Epetra_Export> errorexport_;
    RCP<Epetra_Import> errorimport_;
    RCP<Epetra_Map> rederrormap_;
    RCP<Epetra_Map> errormap_;
    RCP<Epetra_Vector> initerror_;
    RCP<Epetra_Vector> acterror_;


  }; //class
}
#endif /*MPCONSTRAINT_H*/
