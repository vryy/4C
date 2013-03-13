/*!----------------------------------------------------------------------
\file multipointconstraint3.H

\brief Basic constraint class, dealing with multi point constraints

<pre>
Maintainer: Thomas Kloeppel
            kloeppel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>

*----------------------------------------------------------------------*/

#ifndef MPCONSTRAINT3_H
#define MPCONSTRAINT3_H

#include "multipointconstraint.H"

using namespace Teuchos;



namespace UTILS
{
  /*!
  \brief This pure virtual class can handle multi point constraints in 3D.
  It is derived from the basic multipointconstraint class.

  \author Kloeppel (kloeppel@lnm.mw.tum.de)
  */
  class MPConstraint3 : public UTILS::MPConstraint
  {
  public:

    /*!
    \brief Standard Constructor
    */
    MPConstraint3
    (
      RCP<DRT::Discretization> discr,  ///< Discretization constraint lives on
      const string& CondName, ///< Name of condition to create constraint from
      int& minID,  ///< minimum constraint or monitor ID so far
      int& maxID  ///< maximum constraint or monitor ID so far
    );

    /// initialization routine called by the manager ctor to get correct reference base values and activating the right conditions at the beginning
    void Initialize
    (
      Teuchos::ParameterList&        params,  ///< parameter list to communicate between elements and discretization
      RCP<Epetra_Vector>    systemvector3  ///< distributed vector that may be filled by assembly of element contributions
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
    MPConstraint3 operator = (const MPConstraint3& old);
    MPConstraint3(const MPConstraint3& old);

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
    void InitializeConstraint
    (
      RCP<DRT::Discretization> disc,  ///< discretization to evaluate
      Teuchos::ParameterList&        params,  ///< parameter list to communicate between elements and discretization
      RCP<Epetra_Vector>    systemvector3 ///< distributed vector that may be filled by aasembly of element contributions
    );

    //! creating a new discretization based on conditions containing constraint elements
    map<int,RCP<DRT::Discretization> > CreateDiscretizationFromCondition
    (
      RCP<DRT::Discretization>       actdisc,
      std::vector< DRT::Condition* > constrcond, ///< conditions as discretization basis
      const std::string&             discret_name, ///< name of new discretization
      const std::string&             element_name, ///< name of element type to create
      int& startID  ///<ID to start with
    );

    // projected attributes
    map<int,bool> absconstraint_;///<maps condition ID to indicator if absolute values are to use for controlling
    map<int,int> eletocondID_;///<maps element ID to condition ID, to allow use of other maps
    map<int,int> eletocondvecindex_;///<maps element ID to condition index in vector #constrcond_

  }; //class
}
#endif /*MPCONSTRAINT_H*/
