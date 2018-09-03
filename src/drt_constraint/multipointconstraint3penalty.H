/*!----------------------------------------------------------------------
\file multipointconstraint3penalty.H
\brief Basic constraint class, dealing with multi point constraints
\level 2
<pre>

\maintainer Alexander Popp
            popp@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>

*----------------------------------------------------------------------*/

#ifndef MPCONSTRAINT3PENALTY_H
#define MPCONSTRAINT3PENALTY_H

#include "multipointconstraint.H"



namespace UTILS
{
  /*!
  \brief This class can handle multi point constraints in 3D.
  It is derived from the basic multipointconstraint class.
  */
  class MPConstraint3Penalty : public UTILS::MPConstraint
  {
   public:
    /*!
    \brief Standard Constructor
    */
    MPConstraint3Penalty(
        Teuchos::RCP<DRT::Discretization> discr,  ///< Discretization constraint lives on
        const std::string& CondName               ///< Name of condition to create constraint from
    );

    /// unused
    void Initialize(
        Teuchos::ParameterList&
            params,  ///< parameter list to communicate between elements and discretization
        Teuchos::RCP<Epetra_Vector> systemvector3  ///< distributed vector that may be filled by
                                                   ///< assembly of element contributions
    );

    /// initialization routine called by the manager ctor
    void Initialize(Teuchos::ParameterList&
            params  ///< parameter list to communicate between elements and discretization
    );

    /// initialization routine called at restart to activate the right conditions
    void Initialize(const double& time  ///< current time
    );

    //! Evaluate routine to call from outside. In here the right action is determined and the
    //! #EvaluateConstraint routine is called
    void Evaluate(
        Teuchos::ParameterList&
            params,  ///< parameter list to communicate between elements and discretization
        Teuchos::RCP<LINALG::SparseOperator>
            systemmatrix1,  ///< sparse matrix that may be filled by assembly of element
                            ///< contributions
        Teuchos::RCP<LINALG::SparseOperator>
            systemmatrix2,  ///< sparse (rectangular) matrix that may be filled by assembly of
                            ///< element contributions
        Teuchos::RCP<Epetra_Vector> systemvector1,  ///< distributed vector that may be filled by
                                                    ///< assembly of element contributions
        Teuchos::RCP<Epetra_Vector> systemvector2,  ///< distributed vector that may be filled by
                                                    ///< assembly of element contributions
        Teuchos::RCP<Epetra_Vector> systemvector3   ///< distributed vector that may be filled by
                                                    ///< assembly of element contributions
    );

   private:
    // don't want = operator, cctor
    MPConstraint3Penalty operator=(const MPConstraint3Penalty& old);
    MPConstraint3Penalty(const MPConstraint3Penalty& old);

    //! Return the ConstrType based on the condition name
    ConstrType GetConstrType(const std::string& Name);  ///< condition name

    //! Evaluate constraint discretization and assemble the results
    void EvaluateConstraint(Teuchos::RCP<DRT::Discretization> disc,  ///< discretization to evaluate
        Teuchos::ParameterList&
            params,  ///< parameter list to communicate between elements and discretization
        Teuchos::RCP<LINALG::SparseOperator>
            systemmatrix1,  ///< sparse matrix that may be filled by assembly of element
                            ///< contributions
        Teuchos::RCP<LINALG::SparseOperator>
            systemmatrix2,  ///< sparse (rectangular) matrix that may be filled by assembly of
                            ///< element contributions
        Teuchos::RCP<Epetra_Vector> systemvector1,   ///< distributed vector that may be filled by
                                                     ///< assembly of element contributions
        Teuchos::RCP<Epetra_Vector> systemvector2,   ///< distributed vector that may be filled by
                                                     ///< assembly of element contributions
        Teuchos::RCP<Epetra_Vector> systemvector3);  ///< distributed vector that may be filled by
                                                     ///< assembly of element contributions

    //! Initialize constraint discretization and assemble the results to the refbasevector
    void EvaluateError(Teuchos::RCP<DRT::Discretization> disc,  ///< discretization to evaluate
        Teuchos::ParameterList&
            params,  ///< parameter list to communicate between elements and discretization
        Teuchos::RCP<Epetra_Vector> systemvector3,  ///< distributed vector that may be filled by
                                                    ///< aasembly of element contributions
        bool init = false);

    //! creating a new discretization based on conditions containing constraint elements
    std::map<int, Teuchos::RCP<DRT::Discretization>> CreateDiscretizationFromCondition(
        Teuchos::RCP<DRT::Discretization> actdisc,
        std::vector<DRT::Condition*> constrcond,  ///< conditions as discretization basis
        const std::string& discret_name,          ///< name of new discretization
        const std::string& element_name,          ///< name of element type to create
        int& startID);

    // projected attributes
    std::map<int, bool> absconstraint_;  ///< maps condition ID to indicator if absolute values are
                                         ///< to use for controlling
    std::map<int, int>
        eletocondID_;  ///< maps element ID to condition ID, to allow use of other maps
    std::map<int, int>
        eletocondvecindex_;  ///< maps element ID to condition index in vector #constrcond_
    std::map<int, double> penalties_;  ///< maps condition ID to penalty factor
    Teuchos::RCP<Epetra_Export> errorexport_;
    Teuchos::RCP<Epetra_Import> errorimport_;
    Teuchos::RCP<Epetra_Map> rederrormap_;
    Teuchos::RCP<Epetra_Map> errormap_;
    Teuchos::RCP<Epetra_Vector> initerror_;
    Teuchos::RCP<Epetra_Vector> acterror_;


  };  // class
}  // namespace UTILS
#endif /*MPCONSTRAINT_H*/
