/*----------------------------------------------------------------------*/
/*! \file
\brief Basic constraint class, dealing with multi point constraints
\level 2


*----------------------------------------------------------------------*/

#ifndef BACI_CONSTRAINT_MULTIPOINTCONSTRAINT3_H
#define BACI_CONSTRAINT_MULTIPOINTCONSTRAINT3_H

#include "baci_constraint_multipointconstraint.H"



namespace UTILS
{
  /*!
  \brief This pure virtual class can handle multi point constraints in 3D.
  It is derived from the basic multipointconstraint class.
  */
  class MPConstraint3 : public UTILS::MPConstraint
  {
   public:
    /*!
    \brief Standard Constructor
    */
    MPConstraint3(Teuchos::RCP<DRT::Discretization> discr,  ///< Discretization constraint lives on
        const std::string& conditionname,  ///< Name of condition to create constraint from
        int& offsetID,                     ///< minimum constraint or monitor ID so far
        int& maxID                         ///< maximum constraint or monitor ID so far
    );

    /// initialization routine called by the manager ctor to get correct reference base values and
    /// activating the right conditions at the beginning
    void Initialize(
        Teuchos::ParameterList&
            params,  ///< parameter list to communicate between elements and discretization
        Teuchos::RCP<Epetra_Vector> systemvector3  ///< distributed vector that may be filled by
                                                   ///< assembly of element contributions
        ) override;

    /// initialization routine called at restart to activate the right conditions
    void Initialize(const double& time  ///< current time
        ) override;

    //! Evaluate routine to call from outside. In here the right action is determined and the
    //! #EvaluateConstraint routine is called
    void Evaluate(
        Teuchos::ParameterList&
            params,  ///< parameter list to communicate between elements and discretization
        Teuchos::RCP<CORE::LINALG::SparseOperator>
            systemmatrix1,  ///< sparse matrix that may be filled by assembly of element
                            ///< contributions
        Teuchos::RCP<CORE::LINALG::SparseOperator>
            systemmatrix2,  ///< sparse (rectangular) matrix that may be filled by assembly of
                            ///< element contributions
        Teuchos::RCP<Epetra_Vector> systemvector1,  ///< distributed vector that may be filled by
                                                    ///< assembly of element contributions
        Teuchos::RCP<Epetra_Vector> systemvector2,  ///< distributed vector that may be filled by
                                                    ///< assembly of element contributions
        Teuchos::RCP<Epetra_Vector> systemvector3   ///< distributed vector that may be filled by
                                                    ///< assembly of element contributions
        ) override;

   private:
    // don't want = operator, cctor
    MPConstraint3 operator=(const MPConstraint3& old);
    MPConstraint3(const MPConstraint3& old);

    //! Return the ConstrType based on the condition name
    ConstrType GetConstrType(const std::string& Name);  ///< condition name

    //! Evaluate constraint discretization and assemble the results
    void EvaluateConstraint(Teuchos::RCP<DRT::Discretization> disc,  ///< discretization to evaluate
        Teuchos::ParameterList&
            params,  ///< parameter list to communicate between elements and discretization
        Teuchos::RCP<CORE::LINALG::SparseOperator>
            systemmatrix1,  ///< sparse matrix that may be filled by assembly of element
                            ///< contributions
        Teuchos::RCP<CORE::LINALG::SparseOperator>
            systemmatrix2,  ///< sparse (rectangular) matrix that may be filled by assembly of
                            ///< element contributions
        Teuchos::RCP<Epetra_Vector> systemvector1,  ///< distributed vector that may be filled by
                                                    ///< assembly of element contributions
        Teuchos::RCP<Epetra_Vector> systemvector2,  ///< distributed vector that may be filled by
                                                    ///< assembly of element contributions
        Teuchos::RCP<Epetra_Vector> systemvector3)
        override;  ///< distributed vector that may be filled by
                   ///< assembly of element contributions

    //! Initialize constraint discretization and assemble the results to the refbasevector
    void InitializeConstraint(
        Teuchos::RCP<DRT::Discretization> disc,  ///< discretization to evaluate
        Teuchos::ParameterList&
            params,  ///< parameter list to communicate between elements and discretization
        Teuchos::RCP<Epetra_Vector> systemvector3  ///< distributed vector that may be filled by
                                                   ///< aasembly of element contributions
    );

    //! creating a new discretization based on conditions containing constraint elements
    std::map<int, Teuchos::RCP<DRT::Discretization>> CreateDiscretizationFromCondition(
        Teuchos::RCP<DRT::Discretization> actdisc,
        std::vector<DRT::Condition*> constrcond,  ///< conditions as discretization basis
        const std::string& discret_name,          ///< name of new discretization
        const std::string& element_name,          ///< name of element type to create
        int& startID                              ///< ID to start with
        ) override;

    // projected attributes
    std::map<int, bool> absconstraint_;  ///< maps condition ID to indicator if absolute values are
                                         ///< to use for controlling
    std::map<int, int>
        eletocondID_;  ///< maps element ID to condition ID, to allow use of other maps
    std::map<int, int>
        eletocondvecindex_;  ///< maps element ID to condition index in vector #constrcond_

  };  // class
}  // namespace UTILS
#endif  // CONSTRAINT_MULTIPOINTCONSTRAINT3_H
