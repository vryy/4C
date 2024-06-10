/*----------------------------------------------------------------------*/
/*! \file
\brief Basic constraint class, dealing with multi point constraints
\level 2


*----------------------------------------------------------------------*/

#ifndef FOUR_C_CONSTRAINT_MULTIPOINTCONSTRAINT_HPP
#define FOUR_C_CONSTRAINT_MULTIPOINTCONSTRAINT_HPP

#include "4C_config.hpp"

#include "4C_constraint.hpp"

FOUR_C_NAMESPACE_OPEN



namespace CONSTRAINTS
{
  /*!
  \brief This pure virtual class is the common interface for multi point constraints.
  It is derived from the basic constraint class.
  */
  class MPConstraint : public Constraint
  {
   public:
    /*!
    \brief Standard Constructor
    */
    MPConstraint(
        Teuchos::RCP<Core::FE::Discretization> discr,  ///< discretization constraint lives on
        const std::string& conditionname,  ///< Name of condition to create constraint from
        int& minID,                        ///< minimum constraint or monitor ID so far
        int& maxID                         ///< maximum constraint or monitor ID so far
    );

    /*!
        \brief Alternative Constructor
    */
    MPConstraint(
        Teuchos::RCP<Core::FE::Discretization> discr,  ///< discretization constraint lives on
        const std::string& conditionname  ///< Name of condition to create constraint from
    );

    /*!
        \brief Destructor
    */
    virtual ~MPConstraint() { ; };

    /// Set state of the underlying constraint discretization
    void SetConstrState(const std::string& state,  ///< name of state to set
        Teuchos::RCP<const Epetra_Vector> V        ///< values to set
    );

    /// initialization routine called by the manager ctor to get correct reference base values and
    /// activating the right conditions at the beginning
    virtual void Initialize(
        Teuchos::ParameterList&
            params,  ///< parameter list to communicate between elements and discretization
        Teuchos::RCP<Epetra_Vector> systemvector3  ///< distributed vector that may be filled by
                                                   ///< assembly of element contributions
        ) = 0;

    /// initialization routine called at restart to activate the right conditions
    virtual void Initialize(const double& time  ///< current time
        ) = 0;

    //! Evaluate routine to call from outside. In here the right action is determined and the
    //! #evaluate_constraint routine is called
    virtual void Evaluate(
        Teuchos::ParameterList&
            params,  ///< parameter list to communicate between elements and discretization
        Teuchos::RCP<Core::LinAlg::SparseOperator>
            systemmatrix1,  ///< sparse matrix that may be filled by assembly of element
                            ///< contributions
        Teuchos::RCP<Core::LinAlg::SparseOperator>
            systemmatrix2,  ///< sparse (rectangular) matrix that may be filled by assembly of
                            ///< element contributions
        Teuchos::RCP<Epetra_Vector> systemvector1,  ///< distributed vector that may be filled by
                                                    ///< assembly of element contributions
        Teuchos::RCP<Epetra_Vector> systemvector2,  ///< distributed vector that may be filled by
                                                    ///< assembly of element contributions
        Teuchos::RCP<Epetra_Vector> systemvector3   ///< distributed vector that may be filled by
                                                    ///< assembly of element contributions
        ) = 0;

    //! Is there a constraint defined in this class?
    bool HaveConstraint() { return constrtype_ != none; };

   protected:
    // cctor
    MPConstraint(const MPConstraint& old);


    //! additional discretization consisting of constraint elements
    std::map<int, Teuchos::RCP<Core::FE::Discretization>> constraintdis_;

    //! Evaluate constraint discretization and assemble the results
    virtual void evaluate_constraint(
        Teuchos::RCP<Core::FE::Discretization> disc,  ///< discretization to evaluate
        Teuchos::ParameterList&
            params,  ///< parameter list to communicate between elements and discretization
        Teuchos::RCP<Core::LinAlg::SparseOperator>
            systemmatrix1,  ///< sparse matrix that may be filled by assembly of element
                            ///< contributions
        Teuchos::RCP<Core::LinAlg::SparseOperator>
            systemmatrix2,  ///< sparse (rectangular) matrix that may be filled by assembly of
                            ///< element contributions
        Teuchos::RCP<Epetra_Vector> systemvector1,  ///< distributed vector that may be filled by
                                                    ///< assembly of element contributions
        Teuchos::RCP<Epetra_Vector> systemvector2,  ///< distributed vector that may be filled by
                                                    ///< assembly of element contributions
        Teuchos::RCP<Epetra_Vector> systemvector3   ///< distributed vector that may be filled by
                                                    ///< assembly of element contributions
        ) = 0;

    //! creating a new discretization based on conditions containing constraint elements
    virtual std::map<int, Teuchos::RCP<Core::FE::Discretization>>
    create_discretization_from_condition(Teuchos::RCP<Core::FE::Discretization> actdisc,
        std::vector<Core::Conditions::Condition*>
            constrcond,                   ///< conditions as discretization basis
        const std::string& discret_name,  ///< name of new discretization
        const std::string& element_name,  ///< name of element type to create
        int& startID) = 0;


    //    /// find col node map so that we can evaluate the constraint elements
    //    Teuchos::RCP<Epetra_Map> ComputeNodeColMap(
    //             const Teuchos::RCP<Core::FE::Discretization> sourcedis,  ///< standard
    //             discretization we want to redistribute const
    //             Teuchos::RCP<Core::FE::Discretization> constraintdis
    //             ///< constraint discretization prescribing ghosting ) const;

  };  // class
}  // namespace CONSTRAINTS
FOUR_C_NAMESPACE_CLOSE

#endif
