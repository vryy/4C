/*----------------------------------------------------------------------*/
/*! \file

\brief Scatra-scatra interface coupling strategy for electrochemistry problems

\level 2


*----------------------------------------------------------------------*/
#ifndef FOUR_C_SCATRA_TIMINT_MESHTYING_STRATEGY_S2I_ELCH_HPP
#define FOUR_C_SCATRA_TIMINT_MESHTYING_STRATEGY_S2I_ELCH_HPP

#include "4C_config.hpp"

#include "4C_lib_element.hpp"
#include "4C_scatra_timint_elch.hpp"
#include "4C_scatra_timint_meshtying_strategy_s2i.hpp"

FOUR_C_NAMESPACE_OPEN

namespace SCATRA
{
  /*!
  \brief Scatra-scatra interface coupling strategy for electrochemistry problems

  To keep the scalar transport time integrator class and derived classes as plain as possible,
  several algorithmic parts have been encapsulated within separate meshtying strategy classes.
  These algorithmic parts include initializing the system matrix and other relevant objects,
  computing meshtying residual terms and their linearizations, and solving the resulting
  linear system of equations. By introducing a hierarchy of strategies for these algorithmic
  parts, a bunch of unhandy if-else selections within the time integrator classes themselves
  can be circumvented. This class contains the scatra-scatra interface coupling strategy for
  electrochemistry problems.

  */

  class MeshtyingStrategyS2IElch : public MeshtyingStrategyS2I
  {
   public:
    //! constructor
    explicit MeshtyingStrategyS2IElch(
        SCATRA::ScaTraTimIntElch* elchtimint,  //!< elch time integrator
        const Teuchos::ParameterList&
            parameters  //!< input parameters for scatra-scatra interface coupling
    );

    //! compute time step size
    void ComputeTimeStepSize(double& dt) final;

    //! compute meshtying residual terms and their linearizations
    void EvaluateMeshtying() override;

    void evaluate_point_coupling() override;

    //! update solution after convergence of the nonlinear Newton-Raphson iteration
    void Update() const override;

   private:
    //! copy constructor
    MeshtyingStrategyS2IElch(const MeshtyingStrategyS2IElch& old);

    //! return pointer to elch time integrator after cast
    SCATRA::ScaTraTimIntElch* ElchTimInt() const
    {
      return dynamic_cast<SCATRA::ScaTraTimIntElch*>(scatratimint_);
    };

    //! instantiate strategy for Newton-Raphson convergence check
    void init_conv_check_strategy() override;

    //! minimum interfacial overpotential associated with scatra-scatra interface layer growth
    double etagrowthmin_;

    //! time step at the start of scatra-scatra interface layer growth
    int intlayergrowth_startstep_;

    //! flag indicating modified time step size for scatra-scatra interface layer growth
    bool intlayergrowth_timestep_active_;
  };  // class MeshtyingStrategyS2IElch

  class MeshtyingStrategyS2IElchSCL : public MeshtyingStrategyS2IElch
  {
   public:
    explicit MeshtyingStrategyS2IElchSCL(
        SCATRA::ScaTraTimIntElch* elchtimint, const Teuchos::ParameterList& parameters);

    void add_time_integration_specific_vectors() const override{};

    void EvaluateMeshtying() override{};

    void SetupMeshtying() override;

    void Solve(const Teuchos::RCP<CORE::LINALG::Solver>& solver,
        const Teuchos::RCP<CORE::LINALG::SparseOperator>& systemmatrix,
        const Teuchos::RCP<Epetra_Vector>& increment, const Teuchos::RCP<Epetra_Vector>& residual,
        const Teuchos::RCP<Epetra_Vector>& phinp, const int iteration,
        CORE::LINALG::SolverParams& solver_params) const override;
  };

  template <CORE::FE::CellType distypeS, CORE::FE::CellType distypeM>
  class MortarCellCalcElch : public MortarCellCalc<distypeS, distypeM>
  {
   public:
    //! singleton access method
    static MortarCellCalcElch<distypeS, distypeM>* Instance(
        const INPAR::S2I::CouplingType& couplingtype,  //!< flag for meshtying method
        const INPAR::S2I::InterfaceSides&
            lmside,  //!< flag for interface side underlying Lagrange multiplier definition
        const int& numdofpernode_slave,   //!< number of slave-side degrees of freedom per node
        const int& numdofpernode_master,  //!< number of master-side degrees of freedom per node
        const std::string& disname        //!< name of mortar discretization
    );


   protected:
    //! protected constructor for singletons
    MortarCellCalcElch(const INPAR::S2I::CouplingType& couplingtype,  //!< flag for meshtying method
        const INPAR::S2I::InterfaceSides&
            lmside,  //!< flag for interface side underlying Lagrange multiplier definition
        const int& numdofpernode_slave,  //!< number of slave-side degrees of freedom per node
        const int& numdofpernode_master  //!< number of master-side degrees of freedom per node
    );

   private:
    //! abbreviation
    typedef MortarCellCalc<distypeS, distypeM> my;

   protected:
    using my::nen_master_;
    using my::nen_slave_;
    using my::nsd_slave_;

   private:
    //! evaluate and assemble interface linearizations and residuals
    void EvaluateCondition(const DRT::Discretization& idiscret,  //!< interface discretization
        MORTAR::IntCell& cell,                                   //!< mortar integration cell
        MORTAR::Element& slaveelement,                           //!< slave-side mortar element
        MORTAR::Element& masterelement,                          //!< master-side mortar element
        DRT::Element::LocationArray& la_slave,                   //!< slave-side location array
        DRT::Element::LocationArray& la_master,                  //!< master-side location array
        const Teuchos::ParameterList& params,                    //!< parameter list
        CORE::LINALG::SerialDenseMatrix&
            k_ss,  //!< linearizations of slave-side residuals w.r.t. slave-side dofs
        CORE::LINALG::SerialDenseMatrix&
            k_sm,  //!< linearizations of slave-side residuals w.r.t. master-side dofs
        CORE::LINALG::SerialDenseMatrix&
            k_ms,  //!< linearizations of master-side residuals w.r.t. slave-side dofs
        CORE::LINALG::SerialDenseMatrix&
            k_mm,  //!< linearizations of master-side residuals w.r.t. master-side dofs
        CORE::LINALG::SerialDenseVector& r_s,  //!< slave-side residual vector
        CORE::LINALG::SerialDenseVector& r_m   //!< master-side residual vector
        ) override;

    //! evaluate and assemble interface linearizations and residuals for node-to-segment coupling
    void evaluate_condition_nts(
        CORE::Conditions::Condition& condition,  //!< scatra-scatra interface coupling condition
        const MORTAR::Node& slavenode,           //!< slave-side node
        const double&
            lumpedarea,  //!< lumped interface area fraction associated with slave-side node
        MORTAR::Element& slaveelement,   //!< slave-side mortar element
        MORTAR::Element& masterelement,  //!< master-side mortar element
        const std::vector<CORE::LINALG::Matrix<nen_slave_, 1>>&
            ephinp_slave,  //!< state variables at slave-side nodes
        const std::vector<CORE::LINALG::Matrix<nen_master_, 1>>&
            ephinp_master,  //!< state variables at master-side nodes
        CORE::LINALG::SerialDenseMatrix&
            k_ss,  //!< linearizations of slave-side residuals w.r.t. slave-side dofs
        CORE::LINALG::SerialDenseMatrix&
            k_sm,  //!< linearizations of slave-side residuals w.r.t. master-side dofs
        CORE::LINALG::SerialDenseMatrix&
            k_ms,  //!< linearizations of master-side residuals w.r.t. slave-side dofs
        CORE::LINALG::SerialDenseMatrix&
            k_mm,  //!< linearizations of master-side residuals w.r.t. master-side dofs
        CORE::LINALG::SerialDenseVector& r_s,  //!< slave-side residual vector
        CORE::LINALG::SerialDenseVector& r_m   //!< master-side residual vector
        ) override;

    //! evaluate factor F/RT
    virtual double GetFRT() const;
  };  // class MortarCellCalcElch


  template <CORE::FE::CellType distypeS, CORE::FE::CellType distypeM>
  class MortarCellCalcElchSTIThermo : public MortarCellCalcElch<distypeS, distypeM>
  {
   public:
    //! singleton access method
    static MortarCellCalcElchSTIThermo<distypeS, distypeM>* Instance(
        const INPAR::S2I::CouplingType& couplingtype,  //!< flag for meshtying method
        const INPAR::S2I::InterfaceSides&
            lmside,  //!< flag for interface side underlying Lagrange multiplier definition
        const int& numdofpernode_slave,   //!< number of slave-side degrees of freedom per node
        const int& numdofpernode_master,  //!< number of master-side degrees of freedom per node
        const std::string& disname        //!< name of mortar discretization
    );


    //! evaluate single mortar integration cell of particular slave-side and master-side
    //! discretization types
    void Evaluate(const DRT::Discretization& idiscret,  //!< interface discretization
        MORTAR::IntCell& cell,                          //!< mortar integration cell
        MORTAR::Element& slaveelement,                  //!< slave-side mortar element
        MORTAR::Element& masterelement,                 //!< master-side mortar element
        DRT::Element::LocationArray& la_slave,          //!< slave-side location array
        DRT::Element::LocationArray& la_master,         //!< master-side location array
        const Teuchos::ParameterList& params,           //!< parameter list
        CORE::LINALG::SerialDenseMatrix& cellmatrix1,   //!< cell matrix 1
        CORE::LINALG::SerialDenseMatrix& cellmatrix2,   //!< cell matrix 2
        CORE::LINALG::SerialDenseMatrix& cellmatrix3,   //!< cell matrix 3
        CORE::LINALG::SerialDenseMatrix& cellmatrix4,   //!< cell matrix 4
        CORE::LINALG::SerialDenseVector& cellvector1,   //!< cell vector 1
        CORE::LINALG::SerialDenseVector& cellvector2    //!< cell vector 2
        ) override;

   private:
    //! abbreviation
    typedef MortarCellCalc<distypeS, distypeM> my;
    typedef MortarCellCalcElch<distypeS, distypeM> myelch;
    using my::nen_master_;
    using my::nen_slave_;
    using my::nsd_slave_;

    //! private constructor for singletons
    MortarCellCalcElchSTIThermo(
        const INPAR::S2I::CouplingType& couplingtype,  //!< flag for meshtying method
        const INPAR::S2I::InterfaceSides&
            lmside,  //!< flag for interface side underlying Lagrange multiplier definition
        const int& numdofpernode_slave,  //!< number of slave-side degrees of freedom per node
        const int& numdofpernode_master  //!< number of master-side degrees of freedom per node
    );

    //! evaluate and assemble off-diagonal interface linearizations
    void EvaluateConditionOD(const DRT::Discretization& idiscret,  //!< interface discretization
        MORTAR::IntCell& cell,                                     //!< mortar integration cell
        MORTAR::Element& slaveelement,                             //!< slave-side mortar element
        MORTAR::Element& masterelement,                            //!< master-side mortar element
        DRT::Element::LocationArray& la_slave,                     //!< slave-side location array
        DRT::Element::LocationArray& la_master,                    //!< master-side location array
        const Teuchos::ParameterList& params,                      //!< parameter list
        CORE::LINALG::SerialDenseMatrix&
            k_ss,  //!< linearizations of slave-side residuals w.r.t. slave-side dofs
        CORE::LINALG::SerialDenseMatrix&
            k_ms  //!< linearizations of master-side residuals w.r.t. slave-side dofs
    );

    //! extract nodal state variables associated with mortar integration cell
    void ExtractNodeValues(const DRT::Discretization& idiscret,  //!< interface discretization
        DRT::Element::LocationArray& la_slave,                   //!< slave-side location array
        DRT::Element::LocationArray& la_master                   //!< master-side location array
        ) override;

    //! evaluate factor F/RT
    double GetFRT() const override;

    //! nodal, slave-side temperature variables associated with time t_{n+1} or t_{n+alpha_f}
    CORE::LINALG::Matrix<nen_slave_, 1> etempnp_slave_;
  };  // class MortarCellCalcElchSTIThermo


  template <CORE::FE::CellType distypeS, CORE::FE::CellType distypeM>
  class MortarCellCalcSTIElch : public MortarCellCalc<distypeS, distypeM>
  {
   public:
    //! singleton access method
    static MortarCellCalcSTIElch<distypeS, distypeM>* Instance(
        const INPAR::S2I::CouplingType& couplingtype,  //!< flag for meshtying method
        const INPAR::S2I::InterfaceSides&
            lmside,  //!< flag for interface side underlying Lagrange multiplier definition
        const int& numdofpernode_slave,   //!< number of slave-side degrees of freedom per node
        const int& numdofpernode_master,  //!< number of master-side degrees of freedom per node
        const std::string& disname        //!< name of mortar discretization
    );

    //! evaluate single mortar integration cell of particular slave-side and master-side
    //! discretization types
    void Evaluate(const DRT::Discretization& idiscret,  //!< interface discretization
        MORTAR::IntCell& cell,                          //!< mortar integration cell
        MORTAR::Element& slaveelement,                  //!< slave-side mortar element
        MORTAR::Element& masterelement,                 //!< master-side mortar element
        DRT::Element::LocationArray& la_slave,          //!< slave-side location array
        DRT::Element::LocationArray& la_master,         //!< master-side location array
        const Teuchos::ParameterList& params,           //!< parameter list
        CORE::LINALG::SerialDenseMatrix& cellmatrix1,   //!< cell matrix 1
        CORE::LINALG::SerialDenseMatrix& cellmatrix2,   //!< cell matrix 2
        CORE::LINALG::SerialDenseMatrix& cellmatrix3,   //!< cell matrix 3
        CORE::LINALG::SerialDenseMatrix& cellmatrix4,   //!< cell matrix 4
        CORE::LINALG::SerialDenseVector& cellvector1,   //!< cell vector 1
        CORE::LINALG::SerialDenseVector& cellvector2    //!< cell vector 2
        ) override;

   private:
    //! abbreviation
    typedef MortarCellCalc<distypeS, distypeM> my;
    using my::nen_master_;
    using my::nen_slave_;
    using my::nsd_slave_;

    //! private constructor for singletons
    MortarCellCalcSTIElch(
        const INPAR::S2I::CouplingType& couplingtype,  //!< flag for meshtying method
        const INPAR::S2I::InterfaceSides&
            lmside,  //!< flag for interface side underlying Lagrange multiplier definition
        const int& numdofpernode_slave,  //!< number of slave-side degrees of freedom per node
        const int& numdofpernode_master  //!< number of master-side degrees of freedom per node
    );

    //! evaluate and assemble interface linearizations and residuals
    void EvaluateCondition(const DRT::Discretization& idiscret,  //!< interface discretization
        MORTAR::IntCell& cell,                                   //!< mortar integration cell
        MORTAR::Element& slaveelement,                           //!< slave-side mortar element
        MORTAR::Element& masterelement,                          //!< master-side mortar element
        DRT::Element::LocationArray& la_slave,                   //!< slave-side location array
        DRT::Element::LocationArray& la_master,                  //!< master-side location array
        const Teuchos::ParameterList& params,                    //!< parameter list
        CORE::LINALG::SerialDenseMatrix&
            k_ss,  //!< linearizations of slave-side residuals w.r.t. slave-side dofs
        CORE::LINALG::SerialDenseVector& r_s  //!< slave-side residual vector
    );

    //! evaluate and assemble off-diagonal interface linearizations
    void EvaluateConditionOD(const DRT::Discretization& idiscret,  //!< interface discretization
        MORTAR::IntCell& cell,                                     //!< mortar integration cell
        MORTAR::Element& slaveelement,                             //!< slave-side mortar element
        MORTAR::Element& masterelement,                            //!< master-side mortar element
        DRT::Element::LocationArray& la_slave,                     //!< slave-side location array
        DRT::Element::LocationArray& la_master,                    //!< master-side location array
        const Teuchos::ParameterList& params,                      //!< parameter list
        CORE::LINALG::SerialDenseMatrix&
            k_ss,  //!< linearizations of slave-side residuals w.r.t. slave-side dofs
        CORE::LINALG::SerialDenseMatrix&
            k_sm  //!< linearizations of slave-side residuals w.r.t. master-side dofs
    );

    //! extract nodal state variables associated with mortar integration cell
    void ExtractNodeValues(const DRT::Discretization& idiscret,  //!< interface discretization
        DRT::Element::LocationArray& la_slave,                   //!< slave-side location array
        DRT::Element::LocationArray& la_master                   //!< master-side location array
        ) override;

    //! nodal, slave-side electrochemistry variables associated with time t_{n+1} or t_{n+alpha_f}
    std::vector<CORE::LINALG::Matrix<nen_slave_, 1>> eelchnp_slave_;

    //! nodal, master-side electrochemistry variables associated with time t_{n+1} or t_{n+alpha_f}
    std::vector<CORE::LINALG::Matrix<nen_master_, 1>> eelchnp_master_;
  };  // class MortarCellCalcSTIElch
}  // namespace SCATRA
FOUR_C_NAMESPACE_CLOSE

#endif
