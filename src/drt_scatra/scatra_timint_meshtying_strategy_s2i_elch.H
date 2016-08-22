/*!----------------------------------------------------------------------
\file scatra_timint_meshtying_strategy_s2i_elch.H

\brief Scatra-scatra interface coupling strategy for electrochemistry problems

\level 2

<pre>
\maintainer Rui Fang
            fang@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089 - 289-15251
</pre>

*----------------------------------------------------------------------*/
#ifndef SCATRA_TIMINT_MESHTYING_STRATEGY_S2I_ELCH_H
#define SCATRA_TIMINT_MESHTYING_STRATEGY_S2I_ELCH_H

#include "../drt_scatra/scatra_timint_elch.H"

#include "scatra_timint_meshtying_strategy_s2i.H"

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

\author fang (fang@lnm.mw.tum.de)
*/

class MeshtyingStrategyS2IElch : public MeshtyingStrategyS2I
{
  public:

    //! constructor
    explicit MeshtyingStrategyS2IElch(
        SCATRA::ScaTraTimIntElch*       elchtimint,   //!< elch time integrator
        const Teuchos::ParameterList&   parameters    //!< input parameters for scatra-scatra interface coupling
        );

    //! destructor
    virtual ~MeshtyingStrategyS2IElch(){return;};

    //! compute meshtying residual terms and their linearizations
    void EvaluateMeshtying();

  private:

    //! copy constructor
    MeshtyingStrategyS2IElch(const MeshtyingStrategyS2IElch& old);

    //! return pointer to elch time integrator after cast
    SCATRA::ScaTraTimIntElch* ElchTimInt() const {return dynamic_cast<SCATRA::ScaTraTimIntElch*>(scatratimint_);};

    //! build maps associated with blocks of global system matrix
    void BuildBlockMaps(
        const std::vector<Teuchos::RCP<DRT::Condition> >&   partitioningconditions,   //!< domain partitioning conditions
        std::vector<Teuchos::RCP<const Epetra_Map> >&       blockmaps                 //!< empty vector for maps to be built
        ) const;

    //! build null spaces associated with blocks of global system matrix
    void BuildBlockNullSpaces() const;

    //! instantiate strategy for Newton-Raphson convergence check
    void InitConvCheckStrategy();
}; // class MeshtyingStrategyS2IElch


template<DRT::Element::DiscretizationType distypeS,DRT::Element::DiscretizationType distypeM>
class MortarCellCalcElch : public MortarCellCalc<distypeS,distypeM>
{
  public:

    //! singleton access method
    static MortarCellCalcElch<distypeS,distypeM>* Instance(
        const INPAR::S2I::CouplingType&     couplingtype,           //!< flag for meshtying method
        const INPAR::S2I::InterfaceSides&   lmside,                 //!< flag for interface side underlying Lagrange multiplier definition
        const int&                          numdofpernode_slave,    //!< number of slave-side degrees of freedom per node
        const int&                          numdofpernode_master,   //!< number of master-side degrees of freedom per node
        bool                                create = true           //!< creation flag
        );

    //! singleton destruction
    virtual void Done();

  private:

    //! abbreviation
    typedef MortarCellCalc<distypeS,distypeM> my;

    //! private constructor for singletons
    MortarCellCalcElch(
        const INPAR::S2I::CouplingType&     couplingtype,          //!< flag for meshtying method
        const INPAR::S2I::InterfaceSides&   lmside,                //!< flag for interface side underlying Lagrange multiplier definition
        const int&                          numdofpernode_slave,   //!< number of slave-side degrees of freedom per node
        const int&                          numdofpernode_master   //!< number of master-side degrees of freedom per node
        );

    //! private destructor for singletons
    virtual ~MortarCellCalcElch(){return;};

    //! evaluate and assemble interface linearizations and residuals
    void EvaluateCondition(
        DRT::Condition&                                          condition,       //!< scatra-scatra interface coupling condition
        MORTAR::IntCell&                                         cell,            //!< mortar integration cell
        MORTAR::MortarElement&                                   slaveelement,    //!< slave-side mortar element
        MORTAR::MortarElement&                                   masterelement,   //!< master-side mortar element
        const std::vector<LINALG::Matrix<my::nen_slave_,1> >&    ephinp_slave,    //!< state variables at slave-side nodes
        const std::vector<LINALG::Matrix<my::nen_master_,1> >&   ephinp_master,   //!< state variables at master-side nodes
        Epetra_SerialDenseMatrix&                                k_ss,            //!< linearizations of slave-side residuals w.r.t. slave-side dofs
        Epetra_SerialDenseMatrix&                                k_sm,            //!< linearizations of slave-side residuals w.r.t. master-side dofs
        Epetra_SerialDenseMatrix&                                k_ms,            //!< linearizations of master-side residuals w.r.t. slave-side dofs
        Epetra_SerialDenseMatrix&                                k_mm,            //!< linearizations of master-side residuals w.r.t. master-side dofs
        Epetra_SerialDenseVector&                                r_s,             //!< slave-side residual vector
        Epetra_SerialDenseVector&                                r_m              //!< master-side residual vector
        );

    //! evaluate and assemble interface linearizations and residuals for node-to-segment coupling
    void EvaluateConditionNTS(
        DRT::Condition&                                          condition,       //!< scatra-scatra interface coupling condition
        const MORTAR::MortarNode&                                slavenode,       //!< slave-side node
        const double&                                            lumpedarea,      //!< lumped interface area fraction associated with slave-side node
        MORTAR::MortarElement&                                   slaveelement,    //!< slave-side mortar element
        MORTAR::MortarElement&                                   masterelement,   //!< master-side mortar element
        const std::vector<LINALG::Matrix<my::nen_slave_,1> >&    ephinp_slave,    //!< state variables at slave-side nodes
        const std::vector<LINALG::Matrix<my::nen_master_,1> >&   ephinp_master,   //!< state variables at master-side nodes
        Epetra_SerialDenseMatrix&                                k_ss,            //!< linearizations of slave-side residuals w.r.t. slave-side dofs
        Epetra_SerialDenseMatrix&                                k_sm,            //!< linearizations of slave-side residuals w.r.t. master-side dofs
        Epetra_SerialDenseMatrix&                                k_ms,            //!< linearizations of master-side residuals w.r.t. slave-side dofs
        Epetra_SerialDenseMatrix&                                k_mm,            //!< linearizations of master-side residuals w.r.t. master-side dofs
        Epetra_SerialDenseVector&                                r_s,             //!< slave-side residual vector
        Epetra_SerialDenseVector&                                r_m              //!< master-side residual vector
        );
}; // class MortarCellCalcElch
} // namespace SCATRA
#endif // #ifndef SCATRA_TIMINT_MESHTYING_STRATEGY_S2I_ELCH_H
