/*-----------------------------------------------------------------------*/
/*! \file
\brief Contact interface capable of TSI

\maintainer Matthias Mayr

\level 3
*/
/*-----------------------------------------------------------------------*/
#ifndef CONTACT_TSI_INTERFACE_H_
#define CONTACT_TSI_INTERFACE_H_

#include "../drt_inpar/inpar_contact.H"
#include "contact_interface.H"


namespace CONTACT
{
  class CoTSIInterface : public CoInterface
  {
   public:
    /*!
    \brief Constructor

    */
    CoTSIInterface(const Teuchos::RCP<MORTAR::InterfaceDataContainer>& interfaceData_ptr,
        const int id, const Epetra_Comm& comm, const int dim,
        const Teuchos::ParameterList& icontact, bool selfcontact);

    /*!
    \brief Destructor

    */
    virtual ~CoTSIInterface() {}



    enum LinDM_X_mode
    {
      LinDM_Diss,
      LinDM_ThermoLM,
      linDM_ContactLMnormal
    };

    /*!
    \brief Assemble matrix LinStick containing linearizations

    This method builds an algebraic form of the FULL linearization
    of the tangential stick condition.

    */
    virtual void AssembleLinStick(LINALG::SparseMatrix& linstickLMglobal,
        LINALG::SparseMatrix& linstickDISglobal, LINALG::SparseMatrix& linstickTEMPglobal,
        Epetra_Vector& linstickRHSglobal);

    /*!
    \brief Assemble matrix LinSlip containing linearizations

    This method builds an algebraic form of the FULL linearization
    of the tangential slip condition. Concretely, this
    includes assembling the linearizations of the slave side
    nodal tangents and of the Mortar matrices D  and M.

    */
    virtual void AssembleLinSlip(LINALG::SparseMatrix& linslipLMglobal,
        LINALG::SparseMatrix& linslipDISglobal, LINALG::SparseMatrix& linslipTEMPglobal,
        Epetra_Vector& linslipRHSglobal);


    /*!
    \brief Assemble contact heat conduction

    This method assembles the Thermo-Lagrange-Multiplier line
    in the global system of equations, containing linearizations
    wrt displacements, contact forces, temperatures
    */
    virtual void AssembleLinConduct(LINALG::SparseMatrix& linConductDISglobal,
        LINALG::SparseMatrix& linConductTEMPglobal, LINALG::SparseMatrix& linConductThermoLMglobal,
        LINALG::SparseMatrix& linConductContactLMglobal);

    /*!
    \brief Assemble lumped mass matrix of the dual basis
           This is actually the D-Matrix. However, to keep
           the info where it comes from, we stick to that name

    */
    virtual void AssembleDualMassLumped(LINALG::SparseMatrix& dualMassGlobal,  /// dual mass matrix
        LINALG::SparseMatrix& linDualMassGlobal  // derivative of dual mass matrix wrt displacements
                                                 // multiplied with thermo-LM
    );

    /*!
      \brief Assemble the linearization of D and M times some nodal value

      This assembles  D_{jk,c) X_j and
                     -M_{jl,c} X_j (mind the minus sign!)
      where X_j is some (scalar) nodal value determined by "mode"
      */
    virtual void AssembleLinDM_X(LINALG::SparseMatrix* linD_X, LINALG::SparseMatrix* linM_X,
        const double fac, const LinDM_X_mode mode, const Teuchos::RCP<Epetra_Map> node_rowmap);

    /*!
      \brief Assemble D and M times the linearization of the nodal frictional dissipation

      This assembles  D_{jk) Diss_{j,c} and
                      M_{jl} Diss_{j,c}
                     where c is some discrete nodal DISPLACEMENT or LM dof.
      The dissipation is calculated by the tangential slip times the contact Lagrange multiplier
      Diss = \lambda^{contact}  \cdot  (1 - n \otimes n) \cdot jump
      */
    virtual void AssembleDM_linDiss(LINALG::SparseMatrix* d_LinDissDISP,
        LINALG::SparseMatrix* m_LinDissDISP, LINALG::SparseMatrix* d_LinDissContactLM,
        LINALG::SparseMatrix* m_LinDissContactLM, const double fac);

    /*!
      \brief Assemble the linearization of D and M times the temperature

      This assembles the derivative wrt displacements
      linDM_Temp_{jl}= LMj_n* ( D_{jk,c) T_k
                               -M_{jm,c} T_m )
                       + LMj_d * n_{d,l} * ( D_{jk} T_k
                                            -M_{jm} T_m )
      and the derivative wrt the contact Lagrange multiplier
      lin_lm_{jk} = nj_{k} * (D_jk T_k - M_jl T_l)
                                 (mind the minus sign!)
      */
    virtual void AssembleLinLMnDM_Temp(
        const double fac, LINALG::SparseMatrix* lin_disp, LINALG::SparseMatrix* lin_lm);

    /*!
      \brief Assemble the D and M times the normal Lagrange multiplier

      This assembles lambda_n* ( D_{jk) LMn_j
                                -M_{jl} LMn_j ) (mind the minus sign! No sum over j!)
      */
    virtual void AssembleDM_LMn(const double fac, LINALG::SparseMatrix* DM_LMn);

    /*!
      \brief Assemble inactive part of the thermal heat conduction equation

      This assembles lambda(thermo) = 0 for all inactive contact nodes
      */
    virtual void AssembleInactive(LINALG::SparseMatrix* linConductThermoLM);


    /*!
    \brief Initialize / reset interface for contact

    Derived version with some additional TSI related stuff!

    */
    void Initialize() final;

   protected:
  };  // class


}  // namespace CONTACT

#endif /* CONTACT_TSI_INTERFACE_H_ */
