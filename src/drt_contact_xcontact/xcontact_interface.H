/*---------------------------------------------------------------------------*/
/*!
\file xcontact_interface.H

\brief Interface class for the eXtended contact evaluation

\level 3

\maintainer Michael Hiermeier

\date Aug 17, 2016
*/
/*---------------------------------------------------------------------------*/

#ifndef XCONTACT_INTERFACE_H
#define XCONTACT_INTERFACE_H

#include "../drt_contact/contact_interface.H"
#include "../drt_xfem/xfem_enum_lists.H"

namespace DRT
{
  class DiscretizationInterface;
}  // namespace DRT
namespace XCONTACT
{
  class Interface : public CONTACT::CoInterface
  {
   public:
    //! @name Constructors and destructors and related methods
    //! @{

    //! Constructor
    Interface(const Teuchos::RCP<MORTAR::IDataContainer>& idata_ptr, const int id,
        const Epetra_Comm& comm, const int dim, const Teuchos::ParameterList& icontact,
        bool selfcontact, INPAR::MORTAR::RedundantStorage redundant);

    //! Destructor
    virtual ~Interface() {}

    //! @}


    //! @name Accessors
    //! @{

    /*! Access local slave dof row map in normal direction. */
    const Teuchos::RCP<Epetra_Map> SlaveRowNDofs() const
    {
      if (Filled())
      {
        return sndofrowmap_;
      }
      else
      {
        dserror("XContact interface: FillComplete was not called.");
        exit(EXIT_FAILURE);
      }
    }

    /*! Access local slave dof row map in tangential direction. */
    const Teuchos::RCP<Epetra_Map> SlaveRowTDofs() const
    {
      if (Filled())
      {
        return stdofrowmap_;
      }
      else
      {
        dserror("XContact interface: FillComplete was not called.");
        exit(EXIT_FAILURE);
      }
    }

    //! get the field name enumerator of the parent discretization
    const enum XFEM::FieldName& ParentDiscretType() const { return parent_discret_.first; }

    //! access the parent discretization of this interface
    const DRT::DiscretizationInterface& ParentDiscret() const { return *parent_discret_.second; }

    //! @}


    //! @name Initialize and evaluate interface element contributions
    //! @{

    //! Initialize all XContact related quantities
    void Initialize();

    //! @}

    //! @name Assemble partial global vectors
    //! @{

    /*! \brief Assemble weighted gap vector (slave normal dof row map layout)
     *
     *  \author hiermeier
     *  \date 11/16 */
    void AssembleWeightedGap(Epetra_Vector& wgap) const;

    /*! \brief Assemble constraint residual vector and Lagrange multiplier vector
     *  in normal direction
     *
     *  \author hofer
     *  \date 08/16 */
    void AssembleContactRHS(Epetra_Vector& Wc_lm, Epetra_Vector& lm) const;

    //! @}

    //! @name Assemble partial global matrices
    //! @{

    /*! \brief Assemble mortar matrices
     *
     *  D (slave) and M (master)
     *
     *  \author hofer
     *  \date 08/16 */
    void AssembleMortar(LINALG::SparseMatrix& D, LINALG::SparseMatrix& M) const;

    /*! \brief Assemble structural contact tangent matrix
     *
     *  \author hofer
     *  \date 08/16 */
    void AssembleWcUU(LINALG::SparseMatrix& Wc_su_u, LINALG::SparseMatrix& Wc_mu_u) const;

    //! @}

    //! @name Basic setup
    //! @{

    /*! \brief Split sdofrowmap_ into normal and tangential part
     *
     *  \author hofer
     *  \date 08/16 */
    void SplitSlaveDofs();

   protected:
    /// derived
    void PreMortarCoupling(const MORTAR::MortarElement* sele,
        const std::vector<MORTAR::MortarElement*> mele,
        const Teuchos::RCP<MORTAR::ParamsInterface>& mparams_ptr) const;
    /// derived
    void PostMortarCoupling(const MORTAR::MortarElement* sele,
        const std::vector<MORTAR::MortarElement*> mele,
        const Teuchos::RCP<MORTAR::ParamsInterface>& mparams_ptr) const;

    //! [derived]
    virtual void UpdateMasterSlaveSets();

    //! @}

   protected:
    //! @name Member variables
    //! @{

    //! Pointer to map of slave DOFs in normal direction
    Teuchos::RCP<Epetra_Map> sndofrowmap_;

    //! Pointer to map of slave DOFs in tangential direction
    Teuchos::RCP<Epetra_Map> stdofrowmap_;

    //! Parent discretization pair
    std::pair<enum XFEM::FieldName, Teuchos::RCP<const DRT::DiscretizationInterface>>
        parent_discret_;

    //! @}

  };  // class Interface
}  // namespace XCONTACT

#endif  // #ifndef XCONTACT_INTERFACE_H
