/*----------------------------------------------------------------------------*/
/**
\file xcontact_multi_discretization_wrapper.H

\brief wrapper of the discretizations involved in a xcontact simulation

\maintainer Michael Hiermeier

\date Nov 10, 2016

\level 3

*/
/*----------------------------------------------------------------------------*/

#ifndef SRC_DRT_CONTACT_XCONTACT_XCONTACT_MULTI_DISCRETIZATION_WRAPPER_H_
#define SRC_DRT_CONTACT_XCONTACT_XCONTACT_MULTI_DISCRETIZATION_WRAPPER_H_


#include "../drt_structure_xstructure/xstr_multi_discretization_wrapper.H"
#include "../headers/pairedvector.H"

namespace XCONTACT
{
  class MultiDiscretizationWrapper : public ::XSTR::MultiDiscretizationWrapper
  {
   public:
    typedef GEN::pairedvector<XFEM::FieldName, DRT::Discretization*> XDisPairedPtrVector;
    typedef GEN::pairedvector<XFEM::FieldName, Teuchos::RCP<DRT::Discretization>>
        XDisPairedRCPVector;
    typedef GEN::pairedvector<XFEM::FieldName, Teuchos::RCP<Epetra_Map>> XMapPairedRCPVector;

   public:
    /// constructor
    MultiDiscretizationWrapper();

    /// destructor
    virtual ~MultiDiscretizationWrapper(){};

    /// add the mortar contact interface discretizations to the wrapper
    void AddContactIDiscret(
        const GEN::pairedvector<XFEM::FieldName, DRT::Discretization*>& idiscrets);

    /// add the scatra interface discretizations to the wrapper
    void AddScaTraIDiscret(
        const XFEM::FieldName& field, Teuchos::RCP<DRT::Discretization> idiscret);

    /// derived
    int FillComplete(bool assigndegreesoffreedom = true, bool initelements = true,
        bool doboundaryconditions = true, bool buildsystemmaps = true,
        bool setupmapextractor = true);

    void SetupContactScaTraDofExporter();

    const XDisPairedPtrVector& GetContactIDiscret() const;

    inline const bool& IsExporter() const { return is_exporter_; };

    /// @name Transfer vectors with contact map layout to the scatra map layout
    /// @{
    Teuchos::RCP<Epetra_Vector> Contact2ScaTra(const Teuchos::RCP<Epetra_Vector>& contact_vec,
        const enum XFEM::FieldName& field, const bool& do_normalize = false) const;
    Teuchos::RCP<Epetra_Vector> Contact2ScaTra(const Epetra_Vector& contact_vec,
        const enum XFEM::FieldName& field, const bool& do_normalize = false) const;

    Teuchos::RCP<Epetra_MultiVector> Contact2ScaTra(
        const Teuchos::RCP<Epetra_MultiVector>& contact_vec, const enum XFEM::FieldName& field,
        const bool& do_normalize = false) const;
    Teuchos::RCP<Epetra_MultiVector> Contact2ScaTra(const Epetra_MultiVector& contact_vec,
        const enum XFEM::FieldName& field, const bool& do_normalize = false) const;

    void Contact2ScaTra(const Epetra_MultiVector& contact_vec, Epetra_MultiVector& scatra_vec,
        const enum XFEM::FieldName& field, const bool& do_normalize = false) const;
    /// @}

    /// @name Transfer vectors with scaTra map layout to the contact map layout
    /// @{
    Teuchos::RCP<Epetra_Vector> ScaTra2Contact(
        const Teuchos::RCP<Epetra_Vector>& scatra_vec, const enum XFEM::FieldName& field) const;

    Teuchos::RCP<Epetra_Vector> ScaTra2Contact(
        const Epetra_Vector& scatra_vec, const enum XFEM::FieldName& field) const;

    Teuchos::RCP<Epetra_MultiVector> ScaTra2Contact(
        const Teuchos::RCP<Epetra_MultiVector>& scatra_vec,
        const enum XFEM::FieldName& field) const;

    Teuchos::RCP<Epetra_MultiVector> ScaTra2Contact(
        const Epetra_MultiVector& scatra_vec, const enum XFEM::FieldName& field) const;

    void ScaTra2Contact(const Epetra_MultiVector& scatra_vec, Epetra_MultiVector& contact_vec,
        const enum XFEM::FieldName& field) const;
    /// @}

    /** \brief Access one of the wrapped ScaTra discretizations (as reference)
     *
     *  \param field (in) : field enumerator corresponding to the desired discretization
     *
     *  \author hiermeier \date 10/16 */
    inline const DRT::Discretization& ScaTraDiscret(const XFEM::FieldName& field) const
    {
      CheckInit();
      return *(GetScaTraDiscretPtr(field));
    }
    inline DRT::Discretization& ScaTraDiscret(const XFEM::FieldName& field)
    {
      CheckInit();
      return *(GetScaTraDiscretPtr(field));
    }

    /** \brief Access one of the wrapped ScaTra discretizations (as rcp)
     *
     *  \param field (in) : field enumerator corresponding to the desired discretization
     *
     *  \author hiermeier \date 10/16 */
    inline Teuchos::RCP<const DRT::Discretization> ScaTraDiscretPtr(
        const XFEM::FieldName& field) const
    {
      CheckInit();
      return GetScaTraDiscretPtr(field).getConst();
    }
    inline const Teuchos::RCP<DRT::Discretization>& ScaTraDiscretPtr(const XFEM::FieldName& field)
    {
      CheckInit();
      return GetScaTraDiscretPtr(field);
    }

   private:
    void SetupContactScaTraDofExporter(const XFEM::FieldName& field,
        const DRT::Discretization& icontact, const DRT::Discretization& iscatra,
        const bool& resize);

    void BuildSlaveMasterMaps();

    const Teuchos::RCP<DRT::Discretization>& GetScaTraDiscretPtr(
        const enum XFEM::FieldName& field) const;

    void NormalizeInitialScaTraValues(Epetra_MultiVector& scatra_vec) const;

   private:
    XDisPairedPtrVector icontact_;

    XDisPairedRCPVector iscatra_;

    XMapPairedRCPVector slave_dof_row_maps_;
    XMapPairedRCPVector slave_dof_col_maps_;

    /// slave dofs in normal direction
    XMapPairedRCPVector slave_ndof_row_maps_;

    XMapPairedRCPVector master_dof_row_maps_;
    XMapPairedRCPVector master_dof_col_maps_;

    GEN::pairedvector<XFEM::FieldName, Teuchos::RCP<Epetra_Export>> contact2scatra_exporter_;
    GEN::pairedvector<XFEM::FieldName, Teuchos::RCP<Epetra_Export>> scatra2contact_exporter_;

    bool is_exporter_;
  };  // class MultiDiscretizationWrapper
}  // namespace XCONTACT


#endif /* SRC_DRT_CONTACT_XCONTACT_XCONTACT_MULTI_DISCRETIZATION_WRAPPER_H_ */
