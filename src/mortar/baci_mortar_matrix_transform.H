/*----------------------------------------------------------------------------*/
/*! \file
\brief matrix transformation tools: Switch between different parallel
distributions

\level 3

*/
/*----------------------------------------------------------------------------*/

#ifndef BACI_MORTAR_MATRIX_TRANSFORM_H
#define BACI_MORTAR_MATRIX_TRANSFORM_H

#include "baci_contact_utils.H"
#include "baci_utils_exceptions.H"
#include "baci_utils_pairedvector.H"

#include <Teuchos_RCP.hpp>

class Epetra_Map;
class Epetra_Export;
namespace CORE::LINALG
{
  class SparseMatrix;
}  // namespace CORE::LINALG
namespace MORTAR
{
  class MatrixRowColTransformer
  {
    typedef CORE::GEN::pairedvector<CONTACT::MatBlockType, Teuchos::RCP<Epetra_Export>>
        plain_block_export_pairs;

   public:
    typedef CORE::GEN::pairedvector<CONTACT::MatBlockType, Teuchos::RCP<Epetra_Map>*>
        plain_block_map_pairs;

   public:
    /// constructor
    MatrixRowColTransformer(const unsigned num_transformer);


    void Init(const plain_block_map_pairs& redistributed_row,
        const plain_block_map_pairs& redistributed_column,
        const plain_block_map_pairs& unredistributed_row,
        const plain_block_map_pairs& unredistributed_column);

    void Setup();

    Teuchos::RCP<CORE::LINALG::SparseMatrix> RedistributedToUnredistributed(
        const CONTACT::MatBlockType bt, const CORE::LINALG::SparseMatrix& src_mat);

    void RedistributedToUnredistributed(const CONTACT::MatBlockType bt,
        const CORE::LINALG::SparseMatrix& src_mat, CORE::LINALG::SparseMatrix& dst_mat);

    Teuchos::RCP<CORE::LINALG::SparseMatrix> UnredistributedToRedistributed(
        const CONTACT::MatBlockType bt, const CORE::LINALG::SparseMatrix& src_mat);

    void UnredistributedToRedistributed(const CONTACT::MatBlockType bt,
        const CORE::LINALG::SparseMatrix& src_mat, CORE::LINALG::SparseMatrix& dst_mat);

   private:
    void SetSlaveMapPairs(const plain_block_map_pairs& redistributed_row,
        const plain_block_map_pairs& redistributed_column);

    void SetMasterMapPairs(const plain_block_map_pairs& unredistributed_row,
        const plain_block_map_pairs& unredistributed_column);

    // hide empty constructor
    MatrixRowColTransformer();

    inline void ThrowIfNotInit() const
    {
      if (not isinit_) dserror("Call Init() first!");
    }

    inline void ThrowIfNotInitAndSetup() const
    {
      ThrowIfNotInit();
      if (not issetup_) dserror("Call Setup() first!");
    }

    void ResetExporter(Teuchos::RCP<Epetra_Export>& exporter) const;

   private:
    bool isinit_;
    bool issetup_;

    plain_block_export_pairs slave_to_master_;
    plain_block_export_pairs master_to_slave_;

    plain_block_map_pairs slave_row_;
    plain_block_map_pairs slave_col_;

    plain_block_map_pairs master_row_;
    plain_block_map_pairs master_col_;

  };  // class MatrixTransform

}  // namespace MORTAR



#endif  // MORTAR_MATRIX_TRANSFORM_H
