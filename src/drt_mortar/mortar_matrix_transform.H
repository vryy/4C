/*----------------------------------------------------------------------------*/
/*!
\file mortar_matrix_transform.H

\brief matrix transformation tools: Switch between different parallel
distributions

\level 3

\maintainer Michael Hiermeier
\date Apr 4, 2017

*/
/*----------------------------------------------------------------------------*/

#ifndef SRC_DRT_MORTAR_MORTAR_MATRIX_TRANSFORM_H_
#define SRC_DRT_MORTAR_MORTAR_MATRIX_TRANSFORM_H_

#include "../drt_lib/drt_dserror.H"
#include "../headers/pairedvector.H"
#include "../drt_lib/drt_utils_matrix_vector_enums.H"

#include <Teuchos_RCP.hpp>

class Epetra_Map;
class Epetra_Export;
namespace LINALG
{
  class SparseMatrix;
}  // namespace LINALG
namespace MORTAR
{
  class MatrixRowColTransformer
  {
    typedef GEN::pairedvector<DRT::UTILS::MatBlockType, Teuchos::RCP<Epetra_Export>>
        plain_block_export_pairs;

   public:
    typedef GEN::pairedvector<DRT::UTILS::MatBlockType, Teuchos::RCP<Epetra_Map>*>
        plain_block_map_pairs;

   public:
    /// constructor
    MatrixRowColTransformer(const unsigned num_transformer);

    /// destructor
    ~MatrixRowColTransformer(){};

    void Init(const plain_block_map_pairs& redistributed_row,
        const plain_block_map_pairs& redistributed_column,
        const plain_block_map_pairs& unredistributed_row,
        const plain_block_map_pairs& unredistributed_column);

    void Setup();

    Teuchos::RCP<LINALG::SparseMatrix> RedistributedToUnredistributed(
        const DRT::UTILS::MatBlockType bt, const LINALG::SparseMatrix& src_mat);

    void RedistributedToUnredistributed(const DRT::UTILS::MatBlockType bt,
        const LINALG::SparseMatrix& src_mat, LINALG::SparseMatrix& dst_mat);

    Teuchos::RCP<LINALG::SparseMatrix> UnredistributedToRedistributed(
        const DRT::UTILS::MatBlockType bt, const LINALG::SparseMatrix& src_mat);

    void UnredistributedToRedistributed(const DRT::UTILS::MatBlockType bt,
        const LINALG::SparseMatrix& src_mat, LINALG::SparseMatrix& dst_mat);

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



#endif /* SRC_DRT_MORTAR_MORTAR_MATRIX_TRANSFORM_H_ */
