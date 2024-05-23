/*----------------------------------------------------------------------------*/
/*! \file
\brief matrix transformation tools: Switch between different parallel
distributions

\level 3

*/
/*----------------------------------------------------------------------------*/


#include "4C_mortar_matrix_transform.hpp"

#include "4C_linalg_sparsematrix.hpp"

#include <Epetra_Distributor.h>
#include <Epetra_Export.h>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
MORTAR::MatrixRowColTransformer::MatrixRowColTransformer(const unsigned num_transformer)
    : isinit_(false),
      issetup_(false),
      slave_to_master_(num_transformer),
      master_to_slave_(num_transformer),
      slave_row_(num_transformer),
      slave_col_(num_transformer),
      master_row_(num_transformer),
      master_col_(num_transformer)
{
  /* empty */
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void MORTAR::MatrixRowColTransformer::Init(const plain_block_map_pairs& redistributed_row,
    const plain_block_map_pairs& redistributed_column,
    const plain_block_map_pairs& unredistributed_row,
    const plain_block_map_pairs& unredistributed_column)
{
  issetup_ = false;

  SetSlaveMapPairs(redistributed_row, redistributed_column);
  SetMasterMapPairs(unredistributed_row, unredistributed_column);

  isinit_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void MORTAR::MatrixRowColTransformer::SetSlaveMapPairs(
    const plain_block_map_pairs& redistributed_row,
    const plain_block_map_pairs& redistributed_column)
{
  slave_row_.clear();
  slave_row_ = redistributed_row;

  slave_col_.clear();
  slave_col_ = redistributed_column;

  issetup_ = false;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void MORTAR::MatrixRowColTransformer::SetMasterMapPairs(
    const plain_block_map_pairs& unredistributed_row,
    const plain_block_map_pairs& unredistributed_column)
{
  master_row_.clear();
  master_row_ = unredistributed_row;

  master_col_.clear();
  master_col_ = unredistributed_column;

  issetup_ = false;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void MORTAR::MatrixRowColTransformer::Setup()
{
  ThrowIfNotInit();


  for (plain_block_map_pairs::const_iterator cit = slave_row_.begin(); cit != slave_row_.end();
       ++cit)
  {
    const CONTACT::MatBlockType bt = cit->first;

    Teuchos::RCP<Epetra_Export>& slave_to_master = slave_to_master_[bt];
    slave_to_master = Teuchos::null;
    slave_to_master = Teuchos::rcp(new Epetra_Export(**master_row_[bt], **slave_row_[bt]));

    Teuchos::RCP<Epetra_Export>& master_to_slave = master_to_slave_[bt];
    master_to_slave = Teuchos::null;
    master_to_slave = Teuchos::rcp(new Epetra_Export(**slave_row_[bt], **master_row_[bt]));
  }

  issetup_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<CORE::LINALG::SparseMatrix>
MORTAR::MatrixRowColTransformer::redistributed_to_unredistributed(
    const CONTACT::MatBlockType bt, const CORE::LINALG::SparseMatrix& src_mat)
{
  throw_if_not_init_and_setup();

  Teuchos::RCP<CORE::LINALG::SparseMatrix> dst_mat = Teuchos::rcp(new CORE::LINALG::SparseMatrix(
      **master_row_[bt], src_mat.EpetraMatrix()->MaxNumEntries(), false, true));

  redistributed_to_unredistributed(bt, src_mat, *dst_mat);

  dst_mat->Complete(**master_col_[bt], **master_row_[bt]);
  return dst_mat;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void MORTAR::MatrixRowColTransformer::redistributed_to_unredistributed(
    const CONTACT::MatBlockType bt, const CORE::LINALG::SparseMatrix& src_mat,
    CORE::LINALG::SparseMatrix& dst_mat)
{
  throw_if_not_init_and_setup();

  const int err =
      dst_mat.EpetraMatrix()->Import(*src_mat.EpetraMatrix(), *slave_to_master_[bt], Insert);

  // reset the distributor of the exporter after use
  ResetExporter(slave_to_master_[bt]);

  if (err) FOUR_C_THROW("Import failed with err=%d", err);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<CORE::LINALG::SparseMatrix>
MORTAR::MatrixRowColTransformer::unredistributed_to_redistributed(
    const CONTACT::MatBlockType bt, const CORE::LINALG::SparseMatrix& src_mat)
{
  throw_if_not_init_and_setup();

  Teuchos::RCP<CORE::LINALG::SparseMatrix> dst_mat = Teuchos::rcp(new CORE::LINALG::SparseMatrix(
      **slave_row_[bt], src_mat.EpetraMatrix()->MaxNumEntries(), false, true));

  redistributed_to_unredistributed(bt, src_mat, *dst_mat);

  dst_mat->Complete(**slave_col_[bt], **slave_row_[bt]);
  return dst_mat;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void MORTAR::MatrixRowColTransformer::unredistributed_to_redistributed(
    const CONTACT::MatBlockType bt, const CORE::LINALG::SparseMatrix& src_mat,
    CORE::LINALG::SparseMatrix& dst_mat)
{
  throw_if_not_init_and_setup();

  const int err =
      dst_mat.EpetraMatrix()->Import(*src_mat.EpetraMatrix(), *master_to_slave_[bt], Insert);

  // reset the distributor of the exporter after use
  ResetExporter(master_to_slave_[bt]);

  if (err) FOUR_C_THROW("Import failed with err=%d", err);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void MORTAR::MatrixRowColTransformer::ResetExporter(Teuchos::RCP<Epetra_Export>& exporter) const
{
  exporter = Teuchos::rcp(new Epetra_Export(*exporter));
}

FOUR_C_NAMESPACE_CLOSE
