/*----------------------------------------------------------------------*/
/*! \file

\brief prestress functionality in solid elements

\level 2

*----------------------------------------------------------------------*/

#include "4C_so3_prestress.hpp"

#include "4C_comm_pack_helpers.hpp"
#include "4C_utils_exceptions.hpp"

FOUR_C_NAMESPACE_OPEN


Discret::ELEMENTS::PreStressType Discret::ELEMENTS::PreStressType::instance_;


/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 07/08|
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::PreStress::PreStress(const int numnode, const int ngp, const bool istet4)
    : ParObject(), isinit_(false), numnode_(numnode)
{
  // allocate history memory
  fhist_ = Teuchos::rcp(new Core::LinAlg::SerialDenseMatrix(ngp, 9));
  if (!istet4)
    inv_jhist_ = Teuchos::rcp(new Core::LinAlg::SerialDenseMatrix(ngp, 9));
  else
    inv_jhist_ = Teuchos::rcp(new Core::LinAlg::SerialDenseMatrix(ngp, 12));

  // init the deformation gradient history
  Core::LinAlg::Matrix<3, 3> F(true);
  F(0, 0) = F(1, 1) = F(2, 2) = 1.0;
  for (int i = 0; i < num_gp(); ++i) matrixto_storage(i, F, f_history());
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       mwgee 11/06|
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::PreStress::PreStress(const Discret::ELEMENTS::PreStress& old)
    : ParObject(old),
      isinit_(old.isinit_),
      numnode_(old.numnode_),
      fhist_(Teuchos::rcp(new Core::LinAlg::SerialDenseMatrix(old.f_history()))),
      inv_jhist_(Teuchos::rcp(new Core::LinAlg::SerialDenseMatrix(old.j_history())))
{
  return;
}


/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                            gee 02/07 |
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::PreStress::pack(Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);

  // pack type of this instance of ParObject
  int type = unique_par_object_id();
  add_to_pack(data, type);

  // pack isinit_
  add_to_pack(data, isinit_);

  // pack numnode_
  add_to_pack(data, numnode_);

  // pack Fhist_
  add_to_pack(data, *fhist_);

  // pack invJhist_
  add_to_pack(data, *inv_jhist_);

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                            gee 02/07 |
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::PreStress::unpack(Core::Communication::UnpackBuffer& buffer)
{
  Core::Communication::extract_and_assert_id(buffer, unique_par_object_id());

  // extract isinit_
  isinit_ = extract_int(buffer);

  // extract numnode_
  extract_from_pack(buffer, numnode_);

  // extract Fhist_
  extract_from_pack(buffer, *fhist_);

  // extract invJhist_
  extract_from_pack(buffer, *inv_jhist_);

  FOUR_C_THROW_UNLESS(buffer.at_end(), "Buffer not fully consumed.");
  return;
}

int Discret::ELEMENTS::PreStress::unique_par_object_id() const
{
  return PreStressType::instance().unique_par_object_id();
}

FOUR_C_NAMESPACE_CLOSE
