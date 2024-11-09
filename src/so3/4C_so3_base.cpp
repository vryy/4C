// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_so3_base.hpp"

#include "4C_comm_pack_helpers.hpp"
#include "4C_mat_so3_material.hpp"
#include "4C_structure_new_elements_paramsinterface.hpp"

#include <Teuchos_ParameterList.hpp>

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*
 |  ctor (public)                                            vuong 03/15|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
Discret::Elements::SoBase::SoBase(int id, int owner)
    : Core::Elements::Element(id, owner),
      kintype_(Inpar::Solid::KinemType::vague),
      interface_ptr_(nullptr),
      material_post_setup_(false)
{
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       vuong 03/15|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
Discret::Elements::SoBase::SoBase(const Discret::Elements::SoBase& old)
    : Core::Elements::Element(old),
      kintype_(old.kintype_),
      interface_ptr_(old.interface_ptr_),
      material_post_setup_(false)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                           vuong 03/15|
 *----------------------------------------------------------------------*/
void Discret::Elements::SoBase::pack(Core::Communication::PackBuffer& data) const
{
  // pack type of this instance of ParObject
  int type = unique_par_object_id();
  add_to_pack(data, type);
  // add base class Element
  Element::pack(data);
  // kintype_
  add_to_pack(data, kintype_);

  // material post setup routine
  add_to_pack(data, material_post_setup_);
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                           vuong 03/15|
 *----------------------------------------------------------------------*/
void Discret::Elements::SoBase::unpack(Core::Communication::UnpackBuffer& buffer)
{
  Core::Communication::extract_and_assert_id(buffer, unique_par_object_id());

  // extract base class Element
  Element::unpack(buffer);
  // kintype_
  extract_from_pack(buffer, kintype_);

  // material post setup routine
  extract_from_pack(buffer, material_post_setup_);
}

/*----------------------------------------------------------------------*
 |  return solid material                                      (public) |
 |                                                           seitz 03/15|
 *----------------------------------------------------------------------*/
std::shared_ptr<Mat::So3Material> Discret::Elements::SoBase::solid_material(int nummat) const
{
  return std::dynamic_pointer_cast<Mat::So3Material>(Core::Elements::Element::material(nummat));
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Discret::Elements::SoBase::set_params_interface_ptr(const Teuchos::ParameterList& p)
{
  if (p.isParameter("interface"))
    interface_ptr_ = p.get<std::shared_ptr<Core::Elements::ParamsInterface>>("interface");
  else
    interface_ptr_ = nullptr;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
std::shared_ptr<Core::Elements::ParamsInterface> Discret::Elements::SoBase::params_interface_ptr()
{
  return interface_ptr_;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Solid::Elements::ParamsInterface& Discret::Elements::SoBase::str_params_interface()
{
  if (not is_params_interface()) FOUR_C_THROW("The interface ptr is not set!");
  return *(std::dynamic_pointer_cast<Solid::Elements::ParamsInterface>(interface_ptr_));
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Discret::Elements::SoBase::error_handling(const double& det_curr,
    Teuchos::ParameterList& params, const int line_id, const Solid::Elements::EvalErrorFlag flag)
{
  // check, if errors are tolerated or should throw a FOUR_C_THROW
  if (is_params_interface() and str_params_interface().is_tolerate_errors())
  {
    str_params_interface().set_ele_eval_error_flag(flag);
    return;
  }
  else
  {
    // === DEPRECATED (hiermeier, 11/17) ==================================
    bool error_tol = false;
    if (params.isParameter("tolerate_errors")) error_tol = params.get<bool>("tolerate_errors");
    if (error_tol)
    {
      params.set<bool>("eval_error", true);
      return;
    }
    // if the errors are not tolerated, throw a FOUR_C_THROW
    else
    {
      if (det_curr == 0.0)
        FOUR_C_THROW("ZERO DETERMINANT DETECTED in line %d", line_id);
      else if (det_curr < 0.0)
        FOUR_C_THROW("NEGATIVE DETERMINANT DETECTED in line %d (value = %.5e)", line_id, det_curr);
    }
  }
}

// Check, whether the material post setup routine was
void Discret::Elements::SoBase::ensure_material_post_setup(Teuchos::ParameterList& params)
{
  if (!material_post_setup_)
  {
    material_post_setup(params);
  }
}

void Discret::Elements::SoBase::material_post_setup(Teuchos::ParameterList& params)
{
  // This is the minimal implementation. Advanced materials may need extra implementation here.
  solid_material()->post_setup(params, id());
  material_post_setup_ = true;
}

FOUR_C_NAMESPACE_CLOSE
