/*----------------------------------------------------------------------*/
/*! \file

\brief a common base class for all solid elements

\level 2


 *----------------------------------------------------------------------*/


#include "4C_so3_base.hpp"

#include "4C_mat_so3_material.hpp"
#include "4C_structure_new_elements_paramsinterface.hpp"

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*
 |  ctor (public)                                            vuong 03/15|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::SoBase::SoBase(int id, int owner)
    : Core::Elements::Element(id, owner),
      kintype_(Inpar::STR::KinemType::vague),
      interface_ptr_(Teuchos::null),
      material_post_setup_(false)
{
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       vuong 03/15|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::SoBase::SoBase(const Discret::ELEMENTS::SoBase& old)
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
void Discret::ELEMENTS::SoBase::pack(Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  add_to_pack(data, type);
  // add base class Element
  Element::pack(data);
  // kintype_
  add_to_pack(data, kintype_);

  // material post setup routine
  add_to_pack(data, static_cast<int>(material_post_setup_));
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                           vuong 03/15|
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::SoBase::unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  Core::Communication::ExtractAndAssertId(position, data, UniqueParObjectId());

  // extract base class Element
  std::vector<char> basedata(0);
  extract_from_pack(position, data, basedata);
  Element::unpack(basedata);
  // kintype_
  kintype_ = static_cast<Inpar::STR::KinemType>(extract_int(position, data));

  // material post setup routine
  material_post_setup_ = (extract_int(position, data) != 0);
}

/*----------------------------------------------------------------------*
 |  return solid material                                      (public) |
 |                                                           seitz 03/15|
 *----------------------------------------------------------------------*/
Teuchos::RCP<Mat::So3Material> Discret::ELEMENTS::SoBase::SolidMaterial(int nummat) const
{
  return Teuchos::rcp_dynamic_cast<Mat::So3Material>(
      Core::Elements::Element::Material(nummat), true);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::SoBase::set_params_interface_ptr(const Teuchos::ParameterList& p)
{
  if (p.isParameter("interface"))
    interface_ptr_ = p.get<Teuchos::RCP<Core::Elements::ParamsInterface>>("interface");
  else
    interface_ptr_ = Teuchos::null;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<Core::Elements::ParamsInterface> Discret::ELEMENTS::SoBase::ParamsInterfacePtr()
{
  return interface_ptr_;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
STR::ELEMENTS::ParamsInterface& Discret::ELEMENTS::SoBase::str_params_interface()
{
  if (not IsParamsInterface()) FOUR_C_THROW("The interface ptr is not set!");
  return *(Teuchos::rcp_dynamic_cast<STR::ELEMENTS::ParamsInterface>(interface_ptr_, true));
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Discret::ELEMENTS::SoBase::error_handling(const double& det_curr,
    Teuchos::ParameterList& params, const int line_id, const STR::ELEMENTS::EvalErrorFlag flag)
{
  // check, if errors are tolerated or should throw a FOUR_C_THROW
  if (IsParamsInterface() and str_params_interface().is_tolerate_errors())
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
void Discret::ELEMENTS::SoBase::ensure_material_post_setup(Teuchos::ParameterList& params)
{
  if (!material_post_setup_)
  {
    material_post_setup(params);
  }
}

void Discret::ELEMENTS::SoBase::material_post_setup(Teuchos::ParameterList& params)
{
  // This is the minimal implementation. Advanced materials may need extra implementation here.
  SolidMaterial()->post_setup(params, Id());
  material_post_setup_ = true;
}

FOUR_C_NAMESPACE_CLOSE
