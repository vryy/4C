/*----------------------------------------------------------------------*/
/*! \file

\brief a common base class for all solid elements

\level 2

\maintainer Christoph Meier

 *----------------------------------------------------------------------*/


#include "so_base.H"
#include "../drt_mat/so3_material.H"
#include "../drt_structure_new/str_elements_paramsinterface.H"


/*----------------------------------------------------------------------*
 |  ctor (public)                                            vuong 03/15|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::So_base::So_base(int id, int owner)
    : DRT::Element(id, owner),
      kintype_(INPAR::STR::kinem_vague),
      interface_ptr_(Teuchos::null),
      material_post_setup_(false)
{
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       vuong 03/15|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::So_base::So_base(const DRT::ELEMENTS::So_base& old)
    : DRT::Element(old),
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
void DRT::ELEMENTS::So_base::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);
  // add base class Element
  Element::Pack(data);
  // kintype_
  AddtoPack(data, kintype_);

  // material post setup routine
  AddtoPack(data, static_cast<const int>(material_post_setup_));
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                           vuong 03/15|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_base::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position, data, type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");
  // extract base class Element
  std::vector<char> basedata(0);
  ExtractfromPack(position, data, basedata);
  Element::Unpack(basedata);
  // kintype_
  kintype_ = static_cast<INPAR::STR::KinemType>(ExtractInt(position, data));

  // material post setup routine
  material_post_setup_ = (ExtractInt(position, data) != 0);
}

/*----------------------------------------------------------------------*
 |  return solid material                                      (public) |
 |                                                           seitz 03/15|
 *----------------------------------------------------------------------*/
Teuchos::RCP<MAT::So3Material> DRT::ELEMENTS::So_base::SolidMaterial(int nummat) const
{
  return Teuchos::rcp_dynamic_cast<MAT::So3Material>(DRT::Element::Material(nummat), true);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_base::SetParamsInterfacePtr(const Teuchos::ParameterList& p)
{
  if (p.isParameter("interface"))
    interface_ptr_ = p.get<Teuchos::RCP<DRT::ELEMENTS::ParamsInterface>>("interface");
  else
    interface_ptr_ = Teuchos::null;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::ELEMENTS::ParamsInterface> DRT::ELEMENTS::So_base::ParamsInterfacePtr()
{
  return interface_ptr_;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
STR::ELEMENTS::ParamsInterface& DRT::ELEMENTS::So_base::StrParamsInterface()
{
  if (not IsParamsInterface()) dserror("The interface ptr is not set!");
  return *(Teuchos::rcp_dynamic_cast<STR::ELEMENTS::ParamsInterface>(interface_ptr_, true));
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void DRT::ELEMENTS::So_base::ErrorHandling(const double& det_curr, Teuchos::ParameterList& params,
    const int line_id, const STR::ELEMENTS::EvalErrorFlag flag)
{
  // check, if errors are tolerated or should throw a dserror
  if (IsParamsInterface() and StrParamsInterface().IsTolerateErrors())
  {
    StrParamsInterface().SetEleEvalErrorFlag(flag);
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
    // if the errors are not tolerated, throw a dserror
    else
    {
      if (det_curr == 0.0)
        dserror("ZERO DETERMINANT DETECTED in line %d", line_id);
      else if (det_curr < 0.0)
        dserror("NEGATIVE DETERMINANT DETECTED in line %d (value = %.5e)", line_id, det_curr);
    }
  }
}

// Check, whether the material post setup routine was
void DRT::ELEMENTS::So_base::EnsureMaterialPostSetup(Teuchos::ParameterList& params)
{
  if (!material_post_setup_)
  {
    MaterialPostSetup(params);
  }
}

void DRT::ELEMENTS::So_base::MaterialPostSetup(Teuchos::ParameterList& params)
{
  // This is the minimal implementation. Advanced materials may need extra implementation here.
  SolidMaterial()->PostSetup(params, Id());
  material_post_setup_ = true;
}
