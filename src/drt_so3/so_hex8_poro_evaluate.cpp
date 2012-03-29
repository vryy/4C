/*!----------------------------------------------------------------------
\file so_hex8_poro_evaluate.cpp
\brief

<pre>
   Maintainer: Anh-Tu Vuong
               vuong@lnm.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-15264
</pre>

*----------------------------------------------------------------------*/
#include "so_hex8_poro.H"

/*----------------------------------------------------------------------*
 |  evaluate the element (public)                                       |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::So_hex8_poro::Evaluate(ParameterList& params,
                                    DRT::Discretization&      discretization,
                                    DRT::Element::LocationArray& la,
                                    Epetra_SerialDenseMatrix& elemat1_epetra,
                                    Epetra_SerialDenseMatrix& elemat2_epetra,
                                    Epetra_SerialDenseVector& elevec1_epetra,
                                    Epetra_SerialDenseVector& elevec2_epetra,
                                    Epetra_SerialDenseVector& elevec3_epetra)
{
  // start with "none"
  So_hex8::ActionType act = So_hex8::none;

  // get the required action
  string action = params.get<string>("action","none");
  if (action == "none") dserror("No action supplied");
  else if (action=="calc_struct_update_istep")          act = So_hex8::calc_struct_update_istep;
  else if (action=="calc_struct_internalforce")         act = So_hex8::calc_struct_internalforce;
  else if (action=="calc_struct_nlnstiff")              act = So_hex8::calc_struct_nlnstiff;
  else if (action=="calc_struct_nlnstiffmass")          act = So_hex8::calc_struct_nlnstiffmass;
  else if (action=="calc_struct_multidofsetcoupling")   act = So_hex8::calc_struct_multidofsetcoupling;
  else if (action=="postprocess_stress")                act = So_hex8::postprocess_stress;
  else dserror("Unknown type of action for So_hex8_poro: %s",action.c_str());
  // what should the element do
  switch(act)
  {
  //==================================================================================
  // nonlinear stiffness and internal force vector for poroelasticity
  case So_hex8::calc_struct_nlnstiff:
    // nonlinear stiffness, mass matrix and internal force vector for poroelasticity
  case So_hex8::calc_struct_nlnstiffmass:
    // nonlinear stiffness and internal force vector for poroelasticity
  case So_hex8::calc_struct_internalforce:
  {
    So_hex8::Evaluate(params,
                      discretization,
                      la,
                      elemat1_epetra,
                      elemat2_epetra,
                      elevec1_epetra,
                      elevec2_epetra,
                      elevec3_epetra);

    So3_Poro<DRT::Element::hex8>::Evaluate(params,
                      discretization,
                      la,
                      elemat1_epetra,
                      elemat2_epetra,
                      elevec1_epetra,
                      elevec2_epetra,
                      elevec3_epetra);
  }
  break;

  //==================================================================================
  // coupling terms in force-vector and stiffness matrix for poroelasticity
  case So_hex8::calc_struct_multidofsetcoupling:
  {
    So3_Poro<DRT::Element::hex8>::Evaluate(params,
                      discretization,
                      la,
                      elemat1_epetra,
                      elemat2_epetra,
                      elevec1_epetra,
                      elevec2_epetra,
                      elevec3_epetra);
  }
  break;

  //==================================================================================
  // postprocess stresses/strains at gauss points

  // note that in the following, quantities are always referred to as
  // "stresses" etc. although they might also apply to strains
  // (depending on what this routine is called for from the post filter)
  case So3_Poro<DRT::Element::hex8>::postprocess_stress:
  case So_hex8::calc_struct_update_istep:
  {
    So_hex8::Evaluate(params,
                      discretization,
                      la,
                      elemat1_epetra,
                      elemat2_epetra,
                      elevec1_epetra,
                      elevec2_epetra,
                      elevec3_epetra);
  }
  break;

  //==================================================================================
  default:
  dserror("Unknown type of action for So_hex8_poro");
  } // action
  return 0;
}


/*----------------------------------------------------------------------*
 |  init the element (public)                                           |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::So_hex8PoroType::Initialize(DRT::Discretization& dis)
{
  for (int i=0; i<dis.NumMyColElements(); ++i)
  {
    if (dis.lColElement(i)->ElementType() != *this) continue;
    DRT::ELEMENTS::So_hex8_poro* actele = dynamic_cast<DRT::ELEMENTS::So_hex8_poro*>(dis.lColElement(i));
    if (!actele) dserror("cast to So_hex8_poro* failed");
    actele->So_hex8::InitJacobianMapping();
    actele->So3_Poro<DRT::Element::hex8>::InitJacobianMapping();
  }
  return 0;
}

