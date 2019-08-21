/*----------------------------------------------------------------------------*/
/*! \file
\brief Input for wall1 element.

\level 1

\maintainer Christoph Meier

*/
/*---------------------------------------------------------------------------*/

#include "wall1.H"
#include "../drt_lib/drt_linedefinition.H"
#include "../drt_mat/elasthyper.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool DRT::ELEMENTS::Wall1::ReadElement(
    const std::string& eletype, const std::string& distype, DRT::INPUT::LineDefinition* linedef)
{
  // set discretization type
  SetDisType(DRT::StringToDistype(distype));

  linedef->ExtractDouble("THICK", thickness_);
  if (thickness_ <= 0) dserror("WALL element thickness needs to be < 0");

  std::vector<int> ngp;
  linedef->ExtractIntVector("GP", ngp);

  if ((NumNode() == 4) and ((ngp[0] < 2) or (ngp[1] < 2)))
    dserror("Insufficient number of Gauss points");
  else if ((NumNode() == 8) and ((ngp[0] < 3) or (ngp[1] < 3)))
    dserror("Insufficient number of Gauss points");
  else if ((NumNode() == 9) and ((ngp[0] < 3) or (ngp[1] < 3)))
    dserror("Insufficient number of Gauss points");
  else if ((NumNode() == 6) and (ngp[0] < 3))
    dserror("Insufficient number of Gauss points");

  gaussrule_ = getGaussrule(&ngp[0]);

  // read number of material model
  int material = 0;
  linedef->ExtractInt("MAT", material);
  SetMaterial(material);

  Teuchos::RCP<MAT::Material> mat = Material();

  {
    const DRT::UTILS::IntegrationPoints2D intpoints(gaussrule_);
    const int numgp = intpoints.nquad;
    SolidMaterial()->Setup(numgp, linedef);
  }

  std::string buffer;
  // reduced dimension assumption
  linedef->ExtractString("STRESS_STRAIN", buffer);
  if (buffer == "plane_stress")
    wtype_ = plane_stress;
  else if (buffer == "plane_strain")
    wtype_ = plane_strain;
  else
    dserror("Illegal strain/stress type '%s'", buffer.c_str());

  // kinematics type
  linedef->ExtractString("KINEM", buffer);
  // geometrically linear
  if (buffer == "linear") kintype_ = INPAR::STR::kinem_linear;
  // geometrically non-linear with Total Lagrangean approach
  else if (buffer == "nonlinear")
    kintype_ = INPAR::STR::kinem_nonlinearTotLag;
  else
    dserror("Illegal KINEM type '%s'", buffer.c_str());

  // EAS type
  linedef->ExtractString("EAS", buffer);
  if (buffer == "none")
  {
    iseas_ = false;
  }
  else if (buffer == "full")
  {
    iseas_ = true;

    if (NumNode() == 9)
      dserror("eas-technology not necessary with 9 nodes");
    else if (NumNode() == 8)
      dserror("eas-technology not necessary with 8 nodes");
    else if (NumNode() == 3)
      dserror("eas-technology not implemented for tri3 elements");
    else if (NumNode() == 6)
      dserror("eas-technology not implemented for tri6 elements");
    else
    {
      // EAS enhanced deformation gradient parameters
      Epetra_SerialDenseMatrix alpha(
          Wall1::neas_, 1);  // if you change '4' here, then do it for alphao as well
      Epetra_SerialDenseMatrix alphao(Wall1::neas_, 1);

      // EAS portion of internal forces, also called enhacement vector s or Rtilde
      Epetra_SerialDenseVector feas(Wall1::neas_);
      // EAS matrix K_{alpha alpha}, also called Dtilde
      Epetra_SerialDenseMatrix invKaa(Wall1::neas_, Wall1::neas_);
      // EAS matrix K_{d alpha}
      Epetra_SerialDenseMatrix Kda(2 * NumNode(), Wall1::neas_);
      // EAS matrix K_{alpha d} // ONLY NEEDED FOR GENERALISED ENERGY-MOMENTUM METHOD
      Epetra_SerialDenseMatrix Kad(Wall1::neas_, 2 * NumNode());
      // EAS increment over last Newton step
      Epetra_SerialDenseMatrix eas_inc(Wall1::neas_, 1);

      // save EAS data into element container easdata_
      data_.Add("alpha", alpha);
      data_.Add("alphao", alphao);
      data_.Add("feas", feas);
      data_.Add("invKaa", invKaa);
      data_.Add("Kda", Kda);
      data_.Add("Kad", Kad);  // ONLY NEEDED FOR GENERALISED ENERGY-MOMENTUM METHOD
      data_.Add("eas_inc", eas_inc);
    }
  }
  else
  {
    dserror("Illegal EAS model");
  }

  // EAS type
  if (iseas_)
  {
    eastype_ = eas_q1e4;
  }
  else
  {
    eastype_ = eas_vague;
  }

  stresstype_ = w1_xy;

  // check for invalid combinations
  if (kintype_ == INPAR::STR::kinem_linear && iseas_ == true)
    dserror("ERROR: No EAS for geometrically linear WALL element");

  // validate kinematics of solid material
  SolidMaterial()->ValidKinematics(kintype_);

  // Validate that materials doesn't use extended update call.
  if (SolidMaterial()->UsesExtendedUpdate())
    dserror("This element currently does not support the extended update call.");

  return true;
}

/*----------------------------------------------------------------------*
 |  Get gaussrule on dependance of gausspoints                     mgit |
 *----------------------------------------------------------------------*/
DRT::UTILS::GaussRule2D DRT::ELEMENTS::Wall1::getGaussrule(int* ngp)
{
  DRT::UTILS::GaussRule2D rule = DRT::UTILS::intrule2D_undefined;

  switch (Shape())
  {
    case DRT::Element::quad4:
    case DRT::Element::quad8:
    case DRT::Element::quad9:
    {
      if ((ngp[0] == 2) && (ngp[1] == 2))
      {
        rule = DRT::UTILS::intrule_quad_4point;
      }
      else if ((ngp[0] == 3) && (ngp[1] == 3))
      {
        rule = DRT::UTILS::intrule_quad_9point;
      }
      else
        dserror("Unknown number of Gauss points for quad element");
      break;
    }
    case DRT::Element::nurbs4:
    case DRT::Element::nurbs9:
    {
      if ((ngp[0] == 2) && (ngp[1] == 2))
      {
        rule = DRT::UTILS::intrule_quad_4point;
      }
      else if ((ngp[0] == 3) && (ngp[1] == 3))
      {
        rule = DRT::UTILS::intrule_quad_9point;
      }
      else if ((ngp[0] == 4) && (ngp[1] == 4))
      {
        rule = DRT::UTILS::intrule_quad_16point;
      }
      else if ((ngp[0] == 5) && (ngp[1] == 5))
      {
        rule = DRT::UTILS::intrule_quad_25point;
      }
      else if ((ngp[0] == 10) && (ngp[1] == 10))
      {
        rule = DRT::UTILS::intrule_quad_100point;
      }
      else
        dserror("Unknown number of Gauss points for nurbs element");
      break;
    }
    case DRT::Element::tri3:
    case DRT::Element::tri6:
    {
      if ((ngp[0] == 1) && (ngp[1] == 0))
      {
        rule = DRT::UTILS::intrule_tri_1point;
      }
      else if ((ngp[0] == 3) && (ngp[1] == 0))
      {
        rule = DRT::UTILS::intrule_tri_3point;
      }
      else if ((ngp[0] == 6) && (ngp[1] == 0))
      {
        rule = DRT::UTILS::intrule_tri_6point;
      }
      else
        dserror("Unknown number of Gauss points for tri element");
      break;
    }
    default:
      dserror("Unknown distype");
      break;
  }
  return rule;
}
