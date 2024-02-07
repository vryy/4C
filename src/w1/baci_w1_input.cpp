/*----------------------------------------------------------------------------*/
/*! \file
\brief Input for wall1 element.

\level 1


*/
/*---------------------------------------------------------------------------*/

#include "baci_io_linedefinition.H"
#include "baci_mat_elasthyper.H"
#include "baci_w1.H"

BACI_NAMESPACE_OPEN


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool DRT::ELEMENTS::Wall1::ReadElement(
    const std::string& eletype, const std::string& distype, INPUT::LineDefinition* linedef)
{
  // set discretization type
  SetDisType(CORE::FE::StringToCellType(distype));

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

  gaussrule_ = getGaussrule(ngp.data());

  // read number of material model
  int material = 0;
  linedef->ExtractInt("MAT", material);
  SetMaterial(material);

  Teuchos::RCP<MAT::Material> mat = Material();

  {
    const CORE::FE::IntegrationPoints2D intpoints(gaussrule_);
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
  if (buffer == "linear") kintype_ = INPAR::STR::KinemType::linear;
  // geometrically non-linear with Total Lagrangean approach
  else if (buffer == "nonlinear")
    kintype_ = INPAR::STR::KinemType::nonlinearTotLag;
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
      CORE::LINALG::SerialDenseMatrix alpha(
          Wall1::neas_, 1);  // if you change '4' here, then do it for alphao as well
      CORE::LINALG::SerialDenseMatrix alphao(Wall1::neas_, 1);

      // EAS portion of internal forces, also called enhacement vector s or Rtilde
      CORE::LINALG::SerialDenseVector feas(Wall1::neas_);
      // EAS matrix K_{alpha alpha}, also called Dtilde
      CORE::LINALG::SerialDenseMatrix invKaa(Wall1::neas_, Wall1::neas_);
      // EAS matrix K_{d alpha}
      CORE::LINALG::SerialDenseMatrix Kda(2 * NumNode(), Wall1::neas_);
      // EAS matrix K_{alpha d} // ONLY NEEDED FOR GENERALISED ENERGY-MOMENTUM METHOD
      CORE::LINALG::SerialDenseMatrix Kad(Wall1::neas_, 2 * NumNode());
      // EAS increment over last Newton step
      CORE::LINALG::SerialDenseMatrix eas_inc(Wall1::neas_, 1);

      // save EAS data into element container easdata_
      easdata_.alpha = alpha;
      easdata_.alphao = alphao;
      easdata_.feas = feas;
      easdata_.invKaa = invKaa;
      easdata_.Kda = Kda;
      easdata_.Kad = Kad;  // ONLY NEEDED FOR GENERALISED ENERGY-MOMENTUM METHOD
      easdata_.eas_inc = eas_inc;
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
  if (kintype_ == INPAR::STR::KinemType::linear && iseas_ == true)
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
CORE::FE::GaussRule2D DRT::ELEMENTS::Wall1::getGaussrule(int* ngp)
{
  CORE::FE::GaussRule2D rule = CORE::FE::GaussRule2D::undefined;

  switch (Shape())
  {
    case CORE::FE::CellType::quad4:
    case CORE::FE::CellType::quad8:
    case CORE::FE::CellType::quad9:
    {
      if ((ngp[0] == 2) && (ngp[1] == 2))
      {
        rule = CORE::FE::GaussRule2D::quad_4point;
      }
      else if ((ngp[0] == 3) && (ngp[1] == 3))
      {
        rule = CORE::FE::GaussRule2D::quad_9point;
      }
      else
        dserror("Unknown number of Gauss points for quad element");
      break;
    }
    case CORE::FE::CellType::nurbs4:
    case CORE::FE::CellType::nurbs9:
    {
      if ((ngp[0] == 2) && (ngp[1] == 2))
      {
        rule = CORE::FE::GaussRule2D::quad_4point;
      }
      else if ((ngp[0] == 3) && (ngp[1] == 3))
      {
        rule = CORE::FE::GaussRule2D::quad_9point;
      }
      else if ((ngp[0] == 4) && (ngp[1] == 4))
      {
        rule = CORE::FE::GaussRule2D::quad_16point;
      }
      else if ((ngp[0] == 5) && (ngp[1] == 5))
      {
        rule = CORE::FE::GaussRule2D::quad_25point;
      }
      else if ((ngp[0] == 10) && (ngp[1] == 10))
      {
        rule = CORE::FE::GaussRule2D::quad_100point;
      }
      else
        dserror("Unknown number of Gauss points for nurbs element");
      break;
    }
    case CORE::FE::CellType::tri3:
    case CORE::FE::CellType::tri6:
    {
      if ((ngp[0] == 1) && (ngp[1] == 0))
      {
        rule = CORE::FE::GaussRule2D::tri_1point;
      }
      else if ((ngp[0] == 3) && (ngp[1] == 0))
      {
        rule = CORE::FE::GaussRule2D::tri_3point;
      }
      else if ((ngp[0] == 6) && (ngp[1] == 0))
      {
        rule = CORE::FE::GaussRule2D::tri_6point;
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

BACI_NAMESPACE_CLOSE
