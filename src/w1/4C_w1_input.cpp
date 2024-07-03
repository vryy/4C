/*----------------------------------------------------------------------------*/
/*! \file
\brief Input for wall1 element.

\level 1


*/
/*---------------------------------------------------------------------------*/

#include "4C_io_linedefinition.hpp"
#include "4C_mat_elasthyper.hpp"
#include "4C_w1.hpp"

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool Discret::ELEMENTS::Wall1::ReadElement(
    const std::string& eletype, const std::string& distype, Input::LineDefinition* linedef)
{
  // set discretization type
  SetDisType(Core::FE::StringToCellType(distype));

  linedef->extract_double("THICK", thickness_);
  if (thickness_ <= 0) FOUR_C_THROW("WALL element thickness needs to be < 0");

  std::vector<int> ngp;
  linedef->extract_int_vector("GP", ngp);

  if ((num_node() == 4) and ((ngp[0] < 2) or (ngp[1] < 2)))
    FOUR_C_THROW("Insufficient number of Gauss points");
  else if ((num_node() == 8) and ((ngp[0] < 3) or (ngp[1] < 3)))
    FOUR_C_THROW("Insufficient number of Gauss points");
  else if ((num_node() == 9) and ((ngp[0] < 3) or (ngp[1] < 3)))
    FOUR_C_THROW("Insufficient number of Gauss points");
  else if ((num_node() == 6) and (ngp[0] < 3))
    FOUR_C_THROW("Insufficient number of Gauss points");

  gaussrule_ = get_gaussrule(ngp.data());

  // read number of material model
  int material = 0;
  linedef->extract_int("MAT", material);
  SetMaterial(0, Mat::Factory(material));

  Teuchos::RCP<Core::Mat::Material> mat = Material();

  {
    const Core::FE::IntegrationPoints2D intpoints(gaussrule_);
    const int numgp = intpoints.nquad;
    SolidMaterial()->setup(numgp, linedef);
  }

  std::string buffer;
  // reduced dimension assumption
  linedef->extract_string("STRESS_STRAIN", buffer);
  if (buffer == "plane_stress")
    wtype_ = plane_stress;
  else if (buffer == "plane_strain")
    wtype_ = plane_strain;
  else
    FOUR_C_THROW("Illegal strain/stress type '%s'", buffer.c_str());

  // kinematics type
  linedef->extract_string("KINEM", buffer);
  // geometrically linear
  if (buffer == "linear") kintype_ = Inpar::Solid::KinemType::linear;
  // geometrically non-linear with Total Lagrangean approach
  else if (buffer == "nonlinear")
    kintype_ = Inpar::Solid::KinemType::nonlinearTotLag;
  else
    FOUR_C_THROW("Illegal KINEM type '%s'", buffer.c_str());

  // EAS type
  linedef->extract_string("EAS", buffer);
  if (buffer == "none")
  {
    iseas_ = false;
  }
  else if (buffer == "full")
  {
    iseas_ = true;

    if (num_node() == 9)
      FOUR_C_THROW("eas-technology not necessary with 9 nodes");
    else if (num_node() == 8)
      FOUR_C_THROW("eas-technology not necessary with 8 nodes");
    else if (num_node() == 3)
      FOUR_C_THROW("eas-technology not implemented for tri3 elements");
    else if (num_node() == 6)
      FOUR_C_THROW("eas-technology not implemented for tri6 elements");
    else
    {
      // EAS enhanced deformation gradient parameters
      Core::LinAlg::SerialDenseMatrix alpha(
          Wall1::neas_, 1);  // if you change '4' here, then do it for alphao as well
      Core::LinAlg::SerialDenseMatrix alphao(Wall1::neas_, 1);

      // EAS portion of internal forces, also called enhacement vector s or Rtilde
      Core::LinAlg::SerialDenseVector feas(Wall1::neas_);
      // EAS matrix K_{alpha alpha}, also called Dtilde
      Core::LinAlg::SerialDenseMatrix invKaa(Wall1::neas_, Wall1::neas_);
      // EAS matrix K_{d alpha}
      Core::LinAlg::SerialDenseMatrix Kda(2 * num_node(), Wall1::neas_);
      // EAS matrix K_{alpha d} // ONLY NEEDED FOR GENERALISED ENERGY-MOMENTUM METHOD
      Core::LinAlg::SerialDenseMatrix Kad(Wall1::neas_, 2 * num_node());
      // EAS increment over last Newton step
      Core::LinAlg::SerialDenseMatrix eas_inc(Wall1::neas_, 1);

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
    FOUR_C_THROW("Illegal EAS model");
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
  if (kintype_ == Inpar::Solid::KinemType::linear && iseas_ == true)
    FOUR_C_THROW("ERROR: No EAS for geometrically linear WALL element");

  // validate kinematics of solid material
  SolidMaterial()->ValidKinematics(kintype_);

  // Validate that materials doesn't use extended update call.
  if (SolidMaterial()->UsesExtendedUpdate())
    FOUR_C_THROW("This element currently does not support the extended update call.");

  return true;
}

/*----------------------------------------------------------------------*
 |  Get gaussrule on dependance of gausspoints                     mgit |
 *----------------------------------------------------------------------*/
Core::FE::GaussRule2D Discret::ELEMENTS::Wall1::get_gaussrule(int* ngp)
{
  Core::FE::GaussRule2D rule = Core::FE::GaussRule2D::undefined;

  switch (Shape())
  {
    case Core::FE::CellType::quad4:
    case Core::FE::CellType::quad8:
    case Core::FE::CellType::quad9:
    {
      if ((ngp[0] == 2) && (ngp[1] == 2))
      {
        rule = Core::FE::GaussRule2D::quad_4point;
      }
      else if ((ngp[0] == 3) && (ngp[1] == 3))
      {
        rule = Core::FE::GaussRule2D::quad_9point;
      }
      else
        FOUR_C_THROW("Unknown number of Gauss points for quad element");
      break;
    }
    case Core::FE::CellType::nurbs4:
    case Core::FE::CellType::nurbs9:
    {
      if ((ngp[0] == 2) && (ngp[1] == 2))
      {
        rule = Core::FE::GaussRule2D::quad_4point;
      }
      else if ((ngp[0] == 3) && (ngp[1] == 3))
      {
        rule = Core::FE::GaussRule2D::quad_9point;
      }
      else if ((ngp[0] == 4) && (ngp[1] == 4))
      {
        rule = Core::FE::GaussRule2D::quad_16point;
      }
      else if ((ngp[0] == 5) && (ngp[1] == 5))
      {
        rule = Core::FE::GaussRule2D::quad_25point;
      }
      else if ((ngp[0] == 10) && (ngp[1] == 10))
      {
        rule = Core::FE::GaussRule2D::quad_100point;
      }
      else
        FOUR_C_THROW("Unknown number of Gauss points for nurbs element");
      break;
    }
    case Core::FE::CellType::tri3:
    case Core::FE::CellType::tri6:
    {
      if ((ngp[0] == 1) && (ngp[1] == 0))
      {
        rule = Core::FE::GaussRule2D::tri_1point;
      }
      else if ((ngp[0] == 3) && (ngp[1] == 0))
      {
        rule = Core::FE::GaussRule2D::tri_3point;
      }
      else if ((ngp[0] == 6) && (ngp[1] == 0))
      {
        rule = Core::FE::GaussRule2D::tri_6point;
      }
      else
        FOUR_C_THROW("Unknown number of Gauss points for tri element");
      break;
    }
    default:
      FOUR_C_THROW("Unknown distype");
      break;
  }
  return rule;
}

FOUR_C_NAMESPACE_CLOSE
