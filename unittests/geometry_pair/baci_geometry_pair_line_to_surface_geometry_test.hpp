/*----------------------------------------------------------------------*/
/*! \file

\brief Create the geometry for the unit tests.

\level 1
*/
// End doxygen header.


#ifndef BACI_GEOMETRY_PAIR_LINE_TO_SURFACE_GEOMETRY_TEST_HPP
#define BACI_GEOMETRY_PAIR_LINE_TO_SURFACE_GEOMETRY_TEST_HPP


#include "baci_beam3_reissner.hpp"
#include "baci_geometry_pair_element.hpp"
#include "baci_lib_element.hpp"

#include <memory>

namespace
{
  using namespace BACI;

  /**
   * Setup the surface geometry for the tri3 tests.
   */
  GEOMETRYPAIR::ElementData<GEOMETRYPAIR::t_tri3, double> XtestSetupTri3()
  {
    auto element_data = GEOMETRYPAIR::ElementData<GEOMETRYPAIR::t_tri3, double>();

    element_data.element_position_(0) = 0.;
    element_data.element_position_(1) = 0.;
    element_data.element_position_(2) = 0.;
    element_data.element_position_(3) = 1.;
    element_data.element_position_(4) = -0.5;
    element_data.element_position_(5) = 0.5;
    element_data.element_position_(6) = -0.1;
    element_data.element_position_(7) = 0.95;
    element_data.element_position_(8) = 0.;

    element_data.nodal_normals_(0) = -0.2627627396383057;
    element_data.nodal_normals_(1) = 0.0814482045510598;
    element_data.nodal_normals_(2) = 0.961416628019913;
    element_data.nodal_normals_(3) = -0.7190848597139832;
    element_data.nodal_normals_(4) = -0.0866854771055773;
    element_data.nodal_normals_(5) = 0.6894944471053408;
    element_data.nodal_normals_(6) = -0.6025696958211013;
    element_data.nodal_normals_(7) = 0.1931525301890987;
    element_data.nodal_normals_(8) = 0.7743396294647555;

    return element_data;
  }

  /**
   * Setup the surface geometry for the tri6 tests.
   */
  GEOMETRYPAIR::ElementData<GEOMETRYPAIR::t_tri6, double> XtestSetupTri6()
  {
    auto element_data = GEOMETRYPAIR::ElementData<GEOMETRYPAIR::t_tri6, double>();

    element_data.element_position_(0) = 0.;
    element_data.element_position_(1) = 0.;
    element_data.element_position_(2) = 0.;
    element_data.element_position_(3) = 1.;
    element_data.element_position_(4) = -0.5;
    element_data.element_position_(5) = 0.5;
    element_data.element_position_(6) = -0.1;
    element_data.element_position_(7) = 0.95;
    element_data.element_position_(8) = 0.;
    element_data.element_position_(9) = 0.7;
    element_data.element_position_(10) = -0.1;
    element_data.element_position_(11) = 0.;
    element_data.element_position_(12) = 0.5;
    element_data.element_position_(13) = 0.5;
    element_data.element_position_(14) = 0.2;
    element_data.element_position_(15) = 0.;
    element_data.element_position_(16) = 0.4;
    element_data.element_position_(17) = 0.;

    element_data.nodal_normals_(0) = 0.4647142046388283;
    element_data.nodal_normals_(1) = 0.07560665247369116;
    element_data.nodal_normals_(2) = 0.882226922117334;
    element_data.nodal_normals_(3) = -0.940585437337313;
    element_data.nodal_normals_(4) = -0.3017696640225994;
    element_data.nodal_normals_(5) = -0.155673070711231;
    element_data.nodal_normals_(6) = -0.4731211190964157;
    element_data.nodal_normals_(7) = 0.1945363466954595;
    element_data.nodal_normals_(8) = 0.859250846074265;
    element_data.nodal_normals_(9) = -0.840743980012363;
    element_data.nodal_normals_(10) = 0.01531668856144018;
    element_data.nodal_normals_(11) = 0.5412161852018873;
    element_data.nodal_normals_(12) = -0.7282283511923046;
    element_data.nodal_normals_(13) = -0.4445443397695543;
    element_data.nodal_normals_(14) = 0.5215973528485246;
    element_data.nodal_normals_(15) = -0.1079357859495496;
    element_data.nodal_normals_(16) = -0.06091161819004808;
    element_data.nodal_normals_(17) = 0.992290099154941;

    return element_data;
  }

  /**
   * Setup the surface geometry for the quad4 tests.
   */
  GEOMETRYPAIR::ElementData<GEOMETRYPAIR::t_quad4, double> XtestSetupQuad4()
  {
    auto element_data = GEOMETRYPAIR::ElementData<GEOMETRYPAIR::t_quad4, double>();

    element_data.element_position_(0) = 0;
    element_data.element_position_(1) = 0;
    element_data.element_position_(2) = 0;
    element_data.element_position_(3) = 1.0000000000000000000;
    element_data.element_position_(4) = -0.50000000000000000000;
    element_data.element_position_(5) = 0.50000000000000000000;
    element_data.element_position_(6) = 1.2;
    element_data.element_position_(7) = 1.2;
    element_data.element_position_(8) = 0.5;
    element_data.element_position_(9) = -0.1;
    element_data.element_position_(10) = 0.95;
    element_data.element_position_(11) = 0;

    element_data.nodal_normals_(0) = -0.2627627396383057;
    element_data.nodal_normals_(1) = 0.0814482045510598;
    element_data.nodal_normals_(2) = 0.961416628019913;
    element_data.nodal_normals_(3) = -0.696398235712595;
    element_data.nodal_normals_(4) = -0.00250557448550925;
    element_data.nodal_normals_(5) = 0.7176511822556153;
    element_data.nodal_normals_(6) = -0.5381757744868267;
    element_data.nodal_normals_(7) = 0.2523235485556319;
    element_data.nodal_normals_(8) = 0.804176387740773;
    element_data.nodal_normals_(9) = -0.5695583742587585;
    element_data.nodal_normals_(10) = 0.3920519545501438;
    element_data.nodal_normals_(11) = 0.7224254447658469;

    return element_data;
  }

  /**
   * Setup the surface geometry for the quad8 tests.
   */
  GEOMETRYPAIR::ElementData<GEOMETRYPAIR::t_quad8, double> XtestSetupQuad8()
  {
    auto element_data = GEOMETRYPAIR::ElementData<GEOMETRYPAIR::t_quad8, double>();

    element_data.element_position_(0) = 0.;
    element_data.element_position_(1) = 0.;
    element_data.element_position_(2) = 0.;
    element_data.element_position_(3) = 1.;
    element_data.element_position_(4) = -0.5;
    element_data.element_position_(5) = 0.5;
    element_data.element_position_(6) = 1.2;
    element_data.element_position_(7) = 1.2;
    element_data.element_position_(8) = 0.5;
    element_data.element_position_(9) = -0.1;
    element_data.element_position_(10) = 0.95;
    element_data.element_position_(11) = 0.;
    element_data.element_position_(12) = 0.7;
    element_data.element_position_(13) = -0.1;
    element_data.element_position_(14) = 0.;
    element_data.element_position_(15) = 1.5;
    element_data.element_position_(16) = 0.5;
    element_data.element_position_(17) = 0.4285714285714286;
    element_data.element_position_(18) = 0.6;
    element_data.element_position_(19) = 1.;
    element_data.element_position_(20) = 0.;
    element_data.element_position_(21) = 0.;
    element_data.element_position_(22) = 0.4;
    element_data.element_position_(23) = 0.;

    element_data.nodal_normals_(0) = 0.4647142046388283;
    element_data.nodal_normals_(1) = 0.07560665247369116;
    element_data.nodal_normals_(2) = 0.882226922117334;
    element_data.nodal_normals_(3) = -0.836343726511398;
    element_data.nodal_normals_(4) = 0.4199949186410029;
    element_data.nodal_normals_(5) = 0.3523257575607623;
    element_data.nodal_normals_(6) = -0.685419065478165;
    element_data.nodal_normals_(7) = -0.2637809647101298;
    element_data.nodal_normals_(8) = 0.6786901408858333;
    element_data.nodal_normals_(9) = 0.084751268287718;
    element_data.nodal_normals_(10) = 0.5628096661844073;
    element_data.nodal_normals_(11) = 0.822230200231674;
    element_data.nodal_normals_(12) = -0.4965263011922144;
    element_data.nodal_normals_(13) = -0.01413350221776395;
    element_data.nodal_normals_(14) = 0.867906605770136;
    element_data.nodal_normals_(15) = -0.979140049557007;
    element_data.nodal_normals_(16) = 0.05236426002737823;
    element_data.nodal_normals_(17) = 0.1963230695188068;
    element_data.nodal_normals_(18) = 0.2410702660375203;
    element_data.nodal_normals_(19) = -0.349981608608912;
    element_data.nodal_normals_(20) = 0.905206054149064;
    element_data.nodal_normals_(21) = 0.5745894251085984;
    element_data.nodal_normals_(22) = -0.2036034941085485;
    element_data.nodal_normals_(23) = 0.792712185941506;

    return element_data;
  }

  /**
   * Setup the surface geometry for the quad9 tests.
   */
  GEOMETRYPAIR::ElementData<GEOMETRYPAIR::t_quad9, double> XtestSetupQuad9()
  {
    auto element_data = GEOMETRYPAIR::ElementData<GEOMETRYPAIR::t_quad9, double>();

    element_data.element_position_(0) = 0.;
    element_data.element_position_(1) = 0.;
    element_data.element_position_(2) = 0.;
    element_data.element_position_(3) = 1.;
    element_data.element_position_(4) = -0.5;
    element_data.element_position_(5) = 0.5;
    element_data.element_position_(6) = 1.2;
    element_data.element_position_(7) = 1.2;
    element_data.element_position_(8) = 0.5;
    element_data.element_position_(9) = -0.1;
    element_data.element_position_(10) = 0.95;
    element_data.element_position_(11) = 0.;
    element_data.element_position_(12) = 0.7;
    element_data.element_position_(13) = -0.1;
    element_data.element_position_(14) = 0.;
    element_data.element_position_(15) = 1.5;
    element_data.element_position_(16) = 0.5;
    element_data.element_position_(17) = 0.4285714285714286;
    element_data.element_position_(18) = 0.6;
    element_data.element_position_(19) = 1.;
    element_data.element_position_(20) = 0.;
    element_data.element_position_(21) = 0.;
    element_data.element_position_(22) = 0.4;
    element_data.element_position_(23) = 0.;
    element_data.element_position_(24) = 0.5;
    element_data.element_position_(25) = 0.5;
    element_data.element_position_(26) = 0.2;

    element_data.nodal_normals_(0) = 0.4647142046388283;
    element_data.nodal_normals_(1) = 0.07560665247369116;
    element_data.nodal_normals_(2) = 0.882226922117334;
    element_data.nodal_normals_(3) = -0.836343726511398;
    element_data.nodal_normals_(4) = 0.4199949186410029;
    element_data.nodal_normals_(5) = 0.3523257575607623;
    element_data.nodal_normals_(6) = -0.685419065478165;
    element_data.nodal_normals_(7) = -0.2637809647101298;
    element_data.nodal_normals_(8) = 0.6786901408858333;
    element_data.nodal_normals_(9) = 0.084751268287718;
    element_data.nodal_normals_(10) = 0.5628096661844073;
    element_data.nodal_normals_(11) = 0.822230200231674;
    element_data.nodal_normals_(12) = -0.6102068773340782;
    element_data.nodal_normals_(13) = -0.5672458156179746;
    element_data.nodal_normals_(14) = 0.5530639669315768;
    element_data.nodal_normals_(15) = -0.477142251285473;
    element_data.nodal_normals_(16) = -0.03320684296047792;
    element_data.nodal_normals_(17) = 0.878198484181582;
    element_data.nodal_normals_(18) = 0.0853062914327283;
    element_data.nodal_normals_(19) = 0.898922710199504;
    element_data.nodal_normals_(20) = 0.4297217678097922;
    element_data.nodal_normals_(21) = -0.3719814588140768;
    element_data.nodal_normals_(22) = -0.3844673823348876;
    element_data.nodal_normals_(23) = 0.844875509302472;
    element_data.nodal_normals_(24) = 0.2831165159078166;
    element_data.nodal_normals_(25) = 0.1617600702632027;
    element_data.nodal_normals_(26) = 0.945345819310935;

    return element_data;
  }

  /**
   * Setup the beam geometry for the tests.
   */
  std::pair<std::shared_ptr<DRT::Element>,
      GEOMETRYPAIR::ElementData<GEOMETRYPAIR::t_hermite, double>>
  XtestSetupBeam()
  {
    // Set up the beam element.
    const int dummy_node_ids[2] = {0, 1};
    std::shared_ptr<DRT::Element> element = std::make_shared<DRT::ELEMENTS::Beam3r>(0, 0);
    element->SetNodeIds(2, dummy_node_ids);

    // Set up the beam.
    auto element_data = GEOMETRYPAIR::ElementData<GEOMETRYPAIR::t_hermite, double>();

    element_data.shape_function_data_.ref_length_ = 1.807519343263254585;
    element_data.element_position_(0) = -0.1;
    element_data.element_position_(1) = 0.1;
    element_data.element_position_(2) = 0.1;
    element_data.element_position_(3) = 0.993883734673619;
    element_data.element_position_(4) = -0.1104315260748465;
    element_data.element_position_(5) = 0.;
    element_data.element_position_(6) = 1.5;
    element_data.element_position_(7) = 0.5;
    element_data.element_position_(8) = 0.75;
    element_data.element_position_(9) = 0.975900072948533;
    element_data.element_position_(10) = 0.1951800145897066;
    element_data.element_position_(11) = 0.0975900072948533;

    return {element, element_data};
  }

}  // namespace

#endif
