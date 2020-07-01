/*----------------------------------------------------------------------*/
/*! \file

\brief Create the geometry for the unit tests.

\level 1
*/
// End doxygen header.


#ifndef UNIT_GEOMETRY_PAIR_LINE_TO_SURFACE_GEOMETRY_H_
#define UNIT_GEOMETRY_PAIR_LINE_TO_SURFACE_GEOMETRY_H_

#include "unit_geometry_pair_line_to_surface.H"


namespace GEOMETRYPAIR
{
  /**
   * Setup the surface geometry for the tri4 tests.
   */
  void xtest_setup_tri3(
      LINALG::Matrix<9, 1, double>& q_solid, LINALG::Matrix<9, 1, double>* nodal_normals_ptr = NULL)
  {
    q_solid(0) = 0.;
    q_solid(1) = 0.;
    q_solid(2) = 0.;
    q_solid(3) = 1.;
    q_solid(4) = -0.5;
    q_solid(5) = 0.5;
    q_solid(6) = -0.1;
    q_solid(7) = 0.95;
    q_solid(8) = 0.;

    if (nodal_normals_ptr != NULL)
    {
      auto& nodal_normals = *nodal_normals_ptr;
      nodal_normals(0) = -0.2627627396383057;
      nodal_normals(1) = 0.0814482045510598;
      nodal_normals(2) = 0.961416628019913;
      nodal_normals(3) = -0.7190848597139832;
      nodal_normals(4) = -0.0866854771055773;
      nodal_normals(5) = 0.6894944471053408;
      nodal_normals(6) = -0.6025696958211013;
      nodal_normals(7) = 0.1931525301890987;
      nodal_normals(8) = 0.7743396294647555;
    }
  }

  /**
   * Setup the surface geometry for the tri6 tests.
   */
  void xtest_setup_tri6(LINALG::Matrix<18, 1, double>& q_solid,
      LINALG::Matrix<18, 1, double>* nodal_normals_ptr = NULL)
  {
    q_solid(0) = 0.;
    q_solid(1) = 0.;
    q_solid(2) = 0.;
    q_solid(3) = 1.;
    q_solid(4) = -0.5;
    q_solid(5) = 0.5;
    q_solid(6) = -0.1;
    q_solid(7) = 0.95;
    q_solid(8) = 0.;
    q_solid(9) = 0.7;
    q_solid(10) = -0.1;
    q_solid(11) = 0.;
    q_solid(12) = 0.5;
    q_solid(13) = 0.5;
    q_solid(14) = 0.2;
    q_solid(15) = 0.;
    q_solid(16) = 0.4;
    q_solid(17) = 0.;

    if (nodal_normals_ptr != NULL)
    {
      auto& nodal_normals = *nodal_normals_ptr;
      nodal_normals(0) = 0.4647142046388283;
      nodal_normals(1) = 0.07560665247369116;
      nodal_normals(2) = 0.882226922117334;
      nodal_normals(3) = -0.940585437337313;
      nodal_normals(4) = -0.3017696640225994;
      nodal_normals(5) = -0.155673070711231;
      nodal_normals(6) = -0.4731211190964157;
      nodal_normals(7) = 0.1945363466954595;
      nodal_normals(8) = 0.859250846074265;
      nodal_normals(9) = -0.840743980012363;
      nodal_normals(10) = 0.01531668856144018;
      nodal_normals(11) = 0.5412161852018873;
      nodal_normals(12) = -0.7282283511923046;
      nodal_normals(13) = -0.4445443397695543;
      nodal_normals(14) = 0.5215973528485246;
      nodal_normals(15) = -0.1079357859495496;
      nodal_normals(16) = -0.06091161819004808;
      nodal_normals(17) = 0.992290099154941;
    }
  }

  /**
   * Setup the surface geometry for the quad4 tests.
   */
  void xtest_setup_quad4(LINALG::Matrix<12, 1, double>& q_solid,
      LINALG::Matrix<12, 1, double>* nodal_normals_ptr = NULL)
  {
    q_solid(0) = 0;
    q_solid(1) = 0;
    q_solid(2) = 0;
    q_solid(3) = 1.0000000000000000000;
    q_solid(4) = -0.50000000000000000000;
    q_solid(5) = 0.50000000000000000000;
    q_solid(6) = 1.2;
    q_solid(7) = 1.2;
    q_solid(8) = 0.5;
    q_solid(9) = -0.1;
    q_solid(10) = 0.95;
    q_solid(11) = 0;

    if (nodal_normals_ptr != NULL)
    {
      auto& nodal_normals = *nodal_normals_ptr;
      nodal_normals(0) = -0.2627627396383057;
      nodal_normals(1) = 0.0814482045510598;
      nodal_normals(2) = 0.961416628019913;
      nodal_normals(3) = -0.696398235712595;
      nodal_normals(4) = -0.00250557448550925;
      nodal_normals(5) = 0.7176511822556153;
      nodal_normals(6) = -0.5381757744868267;
      nodal_normals(7) = 0.2523235485556319;
      nodal_normals(8) = 0.804176387740773;
      nodal_normals(9) = -0.5695583742587585;
      nodal_normals(10) = 0.3920519545501438;
      nodal_normals(11) = 0.7224254447658469;
    }
  }

  /**
   * Setup the surface geometry for the quad8 tests.
   */
  void xtest_setup_quad8(LINALG::Matrix<24, 1, double>& q_solid,
      LINALG::Matrix<24, 1, double>* nodal_normals_ptr = NULL)
  {
    q_solid(0) = 0.;
    q_solid(1) = 0.;
    q_solid(2) = 0.;
    q_solid(3) = 1.;
    q_solid(4) = -0.5;
    q_solid(5) = 0.5;
    q_solid(6) = 1.2;
    q_solid(7) = 1.2;
    q_solid(8) = 0.5;
    q_solid(9) = -0.1;
    q_solid(10) = 0.95;
    q_solid(11) = 0.;
    q_solid(12) = 0.7;
    q_solid(13) = -0.1;
    q_solid(14) = 0.;
    q_solid(15) = 1.5;
    q_solid(16) = 0.5;
    q_solid(17) = 0.4285714285714286;
    q_solid(18) = 0.6;
    q_solid(19) = 1.;
    q_solid(20) = 0.;
    q_solid(21) = 0.;
    q_solid(22) = 0.4;
    q_solid(23) = 0.;

    if (nodal_normals_ptr != NULL)
    {
      auto& nodal_normals = *nodal_normals_ptr;
      nodal_normals(0) = 0.4647142046388283;
      nodal_normals(1) = 0.07560665247369116;
      nodal_normals(2) = 0.882226922117334;
      nodal_normals(3) = -0.836343726511398;
      nodal_normals(4) = 0.4199949186410029;
      nodal_normals(5) = 0.3523257575607623;
      nodal_normals(6) = -0.685419065478165;
      nodal_normals(7) = -0.2637809647101298;
      nodal_normals(8) = 0.6786901408858333;
      nodal_normals(9) = 0.084751268287718;
      nodal_normals(10) = 0.5628096661844073;
      nodal_normals(11) = 0.822230200231674;
      nodal_normals(12) = -0.4965263011922144;
      nodal_normals(13) = -0.01413350221776395;
      nodal_normals(14) = 0.867906605770136;
      nodal_normals(15) = -0.979140049557007;
      nodal_normals(16) = 0.05236426002737823;
      nodal_normals(17) = 0.1963230695188068;
      nodal_normals(18) = 0.2410702660375203;
      nodal_normals(19) = -0.349981608608912;
      nodal_normals(20) = 0.905206054149064;
      nodal_normals(21) = 0.5745894251085984;
      nodal_normals(22) = -0.2036034941085485;
      nodal_normals(23) = 0.792712185941506;
    }
  }

  /**
   * Setup the surface geometry for the quad9 tests.
   */
  void xtest_setup_quad9(LINALG::Matrix<27, 1, double>& q_solid,
      LINALG::Matrix<27, 1, double>* nodal_normals_ptr = NULL)
  {
    q_solid(0) = 0.;
    q_solid(1) = 0.;
    q_solid(2) = 0.;
    q_solid(3) = 1.;
    q_solid(4) = -0.5;
    q_solid(5) = 0.5;
    q_solid(6) = 1.2;
    q_solid(7) = 1.2;
    q_solid(8) = 0.5;
    q_solid(9) = -0.1;
    q_solid(10) = 0.95;
    q_solid(11) = 0.;
    q_solid(12) = 0.7;
    q_solid(13) = -0.1;
    q_solid(14) = 0.;
    q_solid(15) = 1.5;
    q_solid(16) = 0.5;
    q_solid(17) = 0.4285714285714286;
    q_solid(18) = 0.6;
    q_solid(19) = 1.;
    q_solid(20) = 0.;
    q_solid(21) = 0.;
    q_solid(22) = 0.4;
    q_solid(23) = 0.;
    q_solid(24) = 0.5;
    q_solid(25) = 0.5;
    q_solid(26) = 0.2;

    if (nodal_normals_ptr != NULL)
    {
      auto& nodal_normals = *nodal_normals_ptr;
      nodal_normals(0) = 0.4647142046388283;
      nodal_normals(1) = 0.07560665247369116;
      nodal_normals(2) = 0.882226922117334;
      nodal_normals(3) = -0.836343726511398;
      nodal_normals(4) = 0.4199949186410029;
      nodal_normals(5) = 0.3523257575607623;
      nodal_normals(6) = -0.685419065478165;
      nodal_normals(7) = -0.2637809647101298;
      nodal_normals(8) = 0.6786901408858333;
      nodal_normals(9) = 0.084751268287718;
      nodal_normals(10) = 0.5628096661844073;
      nodal_normals(11) = 0.822230200231674;
      nodal_normals(12) = -0.6102068773340782;
      nodal_normals(13) = -0.5672458156179746;
      nodal_normals(14) = 0.5530639669315768;
      nodal_normals(15) = -0.477142251285473;
      nodal_normals(16) = -0.03320684296047792;
      nodal_normals(17) = 0.878198484181582;
      nodal_normals(18) = 0.0853062914327283;
      nodal_normals(19) = 0.898922710199504;
      nodal_normals(20) = 0.4297217678097922;
      nodal_normals(21) = -0.3719814588140768;
      nodal_normals(22) = -0.3844673823348876;
      nodal_normals(23) = 0.844875509302472;
      nodal_normals(24) = 0.2831165159078166;
      nodal_normals(25) = 0.1617600702632027;
      nodal_normals(26) = 0.945345819310935;
    }
  }

  /**
   * Setup the beam geometry for the tests.
   */
  void xtest_setup_beam(
      Teuchos::RCP<DRT::Element>& element_1, LINALG::Matrix<12, 1, double>& q_beam)
  {
    // Set up the beam.
    const int dummy_node_ids[2] = {0, 1};
    element_1 = Teuchos::rcp(new DRT::ELEMENTS::Beam3r(0, 0));
    element_1->SetNodeIds(2, dummy_node_ids);
    Teuchos::RCP<DRT::ELEMENTS::Beam3r> beam_element =
        Teuchos::rcp_dynamic_cast<DRT::ELEMENTS::Beam3r>(element_1, true);
    beam_element->SetCenterlineHermite(true);

    // Set the reference geometry of the beam.
    std::vector<double> xref(6);
    xref = {-0.1, 0.1, 0.1, 3.0 / 2.0, 1.0 / 2.0, 3.0 / 4.0};
    std::vector<double> rotref(9, 0.0);
    rotref[0] = 0.0;
    rotref[1] = 0.0;
    rotref[2] = -0.1106572211738956;
    rotref[3] = 0.0;
    rotref[4] = -0.04948741374921195;
    rotref[5] = 0.0989748274984239;
    beam_element->SetUpReferenceGeometry<3, 2, 2>(xref, rotref);

    // Define the coordinates for the beam element.
    q_beam(0) = -0.1;
    q_beam(1) = 0.1;
    q_beam(2) = 0.1;
    q_beam(3) = 0.993883734673619;
    q_beam(4) = -0.1104315260748465;
    q_beam(5) = 0.;
    q_beam(6) = 1.5;
    q_beam(7) = 0.5;
    q_beam(8) = 0.75;
    q_beam(9) = 0.975900072948533;
    q_beam(10) = 0.1951800145897066;
    q_beam(11) = 0.0975900072948533;
  }

}  // namespace GEOMETRYPAIR



#endif
