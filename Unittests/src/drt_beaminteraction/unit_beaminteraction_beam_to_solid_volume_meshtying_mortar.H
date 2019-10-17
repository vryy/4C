/*----------------------------------------------------------------------*/
/*! \file

\brief Unit tests for the beam to volume meshtying mortar pair.

\level 3
\maintainer Ivo Steinbrecher
*/
// End doxygen header.


#ifndef UNIT_BEAMINTERACTION_BEAM_TO_SOLIV_VOLUME_MESHTYING_MORTAR_H_
#define UNIT_BEAMINTERACTION_BEAM_TO_SOLIV_VOLUME_MESHTYING_MORTAR_H_


#include <cxxtest/TestSuite.h>

#include "src/drt_beam3/beam3r.H"
#include "src/drt_so3/so_hex8.H"
#include "src/drt_so3/so_hex27.H"
#include "src/drt_so3/so_tet10.H"
#include "src/drt_geometry_pair/geometry_pair_element.H"
#include "src/drt_geometry_pair/geometry_pair_line_to_volume_segmentation.H"
#include "src/drt_geometry_pair/geometry_pair_line_to_3D_evaluation_data.H"
#include "src/drt_geometry_pair/geometry_pair_utility_classes.H"
#include "src/drt_beaminteraction/beam_to_solid_volume_meshtying_pair_mortar.H"
#include "src/drt_lib/drt_globalproblem.H"
#include "src/linalg/linalg_serialdensevector.H"


namespace BEAMINTERACTION
{
  class BeamToSolidVolumeMeshtyingPairMortar_TestSuite;
}

/**
 * Class to test the local mortar matrices calculated by the beam to volume mesh tying mortar pair.
 */
class BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar_TestSuite : public CxxTest::TestSuite
{
 public:
  /**
   * \brief Set up the testing environment.
   */
  void setUp()
  {
    // Set up the evaluation data container for the geometry pairs.
    Teuchos::ParameterList line_to_volume_params_list;
    INPAR::GEOMETRYPAIR::SetValidParametersLineTo3D(line_to_volume_params_list);
    evaluation_data_ =
        Teuchos::rcp(new GEOMETRYPAIR::LineTo3DEvaluationData(line_to_volume_params_list));
  }

  /**
   * \brief Delete pointers and other class variables.
   */
  void tearDown()
  {
    // Dereference the pointers.
    evaluation_data_ = Teuchos::null;
  }

  /**
   * \brief Test a non straight beam in a hex8 element, with line2 mortar shape functions.
   */
  void test_beam_to_solid_meshtying_mortar_hermite3_hex8_line2()
  {
    // Element types.
    typedef GEOMETRYPAIR::t_hermite beam_type;
    typedef GEOMETRYPAIR::t_hex8 solid_type;
    typedef GEOMETRYPAIR::t_line2 lambda_type;

    // Create the mesh tying mortar pair.
    BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<beam_type, solid_type, lambda_type>
        contact_pair = BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<beam_type, solid_type,
            lambda_type>();

    // Definition of variables for this test case.
    LINALG::Matrix<beam_type::n_dof_, 1, double> q_beam;
    LINALG::Matrix<9, 1, double> q_beam_rot;
    LINALG::Matrix<solid_type::n_dof_, 1, double> q_solid;
    LINALG::SerialDenseMatrix local_D;
    LINALG::SerialDenseMatrix local_M;
    LINALG::SerialDenseVector local_kappa;

    // Matrices for the results.
    LINALG::Matrix<lambda_type::n_dof_, beam_type::n_dof_, double> result_local_D(true);
    LINALG::Matrix<lambda_type::n_dof_, solid_type::n_dof_, double> result_local_M(true);
    LINALG::Matrix<lambda_type::n_dof_, 1, double> result_local_kappa(true);

    // Define the geometry of the two elements.
    q_beam(0) = 0.15;
    q_beam(1) = 0.2;
    q_beam(2) = 0.3;
    q_beam(3) = 0.5773502691896255;
    q_beam(4) = 0.5773502691896258;
    q_beam(5) = 0.577350269189626;
    q_beam(6) = 0.65;
    q_beam(7) = 0.1;
    q_beam(8) = 0.1;
    q_beam(9) = 0.8017837257372733;
    q_beam(10) = -0.5345224838248488;
    q_beam(11) = 0.2672612419124244;

    q_beam_rot(0) = 1.674352746442651;
    q_beam_rot(1) = 0.1425949677148126;
    q_beam_rot(2) = 1.0831163124736984;
    q_beam_rot(3) = 1.4331999091513161;
    q_beam_rot(4) = -0.6560404572957742;
    q_beam_rot(5) = -0.2491376152457331;
    q_beam_rot(6) = 0.0;
    q_beam_rot(7) = 0.0;
    q_beam_rot(8) = 0.0;

    // Positional DOFs of the solid.
    q_solid(0) = -0.954193516126594;
    q_solid(1) = -0.975482672114534;
    q_solid(2) = -1.00311316662815;
    q_solid(3) = 0.923762409575691;
    q_solid(4) = -1.010615727466753;
    q_solid(5) = -1.014160992102648;
    q_solid(6) = 0.905563488894572;
    q_solid(7) = 1.069973754602123;
    q_solid(8) = -0.935488788632622;
    q_solid(9) = -1.047838992508937;
    q_solid(10) = 1.08888017769104;
    q_solid(11) = -1.036911970151991;
    q_solid(12) = -1.089167378981223;
    q_solid(13) = -1.063820485514786;
    q_solid(14) = 1.081000671652452;
    q_solid(15) = 0.972887477248942;
    q_solid(16) = -1.008995039455167;
    q_solid(17) = 0.923459231565508;
    q_solid(18) = 1.089192110164642;
    q_solid(19) = 1.026798931538616;
    q_solid(20) = 0.964761788314679;
    q_solid(21) = -0.941489453299989;
    q_solid(22) = 0.954327975326486;
    q_solid(23) = 0.963168067937965;

    // Results for D.
    result_local_D(0, 0) = 0.1945989639134667;
    result_local_D(0, 3) = 0.01783046400866768;
    result_local_D(0, 6) = 0.0989060324956216;
    result_local_D(0, 9) = -0.01321892508537963;
    result_local_D(1, 1) = 0.1945989639134667;
    result_local_D(1, 4) = 0.01783046400866768;
    result_local_D(1, 7) = 0.0989060324956216;
    result_local_D(1, 10) = -0.01321892508537963;
    result_local_D(2, 2) = 0.1945989639134667;
    result_local_D(2, 5) = 0.01783046400866768;
    result_local_D(2, 8) = 0.0989060324956216;
    result_local_D(2, 11) = -0.01321892508537963;
    result_local_D(3, 0) = 0.0947117421336144;
    result_local_D(3, 3) = 0.01321892508537963;
    result_local_D(3, 6) = 0.2287572183685851;
    result_local_D(3, 9) = -0.02037627478307643;
    result_local_D(4, 1) = 0.0947117421336144;
    result_local_D(4, 4) = 0.01321892508537963;
    result_local_D(4, 7) = 0.2287572183685851;
    result_local_D(4, 10) = -0.02037627478307643;
    result_local_D(5, 2) = 0.0947117421336144;
    result_local_D(5, 5) = 0.01321892508537963;
    result_local_D(5, 8) = 0.2287572183685851;
    result_local_D(5, 11) = -0.02037627478307643;

    // Results for M.
    result_local_M(0, 0) = 0.01399734051851155;
    result_local_M(0, 3) = 0.02819625621789035;
    result_local_M(0, 6) = 0.04322611636999435;
    result_local_M(0, 9) = 0.02172860974957874;
    result_local_M(0, 12) = 0.02492453660501876;
    result_local_M(0, 15) = 0.04814578912323388;
    result_local_M(0, 18) = 0.0743949977000356;
    result_local_M(0, 21) = 0.03889135012482514;
    result_local_M(1, 1) = 0.01399734051851155;
    result_local_M(1, 4) = 0.02819625621789035;
    result_local_M(1, 7) = 0.04322611636999435;
    result_local_M(1, 10) = 0.02172860974957874;
    result_local_M(1, 13) = 0.02492453660501876;
    result_local_M(1, 16) = 0.04814578912323388;
    result_local_M(1, 19) = 0.0743949977000356;
    result_local_M(1, 22) = 0.03889135012482514;
    result_local_M(2, 2) = 0.01399734051851155;
    result_local_M(2, 5) = 0.02819625621789035;
    result_local_M(2, 8) = 0.04322611636999435;
    result_local_M(2, 11) = 0.02172860974957874;
    result_local_M(2, 14) = 0.02492453660501876;
    result_local_M(2, 17) = 0.04814578912323388;
    result_local_M(2, 20) = 0.0743949977000356;
    result_local_M(2, 23) = 0.03889135012482514;
    result_local_M(3, 0) = 0.01369921717609965;
    result_local_M(3, 3) = 0.04153192175012023;
    result_local_M(3, 6) = 0.0578615874518087;
    result_local_M(3, 9) = 0.01967492295908918;
    result_local_M(3, 12) = 0.02020829336851376;
    result_local_M(3, 15) = 0.05849674525787201;
    result_local_M(3, 18) = 0.0825640019910821;
    result_local_M(3, 21) = 0.02943227054761388;
    result_local_M(4, 1) = 0.01369921717609965;
    result_local_M(4, 4) = 0.04153192175012023;
    result_local_M(4, 7) = 0.0578615874518087;
    result_local_M(4, 10) = 0.01967492295908918;
    result_local_M(4, 13) = 0.02020829336851376;
    result_local_M(4, 16) = 0.05849674525787201;
    result_local_M(4, 19) = 0.0825640019910821;
    result_local_M(4, 22) = 0.02943227054761388;
    result_local_M(5, 2) = 0.01369921717609965;
    result_local_M(5, 5) = 0.04153192175012023;
    result_local_M(5, 8) = 0.0578615874518087;
    result_local_M(5, 11) = 0.01967492295908918;
    result_local_M(5, 14) = 0.02020829336851376;
    result_local_M(5, 17) = 0.05849674525787201;
    result_local_M(5, 20) = 0.0825640019910821;
    result_local_M(5, 23) = 0.02943227054761388;

    // Results for Kappa.
    result_local_kappa(0) = 0.2935049964090884;
    result_local_kappa(1) = 0.2935049964090884;
    result_local_kappa(2) = 0.2935049964090884;
    result_local_kappa(3) = 0.3234689605021996;
    result_local_kappa(4) = 0.3234689605021996;
    result_local_kappa(5) = 0.3234689605021996;

    // Perform the unit tests.
    PerformMortarPairUnitTest(contact_pair, q_beam, q_beam_rot, q_solid, result_local_D,
        result_local_M, result_local_kappa);
  }

  /**
   * \brief Test a non straight beam in a hex8 element, with line3 mortar shape functions.
   */
  void test_beam_to_solid_meshtying_mortar_hermite3_hex8_line3()
  {
    // Element types.
    typedef GEOMETRYPAIR::t_hermite beam_type;
    typedef GEOMETRYPAIR::t_hex8 solid_type;
    typedef GEOMETRYPAIR::t_line3 lambda_type;

    // Create the mesh tying mortar pair.
    BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<beam_type, solid_type, lambda_type>
        contact_pair = BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<beam_type, solid_type,
            lambda_type>();

    // Definition of variables for this test case.
    LINALG::Matrix<beam_type::n_dof_, 1, double> q_beam;
    LINALG::Matrix<9, 1, double> q_beam_rot;
    LINALG::Matrix<solid_type::n_dof_, 1, double> q_solid;
    LINALG::SerialDenseMatrix local_D;
    LINALG::SerialDenseMatrix local_M;
    LINALG::SerialDenseVector local_kappa;

    // Matrices for the results.
    LINALG::Matrix<lambda_type::n_dof_, beam_type::n_dof_, double> result_local_D(true);
    LINALG::Matrix<lambda_type::n_dof_, solid_type::n_dof_, double> result_local_M(true);
    LINALG::Matrix<lambda_type::n_dof_, 1, double> result_local_kappa(true);

    // Define the geometry of the two elements.
    q_beam(0) = 0.15;
    q_beam(1) = 0.2;
    q_beam(2) = 0.3;
    q_beam(3) = 0.5773502691896255;
    q_beam(4) = 0.5773502691896258;
    q_beam(5) = 0.577350269189626;
    q_beam(6) = 0.65;
    q_beam(7) = 0.1;
    q_beam(8) = 0.1;
    q_beam(9) = 0.8017837257372733;
    q_beam(10) = -0.5345224838248488;
    q_beam(11) = 0.2672612419124244;

    q_beam_rot(0) = 1.674352746442651;
    q_beam_rot(1) = 0.1425949677148126;
    q_beam_rot(2) = 1.0831163124736984;
    q_beam_rot(3) = 1.4331999091513161;
    q_beam_rot(4) = -0.6560404572957742;
    q_beam_rot(5) = -0.2491376152457331;
    q_beam_rot(6) = 0.0;
    q_beam_rot(7) = 0.0;
    q_beam_rot(8) = 0.0;

    // Positional DOFs of the solid.
    q_solid(0) = -0.954193516126594;
    q_solid(1) = -0.975482672114534;
    q_solid(2) = -1.00311316662815;
    q_solid(3) = 0.923762409575691;
    q_solid(4) = -1.010615727466753;
    q_solid(5) = -1.014160992102648;
    q_solid(6) = 0.905563488894572;
    q_solid(7) = 1.069973754602123;
    q_solid(8) = -0.935488788632622;
    q_solid(9) = -1.047838992508937;
    q_solid(10) = 1.08888017769104;
    q_solid(11) = -1.036911970151991;
    q_solid(12) = -1.089167378981223;
    q_solid(13) = -1.063820485514786;
    q_solid(14) = 1.081000671652452;
    q_solid(15) = 0.972887477248942;
    q_solid(16) = -1.008995039455167;
    q_solid(17) = 0.923459231565508;
    q_solid(18) = 1.089192110164642;
    q_solid(19) = 1.026798931538616;
    q_solid(20) = 0.964761788314679;
    q_solid(21) = -0.941489453299989;
    q_solid(22) = 0.954327975326486;
    q_solid(23) = 0.963168067937965;

    // Results for D.
    result_local_D(0, 0) = 0.0936536895927479;
    result_local_D(0, 3) = 0.00502552312033731;
    result_local_D(0, 6) = -0.01315602206049943;
    result_local_D(0, 9) = 0.0004139841970492526;
    result_local_D(1, 1) = 0.0936536895927479;
    result_local_D(1, 4) = 0.00502552312033731;
    result_local_D(1, 7) = -0.01315602206049943;
    result_local_D(1, 10) = 0.0004139841970492526;
    result_local_D(2, 2) = 0.0936536895927479;
    result_local_D(2, 5) = 0.00502552312033731;
    result_local_D(2, 8) = -0.01315602206049943;
    result_local_D(2, 11) = 0.0004139841970492526;
    result_local_D(3, 0) = -0.006233532187104405;
    result_local_D(3, 3) = 0.0004139841970492533;
    result_local_D(3, 6) = 0.1166951638124641;
    result_local_D(3, 9) = -0.006743365500647556;
    result_local_D(4, 1) = -0.006233532187104405;
    result_local_D(4, 4) = 0.0004139841970492533;
    result_local_D(4, 7) = 0.1166951638124641;
    result_local_D(4, 10) = -0.006743365500647556;
    result_local_D(5, 2) = -0.006233532187104405;
    result_local_D(5, 5) = 0.0004139841970492533;
    result_local_D(5, 8) = 0.1166951638124641;
    result_local_D(5, 11) = -0.006743365500647556;
    result_local_D(6, 0) = 0.2018905486414377;
    result_local_D(6, 3) = 0.02560988177666075;
    result_local_D(6, 6) = 0.2241241091122421;
    result_local_D(6, 9) = -0.02726581856485776;
    result_local_D(7, 1) = 0.2018905486414377;
    result_local_D(7, 4) = 0.02560988177666075;
    result_local_D(7, 7) = 0.2241241091122421;
    result_local_D(7, 10) = -0.02726581856485776;
    result_local_D(8, 2) = 0.2018905486414377;
    result_local_D(8, 5) = 0.02560988177666075;
    result_local_D(8, 8) = 0.2241241091122421;
    result_local_D(8, 11) = -0.02726581856485776;

    // Results for M.
    result_local_M(0, 0) = 0.00442272260370199;
    result_local_M(0, 3) = 0.004818495175629803;
    result_local_M(0, 6) = 0.00835605925960412;
    result_local_M(0, 9) = 0.007159180407398056;
    result_local_M(0, 12) = 0.00937813704127027;
    result_local_M(0, 15) = 0.01178809762126675;
    result_local_M(0, 18) = 0.01956565251460088;
    result_local_M(0, 21) = 0.0150093229087766;
    result_local_M(1, 1) = 0.00442272260370199;
    result_local_M(1, 4) = 0.004818495175629803;
    result_local_M(1, 7) = 0.00835605925960412;
    result_local_M(1, 10) = 0.007159180407398056;
    result_local_M(1, 13) = 0.00937813704127027;
    result_local_M(1, 16) = 0.01178809762126675;
    result_local_M(1, 19) = 0.01956565251460088;
    result_local_M(1, 22) = 0.0150093229087766;
    result_local_M(2, 2) = 0.00442272260370199;
    result_local_M(2, 5) = 0.004818495175629803;
    result_local_M(2, 8) = 0.00835605925960412;
    result_local_M(2, 11) = 0.007159180407398056;
    result_local_M(2, 14) = 0.00937813704127027;
    result_local_M(2, 17) = 0.01178809762126675;
    result_local_M(2, 20) = 0.01956565251460088;
    result_local_M(2, 23) = 0.0150093229087766;
    result_local_M(3, 0) = 0.004124599261290101;
    result_local_M(3, 3) = 0.01815416070785969;
    result_local_M(3, 6) = 0.02299153034141846;
    result_local_M(3, 9) = 0.005105493616908499;
    result_local_M(3, 12) = 0.004661893804765268;
    result_local_M(3, 15) = 0.02213905375590488;
    result_local_M(3, 18) = 0.02773465680564743;
    result_local_M(3, 21) = 0.005550243331565347;
    result_local_M(4, 1) = 0.004124599261290101;
    result_local_M(4, 4) = 0.01815416070785969;
    result_local_M(4, 7) = 0.02299153034141846;
    result_local_M(4, 10) = 0.005105493616908499;
    result_local_M(4, 13) = 0.004661893804765268;
    result_local_M(4, 16) = 0.02213905375590488;
    result_local_M(4, 19) = 0.02773465680564743;
    result_local_M(4, 22) = 0.005550243331565347;
    result_local_M(5, 2) = 0.004124599261290101;
    result_local_M(5, 5) = 0.01815416070785969;
    result_local_M(5, 8) = 0.02299153034141846;
    result_local_M(5, 11) = 0.005105493616908499;
    result_local_M(5, 14) = 0.004661893804765268;
    result_local_M(5, 17) = 0.02213905375590488;
    result_local_M(5, 20) = 0.02773465680564743;
    result_local_M(5, 23) = 0.005550243331565347;
    result_local_M(6, 0) = 0.01914923582961911;
    result_local_M(6, 3) = 0.04675552208452109;
    result_local_M(6, 6) = 0.06974011422078046;
    result_local_M(6, 9) = 0.02913885868436137;
    result_local_M(6, 12) = 0.03109279912749698;
    result_local_M(6, 15) = 0.07271538300393426;
    result_local_M(6, 18) = 0.1096586903708694;
    result_local_M(6, 21) = 0.04776405443209706;
    result_local_M(7, 1) = 0.01914923582961911;
    result_local_M(7, 4) = 0.04675552208452109;
    result_local_M(7, 7) = 0.06974011422078046;
    result_local_M(7, 10) = 0.02913885868436137;
    result_local_M(7, 13) = 0.03109279912749698;
    result_local_M(7, 16) = 0.07271538300393426;
    result_local_M(7, 19) = 0.1096586903708694;
    result_local_M(7, 22) = 0.04776405443209706;
    result_local_M(8, 2) = 0.01914923582961911;
    result_local_M(8, 5) = 0.04675552208452109;
    result_local_M(8, 8) = 0.06974011422078046;
    result_local_M(8, 11) = 0.02913885868436137;
    result_local_M(8, 14) = 0.03109279912749698;
    result_local_M(8, 17) = 0.07271538300393426;
    result_local_M(8, 20) = 0.1096586903708694;
    result_local_M(8, 23) = 0.04776405443209706;

    // Results for Kappa.
    result_local_kappa(0) = 0.0804976675322485;
    result_local_kappa(1) = 0.0804976675322485;
    result_local_kappa(2) = 0.0804976675322485;
    result_local_kappa(3) = 0.1104616316253597;
    result_local_kappa(4) = 0.1104616316253597;
    result_local_kappa(5) = 0.1104616316253597;
    result_local_kappa(6) = 0.4260146577536798;
    result_local_kappa(7) = 0.4260146577536798;
    result_local_kappa(8) = 0.4260146577536798;

    // Perform the unit tests.
    PerformMortarPairUnitTest(contact_pair, q_beam, q_beam_rot, q_solid, result_local_D,
        result_local_M, result_local_kappa);
  }

  /**
   * \brief Test a non straight beam in a hex8 element, with line4 mortar shape functions.
   */
  void test_beam_to_solid_meshtying_mortar_hermite3_hex8_line4()
  {
    // Element types.
    typedef GEOMETRYPAIR::t_hermite beam_type;
    typedef GEOMETRYPAIR::t_hex8 solid_type;
    typedef GEOMETRYPAIR::t_line4 lambda_type;

    // Create the mesh tying mortar pair.
    BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<beam_type, solid_type, lambda_type>
        contact_pair = BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<beam_type, solid_type,
            lambda_type>();

    // Definition of variables for this test case.
    LINALG::Matrix<beam_type::n_dof_, 1, double> q_beam;
    LINALG::Matrix<9, 1, double> q_beam_rot;
    LINALG::Matrix<solid_type::n_dof_, 1, double> q_solid;
    LINALG::SerialDenseMatrix local_D;
    LINALG::SerialDenseMatrix local_M;
    LINALG::SerialDenseVector local_kappa;

    // Matrices for the results.
    LINALG::Matrix<lambda_type::n_dof_, beam_type::n_dof_, double> result_local_D(true);
    LINALG::Matrix<lambda_type::n_dof_, solid_type::n_dof_, double> result_local_M(true);
    LINALG::Matrix<lambda_type::n_dof_, 1, double> result_local_kappa(true);

    // Define the geometry of the two elements.
    q_beam(0) = 0.15;
    q_beam(1) = 0.2;
    q_beam(2) = 0.3;
    q_beam(3) = 0.5773502691896255;
    q_beam(4) = 0.5773502691896258;
    q_beam(5) = 0.577350269189626;
    q_beam(6) = 0.65;
    q_beam(7) = 0.1;
    q_beam(8) = 0.1;
    q_beam(9) = 0.8017837257372733;
    q_beam(10) = -0.5345224838248488;
    q_beam(11) = 0.2672612419124244;

    q_beam_rot(0) = 1.674352746442651;
    q_beam_rot(1) = 0.1425949677148126;
    q_beam_rot(2) = 1.0831163124736984;
    q_beam_rot(3) = 1.4331999091513161;
    q_beam_rot(4) = -0.6560404572957742;
    q_beam_rot(5) = -0.2491376152457331;
    q_beam_rot(6) = 0.0;
    q_beam_rot(7) = 0.0;
    q_beam_rot(8) = 0.0;

    // Positional DOFs of the solid.
    q_solid(0) = -0.954193516126594;
    q_solid(1) = -0.975482672114534;
    q_solid(2) = -1.00311316662815;
    q_solid(3) = 0.923762409575691;
    q_solid(4) = -1.010615727466753;
    q_solid(5) = -1.014160992102648;
    q_solid(6) = 0.905563488894572;
    q_solid(7) = 1.069973754602123;
    q_solid(8) = -0.935488788632622;
    q_solid(9) = -1.047838992508937;
    q_solid(10) = 1.08888017769104;
    q_solid(11) = -1.036911970151991;
    q_solid(12) = -1.089167378981223;
    q_solid(13) = -1.063820485514786;
    q_solid(14) = 1.081000671652452;
    q_solid(15) = 0.972887477248942;
    q_solid(16) = -1.008995039455167;
    q_solid(17) = 0.923459231565508;
    q_solid(18) = 1.089192110164642;
    q_solid(19) = 1.026798931538616;
    q_solid(20) = 0.964761788314679;
    q_solid(21) = -0.941489453299989;
    q_solid(22) = 0.954327975326486;
    q_solid(23) = 0.963168067937965;

    // Results for D.
    result_local_D(0, 0) = 0.05756355583405802;
    result_local_D(0, 3) = 0.001955220263108095;
    result_local_D(0, 6) = 0.005745348903101642;
    result_local_D(0, 9) = -0.0002830518321958758;
    result_local_D(1, 1) = 0.05756355583405802;
    result_local_D(1, 4) = 0.001955220263108095;
    result_local_D(1, 7) = 0.005745348903101642;
    result_local_D(1, 10) = -0.0002830518321958758;
    result_local_D(2, 2) = 0.05756355583405802;
    result_local_D(2, 5) = 0.001955220263108095;
    result_local_D(2, 8) = 0.005745348903101642;
    result_local_D(2, 11) = -0.0002830518321958758;
    result_local_D(3, 0) = 0.004620282991405782;
    result_local_D(3, 3) = 0.0002830518321958757;
    result_local_D(3, 6) = 0.06977827920983277;
    result_local_D(3, 9) = -0.002638102150795207;
    result_local_D(4, 1) = 0.004620282991405782;
    result_local_D(4, 4) = 0.0002830518321958757;
    result_local_D(4, 7) = 0.06977827920983277;
    result_local_D(4, 10) = -0.002638102150795207;
    result_local_D(5, 2) = 0.004620282991405782;
    result_local_D(5, 5) = 0.0002830518321958757;
    result_local_D(5, 8) = 0.06977827920983277;
    result_local_D(5, 11) = -0.002638102150795207;
    result_local_D(6, 0) = 0.1839793570166088;
    result_local_D(6, 3) = 0.01881461423793542;
    result_local_D(6, 6) = 0.02734242802628759;
    result_local_D(6, 9) = -0.00813357387408627;
    result_local_D(7, 1) = 0.1839793570166088;
    result_local_D(7, 4) = 0.01881461423793542;
    result_local_D(7, 7) = 0.02734242802628759;
    result_local_D(7, 10) = -0.00813357387408627;
    result_local_D(8, 2) = 0.1839793570166088;
    result_local_D(8, 5) = 0.01881461423793542;
    result_local_D(8, 8) = 0.02734242802628759;
    result_local_D(8, 11) = -0.00813357387408627;
    result_local_D(9, 0) = 0.04314751020500857;
    result_local_D(9, 3) = 0.00999650276080791;
    result_local_D(9, 6) = 0.2247971947249847;
    result_local_D(9, 9) = -0.0225404720113787;
    result_local_D(10, 1) = 0.04314751020500857;
    result_local_D(10, 4) = 0.00999650276080791;
    result_local_D(10, 7) = 0.2247971947249847;
    result_local_D(10, 10) = -0.0225404720113787;
    result_local_D(11, 2) = 0.04314751020500857;
    result_local_D(11, 5) = 0.00999650276080791;
    result_local_D(11, 8) = 0.2247971947249847;
    result_local_D(11, 11) = -0.0225404720113787;

    // Results for M.
    result_local_M(0, 0) = 0.003346229311968168;
    result_local_M(0, 3) = 0.005354343233335789;
    result_local_M(0, 6) = 0.007966326587387424;
    result_local_M(0, 9) = 0.00509144130802595;
    result_local_M(0, 12) = 0.006517401996246052;
    result_local_M(0, 15) = 0.00999150651667994;
    result_local_M(0, 18) = 0.01506543504315179;
    result_local_M(0, 21) = 0.00997622074036454;
    result_local_M(1, 1) = 0.003346229311968168;
    result_local_M(1, 4) = 0.005354343233335789;
    result_local_M(1, 7) = 0.007966326587387424;
    result_local_M(1, 10) = 0.00509144130802595;
    result_local_M(1, 13) = 0.006517401996246052;
    result_local_M(1, 16) = 0.00999150651667994;
    result_local_M(1, 19) = 0.01506543504315179;
    result_local_M(1, 22) = 0.00997622074036454;
    result_local_M(2, 2) = 0.003346229311968168;
    result_local_M(2, 5) = 0.005354343233335789;
    result_local_M(2, 8) = 0.007966326587387424;
    result_local_M(2, 11) = 0.00509144130802595;
    result_local_M(2, 14) = 0.006517401996246052;
    result_local_M(2, 17) = 0.00999150651667994;
    result_local_M(2, 20) = 0.01506543504315179;
    result_local_M(2, 23) = 0.00997622074036454;
    result_local_M(3, 0) = 0.002807438074321535;
    result_local_M(3, 3) = 0.01177387238958858;
    result_local_M(3, 6) = 0.01466374873603761;
    result_local_M(3, 9) = 0.003530875380735433;
    result_local_M(3, 12) = 0.00363602895885236;
    result_local_M(3, 15) = 0.01484622198499991;
    result_local_M(3, 18) = 0.01852753798073785;
    result_local_M(3, 21) = 0.004612838695965277;
    result_local_M(4, 1) = 0.002807438074321535;
    result_local_M(4, 4) = 0.01177387238958858;
    result_local_M(4, 7) = 0.01466374873603761;
    result_local_M(4, 10) = 0.003530875380735433;
    result_local_M(4, 13) = 0.00363602895885236;
    result_local_M(4, 16) = 0.01484622198499991;
    result_local_M(4, 19) = 0.01852753798073785;
    result_local_M(4, 22) = 0.004612838695965277;
    result_local_M(5, 2) = 0.002807438074321535;
    result_local_M(5, 5) = 0.01177387238958858;
    result_local_M(5, 8) = 0.01466374873603761;
    result_local_M(5, 11) = 0.003530875380735433;
    result_local_M(5, 14) = 0.00363602895885236;
    result_local_M(5, 17) = 0.01484622198499991;
    result_local_M(5, 20) = 0.01852753798073785;
    result_local_M(5, 23) = 0.004612838695965277;
    result_local_M(6, 0) = 0.01041044331130863;
    result_local_M(6, 3) = 0.01592577660857746;
    result_local_M(6, 6) = 0.02732174084944276;
    result_local_M(6, 9) = 0.01713028930475183;
    result_local_M(6, 12) = 0.02024200480788401;
    result_local_M(6, 15) = 0.0326580419402358;
    result_local_M(6, 18) = 0.05462266130342334;
    result_local_M(6, 21) = 0.0330108269172726;
    result_local_M(7, 1) = 0.01041044331130863;
    result_local_M(7, 4) = 0.01592577660857746;
    result_local_M(7, 7) = 0.02732174084944276;
    result_local_M(7, 10) = 0.01713028930475183;
    result_local_M(7, 13) = 0.02024200480788401;
    result_local_M(7, 16) = 0.0326580419402358;
    result_local_M(7, 19) = 0.05462266130342334;
    result_local_M(7, 22) = 0.0330108269172726;
    result_local_M(8, 2) = 0.01041044331130863;
    result_local_M(8, 5) = 0.01592577660857746;
    result_local_M(8, 8) = 0.02732174084944276;
    result_local_M(8, 11) = 0.01713028930475183;
    result_local_M(8, 14) = 0.02024200480788401;
    result_local_M(8, 17) = 0.0326580419402358;
    result_local_M(8, 20) = 0.05462266130342334;
    result_local_M(8, 23) = 0.0330108269172726;
    result_local_M(9, 0) = 0.01113244699701287;
    result_local_M(9, 3) = 0.03667418573650877;
    result_local_M(9, 6) = 0.05113588764893526;
    result_local_M(9, 9) = 0.01565092671515472;
    result_local_M(9, 12) = 0.0147373942105501;
    result_local_M(9, 15) = 0.04914676393919025;
    result_local_M(9, 18) = 0.06874336536380476;
    result_local_M(9, 21) = 0.02072373431883661;
    result_local_M(10, 1) = 0.01113244699701287;
    result_local_M(10, 4) = 0.03667418573650877;
    result_local_M(10, 7) = 0.05113588764893526;
    result_local_M(10, 10) = 0.01565092671515472;
    result_local_M(10, 13) = 0.0147373942105501;
    result_local_M(10, 16) = 0.04914676393919025;
    result_local_M(10, 19) = 0.06874336536380476;
    result_local_M(10, 22) = 0.02072373431883661;
    result_local_M(11, 2) = 0.01113244699701287;
    result_local_M(11, 5) = 0.03667418573650877;
    result_local_M(11, 8) = 0.05113588764893526;
    result_local_M(11, 11) = 0.01565092671515472;
    result_local_M(11, 14) = 0.0147373942105501;
    result_local_M(11, 17) = 0.04914676393919025;
    result_local_M(11, 20) = 0.06874336536380476;
    result_local_M(11, 23) = 0.02072373431883661;

    // Results for Kappa.
    result_local_kappa(0) = 0.06330890473715966;
    result_local_kappa(1) = 0.06330890473715966;
    result_local_kappa(2) = 0.06330890473715966;
    result_local_kappa(3) = 0.07439856220123855;
    result_local_kappa(4) = 0.07439856220123855;
    result_local_kappa(5) = 0.07439856220123855;
    result_local_kappa(6) = 0.2113217850428964;
    result_local_kappa(7) = 0.2113217850428964;
    result_local_kappa(8) = 0.2113217850428964;
    result_local_kappa(9) = 0.2679447049299934;
    result_local_kappa(10) = 0.2679447049299934;
    result_local_kappa(11) = 0.2679447049299934;

    // Perform the unit tests.
    PerformMortarPairUnitTest(contact_pair, q_beam, q_beam_rot, q_solid, result_local_D,
        result_local_M, result_local_kappa);
  }

  /**
   * \brief Test a non straight beam in a hex20 element, with line2 mortar shape functions.
   */
  void test_beam_to_solid_meshtying_mortar_hermite3_hex20_line2()
  {
    // Element types.
    typedef GEOMETRYPAIR::t_hermite beam_type;
    typedef GEOMETRYPAIR::t_hex20 solid_type;
    typedef GEOMETRYPAIR::t_line2 lambda_type;

    // Create the mesh tying mortar pair.
    BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<beam_type, solid_type, lambda_type>
        contact_pair = BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<beam_type, solid_type,
            lambda_type>();

    // Definition of variables for this test case.
    LINALG::Matrix<beam_type::n_dof_, 1, double> q_beam;
    LINALG::Matrix<9, 1, double> q_beam_rot;
    LINALG::Matrix<solid_type::n_dof_, 1, double> q_solid;
    LINALG::SerialDenseMatrix local_D;
    LINALG::SerialDenseMatrix local_M;
    LINALG::SerialDenseVector local_kappa;

    // Matrices for the results.
    LINALG::Matrix<lambda_type::n_dof_, beam_type::n_dof_, double> result_local_D(true);
    LINALG::Matrix<lambda_type::n_dof_, solid_type::n_dof_, double> result_local_M(true);
    LINALG::Matrix<lambda_type::n_dof_, 1, double> result_local_kappa(true);

    // Define the geometry of the two elements.
    q_beam(0) = 0.15;
    q_beam(1) = 0.2;
    q_beam(2) = 0.3;
    q_beam(3) = 0.5773502691896255;
    q_beam(4) = 0.5773502691896258;
    q_beam(5) = 0.577350269189626;
    q_beam(6) = 0.65;
    q_beam(7) = 0.1;
    q_beam(8) = 0.1;
    q_beam(9) = 0.8017837257372733;
    q_beam(10) = -0.5345224838248488;
    q_beam(11) = 0.2672612419124244;

    q_beam_rot(0) = 1.674352746442651;
    q_beam_rot(1) = 0.1425949677148126;
    q_beam_rot(2) = 1.0831163124736984;
    q_beam_rot(3) = 1.4331999091513161;
    q_beam_rot(4) = -0.6560404572957742;
    q_beam_rot(5) = -0.2491376152457331;
    q_beam_rot(6) = 0.0;
    q_beam_rot(7) = 0.0;
    q_beam_rot(8) = 0.0;

    // Positional DOFs of the solid.
    q_solid(0) = -0.954193516126594;
    q_solid(1) = -0.975482672114534;
    q_solid(2) = -1.00311316662815;
    q_solid(3) = 0.923762409575691;
    q_solid(4) = -1.010615727466753;
    q_solid(5) = -1.014160992102648;
    q_solid(6) = 0.905563488894572;
    q_solid(7) = 1.069973754602123;
    q_solid(8) = -0.935488788632622;
    q_solid(9) = -1.047838992508937;
    q_solid(10) = 1.08888017769104;
    q_solid(11) = -1.036911970151991;
    q_solid(12) = -1.089167378981223;
    q_solid(13) = -1.063820485514786;
    q_solid(14) = 1.081000671652452;
    q_solid(15) = 0.972887477248942;
    q_solid(16) = -1.008995039455167;
    q_solid(17) = 0.923459231565508;
    q_solid(18) = 1.089192110164642;
    q_solid(19) = 1.026798931538616;
    q_solid(20) = 0.964761788314679;
    q_solid(21) = -0.941489453299989;
    q_solid(22) = 0.954327975326486;
    q_solid(23) = 0.963168067937965;
    q_solid(24) = 0.0925155012591516;
    q_solid(25) = -0.973383316878771;
    q_solid(26) = -1.065800877068667;
    q_solid(27) = 0.993751379524698;
    q_solid(28) = -0.089078411148891;
    q_solid(29) = -1.080589211553427;
    q_solid(30) = 0.0823626211382376;
    q_solid(31) = 1.033905181036807;
    q_solid(32) = -0.952741283444911;
    q_solid(33) = -1.024714325285187;
    q_solid(34) = -0.0951797166281465;
    q_solid(35) = -0.917694055170081;
    q_solid(36) = -0.94955183268814;
    q_solid(37) = -0.936209915275568;
    q_solid(38) = -0.07531489653912402;
    q_solid(39) = 0.977045420134747;
    q_solid(40) = -0.915831507095275;
    q_solid(41) = -0.03604489795829757;
    q_solid(42) = 1.044996295051512;
    q_solid(43) = 0.912268002554735;
    q_solid(44) = 0.026884640088015;
    q_solid(45) = -1.003774079332873;
    q_solid(46) = 1.016364821141699;
    q_solid(47) = -0.04547583423590318;
    q_solid(48) = -0.07606458957105064;
    q_solid(49) = -0.970767004162723;
    q_solid(50) = 1.09197908262736;
    q_solid(51) = 1.000582069178839;
    q_solid(52) = 0.02495199075528687;
    q_solid(53) = 0.942442567644017;
    q_solid(54) = -0.0840408739597993;
    q_solid(55) = 1.082175787959178;
    q_solid(56) = 1.016887685734179;
    q_solid(57) = -1.02950388437141;
    q_solid(58) = -0.03189318909487643;
    q_solid(59) = 0.962850907596607;

    // Results for D.
    result_local_D(0, 0) = 0.1945989639134667;
    result_local_D(0, 3) = 0.01783046400866768;
    result_local_D(0, 6) = 0.0989060324956216;
    result_local_D(0, 9) = -0.01321892508537963;
    result_local_D(1, 1) = 0.1945989639134667;
    result_local_D(1, 4) = 0.01783046400866768;
    result_local_D(1, 7) = 0.0989060324956216;
    result_local_D(1, 10) = -0.01321892508537963;
    result_local_D(2, 2) = 0.1945989639134667;
    result_local_D(2, 5) = 0.01783046400866768;
    result_local_D(2, 8) = 0.0989060324956216;
    result_local_D(2, 11) = -0.01321892508537963;
    result_local_D(3, 0) = 0.0947117421336144;
    result_local_D(3, 3) = 0.01321892508537963;
    result_local_D(3, 6) = 0.2287572183685851;
    result_local_D(3, 9) = -0.02037627478307643;
    result_local_D(4, 1) = 0.0947117421336144;
    result_local_D(4, 4) = 0.01321892508537963;
    result_local_D(4, 7) = 0.2287572183685851;
    result_local_D(4, 10) = -0.02037627478307643;
    result_local_D(5, 2) = 0.0947117421336144;
    result_local_D(5, 5) = 0.01321892508537963;
    result_local_D(5, 8) = 0.2287572183685851;
    result_local_D(5, 11) = -0.02037627478307643;

    // Results for M.
    result_local_M(0, 0) = -0.0398394162819935;
    result_local_M(0, 3) = -0.05895183035967408;
    result_local_M(0, 6) = -0.07401850868836873;
    result_local_M(0, 9) = -0.05306369423679186;
    result_local_M(0, 12) = -0.05604004397848604;
    result_local_M(0, 15) = -0.07707774046242664;
    result_local_M(0, 18) = -0.0888559670693071;
    result_local_M(0, 21) = -0.07099406180040867;
    result_local_M(0, 24) = 0.03694675213865694;
    result_local_M(0, 27) = 0.06707414635639982;
    result_local_M(0, 30) = 0.05879187720521393;
    result_local_M(0, 33) = 0.03501948001941349;
    result_local_M(0, 36) = 0.03623590539737283;
    result_local_M(0, 39) = 0.06858604782414874;
    result_local_M(0, 42) = 0.1089742104092924;
    result_local_M(0, 45) = 0.05771874418059304;
    result_local_M(0, 48) = 0.06437180493669468;
    result_local_M(0, 51) = 0.1143173854354552;
    result_local_M(0, 54) = 0.1026275611405286;
    result_local_M(0, 57) = 0.06168234424277533;
    result_local_M(1, 1) = -0.0398394162819935;
    result_local_M(1, 4) = -0.05895183035967408;
    result_local_M(1, 7) = -0.07401850868836873;
    result_local_M(1, 10) = -0.05306369423679186;
    result_local_M(1, 13) = -0.05604004397848604;
    result_local_M(1, 16) = -0.07707774046242664;
    result_local_M(1, 19) = -0.0888559670693071;
    result_local_M(1, 22) = -0.07099406180040867;
    result_local_M(1, 25) = 0.03694675213865694;
    result_local_M(1, 28) = 0.06707414635639982;
    result_local_M(1, 31) = 0.05879187720521393;
    result_local_M(1, 34) = 0.03501948001941349;
    result_local_M(1, 37) = 0.03623590539737283;
    result_local_M(1, 40) = 0.06858604782414874;
    result_local_M(1, 43) = 0.1089742104092924;
    result_local_M(1, 46) = 0.05771874418059304;
    result_local_M(1, 49) = 0.06437180493669468;
    result_local_M(1, 52) = 0.1143173854354552;
    result_local_M(1, 55) = 0.1026275611405286;
    result_local_M(1, 58) = 0.06168234424277533;
    result_local_M(2, 2) = -0.0398394162819935;
    result_local_M(2, 5) = -0.05895183035967408;
    result_local_M(2, 8) = -0.07401850868836873;
    result_local_M(2, 11) = -0.05306369423679186;
    result_local_M(2, 14) = -0.05604004397848604;
    result_local_M(2, 17) = -0.07707774046242664;
    result_local_M(2, 20) = -0.0888559670693071;
    result_local_M(2, 23) = -0.07099406180040867;
    result_local_M(2, 26) = 0.03694675213865694;
    result_local_M(2, 29) = 0.06707414635639982;
    result_local_M(2, 32) = 0.05879187720521393;
    result_local_M(2, 35) = 0.03501948001941349;
    result_local_M(2, 38) = 0.03623590539737283;
    result_local_M(2, 41) = 0.06858604782414874;
    result_local_M(2, 44) = 0.1089742104092924;
    result_local_M(2, 47) = 0.05771874418059304;
    result_local_M(2, 50) = 0.06437180493669468;
    result_local_M(2, 53) = 0.1143173854354552;
    result_local_M(2, 56) = 0.1026275611405286;
    result_local_M(2, 59) = 0.06168234424277533;
    result_local_M(3, 0) = -0.04016500977377299;
    result_local_M(3, 3) = -0.07331629813012729;
    result_local_M(3, 6) = -0.0873163666948256;
    result_local_M(3, 9) = -0.05215134283154804;
    result_local_M(3, 12) = -0.05108220555349687;
    result_local_M(3, 15) = -0.0862083946833473;
    result_local_M(3, 18) = -0.0966321010449244;
    result_local_M(3, 21) = -0.06450745866195223;
    result_local_M(3, 24) = 0.04088198027638182;
    result_local_M(3, 27) = 0.0925276468350354;
    result_local_M(3, 30) = 0.06190494185127097;
    result_local_M(3, 33) = 0.03414879279601527;
    result_local_M(3, 36) = 0.03363080222143566;
    result_local_M(3, 39) = 0.0906020222359216;
    result_local_M(3, 42) = 0.1356924685275168;
    result_local_M(3, 45) = 0.05137611239685881;
    result_local_M(3, 48) = 0.05975727027917102;
    result_local_M(3, 51) = 0.132391537087825;
    result_local_M(3, 54) = 0.0912954866970948;
    result_local_M(3, 57) = 0.05063907667166707;
    result_local_M(4, 1) = -0.04016500977377299;
    result_local_M(4, 4) = -0.07331629813012729;
    result_local_M(4, 7) = -0.0873163666948256;
    result_local_M(4, 10) = -0.05215134283154804;
    result_local_M(4, 13) = -0.05108220555349687;
    result_local_M(4, 16) = -0.0862083946833473;
    result_local_M(4, 19) = -0.0966321010449244;
    result_local_M(4, 22) = -0.06450745866195223;
    result_local_M(4, 25) = 0.04088198027638182;
    result_local_M(4, 28) = 0.0925276468350354;
    result_local_M(4, 31) = 0.06190494185127097;
    result_local_M(4, 34) = 0.03414879279601527;
    result_local_M(4, 37) = 0.03363080222143566;
    result_local_M(4, 40) = 0.0906020222359216;
    result_local_M(4, 43) = 0.1356924685275168;
    result_local_M(4, 46) = 0.05137611239685881;
    result_local_M(4, 49) = 0.05975727027917102;
    result_local_M(4, 52) = 0.132391537087825;
    result_local_M(4, 55) = 0.0912954866970948;
    result_local_M(4, 58) = 0.05063907667166707;
    result_local_M(5, 2) = -0.04016500977377299;
    result_local_M(5, 5) = -0.07331629813012729;
    result_local_M(5, 8) = -0.0873163666948256;
    result_local_M(5, 11) = -0.05215134283154804;
    result_local_M(5, 14) = -0.05108220555349687;
    result_local_M(5, 17) = -0.0862083946833473;
    result_local_M(5, 20) = -0.0966321010449244;
    result_local_M(5, 23) = -0.06450745866195223;
    result_local_M(5, 26) = 0.04088198027638182;
    result_local_M(5, 29) = 0.0925276468350354;
    result_local_M(5, 32) = 0.06190494185127097;
    result_local_M(5, 35) = 0.03414879279601527;
    result_local_M(5, 38) = 0.03363080222143566;
    result_local_M(5, 41) = 0.0906020222359216;
    result_local_M(5, 44) = 0.1356924685275168;
    result_local_M(5, 47) = 0.05137611239685881;
    result_local_M(5, 50) = 0.05975727027917102;
    result_local_M(5, 53) = 0.132391537087825;
    result_local_M(5, 56) = 0.0912954866970948;
    result_local_M(5, 59) = 0.05063907667166707;

    // Results for Kappa.
    result_local_kappa(0) = 0.2935049964090884;
    result_local_kappa(1) = 0.2935049964090884;
    result_local_kappa(2) = 0.2935049964090884;
    result_local_kappa(3) = 0.3234689605021996;
    result_local_kappa(4) = 0.3234689605021996;
    result_local_kappa(5) = 0.3234689605021996;

    // Perform the unit tests.
    PerformMortarPairUnitTest(contact_pair, q_beam, q_beam_rot, q_solid, result_local_D,
        result_local_M, result_local_kappa);
  }

  /**
   * \brief Test a non straight beam in a hex20 element, with line3 mortar shape functions.
   */
  void test_beam_to_solid_meshtying_mortar_hermite3_hex20_line3()
  {
    // Element types.
    typedef GEOMETRYPAIR::t_hermite beam_type;
    typedef GEOMETRYPAIR::t_hex20 solid_type;
    typedef GEOMETRYPAIR::t_line3 lambda_type;

    // Create the mesh tying mortar pair.
    BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<beam_type, solid_type, lambda_type>
        contact_pair = BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<beam_type, solid_type,
            lambda_type>();

    // Definition of variables for this test case.
    LINALG::Matrix<beam_type::n_dof_, 1, double> q_beam;
    LINALG::Matrix<9, 1, double> q_beam_rot;
    LINALG::Matrix<solid_type::n_dof_, 1, double> q_solid;
    LINALG::SerialDenseMatrix local_D;
    LINALG::SerialDenseMatrix local_M;
    LINALG::SerialDenseVector local_kappa;

    // Matrices for the results.
    LINALG::Matrix<lambda_type::n_dof_, beam_type::n_dof_, double> result_local_D(true);
    LINALG::Matrix<lambda_type::n_dof_, solid_type::n_dof_, double> result_local_M(true);
    LINALG::Matrix<lambda_type::n_dof_, 1, double> result_local_kappa(true);

    // Define the geometry of the two elements.
    q_beam(0) = 0.15;
    q_beam(1) = 0.2;
    q_beam(2) = 0.3;
    q_beam(3) = 0.5773502691896255;
    q_beam(4) = 0.5773502691896258;
    q_beam(5) = 0.577350269189626;
    q_beam(6) = 0.65;
    q_beam(7) = 0.1;
    q_beam(8) = 0.1;
    q_beam(9) = 0.8017837257372733;
    q_beam(10) = -0.5345224838248488;
    q_beam(11) = 0.2672612419124244;

    q_beam_rot(0) = 1.674352746442651;
    q_beam_rot(1) = 0.1425949677148126;
    q_beam_rot(2) = 1.0831163124736984;
    q_beam_rot(3) = 1.4331999091513161;
    q_beam_rot(4) = -0.6560404572957742;
    q_beam_rot(5) = -0.2491376152457331;
    q_beam_rot(6) = 0.0;
    q_beam_rot(7) = 0.0;
    q_beam_rot(8) = 0.0;

    // Positional DOFs of the solid.
    q_solid(0) = -0.954193516126594;
    q_solid(1) = -0.975482672114534;
    q_solid(2) = -1.00311316662815;
    q_solid(3) = 0.923762409575691;
    q_solid(4) = -1.010615727466753;
    q_solid(5) = -1.014160992102648;
    q_solid(6) = 0.905563488894572;
    q_solid(7) = 1.069973754602123;
    q_solid(8) = -0.935488788632622;
    q_solid(9) = -1.047838992508937;
    q_solid(10) = 1.08888017769104;
    q_solid(11) = -1.036911970151991;
    q_solid(12) = -1.089167378981223;
    q_solid(13) = -1.063820485514786;
    q_solid(14) = 1.081000671652452;
    q_solid(15) = 0.972887477248942;
    q_solid(16) = -1.008995039455167;
    q_solid(17) = 0.923459231565508;
    q_solid(18) = 1.089192110164642;
    q_solid(19) = 1.026798931538616;
    q_solid(20) = 0.964761788314679;
    q_solid(21) = -0.941489453299989;
    q_solid(22) = 0.954327975326486;
    q_solid(23) = 0.963168067937965;
    q_solid(24) = 0.0925155012591516;
    q_solid(25) = -0.973383316878771;
    q_solid(26) = -1.065800877068667;
    q_solid(27) = 0.993751379524698;
    q_solid(28) = -0.089078411148891;
    q_solid(29) = -1.080589211553427;
    q_solid(30) = 0.0823626211382376;
    q_solid(31) = 1.033905181036807;
    q_solid(32) = -0.952741283444911;
    q_solid(33) = -1.024714325285187;
    q_solid(34) = -0.0951797166281465;
    q_solid(35) = -0.917694055170081;
    q_solid(36) = -0.94955183268814;
    q_solid(37) = -0.936209915275568;
    q_solid(38) = -0.07531489653912402;
    q_solid(39) = 0.977045420134747;
    q_solid(40) = -0.915831507095275;
    q_solid(41) = -0.03604489795829757;
    q_solid(42) = 1.044996295051512;
    q_solid(43) = 0.912268002554735;
    q_solid(44) = 0.026884640088015;
    q_solid(45) = -1.003774079332873;
    q_solid(46) = 1.016364821141699;
    q_solid(47) = -0.04547583423590318;
    q_solid(48) = -0.07606458957105064;
    q_solid(49) = -0.970767004162723;
    q_solid(50) = 1.09197908262736;
    q_solid(51) = 1.000582069178839;
    q_solid(52) = 0.02495199075528687;
    q_solid(53) = 0.942442567644017;
    q_solid(54) = -0.0840408739597993;
    q_solid(55) = 1.082175787959178;
    q_solid(56) = 1.016887685734179;
    q_solid(57) = -1.02950388437141;
    q_solid(58) = -0.03189318909487643;
    q_solid(59) = 0.962850907596607;

    // Results for D.
    result_local_D(0, 0) = 0.0936536895927479;
    result_local_D(0, 3) = 0.00502552312033731;
    result_local_D(0, 6) = -0.01315602206049943;
    result_local_D(0, 9) = 0.0004139841970492526;
    result_local_D(1, 1) = 0.0936536895927479;
    result_local_D(1, 4) = 0.00502552312033731;
    result_local_D(1, 7) = -0.01315602206049943;
    result_local_D(1, 10) = 0.0004139841970492526;
    result_local_D(2, 2) = 0.0936536895927479;
    result_local_D(2, 5) = 0.00502552312033731;
    result_local_D(2, 8) = -0.01315602206049943;
    result_local_D(2, 11) = 0.0004139841970492526;
    result_local_D(3, 0) = -0.006233532187104405;
    result_local_D(3, 3) = 0.0004139841970492533;
    result_local_D(3, 6) = 0.1166951638124641;
    result_local_D(3, 9) = -0.006743365500647556;
    result_local_D(4, 1) = -0.006233532187104405;
    result_local_D(4, 4) = 0.0004139841970492533;
    result_local_D(4, 7) = 0.1166951638124641;
    result_local_D(4, 10) = -0.006743365500647556;
    result_local_D(5, 2) = -0.006233532187104405;
    result_local_D(5, 5) = 0.0004139841970492533;
    result_local_D(5, 8) = 0.1166951638124641;
    result_local_D(5, 11) = -0.006743365500647556;
    result_local_D(6, 0) = 0.2018905486414377;
    result_local_D(6, 3) = 0.02560988177666075;
    result_local_D(6, 6) = 0.2241241091122421;
    result_local_D(6, 9) = -0.02726581856485776;
    result_local_D(7, 1) = 0.2018905486414377;
    result_local_D(7, 4) = 0.02560988177666075;
    result_local_D(7, 7) = 0.2241241091122421;
    result_local_D(7, 10) = -0.02726581856485776;
    result_local_D(8, 2) = 0.2018905486414377;
    result_local_D(8, 5) = 0.02560988177666075;
    result_local_D(8, 8) = 0.2241241091122421;
    result_local_D(8, 11) = -0.02726581856485776;

    // Results for M.
    result_local_M(0, 0) = -0.0121767194956414;
    result_local_M(0, 3) = -0.01367035141522692;
    result_local_M(0, 6) = -0.01806871882624224;
    result_local_M(0, 9) = -0.01618899309855407;
    result_local_M(0, 12) = -0.0189602562474041;
    result_local_M(0, 15) = -0.02091748354714573;
    result_local_M(0, 18) = -0.02462226584119422;
    result_local_M(0, 21) = -0.0234972843478389;
    result_local_M(0, 24) = 0.01011337773588912;
    result_local_M(0, 27) = 0.01280854017900686;
    result_local_M(0, 30) = 0.0162577987896318;
    result_local_M(0, 33) = 0.01100993504629004;
    result_local_M(0, 36) = 0.01218380006748016;
    result_local_M(0, 39) = 0.01470101100682746;
    result_local_M(0, 42) = 0.02401358397204773;
    result_local_M(0, 45) = 0.0194060837364106;
    result_local_M(0, 48) = 0.02151009162200237;
    result_local_M(0, 51) = 0.02943145177650484;
    result_local_M(0, 54) = 0.03431275584911837;
    result_local_M(0, 57) = 0.02285131057028671;
    result_local_M(1, 1) = -0.0121767194956414;
    result_local_M(1, 4) = -0.01367035141522692;
    result_local_M(1, 7) = -0.01806871882624224;
    result_local_M(1, 10) = -0.01618899309855407;
    result_local_M(1, 13) = -0.0189602562474041;
    result_local_M(1, 16) = -0.02091748354714573;
    result_local_M(1, 19) = -0.02462226584119422;
    result_local_M(1, 22) = -0.0234972843478389;
    result_local_M(1, 25) = 0.01011337773588912;
    result_local_M(1, 28) = 0.01280854017900686;
    result_local_M(1, 31) = 0.0162577987896318;
    result_local_M(1, 34) = 0.01100993504629004;
    result_local_M(1, 37) = 0.01218380006748016;
    result_local_M(1, 40) = 0.01470101100682746;
    result_local_M(1, 43) = 0.02401358397204773;
    result_local_M(1, 46) = 0.0194060837364106;
    result_local_M(1, 49) = 0.02151009162200237;
    result_local_M(1, 52) = 0.02943145177650484;
    result_local_M(1, 55) = 0.03431275584911837;
    result_local_M(1, 58) = 0.02285131057028671;
    result_local_M(2, 2) = -0.0121767194956414;
    result_local_M(2, 5) = -0.01367035141522692;
    result_local_M(2, 8) = -0.01806871882624224;
    result_local_M(2, 11) = -0.01618899309855407;
    result_local_M(2, 14) = -0.0189602562474041;
    result_local_M(2, 17) = -0.02091748354714573;
    result_local_M(2, 20) = -0.02462226584119422;
    result_local_M(2, 23) = -0.0234972843478389;
    result_local_M(2, 26) = 0.01011337773588912;
    result_local_M(2, 29) = 0.01280854017900686;
    result_local_M(2, 32) = 0.0162577987896318;
    result_local_M(2, 35) = 0.01100993504629004;
    result_local_M(2, 38) = 0.01218380006748016;
    result_local_M(2, 41) = 0.01470101100682746;
    result_local_M(2, 44) = 0.02401358397204773;
    result_local_M(2, 47) = 0.0194060837364106;
    result_local_M(2, 50) = 0.02151009162200237;
    result_local_M(2, 53) = 0.02943145177650484;
    result_local_M(2, 56) = 0.03431275584911837;
    result_local_M(2, 59) = 0.02285131057028671;
    result_local_M(3, 0) = -0.01250231298742088;
    result_local_M(3, 3) = -0.02803481918568012;
    result_local_M(3, 6) = -0.03136657683269906;
    result_local_M(3, 9) = -0.01527664169331026;
    result_local_M(3, 12) = -0.01400241782241495;
    result_local_M(3, 15) = -0.03004813776806636;
    result_local_M(3, 18) = -0.03239839981681154;
    result_local_M(3, 21) = -0.01701068120938247;
    result_local_M(3, 24) = 0.01404860587361401;
    result_local_M(3, 27) = 0.03826204065764244;
    result_local_M(3, 30) = 0.01937086343568884;
    result_local_M(3, 33) = 0.01013924782289183;
    result_local_M(3, 36) = 0.00957869689154299;
    result_local_M(3, 39) = 0.03671698541860033;
    result_local_M(3, 42) = 0.05073184209027208;
    result_local_M(3, 45) = 0.01306345195267638;
    result_local_M(3, 48) = 0.01689555696447874;
    result_local_M(3, 51) = 0.04750560342887469;
    result_local_M(3, 54) = 0.02298068140568459;
    result_local_M(3, 57) = 0.01180804299917844;
    result_local_M(4, 1) = -0.01250231298742088;
    result_local_M(4, 4) = -0.02803481918568012;
    result_local_M(4, 7) = -0.03136657683269906;
    result_local_M(4, 10) = -0.01527664169331026;
    result_local_M(4, 13) = -0.01400241782241495;
    result_local_M(4, 16) = -0.03004813776806636;
    result_local_M(4, 19) = -0.03239839981681154;
    result_local_M(4, 22) = -0.01701068120938247;
    result_local_M(4, 25) = 0.01404860587361401;
    result_local_M(4, 28) = 0.03826204065764244;
    result_local_M(4, 31) = 0.01937086343568884;
    result_local_M(4, 34) = 0.01013924782289183;
    result_local_M(4, 37) = 0.00957869689154299;
    result_local_M(4, 40) = 0.03671698541860033;
    result_local_M(4, 43) = 0.05073184209027208;
    result_local_M(4, 46) = 0.01306345195267638;
    result_local_M(4, 49) = 0.01689555696447874;
    result_local_M(4, 52) = 0.04750560342887469;
    result_local_M(4, 55) = 0.02298068140568459;
    result_local_M(4, 58) = 0.01180804299917844;
    result_local_M(5, 2) = -0.01250231298742088;
    result_local_M(5, 5) = -0.02803481918568012;
    result_local_M(5, 8) = -0.03136657683269906;
    result_local_M(5, 11) = -0.01527664169331026;
    result_local_M(5, 14) = -0.01400241782241495;
    result_local_M(5, 17) = -0.03004813776806636;
    result_local_M(5, 20) = -0.03239839981681154;
    result_local_M(5, 23) = -0.01701068120938247;
    result_local_M(5, 26) = 0.01404860587361401;
    result_local_M(5, 29) = 0.03826204065764244;
    result_local_M(5, 32) = 0.01937086343568884;
    result_local_M(5, 35) = 0.01013924782289183;
    result_local_M(5, 38) = 0.00957869689154299;
    result_local_M(5, 41) = 0.03671698541860033;
    result_local_M(5, 44) = 0.05073184209027208;
    result_local_M(5, 47) = 0.01306345195267638;
    result_local_M(5, 50) = 0.01689555696447874;
    result_local_M(5, 53) = 0.04750560342887469;
    result_local_M(5, 56) = 0.02298068140568459;
    result_local_M(5, 59) = 0.01180804299917844;
    result_local_M(6, 0) = -0.05532539357270422;
    result_local_M(6, 3) = -0.0905629578888943;
    result_local_M(6, 6) = -0.111899579724253;
    result_local_M(6, 9) = -0.07374940227647557;
    result_local_M(6, 12) = -0.07415957546216387;
    result_local_M(6, 15) = -0.1123205138305618;
    result_local_M(6, 18) = -0.1284674024562258;
    result_local_M(6, 21) = -0.0949935549051395;
    result_local_M(6, 24) = 0.05366674880553564;
    result_local_M(6, 27) = 0.1085312123547859;
    result_local_M(6, 30) = 0.0850681568311642;
    result_local_M(6, 33) = 0.04801908994624688;
    result_local_M(6, 36) = 0.04810421065978533;
    result_local_M(6, 39) = 0.1077700736346426;
    result_local_M(6, 42) = 0.1699212528744894;
    result_local_M(6, 45) = 0.07662532088836487;
    result_local_M(6, 48) = 0.0857234266293846;
    result_local_M(6, 51) = 0.1697718673179007;
    result_local_M(6, 54) = 0.1366296105828205;
    result_local_M(6, 57) = 0.07766206734497726;
    result_local_M(7, 1) = -0.05532539357270422;
    result_local_M(7, 4) = -0.0905629578888943;
    result_local_M(7, 7) = -0.111899579724253;
    result_local_M(7, 10) = -0.07374940227647557;
    result_local_M(7, 13) = -0.07415957546216387;
    result_local_M(7, 16) = -0.1123205138305618;
    result_local_M(7, 19) = -0.1284674024562258;
    result_local_M(7, 22) = -0.0949935549051395;
    result_local_M(7, 25) = 0.05366674880553564;
    result_local_M(7, 28) = 0.1085312123547859;
    result_local_M(7, 31) = 0.0850681568311642;
    result_local_M(7, 34) = 0.04801908994624688;
    result_local_M(7, 37) = 0.04810421065978533;
    result_local_M(7, 40) = 0.1077700736346426;
    result_local_M(7, 43) = 0.1699212528744894;
    result_local_M(7, 46) = 0.07662532088836487;
    result_local_M(7, 49) = 0.0857234266293846;
    result_local_M(7, 52) = 0.1697718673179007;
    result_local_M(7, 55) = 0.1366296105828205;
    result_local_M(7, 58) = 0.07766206734497726;
    result_local_M(8, 2) = -0.05532539357270422;
    result_local_M(8, 5) = -0.0905629578888943;
    result_local_M(8, 8) = -0.111899579724253;
    result_local_M(8, 11) = -0.07374940227647557;
    result_local_M(8, 14) = -0.07415957546216387;
    result_local_M(8, 17) = -0.1123205138305618;
    result_local_M(8, 20) = -0.1284674024562258;
    result_local_M(8, 23) = -0.0949935549051395;
    result_local_M(8, 26) = 0.05366674880553564;
    result_local_M(8, 29) = 0.1085312123547859;
    result_local_M(8, 32) = 0.0850681568311642;
    result_local_M(8, 35) = 0.04801908994624688;
    result_local_M(8, 38) = 0.04810421065978533;
    result_local_M(8, 41) = 0.1077700736346426;
    result_local_M(8, 44) = 0.1699212528744894;
    result_local_M(8, 47) = 0.07662532088836487;
    result_local_M(8, 50) = 0.0857234266293846;
    result_local_M(8, 53) = 0.1697718673179007;
    result_local_M(8, 56) = 0.1366296105828205;
    result_local_M(8, 59) = 0.07766206734497726;

    // Results for Kappa.
    result_local_kappa(0) = 0.0804976675322485;
    result_local_kappa(1) = 0.0804976675322485;
    result_local_kappa(2) = 0.0804976675322485;
    result_local_kappa(3) = 0.1104616316253597;
    result_local_kappa(4) = 0.1104616316253597;
    result_local_kappa(5) = 0.1104616316253597;
    result_local_kappa(6) = 0.4260146577536798;
    result_local_kappa(7) = 0.4260146577536798;
    result_local_kappa(8) = 0.4260146577536798;

    // Perform the unit tests.
    PerformMortarPairUnitTest(contact_pair, q_beam, q_beam_rot, q_solid, result_local_D,
        result_local_M, result_local_kappa);
  }

  /**
   * \brief Test a non straight beam in a hex20 element, with line4 mortar shape functions.
   */
  void test_beam_to_solid_meshtying_mortar_hermite3_hex20_line4()
  {
    // Element types.
    typedef GEOMETRYPAIR::t_hermite beam_type;
    typedef GEOMETRYPAIR::t_hex20 solid_type;
    typedef GEOMETRYPAIR::t_line4 lambda_type;

    // Create the mesh tying mortar pair.
    BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<beam_type, solid_type, lambda_type>
        contact_pair = BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<beam_type, solid_type,
            lambda_type>();

    // Definition of variables for this test case.
    LINALG::Matrix<beam_type::n_dof_, 1, double> q_beam;
    LINALG::Matrix<9, 1, double> q_beam_rot;
    LINALG::Matrix<solid_type::n_dof_, 1, double> q_solid;
    LINALG::SerialDenseMatrix local_D;
    LINALG::SerialDenseMatrix local_M;
    LINALG::SerialDenseVector local_kappa;

    // Matrices for the results.
    LINALG::Matrix<lambda_type::n_dof_, beam_type::n_dof_, double> result_local_D(true);
    LINALG::Matrix<lambda_type::n_dof_, solid_type::n_dof_, double> result_local_M(true);
    LINALG::Matrix<lambda_type::n_dof_, 1, double> result_local_kappa(true);

    // Define the geometry of the two elements.
    q_beam(0) = 0.15;
    q_beam(1) = 0.2;
    q_beam(2) = 0.3;
    q_beam(3) = 0.5773502691896255;
    q_beam(4) = 0.5773502691896258;
    q_beam(5) = 0.577350269189626;
    q_beam(6) = 0.65;
    q_beam(7) = 0.1;
    q_beam(8) = 0.1;
    q_beam(9) = 0.8017837257372733;
    q_beam(10) = -0.5345224838248488;
    q_beam(11) = 0.2672612419124244;

    q_beam_rot(0) = 1.674352746442651;
    q_beam_rot(1) = 0.1425949677148126;
    q_beam_rot(2) = 1.0831163124736984;
    q_beam_rot(3) = 1.4331999091513161;
    q_beam_rot(4) = -0.6560404572957742;
    q_beam_rot(5) = -0.2491376152457331;
    q_beam_rot(6) = 0.0;
    q_beam_rot(7) = 0.0;
    q_beam_rot(8) = 0.0;

    // Positional DOFs of the solid.
    q_solid(0) = -0.954193516126594;
    q_solid(1) = -0.975482672114534;
    q_solid(2) = -1.00311316662815;
    q_solid(3) = 0.923762409575691;
    q_solid(4) = -1.010615727466753;
    q_solid(5) = -1.014160992102648;
    q_solid(6) = 0.905563488894572;
    q_solid(7) = 1.069973754602123;
    q_solid(8) = -0.935488788632622;
    q_solid(9) = -1.047838992508937;
    q_solid(10) = 1.08888017769104;
    q_solid(11) = -1.036911970151991;
    q_solid(12) = -1.089167378981223;
    q_solid(13) = -1.063820485514786;
    q_solid(14) = 1.081000671652452;
    q_solid(15) = 0.972887477248942;
    q_solid(16) = -1.008995039455167;
    q_solid(17) = 0.923459231565508;
    q_solid(18) = 1.089192110164642;
    q_solid(19) = 1.026798931538616;
    q_solid(20) = 0.964761788314679;
    q_solid(21) = -0.941489453299989;
    q_solid(22) = 0.954327975326486;
    q_solid(23) = 0.963168067937965;
    q_solid(24) = 0.0925155012591516;
    q_solid(25) = -0.973383316878771;
    q_solid(26) = -1.065800877068667;
    q_solid(27) = 0.993751379524698;
    q_solid(28) = -0.089078411148891;
    q_solid(29) = -1.080589211553427;
    q_solid(30) = 0.0823626211382376;
    q_solid(31) = 1.033905181036807;
    q_solid(32) = -0.952741283444911;
    q_solid(33) = -1.024714325285187;
    q_solid(34) = -0.0951797166281465;
    q_solid(35) = -0.917694055170081;
    q_solid(36) = -0.94955183268814;
    q_solid(37) = -0.936209915275568;
    q_solid(38) = -0.07531489653912402;
    q_solid(39) = 0.977045420134747;
    q_solid(40) = -0.915831507095275;
    q_solid(41) = -0.03604489795829757;
    q_solid(42) = 1.044996295051512;
    q_solid(43) = 0.912268002554735;
    q_solid(44) = 0.026884640088015;
    q_solid(45) = -1.003774079332873;
    q_solid(46) = 1.016364821141699;
    q_solid(47) = -0.04547583423590318;
    q_solid(48) = -0.07606458957105064;
    q_solid(49) = -0.970767004162723;
    q_solid(50) = 1.09197908262736;
    q_solid(51) = 1.000582069178839;
    q_solid(52) = 0.02495199075528687;
    q_solid(53) = 0.942442567644017;
    q_solid(54) = -0.0840408739597993;
    q_solid(55) = 1.082175787959178;
    q_solid(56) = 1.016887685734179;
    q_solid(57) = -1.02950388437141;
    q_solid(58) = -0.03189318909487643;
    q_solid(59) = 0.962850907596607;

    // Results for D.
    result_local_D(0, 0) = 0.05756355583405802;
    result_local_D(0, 3) = 0.001955220263108095;
    result_local_D(0, 6) = 0.005745348903101642;
    result_local_D(0, 9) = -0.0002830518321958758;
    result_local_D(1, 1) = 0.05756355583405802;
    result_local_D(1, 4) = 0.001955220263108095;
    result_local_D(1, 7) = 0.005745348903101642;
    result_local_D(1, 10) = -0.0002830518321958758;
    result_local_D(2, 2) = 0.05756355583405802;
    result_local_D(2, 5) = 0.001955220263108095;
    result_local_D(2, 8) = 0.005745348903101642;
    result_local_D(2, 11) = -0.0002830518321958758;
    result_local_D(3, 0) = 0.004620282991405782;
    result_local_D(3, 3) = 0.0002830518321958757;
    result_local_D(3, 6) = 0.06977827920983277;
    result_local_D(3, 9) = -0.002638102150795207;
    result_local_D(4, 1) = 0.004620282991405782;
    result_local_D(4, 4) = 0.0002830518321958757;
    result_local_D(4, 7) = 0.06977827920983277;
    result_local_D(4, 10) = -0.002638102150795207;
    result_local_D(5, 2) = 0.004620282991405782;
    result_local_D(5, 5) = 0.0002830518321958757;
    result_local_D(5, 8) = 0.06977827920983277;
    result_local_D(5, 11) = -0.002638102150795207;
    result_local_D(6, 0) = 0.1839793570166088;
    result_local_D(6, 3) = 0.01881461423793542;
    result_local_D(6, 6) = 0.02734242802628759;
    result_local_D(6, 9) = -0.00813357387408627;
    result_local_D(7, 1) = 0.1839793570166088;
    result_local_D(7, 4) = 0.01881461423793542;
    result_local_D(7, 7) = 0.02734242802628759;
    result_local_D(7, 10) = -0.00813357387408627;
    result_local_D(8, 2) = 0.1839793570166088;
    result_local_D(8, 5) = 0.01881461423793542;
    result_local_D(8, 8) = 0.02734242802628759;
    result_local_D(8, 11) = -0.00813357387408627;
    result_local_D(9, 0) = 0.04314751020500857;
    result_local_D(9, 3) = 0.00999650276080791;
    result_local_D(9, 6) = 0.2247971947249847;
    result_local_D(9, 9) = -0.0225404720113787;
    result_local_D(10, 1) = 0.04314751020500857;
    result_local_D(10, 4) = 0.00999650276080791;
    result_local_D(10, 7) = 0.2247971947249847;
    result_local_D(10, 10) = -0.0225404720113787;
    result_local_D(11, 2) = 0.04314751020500857;
    result_local_D(11, 5) = 0.00999650276080791;
    result_local_D(11, 8) = 0.2247971947249847;
    result_local_D(11, 11) = -0.0225404720113787;

    // Results for M.
    result_local_M(0, 0) = -0.00928183897191738;
    result_local_M(0, 3) = -0.01203084322786457;
    result_local_M(0, 6) = -0.01495153966074208;
    result_local_M(0, 9) = -0.01195508304052108;
    result_local_M(0, 12) = -0.01363179242985657;
    result_local_M(0, 15) = -0.01666797468851883;
    result_local_M(0, 18) = -0.0191569220943395;
    result_local_M(0, 21) = -0.01659665530850902;
    result_local_M(0, 24) = 0.00817950847800659;
    result_local_M(0, 27) = 0.01277459886257019;
    result_local_M(0, 30) = 0.012396756712724;
    result_local_M(0, 33) = 0.0082557317872665;
    result_local_M(0, 36) = 0.00897441588267774;
    result_local_M(0, 39) = 0.01375025433927722;
    result_local_M(0, 42) = 0.02072255253995805;
    result_local_M(0, 45) = 0.0136619286755524;
    result_local_M(0, 48) = 0.01556178557559656;
    result_local_M(0, 51) = 0.02375569780721402;
    result_local_M(0, 54) = 0.02370405588153632;
    result_local_M(0, 57) = 0.01584426761704911;
    result_local_M(1, 1) = -0.00928183897191738;
    result_local_M(1, 4) = -0.01203084322786457;
    result_local_M(1, 7) = -0.01495153966074208;
    result_local_M(1, 10) = -0.01195508304052108;
    result_local_M(1, 13) = -0.01363179242985657;
    result_local_M(1, 16) = -0.01666797468851883;
    result_local_M(1, 19) = -0.0191569220943395;
    result_local_M(1, 22) = -0.01659665530850902;
    result_local_M(1, 25) = 0.00817950847800659;
    result_local_M(1, 28) = 0.01277459886257019;
    result_local_M(1, 31) = 0.012396756712724;
    result_local_M(1, 34) = 0.0082557317872665;
    result_local_M(1, 37) = 0.00897441588267774;
    result_local_M(1, 40) = 0.01375025433927722;
    result_local_M(1, 43) = 0.02072255253995805;
    result_local_M(1, 46) = 0.0136619286755524;
    result_local_M(1, 49) = 0.01556178557559656;
    result_local_M(1, 52) = 0.02375569780721402;
    result_local_M(1, 55) = 0.02370405588153632;
    result_local_M(1, 58) = 0.01584426761704911;
    result_local_M(2, 2) = -0.00928183897191738;
    result_local_M(2, 5) = -0.01203084322786457;
    result_local_M(2, 8) = -0.01495153966074208;
    result_local_M(2, 11) = -0.01195508304052108;
    result_local_M(2, 14) = -0.01363179242985657;
    result_local_M(2, 17) = -0.01666797468851883;
    result_local_M(2, 20) = -0.0191569220943395;
    result_local_M(2, 23) = -0.01659665530850902;
    result_local_M(2, 26) = 0.00817950847800659;
    result_local_M(2, 29) = 0.01277459886257019;
    result_local_M(2, 32) = 0.012396756712724;
    result_local_M(2, 35) = 0.0082557317872665;
    result_local_M(2, 38) = 0.00897441588267774;
    result_local_M(2, 41) = 0.01375025433927722;
    result_local_M(2, 44) = 0.02072255253995805;
    result_local_M(2, 47) = 0.0136619286755524;
    result_local_M(2, 50) = 0.01556178557559656;
    result_local_M(2, 53) = 0.02375569780721402;
    result_local_M(2, 56) = 0.02370405588153632;
    result_local_M(2, 59) = 0.01584426761704911;
    result_local_M(3, 0) = -0.00848151931455688;
    result_local_M(3, 3) = -0.01835395763693071;
    result_local_M(3, 6) = -0.02049630853266761;
    result_local_M(3, 9) = -0.0102918764667838;
    result_local_M(3, 12) = -0.01006093470719201;
    result_local_M(3, 15) = -0.02025758239787303;
    result_local_M(3, 18) = -0.02180531825663805;
    result_local_M(3, 21) = -0.0120371158855699;
    result_local_M(3, 24) = 0.00927413153080459;
    result_local_M(3, 27) = 0.02472958042973087;
    result_local_M(3, 30) = 0.01259838590870111;
    result_local_M(3, 33) = 0.00689106483863452;
    result_local_M(3, 36) = 0.006775054743872253;
    result_local_M(3, 39) = 0.02419648288182024;
    result_local_M(3, 42) = 0.03278271691305059;
    result_local_M(3, 45) = 0.00922944190248898;
    result_local_M(3, 48) = 0.01212843468221148;
    result_local_M(3, 51) = 0.03195987398342793;
    result_local_M(3, 54) = 0.01651068005041408;
    result_local_M(3, 57) = 0.00910732753429389;
    result_local_M(4, 1) = -0.00848151931455688;
    result_local_M(4, 4) = -0.01835395763693071;
    result_local_M(4, 7) = -0.02049630853266761;
    result_local_M(4, 10) = -0.0102918764667838;
    result_local_M(4, 13) = -0.01006093470719201;
    result_local_M(4, 16) = -0.02025758239787303;
    result_local_M(4, 19) = -0.02180531825663805;
    result_local_M(4, 22) = -0.0120371158855699;
    result_local_M(4, 25) = 0.00927413153080459;
    result_local_M(4, 28) = 0.02472958042973087;
    result_local_M(4, 31) = 0.01259838590870111;
    result_local_M(4, 34) = 0.00689106483863452;
    result_local_M(4, 37) = 0.006775054743872253;
    result_local_M(4, 40) = 0.02419648288182024;
    result_local_M(4, 43) = 0.03278271691305059;
    result_local_M(4, 46) = 0.00922944190248898;
    result_local_M(4, 49) = 0.01212843468221148;
    result_local_M(4, 52) = 0.03195987398342793;
    result_local_M(4, 55) = 0.01651068005041408;
    result_local_M(4, 58) = 0.00910732753429389;
    result_local_M(5, 2) = -0.00848151931455688;
    result_local_M(5, 5) = -0.01835395763693071;
    result_local_M(5, 8) = -0.02049630853266761;
    result_local_M(5, 11) = -0.0102918764667838;
    result_local_M(5, 14) = -0.01006093470719201;
    result_local_M(5, 17) = -0.02025758239787303;
    result_local_M(5, 20) = -0.02180531825663805;
    result_local_M(5, 23) = -0.0120371158855699;
    result_local_M(5, 26) = 0.00927413153080459;
    result_local_M(5, 29) = 0.02472958042973087;
    result_local_M(5, 32) = 0.01259838590870111;
    result_local_M(5, 35) = 0.00689106483863452;
    result_local_M(5, 38) = 0.006775054743872253;
    result_local_M(5, 41) = 0.02419648288182024;
    result_local_M(5, 44) = 0.03278271691305059;
    result_local_M(5, 47) = 0.00922944190248898;
    result_local_M(5, 50) = 0.01212843468221148;
    result_local_M(5, 53) = 0.03195987398342793;
    result_local_M(5, 56) = 0.01651068005041408;
    result_local_M(5, 59) = 0.00910732753429389;
    result_local_M(6, 0) = -0.02943166416093614;
    result_local_M(6, 3) = -0.03887963377042245;
    result_local_M(6, 6) = -0.05131387989309534;
    result_local_M(6, 9) = -0.04035775602777731;
    result_local_M(6, 12) = -0.04379523225095406;
    result_local_M(6, 15) = -0.05486871926234138;
    result_local_M(6, 18) = -0.06457130716164883;
    result_local_M(6, 21) = -0.05632447020741695;
    result_local_M(6, 24) = 0.02592663857572348;
    result_local_M(6, 27) = 0.04080102858235475;
    result_local_M(6, 30) = 0.04348368504240999;
    result_local_M(6, 33) = 0.02626976850691323;
    result_local_M(6, 36) = 0.02766723155182676;
    result_local_M(6, 39) = 0.04326604761564169;
    result_local_M(6, 42) = 0.07359356412420257;
    result_local_M(6, 45) = 0.04596696051571144;
    result_local_M(6, 48) = 0.04999120312523665;
    result_local_M(6, 51) = 0.0806917121520852;
    result_local_M(6, 54) = 0.0830622038713038;
    result_local_M(6, 57) = 0.05014440411407929;
    result_local_M(7, 1) = -0.02943166416093614;
    result_local_M(7, 4) = -0.03887963377042245;
    result_local_M(7, 7) = -0.05131387989309534;
    result_local_M(7, 10) = -0.04035775602777731;
    result_local_M(7, 13) = -0.04379523225095406;
    result_local_M(7, 16) = -0.05486871926234138;
    result_local_M(7, 19) = -0.06457130716164883;
    result_local_M(7, 22) = -0.05632447020741695;
    result_local_M(7, 25) = 0.02592663857572348;
    result_local_M(7, 28) = 0.04080102858235475;
    result_local_M(7, 31) = 0.04348368504240999;
    result_local_M(7, 34) = 0.02626976850691323;
    result_local_M(7, 37) = 0.02766723155182676;
    result_local_M(7, 40) = 0.04326604761564169;
    result_local_M(7, 43) = 0.07359356412420257;
    result_local_M(7, 46) = 0.04596696051571144;
    result_local_M(7, 49) = 0.04999120312523665;
    result_local_M(7, 52) = 0.0806917121520852;
    result_local_M(7, 55) = 0.0830622038713038;
    result_local_M(7, 58) = 0.05014440411407929;
    result_local_M(8, 2) = -0.02943166416093614;
    result_local_M(8, 5) = -0.03887963377042245;
    result_local_M(8, 8) = -0.05131387989309534;
    result_local_M(8, 11) = -0.04035775602777731;
    result_local_M(8, 14) = -0.04379523225095406;
    result_local_M(8, 17) = -0.05486871926234138;
    result_local_M(8, 20) = -0.06457130716164883;
    result_local_M(8, 23) = -0.05632447020741695;
    result_local_M(8, 26) = 0.02592663857572348;
    result_local_M(8, 29) = 0.04080102858235475;
    result_local_M(8, 32) = 0.04348368504240999;
    result_local_M(8, 35) = 0.02626976850691323;
    result_local_M(8, 38) = 0.02766723155182676;
    result_local_M(8, 41) = 0.04326604761564169;
    result_local_M(8, 44) = 0.07359356412420257;
    result_local_M(8, 47) = 0.04596696051571144;
    result_local_M(8, 50) = 0.04999120312523665;
    result_local_M(8, 53) = 0.0806917121520852;
    result_local_M(8, 56) = 0.0830622038713038;
    result_local_M(8, 59) = 0.05014440411407929;
    result_local_M(9, 0) = -0.0328094036083561;
    result_local_M(9, 3) = -0.06300369385458366;
    result_local_M(9, 6) = -0.07457314729668925;
    result_local_M(9, 9) = -0.0426103215332577;
    result_local_M(9, 12) = -0.03963429014398028;
    result_local_M(9, 15) = -0.07149185879704066;
    result_local_M(9, 18) = -0.07995452060160515;
    result_local_M(9, 21) = -0.05054327906086503;
    result_local_M(9, 24) = 0.0344484538305041;
    result_local_M(9, 27) = 0.0812965853167794;
    result_local_M(9, 30) = 0.0522179913926498;
    result_local_M(9, 33) = 0.02775170768261452;
    result_local_M(9, 36) = 0.02645000544043174;
    result_local_M(9, 39) = 0.0779752852233312;
    result_local_M(9, 42) = 0.117567845359598;
    result_local_M(9, 45) = 0.04023652548369903;
    result_local_M(9, 48) = 0.04644765183282101;
    result_local_M(9, 51) = 0.110301638580553;
    result_local_M(9, 54) = 0.07064610803436924;
    result_local_M(9, 57) = 0.03722542164902012;
    result_local_M(10, 1) = -0.0328094036083561;
    result_local_M(10, 4) = -0.06300369385458366;
    result_local_M(10, 7) = -0.07457314729668925;
    result_local_M(10, 10) = -0.0426103215332577;
    result_local_M(10, 13) = -0.03963429014398028;
    result_local_M(10, 16) = -0.07149185879704066;
    result_local_M(10, 19) = -0.07995452060160515;
    result_local_M(10, 22) = -0.05054327906086503;
    result_local_M(10, 25) = 0.0344484538305041;
    result_local_M(10, 28) = 0.0812965853167794;
    result_local_M(10, 31) = 0.0522179913926498;
    result_local_M(10, 34) = 0.02775170768261452;
    result_local_M(10, 37) = 0.02645000544043174;
    result_local_M(10, 40) = 0.0779752852233312;
    result_local_M(10, 43) = 0.117567845359598;
    result_local_M(10, 46) = 0.04023652548369903;
    result_local_M(10, 49) = 0.04644765183282101;
    result_local_M(10, 52) = 0.110301638580553;
    result_local_M(10, 55) = 0.07064610803436924;
    result_local_M(10, 58) = 0.03722542164902012;
    result_local_M(11, 2) = -0.0328094036083561;
    result_local_M(11, 5) = -0.06300369385458366;
    result_local_M(11, 8) = -0.07457314729668925;
    result_local_M(11, 11) = -0.0426103215332577;
    result_local_M(11, 14) = -0.03963429014398028;
    result_local_M(11, 17) = -0.07149185879704066;
    result_local_M(11, 20) = -0.07995452060160515;
    result_local_M(11, 23) = -0.05054327906086503;
    result_local_M(11, 26) = 0.0344484538305041;
    result_local_M(11, 29) = 0.0812965853167794;
    result_local_M(11, 32) = 0.0522179913926498;
    result_local_M(11, 35) = 0.02775170768261452;
    result_local_M(11, 38) = 0.02645000544043174;
    result_local_M(11, 41) = 0.0779752852233312;
    result_local_M(11, 44) = 0.117567845359598;
    result_local_M(11, 47) = 0.04023652548369903;
    result_local_M(11, 50) = 0.04644765183282101;
    result_local_M(11, 53) = 0.110301638580553;
    result_local_M(11, 56) = 0.07064610803436924;
    result_local_M(11, 59) = 0.03722542164902012;

    // Results for Kappa.
    result_local_kappa(0) = 0.06330890473715966;
    result_local_kappa(1) = 0.06330890473715966;
    result_local_kappa(2) = 0.06330890473715966;
    result_local_kappa(3) = 0.07439856220123855;
    result_local_kappa(4) = 0.07439856220123855;
    result_local_kappa(5) = 0.07439856220123855;
    result_local_kappa(6) = 0.2113217850428964;
    result_local_kappa(7) = 0.2113217850428964;
    result_local_kappa(8) = 0.2113217850428964;
    result_local_kappa(9) = 0.2679447049299934;
    result_local_kappa(10) = 0.2679447049299934;
    result_local_kappa(11) = 0.2679447049299934;

    // Perform the unit tests.
    PerformMortarPairUnitTest(contact_pair, q_beam, q_beam_rot, q_solid, result_local_D,
        result_local_M, result_local_kappa);
  }

  /**
   * \brief Test a non straight beam in a hex27 element, with line2 mortar shape functions.
   */
  void test_beam_to_solid_meshtying_mortar_hermite3_hex27_line2()
  {
    // Element types.
    typedef GEOMETRYPAIR::t_hermite beam_type;
    typedef GEOMETRYPAIR::t_hex27 solid_type;
    typedef GEOMETRYPAIR::t_line2 lambda_type;

    // Create the mesh tying mortar pair.
    BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<beam_type, solid_type, lambda_type>
        contact_pair = BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<beam_type, solid_type,
            lambda_type>();

    // Definition of variables for this test case.
    LINALG::Matrix<beam_type::n_dof_, 1, double> q_beam;
    LINALG::Matrix<9, 1, double> q_beam_rot;
    LINALG::Matrix<solid_type::n_dof_, 1, double> q_solid;
    LINALG::SerialDenseMatrix local_D;
    LINALG::SerialDenseMatrix local_M;
    LINALG::SerialDenseVector local_kappa;

    // Matrices for the results.
    LINALG::Matrix<lambda_type::n_dof_, beam_type::n_dof_, double> result_local_D(true);
    LINALG::Matrix<lambda_type::n_dof_, solid_type::n_dof_, double> result_local_M(true);
    LINALG::Matrix<lambda_type::n_dof_, 1, double> result_local_kappa(true);

    // Define the geometry of the two elements.
    q_beam(0) = 0.15;
    q_beam(1) = 0.2;
    q_beam(2) = 0.3;
    q_beam(3) = 0.5773502691896255;
    q_beam(4) = 0.5773502691896258;
    q_beam(5) = 0.577350269189626;
    q_beam(6) = 0.65;
    q_beam(7) = 0.1;
    q_beam(8) = 0.1;
    q_beam(9) = 0.8017837257372733;
    q_beam(10) = -0.5345224838248488;
    q_beam(11) = 0.2672612419124244;

    q_beam_rot(0) = 1.674352746442651;
    q_beam_rot(1) = 0.1425949677148126;
    q_beam_rot(2) = 1.0831163124736984;
    q_beam_rot(3) = 1.4331999091513161;
    q_beam_rot(4) = -0.6560404572957742;
    q_beam_rot(5) = -0.2491376152457331;
    q_beam_rot(6) = 0.0;
    q_beam_rot(7) = 0.0;
    q_beam_rot(8) = 0.0;

    // Positional DOFs of the solid.
    q_solid(0) = -0.954193516126594;
    q_solid(1) = -0.975482672114534;
    q_solid(2) = -1.00311316662815;
    q_solid(3) = 0.923762409575691;
    q_solid(4) = -1.010615727466753;
    q_solid(5) = -1.014160992102648;
    q_solid(6) = 0.905563488894572;
    q_solid(7) = 1.069973754602123;
    q_solid(8) = -0.935488788632622;
    q_solid(9) = -1.047838992508937;
    q_solid(10) = 1.08888017769104;
    q_solid(11) = -1.036911970151991;
    q_solid(12) = -1.089167378981223;
    q_solid(13) = -1.063820485514786;
    q_solid(14) = 1.081000671652452;
    q_solid(15) = 0.972887477248942;
    q_solid(16) = -1.008995039455167;
    q_solid(17) = 0.923459231565508;
    q_solid(18) = 1.089192110164642;
    q_solid(19) = 1.026798931538616;
    q_solid(20) = 0.964761788314679;
    q_solid(21) = -0.941489453299989;
    q_solid(22) = 0.954327975326486;
    q_solid(23) = 0.963168067937965;
    q_solid(24) = 0.0925155012591516;
    q_solid(25) = -0.973383316878771;
    q_solid(26) = -1.065800877068667;
    q_solid(27) = 0.993751379524698;
    q_solid(28) = -0.089078411148891;
    q_solid(29) = -1.080589211553427;
    q_solid(30) = 0.0823626211382376;
    q_solid(31) = 1.033905181036807;
    q_solid(32) = -0.952741283444911;
    q_solid(33) = -1.024714325285187;
    q_solid(34) = -0.0951797166281465;
    q_solid(35) = -0.917694055170081;
    q_solid(36) = -0.94955183268814;
    q_solid(37) = -0.936209915275568;
    q_solid(38) = -0.07531489653912402;
    q_solid(39) = 0.977045420134747;
    q_solid(40) = -0.915831507095275;
    q_solid(41) = -0.03604489795829757;
    q_solid(42) = 1.044996295051512;
    q_solid(43) = 0.912268002554735;
    q_solid(44) = 0.026884640088015;
    q_solid(45) = -1.003774079332873;
    q_solid(46) = 1.016364821141699;
    q_solid(47) = -0.04547583423590318;
    q_solid(48) = -0.07606458957105064;
    q_solid(49) = -0.970767004162723;
    q_solid(50) = 1.09197908262736;
    q_solid(51) = 1.000582069178839;
    q_solid(52) = 0.02495199075528687;
    q_solid(53) = 0.942442567644017;
    q_solid(54) = -0.0840408739597993;
    q_solid(55) = 1.082175787959178;
    q_solid(56) = 1.016887685734179;
    q_solid(57) = -1.02950388437141;
    q_solid(58) = -0.03189318909487643;
    q_solid(59) = 0.962850907596607;
    q_solid(60) = 0.06298633006639215;
    q_solid(61) = -0.06929279965773389;
    q_solid(62) = -0.963099037764495;
    q_solid(63) = 0.005860067798445501;
    q_solid(64) = -1.087116820781742;
    q_solid(65) = 0.01143091251576034;
    q_solid(66) = 0.990119207129434;
    q_solid(67) = -0.07554386166147545;
    q_solid(68) = 0.02213195187583161;
    q_solid(69) = -0.0892155216099889;
    q_solid(70) = 1.051578415322994;
    q_solid(71) = 0.01644372121101889;
    q_solid(72) = -1.084952712553014;
    q_solid(73) = -0.06218066414480983;
    q_solid(74) = 0.03589340707026345;
    q_solid(75) = -0.01441522951114998;
    q_solid(76) = 0.03853468741828525;
    q_solid(77) = 0.974176595435944;
    q_solid(78) = 0.07916144913106371;
    q_solid(79) = 0.0891334591996562;
    q_solid(80) = 0.0835795633628921;

    // Results for D.
    result_local_D(0, 0) = 0.1945989639134667;
    result_local_D(0, 3) = 0.01783046400866768;
    result_local_D(0, 6) = 0.0989060324956216;
    result_local_D(0, 9) = -0.01321892508537963;
    result_local_D(1, 1) = 0.1945989639134667;
    result_local_D(1, 4) = 0.01783046400866768;
    result_local_D(1, 7) = 0.0989060324956216;
    result_local_D(1, 10) = -0.01321892508537963;
    result_local_D(2, 2) = 0.1945989639134667;
    result_local_D(2, 5) = 0.01783046400866768;
    result_local_D(2, 8) = 0.0989060324956216;
    result_local_D(2, 11) = -0.01321892508537963;
    result_local_D(3, 0) = 0.0947117421336144;
    result_local_D(3, 3) = 0.01321892508537963;
    result_local_D(3, 6) = 0.2287572183685851;
    result_local_D(3, 9) = -0.02037627478307643;
    result_local_D(4, 1) = 0.0947117421336144;
    result_local_D(4, 4) = 0.01321892508537963;
    result_local_D(4, 7) = 0.2287572183685851;
    result_local_D(4, 10) = -0.02037627478307643;
    result_local_D(5, 2) = 0.0947117421336144;
    result_local_D(5, 5) = 0.01321892508537963;
    result_local_D(5, 8) = 0.2287572183685851;
    result_local_D(5, 11) = -0.02037627478307643;

    // Results for M.
    result_local_M(0, 0) = -0.0001064905472216395;
    result_local_M(0, 3) = 0.0001818988901377662;
    result_local_M(0, 6) = -0.0002435142253126353;
    result_local_M(0, 9) = 0.0001426123588236535;
    result_local_M(0, 12) = 0.0001596813089286875;
    result_local_M(0, 15) = -0.0002662116270352872;
    result_local_M(0, 18) = 0.0003567584926038978;
    result_local_M(0, 21) = -0.000213915019954504;
    result_local_M(0, 24) = 0.001224882389735379;
    result_local_M(0, 27) = -0.002911870867136624;
    result_local_M(0, 30) = -0.001631169196494466;
    result_local_M(0, 33) = 0.00170209102247987;
    result_local_M(0, 36) = 0.001491761744751473;
    result_local_M(0, 39) = -0.002809045248858824;
    result_local_M(0, 42) = 0.003722814305488443;
    result_local_M(0, 45) = -0.00198547658003252;
    result_local_M(0, 48) = -0.001894328217357991;
    result_local_M(0, 51) = 0.004245530680850893;
    result_local_M(0, 54) = 0.002521616391906724;
    result_local_M(0, 57) = -0.00254927234955889;
    result_local_M(0, 60) = -0.01997070126246551;
    result_local_M(0, 63) = -0.01544314185813551;
    result_local_M(0, 66) = 0.04674763235914034;
    result_local_M(0, 69) = 0.02051886952081389;
    result_local_M(0, 72) = -0.02440849891980948;
    result_local_M(0, 75) = 0.03092656083677221;
    result_local_M(0, 78) = 0.253995922026029;
    result_local_M(1, 1) = -0.0001064905472216395;
    result_local_M(1, 4) = 0.0001818988901377662;
    result_local_M(1, 7) = -0.0002435142253126353;
    result_local_M(1, 10) = 0.0001426123588236535;
    result_local_M(1, 13) = 0.0001596813089286875;
    result_local_M(1, 16) = -0.0002662116270352872;
    result_local_M(1, 19) = 0.0003567584926038978;
    result_local_M(1, 22) = -0.000213915019954504;
    result_local_M(1, 25) = 0.001224882389735379;
    result_local_M(1, 28) = -0.002911870867136624;
    result_local_M(1, 31) = -0.001631169196494466;
    result_local_M(1, 34) = 0.00170209102247987;
    result_local_M(1, 37) = 0.001491761744751473;
    result_local_M(1, 40) = -0.002809045248858824;
    result_local_M(1, 43) = 0.003722814305488443;
    result_local_M(1, 46) = -0.00198547658003252;
    result_local_M(1, 49) = -0.001894328217357991;
    result_local_M(1, 52) = 0.004245530680850893;
    result_local_M(1, 55) = 0.002521616391906724;
    result_local_M(1, 58) = -0.00254927234955889;
    result_local_M(1, 61) = -0.01997070126246551;
    result_local_M(1, 64) = -0.01544314185813551;
    result_local_M(1, 67) = 0.04674763235914034;
    result_local_M(1, 70) = 0.02051886952081389;
    result_local_M(1, 73) = -0.02440849891980948;
    result_local_M(1, 76) = 0.03092656083677221;
    result_local_M(1, 79) = 0.253995922026029;
    result_local_M(2, 2) = -0.0001064905472216395;
    result_local_M(2, 5) = 0.0001818988901377662;
    result_local_M(2, 8) = -0.0002435142253126353;
    result_local_M(2, 11) = 0.0001426123588236535;
    result_local_M(2, 14) = 0.0001596813089286875;
    result_local_M(2, 17) = -0.0002662116270352872;
    result_local_M(2, 20) = 0.0003567584926038978;
    result_local_M(2, 23) = -0.000213915019954504;
    result_local_M(2, 26) = 0.001224882389735379;
    result_local_M(2, 29) = -0.002911870867136624;
    result_local_M(2, 32) = -0.001631169196494466;
    result_local_M(2, 35) = 0.00170209102247987;
    result_local_M(2, 38) = 0.001491761744751473;
    result_local_M(2, 41) = -0.002809045248858824;
    result_local_M(2, 44) = 0.003722814305488443;
    result_local_M(2, 47) = -0.00198547658003252;
    result_local_M(2, 50) = -0.001894328217357991;
    result_local_M(2, 53) = 0.004245530680850893;
    result_local_M(2, 56) = 0.002521616391906724;
    result_local_M(2, 59) = -0.00254927234955889;
    result_local_M(2, 62) = -0.01997070126246551;
    result_local_M(2, 65) = -0.01544314185813551;
    result_local_M(2, 68) = 0.04674763235914034;
    result_local_M(2, 71) = 0.02051886952081389;
    result_local_M(2, 74) = -0.02440849891980948;
    result_local_M(2, 77) = 0.03092656083677221;
    result_local_M(2, 80) = 0.253995922026029;
    result_local_M(3, 0) = -0.0000828596829973999;
    result_local_M(3, 3) = 0.000177830963916105;
    result_local_M(3, 6) = -0.0002328536655118305;
    result_local_M(3, 9) = 0.0001095997732972963;
    result_local_M(3, 12) = 0.0001105089469755819;
    result_local_M(3, 15) = -0.0002292158651404607;
    result_local_M(3, 18) = 0.0003017411934063502;
    result_local_M(3, 21) = -0.0001468230251042916;
    result_local_M(3, 24) = 0.0006997780737962264;
    result_local_M(3, 27) = -0.003138347368404004;
    result_local_M(3, 30) = -0.000931503355297912;
    result_local_M(3, 33) = 0.001400891557583339;
    result_local_M(3, 36) = 0.001907580135065632;
    result_local_M(3, 39) = -0.004782444693177633;
    result_local_M(3, 42) = 0.006077701060949465;
    result_local_M(3, 45) = -0.002456488411834579;
    result_local_M(3, 48) = -0.0009670197260005;
    result_local_M(3, 51) = 0.003961006626067142;
    result_local_M(3, 54) = 0.001291796725646598;
    result_local_M(3, 57) = -0.001835808118443902;
    result_local_M(3, 60) = -0.01152291829627671;
    result_local_M(3, 63) = -0.01404904918470517;
    result_local_M(3, 66) = 0.0950185331035515;
    result_local_M(3, 69) = 0.01826709068527935;
    result_local_M(3, 72) = -0.0359326351641694;
    result_local_M(3, 75) = 0.01569927367886591;
    result_local_M(3, 78) = 0.2547535945348629;
    result_local_M(4, 1) = -0.0000828596829973999;
    result_local_M(4, 4) = 0.000177830963916105;
    result_local_M(4, 7) = -0.0002328536655118305;
    result_local_M(4, 10) = 0.0001095997732972963;
    result_local_M(4, 13) = 0.0001105089469755819;
    result_local_M(4, 16) = -0.0002292158651404607;
    result_local_M(4, 19) = 0.0003017411934063502;
    result_local_M(4, 22) = -0.0001468230251042916;
    result_local_M(4, 25) = 0.0006997780737962264;
    result_local_M(4, 28) = -0.003138347368404004;
    result_local_M(4, 31) = -0.000931503355297912;
    result_local_M(4, 34) = 0.001400891557583339;
    result_local_M(4, 37) = 0.001907580135065632;
    result_local_M(4, 40) = -0.004782444693177633;
    result_local_M(4, 43) = 0.006077701060949465;
    result_local_M(4, 46) = -0.002456488411834579;
    result_local_M(4, 49) = -0.0009670197260005;
    result_local_M(4, 52) = 0.003961006626067142;
    result_local_M(4, 55) = 0.001291796725646598;
    result_local_M(4, 58) = -0.001835808118443902;
    result_local_M(4, 61) = -0.01152291829627671;
    result_local_M(4, 64) = -0.01404904918470517;
    result_local_M(4, 67) = 0.0950185331035515;
    result_local_M(4, 70) = 0.01826709068527935;
    result_local_M(4, 73) = -0.0359326351641694;
    result_local_M(4, 76) = 0.01569927367886591;
    result_local_M(4, 79) = 0.2547535945348629;
    result_local_M(5, 2) = -0.0000828596829973999;
    result_local_M(5, 5) = 0.000177830963916105;
    result_local_M(5, 8) = -0.0002328536655118305;
    result_local_M(5, 11) = 0.0001095997732972963;
    result_local_M(5, 14) = 0.0001105089469755819;
    result_local_M(5, 17) = -0.0002292158651404607;
    result_local_M(5, 20) = 0.0003017411934063502;
    result_local_M(5, 23) = -0.0001468230251042916;
    result_local_M(5, 26) = 0.0006997780737962264;
    result_local_M(5, 29) = -0.003138347368404004;
    result_local_M(5, 32) = -0.000931503355297912;
    result_local_M(5, 35) = 0.001400891557583339;
    result_local_M(5, 38) = 0.001907580135065632;
    result_local_M(5, 41) = -0.004782444693177633;
    result_local_M(5, 44) = 0.006077701060949465;
    result_local_M(5, 47) = -0.002456488411834579;
    result_local_M(5, 50) = -0.0009670197260005;
    result_local_M(5, 53) = 0.003961006626067142;
    result_local_M(5, 56) = 0.001291796725646598;
    result_local_M(5, 59) = -0.001835808118443902;
    result_local_M(5, 62) = -0.01152291829627671;
    result_local_M(5, 65) = -0.01404904918470517;
    result_local_M(5, 68) = 0.0950185331035515;
    result_local_M(5, 71) = 0.01826709068527935;
    result_local_M(5, 74) = -0.0359326351641694;
    result_local_M(5, 77) = 0.01569927367886591;
    result_local_M(5, 80) = 0.2547535945348629;

    // Results for Kappa.
    result_local_kappa(0) = 0.2935049964090884;
    result_local_kappa(1) = 0.2935049964090884;
    result_local_kappa(2) = 0.2935049964090884;
    result_local_kappa(3) = 0.3234689605021996;
    result_local_kappa(4) = 0.3234689605021996;
    result_local_kappa(5) = 0.3234689605021996;

    // Perform the unit tests.
    PerformMortarPairUnitTest(contact_pair, q_beam, q_beam_rot, q_solid, result_local_D,
        result_local_M, result_local_kappa);
  }

  /**
   * \brief Test a non straight beam in a hex27 element, with line3 mortar shape functions.
   */
  void test_beam_to_solid_meshtying_mortar_hermite3_hex27_line3()
  {
    // Element types.
    typedef GEOMETRYPAIR::t_hermite beam_type;
    typedef GEOMETRYPAIR::t_hex27 solid_type;
    typedef GEOMETRYPAIR::t_line3 lambda_type;

    // Create the mesh tying mortar pair.
    BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<beam_type, solid_type, lambda_type>
        contact_pair = BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<beam_type, solid_type,
            lambda_type>();

    // Definition of variables for this test case.
    LINALG::Matrix<beam_type::n_dof_, 1, double> q_beam;
    LINALG::Matrix<9, 1, double> q_beam_rot;
    LINALG::Matrix<solid_type::n_dof_, 1, double> q_solid;
    LINALG::SerialDenseMatrix local_D;
    LINALG::SerialDenseMatrix local_M;
    LINALG::SerialDenseVector local_kappa;

    // Matrices for the results.
    LINALG::Matrix<lambda_type::n_dof_, beam_type::n_dof_, double> result_local_D(true);
    LINALG::Matrix<lambda_type::n_dof_, solid_type::n_dof_, double> result_local_M(true);
    LINALG::Matrix<lambda_type::n_dof_, 1, double> result_local_kappa(true);

    // Define the geometry of the two elements.
    q_beam(0) = 0.15;
    q_beam(1) = 0.2;
    q_beam(2) = 0.3;
    q_beam(3) = 0.5773502691896255;
    q_beam(4) = 0.5773502691896258;
    q_beam(5) = 0.577350269189626;
    q_beam(6) = 0.65;
    q_beam(7) = 0.1;
    q_beam(8) = 0.1;
    q_beam(9) = 0.8017837257372733;
    q_beam(10) = -0.5345224838248488;
    q_beam(11) = 0.2672612419124244;

    q_beam_rot(0) = 1.674352746442651;
    q_beam_rot(1) = 0.1425949677148126;
    q_beam_rot(2) = 1.0831163124736984;
    q_beam_rot(3) = 1.4331999091513161;
    q_beam_rot(4) = -0.6560404572957742;
    q_beam_rot(5) = -0.2491376152457331;
    q_beam_rot(6) = 0.0;
    q_beam_rot(7) = 0.0;
    q_beam_rot(8) = 0.0;

    // Positional DOFs of the solid.
    q_solid(0) = -0.954193516126594;
    q_solid(1) = -0.975482672114534;
    q_solid(2) = -1.00311316662815;
    q_solid(3) = 0.923762409575691;
    q_solid(4) = -1.010615727466753;
    q_solid(5) = -1.014160992102648;
    q_solid(6) = 0.905563488894572;
    q_solid(7) = 1.069973754602123;
    q_solid(8) = -0.935488788632622;
    q_solid(9) = -1.047838992508937;
    q_solid(10) = 1.08888017769104;
    q_solid(11) = -1.036911970151991;
    q_solid(12) = -1.089167378981223;
    q_solid(13) = -1.063820485514786;
    q_solid(14) = 1.081000671652452;
    q_solid(15) = 0.972887477248942;
    q_solid(16) = -1.008995039455167;
    q_solid(17) = 0.923459231565508;
    q_solid(18) = 1.089192110164642;
    q_solid(19) = 1.026798931538616;
    q_solid(20) = 0.964761788314679;
    q_solid(21) = -0.941489453299989;
    q_solid(22) = 0.954327975326486;
    q_solid(23) = 0.963168067937965;
    q_solid(24) = 0.0925155012591516;
    q_solid(25) = -0.973383316878771;
    q_solid(26) = -1.065800877068667;
    q_solid(27) = 0.993751379524698;
    q_solid(28) = -0.089078411148891;
    q_solid(29) = -1.080589211553427;
    q_solid(30) = 0.0823626211382376;
    q_solid(31) = 1.033905181036807;
    q_solid(32) = -0.952741283444911;
    q_solid(33) = -1.024714325285187;
    q_solid(34) = -0.0951797166281465;
    q_solid(35) = -0.917694055170081;
    q_solid(36) = -0.94955183268814;
    q_solid(37) = -0.936209915275568;
    q_solid(38) = -0.07531489653912402;
    q_solid(39) = 0.977045420134747;
    q_solid(40) = -0.915831507095275;
    q_solid(41) = -0.03604489795829757;
    q_solid(42) = 1.044996295051512;
    q_solid(43) = 0.912268002554735;
    q_solid(44) = 0.026884640088015;
    q_solid(45) = -1.003774079332873;
    q_solid(46) = 1.016364821141699;
    q_solid(47) = -0.04547583423590318;
    q_solid(48) = -0.07606458957105064;
    q_solid(49) = -0.970767004162723;
    q_solid(50) = 1.09197908262736;
    q_solid(51) = 1.000582069178839;
    q_solid(52) = 0.02495199075528687;
    q_solid(53) = 0.942442567644017;
    q_solid(54) = -0.0840408739597993;
    q_solid(55) = 1.082175787959178;
    q_solid(56) = 1.016887685734179;
    q_solid(57) = -1.02950388437141;
    q_solid(58) = -0.03189318909487643;
    q_solid(59) = 0.962850907596607;
    q_solid(60) = 0.06298633006639215;
    q_solid(61) = -0.06929279965773389;
    q_solid(62) = -0.963099037764495;
    q_solid(63) = 0.005860067798445501;
    q_solid(64) = -1.087116820781742;
    q_solid(65) = 0.01143091251576034;
    q_solid(66) = 0.990119207129434;
    q_solid(67) = -0.07554386166147545;
    q_solid(68) = 0.02213195187583161;
    q_solid(69) = -0.0892155216099889;
    q_solid(70) = 1.051578415322994;
    q_solid(71) = 0.01644372121101889;
    q_solid(72) = -1.084952712553014;
    q_solid(73) = -0.06218066414480983;
    q_solid(74) = 0.03589340707026345;
    q_solid(75) = -0.01441522951114998;
    q_solid(76) = 0.03853468741828525;
    q_solid(77) = 0.974176595435944;
    q_solid(78) = 0.07916144913106371;
    q_solid(79) = 0.0891334591996562;
    q_solid(80) = 0.0835795633628921;

    // Results for D.
    result_local_D(0, 0) = 0.0936536895927479;
    result_local_D(0, 3) = 0.00502552312033731;
    result_local_D(0, 6) = -0.01315602206049943;
    result_local_D(0, 9) = 0.0004139841970492526;
    result_local_D(1, 1) = 0.0936536895927479;
    result_local_D(1, 4) = 0.00502552312033731;
    result_local_D(1, 7) = -0.01315602206049943;
    result_local_D(1, 10) = 0.0004139841970492526;
    result_local_D(2, 2) = 0.0936536895927479;
    result_local_D(2, 5) = 0.00502552312033731;
    result_local_D(2, 8) = -0.01315602206049943;
    result_local_D(2, 11) = 0.0004139841970492526;
    result_local_D(3, 0) = -0.006233532187104405;
    result_local_D(3, 3) = 0.0004139841970492533;
    result_local_D(3, 6) = 0.1166951638124641;
    result_local_D(3, 9) = -0.006743365500647556;
    result_local_D(4, 1) = -0.006233532187104405;
    result_local_D(4, 4) = 0.0004139841970492533;
    result_local_D(4, 7) = 0.1166951638124641;
    result_local_D(4, 10) = -0.006743365500647556;
    result_local_D(5, 2) = -0.006233532187104405;
    result_local_D(5, 5) = 0.0004139841970492533;
    result_local_D(5, 8) = 0.1166951638124641;
    result_local_D(5, 11) = -0.006743365500647556;
    result_local_D(6, 0) = 0.2018905486414377;
    result_local_D(6, 3) = 0.02560988177666075;
    result_local_D(6, 6) = 0.2241241091122421;
    result_local_D(6, 9) = -0.02726581856485776;
    result_local_D(7, 1) = 0.2018905486414377;
    result_local_D(7, 4) = 0.02560988177666075;
    result_local_D(7, 7) = 0.2241241091122421;
    result_local_D(7, 10) = -0.02726581856485776;
    result_local_D(8, 2) = 0.2018905486414377;
    result_local_D(8, 5) = 0.02560988177666075;
    result_local_D(8, 8) = 0.2241241091122421;
    result_local_D(8, 11) = -0.02726581856485776;

    // Results for M.
    result_local_M(0, 0) = -0.00003101038802622241;
    result_local_M(0, 3) = 0.00003797994702859163;
    result_local_M(0, 6) = -0.00005111745466242084;
    result_local_M(0, 9) = 0.00004127743908547414;
    result_local_M(0, 12) = 0.00005299219411985259;
    result_local_M(0, 15) = -0.00006797519239524894;
    result_local_M(0, 18) = 0.0000909948922169276;
    result_local_M(0, 21) = -0.00007037543447043705;
    result_local_M(0, 24) = 0.0005148107717828944;
    result_local_M(0, 27) = -0.0005944182650412531;
    result_local_M(0, 30) = -0.0006760913267109294;
    result_local_M(0, 33) = 0.0005059940807261272;
    result_local_M(0, 36) = 0.0002062692823087575;
    result_local_M(0, 39) = -0.00007908227992062582;
    result_local_M(0, 42) = 0.000135106039729768;
    result_local_M(0, 45) = -0.0002827509054457556;
    result_local_M(0, 48) = -0.000859769709328038;
    result_local_M(0, 51) = 0.001085894149028442;
    result_local_M(0, 54) = 0.001128255476927984;
    result_local_M(0, 57) = -0.000871376532494545;
    result_local_M(0, 60) = -0.00879817729098251;
    result_local_M(0, 63) = -0.004558398006380783;
    result_local_M(0, 66) = -0.0001847852962478212;
    result_local_M(0, 69) = 0.006024899347670719;
    result_local_M(0, 72) = -0.002966758232136449;
    result_local_M(0, 75) = 0.01472540422980197;
    result_local_M(0, 78) = 0.07603987599606401;
    result_local_M(1, 1) = -0.00003101038802622241;
    result_local_M(1, 4) = 0.00003797994702859163;
    result_local_M(1, 7) = -0.00005111745466242084;
    result_local_M(1, 10) = 0.00004127743908547414;
    result_local_M(1, 13) = 0.00005299219411985259;
    result_local_M(1, 16) = -0.00006797519239524894;
    result_local_M(1, 19) = 0.0000909948922169276;
    result_local_M(1, 22) = -0.00007037543447043705;
    result_local_M(1, 25) = 0.0005148107717828944;
    result_local_M(1, 28) = -0.0005944182650412531;
    result_local_M(1, 31) = -0.0006760913267109294;
    result_local_M(1, 34) = 0.0005059940807261272;
    result_local_M(1, 37) = 0.0002062692823087575;
    result_local_M(1, 40) = -0.00007908227992062582;
    result_local_M(1, 43) = 0.000135106039729768;
    result_local_M(1, 46) = -0.0002827509054457556;
    result_local_M(1, 49) = -0.000859769709328038;
    result_local_M(1, 52) = 0.001085894149028442;
    result_local_M(1, 55) = 0.001128255476927984;
    result_local_M(1, 58) = -0.000871376532494545;
    result_local_M(1, 61) = -0.00879817729098251;
    result_local_M(1, 64) = -0.004558398006380783;
    result_local_M(1, 67) = -0.0001847852962478212;
    result_local_M(1, 70) = 0.006024899347670719;
    result_local_M(1, 73) = -0.002966758232136449;
    result_local_M(1, 76) = 0.01472540422980197;
    result_local_M(1, 79) = 0.07603987599606401;
    result_local_M(2, 2) = -0.00003101038802622241;
    result_local_M(2, 5) = 0.00003797994702859163;
    result_local_M(2, 8) = -0.00005111745466242084;
    result_local_M(2, 11) = 0.00004127743908547414;
    result_local_M(2, 14) = 0.00005299219411985259;
    result_local_M(2, 17) = -0.00006797519239524894;
    result_local_M(2, 20) = 0.0000909948922169276;
    result_local_M(2, 23) = -0.00007037543447043705;
    result_local_M(2, 26) = 0.0005148107717828944;
    result_local_M(2, 29) = -0.0005944182650412531;
    result_local_M(2, 32) = -0.0006760913267109294;
    result_local_M(2, 35) = 0.0005059940807261272;
    result_local_M(2, 38) = 0.0002062692823087575;
    result_local_M(2, 41) = -0.00007908227992062582;
    result_local_M(2, 44) = 0.000135106039729768;
    result_local_M(2, 47) = -0.0002827509054457556;
    result_local_M(2, 50) = -0.000859769709328038;
    result_local_M(2, 53) = 0.001085894149028442;
    result_local_M(2, 56) = 0.001128255476927984;
    result_local_M(2, 59) = -0.000871376532494545;
    result_local_M(2, 62) = -0.00879817729098251;
    result_local_M(2, 65) = -0.004558398006380783;
    result_local_M(2, 68) = -0.0001847852962478212;
    result_local_M(2, 71) = 0.006024899347670719;
    result_local_M(2, 74) = -0.002966758232136449;
    result_local_M(2, 77) = 0.01472540422980197;
    result_local_M(2, 80) = 0.07603987599606401;
    result_local_M(3, 0) = -7.379523801982765e-6;
    result_local_M(3, 3) = 0.00003391202080693039;
    result_local_M(3, 6) = -0.00004045689486161603;
    result_local_M(3, 9) = 8.26485355911692e-6;
    result_local_M(3, 12) = 3.819832166746965e-6;
    result_local_M(3, 15) = -0.00003097943050042247;
    result_local_M(3, 18) = 0.00003597759301937998;
    result_local_M(3, 21) = -3.283439620224666e-6;
    result_local_M(3, 24) = -0.0000102935441562574;
    result_local_M(3, 27) = -0.000820894766308633;
    result_local_M(3, 30) = 0.00002357451448562373;
    result_local_M(3, 33) = 0.000204794615829597;
    result_local_M(3, 36) = 0.0006220876726229164;
    result_local_M(3, 39) = -0.002052481724239436;
    result_local_M(3, 42) = 0.00248999279519079;
    result_local_M(3, 45) = -0.000753762737247814;
    result_local_M(3, 48) = 0.00006753878202945282;
    result_local_M(3, 51) = 0.000801370094244692;
    result_local_M(3, 54) = -0.0001015641893321414;
    result_local_M(3, 57) = -0.0001579123013795568;
    result_local_M(3, 60) = -0.0003503943247937131;
    result_local_M(3, 63) = -0.003164305332950449;
    result_local_M(3, 66) = 0.04808611544816329;
    result_local_M(3, 69) = 0.003773120512136181;
    result_local_M(3, 72) = -0.01449089447649638;
    result_local_M(3, 75) = -0.0005018829281043281;
    result_local_M(3, 78) = 0.07679754850489791;
    result_local_M(4, 1) = -7.379523801982765e-6;
    result_local_M(4, 4) = 0.00003391202080693039;
    result_local_M(4, 7) = -0.00004045689486161603;
    result_local_M(4, 10) = 8.26485355911692e-6;
    result_local_M(4, 13) = 3.81983216674696e-6;
    result_local_M(4, 16) = -0.00003097943050042247;
    result_local_M(4, 19) = 0.00003597759301937998;
    result_local_M(4, 22) = -3.283439620224666e-6;
    result_local_M(4, 25) = -0.0000102935441562574;
    result_local_M(4, 28) = -0.000820894766308633;
    result_local_M(4, 31) = 0.00002357451448562373;
    result_local_M(4, 34) = 0.000204794615829597;
    result_local_M(4, 37) = 0.0006220876726229164;
    result_local_M(4, 40) = -0.002052481724239436;
    result_local_M(4, 43) = 0.00248999279519079;
    result_local_M(4, 46) = -0.000753762737247814;
    result_local_M(4, 49) = 0.00006753878202945282;
    result_local_M(4, 52) = 0.000801370094244692;
    result_local_M(4, 55) = -0.0001015641893321414;
    result_local_M(4, 58) = -0.0001579123013795568;
    result_local_M(4, 61) = -0.0003503943247937131;
    result_local_M(4, 64) = -0.003164305332950449;
    result_local_M(4, 67) = 0.04808611544816329;
    result_local_M(4, 70) = 0.003773120512136181;
    result_local_M(4, 73) = -0.01449089447649638;
    result_local_M(4, 76) = -0.0005018829281043281;
    result_local_M(4, 79) = 0.07679754850489791;
    result_local_M(5, 2) = -7.379523801982765e-6;
    result_local_M(5, 5) = 0.00003391202080693039;
    result_local_M(5, 8) = -0.00004045689486161603;
    result_local_M(5, 11) = 8.26485355911692e-6;
    result_local_M(5, 14) = 3.819832166746965e-6;
    result_local_M(5, 17) = -0.00003097943050042247;
    result_local_M(5, 20) = 0.00003597759301937998;
    result_local_M(5, 23) = -3.283439620224666e-6;
    result_local_M(5, 26) = -0.0000102935441562574;
    result_local_M(5, 29) = -0.000820894766308633;
    result_local_M(5, 32) = 0.00002357451448562373;
    result_local_M(5, 35) = 0.000204794615829597;
    result_local_M(5, 38) = 0.0006220876726229164;
    result_local_M(5, 41) = -0.002052481724239436;
    result_local_M(5, 44) = 0.00248999279519079;
    result_local_M(5, 47) = -0.000753762737247814;
    result_local_M(5, 50) = 0.00006753878202945282;
    result_local_M(5, 53) = 0.000801370094244692;
    result_local_M(5, 56) = -0.0001015641893321414;
    result_local_M(5, 59) = -0.0001579123013795568;
    result_local_M(5, 62) = -0.0003503943247937131;
    result_local_M(5, 65) = -0.003164305332950449;
    result_local_M(5, 68) = 0.04808611544816329;
    result_local_M(5, 71) = 0.003773120512136181;
    result_local_M(5, 74) = -0.01449089447649638;
    result_local_M(5, 77) = -0.0005018829281043281;
    result_local_M(5, 80) = 0.07679754850489791;
    result_local_M(6, 0) = -0.0001509603183908343;
    result_local_M(6, 3) = 0.0002878378862183491;
    result_local_M(6, 6) = -0.000384793541300429;
    result_local_M(6, 9) = 0.0002026698394763587;
    result_local_M(6, 12) = 0.0002133782296176699;
    result_local_M(6, 15) = -0.0003964728692800765;
    result_local_M(6, 18) = 0.0005315272007739404;
    result_local_M(6, 21) = -0.0002870791709681339;
    result_local_M(6, 24) = 0.001420143235904968;
    result_local_M(6, 27) = -0.004634905204190742;
    result_local_M(6, 30) = -0.001910155739567073;
    result_local_M(6, 33) = 0.002392193883507485;
    result_local_M(6, 36) = 0.002570984924885431;
    result_local_M(6, 39) = -0.005459925937876395;
    result_local_M(6, 42) = 0.00717541653151735;
    result_local_M(6, 45) = -0.00340545134917353;
    result_local_M(6, 48) = -0.002069117016059906;
    result_local_M(6, 51) = 0.006319273063644902;
    result_local_M(6, 54) = 0.002786721829957479;
    result_local_M(6, 57) = -0.00335579163412869;
    result_local_M(6, 60) = -0.022345047942966;
    result_local_M(6, 63) = -0.02176948770350945;
    result_local_M(6, 66) = 0.0938648353107763;
    result_local_M(6, 69) = 0.02898794034628635;
    result_local_M(6, 72) = -0.04288348137534605;
    result_local_M(6, 75) = 0.03240231321394049;
    result_local_M(6, 78) = 0.35591209205993;
    result_local_M(7, 1) = -0.0001509603183908343;
    result_local_M(7, 4) = 0.0002878378862183491;
    result_local_M(7, 7) = -0.000384793541300429;
    result_local_M(7, 10) = 0.0002026698394763587;
    result_local_M(7, 13) = 0.0002133782296176699;
    result_local_M(7, 16) = -0.0003964728692800765;
    result_local_M(7, 19) = 0.0005315272007739404;
    result_local_M(7, 22) = -0.0002870791709681339;
    result_local_M(7, 25) = 0.001420143235904968;
    result_local_M(7, 28) = -0.004634905204190742;
    result_local_M(7, 31) = -0.001910155739567073;
    result_local_M(7, 34) = 0.002392193883507485;
    result_local_M(7, 37) = 0.002570984924885431;
    result_local_M(7, 40) = -0.005459925937876395;
    result_local_M(7, 43) = 0.00717541653151735;
    result_local_M(7, 46) = -0.00340545134917353;
    result_local_M(7, 49) = -0.002069117016059906;
    result_local_M(7, 52) = 0.006319273063644902;
    result_local_M(7, 55) = 0.002786721829957479;
    result_local_M(7, 58) = -0.00335579163412869;
    result_local_M(7, 61) = -0.022345047942966;
    result_local_M(7, 64) = -0.02176948770350945;
    result_local_M(7, 67) = 0.0938648353107763;
    result_local_M(7, 70) = 0.02898794034628635;
    result_local_M(7, 73) = -0.04288348137534605;
    result_local_M(7, 76) = 0.03240231321394049;
    result_local_M(7, 79) = 0.35591209205993;
    result_local_M(8, 2) = -0.0001509603183908343;
    result_local_M(8, 5) = 0.0002878378862183491;
    result_local_M(8, 8) = -0.000384793541300429;
    result_local_M(8, 11) = 0.0002026698394763587;
    result_local_M(8, 14) = 0.0002133782296176699;
    result_local_M(8, 17) = -0.0003964728692800765;
    result_local_M(8, 20) = 0.0005315272007739404;
    result_local_M(8, 23) = -0.0002870791709681339;
    result_local_M(8, 26) = 0.001420143235904968;
    result_local_M(8, 29) = -0.004634905204190742;
    result_local_M(8, 32) = -0.001910155739567073;
    result_local_M(8, 35) = 0.002392193883507485;
    result_local_M(8, 38) = 0.002570984924885431;
    result_local_M(8, 41) = -0.005459925937876395;
    result_local_M(8, 44) = 0.00717541653151735;
    result_local_M(8, 47) = -0.00340545134917353;
    result_local_M(8, 50) = -0.002069117016059906;
    result_local_M(8, 53) = 0.006319273063644902;
    result_local_M(8, 56) = 0.002786721829957479;
    result_local_M(8, 59) = -0.00335579163412869;
    result_local_M(8, 62) = -0.022345047942966;
    result_local_M(8, 65) = -0.02176948770350945;
    result_local_M(8, 68) = 0.0938648353107763;
    result_local_M(8, 71) = 0.02898794034628635;
    result_local_M(8, 74) = -0.04288348137534605;
    result_local_M(8, 77) = 0.03240231321394049;
    result_local_M(8, 80) = 0.35591209205993;

    // Results for Kappa.
    result_local_kappa(0) = 0.0804976675322485;
    result_local_kappa(1) = 0.0804976675322485;
    result_local_kappa(2) = 0.0804976675322485;
    result_local_kappa(3) = 0.1104616316253597;
    result_local_kappa(4) = 0.1104616316253597;
    result_local_kappa(5) = 0.1104616316253597;
    result_local_kappa(6) = 0.4260146577536798;
    result_local_kappa(7) = 0.4260146577536798;
    result_local_kappa(8) = 0.4260146577536798;

    // Perform the unit tests.
    PerformMortarPairUnitTest(contact_pair, q_beam, q_beam_rot, q_solid, result_local_D,
        result_local_M, result_local_kappa);
  }

  /**
   * \brief Test a non straight beam in a hex27 element, with line4 mortar shape functions.
   */
  void test_beam_to_solid_meshtying_mortar_hermite3_hex27_line4()
  {
    // Element types.
    typedef GEOMETRYPAIR::t_hermite beam_type;
    typedef GEOMETRYPAIR::t_hex27 solid_type;
    typedef GEOMETRYPAIR::t_line4 lambda_type;

    // Create the mesh tying mortar pair.
    BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<beam_type, solid_type, lambda_type>
        contact_pair = BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<beam_type, solid_type,
            lambda_type>();

    // Definition of variables for this test case.
    LINALG::Matrix<beam_type::n_dof_, 1, double> q_beam;
    LINALG::Matrix<9, 1, double> q_beam_rot;
    LINALG::Matrix<solid_type::n_dof_, 1, double> q_solid;
    LINALG::SerialDenseMatrix local_D;
    LINALG::SerialDenseMatrix local_M;
    LINALG::SerialDenseVector local_kappa;

    // Matrices for the results.
    LINALG::Matrix<lambda_type::n_dof_, beam_type::n_dof_, double> result_local_D(true);
    LINALG::Matrix<lambda_type::n_dof_, solid_type::n_dof_, double> result_local_M(true);
    LINALG::Matrix<lambda_type::n_dof_, 1, double> result_local_kappa(true);

    // Define the geometry of the two elements.
    q_beam(0) = 0.15;
    q_beam(1) = 0.2;
    q_beam(2) = 0.3;
    q_beam(3) = 0.5773502691896255;
    q_beam(4) = 0.5773502691896258;
    q_beam(5) = 0.577350269189626;
    q_beam(6) = 0.65;
    q_beam(7) = 0.1;
    q_beam(8) = 0.1;
    q_beam(9) = 0.8017837257372733;
    q_beam(10) = -0.5345224838248488;
    q_beam(11) = 0.2672612419124244;

    q_beam_rot(0) = 1.674352746442651;
    q_beam_rot(1) = 0.1425949677148126;
    q_beam_rot(2) = 1.0831163124736984;
    q_beam_rot(3) = 1.4331999091513161;
    q_beam_rot(4) = -0.6560404572957742;
    q_beam_rot(5) = -0.2491376152457331;
    q_beam_rot(6) = 0.0;
    q_beam_rot(7) = 0.0;
    q_beam_rot(8) = 0.0;

    // Positional DOFs of the solid.
    q_solid(0) = -0.954193516126594;
    q_solid(1) = -0.975482672114534;
    q_solid(2) = -1.00311316662815;
    q_solid(3) = 0.923762409575691;
    q_solid(4) = -1.010615727466753;
    q_solid(5) = -1.014160992102648;
    q_solid(6) = 0.905563488894572;
    q_solid(7) = 1.069973754602123;
    q_solid(8) = -0.935488788632622;
    q_solid(9) = -1.047838992508937;
    q_solid(10) = 1.08888017769104;
    q_solid(11) = -1.036911970151991;
    q_solid(12) = -1.089167378981223;
    q_solid(13) = -1.063820485514786;
    q_solid(14) = 1.081000671652452;
    q_solid(15) = 0.972887477248942;
    q_solid(16) = -1.008995039455167;
    q_solid(17) = 0.923459231565508;
    q_solid(18) = 1.089192110164642;
    q_solid(19) = 1.026798931538616;
    q_solid(20) = 0.964761788314679;
    q_solid(21) = -0.941489453299989;
    q_solid(22) = 0.954327975326486;
    q_solid(23) = 0.963168067937965;
    q_solid(24) = 0.0925155012591516;
    q_solid(25) = -0.973383316878771;
    q_solid(26) = -1.065800877068667;
    q_solid(27) = 0.993751379524698;
    q_solid(28) = -0.089078411148891;
    q_solid(29) = -1.080589211553427;
    q_solid(30) = 0.0823626211382376;
    q_solid(31) = 1.033905181036807;
    q_solid(32) = -0.952741283444911;
    q_solid(33) = -1.024714325285187;
    q_solid(34) = -0.0951797166281465;
    q_solid(35) = -0.917694055170081;
    q_solid(36) = -0.94955183268814;
    q_solid(37) = -0.936209915275568;
    q_solid(38) = -0.07531489653912402;
    q_solid(39) = 0.977045420134747;
    q_solid(40) = -0.915831507095275;
    q_solid(41) = -0.03604489795829757;
    q_solid(42) = 1.044996295051512;
    q_solid(43) = 0.912268002554735;
    q_solid(44) = 0.026884640088015;
    q_solid(45) = -1.003774079332873;
    q_solid(46) = 1.016364821141699;
    q_solid(47) = -0.04547583423590318;
    q_solid(48) = -0.07606458957105064;
    q_solid(49) = -0.970767004162723;
    q_solid(50) = 1.09197908262736;
    q_solid(51) = 1.000582069178839;
    q_solid(52) = 0.02495199075528687;
    q_solid(53) = 0.942442567644017;
    q_solid(54) = -0.0840408739597993;
    q_solid(55) = 1.082175787959178;
    q_solid(56) = 1.016887685734179;
    q_solid(57) = -1.02950388437141;
    q_solid(58) = -0.03189318909487643;
    q_solid(59) = 0.962850907596607;
    q_solid(60) = 0.06298633006639215;
    q_solid(61) = -0.06929279965773389;
    q_solid(62) = -0.963099037764495;
    q_solid(63) = 0.005860067798445501;
    q_solid(64) = -1.087116820781742;
    q_solid(65) = 0.01143091251576034;
    q_solid(66) = 0.990119207129434;
    q_solid(67) = -0.07554386166147545;
    q_solid(68) = 0.02213195187583161;
    q_solid(69) = -0.0892155216099889;
    q_solid(70) = 1.051578415322994;
    q_solid(71) = 0.01644372121101889;
    q_solid(72) = -1.084952712553014;
    q_solid(73) = -0.06218066414480983;
    q_solid(74) = 0.03589340707026345;
    q_solid(75) = -0.01441522951114998;
    q_solid(76) = 0.03853468741828525;
    q_solid(77) = 0.974176595435944;
    q_solid(78) = 0.07916144913106371;
    q_solid(79) = 0.0891334591996562;
    q_solid(80) = 0.0835795633628921;

    // Results for D.
    result_local_D(0, 0) = 0.05756355583405802;
    result_local_D(0, 3) = 0.001955220263108095;
    result_local_D(0, 6) = 0.005745348903101642;
    result_local_D(0, 9) = -0.0002830518321958758;
    result_local_D(1, 1) = 0.05756355583405802;
    result_local_D(1, 4) = 0.001955220263108095;
    result_local_D(1, 7) = 0.005745348903101642;
    result_local_D(1, 10) = -0.0002830518321958758;
    result_local_D(2, 2) = 0.05756355583405802;
    result_local_D(2, 5) = 0.001955220263108095;
    result_local_D(2, 8) = 0.005745348903101642;
    result_local_D(2, 11) = -0.0002830518321958758;
    result_local_D(3, 0) = 0.004620282991405782;
    result_local_D(3, 3) = 0.0002830518321958757;
    result_local_D(3, 6) = 0.06977827920983277;
    result_local_D(3, 9) = -0.002638102150795207;
    result_local_D(4, 1) = 0.004620282991405782;
    result_local_D(4, 4) = 0.0002830518321958757;
    result_local_D(4, 7) = 0.06977827920983277;
    result_local_D(4, 10) = -0.002638102150795207;
    result_local_D(5, 2) = 0.004620282991405782;
    result_local_D(5, 5) = 0.0002830518321958757;
    result_local_D(5, 8) = 0.06977827920983277;
    result_local_D(5, 11) = -0.002638102150795207;
    result_local_D(6, 0) = 0.1839793570166088;
    result_local_D(6, 3) = 0.01881461423793542;
    result_local_D(6, 6) = 0.02734242802628759;
    result_local_D(6, 9) = -0.00813357387408627;
    result_local_D(7, 1) = 0.1839793570166088;
    result_local_D(7, 4) = 0.01881461423793542;
    result_local_D(7, 7) = 0.02734242802628759;
    result_local_D(7, 10) = -0.00813357387408627;
    result_local_D(8, 2) = 0.1839793570166088;
    result_local_D(8, 5) = 0.01881461423793542;
    result_local_D(8, 8) = 0.02734242802628759;
    result_local_D(8, 11) = -0.00813357387408627;
    result_local_D(9, 0) = 0.04314751020500857;
    result_local_D(9, 3) = 0.00999650276080791;
    result_local_D(9, 6) = 0.2247971947249847;
    result_local_D(9, 9) = -0.0225404720113787;
    result_local_D(10, 1) = 0.04314751020500857;
    result_local_D(10, 4) = 0.00999650276080791;
    result_local_D(10, 7) = 0.2247971947249847;
    result_local_D(10, 10) = -0.0225404720113787;
    result_local_D(11, 2) = 0.04314751020500857;
    result_local_D(11, 5) = 0.00999650276080791;
    result_local_D(11, 8) = 0.2247971947249847;
    result_local_D(11, 11) = -0.0225404720113787;

    // Results for M.
    result_local_M(0, 0) = -0.00001431329535864171;
    result_local_M(0, 3) = 0.0000179413754998695;
    result_local_M(0, 6) = -0.00002265073368330377;
    result_local_M(0, 9) = 0.00001822781255367886;
    result_local_M(0, 12) = 0.00002425427108570378;
    result_local_M(0, 15) = -0.00003023068986189826;
    result_local_M(0, 18) = 0.00003840920214125847;
    result_local_M(0, 21) = -0.0000310009546953883;
    result_local_M(0, 24) = 0.0002846465021201302;
    result_local_M(0, 27) = -0.0003477609635750689;
    result_local_M(0, 30) = -0.0003633966575873966;
    result_local_M(0, 33) = 0.0002691352237516533;
    result_local_M(0, 36) = 0.0001636299692264964;
    result_local_M(0, 39) = -0.0002611458860558529;
    result_local_M(0, 42) = 0.0003215000776352248;
    result_local_M(0, 45) = -0.0002050028300008692;
    result_local_M(0, 48) = -0.0004764843600916477;
    result_local_M(0, 51) = 0.0005731310710650396;
    result_local_M(0, 54) = 0.0006093616633066393;
    result_local_M(0, 57) = -0.0004504231404486054;
    result_local_M(0, 60) = -0.005315858196012471;
    result_local_M(0, 63) = -0.002955192732375306;
    result_local_M(0, 66) = 0.005576254385143274;
    result_local_M(0, 69) = 0.003747954050650781;
    result_local_M(0, 72) = -0.003273754777292538;
    result_local_M(0, 75) = 0.00884779150367707;
    result_local_M(0, 78) = 0.05656388284634183;
    result_local_M(1, 1) = -0.00001431329535864171;
    result_local_M(1, 4) = 0.0000179413754998695;
    result_local_M(1, 7) = -0.00002265073368330377;
    result_local_M(1, 10) = 0.00001822781255367886;
    result_local_M(1, 13) = 0.00002425427108570378;
    result_local_M(1, 16) = -0.00003023068986189826;
    result_local_M(1, 19) = 0.00003840920214125847;
    result_local_M(1, 22) = -0.0000310009546953883;
    result_local_M(1, 25) = 0.0002846465021201302;
    result_local_M(1, 28) = -0.0003477609635750689;
    result_local_M(1, 31) = -0.0003633966575873966;
    result_local_M(1, 34) = 0.0002691352237516533;
    result_local_M(1, 37) = 0.0001636299692264964;
    result_local_M(1, 40) = -0.0002611458860558529;
    result_local_M(1, 43) = 0.0003215000776352248;
    result_local_M(1, 46) = -0.0002050028300008692;
    result_local_M(1, 49) = -0.0004764843600916477;
    result_local_M(1, 52) = 0.0005731310710650396;
    result_local_M(1, 55) = 0.0006093616633066393;
    result_local_M(1, 58) = -0.0004504231404486054;
    result_local_M(1, 61) = -0.005315858196012471;
    result_local_M(1, 64) = -0.002955192732375306;
    result_local_M(1, 67) = 0.005576254385143274;
    result_local_M(1, 70) = 0.003747954050650781;
    result_local_M(1, 73) = -0.003273754777292538;
    result_local_M(1, 76) = 0.00884779150367707;
    result_local_M(1, 79) = 0.05656388284634183;
    result_local_M(2, 2) = -0.00001431329535864171;
    result_local_M(2, 5) = 0.0000179413754998695;
    result_local_M(2, 8) = -0.00002265073368330377;
    result_local_M(2, 11) = 0.00001822781255367886;
    result_local_M(2, 14) = 0.00002425427108570378;
    result_local_M(2, 17) = -0.00003023068986189826;
    result_local_M(2, 20) = 0.00003840920214125847;
    result_local_M(2, 23) = -0.0000310009546953883;
    result_local_M(2, 26) = 0.0002846465021201302;
    result_local_M(2, 29) = -0.0003477609635750689;
    result_local_M(2, 32) = -0.0003633966575873966;
    result_local_M(2, 35) = 0.0002691352237516533;
    result_local_M(2, 38) = 0.0001636299692264964;
    result_local_M(2, 41) = -0.0002611458860558529;
    result_local_M(2, 44) = 0.0003215000776352248;
    result_local_M(2, 47) = -0.0002050028300008692;
    result_local_M(2, 50) = -0.0004764843600916477;
    result_local_M(2, 53) = 0.0005731310710650396;
    result_local_M(2, 56) = 0.0006093616633066393;
    result_local_M(2, 59) = -0.0004504231404486054;
    result_local_M(2, 62) = -0.005315858196012471;
    result_local_M(2, 65) = -0.002955192732375306;
    result_local_M(2, 68) = 0.005576254385143274;
    result_local_M(2, 71) = 0.003747954050650781;
    result_local_M(2, 74) = -0.003273754777292538;
    result_local_M(2, 77) = 0.00884779150367707;
    result_local_M(2, 80) = 0.05656388284634183;
    result_local_M(3, 0) = -5.206576670709177e-6;
    result_local_M(3, 3) = 0.00001797085655835887;
    result_local_M(3, 6) = -0.00002082442317817949;
    result_local_M(3, 9) = 5.980750156367378e-6;
    result_local_M(3, 12) = 5.885476498687035e-6;
    result_local_M(3, 15) = -0.00001916482437376361;
    result_local_M(3, 18) = 0.00002212238299830653;
    result_local_M(3, 21) = -6.773023024256672e-6;
    result_local_M(3, 24) = 0.0000423528210183861;
    result_local_M(3, 27) = -0.0004881889172509752;
    result_local_M(3, 30) = -0.00005035068719202513;
    result_local_M(3, 33) = 0.0001426292373656352;
    result_local_M(3, 36) = 0.0003433538700944985;
    result_local_M(3, 39) = -0.001187927375869659;
    result_local_M(3, 42) = 0.001406671690845665;
    result_local_M(3, 45) = -0.0004058293940460092;
    result_local_M(3, 48) = -0.0000571069401994497;
    result_local_M(3, 51) = 0.000524224039252482;
    result_local_M(3, 54) = 0.00006898939554451826;
    result_local_M(3, 57) = -0.0001593917391594104;
    result_local_M(3, 60) = -0.001039582426893005;
    result_local_M(3, 63) = -0.002046324644017246;
    result_local_M(3, 66) = 0.03059197135292516;
    result_local_M(3, 69) = 0.00242657326587033;
    result_local_M(3, 72) = -0.00882346275942203;
    result_local_M(3, 75) = 0.001325440646278008;
    result_local_M(3, 78) = 0.05178453014712886;
    result_local_M(4, 1) = -5.206576670709177e-6;
    result_local_M(4, 4) = 0.00001797085655835887;
    result_local_M(4, 7) = -0.00002082442317817949;
    result_local_M(4, 10) = 5.980750156367378e-6;
    result_local_M(4, 13) = 5.885476498687035e-6;
    result_local_M(4, 16) = -0.00001916482437376361;
    result_local_M(4, 19) = 0.00002212238299830653;
    result_local_M(4, 22) = -6.773023024256672e-6;
    result_local_M(4, 25) = 0.0000423528210183861;
    result_local_M(4, 28) = -0.0004881889172509752;
    result_local_M(4, 31) = -0.00005035068719202513;
    result_local_M(4, 34) = 0.0001426292373656352;
    result_local_M(4, 37) = 0.0003433538700944985;
    result_local_M(4, 40) = -0.001187927375869659;
    result_local_M(4, 43) = 0.001406671690845665;
    result_local_M(4, 46) = -0.0004058293940460092;
    result_local_M(4, 49) = -0.0000571069401994497;
    result_local_M(4, 52) = 0.000524224039252482;
    result_local_M(4, 55) = 0.00006898939554451826;
    result_local_M(4, 58) = -0.0001593917391594104;
    result_local_M(4, 61) = -0.001039582426893005;
    result_local_M(4, 64) = -0.002046324644017246;
    result_local_M(4, 67) = 0.03059197135292516;
    result_local_M(4, 70) = 0.00242657326587033;
    result_local_M(4, 73) = -0.00882346275942203;
    result_local_M(4, 76) = 0.001325440646278008;
    result_local_M(4, 79) = 0.05178453014712886;
    result_local_M(5, 2) = -5.206576670709177e-6;
    result_local_M(5, 5) = 0.00001797085655835887;
    result_local_M(5, 8) = -0.00002082442317817949;
    result_local_M(5, 11) = 5.980750156367378e-6;
    result_local_M(5, 14) = 5.885476498687035e-6;
    result_local_M(5, 17) = -0.00001916482437376361;
    result_local_M(5, 20) = 0.00002212238299830653;
    result_local_M(5, 23) = -6.773023024256672e-6;
    result_local_M(5, 26) = 0.0000423528210183861;
    result_local_M(5, 29) = -0.0004881889172509752;
    result_local_M(5, 32) = -0.00005035068719202513;
    result_local_M(5, 35) = 0.0001426292373656352;
    result_local_M(5, 38) = 0.0003433538700944985;
    result_local_M(5, 41) = -0.001187927375869659;
    result_local_M(5, 44) = 0.001406671690845665;
    result_local_M(5, 47) = -0.0004058293940460092;
    result_local_M(5, 50) = -0.0000571069401994497;
    result_local_M(5, 53) = 0.000524224039252482;
    result_local_M(5, 56) = 0.00006898939554451826;
    result_local_M(5, 59) = -0.0001593917391594104;
    result_local_M(5, 62) = -0.001039582426893005;
    result_local_M(5, 65) = -0.002046324644017246;
    result_local_M(5, 68) = 0.03059197135292516;
    result_local_M(5, 71) = 0.00242657326587033;
    result_local_M(5, 74) = -0.00882346275942203;
    result_local_M(5, 77) = 0.001325440646278008;
    result_local_M(5, 80) = 0.05178453014712886;
    result_local_M(6, 0) = -0.0001067013973993049;
    result_local_M(6, 3) = 0.0001680549219180473;
    result_local_M(6, 6) = -0.000229697740925012;
    result_local_M(6, 9) = 0.0001451500693990204;
    result_local_M(6, 12) = 0.0001662306052090727;
    result_local_M(6, 15) = -0.0002619108335800808;
    result_local_M(6, 18) = 0.000357079770517235;
    result_local_M(6, 21) = -0.0002257781284381965;
    result_local_M(6, 24) = 0.001223046522452656;
    result_local_M(6, 27) = -0.002478061355970081;
    result_local_M(6, 30) = -0.001654392409708251;
    result_local_M(6, 33) = 0.001607649277238728;
    result_local_M(6, 36) = 0.00109203728607882;
    result_local_M(6, 39) = -0.001501281408297968;
    result_local_M(6, 42) = 0.002131599085602636;
    result_local_M(6, 45) = -0.001510288482274733;
    result_local_M(6, 48) = -0.001925774928731636;
    result_local_M(6, 51) = 0.003908016632757044;
    result_local_M(6, 54) = 0.002601702127098088;
    result_local_M(6, 57) = -0.002521282038936078;
    result_local_M(6, 60) = -0.01882635026352237;
    result_local_M(6, 63) = -0.01297317371083248;
    result_local_M(6, 66) = 0.01791619419736784;
    result_local_M(6, 69) = 0.0177013135209172;
    result_local_M(6, 72) = -0.01516031588028651;
    result_local_M(6, 75) = 0.02978370563360237;
    result_local_M(6, 78) = 0.1918950139716403;
    result_local_M(7, 1) = -0.0001067013973993049;
    result_local_M(7, 4) = 0.0001680549219180473;
    result_local_M(7, 7) = -0.000229697740925012;
    result_local_M(7, 10) = 0.0001451500693990204;
    result_local_M(7, 13) = 0.0001662306052090727;
    result_local_M(7, 16) = -0.0002619108335800808;
    result_local_M(7, 19) = 0.000357079770517235;
    result_local_M(7, 22) = -0.0002257781284381965;
    result_local_M(7, 25) = 0.001223046522452656;
    result_local_M(7, 28) = -0.002478061355970081;
    result_local_M(7, 31) = -0.001654392409708251;
    result_local_M(7, 34) = 0.001607649277238728;
    result_local_M(7, 37) = 0.00109203728607882;
    result_local_M(7, 40) = -0.001501281408297968;
    result_local_M(7, 43) = 0.002131599085602636;
    result_local_M(7, 46) = -0.001510288482274733;
    result_local_M(7, 49) = -0.001925774928731636;
    result_local_M(7, 52) = 0.003908016632757044;
    result_local_M(7, 55) = 0.002601702127098088;
    result_local_M(7, 58) = -0.002521282038936078;
    result_local_M(7, 61) = -0.01882635026352237;
    result_local_M(7, 64) = -0.01297317371083248;
    result_local_M(7, 67) = 0.01791619419736784;
    result_local_M(7, 70) = 0.0177013135209172;
    result_local_M(7, 73) = -0.01516031588028651;
    result_local_M(7, 76) = 0.02978370563360237;
    result_local_M(7, 79) = 0.1918950139716403;
    result_local_M(8, 2) = -0.0001067013973993049;
    result_local_M(8, 5) = 0.0001680549219180473;
    result_local_M(8, 8) = -0.000229697740925012;
    result_local_M(8, 11) = 0.0001451500693990204;
    result_local_M(8, 14) = 0.0001662306052090727;
    result_local_M(8, 17) = -0.0002619108335800808;
    result_local_M(8, 20) = 0.000357079770517235;
    result_local_M(8, 23) = -0.0002257781284381965;
    result_local_M(8, 26) = 0.001223046522452656;
    result_local_M(8, 29) = -0.002478061355970081;
    result_local_M(8, 32) = -0.001654392409708251;
    result_local_M(8, 35) = 0.001607649277238728;
    result_local_M(8, 38) = 0.00109203728607882;
    result_local_M(8, 41) = -0.001501281408297968;
    result_local_M(8, 44) = 0.002131599085602636;
    result_local_M(8, 47) = -0.001510288482274733;
    result_local_M(8, 50) = -0.001925774928731636;
    result_local_M(8, 53) = 0.003908016632757044;
    result_local_M(8, 56) = 0.002601702127098088;
    result_local_M(8, 59) = -0.002521282038936078;
    result_local_M(8, 62) = -0.01882635026352237;
    result_local_M(8, 65) = -0.01297317371083248;
    result_local_M(8, 68) = 0.01791619419736784;
    result_local_M(8, 71) = 0.0177013135209172;
    result_local_M(8, 74) = -0.01516031588028651;
    result_local_M(8, 77) = 0.02978370563360237;
    result_local_M(8, 80) = 0.1918950139716403;
    result_local_M(9, 0) = -0.0000631289607903836;
    result_local_M(9, 3) = 0.0001557627000775956;
    result_local_M(9, 6) = -0.0002031949930379705;
    result_local_M(9, 9) = 0.0000828535000118832;
    result_local_M(9, 12) = 0.000073819903110806;
    result_local_M(9, 15) = -0.0001841211443600053;
    result_local_M(9, 18) = 0.0002408883303534481;
    result_local_M(9, 21) = -0.0000971859389009542;
    result_local_M(9, 24) = 0.0003746146179404325;
    result_local_M(9, 27) = -0.002736206998744504;
    result_local_M(9, 30) = -0.0004945327973047058;
    result_local_M(9, 33) = 0.001083568841707192;
    result_local_M(9, 36) = 0.00180032075441729;
    result_local_M(9, 39) = -0.004641135271812979;
    result_local_M(9, 42) = 0.005940744512354383;
    result_local_M(9, 45) = -0.002320844285545488;
    result_local_M(9, 48) = -0.000401981714335757;
    result_local_M(9, 51) = 0.003201165563843469;
    result_local_M(9, 54) = 0.0005333599316040758;
    result_local_M(9, 57) = -0.001253983549458698;
    result_local_M(9, 60) = -0.006311828672314378;
    result_local_M(9, 63) = -0.01151749995561566;
    result_local_M(9, 66) = 0.0876817455272555;
    result_local_M(9, 69) = 0.01491011936865494;
    result_local_M(9, 72) = -0.0330836006669778;
    result_local_M(9, 75) = 0.006668896732080677;
    result_local_M(9, 78) = 0.2085060895957809;
    result_local_M(10, 1) = -0.0000631289607903836;
    result_local_M(10, 4) = 0.0001557627000775956;
    result_local_M(10, 7) = -0.0002031949930379705;
    result_local_M(10, 10) = 0.0000828535000118832;
    result_local_M(10, 13) = 0.000073819903110806;
    result_local_M(10, 16) = -0.0001841211443600053;
    result_local_M(10, 19) = 0.0002408883303534481;
    result_local_M(10, 22) = -0.0000971859389009542;
    result_local_M(10, 25) = 0.0003746146179404325;
    result_local_M(10, 28) = -0.002736206998744504;
    result_local_M(10, 31) = -0.0004945327973047058;
    result_local_M(10, 34) = 0.001083568841707192;
    result_local_M(10, 37) = 0.00180032075441729;
    result_local_M(10, 40) = -0.004641135271812979;
    result_local_M(10, 43) = 0.005940744512354383;
    result_local_M(10, 46) = -0.002320844285545488;
    result_local_M(10, 49) = -0.000401981714335757;
    result_local_M(10, 52) = 0.003201165563843469;
    result_local_M(10, 55) = 0.0005333599316040758;
    result_local_M(10, 58) = -0.001253983549458698;
    result_local_M(10, 61) = -0.006311828672314378;
    result_local_M(10, 64) = -0.01151749995561566;
    result_local_M(10, 67) = 0.0876817455272555;
    result_local_M(10, 70) = 0.01491011936865494;
    result_local_M(10, 73) = -0.0330836006669778;
    result_local_M(10, 76) = 0.006668896732080677;
    result_local_M(10, 79) = 0.2085060895957809;
    result_local_M(11, 2) = -0.0000631289607903836;
    result_local_M(11, 5) = 0.0001557627000775956;
    result_local_M(11, 8) = -0.0002031949930379705;
    result_local_M(11, 11) = 0.0000828535000118832;
    result_local_M(11, 14) = 0.000073819903110806;
    result_local_M(11, 17) = -0.0001841211443600053;
    result_local_M(11, 20) = 0.0002408883303534481;
    result_local_M(11, 23) = -0.0000971859389009542;
    result_local_M(11, 26) = 0.0003746146179404325;
    result_local_M(11, 29) = -0.002736206998744504;
    result_local_M(11, 32) = -0.0004945327973047058;
    result_local_M(11, 35) = 0.001083568841707192;
    result_local_M(11, 38) = 0.00180032075441729;
    result_local_M(11, 41) = -0.004641135271812979;
    result_local_M(11, 44) = 0.005940744512354383;
    result_local_M(11, 47) = -0.002320844285545488;
    result_local_M(11, 50) = -0.000401981714335757;
    result_local_M(11, 53) = 0.003201165563843469;
    result_local_M(11, 56) = 0.0005333599316040758;
    result_local_M(11, 59) = -0.001253983549458698;
    result_local_M(11, 62) = -0.006311828672314378;
    result_local_M(11, 65) = -0.01151749995561566;
    result_local_M(11, 68) = 0.0876817455272555;
    result_local_M(11, 71) = 0.01491011936865494;
    result_local_M(11, 74) = -0.0330836006669778;
    result_local_M(11, 77) = 0.006668896732080677;
    result_local_M(11, 80) = 0.2085060895957809;

    // Results for Kappa.
    result_local_kappa(0) = 0.06330890473715966;
    result_local_kappa(1) = 0.06330890473715966;
    result_local_kappa(2) = 0.06330890473715966;
    result_local_kappa(3) = 0.07439856220123855;
    result_local_kappa(4) = 0.07439856220123855;
    result_local_kappa(5) = 0.07439856220123855;
    result_local_kappa(6) = 0.2113217850428964;
    result_local_kappa(7) = 0.2113217850428964;
    result_local_kappa(8) = 0.2113217850428964;
    result_local_kappa(9) = 0.2679447049299934;
    result_local_kappa(10) = 0.2679447049299934;
    result_local_kappa(11) = 0.2679447049299934;

    // Perform the unit tests.
    PerformMortarPairUnitTest(contact_pair, q_beam, q_beam_rot, q_solid, result_local_D,
        result_local_M, result_local_kappa);
  }

  /**
   * \brief Test a non straight beam in a tet4 element, with line2 mortar shape functions.
   */
  void test_beam_to_solid_meshtying_mortar_hermite3_tet4_line2()
  {
    // Element types.
    typedef GEOMETRYPAIR::t_hermite beam_type;
    typedef GEOMETRYPAIR::t_tet4 solid_type;
    typedef GEOMETRYPAIR::t_line2 lambda_type;

    // Create the mesh tying mortar pair.
    BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<beam_type, solid_type, lambda_type>
        contact_pair = BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<beam_type, solid_type,
            lambda_type>();

    // Definition of variables for this test case.
    LINALG::Matrix<beam_type::n_dof_, 1, double> q_beam;
    LINALG::Matrix<9, 1, double> q_beam_rot;
    LINALG::Matrix<solid_type::n_dof_, 1, double> q_solid;
    LINALG::SerialDenseMatrix local_D;
    LINALG::SerialDenseMatrix local_M;
    LINALG::SerialDenseVector local_kappa;

    // Matrices for the results.
    LINALG::Matrix<lambda_type::n_dof_, beam_type::n_dof_, double> result_local_D(true);
    LINALG::Matrix<lambda_type::n_dof_, solid_type::n_dof_, double> result_local_M(true);
    LINALG::Matrix<lambda_type::n_dof_, 1, double> result_local_kappa(true);

    // Define the geometry of the two elements.
    q_beam(0) = 0.15;
    q_beam(1) = 0.2;
    q_beam(2) = 0.3;
    q_beam(3) = 0.5773502691896255;
    q_beam(4) = 0.5773502691896258;
    q_beam(5) = 0.577350269189626;
    q_beam(6) = 0.65;
    q_beam(7) = 0.1;
    q_beam(8) = 0.1;
    q_beam(9) = 0.8017837257372733;
    q_beam(10) = -0.5345224838248488;
    q_beam(11) = 0.2672612419124244;

    q_beam_rot(0) = 1.674352746442651;
    q_beam_rot(1) = 0.1425949677148126;
    q_beam_rot(2) = 1.0831163124736984;
    q_beam_rot(3) = 1.4331999091513161;
    q_beam_rot(4) = -0.6560404572957742;
    q_beam_rot(5) = -0.2491376152457331;
    q_beam_rot(6) = 0.0;
    q_beam_rot(7) = 0.0;
    q_beam_rot(8) = 0.0;

    // Positional DOFs of the solid.
    q_solid(0) = 0.04580648387340619;
    q_solid(1) = 0.02451732788546579;
    q_solid(2) = -0.0031131666281497;
    q_solid(3) = 0.923762409575691;
    q_solid(4) = -0.0106157274667526;
    q_solid(5) = -0.0141609921026476;
    q_solid(6) = -0.0944365111054278;
    q_solid(7) = 1.069973754602123;
    q_solid(8) = 0.06451121136737775;
    q_solid(9) = -0.04783899250893692;
    q_solid(10) = 0.0888801776910401;
    q_solid(11) = 0.963088029848009;

    // Results for D.
    result_local_D(0, 0) = 0.1945989639134667;
    result_local_D(0, 3) = 0.01783046400866768;
    result_local_D(0, 6) = 0.0989060324956216;
    result_local_D(0, 9) = -0.01321892508537963;
    result_local_D(1, 1) = 0.1945989639134667;
    result_local_D(1, 4) = 0.01783046400866768;
    result_local_D(1, 7) = 0.0989060324956216;
    result_local_D(1, 10) = -0.01321892508537963;
    result_local_D(2, 2) = 0.1945989639134667;
    result_local_D(2, 5) = 0.01783046400866768;
    result_local_D(2, 8) = 0.0989060324956216;
    result_local_D(2, 11) = -0.01321892508537963;
    result_local_D(3, 0) = 0.0947117421336144;
    result_local_D(3, 3) = 0.01321892508537963;
    result_local_D(3, 6) = 0.2287572183685851;
    result_local_D(3, 9) = -0.02037627478307643;
    result_local_D(4, 1) = 0.0947117421336144;
    result_local_D(4, 4) = 0.01321892508537963;
    result_local_D(4, 7) = 0.2287572183685851;
    result_local_D(4, 10) = -0.02037627478307643;
    result_local_D(5, 2) = 0.0947117421336144;
    result_local_D(5, 5) = 0.01321892508537963;
    result_local_D(5, 8) = 0.2287572183685851;
    result_local_D(5, 11) = -0.02037627478307643;

    // Results for M.
    result_local_M(0, 0) = 0.05443196855983395;
    result_local_M(0, 3) = 0.1077577326062556;
    result_local_M(0, 6) = 0.05535527568629751;
    result_local_M(0, 9) = 0.07596001955670131;
    result_local_M(1, 1) = 0.05443196855983395;
    result_local_M(1, 4) = 0.1077577326062556;
    result_local_M(1, 7) = 0.05535527568629751;
    result_local_M(1, 10) = 0.07596001955670131;
    result_local_M(2, 2) = 0.05443196855983395;
    result_local_M(2, 5) = 0.1077577326062556;
    result_local_M(2, 8) = 0.05535527568629751;
    result_local_M(2, 11) = 0.07596001955670131;
    result_local_M(3, 0) = 0.04322108019762009;
    result_local_M(3, 3) = 0.1729830168839813;
    result_local_M(3, 6) = 0.05257862602525798;
    result_local_M(3, 9) = 0.05468623739534022;
    result_local_M(4, 1) = 0.04322108019762009;
    result_local_M(4, 4) = 0.1729830168839813;
    result_local_M(4, 7) = 0.05257862602525798;
    result_local_M(4, 10) = 0.05468623739534022;
    result_local_M(5, 2) = 0.04322108019762009;
    result_local_M(5, 5) = 0.1729830168839813;
    result_local_M(5, 8) = 0.05257862602525798;
    result_local_M(5, 11) = 0.05468623739534022;

    // Results for Kappa.
    result_local_kappa(0) = 0.2935049964090884;
    result_local_kappa(1) = 0.2935049964090884;
    result_local_kappa(2) = 0.2935049964090884;
    result_local_kappa(3) = 0.3234689605021996;
    result_local_kappa(4) = 0.3234689605021996;
    result_local_kappa(5) = 0.3234689605021996;

    // Perform the unit tests.
    PerformMortarPairUnitTest(contact_pair, q_beam, q_beam_rot, q_solid, result_local_D,
        result_local_M, result_local_kappa);
  }

  /**
   * \brief Test a non straight beam in a tet4 element, with line3 mortar shape functions.
   */
  void test_beam_to_solid_meshtying_mortar_hermite3_tet4_line3()
  {
    // Element types.
    typedef GEOMETRYPAIR::t_hermite beam_type;
    typedef GEOMETRYPAIR::t_tet4 solid_type;
    typedef GEOMETRYPAIR::t_line3 lambda_type;

    // Create the mesh tying mortar pair.
    BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<beam_type, solid_type, lambda_type>
        contact_pair = BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<beam_type, solid_type,
            lambda_type>();

    // Definition of variables for this test case.
    LINALG::Matrix<beam_type::n_dof_, 1, double> q_beam;
    LINALG::Matrix<9, 1, double> q_beam_rot;
    LINALG::Matrix<solid_type::n_dof_, 1, double> q_solid;
    LINALG::SerialDenseMatrix local_D;
    LINALG::SerialDenseMatrix local_M;
    LINALG::SerialDenseVector local_kappa;

    // Matrices for the results.
    LINALG::Matrix<lambda_type::n_dof_, beam_type::n_dof_, double> result_local_D(true);
    LINALG::Matrix<lambda_type::n_dof_, solid_type::n_dof_, double> result_local_M(true);
    LINALG::Matrix<lambda_type::n_dof_, 1, double> result_local_kappa(true);

    // Define the geometry of the two elements.
    q_beam(0) = 0.15;
    q_beam(1) = 0.2;
    q_beam(2) = 0.3;
    q_beam(3) = 0.5773502691896255;
    q_beam(4) = 0.5773502691896258;
    q_beam(5) = 0.577350269189626;
    q_beam(6) = 0.65;
    q_beam(7) = 0.1;
    q_beam(8) = 0.1;
    q_beam(9) = 0.8017837257372733;
    q_beam(10) = -0.5345224838248488;
    q_beam(11) = 0.2672612419124244;

    q_beam_rot(0) = 1.674352746442651;
    q_beam_rot(1) = 0.1425949677148126;
    q_beam_rot(2) = 1.0831163124736984;
    q_beam_rot(3) = 1.4331999091513161;
    q_beam_rot(4) = -0.6560404572957742;
    q_beam_rot(5) = -0.2491376152457331;
    q_beam_rot(6) = 0.0;
    q_beam_rot(7) = 0.0;
    q_beam_rot(8) = 0.0;

    // Positional DOFs of the solid.
    q_solid(0) = 0.04580648387340619;
    q_solid(1) = 0.02451732788546579;
    q_solid(2) = -0.0031131666281497;
    q_solid(3) = 0.923762409575691;
    q_solid(4) = -0.0106157274667526;
    q_solid(5) = -0.0141609921026476;
    q_solid(6) = -0.0944365111054278;
    q_solid(7) = 1.069973754602123;
    q_solid(8) = 0.06451121136737775;
    q_solid(9) = -0.04783899250893692;
    q_solid(10) = 0.0888801776910401;
    q_solid(11) = 0.963088029848009;

    // Results for D.
    result_local_D(0, 0) = 0.0936536895927479;
    result_local_D(0, 3) = 0.00502552312033731;
    result_local_D(0, 6) = -0.01315602206049943;
    result_local_D(0, 9) = 0.0004139841970492526;
    result_local_D(1, 1) = 0.0936536895927479;
    result_local_D(1, 4) = 0.00502552312033731;
    result_local_D(1, 7) = -0.01315602206049943;
    result_local_D(1, 10) = 0.0004139841970492526;
    result_local_D(2, 2) = 0.0936536895927479;
    result_local_D(2, 5) = 0.00502552312033731;
    result_local_D(2, 8) = -0.01315602206049943;
    result_local_D(2, 11) = 0.0004139841970492526;
    result_local_D(3, 0) = -0.006233532187104405;
    result_local_D(3, 3) = 0.0004139841970492533;
    result_local_D(3, 6) = 0.1166951638124641;
    result_local_D(3, 9) = -0.006743365500647556;
    result_local_D(4, 1) = -0.006233532187104405;
    result_local_D(4, 4) = 0.0004139841970492533;
    result_local_D(4, 7) = 0.1166951638124641;
    result_local_D(4, 10) = -0.006743365500647556;
    result_local_D(5, 2) = -0.006233532187104405;
    result_local_D(5, 5) = 0.0004139841970492533;
    result_local_D(5, 8) = 0.1166951638124641;
    result_local_D(5, 11) = -0.006743365500647556;
    result_local_D(6, 0) = 0.2018905486414377;
    result_local_D(6, 3) = 0.02560988177666075;
    result_local_D(6, 6) = 0.2241241091122421;
    result_local_D(6, 9) = -0.02726581856485776;
    result_local_D(7, 1) = 0.2018905486414377;
    result_local_D(7, 4) = 0.02560988177666075;
    result_local_D(7, 7) = 0.2241241091122421;
    result_local_D(7, 10) = -0.02726581856485776;
    result_local_D(8, 2) = 0.2018905486414377;
    result_local_D(8, 5) = 0.02560988177666075;
    result_local_D(8, 8) = 0.2241241091122421;
    result_local_D(8, 11) = -0.02726581856485776;

    // Results for M.
    result_local_M(0, 0) = 0.02302633185141577;
    result_local_M(0, 3) = 0.01149028891192587;
    result_local_M(0, 6) = 0.01586593713817261;
    result_local_M(0, 9) = 0.03011510963073423;
    result_local_M(1, 1) = 0.02302633185141577;
    result_local_M(1, 4) = 0.01149028891192587;
    result_local_M(1, 7) = 0.01586593713817261;
    result_local_M(1, 10) = 0.03011510963073423;
    result_local_M(2, 2) = 0.02302633185141577;
    result_local_M(2, 5) = 0.01149028891192587;
    result_local_M(2, 8) = 0.01586593713817261;
    result_local_M(2, 11) = 0.03011510963073423;
    result_local_M(3, 0) = 0.01181544348920192;
    result_local_M(3, 3) = 0.07671557318965152;
    result_local_M(3, 6) = 0.01308928747713308;
    result_local_M(3, 9) = 0.00884132746937315;
    result_local_M(4, 1) = 0.01181544348920192;
    result_local_M(4, 4) = 0.07671557318965152;
    result_local_M(4, 7) = 0.01308928747713308;
    result_local_M(4, 10) = 0.00884132746937315;
    result_local_M(5, 2) = 0.01181544348920192;
    result_local_M(5, 5) = 0.07671557318965152;
    result_local_M(5, 8) = 0.01308928747713308;
    result_local_M(5, 11) = 0.00884132746937315;
    result_local_M(6, 0) = 0.06281127341683633;
    result_local_M(6, 3) = 0.1925348873886595;
    result_local_M(6, 6) = 0.07897867709624979;
    result_local_M(6, 9) = 0.0916898198519342;
    result_local_M(7, 1) = 0.06281127341683633;
    result_local_M(7, 4) = 0.1925348873886595;
    result_local_M(7, 7) = 0.07897867709624979;
    result_local_M(7, 10) = 0.0916898198519342;
    result_local_M(8, 2) = 0.06281127341683633;
    result_local_M(8, 5) = 0.1925348873886595;
    result_local_M(8, 8) = 0.07897867709624979;
    result_local_M(8, 11) = 0.0916898198519342;

    // Results for Kappa.
    result_local_kappa(0) = 0.0804976675322485;
    result_local_kappa(1) = 0.0804976675322485;
    result_local_kappa(2) = 0.0804976675322485;
    result_local_kappa(3) = 0.1104616316253597;
    result_local_kappa(4) = 0.1104616316253597;
    result_local_kappa(5) = 0.1104616316253597;
    result_local_kappa(6) = 0.4260146577536798;
    result_local_kappa(7) = 0.4260146577536798;
    result_local_kappa(8) = 0.4260146577536798;

    // Perform the unit tests.
    PerformMortarPairUnitTest(contact_pair, q_beam, q_beam_rot, q_solid, result_local_D,
        result_local_M, result_local_kappa);
  }

  /**
   * \brief Test a non straight beam in a tet4 element, with line4 mortar shape functions.
   */
  void test_beam_to_solid_meshtying_mortar_hermite3_tet4_line4()
  {
    // Element types.
    typedef GEOMETRYPAIR::t_hermite beam_type;
    typedef GEOMETRYPAIR::t_tet4 solid_type;
    typedef GEOMETRYPAIR::t_line4 lambda_type;

    // Create the mesh tying mortar pair.
    BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<beam_type, solid_type, lambda_type>
        contact_pair = BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<beam_type, solid_type,
            lambda_type>();

    // Definition of variables for this test case.
    LINALG::Matrix<beam_type::n_dof_, 1, double> q_beam;
    LINALG::Matrix<9, 1, double> q_beam_rot;
    LINALG::Matrix<solid_type::n_dof_, 1, double> q_solid;
    LINALG::SerialDenseMatrix local_D;
    LINALG::SerialDenseMatrix local_M;
    LINALG::SerialDenseVector local_kappa;

    // Matrices for the results.
    LINALG::Matrix<lambda_type::n_dof_, beam_type::n_dof_, double> result_local_D(true);
    LINALG::Matrix<lambda_type::n_dof_, solid_type::n_dof_, double> result_local_M(true);
    LINALG::Matrix<lambda_type::n_dof_, 1, double> result_local_kappa(true);

    // Define the geometry of the two elements.
    q_beam(0) = 0.15;
    q_beam(1) = 0.2;
    q_beam(2) = 0.3;
    q_beam(3) = 0.5773502691896255;
    q_beam(4) = 0.5773502691896258;
    q_beam(5) = 0.577350269189626;
    q_beam(6) = 0.65;
    q_beam(7) = 0.1;
    q_beam(8) = 0.1;
    q_beam(9) = 0.8017837257372733;
    q_beam(10) = -0.5345224838248488;
    q_beam(11) = 0.2672612419124244;

    q_beam_rot(0) = 1.674352746442651;
    q_beam_rot(1) = 0.1425949677148126;
    q_beam_rot(2) = 1.0831163124736984;
    q_beam_rot(3) = 1.4331999091513161;
    q_beam_rot(4) = -0.6560404572957742;
    q_beam_rot(5) = -0.2491376152457331;
    q_beam_rot(6) = 0.0;
    q_beam_rot(7) = 0.0;
    q_beam_rot(8) = 0.0;

    // Positional DOFs of the solid.
    q_solid(0) = 0.04580648387340619;
    q_solid(1) = 0.02451732788546579;
    q_solid(2) = -0.0031131666281497;
    q_solid(3) = 0.923762409575691;
    q_solid(4) = -0.0106157274667526;
    q_solid(5) = -0.0141609921026476;
    q_solid(6) = -0.0944365111054278;
    q_solid(7) = 1.069973754602123;
    q_solid(8) = 0.06451121136737775;
    q_solid(9) = -0.04783899250893692;
    q_solid(10) = 0.0888801776910401;
    q_solid(11) = 0.963088029848009;

    // Results for D.
    result_local_D(0, 0) = 0.05756355583405802;
    result_local_D(0, 3) = 0.001955220263108095;
    result_local_D(0, 6) = 0.005745348903101642;
    result_local_D(0, 9) = -0.0002830518321958758;
    result_local_D(1, 1) = 0.05756355583405802;
    result_local_D(1, 4) = 0.001955220263108095;
    result_local_D(1, 7) = 0.005745348903101642;
    result_local_D(1, 10) = -0.0002830518321958758;
    result_local_D(2, 2) = 0.05756355583405802;
    result_local_D(2, 5) = 0.001955220263108095;
    result_local_D(2, 8) = 0.005745348903101642;
    result_local_D(2, 11) = -0.0002830518321958758;
    result_local_D(3, 0) = 0.004620282991405782;
    result_local_D(3, 3) = 0.0002830518321958757;
    result_local_D(3, 6) = 0.06977827920983277;
    result_local_D(3, 9) = -0.002638102150795207;
    result_local_D(4, 1) = 0.004620282991405782;
    result_local_D(4, 4) = 0.0002830518321958757;
    result_local_D(4, 7) = 0.06977827920983277;
    result_local_D(4, 10) = -0.002638102150795207;
    result_local_D(5, 2) = 0.004620282991405782;
    result_local_D(5, 5) = 0.0002830518321958757;
    result_local_D(5, 8) = 0.06977827920983277;
    result_local_D(5, 11) = -0.002638102150795207;
    result_local_D(6, 0) = 0.1839793570166088;
    result_local_D(6, 3) = 0.01881461423793542;
    result_local_D(6, 6) = 0.02734242802628759;
    result_local_D(6, 9) = -0.00813357387408627;
    result_local_D(7, 1) = 0.1839793570166088;
    result_local_D(7, 4) = 0.01881461423793542;
    result_local_D(7, 7) = 0.02734242802628759;
    result_local_D(7, 10) = -0.00813357387408627;
    result_local_D(8, 2) = 0.1839793570166088;
    result_local_D(8, 5) = 0.01881461423793542;
    result_local_D(8, 8) = 0.02734242802628759;
    result_local_D(8, 11) = -0.00813357387408627;
    result_local_D(9, 0) = 0.04314751020500857;
    result_local_D(9, 3) = 0.00999650276080791;
    result_local_D(9, 6) = 0.2247971947249847;
    result_local_D(9, 9) = -0.0225404720113787;
    result_local_D(10, 1) = 0.04314751020500857;
    result_local_D(10, 4) = 0.00999650276080791;
    result_local_D(10, 7) = 0.2247971947249847;
    result_local_D(10, 10) = -0.0225404720113787;
    result_local_D(11, 2) = 0.04314751020500857;
    result_local_D(11, 5) = 0.00999650276080791;
    result_local_D(11, 8) = 0.2247971947249847;
    result_local_D(11, 11) = -0.0225404720113787;

    // Results for M.
    result_local_M(0, 0) = 0.0179113322450487;
    result_local_M(0, 3) = 0.01556003555258317;
    result_local_M(0, 6) = 0.0106426831473678;
    result_local_M(0, 9) = 0.01919485379215998;
    result_local_M(1, 1) = 0.0179113322450487;
    result_local_M(1, 4) = 0.01556003555258317;
    result_local_M(1, 7) = 0.0106426831473678;
    result_local_M(1, 10) = 0.01919485379215998;
    result_local_M(2, 2) = 0.0179113322450487;
    result_local_M(2, 5) = 0.01556003555258317;
    result_local_M(2, 8) = 0.0106426831473678;
    result_local_M(2, 11) = 0.01919485379215998;
    result_local_M(3, 0) = 0.00907903391846446;
    result_local_M(3, 3) = 0.0485789656254118;
    result_local_M(3, 6) = 0.00844020683090666;
    result_local_M(3, 9) = 0.00830035582645563;
    result_local_M(4, 1) = 0.00907903391846446;
    result_local_M(4, 4) = 0.0485789656254118;
    result_local_M(4, 7) = 0.00844020683090666;
    result_local_M(4, 10) = 0.00830035582645563;
    result_local_M(5, 2) = 0.00907903391846446;
    result_local_M(5, 5) = 0.0485789656254118;
    result_local_M(5, 8) = 0.00844020683090666;
    result_local_M(5, 11) = 0.00830035582645563;
    result_local_M(6, 0) = 0.03889922635041484;
    result_local_M(6, 3) = 0.05999134284877538;
    result_local_M(6, 6) = 0.04528676588350811;
    result_local_M(6, 9) = 0.06714444996019806;
    result_local_M(7, 1) = 0.03889922635041484;
    result_local_M(7, 4) = 0.05999134284877538;
    result_local_M(7, 7) = 0.04528676588350811;
    result_local_M(7, 10) = 0.06714444996019806;
    result_local_M(8, 2) = 0.03889922635041484;
    result_local_M(8, 5) = 0.05999134284877538;
    result_local_M(8, 8) = 0.04528676588350811;
    result_local_M(8, 11) = 0.06714444996019806;
    result_local_M(9, 0) = 0.03176345624352602;
    result_local_M(9, 3) = 0.1566104054634665;
    result_local_M(9, 6) = 0.04356424584977292;
    result_local_M(9, 9) = 0.03600659737322787;
    result_local_M(10, 1) = 0.03176345624352602;
    result_local_M(10, 4) = 0.1566104054634665;
    result_local_M(10, 7) = 0.04356424584977292;
    result_local_M(10, 10) = 0.03600659737322787;
    result_local_M(11, 2) = 0.03176345624352602;
    result_local_M(11, 5) = 0.1566104054634665;
    result_local_M(11, 8) = 0.04356424584977292;
    result_local_M(11, 11) = 0.03600659737322787;

    // Results for Kappa.
    result_local_kappa(0) = 0.06330890473715966;
    result_local_kappa(1) = 0.06330890473715966;
    result_local_kappa(2) = 0.06330890473715966;
    result_local_kappa(3) = 0.07439856220123855;
    result_local_kappa(4) = 0.07439856220123855;
    result_local_kappa(5) = 0.07439856220123855;
    result_local_kappa(6) = 0.2113217850428964;
    result_local_kappa(7) = 0.2113217850428964;
    result_local_kappa(8) = 0.2113217850428964;
    result_local_kappa(9) = 0.2679447049299934;
    result_local_kappa(10) = 0.2679447049299934;
    result_local_kappa(11) = 0.2679447049299934;

    // Perform the unit tests.
    PerformMortarPairUnitTest(contact_pair, q_beam, q_beam_rot, q_solid, result_local_D,
        result_local_M, result_local_kappa);
  }

  /**
   * \brief Test a non straight beam in a tet10 element, with line2 mortar shape functions.
   */
  void test_beam_to_solid_meshtying_mortar_hermite3_tet10_line2()
  {
    // Element types.
    typedef GEOMETRYPAIR::t_hermite beam_type;
    typedef GEOMETRYPAIR::t_tet10 solid_type;
    typedef GEOMETRYPAIR::t_line2 lambda_type;

    // Create the mesh tying mortar pair.
    BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<beam_type, solid_type, lambda_type>
        contact_pair = BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<beam_type, solid_type,
            lambda_type>();

    // Definition of variables for this test case.
    LINALG::Matrix<beam_type::n_dof_, 1, double> q_beam;
    LINALG::Matrix<9, 1, double> q_beam_rot;
    LINALG::Matrix<solid_type::n_dof_, 1, double> q_solid;
    LINALG::SerialDenseMatrix local_D;
    LINALG::SerialDenseMatrix local_M;
    LINALG::SerialDenseVector local_kappa;

    // Matrices for the results.
    LINALG::Matrix<lambda_type::n_dof_, beam_type::n_dof_, double> result_local_D(true);
    LINALG::Matrix<lambda_type::n_dof_, solid_type::n_dof_, double> result_local_M(true);
    LINALG::Matrix<lambda_type::n_dof_, 1, double> result_local_kappa(true);

    // Define the geometry of the two elements.
    q_beam(0) = 0.15;
    q_beam(1) = 0.2;
    q_beam(2) = 0.3;
    q_beam(3) = 0.5773502691896255;
    q_beam(4) = 0.5773502691896258;
    q_beam(5) = 0.577350269189626;
    q_beam(6) = 0.65;
    q_beam(7) = 0.1;
    q_beam(8) = 0.1;
    q_beam(9) = 0.8017837257372733;
    q_beam(10) = -0.5345224838248488;
    q_beam(11) = 0.2672612419124244;

    q_beam_rot(0) = 1.674352746442651;
    q_beam_rot(1) = 0.1425949677148126;
    q_beam_rot(2) = 1.0831163124736984;
    q_beam_rot(3) = 1.4331999091513161;
    q_beam_rot(4) = -0.6560404572957742;
    q_beam_rot(5) = -0.2491376152457331;
    q_beam_rot(6) = 0.0;
    q_beam_rot(7) = 0.0;
    q_beam_rot(8) = 0.0;

    // Positional DOFs of the solid.
    q_solid(0) = 0.04580648387340619;
    q_solid(1) = 0.02451732788546579;
    q_solid(2) = -0.0031131666281497;
    q_solid(3) = 0.923762409575691;
    q_solid(4) = -0.0106157274667526;
    q_solid(5) = -0.0141609921026476;
    q_solid(6) = -0.0944365111054278;
    q_solid(7) = 1.069973754602123;
    q_solid(8) = 0.06451121136737775;
    q_solid(9) = -0.04783899250893692;
    q_solid(10) = 0.0888801776910401;
    q_solid(11) = 0.963088029848009;
    q_solid(12) = 0.4108326210187769;
    q_solid(13) = -0.06382048551478623;
    q_solid(14) = 0.0810006716524518;
    q_solid(15) = 0.4728874772489418;
    q_solid(16) = 0.4910049605448334;
    q_solid(17) = -0.07654076843449248;
    q_solid(18) = 0.0891921101646426;
    q_solid(19) = 0.5267989315386162;
    q_solid(20) = -0.03523821168532098;
    q_solid(21) = 0.05851054670001138;
    q_solid(22) = -0.04567202467351366;
    q_solid(23) = 0.4631680679379652;
    q_solid(24) = 0.5925155012591516;
    q_solid(25) = 0.02661668312122867;
    q_solid(26) = 0.4341991229313329;
    q_solid(27) = -0.006248620475302225;
    q_solid(28) = 0.4109215888511089;
    q_solid(29) = 0.4194107884465731;

    // Results for D.
    result_local_D(0, 0) = 0.1945989639134667;
    result_local_D(0, 3) = 0.01783046400866768;
    result_local_D(0, 6) = 0.0989060324956216;
    result_local_D(0, 9) = -0.01321892508537963;
    result_local_D(1, 1) = 0.1945989639134667;
    result_local_D(1, 4) = 0.01783046400866768;
    result_local_D(1, 7) = 0.0989060324956216;
    result_local_D(1, 10) = -0.01321892508537963;
    result_local_D(2, 2) = 0.1945989639134667;
    result_local_D(2, 5) = 0.01783046400866768;
    result_local_D(2, 8) = 0.0989060324956216;
    result_local_D(2, 11) = -0.01321892508537963;
    result_local_D(3, 0) = 0.0947117421336144;
    result_local_D(3, 3) = 0.01321892508537963;
    result_local_D(3, 6) = 0.2287572183685851;
    result_local_D(3, 9) = -0.02037627478307643;
    result_local_D(4, 1) = 0.0947117421336144;
    result_local_D(4, 4) = 0.01321892508537963;
    result_local_D(4, 7) = 0.2287572183685851;
    result_local_D(4, 10) = -0.02037627478307643;
    result_local_D(5, 2) = 0.0947117421336144;
    result_local_D(5, 5) = 0.01321892508537963;
    result_local_D(5, 8) = 0.2287572183685851;
    result_local_D(5, 11) = -0.02037627478307643;

    // Results for M.
    result_local_M(0, 0) = -0.02421769608243993;
    result_local_M(0, 3) = -0.02674852460567383;
    result_local_M(0, 6) = -0.0354468725543303;
    result_local_M(0, 9) = -0.02924557885452574;
    result_local_M(0, 12) = 0.03094871943896918;
    result_local_M(0, 15) = 0.0850084257522522;
    result_local_M(0, 18) = 0.03825407006883004;
    result_local_M(0, 21) = 0.0469304311474056;
    result_local_M(0, 24) = 0.0943181456996416;
    result_local_M(0, 27) = 0.1137038763989595;
    result_local_M(1, 1) = -0.02421769608243993;
    result_local_M(1, 4) = -0.02674852460567383;
    result_local_M(1, 7) = -0.0354468725543303;
    result_local_M(1, 10) = -0.02924557885452574;
    result_local_M(1, 13) = 0.03094871943896918;
    result_local_M(1, 16) = 0.0850084257522522;
    result_local_M(1, 19) = 0.03825407006883004;
    result_local_M(1, 22) = 0.0469304311474056;
    result_local_M(1, 25) = 0.0943181456996416;
    result_local_M(1, 28) = 0.1137038763989595;
    result_local_M(2, 2) = -0.02421769608243993;
    result_local_M(2, 5) = -0.02674852460567383;
    result_local_M(2, 8) = -0.0354468725543303;
    result_local_M(2, 11) = -0.02924557885452574;
    result_local_M(2, 14) = 0.03094871943896918;
    result_local_M(2, 17) = 0.0850084257522522;
    result_local_M(2, 20) = 0.03825407006883004;
    result_local_M(2, 23) = 0.0469304311474056;
    result_local_M(2, 26) = 0.0943181456996416;
    result_local_M(2, 29) = 0.1137038763989595;
    result_local_M(3, 0) = -0.0239386780673502;
    result_local_M(3, 3) = -0.005008211595613448;
    result_local_M(3, 6) = -0.03834379168377978;
    result_local_M(3, 9) = -0.03575905425830305;
    result_local_M(3, 12) = 0.05232905486228125;
    result_local_M(3, 15) = 0.1240363173943513;
    result_local_M(3, 18) = 0.02723466860182212;
    result_local_M(3, 21) = 0.02728616767788518;
    result_local_M(3, 24) = 0.1196261288554766;
    result_local_M(3, 27) = 0.07600635871542957;
    result_local_M(4, 1) = -0.0239386780673502;
    result_local_M(4, 4) = -0.005008211595613448;
    result_local_M(4, 7) = -0.03834379168377978;
    result_local_M(4, 10) = -0.03575905425830305;
    result_local_M(4, 13) = 0.05232905486228125;
    result_local_M(4, 16) = 0.1240363173943513;
    result_local_M(4, 19) = 0.02723466860182212;
    result_local_M(4, 22) = 0.02728616767788518;
    result_local_M(4, 25) = 0.1196261288554766;
    result_local_M(4, 28) = 0.07600635871542957;
    result_local_M(5, 2) = -0.0239386780673502;
    result_local_M(5, 5) = -0.005008211595613448;
    result_local_M(5, 8) = -0.03834379168377978;
    result_local_M(5, 11) = -0.03575905425830305;
    result_local_M(5, 14) = 0.05232905486228125;
    result_local_M(5, 17) = 0.1240363173943513;
    result_local_M(5, 20) = 0.02723466860182212;
    result_local_M(5, 23) = 0.02728616767788518;
    result_local_M(5, 26) = 0.1196261288554766;
    result_local_M(5, 29) = 0.07600635871542957;

    // Results for Kappa.
    result_local_kappa(0) = 0.2935049964090884;
    result_local_kappa(1) = 0.2935049964090884;
    result_local_kappa(2) = 0.2935049964090884;
    result_local_kappa(3) = 0.3234689605021996;
    result_local_kappa(4) = 0.3234689605021996;
    result_local_kappa(5) = 0.3234689605021996;

    // Perform the unit tests.
    PerformMortarPairUnitTest(contact_pair, q_beam, q_beam_rot, q_solid, result_local_D,
        result_local_M, result_local_kappa);
  }

  /**
   * \brief Test a non straight beam in a tet10 element, with line3 mortar shape functions.
   */
  void test_beam_to_solid_meshtying_mortar_hermite3_tet10_line3()
  {
    // Element types.
    typedef GEOMETRYPAIR::t_hermite beam_type;
    typedef GEOMETRYPAIR::t_tet10 solid_type;
    typedef GEOMETRYPAIR::t_line3 lambda_type;

    // Create the mesh tying mortar pair.
    BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<beam_type, solid_type, lambda_type>
        contact_pair = BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<beam_type, solid_type,
            lambda_type>();

    // Definition of variables for this test case.
    LINALG::Matrix<beam_type::n_dof_, 1, double> q_beam;
    LINALG::Matrix<9, 1, double> q_beam_rot;
    LINALG::Matrix<solid_type::n_dof_, 1, double> q_solid;
    LINALG::SerialDenseMatrix local_D;
    LINALG::SerialDenseMatrix local_M;
    LINALG::SerialDenseVector local_kappa;

    // Matrices for the results.
    LINALG::Matrix<lambda_type::n_dof_, beam_type::n_dof_, double> result_local_D(true);
    LINALG::Matrix<lambda_type::n_dof_, solid_type::n_dof_, double> result_local_M(true);
    LINALG::Matrix<lambda_type::n_dof_, 1, double> result_local_kappa(true);

    // Define the geometry of the two elements.
    q_beam(0) = 0.15;
    q_beam(1) = 0.2;
    q_beam(2) = 0.3;
    q_beam(3) = 0.5773502691896255;
    q_beam(4) = 0.5773502691896258;
    q_beam(5) = 0.577350269189626;
    q_beam(6) = 0.65;
    q_beam(7) = 0.1;
    q_beam(8) = 0.1;
    q_beam(9) = 0.8017837257372733;
    q_beam(10) = -0.5345224838248488;
    q_beam(11) = 0.2672612419124244;

    q_beam_rot(0) = 1.674352746442651;
    q_beam_rot(1) = 0.1425949677148126;
    q_beam_rot(2) = 1.0831163124736984;
    q_beam_rot(3) = 1.4331999091513161;
    q_beam_rot(4) = -0.6560404572957742;
    q_beam_rot(5) = -0.2491376152457331;
    q_beam_rot(6) = 0.0;
    q_beam_rot(7) = 0.0;
    q_beam_rot(8) = 0.0;

    // Positional DOFs of the solid.
    q_solid(0) = 0.04580648387340619;
    q_solid(1) = 0.02451732788546579;
    q_solid(2) = -0.0031131666281497;
    q_solid(3) = 0.923762409575691;
    q_solid(4) = -0.0106157274667526;
    q_solid(5) = -0.0141609921026476;
    q_solid(6) = -0.0944365111054278;
    q_solid(7) = 1.069973754602123;
    q_solid(8) = 0.06451121136737775;
    q_solid(9) = -0.04783899250893692;
    q_solid(10) = 0.0888801776910401;
    q_solid(11) = 0.963088029848009;
    q_solid(12) = 0.4108326210187769;
    q_solid(13) = -0.06382048551478623;
    q_solid(14) = 0.0810006716524518;
    q_solid(15) = 0.4728874772489418;
    q_solid(16) = 0.4910049605448334;
    q_solid(17) = -0.07654076843449248;
    q_solid(18) = 0.0891921101646426;
    q_solid(19) = 0.5267989315386162;
    q_solid(20) = -0.03523821168532098;
    q_solid(21) = 0.05851054670001138;
    q_solid(22) = -0.04567202467351366;
    q_solid(23) = 0.4631680679379652;
    q_solid(24) = 0.5925155012591516;
    q_solid(25) = 0.02661668312122867;
    q_solid(26) = 0.4341991229313329;
    q_solid(27) = -0.006248620475302225;
    q_solid(28) = 0.4109215888511089;
    q_solid(29) = 0.4194107884465731;

    // Results for D.
    result_local_D(0, 0) = 0.0936536895927479;
    result_local_D(0, 3) = 0.00502552312033731;
    result_local_D(0, 6) = -0.01315602206049943;
    result_local_D(0, 9) = 0.0004139841970492526;
    result_local_D(1, 1) = 0.0936536895927479;
    result_local_D(1, 4) = 0.00502552312033731;
    result_local_D(1, 7) = -0.01315602206049943;
    result_local_D(1, 10) = 0.0004139841970492526;
    result_local_D(2, 2) = 0.0936536895927479;
    result_local_D(2, 5) = 0.00502552312033731;
    result_local_D(2, 8) = -0.01315602206049943;
    result_local_D(2, 11) = 0.0004139841970492526;
    result_local_D(3, 0) = -0.006233532187104405;
    result_local_D(3, 3) = 0.0004139841970492533;
    result_local_D(3, 6) = 0.1166951638124641;
    result_local_D(3, 9) = -0.006743365500647556;
    result_local_D(4, 1) = -0.006233532187104405;
    result_local_D(4, 4) = 0.0004139841970492533;
    result_local_D(4, 7) = 0.1166951638124641;
    result_local_D(4, 10) = -0.006743365500647556;
    result_local_D(5, 2) = -0.006233532187104405;
    result_local_D(5, 5) = 0.0004139841970492533;
    result_local_D(5, 8) = 0.1166951638124641;
    result_local_D(5, 11) = -0.006743365500647556;
    result_local_D(6, 0) = 0.2018905486414377;
    result_local_D(6, 3) = 0.02560988177666075;
    result_local_D(6, 6) = 0.2241241091122421;
    result_local_D(6, 9) = -0.02726581856485776;
    result_local_D(7, 1) = 0.2018905486414377;
    result_local_D(7, 4) = 0.02560988177666075;
    result_local_D(7, 7) = 0.2241241091122421;
    result_local_D(7, 10) = -0.02726581856485776;
    result_local_D(8, 2) = 0.2018905486414377;
    result_local_D(8, 5) = 0.02560988177666075;
    result_local_D(8, 8) = 0.2241241091122421;
    result_local_D(8, 11) = -0.02726581856485776;

    // Results for M.
    result_local_M(0, 0) = -0.00863699253265766;
    result_local_M(0, 3) = -0.01089260043938194;
    result_local_M(0, 6) = -0.00967306455078242;
    result_local_M(0, 9) = -0.006121527783446023;
    result_local_M(0, 12) = 0.004090063171667997;
    result_local_M(0, 15) = 0.00875086287021659;
    result_local_M(0, 18) = 0.01766213868102817;
    result_local_M(0, 21) = 0.024519405754108;
    result_local_M(0, 24) = 0.01592105929771795;
    result_local_M(0, 27) = 0.0448783230637778;
    result_local_M(1, 1) = -0.00863699253265766;
    result_local_M(1, 4) = -0.01089260043938194;
    result_local_M(1, 7) = -0.00967306455078242;
    result_local_M(1, 10) = -0.006121527783446023;
    result_local_M(1, 13) = 0.004090063171667997;
    result_local_M(1, 16) = 0.00875086287021659;
    result_local_M(1, 19) = 0.01766213868102817;
    result_local_M(1, 22) = 0.024519405754108;
    result_local_M(1, 25) = 0.01592105929771795;
    result_local_M(1, 28) = 0.0448783230637778;
    result_local_M(2, 2) = -0.00863699253265766;
    result_local_M(2, 5) = -0.01089260043938194;
    result_local_M(2, 8) = -0.00967306455078242;
    result_local_M(2, 11) = -0.006121527783446023;
    result_local_M(2, 14) = 0.004090063171667997;
    result_local_M(2, 17) = 0.00875086287021659;
    result_local_M(2, 20) = 0.01766213868102817;
    result_local_M(2, 23) = 0.024519405754108;
    result_local_M(2, 26) = 0.01592105929771795;
    result_local_M(2, 29) = 0.0448783230637778;
    result_local_M(3, 0) = -0.00835797451756794;
    result_local_M(3, 3) = 0.01084771257067845;
    result_local_M(3, 6) = -0.01256998368023192;
    result_local_M(3, 9) = -0.01263500318722335;
    result_local_M(3, 12) = 0.02547039859498008;
    result_local_M(3, 15) = 0.04777875451231568;
    result_local_M(3, 18) = 0.006642737214020257;
    result_local_M(3, 21) = 0.004875142284587578;
    result_local_M(3, 24) = 0.04122904245355294;
    result_local_M(3, 27) = 0.007180805380247922;
    result_local_M(4, 1) = -0.00835797451756794;
    result_local_M(4, 4) = 0.01084771257067845;
    result_local_M(4, 7) = -0.01256998368023192;
    result_local_M(4, 10) = -0.01263500318722335;
    result_local_M(4, 13) = 0.02547039859498008;
    result_local_M(4, 16) = 0.04777875451231568;
    result_local_M(4, 19) = 0.006642737214020257;
    result_local_M(4, 22) = 0.004875142284587578;
    result_local_M(4, 25) = 0.04122904245355294;
    result_local_M(4, 28) = 0.007180805380247922;
    result_local_M(5, 2) = -0.00835797451756794;
    result_local_M(5, 5) = 0.01084771257067845;
    result_local_M(5, 8) = -0.01256998368023192;
    result_local_M(5, 11) = -0.01263500318722335;
    result_local_M(5, 14) = 0.02547039859498008;
    result_local_M(5, 17) = 0.04777875451231568;
    result_local_M(5, 20) = 0.006642737214020257;
    result_local_M(5, 23) = 0.004875142284587578;
    result_local_M(5, 26) = 0.04122904245355294;
    result_local_M(5, 29) = 0.007180805380247922;
    result_local_M(6, 0) = -0.03116140709956454;
    result_local_M(6, 3) = -0.03171184833258379;
    result_local_M(6, 6) = -0.05154761600709574;
    result_local_M(6, 9) = -0.04624810214215941;
    result_local_M(6, 12) = 0.05371731253460236;
    result_local_M(6, 15) = 0.1525151257640712;
    result_local_M(6, 18) = 0.04118386277560374;
    result_local_M(6, 21) = 0.04482205078659521;
    result_local_M(6, 24) = 0.1567941728038473;
    result_local_M(6, 27) = 0.1376511066703633;
    result_local_M(7, 1) = -0.03116140709956454;
    result_local_M(7, 4) = -0.03171184833258379;
    result_local_M(7, 7) = -0.05154761600709574;
    result_local_M(7, 10) = -0.04624810214215941;
    result_local_M(7, 13) = 0.05371731253460236;
    result_local_M(7, 16) = 0.1525151257640712;
    result_local_M(7, 19) = 0.04118386277560374;
    result_local_M(7, 22) = 0.04482205078659521;
    result_local_M(7, 25) = 0.1567941728038473;
    result_local_M(7, 28) = 0.1376511066703633;
    result_local_M(8, 2) = -0.03116140709956454;
    result_local_M(8, 5) = -0.03171184833258379;
    result_local_M(8, 8) = -0.05154761600709574;
    result_local_M(8, 11) = -0.04624810214215941;
    result_local_M(8, 14) = 0.05371731253460236;
    result_local_M(8, 17) = 0.1525151257640712;
    result_local_M(8, 20) = 0.04118386277560374;
    result_local_M(8, 23) = 0.04482205078659521;
    result_local_M(8, 26) = 0.1567941728038473;
    result_local_M(8, 29) = 0.1376511066703633;

    // Results for Kappa.
    result_local_kappa(0) = 0.0804976675322485;
    result_local_kappa(1) = 0.0804976675322485;
    result_local_kappa(2) = 0.0804976675322485;
    result_local_kappa(3) = 0.1104616316253597;
    result_local_kappa(4) = 0.1104616316253597;
    result_local_kappa(5) = 0.1104616316253597;
    result_local_kappa(6) = 0.4260146577536798;
    result_local_kappa(7) = 0.4260146577536798;
    result_local_kappa(8) = 0.4260146577536798;

    // Perform the unit tests.
    PerformMortarPairUnitTest(contact_pair, q_beam, q_beam_rot, q_solid, result_local_D,
        result_local_M, result_local_kappa);
  }

  /**
   * \brief Test a non straight beam in a tet10 element, with line4 mortar shape functions.
   */
  void test_beam_to_solid_meshtying_mortar_hermite3_tet10_line4()
  {
    // Element types.
    typedef GEOMETRYPAIR::t_hermite beam_type;
    typedef GEOMETRYPAIR::t_tet10 solid_type;
    typedef GEOMETRYPAIR::t_line4 lambda_type;

    // Create the mesh tying mortar pair.
    BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<beam_type, solid_type, lambda_type>
        contact_pair = BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<beam_type, solid_type,
            lambda_type>();

    // Definition of variables for this test case.
    LINALG::Matrix<beam_type::n_dof_, 1, double> q_beam;
    LINALG::Matrix<9, 1, double> q_beam_rot;
    LINALG::Matrix<solid_type::n_dof_, 1, double> q_solid;
    LINALG::SerialDenseMatrix local_D;
    LINALG::SerialDenseMatrix local_M;
    LINALG::SerialDenseVector local_kappa;

    // Matrices for the results.
    LINALG::Matrix<lambda_type::n_dof_, beam_type::n_dof_, double> result_local_D(true);
    LINALG::Matrix<lambda_type::n_dof_, solid_type::n_dof_, double> result_local_M(true);
    LINALG::Matrix<lambda_type::n_dof_, 1, double> result_local_kappa(true);

    // Define the geometry of the two elements.
    q_beam(0) = 0.15;
    q_beam(1) = 0.2;
    q_beam(2) = 0.3;
    q_beam(3) = 0.5773502691896255;
    q_beam(4) = 0.5773502691896258;
    q_beam(5) = 0.577350269189626;
    q_beam(6) = 0.65;
    q_beam(7) = 0.1;
    q_beam(8) = 0.1;
    q_beam(9) = 0.8017837257372733;
    q_beam(10) = -0.5345224838248488;
    q_beam(11) = 0.2672612419124244;

    q_beam_rot(0) = 1.674352746442651;
    q_beam_rot(1) = 0.1425949677148126;
    q_beam_rot(2) = 1.0831163124736984;
    q_beam_rot(3) = 1.4331999091513161;
    q_beam_rot(4) = -0.6560404572957742;
    q_beam_rot(5) = -0.2491376152457331;
    q_beam_rot(6) = 0.0;
    q_beam_rot(7) = 0.0;
    q_beam_rot(8) = 0.0;

    // Positional DOFs of the solid.
    q_solid(0) = 0.04580648387340619;
    q_solid(1) = 0.02451732788546579;
    q_solid(2) = -0.0031131666281497;
    q_solid(3) = 0.923762409575691;
    q_solid(4) = -0.0106157274667526;
    q_solid(5) = -0.0141609921026476;
    q_solid(6) = -0.0944365111054278;
    q_solid(7) = 1.069973754602123;
    q_solid(8) = 0.06451121136737775;
    q_solid(9) = -0.04783899250893692;
    q_solid(10) = 0.0888801776910401;
    q_solid(11) = 0.963088029848009;
    q_solid(12) = 0.4108326210187769;
    q_solid(13) = -0.06382048551478623;
    q_solid(14) = 0.0810006716524518;
    q_solid(15) = 0.4728874772489418;
    q_solid(16) = 0.4910049605448334;
    q_solid(17) = -0.07654076843449248;
    q_solid(18) = 0.0891921101646426;
    q_solid(19) = 0.5267989315386162;
    q_solid(20) = -0.03523821168532098;
    q_solid(21) = 0.05851054670001138;
    q_solid(22) = -0.04567202467351366;
    q_solid(23) = 0.4631680679379652;
    q_solid(24) = 0.5925155012591516;
    q_solid(25) = 0.02661668312122867;
    q_solid(26) = 0.4341991229313329;
    q_solid(27) = -0.006248620475302225;
    q_solid(28) = 0.4109215888511089;
    q_solid(29) = 0.4194107884465731;

    // Results for D.
    result_local_D(0, 0) = 0.05756355583405802;
    result_local_D(0, 3) = 0.001955220263108095;
    result_local_D(0, 6) = 0.005745348903101642;
    result_local_D(0, 9) = -0.0002830518321958758;
    result_local_D(1, 1) = 0.05756355583405802;
    result_local_D(1, 4) = 0.001955220263108095;
    result_local_D(1, 7) = 0.005745348903101642;
    result_local_D(1, 10) = -0.0002830518321958758;
    result_local_D(2, 2) = 0.05756355583405802;
    result_local_D(2, 5) = 0.001955220263108095;
    result_local_D(2, 8) = 0.005745348903101642;
    result_local_D(2, 11) = -0.0002830518321958758;
    result_local_D(3, 0) = 0.004620282991405782;
    result_local_D(3, 3) = 0.0002830518321958757;
    result_local_D(3, 6) = 0.06977827920983277;
    result_local_D(3, 9) = -0.002638102150795207;
    result_local_D(4, 1) = 0.004620282991405782;
    result_local_D(4, 4) = 0.0002830518321958757;
    result_local_D(4, 7) = 0.06977827920983277;
    result_local_D(4, 10) = -0.002638102150795207;
    result_local_D(5, 2) = 0.004620282991405782;
    result_local_D(5, 5) = 0.0002830518321958757;
    result_local_D(5, 8) = 0.06977827920983277;
    result_local_D(5, 11) = -0.002638102150795207;
    result_local_D(6, 0) = 0.1839793570166088;
    result_local_D(6, 3) = 0.01881461423793542;
    result_local_D(6, 6) = 0.02734242802628759;
    result_local_D(6, 9) = -0.00813357387408627;
    result_local_D(7, 1) = 0.1839793570166088;
    result_local_D(7, 4) = 0.01881461423793542;
    result_local_D(7, 7) = 0.02734242802628759;
    result_local_D(7, 10) = -0.00813357387408627;
    result_local_D(8, 2) = 0.1839793570166088;
    result_local_D(8, 5) = 0.01881461423793542;
    result_local_D(8, 8) = 0.02734242802628759;
    result_local_D(8, 11) = -0.00813357387408627;
    result_local_D(9, 0) = 0.04314751020500857;
    result_local_D(9, 3) = 0.00999650276080791;
    result_local_D(9, 6) = 0.2247971947249847;
    result_local_D(9, 9) = -0.0225404720113787;
    result_local_D(10, 1) = 0.04314751020500857;
    result_local_D(10, 4) = 0.00999650276080791;
    result_local_D(10, 7) = 0.2247971947249847;
    result_local_D(10, 10) = -0.0225404720113787;
    result_local_D(11, 2) = 0.04314751020500857;
    result_local_D(11, 5) = 0.00999650276080791;
    result_local_D(11, 8) = 0.2247971947249847;
    result_local_D(11, 11) = -0.0225404720113787;

    // Results for M.
    result_local_M(0, 0) = -0.007265864240329137;
    result_local_M(0, 3) = -0.004372962221471238;
    result_local_M(0, 6) = -0.007635186675347714;
    result_local_M(0, 9) = -0.00535769732751274;
    result_local_M(0, 12) = 0.00670419626700892;
    result_local_M(0, 15) = 0.00966849019628169;
    result_local_M(0, 18) = 0.01397432158241253;
    result_local_M(0, 21) = 0.01873931014889322;
    result_local_M(0, 24) = 0.01209072797200314;
    result_local_M(0, 27) = 0.026763569035221;
    result_local_M(1, 1) = -0.007265864240329137;
    result_local_M(1, 4) = -0.004372962221471238;
    result_local_M(1, 7) = -0.007635186675347714;
    result_local_M(1, 10) = -0.00535769732751274;
    result_local_M(1, 13) = 0.00670419626700892;
    result_local_M(1, 16) = 0.00966849019628169;
    result_local_M(1, 19) = 0.01397432158241253;
    result_local_M(1, 22) = 0.01873931014889322;
    result_local_M(1, 25) = 0.01209072797200314;
    result_local_M(1, 28) = 0.026763569035221;
    result_local_M(2, 2) = -0.007265864240329137;
    result_local_M(2, 5) = -0.004372962221471238;
    result_local_M(2, 8) = -0.007635186675347714;
    result_local_M(2, 11) = -0.00535769732751274;
    result_local_M(2, 14) = 0.00670419626700892;
    result_local_M(2, 17) = 0.00966849019628169;
    result_local_M(2, 20) = 0.01397432158241253;
    result_local_M(2, 23) = 0.01873931014889322;
    result_local_M(2, 26) = 0.01209072797200314;
    result_local_M(2, 29) = 0.026763569035221;
    result_local_M(3, 0) = -0.005833926922450891;
    result_local_M(3, 3) = 0.00829205539434072;
    result_local_M(3, 6) = -0.00816440955477965;
    result_local_M(3, 9) = -0.007617820875386703;
    result_local_M(3, 12) = 0.01614160143281385;
    result_local_M(3, 15) = 0.02779673646574168;
    result_local_M(3, 18) = 0.005182571465685431;
    result_local_M(3, 21) = 0.005052481541477967;
    result_local_M(3, 24) = 0.02546010217878683;
    result_local_M(3, 27) = 0.00808917107500932;
    result_local_M(4, 1) = -0.005833926922450891;
    result_local_M(4, 4) = 0.00829205539434072;
    result_local_M(4, 7) = -0.00816440955477965;
    result_local_M(4, 10) = -0.007617820875386703;
    result_local_M(4, 13) = 0.01614160143281385;
    result_local_M(4, 16) = 0.02779673646574168;
    result_local_M(4, 19) = 0.005182571465685431;
    result_local_M(4, 22) = 0.005052481541477967;
    result_local_M(4, 25) = 0.02546010217878683;
    result_local_M(4, 28) = 0.00808917107500932;
    result_local_M(5, 2) = -0.005833926922450891;
    result_local_M(5, 5) = 0.00829205539434072;
    result_local_M(5, 8) = -0.00816440955477965;
    result_local_M(5, 11) = -0.007617820875386703;
    result_local_M(5, 14) = 0.01614160143281385;
    result_local_M(5, 17) = 0.02779673646574168;
    result_local_M(5, 20) = 0.005182571465685431;
    result_local_M(5, 23) = 0.005052481541477967;
    result_local_M(5, 26) = 0.02546010217878683;
    result_local_M(5, 29) = 0.00808917107500932;
    result_local_M(6, 0) = -0.01579891253932227;
    result_local_M(6, 3) = -0.03145085777845103;
    result_local_M(6, 6) = -0.02544398962896501;
    result_local_M(6, 9) = -0.01963452967110964;
    result_local_M(6, 12) = 0.01230159291445312;
    result_local_M(6, 15) = 0.05444029018333145;
    result_local_M(6, 18) = 0.02650739983669832;
    result_local_M(6, 21) = 0.03414855586061757;
    result_local_M(6, 24) = 0.07028880877858722;
    result_local_M(6, 27) = 0.1059634270870567;
    result_local_M(7, 1) = -0.01579891253932227;
    result_local_M(7, 4) = -0.03145085777845103;
    result_local_M(7, 7) = -0.02544398962896501;
    result_local_M(7, 10) = -0.01963452967110964;
    result_local_M(7, 13) = 0.01230159291445312;
    result_local_M(7, 16) = 0.05444029018333145;
    result_local_M(7, 19) = 0.02650739983669832;
    result_local_M(7, 22) = 0.03414855586061757;
    result_local_M(7, 25) = 0.07028880877858722;
    result_local_M(7, 28) = 0.1059634270870567;
    result_local_M(8, 2) = -0.01579891253932227;
    result_local_M(8, 5) = -0.03145085777845103;
    result_local_M(8, 8) = -0.02544398962896501;
    result_local_M(8, 11) = -0.01963452967110964;
    result_local_M(8, 14) = 0.01230159291445312;
    result_local_M(8, 17) = 0.05444029018333145;
    result_local_M(8, 20) = 0.02650739983669832;
    result_local_M(8, 23) = 0.03414855586061757;
    result_local_M(8, 26) = 0.07028880877858722;
    result_local_M(8, 29) = 0.1059634270870567;
    result_local_M(9, 0) = -0.01925767044768784;
    result_local_M(9, 3) = -0.004224971595705737;
    result_local_M(9, 6) = -0.03254707837901769;
    result_local_M(9, 9) = -0.03239458523881971;
    result_local_M(9, 12) = 0.04813038368697454;
    result_local_M(9, 15) = 0.1171392263012487;
    result_local_M(9, 18) = 0.01982444578585588;
    result_local_M(9, 21) = 0.01627625127430204;
    result_local_M(9, 24) = 0.106104635625741;
    result_local_M(9, 27) = 0.04889406791710205;
    result_local_M(10, 1) = -0.01925767044768784;
    result_local_M(10, 4) = -0.004224971595705737;
    result_local_M(10, 7) = -0.03254707837901769;
    result_local_M(10, 10) = -0.03239458523881971;
    result_local_M(10, 13) = 0.04813038368697454;
    result_local_M(10, 16) = 0.1171392263012487;
    result_local_M(10, 19) = 0.01982444578585588;
    result_local_M(10, 22) = 0.01627625127430204;
    result_local_M(10, 25) = 0.106104635625741;
    result_local_M(10, 28) = 0.04889406791710205;
    result_local_M(11, 2) = -0.01925767044768784;
    result_local_M(11, 5) = -0.004224971595705737;
    result_local_M(11, 8) = -0.03254707837901769;
    result_local_M(11, 11) = -0.03239458523881971;
    result_local_M(11, 14) = 0.04813038368697454;
    result_local_M(11, 17) = 0.1171392263012487;
    result_local_M(11, 20) = 0.01982444578585588;
    result_local_M(11, 23) = 0.01627625127430204;
    result_local_M(11, 26) = 0.106104635625741;
    result_local_M(11, 29) = 0.04889406791710205;

    // Results for Kappa.
    result_local_kappa(0) = 0.06330890473715966;
    result_local_kappa(1) = 0.06330890473715966;
    result_local_kappa(2) = 0.06330890473715966;
    result_local_kappa(3) = 0.07439856220123855;
    result_local_kappa(4) = 0.07439856220123855;
    result_local_kappa(5) = 0.07439856220123855;
    result_local_kappa(6) = 0.2113217850428964;
    result_local_kappa(7) = 0.2113217850428964;
    result_local_kappa(8) = 0.2113217850428964;
    result_local_kappa(9) = 0.2679447049299934;
    result_local_kappa(10) = 0.2679447049299934;
    result_local_kappa(11) = 0.2679447049299934;

    // Perform the unit tests.
    PerformMortarPairUnitTest(contact_pair, q_beam, q_beam_rot, q_solid, result_local_D,
        result_local_M, result_local_kappa);
  }

 private:
  /**
   * \brief Set up the contact pair so it can be evaluated and compare the results.
   */
  template <typename beam_type, typename solid_type, typename lambda_type>
  void PerformMortarPairUnitTest(
      BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<beam_type, solid_type, lambda_type>&
          contact_pair,
      const LINALG::Matrix<beam_type::n_dof_, 1, double>& q_beam,
      const LINALG::Matrix<9, 1, double>& q_beam_rot,
      const LINALG::Matrix<solid_type::n_dof_, 1, double>& q_solid,
      const LINALG::Matrix<lambda_type::n_dof_, beam_type::n_dof_, double>& result_local_D,
      const LINALG::Matrix<lambda_type::n_dof_, solid_type::n_dof_, double>& result_local_M,
      const LINALG::Matrix<lambda_type::n_dof_, 1, double>& result_local_kappa)
  {
    // Create the elements.
    const int dummy_node_ids[2] = {0, 1};
    Teuchos::RCP<DRT::Element> beam_element = Teuchos::rcp(new DRT::ELEMENTS::Beam3r(0, 0));
    beam_element->SetNodeIds(2, dummy_node_ids);
    Teuchos::RCP<DRT::Element> solid_element = Teuchos::rcp(new DRT::ELEMENTS::So_hex8(1, 0));

    // Set up the beam element.
    std::vector<double> xrefe(6);
    for (unsigned int j = 0; j < 2; j++)
      for (unsigned int i = 0; i < 3; i++) xrefe[i + 3 * j] = q_beam(i + j * 6);

    // Get the rotational vector.
    std::vector<double> rotrefe(9);
    for (unsigned int i = 0; i < 9; i++) rotrefe[i] = q_beam_rot(i);

    // Cast beam element and set the hermitian interpolation.
    Teuchos::RCP<DRT::ELEMENTS::Beam3r> beam_element_cast =
        Teuchos::rcp_dynamic_cast<DRT::ELEMENTS::Beam3r>(beam_element, true);
    beam_element_cast->SetCenterlineHermite(true);
    beam_element_cast->SetUpReferenceGeometry<3, 2, 2>(xrefe, rotrefe);

    // Call Init on the beam contact pair.
    std::vector<const DRT::Element*> pair_elements;
    pair_elements.push_back(&(*beam_element));
    pair_elements.push_back(&(*solid_element));
    contact_pair.CreateGeometryPair(evaluation_data_);
    contact_pair.Init(Teuchos::null, pair_elements);

    // Evaluate the local matrices.
    LINALG::SerialDenseMatrix local_D;
    LINALG::SerialDenseMatrix local_M;
    LINALG::SerialDenseVector local_kappa;
    contact_pair.CastGeometryPair()->Setup();
    contact_pair.ele1posref_ = q_beam;
    contact_pair.ele2posref_ = q_solid;
    contact_pair.EvaluateDM(local_D, local_M, local_kappa);

    // Check the results for D.
    TS_ASSERT_EQUALS(local_D.M(), lambda_type::n_dof_);
    TS_ASSERT_EQUALS(local_D.N(), beam_type::n_dof_);
    for (int i_row = 0; i_row < local_D.M(); i_row++)
      for (int i_col = 0; i_col < local_D.N(); i_col++)
        TS_ASSERT_DELTA(local_D(i_row, i_col), result_local_D(i_row, i_col), 1e-11);

    // Check the results for M.
    TS_ASSERT_EQUALS(local_M.M(), lambda_type::n_dof_);
    TS_ASSERT_EQUALS(local_M.N(), solid_type::n_dof_);
    for (int i_row = 0; i_row < local_M.M(); i_row++)
      for (int i_col = 0; i_col < local_M.N(); i_col++)
        TS_ASSERT_DELTA(local_M(i_row, i_col), result_local_M(i_row, i_col), 1e-11);

    // Check the results for kappa.
    TS_ASSERT_EQUALS(local_kappa.M(), lambda_type::n_dof_);
    TS_ASSERT_EQUALS(local_kappa.N(), 1);
    for (int i_row = 0; i_row < local_M.M(); i_row++)
      TS_ASSERT_DELTA(local_kappa(i_row), result_local_kappa(i_row), 1e-11);
  }


 private:
  //! Evaluation data container for geometry pairs.
  Teuchos::RCP<GEOMETRYPAIR::LineTo3DEvaluationData> evaluation_data_;
};

#endif
