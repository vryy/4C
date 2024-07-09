/*-----------------------------------------------------------*/
/*! \file

\brief determinant analysis for Solid Hex8 element

\level 3

*/
/*-----------------------------------------------------------*/

#include "4C_so3_hex8_determinant_analysis.hpp"

#include <Teuchos_LAPACK.hpp>

FOUR_C_NAMESPACE_OPEN

const double Discret::ELEMENTS::SoHex8DeterminantAnalysis::bezier_points_[27][3] = {
    {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {1.0, 1.0, 0.0}, {0.0, 1.0, 0.0},  // 1 to 4
    {0.0, 0.0, 1.0}, {1.0, 0.0, 1.0}, {1.0, 1.0, 1.0}, {0.0, 1.0, 1.0},  // 5 to 8
    {0.5, 0.0, 0.0}, {1.0, 0.5, 0.0}, {0.5, 1.0, 0.0}, {0.0, 0.5, 0.0},  // 9 to 12
    {0.0, 0.0, 0.5}, {1.0, 0.0, 0.5}, {1.0, 1.0, 0.5}, {0.0, 1.0, 0.5},  // 13 to 16
    {0.5, 0.0, 1.0}, {1.0, 0.5, 1.0}, {0.5, 1.0, 1.0}, {0.0, 0.5, 1.0},  // 17 to 20
    {0.5, 0.5, 0.0}, {0.5, 0.0, 0.5}, {1.0, 0.5, 0.5}, {0.5, 1.0, 0.5},  // 21 to 24
    {0.0, 0.5, 0.5}, {0.5, 0.5, 1.0}, {0.5, 0.5, 0.5}                    // 25 to 27
};
const unsigned Discret::ELEMENTS::SoHex8DeterminantAnalysis::bezier_indices_[27][3] = {
    {0, 0, 0}, {2, 0, 0}, {2, 2, 0}, {0, 2, 0},  // 1 to 4
    {0, 0, 2}, {2, 0, 2}, {2, 2, 2}, {0, 2, 2},  // 5 to 8
    {1, 0, 0}, {2, 1, 0}, {1, 2, 0}, {0, 1, 0},  // 9 to 12
    {0, 0, 1}, {2, 0, 1}, {2, 2, 1}, {0, 2, 1},  // 13 to 16
    {1, 0, 2}, {2, 1, 2}, {1, 2, 2}, {0, 1, 2},  // 17 to 20
    {1, 1, 0}, {1, 0, 1}, {2, 1, 1}, {1, 2, 1},  // 21 to 24
    {0, 1, 1}, {1, 1, 2}, {1, 1, 1}              // 25 to 27
};
bool Discret::ELEMENTS::SoHex8DeterminantAnalysis::issetup_ = false;
Core::LinAlg::Matrix<27, 20> Discret::ELEMENTS::SoHex8DeterminantAnalysis::map_q_(true);
Core::LinAlg::Matrix<27, 27> Discret::ELEMENTS::SoHex8DeterminantAnalysis::map_l2b_(true);

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<Discret::ELEMENTS::SoHex8DeterminantAnalysis>
Discret::ELEMENTS::SoHex8DeterminantAnalysis::create()
{
  if (not issetup_)
  {
    build_map_lagrange20_to_bezier27();
    build_map_lagrange_to_bezier();
    issetup_ = true;
  }

  return Teuchos::rcp(new SoHex8DeterminantAnalysis);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Discret::ELEMENTS::SoHex8DeterminantAnalysis::build_sub_map_bezier_to_lagrange(
    const double* left, const double* right, Core::LinAlg::Matrix<27, 27>& sub_map_b2l) const
{
  double scale[3] = {0.0, 0.0, 0.0};
  const double* shift = nullptr;

  for (unsigned j = 0; j < 3; ++j) scale[j] = right[j] - left[j];

  shift = left;

  build_map_bezier_to_lagrange(sub_map_b2l, scale, shift);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Discret::ELEMENTS::SoHex8DeterminantAnalysis::build_map_lagrange20_to_bezier27()
{
  map_q_.put_scalar(0.0);

  const double block_template[4][4] = {
      {1.0, 1.0, 0.0, 0.0}, {0.0, 1.0, 1.0, 0.0}, {0.0, 0.0, 1.0, 1.0}, {1.0, 0.0, 0.0, 1.0}};

  for (unsigned i = 0; i < 8; ++i) map_q_(i, i) = 1.0;

  for (unsigned i = 0; i < 4; ++i)
    for (unsigned j = 0; j < 4; ++j) map_q_(i + 8, j) = -0.5 * block_template[i][j];

  for (unsigned i = 0; i < 4; ++i) map_q_(i + 12, i) = -0.5;

  for (unsigned i = 0; i < 4; ++i) map_q_(i + 12, i + 4) = -0.5;

  for (unsigned i = 0; i < 4; ++i)
    for (unsigned j = 0; j < 4; ++j) map_q_(i + 16, j + 4) = -0.5 * block_template[i][j];

  for (unsigned i = 0; i < 12; ++i) map_q_(i + 8, i + 8) = 2.0;

  for (unsigned i = 0; i < 4; ++i) map_q_(20, i) = -0.75;

  for (unsigned i = 0; i < 4; ++i) map_q_(20, i + 8) = 1.0;

  for (unsigned i = 0; i < 4; ++i)
    for (unsigned j = 0; j < 4; ++j) map_q_(i + 21, j) = -0.75 * block_template[i][j];

  for (unsigned i = 0; i < 4; ++i)
    for (unsigned j = 0; j < 4; ++j) map_q_(i + 21, j + 4) = -0.75 * block_template[i][j];

  for (unsigned i = 0; i < 4; ++i) map_q_(i + 21, i + 8) = 1.0;

  for (unsigned i = 0; i < 4; ++i)
    for (unsigned j = 0; j < 4; ++j) map_q_(i + 21, j + 12) = block_template[i][j];

  for (unsigned i = 0; i < 4; ++i) map_q_(i + 21, i + 16) = 1.0;

  for (unsigned i = 0; i < 4; ++i) map_q_(25, i + 4) = -0.75;

  for (unsigned i = 0; i < 4; ++i) map_q_(25, i + 16) = 1.0;

  for (unsigned i = 0; i < 8; ++i) map_q_(26, i) = -0.625;

  for (unsigned i = 8; i < 20; ++i) map_q_(26, i) = 0.5;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Discret::ELEMENTS::SoHex8DeterminantAnalysis::build_map_bezier_to_lagrange(
    Core::LinAlg::Matrix<27, 27>& map_b2l)
{
  const double scale[3] = {1.0, 1.0, 1.0};
  const double shift[3] = {0.0, 0.0, 0.0};
  build_map_bezier_to_lagrange(map_b2l, scale, shift);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Discret::ELEMENTS::SoHex8DeterminantAnalysis::build_map_bezier_to_lagrange(
    Core::LinAlg::Matrix<27, 27>& map_b2l, const double* scale, const double* shift)
{
  // build map from bezier to lagrange
  for (unsigned i = 0; i < 27; ++i)
    for (unsigned j = 0; j < 27; ++j)
    {
      map_b2l(i, j) =
          bezier_func2(shift[0] + scale[0] * bezier_points_[i][0], bezier_indices_[j][0]) *
          bezier_func2(shift[1] + scale[1] * bezier_points_[i][1], bezier_indices_[j][1]) *
          bezier_func2(shift[2] + scale[2] * bezier_points_[i][2], bezier_indices_[j][2]);
    }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Discret::ELEMENTS::SoHex8DeterminantAnalysis::build_map_lagrange_to_bezier()
{
  // build map from bezier to lagrange
  build_map_bezier_to_lagrange(map_l2b_);

  // (0) compute the factorization of the map (LU-decomposition)
  std::vector<int> ipiv(27);
  int info = 0;

  Teuchos::LAPACK<int, double> lapack;
  lapack.GETRF(27, 27, map_l2b_.data(), 27, ipiv.data(), &info);

  if (info) FOUR_C_THROW("Error detected in LAPACK::GETRF. info = %d", info);

  // (1) compute the optimal block size first
  std::vector<double> work(1);
  int lwork = -1;
  lapack.GETRI(27, map_l2b_.data(), 27, ipiv.data(), work.data(), lwork, &info);
  if (info)
    FOUR_C_THROW(
        "Error detected in LAPACK::GETRI during the calculation of "
        "the optimal block size. info = %d",
        info);

  // (2) compute the inverse: map from lagrange to bezier
  lwork = work[0];
  work.resize(lwork);
  lapack.GETRI(27, map_l2b_.data(), 27, ipiv.data(), work.data(), lwork, &info);
  if (info)
    FOUR_C_THROW(
        "Error detected in LAPACK::GETRI during the calculation of "
        "the inverse. info = %d",
        info);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool Discret::ELEMENTS::SoHex8DeterminantAnalysis::is_valid(
    const Core::LinAlg::Matrix<NUMDIM_SOH8, NUMNOD_SOH8>& x_curr, unsigned* rc) const
{
  // placed on the stack and, thus, not static
  Core::LinAlg::Matrix<20, 1> tet4_volumes;
  compute20_tet4_volumes(tet4_volumes, x_curr);

#ifdef DEBUG_SO_HEX8_DET_ANALYSIS
  std::cout << "20 TET4 volumes:\n";
  tet4_volumes.print(std::cout);
#endif

  // Check the TET4 volumes. If one is not positive, return FALSE
  if (has_invalid_entry(tet4_volumes.data(), 20)) return false;

  BezierCube bcube;
  bcube.init();

  // place on stack and thus not static
  Core::LinAlg::Matrix<27, 1> bcoeffs(false);
  bcoeffs.multiply(map_q_, tet4_volumes);

  // Check the remaining Bezier coefficients. If all are positive, return TRUE.
  if (not has_invalid_entry(bcoeffs.data() + 8, 19)) return true;

  // init recursion counter
  unsigned rcount = 0;
  unsigned* rcount_ptr = nullptr;
  if (rc)
    rcount_ptr = rc;
  else
    rcount_ptr = &rcount;

  const bool isvalid = recursive_subdivision(bcoeffs, bcube.left_, bcube.right_, *rcount_ptr);

#ifdef DEBUG_SO_HEX8_DET_ANALYSIS
  std::cout << "The current element is " << (isvalid ? "valid" : "invalid") << ".\n";
#endif

  return isvalid;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Discret::ELEMENTS::SoHex8DeterminantAnalysis::get_sub_cube_borders(
    const double* l, const double* r, std::list<BezierCube>& subcubes) const
{
  const std::array<double, 3> rl = {0.5 * (r[0] + l[0]), 0.5 * (r[1] + l[1]), 0.5 * (r[2] + l[2])};
  auto it = subcubes.begin();

  // sub-cube #0
  double* sl = it->left_;
  double* sr = it->right_;
  sl[0] = l[0];
  sl[1] = l[1];
  sl[2] = l[2];
  sr[0] = rl[0];
  sr[1] = rl[1];
  sr[2] = rl[2];

  // sub-cube #1
  sl = (++it)->left_;
  sr = it->right_;
  sl[0] = rl[0];
  sl[1] = l[1];
  sl[2] = l[2];
  sr[0] = r[0];
  sr[1] = rl[1];
  sr[2] = rl[2];

  // sub-cube #2
  sl = (++it)->left_;
  sr = it->right_;
  sl[0] = rl[0];
  sl[1] = rl[1];
  sl[2] = l[2];
  sr[0] = r[0];
  sr[1] = r[1];
  sr[2] = rl[2];

  // sub-cube #3
  sl = (++it)->left_;
  sr = it->right_;
  sl[0] = l[0];
  sl[1] = rl[1];
  sl[2] = l[2];
  sr[0] = rl[0];
  sr[1] = r[1];
  sr[2] = rl[2];

  // ----------------------------------------------------------

  // sub-cube #4
  sl = (++it)->left_;
  sr = it->right_;
  sl[0] = l[0];
  sl[1] = l[1];
  sl[2] = rl[2];
  sr[0] = rl[0];
  sr[1] = rl[1];
  sr[2] = r[2];

  // sub-cube #5
  sl = (++it)->left_;
  sr = it->right_;
  sl[0] = rl[0];
  sl[1] = l[1];
  sl[2] = rl[2];
  sr[0] = r[0];
  sr[1] = rl[1];
  sr[2] = r[2];

  // sub-cube #6
  sl = (++it)->left_;
  sr = it->right_;
  sl[0] = rl[0];
  sl[1] = rl[1];
  sl[2] = rl[2];
  sr[0] = r[0];
  sr[1] = r[1];
  sr[2] = r[2];

  // sub-cube #7
  sl = (++it)->left_;
  sr = it->right_;
  sl[0] = l[0];
  sl[1] = rl[1];
  sl[2] = rl[2];
  sr[0] = rl[0];
  sr[1] = r[1];
  sr[2] = r[2];
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool Discret::ELEMENTS::SoHex8DeterminantAnalysis::recursive_subdivision(
    const Core::LinAlg::Matrix<27, 1>& bcoeffs, const double* left, const double* right,
    unsigned& rcount) const
{
  // stop the recursion
  if (++rcount >= 1000) return false;

  std::list<BezierCube> subcubes;
  subcubes.assign(8, BezierCube());
  get_sub_cube_borders(left, right, subcubes);

  // don't allocate new storage in a recursive call
  static Core::LinAlg::Matrix<27, 1> sub_bcoeffs(false);
  auto it = subcubes.begin();
  while (it != subcubes.end())
  {
    get_bezier_coeffs_of_subdomain(bcoeffs, *it, sub_bcoeffs);

    // check first 8 entries of each column
    if (has_invalid_entry(sub_bcoeffs.data(), 8)) return false;

    /* If none of the remaining Bezier coefficients is invalid, erase the
     * corresponding sub-cube. */
    if (not has_invalid_entry(&sub_bcoeffs(8, 0), 19))
    {
      it = subcubes.erase(it);
    }
    else
    {
      ++it;
    }
  }

#ifdef DEBUG_SO_HEX8_DET_ANALYSIS
  std::cout << "\n------------------ RECURSION ------------------\n\n";
#endif

  /* By merging this and the previous loop it would be possible to speed-up
   * the algorithm even more. However, I decided to leave it in this way, since
   * it makes the debugging output better readable.          hiermeier 09/18 */
  for (auto& subcube : subcubes)
  {
    /* if there are undetermined cases left, do a recursive call and
     * subdivide the current subdivision once more. */
    if (not recursive_subdivision(bcoeffs, subcube.left_, subcube.right_, rcount)) return false;
  }

  return true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Discret::ELEMENTS::SoHex8DeterminantAnalysis::get_bezier_coeffs_of_subdomain(
    const Core::LinAlg::Matrix<27, 1>& bcoeffs, BezierCube& subcube,
    Core::LinAlg::Matrix<27, 1>& sub_bcoeffs) const
{
  static Core::LinAlg::Matrix<27, 27> sub_map_b2l(false);
  static Core::LinAlg::Matrix<27, 1> tmp(false);

  build_sub_map_bezier_to_lagrange(subcube.left_, subcube.right_, sub_map_b2l);

  tmp.multiply(sub_map_b2l, bcoeffs);
  sub_bcoeffs.multiply(map_l2b_, tmp);

#ifdef DEBUG_SO_HEX8_DET_ANALYSIS
  std::cout << "Bezier coefficients:\n";
  sub_bcoeffs.print(std::cout);
  subcube.print(std::cout);
#endif
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool Discret::ELEMENTS::SoHex8DeterminantAnalysis::has_invalid_entry(
    const double* entries, const unsigned length) const
{
  for (unsigned i = 0; i < length; ++i)
    if (entries[i] <= 0) return true;

  return false;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
unsigned Discret::ELEMENTS::SoHex8DeterminantAnalysis::compute_tet4_vol_at_corners(
    Core::LinAlg::Matrix<20, 1>& tet4_volumes,
    const Core::LinAlg::Matrix<NUMDIM_SOH8, NUMNOD_SOH8>& x_curr,
    const std::function<unsigned(unsigned i)>& f_index0,
    const std::function<unsigned(unsigned i)>& f_index1,
    const std::function<unsigned(unsigned i)>& f_index2,
    const std::function<unsigned(unsigned i)>& f_index3, unsigned offset) const
{
  static Core::LinAlg::Matrix<NUMDIM_SOH8, 4> tet4_coords;

  // --------------------------------------------------------------------------
  // corner TET4s from 0 to 3
  for (unsigned i = 0; i < 4; ++i)
  {
    // coordinates of the first TET4 node
    unsigned j = f_index0(i);  // 0, 1, 2, 3
    std::copy(&x_curr(0, j), &x_curr(0, j) + NUMDIM_SOH8, tet4_coords.data());

    // coordinates of the second TET4 node
    j = f_index1(i);  // 1, 2, 3, 0
    std::copy(&x_curr(0, j), &x_curr(0, j) + NUMDIM_SOH8, tet4_coords.data() + NUMDIM_SOH8);

    // coordinates of the third TET4 node
    j = f_index2(i);  // 3, 0, 1, 2
    std::copy(&x_curr(0, j), &x_curr(0, j) + NUMDIM_SOH8, tet4_coords.data() + 2 * NUMDIM_SOH8);

    // coordinates of the fourth TET4 node
    j = f_index3(i);  // 4, 5, 6, 7
    std::copy(&x_curr(0, j), &x_curr(0, j) + NUMDIM_SOH8, tet4_coords.data() + 3 * NUMDIM_SOH8);

    tet4_volumes(i + offset) = compute_tet4_volume(tet4_coords);
  }

  return offset + 4;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
unsigned Discret::ELEMENTS::SoHex8DeterminantAnalysis::compute_tet4_vol_at_edges(
    Core::LinAlg::Matrix<20, 1>& tet4_volumes,
    const Core::LinAlg::Matrix<NUMDIM_SOH8, NUMNOD_SOH8>& x_curr,
    const std::function<unsigned(unsigned i)>& f_index0,
    const std::function<unsigned(unsigned i)>& f_index1,
    const std::function<unsigned(unsigned i)>& f_index20,
    const std::function<unsigned(unsigned i)>& f_index21,
    const std::function<unsigned(unsigned i)>& f_index30,
    const std::function<unsigned(unsigned i)>& f_index31, unsigned offset) const
{
  static Core::LinAlg::Matrix<NUMDIM_SOH8, 4> tet4_coords;

  Core::LinAlg::Matrix<NUMDIM_SOH8, 1> tmp;
  for (unsigned i = 0; i < 4; ++i)
  {
    unsigned j = f_index0(i);
    std::copy(&x_curr(0, j), &x_curr(0, j) + NUMDIM_SOH8, tet4_coords.data());

    j = f_index1(i);
    std::copy(&x_curr(0, j), &x_curr(0, j) + NUMDIM_SOH8, tet4_coords.data() + NUMDIM_SOH8);

    j = f_index20(i);
    unsigned jj = f_index21(i);
    for (unsigned k = 0; k < NUMDIM_SOH8; ++k) tmp(k, 0) = 0.5 * (x_curr(k, j) + x_curr(k, jj));
    std::copy(tmp.data(), tmp.data() + NUMDIM_SOH8, tet4_coords.data() + 2 * NUMDIM_SOH8);

    j = f_index30(i);
    jj = f_index31(i);
    for (unsigned k = 0; k < NUMDIM_SOH8; ++k) tmp(k, 0) = 0.5 * (x_curr(k, j) + x_curr(k, jj));
    std::copy(tmp.data(), tmp.data() + NUMDIM_SOH8, tet4_coords.data() + 3 * NUMDIM_SOH8);

    tet4_volumes(i + offset) = compute_tet4_volume(tet4_coords);
  }

  return offset + 4;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Discret::ELEMENTS::SoHex8DeterminantAnalysis::compute20_tet4_volumes(
    Core::LinAlg::Matrix<20, 1>& tet4_volumes,
    const Core::LinAlg::Matrix<NUMDIM_SOH8, NUMNOD_SOH8>& x_curr) const
{
  // --------------------------------------------------------------------------
  // corner TET4s from 0 to 3
  unsigned offset = 0;
  offset = compute_tet4_vol_at_corners(
      tet4_volumes, x_curr, [](unsigned i) -> unsigned { return i; },  // 0, 1, 2, 3
      [](unsigned i) -> unsigned { return fast_mod(i + 1, 4); },       // 1, 2, 3, 0
      [](unsigned i) -> unsigned { return fast_mod(i + 3, 4); },       // 3, 0, 1, 2
      [](unsigned i) -> unsigned { return (i + 4); },                  // 4, 5, 6, 7
      offset);

  // --------------------------------------------------------------------------
  // corner TET4s from 4 to 7
  offset = compute_tet4_vol_at_corners(
      tet4_volumes, x_curr, [](unsigned i) -> unsigned { return i; },  // 0, 1, 2, 3
      [](unsigned i) -> unsigned { return fast_mod(i + 1, 4) + 4; },   // 5, 6, 7, 4
      [](unsigned i) -> unsigned { return fast_mod(i + 3, 4) + 4; },   // 7, 4, 5, 6
      [](unsigned i) -> unsigned { return i + 4; },                    // 4, 5, 6, 7
      offset);

  // --------------------------------------------------------------------------
  // edge TET4s from 8 to 11
  offset = compute_tet4_vol_at_edges(
      tet4_volumes, x_curr,
      [](unsigned i) -> unsigned { return fast_mod(i + 1, 4); },      // 1, 2, 3, 0
      [](unsigned i) -> unsigned { return i; },                       // 0, 1, 2, 3
      [](unsigned i) -> unsigned { return i + 4; },                   // 4, 5, 6, 7
      [](unsigned i) -> unsigned { return fast_mod(i + 1, 4) + 4; },  // 5, 6, 7, 4
      [](unsigned i) -> unsigned { return fast_mod(i + 3, 4); },      // 3, 0, 1, 2
      [](unsigned i) -> unsigned { return fast_mod(i + 2, 4); },      // 2, 3, 0, 1
      offset);

  // --------------------------------------------------------------------------
  // edge TET4s from 12 to 15
  offset = compute_tet4_vol_at_edges(
      tet4_volumes, x_curr, [](unsigned i) -> unsigned { return i; },  // 0, 1, 2, 3
      [](unsigned i) -> unsigned { return (i + 4); },                  // 4, 5, 6, 7
      [](unsigned i) -> unsigned { return fast_mod(i + 1, 4); },       // 1, 2, 3, 0
      [](unsigned i) -> unsigned { return fast_mod(i + 1, 4) + 4; },   // 5, 6, 7, 4
      [](unsigned i) -> unsigned { return fast_mod(i + 3, 4); },       // 3, 0, 1, 2
      [](unsigned i) -> unsigned { return fast_mod(i + 3, 4) + 4; },   // 7, 4, 5, 6
      offset);

  // --------------------------------------------------------------------------
  // edge TET4s from 16 to 19
  offset = compute_tet4_vol_at_edges(
      tet4_volumes, x_curr, [](unsigned i) -> unsigned { return i + 4; },  // 4, 5, 6, 7
      [](unsigned i) -> unsigned { return fast_mod(i + 1, 4) + 4; },       // 5, 6, 7, 4
      [](unsigned i) -> unsigned { return i; },                            // 0, 1, 2, 3
      [](unsigned i) -> unsigned { return fast_mod(i + 1, 4); },           // 1, 2, 3, 0
      [](unsigned i) -> unsigned { return fast_mod(i + 2, 4) + 4; },       // 6, 7, 4, 5
      [](unsigned i) -> unsigned { return fast_mod(i + 3, 4) + 4; },       // 7, 4, 5, 6
      offset);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double Discret::ELEMENTS::SoHex8DeterminantAnalysis::compute_tet4_volume(
    const Core::LinAlg::Matrix<NUMDIM_SOH8, 4>& tet4_ncoords) const
{
  const Core::LinAlg::Matrix<3, 1> xref(tet4_ncoords.data(), true);
  Core::LinAlg::Matrix<3, 3> v;

  for (unsigned i = 0; i < 3; ++i)
  {
    Core::LinAlg::Matrix<3, 1> vr(v.data() + i * 3, true);
    const Core::LinAlg::Matrix<3, 1> xother(tet4_ncoords.data() + (i + 1) * 3, true);
    vr.update(1.0, xother, -1.0, xref);
  }

  const Core::LinAlg::Matrix<3, 1> a(v.data(), true);
  const Core::LinAlg::Matrix<3, 1> b(v.data() + 3, true);
  const Core::LinAlg::Matrix<3, 1> c(v.data() + 6, true);

  // compute the triple product
  return c(0) * (a(1) * b(2) - a(2) * b(1)) + c(1) * (a(2) * b(0) - a(0) * b(2)) +
         c(2) * (a(0) * b(1) - a(1) * b(0));
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Discret::ELEMENTS::SoHex8DeterminantAnalysis::BezierCube::print(std::ostream& os) const
{
  os << std::string(20, '*') << "\n";
  os << "Bezier-Cube (addr. = " << this << ")\n";
  os << std::string(20, '-') << "\n";

  os << "Left borders = { ";
  for (unsigned i = 0; i < 3; ++i) os << (i > 0 ? ", " : "") << left_[i];
  os << " };\n";

  os << "Right borders = { ";
  for (unsigned i = 0; i < 3; ++i) os << (i > 0 ? ", " : "") << right_[i];
  os << " };\n";

  os << std::string(20, '*') << "\n";
}

FOUR_C_NAMESPACE_CLOSE
