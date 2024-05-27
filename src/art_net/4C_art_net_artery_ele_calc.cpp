/*----------------------------------------------------------------------*/
/*! \file

\brief main file containing routines for calculation of artery element

\level 3


*----------------------------------------------------------------------*/


#include "4C_art_net_artery_ele_calc.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
DRT::ELEMENTS::ArteryEleCalc<distype>::ArteryEleCalc(
    const int numdofpernode, const std::string& disname)
    : funct_(), deriv_(), tderiv_(), xjm_(), xji_(), derxy_()
{
}
/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
double DRT::ELEMENTS::ArteryEleCalc<distype>::calculate_ele_length(Artery* ele)
{
  // get node coordinates and number of elements per node
  DRT::Node** nodes = ele->Nodes();
  CORE::LINALG::Matrix<3, iel_> xyze;
  // TODO: does this work for line3?
  for (int inode = 0; inode < iel_; inode++)
  {
    const auto& x = nodes[inode]->X();
    xyze(0, inode) = x[0];
    xyze(1, inode) = x[1];
    xyze(2, inode) = x[2];
  }

  // Calculate the length of artery element
  const double L = sqrt(pow(xyze(0, 0) - xyze(0, 1), 2) + pow(xyze(1, 0) - xyze(1, 1), 2) +
                        pow(xyze(2, 0) - xyze(2, 1), 2));

  return L;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
// template classes

// 1D elements
template class DRT::ELEMENTS::ArteryEleCalc<CORE::FE::CellType::line2>;

FOUR_C_NAMESPACE_CLOSE
