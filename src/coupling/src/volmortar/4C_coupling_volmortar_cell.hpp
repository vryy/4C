// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_COUPLING_VOLMORTAR_CELL_HPP
#define FOUR_C_COUPLING_VOLMORTAR_CELL_HPP

/*---------------------------------------------------------------------*
 | headers                                                 farah 01/14 |
 *---------------------------------------------------------------------*/
#include "4C_config.hpp"

#include "4C_fem_general_element.hpp"
#include "4C_linalg_serialdensematrix.hpp"

FOUR_C_NAMESPACE_OPEN

/*---------------------------------------------------------------------*
 | forward declarations                                    farah 01/14 |
 *---------------------------------------------------------------------*/
namespace Core::LinAlg
{
  class SerialDenseVector;
  class SerialDenseMatrix;
}  // namespace Core::LinAlg

namespace Coupling::VolMortar
{
  class Cell
  {
   public:
    /*!
    \brief constructor

    */
    Cell(int id, int nvertices, const Core::LinAlg::SerialDenseMatrix& coords,
        const Core::FE::CellType& shape);

    /*!
    \brief destructor

    */
    virtual ~Cell() = default;

    /*!
    \brief calc jacobian

    */
    double calc_jac(const double* xi);

    /*!
    \brief get cell id

    */
    int id() const { return id_; }

    /*!
    \brief mapping between para space and global space

    */
    void local_to_global(double* local, double* global);

    /*!
    \brief output for coordinates

    */
    void print();

    /*!
    \brief get shape

    */
    virtual Core::FE::CellType shape() const { return shape_; }

    /*!
    \brief get cell volume

    */
    virtual double vol() { return vol_; }

    //@}
   protected:
    int id_;                                  // local ID of this cell
    Core::LinAlg::SerialDenseMatrix coords_;  // coords of cell vertices (dim,vertices)
    Core::FE::CellType shape_;                // shape of this element (always tet4)
    double vol_;                              // integration cell volume
  };

}  // namespace Coupling::VolMortar

FOUR_C_NAMESPACE_CLOSE

#endif
