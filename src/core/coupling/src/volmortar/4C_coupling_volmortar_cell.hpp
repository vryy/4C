/*----------------------------------------------------------------------*/
/*! \file

\level 1


*---------------------------------------------------------------------*/

/*---------------------------------------------------------------------*
 | definitions                                             farah 01/14 |
 *---------------------------------------------------------------------*/
#ifndef FOUR_C_COUPLING_VOLMORTAR_CELL_HPP
#define FOUR_C_COUPLING_VOLMORTAR_CELL_HPP

/*---------------------------------------------------------------------*
 | headers                                                 farah 01/14 |
 *---------------------------------------------------------------------*/
#include "4C_config.hpp"

#include "4C_fem_general_element.hpp"

FOUR_C_NAMESPACE_OPEN

/*---------------------------------------------------------------------*
 | forward declarations                                    farah 01/14 |
 *---------------------------------------------------------------------*/
namespace Core::LinAlg
{
  class SerialDenseVector;
  class SerialDenseMatrix;
}  // namespace Core::LinAlg

namespace Core::VolMortar
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
    double CalcJac(const double* xi);

    /*!
    \brief get cell id

    */
    int Id() const { return id_; }

    /*!
    \brief mapping between para space and global space

    */
    void LocalToGlobal(double* local, double* global);

    /*!
    \brief output for coordinates

    */
    void Print();

    /*!
    \brief get shape

    */
    virtual Core::FE::CellType Shape() const { return shape_; }

    /*!
    \brief get cell volume

    */
    virtual double Vol() { return vol_; }

    //@}
   protected:
    int id_;                                  // local ID of this cell
    Core::LinAlg::SerialDenseMatrix coords_;  // coords of cell vertices (dim,vertices)
    Core::FE::CellType shape_;                // shape of this element (always tet4)
    double vol_;                              // integration cell volume
  };

}  // namespace Core::VolMortar

FOUR_C_NAMESPACE_CLOSE

#endif
