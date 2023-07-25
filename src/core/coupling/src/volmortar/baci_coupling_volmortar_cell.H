/*----------------------------------------------------------------------*/
/*! \file

\level 1


*---------------------------------------------------------------------*/

/*---------------------------------------------------------------------*
 | definitions                                             farah 01/14 |
 *---------------------------------------------------------------------*/
#ifndef BACI_COUPLING_VOLMORTAR_CELL_H
#define BACI_COUPLING_VOLMORTAR_CELL_H

/*---------------------------------------------------------------------*
 | headers                                                 farah 01/14 |
 *---------------------------------------------------------------------*/
#include "baci_lib_element.H"

/*---------------------------------------------------------------------*
 | forward declarations                                    farah 01/14 |
 *---------------------------------------------------------------------*/
namespace CORE::LINALG
{
  class SerialDenseVector;
  class SerialDenseMatrix;
}  // namespace CORE::LINALG

namespace CORE::VOLMORTAR
{
  class Cell
  {
   public:
    /*!
    \brief constructor

    */
    Cell(int id, int nvertices, const Epetra_SerialDenseMatrix& coords,
        const ::DRT::Element::DiscretizationType& shape);

    /*!
    \brief destructor

    */
    virtual ~Cell(){};

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
    virtual ::DRT::Element::DiscretizationType Shape() const { return shape_; }

    /*!
    \brief get cell volume

    */
    virtual double Vol() { return vol_; }

    //@}
   protected:
    int id_;                                    // local ID of this cell
    Epetra_SerialDenseMatrix coords_;           // coords of cell vertices (dim,vertices)
    ::DRT::Element::DiscretizationType shape_;  // shape of this element (always tet4)
    double vol_;                                // integration cell volume
  };

}  // namespace CORE::VOLMORTAR

#endif  // COUPLING_VOLMORTAR_CELL_H
