/*----------------------------------------------------------------------------*/
/** \file
\brief contains mainly debugging methods for the xcontact contact classes


\level 3
*/
/*----------------------------------------------------------------------------*/

#include <fstream>

#include "../drt_lib/drt_node.H"
#include "../drt_lib/drt_element.H"
#include "../drt_contact/contact_node.H"
#include "../drt_contact/contact_element.H"
#include "../drt_contact_xcontact/xcontact_interface.H"

#ifndef XCONTACT_DEBUG_UTILS_H
#define XCONTACT_DEBUG_UTILS_H

namespace XCONTACT
{
  namespace DEBUG
  {
    /** \brief check the row and column dimensions of input-object \c a
     *
     *  If the weak template parameter is TRUE, the given check value is used
     *  as an upper bound for the dimension, otherwise the dimensions have
     *  to fit exactly.
     *
     *  \param row_check (in) : expected row dimension
     *  \param col_check (in) : expected column dimension
     *  \param a         (in) : check this object
     *
     *  \author hiermeier \date 10/16 */
#if __cplusplus < 201103L
    template <bool weak_row_check, bool weak_col_check, class T>
#else
    template <bool weak_row_check = false, bool weak_col_check = false, class T>
#endif
    inline void CheckDimensions(const T& a, unsigned row_check, unsigned col_check)
    {
      // ToDo Use dsassert as soon as everything is running!
      // ------------------------------------------------------------------------
      // row check
      // ------------------------------------------------------------------------
      if (weak_row_check)
      {
        if (static_cast<unsigned>(a.M()) > row_check)
          dserror(
              "::: WEAK DIMENSION CHECK :::\n"
              "The given matrix has the wrong row dimension!\n"
              "a.M() = %d > %d = row_check",
              a.M(), row_check);
      }
      else
      {
        if (static_cast<unsigned>(a.M()) != row_check)
          dserror(
              "::: STRONG DIMENSION CHECK :::\n"
              "The given matrix has the wrong row dimension!\n"
              "a.M() = %d != %d = row_check",
              a.M(), row_check);
      }
      // ------------------------------------------------------------------------
      // column check
      // ------------------------------------------------------------------------
      if (weak_col_check)
      {
        if (static_cast<unsigned>(a.N()) > col_check)
          dserror(
              "::: WEAK DIMENSION CHECK :::\n"
              "The given matrix has the wrong column dimension!\n"
              "a.N() = %d > %d = col_check",
              a.N(), col_check);
      }
      else
      {
        if (static_cast<unsigned>(a.N()) != col_check)
          dserror(
              "::: STRONG DIMENSION CHECK :::\n"
              "The given matrix has the wrong column dimension!\n"
              "a.N() = %d != %d = col_check",
              a.N(), col_check);
      }
    }

    /** \brief Sanity check, that the matrix was correctly zeroized
     *
     *  \author hiermeier
     *  \date 08/16 */
    template <class T>
    inline void CheckIsZero(const T& mat)
    {
#ifdef DEBUG
      for (unsigned r = 0; r < mat.M(); ++r)
        for (unsigned c = 0; c < mat.N(); ++c)
          // Since we set the value to zero, it should be exactly zero!
          if (mat(r, c) != 0.0)
            dserror("The entry (%d, %d) is not zero! (value=%e)", r, c, mat(r, c));
#endif
    }

    /** \brief Sanity check, that the LINALG::Matrix was correctly zeroized
     *
     *  \author hiermeier
     *  \date 08/16 */
    template <unsigned dim_1, unsigned dim_2, unsigned dim_3, unsigned dim_4>
    inline void CheckIsZero(const LINALG::Matrix<dim_3, dim_4, LINALG::Matrix<dim_1, dim_2>>& mat)
    {
#ifdef DEBUG
      for (unsigned r = 0; r < dim_3; ++r)
        for (unsigned c = 0; c < dim_4; ++c) CheckIsZero(mat(r, c));
#endif
    }

    /** \brief Print double values
     *
     *  \author Hofer
     *  \date 08/16 */
    void PrintDoubleValue(const std::string& name, double a)
    {
      std::cout << name << ":" << std::endl;
      std::cout << "|" << a << std::endl;
    }

    /** \brief Print double array
     *
     *  \author Hofer
     *  \date 07/16 */
    void PrintDoubleArray(const std::string& name, const double a[], int length)
    {
      std::cout << name << ":" << std::endl;
      for (int i = 0; i < length; ++i)
      {
        std::cout << "|" << a[i] << std::endl;
      }
    }

    /** \brief Print double vector
     *
     *  \author Hofer
     *  \date 08/16 */
    void PrintDoubleVector(const std::string& name, const std::vector<double>& a)
    {
      std::cout << name << ":" << std::endl;
      for (unsigned i = 0; i < a.size(); ++i)
      {
        std::cout << "|" << a[i] << std::endl;
      }
    }

    /** \brief Output Gauss point data
     *
     *  \author Hofer
     *  \date 08/16 */
#if __cplusplus < 201103L
    template <unsigned probdim, DRT::Element::DiscretizationType eletype>
    void OutputGaussPoint(const LINALG::Matrix<probdim, 1>& sx, LINALG::Matrix<probdim, 1>& mx,
        double& gN, LINALG::Matrix<1, DRT::UTILS::DisTypeToDim<eletype>::dim>& gN_sxi)
#else
    template <unsigned probdim, DRT::Element::DiscretizationType eletype,
        unsigned dim = DRT::UTILS::DisTypeToDim<eletype>::dim>
    void OutputGaussPoint(const LINALG::Matrix<probdim, 1>& sx, LINALG::Matrix<probdim, 1>& mx,
        double& gN, LINALG::Matrix<1, dim>& gN_sxi)
#endif
    {
      std::ofstream outputfile;
      outputfile.open("../o/GP.txt", std::ios::app);
      outputfile << std::setprecision(16) << sx(0) << ";" << sx(1) << ";" << mx(0) << ";" << mx(1)
                 << ";" << gN << ";" << gN_sxi(0) << std::endl;
      outputfile.close();
    }

    /** \brief Output element geometry
     *
     *  \author Hofer
     *  \date 08/16 */
    void OutputElement(MORTAR::MortarElement& sele)
    {
      std::ofstream outputfile;
      outputfile.open("../o/EG.txt", std::ios::app);

      //  DRT::Element* sele = dynamic_cast<DRT::Element*>(&mrtsele);

      DRT::Node** snodes = sele.Nodes();
      const int snnodes = sele.NumNode();

      outputfile << std::setprecision(16);
      outputfile << "seleId=" << sele.Id() << std::endl;

      for (int k = 0; k < snnodes; ++k)
      {
        //    DRT::Node* snode = snodes[k];
        MORTAR::MortarNode* snode = dynamic_cast<MORTAR::MortarNode*>(snodes[k]);
        outputfile << "snodeId=" << snode->Id() << ","
                   << "x1n=" << snode->xspatial()[0] << ","
                   << "x2n=" << snode->xspatial()[1] << std::endl;
      }

      outputfile.close();
    }

  }  // namespace DEBUG
}  // namespace XCONTACT


#endif  // #ifndef XCONTACT_DEBUG_UTILS_H
