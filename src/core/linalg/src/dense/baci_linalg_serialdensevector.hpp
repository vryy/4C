/*----------------------------------------------------------------------*/
/*! \file

\brief Specification of wrapper to a serial vector

\level 0
*/
/*----------------------------------------------------------------------*/
#ifndef BACI_LINALG_SERIALDENSEVECTOR_HPP
#define BACI_LINALG_SERIALDENSEVECTOR_HPP


#include "baci_config.hpp"

#include <Teuchos_SerialDenseVector.hpp>

BACI_NAMESPACE_OPEN

namespace CORE::LINALG
{
  /*!
 \brief A class that wraps Teuchos::SerialDenseVector

      This is done in favor of typedef to allow forward declaration
 */
  class SerialDenseVector : public Teuchos::SerialDenseVector<int, double>
  {
   public:
    /// Base type definition
    using Base = Teuchos::SerialDenseVector<int, double>;

    /// Using the base class constructor
    using Base::SerialDenseVector;
  };

  // type definition for serial integer vector
  typedef Teuchos::SerialDenseVector<int, int> IntSerialDenseVector;

  /*!
    \brief Update vector components with scaled values of a,
           b = alpha*a + beta*b
    */
  void Update(double alpha, const SerialDenseVector& a, double beta, SerialDenseVector& b);

  // wrapper function to compute Norm of vector
  double Norm2(const SerialDenseVector& v);

  // output stream operator
  inline std::ostream& operator<<(std::ostream& out, const SerialDenseVector& vec)
  {
    vec.print(out);
    return out;
  }

  // output stream operator
  inline std::ostream& operator<<(std::ostream& out, const IntSerialDenseVector& vec)
  {
    vec.print(out);
    return out;
  }
}  // namespace CORE::LINALG


BACI_NAMESPACE_CLOSE

#endif  // BACI_LINALG_SERIALDENSEVECTOR_H
