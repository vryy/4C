// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_COMM_PAROBJECT_HPP
#define FOUR_C_COMM_PAROBJECT_HPP

#include "4C_config.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Core::Communication
{
  class PackBuffer;
  class UnpackBuffer;

  /*!
  * \brief A virtual class with functionality to pack, unpack and communicate
   classes in parallel

   This class is used to pack information usually stored in a class in a vector<char>.
   This vector<char> can then be used to communicate the contents of a class and to
   read/write binary io. Every class derived from ParObject must basically implement
   the Pack and the Unpack method. There are several methods (most of them template specializations)
   to ease the work of packing/unpacking<br><br>
   Here is an example:
   \code
   // stuff in a class 'Fool' that needs to be packed
   int                      i;
   double                   b;
   double*                  vec = new double[50];
   vector<char>             bla;
   Core::LinAlg::SerialDenseMatrix matrix;
   \endcode
   This is how it is packed into a vector<char>& data:<br>
   \code
   Fool::Pack (vector< char > &data) const
   {
   data.resize(0);                       // resize data to zero
   int tmp = unique_par_object_id()         // get the unique parobject id
   add_to_pack(data,tmp);                  // pack this id
   add_to_pack(data,i);                    // pack i
   add_to_pack(data,b);                    // pack b
   add_to_pack(data,vec,50*sizeof(double); // pack vec
   add_to_pack(data,bla);                  // pack bla
   add_to_pack(data,matrix);               // pack matrix
   return;
   }
   \endcode
   Here is how this data can be unpacked again:<br>
   \code
   Fool::unpack(const vector< char > &data)
   {
                         // used to mark current reading
   position in data int tmp; extract_from_pack(position,data,tmp);    // unpack the unique id if
  (tmp
   != unique_par_object_id()) FOUR_C_THROW("data does not belong to this class");
   extract_from_pack(position,data,i);      // extract i
   extract_from_pack(position,data,b);      // extract b
   extract_from_pack(position,data,bla);    // extract bla
   extract_from_pack(position,data,matrix); // extract matrix
   if (position != data.size()) FOUR_C_THROW("Mismatch in size of data");
   return;
   }
   \endcode
   <br>
   Some remarks:

   - Data has to be unpacked the order it was packed

   - The first object in every packed data set has to be the unique parobject id, see head of file
   lib_parobject.H

   - The size of data ( data.size() ) must 'fit' exactly

   - A class should pack everything it needs to be exactly recreated on a different processor.
   this specifically holds for classes used in a discretization where data might be shifted around
   processors.

   - Every object that carefully implements ParObject can very easily be communicated using the
   Exporter.

   - Every class that carefully implements ParObject can pretty easily be written/read to/from
   binary io

   - A class derived from or more base classes is responsible to also pack and unpack the base
   classes' data by calling the base class' implementation of Pack/Unpack

   - The intention of this class is to pack and communicate rather 'small' units of data. Though
   possible, it is not meant to be used at the system level to communicate huge data sets such as
   sparse matrices or vectors of system length. It does therefore not support any
  Core::LinAlg::Vector<double> or Epetra_CrsMatrix objects and is not supposed to in the future
  either.

   <br>
   Here is a list of data types that are currently supported by the existing
  add_to_pack and extract_from_pack methods \code bool, bool* char, char* enum,
  enum* int, int* double, double* float, float* string template<typename T, typename U>
   std::vector<T>  // especially useful to pack other packs into a pack, e.g. a class packing its
   std::array<T, n>
   own base class
   std::map<T,U>   // especially useful to pack other packs into a pack, e.g. a class
   std::unordered_map<T,U>
   packing its own base class std::pair<T,U> std::set<T> Core::LinAlg::Matrix<T,U>
   std::vector<Core::LinAlg::Matrix<T,U> >
   Core::LinAlg::SerialDenseMatrix
   Core::LinAlg::SerialDenseVector
   Core::LinAlg::SerialDenseMatrix

   \endcode

   Note that trying to pack an unsupported type of data (such as e.g. std::list<T> ) might compile
   and link but will result in the most (or least) funny behavior. Also, this type of bug might be
   extremely hard to find.... <br><br> Of course, you are welcome to add more specializations to the
   existing add_to_pack and extract_from_pack templates. If you do so, please
  update this documentation.

   */
  class ParObject
  {
   public:
    /*!
     * \brief Standard Constructor
     */
    ParObject() = default;

    /*!
     * \brief Destructor
     */
    virtual ~ParObject() = default;


    //! @name Pure virtual packing and unpacking

    /*!
     * \brief Return unique ParObject id
     *
     * Every class implementing ParObject needs a unique id defined at the top of parobject.H (this
     * file) and should return it in this method.
     */
    virtual int unique_par_object_id() const = 0;

    /*!
     * \brief Pack this class so it can be communicated
     *
     * Resizes the vector data and stores all information of a class in it. The first information to
     * be stored in data has to be the unique parobject id delivered by unique_par_object_id() which
     * will then identify the exact class on the receiving processor.
     *
     * \param[in,out] data char vector to store class information
     */
    virtual void pack(PackBuffer& data) const = 0;

    /*!
     * \brief Unpack data from a char vector into this class
     *
     * The vector data contains all information to rebuild the exact copy of an instance of a class
     * on a different processor. The first entry in data has to be an integer which is the unique
     * parobject id defined at the top of this file and delivered by unique_par_object_id().
     *
     * \param[in] data vector storing all data to be unpacked into this instance.
     */
    virtual void unpack(UnpackBuffer& buffer) = 0;

    //@}
  };
  // class ParObject
}  // namespace Core::Communication

FOUR_C_NAMESPACE_CLOSE

#endif
