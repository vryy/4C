/*---------------------------------------------------------------------*/
/*! \file

\brief functionality to pack, unpack and communicate classes with MPI

\level 0


*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_COMM_PAROBJECT_HPP
#define FOUR_C_COMM_PAROBJECT_HPP

#include "4C_config.hpp"

#include "4C_comm_pack_buffer.hpp"
#include "4C_linalg_fixedsizematrix.hpp"
#include "4C_linalg_serialdensematrix.hpp"
#include "4C_linalg_serialdensevector.hpp"
#include "4C_utils_pairedmatrix.hpp"

#include <array>
#include <map>
#include <set>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <vector>

FOUR_C_NAMESPACE_OPEN

namespace Core::Communication
{
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
   int tmp = UniqueParObjectId()         // get the unique parobject id
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
   std::vector<char>::size_type position = 0;                      // used to mark current reading
   position in data int tmp; extract_from_pack(position,data,tmp);    // unpack the unique id if
  (tmp
   != UniqueParObjectId()) FOUR_C_THROW("data does not belong to this class");
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
   sparse matrices or vectors of system length. It does therefore not support any Epetra_Vector or
   Epetra_CrsMatrix objects and is not supposed to in the future either.

   <br>
   Here is a list of data types that are currently supported by the existing add_to_pack and
   extract_from_pack methods \code bool, bool* char, char* enum, enum* int, int* double, double*
   float, float*
   string
   template<typename T, typename U>
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
   existing add_to_pack and extract_from_pack templates. If you do so, please update this
  documentation.

   */
  class ParObject
  {
   private:
    template <class T>
    using is_enum_class =
        std::integral_constant<bool, !std::is_convertible<T, int>::value && std::is_enum<T>::value>;

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
    virtual int UniqueParObjectId() const = 0;

    /*!
     * \brief Pack this class so it can be communicated
     *
     * Resizes the vector data and stores all information of a class in it. The first information to
     * be stored in data has to be the unique parobject id delivered by UniqueParObjectId() which
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
     * parobject id defined at the top of this file and delivered by UniqueParObjectId().
     *
     * \param[in] data vector storing all data to be unpacked into this instance.
     */
    virtual void unpack(const std::vector<char>& data) = 0;

    //@}

   public:
    //! @name Routines to help pack stuff into a char vector

    /*!
     * \brief Add stuff to the end of a char vector data
     *
     * This method is templated for all basic types int char double enum bool etc.
     *
     * \param[in,out] data char string stuff shall be added to
     * \param[in] stuff basic data type (float int double char etc) that get's added to stuff
     *
     * \note To be more precise, you can use this template for types of data that sizeof(kind) works
     * for. Do not use for classes or structs or stl containers!
     */
    static void add_to_pack(PackBuffer& data, const int& stuff) { data.add_to_pack(stuff); }

    static void add_to_pack(PackBuffer& data, const unsigned& stuff) { data.add_to_pack(stuff); }

    static void add_to_pack(PackBuffer& data, const double& stuff) { data.add_to_pack(stuff); }

    /*!
     * \brief Add scoped enums to the end of a char vector data
     *
     * \param[in,out] data Pack buffer where stuff should be added
     * \param[in] stuff scoped enum to be added to the data
     *
     * \note This method is template for scoped enums only. Unscoped enums are currently usually
     * unpacked as ints and then converted, so they are excluded here.
     */
    template <class T,
        typename Enable = typename std::enable_if<is_enum_class<T>::value, void>::type>
    static void add_to_pack(PackBuffer& data, const T& stuff)
    {
      data.add_to_pack<T>(stuff);
    }

    static void add_to_pack(PackBuffer& data, const ParObject& obj);

    static void add_to_pack(PackBuffer& data, const ParObject* obj);

    /*!
     * \brief Add stuff to the end of a char vector data
     *
     * This method is templated for std::vector<T>
     * \param[in,out] data char string stuff shall be added to
     * \param[in] stuff std::vector<T> that get's added to stuff
     */
    template <typename T>
    static void add_to_pack(PackBuffer& data, const std::vector<T>& stuff)
    {
      int numele = stuff.size();
      add_to_pack(data, numele);

      // If T is trivially copyable, we can just copy the bytes. Otherwise, recursively call the
      // pack function for every element.
      if constexpr (std::is_trivially_copyable_v<T>)
      {
        add_to_pack(data, stuff.data(), numele * sizeof(T));
      }
      else
      {
        for (const auto& elem : stuff) add_to_pack(data, elem);
      }
    }

    /*!
     * \brief Add stuff to the end of a char vector data
     *
     * This method is templated for std::array<T, n>
     * \param[in,out] data char string stuff shall be added to
     * \param[in] stuff std::array<T, n> that get's added to stuff
     */
    template <typename T, std::size_t numentries>
    static void add_to_pack(PackBuffer& data, const std::array<T, numentries>& stuff)
    {
      add_to_pack(data, stuff.data(), numentries * sizeof(T));
    }

    /*!
     * \brief Add stuff to the end of a char vector data
     *
     * This method is templated for std::map<T,U>
     * \param[in,out] data char string stuff shall be added to
     * \param[in] stuff std::map<T,U> that get's added to stuff
     */
    template <typename T, typename U>
    static void add_to_pack(PackBuffer& data, const std::map<T, U>& stuff)
    {
      int numentries = (int)stuff.size();
      add_to_pack(data, numentries);

      for (const auto& entry : stuff)
      {
        add_to_pack(data, entry.first);
        add_to_pack(data, entry.second);
      }
    }

    /*!
     * \brief Add stuff to the end of a char vector data
     *
     * This method is templated for std::unordered_map<T,U>
     *
     * \param[in,out] data char string stuff shall be added to
     * \param[in] stuff std::unordered_map<T,U> that get's added to stuff
     */
    template <typename T, typename U>
    static void add_to_pack(PackBuffer& data, const std::unordered_map<T, U>& stuff)
    {
      int numentries = (int)stuff.size();
      add_to_pack(data, numentries);

      for (const auto& entry : stuff)
      {
        add_to_pack(data, entry.first);
        add_to_pack(data, entry.second);
      }
    }

    /*!
     * \brief Add stuff to the end of a char vector data
     *
     * This method is templated for std::pair<T,U>
     * \param[in,out] data char string stuff shall be added to
     * \param[in] stuff std::pair<T,U> that get's added to stuff
     */
    template <typename T, typename U>
    static void add_to_pack(PackBuffer& data, const std::pair<T, U>& stuff)
    {
      add_to_pack(data, stuff.first);
      add_to_pack(data, stuff.second);
    }

    /*!
     * \brief Add stuff to the end of a char vector data
     *
     * This method is a template for Pairedvector<Ts...>
     * \param[in,out] data char string stuff shall be added to
     * \param[in] stuff Pairedvector<Ts...> that get's added to stuff
     */
    template <typename... Ts>
    static void add_to_pack(PackBuffer& data, const Core::Gen::Pairedvector<Ts...>& stuff)
    {
      int numentries = (int)stuff.size();
      add_to_pack(data, numentries);

      int i = 0;
      for (const auto& colcurr : stuff)
      {
        add_to_pack(data, colcurr.first);
        add_to_pack(data, colcurr.second);
        ++i;
      }

      if (i != numentries) FOUR_C_THROW("Something wrong with number of elements");
    }

    /*!
     * \brief Add stuff to the end of a char vector data first
     *
     * This method is a template for std::vector< Pairedvector<Ts...> >
     * \param[in,out] data char string stuff shall be added to
     * \param[in] stuff std::vector<Pairedvector<Ts...> > that get's added to stuff
     */
    template <typename... Ts>
    static void add_to_pack(
        PackBuffer& data, const std::vector<Core::Gen::Pairedvector<Ts...>>& stuff)
    {
      int numentries = (int)stuff.size();
      add_to_pack(data, numentries);

      int i = 0;
      for (auto& paired_vec : stuff)
      {
        add_to_pack(data, paired_vec);
        ++i;
      }

      if (i != numentries) FOUR_C_THROW("Something wrong with number of elements");
    }

    /*!
     * \brief Add stuff to the end of a char vector data first
     *
     * This method is a template for std::vector< pairedmatrix<Ts...> >
     * \param[in,out] data char string stuff shall be added to
     * \param[in] stuff std::vector<pairedmatrix<Ts...> > that get's added to stuff
     */
    template <typename... Ts>
    static void add_to_pack(
        PackBuffer& data, const std::vector<Core::Gen::Pairedmatrix<Ts...>>& stuff)
    {
      int numentries = (int)stuff.size();
      add_to_pack(data, numentries);

      int i = 0;
      for (const typename Core::Gen::PairedmatrixBase<Ts...>::type& paired_mat : stuff)
      {
        add_to_pack(data, paired_mat);
        ++i;
      }

      if (i != numentries) FOUR_C_THROW("Something wrong with number of elements");
    }

    /*!
     * \brief Add stuff to the end of a char vector data
     *
     * This method is templated for std::set<T, U>
     * \param[in,out] data char string stuff shall be added to
     * \param[in] stuff std::set<T> that get's added to stuff
     */
    template <typename T, typename U>
    static void add_to_pack(PackBuffer& data, const std::set<T, U>& stuff)
    {
      int numentries = (int)stuff.size();
      add_to_pack(data, numentries);

      // iterator
      typename std::set<T, U>::const_iterator colcurr;

      int i = 0;
      for (colcurr = stuff.begin(); colcurr != stuff.end(); ++colcurr)
      {
        add_to_pack(data, *colcurr);
        ++i;
      }

      if (i != numentries) FOUR_C_THROW("Something wrong with number of elements");
    }

    /*!
     * \brief Add stuff to the end of a char vector data
     *
     * This method is an overload of the above template for Core::LinAlg::SerialDenseMatrix
     * \param[in,out] data char string stuff shall be added to
     * \param[in] stuff Core::LinAlg::SerialDenseMatrix that get's added to stuff
     */
    static void add_to_pack(PackBuffer& data, const Core::LinAlg::SerialDenseMatrix& stuff);

    /*!
     * \brief Add stuff to the end of a char vector data
     *
     * This method is an overload of the above template for Core::LinAlg::SerialDenseMatrix
     * \param[in,out] data char string stuff shall be added to
     * \param[in] stuff Core::LinAlg::SerialDenseVector that get's added to stuff
     */
    static void add_to_pack(PackBuffer& data, const Core::LinAlg::SerialDenseVector& stuff);

    /*!
     * \brief Add stuff to the end of a char vector data
     *
     * This method is an overload of the above template for Matrix
     * \param[in,out] data char string stuff shall be added to
     * \param[in] stuff Matrix that get's added to data
     */
    template <unsigned i, unsigned j>
    static void add_to_pack(PackBuffer& data, const Core::LinAlg::Matrix<i, j>& stuff)
    {
      add_to_pack(data, i);
      add_to_pack(data, j);
      add_to_pack(data, stuff.data(), stuff.m() * stuff.n() * sizeof(double));
    }

    /*!
     * \brief Add stuff to the end of a char vector data
     *
     * This method is an overload of the above template for an STL vector containing matrices
     * \param[in,out] data char string stuff shall be added to
     * \param[in] stuff vector of matrices that get's added to data
     */
    template <unsigned i, unsigned j>
    static void add_to_pack(PackBuffer& data, const std::vector<Core::LinAlg::Matrix<i, j>>& stuff)
    {
      // add length of vector to be packed so that later the vector can be restored with correct
      // length when unpacked
      int vectorlength = stuff.size();
      add_to_pack(data, vectorlength);

      for (int p = 0; p < vectorlength; ++p)
      {
        const double* A = stuff[p].data();

        // add all data in vector to pack
        add_to_pack(data, A, i * j * sizeof(double));
      }
    }

    /*!
     * \brief Add stuff to the end of a char vector data
     *
     * This method is an overload of the above template for string
     * \param[in,out] data char string stuff shall be added to
     * \param[in] stuff string that get's added to stuff
     */
    static void add_to_pack(PackBuffer& data, const std::string& stuff);

    /*!
     * \brief Add stuff to a char vector data
     *
     * This method adds an array to data
     * \param[in,out] data  char vector stuff shall be added to
     * \param[in] stuff  ptr to stuff that has length stuffsize (in byte)
     * \param[in] stuffsize length of stuff in byte
     */
    template <typename kind>
    static void add_to_pack(PackBuffer& data, const kind* stuff, const int stuffsize)
    {
      data.add_to_pack(stuff, stuffsize);
    }

    /*!
     * \brief Extract stuff from a char vector data and increment position
     *
     * This method is a template for all basic types like int char double enum. To be precise, it
     * will work for all objects where sizeof(kind) is well defined.
     *
     * \param[in,out] position place in data where to extract stuff. Position will be incremented
     * by sizeof(kind)
     * \param[in] data char vector where stuff is extracted from
     * \param[out] stuff basic data type (float int double char ...) to extract from data
     */
    template <typename T>
    static typename std::enable_if<std::is_pod<T>::value, void>::type extract_from_pack(
        std::vector<char>::size_type& position, const std::vector<char>& data, T& stuff)
    {
      int size = sizeof(T);
      memcpy(&stuff, &data[position], size);
      position += size;
    }

    static void extract_from_pack(
        std::vector<char>::size_type& position, const std::vector<char>& data, int& stuff)
    {
      int size = sizeof(int);
      memcpy(&stuff, &data[position], size);
      position += size;
    }

    static void extract_from_pack(
        std::vector<char>::size_type& position, const std::vector<char>& data, unsigned& stuff)
    {
      int size = sizeof(unsigned);
      memcpy(&stuff, &data[position], size);
      position += size;
    }

    static void extract_from_pack(
        std::vector<char>::size_type& position, const std::vector<char>& data, double& stuff)
    {
      int size = sizeof(double);
      memcpy(&stuff, &data[position], size);
      position += size;
    }

    static int extract_int(std::vector<char>::size_type& position, const std::vector<char>& data)
    {
      int i;
      extract_from_pack(position, data, i);
      return i;
    }

    static double extract_double(
        std::vector<char>::size_type& position, const std::vector<char>& data)
    {
      double f;
      extract_from_pack(position, data, f);
      return f;
    }

    /*!
     * \brief Extract stuff from a char vector data and increment position
     *
     * This method is templated for stuff of type std::vector<T>
     *
     * \param[in,out] position place in data where to extract stuff. Position will be incremented
     * by this method
     * \param[in] data char vector where stuff is extracted from
     * \param[out] stuff std::vector<T> to extract from data
     */
    template <typename T>
    static void extract_from_pack(std::vector<char>::size_type& position,
        const std::vector<char>& data, std::vector<T>& stuff)
    {
      int dim = 0;
      extract_from_pack(position, data, dim);
      stuff.resize(dim);

      if constexpr (std::is_trivially_copyable_v<T>)
      {
        extract_from_pack(position, data, stuff.data(), dim * sizeof(T));
      }
      else
      {
        for (auto& elem : stuff) extract_from_pack(position, data, elem);
      }
    }

    /*!
     * \brief Extract stuff from a char vector data and increment position
     *
     * This method is templated for stuff of type std::array<T, n>
     *
     * \param[in,out] position place in data where to extract stuff. Position will be incremented
     * by this method
     * \param[in] data char vector where stuff is extracted from
     * \param[out] stuff std::map<T,n> to extract from data
     */
    template <typename T, std::size_t numentries>
    static void extract_from_pack(std::vector<char>::size_type& position,
        const std::vector<char>& data, std::array<T, numentries>& stuff)
    {
      int size = numentries * sizeof(T);
      extract_from_pack(position, data, stuff.data(), size);
    }

    /*!
     * \brief Extract stuff from a char vector data and increment position
     *
     * This method is templated for stuff of type std::map<T,U>
     *
     * \param[in,out] position place in data where to extract stuff. Position will be incremented by
     * this method
     * \param[in] data char vector where stuff is extracted from
     * \param[out] stuff std::map<T,U> to extract from data
     */
    template <typename T, typename U>
    static void extract_from_pack(std::vector<char>::size_type& position,
        const std::vector<char>& data, std::map<T, U>& stuff)
    {
      int numentries = 0;
      extract_from_pack(position, data, numentries);

      stuff.clear();

      for (int i = 0; i < numentries; i++)
      {
        T first;
        U second;
        extract_from_pack(position, data, first);
        extract_from_pack(position, data, second);

        // add to map
        stuff.insert(std::pair<T, U>(first, second));
      }
    }

    /*!
     * \brief Extract stuff from a char vector data and increment position
     *
     * This method is templated for stuff of type std::unordered_map<T,U>
     *
     * \param[in,out] position place in data where to extract stuff. Position will be incremented
     * by this method
     * \param[in] data  char vector where stuff is extracted from
     * \param[out] stuff std::unordered_map<T,U> to extract from data
     */
    template <typename T, typename U>
    static void extract_from_pack(std::vector<char>::size_type& position,
        const std::vector<char>& data, std::unordered_map<T, U>& stuff)
    {
      int numentries = 0;
      extract_from_pack(position, data, numentries);

      stuff.clear();

      for (int i = 0; i < numentries; i++)
      {
        T first;
        U second;
        extract_from_pack(position, data, first);
        extract_from_pack(position, data, second);

        // add to map
        stuff.insert({first, second});
      }
    }

    /*!
     * \brief Extract stuff from a char vector data and increment position
     *
     * This method is templated for stuff of type std::pair<T,U>
     *
     * \param[in,out] position place in data where to extract stuff. Position will be incremented
     * by this method
     * \param[in] data char vector where stuff is extracted from
     * \param[out] stuff std::pair<T,U> to extract from data
     */
    template <typename T, typename U>
    static void extract_from_pack(std::vector<char>::size_type& position,
        const std::vector<char>& data, std::pair<T, U>& stuff)
    {
      extract_from_pack(position, data, stuff.first);
      extract_from_pack(position, data, stuff.second);
    }

    /*!
     * \brief Extract stuff from a char vector data and increment position
     *
     * This method is templated for stuff of type std::vector<std::pair<T,U>>*
     *
     * \param[in,out] position place in data where to extract stuff. Position will be incremented
     * by this method
     * \param[in] data char vector where stuff is extracted from
     * \param[out] stuff std::vector<std::pair<T,U>> to extract from data
     */
    template <typename T, typename U>
    static void extract_from_pack(std::vector<char>::size_type& position,
        const std::vector<char>& data, std::vector<std::pair<T, U>>& stuff)
    {
      int numentries = 0;
      extract_from_pack(position, data, numentries);

      stuff.clear();

      for (int i = 0; i < numentries; i++)
      {
        T first;
        U second;
        extract_from_pack(position, data, first);
        extract_from_pack(position, data, second);

        // add to map
        stuff.push_back(std::pair<T, U>(first, second));
      }
    }

    /*!
     * \brief Extract stuff from a char vector data and increment position
     *
     * This method is a template for stuff of type Pairedvector<Key,T0,Ts...>*
     *
     * \param[in,out] position place in data where to extract stuff. Position will be incremented
     * by this method
     * \param[in] data char vector where stuff is extracted from
     * \param[out] stuff Pairedvector<Key,T0,Ts...> to extract from data
     */
    template <typename Key, typename T0, typename... Ts>
    static void extract_from_pack(std::vector<char>::size_type& position,
        const std::vector<char>& data, Core::Gen::Pairedvector<Key, T0, Ts...>& stuff)
    {
      int numentries = 0;
      extract_from_pack(position, data, numentries);

      stuff.clear();
      stuff.resize(numentries);

      for (int i = 0; i < numentries; i++)
      {
        Key first;
        T0 second;
        extract_from_pack(position, data, first);
        extract_from_pack(position, data, second);

        // add to map
        stuff[first] = second;
      }
    }

    /*!
     * \brief Extract stuff from a char vector data and increment position
     *
     * This method is a template for stuff of type std::vector< Pairedvector<Ts...> >*
     *
     * \param[in,out] position place in data where to extract stuff. Position will be incremented
     * by this method
     * \param[in] data char vector where stuff is extracted from
     * \param[out] stuff std::vector< Pairedvector<Ts...> > to extract from data
     */
    template <typename... Ts>
    static void extract_from_pack(std::vector<char>::size_type& position,
        const std::vector<char>& data, std::vector<Core::Gen::Pairedvector<Ts...>>& stuff)
    {
      int numentries = 0;
      extract_from_pack(position, data, numentries);

      stuff.clear();
      stuff.resize(numentries);

      Core::Gen::Pairedvector<Ts...> paired_vec;
      for (int i = 0; i < numentries; i++)
      {
        extract_from_pack(position, data, paired_vec);


        // add to map
        stuff[i] = paired_vec;
      }
    }

    /*!
     * \brief Extract stuff from a char vector data and increment position
     *
     * This method is templated for stuff of type std::vector< pairedmatrix<Ts...> >
     *
     * \param[in,out] position place in data where to extract stuff. Position will be incremented
     * by this method
     * \param[in] data char vector where stuff is extracted from
     * \param[out] stuff std::vector< pairedmatrix<Ts...> > to extract from data
     */
    template <typename... Ts>
    static void extract_from_pack(std::vector<char>::size_type& position,
        const std::vector<char>& data, std::vector<Core::Gen::Pairedmatrix<Ts...>>& stuff)
    {
      int numentries = 0;
      extract_from_pack(position, data, numentries);

      stuff.clear();
      stuff.resize(numentries);

      typename Core::Gen::PairedmatrixBase<Ts...>::type paired_mat;
      for (int i = 0; i < numentries; i++)
      {
        extract_from_pack(position, data, paired_mat);

        // add to map
        stuff[i] = paired_mat;
      }
    }

    /*!
     * \brief Extract stuff from a char vector data and increment position
     *
     * This method is templated for stuff of type std::set<T>
     *
     * \param[in,out] position place in data where to extract stuff. Position will be incremented
     * by this method
     * \param[in] data char vector where stuff is extracted from
     * \param[out] stuff std::set<T> to extract from data
     */
    template <typename T, typename U>
    static void extract_from_pack(std::vector<char>::size_type& position,
        const std::vector<char>& data, std::set<T, U>& stuff)
    {
      int numentries = 0;
      extract_from_pack(position, data, numentries);

      stuff.clear();

      for (int i = 0; i < numentries; i++)
      {
        T value;
        extract_from_pack(position, data, value);

        // add to set
        stuff.insert(value);
      }
    }

    /*!
     * \brief Extract stuff from a char vector data and increment position
     *
     * This method is an overload of the above template for stuff of type
     * Core::LinAlg::SerialDenseMatrix
     *
     * \param[in,out] position place in data where to extract stuff. Position will be incremented
     * by this method
     * \param[in] data char vector where stuff is extracted from
     * \param[out] stuff Core::LinAlg::SerialDenseMatrix to extract from data
     */
    static void extract_from_pack(std::vector<char>::size_type& position,
        const std::vector<char>& data, Core::LinAlg::SerialDenseMatrix& stuff);

    /*!
     * \brief Extract stuff from a char vector data and increment position
     *
     * This method is an overload of the above template for stuff of type
     * Core::LinAlg::SerialDenseVector
     *
     * \param[in,out] position place in data where to extract stuff. Position will be incremented
     * by this method
     * \param[in] data char vector where stuff is extracted from
     * \param[out] stuff Core::LinAlg::SerialDenseVector to extract from data
     */
    static void extract_from_pack(std::vector<char>::size_type& position,
        const std::vector<char>& data, Core::LinAlg::SerialDenseVector& stuff);

    /*!
     * \brief Extract stuff from a char vector data and increment position
     *
     * This method is an overload of the above template for stuff of type Matrix
     *
     * \param[in,out] position place in data where to extract stuff. Position will be incremented
     * by this method
     * \param[in] data char vector where stuff is extracted from
     * \param[out] stuff Matrix to extract from data
     */
    template <unsigned int i, unsigned int j>
    static void extract_from_pack(std::vector<char>::size_type& position,
        const std::vector<char>& data, Core::LinAlg::Matrix<i, j>& stuff)
    {
      int m = 0;
      extract_from_pack(position, data, m);
      if (m != i) FOUR_C_THROW("first dimension mismatch");
      int n = 0;
      extract_from_pack(position, data, n);
      if (n != j) FOUR_C_THROW("second dimension mismatch");
      extract_from_pack(position, data, stuff.data(), stuff.m() * stuff.n() * sizeof(double));
    }

    /*!
     * \brief Extract stuff from a char vector data and increment position
     *
     * This method is an overload of the above template for stuff of type STL vector of matrices
     *
     * \param[in,out] position place in data where to extract stuff. Position will be incremented
     * by this method
     * \param[in] data char vector where stuff is extracted from
     * \param[out] stuff STL vector of matrices to extract from data
     */
    template <unsigned int i, unsigned int j>
    static void extract_from_pack(std::vector<char>::size_type& position,
        const std::vector<char>& data, std::vector<Core::LinAlg::Matrix<i, j>>& stuff)
    {
      // get length of vector to be extracted and allocate according amount of memory for all
      // extracted data
      int vectorlength;
      extract_from_pack(position, data, vectorlength);

      // resize vector stuff appropriately
      stuff.resize(vectorlength);

      for (int p = 0; p < vectorlength; ++p)
      {
        double* A = stuff[p].data();

        // actual extraction of data
        extract_from_pack(position, data, A, i * j * sizeof(double));
      }
    }

    /*!
     * \brief Extract stuff from a char vector data and increment position
     *
     * This method is an overload of the above template for stuff of type string
     *
     * \param[in,out] position place in data where to extract stuff. Position will be incremented
     * by this method
     * \param[in] data char vector where stuff is extracted from
     * \param[out] stuff string to extract from data
     */
    static void extract_from_pack(
        std::vector<char>::size_type& position, const std::vector<char>& data, std::string& stuff);

    /*!
     * \brief Extract stuff from a char vector data and increment position
     *
     * This method extracts an array from data
     *
     * \param[in,out] position place in data where to extract stuff. Position will be incremented
     * by stuffsize
     * \param[in] data char string where stuff is extracted from
     * \param[out] stuff array of total length stuffsize (in byte)
     * \param[in] stuffsize length of stuff in byte
     */
    template <typename kind>
    static void extract_from_pack(std::vector<char>::size_type& position,
        const std::vector<char>& data, kind* stuff, const int stuffsize)
    {
      memcpy(stuff, &data[position], stuffsize);
      position += stuffsize;
    }

    //@}
  };
  // class ParObject
}  // namespace Core::Communication


namespace Core::Communication
{
  /*!
   *  \brief Extract the type id at a given @p position from a @p data vector, assert if the
   * extracted type id matches the @p desired_type_id, and return the extracted type id.
   *
   * During extraction the position is incremented.
   *
   * \param[in,out] position Position in data vector where type id is extracted
   * \param[in] data Data vector where type id is extracted from
   * \param[in] desired_type_id Id of the desired type
   */
  int ExtractAndAssertId(std::vector<char>::size_type& position, const std::vector<char>& data,
      const int desired_type_id);
}  // namespace Core::Communication

FOUR_C_NAMESPACE_CLOSE

#endif
