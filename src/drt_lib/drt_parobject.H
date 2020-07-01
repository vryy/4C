/*---------------------------------------------------------------------*/
/*! \file

\brief functionality to pack, unpack and communicate classes with MPI

\level 0


*/
/*----------------------------------------------------------------------*/

#ifndef PAROBJECT_H
#define PAROBJECT_H

#include <vector>
#include <map>
#include <set>
#include <type_traits>
#include <string>
#include "../linalg/linalg_serialdensematrix.H"
#include <Epetra_SerialDenseVector.h>
#include "../linalg/linalg_fixedsizematrix.H"
#include "drt_pack_buffer.H"
#include "../headers/pairedmatrix.H"

namespace DRT
{
  /*!
   \brief A virtual class with functionality to pack, unpack and communicate
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
   Epetra_SerialDenseMatrix matrix;
   \endcode
   This is how it is packed into a vector<char>& data:<br>
   \code
   Fool::Pack (vector< char > &data) const
   {
   data.resize(0);                       // resize data to zero
   int tmp = UniqueParObjectId()         // get the unique parobject id
   AddtoPack(data,tmp);                  // pack this id
   AddtoPack(data,i);                    // pack i
   AddtoPack(data,b);                    // pack b
   AddtoPack(data,vec,50*sizeof(double); // pack vec
   AddtoPack(data,bla);                  // pack bla
   AddtoPack(data,matrix);               // pack matrix
   return;
   }
   \endcode
   Here is how this data can be unpacked again:<br>
   \code
   Fool::Unpack(const vector< char > &data)
   {
   std::vector<char>::size_type position = 0;                      // used to mark current reading
   position in data int tmp; ExtractfromPack(position,data,tmp);    // unpack the unique id if (tmp
   != UniqueParObjectId()) dserror("data does not belong to this class");
   ExtractfromPack(position,data,i);      // extract i
   ExtractfromPack(position,data,b);      // extract b
   ExtractfromPack(position,data,bla);    // extract bla
   ExtractfromPack(position,data,matrix); // extract matrix
   if (position != data.size()) dserror("Mismatch in size of data");
   return;
   }
   \endcode
   <br>
   Some remarks:

   - Data has to be unpacked the order it was packed

   - The first object in every packed data set has to be the unique parobject id, see head of file
   drt_parobject.H

   - The size of data ( data.size() ) must 'fit' exactly

   - A class should pack everything it needs to be exactly recreated on a different processor.
   this specifically holds for classes used in a Discretization where data might be shifted around
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
   Here is a list of data types that are currently supported by the existing AddtoPack and
   ExtractfromPack methods: \code bool, bool* char, char* enum, enum* int, int* double, double*
   float, float*
   string
   template<typename T, typename U>
   std::vector<T>  // especially useful to pack other packs into a pack, e.g. a class packing its
   own base class std::map<T,U>   // especially useful to pack other packs into a pack, e.g. a class
   packing its own base class std::pair<T,U> std::set<T> LINALG::Matrix<T,U>
   std::vector<LINALG::Matrix<T,U> >
   Epetra_SerialDenseMatrix
   Epetra_SerialDenseVector
   LINALG::SerialDenseMatrix

   \endcode

   Note that trying to pack an unsupported type of data (such as e.g. std::list<T> ) might compile
   and link but will result in the most (or least) funny behavior. Also, this type of bug might be
   extremely hard to find.... <br><br> Of course, you are welcome to add more specializations to the
   existing AddtoPack and ExtractfromPack templates. If you do so, please update this documentation.

   */
  class ParObject
  {
   public:
    /*!
     \brief Standard Constructor

     */
    ParObject();

    /*!
     \brief Copy Constructor

     Makes a deep copy of a ParObject

     */
    ParObject(const DRT::ParObject& old);

    /*!
     \brief Destructor

     */
    virtual ~ParObject();

    //! @name pure virtual Packing and Unpacking

    /*!
     \brief Return unique ParObject id

     every class implementing ParObject needs a unique id defined at the
     top of drt_parobject.H (this file) and should return it in this method.
     */
    virtual int UniqueParObjectId() const = 0;

    /*!
     \brief Pack this class so it can be communicated

     Resizes the vector data and stores all information of a class in it.
     The first information to be stored in data has to be the
     unique parobject id delivered by UniqueParObjectId() which will then
     identify the exact class on the receiving processor.

     \param data (in/out): char vector to store class information

     */
    virtual void Pack(DRT::PackBuffer& data) const = 0;

    /*!
     \brief Unpack data from a char vector into this class

     The vector data contains all information to rebuild the
     exact copy of an instance of a class on a different processor.
     The first entry in data has to be an integer which is the unique
     parobject id defined at the top of this file and delivered by
     UniqueParObjectId().

     \param data (in) : vector storing all data to be unpacked into this
     instance.
     */
    virtual void Unpack(const std::vector<char>& data) = 0;

    //@}

   public:
    //! @name Routines to help pack stuff into a char vector

    /*!
     \brief Add stuff to the end of a char vector data

     This method is templated for all basic types int char double enum bool etc.
     \param data (in/out) : char string stuff shall be added to
     \param stuff (in)    : basic data type (float int double char etc) that get's added to stuff

     \note To be more precise, you can use this template for types of data that sizeof(kind)
     works for. Do not use for classes or structs or stl containers!

     */
    static void AddtoPack(PackBuffer& data, const int& stuff) { data.AddtoPack(stuff); }

    static void AddtoPack(PackBuffer& data, const unsigned& stuff) { data.AddtoPack(stuff); }

    static void AddtoPack(PackBuffer& data, const double& stuff) { data.AddtoPack(stuff); }

    static void AddtoPack(PackBuffer& data, const ParObject& obj);

    static void AddtoPack(PackBuffer& data, const ParObject* obj);

    /*!
     \brief Add stuff to the end of a char vector data

     This method is templated for std::vector<T>
     \param data (in/out) : char string stuff shall be added to
     \param stuff (in)    : std::vector<T> that get's added to stuff

     */
    template <typename T>
    static void AddtoPack(PackBuffer& data, const std::vector<T>& stuff)
    {
      int numele = stuff.size();
      AddtoPack(data, numele);
      AddtoPack(data, &stuff[0], numele * sizeof(T));
      return;
    }

    /*!
     \brief Add stuff to the end of a char vector data

     This method is templated for std::map<T,U>
     \param data (in/out) : char string stuff shall be added to
     \param stuff (in)    : std::map<T,U> that get's added to stuff

     */
    template <typename T, typename U>
    static void AddtoPack(PackBuffer& data, const std::map<T, U>& stuff)
    {
      int numentries = (int)stuff.size();
      AddtoPack(data, numentries);

      // iterator
      typename std::map<T, U>::const_iterator colcurr;

      int i = 0;
      for (colcurr = stuff.begin(); colcurr != stuff.end(); ++colcurr)
      {
        AddtoPack(data, colcurr->first);
        AddtoPack(data, colcurr->second);
        ++i;
      }

      if (i != numentries) dserror("Something wrong with number of elements");

      return;
    }

    /*!
     \brief Add stuff to the end of a char vector data

     This method is templated for std::pair<T,U>
     \param data (in/out) : char string stuff shall be added to
     \param stuff (in)    : std::pair<T,U> that get's added to stuff

     */
    template <typename T, typename U>
    static void AddtoPack(PackBuffer& data, const std::pair<T, U>& stuff)
    {
      AddtoPack(data, stuff.first);
      AddtoPack(data, stuff.second);

      return;
    }

    /*!
     \brief Add stuff to the end of a char vector data

     This method is a template for pairedvector<Ts...>
     \param data (in/out) : char string stuff shall be added to
     \param stuff (in)    : pairedvector<Ts...> that get's added to stuff

     */
    template <typename... Ts>
    static void AddtoPack(PackBuffer& data, const GEN::pairedvector<Ts...>& stuff)
    {
      int numentries = (int)stuff.size();
      AddtoPack(data, numentries);

      int i = 0;
      for (const auto& colcurr : stuff)
      {
        AddtoPack(data, colcurr.first);
        AddtoPack(data, colcurr.second);
        ++i;
      }

      if (i != numentries) dserror("Something wrong with number of elements");

      return;
    }

    /*!
     \brief Add stuff to the end of a char vector data first

     This method is a template for std::vector< pairedvector<Ts...> >
     \param data (in/out) : char string stuff shall be added to
     \param stuff (in)    : std::vector<pairedvector<Ts...> > that get's added to stuff

     */
    template <typename... Ts>
    static void AddtoPack(PackBuffer& data, const std::vector<GEN::pairedvector<Ts...>>& stuff)
    {
      int numentries = (int)stuff.size();
      AddtoPack(data, numentries);

      int i = 0;
      for (auto& paired_vec : stuff)
      {
        AddtoPack(data, paired_vec);
        ++i;
      }

      if (i != numentries) dserror("Something wrong with number of elements");

      return;
    }

    /*!
     \brief Add stuff to the end of a char vector data first

     This method is a template for std::vector< pairedmatrix<Ts...> >
     \param data (in/out) : char string stuff shall be added to
     \param stuff (in)    : std::vector<pairedmatrix<Ts...> > that get's added to stuff

     */
    template <typename... Ts>
    static void AddtoPack(PackBuffer& data, const std::vector<GEN::pairedmatrix<Ts...>>& stuff)
    {
      int numentries = (int)stuff.size();
      AddtoPack(data, numentries);

      int i = 0;
      for (const typename GEN::pairedmatrix_base<Ts...>::type& paired_mat : stuff)
      {
        AddtoPack(data, paired_mat);
        ++i;
      }

      if (i != numentries) dserror("Something wrong with number of elements");

      return;
    }

    /*!
     \brief Add stuff to the end of a char vector data

     This method is templated for std::set<T, U>
     \param data (in/out) : char string stuff shall be added to
     \param stuff (in)    : std::set<T> that get's added to stuff

     */
    template <typename T, typename U>
    static void AddtoPack(PackBuffer& data, const std::set<T, U>& stuff)
    {
      int numentries = (int)stuff.size();
      AddtoPack(data, numentries);

      // iterator
      typename std::set<T, U>::const_iterator colcurr;

      int i = 0;
      for (colcurr = stuff.begin(); colcurr != stuff.end(); ++colcurr)
      {
        AddtoPack(data, *colcurr);
        ++i;
      }

      if (i != numentries) dserror("Something wrong with number of elements");

      return;
    }

    /*!
     \brief Add stuff to the end of a char vector data

     This method is an overload of the above template for Epetra_SerialDenseMatrix
     \param data (in/out) : char string stuff shall be added to
     \param stuff (in)    : Epetra_SerialDenseMatrix that get's added to stuff

     */
    static void AddtoPack(PackBuffer& data, const Epetra_SerialDenseMatrix& stuff);

    /*!
     \brief Add stuff to the end of a char vector data

     This method is an overload of the above template for Epetra_SerialDenseMatrix
     \param data (in/out) : char string stuff shall be added to
     \param stuff (in)    : Epetra_SerialDenseVector that get's added to stuff

     */
    static void AddtoPack(PackBuffer& data, const Epetra_SerialDenseVector& stuff);

    /*!
     \brief Add stuff to the end of a char vector data

     This method is an overload of the above template for LINALG::SerialDenseMatrix
     \param data (in/out) : char string stuff shall be added to
     \param stuff (in)    : LINALG::SerialDenseMatrix that get's added to stuff

     */
    static void AddtoPack(PackBuffer& data, const LINALG::SerialDenseMatrix& stuff);

    /*!
     \brief Add stuff to the end of a char vector data

     This method is an overload of the above template for Matrix
     \param data (in/out) : char string stuff shall be added to
     \param stuff (in)    : Matrix that get's added to data

     */
    template <unsigned i, unsigned j>
    static void AddtoPack(PackBuffer& data, const LINALG::Matrix<i, j>& stuff)
    {
      AddtoPack(data, i);
      AddtoPack(data, j);
      AddtoPack(data, stuff.A(), stuff.M() * stuff.N() * sizeof(double));
    }

    /*!
     \brief Add stuff to the end of a char vector data

     This method is an overload of the above template for an STL vector containing matrices
     \param data (in/out) : char string stuff shall be added to
     \param stuff (in)    : vector of matrices that get's added to data

     */
    template <unsigned i, unsigned j>
    static void AddtoPack(PackBuffer& data, const std::vector<LINALG::Matrix<i, j>>& stuff)
    {
      // add length of vector to be packed so that later the vector can be restored with correct
      // length when unpacked
      int vectorlength = stuff.size();
      AddtoPack(data, vectorlength);

      for (int p = 0; p < vectorlength; ++p)
      {
        const double* A = stuff[p].A();

        // add all data in vector to pack
        AddtoPack(data, A, i * j * sizeof(double));
      }
      return;
    }

    /*!
     \brief Add stuff to the end of a char vector data

     This method is an overload of the above template for string
     \param data (in/out) : char string stuff shall be added to
     \param stuff (in)    : string that get's added to stuff

     */
    static void AddtoPack(PackBuffer& data, const std::string& stuff);

    /*!
     \brief Add stuff to a char vector data

     This method adds an array to data
     \param data (in/out)  : char vector stuff shall be added to
     \param stuff (in)     : ptr to stuff that has length stuffsize (in byte)
     \param stuffsize (in) : length of stuff in byte
     */
    template <typename kind>
    static void AddtoPack(PackBuffer& data, const kind* stuff, const int stuffsize)
    {
      data.AddtoPack(stuff, stuffsize);
    }

    /*!
     \brief Extract stuff from a char vector data and increment position

     This method is a template for all basic types like int char double enum.
     To be precise, it will work for all objects where sizeof(kind) is well defined.

     \param position (in/out) : place in data where to extract stuff. Position will be incremented
     by sizeof(kind) \param data (in)         : char vector where stuff is extracted from \param
     stuff (out)       : basic data type (float int double char ...) to extract from data

     */
    template <typename T>
    static typename std::enable_if<std::is_pod<T>::value, void>::type ExtractfromPack(
        std::vector<char>::size_type& position, const std::vector<char>& data, T& stuff)
    {
      int size = sizeof(T);
      memcpy(&stuff, &data[position], size);
      position += size;
    }

    static void ExtractfromPack(
        std::vector<char>::size_type& position, const std::vector<char>& data, int& stuff)
    {
      int size = sizeof(int);
      memcpy(&stuff, &data[position], size);
      position += size;
    }

    static void ExtractfromPack(
        std::vector<char>::size_type& position, const std::vector<char>& data, unsigned& stuff)
    {
      int size = sizeof(unsigned);
      memcpy(&stuff, &data[position], size);
      position += size;
    }

    static void ExtractfromPack(
        std::vector<char>::size_type& position, const std::vector<char>& data, double& stuff)
    {
      int size = sizeof(double);
      memcpy(&stuff, &data[position], size);
      position += size;
    }

    static int ExtractInt(std::vector<char>::size_type& position, const std::vector<char>& data)
    {
      int i;
      ExtractfromPack(position, data, i);
      return i;
    }

    static unsigned ExtractUnsigned(
        std::vector<char>::size_type& position, const std::vector<char>& data)
    {
      unsigned i;
      ExtractfromPack(position, data, i);
      return i;
    }

    static double ExtractDouble(
        std::vector<char>::size_type& position, const std::vector<char>& data)
    {
      double f;
      ExtractfromPack(position, data, f);
      return f;
    }



    /*!
     \brief Extract stuff from a char vector data and increment position

     This method is templated for stuff of type std::vector<T>

     \param position (in/out) : place in data where to extract stuff. Position will be incremented
     by this method \param data (in)         : char vector where stuff is extracted from \param
     stuff (out)       : std::vector<T> to extract from data

     */
    template <typename T>
    static void ExtractfromPack(std::vector<char>::size_type& position,
        const std::vector<char>& data, std::vector<T>& stuff)
    {
      int dim = 0;
      ExtractfromPack(position, data, dim);
      stuff.resize(dim);
      int size = dim * sizeof(T);
      ExtractfromPack(position, data, &stuff[0], size);
      return;
    }

    /*!
     \brief Extract stuff from a char vector data and increment position

     This method is templated for stuff of type std::map<T,U>

     \param position (in/out) : place in data where to extract stuff. Position will be incremented
     by this method \param data (in)         : char vector where stuff is extracted from \param
     stuff (out)       : std::map<T,U> to extract from data

     */
    template <typename T, typename U>
    static void ExtractfromPack(std::vector<char>::size_type& position,
        const std::vector<char>& data, std::map<T, U>& stuff)
    {
      int numentries = 0;
      ExtractfromPack(position, data, numentries);

      stuff.clear();

      for (int i = 0; i < numentries; i++)
      {
        T first;
        U second;
        ExtractfromPack(position, data, first);
        ExtractfromPack(position, data, second);

        // add to map
        stuff.insert(std::pair<T, U>(first, second));
      }

      return;
    }

    /*!
     \brief Extract stuff from a char vector data and increment position

     This method is templated for stuff of type std::pair<T,U>

     \param position (in/out) : place in data where to extract stuff. Position will be incremented
     by this method \param data (in)         : char vector where stuff is extracted from \param
     stuff (out)       : std::pair<T,U> to extract from data

     */
    template <typename T, typename U>
    static void ExtractfromPack(std::vector<char>::size_type& position,
        const std::vector<char>& data, std::pair<T, U>& stuff)
    {
      ExtractfromPack(position, data, stuff.first);
      ExtractfromPack(position, data, stuff.second);

      return;
    }

    /*!
     \brief Extract stuff from a char vector data and increment position

     This method is a template for stuff of type pairedvector<Key,T0,Ts...>

     \param position (in/out) : place in data where to extract stuff. Position will be incremented
     by this method \param data (in)         : char vector where stuff is extracted from \param
     stuff (out)       : pairedvector<Key,T0,Ts...> to extract from data

     */
    template <typename Key, typename T0, typename... Ts>
    static void ExtractfromPack(std::vector<char>::size_type& position,
        const std::vector<char>& data, GEN::pairedvector<Key, T0, Ts...>& stuff)
    {
      int numentries = 0;
      ExtractfromPack(position, data, numentries);

      stuff.clear();
      stuff.resize(numentries);

      for (int i = 0; i < numentries; i++)
      {
        Key first;
        T0 second;
        ExtractfromPack(position, data, first);
        ExtractfromPack(position, data, second);

        // add to map
        stuff[first] = second;
      }

      return;
    }

    /*!
     \brief Extract stuff from a char vector data and increment position

     This method is a template for stuff of type std::vector< pairedvector<Ts...> >

     \param position (in/out) : place in data where to extract stuff. Position will be incremented
     by this method \param data (in)         : char vector where stuff is extracted from \param
     stuff (out)       : std::vector< pairedvector<Ts...> > to extract from data

     */
    template <typename... Ts>
    static void ExtractfromPack(std::vector<char>::size_type& position,
        const std::vector<char>& data, std::vector<GEN::pairedvector<Ts...>>& stuff)
    {
      int numentries = 0;
      ExtractfromPack(position, data, numentries);

      stuff.clear();
      stuff.resize(numentries);

      GEN::pairedvector<Ts...> paired_vec;
      for (int i = 0; i < numentries; i++)
      {
        ExtractfromPack(position, data, paired_vec);


        // add to map
        stuff[i] = paired_vec;
      }

      return;
    }

    /*!
     \brief Extract stuff from a char vector data and increment position

     This method is templated for stuff of type std::vector< pairedmatrix<Ts...> >

     \param position (in/out) : place in data where to extract stuff. Position will be incremented
     by this method \param data (in)         : char vector where stuff is extracted from \param
     stuff (out)       : std::vector< pairedmatrix<Ts...> > to extract from data

     */
    template <typename... Ts>
    static void ExtractfromPack(std::vector<char>::size_type& position,
        const std::vector<char>& data, std::vector<GEN::pairedmatrix<Ts...>>& stuff)
    {
      int numentries = 0;
      ExtractfromPack(position, data, numentries);

      stuff.clear();
      stuff.resize(numentries);

      typename GEN::pairedmatrix_base<Ts...>::type paired_mat;
      for (int i = 0; i < numentries; i++)
      {
        ExtractfromPack(position, data, paired_mat);

        // add to map
        stuff[i] = paired_mat;
      }

      return;
    }

    /*!
     \brief Extract stuff from a char vector data and increment position

     This method is templated for stuff of type std::set<T>

     \param position (in/out) : place in data where to extract stuff. Position will be incremented
     by this method \param data (in)         : char vector where stuff is extracted from \param
     stuff (out)       : std::set<T> to extract from data

     */
    template <typename T, typename U>
    static void ExtractfromPack(std::vector<char>::size_type& position,
        const std::vector<char>& data, std::set<T, U>& stuff)
    {
      int numentries = 0;
      ExtractfromPack(position, data, numentries);

      stuff.clear();

      for (int i = 0; i < numentries; i++)
      {
        T value;
        ExtractfromPack(position, data, value);

        // add to set
        stuff.insert(value);
      }

      return;
    }

    /*!
     \brief Extract stuff from a char vector data and increment position

     This method is an overload of the above template for stuff of type Epetra_SerialDenseMatrix

     \param position (in/out) : place in data where to extract stuff. Position will be incremented
     by this method \param data (in)         : char vector where stuff is extracted from \param
     stuff (out)       : Epetra_SerialDenseMatrix to extract from data

     */
    static void ExtractfromPack(std::vector<char>::size_type& position,
        const std::vector<char>& data, Epetra_SerialDenseMatrix& stuff);

    /*!
     \brief Extract stuff from a char vector data and increment position

     This method is an overload of the above template for stuff of type Epetra_SerialDenseVector

     \param position (in/out) : place in data where to extract stuff. Position will be incremented
     by this method \param data (in)         : char vector where stuff is extracted from \param
     stuff (out)       : Epetra_SerialDenseVector to extract from data

     */
    static void ExtractfromPack(std::vector<char>::size_type& position,
        const std::vector<char>& data, Epetra_SerialDenseVector& stuff);

    /*!
     \brief Extract stuff from a char vector data and increment position

     This method is an overload of the above template for stuff of type Matrix

     \param position (in/out) : place in data where to extract stuff. Position will be incremented
     by this method \param data (in)         : char vector where stuff is extracted from \param
     stuff (out)       : Matrix to extract from data

     */
    /*!
     \brief Extract stuff from a char vector data and increment position

     This method is an overload of the above template for stuff of type LINALG::SerialDenseMatrix

     \param position (in/out) : place in data where to extract stuff. Position will be incremented
     by this method \param data (in)         : char vector where stuff is extracted from \param
     stuff (out)       : LINALG::SerialDenseMatrix to extract from data

     */
    static void ExtractfromPack(std::vector<char>::size_type& position,
        const std::vector<char>& data, LINALG::SerialDenseMatrix& stuff);

    template <unsigned int i, unsigned int j>
    static void ExtractfromPack(std::vector<char>::size_type& position,
        const std::vector<char>& data, LINALG::Matrix<i, j>& stuff)
    {
      int m = 0;
      ExtractfromPack(position, data, m);
      if (m != i) dserror("first dimension mismatch");
      int n = 0;
      ExtractfromPack(position, data, n);
      if (n != j) dserror("second dimension mismatch");
      ExtractfromPack(position, data, stuff.A(), stuff.M() * stuff.N() * sizeof(double));
    }

    /*!
     \brief Extract stuff from a char vector data and increment position

     This method is an overload of the above template for stuff of type STL vector of matrices

     \param position (in/out) : place in data where to extract stuff. Position will be incremented
     by this method \param data (in)         : char vector where stuff is extracted from \param
     stuff (out)       : STL vector of matrices to extract from data

     */
    template <unsigned int i, unsigned int j>
    static void ExtractfromPack(std::vector<char>::size_type& position,
        const std::vector<char>& data, std::vector<LINALG::Matrix<i, j>>& stuff)
    {
      // get length of vector to be extracted and allocate according amount of memory for all
      // extracted data
      int vectorlength;
      ExtractfromPack(position, data, vectorlength);

      // resize vector stuff appropriately
      stuff.resize(vectorlength);

      for (int p = 0; p < vectorlength; ++p)
      {
        double* A = stuff[p].A();

        // actual extraction of data
        ExtractfromPack(position, data, A, i * j * sizeof(double));
      }

      return;
    }

    /*!
     \brief Extract stuff from a char vector data and increment position

     This method is an overload of the above template for stuff of type string

     \param position (in/out) : place in data where to extract stuff. Position will be incremented
     by this method \param data (in)         : char vector where stuff is extracted from \param
     stuff (out)       : string to extract from data

     */
    static void ExtractfromPack(
        std::vector<char>::size_type& position, const std::vector<char>& data, std::string& stuff);

    /*!
     \brief Extract stuff from a char vector data and increment position

     This method extracts an array from data
     \param position (in/out) : place in data where to extract stuff. Position will be incremented
     by stuffsize \param data (in)         : char string where stuff is extracted from \param stuff
     (out)       : array of total length stuffsize (in byte) \param stuffsize (in)    : length of
     stuff in byte
     */
    template <typename kind>
    static void ExtractfromPack(std::vector<char>::size_type& position,
        const std::vector<char>& data, kind* stuff, const int stuffsize)
    {
      memcpy(stuff, &data[position], stuffsize);
      position += stuffsize;
      return;
    }

    //@}
  };
  // class ParObject
}  // namespace DRT

#endif  // #ifndef PAROBJECT_H
