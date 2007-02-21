/*!----------------------------------------------------------------------
\file parobject.H

\class DRT::ParObject

\brief A virtual class with functionality to pack, unpack and communicate
       classes in parallel

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE
#ifndef PAROBJECT_H
#define PAROBJECT_H

#include "vector"
#include "string"
#include "Epetra_SerialDenseMatrix.h"

using namespace std;


// These defines are sued to have a unique identifier for every class
// that implements the pure virtual ParObject class.
// The unique identifier is stored in a communcation array as an integer
// so the receiving processor knows what type of class has been sent.
// the method DRT::Factory(data) in utils.cpp can then be used to
// create an Unpack the correct instance of ParObject
#define ParObject_Container             (int)0
#define ParObject_Condition             (int)1
#define ParObject_DofSet                (int)2
#define ParObject_Node                  (int)3
#define ParObject_DesignNode            (int)4
#define ParObject_Element               (int)5
#define ParObject_ElementRegister       (int)6
#define ParObject_Shell8                (int)7
#define ParObject_Shell8Register        (int)8
#define ParObject_Shell8Line            (int)9


/*!
\brief DRT: namespace of the ccarat discretization module

*/
namespace DRT
{

/*!
\class ParObject

\brief A virtual class with functionality to pack, unpack and communicate
       classes in parallel

\author gee (gee@lnm.mw.tum.de)
*/
class ParObject
{
public:


  // @{ \name Constructors and destructors

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

  //@}

  // @{ \name virtual Packing and Unpacking

  /*!
  \brief Return unique ParObject id

  every class implementing ParObject needs a unique id defined at the
  top of this file.
  */
  virtual int UniqueParObjectId() const = 0;

  /*!
  \brief Pack this class so it can be communicated

  Resizes the vector data and stores all information of a class in it.
  The first information to be stored in data has to be the
  unique parobject delivered by UniqueParObjectId() which will then
  identify the exact class on the receiving processor.

  \param data (in/out): char vector to store class information

  */
  virtual void Pack(vector<char>& data) const = 0;

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
  virtual void Unpack(const vector<char>& data) = 0;


  //@}

public:

  // @{ \name Routines to help pack stuff into a char string

  /*!
  \brief Add stuff to the end of a char vector data

  This method is templated for all basic types int char double enum bool etc.
  \param data (in/out) : char string stuff shall be added to
  \param stuff (in)    : basic data type (float int double char etc) that get's added to stuff

  \note To be more precise, you can use this template for types of data that sizeof(kind)
        works for. Do not use for classes or structs or stl containers!

  */
  template<typename kind> static void AddtoPack(vector<char>& data, const kind& stuff)
  {
    int size = sizeof(kind);
    int oldsize = data.size();
    data.resize(oldsize+size);
    memcpy(&data[oldsize],&stuff,size);
    return;
  }

  /*!
  \brief Add stuff to a char vector data

  This method adds an array to data
  \param data (in/out)  : char vector stuff shall be added to
  \param stuff (in)     : ptr to stuff that has length stuffsize (in byte)
  \param stuffsize (in) : length of stuff in byte
  */
  template<typename kind> static void AddtoPack(vector<char>& data, const kind* stuff, const int stuffsize)
  {
    int oldsize = data.size();
    data.resize(oldsize+stuffsize);
    memcpy(&data[oldsize],stuff,stuffsize);
    return;
  }

  /*!
  \brief Add a vector of a basic type to a char vector data

  This method adds a vector to data
  \param data (in/out) : char vector stuff shall be added to
  \param stuff (in)    : vector of kind
  */
  template<typename kind> static void AddVectortoPack(vector<char>& data, const vector<kind>& stuff)
  {
    int numele = stuff.size();
    AddtoPack(data,numele);
    AddtoPack(data,&stuff[0],numele*sizeof(kind));
    return;
  }

  /*!
  \brief Add an Epetra_SerialDenseMatrix to data

  This method adds a Epetra_SerialDenseMatrix to data

  \param data (in/out) : char vector stuff shall be added to
  \param stuff (in)    : matrix to add
  */
  static void AddMatrixtoPack(vector<char>& data,
                              const Epetra_SerialDenseMatrix& stuff)
  {
    int m = stuff.M();
    int n = stuff.N();
    AddtoPack(data,m);
    AddtoPack(data,n);
    double* A = stuff.A();
    AddtoPack(data,A,n*m*sizeof(double));
    return;
  }

  /*!
  \brief Add a string to a char vector data

  This method adds a string to data
  \param data (in/out) : char string stuff shall be added to
  \param stuff (in)    : string to add
  */
  static void AddStringtoPack(vector<char>& data, const string& stuff)
  {
    int numele = stuff.size();
    AddtoPack(data,numele);
    AddtoPack(data,&stuff[0],numele*sizeof(char));
    return;
  }

  /*!
  \brief Extract stuff from a char vector data and increment position

  This method is templated for all basic type like int char double enum.
  To be precise, it will work for all objects where sizeof(kind) is well defined.

  \param position (in/out) : place in data where to extract stuff. Position will be incremented by sizeof(kind)
  \param data (in)         : char vector where stuff is extracted from
  \param stuff (out)       : basic data type (float int double char ...) to extract from data

  */
  template<typename kind> static void ExtractfromPack(int& position, const vector<char>& data, kind& stuff)
  {
    int size = sizeof(kind);
    memcpy(&stuff,&data[position],size);
    position += size;
    return;
  }

  /*!
  \brief Extract stuff from a char vector data and increment position

  This method extracts an array from data
  \param position (in/out) : place in data where to extract stuff. Position will be incremented by stuffsize
  \param data (in)         : char string where stuff is extracted from
  \param stuff (out)       : array of total length stuffsize (in byte)
  \param stuffsize (in)    : length of stuff in byte
  */
  template<typename kind> static void ExtractfromPack(int& position, const vector<char>& data, kind* stuff, const int stuffsize)
  {
    memcpy(stuff,&data[position],stuffsize);
    position += stuffsize;
    return;
  }

  /*!
  \brief Extract a vector of a basic type from a char vector data and increment position

  This method extracts a vector<kind> from data
  \param position (in/out): place in data where to extract stuff. Position will be incremented
                            by sizeof(int)+stuff.size()*sizeof(kind)
  \param data (in)        : char string where stuff is extracted from
  \param stuff (out)      : vector of total length stuff.size()
  */
  template<typename kind> static void ExtractVectorfromPack(int& position, const vector<char>& data, vector<kind>& stuff)
  {
    int dim = 0;
    ExtractfromPack(position,data,dim);
    stuff.resize(dim);
    int size = dim*sizeof(kind);
    ExtractfromPack(position,data,&stuff[0],size);
    return;
  }

  /*!
  \brief Extract a Epetra_SerialDenseMatrix from a char vector data and increment position

  \param position (in/out) : place in data where to extract stuff. Position will be incremented by this method by
                             the size of the matrix
  \param data (in)         : char vector where stuff is extracted from
  \param stuff (out)       : Epetra_SerialDenseMatrix to be extracted
  */
  static void ExtractMatrixfromPack(int& position, const vector<char>& data,
                             Epetra_SerialDenseMatrix& stuff)
  {
    int m = 0;
    ExtractfromPack(position,data,m);
    int n = 0;
    ExtractfromPack(position,data,n);
    stuff.Reshape(m,n);
    double* A = stuff.A();
    ExtractfromPack(position,data,A,n*m*sizeof(double));
    return;
  }

  /*!
  \brief Extract a string from a char vector data and increment position

  This method extracts a string
  \param position (in/out) : place in data where to extract stuff. Position will be incremented by
                             the sizeof the string in byte + sizeof(int)
  \param data (in)         : char vector where stuff is extracted from
  \param stuff (out)       : string to be extracted to
  */
  static void ExtractStringfromPack(int& position, const vector<char>& data, string& stuff)
  {
    int dim = 0;
    ExtractfromPack(position,data,dim);
    stuff.resize(dim);
    int size = dim*sizeof(char);
    ExtractfromPack(position,data,&stuff[0],size);
    return;
  }
  //@}



}; // class ParObject
} // namespace DRT











#endif  // #ifndef PAROBJECT_H
#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
