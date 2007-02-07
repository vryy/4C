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
#define ParObject_DesignElement         (int)7
#define ParObject_DesignElementRegister (int)8
#define ParObject_Shell8                (int)9
#define ParObject_Shell8Register        (int)10
#define ParObject_ElementSurface        (int)11
#define ParObject_ElementLine           (int)12
#define ParObject_ElementVolume         (int)13
#define ParObject_Shell8Line            (int)14


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

  Makes a deep copy of a Element
  
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
  
  every class imploementing ParObject needs a unique id defined at the
  top of this file.
  */
  virtual int UniqueParObjectId() const = 0;

  /*!
  \brief Pack this class so it can be communicated
  
  Allocates a char array of size size and returns it as well as size.
  Note that the first entry in the char array has to be the size itself.
  Note that the second entry in the char array has to be the type of
  instance of ParObject, also see top of this file.
  
  Pack and Unpack are used to communicate this element
  
  */
  virtual const char* Pack(int& size) const = 0;

  /*!
  \brief Unpack data from a char vector into this class
  
  \ref Pack and \ref Unpack are used to communicate this element

  */
  virtual bool Unpack(const char* data) = 0;


  /*!
  \brief return size of a char vector used for communication
  
  The size of a char vector used for communication should ALWAYS be stored
  at the beginning of the vector

  */
  inline int SizePack(const char* data) const 
  { int s = ((int*)data)[0]; return s; }
  

  //@}

protected:

  // @{ \name Routines to help pack stuff into a char string

  /*!
  \brief return packing size of a string
  
  Calculate the size it needs to pack a string
  */
  inline int SizeString(const string& data) const 
  { return (sizeof(int)+data.size()*sizeof(char));}

  /*!
  \brief return packing size of a vector
  
  Calculate the size it needs to pack a vector<...>
  
  \note Only vectors of basic types allowed!
  */
  template<typename kind> int SizeVector(const vector<kind>& data) const 
  { return (sizeof(int)+data.size()*sizeof(kind));}

  /*!
  \brief return packing size of a Epetra_SerialDenseMatrix
  
  Calculate the size it needs to pack a Epetra_SerialDenseMatrix
  
  */
  inline int SizeDenseMatrix(const Epetra_SerialDenseMatrix& data) const 
  { return (2*sizeof(int)+data.M()*data.N()*sizeof(double));}

  /*!
  \brief Add stuff to a char vector data and increment position
  
  This method is templated for all basic type int char double enum.
  \param position : place in data where to add stuff. Position will be incremented by this method.
  \param data     : char string stuff shall be added to
  \param stuff    : basic data type (float int double char) that get's added to stuff
  
  */
  template<typename kind> void AddtoPack(int& position, char* data, kind& stuff) const
  {
    int size = sizeof(kind);
    memcpy(&data[position],&stuff,size);
    position += size;
    return;
  }
  
  /*!
  \brief Add stuff to a char vector data and increment position
  
  This method adds an array to data
  \param position : place in data where to add stuff. Position will be incremented by this method.
  \param data     : char string stuff shall be added to
  \param stuff    : array of total length stuffsize (in byte)
  \param stuffsize: length of stuff in byte
  */
  template<typename kind> void AddtoPack(int& position, char* data, kind* stuff, int stuffsize) const
  {
    memcpy(&data[position],stuff,stuffsize);
    position += stuffsize;
    return;
  }

  /*!
  \brief Add a vector of a basic type to a char vector data and increment position
  
  This method adds a vector to data
  \param position : place in data where to add stuff. Position will be incremented by this method.
  \param data     : char string stuff shall be added to
  \param stuff    : vector of things of total length stuffsize
  */
  template<typename kind> void AddVectortoPack(int& position, char* data, const vector<kind>& stuff) const
  {
    int numele = stuff.size();
    AddtoPack(position,data,numele);
    AddtoPack(position,data,&stuff[0],numele*sizeof(kind));
    return;
  }

  /*!
  \brief Add an Epetra_SerialDenseMatrix to data and increment position
  
  This method adds a Epetra_SerialDenseMatrix to data
  \param position : place in data where to add stuff. Position will be incremented by this method.
  \param data     : char string stuff shall be added to
  \param stuff    : matrix to add
  */
  void AddMatrixtoPack(int& position, char* data, 
                       const Epetra_SerialDenseMatrix& stuff) const
  {
    int m = stuff.M();
    int n = stuff.N();
    AddtoPack(position,data,m);
    AddtoPack(position,data,n);
    double* A = stuff.A();
    AddtoPack(position,data,A,n*m*sizeof(double));
    return;
  }

  /*!
  \brief Add a string to a char vector data and increment position
  
  This method adds a string to data
  \param position : place in data where to add stuff. Position will be incremented by this method.
  \param data     : char string stuff shall be added to
  \param stuff    : string to add
  */
  void AddStringtoPack(int& position, char* data, const string& stuff) const
  {
    int numele = stuff.size();
    AddtoPack(position,data,numele);
    AddtoPack(position,data,&stuff[0],numele*sizeof(char));
    return;
  }

  /*!
  \brief Extract stuff from a char vector data and increment position
  
  This method is templated for all basic type like int char double enum.
  \param position : place in data where to extract stuff. Position will be incremented by this method by sizeof(kind).
  \param data     : char string where stuff is extracted from
  \param stuff    : basic data type (float int double char ...) to extract from data
  
  */
  template<typename kind> void ExtractfromPack(int& position, const char* data, kind& stuff) const
  {
    int size = sizeof(kind);
    memcpy(&stuff,&data[position],size);
    position += size;
    return;
  }

  /*!
  \brief Extract stuff from a char vector data and increment position
  
  This method extracts an array from data
  \param position : place in data where to extract stuff. Position will be incremented by this method by stuffsize.
  \param data     : char string where stuff is extracted from
  \param stuff    : array of total length stuffsize (in byte)
  \param stuffsize: length of stuff in byte
  */
  template<typename kind> void ExtractfromPack(int& position, const char* data, kind* stuff, int stuffsize) const
  {
    memcpy(stuff,&data[position],stuffsize);
    position += stuffsize;
    return;
  }

  /*!
  \brief Extract a vector of a basic type from a char vector data and increment position
  
  This method extracts an array from data
  \param position : place in data where to extract stuff. Position will be incremented by this method by stuffsize.
  \param data     : char string where stuff is extracted from
  \param stuff    : array of total length stuff.size() (in byte)
  */
  template<typename kind> void ExtractVectorfromPack(int& position, const char* data, vector<kind>& stuff) const
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
  
  \param position : place in data where to extract stuff. Position will be incremented by this method by stuffsize.
  \param data     : char string where stuff is extracted from
  \param stuff    : Epetra_SerialDenseMatrix to be extracted
  */
  void ExtractMatrixfromPack(int& position, const char* data, 
                             Epetra_SerialDenseMatrix& stuff) const
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
  \brief Extract a vector of a basic type from a char vector data and increment position
  
  This method extracts an array from data
  \param position : place in data where to extract stuff. Position will be incremented by this method by stuffsize.
  \param data     : char string where stuff is extracted from
  \param stuff    : array of total length stuffsize (in byte)
  */
  void ExtractStringfromPack(int& position, const char* data, string& stuff) const
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
