/*!----------------------------------------------------------------------
\file drt_container.H

\class DRT::Container

\brief A data storage container

You can store arrays of integer, double, and strings in this container
and access the data by keys. As the Container implements ParObject
everything stored in here can easily be communicated using the DRT::Exporter.

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef CONTAINER_H
#define CONTAINER_H


#include <string>
#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Epetra_SerialDenseMatrix.h"

using namespace std;
using namespace Teuchos;

#include "drt_parobject.H"
#include "drt_dserror.H"

/*!
\brief DRT: namespace of the ccarat discretization module

*/
namespace DRT
{

// forward declarations

/*!
\class Container

\brief A data storage container

You can store arrays of integer, double, and strings in this container
and access the data by keys. As the Container implements ParObject
everything stored in here can easily be communicated using the DRT::Exporter.

The intention of this class is to store rather 'small' units of data. Though possible,
it is not meant to be used at the system level to store huge data sets such as sparse matrices or
vectors of system length. It does therefore not support any Epetra_Vector
or Epetra_CrsMatrix objects and is not supposed to in the future either.

\author gee (gee@lnm.mw.tum.de)
*/
class Container : public DRT::ParObject
{
public:


  //! @name Constructors and destructors and related methods

  /*!
  \brief Standard Constructor

  */
  Container();

  /*!
  \brief Copy Constructor

  Makes a deep copy of a Node

  */
  Container(const DRT::Container& old);

  /*!
  \brief Destructor

  */
  virtual ~Container();

  /*!
  \brief Return unique ParObject id

  every class imploementing ParObject needs a unique id defined at the
  top of this file.
  */
  virtual int UniqueParObjectId() const { return ParObject_Container;}

  /*!
  \brief Pack this class so it can be communicated

  \ref Pack and \ref Unpack are used to communicate this class

  */
  virtual void Pack(vector<char>& data) const;

  /*!
  \brief Unpack data from a char vector into this class

  \ref Pack and \ref Unpack are used to communicate this class

  */
  virtual void Unpack(const vector<char>& data);

  //@}

  //! @name Query methods

  /*!
  \brief Print this element
  */
  virtual void Print(ostream& os) const ;

  //@}

  //! @name Construction methods

  /*!
  \brief Add an integer to the container

  \param name : Name of data to store data with
  \param data : the integer to add

  \note If record with name name already exists, it will be overwritten
  */
  void Add(const string& name, const int& data)
  { Add(name,&data,1); return; }

  /*!
  \brief Add vector of int to the container

  \param name : Name of data to store data with
  \param data : ptr to beginning of data
  \param num  : Number of objects in data

  \note If record with name name already exists, it will be overwritten
  */
  void Add(const string& name, const int* data, const int num);

  /*!
  \brief Add vector of int to the container

  \param name : Name of data to store data with
  \param data : vector of integers to be added

  \note If record with name name already exists, it will be overwritten
  */
  void Add(const string& name, const vector<int>& data)
  { Add(name,&data[0],data.size()); return;}

  /*!
  \brief Add a double to the container

  \param name : Name of data to store data with
  \param data : one double to be added to the container

  \note If record with name name already exists, it will be overwritten
  */
  void Add(const string& name, const double& data)
  { Add(name,&data,1); return; }

  /*!
  \brief Add vector of double to the container

  \param name : Name of data to store data with
  \param data : ptr to beginning of data
  \param num  : Number of objects in data

  \note If record with name name already exists, it will be overwritten
  */
  void Add(const string& name, const double* data, const int num);

  /*!
  \brief Add vector of double to the container

  \param name : Name of data to store data with
  \param data : vector of doubles to be added

  \note If record with name name already exists, it will be overwritten
  */
  void Add(const string& name, const vector<double>& data)
  { Add(name,&data[0],data.size()); return;}

  /*!
  \brief Add a string to the container

  \param name : Name of data to store data with
  \param data : string to store

  \note If record with name name already exists, it will be overwritten
  */
  void Add(const string& name, const string& data);

  /*!
  \brief Add a string to the container

  \param name : Name of data to store data with
  \param data : string to store

  \note If record with name name already exists, it will be overwritten
  */
  void Add(const string& name, const char data[])
  { Add(name,(string)data); return; }

  /*!
  \brief Add a Epetra_SerialDenseMatrix to the container

  \param name : Name of data to store data with
  \param data : matrix to store

  \note If record with name name already exists, it will be overwritten
  */
  void Add(const string& name, const Epetra_SerialDenseMatrix& matrix);

  //@}

  //! @name Access methods

  /*!
  \brief Get data of type int, double, string or Epetra_SerialDenseMatrix, read only

  Receive <int>, <double>, <string>, <Epetra_SerialDenseMatrix> data of a given name.

  \param name (in): Name of data to receive

  Usage is<br>
    \code
    const vector<int>* ifool = Container::Get<vector<int> >("ifool_name");
    \endcode
    or <br>
    \code
    const vector<double>* dfool = Container::Get<vector<double> >("dfool_name");
    \endcode
    or <br>
    \code
    const string* sfool = Container::Get<string>("sfool_name");
    \endcode
    or <br>
    \code
    const Epetra_SerialDenseMatrix* mfool = Container::Get<Epetra_SerialDenseMatrix>("mfool_name");
    \endcode

  \note This template will let you experience one of the most kryptic error messages
        (at link time) you have ever seen if you try to use it with other data types then
        mentioned above. The same holds for the case where you did not get the
        syntax absolutely correct.

  \note (Advanced users:) This is a template with no general definition but with four
                          specializations only. It will compile for any data type
                          but will discover at link time that it is unable to generate
                          a definition other than one of the four specializations.

  \return const Ptr to data if data with that name exists, NULL otherwise
  */
  template<typename T> const T* Get(const string& name) const;

  /*!
  \brief Get data of type int, double, string or Epetra_SerialDenseMatrix

  Receive <int>, <double>, <string>, <Epetra_SerialDenseMatrix> data of a given name.

  \param name (in): Name of data to receive

  Usage is<br>
    \code
    vector<int>* ifool = Container::Get<vector<int> >("ifool_name");
    \endcode
    or <br>
    \code
    vector<double>* dfool = Container::Get<vector<double> >("dfool_name");
    \endcode
    or <br>
    \code
    string* sfool = Container::Get<string>("sfool_name");
    \endcode
    or <br>
    \code
    Epetra_SerialDenseMatrix* mfool = Container::Get<Epetra_SerialDenseMatrix>("mfool_name");
    \endcode

  \note This template will let you experience one of the most kryptic error messages
        (at link time) you have ever seen if you try to use it with other data types then
        mentioned above. The same holds for the case where you did not get the
        syntax absolutely correct.

  \note (Advanced users:) This is a template with no general definition but with four
                          specializations only. It will compile for any data type
                          but will discover at link time that it is unable to generate
                          a definition other than one of the four specializations.

  \return non-constant Ptr to data if data with that name exists, NULL otherwise
  */
  template<typename T> T* GetMutable(const string& name);

  //@}

  //! @name Erase methods

  /*!
  \brief Delete a record from the container

  \param name : Name of data to be deleted from container

  */
  void Delete(const string& name);

  /*!
  \brief Delete everything in this container
  */
  void Clear()
  { intdata_.clear(); doubledata_.clear(); stringdata_.clear(); matdata_.clear();
    return;}

  //@}

private:  // these methods are for PURE internal use

  // don't want = operator
  Container operator = (const Container& old);

protected:

  //! a map to store integer data in
  map<string,RefCountPtr<vector<int> > >             intdata_;

  //! a map to store double data in
  map<string,RefCountPtr<vector<double> > >          doubledata_;

  //! a map to store string data in
  map<string,string>                                 stringdata_;

  //! a map to store matrices in
  map<string,RefCountPtr<Epetra_SerialDenseMatrix> > matdata_;

}; // class Container
} // namespace DRT

// << operator
ostream& operator << (ostream& os, const DRT::Container& node);



#endif  // #ifndef CONTAINER_H
#endif  // #ifdef CCADISCRET
