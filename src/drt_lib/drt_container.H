/*!----------------------------------------------------------------------
\file drt_container.H

<pre>
-------------------------------------------------------------------------
                 BACI finite element library subsystem
            Copyright (2008) Technical University of Munich

Under terms of contract T004.008.000 there is a non-exclusive license for use
of this work by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library is proprietary software. It must not be published, distributed,
copied or altered in any form or any media without written permission
of the copyright holder. It may be used under terms and conditions of the
above mentioned license by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library may solemnly used in conjunction with the BACI contact library
for purposes described in the above mentioned contract.

This library contains and makes use of software copyrighted by Sandia Corporation
and distributed under LGPL licence. Licensing does not apply to this or any
other third party software used here.

Questions? Contact Dr. Michael W. Gee (gee@lnm.mw.tum.de)
                   or
                   Prof. Dr. Wolfgang A. Wall (wall@lnm.mw.tum.de)

http://www.lnm.mw.tum.de

-------------------------------------------------------------------------
</pre>

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifndef CONTAINER_H
#define CONTAINER_H


#include <string>
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Epetra_SerialDenseMatrix.h"

#include "drt_parobject.H"
#include "drt_dserror.H"
#include "drt_parobjectfactory.H"

class Epetra_MultiVector;
class Epetra_Vector;

namespace DRT
{

class ContainerType : public DRT::ParObjectType
{
public:

  std::string Name() const { return "ContainerType"; }

  static ContainerType & Instance() { return instance_; };

  virtual DRT::ParObject* Create( const std::vector<char> & data );

private:

  static ContainerType instance_;
};

// forward declarations

/*!
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
  virtual int UniqueParObjectId() const { return ContainerType::Instance().UniqueParObjectId(); }

  /*!
  \brief Pack this class so it can be communicated

  \ref Pack and \ref Unpack are used to communicate this class

  */
  virtual void Pack(DRT::PackBuffer& data) const;

  /*!
  \brief Unpack data from a char vector into this class

  \ref Pack and \ref Unpack are used to communicate this class

  */
  virtual void Unpack(const std::vector<char>& data);

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
  void Add(const std::string& name, const std::vector<int>& data)
  { Add(name,&data[0],data.size()); return;}

  /*!
  \brief Add Teuchos::rcp to vector of int to the container

  \param name : Name of data to store data with
  \param data : Teuchos::rcp to vector of integers to be added

  \note If record with name name already exists, it will be overwritten
  */
  void Add(const std::string& name, Teuchos::RCP<std::vector<int> > data);

  /*!
  \brief Add a double to the container

  \param name : Name of data to store data with
  \param data : one double to be added to the container

  \note If record with name name already exists, it will be overwritten
  */
  void Add(const std::string& name, const double& data)
  { Add(name,&data,1); return; }

  /*!
  \brief Add vector of double to the container

  \param name : Name of data to store data with
  \param data : ptr to beginning of data
  \param num  : Number of objects in data

  \note If record with name name already exists, it will be overwritten
  */
  void Add(const std::string& name, const double* data, const int num);

  /*!
  \brief Add vector of double to the container

  \param name : Name of data to store data with
  \param data : vector of doubles to be added

  \note If record with name name already exists, it will be overwritten
  */
  void Add(const std::string& name, const std::vector<double>& data)
  { Add(name,&data[0],data.size()); return;}

   /*!
  \brief Add Teuchos::rcp to vector of double to the container

  \param name : Name of data to store data with
  \param data : Teuchos::rcp to vector of doubles to be added

  \note If record with name name already exists, it will be overwritten
  */
  void Add(const std::string& name, Teuchos::RCP<std::vector<double> > data);

  /*!
  \brief Add a string to the container

  \param name : Name of data to store data with
  \param data : string to store

  \note If record with name name already exists, it will be overwritten
  */
  void Add(const std::string& name, const std::string& data);

  /*!
  \brief Add a string to the container

  \param name : Name of data to store data with
  \param data : string to store

  \note If record with name name already exists, it will be overwritten
  */
  void Add(const std::string& name, const char data[])
  { Add(name,(std::string)data); return; }

  /*!
  \brief Add a Epetra_SerialDenseMatrix to the container

  \param name : Name of data to store data with
  \param data : matrix to store

  \note If record with name name already exists, it will be overwritten
  */
  void Add(const std::string& name, const Epetra_SerialDenseMatrix& matrix);

   /*!
  \brief Add Teuchos::rcp to Epetra_SerialDenseMatrix to the container

  \param name : Name of data to store data with
  \param data : Teuchos::rcp to Epetra_SerialDenseMatrix to be added

  \note If record with name name already exists, it will be overwritten
  */
  void Add(const std::string& name, Teuchos::RCP<Epetra_SerialDenseMatrix> data);

   /*!
  \brief Add Teuchos::rcp to Epetra_MultiVector to the container

  \param name : Name of data to store data with
  \param data : Teuchos::rcp to Epetra_MultiVector to be added

  \note If record with name name already exists, it will be overwritten

  \warning data stored in Epetra_MultiVector will not participate in
           Pack/Unpack functionality of the container but is deep-copied in
           the copy-constructor
  */
  void Add(const string& name, Epetra_MultiVector& data);

   /*!
  \brief Add Teuchos::rcp to Epetra_Vector to the container

  \param name : Name of data to store data with
  \param data : Teuchos::rcp to Epetra_Vector to be added

  \note If record with name name already exists, it will be overwritten

  \warning data stored in Epetra_Vector will not participate in
           Pack/Unpack functionality of the container but is deep-copied in
           the copy-constructor
  */
  void Add(const string& name, Epetra_Vector& data);

  /*!
  \brief Add a LINALG::Matrix to the container, templated on the dimension
         of the LINALG::Matrix

  \param name : Name of data to store data with
  \param data : matrix to store

  \note Internally transformed into a Epetra_SerialDenseMatrix.

  Usage is<br>
  \code
    DRT::Container->Add<6,6>("mfool_name", LINALG::Matrix<6,6>);
  \endcode

  \note Extracting only possible via Epetra_SerialDenseMatrix, which then
        can be used to set a view on the data:

  Usage is<br>
  \code
        LINALG::Matrix<6,6> xyz(*(const_cast<Epetra_SerialDenseMatrix*>(Container::
        Get<Epetra_SerialDenseMatrix>("mfool_name"))), true);
  \endcode
  */
  template<int row, int col>
  void Add(const string& name, const LINALG::Matrix<row,col>& matrix)
  {
    Epetra_SerialDenseMatrix data(row,col);
    for(int i=0; i<row; i++)
      for(int j=0; j<col; j++)
        data(i,j) = matrix(i,j);

    Add(name,data);
    return;
  };

  //@}

  //! @name Access methods

  /*!
  \brief Get data of type int, double, string or Epetra_SerialDenseMatrix, read only

  Receive <int>, <double>, <string>, <Epetra_SerialDenseMatrix> data of a given name.

  \param name (in): Name of data to receive

  Usage is<br>
    \code
    const std::vector<int>* ifool = Container::Get<std::vector<int> >("ifool_name");
    \endcode
    or <br>
    \code
    const std::vector<double>* dfool = Container::Get<std::vector<double> >("dfool_name");
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
  \brief Get a single integer back without all that hassle

  Receive <int>

  \param name (in): Name of data to receive

  Usage is<br>
  \code
  int mynumber = Container::Getint("name_of_my_number")
  \endcode
  */
  int GetInt(const string& name) const;

  /*!
  \brief Get a single double back without all that hassle

  Receive <int>

  \param name (in): Name of data to receive

  Usage is<br>
  \code
  double mynumber = Container::GetDouble("name_of_my_number")
  \endcode
  */
  double GetDouble(const string& name) const;

  /*!
  \brief Get data of type int, double, string or Epetra_SerialDenseMatrix

  Receive <int>, <double>, <string>, <Epetra_SerialDenseMatrix> data of a given name.

  \param name (in): Name of data to receive

  Usage is<br>
    \code
    std::vector<int>* ifool = Container::Get<std::vector<int> >("ifool_name");
    \endcode
    or <br>
    \code
    std::vector<double>* dfool = Container::Get<std::vector<double> >("dfool_name");
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
  {
    intdata_.clear(); doubledata_.clear(); stringdata_.clear(); matdata_.clear();
    return;
  }

  //@}

  //! @name deep copy access methods
  /// to extract all data

  std::pair<std::map<std::string,Teuchos::RCP<std::vector<int> > >::const_iterator,
            std::map<std::string,Teuchos::RCP<std::vector<int> > >::const_iterator>
  IntRange()
  { return std::make_pair( intdata_.begin(), intdata_.end() ); }

  std::pair<std::map<std::string,Teuchos::RCP<std::vector<double> > >::const_iterator,
            std::map<std::string,Teuchos::RCP<std::vector<double> > >::const_iterator>
  DoubleRange()
  { return std::make_pair( doubledata_.begin(), doubledata_.end() ); }

  std::pair<std::map<std::string,std::string>::const_iterator,
            std::map<std::string,std::string>::const_iterator>
  StringRange()
  { return std::make_pair( stringdata_.begin(), stringdata_.end() ); }

  std::pair<std::map<std::string,Teuchos::RCP<Epetra_SerialDenseMatrix> >::const_iterator,
            std::map<std::string,Teuchos::RCP<Epetra_SerialDenseMatrix> >::const_iterator>
  MatRange()
  { return std::make_pair( matdata_.begin(), matdata_.end() ); }

  std::pair<std::map<std::string,Teuchos::RCP<Epetra_MultiVector> >::const_iterator,
            std::map<std::string,Teuchos::RCP<Epetra_MultiVector> >::const_iterator>
  EvecRange()
  { return std::make_pair( evecdata_.begin(), evecdata_.end() ); }

  //@

private:  // these methods are for PURE internal use

  // don't want = operator
  Container operator = (const Container& old);

protected:

  //! a map to store integer data in
  std::map<string,Teuchos::RCP<std::vector<int> > >             intdata_;

  //! a map to store double data in
  std::map<string,Teuchos::RCP<std::vector<double> > >          doubledata_;

  //! a map to store string data in
  std::map<string,string>                         stringdata_;

  //! a map to store matrices in
  std::map<string,Teuchos::RCP<Epetra_SerialDenseMatrix> > matdata_;

  //! a map to store epetra vectors in
  std::map<string,Teuchos::RCP<Epetra_MultiVector> >       evecdata_;

}; // class Container
} // namespace DRT

// << operator
ostream& operator << (ostream& os, const DRT::Container& node);



#endif  // #ifndef CONTAINER_H
