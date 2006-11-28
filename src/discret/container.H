/*!----------------------------------------------------------------------
\file container.H

\class DRT::Container

\brief A data storage container

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE
#ifndef CONTAINER_H
#define CONTAINER_H


#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_ParameterList.hpp"

using namespace std;
using namespace Teuchos;

#include "parobject.H"
#include "dserror.H"

/*!
\brief DRT: namespace of the ccarat discretization module

*/
namespace DRT
{

// forward declarations

/*!
\class Container

\brief A data storage container

\author gee (gee@lnm.mw.tum.de)
*/
class Container : public DRT::ParObject
{
public:


  // @{ \name Constructors and destructors

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
  \brief Pack this class so it can be communicated
  
  \ref Pack and \ref Unpack are used to communicate this class
  
  */
  virtual const char* Pack(int& size) const;

  /*!
  \brief Unpack data from a char vector into this class
  
  \ref Pack and \ref Unpack are used to communicate this class

  */
  virtual bool Unpack(const char* data);

  //@}

  // @{ \name Acess methods

  /*!
  \brief Print this element
  */
  virtual void Print(ostream& os) const ;

  //@}

  // @{ \name Construction

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
  \brief Get vector of type int or double
  
  Receive <int> or <double> data of a given name.
  
  \param name : Name of data to receive

  \note Usage is<br> 
        const vector<int>* ifool = Container::Get<int>("ifool_name");<br>
        or <br> 
        const vector<double>* dfool = Container::Get<double>("dfool_name");
  
  \return const ptr to vector if data with that name exists, NULL otherwise
  */
  template<typename T> const vector<T>* Get(const string& name) const
  {
    T tester;
    if (IsInt(tester))
    {
      map<string,RefCountPtr<vector<int> > >::const_iterator icurr = intdata_.find(name);
      if (icurr != intdata_.end())
        return (const vector<T>*)icurr->second.get();
      else return NULL;
    }
    else if (IsDouble(tester))
    {
      map<string,RefCountPtr<vector<double> > >::const_iterator dcurr = doubledata_.find(name);
      if (dcurr != doubledata_.end())
        return (const vector<T>*)dcurr->second.get();
      else return NULL;
    }
    else dserror("This templated method returns <int> and <double> data only!");
    return NULL;
  }
  
  /*!
  \brief Get a string from the container of type int or double
  
  \param name : Name of data to receive

  \note Usage is<br> 
        const string* ifool = Container::GetString("ifool_name");
  
  \return const ptr to string if data with that name exists, NULL otherwise
  */
  const string* GetString(const string& name) const;

  /*!
  \brief Delete a record from the container
  
  \param name : Name of data

  */
  void Delete(const string& name);

  /*!
  \brief Delete everything in this container
  */
  void Clear() { intdata_.clear(); doubledata_.clear(); return;}

  //@}

protected:

  // don't want = operator
  Container operator = (const Container& old);
  
  /*!
  \brief Detection of integer type
  */
  template<typename T> bool IsInt(const T& data) const { return false; }
  /*!
  \brief Detection of integer type
  */
  bool IsInt(const int& data) const { return true; }
  
  /*!
  \brief Detection of double type
  */
  template<typename T> bool IsDouble(const T& data) const { return false; }
  /*!
  \brief Detection of double type
  */
  bool IsDouble(const double& data) const { return true; }
    
  /*!
  \brief Detection of double type
  */
  template<typename T> bool IsString(const T& data) const { return false; }
  /*!
  \brief Detection of double type
  */
  bool IsString(const string& data) const { return true; }

protected:

  map<string,RefCountPtr<vector<int> > >    intdata_;      // a map to store integer data in
  map<string,RefCountPtr<vector<double> > > doubledata_;   // a map to store double data in
  map<string,string>                        stringdata_;   // a map to store string data in

}; // class Container
} // namespace DRT


// << operator
ostream& operator << (ostream& os, const DRT::Container& node);


#endif  // #ifndef CONTAINER_H
#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
