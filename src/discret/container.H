/*!----------------------------------------------------------------------
\file container.H

\class CCADISCRETIZATION::Container

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

/*!
\brief CCADISCRET: namespace of the ccarat discretization module

*/
namespace CCADISCRETIZATION
{

// forward declarations

/*!
\class Container

\brief A data storage container

\author gee (gee@lnm.mw.tum.de)
*/
class Container : public CCADISCRETIZATION::ParObject
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
  Container(const CCADISCRETIZATION::Container& old);

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
  */
  void Add(const string& name, const int* data, const int num);

  /*!
  \brief Add vector of double to the container
  
  \param name : Name of data to store data with
  \param data : ptr to beginning of data
  \param num  : Number of objects in data
  */
  void Add(const string& name, const double* data, const int num);

  //@}

protected:

  // don't want = operator
  Container operator = (const Container& old);
  
  map<string,RefCountPtr<vector<int> > >    intdata_;      // a map to store integer data in
  map<string,RefCountPtr<vector<double> > > doubledata_;   // a map to store double data in

}; // class Container
} // namespace CCADISCRETIZATION


// << operator
ostream& operator << (ostream& os, const CCADISCRETIZATION::Container& node);


#endif  // #ifndef CONTAINER_H
#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
