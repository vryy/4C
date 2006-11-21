/*!----------------------------------------------------------------------
\file utils.H
\brief A collection of helper methods for namespace CCADISCRETIZATION

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE
#ifndef UTILS_H
#define UTILS_H

#include "parobject.H"

/*!
\brief CCADISCRET: namespace of the ccarat discretization module

*/
namespace CCADISCRETIZATION
{

  /*!
  \brief Create an instance of a ParObject depending on the type stored in data
  
  An instance of ParObject is allocated and returned. The type of instance
  depends on the second entry in data and is an int (defined at the top of
  ParObject.H)
  
  \param data: A char string used for communication
  
  \warning All instances of ParObject have to store the size of data in position 0
           and the type in position 1
           
  \return Allocates and returns the correct instance of ParObject where data is
          already unpacked in the instance. The calling method is responsible for
          freeing this instance!
  */
  ParObject* Factory(const char* data);


} // namespace CCADISCRETIZATION


#endif  // #ifndef UTILS_H
#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
