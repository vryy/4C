/*!
\file element_ansatz.H

\brief defines class ElementAnsatz

\level 2

<pre>
\maintainer Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>

\warning this combustion module related file will be deleted within the next time!!!
*/


#ifndef ELEMENT_ANSATZ_H
#define ELEMENT_ANSATZ_H


#include "../drt_xfem/physics.H"
#include "../drt_lib/drt_element.H"


namespace XFEM
{
  //! Implement this abstract class to tell the dofmanagers,
  //! what stress ansatz you want to use for Distributed Lagrange multipliers
  class ElementAnsatz
  {
    public:
      //! default constructor
      explicit ElementAnsatz(){};
      //! for a given velocity/displacement interpolation, a map is produced which has the number of element unknowns per physical field
      virtual std::map<XFEM::PHYSICS::Field, DRT::Element::DiscretizationType> getElementAnsatz(DRT::Element::DiscretizationType eledistype) const = 0;
    protected:
      //! no data to delete, but this way the compiler does not warn
      virtual ~ElementAnsatz(){};
  };


} // namespace XFEM

#endif  // #ifndef XDOFMAPCREATION_H

