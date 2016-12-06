/*!
\file xdofmapcreation_parallel_utils.H

\brief defines unknowns based on the intersection pattern from the xfem intersection

\level 2

<pre>
\maintainer Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>

\warning this combustion module related file will be deleted within the next time!!!
 *------------------------------------------------------------------------------------------------*/


#ifndef XDOFMAPCREATION_PARALLEL_UTILS_H
#define XDOFMAPCREATION_PARALLEL_UTILS_H


#include <set>
#include <map>
#include <vector>


namespace COMBUST
{
  class InterfaceHandleCombust;
}

namespace DRT
{
  class PackBuffer;
}

namespace XFEM
{
  class DofKey;
  class FieldEnr;

  /*!
   * \brief create dofkeyset from nodal dofmap
   */
  void fillNodalDofKeySet(
      const COMBUST::InterfaceHandleCombust& ih,
      const std::map<int, std::set<XFEM::FieldEnr> >&  nodalDofSet,
      std::set<XFEM::DofKey >&      nodaldofkeyset
  );


  void updateNodalDofMap(
      const COMBUST::InterfaceHandleCombust& ih,
      std::map<int, std::set<XFEM::FieldEnr> >&      nodalDofSet,
      const std::set<XFEM::DofKey >&   nodaldofkeyset
      );


  /*!
   * \brief pack dofkeys from set into char array
   */
  void packDofKeys(
      const std::set<XFEM::DofKey >&     dofkeyset,
      DRT::PackBuffer&                            dataSend
      );


  /*!
   * \brief unpack dofkeys from char array
   */
  void unpackDofKeys(
      const std::vector<char>&                     dataRecv,
      std::set<XFEM::DofKey >&       dofkeyset
      );

#ifdef PARALLEL
  /*!
   * \brief parallel syncronization
   */
  void syncNodalDofs(
      const COMBUST::InterfaceHandleCombust& ih,
      std::map<int, std::set<XFEM::FieldEnr> >&  nodalDofSet
      );
#endif

} // namespace XFEM

#endif  // #ifndef XDOFMAPCREATION_H

