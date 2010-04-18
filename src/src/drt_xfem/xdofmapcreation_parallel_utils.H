/*!
\file xdofmapcreation_parallel_utils.H

\brief defines unknowns based on the intersection pattern from the xfem intersection


<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
*/
#ifdef CCADISCRET
#ifndef XDOFMAPCREATION_PARALLEL_UTILS_H
#define XDOFMAPCREATION_PARALLEL_UTILS_H

#include "../drt_xfem/dof_management.H"
#include "../drt_xfem/interface.H"
#include "../drt_xfem/element_ansatz.H"

namespace XFEM
{

  /*!
   * \brief create dofkeyset from nodal dofmap
   */
  void fillNodalDofKeySet(
      const XFEM::InterfaceHandle& ih,
      const std::map<int, std::set<XFEM::FieldEnr> >&  nodalDofSet,
      std::set<XFEM::DofKey<XFEM::onNode> >&      nodaldofkeyset
  );


  void updateNodalDofMap(
      const XFEM::InterfaceHandle& ih,
      std::map<int, std::set<XFEM::FieldEnr> >&      nodalDofSet,
      const std::set<XFEM::DofKey<XFEM::onNode> >&   nodaldofkeyset
      );


  /*!
   * \brief pack dofkeys from set into char array
   */
  void packDofKeys(
      const set<XFEM::DofKey<XFEM::onNode> >&     dofkeyset,
      vector<char>&                               dataSend
      );


  /*!
   * \brief unpack dofkeys from char array
   */
  void unpackDofKeys(
      const vector<char>&                     dataRecv,
      set<XFEM::DofKey<XFEM::onNode> >&       dofkeyset
      );

#ifdef PARALLEL
  /*!
   * \brief parallel syncronization
   */
  void syncNodalDofs(
      const XFEM::InterfaceHandle& ih,
      std::map<int, std::set<XFEM::FieldEnr> >&  nodalDofSet
      );
#endif

} // namespace XFEM

#endif  // #ifndef XDOFMAPCREATION_H
#endif  // #ifdef CCADISCRET
