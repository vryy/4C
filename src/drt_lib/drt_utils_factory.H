/*!----------------------------------------------------------------------
\file drt_utils_factory.H
\brief Declaration A collection of helper methods for namespace DRT
<pre>
\brief Declaration
\level 0
\maintainer Martin Kronbichler
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>

*----------------------------------------------------------------------*/
#ifndef UTILS_FACTORY_H
#define UTILS_FACTORY_H

#include "drt_discret.H"
#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"

namespace DRT
{
/// collection of useful and important but less central stuff
namespace UTILS
{

  /*!
  \brief Declaration Create an instance of a ParObject depending on the type stored in data

  An instance of ParObject is allocated and returned. The type of instance
  depends on the first entry in data and is an int (defined at the top of
  ParObject.H)

  \param data (in): A char vector used for communication or io

  \warning All instances of ParObject have to store the parobject id
           at the very beginning of the char vector

  \return Allocates and returns the correct instance of ParObject where data is
          already unpacked in the instance. The calling method is responsible for
          freeing this instance!
  */
  DRT::ParObject* Factory(const std::vector<char>& data);

  /*!
  \brief Declaration Create an instance of a finite element depending on the type of element

  A DRT::Element derived class is allocated and returned. The type of element
  allocated depends on the input parameter eletype.

  \param eletype (in): A string containing the type of element
  \param distype (in): A string containing the distype of the element
  \param id      (in): id of the new element to be created
  \param owner   (in): owner of the new element

  */
  Teuchos::RCP<DRT::Element> Factory(const std::string eletype,
                                     const std::string distype,
                                     const int id        ,
                                     const int owner     );

  //! flag, whether surfaces or lines have to be created in the ElementBoundaryFactory
  enum BoundaryBuildType
  {
    buildSurfaces,  ///< build surfaces
    buildLines,     ///< build lines
    buildNothing    ///< build nothing
  };

  /*!
   * create new instances of volume / surface / line elements for a given parent element
   *
   * template version
   *
   * \tparam BoundaryEle class name of desired volume/surface/line element, e.g. FluidSurface
   * \tparam ParentEle   class name of parent element, e.g. Fluid
   *
   * this templated function creates all volume / surface / line elements for a given element
   * and fills handed over vectors with Teuchos::RCPs and raw pointers of these boundary elements
   * This method is a very powerful helper function for implementing the necessary
   * Surfaces() and Lines() methods for every element class (especially in 3D)
   * using the BACI conventions for element connectivity
   *
   * \return boundaryeles   vector filled with Teuchos::RCPs of allocated boundary elements
   *
   * \author gjb
   * \date 05/08
   */
  template <class BoundaryEle, class ParentEle>
  std::vector<Teuchos::RCP<DRT::Element> > ElementBoundaryFactory(
    const BoundaryBuildType  buildtype,  ///< flag, whether volumes, surfaces or lines have to be created
    ParentEle*               ele         ///< pointer on the parent element
    )
  {
    // do we have to build volume, surface or line elements?
    // get node connectivity for specific distype of parent element
    unsigned int nele=0;
    const DRT::Element::DiscretizationType distype = ele->Shape();
    std::vector< std::vector<int> > connectivity;
    switch (buildtype)
    {
      case buildSurfaces:
      {
        nele = ele->NumSurface();
        connectivity = DRT::UTILS::getEleNodeNumberingSurfaces(distype);
        break;
      }
      case buildLines:
      {
        nele = ele->NumLine();
        connectivity = DRT::UTILS::getEleNodeNumberingLines(distype);
        break;
      }
      default: dserror("buildNothing case not handled in ElementBoundaryFactory");
    }
    // create vectors that will contain the volume, surface or line elements
    std::vector<Teuchos::RCP<DRT::Element> > boundaryeles(nele);

    // does DRT::UTILS convention match your implementation of NumSurface() or NumLine()?
    if (nele != connectivity.size()) dserror("number of surfaces or lines does not match!");

    // now, build the new surface/line elements
    for (unsigned int iele = 0; iele < nele; iele++)
    {
      // allocate node vectors
      unsigned int nnode = connectivity[iele].size(); // this number changes for pyramids or wedges
      std::vector<int> nodeids(nnode);
      std::vector<DRT::Node*> nodes(nnode);

      // get connectivity infos
      for (unsigned int inode=0;inode<nnode;inode++)
      {
        nodeids[inode] = ele->PointIds()[connectivity[iele][inode]];
        nodes  [inode] = ele->Points  ()[connectivity[iele][inode]];
      }

      // allocate a new boundary element
      boundaryeles[iele] = Teuchos::rcp(new BoundaryEle(iele          ,
                                                        ele->Owner()  ,
                                                        nodeids.size(),
                                                        &nodeids[0]   ,
                                                        &nodes[0]     ,
                                                        ele           ,
                                                        iele          ));
    }

    return boundaryeles;
  }

  /*!
   * create new instances of volume / surface / line elements for a given parent element
   *
   * template version
   *
   * \tparam IntFaceEle  class name of desired volume/surface/line element, e.g. FluidSurface
   * \tparam ParentEle   class name of parent element, e.g. Fluid
   *
   * this templated function creates an internal face element for two given parent elements
   * and fills an Teuchos::RCP and raw pointers of this internal faces elements
   * This method is a very powerful helper function for implementing the necessary
   * Surfaces() and Lines() methods for every element class (especially in 3D)
   * using the BACI conventions for element connectivity
   *
   * \return intface   Teuchos::RCP of allocated internal face element
   *
   * \author schott
   * \date 03/12
   */
  template <class IntFaceEle, class ParentEle>
  Teuchos::RCP<DRT::Element> ElementIntFaceFactory(
    int id,                                  ///< element id
    int owner,                               ///< owner (= owner of parent element with smallest gid)
    int nnode,                               ///< number of nodes
    const int* nodeids,                      ///< node ids
    DRT::Node** nodes,                       ///< nodes of surface
    ParentEle*  master_ele,                  ///< pointer on the master parent element
    ParentEle*  slave_ele,                   ///< pointer on the slave parent element
    const int   lsurface_master,             ///< local surface index with respect to master parent element
    const int   lsurface_slave,              ///< local surface index with respect to slave parent element
    const std::vector<int> localtrafomap     ///< local trafo map
    )
  {

    // create a new internal face element
    return Teuchos::rcp(new IntFaceEle(id,
                                       owner,
                                       nnode,
                                       nodeids,
                                       nodes,
                                       master_ele,
                                       slave_ele,
                                       lsurface_master,
                                       lsurface_slave,
                                       localtrafomap) );
  }


} // namespace UTILS
} // namespace DRT


#endif  // #ifndef UTILS_FACTORY_H




