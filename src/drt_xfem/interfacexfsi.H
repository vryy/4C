/*!
 \file interfacexfsi.H

 \brief interface handle that transports the intersection related things around for XFSI problems

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
 */

#ifdef CCADISCRET
#ifndef INTERFACEXFSI_H
#define INTERFACEXFSI_H

#include "interface.H"

#include "../drt_geometry/intersection.H"
#include "spacetime_boundary.H"

namespace XFEM
{

  class DRT::Discretization;

  class InterfaceHandleXFSI : public InterfaceHandle
  {
  public:
    //! constructor
    explicit InterfaceHandleXFSI(
             const Teuchos::RCP<DRT::Discretization>  xfemdis,   ///< discretization that will be intersected
             const Teuchos::RCP<DRT::Discretization>  cutterdis, ///< discretization with boundary elements
             const int step                                      ///< counter to control file numbering of gmsh debug output
             );
    //! destructor
    ~InterfaceHandleXFSI();

    //! return string representation of the interface handle
    void toGmsh(const int step) const;

    //! return pointer to the intersecting discretization
    Teuchos::RCP<DRT::Discretization> cutterdis() const {return cutterdis_;};
    
    //! returns pointer to map with (current) nodal positions per node of the cutter discretization 
    const std::map<int,BlitzVec3>* cutterposnp() const {return &cutterposnp_;};
    
    //! returns pointer to map with (current) nodal positions per node of the cutter discretization 
    const std::map<int,BlitzVec3>* cutterposn()  const {return &cutterposn_;};

    /*!
     * get pointer to boundary element
     * \return pointer to boundary element
     */
    const DRT::Element* GetBoundaryEle(
            const int element_gid   ///< unique global element id of a boundary element
            ) const
    {
        return cutterdis_->gElement(element_gid);
    }

    //! get domain information at x_in at timestep n+1
    int PositionWithinConditionNP(
        const BlitzVec3&                  x_in
    ) const;

    //! get domain information at x_in at timestep n
    int PositionWithinConditionN(
        const BlitzVec3&                  x_in
    ) const;
    
    //! get domain information at x_in at timestep n+1
    int PositionWithinConditionNP(
        const BlitzVec3&                  x_in,
        GEO::NearestObject&               nearestobject
    ) const;

    //! get domain information at x_in at timestep n
    int PositionWithinConditionN(
        const BlitzVec3&                  x_in,
        GEO::NearestObject&               nearestobject
    ) const;


    
    bool FindSpaceTimeLayerCell(
        const BlitzVec3&                  querypos,
        XFEM::SpaceTimeBoundaryCell&      slab,
        BlitzVec3&                        rst
    ) const;

  private:

    //! don't want default constructor
    explicit InterfaceHandleXFSI();
    
    //! don't want copy constructor
    explicit InterfaceHandleXFSI(const InterfaceHandleXFSI&);
    
    //! put a position into the map for each nodal id
    void FillCurrentCutterPositionMap(
        const Teuchos::RCP<DRT::Discretization>  cutterdis,              ///< dis containing nodes
        const Epetra_Vector&                     idispcol,               ///< vector with displacements
        std::map<int,BlitzVec3>&                 currentcutterpositions  ///< map to be filled with positions
        ) const;

    //! for difficult cases during intersection, do sanity checks and create list of redundant integration cells 
    std::set<int> FindDoubleCountedIntersectedElements() const;

    //! from surface element, create 3d volume space-time cells
    void GenerateSpaceTimeLayer(
        const Teuchos::RCP<DRT::Discretization>  cutterdis,   ///< discretization with boundary elements
        const std::map<int,BlitzVec3>&           cutterposnp, ///< nodal displacements at n+1
        const std::map<int,BlitzVec3>&           cutterposn   ///< nodal displacements at n
    );

    //! discretization that intersects the xfemdis_ -> should be a boundary discretization (e.g. surface of structure)
    const Teuchos::RCP<DRT::Discretization> cutterdis_;
    
    //! current positions of all cutter nodes at n+1
    std::map<int,BlitzVec3>                 cutterposnp_;
    
    //! current positions of all cutter nodes at n
    std::map<int,BlitzVec3>                 cutterposn_;

    std::map<int,XFEM::SpaceTimeBoundaryCell>     stlayer_;

  }; // class InterfacehandleXFSI;

} // namespace XFEM

#endif // INTERFACEXFSI_H
#endif // CCADISCRET
