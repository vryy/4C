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
#include "spacetime_boundary.H"
#include "../drt_geometry/position_array.H"

namespace XFEM
{

  /// This class describes the interface using an extra discretization with surface element
  class InterfaceHandleXFSI : public InterfaceHandle
  {
  public:
    
    //! constructor
    explicit InterfaceHandleXFSI(
             const Teuchos::RCP<DRT::Discretization>&  xfemdis,   ///< discretization that will be intersected
             const Teuchos::RCP<DRT::Discretization>&  cutterdis  ///< discretization with boundary elements
             );
    
    //! destructor
    ~InterfaceHandleXFSI();

    //! return pointer to the intersecting discretization
    const Teuchos::RCP<DRT::Discretization> cutterdis() const {return cutterdis_;};
    
    //! get pointer to boundary element
    const DRT::Element* GetBoundaryEle(
            const int element_gid   ///< unique global element id of a boundary element
            ) const
    {
        return cutterdis_->gElement(element_gid);
    }
    
    template<class M>
    void fillBoundaryNodalPositionsNP(
        const DRT::Element*                   ele,                ///< element with nodal pointers
        M&                                    xyze                ///< position array (3,numnode) to be filled
        ) const
    {
      GEO::fillCurrentNodalPositions(ele, cutterposnp_, xyze);
    }
    
    template<class M>
    void fillBoundaryNodalPositionsN(
        const DRT::Element*                   ele,                ///< element with nodal pointers
        M&                                    xyze                ///< position array (3,numnode) to be filled
        ) const
    {
      GEO::fillCurrentNodalPositions(ele, cutterposn_, xyze);
    }
    
    //! get domain information at x_in at timestep n+1
    int PositionWithinConditionNP(const LINALG::Matrix<3,1>&    x_in) const;

    //! get domain information at x_in at timestep n
    int PositionWithinConditionN(const LINALG::Matrix<3,1>&     x_in) const;
    
    //! get domain information at x_in at timestep n+1
    int PositionWithinConditionNP(
        const LINALG::Matrix<3,1>&        x_in,
        GEO::NearestObject&               nearestobject
    ) const;

    //! get domain information at x_in at timestep n
    int PositionWithinConditionN(
        const LINALG::Matrix<3,1>&        x_in,
        GEO::NearestObject&               nearestobject
    ) const;

    //! find space time layer cell
    bool FindSpaceTimeLayerCell(
        const LINALG::Matrix<3,1>&        querypos,
        XFEM::SpaceTimeBoundaryCell&      slab,
        LINALG::Matrix<3,1>&              rst
    ) const;
    
    //! test domain integration cells
    void TestDomainIntCells() const;

    //! print statistics (Debug)
    void PrintStatistics() const;
    
    //! return string representation of the interface handle
    void toGmsh(const int step) const;

    
    
  private:

    //! don't want default constructor
    explicit InterfaceHandleXFSI();
    
    //! don't want copy constructor
    explicit InterfaceHandleXFSI(const InterfaceHandleXFSI&);
    
    //! disabled assignment operator
    InterfaceHandleXFSI operator = (const InterfaceHandleXFSI& old);    
    
    //! put a position into the map for each nodal id
    void FillCurrentCutterPositionMap(
        const Teuchos::RCP<DRT::Discretization>& cutterdis,              ///< dis containing nodes
        const Epetra_Vector&                     idispcol,               ///< vector with displacements
        std::map<int,LINALG::Matrix<3,1> >&      currentcutterpositions  ///< map to be filled with positions
        ) const;
        
    //! for difficult cases during intersection, do sanity checks and create list of redundant integration cells 
    void ClassifyIntegrationCells();

    //! from surface element, create 3d volume space-time cells
    void GenerateSpaceTimeLayer(
        const Teuchos::RCP<DRT::Discretization>&        cutterdis,   ///< discretization with boundary elements
        const std::map<int,LINALG::Matrix<3,1> >&       cutterposnp, ///< nodal displacements at n+1
        const std::map<int,LINALG::Matrix<3,1> >&       cutterposn   ///< nodal displacements at n
    );
    
    
    //! discretization that intersects the xfemdis_ -> should be a boundary discretization (e.g. surface of structure)
    const Teuchos::RCP<DRT::Discretization>       cutterdis_;
    
    //! bounding boxes of cutter elements           
    std::map<int,LINALG::Matrix<3,2> >            currentXAABBs_;
    
    //! map between global node id and its current positions at n+1  
    std::map<int,LINALG::Matrix<3,1> >            cutterposnp_;
    
    //! map between global node id and its current positions at n
    std::map<int,LINALG::Matrix<3,1> >            cutterposn_;

    //! space time bounday cell layer
    std::map<int,XFEM::SpaceTimeBoundaryCell>     stlayer_;

  }; // class InterfacehandleXFSI;

} // namespace XFEM

#endif // INTERFACEXFSI_H
#endif // CCADISCRET
