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

namespace XFEM
{

  class DRT::Discretization;

  /// This class describes the interface using an extra discretization with surface element
  class InterfaceHandleXFSI : public InterfaceHandle
  {
  public:
    //! constructor
    explicit InterfaceHandleXFSI(
             const Teuchos::RCP<DRT::Discretization>&  xfemdis,   ///< discretization that will be intersected
             const Teuchos::RCP<DRT::Discretization>&  cutterdis, ///< discretization with boundary elements
             const int step,                                      ///< counter to control file numbering of gmsh debug output
             const bool experimental_intersection                 ///< choose intersection class
             );
    //! destructor
    ~InterfaceHandleXFSI();

    //! return string representation of the interface handle
    void toGmsh(const int step) const;

    //! return pointer to the intersecting discretization
    const Teuchos::RCP<DRT::Discretization> cutterdis() const {return cutterdis_;};
    
    //! returns pointer to map with (current) nodal positions per node of the cutter discretization 
    const std::map<int,LINALG::Matrix<3,1> >& cutterposnp() const {return cutterposnp_;};
    
    //! returns pointer to map with (current) nodal positions per node of the cutter discretization 
    const std::map<int,LINALG::Matrix<3,1> >& cutterposn()  const {return cutterposn_;};
    
    //! get pointer to boundary element
    const DRT::Element* GetBoundaryEle(
            const int element_gid   ///< unique global element id of a boundary element
            ) const
    {
        return cutterdis_->gElement(element_gid);
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
    
    //! print statistics
    void PrintStatistics() const;

  private:

    //! don't want default constructor
    explicit InterfaceHandleXFSI();
    
    //! don't want copy constructor
    explicit InterfaceHandleXFSI(const InterfaceHandleXFSI&);
    
    //! put a position into the map for each nodal id
    void FillCurrentCutterPositionMap(
        const Teuchos::RCP<DRT::Discretization>& cutterdis,              ///< dis containing nodes
        const Epetra_Vector&                     idispcol,               ///< vector with displacements
        std::map<int,LINALG::Matrix<3,1> >&      currentcutterpositions  ///< map to be filled with positions
        ) const;
        
    //! for difficult cases during intersection, do sanity checks and create list of redundant integration cells 
    std::set<int> FindDoubleCountedIntersectedElements() const;

    //! from surface element, create 3d volume space-time cells
    void GenerateSpaceTimeLayer(
        const Teuchos::RCP<DRT::Discretization>&        cutterdis,   ///< discretization with boundary elements
        const std::map<int,LINALG::Matrix<3,1> >&       cutterposnp, ///< nodal displacements at n+1
        const std::map<int,LINALG::Matrix<3,1> >&       cutterposn   ///< nodal displacements at n
    );
    
    //! discretization that intersects the xfemdis_ -> should be a boundary discretization (e.g. surface of structure)
    const Teuchos::RCP<DRT::Discretization>       cutterdis_;
    
    //! map between global node id and its current positions at n+1  
    std::map<int,LINALG::Matrix<3,1> >            cutterposnp_;
    
    //! map between global node id and its current positions at n
    std::map<int,LINALG::Matrix<3,1> >            cutterposn_;

    //! space time bounday cell layer
    std::map<int,XFEM::SpaceTimeBoundaryCell>     stlayer_;

  }; // class InterfacehandleXFSI;
  
  
  /*! 
   * \brief calculates the volume of an element in initial configuration
   * 
   * \return physical volume of the element in initial configuration (without any displacement)
   */
  double ElementVolume(
      const DRT::Element&           ele            ///< the element whose physical volume we want to compute
      );

} // namespace XFEM

#endif // INTERFACEXFSI_H
#endif // CCADISCRET
