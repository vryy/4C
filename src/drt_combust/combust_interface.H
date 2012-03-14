/*!-----------------------------------------------------------------------------------------------*
 \file combust_interface.H

 \brief interface handle that transports the intersection related things around for combustion problems

<pre>
Maintainer: Florian Henke
            henke@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15265
</pre>
 *------------------------------------------------------------------------------------------------*/

#ifndef COMBUST_INTERFACE_H
#define COMBUST_INTERFACE_H

#include "combust_flamefront.H"


namespace COMBUST
{
  class FlameFront;

  /// COMBUST Interface class
  /*!

    detailed description

    \author henke
    \date 10/08
   */
  class InterfaceHandleCombust
  {
  public:
    /// constructor
    explicit InterfaceHandleCombust(
             const Teuchos::RCP<DRT::Discretization> fluiddis,  ///< discretization (combust3 elements) which will be intersected
             const Teuchos::RCP<const DRT::Discretization> gfuncdis,  ///< discretization (condif3 elements) with level set function (G-function)
             const Teuchos::RCP<const COMBUST::FlameFront> flamefront  ///< geometrical interface (flame front) information
             );

    /// destructor
    ~InterfaceHandleCombust();

    //! return pointer to the intersecting discretization (does anyone need this?)
    const Teuchos::RCP<const DRT::Discretization>& GfuncDis() const { return gfuncdis_; }

    //! return pointer to flame front (does anyone need this?)
    const Teuchos::RCP<const COMBUST::FlameFront>& FlameFront() const { return flamefront_; }

    //! update interfacehandle according to current flame front
    void UpdateInterfaceHandle();

    //! compute volume of domain minus
    const double ComputeVolumeMinus();

    //! compute surface of interface
    const double ComputeSurface();

    //! output of domain and boundary integration cells to Gmsh!
    void toGmsh(const int step) const;

    //! tell me, whether the element has a whole touched face and lies in the plus domain or not
    bool ElementTouched(const int xfemeleid) const;

    //! tell me, whether the element is bisected or not
    bool ElementBisected(const int xfemeleid) const;

    //! tell me, whether the element is trisected or not
    bool ElementTrisected(const int xfemeleid) const;

    //! tell me, whether the element is trisected or not
    COMBUST::FlameFront::CutStatus ElementCutStatus(const int xfemeleid) const;

    //! return string representation of the interface handle
    std::string toString() const;

    //! return pointer to the intersected discretization
    Teuchos::RCP<DRT::Discretization> xfemdis() const {return xfemdis_;}

    //! tell me, whether the element is intersected or not
    bool ElementIntersected(const int element_gid) const;

    //! tell me, whether the element is bisected or not
    bool ElementSplit(const DRT::Element* xfemElement) const;

    //! check, whether given element id has been intersected by the given interface (label)
    bool ElementHasLabel(
        const int element_gid,
        const int label
    ) const;

    //! get labels for each element
    std::set<int> LabelsPerElement( const int element_gid) const;

    //! get set of boundary element GIDs that intersect an element
    std::set<int> GetIntersectingBoundaryElementsGID(
        const int element_gid  ///< element GID of intersected element
        ) const;

    //! return list of domain integrationcells for a given global element id, returned list can be empty -> always terminates
    GEO::DomainIntCells GetDomainIntCells( const DRT::Element*  xfemElement) const;

    //! return number of domain integration cells for a given element
    std::size_t GetNumDomainIntCells( const DRT::Element*  xfemElement) const;

    //! return map of elemental boundary integration cells
    std::map<int, GEO::BoundaryIntCells> GetElementalBoundaryIntCells () const;

    //! return number of boundary integration cells for a given element
    std::size_t GetNumBoundaryIntCells( const DRT::Element*  xfemElement) const;

    //! return list of boundary integrationcells for a given global element id, returned list can be empty -> always terminates
    GEO::BoundaryIntCells GetBoundaryIntCells(const int   gid) const;

    //! get domain information at x_in at timestep n+1
    int PositionWithinConditionNP(const LINALG::Matrix<3,1>& x_in) const;

  protected:

    //! fill label per element list
    void InvertElementsPerLabel();

    //! intersected domain
    const Teuchos::RCP<DRT::Discretization>   xfemdis_;

    //! integration cells for domain (global id for element -> list of volume integration cells)
    std::map<int,GEO::DomainIntCells >        elementalDomainIntCells_;

    //! integration cells for boundary
    std::map<int,GEO::BoundaryIntCells >      elementalBoundaryIntCells_;

    //! for each xfem condition label there is a set of boundary elements
    std::map<int,std::set<int> >              boundaryElementsByLabel_;

    //! each boundary element belongs to one label (inverted information of elementsByLabel_)
    //! i.e. its a map like (ElementGID,Label)
    std::map<int,int>                         labelPerBoundaryElementId_;

  private:

    //! disabled assignment operator
    InterfaceHandleCombust operator = (const InterfaceHandleCombust& old);

    /// don't want default constructor, therefore declared private!
    explicit InterfaceHandleCombust();

    /// don't want copy constructor, therefore declared private!
    explicit InterfaceHandleCombust(const InterfaceHandleCombust&);

    /// pointer to level set (G-function) discretization that intersects the fluiddis_
    const Teuchos::RCP<const DRT::Discretization> gfuncdis_;

    /// pointer to the flame front
    const Teuchos::RCP<const COMBUST::FlameFront> flamefront_;

    /// map holding status of elements with respect to the interface
    std::map<int,COMBUST::FlameFront::CutStatus > elementcutstatus_;



  }; // class InterfaceHandleCombust
} // namespace COMBUST

#endif // COMBUST_INTERFACE_H
