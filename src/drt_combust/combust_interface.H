/*!-----------------------------------------------------------------------------------------------*
 \file combust_interface.H

 \brief interface handle that transports the intersection related things around for combustion problems

\level 2

<pre>
\maintainer Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
 *------------------------------------------------------------------------------------------------*/

#ifndef COMBUST_INTERFACE_H
#define COMBUST_INTERFACE_H


#include <Teuchos_RCP.hpp>
#include "../drt_geometry/geo_utils.H"


namespace DRT
{
  class Discretization;
  class Element;
}

namespace COMBUST
{
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
             const Teuchos::RCP<DRT::Discretization> fluiddis,      ///< discretization (combust3 elements) which will be intersected
             const Teuchos::RCP<DRT::Discretization> gfuncdis ///< discretization (condif3 elements) with level set function (G-function)
             );

    /// destructor
    ~InterfaceHandleCombust();

    enum CutStatus
    {
      undefined, ///< element has no idea what happened to it
      uncut,     ///< element is not intersected
      touched,   ///< element in plus domain is touched by the interface on a 2D face (surface)
      bisected,  ///< element is intersected (two subdomains)
      trisected, ///< element is double intersected (three subdomains)
    };

    //! return pointer to the intersecting discretization (does anyone need this?)
    const Teuchos::RCP<const DRT::Discretization>& GfuncDis() const { return gfuncdis_; }

    //! return pointer to the intersected discretization
    Teuchos::RCP<DRT::Discretization> FluidDis() const {return fluiddis_;}

    //! update interfacehandle according to current flame front
    void UpdateInterfaceHandle(
        std::map<int,GEO::DomainIntCells >&                         elementalDomainIntCells,
        std::map<int,GEO::BoundaryIntCells >&                       elementalBoundaryIntCells,
        std::map<int,COMBUST::InterfaceHandleCombust::CutStatus >&  elementcutstatus
    );

    //! compute volume of domain minus
    double ComputeVolumeMinus();

    //! compute surface of interface
    double ComputeSurface();

    //! output of domain and boundary integration cells to Gmsh!
    void toGmsh(const int step) const;

    //! tell me, whether the element has a whole touched face and lies in the plus domain or not
    bool ElementTouched(const int xfemeleid) const;

    //! tell me, whether the element is bisected or not
    bool ElementBisected(const int xfemeleid) const;

    //! tell me, whether the element is trisected or not
    bool ElementTrisected(const int xfemeleid) const;

    //! return std::string representation of the interface handle
    std::string toString() const;

    //! tell me, whether the element is intersected or not
    bool ElementIntersected(const int element_gid) const;

    //! tell me, whether the element is bisected or not
    bool ElementSplit(const int gid) const;


    //! return map of elemental domain integration cells
    std::map<int, GEO::DomainIntCells> DomainIntCells () const;

    //! return map of elemental boundary integration cells
    std::map<int, GEO::BoundaryIntCells> BoundaryIntCells () const;

    //! return map of elemental cut status
    std::map<int, CutStatus> CutState () const;


    //! return list of domain integrationcells for a given global element id, returned list can be empty -> always terminates
    GEO::DomainIntCells ElementDomainIntCells( const int gid) const;

    //! return list of boundary integrationcells for a given global element id, returned list can be empty -> always terminates
    GEO::BoundaryIntCells ElementBoundaryIntCells(const int   gid) const;

    //! tell me, whether the element is trisected or not
    CutStatus ElementCutStatus(const int gid) const;


    //! return number of domain integration cells for a given element
    std::size_t NumDomainIntCells( const int gid) const;

    //! return number of boundary integration cells for a given element
    std::size_t NumBoundaryIntCells( const int gid) const;

  private:

    /// don't want default constructor, therefore declared private!
    explicit InterfaceHandleCombust();

    /// don't want copy constructor, therefore declared private!
    explicit InterfaceHandleCombust(const InterfaceHandleCombust&);

    //! intersected domain
    const Teuchos::RCP<DRT::Discretization>   fluiddis_;

    /// pointer to level set (G-function) discretization that intersects the fluiddis_
    const Teuchos::RCP<const DRT::Discretization> gfuncdis_;

    //! integration cells for domain (global id for element -> list of volume integration cells)
    std::map<int,GEO::DomainIntCells >        elementalDomainIntCells_;

    //! integration cells for boundary
    std::map<int,GEO::BoundaryIntCells >      elementalBoundaryIntCells_;

    /// map holding status of elements with respect to the interface
    std::map<int,CutStatus > elementcutstatus_;


  }; // class InterfaceHandleCombust
} // namespace COMBUST

#endif // COMBUST_INTERFACE_H
