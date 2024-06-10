/*----------------------------------------------------------------------*/
/*! \file

\brief integration cell classes for domain and boundary integration

--> THIS FUNCTIONALITY IS JUST USED IN COMBUST AND WILL LEAVE 4C SOON

\level 3

*----------------------------------------------------------------------*/


#ifndef FOUR_C_FEM_GEOMETRY_INTEGRATIONCELL_HPP
#define FOUR_C_FEM_GEOMETRY_INTEGRATIONCELL_HPP


#include "4C_config.hpp"

#include "4C_fem_general_element.hpp"
#include "4C_fem_general_utils_local_connectivity_matrices.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Core::Geo
{
  /*----------------------------------------------------------------------------*/
  /*!
   * \brief An integration cell is used for specialized element integration routines
   */
  class IntCell
  {
   public:
    //! Standard Constructor
    explicit IntCell(const Core::FE::CellType& distype  ///< distype of the integration cell
    );

    //! Copy Constructor
    explicit IntCell(const IntCell& old);

    //! virtual destructor
    virtual ~IntCell() = default;

    virtual IntCell& operator=(const IntCell& intcell);

    //! brief returns the shape of the integration cell
    Core::FE::CellType Shape() const { return distype_; };

    //! return boolean indicating plus domain
    bool getDomainPlus() const { return indomainplus_; }

   private:
    //! hidden default constructor
    explicit IntCell();

   protected:
    //! shape
    Core::FE::CellType distype_;

    // boolean indicating that cell belongs to "+"-part of the domain
    bool indomainplus_;

    //! get geometric center of the cell in physical coordinates
    Core::LinAlg::Matrix<3, 1> compute_physical_center_position(
        const Core::FE::CellType& distype, const Core::LinAlg::SerialDenseMatrix& xyze) const;
  };

  /*----------------------------------------------------------------------------*/
  /*!
   * \brief An boundary integration cell is used for
   *        integrating at a discontinuity caused by XFEM enrichments
   */
  class BoundaryIntCell : public IntCell
  {
   public:
    /** \brief Create a new boundary integration cell
     *
     *  \param distype                  (in) : shape of the integration cell
     *  \param surface_ele_gid          (in) : global id of the boundary element (cutter)
     *  \param xfemEleDomainCoordinates (in) : coordinates in element parameter space xsi
     *  \param eleBoundaryCoordinates   (in) : coordinates in boundary parameter space eta
     *  \param physDomainCoordinates    (in) : coordinates of the integrationcell in physical domain
     *  \param indomainplus             (in) : domain part of the integration cell
     *
     *  \author hiermeier \date 11/16 */
    static BoundaryIntCell* Create(const Core::FE::CellType& distype, const int& surface_ele_gid,
        const Core::LinAlg::SerialDenseMatrix& xfemEleDomainCoordinates,
        const Core::LinAlg::SerialDenseMatrix* eleBoundaryCoordinates,
        const Core::LinAlg::SerialDenseMatrix& physDomainCoordinates, const bool& indomainplus);

   public:
    //! Standard Constructor
    explicit BoundaryIntCell(const Core::FE::CellType& distype,  ///< shape of the integration cell
        const int surface_ele_gid,  ///< global id of the boundary element (cutter)
        const Core::LinAlg::SerialDenseMatrix&
            xfemEleDomainCoordinates,  ///< coordinates in element parameter space xsi
        const Core::LinAlg::SerialDenseMatrix&
            eleBoundaryCoordinates,  ///< coordinates in boundary parameter space eta
        const Core::LinAlg::SerialDenseMatrix&
            physDomainCoordinates  ///< coordinates of the integrationcell in physical domain
    );

    //! constructor used for combustion problems
    explicit BoundaryIntCell(const Core::FE::CellType& distype,  ///< shape of the integration cell
        const int surface_ele_gid,  ///< global id of the boundary element (cutter)
        const Core::LinAlg::SerialDenseMatrix&
            xfemEleDomainCoordinates,  ///< coordinates in element parameter space xsi
        const Core::LinAlg::SerialDenseMatrix&
            eleBoundaryCoordinates,  ///< coordinates in boundary parameter space eta
        const Core::LinAlg::SerialDenseMatrix&
            physDomainCoordinates,  ///< coordinates of the integrationcell in physical domain
        const bool indomainplus     ///< domain part of the integration cell
    );

    //! Copy Constructor
    BoundaryIntCell(const BoundaryIntCell& old);


    //! assignment operator
    virtual BoundaryIntCell& operator=(const BoundaryIntCell& boundaryintcell);

    //! returns the coordinates of the integration cell in parent element coordinates xsi
    const Core::LinAlg::SerialDenseMatrix& cell_nodal_pos_xi_domain() const
    {
      return nodalpos_xi_domain_;
    };

    //! returns the coordinates of the integration cell in boundary parent space eta
    const Core::LinAlg::SerialDenseMatrix& cell_nodal_pos_xi_boundary() const
    {
      return nodalpos_xi_boundary_;
    };

    //! returns an array with the coordinates of the integration cell in physical coordinates
    const Core::LinAlg::SerialDenseMatrix& CellNodalPosXYZ() const { return nodalpos_xyz_domain_; }

    //! returns an array with the coordinates of the integration cell in physical coordinates
    const Core::LinAlg::Matrix<3, 1>& get_physical_center_position() const { return phys_center_; }

    //! return "parent" cutter element id (global id)
    int GetSurfaceEleGid() const { return surface_ele_gid_; }

   protected:
    //! constructor for derived class only.
    BoundaryIntCell(Core::FE::CellType distype,  ///< shape of the integration cell
        const int& surface_ele_gid               ///< global id of the boundary element (cutter)
    );

    //! the boundaryIntCell should know, to which cutterElement it belongs!?
    int surface_ele_gid_;

    //! coordinates of the nodes of the integration cell in parent element coordinates xsi
    Core::LinAlg::SerialDenseMatrix nodalpos_xi_domain_;

    //! boundary coordinates of the nodes of the integration cell in boundary element coordinates
    //! eta
    Core::LinAlg::SerialDenseMatrix nodalpos_xi_boundary_;

    //! coordinates of the nodes of the integration cell in physical coordinates
    Core::LinAlg::SerialDenseMatrix nodalpos_xyz_domain_;

    //! get geometric center of the cell in physical coordinates
    Core::LinAlg::Matrix<3, 1> phys_center_;

   private:
    //! hide default constructor
    explicit BoundaryIntCell();
  };

}  // namespace Core::Geo

FOUR_C_NAMESPACE_CLOSE

#endif
