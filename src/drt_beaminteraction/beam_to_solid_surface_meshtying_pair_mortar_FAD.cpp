/*----------------------------------------------------------------------*/
/*! \file

\brief Mortar mesh tying element for between a 3D beam and a surface element, coupling terms are
evaluated with FAD.

\level 3
\maintainer Ivo Steinbrecher
*/


#include "beam_to_solid_surface_meshtying_pair_mortar_FAD.H"

#include "beaminteraction_calc_utils.H"
#include "beam_to_solid_mortar_manager.H"
#include "../drt_geometry_pair/geometry_pair_line_to_surface.H"
#include "../drt_geometry_pair/geometry_pair_element_faces.H"
#include "../drt_inpar/inpar_beam_to_solid.H"

#include "Epetra_FEVector.h"


/**
 *
 */
template <typename scalar_type, typename beam, typename surface, typename mortar>
BEAMINTERACTION::BeamToSolidSurfaceMeshtyingPairMortarFAD<scalar_type, beam, surface,
    mortar>::BeamToSolidSurfaceMeshtyingPairMortarFAD()
    : base_class()
{
  // Empty constructor.
}


/**
 *
 */
Teuchos::RCP<BEAMINTERACTION::BeamContactPair>
BEAMINTERACTION::BeamToSolidSurfaceMeshtyingPairMortarFADFactory(
    const DRT::Element::DiscretizationType surface_shape,
    const INPAR::BEAMTOSOLID::BeamToSolidMortarShapefunctions mortar_shapefunction)
{
  using namespace GEOMETRYPAIR;

  switch (mortar_shapefunction)
  {
    case INPAR::BEAMTOSOLID::BeamToSolidMortarShapefunctions::line2:
    {
      switch (surface_shape)
      {
        case DRT::Element::tri3:
          return Teuchos::rcp(
              new BeamToSolidSurfaceMeshtyingPairMortarFAD<line_to_surface_patch_scalar_type,
                  t_hermite, t_tri3, t_line2>());
        case DRT::Element::tri6:
          return Teuchos::rcp(
              new BeamToSolidSurfaceMeshtyingPairMortarFAD<line_to_surface_patch_scalar_type,
                  t_hermite, t_tri6, t_line2>());
        case DRT::Element::quad4:
          return Teuchos::rcp(
              new BeamToSolidSurfaceMeshtyingPairMortarFAD<line_to_surface_patch_scalar_type,
                  t_hermite, t_quad4, t_line2>());
        case DRT::Element::quad8:
          return Teuchos::rcp(
              new BeamToSolidSurfaceMeshtyingPairMortarFAD<line_to_surface_patch_scalar_type,
                  t_hermite, t_quad8, t_line2>());
        case DRT::Element::quad9:
          return Teuchos::rcp(
              new BeamToSolidSurfaceMeshtyingPairMortarFAD<line_to_surface_patch_scalar_type,
                  t_hermite, t_quad9, t_line2>());
        case DRT::Element::nurbs9:
          return Teuchos::rcp(new BeamToSolidSurfaceMeshtyingPairMortarFAD<
              line_to_surface_patch_nurbs_scalar_type<t_hermite, t_nurbs9>, t_hermite, t_nurbs9,
              t_line2>());
        default:
          dserror("Wrong element type for surface element.");
      }
      break;
    }
    case INPAR::BEAMTOSOLID::BeamToSolidMortarShapefunctions::line3:
    {
      switch (surface_shape)
      {
        case DRT::Element::tri3:
          return Teuchos::rcp(
              new BeamToSolidSurfaceMeshtyingPairMortarFAD<line_to_surface_patch_scalar_type,
                  t_hermite, t_tri3, t_line3>());
        case DRT::Element::tri6:
          return Teuchos::rcp(
              new BeamToSolidSurfaceMeshtyingPairMortarFAD<line_to_surface_patch_scalar_type,
                  t_hermite, t_tri6, t_line3>());
        case DRT::Element::quad4:
          return Teuchos::rcp(
              new BeamToSolidSurfaceMeshtyingPairMortarFAD<line_to_surface_patch_scalar_type,
                  t_hermite, t_quad4, t_line3>());
        case DRT::Element::quad8:
          return Teuchos::rcp(
              new BeamToSolidSurfaceMeshtyingPairMortarFAD<line_to_surface_patch_scalar_type,
                  t_hermite, t_quad8, t_line3>());
        case DRT::Element::quad9:
          return Teuchos::rcp(
              new BeamToSolidSurfaceMeshtyingPairMortarFAD<line_to_surface_patch_scalar_type,
                  t_hermite, t_quad9, t_line3>());
        case DRT::Element::nurbs9:
          return Teuchos::rcp(new BeamToSolidSurfaceMeshtyingPairMortarFAD<
              line_to_surface_patch_nurbs_scalar_type<t_hermite, t_nurbs9>, t_hermite, t_nurbs9,
              t_line3>());
        default:
          dserror("Wrong element type for surface element.");
      }
      break;
    }
    case INPAR::BEAMTOSOLID::BeamToSolidMortarShapefunctions::line4:
    {
      switch (surface_shape)
      {
        case DRT::Element::tri3:
          return Teuchos::rcp(
              new BeamToSolidSurfaceMeshtyingPairMortarFAD<line_to_surface_patch_scalar_type,
                  t_hermite, t_tri3, t_line4>());
        case DRT::Element::tri6:
          return Teuchos::rcp(
              new BeamToSolidSurfaceMeshtyingPairMortarFAD<line_to_surface_patch_scalar_type,
                  t_hermite, t_tri6, t_line4>());
        case DRT::Element::quad4:
          return Teuchos::rcp(
              new BeamToSolidSurfaceMeshtyingPairMortarFAD<line_to_surface_patch_scalar_type,
                  t_hermite, t_quad4, t_line4>());
        case DRT::Element::quad8:
          return Teuchos::rcp(
              new BeamToSolidSurfaceMeshtyingPairMortarFAD<line_to_surface_patch_scalar_type,
                  t_hermite, t_quad8, t_line4>());
        case DRT::Element::quad9:
          return Teuchos::rcp(
              new BeamToSolidSurfaceMeshtyingPairMortarFAD<line_to_surface_patch_scalar_type,
                  t_hermite, t_quad9, t_line4>());
        case DRT::Element::nurbs9:
          return Teuchos::rcp(new BeamToSolidSurfaceMeshtyingPairMortarFAD<
              line_to_surface_patch_nurbs_scalar_type<t_hermite, t_nurbs9>, t_hermite, t_nurbs9,
              t_line4>());
        default:
          dserror("Wrong element type for surface element.");
      }
      break;
    }
    default:
      dserror("Wrong mortar shape function.");
  }

  return Teuchos::null;
}
