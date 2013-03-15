/*!----------------------------------------------------------------------
\file so_disp_input.cpp
\brief

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>

*----------------------------------------------------------------------*/


#include "so_disp.H"
#include "../drt_lib/drt_linedefinition.H"
#include "../drt_mat/so3_material.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool DRT::ELEMENTS::SoDisp::ReadElement(const std::string& eletype,
                                        const std::string& distype,
                                        DRT::INPUT::LineDefinition* linedef)
{
  // read number of material model
  int material = 0;
  linedef->ExtractInt("MAT",material);
  SetMaterial(material);

  DiscretizationType shape = StringToDistype(distype);

  switch (shape)
  {
  case hex8: case hex20: case hex27:
  {
    std::vector<int> ngp;
    linedef->ExtractIntVector("GP",ngp);

    switch (ngp[0])
    {
    case 1:
      gaussrule_ = DRT::UTILS::intrule_hex_1point;
      break;
    case 2:
      gaussrule_ = DRT::UTILS::intrule_hex_8point;
      break;
    case 3:
      gaussrule_ = DRT::UTILS::intrule_hex_27point;
      break;
    default:
      dserror("Reading of SOLID3 element failed: Gaussrule for hexaeder not supported!");
    }
    break;
  }
  case pyramid5:
  {
    int ngp;
    linedef->ExtractInt("GP_PYRAMID",ngp);
    switch (ngp)
    {
    case 1:
      gaussrule_ = DRT::UTILS::intrule_pyramid_1point;
      break;
    case 8:
      gaussrule_ = DRT::UTILS::intrule_pyramid_8point;
      break;
    default:
      dserror("Reading of SOLID3 element failed: Gaussrule for pyramid not supported!\n");
    }
    break;
  }
  case tet4: case tet10:
  {
    int ngp;
    linedef->ExtractInt("GP_TET",ngp);

    std::string buffer;
    linedef->ExtractString("GP_ALT",buffer);

    switch(ngp)
    {
    case 1:
      if (buffer=="standard")
        gaussrule_ = DRT::UTILS::intrule_tet_1point;
      else
        dserror("Reading of SOLID3 element failed: GP_ALT: gauss-radau not possible!\n");
      break;
    case 4:
      if (buffer=="standard")
        gaussrule_ = DRT::UTILS::intrule_tet_4point;
      else if (buffer=="gaussrad")
        gaussrule_ = DRT::UTILS::intrule_tet_4point_gauss_radau;
      else
        dserror("Reading of SOLID3 element failed: GP_ALT\n");
      break;
    case 10:
      if (buffer=="standard")
        gaussrule_ = DRT::UTILS::intrule_tet_5point;
      else
        dserror("Reading of SOLID3 element failed: GP_ALT: gauss-radau not possible!\n");
      break;
    default:
      dserror("Reading of SOLID3 element failed: Gaussrule for tetraeder not supported!\n");
    }
    break;
  } // end reading gaussian points for tetrahedral elements
  default:
    dserror("Reading of SOLID3 element failed: integration points\n");
  } // end switch distype

  // set up of materials with GP data (e.g., history variables)
  Teuchos::RCP<MAT::So3Material> so3mat = Teuchos::rcp_dynamic_cast<MAT::So3Material>(Material());
  so3mat->Setup(numgpt_disp_, linedef);

  std::string buffer;
  linedef->ExtractString("KINEM",buffer);

  // geometrically linear
  if      (buffer=="Geolin")    kintype_ = sodisp_geolin;
  // geometrically non-linear with Total Lagrangean approach
  else if (buffer=="Totlag")    kintype_ = sodisp_totlag;
  // geometrically non-linear with Updated Lagrangean approach
  else if (buffer=="Updlag")
  {
    kintype_ = sodisp_updlag;
    dserror("Updated Lagrange for SOLID3 is not implemented!");
  }
  else dserror("Reading of SOLID3 element failed");

  numnod_disp_ = NumNode();      // number of nodes
  numdof_disp_ = NumNode() * NODDOF_DISP;     // total dofs per element
  const DRT::UTILS::IntegrationPoints3D  intpoints(gaussrule_);
  numgpt_disp_ = intpoints.nquad;      // total gauss points per element

  return true;
}
