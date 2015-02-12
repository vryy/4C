/*!----------------------------------------------------------------------
\file xfem_utils.cpp
\brief Basic tools used in XFEM routines

<pre>
Maintainer:  Magnus Winter / Raffaela Kruse
             winter@lnm.mw.tum.de
             kruse@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15236
</pre>
*----------------------------------------------------------------------*/

#include "xfem_utils.H"

#include "../drt_lib/drt_discret_faces.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_utils_parallel.H"
#include "../linalg/linalg_utils.H"

#include "../drt_io/io_gmsh.H"

#include "../drt_lib/drt_element.H"

//Materials supported in XFEM currently
#include "../drt_mat/material.H"
#include "../drt_mat/matlist.H"
#include "../drt_mat/newtonianfluid.H"

void XFEM::UTILS::PrintDiscretizationToStream(
  Teuchos::RCP<DRT::Discretization> dis,
  const std::string& disname,
  bool elements,
  bool elecol,
  bool nodes,
  bool nodecol,
  bool faces,
  bool facecol,
  std::ostream& s,
  std::map<int, LINALG::Matrix<3,1> >* curr_pos
)
{
  if(elements)
  {
    // draw bg elements with associated gid
    s << "View \" " << disname;
    if(elecol)
    {
      s << " col e->Id() \" {\n";
      for (int i=0; i<dis->NumMyColElements(); ++i)
      {
        const DRT::Element* actele = dis->lColElement(i);
        if(curr_pos == NULL)
          IO::GMSH::elementAtInitialPositionToStream(double(actele->Id()), actele, s);
        else
          IO::GMSH::elementAtCurrentPositionToStream(double(actele->Id()), actele, *curr_pos, s);
      };
    }
    else
    {
      s << " row e->Id() \" {\n";
      for (int i=0; i<dis->NumMyRowElements(); ++i)
      {
        const DRT::Element* actele = dis->lRowElement(i);
        if(curr_pos == NULL)
          IO::GMSH::elementAtInitialPositionToStream(double(actele->Id()), actele, s);
        else
          IO::GMSH::elementAtCurrentPositionToStream(double(actele->Id()), actele, *curr_pos, s);
      };
    }
    s << "};\n";
  }

  if(nodes)
  {
    s << "View \" " << disname;
    if(nodecol)
    {
      s << " col n->Id() \" {\n";
      for (int i=0; i<dis->NumMyColNodes(); ++i)
      {
        const DRT::Node* actnode = dis->lColNode(i);
        LINALG::Matrix<3,1> pos(true);

        if(curr_pos != NULL)
        {
          const LINALG::Matrix<3,1>& curr_x = curr_pos->find(actnode->Id())->second;
          pos(0) = curr_x(0);
          pos(1) = curr_x(1);
          pos(2) = curr_x(2);
        }
        else
        {
          const LINALG::Matrix<3,1> x(actnode->X());
          pos(0) = x(0);
          pos(1) = x(1);
          pos(2) = x(2);
        }
        IO::GMSH::cellWithScalarToStream(DRT::Element::point1, actnode->Id(), pos, s);
      }
    }
    else
    {
      s << " row n->Id() \" {\n";
      for (int i=0; i<dis->NumMyRowNodes(); ++i)
      {
        const DRT::Node* actnode = dis->lRowNode(i);
        LINALG::Matrix<3,1> pos(true);

        if(curr_pos != NULL)
        {
          const LINALG::Matrix<3,1>& curr_x = curr_pos->find(actnode->Id())->second;
          pos(0) = curr_x(0);
          pos(1) = curr_x(1);
          pos(2) = curr_x(2);
        }
        else
        {
          const LINALG::Matrix<3,1> x(actnode->X());
          pos(0) = x(0);
          pos(1) = x(1);
          pos(2) = x(2);
        }
        IO::GMSH::cellWithScalarToStream(DRT::Element::point1, actnode->Id(), pos, s);
      }
    }
    s << "};\n";
  }

  if(faces)
  {
    // cast to DiscretizationXFEM
    Teuchos::RCP<DRT::DiscretizationFaces> xdis = Teuchos::rcp_dynamic_cast<DRT::DiscretizationFaces>(dis, true);
    if (xdis == Teuchos::null)
      dserror("Failed to cast DRT::Discretization to DRT::DiscretizationFaces.");

    s << "View \" " << disname;

    if( xdis->FilledExtension() == true )     // faces output
    {

      if(facecol)
      {
        s << " col f->Id() \" {\n";
        for (int i=0; i<xdis->NumMyColFaces(); ++i)
        {
          const DRT::Element* actele = xdis->lColFace(i);
          if(curr_pos == NULL)
            IO::GMSH::elementAtInitialPositionToStream(double(actele->Id()), actele, s);
          else
            IO::GMSH::elementAtCurrentPositionToStream(double(actele->Id()), actele, *curr_pos, s);
        };
      }
      else
      {
        s << " row f->Id() \" {\n";
        for (int i=0; i<xdis->NumMyRowFaces(); ++i)
        {
          const DRT::Element* actele = xdis->lRowFace(i);
          if(curr_pos == NULL)
            IO::GMSH::elementAtInitialPositionToStream(double(actele->Id()), actele, s);
          else
            IO::GMSH::elementAtCurrentPositionToStream(double(actele->Id()), actele, *curr_pos, s);
        };
      }
      s << "};\n";
    }
  }
}

void XFEM::UTILS::ExtractNodeVectors(
  Teuchos::RCP<DRT::Discretization>    dis,
  std::map<int, LINALG::Matrix<3,1> >& nodevecmap)
{
  Teuchos::RCP<Epetra_Vector> dofrowvec     = LINALG::CreateVector(*dis->DofRowMap(),true);
  Teuchos::RCP<const Epetra_Vector> dispcol = DRT::UTILS::GetColVersionOfRowVector(dis,dofrowvec);
  nodevecmap.clear();

  for (int lid = 0; lid < dis->NumMyColNodes(); ++lid)
  {
    const DRT::Node* node = dis->lColNode(lid);
    std::vector<int> lm;
    dis->Dof(node, lm);
    std::vector<double> mydisp;
    DRT::UTILS::ExtractMyValues(*dispcol,mydisp,lm);
    if (mydisp.size() < 3)
      dserror("we need at least 3 dofs here");

    LINALG::Matrix<3,1> currpos;
    currpos(0) = node->X()[0] + mydisp[0];
    currpos(1) = node->X()[1] + mydisp[1];
    currpos(2) = node->X()[2] + mydisp[2];
    nodevecmap.insert(std::make_pair(node->Id(),currpos));
  }
}

// -------------------------------------------------------------------
// set master and slave parameters (winter 01/2015)
// -------------------------------------------------------------------
void XFEM::UTILS::GetVolumeCellMaterial(
  DRT::Element* actele,
  Teuchos::RCP<MAT::Material> & mat,
  GEO::CUT::Point::PointPosition position
)
{
  int position_id = 0;
  if (position == GEO::CUT::Point::inside)
    position_id = 1;
  else if (position != GEO::CUT::Point::outside)
    dserror("Volume cell is either undecided or on surface. That can't be good....");

  Teuchos::RCP<MAT::Material> material = actele->Material();

  if (material->MaterialType() == INPAR::MAT::m_matlist)
  {
    // get material list for this element
    const MAT::MatList* matlist = static_cast<const MAT::MatList*>(material.get());
    int numofmaterials = matlist->NumMat();

    //Error messages
    if(numofmaterials>2)
    {
      dserror("More than two materials is currently not supported.");
    }

    // set default id in list of materials
    int matid = -1;
    matid     = matlist->MatID(position_id);
    mat       = matlist->MaterialById(matid);

  }
  else
  {
    mat = material;
  }

  //std::cout << "mat->MaterialType(): " << mat->MaterialType() << std::endl;
  return;
}

/*----------------------------------------------------------------------*
 | Checks if Materials in parent and neighbor element are identical     |
 |                                                         winter 01/15 |
 *----------------------------------------------------------------------*/
void XFEM::UTILS::SafetyCheckMaterials(
  Teuchos::RCP<MAT::Material> &          pmat,
  Teuchos::RCP<MAT::Material> &          nmat
)
{

  //------------------------------ see whether materials in patch are equal

  if(pmat->MaterialType() != nmat->MaterialType())
    dserror(" not the same material for master and slave parent element");

  if(pmat->MaterialType() == INPAR::MAT::m_matlist)
    dserror("A matlist has been found in edge based stabilization! If you are running XTPF, check calls as this should NOT happen!!!");

  if( pmat->MaterialType() != INPAR::MAT::m_carreauyasuda
   && pmat->MaterialType() != INPAR::MAT::m_modpowerlaw
   && pmat->MaterialType() != INPAR::MAT::m_herschelbulkley
   && pmat->MaterialType() != INPAR::MAT::m_fluid)
    dserror("Material law for parent element is not a fluid");

  if(pmat->MaterialType() == INPAR::MAT::m_fluid)
  {
    {

      const MAT::NewtonianFluid* actmat_p = static_cast<const MAT::NewtonianFluid*>(pmat.get());
      const double pvisc=actmat_p->Viscosity();
      const double pdens = actmat_p->Density();

      const MAT::NewtonianFluid* actmat_m = static_cast<const MAT::NewtonianFluid*>(nmat.get());
      const double nvisc = actmat_m->Viscosity();
      const double ndens = actmat_m->Density();

      if(std::abs(nvisc - pvisc) > 1e-14)
      {
        std::cout << "Parent element viscosity: " << pvisc << " ,neighbor element viscosity: " << nvisc << std::endl;
        dserror("parent and neighbor element do not have the same viscosity!");
      }
      if(std::abs(ndens - pdens) > 1e-14)
      {
        std::cout << "Parent element density: " << pdens << " ,neighbor element density: " << ndens << std::endl;
        dserror("parent and neighbor element do not have the same density!");
      }
    }
  }
  else
  {
    dserror("up to now I expect a FLUID (m_fluid) material for edge stabilization\n");
  }

  return;
}
