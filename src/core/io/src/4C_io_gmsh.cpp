/*----------------------------------------------------------------------*/
/*! \file

\brief simple element print library for Gmsh


\level 2

*/


#include "4C_io_gmsh.hpp"

#include "4C_fem_discretization.hpp"
#include "4C_fem_general_extract_values.hpp"
#include "4C_fem_geometry_position_array.hpp"
#include "4C_io_control.hpp"
#include "4C_linalg_serialdensevector.hpp"
#include "4C_linalg_utils_sparse_algebra_manipulation.hpp"  // Core::LinAlg::Export
#include "4C_rebalance_binning_based.hpp"

FOUR_C_NAMESPACE_OPEN


/*------------------------------------------------------------------------------------------------*
 | write scalar field to Gmsh postprocessing file                                     henke 12/09 |
 *------------------------------------------------------------------------------------------------*/
void Core::IO::Gmsh::ScalarFieldToGmsh(const Teuchos::RCP<Core::FE::Discretization> discret,
    const Teuchos::RCP<const Epetra_Vector> scalarfield_row, std::ostream& s)
{
  // tranform solution vector from dof_row_map to DofColMap
  const Teuchos::RCP<const Epetra_Vector> scalarfield =
      Core::Rebalance::GetColVersionOfRowVector(discret, scalarfield_row);

  // loop all row elements on this processor
  for (int iele = 0; iele < discret->num_my_row_elements(); ++iele)
  {
    const Core::Elements::Element* ele = discret->l_row_element(iele);
    const Core::FE::CellType distype = ele->shape();
    const int numnode = distypeToGmshNumNode(distype);

    Core::LinAlg::SerialDenseMatrix xyze(3, numnode);

    const Core::Nodes::Node* const* nodes = ele->nodes();
    for (int inode = 0; inode < numnode; ++inode)
    {
      xyze(0, inode) = nodes[inode]->x()[0];
      xyze(1, inode) = nodes[inode]->x()[1];
      xyze(2, inode) = nodes[inode]->x()[2];
    }

    s << "S";  // scalar field indicator
    s << distypeToGmshElementHeader(distype);

    // write node coordinates to Gmsh stream
    CoordinatesToStream(xyze, distype, s);

    std::vector<int> lm;
    std::vector<int> lmowner;
    std::vector<int> lmstride;
    ele->location_vector(*discret, lm, lmowner, lmstride);

    // extract local values from the global vector
    Core::LinAlg::SerialDenseVector myscalarfield(lm.size());
    Core::FE::ExtractMyValues(*scalarfield, myscalarfield, lm);

    // write scalar field to Gmsh stream
    ScalarFieldToStream(myscalarfield, distype, s);

    s << "\n";
  }
}


/*------------------------------------------------------------------------------------------------*
 | write scalar field to Gmsh postprocessing file                                     henke 12/09 |
 *------------------------------------------------------------------------------------------------*/
void Core::IO::Gmsh::ScalarFieldDofBasedToGmsh(const Teuchos::RCP<Core::FE::Discretization> discret,
    const Teuchos::RCP<const Epetra_Vector> scalarfield_row, const int nds, std::ostream& s)
{
  // tranform solution vector from dof_row_map to DofColMap
  const Teuchos::RCP<const Epetra_Vector> scalarfield =
      Core::Rebalance::GetColVersionOfRowVector(discret, scalarfield_row, nds);

  // loop all row elements on this processor
  for (int iele = 0; iele < discret->num_my_row_elements(); ++iele)
  {
    const Core::Elements::Element* ele = discret->l_row_element(iele);
    const Core::FE::CellType distype = ele->shape();
    const int numnode = distypeToGmshNumNode(distype);

    Core::LinAlg::SerialDenseMatrix xyze(3, numnode);

    const Core::Nodes::Node* const* nodes = ele->nodes();
    for (int inode = 0; inode < numnode; ++inode)
    {
      xyze(0, inode) = nodes[inode]->x()[0];
      xyze(1, inode) = nodes[inode]->x()[1];
      xyze(2, inode) = nodes[inode]->x()[2];
    }

    s << "S";  // scalar field indicator
    s << distypeToGmshElementHeader(distype);

    // write node coordinates to Gmsh stream
    CoordinatesToStream(xyze, distype, s);

    Core::Elements::Element::LocationArray la(discret->num_dof_sets());
    ele->location_vector(*discret, la, false);

    // extract local values from the global vector
    Core::LinAlg::SerialDenseVector myscalarfield(la[nds].lm_.size());
    Core::FE::ExtractMyValues(*scalarfield, myscalarfield, la[nds].lm_);

    //    // Extract velocity from local velnp_
    //    Core::LinAlg::SerialDenseMatrix myvectorfield(nsd,numnode);
    //    for (int inode = 0; inode < numnode; ++inode)
    //    {
    //      for (int idim = 0; idim < nsd; ++idim)
    //      {
    //        double val = extractmyvectorfield(idim + inode*nsd);
    //        myvectorfield(idim,inode)= val;
    //        if( displacenodes )
    //          xyze(idim,inode) += val;
    //      }
    //    }
    //
    //    std::vector<int> lm;
    //    std::vector<int> lmowner;
    //    std::vector<int> lmstride;
    //    ele->LocationVector(*discret, lm, lmowner, lmstride);
    //
    //    // extract local values from the global vector
    //    Core::LinAlg::SerialDenseVector myscalarfield(lm.size());
    //    Core::FE::ExtractMyValues(*scalarfield, myscalarfield, lm);

    // write scalar field to Gmsh stream
    ScalarFieldToStream(myscalarfield, distype, s);

    s << "\n";
  }
}

/*------------------------------------------------------------------------------------------------*
 | write scalar element-based field to Gmsh postprocessing file                  winklmaier 12/09 |
 *------------------------------------------------------------------------------------------------*/
void Core::IO::Gmsh::ScalarElementFieldToGmsh(const Teuchos::RCP<Core::FE::Discretization> discret,
    const Teuchos::RCP<const Epetra_Vector> scalarfield_ele_row, std::ostream& s)
{
  if (scalarfield_ele_row->Map().SameAs(*discret->element_row_map()) == false)
  {
    std::cout << scalarfield_ele_row->Map() << std::endl
              << *discret->element_row_map() << std::endl;
    FOUR_C_THROW("The written field should be based on the element row map");
  }

  // currently the element value is used for every node of the element so that
  // we result in the correct solution - gmsh writes nodal values for every
  // element, so that the same node can get different values for different
  // elements (similar to discontinuous methods)
  for (int iele = 0; iele < discret->num_my_row_elements(); ++iele)
  {
    const Core::Elements::Element* ele = discret->l_row_element(iele);
    const Core::FE::CellType distype = ele->shape();
    const int numnode = distypeToGmshNumNode(distype);

    Core::LinAlg::SerialDenseMatrix xyze(3, numnode);

    const Core::Nodes::Node* const* nodes = ele->nodes();
    for (int inode = 0; inode < numnode; ++inode)
    {
      xyze(0, inode) = nodes[inode]->x()[0];
      xyze(1, inode) = nodes[inode]->x()[1];
      xyze(2, inode) = nodes[inode]->x()[2];
    }

    s << "S";  // scalar field indicator
    s << distypeToGmshElementHeader(distype);

    // write node coordinates to Gmsh stream
    CoordinatesToStream(xyze, distype, s);

    double eleval = (*scalarfield_ele_row)[iele];

    // constant value for all nodes
    Core::LinAlg::SerialDenseVector myscalarfield(ele->num_node());
    for (int inode = 0; inode < numnode; ++inode) myscalarfield(inode) = eleval;

    // write scalar field to Gmsh stream
    ScalarFieldToStream(myscalarfield, distype, s);

    s << "\n";
  }
}

/*------------------------------------------------------------------------------------------------*
 | write dof-based vector field to Gmsh postprocessing file                           henke 12/09 |
 *------------------------------------------------------------------------------------------------*/
void Core::IO::Gmsh::VectorFieldDofBasedToGmsh(const Teuchos::RCP<Core::FE::Discretization> discret,
    const Teuchos::RCP<const Epetra_Vector> vectorfield_row, std::ostream& s, const int nds,
    bool displacenodes)
{
  // tranform solution vector from dof_row_map to DofColMap
  const Teuchos::RCP<const Epetra_Vector> vectorfield =
      Core::Rebalance::GetColVersionOfRowVector(discret, vectorfield_row, nds);

  // loop all row elements on this processor
  for (int iele = 0; iele < discret->num_my_row_elements(); ++iele)
  {
    const Core::Elements::Element* ele = discret->l_row_element(iele);
    const Core::FE::CellType distype = ele->shape();
    const int numnode = distypeToGmshNumNode(distype);
    const int nsd = Core::FE::getDimension(distype);

    Core::LinAlg::SerialDenseMatrix xyze(nsd, numnode);

    const Core::Nodes::Node* const* nodes = ele->nodes();
    for (int inode = 0; inode < numnode; ++inode)
    {
      for (int idim = 0; idim < nsd; ++idim) xyze(idim, inode) = nodes[inode]->x()[idim];
    }

    Core::Elements::Element::LocationArray la(discret->num_dof_sets());
    ele->location_vector(*discret, la, false);

    // extract local values from the global vector
    Core::LinAlg::SerialDenseVector extractmyvectorfield(la[nds].lm_.size());
    Core::FE::ExtractMyValues(*vectorfield, extractmyvectorfield, la[nds].lm_);

    // Extract velocity from local velnp_
    Core::LinAlg::SerialDenseMatrix myvectorfield(nsd, numnode);
    for (int inode = 0; inode < numnode; ++inode)
    {
      for (int idim = 0; idim < nsd; ++idim)
      {
        double val = extractmyvectorfield(idim + inode * nsd);
        myvectorfield(idim, inode) = val;
        if (displacenodes) xyze(idim, inode) += val;
      }
    }

    s << "V";  // vector field indicator
    s << distypeToGmshElementHeader(distype);

    // write node coordinates to Gmsh stream
    CoordinatesToStream(xyze, distype, s);

    // write vector field to Gmsh stream
    // remark: only the first 3 components are written to the Gmsh postprocessing file
    //         -> e.g. pressure in velnp_ fluid state vector is not written ('myvectorfield':
    //         velx,vely,velz,pressure)
    VectorFieldToStream(myvectorfield, distype, s);

    s << "\n";
  }
}

/*------------------------------------------------------------------------------------------------*
 | write multivector field to Gmsh postprocessing file                               winter 04/17 |
 *------------------------------------------------------------------------------------------------*/
void Core::IO::Gmsh::VectorFieldMultiVectorDofBasedToGmsh(
    const Teuchos::RCP<const Core::FE::Discretization> discret,
    const Teuchos::RCP<const Epetra_MultiVector> vectorfield_row, std::ostream& s, const int nds)
{
  // TODO: Remove dependence on size of Epetra_Multivector!!!

  // tranform solution vector from dof_row_map to DofColMap
  const Teuchos::RCP<Epetra_MultiVector> vectorfield =
      Teuchos::rcp(new Epetra_MultiVector(*discret->dof_col_map(nds), 3, true));
  Core::LinAlg::Export(*vectorfield_row, *vectorfield);

  // loop all row elements on this processor
  for (int iele = 0; iele < discret->num_my_row_elements(); ++iele)
  {
    const Core::Elements::Element* ele = discret->l_row_element(iele);
    const Core::FE::CellType distype = ele->shape();
    const int numnode = distypeToGmshNumNode(distype);
    const int nsd = Core::FE::getDimension(distype);

    Core::LinAlg::SerialDenseMatrix xyze(nsd, numnode);

    const Core::Nodes::Node* const* nodes = ele->nodes();
    for (int inode = 0; inode < numnode; ++inode)
    {
      for (int idim = 0; idim < nsd; ++idim)
      {
        xyze(idim, inode) = nodes[inode]->x()[idim];
      }
    }

    s << "V";  // vector field indicator
    s << distypeToGmshElementHeader(distype);

    // write node coordinates to Gmsh stream
    CoordinatesToStream(xyze, distype, s);

    Core::Elements::Element::LocationArray la(discret->num_dof_sets());
    ele->location_vector(*discret, la, false);

    // extract local values from the global vector
    // Core::LinAlg::SerialDenseVector extractmyvectorfield(la[nds].lm_.size());
    Core::LinAlg::SerialDenseMatrix myvectorfield(nsd, numnode);
    std::vector<double> local_vector;
    Core::FE::ExtractMyValues(*vectorfield, local_vector, la[nds].lm_);

    for (int inode = 0; inode < numnode; ++inode)
    {
      for (int idim = 0; idim < nsd; ++idim)
      {
        myvectorfield(idim, inode) = local_vector[idim + nsd * inode];
      }
    }

    // write vector field to Gmsh stream
    VectorFieldToStream(myvectorfield, distype, s);

    s << "\n";
  }
}

/*------------------------------------------------------------------------------------------------*
 | write dof-based vector field to Gmsh postprocessing file at current position      schott 12/09 |
 *------------------------------------------------------------------------------------------------*/
void Core::IO::Gmsh::SurfaceVectorFieldDofBasedToGmsh(
    const Teuchos::RCP<Core::FE::Discretization> discret,
    const Teuchos::RCP<const Epetra_Vector> vectorfield_row,
    std::map<int, Core::LinAlg::Matrix<3, 1>>& currpos, std::ostream& s, const int nsd,
    const int numdofpernode)
{
  // tranform solution vector from dof_row_map to DofColMap
  const Teuchos::RCP<const Epetra_Vector> vectorfield =
      Core::Rebalance::GetColVersionOfRowVector(discret, vectorfield_row);

  // loop all row elements on this processor
  for (int iele = 0; iele < discret->num_my_row_elements(); ++iele)
  {
    const Core::Elements::Element* ele = discret->l_row_element(iele);
    const Core::FE::CellType distype = ele->shape();
    const int numnode = distypeToGmshNumNode(distype);

    Core::LinAlg::SerialDenseMatrix xyze(nsd, numnode);

    const Core::Nodes::Node* const* nodes = ele->nodes();
    for (int inode = 0; inode < numnode; ++inode)
    {
      int nid = nodes[inode]->id();
      for (int idim = 0; idim < nsd; ++idim)
        xyze(idim, inode) = ((currpos.find(nid))->second)(idim, 0);
    }

    s << "V";  // vector field indicator
    s << distypeToGmshElementHeader(distype);

    // write node coordinates to Gmsh stream
    CoordinatesToStream(xyze, distype, s);

    std::vector<int> lm;
    std::vector<int> lmowner;
    std::vector<int> lmstride;
    ele->location_vector(*discret, lm, lmowner, lmstride);

    // extract local values from the global vector
    Core::LinAlg::SerialDenseVector extractmyvectorfield(lm.size());
    Core::FE::ExtractMyValues(*vectorfield, extractmyvectorfield, lm);

    // Extract velocity from local velnp_
    Core::LinAlg::SerialDenseMatrix myvectorfield(nsd, numnode);
    for (int inode = 0; inode < numnode; ++inode)
      for (int idim = 0; idim < nsd; ++idim)
        myvectorfield(idim, inode) = extractmyvectorfield(idim + inode * numdofpernode);

    // write vector field to Gmsh stream
    // remark: only the first 3 components are written to the Gmsh postprocessing file
    //         -> e.g. pressure in velnp_ fluid state vector is not written ('myvectorfield':
    //         velx,vely,velz,pressure)
    VectorFieldToStream(myvectorfield, distype, s);

    s << "\n";
  }
}

/*------------------------------------------------------------------------------------------------*
 | write scalar / vector from a dof-based vector field (e.g. velocity)                            |
 | to Gmsh postprocessing file                                                         ehrl 05/11 |
 *------------------------------------------------------------------------------------------------*/
void Core::IO::Gmsh::VelocityPressureFieldDofBasedToGmsh(
    const Teuchos::RCP<Core::FE::Discretization> discret,
    const Teuchos::RCP<const Epetra_Vector> vectorfield_row, const std::string field,
    std::ostream& s, const int nds)
{
  // tranform solution vector from dof_row_map to DofColMap
  const Teuchos::RCP<const Epetra_Vector> vectorfield =
      Core::Rebalance::GetColVersionOfRowVector(discret, vectorfield_row, nds);

  // loop all row elements on this processor
  for (int iele = 0; iele < discret->num_my_row_elements(); ++iele)
  {
    const Core::Elements::Element* ele = discret->l_row_element(iele);
    const Core::FE::CellType distype = ele->shape();
    const int numnode = distypeToGmshNumNode(distype);
    const int nsd = Core::FE::getDimension(distype);

    Core::LinAlg::SerialDenseMatrix xyze(nsd, numnode);

    const Core::Nodes::Node* const* nodes = ele->nodes();
    for (int inode = 0; inode < numnode; ++inode)
    {
      for (int idim = 0; idim < nsd; ++idim) xyze(idim, inode) = nodes[inode]->x()[idim];
    }

    Core::Elements::Element::LocationArray la(discret->num_dof_sets());
    ele->location_vector(*discret, la, false);

    // extract local values from the global vector
    Core::LinAlg::SerialDenseVector extractmyvectorfield(la[nds].lm_.size());
    Core::FE::ExtractMyValues(*vectorfield, extractmyvectorfield, la[nds].lm_);

    if (field == "pressure")
    {
      // Extract scalar from local velnp_
      Core::LinAlg::SerialDenseVector myscalarfield(numnode);
      for (int inode = 0; inode < numnode; ++inode)
        myscalarfield(inode) = extractmyvectorfield(nsd + inode * (nsd + 1));

      // replace function: cellWithScalarFieldToStream(distype, myscalarfield, xyze, s);
      {
        s << "S";  // scalar field indicator
        s << distypeToGmshElementHeader(distype);
        if (nsd == 3)
          CoordinatesToStream(xyze, distype, s);
        else if (nsd == 2)
          CoordinatesToStream2D(xyze, distype, s);
        else
          FOUR_C_THROW("only two- and three-dimensional domains are supported");

        ScalarFieldToStream(myscalarfield, distype, s);
        s << "\n";
      }
    }
    else if (field == "velocity")
    {
      // Extract velocity from local velnp_
      Core::LinAlg::SerialDenseMatrix myvectorfield(nsd, numnode);
      for (int inode = 0; inode < numnode; ++inode)
        for (int idim = 0; idim < nsd; ++idim)
          myvectorfield(idim, inode) = extractmyvectorfield(idim + inode * (nsd + 1));

      // replace function : cellWithVectorFieldToStream(distype, myvectorfield, xyze, s);
      {
        s << "V";  // scalar field indicator
        s << distypeToGmshElementHeader(distype);
        if (nsd == 3)
        {
          CoordinatesToStream(xyze, distype, s);
          VectorFieldToStream(myvectorfield, distype, s);
        }
        else if (nsd == 2)
        {
          CoordinatesToStream2D(xyze, distype, s);
          VectorFieldToStream2D(myvectorfield, distype, s);
        }
        else
          FOUR_C_THROW("");
        s << "\n";
      }
    }
    else
      FOUR_C_THROW("The choosen field does not exist (wrong writting, ...)");

    s << "\n";
  }
}

/*------------------------------------------------------------------------------------------------*
 | write node-based vector field to Gmsh postprocessing file                          henke 12/09 |
 *------------------------------------------------------------------------------------------------*/
void Core::IO::Gmsh::VectorFieldNodeBasedToGmsh(
    const Teuchos::RCP<const Core::FE::Discretization> discret,
    const Teuchos::RCP<const Epetra_MultiVector> vectorfield_row, std::ostream& s)
{
  // tranform solution vector from NodeRowMap to NodeColMap
  // remark: Core::Rebalance::GetColVersionOfRowVector() does only work for Epetra_Vectors
  // on dof_row_map
  const Teuchos::RCP<Epetra_MultiVector> vectorfield =
      Teuchos::rcp(new Epetra_MultiVector(*discret->node_col_map(), 3, true));
  Core::LinAlg::Export(*vectorfield_row, *vectorfield);

  // loop all row elements on this processor
  for (int iele = 0; iele < discret->num_my_row_elements(); ++iele)
  {
    const Core::Elements::Element* ele = discret->l_row_element(iele);
    const Core::FE::CellType distype = ele->shape();
    const int numnode = ele->num_node();
    const int nsd = Core::FE::getDimension(distype);

    Core::LinAlg::SerialDenseMatrix xyze(nsd, numnode);

    const Core::Nodes::Node* const* nodes = ele->nodes();
    for (int inode = 0; inode < numnode; ++inode)
    {
      for (int idim = 0; idim < nsd; ++idim) xyze(idim, inode) = nodes[inode]->x()[idim];
    }

    s << "V";  // vector field indicator
    s << distypeToGmshElementHeader(distype);

    // write node coordinates to Gmsh stream
    CoordinatesToStream(xyze, distype, s);

    // extract local values from the global vector
    Core::LinAlg::SerialDenseMatrix myvectorfield(nsd, numnode);
    Core::FE::ExtractMyNodeBasedValues(ele, myvectorfield, vectorfield, nsd);

    // write vector field to Gmsh stream
    VectorFieldToStream(myvectorfield, distype, s);

    s << "\n";
  }
}


/*------------------------------------------------------------------------------------------------*
 | write node-based scalar field to Gmsh postprocessing file                      rasthofer 09/13 |
 *------------------------------------------------------------------------------------------------*/
void Core::IO::Gmsh::ScalarFieldNodeBasedToGmsh(
    const Teuchos::RCP<const Core::FE::Discretization> discret,
    const Teuchos::RCP<const Epetra_Vector> scalarfield_row, std::ostream& s)
{
  // tranform solution vector from NodeRowMap to NodeColMap
  // remark: Core::Rebalance::GetColVersionOfRowVector() does only work for Epetra_Vectors
  // on dof_row_map
  //         something similar is done in COMBUST::FlameFront::ProcessFlameFront, although not for
  //         Epetra_MultiVectors
  const Teuchos::RCP<Epetra_Vector> scalarfield =
      Teuchos::rcp(new Epetra_Vector(*discret->node_col_map(), true));
  Core::LinAlg::Export(*scalarfield_row, *scalarfield);

  // loop all row elements on this processor
  for (int iele = 0; iele < discret->num_my_row_elements(); ++iele)
  {
    const Core::Elements::Element* ele = discret->l_row_element(iele);
    const Core::FE::CellType distype = ele->shape();
    const int numnode = ele->num_node();

    Core::LinAlg::SerialDenseMatrix xyze(3, numnode);

    const Core::Nodes::Node* const* nodes = ele->nodes();
    for (int inode = 0; inode < numnode; ++inode)
    {
      xyze(0, inode) = nodes[inode]->x()[0];
      xyze(1, inode) = nodes[inode]->x()[1];
      xyze(2, inode) = nodes[inode]->x()[2];
    }

    s << "S";  // vector field indicator
    s << distypeToGmshElementHeader(distype);

    // write node coordinates to Gmsh stream
    CoordinatesToStream(xyze, distype, s);

    // extract local values from the global vector
    Core::LinAlg::SerialDenseVector myscalarfield(numnode);
    Core::FE::ExtractMyNodeBasedValues(ele, myscalarfield, scalarfield, 1);

    // write vector field to Gmsh stream
    ScalarFieldToStream(myscalarfield, distype, s);

    s << "\n";
  }
}

void Core::IO::Gmsh::ScalarToStream(
    const Core::LinAlg::Matrix<3, 1>& pointXYZ,  ///< coordinates of point
    const double scalarvalue,                    ///< scalar value at this point
    std::ostream& s                              ///< stream
)
{
  s.setf(std::ios::scientific, std::ios::floatfield);
  s.precision(12);

  s << "SP(";  // scalar field indicator
  s << pointXYZ(0) << ",";
  s << pointXYZ(1) << ",";
  s << pointXYZ(2) << ")";

  s << "{" << scalarvalue << "};";
  s << "\n";
}

void Core::IO::Gmsh::VectorToStream(
    const Core::LinAlg::Matrix<3, 1>& pointXYZ,     ///< coordinates of point
    const Core::LinAlg::Matrix<3, 1>& vectorvalue,  ///< vector at this point
    std::ostream& s                                 ///< stream
)
{
  s.setf(std::ios::scientific, std::ios::floatfield);
  s.precision(12);

  s << "VP(";  // vector field indicator
  s << pointXYZ(0) << ",";
  s << pointXYZ(1) << ",";
  s << pointXYZ(2) << ")";

  s << "{";
  s << vectorvalue(0);
  s << "," << vectorvalue(1);
  s << "," << vectorvalue(2);
  s << "};";  //<< endl;
  s << "\n";
}

void Core::IO::Gmsh::elementAtInitialPositionToStream(
    const double scalar, const Core::Elements::Element* ele, std::ostream& s)
{
  const Core::Nodes::Node* const* nodes = ele->nodes();

  const Core::FE::CellType distype = ele->shape();
  const int numnode = distypeToGmshNumNode(distype);

  s.setf(std::ios::scientific, std::ios::floatfield);
  s.precision(12);

  s << "S" << distypeToGmshElementHeader(distype) << "(";
  for (int i = 0; i < numnode; ++i)
  {
    const Core::Nodes::Node* node = nodes[i];
    const auto& x = node->x();
    s << x[0] << ",";
    s << x[1] << ",";
    s << x[2];
    if (i < numnode - 1)
    {
      s << ",";
    }
  }
  s << ")";
  // values
  ScalarToStream(scalar, distype, s);
  s << "\n";
}


std::string Core::IO::Gmsh::elementAtInitialPositionToString(
    const double scalar, const Core::Elements::Element* ele)
{
  std::ostringstream s;
  elementAtInitialPositionToStream(scalar, ele, s);
  return s.str();
}


void Core::IO::Gmsh::elementAtCurrentPositionToStream(const double scalar,
    const Core::Elements::Element* ele,
    const std::map<int, Core::LinAlg::Matrix<3, 1>>& currentelepositions, std::ostream& s)
{
  Core::IO::Gmsh::cellWithScalarToStream(
      ele->shape(), scalar, Core::Geo::getCurrentNodalPositions(ele, currentelepositions), s);
}


std::string Core::IO::Gmsh::elementAtCurrentPositionToString(const double scalar,
    const Core::Elements::Element* ele,
    const std::map<int, Core::LinAlg::Matrix<3, 1>>& currentelepositions)
{
  std::ostringstream s;
  Core::IO::Gmsh::elementAtCurrentPositionToStream(scalar, ele, currentelepositions, s);
  return s.str();
}


std::string Core::IO::Gmsh::text3dToString(
    const Core::LinAlg::Matrix<3, 1>& xyz,  ///< 3d Position of text
    const std::string& text,                ///< text to be printed
    const int fontsize                      ///< font size
)
{
  std::ostringstream s;

  s << "T3";
  // coordinates
  s << "(";
  s << std::scientific << xyz(0) << ",";
  s << std::scientific << xyz(1) << ",";
  s << std::scientific << xyz(2) << ",";
  s << fontsize << ")";
  s << "{\"" << text << "\"};";
  s << "\n";
  return s.str();
}

void Core::IO::Gmsh::disToStream(const std::string& text, const double scalar,
    const Teuchos::RCP<Core::FE::Discretization> dis, std::ostream& s)
{
  s << "View \" " << text << " Elements \" {\n";
  for (int i = 0; i < dis->num_my_row_elements(); ++i)
  {
    const Core::Elements::Element* actele = dis->l_row_element(i);
    Core::IO::Gmsh::elementAtInitialPositionToStream(scalar, actele, s);
  }
  s << "};\n";
}

std::string Core::IO::Gmsh::disToString(
    const std::string& text, const double scalar, const Teuchos::RCP<Core::FE::Discretization> dis)
{
  std::ostringstream s;
  disToStream(text, scalar, dis, s);
  return s.str();
}

void Core::IO::Gmsh::disToStream(const std::string& text, const double scalar,
    const Teuchos::RCP<Core::FE::Discretization> dis,
    const std::map<int, Core::LinAlg::Matrix<3, 1>>& currentpositions, std::ostream& s)
{
  s << "View \" " << text << " Elements \" {\n";

  for (int i = 0; i < dis->num_my_col_elements(); ++i)
  {
    const Core::Elements::Element* actele = dis->l_col_element(i);
    Core::IO::Gmsh::cellWithScalarToStream(
        actele->shape(), scalar, Core::Geo::getCurrentNodalPositions(actele, currentpositions), s);
  }
  s << "};\n";
}

std::string Core::IO::Gmsh::disToString(const std::string& text, const double scalar,
    const Teuchos::RCP<Core::FE::Discretization> dis,
    const std::map<int, Core::LinAlg::Matrix<3, 1>>& currentpositions)
{
  std::ostringstream s;
  disToStream(text, scalar, dis, currentpositions, s);
  return s.str();
}

std::string Core::IO::Gmsh::GetNewFileNameAndDeleteOldFiles(const std::string& filename_base,
    const std::string& file_name_prefix, const int& actstep, const int& step_diff,
    const bool screen_out, const int pid)
{
  std::ostringstream filename;
  std::ostringstream filenamedel;

  std::ostringstream pid_stream;
  pid_stream << ".p" << std::setw(2) << std::setfill('0') << pid;

  filename << file_name_prefix << "." << filename_base << "_" << std::setw(5) << std::setfill('0')
           << actstep << pid_stream.str() << ".pos";
  filenamedel << file_name_prefix << "." << filename_base << "_" << std::setw(5)
              << std::setfill('0') << actstep - step_diff << pid_stream.str() << ".pos";
  std::remove(filenamedel.str().c_str());
  if (screen_out)
    std::cout << "writing " << std::left << std::setw(60) << filename.str() << "..." << '\n';
  return filename.str();
}

std::string Core::IO::Gmsh::GetFileName(const std::string& filename_base,
    const std::string& file_name_prefix, const int& actstep, const bool screen_out, const int pid)
{
  std::ostringstream filename;
  std::ostringstream pid_stream;
  pid_stream << ".p" << std::setw(2) << std::setfill('0') << pid;

  filename << file_name_prefix << "." << filename_base << "_" << std::setw(5) << std::setfill('0')
           << actstep << pid_stream.str() << ".pos";

  if (screen_out)
    std::cout << "writing " << std::left << std::setw(60) << filename.str() << "..." << '\n';

  return filename.str();
}

FOUR_C_NAMESPACE_CLOSE
