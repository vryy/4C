/*----------------------------------------------------------------------*/
/*! \file
\brief Basic tools used in XFEM routines

\level 3


\warning this file should be cleaned up
*/
/*----------------------------------------------------------------------*/

#include "4C_xfem_utils.hpp"

#include "4C_fem_discretization_faces.hpp"
#include "4C_fem_general_element.hpp"
#include "4C_fem_general_extract_values.hpp"
#include "4C_linalg_utils_sparse_algebra_math.hpp"
#include "4C_mat_list.hpp"
#include "4C_mat_material_factory.hpp"
#include "4C_mat_newtonianfluid.hpp"
#include "4C_material_base.hpp"
#include "4C_rebalance_binning_based.hpp"

FOUR_C_NAMESPACE_OPEN

void XFEM::UTILS::extract_node_vectors(Teuchos::RCP<Core::FE::Discretization> dis,
    std::map<int, Core::LinAlg::Matrix<3, 1>>& nodevecmap, Teuchos::RCP<Epetra_Vector> idispnp)
{
  Teuchos::RCP<const Epetra_Vector> dispcol =
      Core::Rebalance::GetColVersionOfRowVector(dis, idispnp);
  nodevecmap.clear();

  for (int lid = 0; lid < dis->NumMyColNodes(); ++lid)
  {
    const Core::Nodes::Node* node = dis->lColNode(lid);
    std::vector<int> lm;
    dis->Dof(node, lm);
    std::vector<double> mydisp;
    Core::FE::ExtractMyValues(*dispcol, mydisp, lm);
    if (mydisp.size() < 3) FOUR_C_THROW("we need at least 3 dofs here");

    Core::LinAlg::Matrix<3, 1> currpos;
    currpos(0) = node->X()[0] + mydisp[0];
    currpos(1) = node->X()[1] + mydisp[1];
    currpos(2) = node->X()[2] + mydisp[2];
    nodevecmap.insert(std::make_pair(node->Id(), currpos));
  }
}

// -------------------------------------------------------------------
// set master and slave parameters (winter 01/2015)
// -------------------------------------------------------------------
void XFEM::UTILS::get_volume_cell_material(Core::Elements::Element* actele,
    Teuchos::RCP<Core::Mat::Material>& mat, Core::Geo::Cut::Point::PointPosition position)
{
  int position_id = 0;
  if (position == Core::Geo::Cut::Point::inside)  // minus domain, Omega^i with i<j
    position_id = 1;
  else if (position != Core::Geo::Cut::Point::outside)  // plus domain, \Omega^j with j>i
    FOUR_C_THROW("Volume cell is either undecided or on surface. That can't be good....");

  Teuchos::RCP<Core::Mat::Material> material = actele->Material();

  if (material->MaterialType() == Core::Materials::m_matlist)
  {
    // get material list for this element
    const Mat::MatList* matlist = static_cast<const Mat::MatList*>(material.get());
    int numofmaterials = matlist->NumMat();

    // Error messages
    if (numofmaterials > 2)
    {
      FOUR_C_THROW("More than two materials is currently not supported.");
    }

    // set default id in list of materials
    int matid = -1;
    matid = matlist->MatID(position_id);
    mat = matlist->MaterialById(matid);
  }
  else
  {
    mat = material;
  }

  return;
}

/*----------------------------------------------------------------------*
 | Checks if Materials in parent and neighbor element are identical     |
 |                                                         winter 01/15 |
 *----------------------------------------------------------------------*/
void XFEM::UTILS::SafetyCheckMaterials(
    Teuchos::RCP<Core::Mat::Material>& pmat, Teuchos::RCP<Core::Mat::Material>& nmat)
{
  //------------------------------ see whether materials in patch are equal

  if (pmat->MaterialType() != nmat->MaterialType())
    FOUR_C_THROW(" not the same material for master and slave parent element");

  if (pmat->MaterialType() == Core::Materials::m_matlist)
    FOUR_C_THROW(
        "A matlist has been found in edge based stabilization! If you are running XTPF, check "
        "calls as this should NOT happen!!!");

  if (pmat->MaterialType() != Core::Materials::m_carreauyasuda &&
      pmat->MaterialType() != Core::Materials::m_modpowerlaw &&
      pmat->MaterialType() != Core::Materials::m_herschelbulkley &&
      pmat->MaterialType() != Core::Materials::m_fluid)
    FOUR_C_THROW("Material law for parent element is not a fluid");

  if (pmat->MaterialType() == Core::Materials::m_fluid)
  {
    {
      const Mat::NewtonianFluid* actmat_p = static_cast<const Mat::NewtonianFluid*>(pmat.get());
      const double pvisc = actmat_p->Viscosity();
      const double pdens = actmat_p->Density();

      const Mat::NewtonianFluid* actmat_m = static_cast<const Mat::NewtonianFluid*>(nmat.get());
      const double nvisc = actmat_m->Viscosity();
      const double ndens = actmat_m->Density();

      if (std::abs(nvisc - pvisc) > 1e-14)
      {
        std::cout << "Parent element viscosity: " << pvisc
                  << " ,neighbor element viscosity: " << nvisc << std::endl;
        FOUR_C_THROW("parent and neighbor element do not have the same viscosity!");
      }
      if (std::abs(ndens - pdens) > 1e-14)
      {
        std::cout << "Parent element density: " << pdens << " ,neighbor element density: " << ndens
                  << std::endl;
        FOUR_C_THROW("parent and neighbor element do not have the same density!");
      }
    }
  }
  else
  {
    FOUR_C_THROW("up to now I expect a FLUID (m_fluid) material for edge stabilization\n");
  }

  return;
}

//! Extract a quantity for an element
void XFEM::UTILS::ExtractQuantityAtElement(Core::LinAlg::SerialDenseMatrix::Base& element_vector,
    const Core::Elements::Element* element,
    const Teuchos::RCP<const Epetra_MultiVector>& global_col_vector,
    Teuchos::RCP<Core::FE::Discretization>& dis, const int nds_vector, const int nsd)
{
  // get the other nds-set which is connected to the current one via this boundary-cell
  Core::Elements::Element::LocationArray la(dis->NumDofSets());
  element->LocationVector(*dis, la, false);

  const size_t numnode = element->num_node();

  if (la[nds_vector].lm_.size() != numnode)
  {
    std::cout << "la[nds_vector].lm_.size(): " << la[nds_vector].lm_.size() << std::endl;
    FOUR_C_THROW("assume a unique level-set dof in cutterdis-Dofset per node");
  }

  std::vector<double> local_vector(nsd * numnode);
  Core::FE::ExtractMyValues(*global_col_vector, local_vector, la[nds_vector].lm_);

  if (local_vector.size() != nsd * numnode)
    FOUR_C_THROW("wrong size of (potentially resized) local matrix!");

  // copy local to normal....
  Core::LinAlg::copy(local_vector.data(), element_vector);
}


//! Extract a quantity for a node
void XFEM::UTILS::ExtractQuantityAtNode(Core::LinAlg::SerialDenseMatrix::Base& element_vector,
    const Core::Nodes::Node* node, const Teuchos::RCP<const Epetra_MultiVector>& global_col_vector,
    Teuchos::RCP<Core::FE::Discretization>& dis, const int nds_vector, const unsigned int nsd)
{
  const std::vector<int> lm = dis->Dof(nds_vector, node);
  if (lm.size() != 1) FOUR_C_THROW("assume a unique level-set dof in cutterdis-Dofset");

  std::vector<double> local_vector(nsd);
  Core::FE::ExtractMyValues(*global_col_vector, local_vector, lm);

  if (local_vector.size() != nsd) FOUR_C_THROW("wrong size of (potentially resized) local matrix!");

  // copy local to nvec....
  Core::LinAlg::copy(local_vector.data(), element_vector);
}

FOUR_C_NAMESPACE_CLOSE
