/*----------------------------------------------------------------------*/
/*! \file

  \brief postprocessing of thermal heatfluxes


  \level 1
*/

/*----------------------------------------------------------------------*
 | definitions                                               dano 11/09 |
 *----------------------------------------------------------------------*/


/*----------------------------------------------------------------------*
 | headers                                                   dano 11/09 |
 *----------------------------------------------------------------------*/
#include "4C_io_legacy_table.hpp"
#include "4C_lib_discret.hpp"
#include "4C_linalg_utils_sparse_algebra_create.hpp"
#include "4C_post_common.hpp"
#include "4C_post_processor_single_field_writers.hpp"
#include "4C_post_writer_base.hpp"
#include "4C_thermo_ele_action.hpp"

#include <string>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 | constructor                                               dano 11/09 |
 *----------------------------------------------------------------------*/
void ThermoFilter::PostHeatflux(const std::string groupname, const std::string heatfluxtype)
{
  PostField* field = writer_->GetField();
  PostResult result = PostResult(field);
  result.next_result();

  if (!map_has_map(result.group(), groupname.c_str())) return;

  //--------------------------------------------------------------------
  // calculation and output of nodal heatfluxes in xyz-(cartesian)-reference frame
  //--------------------------------------------------------------------
  if (heatfluxtype == "ndxyz")
  {
    WriteHeatflux(groupname, result, nodebased);
  }

  //-------------------------------------------------------------------------
  // calculation and output of element center heatfluxes in xyz-reference frame
  //-------------------------------------------------------------------------
  else if (heatfluxtype == "cxyz")
  {
    WriteHeatflux(groupname, result, elementbased);
  }

  //-----------------------------------------------------------------------------------
  // calculation and output of nodal and element center heatfluxes in xyz-reference frame
  //-----------------------------------------------------------------------------------
  else if (heatfluxtype == "cxyz_ndxyz")
  {
    WriteHeatflux(groupname, result, nodebased);

    // reset result for postprocessing and output of element center heatfluxes
    PostResult resulteleheatflux = PostResult(field);
    resulteleheatflux.next_result();
    WriteHeatflux(groupname, resulteleheatflux, elementbased);
  }
  else
  {
    FOUR_C_THROW("Unknown heatflux/tempgrad type");
  }

}  // ThermoFilter::PostHeatflux



/*----------------------------------------------------------------------*
 |  write nodal heatfluxes                                   dano 11/09 |
 *----------------------------------------------------------------------*/
struct WriteNodalHeatfluxStep : SpecialFieldInterface
{
  WriteNodalHeatfluxStep(ThermoFilter& filter) : filter_(filter) {}

  int Numdf()
  {
    int numdf = -1;
    if (filter_.GetWriter().GetField()->problem()->num_dim() == 3)
      numdf = 3;
    else if (filter_.GetWriter().GetField()->problem()->num_dim() == 2)
      numdf = 2;
    else if (filter_.GetWriter().GetField()->problem()->num_dim() == 1)
      numdf = 1;
    else
      FOUR_C_THROW(
          "Cannot handle dimension %g", filter_.GetWriter().GetField()->problem()->num_dim());
    return numdf;
  }

  std::vector<int> NumDfMap() override { return std::vector<int>(1, Numdf()); }

  void operator()(std::vector<Teuchos::RCP<std::ofstream>>& files, PostResult& result,
      std::map<std::string, std::vector<std::ofstream::pos_type>>& resultfilepos,
      const std::string& groupname, const std::vector<std::string>& name) override
  {
    using namespace FourC;

    FOUR_C_ASSERT(name.size() == 1, "Unexpected number of names");

    int numdf = Numdf();

    //--------------------------------------------------------------------
    // calculate nodal heatfluxes from gauss point heatfluxes
    //--------------------------------------------------------------------
    const Teuchos::RCP<std::map<int, Teuchos::RCP<CORE::LINALG::SerialDenseMatrix>>> data =
        result.read_result_serialdensematrix(groupname);

    const Teuchos::RCP<DRT::Discretization> dis = result.field()->discretization();

    // create the parameters for the discretization
    Teuchos::ParameterList p;
    // other parameters that might be needed by the elements
    p.set<int>("action", THR::postproc_thermo_heatflux);
    p.set("heatfluxtype", "ndxyz");
    p.set("gpheatfluxmap", data);
    p.set("total time", -1.0);

    // create heatfluxes, three scalarvalued vectors
    Teuchos::RCP<Epetra_Vector> heatfluxx = CORE::LINALG::CreateVector(*(dis->DofRowMap()), true);
    Teuchos::RCP<Epetra_Vector> heatfluxy = CORE::LINALG::CreateVector(*(dis->DofRowMap()), true);
    Teuchos::RCP<Epetra_Vector> heatfluxz = CORE::LINALG::CreateVector(*(dis->DofRowMap()), true);
    dis->Evaluate(p, Teuchos::null, Teuchos::null, heatfluxx, heatfluxy, heatfluxz);

    // change the dis from a DofRowMap to a NodeRowMap, because Paraview can only visualize
    // nodebased date
    const Epetra_Map* nodemap = dis->NodeRowMap();
    Teuchos::RCP<Epetra_MultiVector> nodal_heatfluxes =
        Teuchos::rcp(new Epetra_MultiVector(*nodemap, numdf));

    const int numnodes = dis->NumMyRowNodes();
    const unsigned numdofpernode = 1;

    if (numdf == 3)  // 3 heatflux terms per node in 3D
    {
      for (int i = 0; i < numnodes; ++i)
      {
        const DRT::Node* lnode = dis->lRowNode(i);
        const std::vector<int> lnodedofs = dis->Dof(lnode);

        if (lnodedofs.size() < numdofpernode) FOUR_C_THROW("Too few DOFs at node of interest");
        const int adjele = lnode->NumElement();
        // build three scalar valued vectors for the heatflux output
        (*((*nodal_heatfluxes)(0)))[i] = (*heatfluxx)[dis->DofRowMap()->LID(lnodedofs[0])] / adjele;
        (*((*nodal_heatfluxes)(1)))[i] = (*heatfluxy)[dis->DofRowMap()->LID(lnodedofs[0])] / adjele;
        (*((*nodal_heatfluxes)(2)))[i] = (*heatfluxz)[dis->DofRowMap()->LID(lnodedofs[0])] / adjele;
      }
    }
    else if (numdf == 2)  // 2 heatflux entries per node  in 2D
    {
      for (int i = 0; i < numnodes; ++i)
      {
        const DRT::Node* lnode = dis->lRowNode(i);
        const std::vector<int> lnodedofs = dis->Dof(lnode);

        if (lnodedofs.size() < numdofpernode) FOUR_C_THROW("Too few DOFs at node of interest");
        const int adjele = lnode->NumElement();
        // build two scalar valued vectors for the heatflux output
        (*((*nodal_heatfluxes)(0)))[i] = (*heatfluxx)[dis->DofRowMap()->LID(lnodedofs[0])] / adjele;
        (*((*nodal_heatfluxes)(1)))[i] = (*heatfluxy)[dis->DofRowMap()->LID(lnodedofs[0])] / adjele;
      }
    }
    else if (numdf == 1)  // 1 heatflux entry per node  in 1D
    {
      for (int i = 0; i < numnodes; ++i)
      {
        const DRT::Node* lnode = dis->lRowNode(i);
        const std::vector<int> lnodedofs = dis->Dof(lnode);

        if (lnodedofs.size() < numdofpernode) FOUR_C_THROW("Too few DOFs at node of interest");
        const int adjele = lnode->NumElement();
        // build one scalar valued vectors for the heatflux output
        (*((*nodal_heatfluxes)(0)))[i] = (*heatfluxx)[dis->DofRowMap()->LID(lnodedofs[0])] / adjele;
      }
    }
    else
    {
      FOUR_C_THROW("Cannot handle numdf=%g", numdf);
    }

    filter_.GetWriter().WriteNodalResultStep(
        *files[0], nodal_heatfluxes, resultfilepos, groupname, name[0], numdf);
  }

  ThermoFilter& filter_;
};



/*----------------------------------------------------------------------*
 |  write nodal heatfluxes                                   dano 11/09 |
 *----------------------------------------------------------------------*/
struct WriteElementCenterHeatfluxStep : SpecialFieldInterface
{
  WriteElementCenterHeatfluxStep(ThermoFilter& filter) : filter_(filter) {}

  int Numdf()
  {
    int numdf = -1;
    if (filter_.GetWriter().GetField()->problem()->num_dim() == 3)
      numdf = 3;
    else if (filter_.GetWriter().GetField()->problem()->num_dim() == 2)
      numdf = 2;
    else if (filter_.GetWriter().GetField()->problem()->num_dim() == 1)
      numdf = 1;
    else
      FOUR_C_THROW(
          "Cannot handle dimension %g", filter_.GetWriter().GetField()->problem()->num_dim());
    return numdf;
  }

  std::vector<int> NumDfMap() override { return std::vector<int>(1, Numdf()); }

  void operator()(std::vector<Teuchos::RCP<std::ofstream>>& files, PostResult& result,
      std::map<std::string, std::vector<std::ofstream::pos_type>>& resultfilepos,
      const std::string& groupname, const std::vector<std::string>& name) override
  {
    using namespace FourC;

    FOUR_C_ASSERT(name.size() == 1, "Unexpected number of names");

    int numdf = Numdf();

    //--------------------------------------------------------------------
    // calculate element center heatfluxes from gauss point heatfluxes
    //--------------------------------------------------------------------
    const Teuchos::RCP<DRT::Discretization> dis = result.field()->discretization();
    const Teuchos::RCP<std::map<int, Teuchos::RCP<CORE::LINALG::SerialDenseMatrix>>> data =
        result.read_result_serialdensematrix(groupname);
    // create the parameters for the discretization
    Teuchos::ParameterList p;
    // other parameters that might be needed by the elements
    p.set<int>("action", THR::postproc_thermo_heatflux);
    p.set("heatfluxtype", "cxyz");
    p.set("gpheatfluxmap", data);
    p.set("total time", -1.0);
    Teuchos::RCP<Epetra_MultiVector> eleheatflux =
        Teuchos::rcp(new Epetra_MultiVector(*(dis->ElementRowMap()), numdf));
    p.set("eleheatflux", eleheatflux);
    dis->Evaluate(p, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null);
    if (eleheatflux == Teuchos::null)
    {
      FOUR_C_THROW("vector containing element center heatfluxes/tempgradients not available");
    }

    filter_.GetWriter().WriteElementResultStep(
        *files[0], eleheatflux, resultfilepos, groupname, name[0], numdf, 0);
  }

  ThermoFilter& filter_;
};



/*----------------------------------------------------------------------*
 | write nodal heatflux or temperature gradient              dano 11/09 |
 *----------------------------------------------------------------------*/
void ThermoFilter::WriteHeatflux(
    const std::string groupname, PostResult& result, const ResultType kind)
{
  std::string name;
  std::string out;

  if (groupname == "gauss_initial_heatfluxes_xyz")
  {
    name = "initial_heatfluxes_xyz";
    out = "Initial heatfluxes";
  }
  else if (groupname == "gauss_current_heatfluxes_xyz")
  {
    name = "current_heatfluxes_xyz";
    out = "Current heatfluxes";
  }
  else if (groupname == "gauss_initial_tempgrad_xyz")
  {
    name = "initial_tempgrad_xyz";
    out = "Initial temperature gradients";
  }
  else if (groupname == "gauss_current_tempgrad_xyz")
  {
    name = "current_tempgrad_xyz";
    out = "Current temperature gradients";
  }
  else
  {
    FOUR_C_THROW("trying to write something that is not a heatflux or a temperature gradient");
    exit(1);
  }

  if (kind == nodebased)
  {
    name = "nodal_" + name;
    WriteNodalHeatfluxStep heatflux(*this);
    writer_->WriteSpecialField(
        heatflux, result, nodebased, groupname, std::vector<std::string>(1, name), out);
  }
  else if (kind == elementbased)
  {
    name = "element_" + name;
    WriteElementCenterHeatfluxStep heatflux(*this);
    writer_->WriteSpecialField(
        heatflux, result, elementbased, groupname, std::vector<std::string>(1, name), out);
  }
  else
    FOUR_C_THROW("Unknown heatflux type");
}  // ThermoFilter::WriteNodalHeatflux

FOUR_C_NAMESPACE_CLOSE

/*----------------------------------------------------------------------*/
