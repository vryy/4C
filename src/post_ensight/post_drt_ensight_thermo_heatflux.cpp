/*----------------------------------------------------------------------*/
/*!
  \file post_drt_ensight_thermo_heatflux.cpp

  \brief postprocessing of thermal heatfluxes

  <pre>
  Maintainer: Caroline Danowski
              danowski@lnm.mw.tum.de
              http://www.lnm.mw.tum.de
  </pre>
*/

/*----------------------------------------------------------------------*
 | definitions                                               dano 11/09 |
 *----------------------------------------------------------------------*/


/*----------------------------------------------------------------------*
 | headers                                                   dano 11/09 |
 *----------------------------------------------------------------------*/
#include "post_drt_ensight_writer.H"
#include "../post_drt_common/post_drt_common.H"
#include <string>
#include "post_drt_ensight_single_field_writers.H"
#include "../linalg/linalg_utils.H"
#include "../pss_full/pss_cpp.h"
#include "../drt_thermo/thermo_ele_action.H"

/*----------------------------------------------------------------------*
 | constructor                                               dano 11/09 |
 *----------------------------------------------------------------------*/
void ThermoEnsightWriter::PostHeatflux(
  const std::string groupname,
  const std::string heatfluxtype
  )
{
  PostResult result = PostResult(field_);
  result.next_result();

  if (!map_has_map(result.group(), groupname.c_str()))
    return;

  //--------------------------------------------------------------------
  // calculation and output of nodal heatfluxes in xyz-(cartesian)-reference frame
  //--------------------------------------------------------------------
  if (heatfluxtype == "ndxyz")
  {
    WriteNodalHeatflux(groupname, result);
  }

  //-------------------------------------------------------------------------
  // calculation and output of element center heatfluxes in xyz-reference frame
  //-------------------------------------------------------------------------
  else if (heatfluxtype == "cxyz")
  {
    WriteElementCenterHeatflux(groupname, result);
  }

  //-----------------------------------------------------------------------------------
  // calculation and output of nodal and element center heatfluxes in xyz-reference frame
  //-----------------------------------------------------------------------------------
  else if (heatfluxtype == "cxyz_ndxyz")
  {
    WriteNodalHeatflux(groupname, result);

    // reset result for postprocessing and output of element center heatfluxes
    PostResult resulteleheatflux = PostResult(field_);
    resulteleheatflux.next_result();
    WriteElementCenterHeatflux(groupname,resulteleheatflux);
  }
  else
  {
    dserror("Unknown heatflux/tempgrad type");
  }

  return;
} // ThermoEnsightWriter::PostHeatflux


/*----------------------------------------------------------------------*
 | write nodal heatflux or temperature gradient              dano 11/09 |
 *----------------------------------------------------------------------*/
void ThermoEnsightWriter::WriteNodalHeatflux(
  const std::string groupname,
  PostResult& result
  )
{
  std::string name;
  std::string out;

  if (groupname == "gauss_initial_heatfluxes_xyz")
  {
    name = "nodal_initial_heatfluxes_xyz";
    out = "Initial heatfluxes";
  }
  else if (groupname == "gauss_current_heatfluxes_xyz")
  {
    name = "nodal_current_heatfluxes_xyz";
    out = "Current heatfluxes";
  }
  else if (groupname == "gauss_initial_tempgrad_xyz")
  {
    name = "nodal_initial_tempgrad_xyz";
    out = "Initial temperature gradients";
  }
  else if (groupname == "gauss_current_tempgrad_xyz")
  {
    name = "nodal_current_tempgrad_xyz";
    out = "Current temperature gradients";
  }
  else
  {
    dserror("trying to write something that is not a heatflux or a temperature gradient");
    exit(1);
  }

  // new for file continuation
  bool multiple_files = false;

  // open file
  const std::string filename = filename_ + "_"+ field_->name() + "."+ name;
  std::ofstream file;
  int startfilepos = 0;
  if (myrank_ == 0)
  {
    file.open(filename.c_str());
    startfilepos = file.tellp(); // file position should be zero, but we stay flexible
  }

  std::map<std::string, std::vector<std::ofstream::pos_type> > resultfilepos;
  int stepsize = 0;

  if (myrank_ == 0)
    std::cout << "writing node-based " << out << std::endl;

  // store information for later case file creation
  variableresulttypemap_[name] = "node";

  // number of heatflux/tempgrad DOFs
  int numdf = -1;
  if (field_->problem()->num_dim() == 3) numdf = 3;
  else if (field_->problem()->num_dim() == 2) numdf = 2;
  else if (field_->problem()->num_dim() == 1) numdf = 1;
  else dserror("Cannot handle dimension %g", field_->problem()->num_dim());

  WriteNodalHeatfluxStep(file,result,resultfilepos,groupname,name,numdf);
  // how many bits are necessary per time step (we assume a fixed size)?
  if (myrank_ == 0)
  {
    stepsize = ((int) file.tellp()) - startfilepos;
    if (stepsize <= 0) dserror("found invalid step size for result file");
  }
  else
    stepsize = 1; //use dummy value on other procs

  while (result.next_result())
  {
    const int indexsize = 80+2*sizeof(int)+(file.tellp()/stepsize+2)*sizeof(long);
    if (static_cast<long unsigned int>(file.tellp())+stepsize+indexsize>= FILE_SIZE_LIMIT_)
    {
      // multi file output, because its too much for single output file
      FileSwitcher(file, multiple_files, filesetmap_, resultfilepos, stepsize, name, filename);
    }

    WriteNodalHeatfluxStep(file,result,resultfilepos,groupname,name,numdf);
  }
  // store information for later case file creation
  filesetmap_[name].push_back(file.tellp()/stepsize);// has to be done BEFORE writing the index table
  variablenumdfmap_[name] = numdf;
  variablefilenamemap_[name] = filename;
  // store solution times vector for later case file creation
  {
    PostResult res = PostResult(field_); // this is needed!
    std::vector<double> restimes = res.get_result_times(field_->name(),groupname);
    timesetmap_[name] = restimes;
  }

  // append index table
  WriteIndexTable(file, resultfilepos[name]);
  resultfilepos[name].clear();

  // close result file
  if (file.is_open())
    file.close();

  return;
} // ThermoEnsightWriter::WriteNodalHeatflux


/*----------------------------------------------------------------------*
 |  write nodal heatfluxes                                   dano 11/09 |
 *----------------------------------------------------------------------*/
void ThermoEnsightWriter::WriteNodalHeatfluxStep(
  std::ofstream& file,
  PostResult& result,
  std::map<std::string, std::vector<std::ofstream::pos_type> >& resultfilepos,
  const std::string groupname,
  const std::string name,
  const int numdf
  ) const
{
  //--------------------------------------------------------------------
  // calculate nodal heatfluxes from gauss point heatfluxes
  //--------------------------------------------------------------------
  const Teuchos::RCP<std::map<int,Teuchos::RCP<Epetra_SerialDenseMatrix> > > data
    = result.read_result_serialdensematrix(groupname);

  const Teuchos::RCP<DRT::Discretization> dis = field_->discretization();

  // create the parameters for the discretization
  Teuchos::ParameterList p;
  // other parameters that might be needed by the elements
  p.set<int>("action",THR::postproc_thermo_heatflux);
  p.set("heatfluxtype","ndxyz");
  p.set("gpheatfluxmap", data);
  p.set("total time", -1.0);

  // create heatfluxes, three scalarvalued vectors
  Teuchos::RCP<Epetra_Vector> heatfluxx = LINALG::CreateVector(*(dis->DofRowMap()),true);
  Teuchos::RCP<Epetra_Vector> heatfluxy = LINALG::CreateVector(*(dis->DofRowMap()),true);
  Teuchos::RCP<Epetra_Vector> heatfluxz = LINALG::CreateVector(*(dis->DofRowMap()),true);
  dis->Evaluate(p,Teuchos::null,Teuchos::null,heatfluxx,heatfluxy,heatfluxz);

  // change the dis from a DofRowMap to a NodeRowMap, because Paraview can only visualize nodebased date
  const Epetra_Map* nodemap = dis->NodeRowMap();
  Teuchos::RCP<Epetra_MultiVector> nodal_heatfluxes = Teuchos::rcp(new Epetra_MultiVector(*nodemap, numdf));

  const int numnodes = dis->NumMyRowNodes();
  const unsigned numdofpernode = 1;

  if (numdf == 3) // 3 heatflux terms per node in 3D
  {
    for (int i=0; i<numnodes; ++i)
    {
      const DRT::Node* lnode = dis->lRowNode(i);
      const std::vector<int> lnodedofs = dis->Dof(lnode);

      if (lnodedofs.size() < numdofpernode)
        dserror("Too few DOFs at node of interest");
      const int adjele = lnode->NumElement();
      // build three scalar valued vectors for the heatflux output
      (*((*nodal_heatfluxes)(0)))[i] = (*heatfluxx)[dis->DofRowMap()->LID(lnodedofs[0])]/adjele;
      (*((*nodal_heatfluxes)(1)))[i] = (*heatfluxy)[dis->DofRowMap()->LID(lnodedofs[0])]/adjele;
      (*((*nodal_heatfluxes)(2)))[i] = (*heatfluxz)[dis->DofRowMap()->LID(lnodedofs[0])]/adjele;
    }
  }
  else if (numdf == 2) // 2 heatflux entries per node  in 2D
  {
    for (int i=0; i<numnodes; ++i)
    {
      const DRT::Node* lnode = dis->lRowNode(i);
      const std::vector<int> lnodedofs = dis->Dof(lnode);

      if (lnodedofs.size() < numdofpernode)
        dserror("Too few DOFs at node of interest");
      const int adjele = lnode->NumElement();
      // build two scalar valued vectors for the heatflux output
      (*((*nodal_heatfluxes)(0)))[i] = (*heatfluxx)[dis->DofRowMap()->LID(lnodedofs[0])]/adjele;
      (*((*nodal_heatfluxes)(1)))[i] = (*heatfluxy)[dis->DofRowMap()->LID(lnodedofs[0])]/adjele;
    }
  }
  else if (numdf == 1) // 1 heatflux entry per node  in 1D
  {
    for (int i=0; i<numnodes; ++i)
    {
      const DRT::Node* lnode = dis->lRowNode(i);
      const std::vector<int> lnodedofs = dis->Dof(lnode);

      if (lnodedofs.size() < numdofpernode)
        dserror("Too few DOFs at node of interest");
      const int adjele = lnode->NumElement();
      // build one scalar valued vectors for the heatflux output
      (*((*nodal_heatfluxes)(0)))[i] = (*heatfluxx)[dis->DofRowMap()->LID(lnodedofs[0])]/adjele;
    }
  }
  else
  {
    dserror("Cannot handle numdf=%g", numdf);
  }

  const Epetra_BlockMap& datamap = nodal_heatfluxes->Map();
  // contract Epetra_MultiVector on proc0 (proc0 gets everything, other procs empty)
  Teuchos::RCP<Epetra_MultiVector> data_proc0 = Teuchos::rcp(new Epetra_MultiVector(*proc0map_,numdf));
  Epetra_Import proc0dofimporter(*proc0map_,datamap);
  int err = data_proc0->Import(*nodal_heatfluxes,proc0dofimporter,Insert);
  if (err>0)
    dserror("Importing everything to proc 0 went wrong. Import returns %d",err);

  //--------------------------------------------------------------------
  // write some key words
  //--------------------------------------------------------------------
  std::vector<std::ofstream::pos_type>& filepos = resultfilepos[name];
  Write(file, "BEGIN TIME STEP");
  filepos.push_back(file.tellp());
  Write(file, "description");
  Write(file, "part");
  Write(file, field_->field_pos()+1);
  Write(file, "coordinates");

  //--------------------------------------------------------------------
  // write results
  //--------------------------------------------------------------------
  const int finalnumnode = proc0map_->NumGlobalElements();

  if (myrank_ == 0) // ensures pointer dofgids is valid
  {
    for (int idf=0; idf<numdf; ++idf)
    {
      for (int inode=0; inode<finalnumnode; inode++) // inode == lid of node because we use proc0map_
      {
        Write(file, static_cast<float>((*((*data_proc0)(idf)))[inode]));
      }
    }
  } // if (myrank_==0)

  Write(file, "END TIME STEP");
  return;
} // ThermoEnsightWriter::WriteNodalHeatfluxStep


/*----------------------------------------------------------------------*
 | write the output at the element center                    dano 11/09 |
 *----------------------------------------------------------------------*/
void ThermoEnsightWriter::WriteElementCenterHeatflux(
  const std::string groupname,
  PostResult& result
  )
{
  std::string name;
  std::string out;

  // output in Cartesian coordinates
  if (groupname == "gauss_initial_heatfluxes_xyz")
  {
    name = "element_initial_heatfluxes_xyz";
    out = "Initial heatfluxes";
  }
  else if (groupname == "gauss_current_heatfluxes_xyz")
  {
    name = "element_current_heatfluxes_xyz";
    out = "Current heatfluxes";
  }
  else if (groupname == "gauss_initial_tempgrad_xyz")
  {
    name = "element_initial_tempgrad_xyz";
    out = "Initial temperature gradients";
  }
  else if (groupname == "gauss_current_tempgrad_xyz")
  {
    name = "element_current_tempgrad_xyz";
    out = "Current temperature gradients";
  }
  else
  {
    dserror("trying to write something that is not a heatflux or a temperature gradient");
    exit(1);
  }

  // new for file continuation
  bool multiple_files = false;

  // open file
  const std::string filename = filename_ + "_"+ field_->name() + "."+ name;
  std::ofstream file;
  int startfilepos = 0;
  if (myrank_ == 0)
  {
    file.open(filename.c_str());
    startfilepos = file.tellp(); // file position should be zero, but we stay flexible
  }

  std::map<std::string, std::vector<std::ofstream::pos_type> > resultfilepos;
  int stepsize = 0;

  if (myrank_ == 0)
    std::cout << "writing element-based center " << out << std::endl;

  // store information for later case file creation
  variableresulttypemap_[name] = "element";

  // differ the cases of dimensions: 3 entries of heatflux in 3D, 2 in 2D, 1 in 1D
  int numdf = -1;
  if (field_->problem()->num_dim()==3) numdf = 3;
  else if (field_->problem()->num_dim()==2) numdf = 2;
  else if (field_->problem()->num_dim()==1) numdf = 1;
  else dserror("Cannot handle dimension %g", field_->problem()->num_dim());
  WriteElementCenterHeatfluxStep(file,result,resultfilepos,groupname,name,numdf);

  // how many bits are necessary per time step (we assume a fixed size)?
  if (myrank_ == 0)
  {
    stepsize = ((int) file.tellp())-startfilepos;
    if (stepsize <= 0) dserror("found invalid step size for result file");
  }
  else
  {
    stepsize = 1; //use dummy value on other procs
  }

  while (result.next_result())
  {
    const int indexsize = 80+2*sizeof(int)+(file.tellp()/stepsize+2)*sizeof(long);
    if (static_cast<long unsigned int>(file.tellp())+stepsize+indexsize>= FILE_SIZE_LIMIT_)
    {
      FileSwitcher(file, multiple_files, filesetmap_, resultfilepos, stepsize, name, filename);
    }
    WriteElementCenterHeatfluxStep(file,result,resultfilepos,groupname,name,numdf);
  }
  // store information for later case file creation
  filesetmap_[name].push_back(file.tellp()/stepsize);// has to be done BEFORE writing the index table
  variablenumdfmap_[name] = numdf;
  variablefilenamemap_[name] = filename;
  // store solution times vector for later case file creation
  {
    PostResult res = PostResult(field_); // this is needed!
    std::vector<double> restimes = res.get_result_times(field_->name(),groupname);
    timesetmap_[name] = restimes;
  }

  // append index table
  WriteIndexTable(file, resultfilepos[name]);
  resultfilepos[name].clear();

  // close result file
  if (file.is_open())
    file.close();

  return;
} // ThermoEnsightWriter::WriteElementCenterHeatflux


/*----------------------------------------------------------------------*
 |  output at the center                                     dano 11/09 |
 *----------------------------------------------------------------------*/
void ThermoEnsightWriter::WriteElementCenterHeatfluxStep(
  std::ofstream& file,
  PostResult& result,
  std::map<std::string,
  std::vector<std::ofstream::pos_type> >& resultfilepos,
  const std::string groupname,
  const std::string name,
  const int numdf
  ) const
{
  //--------------------------------------------------------------------
  // calculate element center heatfluxes from gauss point heatfluxes
  //--------------------------------------------------------------------
  const Teuchos::RCP<DRT::Discretization> dis = field_->discretization();
  const Teuchos::RCP<std::map<int,Teuchos::RCP<Epetra_SerialDenseMatrix> > > data
    = result.read_result_serialdensematrix(groupname);
  // create the parameters for the discretization
  Teuchos::ParameterList p;
  // other parameters that might be needed by the elements
  p.set<int>("action",THR::postproc_thermo_heatflux);
  p.set("heatfluxtype","cxyz");
  p.set("gpheatfluxmap",data);
  p.set("total time", -1.0);
  Teuchos::RCP<Epetra_MultiVector> eleheatflux = Teuchos::rcp(new Epetra_MultiVector(*(dis->ElementRowMap()),numdf));
  p.set("eleheatflux",eleheatflux);
  dis->Evaluate(p,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null);
  if (eleheatflux == Teuchos::null)
  {
    dserror("vector containing element center heatfluxes/tempgradients not available");
  }

  //--------------------------------------------------------------------
  // write some key words
  //--------------------------------------------------------------------
  std::vector<std::ofstream::pos_type>& filepos = resultfilepos[name];
  Write(file, "BEGIN TIME STEP");
  filepos.push_back(file.tellp());
  Write(file, "description");
  Write(file, "part");
  Write(file, field_->field_pos()+1);

  const Epetra_BlockMap& datamap = eleheatflux->Map();

  // do stupid conversion into Epetra map
  Teuchos::RCP<Epetra_Map> epetradatamap;
  epetradatamap = Teuchos::rcp(new Epetra_Map(datamap.NumGlobalElements(),
                                              datamap.NumMyElements(),
                                              datamap.MyGlobalElements(),
                                              0,
                                              datamap.Comm()));

  Teuchos::RCP<Epetra_Map> proc0datamap;
  proc0datamap = LINALG::AllreduceEMap(*epetradatamap,0);
  // sort proc0datamap so that we can loop it and get nodes in ascending order.
  std::vector<int> sortmap;
  sortmap.reserve(proc0datamap->NumMyElements());
  sortmap.assign(proc0datamap->MyGlobalElements(), proc0datamap->MyGlobalElements()+proc0datamap->NumMyElements());
  std::sort(sortmap.begin(), sortmap.end());
  proc0datamap = Teuchos::rcp(new Epetra_Map(-1, sortmap.size(), &sortmap[0], 0, proc0datamap->Comm()));

  // contract Epetra_MultiVector on proc0 (proc0 gets everything, other procs empty)
  Teuchos::RCP<Epetra_MultiVector> data_proc0 = Teuchos::rcp(new Epetra_MultiVector(*proc0datamap,numdf));
  Epetra_Import proc0dofimporter(*proc0datamap,datamap);
  int err = data_proc0->Import(*eleheatflux,proc0dofimporter,Insert);
  if (err>0) dserror("Importing everything to proc 0 went wrong. Import returns %d",err);

  //--------------------------------------------------------------------
  // specify the element type
  //--------------------------------------------------------------------
  // loop over the different element types present
  if (myrank_ == 0)
  {
    if (eleGidPerDisType_.empty()==true) dserror("no element types available");
  }
  EleGidPerDisType::const_iterator iter;
  for (iter = eleGidPerDisType_.begin(); iter != eleGidPerDisType_.end(); ++iter)
  {
    const std::string ensighteleString = GetEnsightString(iter->first);
    const int numelepertype = (iter->second).size();
    std::vector<int> actelegids(numelepertype);
    actelegids = iter->second;
    // write element type
    Write(file, ensighteleString);

    //------------------------------------------------------------------
    // write results
    //------------------------------------------------------------------
    if (myrank_ == 0) // ensures pointer dofgids is valid
    {
      for (int idf=0; idf<numdf; ++idf)
      {
        for (int iele=0; iele<numelepertype; iele++) // inode == lid of node because we use proc0map_
        {
          // extract element global id
          const int gid = actelegids[iele];
          // get the dof local id w.r.t. the final datamap
          int lid = proc0datamap->LID(gid);
          if (lid > -1)
          {
            Write(file, static_cast<float>((*((*data_proc0)(idf)))[lid]));
          }
        }
      }
    } // if (myrank_ == 0)
  }

  Write(file, "END TIME STEP");

  return;
} // ThermoEnsightWriter::WriteElementCenterHeatfluxStep


/*----------------------------------------------------------------------*/
