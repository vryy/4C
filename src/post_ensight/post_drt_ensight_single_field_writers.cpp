/*!
  \file post_drt_ensight_single_field_writers.cpp

  \brief main routine of the Ensight filter

  <pre>
  Maintainer: Axel Gerstenberger
  gerstenberger@lnm.mw.tum.de
  http://www.lnm.mw.tum.de/Members/gerstenberger
  089 - 289-15236
  </pre>

*/

#ifdef CCADISCRET

#include "post_drt_ensight_writer.H"
#include <string>
#include "post_drt_ensight_single_field_writers.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_xfem/integrationcell.H"
#include "../drt_xfem/physics.H"
#include "../drt_xfem/dof_management.H"
#include "../drt_xfem/enrichment_utils.H"
#include "../drt_f3/xfluid3_interpolation.H" // obviosly, this is fluid element specific and needs more generalization

using namespace std;


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void StructureEnsightWriter::WriteAllResults(PostField* field)
{
  EnsightWriter::WriteResult("displacement", "displacement", dofbased, field->problem()->num_dim());
  EnsightWriter::WriteResult("velocity", "velocity", dofbased, field->problem()->num_dim());
  EnsightWriter::WriteResult("acceleration", "acceleration", dofbased, field->problem()->num_dim());
  EnsightWriter::WriteElementResults(field);
  if (stresstype_!="none")
  {
    // although appearing here twice, only one function call to PostStress
    // is really postprocessing Gauss point stresses, since only _either_
    // Cauchy _or_ 2nd Piola-Kirchhoff stresses are written during simulation!
    PostStress("gauss_cauchy_stresses_xyz", stresstype_);
    PostStress("gauss_2PK_stresses_xyz", stresstype_);
  }
  if (straintype_!="none")
  {
    PostStress("gauss_GL_strains_xyz", straintype_);
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void StructureEnsightWriter::PostStress(const string groupname, const string stresstype)
{
  PostResult result = PostResult(field_);
  result.next_result();

  if (!map_has_map(result.group(), groupname.c_str()))
    return;

  //--------------------------------------------------------------------
  // calculation and output of nodal stresses in xyz-reference frame
  //--------------------------------------------------------------------

  if (stresstype == "ndxyz")
  {
    WriteNodalStress(groupname, result);
  }

  //-------------------------------------------------------------------------
  // calculation and output of element center stresses in xyz-reference frame
  //-------------------------------------------------------------------------

  else if (stresstype == "cxyz")
  {
    WriteElementCenterStress(groupname, result);
  }

  //-----------------------------------------------------------------------------------
  // calculation and output of nodal and element center stresses in xyz-reference frame
  //-----------------------------------------------------------------------------------

  else if (stresstype == "cxyz_ndxyz")
  {
    WriteNodalStress(groupname, result);

    // reset result for postprocessing and output of element center stresses
    PostResult resultelestress = PostResult(field_);
    resultelestress.next_result();
    WriteElementCenterStress(groupname,resultelestress);
  }

  else
  {
    dserror("Unknown stress/strain type");
  }

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void StructureEnsightWriter::WriteNodalStress(const string groupname,
                                              PostResult& result)
{
  string name;
  string out;

  if (groupname=="gauss_2PK_stresses_xyz")
  {
    name="nodal_2PK_stresses_xyz";
    out="2nd Piola-Kirchhoff stresses";
  }
  else if (groupname=="gauss_cauchy_stresses_xyz")
  {
    name="nodal_cauchy_stresses_xyz";
    out="Cauchy stresses";
  }
  else if (groupname=="gauss_GL_strains_xyz")
  {
    name="nodal_GL_strains_xyz";
    out="Green-Lagrange strains";
  }
  else
  {
    dserror("trying to write something that is not a stress or a strain");
    exit(1);            // especially for Axel ;-)
  }

  // new for file continuation
  bool multiple_files = false;

  // open file
  const string filename = filename_ + "_"+ field_->name() + "."+ name;
  ofstream file;
  int startfilepos = 0;
  if (myrank_==0)
  {
    file.open(filename.c_str());
    startfilepos = file.tellp(); // file position should be zero, but we stay flexible
  }

  map<string, vector<ofstream::pos_type> > resultfilepos;
  int stepsize = 0;

  if (myrank_==0)
    cout<<"writing node-based " << out << endl;

  // store information for later case file creation
  variableresulttypemap_[name] = "node";

  WriteNodalStressStep(file,result,resultfilepos,groupname,name,6);
  // how many bits are necessary per time step (we assume a fixed size)?
  if (myrank_==0)
  {
    stepsize = ((int) file.tellp())-startfilepos;
    if (stepsize <= 0) dserror("found invalid step size for result file");
  }
  else
    stepsize = 1; //use dummy value on other procs

  while (result.next_result())
  {
    const int indexsize = 80+2*sizeof(int)+(file.tellp()/stepsize+2)*sizeof(long);
    if (static_cast<long unsigned int>(file.tellp())+stepsize+indexsize>= FILE_SIZE_LIMIT_)
    {
      FileSwitcher(file, multiple_files, filesetmap_, resultfilepos, stepsize, name, filename);
    }

    WriteNodalStressStep(file,result,resultfilepos,groupname,name,6);
  }
  // store information for later case file creation
  filesetmap_[name].push_back(file.tellp()/stepsize);// has to be done BEFORE writing the index table
  variablenumdfmap_[name] = 6;
  variablefilenamemap_[name] = filename;
  // store solution times vector for later case file creation
  {
    PostResult res = PostResult(field_); // this is needed!
    vector<double> restimes = res.get_result_times(field_->name(),groupname);
    timesetmap_[name] = restimes;
  }

  // append index table
  WriteIndexTable(file, resultfilepos[name]);
  resultfilepos[name].clear();

  // close result file
  if (file.is_open())
    file.close();

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void StructureEnsightWriter::WriteNodalStressStep(ofstream& file,
                                                  PostResult& result,
                                                  map<string, vector<ofstream::pos_type> >& resultfilepos,
                                                  const string groupname,
                                                  const string name,
                                                  const int numdf) const
{
  //--------------------------------------------------------------------
  // calculate nodal stresses from gauss point stresses
  //--------------------------------------------------------------------

  const RefCountPtr<std::map<int,RefCountPtr<Epetra_SerialDenseMatrix> > > data =
    result.read_result_serialdensematrix(groupname);

  const RefCountPtr<DRT::Discretization> dis = field_->discretization();

  ParameterList p;
  p.set("action","postprocess_stress");
  p.set("stresstype","ndxyz");
  p.set("gpstressmap", data);
  RefCountPtr<Epetra_Vector> normal_stresses = LINALG::CreateVector(*(dis->DofRowMap()),true);
  RefCountPtr<Epetra_Vector> shear_stresses = LINALG::CreateVector(*(dis->DofRowMap()),true);
  dis->Evaluate(p,null,null,normal_stresses,shear_stresses,null);

  const Epetra_Map* nodemap = dis->NodeRowMap();
  RefCountPtr<Epetra_MultiVector> nodal_stresses = rcp(new Epetra_MultiVector(*nodemap, 6));

  const int numnodes = dis->NumMyRowNodes();

  for (int i=0;i<numnodes;++i)
  {
    (*((*nodal_stresses)(0)))[i] = (*normal_stresses)[3*i];
    (*((*nodal_stresses)(1)))[i] = (*normal_stresses)[3*i+1];
    (*((*nodal_stresses)(2)))[i] = (*normal_stresses)[3*i+2];
    (*((*nodal_stresses)(3)))[i] = (*shear_stresses)[3*i];
    (*((*nodal_stresses)(4)))[i] = (*shear_stresses)[3*i+1];
    (*((*nodal_stresses)(5)))[i] = (*shear_stresses)[3*i+2];
  }

  //--------------------------------------------------------------------
  // write some key words
  //--------------------------------------------------------------------

  vector<ofstream::pos_type>& filepos = resultfilepos[name];
  Write(file, "BEGIN TIME STEP");
  filepos.push_back(file.tellp());
  Write(file, "description");
  Write(file, "part");
  Write(file, field_->field_pos()+1);
  Write(file, "coordinates");

  const Epetra_BlockMap& datamap = nodal_stresses->Map();

  // do stupid conversion into Epetra map
  RefCountPtr<Epetra_Map> epetradatamap;
  epetradatamap = rcp(new Epetra_Map(datamap.NumGlobalElements(),
                                     datamap.NumMyElements(),
                                     datamap.MyGlobalElements(),
                                     0,
                                     datamap.Comm()));



  // contract Epetra_MultiVector on proc0 (proc0 gets everything, other procs empty)
  RefCountPtr<Epetra_MultiVector> data_proc0 = rcp(new Epetra_MultiVector(*proc0map_,numdf));
  Epetra_Import proc0dofimporter(*proc0map_,datamap);
  int err = data_proc0->Import(*nodal_stresses,proc0dofimporter,Insert);
  if (err>0) dserror("Importing everything to proc 0 went wrong. Import returns %d",err);

  //--------------------------------------------------------------------
  // write results
  //--------------------------------------------------------------------

  const int finalnumnode = proc0map_->NumGlobalElements();

  if (myrank_==0) // ensures pointer dofgids is valid
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
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void StructureEnsightWriter::WriteElementCenterStress(const string groupname,
                                                      PostResult& result)
{
  string name;
  string out;

  if (groupname=="gauss_2PK_stresses_xyz")
  {
    name="element_2PK_stresses_xyz";
    out="2nd Piola-Kirchhoff stresses";
  }
  else if (groupname=="gauss_cauchy_stresses_xyz")
  {
    name="element_cauchy_stresses_xyz";
    out="Cauchy stresses";
  }
  else if (groupname=="gauss_GL_strains_xyz")
  {
    name="element_GL_strains_xyz";
    out="Green-Lagrange strains";
  }
  else
  {
    dserror("trying to write something that is not a stress or a strain");
    exit(1);            // especially for Axel ;-)
  }

  // new for file continuation
  bool multiple_files = false;

  // open file
  const string filename = filename_ + "_"+ field_->name() + "."+ name;
  ofstream file;
  int startfilepos = 0;
  if (myrank_==0)
  {
    file.open(filename.c_str());
    startfilepos = file.tellp(); // file position should be zero, but we stay flexible
  }

  map<string, vector<ofstream::pos_type> > resultfilepos;
  int stepsize = 0;

  if (myrank_==0)
    cout<<"writing element-based center " << out << endl;

  // store information for later case file creation
  variableresulttypemap_[name] = "element";

  WriteElementCenterStressStep(file,result,resultfilepos,groupname,name,6);

  // how many bits are necessary per time step (we assume a fixed size)?
  if (myrank_==0)
  {
    stepsize = ((int) file.tellp())-startfilepos;
    if (stepsize <= 0) dserror("found invalid step size for result file");
  }
  else
    stepsize = 1; //use dummy value on other procs

  while (result.next_result())
  {
    const int indexsize = 80+2*sizeof(int)+(file.tellp()/stepsize+2)*sizeof(long);
    if (static_cast<long unsigned int>(file.tellp())+stepsize+indexsize>= FILE_SIZE_LIMIT_)
    {
      FileSwitcher(file, multiple_files, filesetmap_, resultfilepos, stepsize, name, filename);
    }
    WriteElementCenterStressStep(file,result,resultfilepos,groupname,name,6);
  }
  // store information for later case file creation
  filesetmap_[name].push_back(file.tellp()/stepsize);// has to be done BEFORE writing the index table
  variablenumdfmap_[name] = 6;
  variablefilenamemap_[name] = filename;
  // store solution times vector for later case file creation
  {
    PostResult res = PostResult(field_); // this is needed!
    vector<double> restimes = res.get_result_times(field_->name(),groupname);
    timesetmap_[name] = restimes;
  }

  // append index table
  WriteIndexTable(file, resultfilepos[name]);
  resultfilepos[name].clear();

  // close result file
  if (file.is_open())
    file.close();

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void StructureEnsightWriter::WriteElementCenterStressStep(ofstream& file,
                                                          PostResult& result,
                                                          map<string, vector<ofstream::pos_type> >& resultfilepos,
                                                          const string groupname,
                                                          const string name,
                                                          const int numdf) const
{
  //--------------------------------------------------------------------
  // calculate element center stresses from gauss point stresses
  //--------------------------------------------------------------------

  const RefCountPtr<DRT::Discretization> dis = field_->discretization();
  const RefCountPtr<std::map<int,RefCountPtr<Epetra_SerialDenseMatrix> > > data =
    result.read_result_serialdensematrix(groupname);
  ParameterList p;
  p.set("action","postprocess_stress");
  p.set("stresstype","cxyz");
  p.set("gpstressmap", data);
  RCP<Epetra_MultiVector> elestress = rcp(new Epetra_MultiVector(*(dis->ElementRowMap()),6));
  p.set("elestress",elestress);
  dis->Evaluate(p,null,null,null,null,null);
  if (elestress==null)
  {
    dserror("vector containing element center stresses/strains not available");
  }

  //--------------------------------------------------------------------
  // write some key words
  //--------------------------------------------------------------------

  vector<ofstream::pos_type>& filepos = resultfilepos[name];
  Write(file, "BEGIN TIME STEP");
  filepos.push_back(file.tellp());
  Write(file, "description");
  Write(file, "part");
  Write(file, field_->field_pos()+1);

  const Epetra_BlockMap& datamap = elestress->Map();

  // do stupid conversion into Epetra map
  RefCountPtr<Epetra_Map> epetradatamap;
  epetradatamap = rcp(new Epetra_Map(datamap.NumGlobalElements(),
                                     datamap.NumMyElements(),
                                     datamap.MyGlobalElements(),
                                     0,
                                     datamap.Comm()));

  RefCountPtr<Epetra_Map> proc0datamap;
  proc0datamap = LINALG::AllreduceEMap(*epetradatamap,0);
  // sort proc0datamap so that we can loop it and get nodes in ascending order.
  std::vector<int> sortmap;
  sortmap.reserve(proc0datamap->NumMyElements());
  sortmap.assign(proc0datamap->MyGlobalElements(), proc0datamap->MyGlobalElements()+proc0datamap->NumMyElements());
  std::sort(sortmap.begin(), sortmap.end());
  proc0datamap = Teuchos::rcp(new Epetra_Map(-1, sortmap.size(), &sortmap[0], 0, proc0datamap->Comm()));

  // contract Epetra_MultiVector on proc0 (proc0 gets everything, other procs empty)
  RefCountPtr<Epetra_MultiVector> data_proc0 = rcp(new Epetra_MultiVector(*proc0datamap,numdf));
  Epetra_Import proc0dofimporter(*proc0datamap,datamap);
  int err = data_proc0->Import(*elestress,proc0dofimporter,Insert);
  if (err>0) dserror("Importing everything to proc 0 went wrong. Import returns %d",err);

  //--------------------------------------------------------------------
  // specify the element type
  //--------------------------------------------------------------------
  // loop over the different element types present
  if (myrank_==0)
  {
    if (eleGidPerDisType_.empty()==true) dserror("no element types available");
  }
  EleGidPerDisType::const_iterator iter;
  for (iter=eleGidPerDisType_.begin(); iter != eleGidPerDisType_.end(); ++iter)
  {
    const string ensighteleString = GetEnsightString(iter->first);
    const int numelepertype = (iter->second).size();
    vector<int> actelegids(numelepertype);
    actelegids = iter->second;
    // write element type
    Write(file, ensighteleString);

    //------------------------------------------------------------------
    // write results
    //------------------------------------------------------------------

    if (myrank_==0) // ensures pointer dofgids is valid
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
    } // if (myrank_==0)
  }
  Write(file, "END TIME STEP");
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FluidEnsightWriter::WriteAllResults(PostField* field)
{
  EnsightWriter::WriteResult("velnp", "velocity", dofbased, field->problem()->num_dim());
  EnsightWriter::WriteResult("velnp", "pressure", dofbased, 1, field->problem()->num_dim());
  EnsightWriter::WriteResult("residual", "residual", dofbased, field->problem()->num_dim());
  EnsightWriter::WriteResult("dispnp", "displacement", dofbased, field->problem()->num_dim());
  EnsightWriter::WriteResult("traction", "traction", dofbased, field->problem()->num_dim());
  WriteElementResults(field);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DGFEMFluidEnsightWriter::WriteAllResults(PostField* field)
{
  EnsightWriter::WriteResult("velnp", "elemean_velocity", elementdof, field->problem()->num_dim());
  EnsightWriter::WriteResult("velnp", "elemean_pressure", elementdof, 1, field->problem()->num_dim());
  WriteElementResults(field);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void AleEnsightWriter::WriteAllResults(PostField* field)
{
  EnsightWriter::WriteResult("dispnp", "displacement", dofbased, field->problem()->num_dim());
  WriteElementResults(field);
}


/*----------------------------------------------------------------------*
  |                                                             gjb 12/07  |
  \*----------------------------------------------------------------------*/
void ConDifEnsightWriter::WriteAllResults(PostField* field)
{
  //phinp is a scalar result field with ONE dof per node.
  //Therefore it is NOT possible to hand over field->problem()->num_dim()
  // (equals 2 or 3) as a number of dofs
  EnsightWriter::WriteResult("phinp", "phi", dofbased, 1);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void XFluidEnsightWriter::WriteAllResults(PostField* field)
{
//    EnsightWriter::WriteResult("velnp", "velocity", field->problem()->num_dim());
//    EnsightWriter::WriteResult("velnp", "pressure", 1, field->problem()->num_dim());
//    EnsightWriter::WriteResult("residual", "residual", field->problem()->num_dim());

  cout << "now to xfem solutions" << endl;
  set<XFEM::PHYSICS::Field> velocity_fieldset;
  velocity_fieldset.insert(XFEM::PHYSICS::Velx);
  velocity_fieldset.insert(XFEM::PHYSICS::Vely);
  velocity_fieldset.insert(XFEM::PHYSICS::Velz);
  set<XFEM::PHYSICS::Field> pressure_fieldset;
  pressure_fieldset.insert(XFEM::PHYSICS::Pres);

  XFluidEnsightWriter::WriteResult("velnp", "velocity_physical", velocity_fieldset);
  //XFluidEnsightWriter::WriteResult("residual", "residual_physical", velocity_fieldset);
  XFluidEnsightWriter::WriteResult("velnp", "pressure_physical", pressure_fieldset);

}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void XFluidEnsightWriter::WriteFiles()
{
#ifndef PARALLEL
  if (myrank_ > 0) dserror("have serial filter version, but myrank_ = %d",myrank_);
#endif

  PostResult result = PostResult(field_);

  // timesteps when the solution is written
  const vector<double> soltime = result.get_result_times(field_->name());
  unsigned int numsoltimes = soltime.size();

  ///////////////////////////////////
  //  write geometry file          //
  ///////////////////////////////////
  const string geofilename = filename_ + "_"+ field_->name() + ".geo";
  const size_t found_path = geofilename.find_last_of("/\\");
  const string geofilename_nopath = geofilename.substr(found_path+1);
  WriteGeoFile(geofilename);
  vector<int> filesteps;
  filesteps.push_back(1);
  filesetmap_["geo"] = filesteps;
  vector<double> timesteps;
  timesteps.push_back(soltime[0]);
  timesetmap_["geo"] = timesteps;
  // at the moment, we can only print out the first step -> to be changed


  ///////////////////////////////////
  //  write solution fields files  //
  ///////////////////////////////////
  WriteAllResults(field_);

  // prepare the time sets and file sets for case file creation
  int setcounter = 0;
  int allresulttimeset = 0;
  for (map<string,vector<double> >::const_iterator entry = timesetmap_.begin(); entry != timesetmap_.end(); ++entry)
  {
      string key = entry->first;
      if ((entry->second).size()== numsoltimes)
      {
        if (allresulttimeset == 0)
        {
            setcounter++;
            allresulttimeset = setcounter;
        }
        timesetnumbermap_[key] = allresulttimeset; // reuse the default result time set, when possible
      }
      else
      {
        setcounter++;
        timesetnumbermap_[key] = setcounter; // a new time set number is needed
      }
  }

  setcounter = 0;
  for (map<string,vector<int> >::const_iterator entry = filesetmap_.begin(); entry != filesetmap_.end(); ++entry) {
      setcounter++;
      string key = entry->first;
      filesetnumbermap_[key] = setcounter;
  }


  ///////////////////////////////////
  //  now write the case file      //
  ///////////////////////////////////
  if (myrank_ == 0)
  {
    const string casefilename = filename_ + "_"+ field_->name() + ".case";
    ofstream casefile;
    casefile.open(casefilename.c_str());
    casefile << "# created using post_drt_ensight\n"<< "FORMAT\n\n"<< "type:\tensight gold\n";

    casefile << "\nGEOMETRY\n\n";
    casefile << "model:\t"<<timesetnumbermap_["geo"]<<"\t"<<filesetnumbermap_["geo"]<<"\t"<< geofilename_nopath<< "\n";

    casefile << "\nVARIABLE\n\n";
    casefile << GetVariableSection(filesetmap_, variablenumdfmap_, variablefilenamemap_);

    casefile << "\nTIME\n\n";
    casefile << GetTimeSectionStringFromTimesets(timesetmap_);

    casefile << "\nFILE\n\n";
    casefile << GetFileSectionStringFromFilesets(filesetmap_);

    casefile.close();
  }

        return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void XFluidEnsightWriter::WriteGeoFile(
  const string& geofilename)
{
  // open file
  ofstream geofile;
  if (myrank_ == 0)
  {
    geofile.open(geofilename.c_str());
    if (!geofile)
      dserror("failed to open file: %s", geofilename.c_str());
  }

  // header
  Write(geofile, "C Binary");

  // print out one timestep
  // if more are needed, this has to go into a loop
  map<string, vector<ofstream::pos_type> > resultfilepos;


  vector<ofstream::pos_type> fileposition;
  {
    WriteGeoFileOneTimeStep(geofile, resultfilepos, "geo");
  }

  // append index table
  // TODO: ens_checker complains if this is turned!!!! but I can't see, whats wrong here a.ger 11/07
  // it is also correct to ommit WriteIndexTable, however the EnsightGold Format manual says,
  // it would improve performance to have it on...
  // WriteIndexTable(geofile, fileposition);

  if (geofile.is_open())
    geofile.close();

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void XFluidEnsightWriter::WriteGeoFileOneTimeStep(
  ofstream& file,
  map<string, vector<ofstream::pos_type> >& resultfilepos,
  const string name) const
{
  // initial Intersection
  RCP<XFEM::InterfaceHandle> ih = rcp(new XFEM::InterfaceHandle(field_->discretization(),
                                                                cutterfield_->discretization()));
  // apply enrichments
  RCP<XFEM::DofManager> initialdofmanager = rcp(new XFEM::DofManager(ih));

  vector<ofstream::pos_type>& filepos = resultfilepos[name];
  Write(file, "BEGIN TIME STEP");
  filepos.push_back(file.tellp());

  Write(file, field_->name() + " geometry");
  Write(file, "Comment");

  //nodeidgiven_ is set to true inside the class constructor
  if (false)
    Write(file,"node id given");
  else
    Write(file, "node id assign");

  Write(file, "element id off");

  WriteGeoFilePart(file, resultfilepos, name, ih);

  Write(file, "END TIME STEP");
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void XFluidEnsightWriter::WriteGeoFilePart(
  ofstream& file,
  map<string, vector<ofstream::pos_type> >& resultfilepos,
  const string name,
  const RefCountPtr<XFEM::InterfaceHandle> ih ///< interfacehandle
  ) const
{
  // part + partnumber + comment
  // careful! field_->field_pos() returns the position of the ccarat
  // field, ignoring the discretizations. So if there are many
  // discretizations in one field, we have to do something different...
  Write(file, "part");
  Write(file, field_->field_pos()+1);
  Write(file, field_->name() + " field");

  Write(file, "coordinates");
  Write(file, NumNodesPerField(ih));

  cout << "writing " << NumNodesPerField(ih) << " nodes" << endl;

  // write the grid information
  WriteCoordinates(file, field_->discretization(), ih);
  WriteCells(file, field_->discretization(), ih);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void XFluidEnsightWriter::WriteCoordinates(
  ofstream& geofile,                          ///< filestream for the geometry
  const RefCountPtr<DRT::Discretization> dis, ///< discretization where the nodal positions are take from
  const RefCountPtr<XFEM::InterfaceHandle> ih ///< interfacehandle
  ) const
{
  const Epetra_Map* elementmap = dis->ElementRowMap();
  dsassert(elementmap->NumMyElements() == elementmap->NumGlobalElements(),
           "xfem filter cannot be run in parallel");
  const int nsd = 3; ///< number of space dimensions

                  // get the number of elements for each distype
const NumElePerDisType numElePerDisType = GetNumElePerDisType(dis, ih);

for (int isd = 0; isd < nsd; ++isd)
{
  // for each found distype write block of the same typed elements
  for (NumElePerDisType::const_iterator iter=numElePerDisType.begin(); iter != numElePerDisType.end(); ++iter)
  {
    const DRT::Element::DiscretizationType distypeiter = iter->first;

    for (int iele=0; iele<elementmap->NumMyElements(); ++iele)
    {
      const int elegid = elementmap->GID(iele);
      const DRT::Element* const actele = dis->gElement(elegid);
      const XFEM::DomainIntCells& domainintcells = ih->GetDomainIntCells(elegid, actele->Shape());
      for (XFEM::DomainIntCells::const_iterator cell = domainintcells.begin(); cell != domainintcells.end(); ++cell)
      {
        if (cell->Shape() == distypeiter)
        {
          const blitz::Array<double,2> xarray = cell->NodalPosXYZ(*actele);
          int numnode = cell->NumNode();
          if (distypeiter == DRT::Element::hex27)
          {
            numnode = 20;
          }
          for (int inen = 0; inen < numnode; ++inen)
          {
            Write(geofile, static_cast<float>(xarray(isd,inen)));
          }
        }
      }
    }
  }
}
}




/*----------------------------------------------------------------------*
  | write node connectivity for every element                  gjb 12/07 |
  *----------------------------------------------------------------------*/
void XFluidEnsightWriter::WriteCells(
  ofstream& geofile,
  const RefCountPtr<DRT::Discretization> dis,
  const RefCountPtr<XFEM::InterfaceHandle> ih ///< interfacehandle
  ) const
{
  const Epetra_Map* elementmap = dis->ElementRowMap();

  // get the number of elements for each distype
  const NumElePerDisType numElePerDisType_ = GetNumElePerDisType(dis, ih);

  // for each found distype write block of the same typed elements
  int counter = 0;
  NumElePerDisType::const_iterator iter;
  for (iter=numElePerDisType_.begin(); iter != numElePerDisType_.end(); ++iter)
  {
    const DRT::Element::DiscretizationType distypeiter = iter->first;
    const int ne = iter->second;
    const string ensightCellType = GetEnsightString(distypeiter);

    if (myrank_ == 0)
    {
      cout << "writing "<< iter->second<< " "<< DRT::DistypeToString(distypeiter) << " element(s) as "
           << ne << " " << ensightCellType << " ensight cell(s)..." << endl;
      Write(geofile, ensightCellType);
      Write(geofile, ne);
    }


    // loop all available elements
    for (int iele=0; iele<elementmap->NumMyElements(); ++iele)
    {
      const int elegid = elementmap->GID(iele);
      const XFEM::DomainIntCells& domainintcells = ih->GetDomainIntCells(elegid, dis->gElement(elegid)->Shape());
      for (XFEM::DomainIntCells::const_iterator cell = domainintcells.begin(); cell != domainintcells.end(); ++cell)
      {
        if (cell->Shape() == distypeiter)
        {
          switch (cell->Shape())
          {
          case DRT::Element::hex20:
          case DRT::Element::hex8:
          case DRT::Element::quad4:
          case DRT::Element::quad8:
          case DRT::Element::tet4:
          case DRT::Element::tet10:
          case DRT::Element::tri3:
          case DRT::Element::wedge6:
          case DRT::Element::wedge15:
          case DRT::Element::pyramid5:
          {
            // standard case with direct support
            const int numnp = cell->NumNode();
            for (int inode=0; inode<numnp; ++inode)
            {
              Write(geofile, counter+1);
              counter++;
            }
            break;
          }
          case DRT::Element::hex27:
          {
            // standard case with direct support
            const int numnp = 20;
            for (int inode=0; inode<numnp; ++inode)
            {
              Write(geofile, counter+1);
              counter++;
            }
            break;
          }
          default:
            dserror("don't know, how to write this element type as a Cell");
          };
        };
      };
    };
  };
}


/*!
 * \brief parse all elements and get the global(!) number of elements for each distype
 * \author gjb
 * \date 01/08
 */
NumElePerDisType XFluidEnsightWriter::GetNumElePerDisType(
  const RefCountPtr<DRT::Discretization> dis,
  const RefCountPtr<XFEM::InterfaceHandle> ih ///< interfacehandle
  ) const
{
  const Epetra_Map* elementmap = dis->ElementRowMap();

  NumElePerDisType numElePerDisType;
  for (int iele=0; iele<elementmap->NumMyElements(); ++iele)
  {
    const int elegid = elementmap->GID(iele);
    const XFEM::DomainIntCells& domainintcells = ih->GetDomainIntCells(elegid, dis->gElement(elegid)->Shape());
    for (XFEM::DomainIntCells::const_iterator cell = domainintcells.begin(); cell != domainintcells.end(); ++cell)
    {
      // update counter for current distype
      numElePerDisType[cell->Shape()]++;
    }
  }
  return numElePerDisType;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void XFluidEnsightWriter::WriteResult(
  const string groupname,
  const string name,
  const set<XFEM::PHYSICS::Field> fieldset
  )
{
  // Intersection
  RCP<XFEM::InterfaceHandle> ih = rcp(new XFEM::InterfaceHandle(field_->discretization(),
                                                                cutterfield_->discretization()));
  // apply enrichments
  RCP<XFEM::DofManager> dofman = rcp(new XFEM::DofManager(ih));

  // tell elements about the dofs and the integration
  {
    ParameterList eleparams;
    eleparams.set("action","store_xfem_info");
    eleparams.set("dofmanager",dofman);
    eleparams.set("assemble matrix 1",false);
    eleparams.set("assemble matrix 2",false);
    eleparams.set("assemble vector 1",false);
    eleparams.set("assemble vector 2",false);
    eleparams.set("assemble vector 3",false);
    field_->discretization()->Evaluate(eleparams,null,null,null,null,null);
  }

  // ensure that degrees of freedom in the discretization have been set
  field_->discretization()->FillComplete();

  PostResult result = PostResult(field_);
  result.next_result();
  if (!map_has_map(result.group(), groupname.c_str()))
    return;

  // new for file continuation
  bool multiple_files = false;
  const int numdf = fieldset.size();

  // open file
  const string filename = filename_ + "_"+ field_->name() + "."+ name;
  ofstream file;
  int startfilepos = 0;
  if (myrank_==0)
  {
    file.open(filename.c_str());
    startfilepos = file.tellp(); // file position should be zero, but we stay flexible
  }

  map<string, vector<ofstream::pos_type> > resultfilepos;
  //int stepsize = 0;


  cout<<"writing node-based field "<<name<<endl;
  // store information for later case file creation
  variableresulttypemap_[name] = "node";

  WriteNodalResultStep(file, result, resultfilepos, groupname, name, fieldset, ih, dofman);
  // how many bits are necessary per time step (we assume a fixed size)?
  const int stepsize = ((int) file.tellp())-startfilepos;
  if (stepsize <= 0) dserror("found invalid step size for result file");

  while (result.next_result())
  {
    const int indexsize = 80+2*sizeof(int)+(file.tellp()/stepsize+2)*sizeof(long);
    if (static_cast<long unsigned int>(file.tellp())+stepsize+indexsize>= FILE_SIZE_LIMIT_)
    {
      dserror("nope");
      FileSwitcher(file, multiple_files, filesetmap_, resultfilepos, stepsize, name, filename);
    }
    WriteNodalResultStep(file, result, resultfilepos, groupname, name, fieldset, ih, dofman);
  }

  // store information for later case file creation
  filesetmap_[name].push_back(file.tellp()/stepsize);// has to be done BEFORE writing the index table
  variablenumdfmap_[name] = numdf;
  variablefilenamemap_[name] = filename;
  // store solution times vector for later case file creation
  {
    PostResult res = PostResult(field_); // this is needed!
    vector<double> restimes = res.get_result_times(field_->name(),groupname);
    timesetmap_[name] = restimes;
  }

  // append index table
  WriteIndexTable(file, resultfilepos[name]);
  resultfilepos[name].clear();

  // close result file
  if (file.is_open())
    file.close();

  return;
}


vector<double> computeScalarCellNodeValues(
  const DRT::Element&  ele,
  const RCP<XFEM::InterfaceHandle>&  ih,
  const XFEM::ElementDofManager& dofman,
  const XFEM::DomainIntCell& cell,
  const XFEM::PHYSICS::Field field,
  const blitz::Array<double,1> elementvalues
  )
{
  const int nen_cell = DRT::UTILS::getNumberOfElementNodes(cell.Shape());
  const int numparam  = dofman.NumDofPerField(field);

  const blitz::Array<double,2> nodalPosXiDomain(cell.NodalPosXiDomainBlitz());

  // return value
  vector<double> cellvalues(nen_cell);

  //const blitz::Range _  = blitz::Range::all();

  // if cell node is on the interface, the value is not defined for a jump.
  // however, we approach the interface from one particular side and therefore,
  // -> we use the center of the cell to determine, where we come from
  const blitz::Array<double,1> cellcenterpos(cell.GetPhysicalCenterPosition(ele));

  // cell corner nodes
  //const blitz::Array<double,2> cellnodeposvectors = cell.NodalPosXYZ(ele);
  blitz::Array<double,1> enr_funct(numparam);
  for (int inen = 0; inen < nen_cell; ++inen)
  {
    //const blitz::Array<double,1> cellnodepos = cellnodeposvectors(_,inen);

    // shape functions
    const blitz::Array<double,1> funct(DRT::UTILS::shape_function_3D(
      nodalPosXiDomain(0,inen),
      nodalPosXiDomain(1,inen),
      nodalPosXiDomain(2,inen),
      ele.Shape()));

    XFEM::ComputeEnrichedShapefunction(ele, ih, dofman, field, cellcenterpos, XFEM::Enrichment::approachUnknown, funct, enr_funct);
    // interpolate value
    cellvalues[inen] = blitz::sum(elementvalues * enr_funct);
  }
  return cellvalues;
}




/*!
  \brief Write nodal values for one timestep

  No node has to have the same number of dofs.
*/
void XFluidEnsightWriter::WriteNodalResultStep(
  ofstream& file,
  PostResult& result,
  map<string, vector<ofstream::pos_type> >& resultfilepos,
  const string groupname,
  const string name,
  const set<XFEM::PHYSICS::Field> fieldset,
  const RCP<XFEM::InterfaceHandle> ih,
  const RCP<XFEM::DofManager> dofman
  ) const
{
  //-------------------------------------------
  // write some key words and read result data
  //-------------------------------------------

  vector<ofstream::pos_type>& filepos = resultfilepos[name];
  Write(file, "BEGIN TIME STEP");
  filepos.push_back(file.tellp());
  Write(file, "description");
  Write(file, "part");
  Write(file, field_->field_pos()+1);
  Write(file, "coordinates");

  const RefCountPtr<DRT::Discretization> dis = field_->discretization();
  //const Epetra_Map* nodemap = dis->NodeRowMap();
  const RefCountPtr<Epetra_Vector> data = result.read_result(groupname);

  const Epetra_Map* elementmap = dis->ElementRowMap();
  dsassert(elementmap->NumMyElements() == elementmap->NumGlobalElements(),
           "xfem filter cannot be run in parallel");

  const int numdf = fieldset.size();

  bool propernumdf;
  if (numdf == 1 or numdf == 3 or numdf == 6) {
    propernumdf = true;
  } else {
    propernumdf = false;
  }
  if (not propernumdf)
  {
    dserror("number of output dofs is not 1, 3 or 6");
  }



  // get the number of elements for each distype
  const NumElePerDisType numElePerDisType = GetNumElePerDisType(dis, ih);

  //
  const XFEM::Enrichment enr_std(0, XFEM::Enrichment::typeStandard);

  for (set<XFEM::PHYSICS::Field>::const_iterator fielditer=fieldset.begin(); fielditer!=fieldset.end(); ++fielditer)
  {
    const XFEM::PHYSICS::Field field = *fielditer;

    // for each found distype, write block of the same typed elements
    int counter = 0;
    for (NumElePerDisType::const_iterator iter=numElePerDisType.begin(); iter != numElePerDisType.end(); ++iter)
    {
      const DRT::Element::DiscretizationType distypeiter = iter->first;
      for (int iele=0; iele<elementmap->NumMyElements(); ++iele)
      {
        const int elegid = elementmap->GID(iele);
        DRT::Element* const actele = dis->gElement(elegid);

        const DRT::Element::DiscretizationType stressdistype = XFLUID::getStressInterpolationType3D(actele->Shape());
        const int numeleparam = DRT::UTILS::getNumberOfElementNodes(stressdistype);

        // create local copy of information about dofs
        const XFEM::ElementDofManager eledofman = dofman->constructElementDofManager(*actele, numeleparam);

        vector<int> lm;
        vector<int> lmowner;
        actele->LocationVector(*(field_->discretization()),lm,lmowner);

        // extract local values from the global vector
        vector<double> myvelnp(lm.size());

        DRT::UTILS::ExtractMyValues(*data,myvelnp,lm);

        const int numparam = eledofman.NumDofPerField(field);
        const vector<int>& dofpos = eledofman.LocalDofPosPerField(field);
        //cout << XFEM::PHYSICS::physVarToString(field) << ": numparam = " << numparam << ": lm.size() = " << lm.size() << endl;

        blitz::Array<double,1> elementvalues(numparam);
        for (int iparam=0; iparam<numparam; ++iparam)   elementvalues(iparam) = myvelnp[dofpos[iparam]];

        const XFEM::DomainIntCells& domainintcells = ih->GetDomainIntCells(elegid, actele->Shape());
        for (XFEM::DomainIntCells::const_iterator cell = domainintcells.begin(); cell != domainintcells.end(); ++cell)
        {
          if (cell->Shape() == distypeiter)
          {
            const vector<double> cellvalues = computeScalarCellNodeValues(*actele, ih, eledofman, *cell, field, elementvalues);
            int numnode = cell->NumNode();
            if (distypeiter == DRT::Element::hex27)
            {
              numnode = 20;
            }
            for (int inode = 0; inode < numnode; ++inode)
            {
              Write(file, static_cast<float>(cellvalues[inode]));
              counter++;
            }
          }
        }
      }
    }
    cout << "number of entries per field " << counter << endl;
  }

  Write(file, "END TIME STEP");

  return;
}

int XFluidEnsightWriter::NumNodesPerField(
  const RefCountPtr<XFEM::InterfaceHandle> ih ///< interfacehandle
  ) const
{
  RCP<DRT::Discretization> dis = field_->discretization();

  const Epetra_Map* elementmap = dis->ElementRowMap();
  dsassert(elementmap->NumMyElements() == elementmap->NumGlobalElements(),
           "xfem filter cannot be run in parallel");

  // loop all available elements
  int counter = 0;
  for (int iele=0; iele<elementmap->NumMyElements(); ++iele)
  {
    const int elegid = elementmap->GID(iele);
    const XFEM::DomainIntCells& domainintcells = ih->GetDomainIntCells(elegid, dis->gElement(elegid)->Shape());
    for (XFEM::DomainIntCells::const_iterator cell = domainintcells.begin(); cell != domainintcells.end(); ++cell)
    {
      switch (cell->Shape())
      {
      case DRT::Element::hex20:
      case DRT::Element::hex8:
      case DRT::Element::quad4:
      case DRT::Element::quad8:
      case DRT::Element::tet4:
      case DRT::Element::tet10:
      case DRT::Element::tri3:
      case DRT::Element::wedge6:
      case DRT::Element::wedge15:
      case DRT::Element::pyramid5:
      {
        // standard case with direct support
        counter += cell->NumNode();
        break;
      }
      case DRT::Element::hex27:
      {
        counter += 20;
        break;
      }
      default:
        dserror("don't know, how to write this element type as a Cell");
      };
    };
  };
  return counter;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
string XFluidEnsightWriter::GetEnsightString(
  const DRT::Element::DiscretizationType distype) const
{
  map<DRT::Element::DiscretizationType, string>::const_iterator entry;
  switch (distype)
  {
  case DRT::Element::hex27:
    entry = distype2ensightstring_.find(DRT::Element::hex20);
    break;
  case DRT::Element::quad9:
    entry = distype2ensightstring_.find(DRT::Element::quad8);
    break;
  default:
    entry = distype2ensightstring_.find(distype);
  }
  if (entry == distype2ensightstring_.end())
    dserror("no entry in distype2ensightstring_ found");
  return entry->second;
}


#endif
