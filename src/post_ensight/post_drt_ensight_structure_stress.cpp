/*!
  \file post_drt_ensight_structure_stress.cpp

  \brief postprocessing of structural stresses

  <pre>
  Maintainer: Lena Wiechert
              wiechert@lnm.mw.tum.de
              http://www.lnm.mw.tum.de/Members/wiechert
              089-289-15303
  </pre>

*/

#ifdef CCADISCRET

#include "post_drt_ensight_writer.H"
#include <string>
#include "post_drt_ensight_single_field_writers.H"
#include "../drt_lib/drt_utils.H"

using namespace std;


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
    exit(1);
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
  // write some key words
  //--------------------------------------------------------------------

  vector<ofstream::pos_type>& filepos = resultfilepos[name];
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
    exit(1);
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


#endif
