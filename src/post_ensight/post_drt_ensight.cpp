/*!
\file post_drt_ensight.cpp

\brief Ensight filter

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/kuettler
            089 - 289-15238
</pre>

*/

#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE

#include <iostream>
#include <sstream>
#include <cstdio>
#include <vector>
#include <map>
#include <fstream>
#include <string>

#include "../post_drt_common/post_drt_common.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_node.H"
#include "../drt_lib/drt_element.H"

using namespace std;


//! defines how 4 quad4 elements are constructed from a quad9
const int subquadmap[4][4] = {{ 0, 4, 8, 7},
                              { 4, 1, 5, 8},
                              { 8, 5, 2, 6},
                              { 7, 8, 6, 3}};

//! defines how 8 hex8 elements are constructed from a hex27
//  ;-) its symetric for some reason
const int subhexmap[8][8] = {{  0,  8, 20, 11, 12, 21, 26, 24},
                             {  8,  1,  9, 20, 21, 13, 22, 26},
                             { 20,  9,  2, 10, 26, 22, 14, 23},
                             { 11, 20, 10,  3, 24, 26, 23, 15},
                             { 12, 21, 26, 24,  4, 16, 25, 19},
                             { 21, 13, 22, 26, 16,  5, 17, 25},
                             { 26, 22, 14, 23, 25, 17,  6, 18},
                             { 24, 26, 23, 15, 19, 25, 18,  7}};

typedef map<DRT::Element::DiscretizationType, int> NumElePerDisType;

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
class EnsightWriter
{
public:

  EnsightWriter(PostField* field, const string filename);
  virtual ~EnsightWriter() {}

  //! write the whole thing
  void WriteFiles();

protected:

  //! look for problem dependent result entries and write them
  virtual void WriteResults(PostField* field) = 0;

  /*!
    \brief write all time steps of a result

    Write nodal results. The results are taken from a reconstructed
    Epetra_Vector. In many cases this vector will contain just one
    variable (displacements) and thus is easy to write as a whole. At
    other times, however, there is more than one result (velocity,
    pressure) and we want to write just one part of it. So we have to
    specify which part.

    \author u.kue
    \date 03/07
  */
  void WriteResultsAllSteps(
    const string groupname,  ///< name of the result group in the control file
    const string name,       ///< name of the result to be written
    const int numdf,         ///< number of dofs per node to this result
    const int from=0         ///< start position of values in nodes
    );

private:

  template <class T> void Write(ofstream& os,T i) { os.write(reinterpret_cast<const char*>(&i),sizeof(T)); }
  void Write(ofstream& os,const string s) { WriteString(os,s); }
  void Write(ofstream& os,const char* s) { WriteString(os,s); }

  void WriteString(
    ofstream& stream,
    const string str);
    
  void WriteCoordinates(
    ofstream&                              geofile,
    const RefCountPtr<DRT::Discretization> dis);
  void WriteCells(ofstream& geofile);
  void WriteResultStep(
    ofstream& file,
    PostResult& result,
    const string groupname,
    const string name,
    const int numdf,
    const int from);
  void WriteIndexTable(
    ofstream& file,
    const vector<ofstream::pos_type>& filepos);
  /*!
   * \brief create string for the VARIABLE section 
   *        that corresponds to the current field
   */
  string GetVariableEntryForCaseFile(
    const int numdf,             ///< degrees of freedom per node for this field
    const unsigned int fileset,
    const string name,
    const string filename);
    
  NumElePerDisType GetNumElePerDisType(
    const RefCountPtr<DRT::Discretization> dis
  ) const;
    
  string GetEnsightString(
    const DRT::Element::DiscretizationType distype
  );
    
  /*!
   * \brief if files become to big, this routine switches to a new file
   */
  void FileSwitcher(
    ofstream& file,
    unsigned int& fileset,
    const int stepsize,
    const string name,
    const string filename);
  
  int GetNumEleOutput(
    const DRT::Element::DiscretizationType distype,
    const int                              numele
  ) const;

  PostField* field_;
  string filename_;
  ofstream casefile_;
  
  map<string, vector<ofstream::pos_type> > resultfilepos_;
  vector<vector<int> > filesets_;
  
  //! maps a distype to the corresponding Ensight cell type
  map<DRT::Element::DiscretizationType, string> distype2ensightstring_;
};


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
EnsightWriter::EnsightWriter(PostField* field, const string filename)
  : field_(field), filename_(filename)
{
  // map between distype in BACI and Ensight string
  distype2ensightstring_.clear();
  distype2ensightstring_[DRT::Element::hex8]  = "hexa8";
  distype2ensightstring_[DRT::Element::hex20] = "hexa20";
  distype2ensightstring_[DRT::Element::hex27] = "hexa20";
  distype2ensightstring_[DRT::Element::tet4]  = "tetra4";
  distype2ensightstring_[DRT::Element::tet10] = "tetra10";
  distype2ensightstring_[DRT::Element::quad4] = "quad4";
  distype2ensightstring_[DRT::Element::quad8] = "quad8";
  distype2ensightstring_[DRT::Element::quad9] = "quad9";
  distype2ensightstring_[DRT::Element::tri3]  = "tria3";
  distype2ensightstring_[DRT::Element::tri6]  = "tria6";
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void EnsightWriter::WriteFiles()
{
  //
  // Print out the .geo file
  //
  
  // open file
  string geofilename = filename_ + "_" + field_->name() + ".geo";
  ofstream geofile;
  geofile.open(geofilename.c_str());
  if (!geofile)  dserror("failed to open file: %s", geofilename.c_str());
  
  // header
  Write(geofile,"C Binary");
  
  // print out one timestep
  // if more are needed, this has to go into a loop
  vector<ofstream::pos_type> fileposition_;
  {
    Write(geofile,"BEGIN TIME STEP");
    fileposition_.push_back(geofile.tellp());
    Write(geofile,field_->name() + " geometry");
    Write(geofile,"Comment");
    //Write(geofile_,"node id given");
    Write(geofile,"node id assign");
    Write(geofile,"element id off");
  
    // part + partnumber + comment
    // careful! field_->field_pos() returns the position of the ccarat
    // field, ignoring the discretizations. So if there are many
    // discretizations in one field, we have to do something different...
    Write(geofile,"part");
    Write(geofile,field_->field_pos()+1);
    Write(geofile,field_->name() + " field");
  
    Write(geofile,"coordinates");
    Write(geofile,field_->num_nodes());
  
    // write the grid information
    WriteCoordinates(geofile, field_->discretization());
    WriteCells(geofile);
  
    Write(geofile,"END TIME STEP");
  }
  
  // append index table
  WriteIndexTable(geofile,fileposition_);
  
  geofile.close();

  // loop all results to get number of result timesteps
  PostResult result = PostResult(field_);
  vector<double> soltime_;  // timesteps when the solution is written
  if (result.next_result())
    soltime_.push_back(result.time());
  else
    dserror("kann das ueberhaupt sein?");

  while (result.next_result())  soltime_.push_back(result.time());

  //
  // now do the case file
  //
  const string casefilename = filename_ + "_" + field_->name() + ".case";
  casefile_.open(casefilename.c_str());
  casefile_ << "# created using post_drt_ensight_xfem\n"
            << "FORMAT\n\n"
            << "type:\tensight gold\n\n"
            << "GEOMETRY\n\n"
            << "model:\t2\t2\t" << geofilename << "\n\n"
            << "VARIABLE\n\n";

  // whatever result we need, inside the variable entries in the case file
  //  are generated so do not change the order of the calls here
  WriteResults(field_);
  
  // TIME section
  //
  //
  // write time steps for result file (time set 1)
  casefile_ << "\nTIME\n"
            << "time set:\t\t1\n"
            << "number of steps:\t" << soltime_.size() << "\n"
            << "time values: ";
  for (unsigned i=0; i<soltime_.size(); ++i)
  {
    casefile_ << soltime_[i] << " ";
    if (i%8==0 && i!=0)
      casefile_ << "\n";
  }

  // write time steps for geometry file (time set 2)
  // at the moment, we can only print out the first step -> to be changed
  vector<double> geotime_;  // timesteps when the geometry is written
  geotime_.push_back(soltime_[0]);
  casefile_ << "\n\ntime set:\t\t2\n"
            << "number of steps:\t" << geotime_.size() << "\n"
            << "time values: ";
  for (unsigned i=0; i<geotime_.size(); ++i)
  {
    casefile_ << geotime_[i] << " ";
    if (i%8==0 && i!=0)
      casefile_ << "\n";
  }
  
  //
  // FILE section
  //
  casefile_ << "\n\n"
            << "FILE\n"
            << "file set:\t\t1\n"
            << "number of steps:\t" << soltime_.size() << "\n\n"
            << "file set:\t\t2\n"
            << "number of steps:\t" << geotime_.size() << "\n\n"
    ;
  for (unsigned int i = 0; i < filesets_.size(); ++i)
  {
    casefile_ << "\nfile set:\t\t" << 3+i << "\n";
    for (unsigned int j = 0; j < filesets_[i].size(); ++j)
    {
      casefile_ << "filename index:\t" << 1+j << "\n"
                << "number of steps:\t" << filesets_[i][j] << "\n";
    }
  }

  casefile_.close();

}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void EnsightWriter::WriteCoordinates(
ofstream&                              geofile,
const RefCountPtr<DRT::Discretization> dis
)
{
  const Epetra_Map* nodemap = dis->NodeRowMap();
  dsassert(nodemap->NumMyElements() == nodemap->NumGlobalElements(), "filter cannot be run in parallel");

  const int numnp = nodemap->NumMyElements();
  // write node ids
//  for (int inode=0; inode<numnp; ++inode)
//  {
//    Write(geofile,nodemap->GID(inode)+1);
//  }
  
  const int NSD = 3; ///< number of space dimensions
  // ensight format requires x_1 .. x_n, y_1 .. y_n, z_1 ... z_n
  for (int isd=0; isd<NSD; ++isd)
  {
    for (int inode=0; inode<numnp; inode++)
    {
      const int gid = nodemap->GID(inode);
      DRT::Node* actnode = dis->gNode(gid);
      Write(geofile,static_cast<float>(actnode->X()[isd]));
    }
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void EnsightWriter::WriteCells(ofstream& geofile)
{
  const RefCountPtr<DRT::Discretization> dis = field_->discretization();
  const Epetra_Map* nodemap = dis->NodeRowMap();
  dsassert(nodemap->NumMyElements() == nodemap->NumGlobalElements(), "filter cannot be run in parallel");

  const Epetra_Map* elementmap = dis->ElementRowMap();
  dsassert(elementmap->NumMyElements() == elementmap->NumGlobalElements(), "filter cannot be run in parallel");

  // get the number of elements for each distype
  const NumElePerDisType numElePerDisType = GetNumElePerDisType(dis);
  
  // for each found distype write block of the same typed elements
  NumElePerDisType::const_iterator iter;
  for (iter=numElePerDisType.begin(); iter != numElePerDisType.end(); ++iter)
  {
    const DRT::Element::DiscretizationType distypeiter = iter->first;
    const int ne = GetNumEleOutput(distypeiter, iter->second);
    const string ensightString = GetEnsightString(distypeiter);
    
    cout << "writing "<< iter->second << " " 
         << distype2ensightstring_[distypeiter] << " elements"
         << " (" << ne << " cells) per distype."
         << " ensight output celltype: " << ensightString
         << endl;
    
    Write(geofile,ensightString);                                         
    Write(geofile,ne);

    // loop all available elements
    for (int iele=0; iele<elementmap->NumMyElements(); ++iele)
    {
      DRT::Element* actele = dis->gElement(elementmap->GID(iele));
      DRT::Node** nodes = actele->Nodes();
      if (actele->Shape() == distypeiter)
      {
        switch (actele->Shape())
        {
        // standard case with direct support
        case DRT::Element::hex20:
        case DRT::Element::hex8:
        case DRT::Element::quad4:
        case DRT::Element::quad8:
        case DRT::Element::tet4:
        case DRT::Element::tet10:
        case DRT::Element::tri3:
        case DRT::Element::wedge6:
        case DRT::Element::wedge15:
        {
          const int numnp = actele->NumNode();
          for (int inode=0; inode<numnp; ++inode)
            Write(geofile,nodemap->LID(nodes[inode]->Id())+1);
          break;
        }
        // special cases
        case DRT::Element::hex27:
        {
          for (int isubele=0; isubele<8; ++isubele)
            for (int isubnode=0; isubnode<8; ++isubnode)
              Write(geofile,nodemap->LID(nodes[subhexmap[isubele][isubnode]]->Id())+1);
          break;
        }
        case DRT::Element::quad9:
        {
          // write subelements
          for (int isubele=0; isubele<4; ++isubele)
            for (int isubnode=0; isubnode<4; ++isubnode)
              Write(geofile,nodemap->LID(nodes[subquadmap[isubele][isubnode]]->Id())+1);
          break;
        }
        default:
          dserror("don't know, how to write this element type as a Cell");
        };
      };
    };
  };
}


/*!
 * \brief parse all elements and get the number of elements for each distype
 */
NumElePerDisType EnsightWriter::GetNumElePerDisType(
  const RefCountPtr<DRT::Discretization> dis
  ) const
{
  const Epetra_Map* elementmap = dis->ElementRowMap();
  
  NumElePerDisType numElePerDisType;
  for (int iele=0; iele<elementmap->NumMyElements(); ++iele)
  {
    DRT::Element* actele = dis->gElement(elementmap->GID(iele));
    const DRT::Element::DiscretizationType distype = actele->Shape();
    // update counter, if element has been identified
    numElePerDisType[distype]++;
  }
  return numElePerDisType;
}

int EnsightWriter::GetNumEleOutput(
  const DRT::Element::DiscretizationType distype,
  const int                              numele
  ) const
{
  
  int numeleout = 0;
  switch (distype)
  {
  case DRT::Element::hex27:
    numeleout = 8*numele;
    break;
  case DRT::Element::quad9:
    numeleout = 4*numele;
    break;
  default:
    numeleout = numele;
  }
  return numeleout;
}

string EnsightWriter::GetEnsightString(
  const DRT::Element::DiscretizationType distype
  )
{
  string str;
  switch (distype)
  {
  case DRT::Element::hex27:
    str = distype2ensightstring_[DRT::Element::hex8];
    break;
  case DRT::Element::quad9:
    str = distype2ensightstring_[DRT::Element::quad4];
    break;
  default:
    str = distype2ensightstring_[distype];
  }
  return str;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void EnsightWriter::WriteResultsAllSteps(
  const string groupname,
  const string name,
  const int numdf,
  const int from)
{
  // maximum file size
  const unsigned FILE_SIZE_LIMIT = 0x7fffffff; // 2GB
  
  PostResult result = PostResult(field_);
  result.next_result();
  if (!map_has_map(result.group(), const_cast<char*>(groupname.c_str())))
    return;

  // new for file continuation
  unsigned int fileset = 1;
  const Epetra_Map* nodemap = field_->discretization()->NodeRowMap();
  const int numnp = nodemap->NumMyElements();
  const int stepsize = 5*80+sizeof(int)+numdf*numnp*sizeof(float);

  string filename = filename_ + "_" + field_->name() + "." + name;
  ofstream file(filename.c_str());

  WriteResultStep(file, result, groupname, name, numdf, from);
  while (result.next_result())
  {
    const int indexsize = 80+2*sizeof(int)+(file.tellp()/stepsize+2)*sizeof(long);
    if (static_cast<long unsigned int>(file.tellp())+stepsize+indexsize >= FILE_SIZE_LIMIT)
    {
      FileSwitcher(file, fileset, stepsize, name, filename);
    }
    WriteResultStep(file, result, groupname, name, numdf, from);
  }

  // append index table
  WriteIndexTable(file,resultfilepos_[name]);
  resultfilepos_[name].clear();

  if (fileset != 1)
  {
    filesets_[fileset-3].push_back(file.tellp()/stepsize);
    filename += "***";
  }
  
  casefile_ << GetVariableEntryForCaseFile(numdf, fileset, name, filename);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void EnsightWriter::FileSwitcher(
  ofstream& file,
  unsigned int& fileset,
  const int stepsize,
  const string name,
  const string filename
)
{
  if (fileset == 1) {
    fileset = filesets_.size()+3;
    vector<int> numsteps;
    numsteps.push_back(file.tellp()/stepsize);
    filesets_.push_back(numsteps);
    // append index table
    WriteIndexTable(file,resultfilepos_[name]);
    resultfilepos_[name].clear();
    file.close();
    rename(filename.c_str(), (filename+"001").c_str());
    file.open((filename+"002").c_str());
  }
  else 
  {
    filesets_[fileset-3].push_back(file.tellp()/stepsize);
    ostringstream newfilename;
    newfilename << filename;
    newfilename.width(3);
    newfilename.fill('0');
    newfilename << filesets_[fileset-3].size()+1;
    // append index table
    WriteIndexTable(file,resultfilepos_[name]);
    resultfilepos_[name].clear();
    file.close();
    file.open(newfilename.str().c_str());
  }
}


/*
 */
string EnsightWriter::GetVariableEntryForCaseFile(
  const int numdf,
  const unsigned int fileset,
  const string name,
  const string filename
  )
{
  stringstream str;
  // create variable entry
  switch (numdf)
  {
  case 6:
    str << "tensor symm per node:\t1\t" << fileset << "\t" << name << "\t" << filename << "\n";
    break;
  case 3: case 2:
    str << "vector per node:\t1\t"      << fileset << "\t" << name << "\t" << filename << "\n";
    break;
  case 1:
    str << "scalar per node:\t1\t"      << fileset << "\t" << name << "\t" << filename << "\n";
    break;
  default:
    dserror("unknown number of dof per node");
  };
  return str.str();
}

/*!
 \brief Write nodal values for one timestep
 
 Each node has to have the same number of dofs.
 */
void EnsightWriter::WriteResultStep(
  ofstream& file,
  PostResult& result,
  const string groupname, 
  const string name,
  const int numdf, 
  const int from)
{
  vector<ofstream::pos_type>& filepos = resultfilepos_[name];
  Write(file,"BEGIN TIME STEP");
  filepos.push_back(file.tellp());

  Write(file,"description");
  Write(file,"part");
  Write(file,field_->field_pos()+1);
  Write(file,"coordinates");

  const RefCountPtr<DRT::Discretization> dis = field_->discretization();
  const Epetra_Map* nodemap = dis->NodeRowMap();
  const RefCountPtr<Epetra_Vector> data = result.read_result(groupname);
  const Epetra_BlockMap& datamap = data->Map();

  const int numnp = nodemap->NumMyElements();
  for (int idf=0; idf<numdf; ++idf)
  {
    for (int inode=0; inode<numnp; inode++)
    {
      DRT::Node* n = dis->lRowNode(inode);
      const int lid = datamap.LID(dis->Dof(n,from+idf));
      if (lid > -1)
      {
        Write(file,static_cast<float>((*data)[lid]));
      }
      else
      {
        // Assume we have to write a value here.
        Write<float>(file,0.);
      }
    }
  }

  // 2 component vectors in a 3d problem require a row of zeros.
  // do we really need this?
  if (numdf==2)
  {
    for (int inode=0; inode<numnp; inode++)
    {
      Write<float>(file,0.);
    }
  }

  Write(file,"END TIME STEP");
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void EnsightWriter::WriteIndexTable(
  ofstream&                         file,
  const vector<ofstream::pos_type>& filepos)
{
  ofstream::pos_type lastpos = file.tellp();
  const unsigned steps = filepos.size();
  Write(file,steps);
  for (unsigned i=0; i<steps; ++i)
  {
    Write<long>(file,filepos[i]);
  }
  Write(file,0);
  Write<long>(file,lastpos);
  Write(file,"FILE_INDEX");
}


/*----------------------------------------------------------------------*/
/*!
  \brief write strings of exactly 80 chars
 */
/*----------------------------------------------------------------------*/
void EnsightWriter::WriteString(ofstream& stream, const string str)
{
  // we need to write 80 bytes per string
  vector<char> s(str.begin(),str.end());
  while (s.size()<80)
    s.push_back('\0');
  stream.write(&s[0],80);
}



/*
  \brief Writer for structural problems
 */
class StructureEnsightWriter : public EnsightWriter
{
public:
  StructureEnsightWriter(PostField* field, string filename)
    : EnsightWriter(field,filename) {}

protected:

  virtual void WriteResults(PostField* field);
};


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void StructureEnsightWriter::WriteResults(PostField* field)
{
  EnsightWriter::WriteResultsAllSteps("displacement","displacement",field->problem()->num_dim());
  EnsightWriter::WriteResultsAllSteps("velocity","velocity",field->problem()->num_dim());
  EnsightWriter::WriteResultsAllSteps("acceleration","acceleration",field->problem()->num_dim());
}


/*----------------------------------------------------------------------*/
/*
  \brief Writer for fluid problems
 */
/*----------------------------------------------------------------------*/
class FluidEnsightWriter : public EnsightWriter
{
public:
  FluidEnsightWriter(PostField* field, string filename)
    : EnsightWriter(field,filename) {}

protected:

  virtual void WriteResults(PostField* field);
};


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FluidEnsightWriter::WriteResults(PostField* field)
{
  
  EnsightWriter::WriteResultsAllSteps("velnp","velocity",field->problem()->num_dim());
  EnsightWriter::WriteResultsAllSteps("velnp","pressure",1,field->problem()->num_dim());
  EnsightWriter::WriteResultsAllSteps("residual","residual",field->problem()->num_dim());
  EnsightWriter::WriteResultsAllSteps("dispnp","displacement",field->problem()->num_dim());
}


/*----------------------------------------------------------------------*/
/*
  \brief Writer for ale problems
 */
/*----------------------------------------------------------------------*/
class AleEnsightWriter : public EnsightWriter
{
public:
  AleEnsightWriter(PostField* field, string filename)
    : EnsightWriter(field,filename) {}

protected:

  virtual void WriteResults(PostField* field);
};


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void AleEnsightWriter::WriteResults(PostField* field)
{
  EnsightWriter::WriteResultsAllSteps("dispnp","displacement",field->problem()->num_dim());
}


/*!
  \brief filter main routine

  Write binary ensight format.

  The ens_checker that is part of ensight can be used to verify the
  files generated here.

  \author u.kue
  \date 03/07
 */
int main(int argc, char** argv)
{
  Teuchos::CommandLineProcessor My_CLP;
  My_CLP.setDocString(
    "Post DRT ensight Filter\n"
    );

  PostProblem problem = PostProblem(My_CLP,argc,argv);

#if 0
  for (int i = 0; i<problem.num_discr(); ++i)
  {
    PostField* field = problem.get_discretization(i);
    StructureEnsightWriter writer(field, problem.basename());
    writer.WriteFiles();
  }
#endif

  // each problem type is different and writes different results
  switch (problem.Problemtype())
  {
  case prb_fsi:
  {
    string basename = problem.basename();
//     PostField* structfield = problem.get_discretization(0);
//     StructureEnsightWriter structwriter(structfield, basename);
//     structwriter.WriteFiles();

    PostField* fluidfield = problem.get_discretization(1);
    FluidEnsightWriter fluidwriter(fluidfield, basename);
    fluidwriter.WriteFiles();
    break;
  }
  case prb_structure:
  {
    PostField* field = problem.get_discretization(0);
    StructureEnsightWriter writer(field, problem.basename());
    writer.WriteFiles();
    break;
  }
  case prb_fluid:
  {
    PostField* field = problem.get_discretization(0);
    FluidEnsightWriter writer(field, problem.basename());
    writer.WriteFiles();
    break;
  }
  case prb_ale:
  {
    PostField* field = problem.get_discretization(0);
    AleEnsightWriter writer(field, problem.basename());
    writer.WriteFiles();
    break;
  }
  default:
    dserror("problem type %d not yet supported", problem.Problemtype());
  }

  return 0;
}


#endif
#endif
