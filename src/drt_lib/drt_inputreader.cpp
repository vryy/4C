/*----------------------------------------------------------------------*/
/*!
\file drt_inputreader.cpp

\brief Internal classes to read elements and nodes

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>
*/
/*----------------------------------------------------------------------*/

#ifdef CCADISCRET


#include "drt_inputreader.H"
#include "drt_utils.H"
#include "linalg_utils.H"
#include "standardtypes_cpp.H"

#include <Epetra_Time.h>
#include <iterator>

#ifdef PARMETIS
typedef int idxtype;
extern "C"
{
  void ParMETIS_V3_PartKway(idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, 
                            idxtype *adjwgt, int *wgtflag, int *numflag, int *ncon, int *nparts, 
                            float *tpwgts, float *ubvec, int *options, int *edgecut, idxtype *part, 
                            MPI_Comm *comm);
}
#endif

/*!----------------------------------------------------------------------
  \brief file pointers

  <pre>                                                         m.gee 8/00
  This structure struct _FILES allfiles is defined in input_control_global.c
  and the type is in standardtypes.h
  It holds all file pointers and some variables needed for the FRSYSTEM
  </pre>
 *----------------------------------------------------------------------*/
extern struct _FILES  allfiles;


/*----------------------------------------------------------------------*/
/*!
  \brief remove all leading and trailing whitespaces from a string.

  Note: consecutive whitespaces inside the string will be reduced to a
  single space.

  \author u.kue
  \date 03/07
*/
/*----------------------------------------------------------------------*/
static std::string trim(const std::string& line)
{
  std::istringstream t;
  std::string s;
  std::string newline;
  t.str(line);
  while (t >> s)
  {
    newline.append(s);
    newline.append(" ");
  }
  if (newline.size()>0)
    newline.resize(newline.size()-1);
  return newline;
}


namespace DRT
{
namespace INPUT
{


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
DatFileReader::DatFileReader(string filename, Teuchos::RCP<Epetra_Comm> comm, int outflag, bool dumpinput)
  : filename_(filename), comm_(comm), outflag_(outflag)
{
  ReadDat();
  if (dumpinput)
    DumpInput();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
string DatFileReader::MyInputfileName() const
{
  return filename_;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
int DatFileReader::MyOutputFlag() const
{
  return outflag_;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DatFileReader::Activate()
{
  // Publish (some) internal data to the old ccarat reading system
  // Note that these links remain intact even when the reader goes
  // away.

  allfiles.numrows = numrows_;
  allfiles.input_file_hook = &inputfile_[0];
  allfiles.input_file = &lines_[0];
  //allfiles.inputfile_name = const_cast<char*>(filename_.c_str());

  // set fr-system to begin of input_file
  allfiles.actrow = 0;
  allfiles.actplace = &(allfiles.input_file[0][0]);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool DatFileReader::ReadSection(string name, Teuchos::ParameterList& list)
{
  if (name.length() < 3 or name[0]!='-' or name[1]!='-')
    dserror("illegal section name '%s'", name.c_str());

  Teuchos::ParameterList& sublist = FindSublist(name.substr(2), list);

  if (positions_.find(name)==positions_.end())
    return false;

  for (unsigned pos = positions_[name]+1;
       pos < lines_.size();
       ++pos)
  {
    string line = lines_[pos];
    if (line[0]=='-' and line[1]=='-')
    {
      break;
    }

    // we expect a line: key = value
    // The first = in the line will be taken for the
    // separator. Thus we cannot have a = in a key.
    std::string::size_type delim = line.find('=');
    if (delim==std::string::npos)
      dserror("no key=value pair in line %d: %s", pos, line.c_str());

    std::string key   = line.substr(0,delim-1);
    std::string value = line.substr(delim+2);

    // Now parse the value. Find integers and doubles if there are
    // any.
    AddEntry(key, value, sublist);
  }
  return true;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool DatFileReader::ReadGidSection(string name, Teuchos::ParameterList& list)
{
  if (name.length() < 3 or name[0]!='-' or name[1]!='-')
    dserror("illegal section name '%s'", name.c_str());

  Teuchos::ParameterList& sublist = FindSublist(name.substr(2), list);

  if (positions_.find(name)==positions_.end())
    return false;

  for (unsigned pos = positions_[name]+1;
       pos < lines_.size();
       ++pos)
  {
    string line = lines_[pos];
    if (line[0]=='-' and line[1]=='-')
    {
      break;
    }

    string key;
    string value;

    string::size_type loc = line.find(" ");
    if (loc==string::npos)
    {
      //dserror("line '%s' with just one word in GiD parameter section", line.c_str());
      key = line;
    }
    else
    {
      //if (line.find(" ", loc+1)!=string::npos)
      //  dserror("more that two words on line '%s' in GiD parameter section", line.c_str());
      key = line.substr(0,loc);
      value = line.substr(loc+1);
    }

    // Now parse the value. Find integers and doubles if there are
    // any.
    AddEntry(key, value, sublist);
  }

  return true;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
vector<const char*> DatFileReader::Section(string name) const
{
  vector<const char*> sec;

  map<string,unsigned>::const_iterator i = positions_.find(name);
  if (i!=positions_.end())
  {

    for (unsigned pos = i->second+1;
         pos < lines_.size();
         ++pos)
    {
      const char* line = lines_[pos];
      if (line[0]=='-' and line[1]=='-')
      {
        break;
      }
      sec.push_back(line);
    }
  }

  return sec;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DatFileReader::ReadDesign(const std::string& name, std::vector<std::vector<int> >& dobj_fenode) const
{
  std::map<int,std::set<int> > topology;

  std::string sectionname = name + "-NODE TOPOLOGY";
  std::string marker = std::string("--") + sectionname;

  // this is still the old fr* stuff
  if (frfind(marker.c_str()))
  {
    frread();

    // read the whole thing
    while (strncmp(allfiles.actplace,"------",6)!=0)
    {
      int dobj;
      int ierr;
      int nodeid;
      frint(name.c_str(),&dobj,&ierr);
      if (ierr!=1)
        dserror("Cannot read %s", sectionname.c_str());
      frint("NODE",&nodeid,&ierr);
      if (ierr!=1)
        dserror("Cannot read %s", sectionname.c_str());
      topology[dobj-1].insert(nodeid-1);
      frread();
    }

    // copy all design object entries
    for (std::map<int,std::set<int> >::iterator i=topology.begin();
         i!=topology.end();
         ++i)
    {
      if (i->first >= static_cast<int>(dobj_fenode.size()))
      {
        dserror("Illegal design object number %d in section '%s'", i->first+1, sectionname.c_str());
      }

      // we copy from a std::set, thus the gids are sorted
      dobj_fenode[i->first].reserve(i->second.size());
      dobj_fenode[i->first].assign(i->second.begin(),i->second.end());
    }
  }
}

//----------------------------------------------------------------------
/// read a knotvector section (for isogeometric analysis)
//----------------------------------------------------------------------
void DatFileReader::ReadKnots(
  const int                            npatches,
  const int                            dim     ,
  Teuchos::RCP<DRT::NURBS::Knotvector> disknots
  ) const
{
  // make sure that we have some Knotvector obeject to fill
  if (disknots==Teuchos::null)
  {
    dserror("disknots should have been allocated before");
  }

  // this is a pointer to the knots of one patch in one direction
  // we will read them and put them
  vector<Teuchos::RCP<vector<double> > > patch_knots(dim);

  // open input file --- this is done on all procs
  ifstream file;

  file.open(filename_.c_str());

  // temporary strings
  string tmp;


  int filecount=0;

  // start to read something when read is true
  bool read=false;


  bool knotvectorsection=false;

  // index for number of patch
  int            npatch          = 0;
  // index for u/v/w
  int            actdim          =-1;
  // ints for the number of knots
  vector<int>    n_x_m_x_l(dim);
  // ints for patches degrees
  vector<int>    degree(dim);
  // a vector of strings holding the knotvectortypes read
  vector<string> knotvectortype(dim);

  // loop lines in file
  for (; file; ++filecount)
  {
    file >> tmp;

    // if this a new section
    if((tmp[0]=='-'&&tmp[1]=='-'))
    {

      // check whether it is the knotvectorsection
      string::size_type loc = tmp.rfind("KNOTVECTORS");

      if (loc != string::npos)
      {
        // if this is true, we are at the beginning of a knot section
	knotvectorsection=true;
	// there is nothing more to be done in this line
	continue;
      }
      else
      {
	knotvectorsection=false;

	// there is nothing more to be done in this line
	continue;
      }
    }

    if(knotvectorsection)
    {
      // check for a new patch
      string::size_type loc;

      loc = tmp.rfind("BEGIN");
      if (loc != string::npos)
      {
	file >> tmp;
	read=true;

	actdim=-1;

	for(int rr=0;rr<dim;++rr)
	{
	  patch_knots[rr]=rcp(new vector<double>);
	  (*(patch_knots[rr])).clear();
	}

	continue;
      }

      loc = tmp.rfind("ID");
      if (loc != string::npos)
      {
	// get ID of patch we are currently reading
	string str_npatch;
	file >> str_npatch;

	char* endptr = NULL;
	npatch=strtol(str_npatch.c_str(),&endptr,10);
	npatch--;

	continue;
      }

      loc = tmp.rfind("NUMKNOTS");
      if (loc != string::npos)
      {
	string str_numknots;
	file >> str_numknots;

	// new dimesion for knotvector
	actdim++;
	if(actdim>dim)
	{
	  dserror("too many knotvectors for (we only need dim)\n");
	}

	char* endptr = NULL;
	n_x_m_x_l[actdim]=strtol(str_numknots.c_str(),&endptr,10);

	continue;
      }

      loc = tmp.rfind("DEGREE");
      if (loc != string::npos)
      {
	string str_degree;
	file >> str_degree;

	char* endptr = NULL;
	degree[actdim]=strtol(str_degree.c_str(),&endptr,10);

	continue;
      }

      loc = tmp.rfind("TYPE");
      if (loc != string::npos)
      {
	string type;

	file >> type;
        knotvectortype[actdim]=type;

	continue;
      }

      loc = tmp.rfind("END");
      if (loc != string::npos)
      {
	for (int rr=0;rr<dim;++rr)
	{
	  disknots->SetKnots(
	    rr                ,
	    npatch            ,
	    degree[rr]        ,
	    n_x_m_x_l[rr]     ,
	    knotvectortype[rr],
	    patch_knots[rr]   );
	}
	file >> tmp;
	read=false;
	continue;
      }

      if(read)
      {
	char* endptr = NULL;

        double dv = strtod(tmp.c_str(), &endptr);

	(*(patch_knots[actdim])).push_back(dv);
      }
    }

  } // end loop through file

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::ParameterList& DatFileReader::FindSublist(string name, Teuchos::ParameterList& list)
{
  Teuchos::ParameterList* sublist = &list;

  for (string::size_type pos=name.find('/');
       pos!=string::npos;
       pos=name.find('/'))
  {
    sublist = &sublist->sublist(name.substr(0,pos));
    name = name.substr(pos+1);
  }

  return sublist->sublist(name);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DatFileReader::AddEntry(string key, string value, Teuchos::ParameterList& list)
{
  const char* v = value.c_str();
  char* endptr = NULL;

  // Try converging to int first. If the end pointer points to
  // the trailing zero, we are done.
  long int iv = strtol(v, &endptr, 10);
  if (*endptr=='\0')
  {
    list.set(key,static_cast<int>(iv));
  }
  else
  {
    double dv = strtod(v, &endptr);
    if (*endptr=='\0')
    {
      list.set(key,dv);
    }
    else
    {
      list.set(key,value);
    }
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DatFileReader::ReadDat()
{
  vector<string> exclude;

  exclude.push_back("--NODE COORDS");
  exclude.push_back("--STRUCTURE ELEMENTS");
  exclude.push_back("--FLUID ELEMENTS");
  exclude.push_back("--ALE ELEMENTS");
  exclude.push_back("--TRANSPORT ELEMENTS");

  int arraysize = 0;

  if (comm_->MyPID()==0)
  {
    ifstream file(filename_.c_str());
    if (not file)
      dserror("unable to open file: %s", filename_.c_str());

    list<string> content;
    bool ignoreline = false;
    string line;

    // loop all input lines
    while (getline(file, line))
    {
      // remove comments
      string::size_type loc = line.find("//");
      if (loc != string::npos)
      {
        line = line.substr(0,loc);
      }

      // remove trailing and leading whitespaces
      // compact internal whitespaces
      line = trim(line);

      // exclude all special sections
      // this includes the section header and all lines in that section
      if (ignoreline)
      {
        if (line.find("--")==0)
        {
          ignoreline = false;
        }
      }

      // Two sections to be ignored can follow each other. We need
      // independent tests.
      if (!ignoreline)
      {

        // remember all section positions
        if (line.find("--")==0)
        {
          for (vector<int>::size_type i=0; i<exclude.size(); ++i)
          {
            if (line.find(exclude[i]) != string::npos)
            {
              excludepositions_[exclude[i]] = file.tellg();
              ignoreline = true;
              break;
            }
          }
        }

      }

      // remember line
      if (!ignoreline && line.length() > 0)
      {
        content.push_back(line);
        arraysize += line.length()+1;
      }
    }

    // setup global variables
    numrows_ = static_cast<int>(content.size());

    // allocate space for copy of file
    inputfile_.resize(arraysize);

    // CAUTION: We allocate one more row pointer that necessary. This
    // pointer will be used to point to temporary lines when the
    // excluded section are read (on proc 0). Of course that's just
    // another EVIL HACK. Don't tell anybody.
    lines_.reserve(numrows_+1);

    // fill file buffer
    char* ptr = &inputfile_[0];
    for (list<string>::iterator i=content.begin(); i!=content.end(); ++i)
    {
      strcpy(ptr,i->c_str());
      lines_.push_back(ptr);
      ptr += i->length()+1;
    }

    if (ptr - &inputfile_[0] != static_cast<int>(arraysize))
      dserror("internal error in file read");

    // add the slot for the temporary line...
    lines_.push_back(NULL);
  }

  // Now lets do all the parallel setup. Afterwards all processors
  // have to be the same.

  if (comm_->NumProc()>1)
  {
    /* Now that we use a variable number of bytes per line we have to
     * communicate the buffer size as well. */
    comm_->Broadcast(&arraysize,1,0);
    comm_->Broadcast(&numrows_,1,0);

    if (comm_->MyPID()>0)
    {
      /*--------------------------------------allocate space for copy of file */
      inputfile_.resize(arraysize);
      lines_.reserve(numrows_+1);
    }

#ifdef PARALLEL

    // There are no char based functions available! Do it by hand!
    //comm_->Broadcast(&inputfile_[0],arraysize,0);

    const Epetra_MpiComm& mpicomm = dynamic_cast<const Epetra_MpiComm&>(*comm_);

    MPI_Bcast(&inputfile_[0], arraysize, MPI_CHAR, 0, mpicomm.GetMpiComm());
#else
    dserror("How did you get here? Go away!");
#endif

    /* We have not yet set the row pointers on procs > 0. So do it now. */
    if (comm_->MyPID()>0)
    {
      lines_.push_back(&inputfile_[0]);
      for (int i=0; i<arraysize; ++i)
      {
        if (inputfile_[i]=='\0')
        {
          lines_.push_back(&inputfile_[i+1]);
        }
      }

      if (static_cast<int>(lines_.size()) != numrows_+1)
        dserror("line count mismatch: %d lines expected but %d lines received",
                numrows_+1, lines_.size());
    }

#ifdef PARALLEL
    // distribute excluded section positions
    for (vector<int>::size_type i=0; i<exclude.size(); ++i)
      //comm_->Broadcast(&excludepositions_[exclude[i]],1,0);
      MPI_Bcast(&excludepositions_[exclude[i]],1,MPI_INT,0,mpicomm.GetMpiComm());
#endif
  }

  // Now finally find the section names. We have to do this on all
  // processors, so it cannot be done while reading.
  for (vector<char*>::size_type i=0; i<lines_.size()-1; ++i)
  {
    char* l = lines_[i];
    if (l and l[0]=='-' and l[1]=='-')
    {
      std::string line(l);

      // take the last "--" and all that follows as section name
      std::string::size_type loc = line.rfind("--");
      std::string sectionname = line.substr(loc);
      positions_[sectionname] = i;
    }
  }

  if (positions_.find("--END")==positions_.end())
    dserror("end section missing. incomplete dat file?");
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DatFileReader::DumpInput()
{
#if defined(DEBUG) || defined(OUTPUT_INPUT)
  if (comm_->MyPID()==0 and allfiles.out_err!=NULL)
  {
    fprintf(allfiles.out_err,
            "============================================================================\n"
            "broadcasted copy of input file:\n"
            "============================================================================\n"
      );
    for (unsigned i=0; i<lines_.size()-1; ++i)
    {
      fprintf(allfiles.out_err,"%s\n", lines_[i]);
    }
    fprintf(allfiles.out_err,
            "============================================================================\n"
            "end of broadcasted copy of input file\n"
            "============================================================================\n"
      );
    fflush(allfiles.out_err);
  }
#endif
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ElementReader::ElementReader(Teuchos::RCP<Discretization> dis,
                             const DRT::INPUT::DatFileReader& reader,
                             string sectionname)
  : name_(dis->Name()),
    reader_(reader),
    comm_(reader.Comm()),
    sectionname_(sectionname),
    dis_(dis)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ElementReader::Partition()
{
  const int myrank  = comm_->MyPID();
  const int numproc = comm_->NumProc();

  Epetra_Time time(*comm_);

#ifdef PARMETIS

  // - read global ids of elements of this discretization
  //   (this is one fully redundant vector for elements)
  // - determine a preliminary element distribution. The fully redundant 
  //   vector is trashed after the construction.
  // - define blocksizes for blocks of elements we read (not necessarily
  //   the same as it was used to construct the map --- we may have a
  //   smaller blocksize here).
  // - read elements of this discretization and distribute according 
  //   to a linear map. While reading, remember node gids an assemble 
  //   them into a second fully redundant vector (mapping vertex id->gid).
  //   In addition, we keep them in a fully redundant set (required by 
  //   node reader). Construct reverse lookup from gids to vertex ids. 
  //   Again, this is a global, fully redundant map!
  // - define preliminary linear distributed nodal row map
  // - determine adjacency array (i.e. the infos for the node graph)
  //   using the nodal row distribution and a round robin communication
  //   of element connectivity information.
  //   Use adjacency array to build an initial Crsgraph on the linear map.
  // - do partitioning using parmetis
  //   Results are distributed to other procs using two global vectors!
  // - build final nodal row map, export graph to the new map
  // 

  // --------------------------------------------------
  // - read global ids of elements of this discretization

  // vector of all global element ids 
  vector<int> eids;            
  int numele = 0;
  string inputfile_name = reader_.MyInputfileName();

  // all reading is done on proc 0
  if (myrank==0)
  {
    if (!reader_.MyOutputFlag())
    {
      cout << "Entering jumbo reading mode for " << name_ << " discretization ...\n"
           << "Read, create and partition elements      in....";
      fflush(stdout);
    }

    // open input file at the right position
    ifstream file(inputfile_name.c_str());
    file.seekg(reader_.ExcludedSectionPosition(sectionname_));

    // loop all element lines
    // Comments in the element section are not supported!

    // Ok. We do this twice here! The first time we just gather the
    // element numbers. With those we construct a preliminary element
    // row map.
    string line;
    for (int i=0; getline(file, line); ++i)
    {
      if (line.find("--")==0)
      {
        break;
      }
      else
      {
        istringstream t;
        t.str(line);
        int elenumber;
        t >> elenumber;
        elenumber -= 1;
        eids.push_back(elenumber);
      }
    }
    numele = static_cast<int>(eids.size());
  }

  // Simply allreduce the element ids
  comm_->Broadcast(&numele,1,0);
  eids.resize(numele);
  comm_->Broadcast(&eids[0],numele,0);

  // --------------------------------------------------
  // - determine a preliminary element distribution

  // number of element junks to split the reading process in
  // approximate block size (just a guess!)
  int nblock = numproc;
  int bsize = static_cast<int>(numele)/nblock;

  // create a simple (pseudo linear) map
  int mysize = bsize;
  if (myrank==numproc-1)
    mysize = numele-(numproc-1)*bsize;

  // construct the map
  roweles_ = rcp(new Epetra_Map(-1,mysize,&eids[myrank*bsize],0,*comm_));

  // throw away redundant vector of elements
  eids.clear();

  // --------------------------------------------------
  // - define blocksizes for blocks of elements we read

  // for block sizes larger than about 250000 elements (empirical value !)
  // the code sometimes hangs during ExportRowElements call for the
  // second block (block 1).
  // Therefore an upper limit of 200000 for bsize is ensured below.
  int maxblocksize = 200000;

  if (bsize > maxblocksize)
  {
    // without an additional increase of nblock by 1 the last block size
    // could reach a maximum value of (2*maxblocksize)-1, potentially
    // violating the intended upper limit!
    nblock = 1+ static_cast<int>(numele/maxblocksize);
    bsize = maxblocksize;
  }

  // --------------------------------------------------
  // - read elements of this discretization and distribute according 
  //   to a linear map. While reading, remember node gids an assemble 
  //   them into a second fully redundant vector.

  // open input file at correct position,
  // valid on proc 0 only!
  ifstream file;
  if (0==myrank)
  {
    file.open(inputfile_name.c_str());
    file.seekg(reader_.ExcludedSectionPosition(sectionname_));
  }
  string line;
  int filecount=0;
  bool endofsection = false;

  // note that the last block is special....
  for (int block=0; block<nblock; ++block)
  {
    if (not endofsection and 0==myrank)
    {
      int bcount=0;
      for (;getline(file, line); ++filecount)
      {
        if (line.find("--")==0)
        {
          // If we have an empty element section (or fewer elements
          // that processors) we cannot read on. But we cannot exit
          // the block loop beforehand either because the other
          // processors need to syncronize nblock times.
          endofsection = true;
          break;
        }
        else
        {
          istringstream t;
          t.str(line);
          int elenumber;
          string eletype;
          string distype;
          // read element id type and distype
          t >> elenumber >> eletype >> distype;
          elenumber -= 1;

          // Set the current row to the empty slot after the file rows
          // and store the current line. This way the elements can use
          // the normal fr* functions to read the line.
          // Of course this is a hack.
          allfiles.actrow = allfiles.numrows;
          allfiles.actplace = allfiles.input_file[allfiles.actrow] = const_cast<char*>(line.c_str());
          // let the factory create a matching empty element
          Teuchos::RCP<DRT::Element> ele = DRT::UTILS::Factory(eletype,distype,elenumber,0);
          // let this element read its input line
          ele->ReadElement();
          // add element to discretization
          dis_->AddElement(ele);

          // get the node ids of this element
          const int  numnode = ele->NumNode();
          const int* nodeids = ele->NodeIds();

          // all node gids of this element are inserted into a set of 
          // node ids --- it will be used later during reading of nodes
          // to add the node to one or more discretisations
          copy(nodeids, nodeids+numnode, inserter(nodes_, nodes_.begin()));

          ++bcount;
          if (block != nblock-1) // last block is different....
          {
            if (bcount==bsize)
            {
	      filecount++;
              break;
            }
          }
        }
      } // for (;getline(file, line); ++filecount)
    } // if (0==myrank)

    // export junk of elements to other processors as reflected in the linear
    // map roweles
    dis_->ExportRowElements(*roweles_);

  } // end loop blocks

  // global node ids --- this will be a fully redundant vector!
  int numnodes=0;
  vector<int> nids;           

  if (myrank==0)
  {
    // Reset fr* functions. Still required.
    frrewind();

    numnodes = (int)nodes_.size();

    // copy set content into nids vector
    nids.reserve(numnodes);
    copy(nodes_.begin(), nodes_.end(), back_inserter(nids));
  }

  // Simply allreduce the node ids --- the vector is ordered according
  // to the < operator from the set which was used to constuct it
  comm_->Broadcast(&numnodes,1,0);
  nids.resize(numnodes);
  comm_->Broadcast(&nids[0],numnodes,0);


  // We want to be able to read empty fields. If we have such a beast
  // just skip the partitioning.
  if (numnodes>0)
  {
    // construct reverse lookup from gids to vertex ids. Again, 
    // this is a global, fully redundant map!
    // We need this lookup for the construction of the parmetis 
    // adjacency array
    //
    //        
    //                         gidtoidx 
    //                      ------------->
    //                gid_i                vertexid
    //                      <-------------
    //                           nids
    //            
    map<int,int> gidtoidx;
    for (int i=0;i<numnodes;++i)
    {
      gidtoidx[nids[i]]=i;
    }
    
    if (myrank==0)
    {
      
      if (!reader_.MyOutputFlag())
      cout << time.ElapsedTime() << " secs\n";
      time.ResetStartTime();
      
      if (!reader_.MyOutputFlag())
      {
        cout << "Build initial node graph, call PARMETIS  in....";
        fflush(stdout);
      }
    }
    
    // --------------------------------------------------
    // - define distributed nodal row map
    
    // vertices and nodes --- the parmetis call will be based on a 
    // consecutive numbering of vertices. We will use a lookup 
    // for global node ids which will map the gid to its vertex id.
    // We need this for the construction of the adjacency array.
    //
    //         rownode gid   | (row)vertex id            vtxdist
    //                       |                          
    //      -----------------+----------------------------+---+
    //   +-        gid0      |      0.....................| 0 |
    //   |         gid1      |      1                     +---+
    //  p|         gid2      |      2                       ^
    //  r|         gid3      |      3                       |
    //  o|         gid4      |      4                       |
    //  c|         gid5      |      5                       |
    //  0|         gid6      |      6                       |
    //   |         gid7      |      7                       |
    //   |         gid8      |      8                       |
    //   +-        gid9      |      9                       v
    //      -----------------+----------------------------+---+ 
    //  .+-       gid10      |     10.....................| 10|
    //  p|        gid11      |     11                     +---+
    //  r|        gid12      |     12                       ^
    //  o|        gid13      |     13                       |
    //  c|        gid14      |     14                       |
    //  1|        gid15      |     15                       |
    //   +-       gid16      |     16                       v
    //      -----------------+----------------------------+---+
    //  p+-       gid17      |     17.....................| 17|
    //  r|        gid18      |     18                     +---+
    //  o|        gid19      |     19                       ^
    //  c|        gid20      |     20                       |
    //  2|        gid21      |     21                       |
    //   +-       gid22      |     22                       v
    //      ----------------------------------------------+---+
    //      ..............................................| 23|
    //                                                    +---+
    //
    
    // number of node id junks
    int nbsize = numnodes/nblock;
    
    // create a simple (pseudo linear) map for nodes
    int mynsize = nbsize;
    if (myrank==numproc-1)
    mynsize = numnodes-(numproc-1)*nbsize;
    
    // construct the initial linear node rowmap
    RCP<Epetra_Map> lin_noderowmap = rcp(new Epetra_Map(-1,mynsize,&nids[myrank*nbsize],0,*comm_));

    // remember my vertex distribution for the later parmetis call
    vector<int> vtxdist(numproc+1);
    for(int np=0;np<numproc;++np)
    {
      vtxdist[np]=np*nbsize;
    }
    vtxdist[numproc]=numnodes;

    // --------------------------------------------------
    // - determine adjacency array (i.e. the infos for the node graph)
    //   using the nodal row distribution and a round robin communication
    //   of element connectivity information

    // this is a set of gids of all nodes for which we have a 
    // connectivity information on this proc
    set<int> procnodes;
    
    // loop all eles on this proc and determine all gids for 
    // which we have some connectivity information
    for(int lid=0;lid<roweles_->NumMyElements();++lid)
    {
      int gid=roweles_->GID(lid);
      
      DRT::Element* ele=dis_->gElement(gid);
      
      // get the node ids of this element
      const int  numnode = ele->NumNode();
      const int* nodeids = ele->NodeIds();
      
      // all node gids of this element are inserted into a set of 
      // node ids
      copy(nodeids, nodeids+numnode, inserter(procnodes,procnodes.begin()));
    }
  
    // ---------------------
    // build a processor local node connectivity


    //     node gid contained in one of the elements 
    //                    on this proc
    //                         |
    //                         | lcon
    //                         |
    //                         v
    //    set of all gids of adjacent nodes on this proc
    //          
    map<int,set<int> > lcon;

    // construct empty local map
    set<int>::iterator procnode;

    for(procnode=procnodes.begin();procnode!=procnodes.end();++procnode)
    {
      lcon.insert(pair<int, set<int> >(*procnode,set<int> ()));
    }

    // loop all eles on this proc and construct the local 
    // connectivity information
    map<int,set<int> >::iterator gidinlcon;

    for(int lid=0;lid<roweles_->NumMyElements();++lid)
    {
      int gid=roweles_->GID(lid);
      
      DRT::Element* ele=dis_->gElement(gid);

      // get the node ids of this element
      const int  numnode = ele->NumNode();
      const int* nodeids = ele->NodeIds();

      // loop nodeids
      for(int rr=0;rr<numnode;++rr)
      {
        // find gid of rr in map
        gidinlcon=lcon.find(nodeids[rr]);
        
        if(gidinlcon==lcon.end())
        {
          dserror("GID %d should be already contained in the map",nodeids[rr]);
        }

        // add nodeids to local connectivity 
        copy(nodeids, 
             nodeids+numnode, 
             inserter((gidinlcon->second),(gidinlcon->second).begin()));
      }
    }

    //-------------------------
    // (incomplete) round robin loop to build complete node
    //  connectivity for local nodes across all processors
    
    //       node gid of one row node on this proc
    //                         |
    //                         | gcon
    //                         |
    //                         v
    //    set of all gids of adjacent nodes on all procs

    // prepare empty map
    map<int,set<int> > gcon;
    for(int j=0;j<lin_noderowmap->NumMyElements();++j)
    {
      int gid=lin_noderowmap->GID(j);
      
      gcon.insert(pair<int, set<int> >(gid,set<int> ()));
    }

    {
#ifdef PARALLEL
      // create an exporter for point to point comunication
      DRT::Exporter exporter(dis_->Comm());
      
      // necessary variables
      MPI_Request request;

      int         tag    =-1;
      int         frompid=-1;
      int         topid  =-1;
      int         length =-1;
      
      // define send and receive blocks
      vector<char> sblock;
      vector<char> rblock;
 
#endif
      
      for (int np=0;np<numproc;++np)
      {
        // in the first step, we cannot receive anything
        if(np >0) 
        {
#ifdef PARALLEL
          //----------------------
          // Unpack local graph from the receive block from the 
          // last proc

          // make sure that you do not think you received something if
          // you didn't
          if(rblock.empty()==false)
          {
            dserror("rblock not empty");
          }

          // receive from predecessor
          frompid=(myrank+numproc-1)%numproc;
          exporter.ReceiveAny(frompid,tag,rblock,length);
          
          if(tag!=(myrank+numproc-1)%numproc)
          {
            dserror("received wrong message (ReceiveAny)");
          }
          
          exporter.Wait(request);

          // for safety
          exporter.Comm().Barrier();
#endif

          // Unpack received block
          UnpackLocalConnectivity(lcon,rblock);

        }

        // -----------------------
        // add local connectivity passing by to global 
        // connectivity on this proc
      
        // loop this procs global connectivity
        map<int,set<int> >::iterator gidingcon;

        for(gidingcon=gcon.begin();gidingcon!=gcon.end();++gidingcon)
        {

          // search in (other) procs local connectivity
          gidinlcon=lcon.find(gidingcon->first);
          
          if(gidinlcon!=lcon.end())
          {
            // in this case we do have some additional 
            // connectivity info from the proc owning lcon

            copy((gidinlcon->second).begin(),
                 (gidinlcon->second).end(),
                 inserter((gidingcon->second),(gidingcon->second).begin()));
          }
        }

        // in the last step, we keep everything on this proc
        if(np < numproc-1)
        {
          //-------------------------
          // Pack local graph into block to send
          PackLocalConnectivity(lcon,sblock);

#ifdef PARALLEL
          // Send block to next proc.

          tag    =myrank;
          frompid=myrank;
          topid  =(myrank+1)%numproc;

          exporter.ISend(frompid,topid,
                         &(sblock[0]),sblock.size(),
                         tag,request);
#endif
        }
      }
    }
    
    lcon.clear();
  
    // --------------------------------------------------
    // - do partitioning using parmetis
    
    //    gid(i)           gid(i+1)                        index gids
    //      |^                 |^
    //  nids||             nids||
    //      ||gidtoidxd        ||gidtoidx        
    //      ||                 || 
    //      v|                 v|
    //   nodegid(i)       nodegid(i+1)                      node gids
    //      ^                 ^
    //      |                 |  
    //      | lin_noderowmap  | lin_noderowmap       
    //      |                 | 
    //      v                 v
    //      i                i+1                        local equivalent indices
    //      |                 |
    //      | xadj            | xadj
    //      |                 |
    //      v                 v
    //     +-+-+-+-+-+-+-+-+-+-+                -+-+-+
    //     | | | | | | | | | | | ............... | | |      adjncy
    //     +-+-+-+-+-+-+-+-+-+-+                -+-+-+
    //
    //     |     gid(i)'s    |   gid(i+1)'s
    //     |    neighbours   |   neighbours           (numbered by global equivalent 
    //                                                 indices, mapping by nids vector 
    //                                                 to global ids)
    //
  
    // xadj points from vertex index i to the index of the 
    // first adjacent vertex.
    vector<int> xadj(lin_noderowmap->NumMyElements()+1);

    // a list of adjacent nodes, adressed using xadj
    vector<int> adjncy;
    
    int count=0;
    xadj[0] = 0;

    for (int idx=0;idx<lin_noderowmap->NumMyElements();++idx)
    {
      // get global node id of rownode
      int  growid =lin_noderowmap->GID(idx);
    
      map<int,set<int> >::iterator gidingcon;
      set<int>::iterator           nbgid;

      gidingcon=gcon.find(growid);

      if(gidingcon!=gcon.end())
      {

        for(nbgid=(gidingcon->second).begin();nbgid!=(gidingcon->second).end();++nbgid)
        {

          if(*nbgid!=growid)
          {
            // reverse lookup to determine neighbours index
            int nbidx=gidtoidx[*nbgid];

            adjncy.push_back(nbidx);
            ++count;
          }
        }
      }
      else
      {
        dserror("GID %d not found in global connectivity during setup of adjncy",growid);
      }
      xadj[idx+1] = count;
    }

    // define a vector of nodeweights
    vector<int> vwgt(lin_noderowmap->NumMyElements());

    // at the moment, we use a constant weight distribution at this point
    for (int i=0; i<lin_noderowmap->NumMyElements(); ++i)
    {
      vwgt[i] = 1;
    }

    /*
      This is used to indicate if the graph is weighted. 
      wgtflag can take one of four values:
      0  No weights (vwgt and adjwgt are both NULL).
      1  Weights on the edges only (vwgt is NULL).
      2  Weights on the vertices only (adjwgt is NULL).
      3  Weights on both the vertices and edges.
    */
    int wgtflag=2;
    /*
      This is used to indicate the numbering scheme that
      is used for the vtxdist, xadj, adjncy, and part
      arrays. numflag can take one of two values:
      0 C-style numbering that starts from 0.
      1 Fortran-style numbering that starts from 1.
    */
    int numflag=0;
    /*
      This is used to specify the number of weights that 
      each vertex has. It is also the number of balance
      constraints that must be satisfied.
    */
    int ncon=1;
    /*
      This is used to specify the number of sub-domains 
      that are desired. Note that the number of sub-domains 
      is independent of the number of processors that call 
      this routine.
    */
    int npart=numproc;
    /*
      This is an array of integers that is used to pass 
      additional parameters for the routine. If options[0]=0,
      then the default values are used. If options[0]=1, 
      then the remaining two elements of options are 
      interpreted as follows:
      options[1]     This specifies the level of information 
                     to be returned during the execution of 
                     the algorithm. Timing information can be 
                     obtained by setting this to 1. Additional 
                     options for this parameter can be obtained 
                     by looking at the the file defs.h in the 
                     ParMETIS-Lib directory. The numerical values 
                     there should be added to obtain the correct 
                     value. The default value is 0.
      options[2]     This is the random number seed for the routine. 
                     The default value is 15.
    */
    int options[3] = { 0,0,15 };
    /*
      Upon successful completion, the number of edges that are cut 
      by the partitioning is written to this parameter.
    */
    int edgecut=0;
    /*
      This is an array of size equal to the number of locally-stored
      vertices. Upon successful completion the partition vector of 
      the locally-stored vertices is written to this array.
      Note that this vector will not be redistributed by metis, it 
      will just contain the information how to redistribute.
    */
    vector<int> part(lin_noderowmap->NumMyElements());
    /*
      An array of size ncon that is used to specify the imbalance 
      tolerance for each vertex weight, with 1 being perfect balance 
      and nparts being perfect imbalance. A value of 1.05 for each 
      of the ncon weights is recommended.
    */
    float ubvec =1.05;
    /*
      An array of size ncon x nparts that is used to specify 
      the fraction of vertex weight that should be distributed 
      to each sub-domain for each balance constraint. If all 
      of the sub-domains are to be of the same size for every 
      vertex weight, then each of the ncon x nparts elements 
      should be set to a value of 1/nparts. If ncon is greater 
      than one, the target sub-domain weights for each sub-domain
      are stored contiguously (similar to the vwgt array). Note 
      that the sum of all of the tpwgts for a give vertex weight 
      should be one.
    */
    vector<float> tpwgts(npart,1.0/(double)npart);
    /*
      This is a pointer to the MPI communicator of the processes that 
      call PARMETIS. 
    */
    MPI_Comm mpicomm=(dynamic_cast<const Epetra_MpiComm*>(&(dis_->Comm())))->Comm();

    ParMETIS_V3_PartKway(
      &(vtxdist[0]), 
      &(xadj   [0]), 
      &(adjncy [0]), 
      &(vwgt   [0]), 
      NULL         , 
      &wgtflag     , 
      &numflag     , 
      &ncon        , 
      &npart       , 
      &(tpwgts[0]) , 
      &ubvec       ,
      &(options[0]), 
      &edgecut     , 
      &(part[0])   , 
      &mpicomm);

    // for each vertex j, the proc number to which this vertex 
    // belongs was written to part[j]. Note that PARMETIS does 
    // not redistribute the graph according to the new partitioning, 
    // it simply computes the partitioning and writes it to 
    // the part array.

    // construct epetra graph on linear noderowmap
    Teuchos::RCP<Epetra_CrsGraph> graph = rcp(new Epetra_CrsGraph(Copy,*lin_noderowmap,108,false));

    for (map<int,set<int> >::iterator gid=gcon.begin();
         gid!=gcon.end();
         ++gid
      )
    {
      set<int>& rowset = gid->second;
      vector<int> row;
      row.reserve(rowset.size());
      row.assign(rowset.begin(),rowset.end());
      rowset.clear();

      int err = graph->InsertGlobalIndices(gid->first,row.size(),&row[0]);
      if (err<0) dserror("graph->InsertGlobalIndices returned %d",err);
    }
    
    // trash the stl based graph
    gcon.clear();

    // fill graph and optimize storage
    graph->FillComplete();
    graph->OptimizeStorage();

    // redundant, global vectors to broadcast partition results
    vector<int> lglobalpart(numnodes,0);
    vector<int> globalpart (numnodes,0);

    // insert proc ids to distribute to into global vector
    for (unsigned i=0;i<part.size();++i)
    {
      lglobalpart[gidtoidx[lin_noderowmap->GID(i)]]=part[i];
    }

    // finally broadcast partition results
    comm_->SumAll(&(lglobalpart[0]),&(globalpart[0]),numnodes);
    
    lglobalpart.clear();

    // --------------------------------------------------
    // - build final nodal row map 
    
    // resize part array to new size
    count=0;
    for (unsigned i=0; i<globalpart.size(); ++i)
    {
      if (globalpart[i]==myrank)
      {
        ++count;
      }
    }
    part.resize(count);

    // fill it with the new distribution from the partitioning results
    count=0;
    for (unsigned i=0; i<globalpart.size(); ++i)
    {
      if (globalpart[i]==myrank)
      {
        part[count] = nids[i];
        ++count;
      }
    }
  
    nids.clear();

    // create map with new layout
    Epetra_Map newmap(globalpart.size(),count,&part[0],0,*comm_);

    // create the output graph and export to it
    RefCountPtr<Epetra_CrsGraph> outgraph =
      rcp(new Epetra_CrsGraph(Copy,newmap,108,false));
    Epetra_Export exporter(graph->RowMap(),newmap);
    int err = outgraph->Export(*graph,exporter,Add);
    if (err<0) dserror("Graph export returned err=%d",err);

    //trash old graph
    graph=null;

    // call fill complete and optimize storage
    outgraph->FillComplete();
    outgraph->OptimizeStorage();

    // replace rownodes, colnodes with row and column maps from the graph
    // do stupid conversion from Epetra_BlockMap to Epetra_Map
    const Epetra_BlockMap& brow = outgraph->RowMap();
    const Epetra_BlockMap& bcol = outgraph->ColMap();
    rownodes_ = rcp(new Epetra_Map(brow.NumGlobalElements(),
                                   brow.NumMyElements(),
                                   brow.MyGlobalElements(),
                                   0,
                                   *comm_));
    colnodes_ = rcp(new Epetra_Map(bcol.NumGlobalElements(),
                                   bcol.NumMyElements(),
                                   bcol.MyGlobalElements(),
                                   0,
                                   *comm_));
  }
  else
  {
    // create empty node row map
    rownodes_ = rcp(new Epetra_Map(-1,nids.size(),&nids[0],0,*comm_));

    // We are empty. Just a proper initialization.
    colnodes_ = rownodes_;
  }

#else

  // - read global ids of elements of this discretization
  // - determine a preliminary element distribution
  // - read elements of this discretization and distribute
  // - construct nodal row map with everything on proc 0
  // - construct and partition graph
  // - build final nodal row map

  vector<int> eids;             // global element ids
  int numele = 0;
  string inputfile_name = reader_.MyInputfileName();

  // We read all elements at proc 0 and keep them. This way we do not
  // need to read (initialize) them again.
  //
  // If we happen to consume too much memory here, we have to
  // distribute at this point already and use parmetis.
  if (myrank==0)
  {
    if (!reader_.MyOutputFlag())
    {
      cout << "Entering jumbo reading mode for " << name_ << " discretization ...\n"
           << "Read, create and partition elements      in....";
      fflush(stdout);
    }

    // open input file at the right position
    ifstream file(inputfile_name.c_str());
    file.seekg(reader_.ExcludedSectionPosition(sectionname_));

    // loop all element lines
    // Comments in the element section are not supported!

    // Ok. We do this twice here! The first time we just gather the
    // element numbers. With those we construct a preliminary element
    // row map.
    string line;
    for (int i=0; getline(file, line); ++i)
    {
      if (line.find("--")==0)
      {
        break;
      }
      else
      {
        istringstream t;
        t.str(line);
        int elenumber;
        t >> elenumber;
        elenumber -= 1;
        eids.push_back(elenumber);
      }
    }
    numele = static_cast<int>(eids.size());
  }

  // Simply allreduce the element ids
  comm_->Broadcast(&numele,1,0);
  eids.resize(numele);
  comm_->Broadcast(&eids[0],numele,0);

  // number of element junks to split the reading process in
  // approximate block size (just a guess!)
  int nblock = numproc;
  int bsize = static_cast<int>(eids.size())/nblock;

  // create a simple (pseudo linear) map
  int mysize = bsize;
  if (myrank==numproc-1)
    mysize = eids.size()-(numproc-1)*bsize;

  roweles_ = rcp(new Epetra_Map(-1,mysize,&eids[myrank*bsize],0,*comm_));

  // for block sizes larger than about 250000 elements (empirical value !)
  // the code sometimes hangs during ExportRowElements call for the
  // second block (block 1).
  // Therefore an upper limit of 200000 for bsize is ensured below.

  int maxblocksize = 200000;

  if (bsize > maxblocksize)
  {
    // without an additional increase of nblock by 1 the last block size
    // could reach a maximum value of (2*maxblocksize)-1, potentially
    // violating the intended upper limit!
    nblock = 1+ static_cast<int>(eids.size()/maxblocksize);
    bsize = maxblocksize;
  }

  eids.clear();

  // For simplicity we remember all node ids of all elements on
  // processor 0. This way we can create the graph on processor 0 and
  // use serial metis. If this turns out to be too memory consuming,
  // we have to use parmetis and use the distributed elements from the
  // discretization.
  list<vector<int> > elementnodes;

  // open input file at correct position,
  // valid on proc 0 only!
  ifstream file;
  if (0==myrank)
  {
    file.open(inputfile_name.c_str());
    file.seekg(reader_.ExcludedSectionPosition(sectionname_));
  }
  string line;
  int filecount=0;
  bool endofsection = false;

  // note that the last block is special....
  for (int block=0; block<nblock; ++block)
  {
    if (not endofsection and 0==myrank)
    {
      int bcount=0;
      for (;getline(file, line); ++filecount)
      {
        if (line.find("--")==0)
        {
          // If we have an empty element section (or fewer elements
          // that processors) we cannot read on. But we cannot exit
          // the block loop beforehand either because the other
          // processors need to syncronize nblock times.
          endofsection = true;
          break;
        }
        else
        {
          istringstream t;
          t.str(line);
          int elenumber;
          string eletype;
          string distype;
          // read element id type and distype
          t >> elenumber >> eletype >> distype;
          elenumber -= 1;

          // Set the current row to the empty slot after the file rows
          // and store the current line. This way the elements can use
          // the normal fr* functions to read the line.
          // Of course this is a hack.
          allfiles.actrow = allfiles.numrows;
          allfiles.actplace = allfiles.input_file[allfiles.actrow] = const_cast<char*>(line.c_str());
          // let the factory create a matching empty element
          Teuchos::RCP<DRT::Element> ele = DRT::UTILS::Factory(eletype,distype,elenumber,0);
          // let this element read its input line
          ele->ReadElement();
          // add element to discretization
          dis_->AddElement(ele);

          //elements_[elenumber] = ele;

          // get the node ids of this element
          const int  numnode = ele->NumNode();
          const int* nodeids = ele->NodeIds();

          copy(nodeids, nodeids+numnode, inserter(nodes_, nodes_.begin()));
          elementnodes.push_back(vector<int>(nodeids, nodeids+numnode));

          ++bcount;
          if (block != nblock-1) // last block is different....
            if (bcount==bsize)
            {
	      filecount++;
              break;
            }
        }
      } // for (;getline(file, line); ++filecount)
    } // if (0==myrank)
    // export block of elements to other processors as reflected in the linear
    // map roweles

    // export junk of elements to other processors
    dis_->ExportRowElements(*roweles_);
  } 

  vector<int> nids;             // global node ids

  if (myrank==0)
  {
    // Reset fr* functions. Still required.
    frrewind();

    if (!reader_.MyOutputFlag())
      cout << time.ElapsedTime() << " secs\n";
    time.ResetStartTime();
    if (!reader_.MyOutputFlag())
    {
      cout << "Read, create and partition problem graph in....";
      fflush(stdout);
    }

    nids.reserve(nodes_.size());
    copy(nodes_.begin(), nodes_.end(), back_inserter(nids));
  }

  // create preliminary node row map
  rownodes_ = rcp(new Epetra_Map(-1,nids.size(),&nids[0],0,*comm_));
  nids.clear();

  // We want to be able to read empty fields. If we have such a beast
  // just skip the partitioning.
  if (rownodes_->NumGlobalElements()>0)
  {

#if 0
    if (myrank==0)
    {
      cout << "\n\nelementnodes: size=" << elementnodes.size() << endl;
      for (list<vector<int> >::iterator i=elementnodes.begin();
           i!=elementnodes.end();
           ++i)
      {
        copy(i->begin(), i->end(), ostream_iterator<int>(cout, " "));
        cout << endl;
      }
    }
#endif

#ifdef WE_DO_NOT_HAVE_MUCH_MEMORY_BUT_A_LOT_OF_TIME

    // construct graph
    Teuchos::RCP<Epetra_CrsGraph> graph = rcp(new Epetra_CrsGraph(Copy,*rownodes_,81,false));
    if (myrank==0)
    {
      for (list<vector<int> >::iterator i=elementnodes.begin();
           i!=elementnodes.end();
           ++i)
      {
        // get the node ids of this element
        int  numnode = static_cast<int>(i->size());
        int* nodeids = &(*i)[0];

        // loop nodes and add this topology to the row in the graph of every node
        for (int i=0; i<numnode; ++i)
        {
          int err = graph->InsertGlobalIndices(nodeids[i],numnode,nodeids);
          if (err<0) dserror("graph->InsertGlobalIndices returned %d",err);
        }
      }
    }

#else

    // No need to test for myrank==0 as elementnodes is filled on proc 0 only.

    // build the graph ourselves
    vector<set<int> > localgraph(rownodes_->NumMyElements());
    for (list<vector<int> >::iterator i=elementnodes.begin();
         i!=elementnodes.end();
         ++i)
    {
      // get the node ids of this element
      int  numnode = static_cast<int>(i->size());
      int* nodeids = &(*i)[0];

      // loop nodes and add this topology to the row in the graph of every node
      for (int n=0; n<numnode; ++n)
      {
        int nodelid = rownodes_->LID(nodeids[n]);
        copy(nodeids,
             nodeids+numnode,
             inserter(localgraph[nodelid],
                      localgraph[nodelid].begin()));
      }
    }

    elementnodes.clear();

    // fill exact entries per row vector
    // this will really speed things up for long lines
    vector<int> entriesperrow;
    entriesperrow.reserve(rownodes_->NumMyElements());

    transform(localgraph.begin(),
              localgraph.end(),
              back_inserter(entriesperrow),
              mem_fun_ref(&set<int>::size));

    // construct graph
    Teuchos::RCP<Epetra_CrsGraph> graph = rcp(new Epetra_CrsGraph(Copy,*rownodes_,&entriesperrow[0],false));

    entriesperrow.clear();

    for (unsigned i = 0; i<localgraph.size(); ++i)
    {
      set<int>& rowset = localgraph[i];
      vector<int> row;
      row.reserve(rowset.size());
      row.assign(rowset.begin(),rowset.end());
      rowset.clear();

      int err = graph->InsertGlobalIndices(rownodes_->GID(i),row.size(),&row[0]);
      if (err<0) dserror("graph->InsertGlobalIndices returned %d",err);
    }

    localgraph.clear();

#endif

    elementnodes.clear();

    // finalize construction of this graph
    int err = graph->FillComplete(*rownodes_,*rownodes_);
    if (err) dserror("graph->FillComplete returned %d",err);

    // partition graph using metis
    Epetra_Vector weights(graph->RowMap(),false);
    weights.PutScalar(1.0);
    graph = DRT::UTILS::PartGraphUsingMetis(*graph,weights);

    // replace rownodes, colnodes with row and column maps from the graph
    // do stupid conversion from Epetra_BlockMap to Epetra_Map
    const Epetra_BlockMap& brow = graph->RowMap();
    const Epetra_BlockMap& bcol = graph->ColMap();
    rownodes_ = rcp(new Epetra_Map(brow.NumGlobalElements(),
                                   brow.NumMyElements(),
                                   brow.MyGlobalElements(),
                                   0,
                                   *comm_));
    colnodes_ = rcp(new Epetra_Map(bcol.NumGlobalElements(),
                                   bcol.NumMyElements(),
                                   bcol.MyGlobalElements(),
                                   0,
                                   *comm_));

    // At this point we have a good guess about the system matrix bandwidth.
    //graph->MaxNumIndices();

    graph = null;
  }
  else
  {
    // We are empty. Just a proper initialization.
    colnodes_ = rownodes_;
  }

#endif // not PARMETIS

  // now we have all elements in a linear map roweles
  // build reasonable maps for elements from the
  // already valid and final node maps
  // note that nothing is actually redistributed in here
  dis_->BuildElementRowColumn(*rownodes_,*colnodes_,roweles_,coleles_);

  // we can now export elements to resonable row element distribution
  dis_->ExportRowElements(*roweles_);

  // export to the column map / create ghosting of elements
  dis_->ExportColumnElements(*coleles_);

  if (comm_->MyPID()==0 && reader_.MyOutputFlag() == 0)
  {
    cout << time.ElapsedTime() << " secs\n";
    fflush(stdout);
  }

  return;
} // end Partition()


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ElementReader::Complete()
{
  const int myrank  = comm_->MyPID();

  Epetra_Time time(*comm_);

  if (!myrank && !reader_.MyOutputFlag())
  {
    cout << "Complete discretization ";
    printf("%-16s",name_.c_str());
    cout << " in...." << flush;
  }

  int err = dis_->FillComplete(false,false,false);
  if (err)
    dserror("dis_->FillComplete() returned %d",err);

  if (!myrank && !reader_.MyOutputFlag())
  {
    cout << time.ElapsedTime() << " secs" << endl;
  }

  int numproc=comm_->NumProc();
  if(numproc>1)
  {
    vector<int> my_n_nodes     (numproc,0);
    vector<int>    n_nodes     (numproc,0);
    vector<int> my_n_ghostnodes(numproc,0);
    vector<int>    n_ghostnodes(numproc,0);
    vector<int> my_n_elements  (numproc,0);
    vector<int>    n_elements  (numproc,0);
    vector<int> my_n_ghostele  (numproc,0);
    vector<int>    n_ghostele  (numproc,0);

    my_n_nodes     [myrank]=dis_->NumMyRowNodes();
    my_n_ghostnodes[myrank]=dis_->NumMyColNodes()-my_n_nodes[myrank];
    my_n_elements  [myrank]=dis_->NumMyRowElements();
    my_n_ghostele  [myrank]=dis_->NumMyColElements()-my_n_elements[myrank];

    dis_->Comm().SumAll(&my_n_nodes     [0],&n_nodes     [0],numproc);
    dis_->Comm().SumAll(&my_n_ghostnodes[0],&n_ghostnodes[0],numproc);
    dis_->Comm().SumAll(&my_n_elements  [0],&n_elements  [0],numproc);
    dis_->Comm().SumAll(&my_n_ghostele  [0],&n_ghostele  [0],numproc);

    if(myrank==0)
    {
      cout << endl;
      printf("   +-----+---------------+--------------+-----------------+----------------+\n");
      printf("   | PID |  n_rownodes   | n_ghostnodes |  n_rowelements  |   n_ghostele   |\n");
      printf("   +-----+---------------+--------------+-----------------+----------------+\n");
      for(int npid=0;npid<numproc;++npid)
      {
        printf("   | %3d | %13d | %12d | %15d | %14d |\n",npid,n_nodes[npid],n_ghostnodes[npid],n_elements[npid],n_ghostele[npid]);
        printf("   +-----+---------------+--------------+-----------------+----------------+\n");
      }
      cout << endl;
    }
  
  }

}

#ifdef PARMETIS

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ElementReader::PackLocalConnectivity(
  map<int,set<int> > & lcon,
  vector<char>       & sblock
  )
{
  map<int,set<int> >::iterator gidinlcon;
  set<int>::iterator           adjacentgid;
  
  sblock.clear();
  // size (number of nodes we have a connectivity for)
  int size=lcon.size();
  DRT::ParObject::AddtoPack(sblock,size);
    
  for(gidinlcon=lcon.begin();gidinlcon!=lcon.end();++gidinlcon)
  {
    // add node gid we store the connectivity for
    DRT::ParObject::AddtoPack(sblock,gidinlcon->first);

    // add number of nodes adjacent to this one
    DRT::ParObject::AddtoPack(sblock,(gidinlcon->second).size());
   
    // add list of neighbours to this node
    for(adjacentgid =(gidinlcon->second).begin();
        adjacentgid!=(gidinlcon->second).end();
        ++adjacentgid)
    {
      DRT::ParObject::AddtoPack(sblock,*adjacentgid);
    }
  }
  
  lcon.clear();
  
  return;
} // end Pack_lcon

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ElementReader::UnpackLocalConnectivity(
  map<int,set<int> > & lcon  ,
  vector<char>       & rblock
  )
{
  lcon.clear();

  // position to extract
  int position = 0;

  // extract size (number of nodes we have a connectivity for)
  int size=0;
  DRT::ParObject::ExtractfromPack(position,rblock,size);

  for(int i=0;i<size;++i)
  {
    // extract node gid we store the connectivity for
    int gid=-1;
    DRT::ParObject::ExtractfromPack(position,rblock,gid);

    if(gid<0)
    {
      dserror("Unable to unpack a proper gid");
    }

    // extract number of adjacent nodes
    int numnb=0;
    DRT::ParObject::ExtractfromPack(position,rblock,numnb);
    if(numnb<1)
    {
      dserror("Everybody should have at least one unpackable neighbour");
    }

    set<int> neighbourset;

    // extract all adjacent nodes and feed them into the set
    for(int j=0;j<numnb;++j)
    {
      int nbgid=0;
      DRT::ParObject::ExtractfromPack(position,rblock,nbgid);
      neighbourset.insert(nbgid);
    }
    
    // add this node connectivity to local connectivity map
    lcon.insert(pair<int,set<int> >(gid,neighbourset));
  }

  // trash receive block
  rblock.clear();

  return;
}

#endif

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
NodeReader::NodeReader(const DRT::INPUT::DatFileReader& reader, string sectionname)
  : reader_(reader), comm_(reader.Comm()), sectionname_(sectionname)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Discretization> > NodeReader::FindDisNode(int nodeid)
{
  std::vector<Teuchos::RCP<DRT::Discretization> > v;
  for (unsigned i=0; i<ereader_.size(); ++i)
  {
    if (ereader_[i]->HasNode(nodeid))
    {
      v.push_back(ereader_[i]->MyDis());
    }
  }
  return v;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void NodeReader::Read()
{
  const int myrank  = comm_->MyPID();
  const int numproc = comm_->NumProc();
  string inputfile_name = reader_.MyInputfileName();

  int numnodes = 0;
  for (unsigned i=0; i<ereader_.size(); ++i)
  {
    ereader_[i]->Partition();
    numnodes += ereader_[i]->rownodes_->NumGlobalElements();
  }

  Epetra_Time time(*comm_);

  if (!myrank && !reader_.MyOutputFlag())
  {
    cout << "Read, create and partition nodes         in...." << flush;
  }

  // We will read the nodes block wise. we will use one block per processor
  // so the number of blocks is numproc
  // determine a rough blocksize
  const int nblock = numproc;
  int bsize = max(numnodes/nblock, 1);

  // open input file at the right position
  // note that stream is valid on proc 0 only!
  ifstream file;
  if (myrank==0)
  {
    file.open(inputfile_name.c_str());
    file.seekg(reader_.ExcludedSectionPosition(sectionname_));
  }
  string tmp;

  // note that the last block is special....
  int filecount=0;
  for (int block=0; block<nblock; ++block)
  {
    if (0==myrank)
    {
      int bcount=0;
      for (; file; ++filecount)
      {
        file >> tmp;

        if (tmp=="NODE")
        {
          double coords[3];
          int nodeid;
          file >> nodeid >> tmp >> coords[0] >> coords[1] >> coords[2];
          nodeid--;
          if (nodeid != filecount)
            dserror("Reading of nodes failed: Nodes must be numbered consecutive!!");
          if (tmp!="COORD")
            dserror("failed to read node %d",nodeid);
          std::vector<Teuchos::RCP<DRT::Discretization> > diss = FindDisNode(nodeid);
          for (unsigned i=0; i<diss.size(); ++i)
          {
            Teuchos::RCP<DRT::Discretization> dis = diss[i];
            // create node and add to discretization
            Teuchos::RCP<DRT::Node> node = rcp(new DRT::Node(nodeid,coords,myrank));
            dis->AddNode(node);
          }
          ++bcount;
          if (block != nblock-1) // last block takes all the rest
            if (bcount==bsize)   // block is full
            {
              ++filecount;
              break;
            }
        }
        else if (tmp=="CP")
        {
          // read control points for isogeometric analysis
          double coords[3];
          double weight;

          int cpid;
          file >> cpid >> tmp >> coords[0] >> coords[1] >> coords[2] >> weight;
          cpid--;
          if (cpid != filecount)
            dserror("Reading of control points failed: They must be numbered consecutive!!");
          if (tmp!="COORD")
            dserror("failed to read control point %d",cpid);
          std::vector<Teuchos::RCP<DRT::Discretization> > diss = FindDisNode(cpid);
          for (unsigned i=0; i<diss.size(); ++i)
          {
            Teuchos::RCP<DRT::Discretization> dis = diss[i];
            // create node/control point and add to discretization
            Teuchos::RCP<DRT::NURBS::ControlPoint> node = rcp(new DRT::NURBS::ControlPoint(cpid,coords,weight,myrank));
            dis->AddNode(node);
          }
          ++bcount;
          if (block != nblock-1) // last block takes all the rest
            if (bcount==bsize)   // block is full
            {
              ++filecount;
              break;
            }
        }
        else if (tmp.find("--")==0)
          break;
        else
          dserror("unexpected word '%s'",tmp.c_str());
      } // for (filecount; file; ++filecount)
    } // if (0==myrank)

    // export block of nodes to other processors as reflected in rownodes,
    // changes ownership of nodes
    for (unsigned i=0; i<ereader_.size(); ++i)
    {
      ereader_[i]->dis_->ExportRowNodes(*ereader_[i]->rownodes_);
    }

  } // for (int block=0; block<nblock; ++block)

  // last thing to do here is to produce nodal ghosting/overlap
  for (unsigned i=0; i<ereader_.size(); ++i)
  {
    ereader_[i]->dis_->ExportColumnNodes(*ereader_[i]->colnodes_);
  }

  if (!myrank && !reader_.MyOutputFlag())
  {
    cout << time.ElapsedTime() << " secs\n" << endl;
  }
  
  for (unsigned i=0; i<ereader_.size(); ++i)
  {
    ereader_[i]->Complete();
  }
}

}
}

#endif
