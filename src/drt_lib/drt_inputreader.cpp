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
  const int                              dim     ,
  const string                           name    ,
  Teuchos::RCP<DRT::NURBS::Knotvector> & disknots
  ) const
{

  // io to shell
  const int myrank  = comm_->MyPID();

  Epetra_Time time(*comm_);

  if (myrank==0)
  {
    if (!MyOutputFlag())
    {
      cout << "Reading knot vectors for " << name << " discretization :\n";
      fflush(stdout);
    }
  }

  // number of patches to be determined
  int  npatches  = 0;

  //--------------------------------------------------------------------
  //--------------------------------------------------------------------
  //             first, determine number of patches
  //--------------------------------------------------------------------
  //--------------------------------------------------------------------
  {
    // open input file --- this is done on all procs
    ifstream file;
    file.open(filename_.c_str());

    // temporary string
    string tmp;

    // flag indicating knot vector section in input
    bool knotvectorsection=false;

    // loop lines in file
    for (; file;)
    {
      // read piece of file until next seperator (whitespace, newline)
      file >> tmp;

      // if this a new section, i.e. starts like ------
      if((tmp[0]=='-'&&tmp[1]=='-'))
      {
        // check whether it is the knotvectorsection
        string::size_type loc=string::npos;

        // only the knotvector section of this discretisation
        // type is of interest
        if(name=="fluid")
        {
          loc= tmp.rfind("FLUID");
        }
        else if(name=="structure")
        {
          loc= tmp.rfind("STRUCTURE");
        }
        else if(name=="ale")
        {
          loc= tmp.rfind("ALE");
        }
        else if(name=="scatra")
        {
          loc= tmp.rfind("TRANSPORT");
        }
        else
        {
          dserror("Unknown discretization name for knotvector input\n");
        }
        if (loc == string::npos)
        {
          knotvectorsection=false;

          // there is nothing more to be done in this line
          continue;
        }
        else
        {
          // continue reading of second keyword
          file >> tmp;

          // check whether second keyword is knotvector
          loc= tmp.rfind("KNOTVECTORS");

          if (loc != string::npos)
          {
            // if this is true, we are at the beginning of a
            // knot section
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
      }

      // count number of patches in knotvector section of
      // this discretisation
      if(knotvectorsection)
      {
        // check for a new patch
        string::size_type loc;

        loc = tmp.rfind("ID");
        if (loc != string::npos)
        {
          // increase number of patches
          npatches++;

          continue;
        }
      }
    } // end loop through file
  }

  if (myrank==0)
  {
    if (!MyOutputFlag())
    {
      printf("                        %8d patches",npatches);
      fflush(stdout);
    }
  }


  //--------------------------------------------------------------------
  //--------------------------------------------------------------------
  //                alloc knotvector object to fill
  //--------------------------------------------------------------------
  //--------------------------------------------------------------------

  // allocate knotvector for this dis
  disknots=Teuchos::rcp (new DRT::NURBS::Knotvector(dim,npatches));

  // make sure that we have some Knotvector object to fill
  if (disknots==Teuchos::null)
  {
    dserror("disknots should have been allocated before");
  }

  //--------------------------------------------------------------------
  //--------------------------------------------------------------------
  //                finally read knotvector section
  //--------------------------------------------------------------------
  //--------------------------------------------------------------------
  {
    // this is a pointer to the knots of one patch in one direction
    // we will read them and put them
    vector<Teuchos::RCP<vector<double> > > patch_knots(dim);

    // open input file --- this is done on all procs
    ifstream file;

    file.open(filename_.c_str());

    // temporary string
    string tmp;

    // start to read something when read is true
    bool read=false;

    // flag indicating knot vector section in input
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

    // count for sanity check
    int            count_read=0;
    vector<int>    count_vals(dim);

    // loop lines in file
    for (; file;)
    {
      file >> tmp;

      // if this a new section
      if((tmp[0]=='-'&&tmp[1]=='-'))
      {
        // check whether it is the knotvectorsection
        string::size_type loc=string::npos;

        // only the knotvector section of this discretisation
        // type is of interest
        if(name=="fluid")
        {
          loc= tmp.rfind("FLUID");
        }
        else if(name=="structure")
        {
          loc= tmp.rfind("STRUCTURE");
        }
        else if(name=="ale")
        {
          loc= tmp.rfind("ALE");
        }
        else if(name=="scatra")
        {
          loc= tmp.rfind("TRANSPORT");
        }
        else
        {
          dserror("Unknown discretization name for knotvector input\n");
        }
        if (loc == string::npos)
        {
          knotvectorsection=false;

          // there is nothing more to be done in this line
          continue;
        }
        else
        {
          // continue reading of second keyword
          file >> tmp;

          // check whether second keyword is knotvector
          loc= tmp.rfind("KNOTVECTORS");

          if (loc != string::npos)
          {
            // if this is true, we are at the beginning of a
            // knot section
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
      }

      // do reading in knotvecor section
      if(knotvectorsection)
      {
        // check for a new patch
        string::size_type loc=string::npos;

        loc = tmp.rfind("BEGIN");
        if (loc != string::npos)
        {
          file >> tmp;

          // activate reading
          read=true;

          actdim=-1;

          // create vectors for knots in this patch
          for(int rr=0;rr<dim;++rr)
          {
            patch_knots[rr]=rcp(new vector<double>);
            (*(patch_knots[rr])).clear();
          }

          // reset counter for knot values
          for(int rr=0;rr<dim;rr++)
          {
            count_vals[rr]=0;
          }

          continue;
        }

        // get ID of patch we are currently reading
        loc = tmp.rfind("ID");
        if (loc != string::npos)
        {
          string str_npatch;
          file >> str_npatch;

          char* endptr = NULL;
          npatch=strtol(str_npatch.c_str(),&endptr,10);
          npatch--;

          continue;
        }

        // get number of knots in the knotvector direction
        // we are currently reading
        loc = tmp.rfind("NUMKNOTS");
        if (loc != string::npos)
        {
          string str_numknots;
          file >> str_numknots;

          // increase dimesion for knotvector (i.e. next time
          // we'll fill the following knot vector)
          actdim++;
          if(actdim>dim)
          {
            dserror("too many knotvectors (we only need dim)\n");
          }

          char* endptr = NULL;
          n_x_m_x_l[actdim]=strtol(str_numknots.c_str(),&endptr,10);

          continue;
        }

        // get number of bspline polinomial associated with
        // knots in this direction
        loc = tmp.rfind("DEGREE");
        if (loc != string::npos)
        {
          string str_degree;
          file >> str_degree;

          char* endptr = NULL;
          degree[actdim]=strtol(str_degree.c_str(),&endptr,10);

          continue;
        }

        // get type of knotvector (interpolated or periodic)
        loc = tmp.rfind("TYPE");
        if (loc != string::npos)
        {
          string type;

          file >> type;
          knotvectortype[actdim]=type;

          continue;
        }

        // locate end of patch
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
          // stop reading of knot values if we are here
          read=false;

          for(int rr=0;rr<dim;rr++)
          {
            if(n_x_m_x_l[rr]!=count_vals[rr])
            {
              dserror("not enough knots read in dim %d (%d!=NUMKNOTS=%d)\n",
                      rr            ,
                      count_vals[rr],
                      n_x_m_x_l[rr]);
            }
          }

          // count for sanity check
          count_read++;

          continue;
        }

        //  reading of knot values if read is true and no
        // other keyword was found
        if(read)
        {
          char* endptr = NULL;

          double dv = strtod(tmp.c_str(), &endptr);

          // count for sanity check
          count_vals[actdim]++;

          (*(patch_knots[actdim])).push_back(dv);
        }
      }
    } // end loop through file

    if(count_read!=npatches)
    {
      dserror("wasn't able to read enough patches\n");
    }
  }

  if (myrank==0)
  {
    if (!MyOutputFlag())
    {
      cout << " in...." << time.ElapsedTime() << " secs\n";

      time.ResetStartTime();
      fflush(stdout);
    }
  }
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
  exclude.push_back("--ARTERY ELEMENTS");
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
      if (positions_.find(sectionname)!=positions_.end())
        dserror("section '%s' defined more than once", sectionname.c_str());
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

#ifndef PARMETIS
  // For simplicity we remember all node ids of all elements on
  // processor 0. This way we can create the graph on processor 0 and
  // use serial metis. If this turns out to be too memory consuming,
  // we have to use parmetis and use the distributed elements from the
  // discretization.
  list<vector<int> > elementnodes;
#endif

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
#ifndef PARMETIS
          elementnodes.push_back(vector<int>(nodeids, nodeids+numnode));
#endif

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

  // create preliminary node row map
  rownodes_ = rcp(new Epetra_Map(-1,nids.size(),&nids[0],0,*comm_));

  // We want to be able to read empty fields. If we have such a beast
  // just skip the partitioning.
  if (rownodes_->NumGlobalElements()>0)
  {
#ifdef PARMETIS

    // Simply allreduce the node ids --- the vector is ordered according
    // to the < operator from the set which was used to constuct it
    comm_->Broadcast(&numnodes,1,0);
    nids.resize(numnodes);
    comm_->Broadcast(&nids[0],numnodes,0);

    DRT::UTILS::PartUsingParMetis(dis_,roweles_,rownodes_,colnodes_,nids,nblock,
                                  comm_,time,not reader_.MyOutputFlag());

#else
    nids.clear();
    DRT::UTILS::PartUsingMetis(rownodes_,colnodes_,elementnodes,comm_);
#endif
  }
  else
  {
    // We are empty. Just a proper initialization.
    colnodes_ = rownodes_;
  }

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

  DRT::UTILS::PrintParallelDistribution(*dis_);

}

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
    cout << "Read, create and partition nodes " << flush;
#ifdef PARMETIS
    cout << "block " << flush;
#else
    cout << "        " << flush;
#endif
  }

  // We will read the nodes block wise. we will use one block per processor
  // so the number of blocks is numproc
  // determine a rough blocksize
  int nblock = numproc;
  int bsize = max(numnodes/nblock, 1);

  // for block sizes larger than about 50000 elements (empirical value !)
  // the code sometimes hangs during ExportRowElements
  // Therefore an upper limit for bsize is ensured below.
  int maxblocksize = 50000;

  if (bsize > maxblocksize)
  {
    // without an additional increase of nblock by 1 the last block size
    // could reach a maximum value of (2*maxblocksize)-1, potentially
    // violating the intended upper limit!
    nblock = 1+ static_cast<int>(numnodes/maxblocksize);
    bsize = maxblocksize;
  }

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
#ifdef PARMETIS
      if (!reader_.MyOutputFlag())
      {
        printf("%d",block);
        if(block != nblock-1)
        {
          printf(",");
        }
        else
        {
          printf(" ");
        }
        fflush(stdout);
      }
#endif

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
    cout << "in...." << time.ElapsedTime() << " secs\n" << endl;
  }

  for (unsigned i=0; i<ereader_.size(); ++i)
  {
    ereader_[i]->Complete();
  }
} // NodeReader::Read

}
}

#endif
