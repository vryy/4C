/*----------------------------------------------------------------------*/
/*! \file

\brief Internal classes to read elements and nodes

\level 1


*/
/*----------------------------------------------------------------------*/


#include "baci_io_inputreader.hpp"

#include "baci_comm_utils.hpp"
#include "baci_global_data.hpp"
#include "baci_io_control.hpp"
#include "baci_io_linedefinition.hpp"
#include "baci_io_pstream.hpp"
#include "baci_io_utils_reader.hpp"
#include "baci_linalg_utils_densematrix_communication.hpp"
#include "baci_nurbs_discret_knotvector.hpp"

#include <Teuchos_Time.hpp>

#include <sstream>
#include <utility>

#ifdef BACI_TRAP_FE
#include <cfenv>
#endif /* BACI_TRAP_FE */


BACI_NAMESPACE_OPEN

namespace
{
  constexpr double tolerance_n = 1.0e-14;

  /*----------------------------------------------------------------------*/
  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& FindSublist(std::string name, Teuchos::ParameterList& list)
  {
    Teuchos::ParameterList* sublist = &list;

    for (std::string::size_type pos = name.find('/'); pos != std::string::npos;
         pos = name.find('/'))
    {
      sublist = &sublist->sublist(name.substr(0, pos));
      name = name.substr(pos + 1);
    }

    return sublist->sublist(name);
  }


}  // namespace

namespace INPUT
{
  /*----------------------------------------------------------------------*/
  /*----------------------------------------------------------------------*/
  DatFileReader::DatFileReader(std::string filename)
      : filename_(std::move(filename)), comm_(Teuchos::null), numrows_(0), outflag_(0)
  {
  }

  /*----------------------------------------------------------------------*/
  /*----------------------------------------------------------------------*/
  DatFileReader::DatFileReader(std::string filename, Teuchos::RCP<Epetra_Comm> comm, int outflag)
      : filename_(std::move(filename)), comm_(std::move(comm)), outflag_(outflag)
  {
    ReadDat();
  }


  /*----------------------------------------------------------------------*/
  /*----------------------------------------------------------------------*/
  std::string DatFileReader::MyInputfileName() const { return filename_; }


  /*----------------------------------------------------------------------*/
  /*----------------------------------------------------------------------*/
  int DatFileReader::MyOutputFlag() const { return outflag_; }


  /*----------------------------------------------------------------------*/
  /*----------------------------------------------------------------------*/
  std::ifstream::pos_type DatFileReader::ExcludedSectionPosition(const std::string& section) const
  {
    auto i = excludepositions_.find(section);
    if (i == excludepositions_.end())
    {
      return -1;
    }
    return i->second.first;
  }


  /*----------------------------------------------------------------------*/
  /*----------------------------------------------------------------------*/
  unsigned int DatFileReader::ExcludedSectionLength(const std::string& section) const
  {
    auto i = excludepositions_.find(section);
    if (i == excludepositions_.end())
    {
      return 0;
    }
    return i->second.second;
  }


  /*----------------------------------------------------------------------*/
  /*----------------------------------------------------------------------*/
  bool DatFileReader::ReadSection(std::string name, Teuchos::ParameterList& list)
  {
    if (name.length() < 3 or name[0] != '-' or name[1] != '-')
      dserror("Illegal section name '%s'", name.c_str());

    // The section name is desired from outside. Thus, we consider it as valid
    knownsections_[name] = true;

    Teuchos::ParameterList& sublist = FindSublist(name.substr(2), list);

    if (positions_.find(name) == positions_.end()) return false;

    for (size_t pos = positions_[name] + 1; pos < lines_.size(); ++pos)
    {
      const std::string line = lines_[pos];

      // If the line starts with dashes we found the beginning of the next section and terminate
      // the read process of the current section.
      if (line[0] == '-' and line[1] == '-')
      {
        break;
      }

      const auto& [key, value] = ReadKeyValue(line);

      AddEntry(key, value, sublist);
    }

    return true;
  }


  /*----------------------------------------------------------------------*/
  /*----------------------------------------------------------------------*/
  std::vector<const char*> DatFileReader::Section(const std::string& name)
  {
    // The section name is desired from outside. Thus, we consider it as valid
    knownsections_[name] = true;

    std::vector<const char*> sec;

    auto i = positions_.find(name);
    if (i != positions_.end())
    {
      for (size_t pos = i->second + 1; pos < lines_.size(); ++pos)
      {
        const char* line = lines_[pos];
        if (line[0] == '-' and line[1] == '-')
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
  void DatFileReader::ReadDesign(
      const std::string& name, std::vector<std::vector<int>>& dobj_fenode)
  {
    std::map<int, std::set<int>> topology;

    std::string sectionname = name + "-NODE TOPOLOGY";
    std::string marker = std::string("--") + sectionname;

    auto position = positions_.find(marker);
    if (position != positions_.end())
    {
      for (size_t pos = position->second + 1; pos < lines_.size(); ++pos)
      {
        const char* l = lines_[pos];
        if (l[0] == '-' and l[1] == '-')
        {
          break;
        }

        int dobj;
        int nodeid;
        std::string nname;
        std::string dname;
        std::string disname;
        std::array<int, 3> dir = {0, 0, 0};

        std::istringstream stream(l);
        stream >> nname;
        if (not stream) dserror("Illegal line in section '%s': '%s'", marker.c_str(), l);

        if (nname == "NODE")  // plain old reading of the design nodes from the .dat-file
        {
          stream >> nodeid >> dname >> dobj;
          topology[dobj - 1].insert(nodeid - 1);
        }
        else  // fancy specification of the design nodes by specifying min or max of the domain
        {     // works best on rectangular domains ;)
          if (nname == "CORNER" && name == "DNODE")
          {
            std::string tmp;
            stream >> disname;
            for (int i = 0; i < 3; ++i)
            {
              stream >> tmp;
              if (tmp.size() != 2 || tmp[0] < 'x' || tmp[0] > 'z' ||
                  (tmp[1] != '+' && tmp[1] != '-'))
                dserror("Illegal design node definition.");
              dir[tmp[0] - 'x'] = (tmp[1] == '+') ? 1 : -1;
            }
            stream >> dname >> dobj;
          }
          else if (nname == "EDGE" && name == "DLINE")
          {
            std::string tmp;
            stream >> disname;
            for (int i = 0; i < 2; ++i)
            {
              stream >> tmp;
              if (tmp.size() != 2 || tmp[0] < 'x' || tmp[0] > 'z' ||
                  (tmp[1] != '+' && tmp[1] != '-'))
                dserror("Illegal design node definition.");
              dir[tmp[0] - 'x'] = (tmp[1] == '+') ? 1 : -1;
            }
            stream >> dname >> dobj;
          }
          else if (nname == "SIDE" && name == "DSURF")
          {
            std::string tmp;
            stream >> disname;
            stream >> tmp;
            if (tmp.size() != 2 || tmp[0] < 'x' || tmp[0] > 'z' || (tmp[1] != '+' && tmp[1] != '-'))
              dserror("Illegal design node definition.");
            dir[tmp[0] - 'x'] = (tmp[1] == '+') ? 1 : -1;
            stream >> dname >> dobj;
          }
          else if (nname == "VOLUME" && name == "DVOL")
          {
            stream >> disname;
            stream >> dname >> dobj;
          }
          else
          {
            dserror("Illegal line in section '%s': '%s'", marker.c_str(), l);
          }

          Teuchos::RCP<DRT::Discretization> actdis = GLOBAL::Problem::Instance()->GetDis(disname);

          std::vector<double> box_specifications;
          if (cached_box_specifications_.find(disname) != cached_box_specifications_.end())
          {
            box_specifications = cached_box_specifications_[disname];
          }
          else
          {
            for (int init = 0; init < 9; ++init) box_specifications.push_back(0.0);
            if (comm_->MyPID() == 0)  // Reading is done by proc 0
            {
              // get original domain section from the *.dat-file
              std::string dommarker = "--" + disname + " DOMAIN";
              std::transform(dommarker.begin(), dommarker.end(), dommarker.begin(), ::toupper);
              auto di = excludepositions_.find(dommarker);
              if (di != excludepositions_.end())
              {
                std::ifstream tmpfile(filename_.c_str());
                tmpfile.seekg(di->second.first);
                std::string line;
                for (int ii = 0; getline(tmpfile, line); ++ii)
                {
                  // remove comments, trailing and leading whitespaces
                  // compact internal whitespaces
                  line = DRT::UTILS::StripComment(line);

                  // line is now empty
                  if (line.size() == 0) continue;

                  if (line.find("--") == 0)
                  {
                    break;
                  }
                  else
                  {
                    std::istringstream t;
                    t.str(line);
                    std::string key;
                    t >> key;

                    if (key == "LOWER_BOUND")
                    {
                      t >> box_specifications[0] >> box_specifications[1] >> box_specifications[2];
                    }
                    else if (key == "UPPER_BOUND")
                    {
                      t >> box_specifications[3] >> box_specifications[4] >> box_specifications[5];
                    }
                    else if (key == "ROTATION")
                    {
                      t >> box_specifications[6] >> box_specifications[7] >> box_specifications[8];
                    }
                  }
                }
              }
              else
                dserror("Inputreader: Couldn't find domain section for discretization %s !",
                    disname.c_str());
            }
            // All other processors get this info broadcasted
            comm_->Broadcast(
                box_specifications.data(), static_cast<int>(box_specifications.size()), 0);
            cached_box_specifications_[disname] = box_specifications;
          }

          // determine the active discretizations bounding box
          std::array<double, 6> bbox;
          for (size_t i = 0; i < sizeof(bbox) / sizeof(bbox[0]); ++i)
            bbox[i] = box_specifications[i];

          // manipulate the bounding box according to the specified condition
          for (size_t i = 0; i < 3; ++i)
          {
            switch (dir[i])
            {
              case 0:
                bbox[i + 0] = std::numeric_limits<double>::max();
                bbox[i + 3] = -std::numeric_limits<double>::max();
                break;
              case -1:
                bbox[i] += tolerance_n;
                bbox[i + 3] = std::numeric_limits<double>::max();
                break;
              case 1:
                bbox[i] = -std::numeric_limits<double>::max();
                bbox[i + 3] -= tolerance_n;
                break;
              default:
                dserror("Invalid BC specification");
            }
          }

          // collect all nodes which are outside the adapted bounding box
          std::set<int> dnodes;
          for (int lid = 0; lid < actdis->NumMyRowNodes(); ++lid)
          {
            const DRT::Node* node = actdis->lRowNode(lid);
            const auto& coord = node->X();
            std::array<double, 3> coords;
            coords[0] = coord[0];
            coords[1] = coord[1];
            coords[2] = coord[2];
            // rotate back to identify condition, if a rotation is defined
            static const int rotoffset = 6;
            for (int rotaxis = 2; rotaxis > -1; --rotaxis)
            {
              if (box_specifications[rotaxis + rotoffset] != 0.0)
              {
                std::array<double, 3> coordm;
                coordm[0] = (box_specifications[0] + box_specifications[3]) / 2.;
                coordm[1] = (box_specifications[1] + box_specifications[4]) / 2.;
                coordm[2] = (box_specifications[2] + box_specifications[5]) / 2.;
                // add rotation around mitpoint here.
                std::array<double, 3> dx;
                dx[0] = coords[0] - coordm[0];
                dx[1] = coords[1] - coordm[1];
                dx[2] = coords[2] - coordm[2];

                double calpha = cos(-box_specifications[rotaxis + rotoffset] * M_PI / 180);
                double salpha = sin(-box_specifications[rotaxis + rotoffset] * M_PI / 180);

                coords[0] = coordm[0];  //+ calpha*dx[0] + salpha*dx[1];
                coords[1] = coordm[1];  //+ -salpha*dx[0] + calpha*dx[1];
                coords[2] = coordm[2];

                coords[(rotaxis + 1) % 3] +=
                    calpha * dx[(rotaxis + 1) % 3] + salpha * dx[(rotaxis + 2) % 3];
                coords[(rotaxis + 2) % 3] +=
                    calpha * dx[(rotaxis + 2) % 3] - salpha * dx[(rotaxis + 1) % 3];
                coords[rotaxis] += dx[rotaxis];
              }
            }

            if ((coords[0] <= bbox[0] || coords[0] >= bbox[3]) &&
                (coords[1] <= bbox[1] || coords[1] >= bbox[4]) &&
                (coords[2] <= bbox[2] || coords[2] >= bbox[5]))
              dnodes.insert(node->Id());
          }
          CORE::LINALG::GatherAll(dnodes, *comm_);
          topology[dobj - 1].insert(dnodes.begin(), dnodes.end());
        }

        if (dname.substr(0, name.length()) != name)
          dserror("Illegal line in section '%s': '%s'\n%s found, where %s was expected",
              marker.c_str(), l, dname.substr(0, name.length()).c_str(), name.c_str());
      }
      if (topology.size() > 0)
      {
        int max_num_dobj = topology.rbegin()->first;
        if (max_num_dobj >= static_cast<int>(dobj_fenode.size()))
          dobj_fenode.resize(max_num_dobj + 1);
        // copy all design object entries
        for (auto& topo : topology)
        {
          // we copy from a std::set, thus the gids are sorted
          dobj_fenode[topo.first].reserve(topo.second.size());
          dobj_fenode[topo.first].assign(topo.second.begin(), topo.second.end());
        }
      }
    }

    // The section name is desired from outside. Thus, we consider it as valid
    knownsections_[marker] = true;
  }


  //----------------------------------------------------------------------
  /// read a knotvector section (for isogeometric analysis)
  //----------------------------------------------------------------------
  void DatFileReader::ReadKnots(
      const int dim, const std::string& name, Teuchos::RCP<DRT::NURBS::Knotvector>& disknots)
  {
    // io to shell
    const int myrank = comm_->MyPID();

    Teuchos::Time time("", true);

    // only the knotvector section of this discretisation
    // type is of interest
    std::string field;
    if (name == "fluid" or name == "xfluid" or name == "porofluid")
    {
      field = "FLUID";
    }
    else if (name == "structure")
    {
      field = "STRUCTURE";
    }
    else if (name == "ale")
    {
      field = "ALE";
    }
    else if (name == "scatra")
    {
      field = "TRANSPORT";
    }
    else if (name == "thermo")
    {
      field = "THERMO";
    }
    else if (name == "scatra_micro")
    {
      field = "TRANSPORT2";
    }
    else
    {
      dserror("Unknown discretization name for knotvector input\n");
    }

    // another valid section name was found
    const std::string sectionname = "--" + field + " KNOTVECTORS";
    knownsections_[sectionname] = true;

    if (myrank == 0)
    {
      if (!MyOutputFlag())
      {
        IO::cout << "Reading knot vectors for " << name << " discretization :\n";
        fflush(stdout);
      }
    }

    // number of patches to be determined
    int npatches = 0;

    //--------------------------------------------------------------------
    //--------------------------------------------------------------------
    //             first, determine number of patches
    //--------------------------------------------------------------------
    //--------------------------------------------------------------------
    {
      // open input file --- this is done on all procs
      std::ifstream file;
      file.open(filename_.c_str());

      // temporary string
      std::string tmp;

      // flag indicating knot vector section in input
      bool knotvectorsection = false;

      // loop lines in file
      for (; file;)
      {
        // read piece of file until next seperator (whitespace, newline)
        file >> tmp;

        // if this a new section, i.e. starts like ------
        if ((tmp[0] == '-' && tmp[1] == '-'))
        {
          // check whether it is the knotvectorsection
          // only the knotvector section of this discretisation type is of interest
          std::string::size_type loc = tmp.rfind(field);

          if (loc == std::string::npos)
          {
            knotvectorsection = false;

            // there is nothing more to be done in this line
            continue;
          }
          else
          {
            // continue reading of second keyword
            file >> tmp;

            // check whether second keyword is knotvector
            loc = tmp.rfind("KNOTVECTORS");

            if (loc != std::string::npos)
            {
              // if this is true, we are at the beginning of a
              // knot section
              knotvectorsection = true;
              // there is nothing more to be done in this line
              continue;
            }
            else
            {
              knotvectorsection = false;

              // there is nothing more to be done in this line
              continue;
            }
          }
        }

        // count number of patches in knotvector section of
        // this discretisation
        if (knotvectorsection)
        {
          // check for a new patch
          std::string::size_type loc;

          loc = tmp.rfind("ID");
          if (loc != std::string::npos)
          {
            // increase number of patches
            npatches++;

            continue;
          }
        }
      }  // end loop through file
    }

    if (myrank == 0)
    {
      if (!MyOutputFlag())
      {
        printf("                        %8d patches", npatches);
        fflush(stdout);
      }
    }


    //--------------------------------------------------------------------
    //--------------------------------------------------------------------
    //                alloc knotvector object to fill
    //--------------------------------------------------------------------
    //--------------------------------------------------------------------

    // allocate knotvector for this dis
    disknots = Teuchos::rcp(new DRT::NURBS::Knotvector(dim, npatches));

    // make sure that we have some Knotvector object to fill
    if (disknots == Teuchos::null)
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
      std::vector<Teuchos::RCP<std::vector<double>>> patch_knots(dim);

      // open input file --- this is done on all procs
      std::ifstream file;

      file.open(filename_.c_str());

      // temporary string
      std::string tmp;

      // start to read something when read is true
      bool read = false;

      // flag indicating knot vector section in input
      bool knotvectorsection = false;

      // index for number of patch
      int npatch = 0;
      // index for u/v/w
      int actdim = -1;
      // ints for the number of knots
      std::vector<int> n_x_m_x_l(dim);
      // ints for patches degrees
      std::vector<int> degree(dim);
      // a vector of strings holding the knotvectortypes read
      std::vector<std::string> knotvectortype(dim);

      // count for sanity check
      int count_read = 0;
      std::vector<int> count_vals(dim);

      // loop lines in file
      for (; file;)
      {
        file >> tmp;

        // if this a new section
        if ((tmp[0] == '-' && tmp[1] == '-'))
        {
          // check whether it is the knotvectorsection
          // only the knotvector section of this discretisation
          // type is of interest
          std::string::size_type loc = tmp.rfind(field);

          if (loc == std::string::npos)
          {
            knotvectorsection = false;

            // there is nothing more to be done in this line
            continue;
          }
          else
          {
            // continue reading of second keyword
            file >> tmp;

            // check whether second keyword is knotvector
            loc = tmp.rfind("KNOTVECTORS");

            if (loc != std::string::npos)
            {
              // if this is true, we are at the beginning of a
              // knot section
              knotvectorsection = true;
              // there is nothing more to be done in this line
              continue;
            }
            else
            {
              knotvectorsection = false;

              // there is nothing more to be done in this line
              continue;
            }
          }
        }

        // do reading in knotvecor section
        if (knotvectorsection)
        {
          // check for a new patch
          std::string::size_type loc = tmp.rfind("BEGIN");
          if (loc != std::string::npos)
          {
            file >> tmp;

            // activate reading
            read = true;

            actdim = -1;

            // create vectors for knots in this patch
            for (int rr = 0; rr < dim; ++rr)
            {
              patch_knots[rr] = Teuchos::rcp(new std::vector<double>);
              (*(patch_knots[rr])).clear();
            }

            // reset counter for knot values
            for (int rr = 0; rr < dim; rr++)
            {
              count_vals[rr] = 0;
            }

            continue;
          }

          // get ID of patch we are currently reading
          loc = tmp.rfind("ID");
          if (loc != std::string::npos)
          {
            std::string str_npatch;
            file >> str_npatch;

            char* endptr = nullptr;
            npatch = static_cast<int>(strtol(str_npatch.c_str(), &endptr, 10));
            npatch--;

            continue;
          }

          // get number of knots in the knotvector direction
          // we are currently reading
          loc = tmp.rfind("NUMKNOTS");
          if (loc != std::string::npos)
          {
            std::string str_numknots;
            file >> str_numknots;

            // increase dimesion for knotvector (i.e. next time
            // we'll fill the following knot vector)
            actdim++;
            if (actdim > dim)
            {
              dserror("too many knotvectors (we only need dim)\n");
            }

            char* endptr = nullptr;
            n_x_m_x_l[actdim] = static_cast<int>(strtol(str_numknots.c_str(), &endptr, 10));

            continue;
          }

          // get number of bspline polinomial associated with
          // knots in this direction
          loc = tmp.rfind("DEGREE");
          if (loc != std::string::npos)
          {
            std::string str_degree;
            file >> str_degree;

            char* endptr = nullptr;
            degree[actdim] = static_cast<int>(strtol(str_degree.c_str(), &endptr, 10));

            continue;
          }

          // get type of knotvector (interpolated or periodic)
          loc = tmp.rfind("TYPE");
          if (loc != std::string::npos)
          {
            std::string type;

            file >> type;
            knotvectortype[actdim] = type;

            continue;
          }

          // locate end of patch
          loc = tmp.rfind("END");
          if (loc != std::string::npos)
          {
            for (int rr = 0; rr < dim; ++rr)
            {
              disknots->SetKnots(
                  rr, npatch, degree[rr], n_x_m_x_l[rr], knotvectortype[rr], patch_knots[rr]);
            }
            file >> tmp;
            // stop reading of knot values if we are here
            read = false;

            for (int rr = 0; rr < dim; rr++)
            {
              if (n_x_m_x_l[rr] != count_vals[rr])
              {
                dserror("not enough knots read in dim %d (%d!=NUMKNOTS=%d)\n", rr, count_vals[rr],
                    n_x_m_x_l[rr]);
              }
            }

            // count for sanity check
            count_read++;

            continue;
          }

          //  reading of knot values if read is true and no
          // other keyword was found
          if (read)
          {
            char* endptr = nullptr;

            double dv = strtod(tmp.c_str(), &endptr);

            // count for sanity check
            count_vals[actdim]++;

            (*(patch_knots[actdim])).push_back(dv);
          }
        }
      }  // end loop through file

      if (count_read != npatches)
      {
        dserror("wasn't able to read enough patches\n");
      }
    }

    if (myrank == 0)
    {
      if (!MyOutputFlag())
      {
        IO::cout << " in...." << time.totalElapsedTime(true) << " secs\n";

        time.reset();
        fflush(stdout);
      }
    }
  }


  /*----------------------------------------------------------------------*/
  /*----------------------------------------------------------------------*/
  void DatFileReader::AddEntry(
      const std::string& key, const std::string& value, Teuchos::ParameterList& list)
  {
    // safety check: Is there a duplicate of the same parameter?
    if (list.isParameter(key))
      dserror("Duplicate parameter %s in sublist %s", key.c_str(), list.name().c_str());

    if (key.empty()) dserror("Internal error: missing key.", key.c_str());
    // safety check: Is the parameter without any specified value?
    if (value.empty()) dserror("Missing value for parameter %s. Fix your input file!", key.c_str());

    {  // try to find an int
      std::stringstream ssi;
      int iv;

      ssi << value;
      ssi >> iv;

      if (ssi.eof())
      {
        list.set(key, iv);
        return;
      }
    }

#ifdef BACI_TRAP_FE
    // somehow the following test whether we have a double or not
    // creates always an internal floating point exception (FE_INVALID). An alternative
    // implementation using boost::lexical_cast<double> does not solve this problem!
    // Better temporarily disable this floating point exception in the following,
    // so that we can go on.
    feclearexcept(FE_INVALID);
    /*feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_UNDERFLOW | FE_OVERFLOW);*/
    fedisableexcept(FE_INVALID);
#endif

    {  // try to find a double
      std::stringstream ssd;
      double dv;

      ssd << value;
      ssd >> dv;

#ifdef BACI_TRAP_FE
      feclearexcept(FE_INVALID);
      /*feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_UNDERFLOW | FE_OVERFLOW);*/
      feenableexcept(FE_INVALID | FE_DIVBYZERO);
#endif

      if (ssd.eof())
      {
        list.set(key, dv);
        return;
      }
    }

    // if it is not an int or a double it must be a string
    list.set(key, value);
  }


  /*----------------------------------------------------------------------*/
  /*----------------------------------------------------------------------*/
  void DatFileReader::ReadDat()
  {
    std::vector<std::string> exclude;

    exclude.emplace_back("--NODE COORDS");
    exclude.emplace_back("--STRUCTURE ELEMENTS");
    exclude.emplace_back("--STRUCTURE DOMAIN");
    exclude.emplace_back("--FLUID ELEMENTS");
    exclude.emplace_back("--FLUID DOMAIN");
    exclude.emplace_back("--ALE ELEMENTS");
    exclude.emplace_back("--ALE DOMAIN");
    exclude.emplace_back("--ARTERY ELEMENTS");
    exclude.emplace_back("--REDUCED D AIRWAYS ELEMENTS");
    exclude.emplace_back("--LUBRICATION ELEMENTS");
    exclude.emplace_back("--LUBRICATION DOMAIN");
    exclude.emplace_back("--TRANSPORT ELEMENTS");
    exclude.emplace_back("--TRANSPORT2 ELEMENTS");
    exclude.emplace_back("--TRANSPORT DOMAIN");
    exclude.emplace_back("--THERMO ELEMENTS");
    exclude.emplace_back("--THERMO DOMAIN");
    exclude.emplace_back("--ELECTROMAGNETIC ELEMENTS");
    exclude.emplace_back("--PERIODIC BOUNDINGBOX ELEMENTS");
    exclude.emplace_back("--CELL ELEMENTS");
    exclude.emplace_back("--CELL DOMAIN");
    exclude.emplace_back("--CELLSCATRA ELEMENTS");
    exclude.emplace_back("--CELLSCATRA DOMAIN");
    exclude.emplace_back("--PARTICLES");

    Teuchos::RCP<Epetra_Comm> comm = comm_;

    // if we are in copydatafile mode use global comm instead of local comm
    // and only read in input file once instead of npgroup times
    GLOBAL::Problem* problem = GLOBAL::Problem::Instance();
    CORE::COMM::NestedParallelismType npType = problem->GetCommunicators()->NpType();
    if (npType == CORE::COMM::NestedParallelismType::copy_dat_file)
      comm = problem->GetCommunicators()->GlobalComm();

    int arraysize = 0;
    if (comm->MyPID() == 0)
    {
      std::ifstream file(filename_.c_str());
      if (not file) dserror("unable to open file: %s", filename_.c_str());

      std::list<std::string> content;
      bool ignoreline = false;
      std::string line;
      unsigned int* linecount = nullptr;

      // loop all input lines
      while (getline(file, line))
      {
        // remove comments, trailing and leading whitespaces
        // compact internal whitespaces
        line = DRT::UTILS::StripComment(line);

        // line is now empty
        if (line.size() == 0) continue;

        // exclude all special sections
        // this includes the section header and all lines in that section
        if (ignoreline)
        {
          if (line.find("--") == 0)
          {
            ignoreline = false;
            linecount = nullptr;
          }
          else
          {
            *linecount += 1;
          }
        }

        // Two sections to be ignored can follow each other. We need
        // independent tests.
        if (!ignoreline)
        {
          // remember all section positions
          if (line.find("--") == 0)
          {
            for (auto& i : exclude)
            {
              if (line.find(i) != std::string::npos)
              {
                if (excludepositions_.find(i) != excludepositions_.end())
                  dserror("section '%s' defined more than once", i.c_str());
                std::pair<std::ifstream::pos_type, unsigned int>& p = excludepositions_[i];
                p.first = file.tellg();
                p.second = 0;
                linecount = &p.second;
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
          arraysize += static_cast<int>(line.length()) + 1;
        }
      }

      // setup global variables
      numrows_ = static_cast<int>(content.size());

      // allocate space for copy of file
      inputfile_.clear();
      inputfile_.reserve(arraysize);

      // CAUTION: We allocate one more row pointer that necessary. This
      // pointer will be used to point to temporary lines when the
      // excluded section are read (on proc 0). Of course that's just
      // another EVIL HACK. Don't tell anybody.
      lines_.reserve(numrows_ + 1);

      lines_.push_back(inputfile_.data());
      for (auto& i : content)
      {
        inputfile_.insert(inputfile_.end(), i.begin(), i.end());
        inputfile_.push_back('\0');
        lines_.push_back(&(inputfile_.back()) + 1);
      }
      // add the slot for the temporary line...
      lines_.back() = nullptr;

      if (inputfile_.size() != static_cast<size_t>(arraysize))
      {
        dserror(
            "internal error in file read: inputfile has %d chars, but was predicted to be %d "
            "chars long",
            inputfile_.size(), arraysize);
      }
    }

    // Now lets do all the parallel setup. Afterwards all processors
    // have to be the same.

    if (comm->NumProc() > 1)
    {
      /* Now that we use a variable number of bytes per line we have to
       * communicate the buffer size as well. */
      comm->Broadcast(&arraysize, 1, 0);
      comm->Broadcast(&numrows_, 1, 0);

      if (comm->MyPID() > 0)
      {
        /*--------------------------------------allocate space for copy of file */
        inputfile_.resize(arraysize);
        lines_.reserve(numrows_ + 1);
      }

      // There are no char based functions available! Do it by hand!
      // comm->Broadcast(inputfile_.data(),arraysize,0);

      const auto& mpicomm = dynamic_cast<const Epetra_MpiComm&>(*comm);

      MPI_Bcast(inputfile_.data(), arraysize, MPI_CHAR, 0, mpicomm.GetMpiComm());

      /* We have not yet set the row pointers on procs > 0. So do it now. */
      if (comm->MyPID() > 0)
      {
        lines_.push_back(inputfile_.data());
        for (int i = 0; i < arraysize; ++i)
        {
          if (inputfile_[i] == '\0')
          {
            lines_.push_back(&inputfile_[i + 1]);
          }
        }

        if (static_cast<int>(lines_.size()) != numrows_ + 1)
          dserror("line count mismatch: %d lines expected but %d lines received", numrows_ + 1,
              lines_.size());
      }

      // distribute excluded section positions
      for (auto& i : exclude)
      {
        if (comm->MyPID() == 0)
        {
          auto ep = excludepositions_.find(i);
          if (ep == excludepositions_.end())
          {
            excludepositions_[i] = std::pair<std::ifstream::pos_type, unsigned int>(-1, 0);
          }
        }
        std::pair<std::ifstream::pos_type, unsigned int>& p = excludepositions_[i];
        // comm->Broadcast(&p.second,1,0);
        MPI_Bcast(&p.second, 1, MPI_INT, 0, mpicomm.GetMpiComm());
      }
    }

    // Now finally find the section names. We have to do this on all
    // processors, so it cannot be done while reading.
    for (std::vector<char*>::size_type i = 0; i < lines_.size() - 1; ++i)
    {
      char* l = lines_[i];
      if (l and l[0] == '-' and l[1] == '-')
      {
        std::string line(l);

        // take the last "--" and all that follows as section name
        std::string::size_type loc = line.rfind("--");
        std::string sectionname = line.substr(loc);
        if (positions_.find(sectionname) != positions_.end())
          dserror("section '%s' defined more than once", sectionname.c_str());
        positions_[sectionname] = i;
        knownsections_[sectionname] = false;  // initialize as unknown
      }
    }

    if (positions_.find("--END") == positions_.end())
      dserror("end section missing. incomplete dat file?");

    // the following section names are always regarded as valid
    knownsections_["--END"] = true;
    knownsections_["--TITLE"] = true;
    knownsections_["--FUNCT1"] = true;
    knownsections_["--FUNCT2"] = true;
    knownsections_["--FUNCT3"] = true;
    knownsections_["--FUNCT4"] = true;
    knownsections_["--FUNCT5"] = true;
    knownsections_["--FUNCT6"] = true;
    knownsections_["--FUNCT7"] = true;
    knownsections_["--FUNCT8"] = true;
    knownsections_["--FUNCT9"] = true;
    knownsections_["--FUNCT10"] = true;
    knownsections_["--FUNCT11"] = true;
    knownsections_["--FUNCT12"] = true;
    knownsections_["--FUNCT13"] = true;
    knownsections_["--FUNCT14"] = true;
    knownsections_["--FUNCT15"] = true;
    knownsections_["--FUNCT16"] = true;
    knownsections_["--FUNCT17"] = true;
    knownsections_["--FUNCT18"] = true;
    knownsections_["--FUNCT19"] = true;
    knownsections_["--FUNCT20"] = true;
  }


  /*----------------------------------------------------------------------*/
  /*----------------------------------------------------------------------*/
  bool DatFileReader::PrintUnknownSections()
  {
    // This function shell be called only after all reading with DatFileReader
    // is finished. Only then we have a proper protocol of (in-)valid section names

    const bool printout = std::any_of(
        knownsections_.begin(), knownsections_.end(), [](const auto& kv) { return !kv.second; });

    // now it's time to create noise on the screen
    if (printout and (Comm()->MyPID() == 0))
    {
      IO::cout << "\nERROR!"
               << "\n--------"
               << "\nThe following input file sections remained unused (obsolete or typo?):"
               << IO::endl;
      for (const auto& [section_name, known] : knownsections_)
      {
        if (!known) IO::cout << section_name << IO::endl;
      }
      IO::cout << IO::endl;
    }

    // we wait till all procs are here. Otherwise a hang up might occur where
    // one proc ended with dserror but other procs were not finished and waited...
    // we also want to have the printing above being finished.
    Comm()->Barrier();
    if (printout)
      dserror(
          "Unknown sections detected. Correct this! Find hints on these unknown sections above.");

    return printout;
  }

  std::pair<std::string, std::string> ReadKeyValue(const std::string& line)
  {
    std::string::size_type separator_index = line.find('=');
    // The equals sign is only treated as a separator when surrounded by whitespace.
    if (separator_index != std::string::npos &&
        !(std::isspace(line[separator_index - 1]) && std::isspace(line[separator_index + 1])))
      separator_index = std::string::npos;

    // In case we didn't find an "=" separator, look for a space instead
    if (separator_index == std::string::npos)
    {
      separator_index = line.find(' ');

      if (separator_index == std::string::npos)
        dserror("Line '%s' with just one word in parameter section", line.c_str());
    }

    std::string key = DRT::UTILS::Trim(line.substr(0, separator_index));
    std::string value = DRT::UTILS::Trim(line.substr(separator_index + 1));

    if (key.empty()) dserror("Cannot get key from line '%s'", line.c_str());
    if (value.empty()) dserror("Cannot get value from line '%s'", line.c_str());

    return {std::move(key), std::move(value)};
  }
}  // namespace INPUT

BACI_NAMESPACE_CLOSE
