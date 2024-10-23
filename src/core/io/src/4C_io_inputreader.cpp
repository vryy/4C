// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_io_inputreader.hpp"

#include "4C_fem_general_node.hpp"
#include "4C_fem_nurbs_discretization_knotvector.hpp"
#include "4C_io_pstream.hpp"
#include "4C_linalg_utils_densematrix_communication.hpp"
#include "4C_utils_string.hpp"

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_Time.hpp>

#include <sstream>
#include <utility>

#ifdef FOUR_C_ENABLE_FE_TRAPPING
#include <cfenv>
#endif


FOUR_C_NAMESPACE_OPEN

namespace
{
  constexpr double tolerance_n = 1.0e-14;

  /*----------------------------------------------------------------------*/
  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& find_sublist(std::string name, Teuchos::ParameterList& list)
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

  void add_entry(const std::string& key, const std::string& value, Teuchos::ParameterList& list)
  {
    // safety check: Is there a duplicate of the same parameter?
    if (list.isParameter(key))
      FOUR_C_THROW("Duplicate parameter %s in sublist %s", key.c_str(), list.name().c_str());

    if (key.empty()) FOUR_C_THROW("Internal error: missing key.", key.c_str());
    // safety check: Is the parameter without any specified value?
    if (value.empty())
      FOUR_C_THROW("Missing value for parameter %s. Fix your input file!", key.c_str());

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

#ifdef FOUR_C_ENABLE_FE_TRAPPING
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

#ifdef FOUR_C_ENABLE_FE_TRAPPING
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

}  // namespace

namespace Core::IO
{
  namespace Internal
  {

    StreamLineIterator::StreamLineIterator(std::shared_ptr<std::istream> stream)
        : StreamLineIterator(std::move(stream), std::numeric_limits<int>::max())
    {
    }

    StreamLineIterator::StreamLineIterator(std::shared_ptr<std::istream> stream, int max_reads)
        : stream_(std::move(stream)), max_reads_(max_reads)
    {
      ++(*this);
    }


    StreamLineIterator::StreamLineIterator() : line_number_(-1) {}


    StreamLineIterator& StreamLineIterator::operator++()
    {
      if (stream_)
      {
        if (line_number_ < max_reads_ && std::getline(*stream_, line_))
        {
          line_number_++;
        }
        else
        {
          // we hit EOF, set special value
          line_number_ = -1;
        }
      }

      return *this;
    }


    StreamLineIterator::reference StreamLineIterator::operator*() const { return line_; }
    bool StreamLineIterator::operator==(const StreamLineIterator& other) const
    {
      return line_number_ == other.line_number_;
    }


    bool StreamLineIterator::operator!=(const StreamLineIterator& other) const
    {
      return !(*this == other);
    }


    DatFileLineIterator::DatFileLineIterator(
        std::variant<StreamLineIterator, PreReadIterator> iterator)
        : iterator_(std::move(iterator))
    {
    }


    DatFileLineIterator& DatFileLineIterator::operator++()
    {
      std::visit([](auto& it) { ++it; }, iterator_);
      return *this;
    }


    DatFileLineIterator::reference DatFileLineIterator::operator*() const
    {
      return std::visit([](const auto& it) -> reference { return *it; }, iterator_);
    }


    bool DatFileLineIterator::operator==(const DatFileLineIterator& other) const
    {
      if (other.iterator_.index() != iterator_.index()) return false;

      if (iterator_.index() == 0)
      {
        return std::get<0>(iterator_) == std::get<0>(other.iterator_);
      }
      else
      {
        return std::get<1>(iterator_) == std::get<1>(other.iterator_);
      }
    }


    bool DatFileLineIterator::operator!=(const DatFileLineIterator& other) const
    {
      return !(*this == other);
    }
  }  // namespace Internal

  /*----------------------------------------------------------------------*/
  /*----------------------------------------------------------------------*/
  DatFileReader::DatFileReader(std::string filename, const Epetra_Comm& comm, int outflag)
      : filename_(std::move(filename)), comm_(std::move(comm)), outflag_(outflag)
  {
    read_dat();
  }


  /*----------------------------------------------------------------------*/
  /*----------------------------------------------------------------------*/
  std::string DatFileReader::my_inputfile_name() const { return filename_; }


  /*----------------------------------------------------------------------*/
  /*----------------------------------------------------------------------*/
  int DatFileReader::my_output_flag() const { return outflag_; }


  /*----------------------------------------------------------------------*/
  /*----------------------------------------------------------------------*/
  bool DatFileReader::has_section(const std::string& section_name) const
  {
    const auto range = line_range(section_name);
    return range.begin() != range.end();
  }


  /*----------------------------------------------------------------------*/
  /*----------------------------------------------------------------------*/
  bool read_parameters_in_section(
      DatFileReader& reader, const std::string& section_name, Teuchos::ParameterList& list)
  {
    if (section_name.length() < 3 or section_name[0] != '-' or section_name[1] != '-')
      FOUR_C_THROW("Illegal section name '%s'", section_name.c_str());

    Teuchos::ParameterList& sublist = find_sublist(section_name.substr(2), list);

    for (const auto& line : reader.lines_in_section(section_name))
    {
      // If the line starts with dashes we found the beginning of the next section and terminate
      // the read process of the current section.
      if (line[0] == '-' and line[1] == '-')
      {
        break;
      }

      const auto& [key, value] = read_key_value(std::string(line));

      add_entry(key, value, sublist);
    }

    return true;
  }

  /*----------------------------------------------------------------------*/
  /*----------------------------------------------------------------------*/
  void read_design(DatFileReader& reader, const std::string& name,
      std::vector<std::vector<int>>& dobj_fenode,
      const std::function<const Core::FE::Discretization&(const std::string& name)>&
          get_discretization)
  {
    std::map<int, std::set<int>> topology;

    std::string sectionname = name + "-NODE TOPOLOGY";
    std::string marker = std::string("--") + sectionname;

    for (const auto& l : reader.lines_in_section(marker))
    {
      int dobj;
      int nodeid;
      std::string nname;
      std::string dname;
      std::string disname;
      std::array<int, 3> dir = {0, 0, 0};

      std::istringstream stream{std::string(l)};
      stream >> nname;
      if (not stream) FOUR_C_THROW("Illegal line in section '%s': '%s'", marker.c_str(), l.data());

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
            if (tmp.size() != 2 || tmp[0] < 'x' || tmp[0] > 'z' || (tmp[1] != '+' && tmp[1] != '-'))
              FOUR_C_THROW("Illegal design node definition.");
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
            if (tmp.size() != 2 || tmp[0] < 'x' || tmp[0] > 'z' || (tmp[1] != '+' && tmp[1] != '-'))
              FOUR_C_THROW("Illegal design node definition.");
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
            FOUR_C_THROW("Illegal design node definition.");
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
          FOUR_C_THROW("Illegal line in section '%s': '%s'", marker.c_str(), l.data());
        }

        const Core::FE::Discretization& actdis = get_discretization(disname);

        std::vector<double> box_specifications;
        {
          for (int init = 0; init < 9; ++init) box_specifications.push_back(0.0);
          if (reader.get_comm().MyPID() == 0)  // Reading is done by proc 0
          {
            // get original domain section from the *.dat-file
            std::string dommarker = "--" + disname + " DOMAIN";
            std::transform(dommarker.begin(), dommarker.end(), dommarker.begin(), ::toupper);

            FOUR_C_THROW_UNLESS(reader.has_section(dommarker),
                "Inputreader: Couldn't find domain section for discretization %s !",
                disname.c_str());

            for (const auto& line : reader.lines_in_section(dommarker))
            {
              std::istringstream t{std::string{line}};
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
          // All other processors get this info broadcasted
          reader.get_comm().Broadcast(
              box_specifications.data(), static_cast<int>(box_specifications.size()), 0);
        }

        // determine the active discretizations bounding box
        std::array<double, 6> bbox;
        for (size_t i = 0; i < sizeof(bbox) / sizeof(bbox[0]); ++i) bbox[i] = box_specifications[i];

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
              FOUR_C_THROW("Invalid BC specification");
          }
        }

        // collect all nodes which are outside the adapted bounding box
        std::set<int> dnodes;
        for (const auto* node : actdis.my_row_node_range())
        {
          const auto& coord = node->x();
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
            dnodes.insert(node->id());
        }
        Core::LinAlg::gather_all(dnodes, reader.get_comm());
        topology[dobj - 1].insert(dnodes.begin(), dnodes.end());
      }

      if (dname.substr(0, name.length()) != name)
        FOUR_C_THROW("Illegal line in section '%s': '%s'\n%s found, where %s was expected",
            marker.c_str(), l.data(), dname.substr(0, name.length()).c_str(), name.c_str());
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


  //----------------------------------------------------------------------
  /// read a knotvector section (for isogeometric analysis)
  //----------------------------------------------------------------------
  void read_knots(DatFileReader& reader, const std::string& name,
      Teuchos::RCP<Core::FE::Nurbs::Knotvector>& disknots)
  {
    // io to shell
    const int myrank = reader.get_comm().MyPID();

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
      FOUR_C_THROW("Unknown discretization name for knotvector input\n");
    }

    // another valid section name was found
    const std::string sectionname = "--" + field + " KNOTVECTORS";

    if (myrank == 0)
    {
      if (!reader.my_output_flag())
      {
        Core::IO::cout << "Reading knot vectors for " << name << " discretization :\n";
        fflush(stdout);
      }
    }

    // number of patches to be determined
    int npatches = 0;

    // dimension of nurbs patches
    int nurbs_dim = 0;

    //--------------------------------------------------------------------
    //--------------------------------------------------------------------
    //      first, determine number of patches and dimension of nurbs
    //--------------------------------------------------------------------
    //--------------------------------------------------------------------
    {
      // temporary string
      std::string tmp;
      // loop lines in file
      for (const auto& line : reader.lines_in_section(sectionname))
      {
        // count number of patches in knotvector section of
        // this discretisation
        {
          std::string::size_type loc;
          std::istringstream file{std::string{line}};
          file >> tmp;

          // check for the number of dimensions
          loc = tmp.rfind("NURBS_DIMENSION");
          if (loc != std::string::npos)
          {
            // set number of nurbs dimension
            std::string str_nurbs_dim;
            file >> str_nurbs_dim;
            char* endptr = nullptr;
            nurbs_dim = static_cast<int>(strtol(str_nurbs_dim.c_str(), &endptr, 10));

            continue;
          }

          // check for a new patch
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
      if (!reader.my_output_flag())
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
    disknots = Teuchos::make_rcp<Core::FE::Nurbs::Knotvector>(nurbs_dim, npatches);

    // make sure that we have some Knotvector object to fill
    if (disknots == Teuchos::null)
    {
      FOUR_C_THROW("disknots should have been allocated before");
    }

    //--------------------------------------------------------------------
    //--------------------------------------------------------------------
    //                finally read knotvector section
    //--------------------------------------------------------------------
    //--------------------------------------------------------------------
    {
      // this is a pointer to the knots of one patch in one direction
      // we will read them and put them
      std::vector<Teuchos::RCP<std::vector<double>>> patch_knots(nurbs_dim);

      // temporary string
      std::string tmp;

      // start to read something when read is true
      bool read = false;

      // index for number of patch
      int npatch = 0;
      // index for u/v/w
      int actdim = -1;
      // ints for the number of knots
      std::vector<int> n_x_m_x_l(nurbs_dim);
      // ints for patches degrees
      std::vector<int> degree(nurbs_dim);
      // a vector of strings holding the knotvectortypes read
      std::vector<std::string> knotvectortype(nurbs_dim);

      // count for sanity check
      int count_read = 0;
      std::vector<int> count_vals(nurbs_dim);

      // loop lines in file
      for (const auto& line : reader.lines_in_section(sectionname))
      {
        std::istringstream file{std::string{line}};
        file >> tmp;

        // check for a new patch
        std::string::size_type loc = tmp.rfind("BEGIN");
        if (loc != std::string::npos)
        {
          file >> tmp;

          // activate reading
          read = true;

          actdim = -1;

          // create vectors for knots in this patch
          for (int rr = 0; rr < nurbs_dim; ++rr)
          {
            patch_knots[rr] = Teuchos::make_rcp<std::vector<double>>();
            (*(patch_knots[rr])).clear();
          }

          // reset counter for knot values
          for (int rr = 0; rr < nurbs_dim; rr++)
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
          if (actdim > nurbs_dim)
          {
            FOUR_C_THROW(
                "too many knotvectors, we only need one for each dimension (nurbs_dim = %d)\n",
                nurbs_dim);
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
          for (int rr = 0; rr < nurbs_dim; ++rr)
          {
            disknots->set_knots(
                rr, npatch, degree[rr], n_x_m_x_l[rr], knotvectortype[rr], patch_knots[rr]);
          }
          file >> tmp;
          // stop reading of knot values if we are here
          read = false;

          for (int rr = 0; rr < nurbs_dim; rr++)
          {
            if (n_x_m_x_l[rr] != count_vals[rr])
            {
              FOUR_C_THROW("not enough knots read in dim %d (%d!=NUMKNOTS=%d), nurbs_dim=%d\n", rr,
                  count_vals[rr], n_x_m_x_l[rr], nurbs_dim);
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
      }  // end loop through file

      if (count_read != npatches)
      {
        FOUR_C_THROW("wasn't able to read enough patches\n");
      }
    }

    if (myrank == 0)
    {
      if (!reader.my_output_flag())
      {
        Core::IO::cout << " in...." << time.totalElapsedTime(true) << " secs\n";

        time.reset();
        fflush(stdout);
      }
    }
  }



  /*----------------------------------------------------------------------*/
  /*----------------------------------------------------------------------*/
  void DatFileReader::read_dat()
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

    int arraysize = 0;
    if (comm_.MyPID() == 0)
    {
      std::ifstream file(filename_.c_str());
      if (not file) FOUR_C_THROW("unable to open file: %s", filename_.c_str());

      std::list<std::string> content;
      std::string current_excluded_section_name{};
      unsigned current_section_linecount = 0;

      const auto name_of_excluded_section = [&exclude](const std::string& section_header)
      {
        auto it = std::find_if(exclude.begin(), exclude.end(),
            [&section_header](const std::string& section)
            { return section_header.find(section) != std::string::npos; });

        if (it != exclude.end())
          return *it;
        else
          return std::string{};
      };

      std::string line;

      const auto finalize_section_read = [&](int number_of_lines)
      {
        if (!current_excluded_section_name.empty())
        {
          excludepositions_[current_excluded_section_name].second = number_of_lines;
        }

        const auto maybe_excluded_section = name_of_excluded_section(line);
        if (!maybe_excluded_section.empty())
        {
          // Start a new excluded section. This starts at the next line.
          excludepositions_[maybe_excluded_section].first = file.tellg();
          current_excluded_section_name = maybe_excluded_section;
        }
        else
        {
          current_excluded_section_name.clear();
        }

        current_section_linecount = 0;
      };

      // First loop over all input lines. This reads the actual file contents and determines
      // whether a line is to be read immediately or should be excluded because it is in one of
      // the excluded sections.
      while (getline(file, line))
      {
        ++current_section_linecount;

        // remove comments, trailing and leading whitespaces
        // compact internal whitespaces
        line = Core::Utils::strip_comment(line);

        // line is now empty
        if (line.size() == 0) continue;

        // This line starts a new section
        if (line.find("--") == 0)
        {
          // finalize the last excluded section. Subtract the current line which is not part of
          // the section anymore.
          finalize_section_read(current_section_linecount - 1);
        }

        if (!current_excluded_section_name.empty())
        {
          // We are in an excluded section. Skip the line.
          continue;
        }

        // This line contains content that we want to read now.
        content.push_back(line);
        // Count the number of characters in the line + 1 for the null-terminator we will add later.
        arraysize += static_cast<int>(line.length()) + 1;
      }
      // Finalize the last section
      finalize_section_read(current_section_linecount);

      // allocate space for copy of file
      inputfile_.clear();
      inputfile_.reserve(arraysize);
      lines_.reserve(content.size());

      std::size_t pos = 0u;
      for (auto& i : content)
      {
        inputfile_.insert(inputfile_.end(), i.begin(), i.end());
        // This fixup is crucial in the current implementation to get a null-terminated string
        // we can view with string_view without supplying a length.
        inputfile_.push_back('\0');
        lines_.emplace_back(&inputfile_[pos]);
        // Next entry would begin here.
        pos = inputfile_.size();
      }

      if (inputfile_.size() != static_cast<size_t>(arraysize))
      {
        FOUR_C_THROW(
            "internal error in file read: inputfile has %d chars, but was predicted to be %d "
            "chars long",
            inputfile_.size(), arraysize);
      }
    }

    // Now lets do all the parallel setup. Afterwards all processors
    // have to be the same.

    if (comm_.NumProc() > 1)
    {
      int num_lines = lines_.size();
      /* Now that we use a variable number of bytes per line we have to
       * communicate the buffer size as well. */
      comm_.Broadcast(&arraysize, 1, 0);
      comm_.Broadcast(&num_lines, 1, 0);

      if (comm_.MyPID() > 0)
      {
        /*--------------------------------------allocate space for copy of file */
        inputfile_.resize(arraysize);
        lines_.reserve(num_lines);
      }

      // There are no char based functions available! Do it by hand!
      // comm_.Broadcast(inputfile_.data(),arraysize,0);

      const auto& mpicomm = dynamic_cast<const Epetra_MpiComm&>(comm_);

      MPI_Bcast(inputfile_.data(), arraysize, MPI_CHAR, 0, mpicomm.GetMpiComm());

      if (comm_.MyPID() > 0)
      {
        // Chunk the raw input file into lines: whenever encountering a null-terminator,
        // add the line that just ended to the vector of lines.
        std::size_t pos = 0u;
        for (std::size_t i = 0; i < inputfile_.size(); ++i)
        {
          if (inputfile_[i] == '\0')
          {
            lines_.emplace_back(&inputfile_[pos]);
            pos = i + 1;
          }
        }

        FOUR_C_THROW_UNLESS(static_cast<int>(lines_.size()) == num_lines,
            "line count mismatch: %d lines expected but %d lines received", num_lines,
            lines_.size());
      }

      // distribute excluded section positions
      for (auto& i : exclude)
      {
        if (comm_.MyPID() == 0)
        {
          auto ep = excludepositions_.find(i);
          if (ep == excludepositions_.end())
          {
            excludepositions_[i] = std::pair<std::ifstream::pos_type, unsigned int>(-1, 0);
          }
        }
        std::pair<std::ifstream::pos_type, unsigned int>& p = excludepositions_[i];
        // comm_.Broadcast(&p.second,1,0);
        MPI_Bcast(&p.second, 1, MPI_INT, 0, mpicomm.GetMpiComm());
      }
    }

    // Now finally find the section names. We have to do this on all
    // processors, so it cannot be done while reading.
    std::vector<std::pair<std::size_t, std::string>> section_header_lines;
    for (std::vector<char*>::size_type i = 0; i < lines_.size(); ++i)
    {
      const std::string_view l = lines_[i];
      if (l[0] == '-' and l[1] == '-')
      {
        std::string line(l);

        // take the last "--" and all that follows as section name
        std::string::size_type loc = line.rfind("--");
        std::string sectionname = line.substr(loc);
        section_header_lines.emplace_back(i, sectionname);
      }
    }

    // Include a past-the-end pseudo section such that the last actual section can use the
    // past-the-end position for its position
    section_header_lines.emplace_back(lines_.size(), "unused");

    // Leave out the past-the-end section
    for (unsigned i = 0; i < section_header_lines.size() - 1; ++i)
    {
      const auto& [position, name] = section_header_lines[i];
      const auto& [next_position, _] = section_header_lines[i + 1];

      if (positions_.find(name) != positions_.end())
        FOUR_C_THROW("section '%s' defined more than once", name.c_str());

      // Shift by 1 to start after the header itself
      positions_[name] = {position + 1, next_position};
      // Remember the section name to later check if it was ever queried.
      knownsections_[name] = false;
    }

    // the following section names are always regarded as valid
    record_section_used("--TITLE");
    record_section_used("--FUNCT1");
    record_section_used("--FUNCT2");
    record_section_used("--FUNCT3");
    record_section_used("--FUNCT4");
    record_section_used("--FUNCT5");
    record_section_used("--FUNCT6");
    record_section_used("--FUNCT7");
    record_section_used("--FUNCT8");
    record_section_used("--FUNCT9");
    record_section_used("--FUNCT10");
    record_section_used("--FUNCT11");
    record_section_used("--FUNCT12");
    record_section_used("--FUNCT13");
    record_section_used("--FUNCT14");
    record_section_used("--FUNCT15");
    record_section_used("--FUNCT16");
    record_section_used("--FUNCT17");
    record_section_used("--FUNCT18");
    record_section_used("--FUNCT19");
    record_section_used("--FUNCT20");
  }


  /*----------------------------------------------------------------------*/
  /*----------------------------------------------------------------------*/
  bool DatFileReader::print_unknown_sections(std::ostream& out) const
  {
    const bool printout = std::any_of(
        knownsections_.begin(), knownsections_.end(), [](const auto& kv) { return !kv.second; });

    // now it's time to create noise on the screen
    if (printout and (get_comm().MyPID() == 0))
    {
      out << "\nERROR!"
          << "\n--------"
          << "\nThe following input file sections remained unused (obsolete or typo?):\n";
      for (const auto& [section_name, known] : knownsections_)
      {
        if (!known) out << section_name << '\n';
      }
      out << std::endl;
    }

    return printout;
  }

  void DatFileReader::record_section_used(const std::string& section_name)
  {
    knownsections_[section_name] = true;
  }

  std::pair<std::string, std::string> read_key_value(const std::string& line)
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
        FOUR_C_THROW("Line '%s' with just one word in parameter section", line.c_str());
    }

    std::string key = Core::Utils::trim(line.substr(0, separator_index));
    std::string value = Core::Utils::trim(line.substr(separator_index + 1));

    if (key.empty()) FOUR_C_THROW("Cannot get key from line '%s'", line.c_str());
    if (value.empty()) FOUR_C_THROW("Cannot get value from line '%s'", line.c_str());

    return {std::move(key), std::move(value)};
  }
}  // namespace Core::IO

FOUR_C_NAMESPACE_CLOSE
