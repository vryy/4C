/*----------------------------------------------------------------------*/
/*! \file

\brief Read element sections of dat files.

\level 0

*/
/*----------------------------------------------------------------------*/

#include "4C_io_elementreader.hpp"

#include "4C_comm_utils_factory.hpp"
#include "4C_discretization_fem_general_element_definition.hpp"
#include "4C_rebalance_binning_based.hpp"
#include "4C_rebalance_print.hpp"

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
IO::ElementReader::ElementReader(Teuchos::RCP<DRT::Discretization> dis,
    const INPUT::DatFileReader& reader, std::string sectionname)
    : name_(dis->Name()),
      reader_(reader),
      comm_(reader.Comm()),
      sectionname_(sectionname),
      dis_(dis)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
IO::ElementReader::ElementReader(Teuchos::RCP<DRT::Discretization> dis,
    const INPUT::DatFileReader& reader, std::string sectionname, std::string elementtype)
    : name_(dis->Name()),
      reader_(reader),
      comm_(reader.Comm()),
      sectionname_(sectionname),
      dis_(dis)
{
  elementtypes_.insert(elementtype);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
IO::ElementReader::ElementReader(Teuchos::RCP<DRT::Discretization> dis,
    const INPUT::DatFileReader& reader, std::string sectionname,
    const std::set<std::string>& elementtypes)
    : name_(dis->Name()),
      reader_(reader),
      comm_(reader.Comm()),
      sectionname_(sectionname),
      dis_(dis)
{
  std::copy(elementtypes.begin(), elementtypes.end(),
      std::inserter(elementtypes_, elementtypes_.begin()));
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IO::ElementReader::ReadAndDistribute()
{
  const int myrank = comm_->MyPID();
  const int numproc = comm_->NumProc();

  // read global ids of elements of this discretization
  const auto& [numele, eids] = get_element_size_and_i_ds();

  // determine a preliminary element distribution
  int nblock, mysize, bsize;
  {
    if (numele == 0)
    {
      // If the element section is empty, we create an empty reader and return
      coleles_ = roweles_ = colnodes_ = rownodes_ =
          Teuchos::rcp(new Epetra_Map(-1, 0, nullptr, 0, *comm_));

      return;
    }

    // number of element chunks to split the reading process in
    // approximate block size (just a guess!)
    nblock = numproc;
    bsize = numele / nblock;

    // create a simple (pseudo linear) map
    mysize = bsize;
    if (myrank == numproc - 1) mysize = numele - (numproc - 1) * bsize;

    // construct the map
    roweles_ = Teuchos::rcp(new Epetra_Map(-1, mysize, &eids[myrank * bsize], 0, *comm_));
  }

  // define blocksizes for blocks of elements we read
  {
    // for block sizes larger than about 250000 elements (empirical value !) the code sometimes
    // hangs during ExportRowElements call for the second block (block 1). Therefore an upper limit
    // of 100000 for bsize is ensured below.
    const int maxblocksize = 100000;

    if (bsize > maxblocksize)
    {
      // without an additional increase of nblock by 1 the last block size
      // could reach a maximum value of (2*maxblocksize)-1, potentially
      // violating the intended upper limit!
      nblock = 1 + numele / maxblocksize;
      bsize = maxblocksize;
    }
  }

  get_and_distribute_elements(nblock, bsize);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::pair<int, std::vector<int>> IO::ElementReader::get_element_size_and_i_ds() const
{
  // vector of all global element ids
  std::vector<int> eids;
  int numele = 0;
  std::string inputfile_name = reader_.MyInputfileName();

  // all reading is done on proc 0
  if (comm_->MyPID() == 0)
  {
    // open input file at the right position
    std::ifstream file(inputfile_name.c_str());
    std::ifstream::pos_type pos = reader_.excluded_section_position(sectionname_);
    if (pos != std::ifstream::pos_type(-1))
    {
      file.seekg(pos);

      // loop all element lines, comments in the element section are not supported!
      std::string line;
      while (getline(file, line))
      {
        if (line.find("--") == 0)
          break;
        else
        {
          std::istringstream t;
          t.str(line);
          int elenumber;
          std::string eletype;
          t >> elenumber >> eletype;
          elenumber -= 1;

          // only read registered element types or all elements if nothing is registered
          if (elementtypes_.size() == 0 or elementtypes_.count(eletype) > 0)
            eids.push_back(elenumber);
        }
      }
      numele = static_cast<int>(eids.size());
    }
  }

  // Simply allreduce the element ids
  comm_->Broadcast(&numele, 1, 0);

  eids.resize(numele);
  comm_->Broadcast(&eids[0], numele, 0);

  return {numele, eids};
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IO::ElementReader::get_and_distribute_elements(const int nblock, const int bsize)
{
  std::ifstream file;
  std::string inputfile_name = reader_.MyInputfileName();

  if (comm_->MyPID() == 0)
  {
    file.open(inputfile_name.c_str());
    file.seekg(reader_.excluded_section_position(sectionname_));
  }
  std::string line;
  bool endofsection = false;

  CORE::Elements::ElementDefinition ed;
  ed.setup_valid_element_lines();

  for (int block = 0; block < nblock; ++block)
  {
    std::vector<int> gidlist;
    if (!endofsection && comm_->MyPID() == 0)
    {
      gidlist.reserve(bsize);
      int bcount = 0;
      while (getline(file, line))
      {
        if (line.find("--") == 0)
        {
          // If we have an empty element section (or fewer elements
          // than processors) we cannot read on. But we cannot exit
          // the block loop beforehand either because the other
          // processors need to syncronize nblock times.
          endofsection = true;
          break;
        }
        else
        {
          std::istringstream t;
          t.str(line);
          int elenumber;
          std::string eletype;
          std::string distype;
          // read element id type and distype
          t >> elenumber >> eletype >> distype;
          elenumber -= 1;
          gidlist.push_back(elenumber);

          // only read registered element types or all elements if nothing is
          // registered
          if (elementtypes_.size() == 0 or elementtypes_.count(eletype) > 0)
          {
            // let the factory create a matching empty element
            Teuchos::RCP<CORE::Elements::Element> ele =
                CORE::COMM::Factory(eletype, distype, elenumber, 0);
            if (ele.is_null()) FOUR_C_THROW("element creation failed");

            // For the time being we support old and new input facilities. To
            // smooth transition.

            INPUT::LineDefinition* linedef = ed.ElementLines(eletype, distype);
            if (linedef != nullptr)
            {
              if (not linedef->Read(t))
              {
                std::cout << "\n" << elenumber << " " << eletype << " " << distype << " ";
                linedef->Print(std::cout);
                std::cout << "\n";
                std::cout << line << "\n";
                FOUR_C_THROW(
                    "failed to read element %d %s %s", elenumber, eletype.c_str(), distype.c_str());
              }

              ele->SetNodeIds(distype, linedef);
              ele->ReadElement(eletype, distype, linedef);
            }
            else
            {
              FOUR_C_THROW("a matching line definition is needed for %s %s", eletype.c_str(),
                  distype.c_str());
            }

            // add element to discretization
            dis_->add_element(ele);

            // get the node ids of this element
            const int numnode = ele->num_node();
            const int* nodeids = ele->NodeIds();

            // all node gids of this element are inserted into a set of
            // node ids --- it will be used later during reading of nodes
            // to add the node to one or more discretisations
            std::copy(nodeids, nodeids + numnode, std::inserter(nodes_, nodes_.begin()));

            ++bcount;
            if (block != nblock - 1)
            {
              if (bcount == bsize)
              {
                break;
              }
            }
          }
        }
      }
    }
    dis_->proc_zero_distribute_elements_to_all(*roweles_, gidlist);
  }
}

FOUR_C_NAMESPACE_CLOSE
