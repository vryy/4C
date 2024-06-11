/*----------------------------------------------------------------------*/
/*! \file

\brief Internal classes to read elements and nodes

\level 0


*/
/*----------------------------------------------------------------------*/


#ifndef FOUR_C_IO_INPUTREADER_HPP
#define FOUR_C_IO_INPUTREADER_HPP

#include "4C_config.hpp"

#include "4C_fem_discretization.hpp"

#include <Epetra_Comm.h>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <list>
#include <map>
#include <set>
#include <string>

FOUR_C_NAMESPACE_OPEN

namespace Core::FE
{
  namespace Nurbs
  {
    class Knotvector;
  }
}  // namespace Core::FE

namespace Core::IO
{
  /*----------------------------------------------------------------------*/
  /// reading, broadcasting and storing of dat file contents
  /*!
    The normal 4C job needs to read one dat file at the beginning of
    its execution. This file contains all input parameters as well as
    the meshes and so on and needs to be available on all
    processors. However, we cannot rely on a shared file system, so only
    on process (the first one) can do the actual reading. All others
    need to get the data via MPI communication.

    This class manages the reading and broadcasting of the input
    file. To do so processor 0 reads the input file line by line and
    stores all lines in an array. Afterwards these lines are
    communicated to all other processors, so each processor has an
    internal copy of the input file and subsequently parses it.  There
    is one instance of this class during the input phase of 4C. When
    everything is set up properly this instance goes away and the
    internal copy of the input file vanishes.

    There is one more details to be considered: we do not want to
    read all elements and nodes on each processor. So we skip reading
    these sections here. Elements and nodes are read in parallel by the
    ElementReader and MeshReader, respectively.

    The point of having this reader is to allow for special problem
    types that need to read more than one dat file. (Scary!) This can be
    done be creating a new (additional) reader to read the additional
    file.

    \author u.kue
    \date 06/07
   */
  /*----------------------------------------------------------------------*/
  class DatFileReader
  {
   public:
    /// Empty constructor, which avoids directly opening a file for reading
    DatFileReader() = default;

    //! constructs a dat file reader object and stores the file name handed in. The file specified
    //! by filename is not directly read
    explicit DatFileReader(std::string filename);

    /// construct a reader for a given file
    DatFileReader(std::string filename, Teuchos::RCP<Epetra_Comm> comm, int outflag = 0);

    /// return my inputfile name
    std::string MyInputfileName() const;

    /// return my output flag
    int MyOutputFlag() const;

    /// get file position of excluded section
    /*!
      The big sections (nodes and elements) are excluded from general
      read so that we can read them in parallel. This is done using the
      MeshReader and the ElementReader.
     */
    std::ifstream::pos_type excluded_section_position(const std::string& section) const;

    /// get number of lines of excluded section
    unsigned int excluded_section_length(const std::string& section) const;

    /// return my communicator
    Teuchos::RCP<Epetra_Comm> Comm() const { return comm_; }

    /// convert a parameter section in extended format in a parameter list
    bool ReadSection(std::string name, Teuchos::ParameterList& list);

    /// return the content of a section
    std::vector<const char*> Section(const std::string& name);

    /// Read a node-design topology section
    ///
    /// @param name Name of the topology to read
    /// @param dobj_fenode Resulting collection of all nodes that belong to a design.
    /// @param get_discretization Callback to return a discretization by name.
    void ReadDesign(const std::string& name, std::vector<std::vector<int>>& dobj_fenode,
        const std::function<const Core::FE::Discretization&(const std::string& name)>&
            get_discretization);

    /*!
      \brief read the knotvector section (for isogeometric analysis)

      \param  name           (in ): Name/type of discretisation
      \param  disknots       (out): node vector coordinates

    */
    void ReadKnots(const std::string& name, Teuchos::RCP<Core::FE::Nurbs::Knotvector>& disknots);


    /// print unknown section names found in the input file
    bool print_unknown_sections();

   private:
    /// parse the value and add the key,value pair to the list
    void add_entry(const std::string& key, const std::string& value, Teuchos::ParameterList& list);

    /// actual read dat file, store and broadcast general sections
    void read_dat();

    /// input file name
    std::string filename_;

    /// my communicator
    Teuchos::RCP<Epetra_Comm> comm_;

    /// number of lines
    int numrows_{};

    /// flag for output (default: output should be written)
    int outflag_{};

    /// the whole input file, lines separated by 0
    std::vector<char> inputfile_;

    /// pointers to line beginnings
    std::vector<char*> lines_;

    /// file positions of skipped sections
    std::map<std::string, std::pair<std::ifstream::pos_type, unsigned int>> excludepositions_;

    /// section positions of all sections
    std::map<std::string, unsigned int> positions_;

    /// protocol of known and unknown section names
    std::map<std::string, bool> knownsections_;

    // a cache of the box-domain specifications
    std::map<std::string, std::vector<double>> cached_box_specifications_;
  };

  /**
   * Split the given @p line into a key-value pair. Key and value are normally separated by
   * whitespace. In case there are multiple distinct whitespace groups in one line, the first of
   * these is assumed to be the separator and all the other whitespace is assumed to be part of
   * the value. Key and value may also be separated by an equals sign "=" and at least one
   * whitespace character on both sides. In this case, key and value may contain spaces
   * internally. Leading and trailing whitespace is trimmed from both key and value.
   *
   * @throws Core::Exception If the @p line cannot be read.
   *
   * @return A pair of key and value.
   */
  std::pair<std::string, std::string> ReadKeyValue(const std::string& line);


}  // namespace Core::IO

FOUR_C_NAMESPACE_CLOSE

#endif
