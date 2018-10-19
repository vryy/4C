/*----------------------------------------------------------------------*/
/*!
\file drt_inputreader.H

\brief Internal classes to read elements and nodes

<pre>
\level 0

\maintainer Martin Kronbichler
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>
*/
/*----------------------------------------------------------------------*/



#ifndef DRT_INPUTREADER_H
#define DRT_INPUTREADER_H

#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <list>
#include <string>
#include <algorithm>

#include <Epetra_Comm.h>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

namespace DRT
{
  namespace NURBS
  {
    class Knotvector;
  }

  namespace INPUT
  {
    /*----------------------------------------------------------------------*/
    /// reading, broadcasting and storing of dat file contents
    /*!
      The normal baci job needs to read one dat file at the beginning of
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
      is one instance of this class during the input phase of baci. When
      everything is set up properly this instance goes away and the
      internal copy of the input file vanishes.

      There is one more details to be considered: we do not want to
      read all elements and nodes on each processor. So we skip reading
      these sections here. Elements and nodes are read in parallel by the
      ElementReader and NodeReader, respectively.

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
      /// empty constructor, which avoids directly opening a file for reading
      DatFileReader();

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
        NodeReader and the ElementReader.
       */
      std::ifstream::pos_type ExcludedSectionPosition(std::string section) const;

      /// get number of lines of excluded section
      unsigned int ExcludedSectionLength(std::string section) const;

      /// return my communicator
      Teuchos::RCP<Epetra_Comm> Comm() const { return comm_; }

      /// convert a parameter section in extended format in a parameter list
      bool ReadSection(std::string name, Teuchos::ParameterList& list);

      /// convert a parameter section in GiD format in a parameter list
      bool ReadGidSection(std::string name, Teuchos::ParameterList& list);

      /// return the content of a section
      std::vector<const char*> Section(std::string name);

      /// read a node-design topology section
      void ReadDesign(const std::string& name, std::vector<std::vector<int>>& dobj_fenode);

      /*!
        \brief read the knotvector section (for isogeometric analysis)

        \param  dim            (in ): dimension, i.e. number of knotvectors in one patch
        \param  name           (in ): Name/type of discretisation
        \param  disknots       (out): node vector coordinates

      */
      void ReadKnots(
          const int dim, const std::string name, Teuchos::RCP<DRT::NURBS::Knotvector>& disknots);

      /// dump input file contents to error file
      void DumpInput();

      /// print unknown section names found in the input file
      bool PrintUnknownSections();

     private:
      /// parse name string an get sublist of global input list
      Teuchos::ParameterList& FindSublist(std::string name, Teuchos::ParameterList& list);

      /// parse the value and add the key,value pair to the list
      void AddEntry(std::string key, std::string value, Teuchos::ParameterList& list);

      /// actual read dat file, store and broadcast general sections
      void ReadDat();

      /// input file name
      std::string filename_;

      /// my communicator
      Teuchos::RCP<Epetra_Comm> comm_;

      /// number of lines
      int numrows_;

      /// flag for output (default: output should be written)
      int outflag_;

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


  }  // namespace INPUT
}  // namespace DRT

#endif
