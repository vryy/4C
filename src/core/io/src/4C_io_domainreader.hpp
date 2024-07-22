/*----------------------------------------------------------------------*/
/*! \file

\brief Read domain sections of dat files.

\level 0


*/
/*----------------------------------------------------------------------*/



#ifndef FOUR_C_IO_DOMAINREADER_HPP
#define FOUR_C_IO_DOMAINREADER_HPP

#include "4C_config.hpp"

#include "4C_io_inputreader.hpp"

#include <Teuchos_RCP.hpp>

#include <memory>
#include <set>

FOUR_C_NAMESPACE_OPEN

namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE

namespace Core::IO
{
  class LineDefinition;

  namespace GridGenerator
  {
    struct RectangularCuboidInputs;
  }

  /*----------------------------------------------------------------------*/
  /*!
    \brief helper class to create a rectangular cuboid discretization of (currently) hex or wedge
    elements based on a few input parameters
   */
  /*----------------------------------------------------------------------*/
  class DomainReader
  {
    friend class MeshReader;

   public:
    /// construct element reader for a given field that reads a given section
    /*!
      Create empty discretization and append it to given field.

      \param dis (i) the new discretization
      \param reader (i) the dat file reader
      \param sectionname (i) the section that contains the element lines
     */
    DomainReader(Teuchos::RCP<Core::FE::Discretization> dis, const Core::IO::DatFileReader& reader,
        std::string sectionname);

    //! Destructor
    virtual ~DomainReader() = default;
    /// give the discretization this reader fills
    Teuchos::RCP<Core::FE::Discretization> my_dis() const { return dis_; }

   private:
    /*! \brief generate elements, partition node graph, create nodes
     *
     * This method processes DOMAIN SECTIONs. It
     * 1) generates elements
     * 2) partitions the element graph
     * 3) exports the col elements
     * 4) generates the node graph
     * 5) creates the nodes
     * 6) exports the col nodes
     *
     * \param nodeoffset [in] Node GID of first newly created node.
     * \note Node GIDs may not overlap
     */
    void create_partitioned_mesh(int nodeGIdOfFirstNewNode) const;

    /*! \brief read input parameters from input file
       \return class holding all input parameters for rectangular cuboid domain
     */
    Core::IO::GridGenerator::RectangularCuboidInputs read_rectangular_cuboid_input_data() const;

    /// finalize reading. fill_complete(false,false,false), that is, do not
    /// initialize elements. This is done later after reading boundary conditions.
    void complete() const;

    /// discretization name
    std::string name_;

    /// the main dat file reader
    const Core::IO::DatFileReader& reader_;

    /// my comm
    Teuchos::RCP<Epetra_Comm> comm_;

    /// my section to read
    std::string sectionname_;

    /// my discretization
    Teuchos::RCP<Core::FE::Discretization> dis_;
  };

}  // namespace Core::IO
FOUR_C_NAMESPACE_CLOSE

#endif
