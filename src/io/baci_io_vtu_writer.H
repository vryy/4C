/*----------------------------------------------------------------------*/
/*! \file

\brief VTU writer

\level 2

*/
/*----------------------------------------------------------------------*/
#ifndef BACI_IO_VTU_WRITER_H
#define BACI_IO_VTU_WRITER_H


#include "baci_config.H"

#include "baci_io_vtk_writer_base.H"

#include <map>
#include <string>
#include <vector>

BACI_NAMESPACE_OPEN

/*
 \brief class for VTU output generation

 \author kronbichler, grill
 \date 03/14, 03/17
*/
class VtuWriter : public VtkWriterBase
{
 public:
  //! constructor
  VtuWriter();

  //! write the geometry defining this unstructured grid
  void WriteGeometryUnstructuredGrid(const std::vector<double>& point_coordinates,
      const std::vector<IO::index_type>& point_cell_connectivity,
      const std::vector<IO::index_type>& cell_offset, const std::vector<uint8_t>& cell_types,
      const std::vector<IO::index_type>& face_connectivity,
      const std::vector<IO::index_type>& face_offset);


  //! write a data vector with num_component values of type T per point
  void WritePointDataVector(const IO::visualization_vector_type_variant& data,
      unsigned int num_components_per_point, const std::string& name);

  //! write a data vector with num_component values of type T per cell
  void WriteCellDataVector(const IO::visualization_vector_type_variant& data,
      unsigned int num_components_per_cell, const std::string& name);


 protected:
  //! Return the opening xml tag for this writer type
  const std::string& WriterOpeningTag() const override;

  //! Return the parallel opening xml tag for this writer type
  const std::string& WriterPOpeningTag() const override;

  //! Return a vector of parallel piece tags for each file
  const std::vector<std::string>& WriterPPieceTags() const override;

  //! Return the parallel file suffix including the dot for this file type
  const std::string& WriterPSuffix() const override;

  //! Return the string of this writer type
  const std::string& WriterString() const override;

  //! Return the file suffix including the dot for this file type
  const std::string& WriterSuffix() const override;
};

BACI_NAMESPACE_CLOSE

#endif
