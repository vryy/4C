/*----------------------------------------------------------------------*/
/*! \file

\brief VTU writer

\level 2

*/
/*----------------------------------------------------------------------*/
#ifndef BACI_IO_VTU_WRITER_H
#define BACI_IO_VTU_WRITER_H


#include "baci_io_vtk_writer_base.H"

#include <map>
#include <string>
#include <vector>


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

  //! destructor
  virtual ~VtuWriter() {}

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
  virtual const std::string& WriterOpeningTag() const;

  //! Return the parallel opening xml tag for this writer type
  virtual const std::string& WriterPOpeningTag() const;

  //! Return a vector of parallel piece tags for each file
  virtual const std::vector<std::string>& WriterPPieceTags() const;

  //! Return the parallel file suffix including the dot for this file type
  virtual const std::string& WriterPSuffix() const;

  //! Return the string of this writer type
  virtual const std::string& WriterString() const;

  //! Return the file suffix including the dot for this file type
  virtual const std::string& WriterSuffix() const;
};

#endif
