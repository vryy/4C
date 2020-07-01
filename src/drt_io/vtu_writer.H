/*----------------------------------------------------------------------*/
/*! \file

\brief VTU writer

\level 2

*/
/*----------------------------------------------------------------------*/
#ifndef IO_VTU_WRITER_H
#define IO_VTU_WRITER_H


#include "vtk_writer_base.H"

#include <string>
#include <map>
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

  //! write the geometry defining this unstructured grid assuming contiguous points
  void WriteGeometryUnstructuredGridContiguous(const std::vector<double>& point_coordinates,
      const std::vector<int32_t>& cell_offset, const std::vector<uint8_t>& cell_types);

  //! write the geometry defining this unstructured grid
  void WriteGeometryUnstructuredGrid(const std::vector<double>& point_coordinates,
      const std::vector<int32_t>& point_cell_connectivity, const std::vector<int32_t>& cell_offset,
      const std::vector<uint8_t>& cell_types);


  //! write a data vector with num_component values of type T per point
  // Todo template <typename T>
  void WritePointDataVector(const std::vector<double>& data, unsigned int num_components_per_point,
      const std::string& name);

  //! write a data vector with num_component values of type T per cell
  // Todo template <typename T>
  void WriteCellDataVector(const std::vector<double>& data, unsigned int num_components_per_cell,
      const std::string& name);


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
