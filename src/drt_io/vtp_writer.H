/*----------------------------------------------------------------------*/
/*! \file

\brief VTP writer

\level 2

\maintainer Martin Kronbichler
*/
/*----------------------------------------------------------------------*/
#ifndef IO_VTP_WRITER_H
#define IO_VTP_WRITER_H


#include "vtk_writer_base.H"

#include <string>
#include <map>
#include <vector>


/*
 \brief class for VTP output generation

 \author hammerl, grill
 \date 03/16, 04/17
*/
class VtpWriter : public VtkWriterBase
{
 public:
  //! constructor
  VtpWriter();

  //! destructor
  virtual ~VtpWriter() {}

  //! write the geometry defining this poly data set
  void WriteGeometryPolyData(const std::vector<double>& point_coordinates);

  //! write a data vector with num_component values of type T per point
  // Todo template <typename T>
  void WritePointDataVector(const std::vector<double>& data, unsigned int num_components_per_point,
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
