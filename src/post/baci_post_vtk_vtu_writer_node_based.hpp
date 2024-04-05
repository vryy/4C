/*----------------------------------------------------------------------*/
/*! \file

\brief node based VTU filter

\level 3

*/
/*----------------------------------------------------------------------*/
#ifndef FOUR_C_POST_VTK_VTU_WRITER_NODE_BASED_HPP
#define FOUR_C_POST_VTK_VTU_WRITER_NODE_BASED_HPP


#include "baci_config.hpp"

#include "baci_lib_element.hpp"
#include "baci_post_vtk_vtu_writer.hpp"

#include <map>
#include <string>
#include <vector>


// forward declarations
BACI_NAMESPACE_OPEN
class PostField;
class PostResult;

namespace DRT
{
  class Discretization;
  class Node;

  namespace ELEMENTS
  {
    class Beam3Base;
  }
}  // namespace DRT

/*
 \brief Base class for VTU node based output generation

 \author biehler
 \date 12/16
*/
class PostVtuWriterNode : public PostVtuWriter
{
 public:
  //! constructor. Initializes the writer to a certain field.
  PostVtuWriterNode(PostField* field, const std::string& name);

 protected:
  //! Return the opening xml tag for this writer type
  const std::string& WriterOpeningTag() const override;

  //! Return the parallel opening xml tag for this writer type
  const std::string& WriterPOpeningTag() const override;

  //! Return a vector of parallel piece tags for each file
  const std::vector<std::string>& WriterPPieceTags() const override;

  //! Give every writer a chance to do preparations before writing
  void WriterPrepTimestep() override{};

  //! Return the parallel file suffix including the dot for this file type
  const std::string& WriterPSuffix() const override;

  //! Return the string of this writer type
  const std::string& WriterString() const override;

  //! Return the file suffix including the dot for this file type
  const std::string& WriterSuffix() const override;

  //! Write a single result step
  void WriteDofResultStep(std::ofstream& file, const Teuchos::RCP<Epetra_Vector>& data,
      std::map<std::string, std::vector<std::ofstream::pos_type>>& resultfilepos,
      const std::string& groupname, const std::string& name, const int numdf, const int from,
      const bool fillzeros) override;

  //! Write a single result step
  void WriteNodalResultStep(std::ofstream& file, const Teuchos::RCP<Epetra_MultiVector>& data,
      std::map<std::string, std::vector<std::ofstream::pos_type>>& resultfilepos,
      const std::string& groupname, const std::string& name, const int numdf) override;

  //! Write a single result step
  void WriteElementResultStep(std::ofstream& file, const Teuchos::RCP<Epetra_MultiVector>& data,
      std::map<std::string, std::vector<std::ofstream::pos_type>>& resultfilepos,
      const std::string& groupname, const std::string& name, const int numdf,
      const int from) override;

  //! write the geometry of one time step
  void WriteGeo() override;

  //! write the geometry of Nurbs Element
  virtual void WriteGeoNurbsEle(const DRT::Element* ele, std::vector<uint8_t>& celltypes,
      int& outNodeId, std::vector<int32_t>& celloffset, std::vector<double>& coordinates);

  //! write the geometry of beam element (special treatment due to Hermite interpolation)
  void WriteGeoBeamEle(const DRT::ELEMENTS::Beam3Base* beamele, std::vector<uint8_t>& celltypes,
      int& outNodeId, std::vector<int32_t>& celloffset, std::vector<double>& coordinates) override;

  //! Write a single result step for one Nurbs Element
  virtual void WirteDofResultStepNurbsEle(const DRT::Element* ele, int ncomponents, const int numdf,
      std::vector<double>& solution, Teuchos::RCP<Epetra_Vector> ghostedData, const int from,
      const bool fillzeros);

  void WriteDofResultStepBeamEle(const DRT::ELEMENTS::Beam3Base* beamele, const int& ncomponents,
      const int& numdf, std::vector<double>& solution, Teuchos::RCP<Epetra_Vector>& ghostedData,
      const int& from, const bool fillzeros) override;

  //! Write a single result step for one Nurbs Element
  virtual void WriteNodalResultStepNurbsEle(const DRT::Element* ele, int ncomponents,
      const int numdf, std::vector<double>& solution, Teuchos::RCP<Epetra_MultiVector> ghostedData);
};

BACI_NAMESPACE_CLOSE

#endif
