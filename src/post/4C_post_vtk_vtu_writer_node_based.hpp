/*----------------------------------------------------------------------*/
/*! \file

\brief node based VTU filter

\level 3

*/
/*----------------------------------------------------------------------*/
#ifndef FOUR_C_POST_VTK_VTU_WRITER_NODE_BASED_HPP
#define FOUR_C_POST_VTK_VTU_WRITER_NODE_BASED_HPP


#include "4C_config.hpp"

#include "4C_discretization_fem_general_element.hpp"
#include "4C_post_vtk_vtu_writer.hpp"

#include <map>
#include <string>
#include <vector>


// forward declarations
FOUR_C_NAMESPACE_OPEN
class PostField;
class PostResult;

namespace Discret
{
  class Discretization;

  namespace ELEMENTS
  {
    class Beam3Base;
  }
}  // namespace Discret

namespace Core::Nodes
{
  class Node;
}

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
  const std::string& writer_opening_tag() const override;

  //! Return the parallel opening xml tag for this writer type
  const std::string& writer_p_opening_tag() const override;

  //! Return a vector of parallel piece tags for each file
  const std::vector<std::string>& writer_p_piece_tags() const override;

  //! Give every writer a chance to do preparations before writing
  void writer_prep_timestep() override{};

  //! Return the parallel file suffix including the dot for this file type
  const std::string& writer_p_suffix() const override;

  //! Return the string of this writer type
  const std::string& writer_string() const override;

  //! Return the file suffix including the dot for this file type
  const std::string& writer_suffix() const override;

  //! Write a single result step
  void write_dof_result_step(std::ofstream& file, const Teuchos::RCP<Epetra_Vector>& data,
      std::map<std::string, std::vector<std::ofstream::pos_type>>& resultfilepos,
      const std::string& groupname, const std::string& name, const int numdf, const int from,
      const bool fillzeros) override;

  //! Write a single result step
  void write_nodal_result_step(std::ofstream& file, const Teuchos::RCP<Epetra_MultiVector>& data,
      std::map<std::string, std::vector<std::ofstream::pos_type>>& resultfilepos,
      const std::string& groupname, const std::string& name, const int numdf) override;

  //! Write a single result step
  void write_element_result_step(std::ofstream& file, const Teuchos::RCP<Epetra_MultiVector>& data,
      std::map<std::string, std::vector<std::ofstream::pos_type>>& resultfilepos,
      const std::string& groupname, const std::string& name, const int numdf,
      const int from) override;

  //! write the geometry of one time step
  void write_geo() override;

  //! write the geometry of Nurbs Element
  virtual void write_geo_nurbs_ele(const Core::Elements::Element* ele,
      std::vector<uint8_t>& celltypes, int& outNodeId, std::vector<int32_t>& celloffset,
      std::vector<double>& coordinates);

  //! write the geometry of beam element (special treatment due to Hermite interpolation)
  void write_geo_beam_ele(const Discret::ELEMENTS::Beam3Base* beamele,
      std::vector<uint8_t>& celltypes, int& outNodeId, std::vector<int32_t>& celloffset,
      std::vector<double>& coordinates) override;

  //! Write a single result step for one Nurbs Element
  virtual void wirte_dof_result_step_nurbs_ele(const Core::Elements::Element* ele, int ncomponents,
      const int numdf, std::vector<double>& solution, Teuchos::RCP<Epetra_Vector> ghostedData,
      const int from, const bool fillzeros);

  void write_dof_result_step_beam_ele(const Discret::ELEMENTS::Beam3Base* beamele,
      const int& ncomponents, const int& numdf, std::vector<double>& solution,
      Teuchos::RCP<Epetra_Vector>& ghostedData, const int& from, const bool fillzeros) override;

  //! Write a single result step for one Nurbs Element
  virtual void write_nodal_result_step_nurbs_ele(const Core::Elements::Element* ele,
      int ncomponents, const int numdf, std::vector<double>& solution,
      Teuchos::RCP<Epetra_MultiVector> ghostedData);
};

FOUR_C_NAMESPACE_CLOSE

#endif
