/*----------------------------------------------------------------------*/
/*! \file

\brief VTU filter

\level 2

*/
/*----------------------------------------------------------------------*/
#ifndef FOUR_C_POST_VTK_VTU_WRITER_HPP
#define FOUR_C_POST_VTK_VTU_WRITER_HPP


#include "4C_config.hpp"

#include "4C_discretization_fem_general_element.hpp"
#include "4C_post_vtk_writer.hpp"

#include <map>
#include <string>
#include <vector>


FOUR_C_NAMESPACE_OPEN

// forward declarations
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
 \brief Base class for VTU output generation

 \author kronbichler
 \date 03/14
*/
class PostVtuWriter : public PostVtkWriter
{
 public:
  //! constructor. Initializes the writer to a certain field.
  PostVtuWriter(PostField* field, const std::string& name);

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
  virtual void write_geo_nurbs_ele(const CORE::Elements::Element* ele,
      std::vector<uint8_t>& celltypes, int& outNodeId, std::vector<int32_t>& celloffset,
      std::vector<double>& coordinates) const;

  /*! Generalization of the former non-template method for all implemented NURBS
   *  discretization types
   *
   *  \author hiermeier (originally Seitz) \date 10/17 */
  template <CORE::FE::CellType nurbs_type>
  void write_geo_nurbs_ele(const CORE::Elements::Element* ele, std::vector<uint8_t>& celltypes,
      int& outNodeId, std::vector<int32_t>& celloffset, std::vector<double>& coordinates) const;

  CORE::FE::CellType map_nurbs_dis_type_to_lagrange_dis_type(
      const CORE::FE::CellType nurbs_dis_type) const;

  //! write the geometry of beam element (special treatment due to Hermite interpolation)
  virtual void write_geo_beam_ele(const DRT::ELEMENTS::Beam3Base* beamele,
      std::vector<uint8_t>& celltypes, int& outNodeId, std::vector<int32_t>& celloffset,
      std::vector<double>& coordinates);

  //! Write a single result step for one Nurbs Element
  virtual void wirte_dof_result_step_nurbs_ele(const CORE::Elements::Element* ele, int ncomponents,
      const int numdf, std::vector<double>& solution, Teuchos::RCP<Epetra_Vector> ghostedData,
      const int from, const bool fillzeros) const;

  /*! Generalization of the former non-template method for all implemented NURBS
   *  discretization types
   *
   *  \author hiermeier (originally Seitz) \date 10/17 */
  template <CORE::FE::CellType nurbs_type>
  void wirte_dof_result_step_nurbs_ele(const CORE::Elements::Element* ele, int ncomponents,
      const int numdf, std::vector<double>& solution, Teuchos::RCP<Epetra_Vector> ghostedData,
      const int from, const bool fillzeros) const;

  virtual void write_dof_result_step_beam_ele(const DRT::ELEMENTS::Beam3Base* beamele,
      const int& ncomponents, const int& numdf, std::vector<double>& solution,
      Teuchos::RCP<Epetra_Vector>& ghostedData, const int& from, const bool fillzeros);

  //! Write a single result step for one Nurbs Element
  virtual void write_nodal_result_step_nurbs_ele(const CORE::Elements::Element* ele,
      int ncomponents, const int numdf, std::vector<double>& solution,
      Teuchos::RCP<Epetra_MultiVector> ghostedData) const;

  /*! Generalization of the former non-template method for all implemented NURBS
   *  discretization types
   *
   *  \author hiermeier (originally Seitz) \date 10/17 */
  template <CORE::FE::CellType nurbs_type>
  void write_nodal_result_step_nurbs_ele(const CORE::Elements::Element* ele, int ncomponents,
      const int numdf, std::vector<double>& solution,
      Teuchos::RCP<Epetra_MultiVector> ghostedData) const;

  ///! width of the proc identifier number in the processor specific file names
  const int proc_file_padding_;
};

FOUR_C_NAMESPACE_CLOSE

#endif
