/*----------------------------------------------------------------------*/
/*! \file
\brief Container for output data of the gauss point level

\level 3
*/
/*----------------------------------------------------------------------*/
#ifndef FOUR_C_STRUCTURE_NEW_GAUSS_POINT_DATA_OUTPUT_MANAGER_HPP
#define FOUR_C_STRUCTURE_NEW_GAUSS_POINT_DATA_OUTPUT_MANAGER_HPP

#include "4C_config.hpp"

#include "4C_comm_exporter.hpp"
#include "4C_inpar_structure.hpp"

#include <Epetra_IntVector.h>
#include <Teuchos_RCP.hpp>

#include <memory>
#include <unordered_map>
#include <vector>

// forward declarations
class Epetra_MultiVector;

FOUR_C_NAMESPACE_OPEN

namespace DRT
{
  class Exporter;
}

namespace STR
{
  namespace MODELEVALUATOR
  {
    class GaussPointDataOutputManager
    {
     public:
      explicit GaussPointDataOutputManager(INPAR::STR::GaussPointDataOutputType output_type);

      void add_quantity_if_not_existant(const std::string& name, int size);

      void MergeQuantities(const std::unordered_map<std::string, int>& quantities);

      void add_element_number_of_gauss_points(int numgp);

      void PrepareData(const Epetra_Map& node_col_map, const Epetra_Map& element_row_map);

      void post_evaluate();

      /*!
       * \brief Distribute and collect all quantities to and from all other procs to ensure that all
       * data is in this list.
       */
      void distribute_quantities(const Epetra_Comm& comm);

      inline std::unordered_map<std::string, Teuchos::RCP<Epetra_MultiVector>>& GetNodalData()
      {
        return data_nodes_;
      }

      inline std::unordered_map<std::string, Teuchos::RCP<Epetra_IntVector>>& GetNodalDataCount()
      {
        return data_nodes_count_;
      }

      inline std::unordered_map<std::string, Teuchos::RCP<Epetra_MultiVector>>&
      get_element_center_data()
      {
        return data_element_center_;
      }

      inline std::unordered_map<std::string, std::vector<Teuchos::RCP<Epetra_MultiVector>>>&
      GetGaussPointData()
      {
        return data_gauss_point_;
      }

      inline const std::unordered_map<std::string, std::vector<Teuchos::RCP<Epetra_MultiVector>>>&
      GetGaussPointData() const
      {
        return data_gauss_point_;
      }
      inline const std::unordered_map<std::string, Teuchos::RCP<Epetra_MultiVector>>& GetNodalData()
          const
      {
        return data_nodes_;
      }

      inline const std::unordered_map<std::string, Teuchos::RCP<Epetra_IntVector>>&
      GetNodalDataCount() const
      {
        return data_nodes_count_;
      }

      inline const std::unordered_map<std::string, Teuchos::RCP<Epetra_MultiVector>>&
      get_element_center_data() const
      {
        return data_element_center_;
      }

      inline const std::unordered_map<std::string, int>& GetQuantities() const
      {
        return quantities_;
      }

      inline INPAR::STR::GaussPointDataOutputType GetOutputType() const { return output_type_; }


     private:
      static constexpr int MPI_TAG = 545;
      static constexpr char MPI_DELIMITER = '!';

      void send_my_quantities_to_proc(const CORE::COMM::Exporter& exporter, int to_proc) const;

      std::unique_ptr<std::unordered_map<std::string, int>> receive_quantities_from_proc(
          const CORE::COMM::Exporter& exporter, int from_proc) const;

      void broadcast_my_quantitites(const CORE::COMM::Exporter& exporter);

      void pack_my_quantities(std::vector<char>& data) const;

      void unpack_quantities(std::size_t pos, const std::vector<char>& data,
          std::unordered_map<std::string, int>& quantities) const;

      void prepare_nodal_data_vectors(const Epetra_Map& node_col_map);

      void prepare_element_center_data_vectors(const Epetra_Map& element_col_map);

      void prepare_gauss_point_data_vectors(const Epetra_Map& element_col_map);

      //! output type of the data
      INPAR::STR::GaussPointDataOutputType output_type_;

      //! maximum number of Gauss points of all elements
      int max_num_gp_;

      //! map holding element data projected to the nodes
      std::unordered_map<std::string, Teuchos::RCP<Epetra_MultiVector>> data_nodes_;

      //! map holding the number of elements that share a quantity at each node
      std::unordered_map<std::string, Teuchos::RCP<Epetra_IntVector>> data_nodes_count_;

      //! map holding element data averaged to the element center
      std::unordered_map<std::string, Teuchos::RCP<Epetra_MultiVector>> data_element_center_;

      //! map holding element data for each Gauss point
      std::unordered_map<std::string, std::vector<Teuchos::RCP<Epetra_MultiVector>>>
          data_gauss_point_;

      //! unordered map holding the quantities and its sizes
      std::unordered_map<std::string, int> quantities_;
    };
  }  // namespace MODELEVALUATOR
}  // namespace STR


FOUR_C_NAMESPACE_CLOSE

#endif
