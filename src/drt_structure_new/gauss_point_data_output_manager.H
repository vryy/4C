/*----------------------------------------------------------------------*/
/*! \file
\brief Container for output data of the gauss point level

\level 3
*/
/*----------------------------------------------------------------------*/
#ifndef GAUSS_POINT_DATA_OUTPUT_MANAGER_H
#define GAUSS_POINT_DATA_OUTPUT_MANAGER_H

#include <memory>
#include <unordered_map>
#include <vector>
#include <Teuchos_RCP.hpp>
#include <Epetra_IntVector.h>

#include "../drt_inpar/inpar_structure.H"

// forward declarations
class Epetra_MultiVector;

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

      void AddQuantityIfNotExistant(const std::string& name, int size);

      void MergeQuantities(const std::unordered_map<std::string, int>& quantities);

      void AddElementNumberOfGaussPoints(int numgp);

      void PrepareData(const Epetra_Map& node_col_map, const Epetra_Map& element_row_map);

      void PostEvaluate();

      /*!
       * \brief Distribute and collect all quantities to and from all other procs to ensure that all
       * data is in this list.
       */
      void DistributeQuantities(const Epetra_Comm& comm);

      inline std::unordered_map<std::string, Teuchos::RCP<Epetra_MultiVector>>&
      GetMutableNodalData()
      {
        return data_nodes_;
      }

      inline std::unordered_map<std::string, Teuchos::RCP<Epetra_IntVector>>&
      GetMutableNodalDataCount()
      {
        return data_nodes_count_;
      }

      inline std::unordered_map<std::string, Teuchos::RCP<Epetra_MultiVector>>&
      GetMutableElementCenterData()
      {
        return data_element_center_;
      }

      inline std::unordered_map<std::string, std::vector<Teuchos::RCP<Epetra_MultiVector>>>&
      GetMutableGaussPointData()
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
      GetElementCenterData() const
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
      static const char MPI_DELIMITER = '!';

      void SendMyQuantitiesToProc(const DRT::Exporter&, int to_proc) const;

      std::unique_ptr<std::unordered_map<std::string, int>> ReceiveQuantitiesFromProc(
          const DRT::Exporter&, int from_proc) const;

      void BroadcastMyQuantitites(const DRT::Exporter&);

      void PackMyQuantities(std::vector<char>& data) const;

      void UnpackQuantities(std::size_t pos, const std::vector<char>& data,
          std::unordered_map<std::string, int>& quantities) const;

      void PrepareNodalDataVectors(const Epetra_Map& nodal_col_map);

      void PrepareElementCenterDataVectors(const Epetra_Map& nodal_col_map);

      void PrepareGaussPointDataVectors(const Epetra_Map& nodal_col_map);

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


#endif  // GAUSS_POINT_DATA_OUTPUT_MANAGER_H