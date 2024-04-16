/*----------------------------------------------------------------------*/
/*! \file
\brief  Coupling Communication Manager automatically creates all required coupling object to
transform matrixes, vectors, ...

\level 2


*----------------------------------------------------------------------*/

#include "baci_config.hpp"

#include "baci_coupling_adapter.hpp"

#include <Epetra_Comm.h>
#include <Epetra_Vector.h>
#include <Teuchos_ParameterList.hpp>

#define COUP_MANAGER_DEBUG_OUT

#ifndef FOUR_C_FSI_XFEM_COUPLING_COMM_MANAGER_HPP
#define FOUR_C_FSI_XFEM_COUPLING_COMM_MANAGER_HPP

BACI_NAMESPACE_OPEN

// forward declarations
namespace DRT
{
  class Discretization;
}

namespace CORE::LINALG
{
  class MultiMapExtractor;
  class SparseMatrix;
  class MatrixRowTransform;
  class MatrixColTransform;
  class MatrixRowColTransform;
  class MatrixLogicalSplitAndTransform;
}  // namespace CORE::LINALG

namespace CORE::ADAPTER
{
  class Coupling;
  class CouplingConverter;
}  // namespace CORE::ADAPTER

namespace XFEM
{
  class Coupling_Comm_Manager
  {
   public:
    enum transfer_type  // partial is on the interfacee //full is the whole discretisation //global
                        // is for all discretizations
    {
      full_to_full,
      full_to_partial,
      partial_to_full,
      partial_to_partial,
      partial_to_global,
      full_to_global,
    };

    enum matrix_transfer_type
    {
      row,
      col,
      row_and_col
    };

    //! constructor
    explicit Coupling_Comm_Manager(std::map<int, Teuchos::RCP<const DRT::Discretization>> dis,
        std::string cond_name, int startdim = 0, int enddim = 3);

    //! constructor
    explicit Coupling_Comm_Manager(Teuchos::RCP<const DRT::Discretization> dis0,
        std::string cond_name, int startdim = 0, int enddim = 3);

    //! constructor
    explicit Coupling_Comm_Manager(Teuchos::RCP<const DRT::Discretization> dis0,
        Teuchos::RCP<const DRT::Discretization> dis1, std::string cond_name, int startdim = 0,
        int enddim = 3);

    //! virtual destructor to support polymorph destruction
    virtual ~Coupling_Comm_Manager() = default;

    //! Insert a Vector A into vector B (choose type of transfer, add or scaling) - Version vor
    //! RCP<const Epetra_Vector> vecA
    void InsertVector(const int idxA, Teuchos::RCP<const Epetra_Vector> vecA, const int idxB,
        Teuchos::RCP<Epetra_Vector> vecB, const Coupling_Comm_Manager::transfer_type ttype,
        bool add = false, double scale = 1.0);

    //! Insert a Vector A into vector B (choose type of transfer, add or scaling) - Version vor
    //! RCP<Epetra_Vector> vecA
    void InsertVector(const int idxA, Teuchos::RCP<Epetra_Vector> vecA, const int idxB,
        Teuchos::RCP<Epetra_Vector> vecB, const Coupling_Comm_Manager::transfer_type ttype,
        bool add = false, double scale = 1.0)
    {
      InsertVector(
          idxA, Teuchos::rcp_static_cast<const Epetra_Vector>(vecA), idxB, vecB, ttype, add, scale);
    }

    //! Insert a Matrix A (from Discretization A) into Matrix B (from Discretization B) (choose type
    //! of transfer, add or scaling)
    bool InsertMatrix(
        int transform_id,  // Unique Id to be set for this transformation object (to be save use
                           // different one, for different matrix transformation)
        int idxA, const CORE::LINALG::SparseMatrix& matA, int idxB,
        CORE::LINALG::SparseMatrix& matB, const Coupling_Comm_Manager::matrix_transfer_type mttype,
        double scale = 1.0, bool exactmatch = true, bool addmatrix = false);

    Teuchos::RCP<CORE::LINALG::MultiMapExtractor> GetMapExtractor(int idx);

   protected:
    Teuchos::RCP<CORE::ADAPTER::CouplingConverter> GetCouplingConverter(int idxA, int idxB);

    Teuchos::RCP<CORE::ADAPTER::Coupling> GetCoupling(int idxA, int idxB);

    Teuchos::RCP<CORE::LINALG::MultiMapExtractor> GetFullMapExtractor() { return fullextractor_; }

    Teuchos::RCP<CORE::LINALG::MatrixLogicalSplitAndTransform> GetTransform(int transform_id);

    void DebugOut(
        std::string str1, std::string str2 = "", std::string str3 = "", std::string str4 = "");

   private:
    void Setup(std::map<int, Teuchos::RCP<const DRT::Discretization>> dis);

    void SetupMultiMapExtractors(std::map<int, Teuchos::RCP<const DRT::Discretization>> dis);

    void SetupFullMapExtractors(std::map<int, Teuchos::RCP<const DRT::Discretization>> dis);

    void SetupCouplings(std::map<int, Teuchos::RCP<const DRT::Discretization>> dis);

    void SetupFullCouplings(std::map<int, Teuchos::RCP<const DRT::Discretization>> dis);

    void SetupFullExtractor(std::map<int, Teuchos::RCP<const DRT::Discretization>> dis);

    std::string cond_name_;
    int startdim_;
    int enddim_;

    // All MultiMapExtractors
    std::map<int, Teuchos::RCP<CORE::LINALG::MultiMapExtractor>> mme_;

    // Couling Objects will just be initizalized in case we have more discretizations!
    std::map<std::pair<int, int>, Teuchos::RCP<CORE::ADAPTER::Coupling>> coup_;

    // Transformation Objects will just be initizalized in case we use matrix transformations!
    std::map<int, Teuchos::RCP<CORE::LINALG::MatrixLogicalSplitAndTransform>> transform_;

    Teuchos::RCP<CORE::LINALG::MultiMapExtractor> fullextractor_;
  };
}  // namespace XFEM

BACI_NAMESPACE_CLOSE

#endif
