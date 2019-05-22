/*----------------------------------------------------------------------*/
/*!
 \brief helper methods/classes for multiphase porous fluid flow problems

   \level 3

   \maintainer  Johannes Kremheller
 *----------------------------------------------------------------------*/

#ifndef SRC_DRT_POROFLUIDMULTIPHASE_POROFLUIDMULTIPHASE_UTILS_H_
#define SRC_DRT_POROFLUIDMULTIPHASE_POROFLUIDMULTIPHASE_UTILS_H_


#include "Teuchos_RCP.hpp"
#include <Epetra_Vector.h>

#include "../drt_inpar/inpar_porofluidmultiphase.H"
#include "../drt_io/io.H"

namespace DRT
{
  class Discretization;
  class Element;
}  // namespace DRT

namespace ADAPTER
{
  class PoroFluidMultiphase;
}

namespace POROFLUIDMULTIPHASE
{
  /// POROFLUIDMULTIPHASE::UTILS: Random stuff that might be helpful when dealing with
  /// poromultiphase problems
  namespace UTILS
  {
    /// setup second materials for porosity evaluation within solid phase
    void SetupMaterial(const Epetra_Comm& comm, const std::string& struct_disname,
        const std::string& fluid_disname);


    /// convert a dof based vector to a node based multi vector
    /*!
      For postprocessing, only vectors based on the primary dof set
      of the discretization can be used. Hence, all other EpetraVectors
      based on secondary dof sets are copied to EpetraMultiVectors with one
      node based vector for each component.

      This method can be deleted, if the post processors would be adapted to
      handle secondary dof sets.

      \param dis              : discretization, the vector is based on
      \param vector           : vector to convert
      \param nds              : number of the dof set the map of the vector corresponds to
      \param numdofpernode    : number of dofs per node of the vector (assumed to be equal for all
      nodes)
     */
    Teuchos::RCP<Epetra_MultiVector> ConvertDofVectorToNodeBasedMultiVector(
        const DRT::Discretization& dis, const Epetra_Vector& vector, const int nds,
        const int numdofpernode);

    /// create solution algorithm depending on input file
    Teuchos::RCP<ADAPTER::PoroFluidMultiphase> CreateAlgorithm(
        INPAR::POROFLUIDMULTIPHASE::TimeIntegrationScheme
            timintscheme,                              //!< time discretization scheme
        Teuchos::RCP<DRT::Discretization> dis,         //!< discretization
        const int linsolvernumber,                     //!< number of linear solver
        const Teuchos::ParameterList& probparams,      //!< parameter list of global problem
        const Teuchos::ParameterList& poroparams,      //!< paramter list of poro problem
        FILE* errfile,                                 //!< error file
        Teuchos::RCP<IO::DiscretizationWriter> output  //!< output writer
    );

    /**
     * extend ghosting for artery discretization
     * @param contdis [in]: discretization of 2D/3D domain
     * @param artdis [in]: discretization of 1D domain
     * @return set of nearby element pairs as seen from the artery discretization, each artery
     * element with a vector of close 3D elements
     */
    std::map<int, std::set<int>> ExtendedGhostingArteryDiscretization(
        Teuchos::RCP<DRT::Discretization> contdis, Teuchos::RCP<DRT::Discretization> artdis);

    /**
     * get axis-aligned bounding box of element
     * @param ele: compute AABB for this element
     * @return AABB of element as 6x1 vector with format [xmin, ymin, zmin, xmax, ymax, zmax]
     */
    std::vector<double> GetAABB(DRT::Element* ele);

    /// maximum distance between two nodes of an element
    double GetMaxNodalDistance(DRT::Element* ele, Teuchos::RCP<DRT::Discretization> dis);

    /**
     * perform brute force search
     * @param contdis [in]: discretization of 2D/3D domain
     * @param artdis [in]: discretization of 1D domain
     * @param artsearchdis [in]: fully overlapping discretization of 1D domain on which search is
     * performed
     * @param elecolset [out]: additional elements that need to be ghosted are filled into this set
     * @param nodecolset [out]: additional nodes that need to be ghosted are filled into this set
     * @return set of nearby element pairs as seen from the artery discretization, each artery
     * element with a vector of close 3D elements
     */
    std::map<int, std::set<int>> BruteForceSearch(Teuchos::RCP<DRT::Discretization> contdis,
        Teuchos::RCP<DRT::Discretization> artdis, Teuchos::RCP<DRT::Discretization> artsearchdis,
        std::set<int>* elecolset, std::set<int>* nodecolset);

    //! Determine norm of vector
    double CalculateVectorNorm(
        const enum INPAR::POROFLUIDMULTIPHASE::VectorNorm norm,  //!< norm to use
        const Teuchos::RCP<const Epetra_Vector> vect             //!< the vector of interest
    );

    //! create (fully overlapping) search discretization
    Teuchos::RCP<DRT::Discretization> CreateSearchDiscretization(
        Teuchos::RCP<DRT::Discretization> artdis);

  }  // namespace UTILS
  // Print the logo
  void PrintLogo();
}  // namespace POROFLUIDMULTIPHASE


#endif /* SRC_DRT_POROFLUIDMULTIPHASE_POROFLUIDMULTIPHASE_UTILS_H_ */
