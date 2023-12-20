/*----------------------------------------------------------------------*/
/*! \file
 \brief helper methods/classes for multiphase porous fluid flow problems

   \level 3

 *----------------------------------------------------------------------*/

#ifndef BACI_POROFLUIDMULTIPHASE_UTILS_H
#define BACI_POROFLUIDMULTIPHASE_UTILS_H

#include "baci_config.H"

#include "baci_inpar_bio.H"
#include "baci_inpar_porofluidmultiphase.H"
#include "baci_io.H"
#include "baci_linalg_fixedsizematrix.H"

#include <Epetra_Vector.h>
#include <Teuchos_RCP.hpp>

#include <set>

BACI_NAMESPACE_OPEN

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
        Teuchos::RCP<IO::DiscretizationWriter> output  //!< output writer
    );

    /**
     * \brief extend ghosting for artery discretization
     * @param[in] contdis  discretization of 2D/3D domain
     * @param[in] artdis   discretization of 1D domain
     * @param[in] evaluate_on_lateral_surface   is coupling evaluated on lateral surface?
     * @return             set of nearby element pairs as seen from the artery discretization, each
     *                     artery element with a vector of close 3D elements
     */
    std::map<int, std::set<int>> ExtendedGhostingArteryDiscretization(
        Teuchos::RCP<DRT::Discretization> contdis, Teuchos::RCP<DRT::Discretization> artdis,
        const bool evaluate_on_lateral_surface,
        const INPAR::ARTNET::ArteryPoroMultiphaseScatraCouplingMethod couplingmethod);

    /**
     * \brief get axis-aligned bounding box of element
     * @param[in] ele        compute AABB for this element
     * @param[in] positions  nodal positions of discretization
     * @param[in] evaluate_on_lateral_surface   is coupling evaluated on lateral surface?
     * @return               AABB of element as 3x2 CORE::LINALG::Matrix
     */
    CORE::LINALG::Matrix<3, 2> GetAABB(DRT::Element* ele,
        std::map<int, CORE::LINALG::Matrix<3, 1>>& positions,
        const bool evaluate_on_lateral_surface);

    /// maximum distance between two nodes of an element
    double GetMaxNodalDistance(DRT::Element* ele, Teuchos::RCP<DRT::Discretization> dis);

    /**
     * \brief perform octtree search for NTP coupling
     * @param[in] contdis       discretization of 2D/3D domain
     * @param[in] artdis        discretization of 1D domain
     * @param[in] artsearchdis  fully overlapping discretization of 1D domain on which search is
     * performed
     * @param[in] evaluate_on_lateral_surface   is coupling evaluated on lateral surface?
     * @param[in] artEleGIDs   vector of artery coupling element Ids
     * @param[out] elecolset    additional elements that need to be ghosted are filled into this set
     * @param[out] nodecolset   additional nodes that need to be ghosted are filled into this set
     * @return                  set of nearby element pairs as seen from the artery discretization,
     *                          each artery element with a vector of close 3D elements
     */
    std::map<int, std::set<int>> OctTreeSearch(Teuchos::RCP<DRT::Discretization> contdis,
        Teuchos::RCP<DRT::Discretization> artdis, Teuchos::RCP<DRT::Discretization> artsearchdis,
        const bool evaluate_on_lateral_surface, const std::vector<int> artEleGIDs,
        std::set<int>& elecolset, std::set<int>& nodecolset);

    /*!
     * \brief get nodal positions of discretization as std::map
     * @param[in] dis      discretization for which nodal positions will be returned
     * @param[in] nodemap  node-map of the discretization (can be either in row or column format)
     * @return             nodal position as a std::map<int, CORE::LINALG::Matrix<3, 1>>
     */
    std::map<int, CORE::LINALG::Matrix<3, 1>> GetNodalPositions(
        Teuchos::RCP<DRT::Discretization> dis, const Epetra_Map* nodemap);

    //! Determine norm of vector
    double CalculateVectorNorm(
        const enum INPAR::POROFLUIDMULTIPHASE::VectorNorm norm,  //!< norm to use
        const Teuchos::RCP<const Epetra_Vector> vect             //!< the vector of interest
    );

    /*!
     * create (fully overlapping) search discretization
     * @param artdis  1D artery discretization
     * @param disname  name of fully-overlapping artery discretization
     * @param doboundaryconditions  also do boundary conditions in fill-complete call
     * @return  fully-overlapping artery discretization
     */
    Teuchos::RCP<DRT::Discretization> CreateFullyOverlappingArteryDiscretization(
        Teuchos::RCP<DRT::Discretization> artdis, std::string disname, bool doboundaryconditions);

  }  // namespace UTILS
  // Print the logo
  void PrintLogo();
}  // namespace POROFLUIDMULTIPHASE


BACI_NAMESPACE_CLOSE

#endif  // POROFLUIDMULTIPHASE_UTILS_H
