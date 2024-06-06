/*----------------------------------------------------------------------*/
/*! \file

 \brief utility functions for poroelasticity problems

\level 2

 *----------------------------------------------------------------------*/

#ifndef FOUR_C_POROELAST_UTILS_HPP
#define FOUR_C_POROELAST_UTILS_HPP

#include "4C_config.hpp"

#include "4C_coupling_volmortar_utils.hpp"
#include "4C_discretization_fem_general_element.hpp"
#include "4C_inpar_poroelast.hpp"

#include <Epetra_Vector.h>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN

namespace Discret
{
  class Discretization;
}  // namespace Discret

namespace Core::Elements
{
  class Element;
}

namespace Core::LinAlg
{
  class MapExtractor;
}  // namespace Core::LinAlg

namespace Mat
{
  class Material;
}  // namespace Mat

namespace PoroElast
{
  class PoroBase;

  //! coupling types of porous media problem
  enum Coupltype
  {
    fluidfluid,
    fluidstructure,
    undefined
  };

  //! PoroElast::UTILS: Random stuff that might be helpful when dealing with poroelasticity problems
  namespace UTILS
  {
    //! check if element is a poro-element
    bool IsPoroElement(const Core::Elements::Element* actele);

    //! check if element is a poro-p1-element
    bool IsPoroP1Element(const Core::Elements::Element* actele);

    Teuchos::RCP<Core::LinAlg::MapExtractor> BuildPoroSplitter(
        Teuchos::RCP<Discret::Discretization> dis);

    //! create solution algorithm depending on input file
    Teuchos::RCP<PoroElast::PoroBase> CreatePoroAlgorithm(
        const Teuchos::ParameterList& timeparams,  //!< problem parameters (i)
        const Epetra_Comm& comm,                   //!< communicator(i)
        bool setup_solve = true,  //!< setup linear solver for Poroelastic problem (only required if
                                  //!< Solve() is called) (i)
        Teuchos::RCP<Core::LinAlg::MapExtractor> porosity_splitter =
            Teuchos::null  //! select porosity splitter
    );

    //! reset Material pointers after redistribution
    void SetMaterialPointersMatchingGrid(Teuchos::RCP<const Discret::Discretization> sourcedis,
        Teuchos::RCP<const Discret::Discretization> targetdis);

    /*!
     Create volume ghosting:
     To ease the contact search algorithms we'll afford the luxury to ghost all nodes
     on all processors in the general mortar coupling framework.
     Here we additionally ghost volume elements and their nodes
     (required if an evaluation of gradients is required)!
     */
    void create_volume_ghosting(
        Discret::Discretization& idiscret);  // redistributed interface discretization of contact!

    /*! Reconnect Face Element - Parent Element Pointers!
     Parent Element need to be ghosted on the processors where Face Elements
     exist already.
     */
    void reconnect_parent_pointers(Discret::Discretization& idiscret,
        Discret::Discretization& voldiscret, Discret::Discretization* voldiscret2 = nullptr);

    //! Determine norm of vector
    double calculate_vector_norm(const enum Inpar::PoroElast::VectorNorm norm,  //!< norm to use
        const Teuchos::RCP<const Epetra_Vector> vect  //!< the vector of interest
    );

    //! Set the slave and master elements of the face element
    void SetSlaveAndMaster(const Discret::Discretization& voldiscret,
        const Discret::Discretization* voldiscret2, const Epetra_Map* elecolmap,
        Core::Elements::FaceElement* faceele);

    //! strategy for material assignment for non matching meshes with poro

    //! Helper class for assigning materials for volumetric coupling of non conforming meshes (poro)
    /*!
     When coupling two overlapping discretizations, most often one discretization needs access
     to the corresponding element/material on the other side. For conforming meshes this is straight
     forward as there is one unique element on the other side and therefore one unique material,
     which can be accessed. However, for non conforming meshes there are potentially several
     elements overlapping. Therefore, some rule for assigning materials is needed. This class is
     meant to do that. It gets the element to which it shall assign a material and a vector of IDs
     of the overlapping elements of the other discretization.

     In case of poro, we also need the initial porosity type of the structural element to be known
     in the fluid element, which is why there is a special strategy for poro. Note that this is not
     yet working for inhomogeneous material properties.

     \author vuong 10/14
     */
    class PoroMaterialStrategy : public Core::VolMortar::UTILS::DefaultMaterialStrategy
    {
     public:
      //! constructor
      PoroMaterialStrategy() = default;

      //! assignment of fluid material to structure material
      void AssignMaterial2To1(const Core::VolMortar::VolMortarCoupl* volmortar,
          Core::Elements::Element* ele1, const std::vector<int>& ids_2,
          Teuchos::RCP<Discret::Discretization> dis1,
          Teuchos::RCP<Discret::Discretization> dis2) override;

      //! assignment of structure material to fluid material
      void AssignMaterial1To2(const Core::VolMortar::VolMortarCoupl* volmortar,
          Core::Elements::Element* ele2, const std::vector<int>& ids_1,
          Teuchos::RCP<Discret::Discretization> dis1,
          Teuchos::RCP<Discret::Discretization> dis2) override;
    };
  }  // namespace UTILS

  void PrintLogo();

}  // namespace PoroElast

FOUR_C_NAMESPACE_CLOSE

#endif
