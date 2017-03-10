/*----------------------------------------------------------------------*/
/*!
\file particleMeshFree_rendering.H

\brief Handler of the rendering-related stuff

\level 1

\maintainer Cattabiani Alessandro
*/

#ifndef SRC_DRT_PARTICLE_PARTICLEMESHFREE_RENDERING_H_
#define SRC_DRT_PARTICLE_PARTICLEMESHFREE_RENDERING_H_

/*----------------------------------------------------------------------*/
/* headers */
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_StandardParameterEntryValidators.hpp>
#include <Epetra_Vector.h>
#include <Epetra_FEVector.h>
#include <Epetra_Time.h>

#include "../drt_lib/drt_globalproblem.H"
#include "particle_node.H"
#include "particleMeshFree_weightFunction.H"
#include "../drt_inpar/inpar_particle.H"

// forward declarations
namespace DRT
{
  class Discretization;
}
namespace PARTICLE
{
  class Algorithm;
}

namespace PARTICLE
{
  class Rendering
  {

  public:

    //! @name Life
    //@{

    //! Constructor (the writer should already be set)
    Rendering(Teuchos::RCP<PARTICLE::Algorithm> particleAlgorithm,
              INPAR::PARTICLE::WeightFunctionDim wf_dim);

    //! Destructor
    virtual ~Rendering() { ; }

    //@}

    //! update the rendering state vectors
    virtual void UpdateStateVectors(
        Teuchos::RCP<DRT::Discretization> pDiscret,
        Teuchos::RCP<const Epetra_Vector> pDis,
        Teuchos::RCP<const Epetra_Vector> pVel,
        Teuchos::RCP<const Epetra_Vector> pAcc,
        Teuchos::RCP<const Epetra_Vector> pDensity,
        Teuchos::RCP<const Epetra_Vector> pSpecEnthalpy,
        Teuchos::RCP<const Epetra_Vector> pTemperature,
        Teuchos::RCP<const Epetra_Vector> pRadius,
        Teuchos::RCP<const Epetra_Vector> pPressure,
        Teuchos::RCP<const Epetra_Vector> pMass);

    //! Rendering output state
    virtual void OutputState();

  private:

    //! search for the neighbouring rendering nodes
    std::list<DRT::Node*> GetNeighbouringRenderingNodes(const int binId);

    //!< rendering discretisation
    Teuchos::RCP<DRT::Discretization> discret_;

    //!< link to the big particle_algorithm box. From this we can access the particle discretization
    Teuchos::RCP<PARTICLE::Algorithm> particle_algorithm_;

    //! weight function pointer
    Teuchos::RCP<PARTICLE::WeightFunction_Base> weightFunctionHandler_;

    bool trg_writeMesh_;

    ///< bins 2 rendering nodes map
    std::map<int,std::list<DRT::Node* > > bins2renderingNodes_;

    //! rendering vectors
    //! rendering velocity
    Teuchos::RCP<Epetra_FEVector> vel_;
    //! rendering acceleration
    Teuchos::RCP<Epetra_FEVector> acc_;
    //! rendering density
    Teuchos::RCP<Epetra_FEVector> density_;
    //! rendering specific enthalpy
    Teuchos::RCP<Epetra_FEVector> specEnthalpy_;
    //! rendering temeprature
    Teuchos::RCP<Epetra_FEVector> temperature_;
    //! rendering pressure
    Teuchos::RCP<Epetra_FEVector> pressure_;

  };

}



#endif /* SRC_DRT_PARTICLE_PARTICLEMESHFREE_RENDERING_H_ */
