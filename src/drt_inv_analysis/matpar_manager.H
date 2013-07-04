/*----------------------------------------------------------------------*/
/*!
\file matpar_manager.H
\brief Class to handle calls to material parameters from an optimization routine

<pre>
Maintainer: Sebastian Kehl
            kehl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15249
</pre>
*/

/*----------------------------------------------------------------------*/
/* definitions */
#ifndef MATPAR_MANAGER_H
#define MATPAR_MANAGER_H

/*----------------------------------------------------------------------*/
/* headers */
#include "Teuchos_RCP.hpp"
#include "Epetra_Vector.h"
#include "../drt_inpar/inpar_material.H"
#include "Epetra_SerialDenseVector.h"


// forward declarations
namespace DRT
{
  class Discretization;
}

/*----------------------------------------------------------------------*/
//! Material parameters manager
namespace STR
{

namespace INVANA

{
  class MatParManager
  {

  public:

    //! Constructor
    MatParManager(Teuchos::RCP<DRT::Discretization> discret);

    //!destructor
    virtual ~MatParManager() { ; };

    //! put derivative of internal force vector into dfint
    //! numvectors: number of parameters (e.g. numelements)
    //! length of each vector comes from the discretizations dofrowmap
    void Evaluate(Teuchos::RCP<Epetra_MultiVector> dfint, Teuchos::RCP<Epetra_Vector> disdual);

    //! update params
    void UpdateParams(Teuchos::RCP<Epetra_MultiVector> toadd);

    //! reset params
    void ResetParams();

    //! Return the material parameters
    const Teuchos::RCP<Epetra_MultiVector> GetParams(){ return params_; };

    //! Return the material parameters
    const Teuchos::RCP<Epetra_MultiVector> GetInitialGuess(){ return params_init_; };

    //! Return weight of the regularization
    double GetRegWeight(){return reg_weight_;};


  protected:

    //! type of material
    //INPAR::MAT::MaterialType mattype_;

    //! map of names for each material parameter
    //! key: matid, value: vector of names of the parameters
    std::map<int, std::vector<std::string> > paramap_;

    //! map of matids to parameterids for the parameter vector
    std::map<int,std::vector<int> > parapos_;

    //! size of the parameter vector
    int sizeparams_;

    //! size of the parameter vector
    int numparams_;

    //! discretization
    Teuchos::RCP<DRT::Discretization> discret_;

    //! a dofrowmap of the above discretization
    const Epetra_Map* dofrowmap_;

    //! a elecolmap for the distribution of parameters
    const Epetra_Map* elecolmap_;

    //! actual set of parameters
    Teuchos::RCP<Epetra_MultiVector> params_;

    //! old set of parameters
    Teuchos::RCP<Epetra_MultiVector> params_o_;

    //! initial parameters as initial guess:
    Teuchos::RCP<Epetra_MultiVector> params_init_;

    //! weight of the regularization:
    double reg_weight_;

    //! construct map of material parameters to be optimized
    void SetupMatOptMap();

    //! set elementwise material parameters
    virtual void SetParams();

    //! initialze parameter vector
    void InitParams();

  };

}
}

/*----------------------------------------------------------------------*/
#endif
