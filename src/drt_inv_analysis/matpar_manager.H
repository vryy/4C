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

namespace STR
{

namespace INVANA

{



  /*!
  \brief A class to handle material parameters for optimization

  The MatParManager handels the mapping between the elemental material parameters and
  the vector of material parameters held in the optimization process.

  Material Parameters are stored as Epetra_MultiVector where each vector refers to a
  physical material property (e.g. youngs modulus for St. Venant Kirchhoff, alpha "stiffness"
  and beta "stiffness" in the AAA material). The layout of each vector is given by the
  underlying discretization's elecolmap.

  It furthermore handles the evaluation of the gradient of the problem residual w.r.t to the material
  parameters choosen to be identified.
 */
  class MatParManager
  {

  public:

    //! Constructor
    MatParManager(Teuchos::RCP<DRT::Discretization> discret);

    //!destructor
    virtual ~MatParManager() { ; };

    //! evaluate derivative of the residual force vector w.r.t to the material parameters
    void Evaluate(Teuchos::RCP<Epetra_MultiVector> dfint);

    //! update params
    void UpdateParams(Teuchos::RCP<Epetra_MultiVector> toadd);

    //! update params
    void ReplaceParams(Teuchos::RCP<Epetra_MultiVector> toreplace);

    //! reset params
    void ResetParams();

    //! Return the material parameters
    const Teuchos::RCP<Epetra_MultiVector> GetParams(){ return params_; };

    //! Return the material parameters of the last run
    const Teuchos::RCP<Epetra_MultiVector> GetParamsOld(){ return params_o_; };

    //! Return number of physically different parameters
    int NumParams(){ return numparams_; };

    //! Return weight of the regularization
    double GetRegWeight(){return reg_weight_;};

    //! Put vector index of a given material parameter <string> into loc
    int GetParameterLocation(int eleid, std::string name);

    //! brief compute optimization parameter multivector from an SMC Particle position
    virtual void ComputeParamsMultiVectorFromSMCParticlePosition(Teuchos::RCP<Epetra_MultiVector> &params , std::vector<double> my_global_params);
  protected:

    //! map of names for each material parameter
    //! key: matid, value: vector of names of the parameters
    //! i.e  < 1 ;  ( youngs, beta)>
    //!      < 2 ;  ( k1 , k2   )  >
    std::map<int, std::vector<std::string> > paramap_;

    //! map of matids to parameterids for the parameter vector
    //! i.e  < 1 ; (1 , 2)>
    //! i.e  < 2 ; (3 , 4)>
    std::map<int,std::vector<int> > parapos_;

    //! number of physically different material parameters
    //! -> number of vectors in params_
    int numparams_;

    //! discretization
    Teuchos::RCP<DRT::Discretization> discret_;

    //! discret_'s elecolmap for the distribution of parameters
    const Epetra_Map* elecolmap_;

    //! actual set of parameters
    Teuchos::RCP<Epetra_MultiVector> params_;

    //! old set of parameters
    Teuchos::RCP<Epetra_MultiVector> params_o_;

    //! initial parameters as initial guess:
    Teuchos::RCP<Epetra_MultiVector> params_init_;

    //! weight of the regularization:
    double reg_weight_;

    //! inidicates meta parmaetrization of the material parameters k_m = k^2+0.1
    bool metaparams_;

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
