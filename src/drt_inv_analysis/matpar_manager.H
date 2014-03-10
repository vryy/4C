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
#include "../drt_lib/drt_dserror.H"


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

    //! evaluate derivative of the residual force vector w.r.t to the material parameters using duals
    void Evaluate(double time, Teuchos::RCP<Epetra_MultiVector> dfint);

    //! evaluate derivative of the residual force vector w.r.t to the material parameters using FD
    void EvaluateFD(Teuchos::RCP<Epetra_MultiVector> dfint);

    //! update params
    void UpdateParams(Teuchos::RCP<Epetra_MultiVector> toadd);

    //! update params
    void ReplaceParams(Teuchos::RCP<Epetra_MultiVector> toreplace);

    //! reset params
    void ResetParams();

    //! Return material parameters
    const Teuchos::RCP<Epetra_MultiVector> GetMatParams(){SetParams(); return params_;};

    //! Return optimization parameters
    const Teuchos::RCP<Epetra_MultiVector> GetParams(){return optparams_;};

    //! Return optimization parameters  of the last run
    const Teuchos::RCP<Epetra_MultiVector> GetParamsOld(){return optparams_o_;};

    //! Return the map layout for the parameter distribution
    const Epetra_Map* ParamLayoutMap(){return paramlayoutmap_.get();};

    //! Return number of physically different parameters
    int NumParams(){ return numparams_; };

    //! Put vector index of a given material parameter <string> into loc
    int GetParameterLocation(int eleid, std::string name);

  protected:

    //! map of names for each material parameter
    //! key: matid, value:
    //! value is of enum type and can be used to set the respective material parameter directly
    //! i.e  < 1 ;  ( youngs, beta)>
    //!      < 2 ;  ( k1 , k2   )  >
    std::map<int, std::vector<int> > paramap_;

    //! map of matids to parameterids for the parameter vector kept in here
    //! i.e  < 1 ; (1 , 2)>
    //! i.e  < 2 ; (3 , 4)>
    std::map<int,std::vector<int> > parapos_;

    //! number of physically different material parameters
    //! -> number of vectors in params_
    int numparams_;

    //! discretization
    Teuchos::RCP<DRT::Discretization> discret_;

    //! actual set of material parameters in elementwise colmap layout
    Teuchos::RCP<Epetra_MultiVector> params_;

    //! actual set of material parameters in some parametrization layout
    Teuchos::RCP<Epetra_MultiVector> optparams_;

    //! old set of material parameters in some parametrization layout
    Teuchos::RCP<Epetra_MultiVector> optparams_o_;

    //! parametrization layout of the optimization parameters
    Teuchos::RCP<Epetra_Map> paramlayoutmap_;

    //! inidicates meta parmaetrization of the material parameters k_m = k^2+0.1
    bool metaparams_;

    //! construct map of material parameters to be optimized
    void SetupMatOptMap();

    //! set elementwise material parameters
    virtual void SetParams();

    //! initialze parameter vector
    void InitParams();

    //! let the vector of material parameters be filled by the derived classes
    virtual void FillParameters(Teuchos::RCP<Epetra_MultiVector> params) = 0;

    //! contract gradient wrt one physical material parameter (dfintp)
    //! to the optimization parameter space (dfint) for that physical parameter
    //! but in a distributed proc environment -> call Consolidate() finally
    virtual void ContractGradient(Teuchos::RCP<Epetra_MultiVector> dfint,
                                  double val,
                                  int elepos,
                                  int parapos) = 0;

    //! consolidate gradient among all procs.
    virtual void Consolidate(Teuchos::RCP<Epetra_MultiVector> dfint) = 0;

  };

}
}

/*----------------------------------------------------------------------*/
#endif
