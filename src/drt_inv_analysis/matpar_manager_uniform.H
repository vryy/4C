/*----------------------------------------------------------------------*/
/*!
\file matpar_manager_uniform.H
\brief Class to handle calls to material parameters from an optimization routine

<pre>
Maintainer: Sebastian Kehl
            kehl@mhpc.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-10361
</pre>
*/

/*----------------------------------------------------------------------*/
/* definitions */
#ifndef MATPAR_MANAGER_UNIFORM_H
#define MATPAR_MANAGER_UNIFORM_H

/*----------------------------------------------------------------------*/
/* headers */
#include "matpar_manager.H"


// forward declarations
namespace DRT
{
  class Discretization;
}

namespace STR
{

namespace INVANA

{

  class MatParManagerUniform  : public MatParManager
  {

  protected:

    void FillParameters(Teuchos::RCP<Epetra_MultiVector> params);

    void InitParameters(int parapos, double val);

    void ContractGradient(Teuchos::RCP<Epetra_MultiVector> dfint, double val, int elepos, int paraposglobal, int paraposlocal);

    void Consolidate(Teuchos::RCP<Epetra_MultiVector> dfint);

  public:

    //! constructor
    MatParManagerUniform(Teuchos::RCP<DRT::Discretization> discret);

    //! destructor
    ~MatParManagerUniform(){};

    int NumVectors(){return NumParams();};

  };

  class MatParManagerPerElement  : public MatParManager
  {

  private:

    //! Map Extractor to extract Blocks of parameters of the parameter vector according to the value in paramap_
    Teuchos::RCP<LINALG::MultiMapExtractor> paramapextractor_;

    //! the two following 'maps' are processor specific: each proc has its own
    //! map keeping correspondence between params in the parameter vector (parameter LID in optparams_) and
    //! eleids (i.e gids of the discretizations ele layout)
    std::vector<int> paramsLIDtoeleGID_;

    //! map keeping correspondence between eleid (GID) + parameter block id and parameter LID
    std::map< int, std::vector<int> > eleGIDtoparamsLID_;

  protected:

    void FillParameters(Teuchos::RCP<Epetra_MultiVector> params);

    void InitParameters(int parapos, double val);

    void ContractGradient(Teuchos::RCP<Epetra_MultiVector> dfint, double val, int elepos, int paraposglobal, int paraposlocal);

    //! nothing to be done in case of elementwise parameters
    void Consolidate(Teuchos::RCP<Epetra_MultiVector> dfint){return;}

  public:

    //! constructor
    MatParManagerPerElement(Teuchos::RCP<DRT::Discretization> discret);

    //! destructor
    ~MatParManagerPerElement(){};

    int NumVectors(){return 1;};

  };

}
}

/*----------------------------------------------------------------------*/
#endif
