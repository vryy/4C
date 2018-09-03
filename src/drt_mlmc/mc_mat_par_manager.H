/*----------------------------------------------------------------------------*/
/*!
\file mc_mat_par_manager.H
\brief Class to handle calls to material parameters from a uq routine

<pre>
Maintainer: Jonas Biehler
            biehler@lnm.mw.tum.de
            089 - 289-15276
</pre>
*/

/*----------------------------------------------------------------------------*/
#ifdef HAVE_FFTW
#ifndef MC_MATPAR_MANAGER_H
#define MC_MATPAR_MANAGER_H

/*----------------------------------------------------------------------------*/
/* header inclusions */
#include "Epetra_MultiVector.h"
#include "Teuchos_RCP.hpp"
#include "../drt_lib/drt_dserror.H"

// forward declarations
namespace DRT
{
  class Discretization;
}

namespace LINALG
{
  class MultiMapExtractor;
}


namespace UQ
{
  class RandomField;
  class RandomVariable;
  /*----------------------------------------------------------------------------*/
  /// Class to handle calls to material parameters from an uq routine
  class MCMatParManager
  {
   public:
    /*!
    \brief standard constructor

    \param discret : discretization for which MatParManager will be used
    */
    MCMatParManager(Teuchos::RCP<DRT::Discretization> discret);

    /*!
    \brief standard destructor

    */
    ~MCMatParManager(){};


    /*!
    \brief setup all random quantities needed for MC

    This funtion just calls SetupRandomFields and SetupRandomVariables

    \param myseed : seed for random number generator

    */
    void SetupAllRandomQuantities(unsigned int myseed);

    /*!
    \brief setup all random fields needed for MC

    \param myseed : seed for random number generator

    */
    void SetupRandomFields(unsigned int myseed);

    /*!
    \brief setup all random variables needed for MC and store them in map

    \param myseed (in): seed for random number generator

    */
    void SetupRandomVariables(unsigned int myseed);

    /*!
    \brief create new realization of all random fields

    create new realization of all random fields and
    push parameters down to material

    \param myseed              (in): seed for random number generator

    \param para_cont_parameter (in): value between 0 and 1 for parameter
                                     continuation scheme

    \param reuse_rf            (in): in case of parameter continuation we can
                                     reuse the computed realization of the
                                     random field. If we want that set to yes
    */
    void SetUpStochMats(unsigned int myseed, double para_cont_parameter, bool reuse_rf);

    /*!
    \brief write all random variables to file

    write all random variables that are used for the simulation
    of stochastic processes to file

    \param filename (in): file to write data to

    \param numrun   (in): which run in the MC analysis are we at

    */
    void WriteRandomVariablesToFile(std::string filename, int numrun);

    /*!
    \brief return number of physical parameters which are modeled by a random field
    */
    int NumPhysStochParams() { return numstochparams_r_field_; };


   private:
    /*!
    \brief initialze maps that contains <matid> <vector of stochastic params>

    */
    void InitStochParaMaps();

    /*!
    \brief Compute params vector and set vectors of respective structure materials

     \param para_cont_parameter (in): value between 0 and 1 for parameter
                                      continuation scheme
    */
    void SetParamsStructure(double para_cont_parameter);

    /*!
    \brief Compute params vector and set parameters of respective airway materials

     \param para_cont_parameter (in): value between 0 and 1 for parameter
                                      continuation scheme
    */
    void SetParamsRedAirways(double para_cont_parameter);

    /*!
    \brief Compute a new realisation for all random fields

     \param myseed (in): seed for random number generator

    */
    void CreateNewRealizationOfRandomQuantities(unsigned int myseed);

    /*!
    \brief Compute material parameter multivector

    Compute multivector with ele col layout which contains the values
    of the uncertain material parameters computed by random fields

    \param para_cont_parameter (in):  parameter continuation parameter [0-1]
    */
    void ComputeMatParamsMultivectorFromRandomFields(const double para_cont_parameter);

    /*!
    \brief Compute material parameter multivector

    Compute multivector with which contains the values
    of the uncertain material parameters computed by random variables

    \param para_cont_parameter (in):  parameter continuation parameter [0-1]
    */
    void ComputeMatParamsMultivectorFromRandomVariables(const double para_cont_parameter);

    /*!
    \brief Compute median material parameter multivector

    Compute multivector with ele col layout which contains the values
    of the median of the uncertain material parameters computed by random fields

    */
    void ComputeMedianMatParamsMultivectorFromRandomFields();

    /*!
    \brief create random field objects

    create random field objects for MC based on input file specifications

    \param random_field_id (in)    : Id of random variable which identifies the
                                     respective section in the input file

    \param myseed (in)             : seed for random number generator
    */
    Teuchos::RCP<RandomField> CreateRandomField(int random_field_id, unsigned int myseed);

    /*!
    \brief create random variables objects

    create random variables objects for MC based on input file specifications

    \param random_variable_id (in) : Id of random variable which identifies the
                                     respective section in the input file

    \param myseed (in)             : seed for random number generator
    */
    Teuchos::RCP<UQ::RandomVariable> CreateRandomVariable(
        int random_variable_id, unsigned int myseed);

   public:
   private:
    //! \brief discretization for convenience
    Teuchos::RCP<DRT::Discretization> discret_;

    //! \brief  map from stochparaid (int) to the random field which models this parameter
    std::map<int, Teuchos::RCP<RandomField>> randomfields_;

    //! \brief map from stochparaid (int) to the random variable which models this parameter
    std::map<int, Teuchos::RCP<RandomVariable>> randomvariables_;

    //! \brief number of stochastic material parameters modelled as random field
    int numstochparams_r_field_;

    //! \brief number of stochastic material parameters modelled as random variable
    int numstochparams_r_var_;

    /*!
    \brief map of names for each material parameter
     key: matid, value: value-enum
     value is of enum type and can be used to set the respective material
     parameter directly
     e.g. < 1 ; ( youngs, beta)>
          < 2 ; ( k1 , k2 )>
    */
    std::map<int, std::vector<int>> stochparamap_r_field_;

    /*!
    \brief map of matids to the id of the  stochastic parameter
    e.g.  < 1 ; (0 , 1)>
          < 2 ; (2)>
    thereby if say matid 1 represents AAAneohook with two stochastic
    parameters thereby the first is modelled by random field number
    0 and the second by random field number 1
    the second material has only one stochastic parameter which is modelled
    by field number 2
    */
    std::map<int, std::vector<int>> stochparaid_r_field_;

    /*!
    \brief map of names for each material parameter
    key: matid, value: value-enum
    value is of enum type and can be used to set the respective material
    parameter directly
     e.g. < 1 ; ( youngs, beta)>
          < 2 ; ( k1 , k2 )>

     */
    std::map<int, std::vector<int>> stochparamap_r_var_;

    /*!
    \brief map of matids to the id of the  stochastic parameter
     e.g.  < 1 ; (0 , 1)>
           < 2 ; (2)>
    thereby if say matid 1 represents AAAneohook with two stochastic
    parameters thereby the first is modelled by random field number
    0 and the second by random field number 1
    the second material has only one stochastic parameter which is modelled
    by field number 2

     */
    std::map<int, std::vector<int>> stochparaid_r_var_;

    //! actual set of material parameters in elementwise column map layout for random fields
    Teuchos::RCP<Epetra_MultiVector> params_r_field_;

    //! actual set of median material parameters in elementwise column map layout for random fields
    Teuchos::RCP<Epetra_MultiVector> params_r_field_median_;

    //! actual set of material parameters in for random variables. thus each vector has only one
    //! entry
    Teuchos::RCP<Epetra_MultiVector> params_r_var_;
  };
}  // namespace UQ

#endif  /* PAT_MATPAR_MANAGER_H */
#endif  // FFTW
