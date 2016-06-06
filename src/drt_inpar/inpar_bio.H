/*----------------------------------------------------------------------*/
/*!
\file inpar_bio.H

\brief Input parameters for biomedical simulations

\level 3

\maintainer Christian Roth

*/
/*----------------------------------------------------------------------*/
#ifndef INPAR_ARTNET_H
#define INPAR_ARTNET_H

#include "inpar_parameterlist_utils.H"

// forward declaration
namespace DRT { namespace INPUT { class ConditionDefinition; } }

/*----------------------------------------------------------------------*/

// ToDo: move these enums to namespace INPAR::ARTNET etc.
//       is the typedef really needed?

/*!----------------------------------------------------------------------
\brief enum of arterial network dynamic types

<pre>                                                      ismail 11/08
This is the enumeration of all types of different integration schemes
</pre>

*-----------------------------------------------------------------------*/
typedef enum _ARTNET_DYNTYPE
{
  typ_tay_gal=0
} _ARTNET_DYNTYPE;


/*!----------------------------------------------------------------------
\brief enum of reduced dimensional airways dynamic types

<pre>                                                     ismail 01/10
This is the enumeration of all types of different integration schemes
</pre>

*-----------------------------------------------------------------------*/
typedef enum _RED_AIRWAYS_DYNTYPE
{
  one_step_theta,
  linear,
  nonlinear
} _RED_AIRWAYS_DYNTYPE;







namespace INPAR
{
  namespace ARTDYN
  {
    /// set the arterial dynamic parameters
    void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);
  }

  namespace ARTNET
  {
    /*!----------------------------------------------------------------------
    \brief enum of reduced dimensional relaxation type

    <pre>                                                     roth 10/13
    This is the enumeration of all types of different relaxation types
    </pre>

    *-----------------------------------------------------------------------*/
    enum RELAXTYPE_3D_0D
    {
      norelaxation,
      fixedrelaxation,
      Aitken,
      SD
    };

    /// set the artnet parameters
    void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

    /// set specific artnet conditions
    void SetValidConditions(std::vector<Teuchos::RCP<DRT::INPUT::ConditionDefinition> >& condlist);

  }

  namespace BIOFILM
  {
    /// set the biofilm parameters
    void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

    /// set specific biofilm conditions
    void SetValidConditions(std::vector<Teuchos::RCP<DRT::INPUT::ConditionDefinition> >& condlist);
  }

  namespace PATSPEC
  {
    /// set the patient specific parameters
    void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

    /// set specific patspec conditions (AAA)
    void SetValidConditions(std::vector<Teuchos::RCP<DRT::INPUT::ConditionDefinition> >& condlist);
  }

  namespace REDAIRWAYS
  {
    /// set the reduced airways parameters
    void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

    /// set specific reduced airways conditions
    void SetValidConditions(std::vector<Teuchos::RCP<DRT::INPUT::ConditionDefinition> >& condlist);
  }
}  // namespace INPAR

/*----------------------------------------------------------------------*/
#endif  // #ifndef INPAR_ARTNET_H
