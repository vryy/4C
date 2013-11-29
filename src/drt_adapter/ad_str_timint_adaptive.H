/*----------------------------------------------------------------------*/
/*!
\file ad_str_timint_adaptive.H

\brief Structure field adapter for time step size adaptivity

<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15237
</pre>
*/

/*----------------------------------------------------------------------*/

#ifndef ADAPTER_STRUCTURE_TIMINT_ADAPT_H
#define ADAPTER_STRUCTURE_TIMINT_ADAPT_H

/*----------------------------------------------------------------------*/
/* headers */
#include "ad_str_fsiwrapper.H"


// forward declarations:
namespace STR
{
  class TimInt;
  class TimAda;
}


/*----------------------------------------------------------------------*/
/* adapting adapter */
namespace ADAPTER {

  /*====================================================================*/
  /*!
   * \brief Adapter to general implicit structural time integration
   *
   * \date 08/08
   */
  class StructureTimIntAda : public FSIStructureWrapper
  {
  public:

    /// Constructor
    StructureTimIntAda(
      Teuchos::RCP<STR::TimAda> sta,
      Teuchos::RCP<STR::TimInt> sti
    );

    /// @name Time step helpers
    //@{

    /// Integrate from \f$t_1\f$ to \f$t_2\f$
    virtual int Integrate();

    /// prepare output (i.e. calculate stresses, strains, energies)
    virtual void PrepareOutput();

    /// output results
    virtual void Output();

    //@}


  private:

    /// the actual structure algorithm
    Teuchos::RCP<STR::TimAda> structure_;

  };  // class StructureTimIntAda

}  // namespace ADAPTER

/*----------------------------------------------------------------------*/
#endif
