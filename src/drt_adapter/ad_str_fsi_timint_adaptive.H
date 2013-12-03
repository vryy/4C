/*----------------------------------------------------------------------*/
/*!
\file ad_str_fsi_timint_adaptive.H

\brief Structure field adapter for time step size adaptivity within monolithic FSI

<pre>
Maintainer: Matthias Mayr
            mayr@lnm.mw.tum.de
            http://www.mhpc.mw.tum.de
            089 - 289-10362
</pre>
*/

/*----------------------------------------------------------------------*/

#ifndef ADAPTER_STRUCTURE_FSI_TIMINT_ADAPT_H
#define ADAPTER_STRUCTURE_FSI_TIMINT_ADAPT_H

/*----------------------------------------------------------------------*/
/* headers */
#include "ad_str_fsiwrapper.H"
#include "ad_str_timint_adaptive.H"

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
   * \brief Structure field adapter for time step size adaptivity within monolithic FSI
   *
   * Use this adapter in case you want to do monolithic FSI with time step size
   * adaptivity. By inheritance, we combine FSI functionalities with structural
   * time adaptivity. The FSI stuff is inherited from ADAPTER::FSIStructureWrapper
   * and the time adaptivity from ADAPTER::StructureTimIntAda
   *
   * The time loop is implemented in FSI::Monolithic, which requires
   * error estimation and time step size calculation based on the structure field.
   * For error estimation and time step size calculation we want to use the standard
   * structural time adaptivity routines. Though, the decision, whether a time step
   * has to be repeated, has to be made by the FSI algorithm.
   *
   * \sa FSIStructureWrapper
   * \sa StructureTimIntAda
   *
   * \author mayr.mt
   * \date 12/2013
   */
  class StructureFSITimIntAda : virtual public FSIStructureWrapper,
                                virtual public StructureTimIntAda
  {
  public:

    //! Constructor
    StructureFSITimIntAda(
      Teuchos::RCP<STR::TimAda> sta,
      Teuchos::RCP<Structure> sti
    );

    //! Destructor
    virtual ~StructureFSITimIntAda() { ; }

  };  // class StructureTimIntAda

}  // namespace ADAPTER

/*----------------------------------------------------------------------*/
#endif
