/*----------------------------------------------------------------------*/
/*! \file

\brief Structure field adapter for time step size adaptivity

\level 2


*/

/*----------------------------------------------------------------------*/

#ifndef ADAPTER_STR_TIMINT_ADAPTIVE_H
#define ADAPTER_STR_TIMINT_ADAPTIVE_H

/*----------------------------------------------------------------------*/
/* headers */
#include "adapter_str_wrapper.H"

// forward declarations:
namespace STR
{
  class TimInt;
  class TimAda;
}  // namespace STR

/*----------------------------------------------------------------------*/
/* adapting adapter */
namespace ADAPTER
{
  /*====================================================================*/
  /*!
   * \brief Adapter to general implicit structural time integration
   *
   * \date 08/08
   */
  class StructureTimIntAda : public StructureWrapper
  {
   public:
    /// Constructor
    StructureTimIntAda(Teuchos::RCP<STR::TimAda> sta, Teuchos::RCP<Structure> sti);

    /// Destructor
    virtual ~StructureTimIntAda() { ; }

    /// @name Time step helpers
    //@{

    /// Integrate from \f$t_1\f$ to \f$t_2\f$
    virtual int Integrate();

    /// prepare output (i.e. calculate stresses, strains, energies)
    void PrepareOutput(bool force_prepare) override;

    /// output results
    virtual void Output();

    //@}

   protected:
    //! Access routines
    //{@

    Teuchos::RCP<STR::TimAda> StrAda() const { return structure_; }

    //@}

   private:
    /// the actual structure algorithm
    Teuchos::RCP<STR::TimAda> structure_;

  };  // class StructureTimIntAda

}  // namespace ADAPTER

/*----------------------------------------------------------------------*/
#endif
