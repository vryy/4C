/*----------------------------------------------------------------------*/
/*! \file

\brief Structure field adapter for time step size adaptivity

\level 2


*/

/*----------------------------------------------------------------------*/

#ifndef FOUR_C_ADAPTER_STR_TIMINT_ADAPTIVE_HPP
#define FOUR_C_ADAPTER_STR_TIMINT_ADAPTIVE_HPP

/*----------------------------------------------------------------------*/
/* headers */
#include "4C_config.hpp"

#include "4C_adapter_str_wrapper.hpp"

FOUR_C_NAMESPACE_OPEN

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
   * \brief Adapter to adaptive implicit structural time integration
   *
   * \date 08/08
   */
  class StructureTimIntAda : public StructureWrapper
  {
   public:
    /// Constructor
    StructureTimIntAda(Teuchos::RCP<STR::TimAda> sta, Teuchos::RCP<Structure> sti);

    /// @name Time step helpers
    //@{

    /// Integrate from \f$t_1\f$ to \f$t_2\f$
    int Integrate() override;

    /// prepare output (i.e. calculate stresses, strains, energies)
    void prepare_output(bool force_prepare) override;

    /// output results
    virtual void Output();

    //@}

   protected:
    //! Access routines
    //{@

    Teuchos::RCP<STR::TimAda> StrAda() const { return sta_; }

    //@}

   private:
    /// the actual structure algorithm
    Teuchos::RCP<STR::TimAda> sta_;  // STR::TimAda is the old time integration

  };  // class StructureTimIntAda

}  // namespace ADAPTER

/*----------------------------------------------------------------------*/
FOUR_C_NAMESPACE_CLOSE

#endif
