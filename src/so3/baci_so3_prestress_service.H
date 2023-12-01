/*! \file

\level 1

\brief Common service function for prestress


*/

#ifndef BACI_SO3_PRESTRESS_SERVICE_H
#define BACI_SO3_PRESTRESS_SERVICE_H


#include "baci_inpar_structure.H"
#include "baci_lib_globalproblem.H"

#include <Teuchos_ParameterList.hpp>

namespace UTILS
{
  namespace PRESTRESS
  {
    /*!
     * \brief Returns the type of the prestress algorithm stored in the parameters of structural
     * dynamics
     *
     * \return INPAR::STR::PreStress
     */
    static inline INPAR::STR::PreStress GetType()
    {
      static INPAR::STR::PreStress pstype = Teuchos::getIntegralValue<INPAR::STR::PreStress>(
          DRT::Problem::Instance()->StructuralDynamicParams(), "PRESTRESS");

      return pstype;
    }

    /*!
     * \brief Returns the prestress time stored in the parameters of structural dynamics
     *
     * \return double
     */
    static inline double GetPrestressTime()
    {
      static double pstime =
          DRT::Problem::Instance()->StructuralDynamicParams().get<double>("PRESTRESSTIME");

      return pstime;
    }

    /*!
     * \brief Returns whether MULF is set for prestressing in the parameters of structural dynamics.
     * This method does not ensure that MULF is actually active
     *
     * \return true MULF is set in input file
     * \return false MULF is not set in input file
     */
    static inline bool IsMulf() { return GetType() == INPAR::STR::PreStress::mulf; }

    /*!
     * \brief Returns whether material iterative prestressing is set in the parameters of structural
     * dynamics. This method does not ensure that prestressing is actually active
     *
     * \return true material iterative prestressing is set in input file
     * \return false material iterative prestressing is not set in input file
     */
    static inline bool IsMaterialIterative()
    {
      return GetType() == INPAR::STR::PreStress::material_iterative;
    }

    /*!
     * \brief Returns whether MULF is set for prestressing as the given prestress type.
     * This method does not ensure that MULF is actually active
     *
     * \param pstype Prestress type that is used
     * \return true MULF is set in input file
     * \return false MULF is not set in input file
     */
    static inline bool IsMulf(INPAR::STR::PreStress pstype)
    {
      return pstype == INPAR::STR::PreStress::mulf;
    }

    /*!
     * \brief Returns whether material iterative prestressing is set as the given prestress type.
     * This method does not ensure that prestressing is actually active
     *
     * \param pstype Prestress type that is used
     * \return true material iterative prestressing is set in input file
     * \return false material iterative prestressing is not set in input file
     */
    static inline bool IsMaterialIterative(INPAR::STR::PreStress pstype)
    {
      return pstype == INPAR::STR::PreStress::material_iterative;
    }


    /*!
     * \brief Returns whether no prestressing is set in the parameters of
     * structural dynamics.
     *
     * \return true No prestressing is set in the input file
     * \return false Prestressing is set in the input file
     */
    static inline bool IsNone() { return GetType() == INPAR::STR::PreStress::none; }


    /*!
     * \brief Returns whether no prestressing is set in the given parameter.
     *
     * \param pstype Prestress type that is used
     * \return true No prestressing is set in the input parameter
     * \return false Prestressing is set in the input parameter
     */
    static inline bool IsNone(INPAR::STR::PreStress pstype)
    {
      return pstype == INPAR::STR::PreStress::none;
    }

    /*!
     * \brief Returns whether any prestressing is set in the parameters of
     * structural dynamics.
     *
     * \return true Prestressing is set in the input file
     * \return false No prestressing is set in the input file
     */
    static inline bool IsAny()
    {
      return Teuchos::getIntegralValue<INPAR::STR::PreStress>(
                 DRT::Problem::Instance()->StructuralDynamicParams(), "PRESTRESS") !=
             INPAR::STR::PreStress::none;
    }

    /*!
     * \brief Returns whether prestressing is set in the given parameter.
     *
     * \param pstype Prestress type that is used
     * \return true Prestressing is set in the input parameter
     * \return false No prestressing is set in the input parameter
     */
    static inline bool IsAny(INPAR::STR::PreStress pstype)
    {
      return pstype != INPAR::STR::PreStress::none;
    }

    /*!
     * \brief Returns true if any prestressing method is currently active with the parameters of
     * strtuctural dynamics.
     *
     * \param currentTime Current time of the simulation
     * \return true Any prestressing method is active
     * \return false No prestressing method is active
     */
    static inline bool IsActive(const double currentTime)
    {
      INPAR::STR::PreStress pstype = Teuchos::getIntegralValue<INPAR::STR::PreStress>(
          DRT::Problem::Instance()->StructuralDynamicParams(), "PRESTRESS");
      const double pstime =
          DRT::Problem::Instance()->StructuralDynamicParams().get<double>("PRESTRESSTIME");
      return pstype != INPAR::STR::PreStress::none && currentTime <= pstime + 1.0e-15;
    }

    /*!
     * \brief Returns true if any prestressing method is currently active with the given parameters.
     *
     * \param currentTimeCurrent time of the simulation
     * \param pstype Prestress type that is used
     * \param pstime Prestress time that is used
     * \return true Any prestressing method is active
     * \return false No prestressing method is active
     */
    static inline bool IsActive(
        const double currentTime, INPAR::STR::PreStress pstype, const double pstime)
    {
      return pstype != INPAR::STR::PreStress::none && currentTime <= pstime + 1.0e-15;
    }

    /*!
     * \brief Returns true if MULF prestressing method is currently active with the parameters of
     * strtuctural dynamics.
     *
     * \param currentTimeCurrent time of the simulation
     * \return true MULF prestressing method is active
     * \return false MULF prestressing method is active
     */
    static inline bool IsMulfActive(const double currentTime)
    {
      return IsMulf() && currentTime <= GetPrestressTime() + 1.0e-15;
    }

    /*!
     * \brief Returns true if MULF prestressing method is currently active with the given
     * parameters.
     *
     * \param currentTimeCurrent time of the simulation
     * \param pstype Prestress type that is used
     * \param pstime Prestress time that is used
     * \return true MULF prestressing method is active
     * \return false MULF prestressing method is active
     */
    static inline bool IsMulfActive(
        const double currentTime, INPAR::STR::PreStress pstype, const double pstime)
    {
      return IsMulf(pstype) && currentTime <= pstime + 1.0e-15;
    }


  }  // namespace PRESTRESS
}  // namespace UTILS

#endif