/*! \file

\level 1

\brief Common service function for prestress


*/

#ifndef PRESTRESS_SERVICE_H_
#define PRESTRESS_SERVICE_H_

#include <Teuchos_ParameterList.hpp>
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_inpar/inpar_structure.H"

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
    static inline double GetTime()
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
     * \brief Returns whether Inverse Design is set for prestressing in the parameters of
     * structural dynamics. This method does not ensure that Inverse Design is actually
     * active
     *
     * \return true Inverse Design is set in input file
     * \return false Inverse Design is not set in input file
     */
    static inline bool IsInverseDesign() { return GetType() == INPAR::STR::PreStress::id; }

    /*!
     * \brief Returns whether Inverse Design is set for prestressing in the given
     * parameter. This method does not ensure that Inverse Design is actually active
     *
     * \param pstype Prestress type that is used
     * \return true Inverse Design is set in input file
     * \return false Inverse Design is not set in input file
     */
    static inline bool IsInverseDesign(INPAR::STR::PreStress pstype)
    {
      return pstype == INPAR::STR::PreStress::id;
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
    static inline bool IsAny() { return GetType() != INPAR::STR::PreStress::none; }

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
      return IsAny() && currentTime <= GetTime() + 1.0e-15;
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
      return IsAny(pstype) && currentTime <= pstime + 1.0e-15;
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
      return IsMulf() && currentTime <= GetTime() + 1.0e-15;
    }

    /*!
     * \brief Returns true if material iterative prestressing method is currently active with the
     * parameters of strtuctural dynamics.
     *
     * \param currentTimeCurrent time of the simulation
     * \return true Material iterative prestressing method is active
     * \return false Material iterative prestressing method is active
     */
    static inline bool IsMaterialIterativeActive(const double currentTime)
    {
      return IsMaterialIterative() && currentTime <= GetTime() + 1.0e-15;
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

    /*!
     * \brief Returns true if Inverse Design prestressing method is currently active with the
     * parameters of strtuctural dynamics.
     *
     * \param currentTimeCurrent time of the simulation
     * \return true Inverse Design prestressing method is active
     * \return false Inverse Design prestressing method is active
     */
    static inline bool IsInverseDesignActive(const double currentTime)
    {
      return IsInverseDesign() && currentTime <= GetTime() + 1.0e-15;
    }

    /*!
     * \brief Returns true if Inverse Design prestressing method is currently active with the given
     * parameters.
     *
     * \param currentTimeCurrent time of the simulation
     * \param pstype Prestress type that is used
     * \param pstime Prestress time that is used
     * \return true Inverse Design prestressing method is active
     * \return false Inverse Design prestressing method is active
     */
    static inline bool IsInverseDesignActive(
        const double currentTime, INPAR::STR::PreStress pstype, const double pstime)
    {
      return IsInverseDesign(pstype) && currentTime <= pstime + 1.0e-15;
    }

  }  // namespace PRESTRESS
}  // namespace UTILS

#endif  // PRESTRESS_SERVICE_H_