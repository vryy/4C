/*---------------------------------------------------------------------*/
/*! \file

\brief Contains ONLY lists of enumerators and is supposed to be included
       into the header files, if necessary and till the C++11 standard is
       available in BACI.


\level 3

*/
/*---------------------------------------------------------------------*/

#ifndef SRC_DRT_LIB_DRT_UTILS_MATRIX_VECTOR_ENUMS_H_
#define SRC_DRT_LIB_DRT_UTILS_MATRIX_VECTOR_ENUMS_H_

namespace DRT
{
  namespace UTILS
  {
    /// enum, which specifies the desired global vector initialization during creation
    enum class VecInitType
    {
      zero,               ///< fill the vector with zeros
      last_time_step,     ///< use the last converged time step state
      init_current_state  ///< use the current state
    };

    /// enum, which specifies the desired matrix block for the different models
    enum class MatBlockType
    {
      displ_displ,          ///< Kdd block (structural block)
      displ_lm,             ///< Kdz block (of the corresponding model evaluator)
      lm_displ,             ///< Kzd block (of the corresponding model evaluator)
      lm_lm,                ///< Kzz block (of the corresponding model evaluator)
      temp_temp,            ///< Ktt block (thermal block)
      temp_displ,           ///< Ktd block (structure-thermo-coupling)
      displ_temp,           ///< Kdt block (thermo-structure-coupling)
      porofluid_porofluid,  ///< Kpp block (porofluid-porofluid)
      porofluid_displ,      ///< Kpd block (porofluid-structure)
      displ_porofluid       ///< Kdp block (structure-porofluid)
    };

    //! map matrix block type enum to std::string
    inline const std::string MatBlockType2String(const enum DRT::UTILS::MatBlockType type)
    {
      switch (type)
      {
        case DRT::UTILS::MatBlockType::displ_displ:
          return "block_displ_displ";
        case DRT::UTILS::MatBlockType::displ_lm:
          return "block_displ_lm";
        case DRT::UTILS::MatBlockType::lm_displ:
          return "block_lm_displ";
        case DRT::UTILS::MatBlockType::lm_lm:
          return "block_lm_lm";
        case DRT::UTILS::MatBlockType::temp_temp:
          return "block_temp_temp";
        case DRT::UTILS::MatBlockType::temp_displ:
          return "block_temp_displ";
        case DRT::UTILS::MatBlockType::displ_temp:
          return "block_displ_temp";
        case DRT::UTILS::MatBlockType::porofluid_porofluid:
          return "block_porofluid_porofluid";
        case DRT::UTILS::MatBlockType::porofluid_displ:
          return "block_porofluid_displ";
        case DRT::UTILS::MatBlockType::displ_porofluid:
          return "block_displ_porofluid";
        default:
          return "unknown matrix block type";
      }
    }

    //! enum, which specifies the desired vector blocks for the different models
    enum class VecBlockType
    {
      displ,       ///< displacement block (structural block)
      constraint,  ///< lagrange multiplier/constraint block of the corresponding model
      temp,        ///< temperature block (thermal block)
      porofluid,   ///< porofluid block (porofluid block)
    };

    //! map vector block type enum to std::string
    inline const std::string VecBlockType2String(const enum VecBlockType type)
    {
      switch (type)
      {
        case DRT::UTILS::VecBlockType::displ:
          return "block_displ";
        case DRT::UTILS::VecBlockType::constraint:
          return "block_constraint";
        case DRT::UTILS::VecBlockType::temp:
          return "block_temp";
        case DRT::UTILS::VecBlockType::porofluid:
          return "block_porofluid";
        default:
          return "unknown vector block type";
      }
    }

  }  // namespace UTILS
}  // namespace DRT


#endif /* SRC_DRT_LIB_DRT_UTILS_MATRIX_VECTOR_ENUMS_H_ */
