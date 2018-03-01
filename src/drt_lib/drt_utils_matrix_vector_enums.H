/*---------------------------------------------------------------------*/
/*!
\file drt_utils_matrix_vector_enums.H

\brief Contains ONLY lists of enumerators and is supposed to be included
       into the header files, if necessary and till the C++11 standard is
       available in BACI.

\maintainer Michael Hiermeier

\date Apr 21, 2016

\level 3

*/
/*---------------------------------------------------------------------*/

#ifndef SRC_DRT_LIB_DRT_UTILS_MATRIX_VECTOR_ENUMS_H_
#define SRC_DRT_LIB_DRT_UTILS_MATRIX_VECTOR_ENUMS_H_

namespace DRT {
namespace UTILS {
  /// enum, which specifies the desired global vector initialization during creation
  enum VecInitType
  {
    vec_init_zero,                ///< fill the vector with zeros
    vec_init_last_time_step,      ///< use the last converged time step state
    vec_init_current_state        ///< use the current state
  };

  /// enum, which specifies the desired matrix block for the different models
  enum MatBlockType
  {
    block_displ_displ,   ///< Kdd block (structural block)
    block_displ_lm,      ///< Kdz block (of the corresponding model evaluator)
    block_lm_displ,      ///< Kzd block (of the corresponding model evaluator)
    block_lm_lm,         ///< Kzz block (of the corresponding model evaluator)
    block_temp_temp,     ///< Ktt block (thermal block)
    block_temp_displ,    ///< Ktd block (structure-thermo-coupling)
    block_displ_temp,    ///< Kdt block (thermo-structure-coupling)
  };

  //! map matrix block type enum to std::string
  inline const std::string MatBlockType2String(const enum MatBlockType type)
  {
    switch (type)
    {
      case block_displ_displ:
        return "block_displ_displ";
      case block_displ_lm:
        return "block_displ_lm";
      case block_lm_displ:
        return "block_lm_displ";
      case block_lm_lm:
        return "block_lm_lm";
      case block_temp_temp:
        return "block_temp_temp";
      case block_temp_displ:
        return "block_temp_displ";
      case block_displ_temp:
        return "block_displ_temp";
      default:
        return "unknown matrix block type";
    }
  }

  /** enum, which specifies the desired vector blocks for the different models */
  enum VecBlockType
  {
    block_displ,      ///< displacement block (structural block)
    block_constraint, ///< lagrange multiplier/constraint block of the corresponding model
    block_temp,       ///< temperature block (thermal block)
  };

  //! map vector block type enum to std::string
  inline const std::string VecBlockType2String( const enum VecBlockType type )
  {
    switch (type)
    {
      case block_displ:
        return "block_displ";
      case block_constraint:
        return "block_constraint";
      case block_temp:
        return "block_temp";
      default:
        return "unknown vector block type";
    }
  }

} // namespace UTILS
} // namespace DRT


#endif /* SRC_DRT_LIB_DRT_UTILS_MATRIX_VECTOR_ENUMS_H_ */
