/*----------------------------------------------------------------------*/
/*! \file
\brief Bundle holds all read-in materials of a #GLOBAL::Problem

\level 1

*/

/*----------------------------------------------------------------------*/
/* macros */
#ifndef FOUR_C_MAT_PAR_BUNDLE_HPP
#define FOUR_C_MAT_PAR_BUNDLE_HPP

/*----------------------------------------------------------------------*/
/* headers */
#include "4C_config.hpp"

#include "4C_material_parameter_base.hpp"
#include "4C_utils_lazy_ptr.hpp"

#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN


namespace MAT
{
  namespace PAR
  {
    /*----------------------------------------------------------------------*/
    /// bundle holds all read-in materials of a #GLOBAL::Problem
    ///
    /// <h4>About</h4>
    /// The bundle provides an interface between unique material IDs and
    /// associated material parameters/data. The material ID is set via
    /// input file. It has to be unique, larger than zero or equal zero.
    /// Material ID and data are hold in #matmap_.
    ///
    /// <h4>Special issues for multi-problem-instance applications</h4>
    /// We have for each GLOBAL::Problem instance an individual material bundle.
    /// However, this fact is not transparanet at the read time of the elements
    /// which are only aware of their material ID. The variable #materialreadfromproblem_
    /// of the material bundle of the 0th GLOBAL::Problem instance make it possible to switch
    /// among different GLOBAL::Problem. (The variable #materialreadfromproblem_ is redundant
    /// in the material bundles of all non-0th GLOBAL::Problem instances.)
    ///
    /// \author bborn
    /// \date 02/09
    class Bundle
    {
     public:
      /**
       * Insert new pair of material ID and the input data. The input data is set up for lazy
       * construction the first time the material is accessed.
       */
      void insert(int matid, CORE::UTILS::LazyPtr<CORE::MAT::PAR::Parameter> mat);

      /**
       * Check whether material parameters exist for provided @p id.
       *
       * @note This call does not check whether material parameters are already constructed. It only
       * checks whether the material ID is known.
       */
      [[nodiscard]] bool id_exists(int id) const;

      /// provide access to material map (a li'l dirty)
      [[nodiscard]] const std::map<int, CORE::UTILS::LazyPtr<CORE::MAT::PAR::Parameter>>& Map()
          const
      {
        return matmap_;
      }

      /// return number of defined materials
      int Num() const { return matmap_.size(); }

      /// return material parameters
      CORE::MAT::PAR::Parameter* ParameterById(
          const int num  ///< request is made for this material ID
      ) const;

      /// return (first) ID by material type
      ///
      /// \return The ID of seached for material type.
      ///         If the search is unsuccessful -1 is returned
      int FirstIdByType(const CORE::Materials::MaterialType type) const;

      /// return problem index to read from
      int GetReadFromProblem() const { return materialreadfromproblem_; }

      /// set problem index to read from
      void SetReadFromProblem(const int p  ///< index of GLOBAL::Problem instance to read for
      )
      {
        materialreadfromproblem_ = p;
      }

      /// reset problem index to read from, i.e. to index 0
      void reset_read_from_problem() { materialreadfromproblem_ = 0; }

     private:
      /// The map linking material IDs to input paramters. The data is stored as a lazy pointer to
      /// allow for lazy construction of material parameters in arbitrary order.
      std::map<int, CORE::UTILS::LazyPtr<CORE::MAT::PAR::Parameter>> matmap_;

      /// the index of problem instance of which material read-in shall be performed
      int materialreadfromproblem_{};
    };

  }  // namespace PAR

}  // namespace MAT


FOUR_C_NAMESPACE_CLOSE

#endif
