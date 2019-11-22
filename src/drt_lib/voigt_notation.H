/*! \file
\level 1
\brief Voigt notation definition and utilities
\maintainer Sebastian Proell
*/

#ifndef VOIGT_NOTATION_H
#define VOIGT_NOTATION_H

#include "../linalg/linalg_fixedsizematrix.H"


namespace UTILS
{
  namespace VOIGT
  {
    /*! Voigt notation types
     *
     * This enum can be used whenever the distinction between stress- and strain-like Voigt notation
     * is important. The classical Voigt notation is only meaningful for symmetric tensors.
     */
    enum class NotationType
    {
      /// stress-like Voigt notation
      ///
      /// A vector in stress-like Voigt notation, contains the off-diagonal values only once.
      stress,
      /// strain-like Voigt notation
      ///
      /// A vector in strain-like Voigt notation, contains the sum of corresponding off-diagonal
      /// values.
      strain
    };

    /*--------------------------------------------------------------------------*/
    /** \brief Utility routines for the perturbed Voigt tensor notation
     *
     * \tparam type specific NotationType this class operates on
     */
    template <NotationType type>
    class VoigtUtils
    {
     public:
      /// instantiation is forbidden
      VoigtUtils() = delete;

      /** \brief compute power of a symmetric 3x3 matrix in perturbed Voigt notation
       *
       *  \f[
       *  [vtensorpow]_{AB} = [vtensor^{pow}]_{AB}
       *  \f]
       *
       *  \param[in] pow          positive integer exponent
       *  \param[in] vtensor      input tensor in Voigt <type> notation
       *  \param[in] vtensor_pow  result, i.e. input tensor to the given power
       *
       *  \author hiermeier \date 09/17 */
      static void PowerOfSymmetricTensor(
          unsigned pow, const LINALG::Matrix<6, 1>& vtensor, LINALG::Matrix<6, 1>& vtensor_pow);

      /** \brief Compute the inverse tensor in perturbed Voigt notation
       *
       *  \param[in]  vtensor      tensor in Voigt <type> notation
       *  \param[out] vtensor_inv  inverse tensor in Voigt <type> notation
       *
       *  \author hiermeier \date 09/17 */
      static void InverseTensor(
          const LINALG::Matrix<6, 1>& vtensor, LINALG::Matrix<6, 1>& vtensor_inv);

      /**
       * \brief Compute the determinant of a matrix in Voigt <type> notation
       *
       * @param vtensor Tensor in Voigt <type> notation
       */
      static inline double Determinant(const LINALG::Matrix<6, 1>& vtensor)
      {
        return TripleEntryProduct(vtensor, 0, 1, 2) + 2 * TripleEntryProduct(vtensor, 3, 4, 5) -
               TripleEntryProduct(vtensor, 1, 5, 5) - TripleEntryProduct(vtensor, 2, 3, 3) -
               TripleEntryProduct(vtensor, 0, 4, 4);
      }

      /** \brief Compute the three principal invariants of a matrix in Voigt <type> notation
       *
       * @param[out] prinv the three principal invariants
       * @param vtensor tensor in Voigt <type> notation
       */
      static inline void InvariantsPrincipal(
          LINALG::Matrix<3, 1>& prinv, const LINALG::Matrix<6, 1>& vtensor)
      {
        // 1st invariant, trace tens
        prinv(0) = vtensor(0) + vtensor(1) + vtensor(2);
        // 2nd invariant, 0.5( (trace(tens))^2 - trace(tens^2))
        prinv(1) = 0.5 * (prinv(0) * prinv(0) - vtensor(0) * vtensor(0) - vtensor(1) * vtensor(1) -
                             vtensor(2) * vtensor(2)) -
                   vtensor(3) * vtensor(3) * UnscaleFactor(3) * UnscaleFactor(3) -
                   vtensor(4) * vtensor(4) * UnscaleFactor(4) * UnscaleFactor(4) -
                   vtensor(5) * vtensor(5) * UnscaleFactor(5) * UnscaleFactor(5);
        // 3rd invariant, determinant tens
        prinv(2) = Determinant(vtensor);
      }

      /** \brief Compute the product of a tensor in perturbed Voigt notation
       *  and a vector
       *
       *  \f$ [vecres]_{A} = vtensor_{AB} vec^{B} \f$
       *
       *  \param[in]  vtensor      tensor in Voigt <type> notation
       *  \param[out] vtensor_inv  inverser tensor in Voigt <type> notation
       *
       *  \author hiermeier \date 09/17 */
      static void MultiplyTensorVector(const LINALG::Matrix<6, 1>& vtensor,
          const LINALG::Matrix<3, 1>& vec, LINALG::Matrix<3, 1>& vec_res);

      /** \brief Compute the symmetric outer product of two vectors
       *
       *  \f$ [abba]_{AB} = [veca]_{A} [vecb]_{B} + [veca]_{B} [vecb]_{A} \f$
       *
       *  \param[in]  vec_a  first vector
       *  \param[in]  vec_b  second vector
       *  \param[out] ab_ba  symmetric outer product of the two input vectors
       *                     in the Voigt <type> notation
       *
       *  \author hiermeier \date 09/17 */
      static void SymmetricOuterProduct(const LINALG::Matrix<3, 1>& vec_a,
          const LINALG::Matrix<3, 1>& vec_b, LINALG::Matrix<6, 1>& ab_ba);

      /*!
       * Converts a <type>-like tensor to stress-like Voigt notation
       * @param vtensor_in tensor in <type>-like Voigt notation
       * @param vtensor_out tensor in stress-like Voigt notation
       */
      static void ToStressLike(
          const LINALG::Matrix<6, 1>& vtensor_in, LINALG::Matrix<6, 1>& vtensor_out);

      /*!
       * Converts a <type>-like tensor to strain-like Voigt notation
       * @param vtensor_in tensor in <type>-like Voigt notation
       * @param vtensor_out tensor in strain-like Voigt notation
       */
      static void ToStrainLike(
          const LINALG::Matrix<6, 1>& vtensor_in, LINALG::Matrix<6, 1>& vtensor_out);

      /*!
       * Converts a <type>-like tensor in Voigt notation to a matrix
       * @param vtensor_in tensor in <type>-like Voigt notation
       * @param vtensor_out tensor as a matrix
       */
      static void ToMatrix(
          const LINALG::Matrix<6, 1>& vtensor_in, LINALG::Matrix<3, 3>& tensor_out);


      /*! Copy matrix contents to <type>-like Voigt notation
       *
       *          A_00 A_01 A_02
       * Matrix   A_10 A_11 A_12   is converted to the Voigt vector below:
       *          A_20 A_21 A_22
       *
       * Vector   V_0 = A_00; V_1 = A_11; V_2 = A_22; V_3 = 0.5*k*(A_01 + A_10);
       *          V_4 = 0.5*k*(A_12 + A_21); V_5 = 0.5*k*(A_02 + A_20)
       *
       * where k is the scale factor for <type>-like Voigt notation.
       *
       * @param tensor_in the matrix to copy from
       * @param[out] vtensor_out target tensor in <type>-like Voigt notation
       */
      template <typename T>
      static void MatrixToVector(
          const LINALG::Matrix<3, 3, T>& tensor_in, LINALG::Matrix<6, 1, T>& vtensor_out);

      /// access scaling factors
      static inline double ScaleFactor(unsigned i) { return scale_fac_[i]; };

      /// access unscaling factors
      static inline double UnscaleFactor(unsigned i) { return unscale_fac_[i]; };

     private:
      /** \brief scale off diagonal values
       *
       *  \note This function changes the values only if the strain notation is used.
       *
       *  \param[out] tensor  scale the off-diagonal values of this tensor
       *
       *  \author hiermeier \date 09/17 */
      static void ScaleOffDiagonalVals(LINALG::Matrix<6, 1>& tensor);

      /** \brief unscale off diagonal values
       *
       *  \note This function changes the values only if the strain notation is used.
       *
       *  \param[out] tensor  unscale the off-diagonal values of this tensor
       *
       *  \author hiermeier \date 09/17 */
      static void UnscaleOffDiagonalVals(LINALG::Matrix<6, 1>& strain);


      /** \brief unscale factors for the perturbed Voigt strain notation
       *
       *  These factors are meaningful if the strain convention is followed. */
      static const double unscale_fac_[6];

      /** \brief scale factors for the perturbed Voigt stress notation
       *
       *  These factors are meaningful if the strain convention is followed. */
      static const double scale_fac_[6];

      /**
       * \brief Compute the product of three matrix entries
       *
       * The entries are correctly scaled depending on the VoigtNotation type of the tensor.
       * @param vtensor the tensor in voigt notaion
       * @param i first entry's Voigt index
       * @param j second entry's Voigt inde
       * @param k third entry's Voigt index
       * @return product of the three entries
       */
      static inline double TripleEntryProduct(
          const LINALG::Matrix<6, 1>& vtensor, unsigned i, unsigned j, unsigned k)
      {
        return vtensor(i) * UnscaleFactor(i) * vtensor(j) * UnscaleFactor(j) * vtensor(k) *
               UnscaleFactor(k);
      }
    };


    /// convert non-symmetric 2-tensor to 9x1 vector
    ///          A_00 A_01 A_02
    /// Matrix   A_10 A_11 A_12  is converted to the vector below:
    ///          A_20 A_21 A_22
    ///
    /// Vector   V_0 = A_00; V_1 = A_11; V_2 = A_22; V_3 = A_01; V_4 = A_12; V_5 = A_02; V_6 = A_10;
    /// V_7 = A_21; V_8 = A_20
    void Matrix3x3to9x1(LINALG::Matrix<3, 3> const& in, LINALG::Matrix<9, 1>& out);

    /**
     * \brief Identity matrix in stress/strain-like Voigt notation
     * @param id (out) : 2nd order identity tensor in stress/strain-like Voigt notation
     */
    inline void IdentityMatrix(LINALG::Matrix<6, 1>& id)
    {
      id.Clear();
      for (unsigned i = 0; i < 3; ++i) id(i) = 1.0;
    }

    /*!
     * \brief Constructs a 4th order identity matrix with rows in <rows_notation>-type Voigt
     * notation and columns in <cols_notation>-type Voigt notation
     * @tparam rows_notation Voigt notation used for the rows
     * @tparam cols_notation  Voigt notation used for the columns
     * @param id (out) : Voigt-Matrix
     */
    template <NotationType rows_notation, NotationType cols_notation>
    void FourthOrderIdentityMatrix(LINALG::Matrix<6, 6>& id);

    /// collection of index mappings from matrix to Voigt-notation or vice versa
    struct IndexMappings
    {
     public:
      /**
       * from 6-Voigt index to corresponding 2-tensor row
       * @param i the index of a 6x1 vector in Voigt notation
       * @return row index of the corresponding 3x3 matrix
       */
      static inline int Voigt6ToRow(unsigned int i)
      {
        assertRangeVoigtIndex(i);
        static constexpr int VOIGT6ROW[6] = {0, 1, 2, 0, 1, 2};
        return VOIGT6ROW[i];
      };

      /// from 6-Voigt index to corresponding 2-tensor col index
      static inline int Voigt6ToCol(unsigned int i)
      {
        assertRangeVoigtIndex(i);
        static constexpr int VOIGT6COL[6] = {0, 1, 2, 1, 2, 0};
        return VOIGT6COL[i];
      };

      /// from symmetric 2-tensor index pair to 6-Voigt index
      // TODO doesn't this need to check if strain or stress notation?
      static inline int SymToVoigt6(unsigned int row, unsigned int col)
      {
        assertRangeMatrixIndex(row, col);
        static constexpr int VOIGT3X3SYM[3][3] = {{0, 3, 5}, {3, 1, 4}, {5, 4, 2}};
        return VOIGT3X3SYM[row][col];
      };

      /// from non-symmetric 2-tensor index pair to 9-Voigt index
      static inline int NonSymToVoigt9(unsigned int row, unsigned int col)
      {
        assertRangeMatrixIndex(row, col);
        static constexpr int VOIGT3X3NONSYM[3][3] = {{0, 3, 5}, {6, 1, 4}, {8, 7, 2}};
        return VOIGT3X3NONSYM[row][col];
      };

      /// from symmetric 6x6 Voigt notation matrix indices (e.g. constitutive tensor) to one of the
      /// four indices of a four tensor
      static inline int Voigt6x6To4Tensor(
          unsigned int voigt_row, unsigned int voigt_col, unsigned int target_index)
      {
        assertRangeVoigtIndex(voigt_row);
        assertRangeVoigtIndex(voigt_col);
        dsassert(target_index < 4, "target index for fourth order tensor out of range");
        static constexpr int FOURTH[6][6][4] = {
            {{0, 0, 0, 0}, {0, 0, 1, 1}, {0, 0, 2, 2}, {0, 0, 0, 1}, {0, 0, 1, 2}, {0, 0, 0, 2}},
            {{1, 1, 0, 0}, {1, 1, 1, 1}, {1, 1, 2, 2}, {1, 1, 0, 1}, {1, 1, 1, 2}, {1, 1, 0, 2}},
            {{2, 2, 0, 0}, {2, 2, 1, 1}, {2, 2, 2, 2}, {2, 2, 0, 1}, {2, 2, 1, 2}, {2, 2, 0, 2}},
            {{0, 1, 0, 0}, {0, 1, 1, 1}, {0, 1, 2, 2}, {0, 1, 0, 1}, {0, 1, 1, 2}, {0, 1, 0, 2}},
            {{1, 2, 0, 0}, {1, 2, 1, 1}, {1, 2, 2, 2}, {1, 2, 0, 1}, {1, 2, 1, 2}, {1, 2, 0, 2}},
            {{0, 2, 0, 0}, {0, 2, 1, 1}, {0, 2, 2, 2}, {0, 2, 0, 1}, {0, 2, 1, 2}, {0, 2, 0, 2}}};
        return FOURTH[voigt_row][voigt_col][target_index];
      };

      /// instancing of this class is forbidden
      IndexMappings() = delete;

     private:
      static inline void assertRangeMatrixIndex(unsigned int row, unsigned int col)
      {
        dsassert(row < 3, "given row index out of range [0,2]");
        dsassert(col < 3, "given col index out of range [0,2]");
      }

      static inline void assertRangeVoigtIndex(unsigned int index)
      {
        dsassert(index < 6, "given index out of range [0,5]");
      }
    };
  }  // namespace VOIGT
}  // namespace UTILS

using VStrainUtils = UTILS::VOIGT::VoigtUtils<UTILS::VOIGT::NotationType::strain>;
using VStressUtils = UTILS::VOIGT::VoigtUtils<UTILS::VOIGT::NotationType::stress>;

#endif  // VOIGT_NOTATION_H