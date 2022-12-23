/*---------------------------------------------------------------------------*/
/*! \file

\brief Unittests for tensor transforations

\level 3
*/
/*----------------------------------------------------------------------*/

#include <gtest/gtest.h>
#include "linalg_fixedsizematrix.H"
#include "tensor_transformation.H"
#include "unittests_assertions.h"

namespace
{
  class TensorTransformationTest : public testing::Test
  {
   public:
    LINALG::Matrix<3, 3> tens1_;
    LINALG::Matrix<3, 3> tens2_;
    LINALG::Matrix<3, 3> rotatedTens1_;
    LINALG::Matrix<3, 3> rotatedTens2_;
    LINALG::Matrix<3, 3> rotationMatrix1_;
    LINALG::Matrix<3, 3> rotationMatrix2_;

    TensorTransformationTest()
    {
      tens1_(0, 0) = 0.771320643266746;
      tens1_(0, 1) = 0.0207519493594015;
      tens1_(0, 2) = 0.6336482349262754;
      tens1_(1, 0) = 0.7488038825386119;
      tens1_(1, 1) = 0.4985070123025904;
      tens1_(1, 2) = 0.22479664553084766;
      tens1_(2, 0) = 0.19806286475962398;
      tens1_(2, 1) = 0.7605307121989587;
      tens1_(2, 2) = 0.16911083656253545;

      tens2_(0, 0) = 0.08833981417401027;
      tens2_(0, 1) = 0.6853598183677972;
      tens2_(0, 2) = 0.9533933461949365;
      tens2_(1, 0) = 0.003948266327914451;
      tens2_(1, 1) = 0.5121922633857766;
      tens2_(1, 2) = 0.8126209616521135;
      tens2_(2, 0) = 0.6125260668293881;
      tens2_(2, 1) = 0.7217553174317995;
      tens2_(2, 2) = 0.29187606817063316;

      rotationMatrix1_(0, 0) = 0.6498181078537578;
      rotationMatrix1_(0, 1) = -0.12009199850537022;
      rotationMatrix1_(0, 2) = 0.7505426960542272;
      rotationMatrix1_(1, 0) = 0.6906850791495796;
      rotationMatrix1_(1, 1) = 0.5055249243539574;
      rotationMatrix1_(1, 2) = -0.517106055173467;
      rotationMatrix1_(2, 0) = -0.31731774004221847;
      rotationMatrix1_(2, 1) = 0.8544135197618961;
      rotationMatrix1_(2, 2) = 0.41144500130951545;

      rotationMatrix2_(0, 0) = 0.7179663436165873;
      rotationMatrix2_(0, 1) = -0.5820205817986273;
      rotationMatrix2_(0, 2) = 0.3818067204707322;
      rotationMatrix2_(1, 0) = 0.6324373962072567;
      rotationMatrix2_(1, 1) = 0.7745645248692012;
      rotationMatrix2_(1, 2) = -0.008528580933186578;
      rotationMatrix2_(2, 0) = -0.2907701313966307;
      rotationMatrix2_(2, 1) = 0.2475920822177799;
      rotationMatrix2_(2, 2) = 0.9242028411072161;

      rotatedTens1_(0, 0) = 0.6849254584795603;
      rotatedTens1_(0, 1) = 0.3872640091723427;
      rotatedTens1_(0, 2) = 0.4809281559181807;
      rotatedTens1_(1, 0) = 0.8889894437547831;
      rotatedTens1_(1, 1) = 0.25464227926695077;
      rotatedTens1_(1, 2) = -0.17427364095084866;
      rotatedTens1_(2, 0) = 0.26718543150292134;
      rotatedTens1_(2, 1) = 0.6679797997370243;
      rotatedTens1_(2, 2) = 0.4993707543853608;

      rotatedTens2_(0, 0) = 0.061836304421077744;
      rotatedTens2_(0, 1) = 0.5474915926576807;
      rotatedTens2_(0, 2) = 0.32897785044946315;
      rotatedTens2_(1, 0) = 0.028546155066549384;
      rotatedTens2_(1, 1) = 0.6617296640042815;
      rotatedTens2_(1, 2) = 1.3250816081673256;
      rotatedTens2_(2, 0) = 0.11660950602923983;
      rotatedTens2_(2, 1) = 0.80128385159843;
      rotatedTens2_(2, 2) = 0.16884217730506076;
    };
  };

  TEST_F(TensorTransformationTest, TensorRotation)
  {
    LINALG::Matrix<3, 3> tmp;
    UTILS::TENSOR::TensorRotation(rotationMatrix1_, tens1_, tmp);
    BACI_EXPECT_NEAR(tmp, rotatedTens1_, 1e-9);

    UTILS::TENSOR::TensorRotation(rotationMatrix2_, tens2_, tmp);
    BACI_EXPECT_NEAR(tmp, rotatedTens2_, 1e-9);
  }

  TEST_F(TensorTransformationTest, InverseTensorRotation)
  {
    LINALG::Matrix<3, 3> tmp;
    UTILS::TENSOR::InverseTensorRotation(rotationMatrix1_, rotatedTens1_, tmp);
    BACI_EXPECT_NEAR(tmp, tens1_, 1e-9);

    UTILS::TENSOR::InverseTensorRotation(rotationMatrix2_, rotatedTens2_, tmp);
    BACI_EXPECT_NEAR(tmp, tens2_, 1e-9);
  }
}  // namespace