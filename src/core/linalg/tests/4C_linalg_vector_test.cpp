/*----------------------------------------------------------------------*/
/*! \file

\brief Unit tests for the Vector wrapper

\level 0
*/


#include <gtest/gtest.h>

#include "4C_linalg_vector.hpp"


// Epetra related headers
#include <Epetra_Map.h>
#include <Epetra_MpiComm.h>
#include <Teuchos_RCPDecl.hpp>

FOUR_C_NAMESPACE_OPEN

namespace
{
  class VectorTest : public testing::Test
  {
   public:
    Teuchos::RCP<Epetra_Comm> comm_;
    Teuchos::RCP<Epetra_Map> map_;
    int NumGlobalElements = 10;

   protected:
    VectorTest()
    {
      // set up communicator
      comm_ = Teuchos::make_rcp<Epetra_MpiComm>(MPI_COMM_WORLD);

      // set up a map
      map_ = Teuchos::make_rcp<Epetra_Map>(NumGlobalElements, 0, *comm_);
    }
  };

  TEST_F(VectorTest, ConstructorsAndNorms)
  {
    // create an epetra vector
    Epetra_Vector my_epetra_vector = Epetra_Vector(*map_, true);

    // try to copy zero vector into wrapper
    Core::LinAlg::Vector<double> epetra_based_test_vector =
        Core::LinAlg::Vector<double>(my_epetra_vector);

    // create vector
    Core::LinAlg::Vector<double> test_vector = Core::LinAlg::Vector<double>(*map_, true);

    // initialize with wrong value
    double norm_of_test_vector = 1;

    test_vector.Print(std::cout);
    test_vector.Norm2(&norm_of_test_vector);
    // test norm2 and success of both vectors
    ASSERT_FLOAT_EQ(0.0, norm_of_test_vector);

    // reset value
    norm_of_test_vector = 1;

    // check result of Norm2
    epetra_based_test_vector.Norm2(&norm_of_test_vector);
    ASSERT_FLOAT_EQ(0.0, norm_of_test_vector);

    // test element access function for proc 0
    if (comm_->MyPID() == 0) test_vector[1] = 1;

    // check result of Norm1
    test_vector.Norm1(&norm_of_test_vector);
    ASSERT_FLOAT_EQ(1.0, norm_of_test_vector);

    test_vector[1] = 100.0;

    // check result of NormInf
    test_vector.NormInf(&norm_of_test_vector);
    ASSERT_FLOAT_EQ(100.0, norm_of_test_vector);
  }

  TEST_F(VectorTest, DeepCopying)
  {
    Core::LinAlg::Vector<double> a(*map_, true);
    a.PutScalar(1.0);

    Core::LinAlg::Vector<double> b(*map_, true);
    // copy assign
    b = a;
    b.PutScalar(2.0);
    double norm_a = 0.0;
    double norm_b = 0.0;
    a.Norm2(&norm_a);
    b.Norm2(&norm_b);

    EXPECT_FLOAT_EQ(norm_a, 1.0 * std::sqrt(NumGlobalElements));
    EXPECT_FLOAT_EQ(norm_b, 2.0 * std::sqrt(NumGlobalElements));

    // copy constructor
    Core::LinAlg::Vector<double> c(a);
    c.PutScalar(3.0);
    double norm_c = 0.0;
    c.Norm2(&norm_c);
    EXPECT_FLOAT_EQ(norm_c, 3.0 * std::sqrt(NumGlobalElements));
  }

  TEST_F(VectorTest, PutScalar)
  {
    // initialize with false value
    double norm_of_test_vector = 0.0;

    // copy zero vector into new interface
    Core::LinAlg::Vector<double> test_vector = Core::LinAlg::Vector<double>(*map_, true);

    test_vector.PutScalar(2.0);

    // check result
    test_vector.Norm2(&norm_of_test_vector);
    ASSERT_FLOAT_EQ(NumGlobalElements * 2.0 * 2.0, norm_of_test_vector * norm_of_test_vector);
  }

  TEST_F(VectorTest, Update)
  {
    Core::LinAlg::Vector<double> a = Core::LinAlg::Vector<double>(*map_, true);
    a.PutScalar(1.0);

    Core::LinAlg::Vector<double> b = Core::LinAlg::Vector<double>(*map_, true);
    b.PutScalar(1.0);

    // update the vector
    b.Update(2.0, a, 3.0);

    // initialize with false value
    double b_norm = 0.0;

    // check norm of vector
    b.Norm2(&b_norm);
    ASSERT_FLOAT_EQ(NumGlobalElements * (2.0 + 3.0) * (2.0 + 3.0), b_norm * b_norm);

    Core::LinAlg::Vector<double> c = Core::LinAlg::Vector<double>(*map_, true);
    c.Update(1, a, -1, b, 0);

    // initialize with false value
    double c_norm = 0.0;

    // check norm of vector
    c.Norm1(&c_norm);
    ASSERT_FLOAT_EQ(4 * NumGlobalElements, c_norm);
  }


  TEST_F(VectorTest, View)
  {
    Epetra_Vector a(*map_, true);
    a.PutScalar(1.0);
    // Scope in which a is modified by the view
    {
      Core::LinAlg::VectorView a_view(a);

      double norm = 0.0;
      ((Core::LinAlg::Vector<double>&)a_view).Norm2(&norm);
      EXPECT_EQ(norm, std::sqrt(NumGlobalElements));

      ((Core::LinAlg::Vector<double>&)a_view).PutScalar(2.0);
    }
    const Epetra_Vector& a_const = a;
    Core::LinAlg::VectorView a_view_const(a_const);
    // Change must be reflected in a
    double norm = 0.0;
    static_cast<const Core::LinAlg::Vector<double>&>(a_view_const).Norm2(&norm);
    EXPECT_EQ(norm, 2.0 * std::sqrt(NumGlobalElements));
  }

  std::vector<double> means_multi_vector(const Core::LinAlg::MultiVector<double>& mv)
  {
    std::vector<double> means(mv.NumVectors());
    mv.MeanValue(means.data());
    return means;
  }


  TEST_F(VectorTest, MultiVectorImplicitConversionView)
  {
    Core::LinAlg::Vector<double> a(*map_, true);
    a.PutScalar(1.0);

    // This views the data that is in a. It does not copy the data.
    // This results in the same behavior as inheritance would give.
    EXPECT_EQ(means_multi_vector(a)[0], 1.0);

    // This copies the data.
    Core::LinAlg::MultiVector<double> mv = a;
    a.PutScalar(2.0);

    // mv should still be 1.0 because we only modified a.
    EXPECT_EQ(means_multi_vector(mv)[0], 1.0);
  }

  TEST_F(VectorTest, MultiVectorImplicitConversionCopy)
  {
    auto a = Teuchos::make_rcp<Core::LinAlg::Vector<double>>(*map_, true);
    a->PutScalar(1.0);

    // This copies the data.
    Core::LinAlg::MultiVector<double> mv = *a;
    a->PutScalar(2.0);
    // Explicitly deallocate a to make sure that mv is not a view.
    a = Teuchos::null;

    // mv should still be 1.0 because we only modified a.
    EXPECT_EQ(means_multi_vector(mv)[0], 1.0);
  }

  TEST_F(VectorTest, MultiVectorImplicitConversionRef)
  {
    Core::LinAlg::Vector<double> a(*map_, true);
    a.PutScalar(1.0);

    Core::LinAlg::MultiVector<double>& mv = a;
    mv.PutScalar(2.0);
    EXPECT_EQ(means_multi_vector(a)[0], 2.0);

    // Reassigning to a must keep mv valid: move assign
    a = Core::LinAlg::Vector<double>(*map_, true);
    EXPECT_EQ(means_multi_vector(mv)[0], 0.0);
    a.PutScalar(3.0);
    EXPECT_EQ(means_multi_vector(mv)[0], 3.0);

    // Reassigning to a must keep mv valid: copy assign
    Core::LinAlg::Vector<double> b(*map_, true);
    a = b;
    EXPECT_EQ(means_multi_vector(mv)[0], 0.0);
    a.PutScalar(4.0);
    EXPECT_EQ(means_multi_vector(mv)[0], 4.0);
  }

  TEST_F(VectorTest, AssignToRef)
  {
    Core::LinAlg::Vector<double> a(*map_, true);
    a.PutScalar(1.0);
    EXPECT_EQ(means_multi_vector(a)[0], 1.0);
    Core::LinAlg::MultiVector<double>& mv = a;
    // Actually assign an MV to a via the ref. Note that this would throw in Trilinos if not using a
    // single column.
    mv = Core::LinAlg::MultiVector<double>(*map_, 1, true);
    EXPECT_EQ(means_multi_vector(mv)[0], 0.0);
  }

  TEST_F(VectorTest, VectorFromMultiVector)
  {
    Core::LinAlg::MultiVector<double> mv(*map_, 3, true);
    mv.PutScalar(1.0);

    const int index = 1;

    Core::LinAlg::Vector<double>& a = mv(index);
    EXPECT_EQ(means_multi_vector(a)[0], 1.0);

    a.PutScalar(2.0);

    // Check that the change is reflected in the MultiVector
    EXPECT_EQ(means_multi_vector(mv), (std::vector{1., 2., 1.}));

    // Another MultiVector conversion
    Core::LinAlg::MultiVector<double>& mv2 = a;
    mv2.PutScalar(3.0);
    EXPECT_EQ(means_multi_vector(mv), (std::vector{1., 3., 1.}));

    // Combine with taking a view
    {
      const auto put_scalar = [](Core::LinAlg::MultiVector<double>& v, double s)
      { v.PutScalar(s); };
      Core::LinAlg::VectorView view_mv2((Epetra_MultiVector&)mv2);
      put_scalar(view_mv2, 4.0);
    }
    EXPECT_EQ(means_multi_vector(mv), (std::vector{1., 4., 1.}));
  }

}  // namespace

FOUR_C_NAMESPACE_CLOSE
