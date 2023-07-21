/*----------------------------------------------------------------------*/
/*! \file

\brief Unittests for the csv reader

\level 1

*-----------------------------------------------------------------------*/
#include <gtest/gtest.h>

#include <fstream>

#include "baci_io_csv_reader.H"

#include "baci_unittest_utils_assertions.h"

namespace
{
  TEST(CsvReaderTest, DataProcessingCsvStream)
  {
    const std::vector<double> x = {0.2, 0.4, 0.45};
    const std::vector<double> y = {4.3, 4.1, 4.15};
    const std::vector<double> z = {-1.0, 0.1, 1.3};

    std::stringstream test_csv_file_stream;
    test_csv_file_stream << "#x,y,z" << std::endl;
    for (std::size_t i = 0; i < x.size(); ++i)
      test_csv_file_stream << std::to_string(x[i]) << "," << std::to_string(y[i]) << ","
                           << std::to_string(z[i]) << std::endl;

    auto csv_values = IO::ReadCsv(3, test_csv_file_stream);

    EXPECT_EQ(csv_values[0], x);
    EXPECT_EQ(csv_values[1], y);
    EXPECT_EQ(csv_values[2], z);
  }

  TEST(CsvReaderTest, DataProcessingCsvFile)
  {
    const std::vector<double> x = {0.3, 0.4, 0.45};
    const std::vector<double> y = {4.3, 4.1, 4.15};
    const std::vector<double> z = {-1.0, 0.1, 1.3};

    const std::string csv_template_file_name = "test.csv";
    std::ofstream test_csv_file(csv_template_file_name);

    // include header line
    test_csv_file << "#x,y,z" << std::endl;
    for (std::size_t i = 0; i < x.size(); ++i)
      test_csv_file << std::to_string(x[i]) << "," << std::to_string(y[i]) << ","
                    << std::to_string(z[i]) << std::endl;
    // close template file
    test_csv_file.close();

    auto csv_values = IO::ReadCsv(3, csv_template_file_name);

    EXPECT_EQ(csv_values[0], x);
    EXPECT_EQ(csv_values[1], y);
    EXPECT_EQ(csv_values[2], z);
  }


  TEST(CsvReaderTest, DifferentColumnLengthThrows)
  {
    std::stringstream test_csv_file;
    test_csv_file << "#x,y" << std::endl;
    test_csv_file << "0.30,4.40" << std::endl;
    test_csv_file << "0.30," << std::endl;

    BACI_EXPECT_THROW_WITH_MESSAGE(
        IO::ReadCsv(3, test_csv_file), std::runtime_error, "same length");
  }

  TEST(CsvReaderTest, TrailingCommaThrows)
  {
    std::stringstream test_csv_file;
    test_csv_file << "#x,y" << std::endl;
    test_csv_file << "0.30,4.40," << std::endl;

    BACI_EXPECT_THROW_WITH_MESSAGE(
        IO::ReadCsv(2, test_csv_file), std::runtime_error, "trailing comma");
  }

  TEST(CsvReaderTest, WrongColumnNumberThrows)
  {
    std::stringstream test_csv_file;
    test_csv_file << "#x,y" << std::endl;
    test_csv_file << "0.30,4.40" << std::endl;

    BACI_EXPECT_THROW_WITH_MESSAGE(IO::ReadCsv(3, test_csv_file), std::runtime_error, "");
  }

  TEST(CsvReaderTest, WrongHeaderStyleThrows)
  {
    std::stringstream test_csv_file;
    test_csv_file << "x,y" << std::endl;
    test_csv_file << "0.30,4.40" << std::endl;

    BACI_EXPECT_THROW_WITH_MESSAGE(IO::ReadCsv(2, test_csv_file), std::runtime_error, "header");
  }

  TEST(CsvReaderTest, WrongInputDataTypeThrows)
  {
    std::stringstream test_csv_file;
    test_csv_file << "x,y" << std::endl;
    test_csv_file << "0.30,a" << std::endl;

    BACI_EXPECT_THROW_WITH_MESSAGE(IO::ReadCsv(2, test_csv_file), std::runtime_error, "numbers");
  }

  TEST(CsvReaderTest, WrongSeparatorThrows)
  {
    std::stringstream test_csv_file;
    test_csv_file << "x;y" << std::endl;
    test_csv_file << "0.30;4.40" << std::endl;

    BACI_EXPECT_THROW_WITH_MESSAGE(
        IO::ReadCsv(2, test_csv_file), std::runtime_error, "separated by commas");
  }

}  // namespace