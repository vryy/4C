/*----------------------------------------------------------------------*/
/*! \file

\brief Unittests for the file reader

\level 1

*-----------------------------------------------------------------------*/
#include <gtest/gtest.h>

#include "baci_io_file_reader.hpp"
#include "baci_unittest_utils_assertions_test.hpp"

#include <fstream>
#include <vector>

namespace
{

  using namespace FourC;

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

    auto csv_values = IO::ReadCsvAsColumns(3, test_csv_file_stream);

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

    auto csv_values = IO::ReadCsvAsColumns(3, csv_template_file_name);

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

    FOUR_C_EXPECT_THROW_WITH_MESSAGE(
        IO::ReadCsvAsColumns(3, test_csv_file), CORE::Exception, "same length");
  }

  TEST(CsvReaderTest, TrailingCommaThrows)
  {
    std::stringstream test_csv_file;
    test_csv_file << "#x,y" << std::endl;
    test_csv_file << "0.30,4.40," << std::endl;

    FOUR_C_EXPECT_THROW_WITH_MESSAGE(
        IO::ReadCsvAsColumns(2, test_csv_file), CORE::Exception, "trailing comma");
  }

  TEST(CsvReaderTest, WrongColumnNumberThrows)
  {
    std::stringstream test_csv_file;
    test_csv_file << "#x,y" << std::endl;
    test_csv_file << "0.30,4.40" << std::endl;

    FOUR_C_EXPECT_THROW_WITH_MESSAGE(IO::ReadCsvAsColumns(3, test_csv_file), CORE::Exception, "");
  }

  TEST(CsvReaderTest, WrongHeaderStyleThrows)
  {
    std::stringstream test_csv_file;
    test_csv_file << "x,y" << std::endl;
    test_csv_file << "0.30,4.40" << std::endl;

    FOUR_C_EXPECT_THROW_WITH_MESSAGE(
        IO::ReadCsvAsColumns(2, test_csv_file), CORE::Exception, "header");
  }

  TEST(CsvReaderTest, WrongInputDataTypeThrows)
  {
    std::stringstream test_csv_file;
    test_csv_file << "x,y" << std::endl;
    test_csv_file << "0.30,a" << std::endl;

    FOUR_C_EXPECT_THROW_WITH_MESSAGE(
        IO::ReadCsvAsColumns(2, test_csv_file), CORE::Exception, "numbers");
  }

  TEST(CsvReaderTest, WrongSeparatorThrows)
  {
    std::stringstream test_csv_file;
    test_csv_file << "x;y" << std::endl;
    test_csv_file << "0.30;4.40" << std::endl;

    FOUR_C_EXPECT_THROW_WITH_MESSAGE(
        IO::ReadCsvAsColumns(2, test_csv_file), CORE::Exception, "separated by commas");
  }

  TEST(ReadFileAsLines, ValidFile)
  {
    using T = std::map<int, std::array<int, 3>>;
    std::stringstream test_file;
    test_file << "1:1,2,3" << std::endl;
    test_file << "2:4,5,6" << std::endl;
    test_file << "3:7,8,9" << std::endl;

    std::vector<T> expected_data = {{{1, {1, 2, 3}}}, {{2, {4, 5, 6}}}, {{3, {7, 8, 9}}}};

    std::vector<T> read_data = IO::ReadFileAsLines<T>(test_file);

    EXPECT_EQ(read_data, expected_data);
  }

  TEST(FoldLines, ReduceMapToVector)
  {
    using T = std::map<int, std::array<int, 3>>;
    std::stringstream test_file;
    test_file << "1:1,2,3" << std::endl;
    test_file << "2:4,5,6" << std::endl;
    test_file << "3:7,8,9" << std::endl;

    std::vector<T> expected_data = {{{1, {1, 2, 3}}}, {{3, {7, 8, 9}}}};

    T read_data = IO::ReadFileAsLines<T, T>(test_file,
        [](T&& acc, T&& next)
        {
          // add only maps that have a key != 2
          if (acc.size() == 1 && acc.count(2) == 0) acc.merge(std::move(next));
          return acc;
        });
  }

  TEST(FoldLines, ReduceMapSumValues)
  {
    using T = std::map<int, std::array<int, 3>>;
    std::stringstream test_file;
    test_file << "1:1,2,3" << std::endl;
    test_file << "2:4,5,6" << std::endl;
    test_file << "3:7,8,9" << std::endl;

    std::map<int, int> expected_data{{1, 6}, {2, 15}, {3, 24}};

    using ReducedType = std::map<int, int>;
    ReducedType read_data = IO::ReadFileAsLines<T, ReducedType>(test_file,
        [](ReducedType&& acc, T&& next)
        {
          for (const auto& [key, value] : next)
          {
            acc[key] = value[0] + value[1] + value[2];
          }
          return acc;
        });

    EXPECT_EQ(read_data, expected_data);
  }

  TEST(FoldLines, ReduceMapToMap)
  {
    using T = std::map<int, std::vector<std::pair<double, double>>>;
    std::stringstream test_file;
    test_file << "1:0.0,0.0;0.1,0.1;0.2,0.2" << std::endl;
    test_file << "2:0.0,0.0;0.1,0.2;0.2,0.4" << std::endl;

    T expected_data = {
        {1, {{0.0, 0.0}, {0.1, 0.1}, {0.2, 0.2}}}, {2, {{0.0, 0.0}, {0.1, 0.2}, {0.2, 0.4}}}};

    T read_data = IO::ReadFileAsLines<T, T>(test_file,
        [](T&& acc, T&& next)
        {
          for (const auto& [key, value] : next)
          {
            acc[key] = value;
          }
          return acc;
        });

    EXPECT_EQ(read_data, expected_data);
  }
}  // namespace
