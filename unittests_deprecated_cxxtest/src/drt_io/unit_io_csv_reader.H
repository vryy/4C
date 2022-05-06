/*----------------------------------------------------------------------*/
/*! \file

\brief Unittests for the csv reader

\level 1

*-----------------------------------------------------------------------*/

#ifndef UNIT_IO_CSVREADER_H_
#define UNIT_IO_CSVREADER_H_

#include "src/common/unit_cxx_test_wrapper.H"
#include "src/common/special_assertions.H"
#include "src/drt_io/csv_reader.H"

namespace IO
{
  class CsvReader_TestSuite;
}

class IO::CsvReader_TestSuite : public BACICxxTestWrapper
{
 public:
  void TestDataProcessingCsvStream()
  {
    const std::vector<double> x = {0.3, 0.4, 0.45};
    const std::vector<double> y = {4.3, 4.1, 4.15};
    const std::vector<double> z = {-1.0, 0.1, 1.3};

    std::stringstream test_csv_file_stream;
    test_csv_file_stream << "#x,y,z" << std::endl;
    for (std::size_t i = 0; i < x.size(); ++i)
      test_csv_file_stream << std::to_string(x[i]) << "," << std::to_string(y[i]) << ","
                           << std::to_string(z[i]) << std::endl;

    auto csv_values = IO::ReadCsv(3, test_csv_file_stream);

    TESTING::AssertDelta(csv_values[0], x, 1.0e-16);
    TESTING::AssertDelta(csv_values[1], y, 1.0e-16);
    TESTING::AssertDelta(csv_values[2], z, 1.0e-16);
  }

  void TestDataProcessingCsvFile()
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

    TESTING::AssertDelta(csv_values[0], x, 1.0e-16);
    TESTING::AssertDelta(csv_values[1], y, 1.0e-16);
    TESTING::AssertDelta(csv_values[2], z, 1.0e-16);
  }

  void TestDifferentColumnLength()
  {
    std::stringstream test_csv_file;
    test_csv_file << "#x,y" << std::endl;
    test_csv_file << "0.30,4.40" << std::endl;
    test_csv_file << "0.30," << std::endl;

    TS_ASSERT_THROWS(IO::ReadCsv(3, test_csv_file), const std::runtime_error &)
  }

  void TestTrailingComma()
  {
    std::stringstream test_csv_file;
    test_csv_file << "#x,y" << std::endl;
    test_csv_file << "0.30,4.40," << std::endl;

    TS_ASSERT_THROWS(IO::ReadCsv(2, test_csv_file), const std::runtime_error &)
  }

  void TestWrongColumnNumber()
  {
    std::stringstream test_csv_file;
    test_csv_file << "#x,y" << std::endl;
    test_csv_file << "0.30,4.40" << std::endl;

    TS_ASSERT_THROWS(IO::ReadCsv(3, test_csv_file), const std::runtime_error &)
  }

  void TestWrongHeaderStyle()
  {
    std::stringstream test_csv_file;
    test_csv_file << "x,y" << std::endl;
    test_csv_file << "0.30,4.40" << std::endl;

    TS_ASSERT_THROWS(IO::ReadCsv(2, test_csv_file), const std::runtime_error &)
  }

  void TestWrongInputDataType()
  {
    std::stringstream test_csv_file;
    test_csv_file << "x,y" << std::endl;
    test_csv_file << "0.30,a" << std::endl;

    TS_ASSERT_THROWS(IO::ReadCsv(2, test_csv_file), const std::runtime_error &)
  }

  void TestWrongSeparator()
  {
    std::stringstream test_csv_file;
    test_csv_file << "x;y" << std::endl;
    test_csv_file << "0.30;4.40" << std::endl;

    TS_ASSERT_THROWS(IO::ReadCsv(2, test_csv_file), const std::runtime_error &)
  }
};

#endif
