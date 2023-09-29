/*----------------------------------------------------------------------*/
/*! \file

\brief Unittests for functions in IO namespace

\level 1

*-----------------------------------------------------------------------*/
#include <gtest/gtest.h>

#include "baci_io_string_converter.H"

#include <map>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace
{
  using namespace BACI;

  TEST(ParseTests, ParseInt) { EXPECT_EQ(IO::StringConverter<int>::Parse("5"), 5); }

  TEST(ParseTests, ParseDouble) { EXPECT_EQ(IO::StringConverter<double>::Parse("5.5"), 5.5); }

  TEST(ParseTests, ParseChar) { EXPECT_EQ(IO::StringConverter<char>::Parse("a"), 'a'); }

  TEST(ParseTests, ParseBool)
  {
    EXPECT_EQ(IO::StringConverter<bool>::Parse("true"), true);
    EXPECT_EQ(IO::StringConverter<bool>::Parse("false"), false);
    EXPECT_THROW(IO::StringConverter<bool>::Parse("abc"), CORE::Exception);
  }

  TEST(ParseTests, ParseArray)
  {
    const std::string str = "1.1,2.2,3.3, 4.4,5.5";
    std::array<double, 5> arr = IO::StringConverter<std::array<double, 5>>::Parse(str);
    std::array<double, 5> expected_arr = {1.1, 2.2, 3.3, 4.4, 5.5};
    EXPECT_EQ(arr, expected_arr);
  }

  TEST(ParseTests, ParseArrayWrongArrayDimension)
  {
    const std::string str = "1.1,2.2,  3.3,4.4,5.5,6.6";
    EXPECT_THROW((IO::StringConverter<std::array<double, 5>>::Parse(str)), CORE::Exception);
  }

  TEST(ParseTests, ParseArrayOfArrayOfArrayOfArrays)
  {
    const std::string str = "1.1,2.2;1.1,2.2 |1.1,2.2;1.1,2.2% 1.1,2.2; 1.1,2.2| 1.1,2.2;1.1,2.2";
    auto arr = IO::StringConverter<
        std::array<std::array<std::array<std::array<double, 2>, 2>, 2>, 2>>::Parse(str);
    std::array<std::array<std::array<std::array<double, 2>, 2>, 2>, 2> expected_arr;
    for (auto& arr1 : expected_arr)
    {
      for (auto& arr2 : arr1)
      {
        for (auto& arr3 : arr2)
        {
          arr3 = {1.1, 2.2};
        }
      }
    }
    EXPECT_EQ(arr, expected_arr);
  }

  TEST(ParseTests, ParseVector)
  {
    const std::string str = "1,2,3";
    auto vec = IO::StringConverter<std::vector<double>>::Parse(str);
    std::vector<double> expected_vec = {1, 2, 3};
    EXPECT_EQ(vec, expected_vec);
  }

  TEST(ParseTests, ParseVectorOfVectorsOfVectors)
  {
    const std::string str = "1,2,3;4,5,6;7,8,9|1,2,3;4,5 ,6;7,8,9|1,2,3;4,5,6;7,8,9";
    auto vec = IO::StringConverter<std::vector<std::vector<std::vector<double>>>>::Parse(str);
    std::vector<std::vector<std::vector<double>>> expected_vec = {{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}},
        {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}, {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}};
    EXPECT_EQ(vec, expected_vec);
  }

  TEST(ParseTests, ParseMap)
  {
    const std::string str = "100 :2.2,200 :4.4";
    std::map<int, double> map = IO::StringConverter<std::map<int, double>>::Parse(str);
    std::map<int, double> expected_map;
    expected_map[100] = 2.2;
    expected_map[200] = 4.4;

    EXPECT_EQ(map, expected_map);
  }

  TEST(ParseTests, ParseMapArrayValues)
  {
    const std::string str = "1.1:10,20,30; 2.2:1,2,3";
    std::map<double, std::array<int, 3>> map =
        IO::StringConverter<std::map<double, std::array<int, 3>>>::Parse(str);
    std::map<double, std::array<int, 3>> expected_map;
    std::array<int, 3> array_val_1 = {10, 20, 30};
    std::array<int, 3> array_val_2 = {1, 2, 3};
    expected_map[1.1] = array_val_1;
    expected_map[2.2] = array_val_2;
    EXPECT_EQ(map, expected_map);
  }

  TEST(ParseTests, ParseMapArrayKeys)
  {
    const std::string str = "1,2,3: 2.2;4,5,6:7.7";
    std::map<std::array<int, 3>, double> map =
        IO::StringConverter<std::map<std::array<int, 3>, double>>::Parse(str);
    std::map<std::array<int, 3>, double> expected_map;
    std::array<int, 3> array_key = {1, 2, 3};
    std::array<int, 3> array_key_2 = {4, 5, 6};

    expected_map[array_key] = 2.2;
    expected_map[array_key_2] = 7.7;

    EXPECT_EQ(map, expected_map);
  }

  TEST(ParseTests, ParseMapPairValues)
  {
    const std::string str = "100:2.2,200:4.4";
    std::map<int, double> map = IO::StringConverter<std::map<int, double>>::Parse(str);
    std::map<int, double> expected_map;
    expected_map[100] = 2.2;
    expected_map[200] = 4.4;

    EXPECT_EQ(map, expected_map);
  }

  TEST(ParseTests, ParseTuple)
  {
    const std::string str = "1,2,3.3";
    std::tuple<int, int, double> tuple =
        IO::StringConverter<std::tuple<int, int, double>>::Parse(str);
    std::tuple<int, int, double> expected_tuple = std::make_tuple(1, 2, 3.3);
    EXPECT_EQ(tuple, expected_tuple);
  }

  TEST(ParseTests, ParseTupleWrongDimension)
  {
    EXPECT_THROW(
        (IO::StringConverter<std::tuple<int, int, double>>::Parse("1,2,3.3,9.0")), CORE::Exception);
    EXPECT_THROW(
        (IO::StringConverter<std::tuple<int, int, double>>::Parse("1,2")), CORE::Exception);
  }

  TEST(ParseTests, ParsePair)
  {
    const std::string str = "1,2.2";
    std::pair<int, double> pair = IO::StringConverter<std::pair<int, double>>::Parse(str);
    std::pair<int, double> expected_pair = std::make_pair(1, 2.2);
    EXPECT_EQ(pair, expected_pair);
  }

  TEST(ParseTests, ParsePairWrongDimension)
  {
    EXPECT_THROW(
        (IO::StringConverter<std::pair<int, double>>::Parse("1,2.2,3.3")), CORE::Exception);
    EXPECT_THROW((IO::StringConverter<std::pair<int, double>>::Parse("1;7")), CORE::Exception);
  }

  TEST(ParseTests, ComplexExample)
  {
    const std::string str = "0:0.0,0.1,0.2;1:1.0, 1.1,1.2 % 0.8;a,b,c%1,2;3,4|5,6;7,8";

    // 0:0.0,0.1,0.2;1:1.0,1.1,1.2
    std::vector<double> vec1 = {0.0, 0.1, 0.2};
    std::vector<double> vec2 = {1.0, 1.1, 1.2};
    std::map<int, std::vector<double>> map;
    map[0] = vec1;
    map[1] = vec2;

    // 0.8;a,b,c
    std::list<char> list = {'a', 'b', 'c'};
    std::pair<double, std::list<char>> pair = std::make_pair(0.8, list);

    // 1,2;3,4|5,6;7,8
    std::vector<std::vector<std::vector<int>>> vec = {{{1, 2}, {3, 4}}, {{5, 6}, {7, 8}}};

    std::tuple<std::map<int, std::vector<double>>, std::pair<std::string, std::list<char>>,
        std::vector<std::vector<std::vector<int>>>>
        type;
    auto tuple = IO::StringConverter<
        std::tuple<std::map<int, std::vector<double>>, std::pair<double, std::list<char>>,
            std::vector<std::vector<std::vector<int>>>>>::Parse(str);
    std::tuple<std::map<int, std::vector<double>>, std::pair<double, std::list<char>>,
        std::vector<std::vector<std::vector<int>>>>
        expected_tuple = std::make_tuple(map, pair, vec);

    EXPECT_EQ(std::get<0>(tuple), std::get<0>(expected_tuple));
    EXPECT_EQ(std::get<1>(tuple), std::get<1>(expected_tuple));
    EXPECT_EQ(std::get<2>(tuple), std::get<2>(expected_tuple));
  }
}  // namespace