/*----------------------------------------------------------------------*/
/*! \file

\brief Unittests for functions in IO namespace

\level 1

*-----------------------------------------------------------------------*/
#include <gtest/gtest.h>

#include "4C_io_string_converter.hpp"

#include <map>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace
{
  using namespace FourC;

  TEST(ParseTests, ParseInt) { EXPECT_EQ(Core::IO::StringConverter<int>::parse("5"), 5); }

  TEST(ParseTests, ParseDouble) { EXPECT_EQ(Core::IO::StringConverter<double>::parse("5.5"), 5.5); }

  TEST(ParseTests, ParseChar) { EXPECT_EQ(Core::IO::StringConverter<char>::parse("a"), 'a'); }

  TEST(ParseTests, ParseBool)
  {
    EXPECT_EQ(Core::IO::StringConverter<bool>::parse("true"), true);
    EXPECT_EQ(Core::IO::StringConverter<bool>::parse("false"), false);
    EXPECT_THROW(Core::IO::StringConverter<bool>::parse("abc"), Core::Exception);
  }

  TEST(ParseTests, ParseArray)
  {
    const std::string str = "1.1,2.2,3.3, 4.4,5.5";
    std::array<double, 5> arr = Core::IO::StringConverter<std::array<double, 5>>::parse(str);
    std::array<double, 5> expected_arr = {1.1, 2.2, 3.3, 4.4, 5.5};
    EXPECT_EQ(arr, expected_arr);
  }

  TEST(ParseTests, ParseArrayWrongArrayDimension)
  {
    const std::string str = "1.1,2.2,  3.3,4.4,5.5,6.6";
    EXPECT_THROW((Core::IO::StringConverter<std::array<double, 5>>::parse(str)), Core::Exception);
  }

  TEST(ParseTests, ParseArrayOfArrayOfArrayOfArrays)
  {
    const std::string str = "1.1,2.2;1.1,2.2 |1.1,2.2;1.1,2.2% 1.1,2.2; 1.1,2.2| 1.1,2.2;1.1,2.2";
    auto arr = Core::IO::StringConverter<
        std::array<std::array<std::array<std::array<double, 2>, 2>, 2>, 2>>::parse(str);
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
    auto vec = Core::IO::StringConverter<std::vector<double>>::parse(str);
    std::vector<double> expected_vec = {1, 2, 3};
    EXPECT_EQ(vec, expected_vec);
  }

  TEST(ParseTests, ParseVectorOfVectorsOfVectors)
  {
    const std::string str = "1,2,3;4,5,6;7,8,9|1,2,3;4,5 ,6;7,8,9|1,2,3;4,5,6;7,8,9";
    auto vec = Core::IO::StringConverter<std::vector<std::vector<std::vector<double>>>>::parse(str);
    std::vector<std::vector<std::vector<double>>> expected_vec = {{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}},
        {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}, {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}};
    EXPECT_EQ(vec, expected_vec);
  }

  TEST(ParseTests, ParseMap)
  {
    const std::string str = "100 :2.2,200 :4.4";
    std::map<int, double> map = Core::IO::StringConverter<std::map<int, double>>::parse(str);
    std::map<int, double> expected_map;
    expected_map[100] = 2.2;
    expected_map[200] = 4.4;

    EXPECT_EQ(map, expected_map);
  }

  TEST(ParseTests, ParseMapArrayValues)
  {
    const std::string str = "1.1:10,20,30; 2.2:1,2,3";
    std::map<double, std::array<int, 3>> map =
        Core::IO::StringConverter<std::map<double, std::array<int, 3>>>::parse(str);
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
        Core::IO::StringConverter<std::map<std::array<int, 3>, double>>::parse(str);
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
    std::map<int, double> map = Core::IO::StringConverter<std::map<int, double>>::parse(str);
    std::map<int, double> expected_map;
    expected_map[100] = 2.2;
    expected_map[200] = 4.4;

    EXPECT_EQ(map, expected_map);
  }

  TEST(ParseTests, ParseTuple)
  {
    const std::string str = "1,2,3.3";
    std::tuple<int, int, double> tuple =
        Core::IO::StringConverter<std::tuple<int, int, double>>::parse(str);
    std::tuple<int, int, double> expected_tuple = std::make_tuple(1, 2, 3.3);
    EXPECT_EQ(tuple, expected_tuple);
  }

  TEST(ParseTests, ParseTupleWrongDimension)
  {
    EXPECT_THROW((Core::IO::StringConverter<std::tuple<int, int, double>>::parse("1,2,3.3,9.0")),
        Core::Exception);
    EXPECT_THROW(
        (Core::IO::StringConverter<std::tuple<int, int, double>>::parse("1,2")), Core::Exception);
  }

  TEST(ParseTests, ParsePair)
  {
    const std::string str = "1,2.2";
    std::pair<int, double> pair = Core::IO::StringConverter<std::pair<int, double>>::parse(str);
    std::pair<int, double> expected_pair = std::make_pair(1, 2.2);
    EXPECT_EQ(pair, expected_pair);
  }

  TEST(ParseTests, ParsePairWrongDimension)
  {
    EXPECT_THROW(
        (Core::IO::StringConverter<std::pair<int, double>>::parse("1,2.2,3.3")), Core::Exception);
    EXPECT_THROW(
        (Core::IO::StringConverter<std::pair<int, double>>::parse("1;7")), Core::Exception);
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
    auto tuple = Core::IO::StringConverter<
        std::tuple<std::map<int, std::vector<double>>, std::pair<double, std::list<char>>,
            std::vector<std::vector<std::vector<int>>>>>::parse(str);
    std::tuple<std::map<int, std::vector<double>>, std::pair<double, std::list<char>>,
        std::vector<std::vector<std::vector<int>>>>
        expected_tuple = std::make_tuple(map, pair, vec);

    EXPECT_EQ(std::get<0>(tuple), std::get<0>(expected_tuple));
    EXPECT_EQ(std::get<1>(tuple), std::get<1>(expected_tuple));
    EXPECT_EQ(std::get<2>(tuple), std::get<2>(expected_tuple));
  }
}  // namespace