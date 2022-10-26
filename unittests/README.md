# BACI's unit test suite

For unit testing we use the [GoogleTest](https://github.com/google/googletest) library.

**Before starting, read the GoogleTest [primer](http://google.github.io/googletest/primer.html)
sections up to section
"Invoking the Tests".
Approximate reading time 10 mins.**

Writing a new unit test roughly works like this:

1. Create a new .cpp file in the directory `unittests` in a subdirectory that matches the
   relative location of the tested file in the `src` directory. You may need to create the
   directory if it does not exist yet.
2. Add the new file to the `CMakeLists.txt` inside the directory (again you may need to create it).
   Consider the existing test directories as a template. In case you have problems listing all
   necessary
   libraries for your test, you can use
   `baci_link_google_test_necessary_libraries(${TESTNAME} -Wl,--start-group ${BACI_LIBRARIES} -Wl,--end-group)`
   as a temporary solution to include all libraries in the correct order.
3. Include the gtest header.
4. Open an anonymous namespace inside the file.
5. Write your tests. Use existing tests and the GoogleTest documentation as a guideline and
   inspiration.

## Special topics

### Testing for dserror

- In unit test executables, `dserror` is automatically replaced by a version that throws a
  `std::runtime_error`.
- You can test for errors with `EXPECT_THROW(<code_with_error>, std::runtime_error);`

### Testing code that runs in parallel

You may execute a unit test executable in parallel by
using `baci_add_google_test_executable(<name> NP <number of processes> SOURCE source1 [source2 ..
.])`. The resulting test executable will then be called with the correct mpi flags by `ctest`.