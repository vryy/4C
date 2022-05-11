# BACI's unit test suite

For unit testing we use the [GoogleTest](https://github.com/google/googletest) library.

## Migrating unit tests from CxxTest

**Before starting, read the GoogleTest [primer](http://google.github.io/googletest/primer.html)
sections up to section
"Invoking the Tests".
Approximate reading time 10 mins.**

1. Create a new .cpp file in the new directory `unittests` in a subdirectory that matches the
   relative location of the tested file in the `src` directory. You may need to create it.
2. Add the new file to the `CMakeLists.txt` inside the directory (again you may need to create it).
   Consider the existing test directories as a template.
3. Include the gtest header.
4. Open an anonymous namespace inside the file.
5. Copy the content from the old file into the anonymous namespace
6. CxxTest always required a test fixture class (the one we derived from BACICxxTestWrapper). If
   the fixture doesn't contain any setup, you can remove it since GoogleTest doesn't need one.
   Only use a fixture if there is some shared setup (see
   the [primer](http://google.github.io/googletest/primer.html) in item 0.).
7. Replace the functions `Test...` with the GoogleTest macro `TEST` or `TEST_F` if you need a
   test fixture (see the [primer](http://google.github.io/googletest/primer.html)). Note that
   the tests must reside _outside_ of the fixture class!
8. Replace the `TS_` assertions with the respective expectations `EXPECT_` from GoogleTest. See
   the [list of assertions](http://google.github.io/googletest/reference/assertions.html).

## Special topics

### Testing for dserror

- In unit test executables, `dserror` is automatically replaced by a version that throws a
  `std::runtime_error`.
- You can test for errors with `EXPECT_THROW(<code_with_error>, std::runtime_error);`

### Testing code that runs in parallel

- This was not possible with the old CxxTest setup and should not be encountered when porting tests.
- If you think a piece of code should run in parallel contact the BACI community and we can
  provide an executable for it.