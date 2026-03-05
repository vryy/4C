**During code reviews, always verify the following:**
* **Testing:** New functionality is covered by either a unit or end-to-end test.
* **Documentation:** Docs are added for new code and updated for modified code.
* **Modern C++:** No new `Teuchos::RCP` instances are introduced.
* **CMake:** Newly added tests in `tests/list_of_tests.cmake` set `REQUIRED_DEPENDENCIES` when necessary.
* **fourcipp:** Any changes made in `utilities/four_c_python/src/four_c_metadata` should be tested for compatibility with `fourcipp`. Some changes might need to be ported to `fourcipp`.