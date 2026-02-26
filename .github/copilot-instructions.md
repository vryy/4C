**During code reviews, always verify the following:**
* **Testing:** New functionality is covered by either a unit or end-to-end test.
* **Documentation:** Docs are added for new code and updated for modified code.
* **Modern C++:** No new `Teuchos::RCP` instances are introduced.
* **CMake:** Newly added tests in `tests/list_of_tests.cmake` set `REQUIRED_DEPENDENCIES` when necessary.