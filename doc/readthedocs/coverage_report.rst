.. _coveragereport:

Coverage report
-----------------

The BACI coverage report can be used to verify existing tests and to detect untested parts of the code.


You can find it by clicking on the coverage badge `here <https://baci.pages.gitlab.lrz.de/baci/coverage_report/index.html>`_
or at the BACI `frontpage <https://gitlab.lrz.de/baci/baci/>`.

The report contains the source code annotated with coverage information.
Coverage in its simplest form gives the ratio of executed lines of code over total lines of code of a testsuite run
(for more info see :ref:`here <code_coverage>`).
Coverage is also determined for individual modules, files and lines.
Code coverage is **one** quantitative measure for the quality of a code's testsuite.

The coverage report supplies the following information about the BACI testsuite:

- the lines of code that are actually executed (statement coverage)
- the number of times each line of code is executed
- the amount of executed functions (function coverage)

For easy navigation, it contains overview pages replicating the file structure.

You may find a **link to the commit** that the current report is based on at the top of each report page.

Best practices
~~~~~~~~~~~~~~~~~~

Remember, code coverage is **one** measure for code quality out of many.
As with any measure that boils complex information down to a single numeric value caution is advised as it may only reflect a small part of the truth. 

1. Coverage comes after good tests
""""""""""""""""""""""""""""""""""""""

High coverage is a result of good tests. 
The opposite is not necessarily true.  

Think for example of a test without an assert statement (a result description).
Such a test gives a relatively high coverage but almost no information about the status of the code
(except for "code is running").

Use coverage to increase the quality of tests:

#. To ensure that they execute the code you intend to test.
#. To identify untested code.

It is highly recommended to check and study the coverage report regularly.
Where coverage is missing, design a :ref:`suitable <firstprinciples>` test.
But try to avoid tests that just execute lines of code.

.. note::

    You do not need to increase coverage for the sake of increasing coverage.
    Sometimes you might also find legacy code that is untested and no longer needed.
    In this case, you could also think about removing the code.

.. note::

    Ideally, coverage is determined with :ref:`unit tests <unittesting>`,
    because unit tests make it particularly easy to relate test execution with executed lines of code.

2. Coverage should increase over time
"""""""""""""""""""""""""""""""""""""""""

An ideal code project has 100% coverage, where all existing lines of code are tested with suitable tests.

Under consideration of "[coverage comes after good tests](#1-coverage-comes-after-good-tests)" 
one development goal should therefore be to try to converge to this limit over time.
In fact, many projects allow only merge requests that do not decrease the coverage. 
With BACI this it is not that simple, because it is unfeasible to create a full coverage report for each merge request.

Nevertheless, it is considered **good practice** to check how your merge requests influence the coverage.
To do so you may study the report a few days after your changes have been merged, i.e.,
when the coverage report has been updated to include the changes.

Together with "[coverage comes after good tests](#1-coverage-comes-after-good-tests)", 
this will bring BACI closer and closer towards the 100% limit.  

To give you an idea on how good we are doing, the following estimates are used to color the coverage report:

.. role:: red
.. role:: green

- :red:`low: < 75 %`
- medium: >= 75 %  (yellow)
- :green:`high: >= 90 %`

**But: there is an exception to every rule!**

BACI is **not** an ideal code project.

Take for example ``dserror`` statements. ``dserror`` statements can only be tested with unit tests but not all code in BACI is testable with unit tests.
Therefore, ``dserror`` statements might decrease code coverage because they can add untestable lines of code.
However, adding ``dserror`` statements as safety checks for valid parameter choices in
fact increases code quality. 

.. _code_coverage:

Code coverage
~~~~~~~~~~~~~~~~~

**Statement coverage** is defined as the ratio of executed lines of code over total lines of code during a program run:

.. math::

    \textit{coverage} = \frac{\textit{executed lines of code}}{\textit{total lines of code}}

Note that we assume that each code statement equals one line of code. Therefore, statement coverage 
is often also called **line coverage**.  
There is also **function coverage**, which gives the amount of functions executed over the total 
amount of functions in the code, and **branch coverage**, which we currently do not track in BACI
(see e.g. the `section of the lecture by S.J. Zeil <https://www.cs.odu.edu/~cs252/Book/branchcov.html>`_ for an explanation).
Coverage can be measured for the whole code or per module, file or even per line. 
On an individual line basis, it breaks down to a binary measure of whether a line was executed or not during a program run.

Some technical details
~~~~~~~~~~~~~~~~~~~~~~~~~~

We use `gcov <https://gcc.gnu.org/onlinedocs/gcc/Gcov.htm>`_ to measure the coverage.
gcov is a tool that can be used together with GCC to test code coverage of programs.

We use `lcov <http://ltp.sourceforge.net/coverage/lcov.php>`_ as a graphical front-end for gcov.
lcov collects gcov data for multiple source files and creates HTML pages containing the source code annotated with coverage information.

