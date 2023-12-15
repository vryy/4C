.. _coding-guidelines:

Coding Guidelines
==================

General Guidelines regarding Coding in C++
--------------------------------------------

Current C++ standard in BACI is C++17. Hence, use of C++17 features is encouraged.
Avoid define flags as much as possible, because they complicate testing and, thus, lead to untested code.
Avoid and actively resolve header-in-header inclusion to speed up compilation time. Use forward declarations instead.
Do not use ``using ... / typedef`` statements at unscoped level.

For more general advice on C++ coding, see for example

- this `C++ tutorial <http://www.cplusplus.com/doc/tutorial/>`_ for new developers
- the `C++ core guidelines <https://github.com/isocpp/CppCoreGuidelines/blob/master/CppCoreGuidelines.md>`_ proposed by Bjarne Stroustrup
- the blog `Modernes C++ <https://www.modernescpp.com/index.php/der-einstieg-in-modernes-c>`_
- recommendations on `how to choose the right data structure that fits best to your algorithm <https://github.com/gibsjose/cpp-cheat-sheet/blob/master/Data%20Structures%20and%20Algorithms.md>`_.
- Useful books might be [Will20]_, [Bancila20]_, [Swidzinski22]_


BACI-specific Design Guidelines
---------------------------------

**Preamble:** BACI is a legacy code.
Large parts were written pre-C++11 and are not modern C++.
In addition, code review was not widely practiced back in these days.
Thus, it is recommended to critically examine the code and techniques while working on the code base and refer to the general guidelines above.

The following guidelines are especially relevant for the current state of BACI:

- If necessary, use smart pointers for their memory management capabilities, e.g., std::shared_ptr, std::unique_ptr

    - more information on passing smart pointers
    - by default, pass parameters by const reference and only expose smart pointers if memory management is necessary

- Prefer parameter container classes over Teuchos::ParameterList for passing parameters among element routines
- Make new code const-correct and fix old code in that regard when working on it.
- Implement new features in the new structural time integration. Work towards migrating existing capabilities to the new structural time integration.
- Test your code! Refer to our documentation about Testing.


BACI-specific Naming Conventions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Namespaces are spelled in all caps: ``LINALG::``

Class names / function names / enum types use camel case, starting with a capital letter:
``ClassToDoSomething / FunctionToPerformAnOperation()``

Variables / enum values use


snake_case, i.e. all small letters separated by underscores, e.g. ``variable_with_descriptive_name``


camelCase staring with a small letter, e.g. ``variableWithDescriptiveName``



Class members end with an underscore: ``variable_``

Define flags are in all caps: ``DEBUG``

Variable names must not be just a single letter, because they are impossible to find in a global search operation.
(Exception: loop indices such as i, j, but remember that even loop indices could/should have descriptive names.)


BACI file levels
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Each file has a header flag \level. This flag should indicate how well tested and applicable for others that implementation is.

- Level 0: Fundamental code infrastructure.
- Level 1: Well-tested and well-documented code with applications in several areas with clearly defined interfaces (e.g. basic fluid or structure implementations).
- Level 2: Well-tested application-specific code (e.g. ALE-FSI implementation); good for use by others
- Level 3: Experimental code, may still contain bugs or missing features that will be fixed by later revisions


Guidelines for BACI input files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``*.dat``-file should include:

- A short summary of what is tested in header
- No unused sections
- No unused parameters (as far as possible)
- No alternative input parameters/values as comment
- No empty lines

In general a "clean" format of the file is demanded (alignment of parameters/values).

The ``*.xml``-file needs to be formatted using the /utilities/bacixmlformat script.
It calls the Open Source Python package xmlformatter 0.2.4 with the following specifications:

- nested indentation with two whitespace indents for a child element
- ``<element></element>`` is collapsed to ``<element/>``
- formatting rules as described in the documentation of the xmlformatter here


To format an ``*.xml``-file run ``./utilities/bacixmlformat <path/to/file>`` from BACI's top-level directory.

.. _firstprinciples:

F.I.R.S.T. principles for writing clean tests
---------------------------------------------

The F.I.R.S.T. principles serve as a general guideline for writing framework tests or unit tests in Baci.
Being an acronym F.I.R.S.T. stands for the following principles:

Fast
~~~~~

- developers should run tests frequently
- slow tests also slow down code development

.. Note::

    unit testing benefits from extremely fast testing times
    write meaningful framework tests while reducing the run time and number of processors of a single test to keep overall testing time and resources in a reasonable scale.

Independent/Isolated
~~~~~~~~~~~~~~~~~~~~~~

- tests should not depend on each other
- tests should pass independently of each other

.. Note::

    in Baci there are some framework tests concerning mesh generation, pre-processing, or post-processing that may depend on a specific order of execution

Repeatable
~~~~~~~~~~~~

- tests produce the same result each time
- tests should be repeatable in any configuration/environment

.. Note::

    Baci is developed by contributors distributed among several institutes working on different configurations
    (including cluster configurations)

Self-Validating
~~~~~~~~~~~~~~~~

- no manual interpretation of results
- a test fails or passes

.. Note::

    manually checking results is time consuming and prone to errors

Thorough
~~~~~~~~~~~~

    cover every use case scenario including corner/edge/boundary values
    test for illegal arguments or bad inputs, exceptions and errors

.. Note::

    following test driven developement (TDD) (refer to Robert C. Martin [Martin08]_ ) this can also be interpreted as:

    **Timely**

        - following TDD write tests just before writing code that makes them pass
        - helps designing code to be testable

    However, TDD is discussed controversial in the community!

