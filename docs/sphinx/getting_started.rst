.. _getting_started:

===============
Getting Started
===============

This page will show you how to quickly get up and running with the RAJA proxy
applications.

------------
Requirements
------------

The primary requirement for building the proxy applications is a C++ compiler
that supports C++11. Additional features and programming model backends have
extra requirements that are detailed in :doc:`advanced_config`. Specific
software versions required to configure and build:

- `CMake <https://cmake.org/>`_ 3.4 or greater.

----------------------
Downloading & Building
----------------------

This repository is hosted on `GitHub <https://github.com/LLNL/RAJA-Proxies>`_.
To clone the repo into your local working space, use the command:

.. code-block:: bash

   $ git clone --recursive https://github.com/LLNL/RAJA-Proxies.git 

The ``--recursive`` argument is used to download the repository's submodules,
`RAJA <https://github.com/LLNL/RAJA>`_ and `BLT
<https://github.com/LLNL/blt>`_.  Once this command has run, the ``master``
branch will be cloned into the ``RAJA-Proxies`` directory.


^^^^^^^^^^^^^
Building RAJA
^^^^^^^^^^^^^

RAJA uses CMake and BLT to handle builds. From the root of the repository
configuration looks like this:

.. code-block:: bash

  $ mkdir build && cd build
  $ cmake -DCMAKE_INSTALL_PREFIX=/path/to/install ../

.. warning:: Builds must be out-of-source.  BLT does not allow building in
             source directory, so you must create a build directory.

CMake will provide output about the compiler that has been detected, and
what features are discovered. Some features, like OpenMP, will be enabled
if they are discovered. For more advanced configuration options, please
see the :doc:`advanced_config`.

Once CMake has completed, the applications can be compiled using Make:

.. code-block:: bash

  $ make


------------------------------
Running the Proxy Applications
------------------------------

Each application will be built in the ``bin`` directory. The binary names will
include the version and features of each application. To run the application,
simply run the relevant executable. For example the following command would
run the RAJA version of LULESH v1.0 with the OpenMP backend: 

.. code-block:: bash

  $ ./bin/lulesh-v1.0-RAJA-OMP.exe

For more detail on the input parameters supported by each application, take a
look at :doc:`applications`.
