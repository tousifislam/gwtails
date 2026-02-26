Installation
============

From PyPI
---------

.. code-block:: bash

   pip install gwtails

From source
-----------

.. code-block:: bash

   git clone https://github.com/tousifislam/gwtails
   cd gwtails
   pip install -e .

Dependencies
------------

The following packages are required and will be installed automatically:

* `numpy <https://numpy.org/>`_
* `scipy <https://scipy.org/>`_
* `matplotlib <https://matplotlib.org/>`_
* `gwtools <https://pypi.org/project/gwtools/>`_
* `emcee <https://emcee.readthedocs.io/>`_
* `corner <https://corner.readthedocs.io/>`_

Building the documentation
--------------------------

To build the documentation locally:

.. code-block:: bash

   pip install gwtails[docs]
   cd docs
   make html

The built documentation will be in ``docs/build/html/``.
