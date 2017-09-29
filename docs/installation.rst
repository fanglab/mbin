.. highlight:: shell

============
Installation
============


Fundamental dependencies
------------------------

.. code-block:: console

	python v2.7.*
	gcc
	hdf5

Python packages
---------------

.. code-block:: console

	numpy>=1.7.1
	pysam == 0.10.0
	h5py >= 2.0.1
	pbcore >= 0.9.4
	scipy >= 0.12.0
	biopython >= 1.6.1
	matplotlib >= 1.5.0

All but Numpy will be installed automatically during the standard installation as described below. Numpy, however, must be installed prior to mBin installation (see below):



Setting up virtualenv
---------------------

mBin should be installed in a Python virtual environment using `virtualenv`_, which creates a clean and isolated Python environment in which to install packages and their dependencies.

.. _virtualenv: https://virtualenv.pypa.io/en/stable/

Virtualenv can be installed using ``pip``:

.. code-block:: console
	
	$ pip install virtualenv

Once installed, navigate to the directory where you would like to keep the virtual environment and create a virtual environment called ``venv``:

.. code-block:: console

	$ virtualenv venv

Finally, activate this virtual environment ``venv``:

.. code-block:: console

	$ . venv/bin/activate

Once activated, you are now operating inside the ``venv`` and should see the following on you command line:

.. code-block:: console

	(venv)<COMMAND LINE>


Installing mBin
---------------

With the virtual environment activated, install numpy using ``pip``:

.. code-block:: console
	
	$ pip install numpy

Next, clone the mBin repository to your desired location:

.. code-block:: console

	$ git clone https://github.com/fanglab/mbin.git

Navigate to the cloned ``mbin`` directory and install using the ``setup.py`` file:

.. code-block:: console

	$ cd mbin
	$ python setup.py install

Testing installation
^^^^^^^^^^^^^^^^^^^^
IN PROGRESS







.. Stable release
.. --------------

.. To install mbin, run this command in your terminal:

.. .. code-block:: console

.. 	$ pip install mbin

.. This is the preferred method to install mbin, as it will always install the most recent stable release. 

.. If you don't have `pip`_ installed, this `Python installation guide`_ can guide
.. you through the process.

.. .. _pip: https://pip.pypa.io
.. .. _Python installation guide: http://docs.python-guide.org/en/latest/starting/installation/


.. From sources
.. ------------

.. The sources for mbin can be downloaded from the `Github repo`_.

.. You can either clone the public repository:

.. .. code-block:: console

.. 	$ git clone git://github.com/fanglab/mbin

.. Or download the `tarball`_:

.. .. code-block:: console

.. 	$ curl  -OL https://github.com/fanglab/mbin/tarball/master

.. Once you have a copy of the source, you can install it with:

.. .. code-block:: console

.. 	$ python setup.py install


.. .. _Github repo: https://github.com/fanglab/mbin
.. .. _tarball: https://github.com/fanglab/mbin/tarball/master
