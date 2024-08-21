Welcome to Virtual Engineering's documentation!
===============================================

.. The Virtual Engineering (VE) project provides a set of Python modules to enable users to connect mathematical models of unit operations to predict outcomes for an entire chemical process.  This VE approach is currently being developed to support the beginning-to-end simulation of the low-temperature conversion of lignocellulosic biomass to a fuel precursor.

Virtual Engineering (VE) is a Python software framework that accelerates research and development of engineering processes. 
It supports multi-physics models of unit operations and joins them to simulate the entire end-to-end process. 
To automate the execution of this model sequence, `VE` provides:

* a robust method to communicate between models, 
* a high-level user-friendly interface to set model parameters and enable optimization, and 
* an overall model-agnostic approach to enable new computational units to be swapped in and out of workflows. 

Although the `VE` approach was developed to support a low-temperature conversion of biomass to fuel workflow, 
we have designed each component to easily support new domains and unit models.

Organization
------------

Documentation is currently organized into three main categories:

* :ref:`How to Guides`: User guides covering basic topics and use cases for the VE software
* :ref:`Technical Reference`: Programming details on the VE API and functions
* :ref:`Applications: Bioconversion`: Information and research sources for the low-temperature conversion unit operations

New users may find it helpful to review the :ref:`Getting Started` materials first.

Contents
--------

.. toctree::
   :maxdepth: 2

   how_to_guides/index
   technical_reference/index
   applications/bioconversion/index
   

