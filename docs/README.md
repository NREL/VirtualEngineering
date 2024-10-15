# Virtual Engineering Documentation

## Organization

The documentation is broken up into three separate folders:

1. How-To Guides: These are walkthroughs intended for experienced programmers to help understand key steps in the setup and simulation process.
2. Background: Information that provides theory and background for unit models and other topics where appropriate.
3. Technical Reference: Programming details for the VE functions, principally this is auto-generated from module docstrings.

## Building Documentation

For developers interested in building the documentation locally for testing purposes, carry out the following steps:

1. In order to build the documentation, you will need a version of Sphinx and some other packages installed in your VE Conda environment. You can do it by running `pip install -r ./docs/requirements.txt`. 
2. Change directory to `VirtualEngineering/docs`
3. Run `make html`
4. Open the file `_build/html/index.html`

The philisophical organization of the documentation is shown above.  From a file system perspective, this means storing documentation using the following directory structure:

```
docs
|  README.md
|  index.rst
|  conf.py
|  make.bat
|  Makefile
|  ...
|
|--background
|  |  index.rst
|
|--how_to_guides
|  |  index.rst
|  |  guide_1.rst
|  |  guide_2.rst
|  |  ...
|
|--technical_reference
|  |  index.rst
|  |  module_1.rst
|  |  module_1.rst
|  |  ...
|
|--_build
|  |  doctrees
|  |  html
|  |  |  index.html
|  |  |  ...
|

```

where the `_build` directory and the resulting html files are the result of Step 3.  Each `index.rst` file enumerates the guides/pages in that section, while the parent `index.rst` file at the `docs/` level contains information to be shown on the home page and a higher-level presentation of the individual sections.
