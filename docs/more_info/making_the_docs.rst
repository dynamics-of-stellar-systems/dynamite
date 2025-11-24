.. _making_the_docs:

*************************
Making the Documentation
*************************

This documentation is made using the package `Sphinx <https://www.sphinx-doc.org/>`_ and the associated packages `nbsphinx <https://nbsphinx.readthedocs.io/>`_, `pandoc <https://pandoc.org/>`_ and `autodocsumm <https://pypi.org/project/autodocsumm/>`_. To build this documentation
yourself, you will need to first install these. If you're using conda, for example, this can be done as follows (depends on your existing setup)::

  conda install -c conda-forge pandoc nbsphinx sphinx autodocsumm

If you are using pip (depending on your Python installation, IPython and ipykernel may not be necessary)::

  python -m pip install -U sphinx nbsphinx numpydoc IPython autodocsumm ipykernel

To make the documentation, change to the ``docs`` directory and run the command::

  make html

to make the HTML version, or::

  make latex

for the PDF.

In case Sphinx reports an error related to ``pandoc`` (containing ``Pandoc wasn't found.``), the ``pandoc`` executable is missing in its Python package.
Please install ``pandoc`` using your package manager and try again (on Mac ``brew install pandoc`` or ``sudo port install pandoc``).
