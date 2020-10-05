.. _making_the_docs:

*************************
Making the Documentation
*************************

This documentation is made using the package `Sphinx <https://www.sphinx-doc.org/>`_ and the associated packages `nbsphinx <https://nbsphinx.readthedocs.io/>`_ and `pandoc <https://pandoc.org/>`_. To build this documentation
yourself, you will need to first install these. If you're using conda, for example, this can be done as follows::

  conda install -c conda-forge pandoc nbsphinx sphinx

Then to make the  documentation run, in the ``docs`` directory run the command::

  make html

to make the HTML version, or::

  make latex

for the PDF.
