Writing Docs
============

This is a rough guide to contributing to blimpy documentation.

TLDR
----
* Write docs in ``.rst`` or ``.md``
* Add the name of the new docs file to ``index.rst``
* Preview docs by installing ``sphinx`` via ``pip`` and run ``make html`` in the ``docs/`` directory
* Site will update automatically after pushing to the blimpy repo

Creating a Page
--------------------
Currently, ``readthedocs`` is able to process two kinds of files: ``reStructuredText (.rst)`` and ``Markdown (.md)``.
You can find a brief guide for reStructuredText `here <http://www.sphinx-doc.org/en/1.8/usage/restructuredtext/basics.html>`_
and a guide for Markdown `here <https://github.com/adam-p/markdown-here/wiki/Markdown-Cheatsheet>`_.

The two file types are rendered equally by sphinx, so feel free to use whichever one you're more comfortable with.

To create a new page, you can create a new file in the same directory as ``index.rst``. After creating the file,
add the filename of the new file to ``index.rst`` to add a link it on the index page. The new file will also show up in the sidebar
after this.

Previewing Docs
---------------
The docs are rendered using ``sphinx``, a python package. Use ``pip install sphinx`` to install the package.

After ``sphinx`` is installed, you can preview your changes by running ``make html`` in the ``docs/`` directory.
The rendered html files will be stored in ``docs/_build/html``. The actual site will look exactly like the rendered
files when built.

Automatic Documentation
-----------------------
You can run ``sphinx-apidoc -o . ../blimpy/ -f`` in ``blimpy/docs`` to generate autodoc pages from all the python modules in blimpy.
Make sure to run this command every time a new file is added to blimpy.

Updating the Site
-----------------
The blimpy Github repo is connected to `readthedocs.org` with a webhook. `readthedocs` will automatically update the site
whenever a new commit is added to the repo.
