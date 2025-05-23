#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Snakemake Wrappers documentation build configuration file, created by
# sphinx-quickstart on Wed Dec 14 09:55:46 2016.
#
# This file is execfile()d with the current directory set to its
# containing dir.
#
# Note that not all possible configuration values are present in this
# autogenerated file.
#
# All configuration values have a default; values that are commented out
# serve to show the default.

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
import subprocess
from sphinxawesome_theme.postprocess import Icons

sys.path.insert(0, os.path.abspath("."))


# -- General configuration ------------------------------------------------

# If your documentation needs a minimal Sphinx version, state it here.
#
# needs_sphinx = '1.0'

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "sphinx.ext.intersphinx",
    "sphinx.ext.todo",
    "generate_docs",
    "myst_parser",
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# The suffix(es) of source filenames.
# You can specify multiple suffix as a list of string:
#
# source_suffix = ['.rst', '.md']
source_suffix = ".rst"

# The master toctree document.
master_doc = "index"

# General information about the project.
project = "Snakemake Wrappers by yangqun"
copyright = "2025, yangqun"
author = "yangqun"

# The version info for the project you're documenting, acts as replacement for
# |version| and |release|, also used in various other places throughout the
# built documents.
#
# The full version, including alpha/beta/rc tags.
release = subprocess.check_output(["git", "describe", "--all"]).strip().decode()
if "master" in release:
    release = "master"
# The short X.Y version.
version = release

# The language for content autogenerated by Sphinx. Refer to documentation
# for a list of supported languages.
#
# This is also used if you do content translation via gettext catalogs.
# Usually you set "language" from the command line for these cases.
language = None

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This patterns also effect to html_static_path and html_extra_path
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# If true, `todo` and `todoList` produce output, else they produce nothing.
todo_include_todos = True


# -- Options for HTML output ----------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
pygments_style = "sphinx"
html_theme = "sphinxawesome_theme"
html_theme_options = {
    "logo_light": "logo-snake.svg",
    "logo_dark": "logo-snake.svg",
    "main_nav_links": {
        "Homepage": "https://conanyangqun.github.io/snakemake-wrappers-yq/",
        "Plugin catalog": "https://conanyangqun.github.io/snakemake-wrappers-yq/",
        "Workflow catalog": "https://conanyangqun.github.io/snakemake-wrappers-yq/",
    },
    "awesome_external_links": True,
    "awesome_headerlinks": True,
    "show_prev_next": False,
}
html_permalinks_icon = Icons.permalinks_icon
html_css_files = ["custom.css"]


# Theme options are theme-specific and customize the look and feel of a theme
# further.  For a list of options available for each theme, see the
# documentation.
#
# html_theme_options = {}

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["_static"]

html_title = "Snakemake wrappers by yangqun"
