# Configuration file for the Sphinx documentation builder.
#
# For the full list of options see:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

import os.path
import sys

# Fetch the version
exec(open('../../wfcommons/version.py').read())

sys.path.insert(0, os.path.abspath('../..'))
sys.setrecursionlimit(1500)

# -- Project information -----------------------------------------------------

project = 'WfCommons'
copyright = '2020-2026, WfCommons Team'
author = 'WfCommons Team'

# The short X.Y version
version = str(__version__)
# The full version, including alpha/beta/rc tags
release = str(__version__)

# -- General configuration ---------------------------------------------------

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.viewcode',
    'sphinx_copybutton',
    'sphinx_design',
]

templates_path = ['_templates']

source_suffix = {".rst": "restructuredtext"}

master_doc = "index"

exclude_patterns = ["_build", "Thumbs.db", ".DS_Store", "_*.rst"]

autodoc_member_order = 'bysource'

# Strip prompts when copying code blocks so pasted commands run as-is.
copybutton_prompt_text = r'\$ |>>> |\.\.\. '
copybutton_prompt_is_regexp = True

# -- Options for HTML output -------------------------------------------------

html_theme = 'furo'
html_favicon = 'favicon.png'
html_logo = 'images/wfcommons-horizontal.png'
html_title = 'WfCommons'

# The documentation uses a light background only: both furo palettes are set
# to the same light values so system/browser dark mode never produces a dark
# page, and custom.css hides the theme toggle.
_light_palette = {
    "color-brand-primary": "#1c53b0",
    "color-brand-content": "#1c53b0",
    "color-background-primary": "#ffffff",
    "color-background-secondary": "#f6f8fb",
    "color-background-hover": "#eef2f8",
    "color-background-border": "#e3e8f0",
    "color-sidebar-background": "#f6f8fb",
    "color-sidebar-background-border": "#e3e8f0",
    "color-foreground-primary": "#1a2233",
    "color-foreground-secondary": "#49556b",
    "color-foreground-muted": "#6b7689",
    "color-inline-code-background": "#f2f5fa",
    "color-highlight-on-target": "#fff8dd",
    "color-admonition-background": "#f6f8fb",
}

html_theme_options = {
    "light_css_variables": _light_palette,
    "dark_css_variables": _light_palette,
    "sidebar_hide_name": True,
}

pygments_style = "friendly"
pygments_dark_style = "friendly"

html_static_path = ['_static']
html_css_files = ['css/custom.css']
