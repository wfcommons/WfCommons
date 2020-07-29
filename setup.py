from setuptools import setup, find_packages

with open("README.md", "r") as fh:
	long_description = fh.read()

setup(
	name="workflowhub",
	version="0.1",
	author="WorkflowHub team",
	author_email="support@workflowhub.org",
	description="Community Framework for Enabling Scientific Workflow Research and Education",
	long_description=long_description,
	long_description_content_type="text/markdown",
	url="https://github.com/workflowhub/workflowhub",
	packages=find_packages(),
	classifiers=[
		"License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
		"Operating System :: OS Independent",
		"Programming Language :: Python :: 3",
		"Programming Language :: Python :: 3.3",
		"Programming Language :: Python :: 3.4",
		"Programming Language :: Python :: 3.5",
		"Programming Language :: Python :: 3.6",
		"Programming Language :: Python :: 3.7",
		"Programming Language :: Python :: 3.8",
		"Programming Language :: Python :: 3.9",
		"Intended Audience :: Developers",
		"Intended Audience :: Education",
		"Intended Audience :: Science/Research",
		"Natural Language :: English",
		"Topic :: Documentation :: Sphinx",
		"Topic :: System :: Distributed Computing"
	],
	python_requires='>=3.3',
)
