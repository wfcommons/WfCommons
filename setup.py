import setuptools

with open("README.md", "r") as fh:
	long_description = fh.read()

setuptools.setup(
	name="workflowhub",
	version="0.0.1",
	author="WorkflowHub team",
	author_email="lpottier@isi.edu",
	description="A package to manage scientific workflow traces",
	long_description=long_description,
	long_description_content_type="text/markdown",
	url="https://github.com/workflowhub/workflowhub",
	packages=setuptools.find_packages(),
	classifiers=[
		"Programming Language :: Python :: 3",
		"License :: OSI Approved :: MIT License",
		"Operating System :: OS Independent",
	],
	python_requires='>=3.3',
)
