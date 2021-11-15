import setuptools
import os

# Utility function to read the README file.
# Used for the long_description.  It's nice, because now 1) we have a top level
# README file and 2) it's easier to type in the README file than to put a raw
# string in below ...
def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()


setuptools.setup(
    name="yulewalker",
    version="0.1.1",
    author="Emmanouil Theofanis Chourdakis",
    author_email="emmanouil.chourdakis@nomono.co",
    description="IIR Filter design using the modified Yule-Walker method",
    classifiers=[
        "Topic :: Scientific/Engineering",
        "Development Status :: 3 - Alpha",
        "License :: OSI Approved :: CEA CNRS Inria Logiciel Libre License, version 2.1 (CeCILL-2.1)",
    ],
    packages=["yulewalker"],
    long_description=read("README.md"),
    license="CeCILL-2.1",
    keywords="yulewalk yulewalker iir signal-processing dsp",
    install_requires=[
        "numpy>=1.19.2",
        "scipy>=1.5.2",
    ],
    url="https://github.com/mmxgn/yulewalker",
)
