#!/usr/bin/env python

#from distutils.core import setup
from setuptools import setup

setup(name='BNfinder',
      version='0.2.1',
      author='Bartek Wilczynski and Norbert Dojer',
      author_email='bartek@mimuw.edu.pl',
      url='http://bioputer.mimuw.edu.pl/software/BNfinder/',
      packages=['BNfinder'],
      scripts=['bnf'],
      license='GNU GPL v.2',
      description='Bayes Net Finder: an efficient and exact method for learning bayesian networks'
     )


