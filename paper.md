---
title: 'blimpy: Breakthrough Listen I/O Methods for Python'
tags:
  - Python
  - astronomy
  - radio astronomy
  - technosignatures
  - SETI
authors:
  - name: Danny C. Price
    orcid: 0000-0003-2783-1608
    affiliation: "1, 2" # (Multiple affiliations must be quoted)
  - name: J. Emilio Enriquez
    orcid: 0000-0003-2516-3546
    affiliation: "1, 3"
affiliations:
  - name: Department of Astronomy,  University of California Berkeley, Berkeley CA 94720
    index: 1
  - name: Centre for Astrophysics & Supercomputing, Swinburne University of Technology, Hawthorn, VIC 3122, Australia
    index: 2
  - name: Radboud University Nijmegen, Nimegen, NL
    index: 3
date: 27 June 2019
bibliography: paper.bib
---

# Summary

The search for extraterrestrial intelligence (SETI) has historically used radio astronomy data as
the main venue to search for artificial signals of extraterrestrial origin. The Breakthrough Listen program
is the latest large scale project for the search of technosignatures,  and thanks to modern telescopes
and instrumentation, as well as significant amounts of dedicated observing time, the program
has become the largest SETI endeavour in history. This has also resulted in an unprecedented amount of 
publicly-avaiable data`[@Lebofsky:2019]`.

The ``Blimpy``--Breakthrough Listen I/O Methods for Python--package provides Python 2.7+/3.6+ utilities
for viewing and interacting with the data formats used within the Breakthrough Listen program.
This includes Sigproc filterbank (.fil) and HDF5 (.h5) files that contain dynamic spectra (aka 'waterfalls'),
and guppi raw (.raw) files that contain voltage-level data. Python methods for data extraction,
calibration, and visualization are provided. A suite of command-line utilities are also available.

``Blimpy`` was designed to be used by both radio astromers, students and anyone interested in accessing
 Breakthrough Listen data. It has already been used in a number of scientific publications
`[@Croft:2016; @Enriquez:2017; @Enriquez2018; @Enriquez:2019; @Gajjar:2018; @Price:2019a;  @Price:2019b]`.


# References






