drugscan: Software for in silico drug repositioning, based on analysis of
Connectivity Map 2.0 data and disease gene sets. This software accompanies the
manuscript "A computational systems biology approach for identifying candidate
drugs for repositioning for cardiovascular disease", by Alvin Yu and Stephen
Ramsey, which has been submitted to the journal Interdisciplinary Life Sciences:
Computational Life Science".

Author:  Alvin Yu, Oregon State University

From the lab of:  Stephen Ramsey, Oregon State University (lab.saramsey.org)

Date:  February 17, 2016

This software is distributed under the Apache Software License 2.0.
Please see the file LICENSE for details on the software licensing
agreement.

Usage notes: There are three R scripts in the subdirectory "R", that comprise
this software release. The scripts load tab-delimited text datafiles, examples
of which are given in the subdirectory "data".  Each of the example datafiles
contains the first ten lines of the actual datafile used in the analysis. All of
the example datafiles have Unix line termination. The three R scripts are used
in this order:

(1) statistical_test.R
(2) ks_drugbank_test.R

The script "synthetic_data_analysis.R" creates the synthetic dataset used for
the analysis of optimal weight values for the gene set scoring (see Table 1 in the
above-referenced manuscript).
