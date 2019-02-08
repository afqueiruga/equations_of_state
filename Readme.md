# afq's Equations of State implementations

Alejandro F Queiruga  
Lawrence Berkeley National Lab  
2019

These are not quite effecient implementations, but I tried to make them as sanitary as possible to input arbitrary types for symbolic manipulations. I use the implementations as my one-off calculator and to build databases to *automatically search for better parameterizations*.

Currently, only water is implemented following the implementations of [IAPWS](http://www.iapws.org). The file [iapws97.py](iapws97.py) spans the gas, liquid, and supercritical regions. Ice Ih is included in the file ... A dataset spanning all of these regions *including the equilibria* is included in [water_data.csv](water_data.csv).

