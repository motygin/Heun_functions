# Heun_functions
Matlab/Octave code for evaluation of the Heun functions

The basic functions to use:

HeunL0:  local Heun function, equal to 1 at z=0

HeunLmv:  multi-valued function, integration along an arbitrary path connecting 0 and z

HeunS0:  second local solution at z=0

HeunSmv:  multi-valued function, integration along an arbitrary path connecting 0 and z

HeunL, HeunS: with improvement near singularities z=1, a, infinity

HeunOpts: tuning of internal parameters 

Auxiliary functions:
HeunGfromZ0
HeunLS
HeunG_near_1
HeunG_near_a
HeunG_near_infty
HeunL00
HeunL00log
HeunS0gamma1
HeunS00gamma1
HeunSmvgamma1

tested in the Octave numerical environment,
intended to be Matlab compatible 

Ref.: http://arxiv.org/abs/1506.03848


