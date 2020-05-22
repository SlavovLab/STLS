# Convex Structured Total Least Squares


This repository contains code for solving structured total least squares problems, i.e., problems for which both dependent and independent variables contain variable amounts of noise. The method was developed and reported by [Malioutov & Slavov (2014)](http://proceedings.mlr.press/v32/malioutov14.html).


<img src="STLS_Comparision.png.png" width="50%">


## Reproducing the analysis reported by [Malioutov & Slavov (2014)](http://proceedings.mlr.press/v32/malioutov14.html)

- The figure comparing the performance of different methods can be generated by running the code **compare_convex_STLS.m**

## Functions
- The main function is **struct_TLS_SDP_Aonly.m**
- The function **RWNN_tls.m** is a wrapper around the main function

## Dependancies
- The convex structured least squares algorithm uses several external packages listed below and included in this repository. The main function **struct_TLS_SDP_Aonly.m** assumes that these packages are in the same directory and adds them to the path. The location of these packages can be changed and the change indicated in the subfunction **addpaths** in the body of the main function **struct_TLS_SDP_Aonly.m**   
  * Yalmip
  * SDPT3-4.0
  * slra-0.5
