# Distributed Trace Ratio (DTRO) optimization

This repository contains the algorithms implementations and examples to solve the TRO problem in a distributed fashion [1], i.e., the following problem:

            minimize       trace(X'AX)/trace(X'BX)
            subject to     X'X=I.


# *Example_DTRO*

The folder *Example_DTRO* contains the implementations to run the TI-DTRO algorithm.  

* `create_data.m`: Functions to create synthetic signals following the models described in [1].  
* `make_sym.m`: Function to force a matrix to be symmetric.  
* `trace_ratio.m`: Function to solve the TRO problem in a centralized way. Results are used as ground truth.  
* `ti_dtro.m`: Function implementing the TI-DTRO algorithm. *Note: This algorithm works for any connected graph, including trees and fully-connected networks.*  
* `shortest_path.m`: Function to compute the shortest paths between a source node and other nodes in the network using Dijkstra's method. *Note: This function can be replaced by other algorithms pruning the graph to a tree.*  
* `example_run.m`: Script to run the TI-DTRO algorithm. Plots the MSE.  


# *dtro_full*

The folder *dtro_full* contains the codes to obtain the graphs shown in [1]. More details can be found within the folder.

## References ##

[1] C. A. Musluoglu and A. Bertrand "Distributed Adaptive Trace Ratio Optimization in Wireless Sensor Networks", Internal Report KU Leuven, 2020.
