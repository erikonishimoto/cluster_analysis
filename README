= cluster_analysis

2012/08/24 Eriko Nishimoto (eriko _at_ gfd-dennou.org )

== Summary
This program performs Cluster Analysis with Ward's method, 
which is one of the hierarchical clustering methods and 
defined based on the notion of square error.

In Ward's method, at each clustering step, the value of the 
sum-of-square error is computed for every possible merger of 
two clusters. The merger which produces the smallest increase 
in the value of sum-of-square error is taken to be the clustering 
of this step. Initially, each cluster contains only one object; 
hence, the value of the sum-of-square error at the beginning is zero.

== Runtime Dependency
* Ruby (>=1.8)
* NArray

== How to Use
  load this file in a ruby program

== Functions
--- cluster_analysis( data, dim, nclass, renumber )

    performing cluster analysis with Ward's method

    ARGUMENTS
    * data  : NArray whose data is calculated
    * dim   : cluster analysis is conducted along this dimention (default=0)
    * nclass: when the number of cluster becomes nclass, calculating
        is stopped (default=10)
    * renumber: renumbering the final clusters (default=true)

    RETURN VALUE
    * element_class: result of cluster analysis
    * dmin_list    : mergers among final clusters

--- euclidean_distance( data, dim=0 )

    called this in cluster_analysis for calculating the value of 
    the sum-of-square error for every possible merger of two clusters


