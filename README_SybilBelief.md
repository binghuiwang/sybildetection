
######### Terms of use #########

The code is provided for research purposes only and without any warranty. Any commercial use is prohibited.
The code implements the SybilBelief algorithm. 
When using it in your research work, you should cite the following paper:

Neil Zhenqiang Gong, Mario Frank, and Prateek Mittal. "SybilBelief: A Semi-supervised Learning Approach for Structure-based Sybil Detection". In IEEE Transactions on Information Forensics and Security (TIFS), 9(6), 2014. 

For any question, please contact Neil Zhenqiang Gong (neilgong@iastate.edu), Le Zhang (lezhang@iastate.edu), or Binghui Wang (binghuiw@iastate.edu).
############################

######### INPUT #########

-graphfile GRAPHFILE
GRAPHFILE stores the edges and weights of an undirected social graph. The format of GRAPHFILE is as follows

0 1 0.8

0 2 0.6

... 

1 0 0.8

...

2 0 0.6

...

It means that node 0 and node 1 are connected with edge weight 0.8, etc.

Note that each edge in the GRAPHFILE appears twice, e.g., 0 1 0.8 and 1 0 0.8, and nodes are consecutive integers starting from 0. 


-priorfile PRIORFILE

PRIORFILE stores the prior probabilities of all nodes of being benign. The format of PRIORFILE is as follows

0 0.9

1 0.5

2 0.1

3 0.8

...  

It means that node 0 has a prior probability 0.9 of being benign; node 2 has a prior probability 0.1 of being benign (or to say, prior probability 0.9 of being Sybil), etc. 

Note that these prior probabilities can be user-defined for labeled benign nodes or/and labled Sybil nodes or/and unlabeled nodes. Or, they can be also learnt via a machine learning classifier. For example, we can extract local node features to train a binary classifier, which produces the probability of being benign for each node. Then, such probabilities can be treated as nodes' priors.


-trainfile TRAINFILE 

TRAINFILE consists of two lines, where the first line includes a set of labeled benign nodes, separated by space; and the second line includes a set of Sybil nodes. The format of TRAINFILE is as follows

1 2 4 6 7 8 9 10 ...

0 3 5 ...

It means that "1 2 4 6 7 8 9 10 ..." are labeled benign nodes and "0 3 5 ..." are labeled Sybil nodes.


-mIter MAXITER 

MAXITER sets the number of maximum iterations. By default, MAXITER=5.


-tp THETA_POS

THETA_POS sets the prior probability of being benign for labeled benign nodes. By default, THETA_POS=0.9.


-tn THETA_NEG

THETA_NEG sets the prior probability of being benign for labeled Sybil nodes. By default, THETA_NEG=0.1.


-tu THETA_UNL

THETA_UNL sets the prior probability of unlabeled nodes. By default, THETA_UNL=0.5.


-wg  WEIGHTED_GRAPH

WEIGHTED_GRAPH indicates whether the considered graph is weighted (WEIGHTED_GRAPH=1) or not (WEIGHTED_GRAPH=0). 

If the graph is weighted, then the weights of all edges can be user defined and are stored in the third column of the GRAPHFILE. 

Otherwise, if it is unweighted, then the parameter

-wei WEIGHT

can be used to set the SAME weight for all edges. 

By default, WEIGHTED_GRAPH=0 and WEIGHT=0.9. 


######### OUTPUT #########

-postfile POSTFILE

PRIORFILE stores the final posterior probabilities of all nodes after running SybilBelief. The format of POSTFILE is as follows 

0 1.0

1 0.8

2 0.1

...  

It means that node 0 has a posterior probability 1.0 of being benign; node 1 has a posterior probability 0.8 of being benign; node 2 has a posterior probability 0.1 of being benign (or to say, probability 0.9 of being Sybil), etc.
###########################


######## USAGE #########
Compile: g++ sybilbelief.cpp -O3 -o sybilbelief

Execute: ./sybilbelief -graphfile GRAPHFILE -trainfile TRAINFILE -postfile POSTFILE [-priorfile PRIORFILE] 
        [-mIter MAXITER] [-tp THETA_POS] [-tn THETA_NEG] [-tu THETA_UNL] [-wg WG] [-wei WEI]  
###########################

