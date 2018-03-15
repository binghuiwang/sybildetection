
######### Terms of use #########

The code is provided for research purposes only and without any warranty. The code implements the SybilSCAR algorithm. 
When using it in your research work, please cite the following paper:

Binghui Wang, Le Zhang, and Neil Zhenqiang Gong. "SybilSCAR: Sybil Detection in Online Social Networks via Local Rule based Propagation", in INFOCOM, 2017. 

For any question, please contact Binghui Wang (binghuiw@iastate.com), Le Zhang (lezhang@iastate.edu), or Neil Zhenqiang Gong (neilgong@iastate.edu).

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

-nt NUM_THREADS
NUM_THREADS is the number of threads. By default, NUM_THREADS=1.

-wg  WEIGHTED_GRAPH
WEIGHTED_GRAPH indicates whether the considered graph is weighted (WEIGHTED_GRAPH=1) or not (WEIGHTED_GRAPH=0). 
If the graph is weighted, then the weights of all edges can be user defined and are stored in the third column of the GRAPHFILE. 
Otherwise, the parameter -wei WEIGHT can be used to set a SAME weight for all edges. 
By default, WEIGHTED_GRAPH=0 and WEIGHT=1/(2*average degree).

############################

######### OUTPUT #########

-postfile POSTFILE
PRIORFILE stores the final posterior probabilities of all nodes after running SybilSCAR. The format of POSTFILE is as follows 
0 1.0
1 0.8
2 0.1
...  
It means that node 0 has a posterior probability 1.0 of being benign; node 1 has a posterior probability 0.8 of being benign; node 2 has a posterior probability 0.1 of being benign (or to say, probability 0.9 of being Sybil), etc.

###########################


######## USAGE #########

Compile: g++ sybilscar.cpp -pthread -O3 -o sybilscar

Execute: ./sybilscar -graphfile GRAPHFILE -trainfile TRAINFILE -postfile POSTFILE [-priorfile PRIORFILE] [-nt NUM_THREADS]
        [-mIter MAXITER] [-tp THETA_POS] [-tn THETA_NEG] [-tu THETA_UNL] [-wg WEIGHTED_GRAPH] [-wei WEIGHT]  

###########################

######## REMARK ############

MAXITER: Different datasets could have different MAXITERs. Usually, setting MAXITER=6 would produce a promising performance (in terms of AUC and top-interval ranking) in both synthetic graphs and real-world OSNs that we have tested. 

WEIGHT: To guarantee the convergence of SybilSCAR, it can be empirically set as  WEIGHT = 1 / (2* average degree). To speed up the computation, its value can be set larger, e.g., WEIGHT=0.51 or WEIGHT=0.6.

TRAINFILE: SybilSCAR is also applicable when only labeled benign nodes or only labeled Sybils are available. In these cases, the posterior probabilities produced by SybilSCAR can be used to rank nodes. When only labeled benign nodes are available, the second line in the TRAINFILE should be an empty line. When only labeled Sybils are available, the first line in the TRAINFILE should be an empty line. 

Prior and edge weight: The node priors and edge weights can be learnt through the following paper: 

Peng Gao, Binghui Wang, Neil Zhenqiang Gong, Sanjeev R. Kulkarni, Kurt Thomas, Prateek Mittal. "SybilFuse: Combining Local Attributes with Global Structure to Perform Robust Sybil Detection". In IEEE Conference on Communications and Network Security (CNS), 2018.

###########################
