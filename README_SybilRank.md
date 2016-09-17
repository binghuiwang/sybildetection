
######### Terms of use #########

The code is provided for research purposes only and without any warranty. Any commercial use is prohibited.
The code implements the SybilRank algorithm and its multi-thread version.

When using the code in your research work, you should cite the following paper:

Neil Zhenqiang Gong, Mario Frank, and Prateek Mittal. "SybilBelief: A Semi-supervised Learning Approach for Structure-based Sybil Detection". In IEEE Transactions on Information Forensics and Security (TIFS), 9(6), 2014. 

For any question, please contact Neil Zhenqiang Gong (neilgong@iastate.edu), Le Zhang (lezhang@iastate.edu), or Binghui Wang (binghuiw@iastate.edu).
############################

######### INPUT #########

-graphfile GRAPHFILE

GRAPHFILE stores the edges of an undirected social graph. The format of GRAPHFILE is as follows

0 1

0 2

... 

1 0

...

2 0

...

It means that node 0 and node 1 are connected; and node 0 and 2 are connected, etc.

Note that each edge in the GRAPHFILE appears twice, e.g., 0 1 and 1 0, and nodes are consecutive integers starting from 0. 


-priorfile PRIORFILE

PRIORFILE stores the trust scores of all nodes of being benign. The format of PRIORFILE is as follows

0 0.9

1 0.5

2 0.1

...  

It means that node 0 has a trust score 0.9 of being benign; node 1 trust score 0.5 of being benign; node 2 trust score 0.1 of being benign, etc.

Note that these trust scores can be user-defined for labeled benign nodes (since SybilRank can only leverage labeled benign nodes). Or, they can be also learnt via a machine learning classifier. For example, we can extract local node features to train a binary classifier, which produces the probability of being benign for each node. Then, such probabilities can be treated as nodes' trust scores.


-trainfile TRAINFILE 

TRAINFILE consists of two lines, where the first line includes a set of labeled benign nodes, separated by space; and the second line includes a set of Sybil nodes. The format of TRAINFILE is as follows

1 2 4 6 7 8 9 10 ...

0 3 5 ...

It means that "1 2 4 6 7 8 9 10 ..." are labeled benign nodes and "0 3 5 ..." are labeled Sybil nodes.

Note that SybilRank only uses labeled benign nodes. For convenient comparison with other algorithms, we also include Sybil nodes in the TRAINFILE. 


-mIter MAXITER 

MAXITER sets the number of iterations. By default, MAXITER=log2(N), where N is the number of nodes.


-alpha ALPHA 

ALPHA is the restart probability, whose value reflects to what extent the influence of labeled benign nodes. 
For original SybilRank, ALPHA=0. However, one can define different ALPHA considering different scenarios.

#############
Additional argument for multi-thread Sybilrank

-nt NT

NT is the number of threads. By default, it is 1.
#############

######### OUTPUT #########

-postfile POSTFILE

PRIORFILE stores the final scores of all nodes of benign benign after running SybilRank. The format of POSTFILE is as follows

0 0.05

1 0.03

2 0.01

...  

It means that node 0 has a score 0.05 of being benign; node 1 has a score 0.03 of being benign, etc.
###########################

######## USAGE #########

Compile: g++ sybilrank.cpp -O3 -o sybilrank

	g++ sybilrank_multi_thread.cpp -pthread -O3 -o sybilrank

Execute: ./sybilrank -graphfile GRAPHFILE -trainfile TRAINFILE -postfile POSTFILE 
		[-priorfile PRIORFILE] [-mIter MAXITER] [-alpha ALPHA]  
	
	./sybilrank_multi_thread -graphfile GRAPHFILE -trainfile TRAINFILE -postfile POSTFILE
		[-priorfile PRIORFILE] [-mIter MAXITER] [-alpha ALPHA] [-nt NT] 
###########################

