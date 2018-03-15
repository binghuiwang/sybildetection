
######### Terms of use #########

The code is used to compute commonly used metrics: 

Area Under the ROC Curve (AUC);
False Positive Rate (FPR); 
False Negative Rate (FNR); 
True Positive Rate (TPR); 
True Negative Rate (TNR); 
Accuracy (ACC) = (TP+TN) / (TP+TN+FP+RN);

AUC is used to evaluate ranking, while the remaining metrics are used to evaluate classification.

Note: in the code, we define benign as positive and Sybil as negative. Therefore, FPR means the fraction of Sybils that are predicted as benign, while FNR is the fraction of benign nodes that are predicted as Sybils.  


The code is provided for research purposes. When using it in your research work, please cite the following paper:

Binghui Wang, Le Zhang, and Neil Zhenqiang Gong. "SybilSCAR: Sybil Detection in Online Social Networks via Local Rule based Propagation", in INFOCOM, 2017. 

For any question, please contact Binghui Wang (binghuiw@iastate.com).

############################

######### INPUT #########

-testfile TESTFILE 
TESTFILE consists of two lines, where the first line includes a set of benign nodes, separated by space; and the second line includes a set of Sybil nodes. The format of TESTFILE is as follows
1 2 4 6 7 8 9 10 ...
0 3 5 ...
It means that "1 2 4 6 7 8 9 10 ..." are benign nodes and "0 3 5 ..." are Sybil nodes.


-postfile POSTFILE
POSTFILE stores the final posterior probabilities of all nodes after running a Sybil detecion algorithm (e.g., SybilSCAR). The format of POSTFILE is as follows 
0 1.0
1 0.8
2 0.1
...  
It means that node 0 has a posterior probability 1.0 of being benign; node 1 has a posterior probability 0.8 of being benign; node 2 has a posterior probability 0.1 of being benign (or to say, probability 0.9 of being Sybil), etc.

-t THRESHOLD

Defines classification threshold when classifying nodes to benign or Sybils. By default, THRESHOLD=0.5. 

###########################


######### OUTPUT #########

The output are the values (between 0 and 1) of those metrics: 
AUC XXX
ACC XXX
FPR XXX
TPR XXX
FNR XXX
TNR XXX

###########################


######## USAGE #########

Compile: g++ metric.cpp -O3 -o metric

Execute: ./metric -testfile TESTFILE -postfile POSTFILE

###########################
