
/*
 This code was written by Neil Zhenqiang Gong, in September 2012.

 For any question, please contact Neil Zhenqiang Gong (neilz.gong@gmail.com)
 
 */

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <set>
#include <map>
#include <list>
#include <vector>
#include <assert.h>
#include <string.h>
#include <algorithm>
#include <math.h>
#include <sstream>
#include <iostream>
#include <time.h>

#define GCC_VERSION (__GNUC__ * 10000 \
+ __GNUC_MINOR__ * 100 + __GNUC_PATCHLEVEL__)

#if GCC_VERSION >= 40300
#include <tr1/unordered_map>
using namespace std::tr1;
#define hash_map unordered_map

#else
#include <unordered_map>
#endif
using namespace std;

// random generator function:
ptrdiff_t myrandom (ptrdiff_t i) { return rand()%i;}

// pointer object to it:
ptrdiff_t (*p_myrandom)(ptrdiff_t) = myrandom;


bool my_comparator ( const pair<double, int>& l, const pair<double, int>& r)
{ 
	return l.first < r.first; 
}

double* rankdata(vector<double> predictions){
	//tiedRank of predictions
	vector< pair<double, int> > new_predictions;
	
	vector<double>::iterator iter;
	int i=0;
	for (iter=predictions.begin(), i = 0; iter != predictions.end(); iter++, i++) {
		new_predictions.push_back(make_pair(*iter, i));
	}
	
	sort(new_predictions.begin(), new_predictions.end(), my_comparator);
	
	int *ivec;
	double *svec;
	
	ivec = (int*) malloc(sizeof(int) * new_predictions.size());
	svec = (double*) malloc(sizeof(double) * new_predictions.size());
	
	
	for (i=0; i < new_predictions.size(); i++) {
		ivec[i] = new_predictions[i].second;
		svec[i] = new_predictions[i].first;
	}
	
	double sumranks = 0, dupcount = 0, averank = 0;
	int n = new_predictions.size();
	
	double *newarray;
	newarray = (double*) malloc(sizeof(double) * new_predictions.size());
	
	for(i = 0; i < n; i++) {
		sumranks += i;
        dupcount += 1;
		if (i == n - 1 || svec[i] != svec[i + 1]){
			averank = sumranks / dupcount + 1;
			for (int j = i - dupcount + 1; j < i + 1; j++) {
				newarray[ivec[j]] = averank;
			}
			sumranks = 0;
            dupcount = 0;
		}
	}
	
	free(ivec);
	free(svec);
	
	return newarray;
    
}

double AUC(vector<double> predictions, vector<int> labels){
    
	double auc = 0, num_P = 0, num_N = 0;
	
	for (int i = 0; i < labels.size(); i++) {
		if (labels[i] == 1){
			num_P += 1;
		}
		else{
			num_N += 1;
		}
	}
	
	double *tiedRank;
	tiedRank = rankdata(predictions);
	double sum_rank = 0.0;
	
	for (int i = 0; i < predictions.size(); i++) {
		if (labels[i] == 1){
			sum_rank += tiedRank[i];
		}
	}
	auc = (sum_rank - (num_P * num_P + num_P) / 2.0) / (num_P * num_N);
	
	free(tiedRank);
	
	return auc;
}



long randint(long min, long max){
    
    return long(rand() / (RAND_MAX + 0.0) * (max - min) + min);
    
}




class Data{
public:
	typedef unsigned long int vertex;
	
   	
	//scores
	unordered_map<vertex, double> post;
	
	//test positive examples, i.e., honest nodes
	set<vertex> pos_test_set;
	
	//test negative examples, i.e., sybil nodes
	set<vertex> neg_test_set;
	
    //test data file. 
	//first line is the ids of positive test nodes
	//second line is the ids of negative test nodes 
    char *test_set_file;
	
	//score file
	char *post_file;
	
    //threshold for classification
	double thresh;
    double thresh_pos;
    double thresh_neg;
    int type;
    
    
	Data(){
		
	}
    
	
	
	/* Read in scores */
	void read_scores(){
		ifstream in(post_file,ifstream::in);
		assert(in);
		
		string line;
		vertex node;
		double score;
        
		while(getline(in,line)!=NULL){
			node=(vertex)atol(strtok((char *)line.c_str()," \n\t\r"));
			score=(double)atof(strtok(NULL," \n\t\r"));
			post[node] = score;
		}
		
		in.close();
        
	}
	
	void read_test_data(){
        ifstream in(test_set_file,ifstream::in);
        assert(in);
        
        string line;
        getline(in,line);
        istringstream pos_test_str(line);
        vertex sub;
        while (pos_test_str){
            pos_test_str >> sub;
            pos_test_set.insert(sub);
        }
        
        getline(in,line);
        istringstream neg_test_str(line);
        while (neg_test_str){
            neg_test_str >> sub;
            neg_test_set.insert(sub);
        }		
        
        in.close();		
    }
    
    
    unordered_map<string, double> test_error(){
        
        
        unordered_map<string, double> result;
        
        vector<int> label_vec;
        
        vector<double> predictions;
        set<vertex>::iterator iter;
        
        for (iter = pos_test_set.begin(); iter != pos_test_set.end(); iter++) {
            label_vec.push_back(1);
            predictions.push_back(post[*iter]);
        }
        
        
        for (iter = neg_test_set.begin(); iter != neg_test_set.end(); iter++) {
            label_vec.push_back(0);
            predictions.push_back(post[*iter]);
        }
        
        double auc;
        auc = AUC(predictions, label_vec);
        result["AUC"] = auc;
        
        //confusion matrix
        double fp = 0;
        double tp = 0;
        double fn = 0;
        double tn = 0;        
        double eq_pos = 0;
        double eq_neg = 0;

        // Regard score=0.5 as half tp/tn to half fn/fp.
        for (iter = pos_test_set.begin(); iter != pos_test_set.end(); iter++) {
            if (post[*iter] > thresh) {
                tp += 1;
            }
            if(post[*iter] == thresh) {
                eq_pos += 1;
            }
            if(post[*iter] < thresh) {
                fn += 1;
            }
            
        }
        tp += int(eq_pos/2);
        fn += eq_pos - int(eq_pos/2);

        for (iter = neg_test_set.begin(); iter != neg_test_set.end(); iter++) {
            if (post[*iter] > thresh) {
                fp += 1;
            }
            if (post[*iter] == thresh) {
                eq_neg += 1;
            }
            if (post[*iter] < thresh) {
                tn += 1;
            }
            
        }
        tn += int(eq_neg/2);
        fp += eq_neg - int(eq_neg/2);
        
        result["ACC"] =(tp + tn) / (neg_test_set.size() + pos_test_set.size());	
        result["FPR"] = fp / (fp + tn);
        result["TPR"] = tp / (fn + tp);
        result["FNR"] = fn / (fn + tp);
        result["TNR"] = tn / (fp + tn);
        return result;
        
    }
    
    
    void parse_par(int argc, char **argv){
        
        //default setting
        test_set_file = "test.txt";
        post_file = "post.txt";
        thresh = 0.5;
        thresh_pos = 0.5;
        thresh_neg = 0.5;
        type = 0;
        
        int i = 1;
        
        while (i < argc) {
            if (strcmp(argv[i],"-t") == 0){
                thresh = (double)atof(argv[i + 1]);
            }
            else if (strcmp(argv[i],"-postfile") == 0){
                post_file = argv[i + 1];
            }
            else if (strcmp(argv[i],"-testfile") == 0){
                test_set_file = argv[i + 1];
            }
            
            else{
                cout << "undefined inputs: " << argv[i] <<endl;
                exit(0);
            }
            
            i += 2;
        }
        
    }
    
};


int main (int argc, char **argv)
{	
    
    srand ( time(NULL) );
    unordered_map<string, double> result;
    
    Data data;
    
    data.parse_par(argc, argv);
    
    data.read_scores();
    
    data.read_test_data();
    
    result = data.test_error();
    
    unordered_map<string, double>::iterator iter;
    for (iter = result.begin(); iter != result.end(); iter++) {
        cout << iter->first << " " << iter->second << endl;
    }
    
    return 0;
    
}
