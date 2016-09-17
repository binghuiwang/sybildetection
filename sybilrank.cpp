
/*
 This code was written by Neil Zhenqiang Gong, in September, 2012, and modified by Binghui Wang in August, 2016. 

 It implements SybilRank algorithm.

 For any question, please contact Neil Zhenqiang Gong (neilz.gong@gmail.com) or Binghui Wang (binghui.wang89@gmail.com)
 
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
#include <iomanip>
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


class SparseMatrix{
public:
    
    long num_row;
    long num_col;
    long num_non_zero;
    long *csr_row;
    long *csr_col;
    double * csr_value;
    
    
    SparseMatrix(){
        num_row = 0;
        num_col = 0;
        num_non_zero = 0;
    }

    //It is fine to not worry about this function for now. Understand the function "construct_transition_matrix". 
    //It helps you know how to transform a matrix into a sparse representation
    void ConSparseMatrix(long *I, long *J, double *V, long num_row_init, long num_col_init, long num_non_zero_init){
        //I must be sorted by row number.
        
        num_row = num_row_init;
        num_col = num_col_init;
        num_non_zero = num_non_zero_init;
        
        csr_row = (long*)malloc(sizeof(long)* (num_row + 1));
        csr_col = (long*)malloc(sizeof(long)* num_non_zero);
        csr_value = (double*)malloc(sizeof(double)* num_non_zero);
        
        long cur_row = 0, k =0, r;
        
        csr_row[0] = 0;
        csr_row[num_row] = num_non_zero;
        
        while (k < num_non_zero) {
            if (I[k] == cur_row){
                k++;
            }
            else if (I[k] > cur_row){
                for (r  = cur_row + 1; r < I[k] + 1; r++){
                    csr_row[r] = k;
                }
                cur_row = I[k];
                k++;
            }
            else{
                cout << "I is not in ascending order" << endl;
            }
            
        }
        
        for (r = cur_row + 1; r < num_row + 1; r++){
            csr_row[r] = num_non_zero;
        }
        
        for( k = 0; k < num_non_zero; k++){
            csr_col[k] = J[k];
            csr_value[k] = V[k];
        }
        
    }
    
    void matrix_product(double *x, double *y){
        //sparse matrix multiplication. y = A*x
        long i, j;
        for(i = 0; i < num_row; i++)
            y[i] = 0.0;
        for(i = 0; i < num_row; i++)
        {
            for(j = csr_row[i]; j < csr_row[i + 1]; j++)
            {
                y[i] += x[csr_col[j]] * csr_value[j];
            }
        }
        
    }
    
};


class Data{
public:
	typedef unsigned long int vertex;
	
    SparseMatrix matrix;
        
    //number of nodes
	long N;
	
	//adjacency list. this is a dictionary. 
    //adjacency list is a popular way to store graphs
	unordered_map<vertex,list<vertex> > network_map;
		
			
	//train positive examples, labeled benign/honest nodes,
	set<vertex> pos_train_set;
		
	//train negative examples, labeled Sybil nodes. 
    //Not used for SybilRank
	set<vertex> neg_train_set;
	    
    //network file
	char *network_file;
    
    //train set file
    char *train_set_file;
	
    //prior score file
    char *prior_file;
    
    //final score file
    char *post_file;
    
    //scores
	double *post;
    double *prior;
	
	//these are damping factors of the random walk;
	double alpha;
    
    //maximum number of iterations
    int max_iter;
    	
	Data(){
		
	}
	
    /* Store edges and weights in the network_map */
	void add_edge(vertex node1, vertex node2){
		// add edge (node1, node2)
		
		// no self loops
		if(node1==node2){
			return;
		}
		
		//add node2 to the adjacency list of node1
		network_map[node1].push_back(node2);
	}
	
	/* Read the social graph */
    //the format for the social graph is
    //each line corresponds to an edge, e.g, 3 2
    //each edge in the graph appears twice, e.g., 
    //3 2
    //2 3
	void read_network(){
		ifstream in(network_file,ifstream::in);
		assert(in);
		
		string line;
		vertex node1,node2;
		
        //read edges
        while(getline(in,line)!=NULL){
			node1=(vertex)atol(strtok((char *)line.c_str()," \n\t\r"));
			node2=(vertex)atol(strtok(NULL," \n\t\r"));
			add_edge(node1, node2);
		}
		
        //number of nodes in the graph
        N = network_map.size();
        
        //allocate space for final scores
		post = (double*) malloc(sizeof(double)*(N));
        
        //allocate space for final scores
        prior = (double*) malloc(sizeof(double)*(N));
        
        in.close();
        
    }
    

    /* Read priorfile or/and training file */
    void read_prior(){
        
        //initialize priors as theta_unl
        vertex node;
        for (node = 0; node < N; node++) {
            prior[node] = 0;
        }
        
        if (prior_file != "") {
            
            ifstream in(prior_file,ifstream::in);
            assert(in);
            
            string line;
            
            double score;
            while(getline(in,line)!=NULL){
                node=(vertex)atol(strtok((char *)line.c_str()," \n\t\r"));
                score=(double)atof(strtok(NULL," \n\t\r"));
                prior[node] = score;
            }
            in.close();
        }
        

        //the training dataset file includes two lines.
        //the first line includes a set of labeled benign nodes. 
        //the second line includes a set of labeled Sybil nodes. 
        //SybilRank only uses labeled benign nodes. 
        //For convenient comparison with other algorithms, we also include Sybil nodes in the training file. 
        if (train_set_file != ""){
            
            ifstream in(train_set_file,ifstream::in);
            assert(in);
            
            string line;
            
            //reading labeled benign nodes.
            getline(in,line);
            istringstream pos_train_str(line);
            vertex sub;
            while (pos_train_str){
                pos_train_str >> sub;
                pos_train_set.insert(sub);
            }
            
            set<vertex>::iterator iter;
            for (iter = pos_train_set.begin(); iter != pos_train_set.end(); iter++) {
                prior[*iter] = 1.0;
            }
            
            in.close();
            
        }
        
    }

    
	/* Construct sparse transition matrix for efficient power iteration */
    void construct_transition_matrix(){
        
        //a sparse matrix can be represented as (row1, column1, value1), (row2, column2, value2)
        //I : store the rows of these triples. It should be a sorted array
        //J : store the columns of the these triples.
        //V : store the values
        long *I;
        long *J, cur_row;
        double *V;
        
        //number of nonzero entries in a sparse matrix
        long num_non_zero;
        
        list<vertex>::iterator iter;
        
        num_non_zero = 0;
        
        for (cur_row = 0; cur_row < network_map.size(); cur_row++) {
            num_non_zero += network_map[cur_row].size();
        }
        
        //allocate space
        
        I = (long*) malloc(sizeof(long) * num_non_zero);
        J = (long*) malloc(sizeof(long) * num_non_zero);
        V = (double*) malloc(sizeof(double) * num_non_zero);
        long k = 0;
        
        for (cur_row = 0; cur_row < network_map.size(); cur_row++) {
            for (iter = network_map[cur_row].begin(); iter != network_map[cur_row].end(); iter++) {
                I[k] = cur_row;
                J[k] = *iter;
                V[k] = 1.0 / network_map[*iter].size();
                
                k ++;
            }
        
        }
        
        matrix.ConSparseMatrix(I, J, V, network_map.size(), network_map.size(), num_non_zero);
        
        free(I);
        free(J);
        free(V);
        
    }
    
    
    /* Compute: z = a*x + b*y + c */
	void vector_linear_operator(double * x, double *y, double* z, long n, double a, double b, double c = 0){
		
		for (long i = 0; i < n; i++) {
			z[i] = a * x[i] + b * y[i] + c;
		}
		
	}
	

    /* Power iteration method */
    void power_iteration(){
        
        //this algorithm computes the final scores, given the set of pos_train_set or prior.
        double* x;
        
        x = (double*)malloc(sizeof(double)*network_map.size());
        
        memcpy(post, prior, sizeof(double) * network_map.size());
        
        int j = 1;
        
        if (log(N) > max_iter) {
            max_iter = (int)log(N);
        }
        
        while (j <= max_iter) {
                        
            //x = A * post;
            matrix.matrix_product(post, x);
                        
            //post = (1 - alpha) * x + alpha * prior;
            vector_linear_operator(x, prior, post, N, 1 - alpha, alpha, 0);
				
            j++;
        }
        
        
    }
    

    /* Normalize scores by node degree */
    void normalize_post(){
            
        for (int i = 0; i < N; i++) {
            post[i] /= network_map[i].size();
        }
        
    }


    /* Write final posteriors to the output file */
    void write_posterior(){
        
        ofstream out(post_file, ofstream::out);
        
        for (int i = 0; i < N; i++) {
            out << i << " " << setprecision(10) << post[i] << endl;
        }
        
        out.close();
        
    }
    
    
    void parse_par(int argc, char **argv){
        
        //default setting
        network_file ="";
		train_set_file = "";
        post_file = "";
        prior_file = "";

        //restart probability
        alpha = 0;
        
        max_iter = 10;
        
        int i = 1;
        
        while (i < argc) {
            if (strcmp(argv[i],"-graphfile") == 0){
                network_file = argv[i + 1];
            }
            else if (strcmp(argv[i],"-trainfile") == 0){
                train_set_file = argv[i + 1];
            }
            else if (strcmp(argv[i],"-priorfile") == 0){
                prior_file = argv[i + 1];
            }
            else if (strcmp(argv[i],"-postfile") == 0){
                post_file = argv[i + 1];
            }
			else if (strcmp(argv[i],"-alpha") == 0){
                alpha = (double)atof(argv[i + 1]);
            }
            else if (strcmp(argv[i],"-mIter") == 0){
                max_iter = atoi(argv[i + 1]);
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
	
	Data data;
    
    data.parse_par(argc, argv);
	
	data.read_network();
    
    data.read_prior();
    
    data.construct_transition_matrix();
        
    data.power_iteration();
    
    data.normalize_post();
    
    data.write_posterior();
	
    return 0;
	
}
