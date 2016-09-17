
/*
 This code was written by Le Zhang in March, 2016, and modified by Binghui Wang in August, 2016.
 
 It implements SybilRank algorithm using multi-thread.

Compile: g++ sybilrank_multi_thread.cpp -pthread -O3 -o sybilrank_multi_thread
Command: ./sybilrank_multi_thread -graphfile GRAPHFILE [-priorfile PRIORFILE] -trainfile TRAINFILE
        -mIter MAXITER -alpha ALPHA -nt NT -wei WEI -postfile POSTFILE 
(Please set -wei 1.0 for unweighted graph)

For any question, please contact Le Zhang (lezhang@iastate.edu) or Binghui Wang (binghuiw@iastate.edu)
  
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
#include <pthread.h>
#include <cmath>
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

long randint(long min, long max){

    return long(rand() / (RAND_MAX + 0.0) * (max - min) + min);

}

bool my_comparator ( const pair<double, int>& l, const pair<double, int>& r)
{
    return l.first > r.first;
}

class Data{
public:
    typedef unsigned long int vertex;

    long N;

    //adjacency list
    //weighted graph
    unordered_map<vertex,set<pair<vertex, double> > > network_map;

    //1 : weighted, 0: unweighted
    int weighted_graph;

    //priors
    //probability of being benign for each node
    double *prior;

    //network file
    char *network_file;

    char *post_file;
    char *prior_file;
    
    //char *post_file_iter;

    //during computation, post is treated as the normalized multiplication of incoming messages of each node.
    //
    double *post;
    double *post_pre;

    //to support different input format
    char *train_set_file;
    
    double alpha;

    //this weight is used when the input is an unweighted graph
    double weight;

    double max_iter;

    int* ordering_array;

    int num_threads;

	set<vertex> pos_train_set;
    

    class RW_arg{
    public:
        Data * data_pointer;
        int current_thread;
        RW_arg(){}
    };

    Data(){


    }

    void add_edge(vertex node1, vertex node2, double w){
        // add edge (node1, node2)
        
        // no self loops
        if(node1==node2){
            return;
        }

		//add node2 (node1) to the adjacency list of node1 (node2)
		network_map[node1].insert(make_pair(node2, w));
        network_map[node2].insert(make_pair(node1, w));
        
    }
    
    /* Read the social graph */
    //the format for the social graph is
    //each line corresponds to an edge, e.g, 3 2.
    //each edge in the graph appears twice, e.g.,
    //3 2
    //2 3 
    void read_network(){
        
        ifstream in(network_file,ifstream::in);
        assert(in);
        
        string line;
        vertex node1,node2;
        double w;
        
        //read edges
        while(getline(in,line)!=NULL){

            node1=(vertex)atol(strtok((char *)line.c_str()," \n\t\r"));
            node2=(vertex)atol(strtok(NULL," \n\t\r"));
            if (weighted_graph == 1) {
                w=(double)atof(strtok(NULL," \n\t\r"));
            }
            else{
                w = weight;
            }

            add_edge(node1, node2, w);
            
        }
        
        in.close();
        

        N = network_map.size();
        
        //allocate space for final scores
        post = (double*) malloc(sizeof(double)*(N));
        post_pre = (double*) malloc(sizeof(double)*(N));
        
        //allocate space for final scores
        prior = (double*) malloc(sizeof(double)*(N));
        
    }


    /* Read priorfile or/and training file */
	void read_prior(){

		//initialize priors as theta_unl
		vertex node;
		for (node = 0; node < N; node++) {
			prior[node] = 0;
		}

		if (prior_file != "") {

			ifstream in(prior_file, ifstream::in);
			assert(in);

			string line;

			double score;
			while (getline(in, line) != NULL){
				node = (vertex)atol(strtok((char *)line.c_str(), " \n\t\r"));
				score = (double)atof(strtok(NULL, " \n\t\r"));
				prior[node] = score;
			}
			in.close();
		}

        //The training dataset file includes two lines.
        //the first line includes a set of labeled benign nodes. 
        //the second line includes a set of labeled Sybil nodes. 
        //SybilRank only uses labeled benign nodes. 
        //For convenient comparison with other algorithms, we also include Sybil nodes in the training file. 
		if (train_set_file != ""){

			ifstream in(train_set_file, ifstream::in);
			assert(in);

			string line;

			//reading labeled benign nodes.
			getline(in, line);
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


    static void * RW_thread(void *arg_pointer){

        Data * pointer = ((RW_arg *)arg_pointer)->data_pointer;
        int current_thread = ((RW_arg *)arg_pointer)->current_thread;

        int num_nodes = ceil(((float)pointer->N) / pointer->num_threads);
        int start = current_thread * num_nodes;
        int end = current_thread * num_nodes + num_nodes;


        if (end > pointer->N) {
            end = pointer->N;
        }
    
        vertex node;
        double message;
        
        set<pair<vertex, double> >::iterator iter;
        double nei_weight;
        double p;
        
        unordered_map<vertex,set<pair<vertex, double> > >::iterator iter_node, find_nei_wei;

        for (vertex index = start; index < end; index++) {
            node = pointer->ordering_array[index];
          
            message = 0;
            
			if (pointer->network_map.find(node) != pointer->network_map.end()){
				for (iter = pointer->network_map[node].begin(); iter != pointer->network_map[node].end(); iter++){
                    
                    double sum_wei = 0;

					find_nei_wei = pointer->network_map.find((*iter).first);
					if (find_nei_wei != pointer->network_map.end()){
                        // unweighted graph: sum weight is the number of neighbors.
                        sum_wei = find_nei_wei->second.size();
                    }else{
                        cout<<"network outputs error at node "<<(*iter).first<<endl;
                    }
                
                    if(sum_wei == 0){
                        message += 0;
                    }else{
                        message += pointer -> post_pre[(*iter).first] * (*iter).second / sum_wei;
                    }
                }
            }
            
            pointer->post[node] = (1 - pointer->alpha) * message + pointer->alpha * pointer->prior[node];

        }

    }

 
    /* Random Walk */
    void RW(){
        
        ordering_array = (int*) malloc(sizeof(int) * (N) );

        //initialize posts
        memcpy(post, prior, sizeof(double) * (N));

        //random ordering
        vector<vertex> ordering;
        for (vertex i = 0; i < N; i++) {
            ordering.push_back(i);
        }

        int iter = 1;
        vector<vertex>::iterator iter_order;

        RW_arg * arg_pointer;
        int current_thread;
        pthread_t thread;
        vector<pthread_t> threads;

        vertex i = 0;

        max_iter = (int)log(N);
        
        while (iter <= max_iter) {

            threads.clear();

            memcpy(post_pre, post, sizeof(double) * (N));

            random_shuffle(ordering.begin(), ordering.end(), p_myrandom);
            for (iter_order = ordering.begin(), i = 0; iter_order != ordering.end(); iter_order++, i++){
                ordering_array[i] = *iter_order;
            }

            for (current_thread = 0; current_thread < num_threads; current_thread++) {
                arg_pointer = new RW_arg();
                arg_pointer->data_pointer = this;
                arg_pointer->current_thread = current_thread;

                pthread_create(&thread, NULL, RW_thread, (void*)arg_pointer);
                threads.push_back(thread);
            }

            for (i = 0; i < threads.size(); i++) {
                pthread_join(threads[i], NULL);
            }
            
            iter += 1;

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
        network_file = "";
        train_set_file = "";
        post_file = "";
        prior_file = "";

        alpha = 0.15;
        max_iter = 10;
        num_threads = 1;

		//by default, unweighted graph
		weighted_graph = 0;
		weight = 1;


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
            else if (strcmp(argv[i], "-postfile") == 0){
                post_file = argv[i + 1];
            }
            else if (strcmp(argv[i],"-mIter") == 0){
                max_iter = atoi(argv[i + 1]);
            }
			else if (strcmp(argv[i],"-alpha") ==0){
          	alpha = atof(argv[i + 1]);
			}
            else if (strcmp(argv[i],"-nt") == 0){
                num_threads = atoi(argv[i + 1]);
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
 
    data.RW();

	data.normalize_post();
   
    data.write_posterior();
  
    return 0;
}
