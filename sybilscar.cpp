
/*
 This code was written by Binghui Wang in December, 2016. 
 
 It implements SybilSCAR algorithm proposed in the following paper:
 Binghui Wang, Le Zhang, and Neil Zhenqiang Gong. "SybilSCAR: Sybil Detection in Online Social Networks via Local Rule based Propagation", in INFOCOM, 2017. 

 For any question, please contact Binghui Wang (binghuiw@iastate.com), Le Zhang (lezhang@iastate.edu), or Neil Zhenqiang Gong (neilgong@iastate.edu).
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

// random generator function
ptrdiff_t myrandom (ptrdiff_t i) {return rand() % i;}

// pointer object to it
ptrdiff_t (*p_myrandom)(ptrdiff_t) = myrandom;


class Data{
public:
    typedef unsigned long int vertex;
    
    long N;
    
    //adjacency list
    //weighted graph
    unordered_map<vertex,list<pair<vertex, double> > > network_map;
    
    //1: weighted, 0: unweighted
    int weighted_graph;
    
    //priors: initial probability of being benign for each node
    double *prior;
    
    //network file
    char *network_file;
    
    //posterior file
    char *post_file;

    //prior file 
    char *prior_file;
    
    //during computation, post is treated as the normalized multiplication of incoming messages of each node.
    double *post;
    double *post_pre;
    
    //to support different input format
    char *train_set_file;
    
    double theta_pos; // positive nodes
    double theta_neg; // negative nodes
    double theta_unl; // unlabeled nodes
    
    //this weight is used when the input is an unweighted graph
    double weight;
    
    //number of maximum iterations
    double max_iter;
    
    //shuffle  
    int* ordering_array;
    
    //number of threads
    int num_threads;
    
    
    class lbp_arg{
    public:
        Data * data_pointer;
        int current_thread;
        lbp_arg(){}
    };
    
    Data(){
        
        
    }
    

    /* Store edges and weights in the network_map */
    void add_edge(vertex node1, vertex node2, double w){
        
        // no self loops
        if(node1==node2){
            return;
        }
        
        //add node2 to the adjacency list of node1
        pair <vertex, double> nei (node2, w);
        network_map[node1].push_back(nei);
    }
    

    /* Read the social graph */
    //the format for the social graph is
    //each line corresponds to an edge, e.g, 3 2 0.8
    //each edge in the graph appears twice, e.g.,
    //3 2 0.8
    //2 3 0.9
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
                w=(double)atof(strtok(NULL," \n\t\r")) - 0.5;
            }
            else{
                w = weight - 0.5;
            }
            add_edge(node1, node2, w);
        }
        
        //number of nodes in the graph
        N = network_map.size();
        
        //allocate space for final scores
        post = (double*) malloc(sizeof(double)*(N));
        post_pre = (double*) malloc(sizeof(double)*(N));
        
        //allocate space for prior scores
        prior = (double*) malloc(sizeof(double)*(N));
        
        in.close();
        
    }
    

    /* Read priorfile or/and training file */
    void read_prior(){
        
        //initialize priors as theta_unl
        vertex node;
        for (node = 0; node < N; node++) {
            prior[node] = theta_unl - 0.5;
        }
        
        if (prior_file != "") {
            
            ifstream in(prior_file,ifstream::in);
            assert(in);
            
            string line;
            
            double score;
            while(getline(in,line)!=NULL){
                node=(vertex)atol(strtok((char *)line.c_str()," \n\t\r"));
                score=(double)atof(strtok(NULL," \n\t\r"));
                prior[node] = score - 0.5;
            }
            
            in.close();
        }
        
        //reassign node priors for labeled benign (positive) nodes or/and Sybils (negative nodes) 
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
                prior[sub] = theta_pos - 0.5;
            }
            
            //reading labeled Sybils
            getline(in,line);
            istringstream neg_train_str(line);
            while (neg_train_str){
                neg_train_str >> sub;
                prior[sub] = theta_neg - 0.5;
            }
            
            in.close();
            
        }
        
    }
    

    /* Write final posteriors to the output file */
    void write_posterior(){
        
        ofstream out(post_file, ofstream::out);
        
        for (vertex i = 0; i < N; i++) {
            out << i << " " << setprecision(10) << post[i] + 0.5 << endl;
        }
        
        out.close();
        
    }
    

    static void * lbp_thread(void *arg_pointer){
        
        Data * pointer = ((lbp_arg *)arg_pointer)->data_pointer;
        int current_thread = ((lbp_arg *)arg_pointer)->current_thread;
        
        int num_nodes = ceil(((float)pointer->N) / pointer->num_threads);
        int start = current_thread * num_nodes;
        int end = current_thread * num_nodes + num_nodes;
        
        
        if (end > pointer->N) {
            end = pointer->N;
        }
        
        vertex node;
        list<pair<vertex, double> >::iterator nei_iter;
        
        for (vertex index = start; index < end; index++) {
            node = pointer->ordering_array[index];
            
            //update the the post for node
            for (nei_iter = pointer->network_map[node].begin(); nei_iter != pointer->network_map[node].end(); nei_iter++) {

                pointer->post[node] += 2 * pointer->post_pre[(*nei_iter).first] * (*nei_iter).second;
                
            }
        
            pointer->post[node] += pointer->prior[node];

            if (pointer->post[node] > 0.5){
                pointer->post[node] = 0.5;
            }

            if (pointer->post[node] < -0.5){
                pointer->post[node] = -0.5;
            }
            
        }
        
    }
     

    void lbp(){
    
        ordering_array = (int*) malloc(sizeof(int) * N);
        
        //initialize posts
        memcpy(post, prior, sizeof(double) * N);
        
        //random ordering
        vector<vertex> ordering;
        for (vertex i = 0; i < N; i++) {
            ordering.push_back(i);
        }
        
        int iter = 1;
        vector<vertex>::iterator iter_order;
        
        lbp_arg * arg_pointer;
        int current_thread;
        pthread_t thread;
        vector<pthread_t> threads;
        
        vertex i = 0;
        
        
        while (iter <= max_iter) {
            
            threads.clear();
            
            memcpy(post_pre, post, sizeof(double) * N);
            
            random_shuffle(ordering.begin(), ordering.end(), p_myrandom);
            for (iter_order = ordering.begin(), i = 0; iter_order != ordering.end(); iter_order++, i++){
                ordering_array[i] = *iter_order;
            }
            
            for (current_thread = 0; current_thread < num_threads; current_thread++) {
                arg_pointer = new lbp_arg();
                arg_pointer->data_pointer = this;
                arg_pointer->current_thread = current_thread;
                
                pthread_create(&thread, NULL, lbp_thread, (void*)arg_pointer);
                threads.push_back(thread);
            }
            
            for (i = 0; i < threads.size(); i++) {
                pthread_join(threads[i], NULL);
            }
            
            iter += 1;
            
        }
    }
    
    
    void parse_par(int argc, char **argv){
        
        //default setting
        network_file = "";
        train_set_file = "";
        post_file = "";
        prior_file = "";
        

        theta_pos = 0.6;
        theta_neg = 0.4;
        theta_unl = 0.5;
        
        weighted_graph = 0;
		weight = 0.6;
		
        num_threads = 1;
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
            else if (strcmp(argv[i], "-postfile") == 0){
                post_file = argv[i + 1];
            }
            else if (strcmp(argv[i],"-mIter") == 0){
                max_iter = atoi(argv[i + 1]);
            }
            else if (strcmp(argv[i],"-tp") == 0){
                theta_pos = atof(argv[i + 1]);
            }
            else if (strcmp(argv[i],"-tn") == 0){
                theta_neg = atof(argv[i + 1]);
            }
            else if (strcmp(argv[i],"-tu") == 0){
                theta_unl = atof(argv[i + 1]);
            }
            else if (strcmp(argv[i],"-wei") == 0){
                weight = atof(argv[i + 1]);
            }
            else if (strcmp(argv[i],"-wg") == 0){
                weighted_graph = atoi(argv[i + 1]);
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
    
    data.lbp();
    
    data.write_posterior();

    return 0;
}
