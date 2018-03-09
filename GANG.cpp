
/*
 This code was written by Binghui in December 2016.

 For any question, please contact Binghui Wang (binghuiw@iastate.edu).
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

 
class Data{
public:
    typedef unsigned long int vertex;

    long N;

    //adjacency list
    //weighted graph
    unordered_map<vertex,list<pair<vertex, double> > > network_map_out, network_map_in, network_map_doub;

    //1 : weighted, 0: unweighted
    int weighted_graph;

    //priors: probability of being benign for each node
    double *prior;

    //network file
    char *network_file;

    // node posterior file
    char *post_file;
    // node prior file
    char *prior_file;

    double *post;
    double *post_pre;

    // labeled training set file
    char *train_set_file;

    //parameters of PMR model when input is a set of labeled nodes, instead of priors
    double theta_pos;
    double theta_neg;
    double theta_unl;

    //this weight is used when the input is an unweighted graph
    double weight;

    // maximal iteration
    double max_iter;

    // shuffle
    int* ordering_array;

    int num_threads;


    class lbp_arg{
    public:
        Data * data_pointer;
        int current_thread;
        lbp_arg(){}
    };

    Data(){


    }
    
    // Outgoing or incoming graph
    void add_edge1(vertex node1, vertex node2, double w){
        // add edge (node1, node2)

        // no self loops
        if(node1==node2){
            return;
        }

        //add node2 to the adjacency list of node1
        pair <vertex, double> nei (node2, w);
        
        network_map_out[node1].push_back(nei);
        network_map_in[node2].push_back(make_pair(node1, w));
        
        
    }

    // Bidirectional graph
    void add_edge2(vertex node1, vertex node2, double w){
        
        // no self loops
        if(node1==node2){
            return;
        }
        
        //add node2 to the adjacency list of node1
        pair <vertex, double> nei (node2, w);

        network_map_doub[node1].push_back(nei);
        network_map_doub[node2].push_back(make_pair(node1, w));
        
    }
    
    // Read outgoing graph and bidirectional graph 
    void read_network(){
 
        string line;
        vertex node1,node2;
        vertex max_node = 0;
        double w;
        

        //outgoing graph
        char* network_type_out = "outgoing";
    	
        /* Certain platforms do not support this way of reading the filename*/
        //  const char *addr_out;
        //  stringstream network_addr;
        //  network_addr << network_file << network_type_out << ".txt";
        //  addr_out = network_addr.str().c_str();
     
    	char* addr_out = (char*)malloc(strlen(network_file) * 4);
    	strcpy(addr_out, network_file);
    	strcat(addr_out, network_type_out);
    	strcat(addr_out, ".txt");

        ifstream in1(addr_out,ifstream::in);
        assert(in1);

        //read edges
        while(getline(in1,line) != NULL){

            node1 = (vertex)atol(strtok((char *)line.c_str()," \n\t\r"));
            node2 = (vertex)atol(strtok(NULL," \n\t\r"));
            if (weighted_graph == 1) {
                w = (double)atof(strtok(NULL," \n\t\r")) - 0.5;
            }
            else{
                w = weight - 0.5;
            }
            
            add_edge1(node1, node2, w);
            if(node1 > max_node){
                max_node = node1;
            }
            if(node2 > max_node){
                max_node = node2;
            }
            
        }

        in1.close();
        
        // bidirectional grph
        char* network_type_bi = "bidirection";

        /* Certain platforms do not support this way of reading the filename*/
        //  const char *addr_bi;
        //  stringstream network_addr;
        //  network_addr << network_file << network_type << ".txt";
        //  addr_bi = network_addr.str().c_str();
        
    	char* addr_bi = (char*)malloc(strlen(network_file) * 4);
    	strcpy(addr_bi, network_file);
    	strcat(addr_bi, network_type_bi);
    	strcat(addr_bi, ".txt");

    	ifstream in2(addr_bi,ifstream::in);
        assert(in2);
        
        //read edges
        while(getline(in2,line) != NULL){

            node1 = (vertex)atol(strtok((char *)line.c_str()," \n\t\r"));
            node2 = (vertex)atol(strtok(NULL," \n\t\r"));
            
            if (weighted_graph == 1) {
                w = (double)atof(strtok(NULL," \n\t\r")) - 0.5;
            }
            else{
                w = weight - 0.5;
            }
            
            add_edge2(node1, node2, w);
            if(node1 > max_node){
                max_node = node1;
            }
            if(node2 > max_node){
                max_node = node2;
            }
            
        }
        
        in2.close();
        
        //number of nodes in the graph
        N = max_node + 1;

        //allocate space for final scores
        post = (double*) malloc(sizeof(double)*(N));
        post_pre = (double*) malloc(sizeof(double)*(N));

        //allocate space for final scores
        prior = (double*) malloc(sizeof(double)*(N));

    }

     // Note that initialization of prior is q, but it changes to q - 0.5 (in the residual form) to perform the calculation. 
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

    
        if (train_set_file != ""){

            ifstream in(train_set_file,ifstream::in);
            assert(in);

            string line;

            //read labeled benign nodes.
            getline(in,line);
            istringstream pos_train_str(line);
            vertex sub;
            while (pos_train_str){
                pos_train_str >> sub;
                prior[sub] = theta_pos - 0.5;
            }

            //read labeled Sybils.
            getline(in,line);
            istringstream neg_train_str(line);
            while (neg_train_str){
                neg_train_str >> sub;
                prior[sub] = theta_neg - 0.5;
            }

            in.close();

        }
       
    }

    // Mainloop of GNAG
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

        unordered_map<vertex,list<pair<vertex, double> > >::iterator iter_node;
        list<pair<vertex, double> >::iterator nei_iter;
       
        for (vertex index = start; index < end; index++) {
            
            node = pointer->ordering_array[index];

            pointer->post[node] += pointer->prior[node];


            iter_node = pointer->network_map_out.find(node);
            if(iter_node != pointer->network_map_out.end()){
                for(nei_iter = iter_node->second.begin(); nei_iter != iter_node->second.end(); nei_iter++){

                    if (pointer->post_pre[(*nei_iter).first] <= 0){
                        pointer->post[node] += 2 * pointer->post_pre[(*nei_iter).first] * (*nei_iter).second;
                    }

                }
            }
            
            iter_node = pointer->network_map_in.find(node);
            if(iter_node != pointer->network_map_in.end()){
                for(nei_iter = iter_node->second.begin(); nei_iter != iter_node->second.end(); nei_iter++){

                    if (pointer->post_pre[(*nei_iter).first] >= 0){
                        pointer->post[node] += 2 * pointer->post_pre[(*nei_iter).first] * (*nei_iter).second;
                    }
                }
            }
            
            iter_node = pointer->network_map_doub.find(node);
            if(iter_node != pointer->network_map_doub.end()){
                for(nei_iter = iter_node->second.begin(); nei_iter != iter_node->second.end(); nei_iter++){

                    pointer->post[node] += 2 * pointer->post_pre[(*nei_iter).first] * (*nei_iter).second;
    
                }
            }

            if (pointer->post[node] > 0.5){
                pointer->post[node] = 0.5;
            }
                    
            if (pointer->post[node] < -0.5){
                pointer->post[node] = -0.5;
            }
        }
    }

    // Multithread to speed up the calculation
    void lbp(){
        
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

        lbp_arg * arg_pointer;
        int current_thread;
        pthread_t thread;
        vector<pthread_t> threads;

        vertex i = 0;


        while (iter <= max_iter) {
           
            threads.clear();

            memcpy(post_pre, post, sizeof(double) * (N));

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

    /* Write final posterior probabilities of nodes to the output file */
    // The final posterior probability is changed from p (in the residual form) to p + 0.5.
    void write_posterior(){

        ofstream out_post(post_file, ofstream::out);

        for (vertex i = 0; i < N; i++) {
            out_post << i << " " << setprecision(10) << post[i] + 0.5 << endl;
        }

        out_post.close();
    }

    void parse_par(int argc, char **argv){

        //default setting
        network_file = "";
        train_set_file = "";
        post_file = "";
        prior_file = "";
        

        max_iter = 10;

        theta_pos = 0.9;
        theta_neg = 0.1;
        theta_unl = 0.5;

        //by default, weighted graph
        weighted_graph = 1;

        //depends on the average degree of the directed graph. 
        weight = 0.51; 

        num_threads = 1;

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
    clock_t start, end;
    start = clock();

    srand ( time(NULL) );
    
    Data data;
    
    data.parse_par(argc, argv);
    
    data.read_network();
    
    data.read_prior();
    
    data.lbp();

    data.write_posterior();

    end = clock();
    
    cout<<endl<<"Total time taken: "<<(double)(end-start)/CLOCKS_PER_SEC*1000<<" ms"<<endl;

    return 0;
}
