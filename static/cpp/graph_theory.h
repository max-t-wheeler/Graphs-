#include "numerical_integration.h"

template <class T>
class vertex {

    private:

    public:

        unsigned int sid;               //vertex species id
        unsigned int id;                //vertex id
        unsigned int index;             //vertex index
        double r;                       //vertex radius
        arma::urowvec degree;           //vertex degree
        arma::vec x;                    //vertex position
        std::string label;              //vertex label

        //vertex specification for directed graphs
        bool is_head;                   //equal to 1 if vertex is a head of an edge
        bool is_tail;                   //equal to 1 if vertex is a tail of an edge

        vertex(

            //exogenous
            unsigned int species,
            unsigned int ID, 
            double radius

            //endogenous

        )
            : sid(species), id(ID), r(radius) {}

        //vertex initialization 

        //assign label to vertex
        void set_label(const std::string entry) {
            
            label = entry;
        
        }
        
        //set coordinates for 2-dimensional graph
        void set_coords(T x1, T x2) {

            x << x1 << x2 << arma::endr;
        
        }

        //set coordinates for 3-dimensional graph
        void set_coords(T x1, T x2, T x3) {

            x << x1 << x2 << x3 << arma::endr;
        
        }

        //vertex qualifiers

        //determine whether or not a vertex is isolated (West Chapter 1.2.8)
        bool is_isolated(const std::string query = "all", const unsigned int layer = 0) {

            if (degree.is_empty()) {

                if(logging) {
                    output("Degree sequence in vertex::is_isolated has not been intialized");
                }

                return 0;
            
            }
            if (query == "all") {

                if (arma::all(degree == 0)) {
                    return 1;
                }
                else {
                    return 0;
                }

            }
            else if (query == "layer"){

                if (degree(layer) == 0) {
                    return 1;
                }
                else {
                    return 0;
                }
                
            }
            else {

                if(logging) {
                    output("Query in vertex::is_isolated must be all or layer");
                }
                
                return 0;
            }

        }

        //determine whether or not a vertex is even (West Chapter 1.2.24)
        bool is_even(const std::string query = "all", const unsigned int layer = 0) {

            if (degree.is_empty()) {

                if(logging) {
                    output("Degree sequence in vertex::is_even has not been intialized");
                }
                
                return 0;
            
            }

            if (query == "all") {

                for (int i = 0; i < degree.size(); ++i) {
                    if (degree(i)%2 == 1) {
                        return 0;
                    }
                }

                return 1;

            }
            else if (query == "layer"){

                if (degree(layer)%2 == 0) {
                    return 1;
                }
                else {
                    return 0;
                }
                    
            }
            else {

                if(logging) {
                    output("Query in vertex::is_even must be all or layer");
                }
                
                return 0;
            }

        }

        //determine whether or not a vertex is odd (West Chapter 1.2.24)
        bool is_odd(const std::string query = "all", const unsigned int layer = 0) {

            if (degree.is_empty()) {

                if(logging) {
                    output("Degree sequence in vertex::is_odd has not been intialized");
                }
                
                return 0;
            
            }

            if (query == "all") {

                for (int i = 0; i < degree.size(); ++i) {

                    if (degree(i)%2 == 0) {
                        return 0;
                    }

                }

                return 1;

            }
            else if (query == "layer"){

                if (degree(layer)%2 == 1) {
                    return 1;
                }
                else {
                    return 0;
                }
                
            }
            else {

                if(logging) {
                    output("Query in vertex::is_odd must be all or layer");
                }
                
                return 0;
            }
            
        }

};

class graph {

    public:
        bool loops_permitted = 0;           //equal to 1 when loops are allowed
        bool is_directed = 0;               //equal to 1 when graph is a directed graph
        int b;                              //maximum allowable distance between vertices
        unsigned int id;                    //graph ID
        unsigned int num_layers;            //number of layers
        double r;                           //graph radius
        std::string type;                   //graph type
        arma::vec x;                        //position of graph center
        arma::uvec v_sizes;                 //list of vertex set sizes
        std::vector< vertex<double> * > v;  //vertex set
        arma::umat V;                       //vertex array
        arma::umat E;                       //edge set
        arma::cube adjacency_tensor;        //adjacency tensor
        arma::imat incidence_matrix;        //incidence tensor
        arma::umat deg_mat;                 //degree matrix
        arma::vec edge_weights;             //edge weights

        graph(

            //graph attributes
            unsigned int ID, 
            arma::uvec vertex_set_sizes,
            double radius,
            arma::vec center,

            //vertex attributes
            int bound,

            //edge attributes
            unsigned int number_of_layers
            )

                : id(ID), v_sizes(vertex_set_sizes), r(radius), x(center), b(bound), num_layers(number_of_layers) {}

        //graph attributes

        void set_type() {

            if (num_layers > 1) {
                type = "multi";
            }
            else {
                type = "simple";
            }

        }

        void contains_loops() {

            loops_permitted = 1;

        }

        //return the girth of graph
        unsigned int girth() {

            unsigned int girth = 0;

            return girth;
        
        }

        //return chromatic number of graph
        unsigned int chromatic_number() {

            unsigned int chi = 0;
            
            return chi;

        }

        //vertex set attributes and algorithms

        //return order 
        unsigned int order() {

            unsigned int sum = 0;
            
            for (int i = 0; i < v_sizes.size(); ++i) {
                sum += v_sizes(i);
            }

            return sum;

        }

        //create vertex set
        void create_vertex_set(const arma::vec rad_vec, const std::string type = "") {

            unsigned int count = 0;
            unsigned int n = order();

            V = arma::zeros<arma::umat>(n, 2);

            if (n != rad_vec.size()) {

                if(logging) {
                    output("Warning: Vertex attributes ill-defined in graph::create_vertex_set");
                }

                return;
            }

            if (!v.empty()) {

                for (int i = 0; i < v.size(); i++){
                    delete v[i];
                }

                v.clear();

            }

            if (type == "random") {

                for (int i = 0; i < v_sizes.size(); ++i) {  

                    for (int j = 0; j < v_sizes(i); ++j) {

                        V(count + j, 0) = i;
                        V(count + j, 1) = j;

                        vertex<double> * node = new vertex<double>(i, j, rad_vec(count+j));
                        node->index = count + j;
                        node->set_coords(Rand(-b, b),Rand(-b, b),Rand(-b, b));

                        v.push_back(node);

                    }

                    count += v_sizes(i);

                }

            }
            else {

                for (int i = 0; i < v_sizes.size(); ++i) {  

                    for (int j = 0; j < v_sizes(i); ++j) {

                        V(count + j, 0) = i;
                        V(count + j, 1) = j;

                        vertex<double> * node = new vertex<double>(i, j, rad_vec(count+j));
                        node->index = count + j;
                        node->set_coords(0,0,0);

                        v.push_back(node);

                    }

                    count += v_sizes(i);
                
                }

            }

        }

        //create vertex array using a file
        void create_vertex_set_w_file(const arma::vec rad_vec, const std::string file) {

            unsigned int count = 0;

            V.load(file);

            if (V.n_rows != rad_vec.size()) {

                if(logging) {
                    output("Warning: Vertex attributes ill-defined in graph::create_vertex_set_w_file");
                }

                return;

            }

            insertion_sort(V, "forward", 1);

            arma::uvec u_V = arma::unique(V.col(0));

            if (!v_sizes.is_empty()) {
                v_sizes.reset();
            }

            arma::uvec temp = arma::zeros<arma::uvec>(u_V.size());
            
            v_sizes = temp;

            for (int i = 0; i < V.n_rows; ++i) {

                if (count < u_V.size() && V(i, 0) == u_V(count)) {
                    v_sizes(count)++;
                }
                else {
                    count++;
                    v_sizes(count)++;
                }

            }

            for (int i = 0; i < v.size(); i++){
                delete v[i];
            }

            v.clear();

            for (int i = 0; i < V.n_rows; ++i) {  

                vertex<double> * node = new vertex<double>(V(i,0), V(i,1), rad_vec(i));
                node->index = i;
                node->set_coords(0,0,0);

                v.push_back(node);

            }

        }

        //create vertex array using a matrix
        void create_vertex_set_w_matrix(const arma::vec rad_vec, const arma::umat vertices) {

            unsigned int count = 0;

            V = vertices;

            if (V.n_rows != rad_vec.size()) {

                if(logging) {
                    output("Warning: Vertex attributes ill-defined in graph::create_vertex_set_w_matrix");
                }

                return;
            }

            insertion_sort(V, "forward", 1);

            arma::uvec u_V = arma::unique(V.col(0)); 

            arma::uvec temp = arma::zeros<arma::uvec>(u_V.size());
            
            v_sizes = temp;

            for (int i = 0; i < V.n_rows; ++i) {

                if (count < u_V.size() && V(i, 0) == u_V(count)) {
                    v_sizes(count)++;
                }
                else {

                    count++;
                    v_sizes(count)++;

                }

            }

            for (int i = 0; i < v.size(); i++){
                delete v[i];
            }
            
            v.clear();

            for (int i = 0; i < V.n_rows; ++i) {  

                vertex<double> * node = new vertex<double>(V(i,0), V(i,1), rad_vec(i));
                node->index = i;
                node->set_coords(0,0,0);

                v.push_back(node);

            }

        }

        //assign coords to vertices in graph
        void assign_coords(const std::string layout = "polar", const double r = 1, const double angle = 0) {

            if (v.empty()) {

                if(logging) {
                    output("Warning: Vertex set used in graph::assign_coords has not been initialized");
                }

                return;
            }

            unsigned int n = order();

            if (layout == "polar") {
                for (int i = 0; i < n; ++i) {
                    v[i]->x(0) = x(0) + r*cos(2*M_PI*i/n + angle);
                    v[i]->x(1) = x(1) + r*sin(2*M_PI*i/n + angle);
                }
            }
            else if (layout == "partite") {
                unsigned int count = 0;
                for (int i = 0; i < v_sizes.size(); ++i) {
                    if (v_sizes(i)%2 == 1) {
                        for (int j = 0; j <  v_sizes(i); ++j) {    
                            v[count + j]->x(0) = (((1.0-v_sizes(i))/2 + j)*r/v_sizes(i))*cos(2*(1-i)*M_PI/v_sizes.size() + angle) + 0.8*cos(2*(1-i)*M_PI/v_sizes.size() - (M_PI/2) + angle);
                            v[count + j]->x(1) = (((1.0-v_sizes(i))/2 + j)*r/v_sizes(i))*sin(2*(1-i)*M_PI/v_sizes.size() + angle) + 0.8*sin(2*(1-i)*M_PI/v_sizes.size() - (M_PI/2) + angle);     
                        }
                    }
                    else {
                        for (int j = 0; j <  v_sizes(i); ++j) {
                            if (j != (v_sizes(i)-1.0)/2) {
                                v[count + j]->x(0) = (((1.0-v_sizes(i))/2 + j)*r/v_sizes(i))*cos(2*(1-i)*M_PI/v_sizes.size() + angle) + 0.8*cos(2*(1-i)*M_PI/v_sizes.size() - (M_PI/2) + angle);
                                v[count + j]->x(1) = (((1.0-v_sizes(i))/2 + j)*r/v_sizes(i))*sin(2*(1-i)*M_PI/v_sizes.size() + angle) + 0.8*sin(2*(1-i)*M_PI/v_sizes.size() - (M_PI/2) + angle);     
                            }
                        }
                    }
                    count += v_sizes(i);
                }

            }
            else if (layout == "concentric") {
                unsigned int count = 0;
                for (int i = 0; i < v_sizes.size(); ++i) {
                    for (int j = 0; j < v_sizes(i); ++j) {
                        if (v_sizes(i) == 1) {
                            v[count + j]->x(0) = 0;
                            v[count + j]->x(1) = 0;
                        }
                        else {
                            v[count + j]->x(0) = x(0) + ((r+i)/(v_sizes.size()))*cos(2*M_PI*j/v_sizes(i) + angle);
                            v[count + j]->x(1) = x(1) + ((r+i)/(v_sizes.size()))*sin(2*M_PI*j/v_sizes(i) + angle);

                        }
                    }
                    count += v_sizes(i);
                }
            }
            //needs work
            else if (layout == "grid") {
                unsigned int count = 0;
                for (int i = 0; i < v_sizes.size(); ++i) {
                    if (v_sizes(i)%2 == 1) {
                        for (int j = 0; j <  v_sizes(i); ++j) {    
                            v[count + j]->x(0) = (((1.0-v_sizes(i))/2 + j)*r/v_sizes(i));
                            v[count + j]->x(1) = (((1.0-v_sizes.size())/2 + i)*r/v_sizes(i));
                        }
                    }
                    else {
                        for (int j = 0; j <  v_sizes(i); ++j) {
                            if (j == (v_sizes(i)-1.0)/2) {

                            }
                            else {
                                v[count + j]->x(0) = (((1.0-v_sizes(i))/2 + j)*r/v_sizes(i));
                                v[count + j]->x(1) = (((1.0-v_sizes.size())/2 + i)*r/v_sizes(i));
                            }
                        }
                    }
                    count += v_sizes(i);
                }

            }

        }

        //access vertex by IDs
        vertex<double> * get_vertex(const unsigned int species, const unsigned int ID) {

            if (v.empty()) {

                if(logging) {
                    output("Warning: Vertex set used in graph::get_vertex is ill-defined");
                }

                return nullptr;
            }

            for (int i = 0; i < v.size(); ++i) {
                if (v[i]->sid == species && v[i]->id == ID) {
                    return v[i];
                }
            }
            
            return nullptr;
        }


        //edge set attributes

        //Degree-Sum Formula  (West Chapter 1.3.3)
        unsigned int size(const std::string query = "all", const unsigned int layer = 0) {

            if (deg_mat.is_empty()) {

                if(logging) {
                    output("Error: Degree matrix in graph::size has not been initialized");
                }
                return 0;

            }

            if (query == "all") {

                unsigned int sum = arma::accu(deg_mat);

                return sum/2;

            }
            else if (query == "layer"){

                if (layer > num_layers-1) {

                    if(logging) {
                        output("Layer specified in graph::size is out of bounds");
                    }

                    return 0;
                }

                unsigned int sum = arma::accu(deg_mat.col(layer));

                return sum/2;
                
            }
            else {

                if(logging) {
                    output("Query in graph::size must be all or layer");
                }
                
                return 0;
            
            }

        }

        //create edge set with file
        void create_edge_set_w_file(const std::string file) {
            E.load(file);
        }

        //create graph's adjacency tensor
        void create_adjacency_tensor(const std::string type = "complete", const std::string query = "layer", const unsigned int layer = 0) {

            unsigned int n = order();

            adjacency_tensor = arma::zeros(n, n, num_layers);

            if (query == "all") {

                if(type == "custom") {

                    if (E.n_rows < 3) {

                        if (logging) {
                            output("Edge set in create_adjacency_tensor is not multilayer");
                        }

                        return;
                    }

                    //needs work
                    for (int i = 0; i < E.n_rows; ++i) {

                        adjacency_tensor(E(i, E.n_cols-2), E(i, E.n_cols-1), E(i, 0))++;
                        adjacency_tensor(E(i, E.n_cols-1), E(i, E.n_cols-2), E(i, 0)) = adjacency_tensor(E(i, E.n_cols-2), E(i, E.n_cols-1), 0);  

                        if (E(i, E.n_cols-2) == E(i, E.n_cols-1)) {
                            adjacency_tensor(E(i, E.n_cols-2), E(i, E.n_cols-1), E(i, 0))++;
                        }

                    }

                }
                else if (type == "complete") {

                    for (int k = 0; k < num_layers; ++k) {

                        for (int i = 0; i < n - 1; ++i) {

                            for (int j = i + 1; j < n; ++j) {

                                adjacency_tensor(i, j, k) = 1;
                                adjacency_tensor(j, i, k) = adjacency_tensor(i, j, k); 

                            }

                            if(loops_permitted) {
                                adjacency_tensor(i, i, k) = 2;
                            }

                        }

                    }

                    if(loops_permitted) {

                        for (int k = 0; k < num_layers; ++k) {
                            adjacency_tensor(n - 1, n - 1, k) = 2;
                        }

                    }

                }
                else if (type == "complete k-partite") {

                    if (v_sizes.is_empty()) {

                        if (logging) {
                            output("Vertex set in create_adjacency_tensor has not been initialized");
                        }

                        return;

                    }

                    unsigned int count = 0;
                    unsigned int l = 0;

                    for (int k = 0; k < num_layers; ++k) {

                        for (int i = 0; i < n-1; ++i) {

                            for (int j = i+1; j < n; ++j) {

                                if (j < (count + v_sizes(l))) {
                                    adjacency_tensor(i, j, k) = 0;
                                }
                                else {
                                    adjacency_tensor(i, j, k) = 1;
                                }

                                adjacency_tensor(j, i, k) = adjacency_tensor(i, j, k); 

                            }

                            if (i >= count + v_sizes(l)-1) {

                                count += v_sizes(l);
                                l++;

                            }

                        }
                    }

                }
                //needs to be retooled
                else if (type == "random") { 

                    for (int k = 0; k < num_layers; ++k) {

                        for (int i = 0; i < n - 1; ++i) {

                            for (int j = i+1; j < n; ++j) {

                                if (arma::norm(v[i]->x - v[j]->x) < 1000) {
                                    
                                    adjacency_tensor(i, j, k) = 1;
                                    adjacency_tensor(j, i, k) = adjacency_tensor(i, j, k); 
                                
                                }

                            }

                        }

                    }

                }
                else if (type == "null") {

                }
                else {

                    if(logging) {
                        output("Type entered in graph::create_adjacency_tensor is invalid");
                    }

                    return;
                }

                //adjust degree matrix
                update_degrees();

                //adjust edge set
                if (type != "custom") {
                    update_edge_set();
                }   

            }
            else if (query == "layer"){

                if (layer > num_layers-1) {

                    if(logging) {
                        output("Layer specified in graph::create_adjacency_tensor is out of bounds");
                    }
                    
                    return;
                
                }
                        //construct adjacency matrix
                if(type == "custom") {
                    for (int i = 0; i < E.n_rows; ++i) {

                        adjacency_tensor(E(i, E.n_cols-2), E(i, E.n_cols-1), layer)++;
                        adjacency_tensor(E(i, E.n_cols-1), E(i, E.n_cols-2), layer) = adjacency_tensor(E(i, E.n_cols-2), E(i, E.n_cols-1), layer);  

                        if (E(i, E.n_cols-2) == E(i, E.n_cols-1)) {
                            adjacency_tensor(E(i, E.n_cols-2), E(i, E.n_cols-1), layer)++;
                        }

                    }

                }
                else if (type == "complete") {

                    for (int i = 0; i < n - 1; ++i) {

                        for (int j = i + 1; j < n; ++j) {

                            adjacency_tensor(i, j, layer) = 1;
                            adjacency_tensor(j, i, layer) = adjacency_tensor(i, j, layer); 
                       
                        }

                        if(loops_permitted) {
                            adjacency_tensor(i, i, layer) = 2;
                        }

                    }

                    if(loops_permitted) {
                        adjacency_tensor(n - 1, n - 1, layer) = 2;
                    }

                }
                else if (type == "complete k-partite") {

                    unsigned int count = 0;
                    unsigned int k = 0;

                    for (int i = 0; i < n-1; ++i) {

                        for (int j = i+1; j < n; ++j) {

                            if (j < (count + v_sizes(k))) {
                                adjacency_tensor(i, j, layer) = 0;
                            }
                            else {
                                adjacency_tensor(i, j, layer) = 1;
                            }

                            adjacency_tensor(j, i, layer) = adjacency_tensor(i, j, layer); 

                        }

                        if (i >= count + v_sizes(k)-1) {

                            count += v_sizes(k);
                            k++;

                        }

                    }

                }
                else if (type == "cycle") {

                    for (int i = 0; i < n-1; ++i) {

                        adjacency_tensor(i, i+1, layer) = 1;
                        adjacency_tensor(i+1, i, layer) = adjacency_tensor(i, i+1, layer); 

                    }

                    adjacency_tensor(0, n-1, layer) = 1;
                    adjacency_tensor(n-1, 0, layer) = adjacency_tensor(0, n-1, layer); 

                }
                else if (type == "path") {

                    for (int i = 0; i < n-1; ++i) {

                        adjacency_tensor(i, i+1, layer) = 1;
                        adjacency_tensor(i+1, i, layer) = adjacency_tensor(i, i+1, layer); 

                    }

                }
                else if (type == "wheel") {

                    for (int i = 1; i < n-1; ++i) {

                        adjacency_tensor(0, i, layer) = 1;
                        adjacency_tensor(i, 0, layer) = adjacency_tensor(0, i, layer); 

                        adjacency_tensor(i, i+1, layer) = 1;
                        adjacency_tensor(i+1, i, layer) = adjacency_tensor(i, i+1, layer); 

                    }

                    adjacency_tensor(0, n-1, layer) = 1;
                    adjacency_tensor(n-1, 0, layer) = adjacency_tensor(0, n-1, layer); 

                    adjacency_tensor(1, n-1, layer) = 1;
                    adjacency_tensor(n-1, 1, layer) = adjacency_tensor(1, n-1, layer); 

                }
                //needs to be retooled
                else if (type == "random") { 

                    for (int i = 0; i < n-1; ++i) {

                        for (int j = i+1; j < n; ++j) {

                            if (arma::norm(v[i]->x - v[j]->x) < 1000) {

                                adjacency_tensor(i, j, layer) = 1;
                                adjacency_tensor(j, i, layer) = adjacency_tensor(i, j, layer); 
                            
                            }

                        }

                    }

                }
                else if (type == "null") {

                }
                else {

                    if(logging) {
                        output("Type entered in graph::create_adjacency_tensor is invalid");
                    }

                    return;

                }

                //adjust degree matrix
                update_degrees();

                //adjust edge set
                if (type != "custom") {
                    update_edge_set();
                }

            }
            else {

                if (logging) {
                    output("Query in graph::create_adjacency_tensor must be all or layer");
                }

                return;

            }

        }

        //create graph's incidence matrix
        void create_incidence_matrix(const unsigned int layer) {

            if (adjacency_tensor.is_empty()) {
                output("Adjacency tensor is ill-defined in create_incidence_matrix");
                return;
            }

            if (size() < 1) {
                output("Degree matrix is ill-defined in create_incidence_matrix");
                return;
            }

            unsigned int k = 0;
            unsigned int increment = 1;
            unsigned int n = order();
            unsigned int m = size("layer", layer);

            incidence_matrix = arma::zeros<arma::imat>(n, m);

            if (loops_permitted) {
                increment = 0;
            }

            for (int i = 0; i < n - increment; ++i) {

                for (int j = i + increment; j < n; ++j) {
                    
                    if (adjacency_tensor(i, j, layer) > 0) {

                        incidence_matrix(i, k)++;
                        incidence_matrix(j, k)++;
                        k++;

                    }
                    
                }

            }

        }

        //update edge set using adjacency tensor
        void update_edge_set() {

            if (type.empty()) {
                
                if(logging) {
                    output("Type in graph::update_edge_set has not been initialized");
                }

                return;

            }

            if (adjacency_tensor.is_empty()) {
                
                if(logging) {
                    output("Adjacency tensor in graph::update_edge_set has not been initialized");
                }

                return;

            }

            if (size() < 1) {

                if(logging) {
                    output("No edges to update in graph::update_edge_set");
                }

                E.reset();

                return;

            }

            if (!E.is_empty()) {
                E.reset();
            }

            unsigned int increment = 1;
            unsigned int m = size();
            unsigned int n = order();

            if (loops_permitted) {
                increment = 0;
            }   


            if (num_layers > 1 && weighted()) {
                E = arma::zeros<arma::umat>(m, 4);
            }
            else if (num_layers > 1 || weighted()) {
                E = arma::zeros<arma::umat>(m, 3);
            }
            else {
                E = arma::zeros<arma::umat>(m, 2);
            }

            if (type == "simple") {

                unsigned int count = 0;

                for (int i = 0; i < n - increment; ++i) {

                    for (int j = i + increment; j < n; ++j)  {

                        if (adjacency_tensor(i, j, 0) > 0) {   

                            E(count,E.n_cols-2) = i;
                            E(count, E.n_cols-1) = j;

                            count++;

                        }

                    }

                }

                insertion_sort(E, "F");

            }
            else if (type == "multi") {

                for (int k = 0; k < num_layers; ++k) {

                    unsigned int count = 0;

                    for (int i = 0; i < n - increment; ++i) {

                        for (int j = i + increment; j < n; ++j)  {

                            if (adjacency_tensor(i, j, k) > 0) {

                                E(count,E.n_cols-2) = i;
                                E(count, E.n_cols-1) = j;

                                count++;

                            }

                        }

                    }

                }

                insertion_sort(E, "F");

            }
            else {

                if(logging) {
                    output("Type in graph::update_edge_set must be simple or multi");
                }

                return;

            }

        }

        //update degree matrix using adjacency tensor (West Chapter 1.1.18 / West Exercise 1.1.30)
        void update_degrees() {

            if (adjacency_tensor.is_empty()) {

                if(logging) {
                    output("Adjacency tensor in graph::update_degrees has not been initialized");
                }

                return;

            }

            //assign new vertex degrees
            if (loops_permitted) {

                deg_mat.reset();

                deg_mat = arma::zeros<arma::umat>(v.size(), num_layers);

                for (int i = 0; i < adjacency_tensor.n_rows; ++i) {

                    v[i]->degree.reset();

                    for (int j = 0; j < adjacency_tensor.n_cols; ++j) {

                        for (int k = 0; k  < num_layers; ++k) {

                            if (adjacency_tensor(i, j, k) > 0) {
                                deg_mat(i, k) = deg_mat(i, k) + adjacency_tensor(i, j, k);
                            }
                        
                        }
                    
                    }

                    v[i]->degree = deg_mat.row(i);
                
                }

            }
            else {

                deg_mat.reset();

                deg_mat.set_size(adjacency_tensor.n_rows, num_layers);
                
                for (int i = 0; i < adjacency_tensor.n_rows; ++i) {

                    v[i]->degree.reset();

                    for (int k = 0; k  < num_layers; ++k) {

                        arma::vec v = adjacency_tensor.slice(k).row(i)*adjacency_tensor.slice(k).col(i);
                        
                        deg_mat(i, k) = arma::accu(v);

                    }

                    v[i]->degree = deg_mat.row(i);

                }

            }

        }

        //general vertex degree attributes
        
        //return minimum degree in graph's degree sequence (West Chapter 1.3.1)
        unsigned int minimum_degree(const std::string query = "all", const unsigned int layer = 0) {

            if (deg_mat.is_empty()) {

                if(logging) {
                    output("Error: Degree matrix in graph::minimum_degree has not been initialized");
                }

                return 0;

            }

            if (query == "all") {

                arma::uvec deg_sum = arma::sum(deg_mat, 1); 

                return deg_sum.min();

            }
            else if (query == "layer"){

                if (layer > num_layers-1) {

                    if(logging) {
                        output("Layer specified in graph::minimum_degree is out of bounds");
                    }
                    
                    return 0;
                }

                return deg_mat.col(layer).min();
                
            }
            else {

                if(logging) {
                    output("Query in graph::minimum_degree must be all or layer");
                }

                return 0;
            
            }

        }

        //return maximum degree in graph's degree sequence (West Chapter 1.3.1)
        unsigned int maximum_degree(const std::string query = "all", const unsigned int layer = 0) {

            if (deg_mat.is_empty()) {

                if(logging) {
                    output("Error: Degree matrix in graph::maximum_degree has not been initialized");
                }

                return 0;
            }

            if (query == "all") {

                arma::uvec deg_sum = arma::sum(deg_mat, 1); 

                return deg_sum.max();

            }
            else if (query == "layer"){

                if (layer > num_layers-1) {

                    if(logging) {
                        output("Layer specified in graph::maximum_degree is out of bounds");
                    }
                    
                    return 0;
                }
            
                return deg_mat.col(layer).max();
                
            }
            else {

                if(logging) {
                    output("Query in graph::maximum_degree must be all or layer");
                }

                return 0;
            
            }

        }

        //return average degree of vertex in layer of graph (West Chapter 1.3.1)
        double average_degree(const std::string query = "all", const unsigned int layer = 0) {

            if (deg_mat.is_empty()) {

                if(logging) {
                    output("Error: Degree matrix in graph::average_degree has not been initialized");
                }

                return 0;
            }

            if (query == "all") {

                return (double)2*size()/order();

            }
            else if (query == "layer"){

                if (layer > num_layers-1) {

                    if(logging) {
                        output("Layer specified in graph::average_degree is out of bounds");
                    }
                    
                    return 0;
                }

                return (double)2*size(query, layer)/order();

            }
            else {

                if(logging) {
                    output("Query in graph::average_degree must be all or layer");
                }

                return 0;
            
            }

        }

        //graph adjustment algorithms 

        //remove a set of vertices from graph
        void delete_vertices(arma::umat vertices) {

            if (vertices.is_empty()) {

                if(logging) {
                    output("Error: Vertex set in graph::delete_vertices has not been initialized");
                }

                return;
            }

            insertion_sort(vertices, "forward", 1);

            unsigned int count = 0;
            unsigned int new_order = order() - vertices.n_rows;
            std::vector < vertex<double> *> w;
            arma::cube adj_ten_temp(new_order, new_order, num_layers);

            //adjust vertex set
            for (int i = 0; i <  v.size(); ++i) {

                if (count < vertices.n_rows && vertices(count, 0) == v[i]->sid && vertices(count, 1) == v[i]->id) {
                    
                    v_sizes(v[i]->sid) -= 1;
                    V.shed_row(i-count);
                    count++;

                }
                /**
                else if (count < vertices.n_rows && v[i]->sid > vertices(count, 0)){
                    w.push_back(v[i]);
                    count++;
                }
                /**/
                else {

                    v[i]->index = i - count;
                    w.push_back(v[i]);
                
                }

            }

            count = 0;

            //update adjacency tensor
            for (int k = 0; k < num_layers; ++k){

                arma::mat adj_ten_slice = adjacency_tensor.slice(k);

                for (int i = 0; i <  v.size(); ++i) {

                    if (count < vertices.n_rows && vertices(count, 0) == v[i]->sid && vertices(count, 1) == v[i]->id) {
                        
                        adj_ten_slice.shed_row(i-count);
                        adj_ten_slice.shed_col(i-count);
                        count++;

                    }

                }

                adj_ten_temp.slice(k) = adj_ten_slice;

            }

            adjacency_tensor.reset();

            adjacency_tensor = adj_ten_temp; 

            v.clear();
            v = w;

            update_degrees();

            //create_incidence_matrix();

            update_edge_set();

        }

        //remove a set of edges from graph
        void delete_edges(arma::umat edges) {

            if (edges.is_empty()) {

                if(logging) {
                    output("Error: Edge set in graph::delete_edges has not been initialized");
                }

                return;
            }

            insertion_sort(edges, "F", 1);
            
            unsigned int count = 0;
            unsigned int increment = 1;

            if (loops_permitted) {
                increment = 0;
            }

            for (int k = 0; k < num_layers; ++k){

                for (int i = 0; i < adjacency_tensor.n_rows-increment; ++i) {
                    
                    for (int j = i + increment; j < adjacency_tensor.n_cols; ++j) {
                        
                        if (count < edges.n_rows && edges(count, edges.n_cols-2) == i && edges(count, edges.n_cols-1) == j) {
                            
                            adjacency_tensor(i, j, k) = 0;
                            adjacency_tensor(j, i, k) = 0;
                            count++;

                        }
                        else if (count < edges.n_rows && i > edges(count, 0)){
                            count++;
                        }

                    }

                }

            }

            update_degrees();

            //create_incidence_matrix();

            update_edge_set();

        }

        //take complement of G (West Chapter 1.1.18a)
        void complement() {

            if (type != "simple") {

                if(logging) {
                    output("Type in graph::complement must be simple");
                }

                return;

            }

            if (adjacency_tensor.is_empty()) {

                if(logging) {
                    output("Adjacency tensor is ill-defined in graph::complement");
                }

                return;

            }

            for (int i = 0; i < adjacency_tensor.n_rows-1; ++i) {

                for (int j = i+1; j < adjacency_tensor.n_cols; ++j) {

                    if (adjacency_tensor(i, j, 0) == 1) {
                        adjacency_tensor(i, j, 0) = 0;
                    }
                    else {
                        adjacency_tensor(i, j, 0) = 1;
                    }
                }
            }

            for (int i = 0; i < adjacency_tensor.n_rows-1; ++i){

                for (int j = i+1; j < adjacency_tensor.n_cols; ++j){
                    adjacency_tensor(j, i, 0) = adjacency_tensor(i, j, 0);                            //make matrix symetric
                }

            }

            update_degrees();

            update_edge_set();

        }

        //qualifiers

        //determine whether or not a graph layer is regular (West Chapter 1.3.1)
        bool is_regular(const std::string query = "all", const unsigned int layer = 0) {

            if (deg_mat.is_empty()) {

                if(logging) {
                    output("Error: Degree matrix in graph::is_regular has not been initialized");
                }

                return 0;
            }

            if (query == "all") {

                if (minimum_degree() == maximum_degree()) {
                    return 1;
                }
                else {
                    return 0;
                }

            }
            else if (query == "layer"){

                if (layer > num_layers-1) {

                    if(logging) {
                        output("Layer specified in graph::is_regular is out of bounds");
                    }
                    
                    return 0;

                }

                if (deg_mat.col(layer).min() == deg_mat.col(layer).max()) {
                    return 1;
                }
                else {
                    return 0;
                }

            }
            else {

                if(logging) {
                    output("Query in graph::is_regular must be all or layer");
                }

                return 0;

            }
        }

        //determine whether or not a graph layer is k-regular (West Chapter 1.3.1)
        bool is_k_regular(const unsigned int k, const std::string query = "all", const unsigned int layer = 0) {


            if (deg_mat.is_empty()) {

                if(logging) {
                    output("Error: Degree matrix in graph::is_k_regular has not been initialized");
                }

                return 0;

            }

            if (query == "all") {

                if (is_regular() && k == maximum_degree()) {
                    return 1;
                }
                else {
                    return 0;
                }

            }
            else if (query == "layer"){

                if (layer > num_layers-1) {

                    if(logging) {
                        output("Layer specified in graph::is_k_regular is out of bounds");
                    }
                    
                    return 0;
                }

                if (is_regular(query, layer) && k == deg_mat.col(layer).max()) {
                    return 1;
                }
                else {
                    return 0;
                }

            }
            else {

                if(logging) {
                    output("Query in graph::is_k_regular must be all or layer");
                }

                return 0;

            }

        }

        bool is_connected(const std::string query = "all", const unsigned int layer = 0) {

            if (deg_mat.is_empty()) {

                if(logging) {
                    output("Error: Degree matrix in graph::is_connected has not been initialized");
                }

                return 0;

            }

            if (query == "all") {

                unsigned int e = size();
                unsigned int n = order();

                //West Chapter 1.3.13
                if (e < n-1) {
                    return 0;
                }


                if(logging) {
                    output("The methods employed in graph::is_connected were not capable of determining whether or not graph was connected");
                }

                return 1;

            }
            else if (query == "layer"){

                if (layer > num_layers-1) {

                    if(logging) {
                        output("Layer specified in graph::is_connected is out of bounds");
                    }
                    
                    return 0;

                }

                unsigned int e = size(query, layer);
                unsigned int n = order();
                unsigned int Delta = maximum_degree(query, layer);
                unsigned int delta = minimum_degree(query, layer);

                //true for all graphs

                //West Chapter 1.3.13
                if (e < n-1) {
                    return 0;
                }

                //true for simple graphs 

                //West Exercise 1.3.40c
                if (n > 2 && e > (n-1)*(n-2)/2) {
                    return 1;
                }

                //West Chapter 1.3.15
                if (delta >= (double)(n-1)/2) {
                    return 1;
                }

                //West Exercise 1.3.41
                if (n > 3 && Delta == ceil(n/2) && delta == floor(n/2)-1) {
                    return 1;
                }

                //West Exercise 1.3.64
                unsigned int index = 0;
                unsigned int count = 0;

                arma::uvec sorted_deg_seq = arma::sort(deg_mat.col(layer), "ascend");

                while (index <= n - 1 - Delta) {

                    if (deg_mat(index, layer) < index) {
                        count++;
                    }

                    index++;

                }

                if (count == index) {
                    return 1;
                }

                if(logging) {
                    output("The methods employed in graph::is_connected were not capable of determining whether or not graph was connected");
                }

                return 1;

            }
            else {

                if(logging) {
                    output("Query in graph::is_connected must be all or layer");
                }

                return 0;

            }
        }

        //determine whether or not a graph is trivial (West Chapter 1.2.8)
        bool is_trivial(const std::string query = "all", const unsigned int layer = 0) {

            if (deg_mat.is_empty()) {

                if(logging) {
                    output("Error: Degree matrix in graph::is_trivial has not been initialized");
                }

                return 0;

            }

            if (query == "all") {

                if (size() == 0) {
                    return 1;
                }
                else {
                    return 0;
                }

            }
            else if (query == "layer"){

                if (layer > num_layers-1) {

                    if(logging) {
                        output("Layer specified in graph::is_trivial is out of bounds");
                    }
                    
                    return 0;

                }

                if (size(query, layer) == 0) {
                    return 1;
                }
                else {
                    return 0;
                }
                
            }
            else {

                if(logging) {
                    output("Query in graph::is_trivial must be all or layer");
                }

                return 0;

            }
        }

        //determine whether or not a graph layer is even (West Chapter 1.2.24)
        bool is_even(const std::string query = "all", const unsigned int layer = 0) {

            if (deg_mat.is_empty()) {

                if(logging) {
                    output("Error: Degree matrix in graph::is_even has not been initialized");
                }

                return 0;

            }

            if (query == "all") {

                for (int j = 0; j < deg_mat.n_cols; ++j) {

                    for (int i = 0; i < deg_mat.n_rows; ++i) {

                       if (deg_mat(i,j)%2 == 1) {
                            return 0;
                        }

                    }
                }

                return 1;

            }
            else if (query == "layer"){

                if (layer > num_layers-1) {

                    if(logging) {
                        output("Layer specified in graph::is_even is out of bounds");
                    }
                    
                    return 0;

                }

                for (int i = 0; i < deg_mat.n_rows; ++i) {

                   if (deg_mat(i,layer)%2 == 1) {
                        return 0;
                    }

                }

                return 1;
                
            }
            else {

                if(logging) {
                    output("Query in graph::is_even must be all or layer");
                }

                return 0;

            }

        }

        //determine whether or not a graph layer contains a cut-edge (West Exercise 1.3.12)
        bool contains_cut_edge(const std::string query = "all", const unsigned int layer = 0) {
            
            if (deg_mat.is_empty()) {

                if(logging) {
                    output("Error: Degree matrix in graph::contains_cut_edge has not been initialized");
                }

                return 0;

            }

            if (query == "all") {

                if (is_even()) {

                    return 0;

                }
                else {

                    if(logging) {
                        output("Unable to determine whether or not graph contains cut-edge");
                    }   

                    return 1;   

                }

            }
            else if (query == "layer"){

                if (layer > num_layers-1) {

                    if(logging) {
                        output("Layer specified in graph::contains_cut_edge is out of bounds");
                    }
                    
                    return 0;

                }

                if (is_even(query, layer)) {
                    return 0;
                }
                else {

                    if(logging) {
                        output("Unable to determine whether or not graph contains cut-edge");
                    }   

                    return 1;   

                }
                
            }
            else {

                if(logging) {
                    output("Query in graph::contains_cut_edge must be all or layer");
                }

                return 0;

            }
        }

        //determine whether or not a graph layer contains a cycle
        bool contains_cycle(const std::string query = "all", const unsigned int layer = 0) {

            if (deg_mat.is_empty()) {

                if(logging) {
                    output("Error: Degree matrix in graph::contains_cycle has not been initialized");
                }

                return 0;

            }

            if (query == "all") {

                //West Exercise 1.2.38
                if (size() >= order()) {
                    return 1;
                }

                //West Chapter 1.2.25
                unsigned int count = 0;

                //loop layer-first
                for (int j = 0; j < deg_mat.n_cols; ++j) {

                    for (int i = 0; i < deg_mat.n_rows; ++i) {

                        if ((int)deg_mat(i, j) < 2) {
                            count++;
                        }

                    }

                }

                if (count == 0) {
                    return 1;
                }
                else {

                    if(logging) {
                        output("Note: Graph may still contain cycle.");
                    }
                    
                    return 0;

                }

            }
            else if (query == "layer"){
                
                if (layer > num_layers-1) {

                    if(logging) {
                        output("Layer specified in graph::contains_cycle is out of bounds");
                    }
                    
                    return 0;

                }

                //West Exercise 1.2.38
                if (size("layer", layer) >= order()) {
                    return 1;
                }

                //West Chapter 1.2.25
                unsigned int count = 0;

                for (int i = 0; i < deg_mat.n_rows; ++i) {

                    if ((int)deg_mat(i, layer) < 2) {
                        count++;
                    }

                }

                if (count == 0) {
                    return 1;
                }
                else {

                    if(logging) {
                        output("Note: Graph may still contain cycle.");
                    }

                    return 0;
                }
            }
            else {

                if(logging) {
                    output("Query in graph::contains_cycle must be all or layer");
                }

                return 0;

            }

        }

        //determine whether or not a graph layer contains a closed euler walk 
        //Trudeau (173)
        bool contains_closed_euler_walk(const std::string query = "all", const unsigned int layer = 0) {

            if (deg_mat.is_empty()) {

                if(logging) {
                    output("Error: Degree matrix in graph::contains_closed_euler_walk has not been initialized");
                }

                return 0;

            }

            if (query == "all") {

                unsigned int count = 0;

                for (int i = 0; i < deg_mat.n_rows; ++i) {

                    if ((int)deg_mat(i, layer)%2 == 1) {
                        return 0;
                    }

                }

                //loop layer-first
                for (int j = 0; j < deg_mat.n_cols; ++j) {

                    for (int i = 0; i < deg_mat.n_rows; ++i) {

                        if ((int)deg_mat(i, j)%2 == 1) {
                            return 0;
                        }

                    }

                }

                if (is_connected("layer", layer)) {
                    return 1;
                }
                else {
                    return 0;
                }

            }
            else if (query == "layer"){

                if (layer > num_layers-1) {

                    if(logging) {
                        output("Layer specified in graph::contains_closed_euler_walk is out of bounds");
                    }
                    
                    return 0;

                }

                unsigned int count = 0;

                for (int i = 0; i < deg_mat.n_rows; ++i) {

                    if ((int)deg_mat(i, layer)%2 == 1) {
                        return 0;
                    }

                }

                if (is_connected("layer", layer)) {
                    return 1;
                }
                else {
                    return 0;
                }

            }
            else {

                if(logging) {
                    output("Query in graph::contains_closed_euler_walk must be all or layer");
                }

                return 0;

            }

        }

        //determine whether or not a graph layer contains an open euler walk
        //Trudeau (175)
        bool contains_open_euler_walk(const std::string query = "all", const unsigned int layer = 0) {

            if (deg_mat.is_empty()) {

                if(logging) {
                    output("Error: Degree matrix in graph::contains_open_euler_walk has not been initialized");
                }

                return 0;

            }

            if (query == "all") {

                unsigned int count = 0;

                for (int j = 0; j < deg_mat.n_cols; ++j) {

                    for (int i = 0; i < deg_mat.n_rows; ++i) {

                        if ((int)deg_mat(i, j)%2 == 1) {
                            count++;
                        }

                        if (count > 2) {
                            return 0;
                        }

                    }

                }

                if (is_connected() && count == 2) {
                    return 1;
                }
                else {
                    return 0;
                }

            }
            else if (query == "layer"){
            
                if (layer > num_layers-1) {

                    if(logging) {
                        output("Layer specified in graph::contains_open_euler_walk is out of bounds");
                    }
                    
                    return 0;

                }

                unsigned int count = 0;

                for (int i = 0; i < deg_mat.n_rows; ++i) {

                    if ((int)deg_mat(i, layer)%2 == 1) {
                        count++;
                    }

                    if (count > 2) {
                        return 0;
                    }

                }

                if (is_connected(query, layer) && count == 2) {
                    return 1;
                }
                else {
                    return 0;
                }

            }
            else {

                if(logging) {
                    output("Query in graph::contains_open_euler_walk must be all or layer");
                }

                return 0;

            }

        }

        //determine whether or not a graph layer contains an euler walk
        bool contains_euler_walk(const std::string query = "all", const unsigned int layer = 0) {

            if (deg_mat.is_empty()) {

                if(logging) {
                    output("Error: Degree matrix in graph::contains_euler_walk has not been initialized");
                }

                return 0;
            }

            if (query == "all") {

                if (contains_closed_euler_walk() == 1) {

                    if(logging) {
                        output("Graph contains closed euler walk");
                    }
                    
                    return 1;

                }
                else if (contains_open_euler_walk() == 1) {

                    if(logging) {
                        output("Graph contains open euler walk");
                    }

                    return 1;

                }
                else {

                    if(logging) {
                        output("Graph contains euler walk");
                    }

                    return 0;

                }
            }
            else if (query == "layer"){

                if (layer > num_layers-1) {

                    if (logging) {
                        output("Layer specified in contains_euler_walk is out of bounds");
                    }
                    return 0;
                }

                if (contains_closed_euler_walk(query, layer) == 1) {

                    if(logging) {
                        output("Graph contains closed euler walk");
                    }

                    return 1;

                }
                else if (contains_open_euler_walk(query, layer) == 1) {

                    if(logging) {
                        output("Graph contains open euler walk");
                    }

                    return 1;

                }
                else {

                    if(logging) {
                        output("Graph contains euler walk");
                    }

                    return 0;

                }

            }
            else {
                
                if (logging) {
                    output("Query in contains_euler_walk must be all or layer");
                }

                return 0;
            
            }

        }

        bool weighted() {

            if (edge_weights.is_empty()) {
                return 0;
            }
            else {
                return 1;
            }

        }

        //initialize graph with preset properties
        void initialize(const std::string type = "complete", const std::string layout = "polar", const double angle = 0) {
            
            unsigned int n = order();
            arma::vec rad_vec(n);

            rad_vec.fill(0.05);

            set_type();
            create_vertex_set(rad_vec);
            assign_coords(layout, r, angle);
            create_adjacency_tensor(type);

        }

};

class digraph: public graph {

    public:

        bool is_directed = 1;                    //ensure graph is directed

        digraph(

            //graph attributes
            const int ID, 
            arma::uvec vertex_set_sizes,
            double radius,
            arma::vec center,

            //vertex attributes
            const int bound,

            //edge attributes
            unsigned int number_of_layers
            )

                : graph(ID, vertex_set_sizes, radius, center, bound, number_of_layers) {} 

};

//graph attributes

//determine the minimum number of components contained in a graph (West Chapter 1.2.11)
int min_num_components(graph G) {

    if (G.v_sizes.is_empty() || G.deg_mat.is_empty()) {

        if (logging) {
            output("Graph used in min_num_components has not been properly initialized");
        }

        return 0;

    }
    
    int num = G.order() - G.size();
    
    if(num < 0) {
        return 0;
    }
    else {
        return num;
    }

}

//return the minimum number of decomposing trails of a graph with exactly 2k odd vertices (West Chapter 1.2.33)
int min_num_decomp_trails(graph G) {

    int count = 0;

    for (int i = 0; i < G.v.size(); ++i) {

        if (G.v[i]->is_odd()) {
            count++;
        }

    }

    if (G.is_connected() && count%2 == 0){
        return std::max(1, count/2);
    }
    else {

        if (logging) {
            output("Criteria not met in min_num_decomp_trails");\
        }

        return 0;

    }

}

unsigned int num_uv_walks(const graph G, const vertex<double>* u, const vertex<double>* v, const unsigned int layer, const unsigned int length) {

    if (u == nullptr || v == nullptr) {

        if (logging) {
            output("One or more vertices in num_uv_walks has not been initialized");
        }

        return 0;

    }

    if (G.adjacency_tensor.is_empty()) {

        if (logging) {
            output("Adjacency tensor of graph in num_uv_walks has not been initialized");
        }

        return 0;

    }

    if (length < 1) {
        
        if (logging) {
            output("Length in num_uv_walks must be a positive integer");
        }

        return 0;

    }

    if (length > 1) {

        arma::mat prod = arma::eye(G.adjacency_tensor.n_rows, G.adjacency_tensor.n_cols);
        
        for (int i = 0; i < length; ++i) {
            prod *= G.adjacency_tensor.slice(layer); 
        }

        return prod(u->index, v->index);
    
    }
    else {
        return G.adjacency_tensor.slice(layer)(u->index, v->index);
    }

}

//return the total number of graphs of order n
long long num_graphs_of_order_n(const unsigned int order, const unsigned int num_layers) {

    long long e = num_layers*order*(order-1)/2;
    long long N = pow(2, e);

    return N;
}

//graph operators

graph graph_union(graph F, graph G) {

    if (F.v.empty() || F.adjacency_tensor.is_empty()) {

        if(logging) {
            output("Warning: Graph F used in graph_union has not been properly initialized");
        }

        return G;
    }

    if (G.v.empty() || G.adjacency_tensor.is_empty()) {

        if(logging) {
            output("Warning: Graph G used in graph_union has not been properly initialized");
        }

        return F;
    }
    
    //merge vertex sets
    arma::umat vertices = arma::join_cols(F.V,G.V);
    arma::uvec u_vertices = arma::unique(vertices.col(0));

    //sort and reduce vertex set
    insertion_sort(vertices, "F", true);
    
    //define inputs
    unsigned int num_layers = std::max(F.num_layers, G.num_layers);
    arma::uvec vertex_set_sizes = arma::zeros<arma::uvec>(u_vertices.size());
    arma::vec center = (F.x+G.x)/2;

    //create graph
    graph H(0, vertex_set_sizes, 0.4, center, 100, num_layers);

    //set graph type
    H.set_type();

    //specify whether or not loops are permitted
    if (F.loops_permitted || G.loops_permitted) {
        H.contains_loops();
    }

    //ensure vertices with the same species/IDs are identical
    for (int i = 0; i <  vertices.n_rows; ++i) {
        
        unsigned int S = vertices(i, 0);
        unsigned int I = vertices(i, 1);
        
        if(F.get_vertex(S,I) != nullptr && G.get_vertex(S,I) != nullptr) {

            if (F.get_vertex(S,I)->r != G.get_vertex(S,I)->r) {

                if(logging) {
                    output("Warning: Conflict in graph_union; vertices with matching identifiers have mismatching radii.");
                }

                H.v_sizes(0)++;
                
                H.initialize("null");
                
                return H;
            }
            /**
            if (arma::sum(F.get_vertex(S,I)->x != G.get_vertex(S,I)->x) > 0) {

                if (logging) {
                    output("Warning: Conflict in graph_union; vertices with matching identifiers have mismatching coordinates.");
                }

                H.v_sizes(0)++;
                
                H.initialize("null");
                
                return H;
            
            }
            /**/

        }

    }

    //determine size of partite vertex sets

    unsigned int index = 0;

    for (int i = 0; i < vertices.n_rows; ++i) {

        if (vertices(i, 0) == u_vertices(index)) {
            H.v_sizes(index)++;
        }
        else {
            index++;
            i -= 1;
        }

    }

    //create vertex set
    H.V = vertices;

    //adjust vertex set
    for (int i = 0; i < H.V.n_rows; ++i) {

        unsigned int S = vertices(i, 0);
        unsigned int I = vertices(i, 1);

        if (F.get_vertex(S,I) != nullptr) {
            H.v.push_back(F.get_vertex(S,I));
            H.v[i]->index = i;
            H.v[i]->degree = arma::zeros<arma::urowvec>(H.num_layers);
        }
        else {
            H.v.push_back(G.get_vertex(S,I));
            H.v[i]->index = i;
            H.v[i]->degree = arma::zeros<arma::urowvec>(H.num_layers);
        }

    }

    //allocate space within adjacency tensor
    H.adjacency_tensor = F.adjacency_tensor;

    unsigned int count = H.adjacency_tensor.n_rows;
    arma::cube adj_ten_temp = arma::zeros(count, count, H.num_layers);

    for (int i = 0; i <  H.v.size(); ++i) {

        unsigned int S = H.v[i]->sid;
        unsigned int I = H.v[i]->id;
        arma::mat adj_ten_slice;

        if (F.get_vertex(S,I) == nullptr) {
            
            count++;

            adj_ten_temp.resize(count, count, num_layers);

            for (int k = 0; k < H.num_layers; ++k){

                adj_ten_slice = H.adjacency_tensor.slice(k);

                if (G.get_vertex(S,I)->sid < G.v_sizes.size()) {

                    adj_ten_slice.insert_rows(i, arma::zeros<arma::rowvec>(adj_ten_slice.n_cols));
                    adj_ten_slice.insert_cols(i, arma::zeros(adj_ten_slice.n_rows));

                    adj_ten_temp.slice(k) = adj_ten_slice;

                }
                else {
                        
                    adj_ten_slice = arma::join_cols(adj_ten_slice, arma::zeros<arma::rowvec>(adj_ten_slice.n_cols));
                    adj_ten_slice = arma::join_rows(adj_ten_slice, arma::zeros(adj_ten_slice.n_rows)); 
                    
                    adj_ten_temp.slice(k) = adj_ten_slice;
                
                }

            }

        H.adjacency_tensor.reset();

        H.adjacency_tensor = adj_ten_temp; 

        }
    
    }

    //fill in adjacency tensor with missing edges
    unsigned int k = 0;
    unsigned int increment = 1;

    if (H.loops_permitted) {
        increment = 0;
    }

    for (int layer = 0; layer < H.num_layers; ++layer) {

        for (int i = 0; i < H.adjacency_tensor.n_rows-increment; ++i) {

            unsigned int l = k+1;

            for (int j = i+increment; j < H.adjacency_tensor.n_cols; ++j) {

                if (G.get_vertex(vertices(i, 0),vertices(j, 1)) != nullptr) {
                    
                    if (G.get_vertex(vertices(j, 0),vertices(j, 1)) != nullptr && l < G.order()) {
                        
                        if (G.adjacency_tensor(k, l, layer) > 0 && H.adjacency_tensor(i,j, layer) == 0) {
                            H.adjacency_tensor(i, j, layer) = G.adjacency_tensor(k, l, layer);
                            H.adjacency_tensor(j, i, layer) = H.adjacency_tensor(i, j, layer);
                        }

                        l++;

                    }

                }

            }

            if (G.get_vertex(vertices(i,0), vertices(i, 1)) != nullptr) {
                k++;
            }

        }

    }

    //assign vertex degrees

    H.update_degrees();

    //create edge set
    H.update_edge_set();

    return H;

}

graph graph_sum(graph F, graph G) {

    if (F.v.empty() || F.adjacency_tensor.is_empty()) {

        if(logging) {
            output("Warning: Graph F used in graph_sum has not been properly initialized");
        }

        return G;

    }

    if (G.v.empty() || G.adjacency_tensor.is_empty()) {

        if(logging) {
            output("Warning: Graph G used in graph_sum has not been properly initialized");
        }

        return F;

    }

    //merge vertex sets
    arma::umat vertices = arma::join_cols(F.V,G.V);
    arma::uvec u_vertices = arma::unique(vertices.col(0));

    //sort vertex set
    insertion_sort(vertices, "F", false);

    if (vertices.n_rows != F.order() + G.order()) {

        if(logging) {
            output("Warning: Vertex sets used in graph_sum are not disjoint");
        }

        return F;

    }

    //define inputs
    unsigned int num_layers = std::max(F.num_layers, G.num_layers);
    arma::uvec vertex_set_sizes = arma::zeros<arma::uvec>(u_vertices.size());
    arma::vec center = (F.x+G.x)/2;

    //create graph
    graph H(0, vertex_set_sizes, 0.4, center, 100, num_layers);

    //set graph type
    H.set_type();

    //specify whether or not loops are permitted
    if (F.loops_permitted || G.loops_permitted) {
        H.contains_loops();
    }

    //determine size of partite vertex sets
    unsigned int index = 0;

    for (int i = 0; i < vertices.n_rows; ++i) {

        if (vertices(i, 0) == u_vertices(index)) {
            H.v_sizes(index)++;
        }
        else {
            index++;
            i -= 1;
        }

    }

    //create vertex set
    H.V = vertices;

    //adjust vertex set
    for (int i = 0; i < H.V.n_rows; ++i) {

        unsigned int S = vertices(i, 0);
        unsigned int I = vertices(i, 1);

        if (F.get_vertex(S,I) != nullptr) {
            H.v.push_back(F.get_vertex(S,I));
            H.v[i]->index = i;
        }
        else {
            H.v.push_back(G.get_vertex(S,I));
            H.v[i]->index = i;
        }

    }

    //allocate space within adjacency tensor
    H.adjacency_tensor = F.adjacency_tensor;

    unsigned int count = H.adjacency_tensor.n_rows;
    arma::cube adj_ten_temp = arma::zeros(count, count, H.num_layers);

    for (int i = 0; i <  H.v.size(); ++i) {

        unsigned int S = H.v[i]->sid;
        unsigned int I = H.v[i]->id;
        arma::mat adj_ten_slice;

        if (F.get_vertex(S,I) == nullptr) {
            
            count++;

            adj_ten_temp.resize(count, count, num_layers);

            for (int k = 0; k < H.num_layers; ++k){

                adj_ten_slice = H.adjacency_tensor.slice(k);

                if (G.get_vertex(S,I)->sid < G.v_sizes.size()) {

                    adj_ten_slice.insert_rows(i, arma::zeros<arma::rowvec>(adj_ten_slice.n_cols));
                    adj_ten_slice.insert_cols(i, arma::zeros(adj_ten_slice.n_rows));

                    adj_ten_temp.slice(k) = adj_ten_slice;

                }
                else {
                        
                    adj_ten_slice = arma::join_cols(adj_ten_slice, arma::zeros<arma::rowvec>(adj_ten_slice.n_cols));
                    adj_ten_slice = arma::join_rows(adj_ten_slice, arma::zeros(adj_ten_slice.n_rows)); 
                    
                    adj_ten_temp.slice(k) = adj_ten_slice;
                
                }

            }

        H.adjacency_tensor.reset();

        H.adjacency_tensor = adj_ten_temp; 

        }
    
    }

    //fill in adjacency tensor with missing edges
    unsigned int k = 0;
    unsigned int increment = 1;

    if (H.loops_permitted) {
        increment = 0;
    }

    for (int layer = 0; layer < H.num_layers; ++layer) {

        for (int i = 0; i < H.adjacency_tensor.n_rows-increment; ++i) {

            unsigned int l = k+1;

            for (int j = i+increment; j < H.adjacency_tensor.n_cols; ++j) {

                if (G.get_vertex(vertices(i, 0),vertices(j, 1)) != nullptr) {
                   
                    if (G.get_vertex(vertices(j, 0),vertices(j, 1)) != nullptr && l < G.order()) {
                        
                        if (G.adjacency_tensor(k, l, layer) > 0 && H.adjacency_tensor(i,j, layer) == 0) {
                            H.adjacency_tensor(i, j, layer) = G.adjacency_tensor(k, l, layer);
                            H.adjacency_tensor(j, i, layer) = H.adjacency_tensor(i, j, layer);
                        }

                        l++;

                    }

                }

            }

            if (G.get_vertex(vertices(i,0), vertices(i, 1)) != nullptr) {
                k++;
            }

        }

    }


    //assign vertex degrees
    H.update_degrees();

    //create edge set
    H.update_edge_set();

    return H;

}

//vertex attributes

//edge attributes

//vertex degree attributes

//connectivity attributes and algorithms

//determine whether or not two vertices are adjacent
bool adjacent(const graph G, const vertex<double> * u, const vertex<double> * v, const unsigned int layer = 0) {

    if (G.adjacency_tensor.is_empty()) {

        if(logging) {
            output("Adjacency tensor in adjacent has not been initialized");
        }

        return 0;
    }

    if (layer > G.num_layers-1) {

        if(logging) {
            output("Layer specified in adjacent is out of bounds");
        }
        
        return 0;
    }

    unsigned int i = u->index;
    unsigned int j = v->index;
    
    if(G.adjacency_tensor(i, j, layer) == 1) {
        return 1;
    }
    else {
        return 0;
    }
    
}

std::vector< vertex<double> * > neighborhood(const graph G, const vertex<double> * v, const unsigned int layer = 0) {

    if (G.v.empty()) {

        if(logging) {
            output("Warning: Vertex set used in neighborhood has not been initialized");
        }

        std::vector< vertex<double> * > exit = {nullptr};

        return exit;
    }

    std::vector< vertex<double> * > N;

    for (int i = 0; i <  G.v.size(); ++i) {

        if (G.v[i] != v && adjacent(G, G.v[i], v, layer)) {
            N.push_back(G.v[i]);
        }

    }

    return N;
}

//graph adjustment algorithms

//return the complement of graph G
graph complement(const graph G) {

    if (G.adjacency_tensor.is_empty()) {

        if(logging) {
            output("Adjacency tensor in complement has not been initialized");
        }

        return G;

    }

    if (G.type != "simple") {

        if(logging) {
            output("Graph type in complement must be simple");
        }

        return G;

    }

    graph G_c = G;
    G_c.complement();

    return G_c;
}

//return graph induced by T
graph induced_graph(const graph G, arma::umat T) {

    if (G.v.empty() || G.adjacency_tensor.is_empty()) {

        if(logging) {
            output("Warning: Graph set used in induced_graph has not been properly initialized");
        }

        return G;
    }

    if (T.is_empty()) {

        if(logging) {
            output("Warning: Vertex set entered in induced_graph has not been initialized");
        }

        return G;
    }

    graph G_T = G;

    insertion_sort(T, "forward", 1);

    unsigned int count_T = 0;
    unsigned int count_remove = 0;
    unsigned int new_order = T.n_rows;
    std::vector < vertex<double> *> w;
    arma::cube adj_ten_temp(new_order, new_order, G_T.num_layers);

    //adjust vertex set
    for (int i = 0; i <  G_T.v.size(); ++i) {

        if (count_T < new_order && T(count_T, 0) == G_T.v[i]->sid && T(count_T, 1) == G_T.v[i]->id) {
            
            G_T.v[i]->index  = i - count_remove;
            w.push_back(G_T.v[i]);
            count_T++;

        }
        else {

            G_T.v_sizes(G_T.v[i]->sid) -= 1;
            G_T.V.shed_row(i - count_remove);
            count_remove++;

        }

    }

    //adjust adjacency tensor
    for (int k = 0; k < G_T.num_layers; ++k){

        count_T = 0;
        count_remove = 0;

        arma::mat adj_ten_slice = G_T.adjacency_tensor.slice(k);

        for (int i = 0; i <  G_T.v.size(); ++i) {

            if (count_T < new_order && T(count_T, 0) == G_T.v[i]->sid && T(count_T, 1) == G_T.v[i]->id) {
                count_T++;
            }
            else {

                adj_ten_slice.shed_row(i-count_remove);
                adj_ten_slice.shed_col(i-count_remove);
                count_remove++;

            }

        }

        adj_ten_temp.slice(k) = adj_ten_slice;

    }

    G_T.adjacency_tensor.reset();

    G_T.adjacency_tensor = adj_ten_temp; 

    G_T.v.clear();
    G_T.v = w;

    G_T.update_degrees();

    //G_T.create_incidence_matrix();

    G_T.update_edge_set();

    return G_T;

}

//perform 2-switch on graph layer (West Chapter 1.3.32)
void two_switch(graph * G, const vertex<double> * u, const vertex<double> * v, const vertex<double> * x, const vertex<double> * y, const unsigned int layer = 0) {
    
    //ensure conditions for 2-switch are met
    if (!adjacent(*G, u, v, layer)) {
        
        if (logging) {
            std::cout << "Graph used in two_switch does not contain edge from vertex " << u->sid <<  ' ' << u->id << " to  vertex " << v->sid <<  ' ' << v->id << '\n';
        }
        
        return;

    }

    if (!adjacent(*G, x, y, layer)) {
        
        if (logging) {
            std::cout << "Graph used in two_switch does not contain edge from vertex " << x->sid <<  ' ' << x->id << " to  vertex " << y->sid <<  ' ' << y->id << '\n';
        }

        return;
    
    }

    if (adjacent(*G, v, x, layer)) {

        if (logging) {
            std::cout << "Graph used in two_switch contains edge from vertex " << v->sid <<  ' ' << v->id << " to  vertex " << x->sid <<  ' ' << x->id << '\n';
        }

        return;
    
    }

    if (adjacent(*G, y, u, layer)) {

        if (logging) {
           std::cout << "Graph used in two_switch contains edge from vertex " << y->sid <<  ' ' << y->id << " to  vertex " << u->sid <<  ' ' << u->id << '\n'; 
        }
        
        return;

    }

    //gather indices of input vertices
    unsigned int i = u->index;
    unsigned int j = v->index;

    unsigned int k = x->index;
    unsigned int l = y->index;

    //adjust adjacency tensor
    G->adjacency_tensor(i, j, layer) = 0;
    G->adjacency_tensor(j, i, layer) = 0;

    G->adjacency_tensor(k, l, layer) = 0;
    G->adjacency_tensor(l, k, layer) = 0;

    //remove edges from edge set
    bool check1 = 0;
    bool check2 = 0;

    for (int row = 0; row < G->E.n_rows; ++row) {
        
        if (G->E(row, 0) == i && G->E(row, 1) == j) {
            G->E.shed_row(row);
            check1 = 1;
        }

        if (G->E(row, 0) == k && G->E(row, 1) == l) {
            G->E.shed_row(row);
            check2 = 1;
        }

        if (check1 && check2) {
            break;
        }

    }

    //adjust adjacency tensor
    G->adjacency_tensor(j, k, layer) = 1;
    G->adjacency_tensor(k, j, layer) = 1;

    G->adjacency_tensor(l, i, layer) = 1;
    G->adjacency_tensor(i, l, layer) = 1;

    //create new nedges
    arma::urowvec vx = {j, k};
    arma::urowvec yu = {l, i};

    //add edges to edge set
    G->E = arma::join_cols(G->E, vx);
    G->E = arma::join_cols(G->E, yu);

    insertion_sort(G->E, "F");

}

graph underlying_graph(graph G) {

    if (!G.is_directed) {

        if (logging) {
            output("Graph entered in underlying_graph must be directed");
        }

        return G;

    }

    return G;

}

//graph comparators

//determine whether or not two graphs are not isomorphic
bool non_isomorphic(graph G, graph H) {

    if (G.v_sizes.is_empty() || G.deg_mat.is_empty()) {

        if (logging) {
            output("Graph G used in non_isomorphic has not been properly initialized");
        }

        return 0;

    }    

    if (H.v_sizes.is_empty() || H.deg_mat.is_empty()) {

        if (logging) {
            output("Graph H used in non_isomorphic has not been properly initialized");
        }

        return 0;

    } 

    arma::umat G_degs = G.deg_mat;
    arma::umat H_degs = H.deg_mat;

    insertion_sort(G_degs, "F");
    insertion_sort(H_degs, "F");

    if (G.order() != H.order()) {
        return 1;
    }
    else if (G.size() != H.size()) {
        return 1;
    }
    else if (arma::accu(G_degs != H_degs) != 0) {
        return 1;
    }
    else if (G.type == "simple" && H.type == "simple") {

        graph G_c = complement(G);
        graph H_c = complement(H);

        arma::umat G_c_degs = G_c.deg_mat;
        arma::umat H_c_degs = H_c.deg_mat;

        insertion_sort(G_c_degs, "F");
        insertion_sort(H_c_degs, "F");

        if (G_c.size() != H_c.size()) {
            return 1;
        }
        else if (arma::accu(G_c_degs != H_c_degs) != 0) {
            return 1;
        }
        else {

            if (logging) {
                output("Unable to determine whether or not graphs are nonisomorphic");
            }

            return 0;
        }
    }
    else {

        if (logging) {
            output("Unable to determine whether or not graphs are nonisomorphic");
        }

        return 0;

    }

}

//graph qualifiers

//determine whether or not a graph is not self-complementary (West Exercise 1.1.31)
bool non_self_complementary(graph G) {

    if (G.v_sizes.is_empty()) {

        if (logging) {
            output("Graph in non_self_complementary has not been intialized");
        }

        return 1;

    }

    if (G.type != "simple") {

        if (logging) {
            output("Graph in non_self_complementary must be simple");
        }

        return 1;

    }

    unsigned int n = G.order();

    if (n%4 != 0 && (n-1)%4 != 0) {
        return 1;
    }

    graph G_c = complement(G);

    if (non_isomorphic(G, G_c)) {
        return 1;
    }
    else {
        return 0;
    }

}

//determine if integer sequence is graphic (West Exercise 1.3.63/West Exercise 1.3.57)
bool is_graphic(const arma::umat deg_mat, const std::string size, const bool contains_loops = 0) {

    if (deg_mat.is_empty()) {

        if (logging) {
            output("Degree matrix in is_graphic has not been initialized");
        }

        return 0;

    }

    unsigned int tally1 = 0;
    unsigned int tally2 = 0;
    arma::uvec deg_seq = arma::sum(deg_mat, 1);

    if (size == "small" || size == "S") {

        if (arma::accu(deg_mat)%2 == 0) {
            tally1++;
            tally2++;
        }
        else {
            return 0;
        }

    }
    else if (size == "large" || size == "L") {

        unsigned int count = 0;
        
        for (int i = 0; i < deg_seq.size(); ++i) {
        
            if ((int)deg_seq(i)%2 == 1) {
                count++;
        
            }
        }

        if (count%2 == 0) {
            tally1++;
            tally2++;
        }
        else {
            return 0;
        }

    }
    else {

        if (logging) {
            output("Must enter a valid qualifier in is_graphic");
        }

        return 0;

    }

    if (contains_loops && tally1 == 1) {
        return 1;
    }

    unsigned int count = 0;
    unsigned int n = deg_seq.size();
    unsigned int Delta = deg_seq.max();

    if (Delta < n) {
        tally2++;
    }

    for (int i = 0; i < n; ++i) {

        if (deg_seq(i) == Delta || deg_seq(i) == Delta - 1) {
            count++;
        }

    }

    if (count == n && tally2 == 2) {
        return 1;
    }

    insertion_sort(deg_seq, "B");

    if (arma::accu(deg_seq) - (2*deg_seq(0)) >= 0) {
        tally1++;
    }

    if (tally1 == 2) {
        return 1;
    }
    else {
        return 0;
    }

}

//Havel-Hakimi algorithm (West Chapter 1.3.31)
//determine if integer sequence is graphic
bool havel_hakimi(const arma::uvec deg_seq, const unsigned int n) {

    if (deg_seq.is_empty()) {

        if (logging) {
            output("Degree sequence in havel_hakimi has not been initialized");
        }

        return 0;

    }
    
    //ensure inputs realistic
    if (n > deg_seq.size() || n < 1) {

        if (logging) {
            output("Error: Second argument in havel_hakimi must be between 1 and the length of the first argument");
        }

        return 0;
    
    }

    //initialize algorithm components
    arma::uvec sds = arma::sort(deg_seq, "descend");

    double deg_max = sds(0);

    //ensure input vector will not cause errors
    if (deg_max > sds.size()-1) {
        return 0;
    }

    //choose algorithm within is_graphic
    std::string qualifier = "small";
    
    if (deg_max > 1e6 || deg_seq.size() > 1e3) {
        qualifier = "large";
    }

    bool cond = 1;

    while (cond) {

        if (sds.size() <= n) {

            if (is_graphic(sds, qualifier)) {
                cond = 0;
                return 1;
            }
            else {
                cond = 0;
                return 0;
            }
        
        }

        sds.shed_row(0);

        for (int i = 0; i < deg_max; ++i) {
            sds(i) -= 1;
        }
        
        sds = arma::sort(sds, "descend");

        deg_max = sds(0);
    
    }

}
