#include "static/cpp/graph_theory.h"
#include <string.h>
#include <jsoncpp/json/json.h>

int string_to_int(std::string str) {
	
	std::string::size_type size;
	
	return std::stoi(str, &size);

}

int string_to_double(std::string str) {
	
	std::string::size_type size;
	
	return std::stod(str, &size);

}

int main(int argc, char * argv[]) {

	srand(time(NULL));

	if (argc < 2) {

		std::cout << "Error. Please provide more input." << '\n';

		return 0;
	}

	//vertex settings

	char * n_sp = strtok(argv[1], "+");
	unsigned int num_species = string_to_int(n_sp);

	//output(num_species);


	arma::uvec vertex_set_sizes(num_species);

	for (int i = 0; i < num_species; ++i) {
		char * v_s_s = strtok(NULL, "+");
		//std::cout << *v_s_s << '\n';
		vertex_set_sizes[i] = string_to_int(v_s_s);
	}

	//vertex_set_sizes.print();

/**/
	//display settings

	char * dim = strtok(NULL, "+");
	unsigned int dimension = string_to_int(dim);
	std::string position = strtok(NULL, "+");

	//output(dimension);
	//output(position);

	arma::vec center(dimension);

	if (position == "origin") {
		center = arma::zeros(dimension);
	}
	else {

		for (int i = 0; i < dimension; ++i) {
			char * cen = strtok(NULL, "+");
			center[i] = string_to_int(cen);
		}

	}

	//center.print();

	char * rad = strtok(NULL, "+");
	double radius = string_to_double(rad);
	std::string layout = strtok(NULL, "+");
	std::string type = strtok(NULL, "+");

	//output(radius);
	//output(layout);
	//output(type);
/**/
	graph g(0, vertex_set_sizes, radius, center, 100, 1);

	g.initialize(type, layout);

	unsigned int n = g.order();

	Json::Value graph;
	
	Json::Value& adjacency_tensor = graph["adjacency_tensor"];

	graph["dimension"] = dimension;
	graph["order"] = n;
	graph["size"] = g.size();
	graph["radius"] = radius;

/**/
	for (int i = 0; i < n; ++i) {
		
		for (int j = 0; j < n; ++j) {

			adjacency_tensor[i]["degree_sequence"].append(g.adjacency_tensor(i, j, 0));

		}

	}
/**/
	Json::Value& vertices = graph["vertices"];

	for (int i = 0; i < n; ++i) {

		vertices[i]["sid"] = g.v[i]->sid;
		vertices[i]["id"] = g.v[i]->id;
		vertices[i]["r"] = g.v[i]->r;
		vertices[i]["x"] = g.v[i]->x(0);
		vertices[i]["y"] = g.v[i]->x(1);
		vertices[i]["z"] = g.v[i]->x(2);
		//g.v[i]->x.print();
	}
/**/
	output(graph);
/**/	
	return 0;

}