#include "qif"

namespace qif {

Graph::Graph(uint vertex_num, std::string& edges) {
	//traslate the string to a vector of integer 1pairs
	std::vector < std::pair<int, int> > new_edges;
	std::string copy = edges;
	std::string current;
	uint pos, pairs = 0;
	uint i = 0;
	while(i < copy.size()) {
		pos = copy.find(";");
		current = copy.substr(i, pos);
		i = pos + 1;
		//vertex_num++;
		while(current[0] == ' ') {
			current = current.substr(1);
		}
		pos = current.find(" ");
		int left = std::atoi(current.substr(0, pos).c_str());
		int rigth = std::atoi(current.substr(pos).c_str());
		std::pair<int, int> new_pair = std::pair<int, int>(left, rigth);
		new_edges.push_back(new_pair);
		pairs++;
		if(i == 0) {
			copy = "";
		} else {
			copy = copy.substr(i);
		}
		i = 0;
	}

	Graph(vertex_num, new_edges);
}

Graph::Graph(uint vertex_num, std::vector< std::pair<int, int> >& edges) {
	//std::cerr << "HOME is not defined." << std::endl;
	V = vertex_num;
	adjacency = arma::mat(V, V);
	adjacency.zeros();
	//std::cerr << "HOME is not defined." << std::endl;
	/*
	for(int i=0; i<V; i++) {
	    for(int j=0; j<V; j++) {
	        adjacency(i,j) = 0;
	    }
	}*/
	//std::cerr << "HOME is not defined." << std::endl;
	for(uint i = 0; i < edges.size(); i++) {
		uint v1 = edges[i].first;
		uint v2 = edges[i].second;
		if(v1 < 1 || v1 > V || v2 < 1 || v2 > V) throw std::runtime_error("error");

		adjacency(v1 - 1, v2 - 1) = 1;
		adjacency(v2 - 1, v1 - 1) = 1;
	}

	distances = arma::mat(V, V);
	// Just use Floyd-Warshall

	// Initialization (we use V+1 as infinity)
	for(uint i = 0; i < V; i++) {
		for(uint j = 0; j < V; j++) {
			if(i == j) distances(i, j) = 0;
			else if(adjacency(i, j) == 1) distances(i, j) = 1;
			else distances(i, j) = V + 1;
		}
	}

	// Main algorithm
	for(uint k = 0; k < V; k++) {
		for(uint i = 0; i < V; i++) {
			for(uint j = 0; j < V; j++) {
				if(distances(i, j) > distances(i, k) + distances(k, j)) {
					distances(i, j) = distances(i, k) + distances(k, j);
				}
			}
		}
	}

	// We set the infinity values to -1
	for(uint i = 0; i < V; i++) {
		for(uint j = 0; j < V; j++) {
			if(distances(i, j) > V) distances(i, j) = -1;
		}
	}

}

//Graph::~Graph()
//{
//    adjacency.~mat();
//    distances.~mat();
//}

uint Graph::vertex_number() {
	return V;
}

uint Graph::get_distance(uint v1, uint v2) {
	return distances(v1 - 1, v2 - 1);
}

bool Graph::is_an_edge(uint v1, uint v2) {
	return adjacency(v1 - 1, v2 - 1) == 1;
}

}
