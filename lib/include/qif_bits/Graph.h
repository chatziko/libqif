class Graph {
	public:
		Graph(uint vertex_num, std::string& edges);

		Graph(uint vertex_num, std::vector< std::pair<int, int> >& edges);

//		~Graph();

		uint vertex_number();
		bool is_an_edge(uint v1, uint v2);
		uint get_distance(uint v1, uint v2);

	protected:
		uint V;
		mat adjacency;
		mat distances;
};
