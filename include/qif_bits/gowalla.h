using ::std::string;

namespace gowalla {

struct Entry {
	uint userId;
	latlon location;
	uint locationId;
};

std::vector<Entry> read_dataset(string filename);
std::vector<uint> to_grid(const std::vector<Entry>& dataset, latlon center, uint width, uint height, double cell_size);
prob to_grid_prior(const std::vector<Entry>& dataset, latlon center, uint width, uint height, double cell_size);

// TODO: extract the main projection functionality to geo::project_dummy
std::vector<point> project_dummy(const std::vector<Entry>& dataset, latlon center, double width, double height);

} // namespace gowalla
