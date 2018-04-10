using ::std::string;

namespace gowalla {

struct Entry {
	uint userId;
	latlon location;
	uint locationId;
};

std::vector<Entry> read_dataset(string filename);
std::vector<uint> to_grid(const std::vector<Entry>& dataset, latlon center, uint width, uint height, double cell_size);

} // namespace gowalla
