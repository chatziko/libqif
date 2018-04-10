#include "qif"

using ::std::string;

namespace qif {

std::map<std::string,latlon> locations = {
	{ "paris",			latlon(48.858072, 2.348050) },
	{ "san_francisco",	latlon(37.755351, -122.440288) }
};

namespace gowalla {

std::vector<Entry> read_dataset(string filename) {
	std::vector<Entry> list;
	string line, date;

	std::ifstream file;
	file.open(filename);
	if(!file)
		throw std::runtime_error("cannot open " + filename);

    while(getline(file, line) ) {
		Entry entry;
		std::istringstream line_s(line);
		line_s >> entry.userId >> date >> entry.location.lat >> entry.location.lon >> entry.locationId;
 		list.push_back(entry);
    }
	file.close();

	return list;
}

std::vector<uint> to_grid(const std::vector<Entry>& dataset, latlon center, uint width, uint height, double cell_size) {

	double lon_d = std::abs(center.lon - center.add_vector(cell_size, pi<double>()/2).lon);
	double lat_d = std::abs(center.lat - center.add_vector(cell_size, 0).lat);

	double lon_min = center.lon - (lon_d * width) / 2;
	double lon_max = center.lon + (lon_d * width) / 2;
	double lat_min = center.lat - (lat_d * height) / 2;
	double lat_max = center.lat + (lat_d * height) / 2;

	std::vector<uint> points;
	for(auto e : dataset) {
		if(!(e.location.lat >= lat_min && e.location.lat < lat_max && e.location.lon >= lon_min && e.location.lon < lon_max))
			continue;

		uint x = std::floor((e.location.lon - lon_min) / lon_d);
		uint y = std::floor((e.location.lat - lat_min) / lat_d);

		points.push_back(y * width + x);
	}

	return points;
}

} // namespace gowalla
} // namespace qif
