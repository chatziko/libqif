
extern std::map<std::string,latlon> locations;

// Point

template<typename eT>
std::ostream& operator<<(std::ostream& os, const Point<eT>& p) {
	return os << '(' << p.x << ',' << p.y << ')';
}


// LatLon

template<typename eT>
eT rad_of_deg(eT ang) { return ang * pi<eT>() / 180; }

template<typename eT>
eT deg_of_rad(eT ang) { return ang * 180 / pi<eT>(); }


template<typename eT>
std::ostream& operator<<(std::ostream& os, const LatLon<eT>& p) {
	return os << '(' << p.lat << ',' << p.lon << ')';
}

// Returns the location obtained after moving distance (meters) in the direction of angle (rads)
// from the current location. 
//
template<typename eT>
LatLon<eT> LatLon<eT>::add_vector(eT distance, eT angle) {

	eT earth_radius(6378137); //const, in meters
	eT pi = qif::pi<eT>();

	eT ang_distance = distance / earth_radius;
	eT lat1 = rad_of_deg(lat);
	eT lon1 = rad_of_deg(lon);

	eT lat2 =	std::asin(
					std::sin(lat1) * std::cos(ang_distance) + 
					std::cos(lat1) * std::sin(ang_distance) * std::cos(angle)
			  	);
	eT lon2 =	lon1 +
			   	std::atan2(
					std::sin(angle) * std::sin(ang_distance) * std::cos(lat1), 
					std::cos(ang_distance) - std::sin(lat1) * std::sin(lat2)
				);
	lon2 = std::fmod(lon2 + 3 * pi, 2 * pi) - pi;		// normalise to -180..+180
	return LatLon<eT>(deg_of_rad(lat2), deg_of_rad(lon2));
}



// TODO: include other geo stuff here
namespace geo {


// returns a function that convers grid cell-id (uint) to Point
//
template<typename eT = eT_def>
std::function<Point<eT>(uint)>
cell_to_point(uint grid_width, eT cell_size = 1.0, Point<eT> corner = Point<eT>(eT(0), eT(0))) {
	return [=](uint i) -> Point<eT> {
		return Point<eT>(
			(i%grid_width) * cell_size + corner.x,
			(i/grid_width) * cell_size + corner.y
		);
	};
}

// returns a function that convers Points to grid cell-id.
// 0 is the bottom-left cell ((0,0), overridable), indexes increase left-to-right and bottom-to-top
//
template<typename eT = eT_def>
std::function<uint(const Point<eT>&)>
point_to_cell(uint grid_width, eT cell_size = eT(1), Point<eT> corner = Point<eT>(eT(0), eT(0))) {
	return [=](Point<eT> p) -> uint {
		// make sure we're within the grid
		assert(
			p.x >= corner.x &&
			p.y >= corner.y &&
			p.x < corner.x + grid_width * cell_size
		);

		return 
			floor((p.x - corner.x) / cell_size) +
			floor((p.y - corner.y) / cell_size) * grid_width;
	};
}

} // namespace geo
