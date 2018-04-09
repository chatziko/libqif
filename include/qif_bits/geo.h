

std::map<std::string,latlon> locations = {
	{ "paris",			latlon(48.858072, 2.348050) },
	{ "san_francisco",	latlon(37.755351, -122.440288) }
};


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
