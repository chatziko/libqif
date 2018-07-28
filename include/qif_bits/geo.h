
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
// 0 is the bottom-left cell, indexes increase left-to-right and bottom-to-top
// corner is the _CENTER_ of cell 0 (default 0,0)
//
template<typename eT = eT_def>
std::function<uint(const Point<eT>&)>
point_to_cell(uint grid_width, eT cell_size = eT(1), Point<eT> corner = Point<eT>(eT(0), eT(0))) {
	// corner is the center of cell 0, make it the bottom-left corder of cell 0
	corner.x -= cell_size/2;
	corner.y -= cell_size/2;

	return [=](Point<eT> p) -> uint {
		// make sure we're within the grid
		if(!(
			p.x >= corner.x &&
			p.y >= corner.y &&
			p.x < corner.x + grid_width * cell_size
		))
			throw std::runtime_error("out of grid area");

		return 
			floor((p.x - corner.x) / cell_size) +
			floor((p.y - corner.y) / cell_size) * grid_width;
	};
}


// iterate points on infinite grid of given cell_size, starting at (0,0)
// iterate per ring r, each ring is defined by max{|xd|,|yd|} == r
//
template<typename eT = eT_def>
class GridWalk {
  private:
	struct generator : public std::iterator<std::forward_iterator_tag, Point<eT>> {
	  private:
		eT cell_size;
		int r = 0, xd = 0, yd = 0;
		Point<eT> cur;

	  public:
		generator(eT _cell_size, uint _r)
			: cell_size(_cell_size), r(_r), cur(0,0) {}

		Point<eT>& operator*() {
			return cur;
		}

		bool operator != (const generator& other) const {
			return !(r == other.r && xd == other.xd && yd == other.yd);
		}

		generator& operator++() {
			// iteration steps:
			//  A. xd =  r, yd = -r .. r
			//  B. xd = -r, yd = -r .. r
			//  C. xd = -r+1 .. r-1, yd = r
			//  D. xd = -r+1 .. r-1, yd = -r
			//
			if(r == -1) { // r == -1 means past the end, cannot advance

			} if(r == 0) { // edge case, the first ring has only one step
				r++;
				xd = r;
				yd = -r;
			} else if(xd == r) {	// step A
				if(++yd > r) {		// goto step B
					xd = -r;
					yd = -r;
				}
			} else if(xd == -r) {	// step B
				if(++yd > r) {		// goto step C
					xd = -r+1;
					yd = r;
				}
			} else if(yd == r) {	// step C
				if(++xd == r) {		// goto step D
					xd = -r+1;
					yd = -r;
				}
			} else {				// step D
				if(++xd == r) {		// inc r, goto step A
					r++;
					xd = r;
					yd = -r;
				}
			}

			cur = Point<eT>(
				xd * cell_size,
				yd * cell_size
			);
			return *this;
		}
	};

  private:
	eT cell_size;

  public:
	GridWalk(eT _cell_size = eT(1))
		: cell_size(_cell_size) {}

	generator begin() const {
		return generator(cell_size, 0);
	}
	generator end() const {
		return generator(cell_size, -1);
	}
};


} // namespace geo
