#define NONCOORD -666666666
#define NONDIST 666666666

#define UNKWN_SHROUD -1
#define FIBER_THICKNESS 0.06

Double_t EPS = 0.00001;
Double_t INF = 1000000;

// Struct 3d point
struct Point{
    Double_t x;
    Double_t y;
    Double_t z;
    Int_t shroud;           // Optional
    Double_t distance;      // Optional
    Point(Double_t xx=NONCOORD, Double_t yy=NONCOORD, Double_t zz=NONCOORD,
          Int_t shrd=UNKWN_SHROUD, Double_t dist=NONDIST){
        x = xx;
        y = yy;
        z = zz;
        shroud = shrd;
        distance = dist;
    }
    void reset(){
        x = NONCOORD;
        y = NONCOORD;
        z = NONCOORD;
        shroud = UNKWN_SHROUD;
        distance = NONDIST;
    }
    bool isDefined(){
        return (x!=NONCOORD || y!=NONCOORD || z!=NONCOORD);
    }
    bool checkSameDirection(Double_t x2, Double_t y2, struct Point p){
        return (x-x2)*(x-p.x)>=0 && (y-y2)*(y-p.y)>=0;
    }
    bool operator < (const Point& point){
        return distance < point.distance;
    }
    bool operator == (const Point& point){
        return abs(x - point.x)<=EPS || abs(y - point.y)<=EPS || abs(z - point.z)<=EPS;
    }
};

std::ostream& operator<<(std::ostream& os, const Point& point){
    os << "X: " << point.x << ", ";
    os << "Y: " << point.y << ", ";
    os << "Z: " << point.z << " | ";
    os << "Shroud: " << point.shroud << ", ";
    os << "Dist_prod: " << point.distance;
    return os;
}

Double_t getPointDistance(struct Point a, struct Point b){
    return sqrt((a.x-b.x)*(a.x-b.x) + (a.y-b.y)*(a.y-b.y) + (a.z-b.z)*(a.z-b.z));
}

vector<Point> computeIntersection(Point prodPoint, Point dirPoint, Double_t r, Double_t maxDistance=INF) {
    vector<Point> hits;
    if (abs(sqrt(prodPoint.x * prodPoint.x + prodPoint.y * prodPoint.y) - r) > FIBER_THICKNESS/2) {    // check if point stay on circle boundary
        // Compute line parameters
        Double_t a = prodPoint.y - dirPoint.y;
        Double_t b = dirPoint.x - prodPoint.x;
        Double_t c = prodPoint.x * dirPoint.y - dirPoint.x * prodPoint.y;
        // Compute the point in line closest to the origin
        Point hit = Point();
        hit.x = -a * c / (a * a + b * b), hit.y = -b * c / (a * a + b * b);
        // Based on dist from origin and radius, we know multiplicity
        if (c * c > r * r * (a * a + b * b) + EPS) {
            hit.reset();
        } else if (abs(c * c - r * r * (a * a + b * b)) >= EPS) {
            double d = r * r - c * c / (a * a + b * b);
            double mult = sqrt(d / (a * a + b * b));
            double ax, ay, bx, by;
            ax = hit.x + b * mult;
            bx = hit.x - b * mult;
            ay = hit.y - a * mult;
            by = hit.y + a * mult;
            // Consider only interesection within maxDistance
            double distance_a = sqrt((ax - prodPoint.x) * (ax - prodPoint.x) + (ay - prodPoint.y) * (ay - prodPoint.y));
            double distance_b = sqrt((bx - prodPoint.x) * (bx - prodPoint.x) + (by - prodPoint.y) * (by - prodPoint.y));
            // assert: a valid intersection has the same direction(both x, y)
            // if 2 intersections can be 1) only 1 in the same direction, 2) both of them (line from outside)
            Point hit = Point(ax, ay);
            if (prodPoint.checkSameDirection(dirPoint.x, dirPoint.y, hit) && (distance_a <= maxDistance))
                hits.push_back(hit);
            hit = Point(bx, by);
            if (prodPoint.checkSameDirection(dirPoint.x, dirPoint.y, hit) && (distance_b <= maxDistance))
                hits.push_back(hit);
        } else { // else only 1 (xf, yf)
            double distance = sqrt((hit.x - prodPoint.x) * (hit.x - prodPoint.x) + (hit.y - prodPoint.y) * (hit.y - prodPoint.y));
            if (prodPoint.checkSameDirection(dirPoint.x, dirPoint.y, hit) && (distance <= maxDistance))     // the only intersection if farther than limit
                hits.push_back(hit);
        }
    } else {      // point on the circle boundary
        hits.push_back(prodPoint);
    }
    for (auto &hit: hits) {
        if (hit.isDefined()) {
            // If there is an interesection, look at where it occurs in Z
            // Use parametrization of line: line = P + t D, where P is a vector (point), D direction vector
            Double_t t = (hit.x - prodPoint.x) / (dirPoint.x - prodPoint.x);     // From param: t = (x-x1)/dx
            hit.z = prodPoint.z + t * (dirPoint.z - prodPoint.z);              // z = z1 + t * dz, where z1 point, dz direction vector
            hit.distance = getPointDistance(prodPoint, hit);
        }
    }
    std::sort(hits.begin(), hits.end());    // sort ascending distance
    return hits;
}

Int_t getPointPlacement(Double_t x, Double_t y, Double_t inner_r, Double_t outer_r){
    Double_t dist0 = sqrt(x*x + y*y);
    if(dist0 <= inner_r)	// point in [0,inner_r]
        return -1;
    if(dist0 <= outer_r)	// point in (inner_r, outer_r]
        return 0;
    return 1;		// point in (outer_r, *)
}


Double_t getAngleBtw3Points(Double_t cx, Double_t cy, Double_t ax, Double_t ay, Double_t bx, Double_t by){
    // Compute the angle between point A, B w.r.t. C (center)
    Double_t angle_a = atan2(ay - cy, ax - cx);
    if(angle_a < 0)
        angle_a += 2 * TMath::Pi();
    Double_t angle_b = atan2(by - cy, bx - cx);
    if(angle_b < 0)
        angle_b += 2 * TMath::Pi();
    return angle_b - angle_a;
}