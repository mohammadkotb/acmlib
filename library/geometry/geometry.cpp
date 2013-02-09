/*
 * Computational Geometry library based on Coach Mohamed Mahmoud Abdel-Wahab's
 * library (fegla on Topcoder)
 */


#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <complex>
#include <vector>
#include <algorithm>
#include <stack>

using namespace std;
typedef complex<double> point;
#define sz(a) ((int)(a).size())
#define all(a) (a).begin(),(a).end()
#define EPS 1.0e-8 //take care of EPS value
#define INF 1.0e9
#define PI acos(-1)
#define X real()
#define Y imag()
#define vec(a,b) ((b)-(a))
#define polar(r,t) ((r)*exp(point(0,(t))))
#define length(a) hypot((a).X,(a).Y)
#define angle(a) atan2((a).Y,(a).X)
#define lengthSqr(a) dot(a,a)
#define dot(a,b) ((conj(a) * (b)).real())
#define cross(a,b) ((conj(a) * (b)).imag())
#define rotate(v,t) (polar((v),(t)))
#define rotateabout(v,t,a)  (rotate(vec(a,v),t) + (a))
#define reflect(p,m) ((conj((p)/(m))) * (m))
#define normalize(a) ((a)/length(a))
#define same(a,b) (lengthSqr(vec(a,b))<EPS)
#define mid(a,b) (((a)+(b))/point(2,0))
#define perp(a) (point(-(a).Y,(a).X))
#define dist(a,b) (length(vec(a,b)))
#define ccw(a,b,c) (cross(vec(a,b),vec(a,c)) > EPS)
#define collinear(a,b,c) (fabs(cross(vec(a,b),vec(a,c))) < EPS)
enum STATE { IN, OUT, BOUNDRY };

bool intersect(const point& a, const point& b, const point& p, const point& q, point& ret) {
	double d1 = cross(p - a, b - a);
	double d2 = cross(q - a, b - a);
	ret = (d1 * q - d2 * p) / (d1 - d2); // TODO handle special cases (segment intersections/end-points)
	return fabs(d1 - d2) > EPS;
}

bool pointOnLine(const point& a, const point& b, const point& p) {
	return fabs(cross(vec(a,b),vec(a,p))) < EPS;
}

bool pointOnRay(const point& a, const point& b, const point& p) {
	return pointOnLine(a, b, p) && dot(vec(a,b),vec(a,p)) > -EPS;
}

bool pointOnSegment(const point& a, const point& b, const point& p) {
	if (same(a,b)) return same(a,p);
	return pointOnRay(a,b,p) && pointOnRay(b,a,p);
}

double pointLineDist(const point& a, const point& b, const point& p) {
	// handle segment case with dot(vec(a,b),vec(a,p)) and dot(vec(b,a),vec(b,p)) checks
	if (same(a,b)) return dist(a,p);
	return fabs(cross(vec(a,b),vec(a,p)) / dist(a,b));
}

int circleLineIntersection(const point& p0, const point& p1, const point& cen, double rad, point& r1, point & r2) {
	if (same(p0,p1)) {
		if(fabs(lengthSqr(vec(p0,cen))-(rad*rad)) < EPS) {
			r1 = r2 = p0;
			return 1;
		}
		return 0;
	}
	double a, b, c, t1, t2;
	a = dot(p1-p0,p1-p0);
	b = 2 * dot(p1-p0,p0-cen);
	c = dot(p0-cen,p0-cen) - rad * rad;
	double det = b * b - 4 * a * c;
	int res;
	if (fabs(det) < EPS) det = 0, res = 1;
	else if (det < 0) res = 0;
	else res = 2;
	det = sqrt(det);
	t1 = (-b + det) / (2 * a);
	t2 = (-b - det) / (2 * a);
	r1 = p0 + t1 * (p1 - p0);
	r2 = p0 + t2 * (p1 - p0);
	return res;
}

int circleSegmentIntersection(const point& p0, const point& p1, const point& cen, double rad, point& r1, point & r2) {
	int res = circleLineIntersection(p0, p1, cen, rad, r1, r2);
	if (res == 2) {
		if (!pointOnSegment(p0, p1, r1)) {
			res--;
			r1 = r2;
		}
		if (!pointOnSegment(p0, p1, r2)) {
			res--;
		}
	} else if (res == 1) {
		if (!pointOnSegment(p0, p1, r1)) {
			return 0;
		}
	}
	return res;
}

double cosRule(long double a, long double b, long double c) {
   double res = (b * b + c * c - a * a) / (2 * b * c);
   if (res > 1) res = 1;
   if (res < -1) res = -1;
   return acos(res);
}

int circleCircleIntersection(const point &c1, const double&r1, const point &c2, const double&r2, point &res1, point &res2) {
	if (same(c1,c2) && fabs(r1 - r2) < EPS) {
			res1 = res2 = c1;
			return fabs(r1) < EPS ? 1 : INF;
	}
	double len = length(vec(c1,c2));
	if (fabs(len - (r1 + r2)) < EPS || fabs(fabs(r1 - r2) - len) < EPS) {
			point d, c; double r;
			if (r1 > r2)  d = vec(c1,c2), c = c1, r = r1;
			else  d = vec(c2,c1), c = c2, r = r2;
			res1 = res2 = normalize(d) * r + c;
			return 1;
	}
	if (len > r1 + r2 || len < fabs(r1 - r2))  return 0;
	long double a = cosRule(r2, r1, len);
	point c1c2 = normalize(vec(c1,c2)) * r1;
	res1 = rotate(c1c2,a) + c1;
	res2 = rotate(c1c2,-a) + c1;
	return 2;
}

void circle2(const point& p1, const point& p2, point& cen, double& r) {
	cen = mid(p1,p2);
	r = length(vec(p1,p2)) / 2;
}

bool circle3(const point& p1, const point& p2, const point& p3, point& cen, double& r) {
	point m1 = mid(p1,p2);
	point m2 = mid(p2,p3);
	point perp1 = perp(vec(p1,p2));
	point perp2 = perp(vec(p2,p3));
	bool res = intersect(m1, m1 + perp1, m2, m2 + perp2, cen);
	r = length(vec(cen,p1));
	return res;
}

STATE circlePoint(const point & cen, const long double & r, const point& p) {
	long double lensqr = lengthSqr(vec(cen,p));
	if (fabs(lensqr - r * r) < EPS) return BOUNDRY;
	if (lensqr < r * r) return IN;
	return OUT;
}

int tangentPoints(const point & cen, const double & r, const point& p, point &r1, point &r2) {
	STATE s = circlePoint(cen, r, p);
	if (s != OUT) {
		r1 = r2 = p;
		return s == BOUNDRY;
	}
	point cp = vec(cen,p);
	long double h = length(cp);
	long double a = acos(r / h);
	cp = normalize(cp) * r;
	r1 = rotate(cp,a) + cen;
	r2 = rotate(cp,-a) + cen;
	return 2;
}

void circleCircleTangencyPoints(const point& c1, const double& r1, const point& c2, const double& r2, vector<point>& t) {
	if (fabs(r1-r2) <= EPS) { // special case
		point unit = normalize(vec(c1,c2));
		point p1 = point(-unit.Y, unit.X) * r1, p2 = point(unit.Y, -unit.X) * r1;
		t.push_back(p1 + c1); t.push_back(p1 + c2);
		t.push_back(p2 + c1); t.push_back(p2 + c2);
		return;
	}
	point t1, t2;
	tangentPoints(c2, r2 - r1, c1, t1, t2);
	point p1 = normalize(vec(c2, t1)), p2 = normalize(vec(c2, t2));
	t.push_back(p1 * r1 + c1); t.push_back(p1 * r2 + c2);
	t.push_back(p2 * r1 + c1); t.push_back(p2 * r2 + c2);
}

double arcLength(const point& p1, const point& p2, const point& c, const double& r, bool ordered) {
	double cross = cross(vec(c, p1), vec(c, p2));
	double dot = dot(vec(c, p1), vec(c, p2));
	dot = dot / r / r; // normalizing to carry only the value of cos(theta)
	if (!ordered) { // order or p1 and p2 doesn't matter
		return acos(dot) * r;
	}

	if (fabs(cross) <= EPS) { // either 0 or PI
		if (dot < 0) // angle = PI
			return PI * r;
	} else if (cross > 0) { // convex angle
		return acos(dot) * r;
	} else if (cross < 0) { // concave angle
		return (2 * PI - acos(dot)) * r;
	}
	return 0.0;
}


// minimum enclosing circle
// init p array with the points and ps with the number of points
// cen and rad are result circle
// you must call random_shuffle(p,p+ps); before you call mec
#define MAXPOINTS 100000
point p[MAXPOINTS], r[3], cen;
int ps, rs;
double rad;
void mec() {
	if (rs == 3) {
		circle3(r[0], r[1], r[2], cen, rad);
		return;
	}
	if (rs == 2 && ps == 0) {
		circle2(r[0], r[1], cen, rad);
		return;
	}
	if (!ps) {
		cen = r[0];
		rad = 0;
		return;
	}
	ps--;
	mec();
	if (circlePoint(cen, rad, p[ps]) == OUT) {
		r[rs++] = p[ps];
		mec();
		rs--;
	}
	ps++;
}

// return the centroid point of the polygon
// The centroid is also known as the "centre of gravity" or the "center of mass". The position of the centroid
// assuming the polygon to be made of a material of uniform density.
point polyginCentroid(vector<point> &polygon) {
	point res(0, 0);
	double a = 0;
	for (int i = 0; i < (int) polygon.size(); i++) {
		int j = (i + 1) % polygon.size();
		res.X += (polygon[i].X + polygon[j].X) * (polygon[i].X * polygon[j].Y
						- polygon[j].X * polygon[i].Y);

		res.Y += (polygon[i].Y + polygon[j].Y) * (polygon[i].X * polygon[j].Y
						- polygon[j].X * polygon[i].Y);
		a += polygon[i].X * polygon[j].Y - polygon[i].Y * polygon[j].X;
	}
	a *= 0.5;
	res.X /= 6 * a;
	res.Y /= 6 * a;
	return res;
}

double polygonArea(vector<point>& pts) {
	double area = 0;
	for (int i = 0; i < sz(pts); i++)
		area += cross(pts[i], pts[(i + 1) % sz(pts)]);
	area = fabs(area) * 0.5;
	return area;
}

void polygonCut(const vector<point>& p, const point&a, const point&b, vector<point>& res) {
	res.clear();
	for (int i = 0; i < sz(p); i++) {
		int j = (i + 1) % sz(p);
		bool in1 = cross(vec(a,b),vec(a,p[i])) > EPS;
		bool in2 = cross(vec(a,b),vec(a,p[j])) > EPS;
		if (in1) res.push_back(p[i]);
		if (in1 ^ in2) {
			point r;
			intersect(a, b, p[i], p[j], r);
			res.push_back(r);
		}
	}
}

//assume that both are anti-clockwise
void convexPolygonIntersect(const vector<point>& p, const vector<point>& q, vector<point>& res) {
	res = q;
	for (int i = 0; i < sz(p); i++) {
		int j = (i + 1) % sz(p);
		vector<point> temp;
		polygonCut(res, p[i], p[j], temp);
		res = temp;
		if (res.empty()) return;
	}
}

void voronoi(const vector<point> &pnts, const vector<point>& rect, vector<vector<point> > &res) {
	res.clear();
	for (int i = 0; i < sz(pnts); i++) {
		res.push_back(rect);
		for (int j = 0; j < sz(pnts); j++) {
			if (j == i) continue;
			point p = perp(vec(pnts[i],pnts[j]));
			point m = mid(pnts[i],pnts[j]);
			vector<point> temp;
			polygonCut(res.back(), m, m + p, temp);
			res.back() = temp;
		}
	}
}

STATE pointInPolygon(const vector<point>& p, const point &pnt) {
	point p2 = pnt + point(1, 0);
	int cnt = 0;
	for (int i = 0; i < sz(p); i++) {
		int j = (i + 1) % sz(p);
		if (pointOnSegment(p[i], p[j], pnt)) return BOUNDRY;
		point r;
		intersect(pnt, p2, p[i], p[j], r);
		if (!pointOnRay(pnt, p2, r)) continue;
		if (same(r,p[i]) || same(r,p[j]))
			if (fabs(r.Y - min(p[i].Y, p[j].Y)) < EPS) continue;
		if (!pointOnSegment(p[i], p[j], r)) continue;
		cnt++;
	}
	return cnt & 1 ? IN : OUT;
}

point pivotCH;
bool angleCmp(const point& a, const point& b) {
	if (collinear(pivotCH, a, b)) {
		return dist(pivotCH, a) < dist(pivotCH, b);
	}
	double cross = cross(vec(pivotCH, a), vec(pivotCH, b));
	return cross > 0;
}

void convexHull(vector<point>& P, vector<point>& res) {
	int i, N = sz(P);
	if (N <= 3) {
		for (int i = 0; i < sz(P); i++) res.push_back(P[i]);
		return;
	}
	int P0 = 0;
	for (i = 1; i < N; i++)
		if (P[i].Y < P[P0].Y || (P[i].Y == P[P0].Y && P[i].X > P[P0].X))
			P0 = i;
	point temp = P[0]; P[0] = P[P0]; P[P0] = temp;
	pivotCH = P[0];
	sort(P.begin()+1, P.end(), angleCmp);
	point prev(0, 0), now(0, 0);
	stack<point> S; S.push(P[N-1]); S.push(P[0]);
	i = 1;
	while (i < N) {
		now = S.top(); S.pop();
		prev = S.top(); S.push(now);
		if (ccw(prev, now, P[i])) S.push(P[i++]);
		else S.pop();
	}
	while (!S.empty()) { res.push_back(S.top()); S.pop(); }
}


int main() {
	printf("Hello world!\n");
	point p1(-1, 2), p2(5, -4), o(0,0);
	cout << pointOnSegment(p1, p2, o) << endl;;
	cout << collinear(p1, p2, o) << endl;
	return 0;
}

