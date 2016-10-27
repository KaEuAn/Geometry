#include <cstdio>
#include <vector>
#include <string>
#include <cmath>

bool approximatelyEqual (const double& a, const double& b) {
    double eps = 2e-10;
    double difference = a - b;
    if (difference < 0)
        difference *= -1;
    return difference < eps;
}

struct Point {
    double x, y;

    Point (double a, double b) : x(a), y(b) {};

    Point () : x(0), y(0) {};

    bool operator == (const Point& a) const {
        return (approximatelyEqual(x, a.x) && approximatelyEqual(y, a.y));
    }

    bool operator != (const Point& a) {
        return ! (*this == a);
    }

    double lengthToPoint (Point b) {
        return sqrt((x - b.x) * (x - b.x) + (y - b.y) * (y - b.y));
    }
};

class Line {

    // A*x + B*y + C = 0
    double A;
    double B;
    double C;

public:

    Line (Point one, Point two) : A(one.y - two.y), B(two.x - one.x), C(one.x * two.y - two.x * one.y) {};

    Line (double k, double b) :  A(-k), B(1), C(-b) {};

    Line (Point p, double k) :  A(-k), B(1), C(k * p.x - p.y) {};

    Line () : A(1), B(1), C(0) {};

    bool operator == (const Line& a) const {
        return (approximatelyEqual(A * a.B, a.A * B) && approximatelyEqual(A * a.C, C * a.A));
    }
    bool operator != (const Line& a) const {
        return ! (a == *this);
    }
    
    Line makePerpendicularLine (Point incoming) {
        Line outcoming;
        outcoming.A = - B * B;
        outcoming.B = A * B;
        outcoming.C = - outcoming.A * incoming.x - outcoming.B * incoming.y;
        return outcoming;
    }

};

class Shape {
    /*virtual reflex(Point center);
    virtual reflex(Line axis);
    virtual scale(Point center, double coefficient);
    virtual rotate(Point center, double angle);

    virtual double perimeter() const;
    virtual double area() const;
    virtual operator==(const Shape& another) const;
    virtual isCongruentTo(const Shape& another) const;
    virtual isSimilarTo(const Shape& another) const;
    virtual containsPoint(Point point) const;*/
};


class Polygon: public Shape {
    unsigned int quantityVertices;
    std::vector <Point> vertices;



public:
    template <class... Args>
    Polygon (Args... args) {

    }


    const std::vector<Point> getVertices () const {
        return vertices;
    }

    unsigned int verticesCount() const {
        return quantityVertices;
    }

    bool isConvex() const {
        bool ans;
        
        return ans;
    }
};


class Ellipse: public Shape {
    std::pair <Point, Point> _focuses;
    double distance;


public:
    Ellipse (Point a, Point b, double c): _focuses(a, b), distance(c) {}
    Ellipse () {}

    const std::pair<Point, Point> focuses() const {
        return _focuses;
    };

    Point center () const {
        Point x(((_focuses.first).x + (_focuses.second).x) / 2, ((_focuses.first).y + (_focuses.second).y) / 2  );
        return x;
    }

    std::pair<Line, Line> directrixes() const {
        double e = eccentricity();
        Point cent = center();
        Line central(_focuses.first, _focuses.second);
        std::pair<Line, Line> c();
        return c;
    };

    double eccentricity() const {
        double c = (_focuses.first).lengthToPoint(_focuses.second);
        return c / distance;
    }
};

class Circle: public Ellipse {
    Point Center;
    double radius


public:
    Circle (Point a, double b) {}
};

