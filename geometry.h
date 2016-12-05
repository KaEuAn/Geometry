#include <cstdio>
#include <vector>
#include <string>
#include <cmath>
#define PI 3.1415926535897932384626433832795


bool approximatelyEqual (const double& a, const double& b) {
    double eps = 2e-10;
    double difference = a - b;
    if (difference < 0)
        difference *= -1;
    return difference < eps;
}

struct Point {
    double x, y;

    Point(double a, double b) : x(a), y(b) {};

    Point() : x(0), y(0) {};
    Point(const Point& a) : x(a.x), y(a.y) {}

    bool operator ==(const Point& a) const {
        return (approximatelyEqual(x, a.x) && approximatelyEqual(y, a.y));
    }

    bool operator != (const Point& a) {
        return ! (*this == a);
    }

    Point operator - (const Point& a) const {
        Point c(*this);
        c.x -= a.x;
        c.y -= a.y;
    }
    Point operator + (const Point& a) const {
        Point c(*this);
        c.x += a.x;
        c.y += a.y;
    }
    Point operator / (double a) const {
        Point c(*this);
        c.x /= a;
        c.y /= a;
        return c;
    }

    double angle() const {
        return y / sqrt(x * x + y * y);
    }
    double length() const {
        return lengthToPoint({0, 0});
    }

    double lengthToPoint(const Point& b) const {
        return sqrt((x - b.x) * (x - b.x) + (y - b.y) * (y - b.y));
    }

    double rotate(Point center, double angle) {
        Point differ = *this - center;
        double old_angle = std::asin(this->angle());
        angle *= PI/180;
        old_angle += angle;
        x = std::sin(old_angle) * differ.length();
        y = std::cos(old_angle) * differ.length();
    }

    double scalarProduct(const Point& a) const {
        return a.length() * length() * std::cos(a.angle() - angle());
    }
    double pseudoScalarProduct(const Point& a) const {
        return a.length() * length() * std::sin(a.angle() - angle());
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

    explicit Line(double a, double b, double c) : A(a), B(b), C(c) {}

    bool isParallel (const Line& a) const {
        return approximatelyEqual(A * a.B, a.A * B);
    }

    bool operator == (const Line& a) const {
        return (isParallel(a) && approximatelyEqual(A * a.C, C * a.A));
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

    const double operator [] (char m) const {
        if (m == 'a')
            return A;
        if (m == 'b')
            return B;
        if (m == 'c')
            return C;
    }

    double getY (double x) {
        return -(C + A * x) / B;
    }

};

class Shape {
    /*virtual reflex(Point center);
    virtual reflex(Line axis);
    virtual scale(Point center, double coefficient);
    virtual rotate(Point center, double angle);

    virtual double perimeter() const;
    virtual double area() const;
    virtual operator ==(const Shape& another) const;
    virtual isCongruentTo(const Shape& another) const;
    virtual isSimilarTo(const Shape& another) const;
    virtual containsPoint(Point point) const;*/
};


class Polygon: public Shape {

protected:
    std::vector<Point> vertices;

    void pushEl(Point a) {
        vertices.push_back(a);
    }
    template<class... Args>
    void pushEl(Point a, Args... args) {
        vertices.push_back(a);
        pushEl(args...);
    }

public:
    //constructors
    template<class... Args>
    Polygon(Args... args) {
        pushEl(args...);
    }
    Polygon (std::vector<Point>& vec): vertices(vec){}
    //

    virtual ~Polygon() {};

    const std::vector<Point> getVertices () const {
        return vertices;
    }
    uint32_t verticesCount() const {
        return vertices.size();
    }
    bool isConvex() const {
        if (verticesCount() == 3)
            return true;
        uint32_t vc = verticesCount();
        for (uint32_t i = 0; i < vc; ++i) {
            Point i1 = vertices[i];
            Point i2 = vertices[(i + 1) % vc];
            Point i3 = vertices[(i + 2) % vc];
            Point i4 = vertices[(i + 3) % vc];
            double psc1 = (i2-i1).pseudoScalarProduct(i3-i2);
            double psc2 = (i3-i2).pseudoScalarProduct(i4-i3);
            if ((psc1 < 0 && psc2 > 0) || (psc1 > 0 && psc2 < 0))
                return false;
        }
        return true;
    }

    double area() const {
        double answer = 0;
        for (uint32_t i = 0; i < vertices.size() - 2; ++i) {
            Point a = vertices[i];
            Point b = vertices[i + 1];
            Point c = vertices[i + 2];
            answer += (b-a).pseudoScalarProduct(c-b);
        }
        return std::abs(answer);
    }

    double perimeter() const {
        double answer = 0;
        for (uint32_t i = 0; i < vertices.size() - 1; ++i) {
            Point a = vertices[i];
            Point b = vertices[i + 1];
            answer += b.lengthToPoint(a);
        }
        return answer;
    }

};

class Rectangle: public Polygon {

public:
    Rectangle (Point a, Point b, double ratio) {
        double diagonal = a.lengthToPoint(b);
        double minSide = diagonal / sqrt(1 + ratio * ratio);
        double maxSide = minSide * ratio;
        if (minSide > maxSide)
            std::swap(minSide, maxSide);

    }

    Point center() {
        Point answer((vertices[0].x + vertices[2].x) / 2, (vertices[0].y + vertices[2].y) / 2);
        return answer;
    }

};

class Triangle: public Polygon {

public:
    Triangle(Point a, Point b, Point c): vertices(a,b,c){}
    Point centroid() const {
        Point sideCenter((vertices[0].x + vertices[1].x) / 2, (vertices[0].y + vertices[1].y) /2);
        Point answer((2 * sideCenter.x + vertices[2].x) / 3, (2 * sideCenter.y + vertices[2].y) / 3);
        return answer;
    }

    Circle ninePointsCircle() {}
    Line EulerLine() {}
    Point orthocenter() {}
    Circle inscribedCircle() {}
    Circle circumscribedCircle() {}

};


class Ellipse: public Shape {

protected:
    std::pair <Point, Point> _focuses;
    double semiMajorAxis;


public:
    Ellipse (Point a, Point b, double c): _focuses(a, b), semiMajorAxis(c / 2) {}
    Ellipse () {}
    virtual ~Ellipse() {};

    const std::pair<Point, Point> focuses() const {
        return _focuses;
    };

    Point center () const {
        Point x = (_focuses.first + _focuses.second) / 2;
        return x;
    }

    std::pair<Line, Line> directrixes() const {
        double e = eccentricity();
        Point cent = center();
        Line central(_focuses.first, _focuses.second);
        std::pair<Line, Line> c;
        /**/
        return c;
    };

    double FocalLength() const {
        return (_focuses.first).lengthToPoint(_focuses.second);
    }
    double semiMinorAxis() const {
        return sqrt(semiMajorAxis * semiMajorAxis - FocalLength() * FocalLength());
    }
    double eccentricity() const {
        return FocalLength() / semiMajorAxis / 2;
    }
    virtual double area() const {
        return semiMajorAxis * PI * semiMinorAxis();
    }
    virtual double perimeter() const {
        double a = semiMajorAxis;
        double b = semiMinorAxis();
        return PI * (3 * (a + b) - sqrt( (3 * a + b) * (3 * b + a) ) );
    }

};

class Circle: public Ellipse {

public:
    Circle (Point a, double b): Ellipse(a, a, b) {}

    double radius() {
        return semiMajorAxis;
    }
    double area() const {
        return radius() * radius() * PI;
    }
    double perimeter() const {
        return radius() * 2 * PI;
    }


};

