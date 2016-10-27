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

std::pair<short, std::pair<double, double> > solveQuadratic (double a, double b, double c) {
    std::pair<short, std::pair<double, double> > answer;
    double discriminant = b * b - 4 * a * c;
    if (approximatelyEqual(discriminant, 0)) {
        answer.first = 1;
        answer.second.first = -b / 2 / a;
    } else if (discriminant < 0) {
        answer.first = 0;
    } else {
        discriminant = sqrt(discriminant);
        answer.first = 2;
        answer.second.first = (discriminant - b) / 2 / a;
        answer.second.second = (-discriminant - b) / 2 / a;
    }
    return answer;
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
    virtual operator==(const Shape& another) const;
    virtual isCongruentTo(const Shape& another) const;
    virtual isSimilarTo(const Shape& another) const;
    virtual containsPoint(Point point) const;*/
};


class Polygon: public Shape {

protected:
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

class Rectangle: private Polygon {


public:
    Rectangle (Point a, Point b, double ratio) {
        double diagonal = a.lengthToPoint(b);
        double minSide = diagonal / sqrt(1 + ratio * ratio);
        double maxSide = minSide * ratio;
        if (minSide > maxSide)
            std::swap(minSide, maxSide);
        Circle one(a, minSide);
        Circle two(b, maxSide);

    }

    Point center() {

        return ;
    }
};


class Ellipse: public Shape {

protected:
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
        std::pair<Line, Line> c;
        /**/
        return c;
    };

    double eccentricity() const {
        double c = (_focuses.first).lengthToPoint(_focuses.second);
        return c / distance;
    }
};

class Circle: public Ellipse {

public:
    Circle (Point a, double b): Ellipse(a, a, b) {}

    double radius() {
        return distance;
    }

    const std::pair<short, std::pair<Point, Point> > intersectionWithLine (Line line) {
        short count = 0;
        // quadratic equation
        if (! approximatelyEqual(line['b'], 0)) {
            double quadratic_a = line['a'] * line['a'] + line['b'] * line['b'];
            double quadratic_b = 2 * line['a'] * line['c'] - 2 * line['b'] * line['b'] * _focuses.first.x +
                                 2 * line['b'] * _focuses.first.y;
            double quadratic_c = line['b'] * line['b'] *
                                 (_focuses.first.x * _focuses.first.x + _focuses.first.y * _focuses.first.y -
                                  distance * distance) + 2 * line['b'] * +line['c'] * line['c'];
            const std::pair<short, std::pair<double, double>> solution = solveQuadratic (quadratic_a, quadratic_b, quadratic_c);
            std::pair<short, std::pair<Point, Point> > answer;
            answer.first = solution.first;
            answer.second.first.x = solution.second.first;
            answer.second.first.y = line.getY(answer.second.first.x);
            answer.second.second.x = solution.second.second;
            answer.second.second.y = line.getY(answer.second.second.x);
            return answer;
        } else {
            double answer_x = -line['c'] / line['a'];
            double quadratic_a = 1;
            double quadratic_b = _focuses.first.y * 2;
            double quadratic_c = (answer_x - _focuses.first.x) * (answer_x - _focuses.first.x) + _focuses.first.y * _focuses.first.y - distance * distance;
            const std::pair<short, std::pair<double, double>> solution = solveQuadratic (quadratic_a, quadratic_b, quadratic_c);
            std::pair<short, std::pair<Point, Point> > answer;
            answer.first = solution.first;
            answer.second.first.x = answer_x;
            answer.second.first.y = solution.second.first;
            answer.second.second.x = answer_x;
            answer.second.second.y = solution.second.second;
            return answer;
        }
    };

    const std::pair<short, std::pair<Point, Point> > intersectionWithCircle (Circle circle) {
        Line line(2 * (_focuses.first.x - circle._focuses.first.x), 2 * (_focuses.first.y - circle._focuses.first.y),
                  circle._focuses.first.x *  circle._focuses.first.x - _focuses.first.x * _focuses.first.x + circle._focuses.first.y *  circle._focuses.first.y - _focuses.first.y * _focuses.first.y + distance * distance - circle.distance * circle.distance);
        return intersectionWithLine(line);
    }
};

