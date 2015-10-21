#include <iostream>
#include <cstdlib>
#include <cassert>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <cmath>
using namespace std;

class UnifGen
{
private:


public:
    static int seed[2];
    static void setSeed(int x0, int x1)
    {
        seed[0] = x0;
        seed[1] = x1;
    }

    static float getUnifNum()
    {

        long long  j,k,a=1999,b=4444,c=2147483647;

        j=a*seed[0]+b*seed[1];
        k=j-c*(int)(j/c);
        seed[1]=seed[0];
        seed[0]=(int)k;
        float res =  (float)(k)/2147483647.f;

        assert (res > 0 && res < 1);
        return res;
    }

};

int UnifGen::seed[2];

// N is the dimension of the space
// K is the number of points

enum DistanceType
{
    Euclidean
};


enum RandomGenType
{

};


const float REL_EPSILON = 0.00001;

bool AlmostEqualRelative(float A, float B, float maxRelativeError)
{
    if (A == B)
        return true;
    float relativeError = fabs((A - B) / B);
    if (relativeError <= maxRelativeError)
        return true;
    return false;
}

class Point
{
private:
    int N;

public:
    vector <float> components;
    int label;

    friend ostream& operator<<(ostream& os, const Point& pt);

    bool operator== (const Point &cP2) const
    {
        for (int i = 0; i < N; i++) {
            if (!AlmostEqualRelative(components[i], cP2.components[i], REL_EPSILON))
                return false;
        }
        return true;
    }

    Point(int _N) : N(_N), label(-1) { }


    // XXX: Who deletes this?
    static const Point *mkRandomPoint(int N)
    {
        Point *newPoint = new Point(N);
        for (int i = 0; i < N; i++) {
            float value = UnifGen::getUnifNum();
            newPoint->components.push_back(value);  // XXX: Should be arbitrary random function here (instead of unif. only)
        }

        newPoint->label = (UnifGen::getUnifNum() < .5);

        return newPoint;
    }

    float getDistance(const Point &p, DistanceType distType) const
    {
        switch(distType) {
            case DistanceType::Euclidean:
                return getEuclDistance(p);
                break;

        };
    }

    float getEuclDistance(const Point &p) const
    {
        float sumOfSquareDiffs = 0.0;
        for (int i = 0; i < N; i++)
        {
            sumOfSquareDiffs += (components[i]-p.components[i])*(components[i]-p.components[i]);
        }
        return sqrt(sumOfSquareDiffs);
    }

};

ostream& operator<<(ostream& os, const Point& pt)
{
    os << "Point(";
    for (int i = 0; i < pt.N; i++) {
        if (i > 0)
            os << ", ";
        os << pt.components[i] ;
    }
    os << ") - Class: " << pt.label << ".";
    return os;
}

typedef vector<const Point *> PointSeq;


// We compare (distance, point) pairs by distance
bool compare(const pair<float,const Point *>&i, const pair<float,const Point *>&j)
{
    return i.first > j.first;
}

class RandomPointSet
{
    private:
    PointSeq points;
    int N, K;

    public:

    friend ostream& operator<<(ostream& os, const RandomPointSet& pt);

    RandomPointSet(int _K, int _N)
    {
        K = _K;
        N = _N;

        for (int i = 0; i < K; i++)
        {
            auto pt = Point::mkRandomPoint(N);
            points.push_back(pt);
        }

        for (int i = 0; i < K; i++) {
            cout << "Closest point to " << *points[i] << ":\n";
            auto closestPoint = getClosestPoint(*points[i], DistanceType::Euclidean);
            cout << closestPoint << endl;

            cout << "Farthest point to " << *points[i] << ":\n";
            auto farthestPoint = getFarthestPoint(*points[i], DistanceType::Euclidean);
            cout << farthestPoint << endl;

        }
    }

    const Point &getClosestPoint(const Point &p, DistanceType distType) const
    {
        // For all points different from p, save point and distance in a pair and then grab the min
        vector<pair<float, const Point *>> distancePoints;

        // XXX: How do we check that it's not the same point? (we check for equality)
        for (const Point *p1 : points) {
            if (*p1 == p) {
                // Do nothing
            } else {
                distancePoints.push_back(make_pair(p.getDistance(*p1, distType), p1));
            }
        }

        auto closestPoint = *min_element(distancePoints.begin(),distancePoints.end());
        return *closestPoint.second;
    }

    const Point &getFarthestPoint(const Point &p, DistanceType distType) const
    {
        // For all points different from p, save point and distance in a pair and then grab the min
        vector<pair<float, const Point *>> distancePoints;

        // XXX: How do we check that it's not the same point? (we check for equality)
        for (const Point *p1 : points) {
            if (*p1 == p) {
                // Do nothing
            } else {
                distancePoints.push_back(make_pair(p.getDistance(*p1, distType), p1));
            }
        }

        auto farthestPoint = *max_element(distancePoints.begin(),distancePoints.end());
        return *farthestPoint.second;
    }

    const PointSeq &getPoints()
    {
        return points;
    }

};

ostream& operator<<(ostream& os, const RandomPointSet& ptSet)
{
    os << "PointSeq:" << endl;
    for (int i = 0; i < ptSet.K; i++) {
        os << "\tP" << i << ": " << *(ptSet.points[i]) << endl;
    }
    return os;
}

class Exp1Simulation
{
private:
    RandomPointSet *pointSet = nullptr;
public:
    // NOTE: At some point you may want to variate on K and N and distance function
    Exp1Simulation(int seed0, int seed1, int K, int N)
    {
        UnifGen::setSeed(seed0,seed1); // XXX
        pointSet = new RandomPointSet(K, N);

        cout << *pointSet;

    }

    ~Exp1Simulation()
    {
        delete pointSet;
    }

    float getAverageRatioDist(DistanceType distType)
    {
        auto points = pointSet->getPoints();
        for (const Point *pPoint : points)
        {
            auto closestPoint = pointSet->getClosestPoint(*pPoint, distType);
            auto farthestPoint = pointSet->getFarthestPoint(*pPoint, distType);

            float distClosest = pPoint->getDistance(closestPoint, distType);
            float distFarthest = pPoint->getDistance(farthestPoint, distType);

            return distClosest / distFarthest;
        }
    }


};

int main(int argc, char **argv)
{
    int N, K;
    if (argc < 3) { // Use default
        N = 1;
        K = 3;
    } else {
        cerr << "Unimplemented\n!";
    }

    // Set output params
    cout << fixed;
    cout << setprecision(4);

    // Set seed


    Exp1Simulation exp1(1, 1, K, N);


    return 0;
}
