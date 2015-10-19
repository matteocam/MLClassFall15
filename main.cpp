#include <iostream>
#include <cstdlib>
#include <cassert>
#include <vector>
using namespace std;

class UnifGen
{
private:


public:
    static int seed[2];
    static float setSeed(int x0, int x1)
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

class Point
{
private:
    int N;

public:
    vector <float> components;

    friend ostream& operator<<(ostream& os, const Point& pt);

    Point(int _N) : N(_N) { }


    // XXX: Who deletes this?
    static const Point *mkRandomPoint(int N)
    {
        Point *newPoint = new Point(N);
        for (int i = 0; i < N; i++)
            newPoint->components.push_back(UnifGen::getUnifNum());  // XXX: Should be arbitrary random fun here
    }

    float getDistance(const Point &p) const
    {
        assert(0); // XXX
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
    os << ")";
    return os;
}

typedef vector<const Point *> PointSeq;

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
    }

    const Point &getClosestPoint(const Point &p) const
    {

    }

    const Point &getFarthestPoint(const Point &p) const
    {

    }

    const PointSeq &getPoints()
    {

    }

};

ostream& operator<<(ostream& os, const RandomPointSet& ptSet)
{
    os << "PointSeq:" << endl;
    for (int i = 0; i < ptSet.K; i++) {
        os << "P" << i << ": " << *(ptSet.points[i]) << endl;
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

        PointSeq points = pointSet->getPoints();

    }

    ~Exp1Simulation()
    {
        delete pointSet;
    }

    float getAverageRatioDist()
    {
        auto points = pointSet->getPoints();
        for (const Point *pPoint : points)
        {
            auto closestPoint = pointSet->getClosestPoint(*pPoint);
            auto farthestPoint = pointSet->getFarthestPoint(*pPoint);

            float distClosest = pPoint->getDistance(closestPoint);
            float distFarthest = pPoint->getDistance(farthestPoint);

            return distClosest / distFarthest;
        }
    }


};

int main(int argc, char **argv)
{
    int N, K;
    if (argc < 3) { // Use default
        N = 2;
        K = 2;
    } else {
        cerr << "Unimplemented\n!";
    }

    // Set seed


    Exp1Simulation exp1(1, 1, K, N);


    return 0;
}
