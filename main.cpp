#include <iostream>
#include <string>
#include <fstream>
#include <cstdlib>
#include <ctime>
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

string distanceType2str(DistanceType dt)
{
    switch(dt)
    {
        case DistanceType::Euclidean:
            return "Euclidean";
            break;
        default:
            assert(0);

    };
}


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
            default:
                assert(0); // You should never get here
                return -1;

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

    int getK() const { return K; }

    RandomPointSet(int _K, int _N)
    {
        K = _K;
        N = _N;

        for (int i = 0; i < K; i++)
        {
            auto pt = Point::mkRandomPoint(N);
            points.push_back(pt);
        }

    /*
        for (int i = 0; i < K; i++) {
            cout << "Closest point to " << *points[i] << ":\n";
            auto closestPoint = getClosestPoint(*points[i], DistanceType::Euclidean);
            cout << closestPoint << endl;

            cout << "Farthest point to " << *points[i] << ":\n";
            auto farthestPoint = getFarthestPoint(*points[i], DistanceType::Euclidean);
            cout << farthestPoint << endl;

        }
    */
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
    static int stdSeed0, stdSeed1;

    // NOTE: At some point you may want to variate on K and N and distance function
    Exp1Simulation(int seed0, int seed1, int K, int N)
    {
        UnifGen::setSeed(seed0,seed1); // XXX
        pointSet = new RandomPointSet(K, N);

        //cout << *pointSet;

    }

    ~Exp1Simulation()
    {
        delete pointSet;
    }

    float getAverageRatioDist(DistanceType distType)
    {
        auto points = pointSet->getPoints();
        float sumDists = 0.0;

        for (const Point *pPoint : points)
        {
            auto closestPoint = pointSet->getClosestPoint(*pPoint, distType);
            auto farthestPoint = pointSet->getFarthestPoint(*pPoint, distType);

            float distClosest = pPoint->getDistance(closestPoint, distType);
            float distFarthest = pPoint->getDistance(farthestPoint, distType);

            sumDists += (distClosest / distFarthest);
        }
        return sumDists / pointSet->getK();
    }

    void logResults(ofstream &outFile)
    {

    }

    void static runFullExperiment()
    {
        // Setup seed
        Exp1Simulation::stdSeed0 = Exp1Simulation::stdSeed1 = time(NULL);

        // Setup sets of variables
        // K, DistanceType
        vector<int> Ks;
        Ks.push_back(100);
        Ks.push_back(1000);
        vector<DistanceType> distanceTypes;
        distanceTypes.push_back(DistanceType::Euclidean);
        // XXX: Set them up

        // Main Loop
        for (int K : Ks) {
            for (DistanceType dType : distanceTypes) {
                runOneExperiment(K, dType);
            }
        }
    }

    void static runOneExperiment(int K, DistanceType dType)
    {
        // Make new file
        string fn = "exp1-data-" + to_string(K) + "-" + distanceType2str(dType);
        ofstream out(fn);

        // Prepare Ns
        vector<int> Ns;
        for (int i = 1; i <= 10; i++)
            Ns.push_back(i);
        for (int i = 20; i <= 100; i += 10)
            Ns.push_back(i);

        // Produce and log results for each N
        for (int N : Ns)
        {
            Exp1Simulation exp1(Exp1Simulation::stdSeed0, Exp1Simulation::stdSeed1, K, N);
            out << exp1.getAverageRatioDist(dType) << endl;
        }

        // Close the file
        out.close();
    }


};

int Exp1Simulation::stdSeed0, Exp1Simulation::stdSeed1;

int main(int argc, char **argv)
{

    // Set output params
    cout << fixed;
    cout << setprecision(4);

    // Set seed
    Exp1Simulation::runFullExperiment();

    return 0;
}
