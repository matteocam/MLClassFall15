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
#include <random>
using namespace std;

#define EXP2


mt19937 *mt;


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
        uniform_real_distribution<float> unifDist(0.0, 1.0);
        float res = unifDist(*mt);


        /*
        long long  j,k,a=1999,b=4444,c=2147483647;

        j=a*seed[0]+b*seed[1];
        k=j-c*(int)(j/c);
        seed[1]=seed[0];
        seed[0]=(int)k;
        float res =  (float)(k)/2147483647.f;
        * */

        assert (res >= 0 && res <= 1);
        return res;
    }

    static float getNormalNum()
    {
        std::normal_distribution<> normalDist(0,1);
        float res = normalDist(*mt);
        return res;
    }

};

int UnifGen::seed[2];

// N is the dimension of the space
// K is the number of points

enum DistanceType
{
    Euclidean,
    CityBlock,
    Max
};

string distanceType2str(DistanceType dt)
{
    switch(dt)
    {
        case DistanceType::Euclidean:
            return "Euclidean";
            break;
        case DistanceType::CityBlock:
            return "CityBlock";
            break;
        case DistanceType::Max:
            return "Max";
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

    float length() const
    {
        return *inner_product(components.begin(), components.end(), components.begin(), components.end());
    }

    void normalize()
    {
        float len = length();
        for (int i = 0; i < N; i++) {
            components[i] /= len;
        }
    }


    // XXX: Who deletes this?
    static const Point *mkRandomPoint(int N)
    {
        Point *newPoint = new Point(N);
        for (int i = 0; i < N; i++) {
            #ifdef EXP1
            float value = UnifGen::getUnifNum();
            #endif
            #ifdef EXP2
            float value = UnifGen::getNormalNum();
            #endif
            newPoint->components.push_back(value);
        }

        #ifdef EXP2
        // Normalize value

        #endif

        newPoint->label = (UnifGen::getUnifNum() < .5);

        return newPoint;
    }

    float getDistance(const Point &p, DistanceType distType) const
    {
        switch(distType) {
            case DistanceType::Euclidean:
                return getEuclDistance(p);
                break;
            case DistanceType::CityBlock:
                return getCityBlockDistance(p);
                break;
            case DistanceType::Max:
                return getMaxDistance(p);
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

    float getCityBlockDistance(const Point &p) const
    {
        float sumOfAbsDiffs = 0.0;
        for (int i = 0; i < N; i++)
        {
            sumOfAbsDiffs += fabs(components[i]-p.components[i]);
        }
        return sumOfAbsDiffs;
    }

    float getMaxDistance(const Point &p) const
    {
        vector<float> absDiffs;
        for (int i = 0; i < N; i++)
        {
            absDiffs.push_back(fabs(components[i]-p.components[i]));
        }
        float maxDiff = *max_element(absDiffs.begin(), absDiffs.end());
        return maxDiff;

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

    ~RandomPointSet()
    {
        for (const Point *p : points)
            delete p;
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

struct Observation
{
    int K, N;
    DistanceType distType;
    float obsVal;

    Observation(int _K, int _N, DistanceType _distType, float _obsVal) :
        K(_K), N(_N), distType(_distType), obsVal(_obsVal)
    {
    }

    static string header()
    {
        return "K\tN\tDistanceType\tr\n";
    }

    friend ostream& operator<<(ostream& os, const Observation& obs);

};

ostream& operator<<(ostream& os, const Observation& obs)
{
    os << obs.K << "\t" << obs.N << "\t" <<
        distanceType2str(obs.distType) << "\t" <<
        obs.obsVal << endl;
    return os;
}

class Exp1Simulation
{
private:
    RandomPointSet *pointSet = nullptr;
public:
    static int stdSeed0, stdSeed1;
    int K, N;
    ofstream &out;

    // NOTE: At some point you may want to variate on K and N and distance function
    Exp1Simulation(int _K, int _N, ofstream &_out) :
        K(_K), N(_N), out(_out)
    {
        pointSet = new RandomPointSet(K, N);

        //cout << *pointSet;

    }

    ~Exp1Simulation()
    {
        delete pointSet;
    }

    void gatherObservations(DistanceType distType)
    {
        auto points = pointSet->getPoints();
        //float sumDists = 0.0;

        // XXX: O(K^2) algorithm!
        for (const Point *pPoint : points)
        {
            auto closestPoint = pointSet->getClosestPoint(*pPoint, distType);
            auto farthestPoint = pointSet->getFarthestPoint(*pPoint, distType);

            float distClosest = pPoint->getDistance(closestPoint, distType);
            float distFarthest = pPoint->getDistance(farthestPoint, distType);

            out << Observation(K, N, distType, distClosest / distFarthest);
            //sumDists += (distClosest / distFarthest);
        }
        //return sumDists / pointSet->getK();
    }

    void logResults(ofstream &outFile)
    {

    }

    void static runFullExperiment()
    {
        // Setup seed
        Exp1Simulation::stdSeed0 = Exp1Simulation::stdSeed1 = time(NULL);
        UnifGen::setSeed(stdSeed0,stdSeed1); // XXX

        // Setup sets of variables
        // K, DistanceType
        vector<int> Ks;
        Ks.push_back(1000);
        //Ks.push_back(2000);
        Ks.push_back(5000);
        Ks.push_back(10000);
        Ks.push_back(20000);
        //Ks.push_back(1000);
        vector<DistanceType> distanceTypes;
        distanceTypes.push_back(DistanceType::Euclidean);
        distanceTypes.push_back(DistanceType::CityBlock);
        distanceTypes.push_back(DistanceType::Max);
        // XXX: Set them up

        // Set logging file
        #ifdef EXP1
        string fn = "exp1-analysis.dat";
        #endif
        #ifdef EXP2
        string fn = "exp2-analysis.dat";
        #endif
        ofstream out(fn);
        out << Observation::header();

        // Main Loop
        for (int K : Ks) {
            for (DistanceType dType : distanceTypes) {
                runOneExperiment(K, dType, out);
            }
        }

        out.close();
    }

    void static runOneExperiment(int K, DistanceType dType, ofstream &out)
    {
        // Make new file
        //string fn = "exp1-data-" + to_string(K) + "-" + distanceType2str(dType);
        //ofstream out(fn);

        // Prepare Ns
        vector<int> Ns;
        for (int i = 1; i <= 10; i++)
            Ns.push_back(i);
        for (int i = 20; i <= 100; i += 10)
            Ns.push_back(i);

        // Produce and log results for each N
        for (int N : Ns)
        {
            cout << "Working on experiment for K = " << K << ", N = "
                 << N << " and DistanceType = " << distanceType2str(dType) << endl;
            Exp1Simulation exp1(K, N, out);
            exp1.gatherObservations(dType);
        }

    }

};

int Exp1Simulation::stdSeed0, Exp1Simulation::stdSeed1;

int main(int argc, char **argv)
{
    // Set randomness
    random_device rd;
    mt = new mt19937(rd());

    // Set output params
    cout << fixed;
    cout << setprecision(4);

    // Set seed
    Exp1Simulation::runFullExperiment();

    return 0;
}
