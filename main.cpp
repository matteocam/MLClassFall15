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

// To use uniformly/normally generated points (EXP1/EXP2)
#define EXP2
// To do classification
//#define CLASS_EXP

// Fractional Distance only
#define FRAC_DIST

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
    Max,
    Fractional
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
        case DistanceType::Fractional:
            return "Fractional";
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
        if (N > 1) { // Something wrong here
          int i;
          float sumSquares = 0.0;
          for (i = 0; i < N; i++)
            sumSquares += components[i]*components[i];
          return sqrt(sumSquares);
        } else
          return components[0];
    }

    void normalize()
    {
        float len = length();
        //cout << *this << endl;
        //cout << "Dividing by " << len << endl;
        for (int i = 0; i < N; i++) {
            components[i] /= len;
        }
        //cout << "NewLength = " << length() << endl;
    }

    const Point *getPerturbedPoint(float sigma) const
    {
        Point *newPt = new Point(N);
        std::normal_distribution<> normalDist(0,sigma);

        for (int i = 0; i < N; i++) {
            float epsilon = normalDist(*mt);
            newPt->components.push_back(components[i] +  epsilon);
        }

        return newPt;
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

        #ifndef CLASS_EXP
        #ifdef EXP2
        // Normalize value
        newPoint->normalize();
        #endif
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
            case DistanceType::Fractional:
                return getFractionalDistance(p);
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


    float getFractionalDistance(const Point &pt) const
    {
        const float p = 1.0/3.0;
        float sumOfPowDiffs = 0.0;
        for (int i = 0; i < N; i++)
        {
            sumOfPowDiffs += pow(fabs(components[i]-pt.components[i]), p);
            //cout << components[i] << " " << pt.components[i] << " -> " << sumOfPowDiffs << endl;
        }
        return sumOfPowDiffs;  // NOTE: It could be raised to the 1/p
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

        int nEqualPoints = 0;

        // XXX: How do we check that it's not the same point? (we check for equality)
        for (const Point *p1 : points) {
            if (*p1 == p) {
                // Do nothing
                nEqualPoints++;
            } else {
                distancePoints.push_back(make_pair(p.getDistance(*p1, distType), p1));
            }
        }

        // Sanity check if we have normalization
        if (distancePoints.empty())
          return p;

        // cout << "I stumbled upon " << nEqualPoints << " equal points\n";

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

        if (distancePoints.empty())
          return p;

        auto farthestPoint = *max_element(distancePoints.begin(),distancePoints.end());
        return *farthestPoint.second;
    }

    const PointSeq &getPoints()
    {
        return points;
    }

    RandomPointSet *mkPerturbedSet(float sigma)
    {
        auto newSet = new RandomPointSet(K, N);
        for (int i = 0; i < K; i++) {
            delete newSet->points[i];
            newSet->points[i] = points[i]->getPerturbedPoint(sigma);
        }
        return newSet;
    }

    void assignFrom(RandomPointSet *classSource)
    {
        for (const Point *pt : points)
        {
            auto closestPoint = classSource->getClosestPoint(*pt, DistanceType::Euclidean);
            const_cast<Point *>(pt)->label = closestPoint.label;
        }
    }

    int getLabel(int i) const
    {
        return points[i]->label;
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

            float r = (distFarthest > 0 ? distClosest / distFarthest : 0.0);

            out << Observation(K, N, distType, r);
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
        //Ks.push_back(1000);
        //Ks.push_back(2000);
        //Ks.push_back(5000);
        //Ks.push_back(10000);
        Ks.push_back(20000);
        //Ks.push_back(1000);
        vector<DistanceType> distanceTypes;
        #ifdef FRAC_DIST
        distanceTypes.push_back(DistanceType::Fractional);
        #else
        distanceTypes.push_back(DistanceType::Euclidean);
        distanceTypes.push_back(DistanceType::CityBlock);
        distanceTypes.push_back(DistanceType::Max);
        #endif
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
        #ifndef FRAC_DIST
        for (int i = 1; i <= 10; i++)
            Ns.push_back(i);
        #else
        for (int i = 8; i <= 10; i++)
            Ns.push_back(i);

        #endif
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

struct ClassificationObservation
{
    int K,N;

    float sigma;

    int correctOnes;
    int onesAsZeros;
    int correctZeros;
    int zerosAsOnes;



    ClassificationObservation(int _K, int _N, float _sigma, int _correctOnes, int _onesAsZeros, int _correctZeros, int _zerosAsOnes) :
        K(_K), N(_N), sigma(_sigma),
        correctOnes(_correctOnes), onesAsZeros(_onesAsZeros), correctZeros(_correctZeros), zerosAsOnes(_zerosAsOnes)

    {
    }

    static string header()
    {
        return "K\tN\tsigma\tcorrectOnes\tonesAsZeros\tcorrectZeros\tzerosAsOnes\n";
    }

    friend ostream& operator<<(ostream& os, const ClassificationObservation& obs);

};

ostream& operator<<(ostream& os, const ClassificationObservation& obs)
{
    os << obs.K << "\t" << obs.N << "\t"
        << obs.sigma << "\t" << obs.correctOnes
        << "\t" << obs.onesAsZeros << "\t" <<
          obs.correctZeros << "\t" <<
          obs.zerosAsOnes << endl;
    return os;
}

class ClassificationSimulation
{
private:
    int K, N;
    float sigma;
    RandomPointSet *X = nullptr;
    RandomPointSet *Y = nullptr;
    RandomPointSet *Y1 = nullptr;
    ofstream &out;

public:

    ClassificationSimulation(int _K, int _N, float _sigma, ofstream &_out) :
        K(_K), N(_N), sigma(_sigma), out(_out)
    {
        // Make X set
        X = new RandomPointSet(K, N);

        // Make Y Set
        Y = new RandomPointSet(K, N);
        // Find class of y from x
        Y->assignFrom(X);

        // Perturb Y
        Y1 = Y->mkPerturbedSet(sigma);

        // Find the closest point
        Y1->assignFrom(X);

        cout << "X:\n" << *X;
        cout << "Y:\n" << *Y;
        cout << "Y1:\n" << *Y1;

    }

    ~ClassificationSimulation()
    {
        delete X;
        delete Y;
        delete Y1;
    }

    static void runFullExperiment()
    {
        // Setup sets of variables
        // K, DistanceType
        vector<int> Ks;
        Ks.push_back(1000);
        //Ks.push_back(2000);
        Ks.push_back(5000);
        //Ks.push_back(10000);
        //Ks.push_back(20000);

        vector<float> sigmas;
        sigmas.push_back(.1);
        sigmas.push_back(.01);
        // XXX: Set them up

        // Set logging file
        string fn = "exp-class-analysis.dat";

        ofstream out(fn);
        out << ClassificationObservation::header();

        // Main Loop
        for (int K : Ks) {
            for (float sigma : sigmas) {
                runOneExperiment(K, sigma, out);
            }
        }

        out.close();

    }

    static void runOneExperiment(int K, float sigma, ofstream &out)
    {
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
                 << N << " and sigma = " << sigma << endl;
            ClassificationSimulation classSim(K, N, sigma, out);
            classSim.gatherObservations();
        }

    }

    void gatherObservations()
    {
        int oneClassifiedAsZero = 0;
        int zeroClassifiedAsOne = 0;
        int correctOnes = 0;
        int correctZeros = 0;
        // Compare classes in Y and Y1
        for (int i = 0; i < K; i++)
        {
            if (Y->getLabel(i) == 0)
            {
                if (Y1->getLabel(i) == 0)
                    correctZeros++;
                else
                    zeroClassifiedAsOne++;
            }

            if (Y->getLabel(i) == 1)
            {
                if (Y1->getLabel(i) == 1)
                    correctOnes++;
                else
                    oneClassifiedAsZero++;
            }
        }

        out << ClassificationObservation(K, N, sigma, correctOnes, oneClassifiedAsZero, correctZeros, zeroClassifiedAsOne);
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

    //cout << pow(0.156, 1.0/3.0) << endl;
    //return 0;

    #ifdef CLASS_EXP
      ClassificationSimulation::runFullExperiment();
    #else
      Exp1Simulation::runFullExperiment();
    #endif

    return 0;
}
