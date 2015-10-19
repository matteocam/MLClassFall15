#include <iostream>
#include <cstdlib>
#include <cassert>
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

class Point
{
public:
    Point(int N) { }
    // XXX: Who deletes this?
    static Point *mkRandomPoint(int N)
    {
        return nullptr; // XXX
    }

    float getDistance(const Point &p)
    {
        return 0; // XXX
    }

};

class RandomPointSet
{
    public:
    RandomPointSet(int K, int N)
    {
    }

    const Point &getClosestPoint(const Point &p)
    {

    }

    const Point &getFarthestPoint(const Point &p)
    {

    }

};

class Exp1Program
{
private:
    RandomPointSet *pointSet = nullptr;
public:
    // NOTE: At some point you may want to variate on K and N and distance function
    Exp1Program(int seed0, int seed1, int K, int N)
    {
        UnifGen::setSeed(seed0,seed1);
        pointSet = new RandomPointSet(K, N);
    }

    // Some plotting methods here

};

int main(int argc, char **argv)
{
    int N, K;
    if (argc < 3) { // Use default
        N = 1;
        K = 2;
    } else {
        cerr << "Unimplemented\n!";
    }

    Exp1Program exp1(1, 1, K, N);


    return 0;
}
