#include <iostream>
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

    static float getUnifNum(int x[])
    {

        long long  j,k,a=1999,b=4444,c=2147483647;

        j=a*x[0]+b*x[1];
        k=j-c*(int)(j/c);
        x[1]=x[0];
        x[0]=(int)k;
        return (float)(k)/2147483647.f;
    }

};

int UnifGen::seed[2];

// N is the dimension of the space
// K is the number of points

class Point
{
    public:
    Point(int N);
    // XXX: Who deletes this?
    static Point *mkRandomPoint(int N)
    {
    }
};

class RandomPointSet
{
    public:
    RandomPointSet(int K, int N)
    {
    }

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

    UnifGen::setSeed(1,1);
    RandomPointSet pointSet(K, N);


    return 0;
}
