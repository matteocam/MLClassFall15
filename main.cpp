#include <iostream>
using namespace std;

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

    RandomPointSet pointSet(K, N);


    return 0;
}
