#include <iostream>
#include <sstream>
#include <fstream>
#include <random>
#include <ctime>
#include <time.h>
#include <chrono>
#include <iomanip>

#include <boost/algorithm/string.hpp>

#include "src/point/l2metric.h"
#include "src/point/squaredl2metric.h"
#include "src/point/point.h"
#include "src/point/pointweightmodifier.h"
#include "src/clustering/bico.h"
#include "src/misc/randomness.h"
#include "src/misc/randomgenerator.h"
#include "src/misc/stopwatch.h"
#include "src/datastructure/proxysolution.h"
#include "src/point/pointcentroid.h"
#include "src/point/pointweightmodifier.h"
#include "src/point/realspaceprovider.h"

using namespace CluE;
using namespace std::chrono;

class Experiment
{
protected:
    const size_t D;
    const size_t N;
    const size_t K;
    const size_t P;
    const size_t T;
    const std::string InputFilePath;
    const bool HasHeader;
    const std::string OutputFilePath;
    Bico<Point> bico;

public:
    Experiment(size_t d, size_t n, size_t k, size_t p, size_t t, std::string inputFilePath, bool hasHeader, std::string outputFilePath) : D(d), N(n), K(k), P(p), T(t), InputFilePath(inputFilePath), HasHeader(hasHeader), OutputFilePath(outputFilePath), bico(d, n, k, p, t, new SquaredL2Metric(), new PointWeightModifier())
    {
    }

    void outputResultsToFile()
    {
        printf("Write results to %s...\n", OutputFilePath.c_str());

        // Retrieve coreset
        ProxySolution<Point> *sol = bico.compute();

        std::ofstream outData(OutputFilePath, std::ifstream::out);

        // Output coreset size
        outData << sol->proxysets[0].size() << "\n";

        // Output coreset points
        for (size_t i = 0; i < sol->proxysets[0].size(); ++i)
        {
            // Output weight
            outData << sol->proxysets[0][i].getWeight() << " ";
            // Output center of gravity
            for (size_t j = 0; j < sol->proxysets[0][i].dimension(); ++j)
            {
                outData << sol->proxysets[0][i][j];
                if (j < sol->proxysets[0][i].dimension() - 1)
                    outData << " ";
            }
            outData << "\n";
        }
        outData.close();
    }

    virtual void parseLine(std::vector<double> &result, const std::string &line)
    {
        throw "ParseLine not implemented!";
    }

    void run()
    {
        printf("Opening input file %s...\n", InputFilePath.c_str());
        std::ifstream inData(InputFilePath, std::ifstream::in);

        std::string line;

        if (HasHeader)
        {
            // Skip the first line because it is the header.
            std::getline(inData, line);
        }

            size_t pointCount = 0;

        StopWatch sw(true);

        while (inData.good())
        {
            // Read line and construct point
            std::getline(inData, line);
            std::vector<double> coords;
            parseLine(coords, line);
            CluE::Point p(coords);

            if (p.dimension() != D)
            {
                std::clog << "Line skipped because line dimension is " << p.dimension() << " instead of " << D << std::endl;
                continue;
            }

            pointCount++;

            if (pointCount % 10000 == 0)
            {
                std::cout << "Read " << pointCount << " points. Run time: " << sw.elapsedStr() << std::endl;
            }

            // Call BICO point update
            bico << p;
        }

        std::cout << "Processed " << pointCount << " points. Run time: " << sw.elapsedStr() << "s" << std::endl;

        outputResultsToFile();

    }
};

class CensusExperiment : public Experiment
{
public:
    CensusExperiment() : Experiment(
        68UL,       // Number of dimensions
        2458285UL,  // Number of points in the dataset.
        200UL,      // Number of clusters.
        50UL,       // Number of random projections
        40000UL,    // Number of target points in the coreset.
        "data/raw/USCensus1990.data.txt",
        true,       // Whether the data contains a header.
        "data/results/USCensus1990.data.txt"
    )
    {
    }

    void parseLine(std::vector<double> &result, const std::string &line)
    {
        std::vector<std::string> stringcoords;
        boost::split(stringcoords, line, boost::is_any_of(","));

        result.reserve(stringcoords.size());

        // Skip the first attribute which is `caseid`
        for (size_t i = 1; i < stringcoords.size(); ++i)
            result.push_back(atof(stringcoords[i].c_str()));
    }
};

int main(int argc, char **argv)
{
    CensusExperiment census;
    census.run();
}

int main_old(int argc, char **argv)
{
    using namespace CluE;

    time_t starttime, endtime;
    double difference;

    if (argc < 8)
    {
        std::cout << "Usage: input n k d space output projections [splitchar [seed]]" << std::endl;
        std::cout << "  input       = path to input file" << std::endl;
        std::cout << "  n           = number of input points" << std::endl;
        std::cout << "  k           = number of desired centers" << std::endl;
        std::cout << "  d           = dimension of an input point" << std::endl;
        std::cout << "  space       = coreset size" << std::endl;
        std::cout << "  output      = path to the output file" << std::endl;
        std::cout << "  projections = number of random projections used for nearest neighbour search" << std::endl;
        std::cout << "                in first level" << std::endl;
        std::cout << "  splitchar   = input CSV split character" << std::endl;
        std::cout << "  seed        = random seed (optional)" << std::endl;
        std::cout << std::endl;
        std::cout << "7 arguments expected, got " << argc - 1 << ":" << std::endl;
        for (int i = 1; i < argc; ++i)
            std::cout << i << ".: " << argv[i] << std::endl;
        return 1;
    }

    // Read arguments
    std::ifstream filestream(argv[1], std::ifstream::in);
    int n = atoi(argv[2]);
    int k = atoi(argv[3]);
    int d = atoi(argv[4]);
    int space = atoi(argv[5]);
    std::ofstream outputstream(argv[6], std::ifstream::out);
    int p = atoi(argv[7]);
    std::string splitchar(",");
    if (argc >= 9)
        splitchar = std::string(argv[8], 1);
    if (argc >= 10)
        Randomness::initialize(atoi(argv[9]));

    time(&starttime);

    // Initialize BICO
    Bico<Point> bico(d, n, k, p, space, new SquaredL2Metric(), new PointWeightModifier());

    int pos = 0;
    while (filestream.good())
    {
        // Read line and construct point
        std::string line;
        std::getline(filestream, line);
        std::vector<std::string> stringcoords;
        boost::split(stringcoords, line, boost::is_any_of(splitchar));

        std::vector<double> coords;
        coords.reserve(stringcoords.size());
        for (size_t i = 0; i < stringcoords.size(); ++i)
            coords.push_back(atof(stringcoords[i].c_str()));
        Point p(coords);

        if (p.dimension() != d)
        {
            std::clog << "Line skipped because line dimension is " << p.dimension() << " instead of " << d << std::endl;
            continue;
        }

        // Call BICO point update
        bico << p;
    }

    // Retrieve coreset
    ProxySolution<Point> *sol = bico.compute();

    // Output coreset size
    outputstream << sol->proxysets[0].size() << "\n";

    // Output coreset points
    for (size_t i = 0; i < sol->proxysets[0].size(); ++i)
    {
        // Output weight
        outputstream << sol->proxysets[0][i].getWeight() << " ";
        // Output center of gravity
        for (size_t j = 0; j < sol->proxysets[0][i].dimension(); ++j)
        {
            outputstream << sol->proxysets[0][i][j];
            if (j < sol->proxysets[0][i].dimension() - 1)
                outputstream << " ";
        }
        outputstream << "\n";
    }
    outputstream.close();

    return 0;
}
