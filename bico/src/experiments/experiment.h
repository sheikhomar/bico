#ifndef EXPERIMENT_H
#define EXPERIMENT_H

#include <iostream>
#include <sstream>
#include <fstream>
#include <random>
#include <ctime>
#include <time.h>
#include <chrono>
#include <iomanip>

#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filter/bzip2.hpp>

#include "../point/l2metric.h"
#include "../point/squaredl2metric.h"
#include "../point/point.h"
#include "../point/pointweightmodifier.h"
#include "../clustering/bico.h"
#include "../misc/randomness.h"
#include "../misc/randomgenerator.h"
#include "../misc/stopwatch.h"
#include "../datastructure/proxysolution.h"
#include "../point/pointcentroid.h"
#include "../point/pointweightmodifier.h"
#include "../point/realspaceprovider.h"

using namespace CluE;

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

        namespace io = boost::iostreams;
        std::ifstream fileStream(InputFilePath, std::ios_base::in | std::ios_base::binary);
        io::filtering_streambuf<io::input> filteredInputStream;
        if (boost::ends_with(InputFilePath, ".gz"))
        {
            filteredInputStream.push(io::gzip_decompressor());
        }
        filteredInputStream.push(fileStream);
        std::istream inData(&filteredInputStream);

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
                             68UL,      // Number of dimensions
                             2458285UL, // Number of points in the dataset.
                             200UL,     // Number of clusters.
                             50UL,      // Number of random projections
                             40000UL,   // Number of target points in the coreset.
                             "data/raw/USCensus1990.data.txt",
                             true, // Whether the data contains a header.
                             "data/results/USCensus1990.data.txt")
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

class CovertypeExperiment : public Experiment
{
public:
    CovertypeExperiment() : Experiment(
                                55UL,     // Number of dimensions, 54 variables and 1 label
                                581012UL, // Number of points in the dataset.
                                200UL,    // Number of clusters.
                                50UL,     // Number of random projections
                                40000UL,  // Number of target points in the coreset.
                                "data/raw/covtype.data.gz",
                                false, // Whether the data contains a header.
                                "data/results/covtype.txt")
    {
    }

    void parseLine(std::vector<double> &result, const std::string &line)
    {
        std::vector<std::string> stringcoords;
        boost::split(stringcoords, line, boost::is_any_of(","));

        result.reserve(stringcoords.size());

        for (size_t i = 0; i < stringcoords.size(); ++i)
            result.push_back(atof(stringcoords[i].c_str()));
    }
};

#endif
