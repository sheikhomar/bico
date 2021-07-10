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
    size_t DimSize;
    size_t DataSize;
    size_t ClusterSize;
    size_t LowDimSize;
    size_t TargetCoresetSize;
    std::string InputFilePath;
    std::string OutputFilePath;

public:
    void outputResultsToFile(ProxySolution<Point> *sol)
    {
        printf("Write results to %s...\n", OutputFilePath.c_str());

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

    virtual void prepareFileStream(std::istream inData)
    {
        throw "prepareFileStream not implemented!";
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
        size_t pointCount = 0;
        StopWatch sw(true);
        Bico<Point> bico(DimSize, DataSize, ClusterSize, LowDimSize, TargetCoresetSize, new SquaredL2Metric(), new PointWeightModifier());

        while (inData.good())
        {
            // Read line and construct point
            std::getline(inData, line);
            std::vector<double> coords;
            parseLine(coords, line);
            CluE::Point p(coords);

            if (p.dimension() != DimSize)
            {
                std::clog << "Line skipped because line dimension is " << p.dimension() << " instead of " << DimSize << std::endl;
                continue;
            }

            pointCount++;

            if (pointCount % 10000 == 0)
            {
                std::cout << "Read " << pointCount << " points. Run time: " << sw.elapsedStr() << std::endl;
            }

            // Call BICO point update
            bico << p;

            // p.debug(pointCount, "%5.0f", 15);
            // if (pointCount > 5) {
            //     break;
            // }
        }

        std::cout << "Processed " << pointCount << " points. Run time: " << sw.elapsedStr() << "s" << std::endl;

        outputResultsToFile(bico.compute());
    }
};

class CensusExperiment : public Experiment
{
public:
    CensusExperiment()
    {
        this->DimSize = 68UL;
        this->DataSize = 2458285UL;
        this->ClusterSize = 200UL;
        this->LowDimSize = 50UL;
        this->TargetCoresetSize = 40000UL;
        this->InputFilePath = "data/raw/USCensus1990.data.txt";
        this->OutputFilePath = "data/results/USCensus1990.data.txt";
    }

    void prepareFileStream(std::istream inData)
    {
        std::string line;
        std::getline(inData, line); // Ignore the header line.
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
    CovertypeExperiment()
    {
        this->DimSize = 55UL;
        this->DataSize = 581012UL;
        this->ClusterSize = 200UL;
        this->LowDimSize = 50UL;
        this->TargetCoresetSize = 40000UL;
        this->InputFilePath = "data/raw/covtype.data.gz";
        this->OutputFilePath = "data/results/covtype.txt";
    }

    void prepareFileStream(std::istream inData)
    {
        std::string line;
        std::getline(inData, line); // Ignore the header line.
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
