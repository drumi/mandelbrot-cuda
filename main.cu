#include "bmp.h"

#include <iostream>
#include <string>
#include <complex>
#include <thread>
#include <vector>
#include <atomic>

struct Complex
{
    float real;
    float imag;
};

namespace 
{
    namespace DEFAULT
    {
        namespace FLAG
        {
            char const* const FILE_NAME     = "-o";
            char const* const IMAGE_WIDTH   = "-w";
            char const* const IMAGE_HEIGHT  = "-h";
            char const* const GRANULARITY   = "-g";
            char const* const THREADS_COUNT = "-t";
            char const* const ZOOM_LEVEL    = "-z";
            char const* const POINT_ORIGIN  = "-p";
            char const* const ITERATIONS    = "-c";
        }   

        namespace IMAGE
        {
            char const* NAME = "mandelbrot.bmp";
            int const WIDTH = 3840;
            int const HEIGHT = 2160;
            float const ZOOM_LEVEL = 1.0;
            int const BYTES_PER_PIXEL = 3;
            Complex const POINT_ORIGIN = {0, 0};

            uint8_t const TINT_ON_ESCAPE = 32;
        }

        namespace THREADS
        {
            int const GRANULARITY = 1;
            int const ITERATIONS = 256;
            int const COUNT = std::thread::hardware_concurrency();
        }

        float const INFINITY_THRESHOLD = 4.0;
    }
}

#include <chrono>
namespace
{
    class Clock
    {
        std::chrono::time_point<std::chrono::high_resolution_clock> start = std::chrono::high_resolution_clock::now();

        public:
        int getElapsedMilliseconds() const
        {
            std::chrono::time_point<std::chrono::high_resolution_clock> now = std::chrono::high_resolution_clock::now();
            return std::chrono::duration_cast<std::chrono::milliseconds>(now - start).count();
        }
    };
}


struct ProgramParameters
{
    int    imageWidth      = DEFAULT::IMAGE::WIDTH;
    int    imageHeight     = DEFAULT::IMAGE::HEIGHT;
    int    threadsCount    = DEFAULT::THREADS::COUNT;
    int    granularity     = DEFAULT::THREADS::GRANULARITY;
    int    iterationsCount = DEFAULT::THREADS::ITERATIONS;
    float zoomLevel       = DEFAULT::IMAGE::ZOOM_LEVEL;

    char const* imageOutputName = DEFAULT::IMAGE::NAME;
    Complex pointOrigin = DEFAULT::IMAGE::POINT_ORIGIN;
};

struct ThreadParameters
{
    int chunkSize;
    int chunksCount;
    int remainderChunkSize;
    int imageTotalSize;
    int bytesPerWidth;
    int bytesPerHeight;

    float dx;
    float dy;

    Complex bottomLeftCoordinates;
    Complex upperRightCoordinates;
};

void printExecutingParameters(ProgramParameters const p)
{
    std::cout << '\n'
              << DEFAULT::FLAG::FILE_NAME     << " for file name in bmp format.  " << "Executing: " << DEFAULT::FLAG::FILE_NAME     <<" " << p.imageOutputName    << "\n"
              << DEFAULT::FLAG::IMAGE_WIDTH   << " for image width.              " << "Executing: " << DEFAULT::FLAG::IMAGE_WIDTH   <<" " << p.imageWidth         << "\n"
              << DEFAULT::FLAG::IMAGE_HEIGHT  << " for image height.             " << "Executing: " << DEFAULT::FLAG::IMAGE_HEIGHT  <<" " << p.imageHeight        << "\n"
              << DEFAULT::FLAG::GRANULARITY   << " for granularity.              " << "Executing: " << DEFAULT::FLAG::GRANULARITY   <<" " << p.granularity        << "\n"
              << DEFAULT::FLAG::THREADS_COUNT << " for thread count.             " << "Executing: " << DEFAULT::FLAG::THREADS_COUNT <<" " << p.threadsCount       << "\n"
              << DEFAULT::FLAG::ZOOM_LEVEL    << " for image zoom.               " << "Executing: " << DEFAULT::FLAG::ZOOM_LEVEL    <<" " << p.zoomLevel          << "\n"
              << DEFAULT::FLAG::ITERATIONS    << " for complex iterations count. " << "Executing: " << DEFAULT::FLAG::ITERATIONS    <<" " << p.iterationsCount    << "\n"
              << DEFAULT::FLAG::POINT_ORIGIN  << " for image center point.       " << "Executing: " << DEFAULT::FLAG::POINT_ORIGIN  <<" " << p.pointOrigin.real << " " << p.pointOrigin.imag << "\n"
              << '\n';
}

ProgramParameters handleInput(int argc, const char** argv)
{
    ProgramParameters result;

    int i = 1;
    while(i < argc)
    {
        std::string inputFlag(argv[i]);

        if(inputFlag == DEFAULT::FLAG::FILE_NAME)
        {
            //result.imageOutputName = argv[i+1].c_str();
            i += 2;
        }
        else if(inputFlag == DEFAULT::FLAG::IMAGE_WIDTH)
        {
            result.imageWidth = atoi(argv[i+1]);
            i += 2;
        }
        else if(inputFlag == DEFAULT::FLAG::IMAGE_HEIGHT)
        {
            result.imageHeight = atoi(argv[i+1]);
            i += 2;
        }
        else if(inputFlag == DEFAULT::FLAG::GRANULARITY)
        {
            result.granularity = atoi(argv[i+1]);
            i += 2;
        }
        else if(inputFlag == DEFAULT::FLAG::THREADS_COUNT) 
        {
           // result.threadsCount = atoi(argv[i+1]);
            i += 2;;
        }
        else if(inputFlag == DEFAULT::FLAG::ZOOM_LEVEL) 
        {
            result.zoomLevel = atof(argv[i+1]);
            i += 2;
        }
        else if(inputFlag == DEFAULT::FLAG::POINT_ORIGIN)
        {
            result.pointOrigin = {(float)atof(argv[i+1]), (float)atof(argv[i+2])};
            i += 3;
        }
        else if(inputFlag == DEFAULT::FLAG::ITERATIONS)
        {
            result.iterationsCount = atoi(argv[i+1]);
            i += 2;
        }
        else
        {
            std::cerr << "Invalid parameter supplied: " << argv[i] << '\n';
            exit(-1);
        }
    }

    return result;
}

ThreadParameters generateThreadParameters(ProgramParameters const p)
{
    ThreadParameters result;

    int const totalPixels = p.imageHeight * p.imageWidth;

    result.imageTotalSize     = totalPixels * DEFAULT::IMAGE::BYTES_PER_PIXEL; 
    result.chunkSize          = (totalPixels / (p.granularity * p.threadsCount)) * DEFAULT::IMAGE::BYTES_PER_PIXEL;
    result.chunksCount        = result.imageTotalSize / result.chunkSize;
    result.remainderChunkSize = result.imageTotalSize % result.chunkSize;
    
    float const zoom = 2.0 / p.zoomLevel;
    float const aspectRatio = p.imageHeight / (float) p.imageWidth;

    result.bottomLeftCoordinates = {-zoom + p.pointOrigin.real, -zoom * aspectRatio + p.pointOrigin.imag};
    result.upperRightCoordinates = { zoom + p.pointOrigin.real,  zoom * aspectRatio + p.pointOrigin.imag};
    
    result.bytesPerWidth = p.imageWidth * DEFAULT::IMAGE::BYTES_PER_PIXEL;
    result.bytesPerHeight = p.imageHeight;

    result.dx = (result.upperRightCoordinates.real - result.bottomLeftCoordinates.real);
    result.dy = (result.upperRightCoordinates.imag - result.bottomLeftCoordinates.imag);

    return result;
}

__device__
int computeSteps(int const iterations, float real, float imag)
{
    float currR = real;
    float currI = imag;

    float squaredR = currR * currR;
    float squaredI = currI * currI;

    for (int i = 1; i <= iterations; ++i)
    {
        float const growthIndex = squaredR + squaredR;

        if(growthIndex > 4.0)
            return i;

        currI = 2 * currI * currR + imag;
        currR = squaredR - squaredI + real;

        squaredR = currR * currR;
        squaredI = currI * currI;
    }

    return 0;
}

__device__
void computePortionOfImage(int const imageStartIndex, int const imageEndIndex, int const iterations, ThreadParameters const t, uint8_t* rawImage)
{
    for(int i = imageStartIndex; i < imageEndIndex; i += 3)
    {
        int const y = i / t.bytesPerWidth;
        int const x = i % t.bytesPerWidth;

        float const realFraction = (x / (float)t.bytesPerWidth);
        float const imagFraction = (y / (float)t.bytesPerHeight);

        float const real = realFraction * t.dx + t.bottomLeftCoordinates.real;
        float const imag = imagFraction * t.dy + t.bottomLeftCoordinates.imag;

        int const steps = computeSteps(iterations, real, imag);
        uint8_t const color = (255 * steps) /iterations;

        rawImage[i  ] = 32 * (steps != 0); // b
        rawImage[i+1] = color; // g
        rawImage[i+2] = 0;     // r; 0 by default
    }
}

__global__
void computeImage(ProgramParameters const p, ThreadParameters const t, uint8_t* rawImage)
{
    int const threadId = blockIdx.x * blockDim.x + threadIdx.x;
    int currentChunkNumber = threadId - p.threadsCount;

    // Handle normal chunks
    while ((currentChunkNumber += p.threadsCount) < t.chunksCount)
    {
        int const imageStartIndex = currentChunkNumber * t.chunkSize;
        int const imageEndIndex = (currentChunkNumber + 1) * t.chunkSize - 1;

        computePortionOfImage(imageStartIndex, imageEndIndex, p.iterationsCount, t, rawImage);
    }

    // Handle remainder
    if(t.remainderChunkSize != 0 && (currentChunkNumber == t.chunksCount))
    {
        int const imageStartIndex = currentChunkNumber * t.chunkSize;
        int const imageEndIndex = t.imageTotalSize - 1;

        computePortionOfImage(imageStartIndex, imageEndIndex, p.iterationsCount, t, rawImage);
    }
}

int main(int const argc, const char** argv) 
{
    Clock const programClock;

    ProgramParameters programParameters = handleInput(argc, argv);
    programParameters.threadsCount = 16384;
    
    printExecutingParameters(programParameters);
    
    ThreadParameters const threadParameters = generateThreadParameters(programParameters);
    ProgramParameters p = programParameters;

    uint8_t* rawImage = nullptr;
    uint8_t* rawImageHost = nullptr;

    rawImageHost = new uint8_t[p.imageWidth * p.imageHeight * DEFAULT::IMAGE::BYTES_PER_PIXEL]();
    cudaMalloc(&rawImage, p.imageWidth * p.imageHeight * DEFAULT::IMAGE::BYTES_PER_PIXEL);

    computeImage<<<128,128>>>(programParameters, threadParameters, rawImage);

    cudaMemcpy(rawImageHost, rawImage, p.imageWidth * p.imageHeight * DEFAULT::IMAGE::BYTES_PER_PIXEL, cudaMemcpyDeviceToHost);

    BMPImage::save(programParameters.imageOutputName, programParameters.imageHeight, programParameters.imageWidth, rawImageHost);
    
    cudaFree(rawImage);
    delete[] rawImageHost;

    std::cout << "Total time for program execution: " << programClock.getElapsedMilliseconds() << "ms\n";

    return 0;
}
