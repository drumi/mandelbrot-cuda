#pragma once

#include <cstdint>

namespace BMPImage 
{
    void save(char const* path, int const height, int const width, uint8_t const* rawImage);
}