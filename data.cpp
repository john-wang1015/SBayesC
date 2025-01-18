#include "data.hpp"

BinaryFileReader::BinaryFileReader(const std::string& filePath) : filePath(filePath) {}

std::vector<char> BinaryFileReader::readAll() {
    std::ifstream file(filePath, std::ios::binary | std::ios::ate);
    if (!file) {
        throw std::runtime_error("Failed to open file: " + filePath);
    }

    std::streamsize size = file.tellg();
    file.seekg(0, std::ios::beg);

    std::vector<char> buffer(size);
    if (!file.read(buffer.data(), size)) {
        throw std::runtime_error("Failed to read file: " + filePath);
    }

    return buffer;
}

void readBinFullLD(const std::string& filePath) {
    std::ifstream file(filePath, std::ios::binary | std::ios::in);
}

void readSummary(const std::string& filePath) {

}

