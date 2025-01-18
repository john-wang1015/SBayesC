#ifndef DATA_HPP
#define DATA_HPP

#include <string>
#include <vector>
#include <fstream>
#include <stdexcept>

class readFile {
public:
    void readBinFullLD(const std::string& filePath);
    void readSummary(const std::string& filePath);
    void readBinSparseLD();
    void readBinShrukLD();
};

class BinaryFileReader {
public:
    BinaryFileReader(const std::string& filePath);

    std::vector<char> readAll(); // Reads all data into a buffer

    template <typename T>
    std::vector<T> readAsType(); // Reads binary data as a specific type

private:
    std::string filePath;
};

// Inline implementation of readAsType<T>()
template <typename T>
std::vector<T> BinaryFileReader::readAsType() {
    std::ifstream file(filePath, std::ios::binary | std::ios::in);
    if (!file) {
        throw std::runtime_error("Failed to open file: " + filePath);
    }

    std::streamsize size = file.tellg();
    file.seekg(0, std::ios::beg);

    if (size % sizeof(T) != 0) {
        throw std::runtime_error("File size is not a multiple of type size");
    }

    std::vector<T> buffer(size / sizeof(T));
    if (!file.read(reinterpret_cast<char*>(buffer.data()), size)) {
        throw std::runtime_error("Failed to read file as type: " + filePath);
    }

    return buffer;
}

#endif // DATA_HPP
