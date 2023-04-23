/**
 * @file FileWriter.hpp
 * @author Sam (ge96daf)
 * @brief Main class to construct ofstream instances
 * @version 0.1
 * @date 2023-01-06
 *
 * @copyright Copyright (c) 2023
 *
 */

#ifndef _FILE_WRITER_H_
#define _FILE_WRITER_H_

#include "HelperFunctions.hpp"
#include <fstream>
#include <sstream>

namespace utilities
{

class FileWriter final
{
  public:
    FileWriter() = default;
    ~FileWriter() { closeFile(); }

    inline bool writeData(std::string&& data)
    {
        if (!file_.is_open() || !file_.good())
        {
            LOG_CRI << "Failed to write to: " << file_name_;
            file_.clear(); // REVIEW needed here?
            return false;
        }

        file_ << std::move(data);

        return file_.good();
    }

    inline bool openFile(const std::string& file_path)
    {
        file_name_ = file_path;
        if (file_.is_open())
        {
            LOG_DBG("File is already open: " << file_path);
            return file_.good();
        }

        file_.open(file_path, std::ios_base::out);
        if (!file_.is_open())
        {
            LOG_CRI << "Failed to open file: " << file_path;
        }

        return file_.good();
    }

    inline void closeFile()
    {
        data_.flush();
        file_ << std::move(data_.str()); // REVIEW is this right place?
        file_.flush();
        file_.close();
    }

    inline void openWriteCloseFileInPlace(const std::string& file_path,
                                          std::string&&      data)
    {
        openFile(std::move(file_path));

        writeData(std::move(data));

        closeFile();
    }

    void reset()
    {
        file_.clear();

        if (file_.is_open())
        {
            file_.flush();
            file_.close();
        }

        data_.clear();
    }

    // Operator for stream overloads
    friend FileWriter& operator<<(FileWriter& file, const std::string& data);
    friend FileWriter& operator<<(FileWriter& file, const std::int32_t& data);
    friend FileWriter& operator<<(FileWriter& file, const std::uint32_t& data);
    friend FileWriter& operator<<(FileWriter& file, const std::int64_t& data);
    friend FileWriter& operator<<(FileWriter& file, const std::uint64_t& data);
    friend FileWriter& operator<<(FileWriter& file, const float& data);
    friend FileWriter& operator<<(FileWriter& file, const double& data);

  private:
    std::ofstream      file_;
    std::ostringstream data_;
    std::string        file_name_;
};

inline FileWriter& operator<<(FileWriter& file, const std::string& data)
{
    file.data_ << std::move(data);
    return file;
}

inline FileWriter& operator<<(FileWriter& file, const std::int32_t& data)
{
    file.data_ << std::move(data);
    return file;
}

inline FileWriter& operator<<(FileWriter& file, const std::uint32_t& data)
{
    file.data_ << std::move(data);
    return file;
}

inline FileWriter& operator<<(FileWriter& file, const std::int64_t& data)
{
    file.data_ << std::move(data);
    return file;
}

inline FileWriter& operator<<(FileWriter& file, const std::uint64_t& data)
{
    file.data_ << std::move(data);
    return file;
}

inline FileWriter& operator<<(FileWriter& file, const float& data)
{
    file.data_ << std::move(data);
    return file;
}

inline FileWriter& operator<<(FileWriter& file, const double& data)
{
    file.data_ << std::move(data);
    return file;
}

} // namespace utilities

#endif // _FILE_WRITER_H_
