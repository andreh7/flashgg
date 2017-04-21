#ifndef TorchUtils_TorchWriter_h
#define TorchUtils_TorchWriter_h

#include <iostream>
#include <vector>
#include <cassert>

namespace flashgg
{

  /** utility class to write data to a Torch 7 file */
  class TorchWriter
  {
  protected:
    unsigned objectIndex = 1;

    std::ostream &os;
    
  public:
    const unsigned MAGIC_NUMBER = 1;
    const unsigned MAGIC_STRING = 2;
    const unsigned MAGIC_TABLE  = 3;
    const unsigned MAGIC_TORCH  = 4;


    TorchWriter(std::ostream &os_) : 
      os(os_)
      {
      }
    
    template<typename DataType>
    inline void writeType(const DataType &value)
    {
        os.write((const char*) &value, sizeof(value));
    }

    inline void writeInt(int32_t value)
    {
        os.write((const char*) &value, sizeof(value));
    }

    /** assumes that double is 8 bytes */
    inline void writeDouble(double value)
    {
        os.write((const char*) &value, sizeof(value));
    }

    inline void writeLong(int64_t value)
    {
        os.write((const char*)&value, sizeof(value));
    }

    inline void writeFloat(float value)
    {
        os.write((const char*)&value, sizeof(value));
    }


    inline void writeString(const std::string &str)
    {
        size_t len = str.size();
        writeInt(len);

        os.write((const char*) &str[0], len);
    }

    template<typename DataType>
    void writeTypeTensorHelper(const std::vector<unsigned> &sizes, const std::vector<DataType> &data, const std::string &tensorTypeName,
                         const std::string &storageTypeName)
    {
        // data must be stored in the order such the an increase of the index in the last dimension
        // by one corresponds to an indcreas of the index into data by one etc.

        writeInt(MAGIC_TORCH);

        // write object index (we do NOT reuse objects here, we assume that we only write one
        // object per file)
        writeInt(objectIndex++);

        // version string
        writeString("V 1");

        // class name
        writeString(tensorTypeName);

        //----------

        // write the number of coordinates
        writeInt(sizes.size());

        // write the sizes of each dimension
        for (unsigned dimsize : sizes)
            writeLong(dimsize);

        // calculate strides
        std::vector<uint64_t> strides(sizes.size());
        {
            uint64_t product = 1;
            for (unsigned i = sizes.size(); i > 0; --i)
                {
                    strides[i-1] = product;
                    product *= sizes[i-1];
                }

            assert(product == data.size());
        }
        for (uint64_t stride : strides)
            writeLong(stride);

        // write storage offset
        writeLong(1);

        //----------
        // write the FloatStorage object: the actual data
        //----------
        writeInt(MAGIC_TORCH);

        // write object index (we do NOT reuse objects here, we assume that we only write one
        // object per file)
        writeInt(objectIndex++);

        // version string
        writeString("V 1");

        // class name
        writeString(storageTypeName);

        // size
        writeLong(data.size());

        // the actual content
        for (unsigned i = 0 ; i < data.size(); ++i)
            writeType<DataType>(data[i]);

    }

    template<typename DataType>
    void writeTypeTensor(const std::vector<unsigned> &sizes, const std::vector<DataType> &data);

    inline int getNextObjectIndex()
    {
      return objectIndex++;
    }

    template<typename DataType>
    void writeTypeVector(const std::vector<DataType> &data)
    {
      std::vector<unsigned> sizes;
      sizes.push_back(data.size());

      writeTypeTensor<DataType>(sizes, data);
    }



  };

  // template specializations
  template<>
  inline void TorchWriter::writeTypeTensor<float>(const std::vector<unsigned> &sizes, const std::vector<float> &data) 
  {
    writeTypeTensorHelper<float>(sizes, data, "torch.FloatTensor", "torch.FloatStorage");
  }

  template<>
  inline void TorchWriter::writeTypeTensor<int32_t>(const std::vector<unsigned> &sizes, const std::vector<int32_t> &data) 
  {
    writeTypeTensorHelper<int32_t>(sizes, data, "torch.IntTensor", "torch.IntStorage");
  }

  template<>
  inline void TorchWriter::writeTypeTensor<long long>(const std::vector<unsigned> &sizes, const std::vector<long long> &data) 
  {
    writeTypeTensorHelper<long long>(sizes, data, "torch.LongTensor", "torch.LongStorage");
  }



} // namespace flashgg

#endif
