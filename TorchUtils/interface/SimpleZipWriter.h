#ifndef TorchUtils_SimpleZipWriter_h
#define TorchUtils_SimpleZipWriter_h

#include <fstream>
#include <string>
#include <vector>

class SimpleZipWriter 
{
 private:
  class ZipEntryData {
  public:
    std::string fname;
  
    unsigned dataLen;

    unsigned crc32;

    // offset of local header
    unsigned headerOffset;
  
  };

  /** for keeping track of added files
      when writing the central directory at the end */
  std::vector<ZipEntryData> entries;

  /** current position within the file */
  unsigned currentPos = 0;
  

  /** the output file */
  std::ostream &zip;
  
 public:  
  SimpleZipWriter(std::ostream &zip_);

  ~SimpleZipWriter();

 protected:
  /** @param local if true, writes a local file header, if false writes a central directory header */
  void writeZipEntry(const ZipEntryData &entry, bool local);
 
 public:
  void addFile(const std::string &fname, const std::string &data);

 protected: 
  /** writes the directory record at the end of the zip file. Automatically called
      in the destructor.*/
  void writeDirectory();


};

#endif
