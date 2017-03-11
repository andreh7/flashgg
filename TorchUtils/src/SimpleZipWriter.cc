#include "flashgg/TorchUtils/interface/SimpleZipWriter.h"

#include <boost/crc.hpp>

using namespace std;

//----------------------------------------------------------------------

SimpleZipWriter::SimpleZipWriter(std::ostream &zip_) : 
  zip(zip_)
{
}

//----------------------------------------------------------------------

void SimpleZipWriter::writeZipEntry(const ZipEntryData &entry, bool local) 
{
  string signature;
  if (local)
    signature = "\x50\x4b\x03\x04";
  else
    // central directory record
    signature = "\x50\x4b\x01\x02";
  
  const unsigned short createVersion = 0x314;
  const unsigned short minVersion = 20;
  const unsigned short flags = 0;
  const unsigned short compression = 0; // no compression for the moment (easier to write)
  const unsigned short modTime = 0;
  const unsigned short modDate = 0;
  
  /*
    see also https://users.cs.jmu.edu/buchhofp/forensics/formats/pkzip.html

    local header example:

00000000  50 4b 03 04                                       |PK..               signature
                      14 00                                      ..          |  version
                            00 00                                  ..        |  flags
                                   00 00                             ..      |  compression
                                         aa 85                         ..    |  mod time
                                               43 4a                     CJ..|  mod date
                                                     fd 83                 ..|  CRC32
00000010  af 47                                             |.G              |  CRC32
                70 00 00 00                                 |  p...          |  compressed size
                            70 00  00 00 08 00 00 00 64 61  |      p...      |  uncompressed size
                                         08 00              |          ..    |  file name length
                                               00 00        |            ..  |  extra length
                                                     64 61  |              da|  file name
00000020  74 61 2e 6e 70 79                                 |ta.npy          |  file name (continued)
00000020                    93 4e  55 4d 50 59 01 00 46 00  |      .NUMPY..F.|  file content


    central directory record example:
    00000090                    50 4b  01 02                    |      PK..      | signature
                                             14 03              |          ..    | version
                                                   14 00        |            ....| version needed
                                                         00 00  |              ..| flags
    000000a0  00 00                                             |..              | compression
                    aa 85                                       |  ..            | mod time
                          43 4a                                 |    CJ...Gp...p.| mod date
                                fd 83  af 47                    |      ...G      | CRC32
                                             70 00 00 00        |          p...  | compressed size
                                                         70 00  |              p.| uncompressed size
    000000b0  00 00                                             |..              | uncompressed size
                    08 00                                       |  ..            | file name length
                          00 00                                 |    ..          | extra field length
                                00 00                           |      ..        | file comment length
                                       00 00                    |        ..      | number of disk
                                             00 00              |          ..    | internal attributes
                                                   00 00 80 81  |          ......| external attributes 
    000000c0  00 00 00 00                                       |....            | offset of local header
                          64 61 74 61  2e 6e 70 79              |    data.npy    | file name


  */
  
  zip.write(signature.c_str(), signature.size());

  if (!local)
    zip.write((char*)&createVersion, 2);

  zip.write((char*)&minVersion, 2);
  zip.write((char*)&flags, 2);
  zip.write((char*)&compression, 2);
  zip.write((char*)&modTime, 2);
  zip.write((char*)&modDate, 2);
  
  // calculate CRC32
  zip.write((const char*)&entry.crc32, 4);

  zip.write((const char*)&entry.dataLen, 4); // compressed size
  zip.write((const char*)&entry.dataLen, 4); // uncompressed size
    
  unsigned fnameLength = entry.fname.size();
  zip.write((char*)&fnameLength, 2);
    
  unsigned short extraSize = 0;
  zip.write((char*)&extraSize, 2);
    
  if (! local) {
    unsigned short fileCommentLength = 0;
    zip.write((char*)&fileCommentLength, 2);
    
    unsigned short diskNumber = 0;
    zip.write((char*)&diskNumber, 2);
    
    unsigned short internalAttrs = 0;
    zip.write((char*)&internalAttrs, 2);
    
    unsigned externalAttrs = 0x81800000;
    zip.write((char*)&externalAttrs, 4);
    
    // offset of local header
    zip.write((const char*)&entry.headerOffset, 4);

  } // if not local

  zip.write(entry.fname.c_str(), entry.fname.size());

  if (local)
    currentPos += 30 + entry.fname.size();
  else
    currentPos += 46 + entry.fname.size();
}

//----------------------------------------------------------------------

void SimpleZipWriter::addFile(const string &fname, 
			      const string &data)
{

  unsigned dataLen = data.size();
  
  ZipEntryData entry;
  entry.fname = fname;
  entry.dataLen = dataLen;
  entry.headerOffset = currentPos;
  
  boost::crc_32_type calculator;
  calculator.process_bytes(data.data(), dataLen);
  entry.crc32 = calculator.checksum();
  
  writeZipEntry(entry, true);
  
  // write the local header
  // write the actual data
  zip.write(data.data(), dataLen);
  
  entries.push_back(entry);
  currentPos += dataLen;
}

//----------------------------------------------------------------------

void SimpleZipWriter::writeDirectory() 
{

  unsigned centralDirStart = currentPos;
  
  // write per file directory records
  for (unsigned index = 0; index < entries.size(); ++index) {
    writeZipEntry(entries[index], false);
  }
  
  unsigned centralDirSize = currentPos - centralDirStart;
  
  /*
    
    // end of central directory record
                                                   50 4b 05 06  |            PK..| signature
    000000d0  00 00                                             |..              | disk number
                    00 00                                       |  ..            | disk # w/cd
                          01 00                                 |    ..          | number of files in this archive (?)
                                01 00                           |      ..        | number of files in this archive (?)
                                       36 00 00 00              |        6...    | size of central directory in bytes
                                                   96 00 00 00  |            ....| offset of start of central directory
    000000e0  00 00                                             |..              | comment length
    */

  const string signature = "\x50\x4b\x05\x06";
  zip.write(signature.c_str(), signature.size());
  
  unsigned short diskNumber = 0;
  zip.write((char*)&diskNumber, 2);
  zip.write((char*)&diskNumber, 2);
  
  unsigned short numFiles = entries.size();
  zip.write((char*)&numFiles, 2);
  zip.write((char*)&numFiles, 2);
  
  zip.write((char*)&centralDirSize, 4);
  zip.write((char*)&centralDirStart, 4);
  
  unsigned short commentLength = 0;
  zip.write((char*)&commentLength, 2);
  
  currentPos += 22;
}

//----------------------------------------------------------------------

SimpleZipWriter::~SimpleZipWriter() {
  writeDirectory();
}

//----------------------------------------------------------------------
