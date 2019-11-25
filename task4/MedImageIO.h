// filename: MedImageIO.h
// author: Robert Hess
// version: 1.1
// date: February 19th 2019
// description: File in- and output for medical images with class cMedImage.
//              The images are stored with index (row,col) with (0,0) top left.

#pragma once
#ifndef HEADER_MEDICAL_IMAGE_IO_INCLUDED
#define HEADER_MEDICAL_IMAGE_IO_INCLUDED

//#define IMPLEMENT_PNG

#include "MedImage.h"
#include <fstream>
#ifdef IMPLEMENT_PNG
#include "png.h"
#endif

//=============================================================================
template<class T> class cMedImageIO {
private:
    // some elementary types
    typedef uint8_t BYTE;
    typedef uint16_t WORD;
    typedef uint32_t DWORD;
    typedef int32_t LONG;

// align to 2 byte blocks
#pragma pack(push, 2)

    // the bitmap file header
    typedef struct sFileHeader {
        WORD type;				// file type (always the letters BM)
        DWORD fileSize;			// size of file including headers and LUT
        DWORD reserved;			// always zero
        DWORD offBytes;			// position of pixel data within file
    } tFileHeader;

    // the bitmap info header
    typedef struct sInfoHeader {
        DWORD infoSize;			// info header size (always 40)
        LONG width;				// number of columns in image
        LONG height;			// number of rows in image
        WORD planes;			// number of image planes (here always 1)
        WORD bitCount;			// bits per pixel
        DWORD compression;		// type of compression (here always 0 for no compression)
        DWORD imageSize;		// total size of pixel data (0 = automatic)
        LONG xPixelsPerMeter;	// x-resolution in pixels per meter
        LONG yPixelsPerMeter;	// y-resolution in pixels per meter
        DWORD colorUsed;		// number of colours in LUT
        DWORD colorImportant;	// number of important colours in LUT
    } tInfoHeader;

    // RGB pixel in bitmap
    typedef struct {
        BYTE b;
        BYTE g;
        BYTE r;
        BYTE reserve;
    } tPixel;

// align elements as before
#pragma pack(pop)

private:
    // some support methods
    static void convertToLittleEndian(char *data, unsigned size);
    static void writeUS(unsigned short group, unsigned short element, unsigned short value, std::ofstream &out);
    static void writeDS(unsigned short group, unsigned short element, double value, std::ofstream &out);

    // elementary methods for DICOM and RAW-files
    static void writeDicom(const cMedImage<T> &image, std::ofstream &out, double tubeVoltage, double exposure);
    static void readRaw(cMedImage<T> &image, std::ifstream &inp, unsigned header, unsigned width, unsigned height, int pixelType,
                        bool isUnsigned, bool isBigEndian);

    // elementary methods for BMP-files
    static void writeBmp(const cMedImage<T> &image, std::ofstream &out);
    static void readBmp(cMedImage<T> &image, std::ifstream &inp);

    // elementary methods for PNG-files
#ifdef IMPLEMENT_PNG
    static void writePng(const cMedImage<T> &image, FILE *out);
	static void readPng(cMedImage<T> &image, FILE *inp);
#endif

public:
    // interfaces for DICOM and RAW-files
    static void writeDicom(const cMedImage<T> &image, const std::string &filename, double tubeVoltage, double exposure);
    static void writeDicom(const cMedImage<T> &image, const std::wstring &filename, double tubeVoltage, double exposure);
    static void readRaw(cMedImage<T> &image, const std::string &filename, unsigned header, unsigned width, unsigned height, int pixelType,
                        bool isUnsigned, bool isBigEndian);
    static void readRaw(cMedImage<T> &image, const std::wstring &filename, unsigned header, unsigned width, unsigned height, int pixelType,
                        bool isUnsigned, bool isBigEndian);

    // interfaces for BMP-files
    static void writeBmp(const cMedImage<T> &image, const std::string &filename);
    static void writeBmp(const cMedImage<T> &image, const std::wstring &filename);
    static void readBmp(cMedImage<T> &image, const std::string &filename);
    static void readBmp(cMedImage<T> &image, const std::wstring &filename);

    // interfaces for PNG-files
#ifdef IMPLEMENT_PNG
    static void writePng(const cMedImage<T> &image, const std::string &filename);
	static void writePng(const cMedImage<T> &image, const std::wstring &filename);
	static void readPng(cMedImage<T> &image, const std::string &filename);
	static void readPng(cMedImage<T> &image, const std::wstring &filename);
#endif
};


//=============================================================================
template<class T> void cMedImageIO<T>::cMedImageIO::convertToLittleEndian(
        // converts big endian to little endian and vice versa
        char *data,		// [inout]pointer to data to be swaped
        unsigned size)	// [in] size of data to be swaped
{
    unsigned i;	// local index variable
    char tmp;	// tamporary value to swap data

    // swap bytes
    for (i = 0; i < size / 2; i++) {
        tmp = data[size - i - 1];
        data[size - i - 1] = data[i];
        data[i] = tmp;
    }
}


//=============================================================================
template<class T> void cMedImageIO<T>::cMedImageIO::writeUS(
        unsigned short group,	// [in] group if tag
        unsigned short element,	// [in] element of tag
        unsigned short value,	// [in] value of tag
        std::ofstream &out)		// [inout] output stream
{
    unsigned short size = 2;	// size of tag

    out.write((char*)&group, sizeof(group));
    out.write((char*)&element, sizeof(element));
    out.write("US", 2);
    out.write((char*)&size, sizeof(size));
    out.write((char*)&value, sizeof(value));
}


//=============================================================================
template<class T> void cMedImageIO<T>::cMedImageIO::writeDS(
        unsigned short group,	// [in] group if tag
        unsigned short element,	// [in] element of tag
        double value,			// [in] value of tag
        std::ofstream &out)		// [inout] output stream
{
    unsigned short size;	// size of tag
    char text[99];			// number as string

    sprintf(text, "%lg", value);
    size = (unsigned short)strlen(text);

    out.write((char*)&group, sizeof(group));
    out.write((char*)&element, sizeof(element));
    out.write("DS", 2);
    out.write((char*)&size, sizeof(size));
    out.write(text, size);
}


//=============================================================================
template<class T> void cMedImageIO<T>::cMedImageIO::writeDicom(
        const cMedImage<T> &image,	// [in] image to be written in file
        std::ofstream &out,			// [inout] stream to be written into
        double tubeVoltage,			// [in] tube voltage of image
        double exposure)			// [in] current time product of image
{
    int i;
    unsigned r, c;
    int a = 0;
    bool explicitVR = true;
    T pvBlack;				// pixel value for black
    T pvWhite;				// pixel value for white

    // check if file is open
    if (!out.is_open())
        throw std::runtime_error("Could not open DICOM file for writing.");

    // leading blanks for explicit VR
    for (i = 0; explicitVR && i<32; i++)
        out.write((char*)&a, 4);

    // write "DICM"
    out.write("DICM", 4);

    // write x-ray data
    if (tubeVoltage>0) writeDS(0x18, 0x60, tubeVoltage, out);	// tube voltage
    if (exposure>0) writeDS(0x18, 0x1152, exposure, out);		// exposure

    // write size
    writeUS(0x28, 0x10, (unsigned short)image.nRow, out);	// rows
    writeUS(0x28, 0x11, (unsigned short)image.nCol, out);	// columns

    // image information
    pvBlack = image.pixelValueForBlack;
    pvWhite = image.pixelValueForWhite;
    writeUS(0x28, 0x100, 16, out);							// bits per pixel
    writeUS(0x28, 0x101, 16, out);							// valid bits per pixel
    writeDS(0x28, 0x1050, 0.5 * (pvBlack + pvWhite), out);	// window center
    writeDS(0x28, 0x1051, (double)(pvWhite - pvBlack), out);// window width
    //writeDS(0x28, 0x1052, -(1 << 15), out);				// rescale intercept

    // write pixel data
    out.write((char*)&(a = 0x7fe0), 2);
    out.write((char*)&(a = 0x10), 2);
    out.write("OW", 2);
    out.write((char*)&(a = 0), 2);
    out.write((char*)&(a = image.nRow * image.nCol), 4);
    for (r = 0; r < image.nRow; r++) {
        for (c = 0; c < image.nCol; c++) {
            out.write((char*)&(a = (int)(image[r][c]/* + (1 << 15)*/)), 2);
        }
    }
}


//=============================================================================
template<class T> void cMedImageIO<T>::cMedImageIO::readRaw(
        cMedImage<T> &image,// [out] image read from file
        std::ifstream &inp,	// [inout] input file stream
        unsigned header,	// [in] header size
        unsigned width,		// [in] number of pixels per row in image
        unsigned height,	// [in] numer of rows in image
        int pixelType,		// [in] pixel type (0: 1-byte int, 1: 2-byte int,
        //   2: 4-byte int, 3: 4-byte float, 4: 8-byte float)
        bool isUnsigned,	// [in] true if unsigned pixel type
        bool isBigEndian)	// [in] true if big endian (high byte first)
{
    unsigned r, c;		// index for row and column
    unsigned char u1;	// pixel of type 1-byte unsigned
    unsigned short u2;	// pixel of type 2-byte unsigned
    unsigned u4;		// pixel of type 4-byte unsigned
    float f4;			// pixel of type 4-byte float
    double f8;			// pixel of type 8-byte float

    // check if file is open
    if (!inp.is_open())
        throw std::runtime_error("Could not open RAW file for reading.");

    try {
        // set size of image
        image.create(height, width);

        // skip header
        inp.seekg(header, std::ios_base::cur);

        // read image
        for (r = 0; r < height; r++) {
            for (c = 0; c < width; c++) {
                switch (pixelType) {
                    case 0:	// 1 byte integer
                        // read value
                        inp.read((char*)&u1, sizeof(u1));
                        // store value
                        if (isUnsigned) f8 = u1;
                        else f8 = (char)u1;
                        break;
                    case 1:	// 2 byte integer
                        // read value
                        inp.read((char*)&u2, sizeof(u2));
                        // if big endian swap bytes
                        if (isBigEndian) convertToLittleEndian((char*)&u2, 2);
                        // store value
                        if (isUnsigned) f8 = u2;
                        else f8 = (short)u2;
                        break;
                    case 2:	// 4 byte integer
                        // read value
                        inp.read((char*)&u4, sizeof(u4));
                        // if big endian swap bytes
                        if (isBigEndian) convertToLittleEndian((char*)&u4, 4);
                        // store value
                        if (isUnsigned) f8 = u4;
                        else f8 = (int)u4;
                        break;
                    case 3:	// 4 byte float
                        // read value
                        inp.read((char*)&f4, sizeof(f4));
                        // if big endian swap bytes
                        if (isBigEndian) convertToLittleEndian((char*)&f4, 4);
                        // store value
                        f8 = f4;
                        break;
                    case 4:	// 8 byte float
                        // read and store value
                        inp.read((char*)&f8, sizeof(f8));
                        // if big endian swap bytes
                        if (isBigEndian) convertToLittleEndian((char*)&f8, 8);
                        break;
                }
                // store pixel in vector
                image[r][c] = (T)f8;
            }
        }

        // find minimum and maximum pixel value
        image.findMinMax(image.pixelValueForWhite, image.pixelValueForBlack);
    }
    catch (...) {
        image.destroy();
        throw std::runtime_error("Could not read raw image file.");
    }
}


//=============================================================================
template<class T> void cMedImageIO<T>::writeBmp(
        const cMedImage<T> &image,	// [in] image to be written in file
        std::ofstream &out)			// [inout] stream to be written into
{
    tFileHeader fileHeader;	// bitmap file header
    tInfoHeader infoHeader;	// bitmap info header
    unsigned r, c;			// index for row and column
    DWORD colour;			// actual colour for look up table
    T dbl;					// actual pixel value
    BYTE byte;				// actual pixel
    T pvBlack;				// pixel value for black
    T pvWhite;				// pixel value for white

    // check if file is open
    if (!out.is_open())
        throw std::runtime_error("Could not open bitmap file for writing.");

    try {
        // prepare and write file header
        fileHeader.type = 'B' + 256 * 'M';
        fileHeader.fileSize = image.nCol + (image.nCol % 4 ? 4 - image.nCol % 4 : 0);
        fileHeader.fileSize = 14 + 40 + 1024 + image.nRow*fileHeader.fileSize;
        fileHeader.reserved = 0;
        fileHeader.offBytes = 54 + 1024;
        out.write((char*)&fileHeader, sizeof(fileHeader));

        // prepare and write info header
        infoHeader.infoSize = 40;			// Größe des info headers
        infoHeader.width = image.nCol;		// Breite des Bilds in Pixeln
        infoHeader.height = image.nRow;		// Höhe des Bilds in Pixeln
        infoHeader.height *= -1;			// Bild vertikal spiegeln
        infoHeader.planes = 1;				// Anzahl der Bildebenen(immer 1)
        infoHeader.bitCount = 8;			// Bits pro Pixel(hier jeweils 8 Bit für die Farben)
        infoHeader.compression = 0;			// Art der Kompression(0 für keine Kompression)
        infoHeader.imageSize = 0;			// Größe der Pixeldaten in Byte(0 = automatisch)
        infoHeader.xPixelsPerMeter = 7000;	// Pixel/Meter in x-Richtung (hier 1 Pixel = 143 µm)
        infoHeader.yPixelsPerMeter = 7000;	// Pixel/Meter in y-Richtung (hier 1 Pixel = 143 µm)
        infoHeader.colorUsed = 256;			// Anzahl der verwendeten Farben in LUT
        infoHeader.colorImportant = 256;	// Anzahl der wichtigen Farben in LUT
        out.write((char*)&infoHeader, sizeof(infoHeader));

        // write look up table
        for (c = 0; c < 256; c++) {
            colour = c | c << 8 | c << 16;
            out.write((char*)&colour, sizeof(colour));
        }

        // write pixel data
        colour = 0;
        pvBlack = image.pixelValueForBlack;
        pvWhite = image.pixelValueForWhite;
        for (r = 0; r < image.nRow; r++) {
            for (c = 0; c < image.nCol; c++) {
                dbl = image[r][c];
                if (pvBlack <= pvWhite && dbl <= pvBlack) byte = 0;
                else if (pvWhite < pvBlack && dbl <= pvWhite) byte = 255;
                else if (pvBlack <= pvWhite && dbl >= pvWhite) byte = 255;
                else if (pvWhite < pvBlack && dbl >= pvBlack) byte = 0;
                else byte = (BYTE)(256 * (dbl - pvBlack) / (pvWhite - pvBlack));
                out.write((char*)&byte, sizeof(byte));
            }
            if (image.nCol % 4) out.write((char*)&colour, 4 - image.nCol % 4);
        }
    }
    catch (...) {
        throw std::runtime_error("Could not write bitmap file (BMP).");
    }
}


//=============================================================================
template<class T> void cMedImageIO<T>::cMedImageIO::readBmp(
        cMedImage<T> &image,	// [out] image read from file
        std::ifstream &inp)		// [inout] stream to read image from
{
    tFileHeader fileHeader;	// file header of bitmap
    tInfoHeader infoHeader;	// info header of bitmap
    unsigned r, c;			// index for row and column in bitmap
    tPixel pixel;			// actual pixel value
    tPixel *lut = NULL;		// look up table for some bitmap types
    tPixel *pPixel;			// pointer to actual colour in LUT
    BYTE index = 0;			// actual index for some bitmap types
    unsigned short us;		// actual pixel in 16 bit bitmap type
    T dbl;					// actual pixel value in image

    // check if file is open
    if (!inp.is_open())
        throw std::runtime_error("Could not open bitmap file for reading.");

    try {
        // read file- and info-header
        inp.read((char*)&fileHeader, sizeof(fileHeader));
        inp.read((char*)&infoHeader, sizeof(infoHeader));

        // check info header
        if (infoHeader.planes != 1) throw - 1;
        if (infoHeader.compression != 0) throw - 1;	// no compression supported

        // prepare image
        image.create(abs(infoHeader.height), infoHeader.width);

        switch (infoHeader.bitCount) {
            case 1:
                // read LUT
                lut = new tPixel[2];
                if (!lut) throw - 1;
                inp.read((char*)lut, 4 * (infoHeader.colorUsed>0 ? infoHeader.colorUsed : 2));
                // loop over all rows
                for (r = 0; r < image.nRow; r++) {
                    // loop over all pixels in a row
                    for (c = 0; c < image.nCol; c++) {
                        if (c % 8 == 0) inp.read((char*)&index, 1);
                        pPixel = lut + ((index >> (7 - c % 8)) & 1);
                        image[infoHeader.height > 0 ? image.nRow - r - 1 : r][c] = (T)((0.0 + pPixel->r + pPixel->g + pPixel->b) / 3.0);
                    }
                    // complement to multiples of 4 byte
                    if (((image.nCol + 7) / 8) % 4) inp.read((char*)&pixel, 4 - ((image.nCol + 7) / 8) % 4);
                }
                // free LUT
                delete[] lut;
                break;
            case 4:
                // read LUT
                lut = new tPixel[16];
                if (!lut) throw - 1;
                inp.read((char*)lut, 4 * (infoHeader.colorUsed>0 ? infoHeader.colorUsed : 16));
                // loop over all rows
                for (r = 0; r < image.nRow; r++) {
                    // loop over all pixels in a row
                    for (c = 0; c < image.nCol; c++) {
                        if (c % 2 == 0) inp.read((char*)&index, 1);
                        pPixel = lut + ((index >> (4 * (1 - c % 2))) & 15);
                        image[infoHeader.height > 0 ? image.nRow - r - 1 : r][c] = (T)((0.0 + pPixel->r + pPixel->g + pPixel->b) / 3.0);
                    }
                    // complement to multiples of 4 byte
                    if (((image.nCol + 1) / 2) % 4) inp.read((char*)&pixel, 4 - ((image.nCol + 1) / 2) % 4);
                }
                // free LUT
                delete[] lut;
                break;
            case 8:
                // read LUT
                lut = new tPixel[256];
                if (!lut) throw - 1;
                inp.read((char*)lut, 4 * (infoHeader.colorUsed>0 ? infoHeader.colorUsed : 256));
                // loop over all rows
                for (r = 0; r < image.nRow; r++) {
                    // loop over all pixels in a row
                    for (c = 0; c < image.nCol; c++) {
                        inp.read((char*)&index, 1);
                        image[infoHeader.height > 0 ? image.nRow - r - 1 : r][c] = (T)((0.0 + lut[index].r + lut[index].g + lut[index].b) / 3.0);
                    }
                    // complement to multiples of 4 byte
                    if (image.nCol % 4) inp.read((char*)&pixel, 4 - image.nCol % 4);
                }
                // free LUT
                delete[] lut;
                break;
            case 16:
                // skip LUT
                inp.seekg(4 * infoHeader.colorUsed, std::ios_base::cur);
                // loop over all rows
                for (r = 0; r < image.nRow; r++) {
                    // loop over all pixels in a row
                    for (c = 0; c < image.nCol; c++) {
                        inp.read((char*)&us, 2);
                        dbl = (T)(((us >> 7 & 0xf8) + (us >> 2 & 0xf8) + (us << 3 & 0xf8)) / 3.0);
                        image[infoHeader.height > 0 ? image.nRow - r - 1 : r][c] = dbl;
                    }
                    // complement to multiples of 4 byte
                    if (image.nCol % 2) inp.read((char*)&us, 2);
                }
                break;
            case 24:
                // skip LUT
                inp.seekg(4 * infoHeader.colorUsed, std::ios_base::cur);
                // loop over all rows
                for (r = 0; r < image.nRow; r++) {
                    // loop over all pixels in a row
                    for (c = 0; c < image.nCol; c++) {
                        inp.read((char*)&pixel, 3);
                        image[infoHeader.height > 0 ? image.nRow - r - 1 : r][c] = (T)((0.0 + pixel.r + pixel.g + pixel.b) / 3.0);
                    }
                    // complement to multiples of 4 byte
                    if (image.nCol % 4) inp.read((char*)&pixel, image.nCol % 4);
                }
                break;
            case 32:
                // skip LUT
                inp.seekg(4 * infoHeader.colorUsed, std::ios_base::cur);
                // loop over all rows
                for (r = 0; r < image.nRow; r++) {
                    // loop over all pixels in a row
                    for (c = 0; c < image.nCol; c++) {
                        inp.read((char*)&pixel, 4);
                        image[infoHeader.height > 0 ? image.nRow - r - 1 : r][c] = (T)((0.0 + pixel.r + pixel.g + pixel.b) / 3.0);
                    }
                }
                break;
            default:
                throw - 1;
        }
        image.pixelValueForBlack = (T)255.999;
        image.pixelValueForWhite = 0;
    }
    catch (...) {
        // free mamory
        if (lut) delete[] lut;
        image.destroy();

        // throw error
        throw std::runtime_error("Could not read bitmap file (BMP).");
    }
}

#ifdef IMPLEMENT_PNG

//=============================================================================
template<class T> void cMedImageIO<T>::cMedImageIO::writePng(
	const cMedImage<T> &image,	// [in] image to be written in file
	FILE *out)					// [inout] file pointer to be written to
{
	png_structp png_ptr = NULL;	// PNG structure
	png_infop info_ptr = NULL;	// PNG info structure
	unsigned short *row = NULL;	// pointer to actual row in PNG
	unsigned r, c;				// index for row and column
	T dbl;						// current pixel value
	T pvBlack;					// pixel value for black
	T pvWhite;					// pixel value for white

	// check if file is open
	if (!out)
		throw std::runtime_error("Could not open PNG file for writing.");

	try {
		// get memory for row
		if ((row = (unsigned short*)malloc(image.nCol*sizeof(unsigned short))) == NULL) throw - 1;

		// create structure to write PNG file
		png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING,
			(png_voidp)NULL/*user_error_ptr*/, NULL/*user_error_fn*/,
			NULL/*user_warning_fn*/);
		if (png_ptr == NULL) throw - 1;
		info_ptr = png_create_info_struct(png_ptr);
		if (info_ptr == NULL) throw - 1;

		// define long jump destination which is called after error handling
		if (setjmp(png_jmpbuf(png_ptr))) throw - 1;

		// give PNG file pointer
		png_init_io(png_ptr, out);

		// fill png_info
		png_set_IHDR(png_ptr, info_ptr, image.nCol, image.nRow,
			16/*bit_depth*/, PNG_COLOR_TYPE_GRAY/*color_type*/,
			PNG_INTERLACE_NONE/*interlace_type*/,
			PNG_COMPRESSION_TYPE_DEFAULT/*compression_type*/,
			PNG_FILTER_TYPE_DEFAULT/*filter_method*/);

		// gamma value at which the image was created
		png_set_gAMA(png_ptr, info_ptr, 1.0);

		// write file heading
		png_write_info(png_ptr, info_ptr);
		// write pixel data
		png_set_swap(png_ptr);
		pvBlack = image.pixelValueForBlack;
		pvWhite = image.pixelValueForWhite;
		for (r = 0; r < image.nRow; r++) {
			for (c = 0; c < image.nCol; c++) {
				dbl = image[r][c];
				if (pvBlack <= pvWhite && dbl <= pvBlack) row[c] = 0;
				else if (pvWhite < pvBlack && dbl <= pvWhite) row[c] = 65535;
				else if (pvBlack <= pvWhite && dbl >= pvWhite) row[c] = 65535;
				else if (pvWhite < pvBlack && dbl >= pvBlack) row[c] = 0;
				else row[c] = (unsigned short)(65536 * (dbl - pvBlack) / (pvWhite - pvBlack));
			}
			png_write_row(png_ptr, (png_bytep)row);
		}
		// finish file
		png_write_end(png_ptr, info_ptr);

		// free memory
		free(row);
		row = NULL;
		png_destroy_write_struct(&png_ptr, &info_ptr);
		png_ptr = NULL;
		info_ptr = NULL;
	}
	catch (...) {
		// free memory
		if (row) free(row);
		row = NULL;
		if (png_ptr) {
			if (info_ptr) png_destroy_write_struct(&png_ptr, &info_ptr);
			else png_destroy_write_struct(&png_ptr, (png_infopp)NULL);
		}
		png_ptr = NULL;
		info_ptr = NULL;

		// trow error
		throw std::runtime_error("Could not write protable network graphics file (PNG).");
	}
}


//=============================================================================
template<class T> void cMedImageIO<T>::cMedImageIO::readPng(
	cMedImage<T> &image,	// [out] image read from file
	FILE *inp)				// [inout] file pointer to read image
{
	png_image img;				// The control structure used by libpng
	png_byte header[8];			// header to check for valid PNG
	png_bytep buffer = NULL;	// read pixel data
	unsigned short *pData;		// pointer to actual pixel in read data
	unsigned r, c;				// index for row and column
	png_structp png_ptr = NULL;
	png_infop info_ptr = NULL;
	png_infop end_info = NULL;

	// check if file is open
	if (!inp)
		throw std::runtime_error("Could not open PNG file for reading.");

	try {
		// check for valid PNG
		if (fread(header, 8, 1, inp) != 1) throw -1;
		if (png_sig_cmp(header, 0, 8)) throw -1;
		fseek(inp, 0, SEEK_SET);

		// initialize ...
		png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING, (png_voidp)NULL, NULL, NULL);
		if (!png_ptr) throw - 1;
		info_ptr = png_create_info_struct(png_ptr);
		if (!info_ptr) {
			png_destroy_read_struct(&png_ptr, (png_infopp)NULL, (png_infopp)NULL);
			throw -1;
		}
		end_info = png_create_info_struct(png_ptr);
		if (!end_info) {
			png_destroy_read_struct(&png_ptr, &info_ptr,
				(png_infopp)NULL);
			throw -1;
		}

		// Initialize the 'png_image' structure.
		memset(&img, 0, (sizeof img));
		img.version = PNG_IMAGE_VERSION;

		// read image header
		if (png_image_begin_read_from_stdio(&img, inp) == 0) throw - 1;

		// define pixel type as 16 bit gray scale
		img.format = PNG_FORMAT_GRAY | PNG_FORMAT_FLAG_LINEAR;

		// get memory for image
		if ((buffer = (png_bytep)malloc(PNG_IMAGE_SIZE(img))) == NULL) throw - 1;

		// read content of PNG image
		if (png_image_finish_read(&img, NULL/*background*/, buffer,
			0/*row_stride*/, NULL/*colormap*/) == 0) throw - 1;

		// prepare medical image
		image.create(img.height, img.width);

		// copy image
		pData = (unsigned short*)buffer;
		for (r = 0; r < image.nRow; r++)
			for (c = 0; c < image.nCol; c++)
				image[r][c] = *(pData++);

		// find minimum and maximum pixel value
		image.findMinMax(image.pixelValueForWhite, image.pixelValueForBlack);
	}
	catch (...) {
		// clear up

		// throw with error
		throw std::runtime_error("Could not read protable network graphics file (PNG).");
	}
}

#endif


//=============================================================================
template<class T> void cMedImageIO<T>::cMedImageIO::writeDicom(
        const cMedImage<T> &image,		// [in] image to be written in file
        const std::string &filename,	// [in] name of file to write image to
        double tubeVoltage,				// [in] tube voltage of image
        double exposure)				// [in] current time product of image
{
    // open output file
    std::ofstream out(filename, std::ios::binary);

    // write bitmap
    writeDicom(image, out, tubeVoltage, exposure);
}


//=============================================================================
template<class T> void cMedImageIO<T>::cMedImageIO::writeDicom(
        const cMedImage<T> &image,		// [in] image to be written in file
        const std::wstring &filename,	// [in] name of file to write image to
        double tubeVoltage,				// [in] tube voltage of image
        double exposure)				// [in] current time product of image
{
    std::string tmp(filename.begin(), filename.end());	// filename with chars

    // open output file
    std::ofstream out(tmp, std::ios::binary);

    // write bitmap
    writeDicom(image, out, tubeVoltage, exposure);
}


//=============================================================================
template<class T> void cMedImageIO<T>::cMedImageIO::readRaw(
        cMedImage<T> &image,		// [out] image read from file
        const std::string &filename,// [in] name of file to read image from
        unsigned header,			// [in] header size
        unsigned width,				// [in] number of pixels per row in image
        unsigned height,			// [in] numer of rows in image
        int pixelType,				// [in] pixel type (0: 1-byte int, 1: 2-byte int,
        //   2: 4-byte int, 3: 4-byte float, 4: 8-byte float)
        bool isUnsigned,			// [in] true if unsigned pixel type
        bool isBigEndian)			// [in] true if big endian (high byte first)
{
    // open input file
    std::ifstream inp(filename, std::ios::binary);

    // read image
    readRaw(image, inp, header, width, height, pixelType, isUnsigned, isBigEndian);
}


//=============================================================================
template<class T> void cMedImageIO<T>::cMedImageIO::readRaw(
        cMedImage<T> &image,			// [out] image read from file
        const std::wstring &filename,	// [in] name of file to read image from
        unsigned header,				// [in] header size
        unsigned width,					// [in] number of pixels per row in image
        unsigned height,				// [in] numer of rows in image
        int pixelType,					// [in] pixel type (0: 1-byte int, 1: 2-byte int,
        //   2: 4-byte int, 3: 4-byte float, 4: 8-byte float)
        bool isUnsigned,				// [in] true if unsigned pixel type
        bool isBigEndian)				// [in] true if big endian (high byte first)
{
    std::string tmp(filename.begin(), filename.end());	// filename with chars

    // open input file
    std::ifstream inp(tmp, std::ios::binary);

    // read image
    readRaw(image, inp, header, width, height, pixelType, isUnsigned, isBigEndian);
}


//=============================================================================
template<class T> void cMedImageIO<T>::cMedImageIO::writeBmp(
        const cMedImage<T> &image,	// [in] image to be written in file
        const std::string &filename)		// [in] name of file to write image to
{
    // open output file
    std::ofstream out(filename, std::ios::binary);

    // write bitmap
    writeBmp(image, out);
}


//=============================================================================
template<class T> void cMedImageIO<T>::cMedImageIO::writeBmp(
        const cMedImage<T> &image,	// [in] image to be written in file
        const std::wstring &filename)	// [in] name of file to write image to
{
    std::string tmp(filename.begin(), filename.end());	// filename with chars

    // open output file
    std::ofstream out(tmp, std::ios::binary);

    // write bitmap
    writeBmp(image, out);
}


//=============================================================================
template<class T> void cMedImageIO<T>::cMedImageIO::readBmp(
        cMedImage<T> &image,			// [out] image read from file
        const std::string &filename)	// [in] name of file to read image from
{
    // open input file
    std::ifstream inp(filename, std::ios::binary);

    // read image
    readBmp(image, inp);
}

//=============================================================================
template<class T> void cMedImageIO<T>::cMedImageIO::readBmp(
        cMedImage<T> &image,			// [out] image read from file
        const std::wstring &filename)	// [in] name of file to read image from
{
    std::string tmp(filename.begin(), filename.end());	// filename with chars

    // open input file
    std::ifstream inp(tmp, std::ios::binary);

    // read image
    readBmp(image, inp);
}

#ifdef IMPLEMENT_PNG
//=============================================================================
template<class T> void cMedImageIO<T>::cMedImageIO::writePng(
	const cMedImage<T> &image,		// [in] image to be written in file
	const std::string &filename)	// [in] name of file to write image to
{
	FILE *out = NULL;	// pointer for output file

	try {
		// open file
		if ((out = fopen(filename, "wb")) == NULL)
			throw std::runtime_error("Could not open PNG-file for writing.");

		// write PNG
		writePng(image, out);

		// close file
		fclose(out);

	}
	catch (std::exception &e) {
		if (out) fclose(out);
		throw e;
	}
	catch (char *text) {
		if (out) fclose(out);
		throw text;
	}
	catch (...) {
		if (out) fclose(out);
		throw std::runtime_error("Could not write PNG-file.");
	}
}

//=============================================================================
template<class T> void cMedImageIO<T>::cMedImageIO::writePng(
	const cMedImage<T> &image,		// [in] image to be written in file
	const std::wstring &filename)	// [in] name of file to write image to
{
	FILE *out = NULL;	// pointer for output file

	try {
		// open file
		if ((out = _wfopen(filename.c_str(), L"wb")) == NULL)
			throw std::runtime_error("Could not open PNG-file for writing.");

		// write PNG
		writePng(image, out);

		// close file
		fclose(out);

	}
	catch (std::exception &e) {
		if (out) fclose(out);
		throw e;
	}
	catch (char *text) {
		if (out) fclose(out);
		throw text;
	}
	catch (...) {
		if (out) fclose(out);
		throw std::runtime_error("Could not write PNG-file.");
	}
}


//=============================================================================
template<class T> void cMedImageIO<T>::cMedImageIO::readPng(
	cMedImage<T> &image,			// [out] image read from file
	const std::string &filename)	// [in] name of file to read image from
{
	FILE *inp = NULL;	// pointer to input file structure

	try {
		// open file
		if ((inp = fopen(filename, "rb")) == NULL)
			throw std::runtime_error("Could not open PNG-file for reading.");

		// read image
		readPng(image, inp);

		// close file
		fclose(inp);

	}
	catch (std::exception &e) {
		if (inp) fclose(inp);
		throw e;
	}
	catch (char *text) {
		if (inp) fclose(inp);
		throw text;
	}
	catch (...) {
		if (inp) fclose(inp);
		throw std::runtime_error("Could not read PNG-file.");
	}
}

//=============================================================================
template<class T> void cMedImageIO<T>::cMedImageIO::readPng(
	cMedImage<T> &image,		// [out] image read from file
	const std::wstring &filename)	// [in] name of file to read image from
{
	FILE *inp = NULL;	// pointer to input file structure

	try {
		// open file
		if((inp = _wfopen(filename.c_str(), L"rb"))==NULL)
			throw std::runtime_error("Could not open PNG-file for reading.");

		// read image
		readPng(image, inp);

		// close file
		fclose(inp);
	}
	catch (std::exception &e) {
		if (inp) fclose(inp);
		throw e;
	}
	catch (char *text) {
		if (inp) fclose(inp);
		throw text;
	}
	catch (...) {
		if (inp) fclose(inp);
		throw std::runtime_error("Could not read PNG-file.");
	}
}
#endif

#endif // #ifndef HEADER_MEDICAL_IMAGE_IO_INCLUDED
