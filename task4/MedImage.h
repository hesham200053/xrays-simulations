// filename: MedImage.h
// author: Robert Hess
// version: 1.1
// date: March 19th 2019
// description: Template class for medical images. Pixel type as template.
//              Coordinates as (row,col) with (0,0) at top left.

#pragma once
#ifndef HEADER_MEDICAL_IMAGE_INCLUDED
#define HEADER_MEDICAL_IMAGE_INCLUDED

#include <iostream>	// for pointer NULL
#include <cstring>


//=========================================================================
template<class T> class cMedImage
{
public:
    unsigned nRow;			// number of rows, i.e. image height
    unsigned nCol;			// number of columns, i.e. image width
    T **data;				// image data: first index for row, second index for column
    T pixelValueForBlack;	// pixel value to be displayed black
    T pixelValueForWhite;	// pixel value to be displayed white

    // construction and destruction
    cMedImage();
    cMedImage(unsigned rows, unsigned cols);
    cMedImage(const cMedImage &image);
    ~cMedImage();

    T *operator[](unsigned r) const;
    bool imagePresent();
    void create(unsigned rows, unsigned cols);
    void create(unsigned rows, unsigned cols, const T &value);
    void destroy();
    void set(const cMedImage<T> &image);
    void set(const T &value);
    void add(const T &value);
    void add(const cMedImage<T> &image);
    void subtract(const T &value);
    void subtract(const cMedImage<T> &image);
    void multiply(const T &value);
    void multiply(const cMedImage<T> &image);
    cMedImage &operator=(const cMedImage &image);
    cMedImage &operator=(const T &value);
    cMedImage &operator+=(const T &value);
    cMedImage &operator-=(const T &value);
    cMedImage &operator*=(const T &value);
    cMedImage &operator+=(const cMedImage<T> &image);
    cMedImage &operator-=(const cMedImage<T> &image);
    cMedImage &operator*=(const cMedImage<T> &image);
    void findMinMax(T &min, T &max);
    void cut(unsigned r1, unsigned c1, unsigned r2, unsigned c2);
    void shift(int r0, int c0, const T &fill);
    void mirror(bool vertical);
    void mirror(const cMedImage<T> &image, bool vertical);
    void rotate(unsigned dir);
    void rotate(const cMedImage<T> &image, unsigned dir);
    void roll(int hori, int vert);
    void roll(const cMedImage<T> &image, int hori, int vert);
    void autoWindowing(bool rad=true);
};


//=========================================================================
template<class T> cMedImage<T>::cMedImage()
{
    // initialize data
    nRow = nCol = 0;
    data = NULL;
    pixelValueForBlack = (T)0;
    pixelValueForWhite = (T)1;
}


//=========================================================================
template<class T> cMedImage<T>::cMedImage(
        unsigned rows,	// [in] number of rows in image
        unsigned cols)	// [in] number of columns in image
{
    // initialize data and vcreate iimage of given size
    nRow = nCol = 0;
    data = NULL;
    pixelValueForBlack = (T)0;
    pixelValueForWhite = (T)1;
    create(rows, cols);
}


//=========================================================================
template<class T> cMedImage<T>::cMedImage(
        const cMedImage &image)	// [in] image to be copied from
{
    // initialize data and copy given image
    nRow = nCol = 0;
    data = NULL;
    set(image);
}


//=========================================================================
template<class T> cMedImage<T>::~cMedImage()
{
    // free memory
    destroy();
}


//=========================================================================
template<class T> T *cMedImage<T>::operator[](
        unsigned r) const	// [in] row index
{
    // get requested row
    return data[r];
}


//=========================================================================
template<class T> bool cMedImage<T>::imagePresent()
{
    // returns ture if valid image present
    return data != NULL && nRow > 0 && nCol > 0;
}


//=========================================================================
template<class T> void cMedImage<T>::create(
        unsigned rows,	// [in] number of rows in image
        unsigned cols)	// [in] number of columns in image
// Creates a new image of given size. The data is not initialized.
// With with rows or colss being zero it only deletes memory.
{
    unsigned i;	// local index varaible

    try {
        if (rows != nRow || cols != nCol) {
            // delete previous memory
            destroy();

            // get new memory
            if (rows > 0 && cols > 0) {
                data = new T*[rows];
                data[0] = new T[rows*cols];
            }

            // initialize pointer array
            for (i = 1; i < rows && cols > 0; i++)
                data[i] = data[i - 1] + cols;

            // keep size
            nRow = rows;
            nCol = cols;
        }
    }
    catch (...) {
        destroy();
        throw std::runtime_error("Could not create medical image.");
    }
}


//=========================================================================
template<class T> void cMedImage<T>::create(
        unsigned rows,	// [in] number of rows in image
        unsigned cols,	// [in] number of columns in image
        const T &value)	// [in] default pixel value
// Creates a new image of given size. The data is not initialized.
// With with rows or colss being zero it only deletes memory.
{
    create(rows, cols);
    set(value);
}


//=========================================================================
template<class T> void cMedImage<T>::destroy()
{
    // clean up
    if (data) {
        if (data[0]) delete[] data[0];
        delete[] data;
        data = NULL;
    }
    nRow = nCol = 0;
    pixelValueForBlack = 0;
    pixelValueForWhite = 1;
}


//=========================================================================
template<class T> void cMedImage<T>::set(
        const cMedImage<T> &image)	// [in] image to be copied from
{
    try {
        // get memory
        create(image.nRow, image.nCol);

        // copy content
        if (nCol > 0 && nRow > 0) {
            std::memcpy(data[0], image.data[0], nRow*nCol*sizeof(T));
            pixelValueForBlack = image.pixelValueForBlack;
            pixelValueForWhite = image.pixelValueForWhite;
        }

    }
    catch (...) {
        destroy();
        throw std::runtime_error("Could not copy medical image");
    }
}


//=========================================================================
template<class T> void cMedImage<T>::set(
        const T &value)	// [in] value to which all pixels will be set to
{
    unsigned r, c;	// index for row and column

    // loop over all pixels
    for (r = 0; r < nRow; r++)
        for (c = 0; c < nCol; c++)
            data[r][c] = value;
}


//=========================================================================
template<class T> void cMedImage<T>::add(
        const T &value)	// [in] value to be added to all pixels
{
    // check if data present
    if (nRow == 0 || nCol == 0 || data == NULL) return;

    try {
        // add value to all elements
        for (unsigned r = 0; r < nRow; r++)
            for (unsigned c = 0; c < nCol; c++)
                data[r][c] += value;
    }
    catch (...) {
        throw std::runtime_error("Could not add constant value to medical image.");
    }
}


//=========================================================================
template<class T> void cMedImage<T>::add(
        const cMedImage<T> &image)	// [in] image to be added
{
    unsigned rmax = nRow<image.nRow?nRow:image.nRow;	// min row number
    unsigned cmax = nCol<image.nCol?nCol:image.nCol;	// min column number

    try {
        // subtract pixel by pixel
        for(unsigned r=0; r<rmax; r++)
            for(unsigned c=0; c<cmax; c++)
                data[r][c] += image[r][c];
    }
    catch (...) {
        throw std::runtime_error("Could not add medical image.");
    }
}


//=========================================================================
template<class T> void cMedImage<T>::subtract(
        const T &value)	// [in] value to be subtracted from all pixels
{
    // check if data present
    if (nRow == 0 || nCol == 0 || data == NULL) return;

    try {
        // subtract value from all elements
        for (unsigned r = 0; r < nRow; r++)
            for (unsigned c = 0; c < nCol; c++)
                data[r][c] -= value;
    }
    catch (...) {
        throw std::runtime_error("Could not subtract constant value from medical image.");
    }
}


//=========================================================================
template<class T> void cMedImage<T>::subtract(
        const cMedImage<T> &image)	// [in] image to be added
{
    unsigned rmax = nRow<image.nRow?nRow:image.nRow;	// min row number
    unsigned cmax = nCol<image.nCol?nCol:image.nCol;	// min column number

    try {
        // subtract pixel by pixel
        for(unsigned r=0; r<rmax; r++)
            for(unsigned c=0; c<cmax; c++)
                data[r][c] -= image[r][c];
    }
    catch (...) {
        throw std::runtime_error("Could not subtract medical image.");
    }
}


//=========================================================================
template<class T> void cMedImage<T>::multiply(
        const T &value)	// [in] factor by which all pixels will be multiplied
{
    // check if data present
    if (nRow == 0 || nCol == 0 || data == NULL) return;

    try {
        // multiply all ellements by value
        for (unsigned r = 0; r < nRow; r++)
            for (unsigned c = 0; c < nCol; c++)
                data[r][c] *= value;
    }
    catch (...) {
        throw std::runtime_error("Could not multiply medical image by constant value.");
    }
}


//=========================================================================
template<class T> void cMedImage<T>::multiply(
        const cMedImage<T> &image)	// [in] image to be added
{
    unsigned rmax = nRow<image.nRow?nRow:image.nRow;	// min row number
    unsigned cmax = nCol<image.nCol?nCol:image.nCol;	// min column number

    try {
        // subtract pixel by pixel
        for(unsigned r=0; r<rmax; r++)
            for(unsigned c=0; c<cmax; c++)
                data[r][c] *= image[r][c];
    }
    catch (...) {
        throw std::runtime_error("Could not multiply medical image pixelwise.");
    }
}


//=========================================================================
template<class T> cMedImage<T> &cMedImage<T>::operator=(
        const cMedImage &image)	// [in] image to be copied from
{
    // copy image
    set(image);
    return *this;
}


//=========================================================================
template<class T> cMedImage<T> &cMedImage<T>::operator=(
        const T &value)	// [in] value to which all pixels will be set to
{
    // set value
    set(value);
    return *this;
}


//=========================================================================
template<class T> cMedImage<T> &cMedImage<T>::operator+=(
        const T &value)	// [in] value to be added to all pixels
{
    add(value);
    return *this;
}


//=========================================================================
template<class T> cMedImage<T> &cMedImage<T>::operator+=(
        const cMedImage<T> &image)	// [in] image to be added
{
    add(image);
    return *this;
}


//=========================================================================
template<class T> cMedImage<T> &cMedImage<T>::operator-=(
        const T &value)	// [in] value to be added to all pixels
{
    subtract(value);
    return *this;
}


//=========================================================================
template<class T> cMedImage<T> &cMedImage<T>::operator-=(
        const cMedImage<T> &image)	// [in] image to be added
{
    subtract(image);
    return *this;
}


//=========================================================================
template<class T> cMedImage<T> &cMedImage<T>::operator*=(
        const T &value)	// [in] factor by which all pixels will be multiplied
{
    multiply(value);
    return *this;
}


//=========================================================================
template<class T> cMedImage<T> &cMedImage<T>::operator*=(
        const cMedImage<T> &image)	// [in] image to be added
{
    multiply(image);
    return *this;
}


//=========================================================================
template<class T> void cMedImage<T>::findMinMax(
        T &min,	// [out] minimum pixel value
        T &max)	// [out] maximum pixel value
{
    unsigned r, c;	// index for row and column

    try {
        // check for valid image
        if (data == NULL || nRow == 0 || nCol == 0) throw - 1;

        // find min and max value
        min = max = data[0][0];
        for (r = 0; r < nRow; r++) {
            for (c = 0; c < nCol; c++) {
                if (min > data[r][c]) min = data[r][c];
                if (max < data[r][c]) max = data[r][c];
            }
        }
    }
    catch (...) {
        min = max = 0;
        throw std::runtime_error("Could not find minimum and maximum pixel value in medical image.");
    }
}


//=========================================================================
template<class T> void cMedImage<T>::cut(
        unsigned r1,	// [in] first row from above which will not be cut away
        unsigned c1,	// [in] first column from left which will not be cut away
        unsigned r2,	// [in] first row on the bottom that will be cut away
        unsigned c2)	// [in] first column on the right that will be cut away
{
    cMedImage<T> reduced;	// new image
    unsigned r, c;			// index for row and column

    try {
        // swap limits if required
        if (r1 > r2) {
            r = r1;
            r1 = r2;
            r2 = r;
        }
        if (c1 > c2) {
            c = c1;
            c1 = c2;
            c2 = c;
        }

        // adjust limits
        if (r1 > nRow) r1 = nRow;
        if (r2 > nRow) r2 = nRow;
        if (c1 > nCol) c1 = nCol;
        if (c2 > nCol) c2 = nCol;

        // if zero size
        if (r1 == r2 || c1 == c2) destroy();

            // if size larger than zero
        else {
            // set size
            reduced.create(r2 - r1, c2 - c1);

            // copy windowing values
            reduced.pixelValueForBlack = pixelValueForBlack;
            reduced.pixelValueForWhite = pixelValueForWhite;

            // copy pixels
            for (r = 0; r < reduced.nRow; r++)
                for (c = 0; c < reduced.nCol; c++)
                    reduced.data[r][c] = data[r + r1][c + c1];

            // copy image over old one
            set(reduced);
        }
    }
    catch (...) {
        throw std::runtime_error("Could not cut medical image.");
    }
}


//=========================================================================
template<class T> void cMedImage<T>::shift(
        int r0,			// [in] vertical shift (positive=down)
        int c0,			// [in] horizontal shift (positive=right)
        const T &fill)	// [in] default values for pixels shifted into the image
{
    int r, c;	// index for row and column

    try {
        // if shift greater than image size
        if (r0 >= (int)nRow || r0 <= -(int)nRow || c0 >= (int)nCol || c0 <= -(int)nCol) set(fill);

            // shift negative row and negative column
        else if (r0 < 0 && c0 < 0) {
            for(r=0; r<(int)nRow; r++) {
                for(c=0; c<(int)nCol; c++) {
                    if(r-r0<(int)nRow && c-c0<(int)nCol) data[r][c] = data[r-r0][c-c0];
                    else data[r][c] = fill;
                }
            }
        }

            // shift negative row and positive column
        else if (r0 < 0 && c0 >= 0) {
            for(r=0; r<(int)nRow; r++) {
                for(c=nCol-1; c>=0; c--) {
                    if(r-r0<(int)nRow && c>=c0) data[r][c] = data[r-r0][c-c0];
                    else data[r][c] = fill;
                }
            }
        }

            // shift positive row and negative column
        else if (r0 >= 0 && c0 < 0) {
            for(r=nRow-1; r>=0; r--) {
                for(c=0; c<(int)nCol; c++) {
                    if(r>=r0 && c-c0<(int)nCol) data[r][c] = data[r-r0][c-c0];
                    else data[r][c] = fill;
                }
            }
        }

            // shift positive row and positive column
        else if (r0 >= 0 && c0 >= 0) {
            for(r=nRow-1; r>=0; r--) {
                for(c=nCol-1; c>=0; c--) {
                    if(r>=r0 && c>=c0) data[r][c] = data[r-r0][c-c0];
                    else data[r][c] = fill;
                }
            }
        }
    }
    catch (...) {
        throw std::runtime_error("Could not shift medical image.");
    }
}


//=========================================================================
template<class T> void cMedImage<T>::mirror(
        bool vertical)	// [in] reflection direction (true=exchange of upper and lower part)
{
    unsigned r, c;	// index for row and column
    T tmp;

    if (vertical) {
        for (c = 0; c<nCol; c++) {
            for (r = 0; r<nRow/2; r++) {
                tmp = data[r][c];
                data[r][c] = data[nRow-r-1][c];
                data[nRow-r-1][c] = tmp;
            }
        }
    } else {
        for (r = 0; r<nRow; r++) {
            for (c = 0; c<nCol/2; c++) {
                tmp = data[r][c];
                data[r][c] = data[r][nCol - c - 1];
                data[r][nCol - c - 1] = tmp;
            }
        }
    }
}


//=========================================================================
template<class T> void cMedImage<T>::mirror(
        const cMedImage<T> &image,	// [in] image to be copied from
        bool vertical)	// [in] reflection direction (true=exchange of upper and lower part)
{
    unsigned r, c;	// index for row and column

    create(image.nRow, image.nCol);

    if (vertical)
        for (r = 0; r < nRow; r++)
            for (c = 0; c < nCol; c++)
                data[nRow - r - 1][c] = image[r][c];
    else
        for (r = 0; r < nRow; r++)
            for (c = 0; c < nCol; c++)
                data[r][nCol - c - 1] = image[r][c];
}


//=========================================================================
template<class T> void cMedImage<T>::rotate(unsigned dir)
{
    cMedImage<T> image=*this;

    rotate(image, dir);
}


//=========================================================================
template<class T> void cMedImage<T>::rotate(
        const cMedImage<T> &image,
        unsigned dir)	// counter clockwise: 1=90deg, 2=180deg, 3=270deg
{
    if(dir%2==0) create(image.nRow, image.nCol);
    else create(image.nCol, image.nRow);

    switch(dir%4) {
        case 0:
            for(unsigned r=0; r<nRow; r++)
                for(unsigned c=0; c<nCol; c++)
                    data[r][c] = image[r][c];
            break;
        case 1:
            for(unsigned r=0; r<nRow; r++)
                for(unsigned c=0; c<nCol; c++)
                    data[r][c] = image[c][nRow-r-1];
            break;
        case 2:
            for(unsigned r=0; r<nRow; r++)
                for(unsigned c=0; c<nCol; c++)
                    data[r][c] = image[nRow-r-1][nCol-c-1];
            break;
        case 3:
            for(unsigned r=0; r<nRow; r++)
                for(unsigned c=0; c<nCol; c++)
                    data[r][c] = image[nCol-c-1][r];
            break;
    }
}


//=========================================================================
template<class T> void cMedImage<T>::roll(int vert, int hori)
{
    cMedImage<T> image=*this;

    roll(image, vert, hori);
}


//=========================================================================
template<class T> void cMedImage<T>::roll(const cMedImage<T> &image, int vert, int hori)
{
    // set size of new image
    create(image.nRow, image.nCol);

    // adjust horizontal and vertical shift
    if(hori<0) hori = nCol-(-hori)%nCol;
    hori %= nCol;
    if(vert<0) vert = nRow-(-vert)%nRow;
    vert %= nRow;

    // copy pixels
    for(unsigned r=0; r<nRow; r++)
        for(unsigned c=0; c<nCol; c++)
            data[(r+vert)%nRow][(c+hori)%nCol] = image[r][c];
}


//=========================================================================
template<class T> void cMedImage<T>::autoWindowing(
        bool rad/*=true*/)	// [in] true if in RAD-mode, i.e. inverted grey scale
{
    if(rad) findMinMax(pixelValueForWhite, pixelValueForBlack);
    else findMinMax(pixelValueForBlack, pixelValueForWhite);
}

#endif // #ifndef HEADER_MEDICAL_IMAGE_INCLUDED
