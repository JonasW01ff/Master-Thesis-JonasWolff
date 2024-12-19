
#ifndef CODE_LIBRARY_H
#define CODE_LIBRARY_H

#include <vector>
#include <iostream>
#include <functional>
#include <random>
#include "numbers.h"
#include "adjointNumeric.h"
#include "settings.h"


using namespace std;

template <class T1>
class mat {
private:
    // Matrix object, it take the form row1, row2, ...
    vector<T1> data;

    // Number of rows
    int N;

    // Number of columns
    int M;

    // Transpose
    bool isTransposed = false;

public:

    // Constructor
    mat<T1>();

    // Overloaded Constructor
    mat<T1>(int &N, int &M);

    // Overloaded Constructor
    mat<T1>(const int &N, const int &M);

    // Overloaded Constructor
    mat<T1>(const int &N, const int &M, const T1& val);

    mat<T1>(const mat<T1> &rhs)=default;

    mat<T1>(mat<T1> &&rhs)=default;

    mat<T1>& operator=(const mat<T1> &) = default;
    mat<T1>& operator=(mat<T1>&&) = default;

    // mat<num> -> mat<adjointIntegral>
    template<typename  = enable_if<is_same_v<num, T1>>>
    explicit mat<T1>(const mat<num> &rhs){
        data.resize(rhs.nRows()* rhs.nCols());
        N = rhs.nRows();
        M = rhs.nCols();
        const auto rhsdata = rhs.begin();
        for (int i=0;i<data.size(); i++)
            data[i] = rhsdata[i]; // Implicit

    }

    template<class T2,
            typename = enable_if<!is_same_v<T2,adjointIntegral> &&
                                 is_same_v<T1,adjointIntegral> &&
                                 !is_same_v<T2,dual<num>>>>
    explicit mat<T1>(const mat<T2> &rhs){
        data.resize(rhs.nRows()* rhs.nCols());
        N = rhs.nRows();
        M = rhs.nCols();
        const auto rhsdata = rhs.begin();
        for (int i=0;i<data.size(); i++)
            data[i] = rhsdata[i];
    };

    template<class T2,
            typename = enable_if<!is_same_v<T2,adjointIntegral> &&
                                is_same_v<T1,adjointIntegral> &&
                                !is_same_v<T2,dual<num>>>>
    mat<T1>&operator=(const mat<T2> &rhs){
        static_assert(!is_same_v<T2,dual<num>>);
        data.resize(rhs.nRows()* rhs.nCols());
        N = rhs.nRows();
        M = rhs.nCols();
        const auto rhsdata = rhs.begin();
        for (int i=0;i<data.size(); i++)
            data[i] = adjointIntegral(rhsdata[i]);
        return *this;
    };


    // Destructor
    ~mat<T1>();

    // Transpose
    void transpose() {
        isTransposed = !isTransposed;
    }

    // In math notation, this is just temporary
    mat<T1> T() const {mat<T1> myMat(*this); myMat.transpose(); return myMat;};

    const bool getIsTransposed() const {return isTransposed;};

    void print() const;

    // Get Entry
    const T1 getEntry(int i, int j) const;

    // Get Iterator starting at entry
    vector<T1>::iterator getIterator(const int&i, const int&j);

    // Get Const Iterator starting at entry
    const vector<T1>::const_iterator getConstIterator(const int&i, const int&j) const;


    // Get entry reference
    const T1& getReference(const int &i,const int &j) const;

    // Get entry mutable pointer
    T1 * getPointer(const int&i, const int &j);

    // Get entry reference (be very careful)
    T1& getMutableReference(const int &i,const int &j);

    // Resize matrix and preserve surviving indexes
    void resize(const int& N, const int& M);

    // Get submatrix
    mat getSubMatrix(const int& rowStart, const int& rowEnd,
                     const int& columnStart, const int& columnEnd) const;

    // Remove row
    void removeRow(const int& i);

    // Remove column
    void removeColumn(const int& j);

    // Set Entry
    void setEntry(int i, int j, T1 value);

    // Perform binary operator entrywise
    void setEntries(function<const T1(const T1, const T1)> op, const T1 &rhs);

    // Perform binary operator entrywise
    void setEntries(function<const T1(const T1, const T1)> op, const mat &rhs);

    // Perform binary operator entrywise
    void setEntries(function<const T1(const T1)> op);

    // Get number of rows
    int nRows() const;

    // Get number of columns
    int nCols() const;

    // get size
    size_t size() const {return data.size();}

    // Get Row
    void getRow(int i, vector<T1> &res);

    // Get Column
    void getCol(int j, vector<T1> &res);

    // Get start of matrix iterator
    const vector<T1>::const_iterator begin() const {return data.begin();};

    // Get end of matrix iterator
    const vector<T1>::const_iterator end() const {return data.end();};


    // Get mutable start of matrix iterator
    vector<T1>::iterator beginMutable() {return data.begin();};

    // Get mutable end of matrix iterator
    vector<T1>::iterator endMutable() {return data.end();};

    // Matrix multiplication
    void matmul(const mat<T1> &xmat, mat<T1> &res) const;

    // Get the orthonormal basis for this' image
    mat<T1> getOrthonormalBasis() const;

    // Get the orthogonal projection matrix
    mat<T1> getOrthogonalProjectionMatrix() const;

    // Project on to image of basis matrix' columns
    mat<T1> project(const mat<T1> &myVector) const;

    mat<T1> matmul(const mat<T1>& rhs) const {
        mat<T1> res(this->nRows(), rhs.nCols());
        matmul(rhs, res);
        return res;};

    template<typename = enable_if<!is_same_v<T1, dual<num>>>>
    mat<adjointIntegral> matmul(const mat<adjointIntegral> &rhs) const{
        mat<adjointIntegral> myMat(nRows(),nCols());
        myMat = *this;
        mat<adjointIntegral> res(nRows(),rhs.nCols());
        for (int i=0; i< res.nRows(); i++){
            for (int j=0;j<res.nCols();j++)
                res[i,j] = 0.;
        }
        myMat.matmul(rhs,res);
        return res;

    }

    T1 norm() const;

    int idxFindNonZeroInColumnAfterRow(const int &colidx, const int &rowidx) const;

    void swapRows(const int &rowidx1, const int &rowidx2);

    void addRow(const int &readRow, const int &writeRow, const T1& multiplier);

    void scaleRow(const int &rowidx, const T1& multiplier);

    mat<T1> makeEchelon();

    mat<T1> makeReducedEchelon();

    vector<int> makeColumnIdependent();

    void inv();

    T1 det() const;

    // Matrix addition
    mat operator+(const mat& rhs) const;

    // Matrix addition
    mat operator+=(const mat& rhs);

    // Matrix unary minus
    mat operator-() const &;

    mat operator-() &&;

    // Matrix subtraction
    mat<T1> operator-(const mat<T1>& rhs) const;

    template<typename = enable_if<!is_same_v<T1, dual<num>>>>
    mat<adjointIntegral> operator-(mat<adjointIntegral>& rhs) const{

        // Check thats size equals
        if (this->nCols() != rhs.nCols() or this->nRows() != rhs.nRows())
            throw runtime_error("Matrix addition needs equal sizes");

        if (isTransposed or rhs.getIsTransposed())
            throw runtime_error("This operation is not implemented yet");

        mat<adjointIntegral> temp(nRows(), nCols());
        for (int i=0;i<size();i++)
            *(temp.beginMutable() + i) =  *(begin()+i) - *(rhs.begin()+i);

        return temp;

    };

    template<typename = enable_if<!is_same_v<T1, dual<num>>>>
    mat<adjointIntegral> operator-(mat<adjointIntegral>&& rhs) const{

        // Check thats size equals
        if (this->nCols() != rhs.nCols() or this->nRows() != rhs.nRows())
            throw runtime_error("Matrix addition needs equal sizes");

        if (isTransposed or rhs.getIsTransposed())
            throw runtime_error("This operation is not implemented yet");

        mat<adjointIntegral> temp(nRows(), nCols());
        for (int i=0;i<size();i++)
            *(temp.beginMutable() + i) =  *(begin()+i) - *(rhs.begin()+i);

        return temp;

    };

    // Matrix suctraction
    mat operator-=(const mat& rhs);

    // Matrix multiplaction
    mat operator*(const mat& rhs) const;

    // Matrix multiplaction
    mat operator*=(const mat& rhs);

    // Scalar multiplication
    mat operator*(const T1& rhs) const;

    // Scalar division
    mat operator/(const T1& rhs) const;

    // Scalar multiplication
    mat operator+(const T1& rhs) const;

    // Get mutable refrence
    T1 operator[](const size_t& row_idx, const size_t& col_idx) && {
        return getEntry(row_idx,col_idx);
    };

    // Get mutable refrence
    T1& operator[](const size_t& row_idx, const size_t& col_idx) & {
        return getMutableReference(row_idx,col_idx);
    };

    //Get constant refrence
    const T1& operator[](const size_t& row_idx, const size_t& col_idx) const & {
        return getReference(row_idx,col_idx);
    };


    // Get mutable row refrance
    T1* operator[](const size_t& row_idx){
        if (isTransposed) throw runtime_error("You are reading memory wrong");
        return &getMutableReference(row_idx,0);
    }

    //Get constant refrence
    const T1 operator[](const size_t& row_idx, const size_t& col_idx) const && {
        return getEntry(row_idx,col_idx);
    };

    //Get constant row refrence
    const T1* operator[](const size_t& row_idx) const {
        if (isTransposed) throw runtime_error("You are reading memory wrong");
        return &getReference(row_idx,0);
    };

    void setAdjointVal(const mat<adjointIntegral>&);
    void setAdjointAdj(const mat<adjointIntegral>&);


};


// Vector on Vector projection
void project(const auto vecStartIt, const auto vecEndIt,
             const auto baseStartIt, auto resStartIt);

// Tests
void matmultest();

void projectionMatrixTest();

#endif //CODE_LIBRARY_H
