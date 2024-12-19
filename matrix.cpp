
/**
 * @file matrix.cpp
 * @brief Implementation file for the matrix class.
 */

#if __APPLE__
#include <Accelerate/Accelerate.h>
//#include "matrixMetal.h"
#endif

#include "matrix.h"
#include <numeric>
#include <execution>
#include <algorithm>
#include <memory_resource>
#include "numbers.h"
#include <array>
#include "ThreadPool.h"
#include "settings.h"

ThreadPool *myPool = &ThreadPool::getThreadPool();

/**
 * @brief Default constructor for the matrix class.
 */
template<class T1>
mat<T1>::mat() {

};

/**
 * @brief Overloaded constructor for the matrix class.
 * @param N The number of rows in the matrix.
 * @param M The number of columns in the matrix.
 */
template<class T1>
mat<T1>::mat(int &myN, int &myM) {

	// Initialize matrix with size NxM
	data.resize(myN*myM, T1(0.));

	this->N = myN;
	this->M = myM;
};

template<class T1>
mat<T1>::mat(const int &myN,const int &myM) {

    N = myN;
    M = myM;

    size_t mySize = myN*myM;
    // Initialize matrix with size NxM
    data.resize(mySize);
};

template<class T1>
mat<T1>::mat(const int &myN,const int &myM, const T1& val) {

    // Initialize matrix with size NxM
    data = vector<T1>(myN*myM, val);

    this->N = myN;
    this->M = myM;
};


/*template<class T1>
template<class T2>
mat<T1>::mat(const mat<T2> &rhs){
    data.resize(rhs.begin() - rhs.end());
    transform(rhs.begin(), rhs.end(), data.begin(),
              [](const T2 &a){return T1(a);});
    N = (rhs.getIsTransposed()) ? nCols() : nRows();
    M = (rhs.getIsTransposed()) ? nRows() : nCols();
    isTransposed = rhs.getIsTransposed();
};*/


/**
 * @brief Destructor for the matrix class.
 */
template<class T1>
mat<T1>::~mat<T1>() = default;

/**
 * @brief Get the value of a specific entry in the matrix.
 * @param i The row index of the entry.
 * @param j The column index of the entry.
 * @return The value of the entry at position (i, j).
 * @throws std::runtime_error if the entry is outside the matrix bounds.
 */
template<class T1>
const T1 mat<T1>::getEntry(int i, int j) const {

    // Get vector idx
    size_t idx;
    if (isTransposed)
        idx = j*M+i;
    else
        idx = i*M+j;

	// Check that entry is inside matrix bounds.
	if (idx > data.size()){
		throw std::runtime_error("The matrix has too few entries.");
	}

	// Return entry.
	return data[idx];
}

template<class T1>
const T1& mat<T1>::getReference(const int &i,const int &j) const{

    // Get vector idx
    size_t idx;
    if (isTransposed)
        idx = j*M+i;
    else
        idx = i*M+j;

	// Check that entry is inside matrix bounds.
	if (idx >= data.size()){
		throw std::runtime_error("The matrix has too few entries.");
	}

	// Return entry.
    return data[idx];
}

template<class T1>
T1& mat<T1>::getMutableReference(const int &i,const int &j){

    // Get vector idx
    size_t idx;
    if (isTransposed)
        idx = j*M+i;
    else
        idx = i*M+j;

    // Check that entry is inside matrix bounds.
    if (idx >= data.size()){
        throw std::runtime_error("The matrix has too few entries.");
    }

    // Return entry.
    return data[idx];
}

template<class T1>
void mat<T1>::setEntry(int i, int j, T1 value) {

    // Get vector idx
    size_t idx;
    if (isTransposed)
        idx = j*M+i;
    else
        idx = i*M+j;

	// Check that entry is inside matrix bounds.
	if (idx >= data.size()){
		throw std::runtime_error("The matrix has too few entries.");
	}

	// Set entry value
	data[idx] = value;

}

template<class T1>
int mat<T1>::nRows() const{
	return (isTransposed) ? M : N;
}

template<class T1>
int mat<T1>::nCols() const{
	return (isTransposed) ? N : M;
}

template<class T1>
void mat<T1>::getRow(int i,vector<T1> &res) {

	// Check that the dimensions match
	if (res.size() != ((isTransposed) ? N : M)){
		throw runtime_error("Vector needs to be the same size as the number of columns");
	}

	// Allocate row
    for (int k = 0; k < res.size(); k++) {
        res[k] = getEntry(i,k);
    }

}

template<class T1>
void mat<T1>::getCol(int j, vector<T1> &res) {

	//Check that dimension match
	if (res.size() != ((isTransposed) ? M : N)){
		throw runtime_error("Vector needs to be the same size as the number of rows");
	}

	// Allocate Column
    for (int k=0; k<res.size(); k++){
        res[k] = getEntry(k,j);
    }

}

#if __APPLE__
/*template<>
void mat<num>::matmul(const mat<num> &xmat, mat<num> &res) const {

    cblas_dgemm(
            CblasRowMajor,
            (getIsTransposed()) ? CblasTrans : CblasNoTrans,
            (xmat.getIsTransposed()) ? CblasTrans : CblasNoTrans,
            nRows(),
            xmat.nCols(),
            nCols(),
            1.,
            &data[0],
            M,
            &xmat.data[0],
            xmat.M,
            0.,
            &res.data[0],
            res.M
    );
    return;
}*/
#endif

// Loop unroling
constexpr size_t unroolsize = 12;

// constant times array stored on array
template<class T1, size_t N>
class vcvClass {
public:
    constexpr static void vcv(const T1&c, const T1* vIn, T1* vOut){
        vOut[N] += c * vIn[N];
        vcvClass<T1, N - 1>::vcv(c, vIn, vOut);
    }
};

template<class T1>
class vcvClass<T1, 0>{
public:
    constexpr static void vcv(const T1& c,const T1* vIn, T1* vOut){
        vOut[0] += c * vIn[0];
    }
};

template<class T1, size_t... Is>
constexpr array<void(*)(const T1&, const T1*, T1*), unroolsize> vcvCompile(index_sequence<Is...>){
    array<void(*)(const T1&, const T1*, T1*), sizeof...(Is)> myVec {&vcvClass<T1, Is>::vcv...};
    /*for (size_t i= 0; i<unroolsize; i++) {
        *myVec[i] = &vcvClass<T1, i>::vcv;
    }*/

    return myVec;
};

template<class T1>
constexpr auto vcvCompiled = vcvCompile<T1>(make_index_sequence<unroolsize>{});

template<class T1>
void mat<T1>::matmul(const mat<T1> &xmat, mat<T1> &res) const {

	// Check that dimensions match
	if(this->nCols() != xmat.nRows()){
		throw std::runtime_error("Matrix multiplication can only take the format MxN * N*K");
	}

    // check res
    if (this->nRows()!=res.nRows() or res.nCols() != xmat.nCols())
        throw runtime_error("Result matrix has wrong dimensions");

    if (res.isTransposed)
        throw runtime_error("Res transposed");

    if (not (isTransposed or xmat.isTransposed)){

        // perform matrix multiplication
        const T1 *this_i, *xmat_k;
        T1 this_ik;

        for (int i=0; i<nRows(); i++){
            T1 *res_i = res[i];
            this_i = (*this)[i];
            for (int k=0; k<nCols(); k++){
                this_ik = this_i[k];
                xmat_k = xmat[k];
                for (int j = 0; j < xmat.nCols(); j++)
                    res_i[j] += this_ik * xmat_k[j];
            }
        }

        /*cout << "I got called" << endl;*/
        /*for (int j=0; j<xmat.nCols(); j++){
            for (int k=0; k<nCols(); k++){
                for (int i=0; i<nRows(); i++)
                    res[i,j] += (*this)[i,k]*xmat[k,j];
            }
        }*/

        return;

        // perform matrix multiplication

        int lhsNCols = nCols();
        int rhsNCols = xmat.nCols();

        // Get current lhs entry
        auto lhsEntry = (*this).begin();

        // Get row start for result matrix
        auto resRowStart = res.getIterator(0,0);

        // Get last entry of row lhs
        auto lhsRowEnd = lhsEntry+M;

        for (; lhsRowEnd<data.end()+1; lhsRowEnd += M){

            // Get xmat start and end iterators
            auto rhsRowStart = xmat.begin();
            auto rhsRowEnd = rhsRowStart + lhsNCols;

            for (; lhsEntry<lhsRowEnd; lhsEntry++){

                // Perform vector sum product
                transform(rhsRowStart,
                          rhsRowEnd,
                          resRowStart,
                          resRowStart,
                          [&lhsEntry](const T1 &a, const T1 &b){return b+(*lhsEntry)*a;}
                );

                // Move forward
                rhsRowStart = rhsRowEnd;
                rhsRowEnd += rhsNCols;

            }

            // Get result start iterator
            resRowStart += rhsNCols;
        }
    } else if (isTransposed and not xmat.isTransposed){

        if (nCols()*xmat.nCols() < 100000000 or nRows()==1){
            const T1 *xmat_k;
            T1 * res_i;
            T1 this_ik;
            for (int i=0; i<nRows(); i++){
                res_i = res[i];
                for (int k=0; k<nCols(); k++){
                    xmat_k = xmat[k];
                    this_ik = (*this)[i,k]; // traverse column for extra speed if needed
                    for (int j=0; j<xmat.nCols(); j++) {
                        res_i[j] += this_ik * xmat_k[j];
                        //res[i, j] += (*this)[i, k] * xmat[k, j];
                    }
                }
            }

        } else {

            vector<future<void>> myfuts(nRows());

            // perform matrix multiplication
            for (int i = 0; i < nRows(); i++) {
                auto calc_i = [this, &xmat, &res, i]() {
                    for (int k = 0; k < nCols(); k++) {
                        for (int j = 0; j < xmat.nCols(); j++)
                            res[i, j] += (*this)[i, k] * xmat[k, j];
                    }
                };

                myfuts[i] = std::move(myPool->spawnTask(calc_i));
                //for (int k=0; k<nCols(); k++){
                //    for (int j=0; j<xmat.nCols(); j++)
                //        res[i,j] += (*this)[i,k]*xmat[k,j];
                //}
            }
            for (auto &fut: myfuts)
                myPool->helpQueueWhile(fut);
        }

    } else if (xmat.isTransposed and not isTransposed) {

        // perform matrix multiplication
        T1 *res_i;
        const T1 *this_i, *xmat_kj;
        for (int i=0; i<nRows(); i++){
            this_i = (*this)[i];
            res_i = res[i];
            for (int j=0; j<xmat.nCols(); j++){
                T1& res_ij = res_i[j];
                xmat_kj = &xmat[0,j];
                for (int k=0; k<nCols(); k++) {
                    //xmat_kj++;
                    res_ij += this_i[k] * xmat_kj[k];//(*xmat_kj);
                    //res[i, j] += (*this)[i, k] * xmat[k, j];
                }
            }
        }



    } else if (isTransposed && xmat.isTransposed){

        // perform matrix multiplication
        /*mat<T1> thiscopy = this->T();
        xmat.matmul(move(thiscopy), res);

        res.transpose();
        return;*/

        for (int k=0; k<nCols(); k++){
            for (int j=0; j<xmat.nCols(); j++){
                for (int i=0; i<nRows(); i++)
                    res[i,j] += (*this)[i,k]*xmat[k,j];

            }
        }
    } else  {
        throw runtime_error("Matrix multiplication has made an error");
    }
}

template<class T1>
T1 mat<T1>::norm() const {
    T1 myNorm = reduce(data.begin(), data.end(), T1(0.),
                       [](const T1 &mySum,const  T1 &rhs){return mySum + rhs*rhs;});
    return sqrt(myNorm);
}

#if __APPLE__XXX
template<>
mat<num> mat<num>::operator+(const mat &rhs) const {

    // Check thats size equals
    if (this->nCols() != rhs.nCols() or this->nRows() != rhs.nRows())
        throw runtime_error("Matrix addition needs equal sizes");

    if (isTransposed or rhs.isTransposed)
        throw runtime_error("This operation is not implemented yet");

    mat temp(nRows(), nCols());

    myMatrixMetal.madd(&data[0], &rhs.data[0], &temp.data[0], data.size());
    return temp;

}
#endif

template<class T1>
mat<T1> mat<T1>::operator+(const mat &rhs) const {

    // Check thats size equals
    if (this->nCols() != rhs.nCols() or this->nRows() != rhs.nRows())
        throw runtime_error("Matrix addition needs equal sizes");

    if (isTransposed or rhs.isTransposed)
        throw runtime_error("This operation is not implemented yet");

    mat temp(nRows(), nCols());

    transform(this->begin(), this->end(),
              rhs.begin(), temp.beginMutable(),
              plus<>());

    return temp;

}

template<class T1>
mat<T1> mat<T1>::operator+=(const mat &rhs) {

    // Check thats size equals
    if (this->nCols() != rhs.nCols() or this->nRows() != rhs.nRows())
        throw runtime_error("Matrix addition needs equal sizes");

    if (isTransposed or rhs.isTransposed)
        throw runtime_error("This operation is not implemented yet");

    transform(this->begin(), this->end(),
              rhs.begin(), this->beginMutable(),
              plus<>());

    return *this;

}

template<class T1>
mat<T1> mat<T1>::operator-() const & {
    mat temp(nRows(), nCols());
    transform(begin(), end(), temp.beginMutable(),
              [](const T1 &x){return -x;});

    return temp;
}

template<class T1>
mat<T1> mat<T1>::operator-() && {
    transform(begin(), end(), beginMutable(),
              [](const T1 &x){return -x;});

    return *this;
}

template<class T1>
mat<T1> mat<T1>::operator-(const mat<T1> &rhs) const {

    // Check thats size equals
    if (this->nCols() != rhs.nCols() or this->nRows() != rhs.nRows())
        throw runtime_error("Matrix addition needs equal sizes");

    if (isTransposed or rhs.isTransposed)
        throw runtime_error("This operation is not implemented yet");

    mat temp(nRows(), nCols());
    transform(this->begin(), this->end(),
              rhs.begin(), temp.beginMutable(),
              minus<>());

    return temp;

}

template<class T1>
mat<T1> mat<T1>::operator-=(const mat &rhs) {

    // Check thats size equals
    if (this->nCols() != rhs.nCols() or this->nRows() != rhs.nRows())
        throw runtime_error("Matrix addition needs equal sizes");

    if (isTransposed or rhs.isTransposed)
        throw runtime_error("This operation is not implemented yet");

    transform(this->begin(), this->end(),
              rhs.begin(), this->beginMutable(),
              minus<>());

    return *this;

}

template<class T1>
mat<T1> mat<T1>::operator*(const mat &rhs) const {

    // Check thats size equals
    if (this->nCols() != rhs.nCols() or this->nRows() != rhs.nRows())
        throw runtime_error("Matrix addition needs equal sizes");

    if (isTransposed or rhs.isTransposed)
        throw runtime_error("This operation is not implemented yet");

    mat temp(nRows(), nCols());
    transform(this->begin(), this->end(),
              rhs.begin(), temp.beginMutable(),
              multiplies<>());

    return temp;

}

template<class T1>
mat<T1> mat<T1>::operator*=(const mat &rhs) {

    // Check thats size equals
    if (this->nCols() != rhs.nCols() or this->nRows() != rhs.nRows())
        throw runtime_error("Matrix addition needs equal sizes");

    if (isTransposed or rhs.isTransposed)
        throw runtime_error("This operation is not implemented yet");

    transform(this->begin(), this->end(),
              rhs.begin(), this->beginMutable(),
              multiplies<>());

    return *this;

}

template<class T1>
mat<T1> mat<T1>::operator*(const T1 &rhs) const {

    mat myMat = (*this);
    transform( this->begin(), this->end(),
               myMat.beginMutable(),
              [&rhs](const T1 &a)  {return a*rhs;});

    return myMat;

}

template<class T1>
mat<T1> mat<T1>::operator/(const T1& rhs) const{

    mat myMat = (*this);
    transform( this->begin(), this->end(),
               myMat.beginMutable(),
               [&rhs](const T1 &a)  {return a/rhs;});

    return myMat;


};

template<class T1>
mat<T1> mat<T1>::operator+(const T1 &rhs) const {

    mat myMat = (*this);
    transform( this->begin(), this->end(),
               myMat.beginMutable(),
               [&rhs](const T1 &a)  {return a+rhs;});


    return myMat;

}

template<class T1>
void mat<T1>::setEntries(function<const T1(const T1, const T1)> op, const T1 &rhs) {

    transform( this->begin(), this->end(),
               this->beginMutable(),
               [&rhs, &op](const T1 &a){return op(a,rhs);});
}

template<class T1>
void mat<T1>::setEntries(function<const T1(const T1, const T1)> op, const mat &rhs) {

    // Check thats size equals
    if (this->nCols() != rhs.nCols() or this->nRows() != rhs.nRows())
        throw runtime_error("Matrix addition needs equal sizes");

    transform( this->begin(), this->end(), rhs.begin(),
               this->beginMutable(), op);
}

template<class T1>
void mat<T1>::setEntries(function<const T1(const T1)> op) {

    transform( this->begin(), this->end(), this->beginMutable(), op);
}

template<class T1>
void mat<T1>::swapRows(const int &rowidx1, const int &rowidx2) {

    if (not isTransposed)
        swap_ranges(this->getIterator(rowidx1, 0),
                    this->getIterator(rowidx1, nCols()),
                    this->getIterator(rowidx2, 0));
    else {
        for (int j=0; j<this->nCols();j++)
            swap((*this)[rowidx1, j], (*this)[rowidx2, j]);
    }
}

template<class T1>
int mat<T1>::idxFindNonZeroInColumnAfterRow(const int &colidx, const int &rowidx) const {

    for (int i=rowidx+1; i< nRows(); i++){
        if (abs((*this)[i,colidx]) > 0.00001) return i;}

    return -1;

}

// Gaussian elimination to echelon form
template<class T1>
mat<T1> mat<T1>::makeEchelon() {

    // res[0][0] = RowRank(this)
    // res[1][0] = det(this)
    mat res(2,1);

    int rowsSolved = 0;
    T1 &det = res.getMutableReference(1,0);

    // Set start value for multiplication
    det = 1.0;

    // Solve columns one at a time
    for (int j = 0; j<nCols() && j<nRows(); j++){

        // Entry to make leading entry
        T1 &futureLeadingEntry = this->getMutableReference(rowsSolved, j);

        // If top entry is 0; swap with first non 0
        if (futureLeadingEntry == 0.) {
            int idx0 = this->idxFindNonZeroInColumnAfterRow(j, rowsSolved);

            // If all entries in column is 0 then continue to next column
            if (idx0 == -1) continue;
            this->swapRows( rowsSolved, idx0);
            futureLeadingEntry = this->getMutableReference(rowsSolved, j);
            det *= -futureLeadingEntry;


        } else {
            det *= futureLeadingEntry;
        }

        // Now first entry in column has a non-zero value
        // Subtract all other rows
        for (int i = rowsSolved+1; i<this->nRows(); i++){
            const T1 &entryTo0 = this->getReference(i, j);
            if (entryTo0 == 0) continue;
            this->addRow(rowsSolved, i, -entryTo0/futureLeadingEntry);
        }

        // Success fully created a new leading row
        rowsSolved++;


    }

    // Return rank and determinant
    res[0][0] = rowsSolved;
    det = (rowsSolved == nRows() && rowsSolved == nCols()) ? det : T1(0.);
    return res;

}

template<class T1>
mat<T1> mat<T1>::makeReducedEchelon() {


    if(isTransposed)
        throw runtime_error("makeReducedEchelon only works on non-transposed matrices");

    // res[0][0] = RowRank(this)
    // res[1][0] = det(this)
    mat res(2,1);

    int rowsSolved = 0;
    T1 &det = res.getMutableReference(1,0);

    // Set start value for multiplication
    det = 1.0;

    // Solve columns one at a time
    for (int j = 0; j<nCols() and rowsSolved<nRows(); j++){

        // Entry to make leading entry
        T1 &futureLeadingEntry = this->getMutableReference(rowsSolved, j);

        // If top entry is 0; swap with first non 0
        if (abs(futureLeadingEntry) < 0.00001) {
            int idx0 = this->idxFindNonZeroInColumnAfterRow(j, rowsSolved);

            // If all entries in column is 0 then continue to next column
            if (idx0 == -1) continue;
            if (rowsSolved != idx0)
                swapRows( rowsSolved, idx0);
            futureLeadingEntry = this->getMutableReference(rowsSolved, j);
            if (abs(futureLeadingEntry) < 0.00001)
                cout << "WTF" << endl;
            det *= -futureLeadingEntry;

            // Be careful, see that multiplier is indeed not a reference to futureLeadingEntry
            scaleRow(rowsSolved, 1./futureLeadingEntry);

        } else {
            det *= futureLeadingEntry;
            // Be careful, see that multiplier is indeed not a reference to futureLeadingEntry
            scaleRow(rowsSolved, 1./futureLeadingEntry);
        }

        // Now first entry in column has a non-zero value
        // Subtract all other rows
        for (int i = 0; i<this->nRows(); i++){
            if (i == rowsSolved) continue;
            const T1 &entryTo0 = this->getReference(i, j);
            if (entryTo0 == 0) continue;
            this->addRow(rowsSolved, i, -entryTo0);
        }

        // Successfully created a new leading row
        rowsSolved++;


    }

    // Return rank and determinant
    res[0][0] = rowsSolved;
    det = (rowsSolved == nRows() && rowsSolved == nCols()) ? det : T1(0.);
    return res;

}

// Make matrix column independent. Returns index of removed columns.
template<class T1>
vector<int> mat<T1>::makeColumnIdependent() {

    // Reduce
    //mat myCopy = *this;
    mat myCopy(nRows(), nCols());
    for (int i = 0; i<nRows(); i++){
        for (int j=0; j<nCols(); j++)
            myCopy[i,j] = (*this)[i,j];
    }
    myCopy.makeReducedEchelon();

    // Remove column index
    vector<int> colIdxRemoved(myCopy.nCols());

    // Remove all non-leading entries
    int leadingRowsFound = 0;
    for (int j=0; j<myCopy.nCols() and leadingRowsFound < myCopy.nRows(); j++){

        // Check current entry
        if (abs(myCopy.getReference(leadingRowsFound,j) -1.) < 0.00001)
            leadingRowsFound++;
        else {
            removeColumn(leadingRowsFound);
            colIdxRemoved[j-leadingRowsFound] = j;
        }
    }

    // Resize index vector to proper size
    colIdxRemoved.resize(myCopy.nCols()-leadingRowsFound);

    return colIdxRemoved;


}


// Gaussian determinant
template<class T1>
T1 mat<T1>::det() const {

    mat temp(*this);
    mat res = temp.makeEchelon();
    return res[1][0];

}

// Invert the matrix
template<class T1>
void mat<T1>::inv() {

    // Check for square matrix
    size_t myN = nRows();
    if (myN != nCols())
        throw runtime_error("Cannot invert rectangular matrix");

    // Append identity matrix
    mat submat = *this;
    submat.resize(myN, myN*2);
    for (int i=0; i< myN; i++)
        submat[i, myN+i] = 1;

    // Get reduced echelon Form
    mat res = submat.makeReducedEchelon();

    // Check if matrix is identity
    if (round(submat[myN-1,myN-1]*10000)/10000 != T1(1))
        throw runtime_error("Cannot invert singular matrix");

    // return inverse
    (*this) = submat.getSubMatrix(0, myN-1, myN, myN*2-1);

}

template<class T1>
void mat<T1>::addRow(const int &readRow, const int &writeRow, const T1& multiplier) {
    if (not isTransposed) {
        auto readStartConst = getIterator(readRow,0);;
        auto readEndConst = readStartConst + nCols();
        auto writeStart = getIterator(writeRow,0);

        transform(readStartConst,
                  readEndConst,
                  writeStart,
                  writeStart,
                  [&multiplier](const T1 &a,const T1 &b)
                            {return b + a*multiplier;});
    }
    else{
        for (int j=0; j< this->nCols(); j++)
            (*this)[writeRow, j] += (*this)[readRow, j]*multiplier;
    }

}

template<class T1>
void mat<T1>::resize(const int& myN, const int& myM) {

    if (isTransposed)
        throw runtime_error("There is no reason to do resizing in transposed form");

    if (myN<=0 or myM <=0)
        throw runtime_error("You cannot resize matrix to less than 1x1");

    int coldiff = M-myM;

    // Make sure we have enough memory
    if (coldiff<0) {
        data.resize(N * myM);

        // Do swap on all old rows
        for (int i=N-1; i>0; i--)
            swap_ranges(getIterator(i,0),
                        getIterator(i, M),
                        &data[i*myM]);

        data.resize(myN*myM);

    } else if (coldiff>0) {

        // Do shift on all old
        for (int i=1;i<N; i++)
            shift_left(getIterator(i,0), getIterator(i,M-1), coldiff*i);

        // Remove memory
        data.resize(myN*myM);

    }

    // Save size
    N = myN;
    M = myM;

}

template<class T1>
void mat<T1>::removeRow(const int &i) {

    if (isTransposed){
        transpose();
        removeColumn(i);
        transpose();
        return;
    }

    // Move end data over the row
    shift_left(getIterator(i,0), endMutable(), M);

    // Deallocate memory
    data.resize((N-1)*M);

    N--;

}

template<class T1>
void mat<T1>::removeColumn(const int &j) {

    if (isTransposed){
        transpose();
        removeRow(j);
        transpose();
        return;
    }

    // Shift all rows to overwrite column

    // Obvious implementation
    for (int i = 0; i<nRows(); i++){
        shift_left(getIterator(i,j)-i, endMutable(), 1);
    }

    /* Might be faster
    for (int i = 0; i<nRows()-1; i++){
        shift_left(getIterator(i,j-i), getIterator(i+1, j), i+1);
    }
    shift_left(getIterator(nRows()-1,j), endMutable(), nRows());*/


    // Deallocate memory
    data.resize(N*(M-1));

    M--;

}

template<class T1>
T1* mat<T1>::getPointer(const int &i, const int &j) {

    // Get vector idx
    size_t idx;
    if (isTransposed)
        idx = j*M+i;
    else
        idx = i*M+j;

    // Check that entry is inside matrix bounds.
    if (i*M+j > data.size()){
        throw std::runtime_error("The matrix has too few entries.");
    }

    // Return entry.
    return &data[idx];

}

template<class T1>
void mat<T1>::scaleRow(const int &rowidx, const T1 &multiplier) {

    if (isTransposed){
        auto ItEnd = getIterator(rowidx,nCols()-1);
        for (auto It = getIterator(rowidx,0); It<ItEnd; It++)
            *It *= multiplier;
    } else {
        transform(this->getIterator(rowidx,0),
                  this->getIterator(rowidx,M),
                  this->getIterator(rowidx,0),
                  [&multiplier](const T1 &a){
                    return a*multiplier;
                });

    }

}

template<class T1>
void mat<T1>::print() const {

    cout << endl;

    for (int i = 0; i< nRows(); i++){
        for (int j = 0 ; j<nCols(); j++){
            cout << (*this)[i,j] << ",";
        }
        cout << endl;
    }

}

template<class T1>
vector<T1>::iterator mat<T1>::getIterator(const int &i, const int &j) {

    // Get vector idx
    size_t idx;
    if (isTransposed)
        idx = j*M+i;
    else
        idx = i*M+j;

    // Check that entry is inside matrix bounds. The end point may be needed
    if (idx > data.size()){
        throw std::runtime_error("The matrix has too few entries.");
    }

    // Return entry.
    return data.begin() + idx;
}

template<class T1>
const vector<T1>::const_iterator mat<T1>::getConstIterator(const int &i, const int &j) const {

    // Get vector idx
    size_t idx;
    if (isTransposed)
        idx = j*M+i;
    else
        idx = i*M+j;

    // Check that entry is inside matrix bounds.
    if (i*M+j > data.size()){
        throw std::runtime_error("The matrix has too few entries.");
    }

    // Return entry.
    return data.begin() + idx;
}

template<class T1>
mat<T1> mat<T1>::getSubMatrix(const int &rowStart, const int &rowEnd,
                       const int &columnStart, const int &columnEnd) const {

    // Check that indexes valid
    const int rowdiff = rowEnd - rowStart;
    const int columndiff = columnEnd - columnStart;
    const int myRows = nRows();
    const int myCols = nCols();
    if (rowdiff<=0 or columndiff<=0 or
        rowStart<0 or columnStart<0 or
        rowEnd>=myRows or columnEnd >= myCols)
        throw runtime_error("Sub matrix is out of bounds");

    // Allocate result
    size_t resN = rowdiff +1;
    size_t resM = columndiff +1;
    mat res(resN, resM);

    if (isTransposed){

        // Fast copy columns
        for (int j = 0; j< resM; j++)
            copy_n(getConstIterator(rowStart, j), resN, res.getIterator(0,j));

    } else {

        // Fast copy rows
        for (int i =0; i<resN; i++)
            copy_n(getConstIterator(i, columnStart), resM, res.getIterator(i, 0));

    }

    return res;
}

template<class T1>
mat<T1> mat<T1>::getOrthonormalBasis() const {

    // Make copy
    mat<T1> myCopy = *this;
    myCopy.makeColumnIdependent();
    //myCopy.setEntries([](const num &a){return round(10000.*abs(a))/10000.;});

    // Allocate memory to transposed orthonormal basis: Reading rows is easier
    mat<T1> basis(myCopy.nCols(), myCopy.nRows());

    // Keep track of all the norms for each of the basis vectors
    vector<T1> norms(myCopy.nCols());

    // Allocate a temporary vector
    mat<T1> temp(1,myCopy.nRows());

    for (int i=0; i<myCopy.nCols(); i++){

        // Read column
        for (int j=0; j<myCopy.nRows();j++)
            basis[i,j] = myCopy[j,i];

        // Subtract projection on to prior columns
        for (int k=0; k< i; k++) {

            // Find Projection
            ::project(basis.getConstIterator(i, 0),
                    basis.getConstIterator(i, 0) +basis.nCols(),
                    basis.getConstIterator(k, 0),
                    temp.beginMutable());

            // Subtract projection
            transform(basis.getConstIterator(i,0),
                      basis.getConstIterator(i,0) + basis.nCols(),
                      temp.begin(),
                      basis.getIterator(i,0), std::minus<>());

        }

        // Find basis length
        T1 norm = inner_product(basis.getConstIterator(i,0),
                                    basis.getConstIterator(i,0)+basis.nCols(),
                                    basis.getConstIterator(i,0),
                                    T1(0.),std::plus<>(), std::multiplies<>());

        // Write down the norm
        norms[i] = sqrt(norm);

    }


    // Get the rank largest norm
    vector<T1> normsCopy = norms;
    sort(normsCopy.begin(), normsCopy.end());
    //T1 normLowerBound = normsCopy[normsCopy.size()-basis.nCols()];

    // Remove 0 basis vectors
    int rowsRemoved = 0;
    for (int i =0; i< basis.nRows(); i++){
        T1 norm = norms[i];
        int idx = i - rowsRemoved;
        if (norms[i] < 0.0001) {
            basis.removeRow(idx);
            rowsRemoved++;
        } else {
            // Scale basis
            transform(basis.getConstIterator(idx,0),
                      basis.getConstIterator(idx,0)+ basis.nCols(),
                      basis.getIterator(idx,0),
                      [&norm](const T1 &a){return a/norm;});
        }
    }

    // Transpose back basis matrix
    basis.transpose();

    return basis;
}

template<class T1>
mat<T1> mat<T1>::getOrthogonalProjectionMatrix() const {
    // This is a very unstable calculation

    // get Orthonormal basis
    mat<T1> basis = getOrthonormalBasis();

    // get Square inverse of basis
    mat<T1> sqrInv = basis.T().matmul(basis);
    if (sqrInv.det() == 0)
        throw runtime_error("Orthonormal basis is too close to be singular.");
    sqrInv.inv();

    // return Orthogonal projection matrix
    return basis.matmul(sqrInv).matmul(basis.T());

}

template<class T1>
mat<T1> mat<T1>::project(const mat<T1> &myVector) const {

    mat<T1> basis = getOrthonormalBasis();
    vector<T1> basisVector(basis.nRows());
    vector<T1> basisVectorProjection(basis.nRows());
    mat<T1> projectedVector(basis.nRows(), 1, T1(0.));

    for (int i=0; i<basis.nCols(); i++){

        // Make memory linear
        for (int j=0; j<basis.nRows(); j++)
            basisVector[j] = basis[j,i];

        // Get projection on to basis vector
        ::project(myVector.begin(), myVector.end(),
                basisVector.begin(),
                basisVectorProjection.begin());

        // Add projection
        transform(projectedVector.begin(), projectedVector.end(),
                  basisVectorProjection.begin(),
                  projectedVector.beginMutable(),
                  std::plus<>());

    }

    return projectedVector;
}

template<class T1>
void mat<T1>::setAdjointVal(const mat<adjointIntegral> &rhs) {

    // Check size
    if (nRows()!=rhs.nRows() or nCols()!=rhs.nCols())
        throw runtime_error("Dimensions dont match");

    auto rhsdata = rhs.begin();
    for (int i=0;i<data.size();i++)
        data[i] = rhs[i]->getVal();

}

template<class T1>
void mat<T1>::setAdjointAdj(const mat<adjointIntegral> &rhs) {

    // Check size
    if (nRows()!=rhs.nRows() or nCols()!=rhs.nCols())
        throw runtime_error("Dimensions dont match");

    auto rhsdata = rhs.begin();
    for (int i=0;i<data.size();i++)
        data[i] = rhs[i]->getAdjoint();

}

template class mat<float>;
template class mat<double>;
template class mat<int>;
template class mat<dual<float>>;
template class mat<dual<double>>;
template class mat<adjointIntegral>;


void matmultest() {
    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_int_distribution<std::mt19937::result_type> dist(0, 100); // distribution in range [1, 6]

    size_t N = 200;
    mat<num> temp1(N, N);
    for (int t = 0; t < 1000; t++) {


        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++)
                temp1[i, j] = 1. * dist(rng);
        }
        mat<num> tempold = temp1;
        temp1.inv();
        mat<num> tempnew = tempold.matmul(temp1);
        if (round(tempnew[3, 3]) != 1) {
            mat<num> myRes = tempold.makeReducedEchelon();
            cout << myRes[1][0] << endl;

        }

    }
    return;
}

void projectionMatrixTest(){

    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_int_distribution<std::mt19937::result_type> dist(0, 5); // distribution in range [1, 6]

    size_t N = 50;
    mat<num> temp1(N, 180);
    int errorCount = 0;
    for (int t = 0; t < 1000; t++) {


        for (int i = 0; i < temp1.nRows(); i++) {
            for (int j = 0; j < temp1.nCols(); j++)
                temp1[i, j] = 1. * dist(rng) ;
        }
        mat<num> projectionMatrix = temp1.getOrthogonalProjectionMatrix();
        mat<num> sameMat = projectionMatrix.matmul(projectionMatrix);
        mat<num> diff = projectionMatrix-sameMat;
        diff.setEntries([](const num &x){return abs(x);});
        num err = reduce(diff.begin(), diff.end())/temp1.nCols()/temp1.nRows();
        if (err> 0.0001) {
            errorCount++;
            cout << "Error " << errorCount << " out of " << t+1 << ": " << err << endl;

        }

    }
    return;

}


// Vector on Vector projection
void project(const auto vecStartIt, const auto vecEndIt,
             const auto baseStartIt, auto resStartIt) {

    // Find vector length
    size_t N = vecEndIt - vecStartIt;

    // Find projection factor
    remove_reference_t<decltype(*resStartIt)> startval = (*vecStartIt)*0.;
    auto angleScale = inner_product(vecStartIt, vecEndIt, baseStartIt, startval, std::plus<>(), std::multiplies<>());

    // Find vector squared norm
    startval = (*vecStartIt)*0.;
    auto squaredNorm = inner_product(baseStartIt, baseStartIt + N,
                                      baseStartIt, startval,
                                      std::plus<>(), std::multiplies<>());

    // If base is zero vector, then return zero vector
    if (squaredNorm == 0) {
        fill(resStartIt, resStartIt + N, 0.);
        return;
    }

    // Find scaling factor
    auto scale = angleScale/squaredNorm;

    // Scale vector
    transform(baseStartIt, baseStartIt + N, resStartIt,
              [&scale](const auto &a){return a*scale;});


}

