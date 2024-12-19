//
// Created by Jonas Wolff on 24/09/2024.
//

#ifndef CODE_OPTIMIZER_H
#define CODE_OPTIMIZER_H
#include <cfloat>
#include <functional>
#include "matrix.h"
#include "type_traits"
#include "adjointNumeric.h"
#include <limits>

template<class T1, class T2>
void func_fd(
        function<void(const T1 &xBeginIt,const T1 &xEndIt, T2 &yBeginIt, T2 &yEndIt)> f,
        T1 &xBeginIt,
        T1 &xEndIt,
        const T2 &yBeginIt, // Pre-calculated y value at x
        const T2 &yEndIt,
        T2 &yDeltaBegin, // Row denominated transposed Jacobian matrix
        const num &sqrteps = sqrt(numeric_limits<num>::epsilon())
);

template<class TX>
void func_fd(
        function<void(const mat<TX> &x, mat<TX> &y)> &f,
        mat<TX> &x,
        mat<TX> &y,
        mat<TX> &jacobian,
        const num &sqrteps = sqrt(numeric_limits<num>::epsilon())){

    // Check dimensions
    if (x.nRows() != jacobian.nCols() or jacobian.nRows() != y.nRows())
        throw runtime_error("Incorrect dimensions");

    int jacobianCol = 0;
    mat<TX> tempY(y.size(),1);
    for (auto xIt = x.beginMutable(); xIt < x.end(); advance(xIt,1)){

        // Increment of x
        TX delta_x = (sqrteps*abs(*xIt)>sqrteps) ? sqrteps*abs(*xIt) : TX(sqrteps);

        // Check if we hit machine epsilon
        if (delta_x == TX(0.)) delta_x = sqrteps;

        // Save the old x
        TX myX = *xIt;

        // Increment x
        *xIt += delta_x;

        // Calculate increment
        f(x, tempY);

        // Calculate secant
        tempY = (tempY-y)/delta_x;

        for (int i=0; i<tempY.size(); i++)
            jacobian[i, jacobianCol] = tempY[i,0];

        // Go to next row
        jacobianCol++;

        // Reset x
        *xIt = myX;

    }

};

// Line search inspired by Numerical Recipies
template<class T>
void lnsrch(mat<T> &xold, const T &fold, mat<T> &grad, mat<T> &p,
            mat<T> &x, T &f, const T stpmax, bool &check, function<T(const mat<T> &X)> func){

    // Contract: xold, fold, grad and p (direction) has to be pre-calculated
    //          x, f and check just have to be allocated
    //          stpmax needs to be provided
    //          func needs to be provided

    const num alfa = 1.e-4, TOLX= numeric_limits<num>::epsilon();
    check = false;
    T sum = *max_element(p.begin(), p.end());

    // Cap step size
    if (sum>stpmax)
        transform(p.begin(), p.end(), p.beginMutable(),[&stpmax, &sum](const T &pi){return pi*stpmax/sum;});

    T slope = inner_product(grad.begin(), grad.end(),p.begin(), T(0.));

    if (slope >= 0.)
        p = -grad;//throw runtime_error("Roundoff error in lnsrch");

    T test(0.);
    for (auto xi = xold.begin(), pi = p.begin(); xi<xold.end(); ++xi, ++pi) {
        T temp = abs(*pi) / max(*xi, T(1.));
        test = (temp > test) ? temp : test;
    }
    T alamin = TOLX/test; // lambda_min
    T alam(1.0); // lambda_a

    // Do newtons steps
    T tempalam, alam2;
    T f2;
    mat<T> matL(2,2), matR(2,1), coef(2,1);
    int its = 0;
    while(true){
        its++;
        if (its>100) {
            // Check if step is too small
            if ((alam < alamin and alamin < 1.) or alam < numEPS*2*2){
                x = xold;
                check = true;
                return;
            }

            throw runtime_error("Too small lambda");
        }


        // Try full newton step (lambda = 1.)
        x = xold +p*alam;
        f = func(x);

        // Do a sanity check
        if (f > fold or f!=f){ // f!=f means NaN
            alam *= .5;
            continue;
        }

        // Check if step is too small
        if ((alam < alamin and alamin < 1.) or alam < numEPS*2*2){
            x = xold;
            check = true;
            return;
        } else if (f <= fold+alfa*alam*slope)
            return; // Check if functions is decreasing enough, and end this step if true.
        else {
            if (alam == 1.)
                tempalam = -slope/(2.*(f - fold- slope));
            else {
                // First step has been calculated to alam2 exists
                // Calculate a and b in 3rd degree polynomial
                matR[0,0] = f - slope*alam - fold;
                matR[0,1] = f2 - slope*alam2 - fold;
                matL[0,0] = 1./(alam*alam);
                matL[0,1] = -1./(alam2*alam2);
                matL[1,0] = -alam2*matL[0,0];
                matL[1,1] = -alam*matL[0,1];
                coef = matL.matmul(matR)/(alam - alam2);
                if(coef[0,0] == 0) // g(1) cancels out g'(0) and g(0)
                    tempalam = -slope/(2.*coef[1,0]);
                else {
                    // Discriminant
                    T d = coef[1,0]*coef[1,0]-3.0*coef[0,0]*slope;

                    // If no minimum exists in approximation, go max distance
                    // n.b. it starts decreasing so it also end decreasing
                    if (d<0) tempalam = .5*alam;
                    else if (coef[1,0] <= 0) tempalam = (-coef[1,0]+sqrt(d))/(3.*coef[0,0]);
                    else tempalam = -slope/(coef[1,0]+sqrt(d));
                    // In the last formula we rearrange as to avoid numbers of similar canceling out
                }
                if (tempalam > .5*alam)
                    tempalam = .5*alam;
            }
        }

        alam2 = alam;
        f2 = f;
        alam = (tempalam == tempalam) ? T(max(tempalam, T(.1)*alam)) : T(.1)*alam;

    }

}


// BFGS minimization of function inspired by Numerical Recipes
template<class T>
void dfpmin(mat<T> &p, const num gtol, int &its, T &fret, function<T(const mat<T>&)> &&func,
            function<void(mat<T>&, const T&, mat<T>&)> &&gradfunc, const int ITMAX =2000){

    // gradfunc takes x, function value at x, and places the gradient in the lat argument
    const num EPS= numEPS; // I intent to use variations of doubles
    const num TOLX=4.*EPS,STPMAX=100.;
    bool check;
    T den,fac,fad,fae,fp,stpmax,sum=T(0.),sumdg,sumxi,temp,test;
    int n=p.size();
    mat<T> dg(n,1), g(n,1), hdg(n,1), pnew(n,1), xi(n,1);
    mat<T> hessin(n,n, T(0.));
    mat<T> Bx(n,1), BxxB(n,n), r(n,1);
    T myTempSlope;

    fp=func(p);
    gradfunc(p,fp,g);
    for (auto hessinIt = hessin.beginMutable(); hessinIt<hessin.end(); hessinIt += hessin.nCols()+1)
        *hessinIt = 1.;
    xi = -g;
    stpmax = STPMAX*max(p.norm(), T(n));
    for (its=0;its<ITMAX;its++){
        //cout << its << ": " << fp << endl;
        //p.print();
        /*myTempSlope = inner_product(xi.begin(), xi.end(), g.begin(), T(0.));
        if (myTempSlope > 0)
            xi = -g;*/
        lnsrch(p,fp,g,xi,pnew,fret,stpmax, check, func);
        fp=fret;
        xi = pnew -p;
        p = pnew;

        // Convergence test for x step-size
        static auto myreducefunc = [](const T &myMax, const T &rhs){return max(myMax, abs(rhs));};
        test = reduce(xi.begin(), xi.end(), T(0.), myreducefunc);
        if (test<TOLX) {
            //cout << "xtol hit: " << test  << endl;
            return;
        }

        dg = g;

        // update gradient
        gradfunc(p,fp,g);

        // Test for small gradient
        test = reduce(g.begin(), g.end(), T(0.), myreducefunc);
        //test = g.norm()/g.size();
        if (test<gtol) {
            //cout << "gtol hit" << endl;
            return;
        }

        // get deltas
        dg = g-dg;
        hdg = hessin.matmul(dg);
        fac = dg.T().matmul(xi)[0,0];
        fae = dg.T().matmul(hdg)[0,0];
        sumdg = pow(dg.norm(),2.);
        sumxi = pow(xi.norm(),2.);
        // Curvature conditions
        if (fac > sqrt(EPS*sumdg*sumxi) and fac > numSQRTEPS){
            // These should not explode in size
            fac = 1./fac;
            fae = 1./fae;

            // u vector
            dg = xi*fac - hdg*fae;

            // update formula
            hessin += xi.matmul(xi.T())*fac;
            hessin -= hdg.matmul(hdg.T())*fae;
            hessin += dg.matmul(dg.T())*fae;

        } else {
            // Damped BFGS
            T theta(1);
            T xBx = xi.T().matmul(hessin).matmul(xi)[0,0];
            if (fac < 0.2*xBx){
                theta = 0.8*xBx/(xBx-fac);
            }
            Bx = hessin.matmul(xi);
            BxxB = Bx.matmul(Bx.T());
            hessin = hessin +(-BxxB/xBx);
            r = dg*theta + hessin.matmul(xi)*(1.-theta);
            hessin = hessin + (r.matmul(r.T()))/(xi.T().matmul(r)[0,0]);

        }

        if (xi.T().matmul(g)[0,0] >= 0){
            fill(hessin.beginMutable(), hessin.endMutable(), T(0.));
            for (auto hessinIt = hessin.beginMutable(); hessinIt<hessin.end(); hessinIt += hessin.nCols()+1)
                *hessinIt = T(1.);
        }


        // update direction
        xi = -hessin.matmul(g);

    }

    // Loop should never finish
    //throw runtime_error("too many iterations in bfgs: in function dfpmin");

}

// BFGS with AAD optimization minimization of function inspired by Numerical Recipes
inline void dfpminAAD(mat<num> &ptemp, const num gtol, int &its, num &fret,
               function<adjointIntegral(mat<adjointIntegral>&)> func,
               function<num(const mat<num>&)> funcD, const int ITMAX =2000){

    // gradfunc takes x, function value at x, and places the gradient in the lat argument
    const num EPS= numEPS; // I intent to use variations of doubles
    const num TOLX=4.*EPS,STPMAX=100.;
    bool check;
    num den,fac,fad,fae,sum,sumdg,sumxi,temp;
    num test, stpmax;
    int n=ptemp.size();
    mat<num> dg(n,1), hdg(n,1), pnew(n,1), xi(n,1);
    mat<num> g(n,1);
    adjointIntegral fp;
    mat<num> hessin(n,n);
    mat<num> Bx(n,1), BxxB(n,n), r(n,1);
    num myTempSlope;
    const adjointNode* myStartNode;
    adjointIntegral* myStartNum;
    mat<adjointIntegral> p(ptemp.nRows(), ptemp.nCols());
    for (int i=0; i<ptemp.nRows();i++){
        for (int j=0; j<ptemp.nCols();j++){
            p[i,j] = ptemp[i,j];
        }
    }
    myStartNode = p.beginMutable()->getNode();
    for (auto pi=p.beginMutable();pi<p.end();pi++){
        if (pi->getNode()<=myStartNode) {
            myStartNode = pi->getNode();
            myStartNum = &*pi;
        }
    }
    myStartNum->setPause();
    fp=func(p);
    fp.setAdjoint(1);
    fp.propagateToPauseInc();
    g.setAdjointAdj(p);

    for (auto hessinIt = hessin.beginMutable(); hessinIt<hessin.end(); hessinIt += hessin.nCols()+1)
        *hessinIt = 1.;

    xi = -g;
    stpmax = STPMAX*max(p.norm().getVal(), num(n));
    for (its=0;its<ITMAX;its++){
        //cout << its << ": " << fp.getVal() << endl;
        /*myTempSlope = inner_product(xi.begin(), xi.end(), g.begin(), T(0.));
        if (myTempSlope > 0)
            xi = -g;*/
        ptemp.setAdjointVal(p);
        lnsrch(ptemp,fp.getVal(),g,xi,pnew,fret,stpmax, check, funcD);
        fp=fret;
        xi = pnew -ptemp;
        p = pnew;

        // Convergence test for x step-size
        test = reduce(xi.begin(), xi.end(), num(0.), [](const num &myMax, const num &rhs)
                    {return max(myMax, abs(rhs));});
        if (test<TOLX) {
            cout << "xtol hit: " << test  << endl;
            return;
        }

        dg = g;

        for (auto pi=p.beginMutable();pi<p.end();pi++){
            if (pi->getNode()<=myStartNode) {
                myStartNode = pi->getNode();
                myStartNum = &*pi;
                myStartNum->setAdjoint(num(0.));
            }
        }
        myStartNum->setPause();
        fp=func(p);
        fp.setAdjoint(1);
        fp.propagateToPauseInc();
        g.setAdjointAdj(p);

        // Test for small gradient
        test = reduce(g.begin(), g.end(), num(0.), [](const num &myMax, const num &rhs)
                    {return max(myMax, abs(rhs));});
        //test = g.norm()/g.size();
        if (test<gtol) {
            cout << "gtol hit" << endl;
            return;
        }

        // get deltas
        dg = g-dg;
        hdg = hessin.matmul(dg);
        fac = dg.T().matmul(xi)[0,0];
        fae = dg.T().matmul(hdg)[0,0];
        sumdg = dg.norm();
        sumdg *= sumdg;
        sumxi = xi.norm();
        sumxi *= sumxi;
        // Curvature conditions
        if (fac > sqrt(EPS*sumdg*sumxi) and fac > numSQRTEPS){
            // These should not explode in size
            fac = 1./fac;
            fae = 1./fae;

            // u vector
            dg = xi*fac - hdg*fae;

            // update formula
            hessin += xi.matmul(xi.T())*fac;
            hessin -= hdg.matmul(hdg.T())*fae;
            hessin += dg.matmul(dg.T())*fae;

        } else {
            // Damped BFGS
            num theta = 1;
            num xBx = xi.T().matmul(hessin).matmul(xi)[0,0];
            if (fac < 0.2*xBx){
                theta = 0.8*xBx/(xBx-fac);
            }
            Bx = hessin.matmul(xi);
            BxxB = Bx.matmul(Bx.T());
            hessin = hessin +(-BxxB/xBx);
            r = dg*theta + hessin.matmul(xi)*(1.-theta);
            hessin = hessin + (r.matmul(r.T()))/(xi.T().matmul(r)[0,0]);

        }

        if (xi.T().matmul(g)[0,0] >= 0){
            for (auto hessinIt = hessin.beginMutable(); hessinIt<hessin.end(); hessinIt += hessin.nCols()+1)
                *hessinIt = 1.;
        }


        // update direction
        xi = -hessin.matmul(g);

    }

    // Loop should never finish
    //throw runtime_error("too many iterations in bfgs: in function dfpmin");

}

template<class T, typename = enable_if_t<is_floating_point_v<T>>>
inline T bisection(function<T(const T&)> &&f, T low, T high, T TSQRTEPS = numeric_limits<T>::epsilon()){

    // Calculate end points
    T flow = f(low);
    T fhigh = f(high);

    // Check order
    if (flow > fhigh){
        swap(fhigh, flow);
        swap(high, low);
    }

    // Check 0 is in interval
    if (fhigh<0 or flow>0){
        throw runtime_error("0 is not in inteval");
    }

    // Initialize
    T mid = (high+low)/2.;
    T fmid;

    int idx = 0;
    while (abs(high - low) > TSQRTEPS) {

        // Max iterations
        if (idx>1e6)
            throw runtime_error("Bisection took 1e6 iterations. Terminating process.");
        idx++;

        // New point
        mid = (low + high) / 2.;
        fmid = f(mid);

        // Bisect
        if (!signbit(fmid)) {
            swap(high, mid);
            swap(fhigh, fmid);
        } else {
            swap(low, mid);
            swap(flow, fmid);
        }

        // Round off error
        if (signbit(fhigh) == signbit(flow))
            throw runtime_error("Round off error in bisection.");

    }

    return mid;
}

inline adjointIntegral bisection(function<adjointIntegral(const adjointIntegral&)> &f, adjointIntegral low,
                                 adjointIntegral high, num TSQRTEPS = numeric_limits<num>::epsilon()){

    // Calculate end points
    adjointIntegral flow = f(low);
    adjointIntegral fhigh = f(high);

    // Check order
    if (flow > fhigh){
        swap(fhigh, flow);
        swap(high, low);
    }

    // Check 0 is in interval
    if (fhigh<0 or flow>0){
        throw runtime_error("0 is not in inteval");
    }

    // Initialize
    adjointIntegral mid = (high+low)/2.;
    adjointIntegral fmid;

    int idx = 0;
    while (abs(fhigh - flow) > TSQRTEPS) {

        // Max iterations
        if (idx>1e6)
            throw runtime_error("Bisection took 1e6 iterations. Terminating process.");
        idx++;

        // New point
        mid = (low + high) / 2.;
        fmid = f(mid);

        // Bisect
        if (fmid > 0) {
            high = mid;
            fhigh = fmid;
        } else if (fmid < 0) {
            low = mid;
            flow = fmid;
        } else {
            break;
        }

        // Round off error
        if (abs(fhigh.getVal()) / (fhigh.getVal()) == abs(flow.getVal()) / (flow.getVal()))
            throw runtime_error("Round off error in bisection.");

    }

    return mid;
}

template<class T1>
class optimizer {

    vector<T1> lastoptimum;

public:

    virtual void minimize(function<T1(const mat<T1> &x)> &&myFunc,
                          mat<T1> &x0, const num tol) = 0;

    // maximization by minimization
    T1 maximize(function<T1(typename vector<T1>::iterator,
                            typename vector<T1>::const_iterator)> myFunc,
                typename vector<T1>::iterator xStartIt,
                typename vector<T1>::const_iterator xEndIt){
        return minimize([&myFunc](auto startIt, auto endIt){
            return -myFunc(startIt, endIt);
        }, xStartIt, xEndIt);
    };


};

template<class T1>
class bfgs : optimizer<T1>{

public:

    bfgs(){};
    ~bfgs(){};

    void minimize(function<T1(const mat<T1>&)> &&myFunc,
                mat<T1> &x, const num tol = numEPS*2*2){


        // give myFunc the right form
        function<void(const mat<T1> &x, mat<T1> &y)> myFdFunc = [&myFunc](const mat<T1> &x, mat<T1> &y)->void{
            y[0,0] = myFunc(x);
        };

        // Get gradient function by finite difference
        function<void(mat<T1>&, const T1&, mat<T1>&)> gradfunc =
                [&myFdFunc](mat<T1>& x, const T1& y, mat<T1>&grad){
            mat<T1> yMat(1,1);
            yMat[0,0] = y;
            grad.transpose();
            func_fd(myFdFunc,x,yMat,grad);
            grad.transpose();
        };

        T1 ysol;
        int iter;
        dfpmin(x, tol, iter, ysol,move(myFunc), move(gradfunc));

    };

};

// Simulated annealing
/*template<class T1>
class heat : optimizer<T1>{

    T1 minimize(auto startIt, auto endIt);

};*/

#endif //CODE_OPTIMIZER_H
