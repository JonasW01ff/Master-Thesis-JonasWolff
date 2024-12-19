//
// Created by Jonas Wolff on 17/11/2024.
// Inspired (strongly) from Antoine Savines AAD library

#ifndef CODE_ADJOINTNUMERIC_H
#define CODE_ADJOINTNUMERIC_H

#include "adjointNode.h"
#include "adjointTape.h"
#include <iostream>
#include "distributionFunctions.h"

// Expresions can have different number of arguments
template<class ARGTYPE>
class adjointExpr{
private:

public:

    const num& getVal() const {return static_cast<const ARGTYPE&>(*this).getVal();} ;

};

template<class LEXPR, class REXPR, class OP>
class binaryExpr: public adjointExpr<binaryExpr<LEXPR, REXPR, OP>>{
private:

    const num val;
    const LEXPR lexpr;
    const REXPR rexpr;

public:

    enum {nActives = LEXPR::nActives + REXPR::nActives};

    // We force it to only take expression derived types
    binaryExpr(const adjointExpr<LEXPR> &mylexpr, const adjointExpr<REXPR> &myrexpr):
            lexpr(static_cast<const LEXPR&>(mylexpr)),
            rexpr(static_cast<const REXPR&>(myrexpr)),
            val(OP::f(mylexpr.getVal(), myrexpr.getVal())){};

    const num& getVal() const {return val;};

    template<size_t idxActives>
    void propagatePartial(adjointNode* exprNode, const num adjoint) const{
        
        if (LEXPR::nActives > 0)
            lexpr.template propagatePartial<idxActives>(exprNode,
                                               adjoint * OP::fPartialL(lexpr.getVal(),
                                                        rexpr.getVal(), val));


        if (REXPR::nActives > 0)
            rexpr.template propagatePartial<idxActives+LEXPR::nActives>(exprNode,
                                                               OP::fPartialR(lexpr.getVal(),
                                                                        rexpr.getVal(), val) * adjoint);

    }

};

struct OpMult{

    static const num f(const num l, const num r){
        return l*r;
    }

    static const num fPartialL(const num l, const num r, const num v){
        return r;
    }

    static const num fPartialR(const num l, const num r, const num v){
        return l;
    }

};

struct OpAdd{

    static const num f(const num l, const num r){
        return l+r;
    }

    static const num fPartialL(const num l, const num r, const num v){
        return 1;
    }

    static const num fPartialR(const num l, const num r, const num v){
        return 1;
    }

};

struct OpSub{

    static const num f(const num l, const num r){
        return l-r;
    }

    static const num fPartialL(const num l, const num r, const num v){
        return 1;
    }

    static const num fPartialR(const num l, const num r, const num v){
        return -1;
    }

};

struct OpDiv{

    static const num f(const num l, const num r){
        return l/r;
    }

    static const num fPartialL(const num l, const num r, const num v){
        return 1./r;
    }

    static const num fPartialR(const num l, const num r, const num v){
        return -v/r; //   -l/(r*r);
    }

};

struct OpPow{

    static const num f(const num l, const num r){
        return pow(l, r);
    }

    static const num fPartialL(const num l, const num r, const num v){
        return r*v/l;
    }

    static const num fPartialR(const num l, const num r, const num v){
        return v*log(l);
    }

};

struct OpMax{

    static const num f(const num l, const num r){
        return max(l, r);
    }

    static const num fPartialL(const num l, const num r, const num v){
        return (l>r) ? 1 : 0;
    }

    static const num fPartialR(const num l, const num r, const num v){
        return (l>r)? 0 : 1;
    }

};

struct OpMin{

    static const num f(const num l, const num r){
        return min(l, r);
    }

    static const num fPartialL(const num l, const num r, const num v){
        return (l>r) ? 0 : 1;
    }

    static const num fPartialR(const num l, const num r, const num v){
        return (l>r)? 1 : 0;
    }

};

template<class LEXPR, class REXPR>
binaryExpr<LEXPR, REXPR, OpMult> operator*(const adjointExpr<LEXPR> &l, const adjointExpr<REXPR> &r){
    return binaryExpr<LEXPR, REXPR, OpMult>(l, r);
}

template<class LEXPR, class REXPR>
binaryExpr<LEXPR, REXPR, OpDiv> operator/(const adjointExpr<LEXPR> &l, const adjointExpr<REXPR> &r){
    return binaryExpr<LEXPR, REXPR, OpDiv>(l, r);
}

template<class LEXPR, class REXPR>
binaryExpr<LEXPR, REXPR, OpAdd> operator+(const adjointExpr<LEXPR> &l, const adjointExpr<REXPR> &r){
    return binaryExpr<LEXPR, REXPR, OpAdd>(l, r);
}

template<class LEXPR, class REXPR>
binaryExpr<LEXPR, REXPR, OpSub> operator-(const adjointExpr<LEXPR> &l, const adjointExpr<REXPR> &r){
    return binaryExpr<LEXPR, REXPR, OpSub>(l, r);
}

template<class LEXPR, class REXPR>
binaryExpr<LEXPR, REXPR, OpPow> pow(const adjointExpr<LEXPR> &l, const adjointExpr<REXPR> &r){
    return binaryExpr<LEXPR, REXPR, OpPow>(l, r);
}

template<class LEXPR, class REXPR>
binaryExpr<LEXPR, REXPR, OpMax> max(const adjointExpr<LEXPR> &l, const adjointExpr<REXPR> &r){
    return binaryExpr<LEXPR, REXPR, OpMax>(l, r);
}

template<class LEXPR, class REXPR>
binaryExpr<LEXPR, REXPR, OpMin> min(const adjointExpr<LEXPR> &l, const adjointExpr<REXPR> &r){
    return binaryExpr<LEXPR, REXPR, OpMin>(l, r);
}

template<class EXPR>
class quadratureExpr: public adjointExpr<quadratureExpr<EXPR>>{
private:

        const EXPR expr;
        const num quadratureValue;
        const num quadraturePartial;

public:

    enum {nActives=EXPR::nActives};

    quadratureExpr(const adjointExpr<EXPR> &myexpr, const num& quadVal, const num& quadPartial): expr(static_cast<const EXPR&>(myexpr)),
        quadratureValue(quadVal), quadraturePartial(quadPartial){};

    const num& getVal() const {return quadratureValue;};

    template<size_t idxActives>
    void propagatePartial(adjointNode* exprNode, const num& adjoint) const{

        // Let the user calculate the partial on their own as this is known at calc time.
        if (EXPR::nActives > 0)
            expr.template propagatePartial<idxActives>(exprNode,
                                                       adjoint * quadraturePartial);
    }

};

template<class EXPR>
quadratureExpr<EXPR> setupAdjointQuadrature(const adjointExpr<EXPR> &myexpr, const num &myQuadVal, const num& myQuadPartial){
    return quadratureExpr<EXPR>(myexpr, myQuadVal, myQuadPartial);
}




template<class EXPR, class OP>
class unaryExpr: public adjointExpr<unaryExpr<EXPR, OP>>{
private:

    num val;
    //adjointNode *node;
    const EXPR expr;
    const num binaryDoubleArg;

public:

    enum {nActives=EXPR::nActives};

    unaryExpr(const adjointExpr<EXPR> &myexpr): expr(static_cast<const EXPR&>(myexpr)), val(OP::f(myexpr.getVal(), 0.)),
                                                binaryDoubleArg(0.){};

    // Technically binary, but one of the arguments is double so it uses unary differentiation
    unaryExpr(const adjointExpr<EXPR> &myexpr, const num doublearg): expr(static_cast<const EXPR&>(myexpr)),
                val(OP::f(myexpr.getVal(), doublearg)), binaryDoubleArg(doublearg){};

    const num& getVal() const {return val;};

    template<size_t idxActives>
    void propagatePartial(adjointNode* exprNode, const num& adjoint) const{

        // We need the initiated double since it becomes part of the expression.
        if (EXPR::nActives > 0)
            expr.template propagatePartial<idxActives>(exprNode,
                                               adjoint * OP::fPartial(expr.getVal(), val, binaryDoubleArg));
    }



};

struct OpNormCDF{

    static const num f(const num x, const num c){
        return NormCDF(x);
    }

    static const num fPartial(const num x, const num y, const num c){
        return NormPDF(x);
    }

};

struct OpExp{

    static const num f(const num x, const num c){
        return exp(x);
    }

    static const num fPartial(const num x, const num y, const num c){
        return y;
    }


};

struct OpLog{

    static const num f(const num x, const num c){
        return log(x);
    }

    static const num fPartial(const num x, const num y, const num c){
        return 1./x;
    }

};

struct OpSqrt{

    static const num f(const num x, const num c){
        return sqrt(x);
    }

    static const num fPartial(const num x, const num y, const num c){
        return .5/y;
    }

};

struct OpAbs{

    static const num f(const num x, const num c){
        return abs(x);
    }

    static const num fPartial(const num x, const num y, const num c){
        return (x>=0) ? 1 : -1;
    }

};

struct OpRound{

    static const num f(const num x, const num c){
        return round(x);
    }

    static const num fPartial(const num x, const num y, const num c){
        return 1; // probably better for convergence
    }

};

struct OpMultConstant{

    static const num f(const num x, const num  c){
        return x*c;
    }
    static const num fPartial(const num x, const num y, const num c){
        return c;
    }

};

struct OpAddConstant{

    static const num f(const num x, const num c){
        return x+c;
    }

    static const num fPartial(const num x, const num y, const num c){
        return 1;
    }

};

struct OpSubLeftConstant{

    static const num f(const num x, const num c){
        return c - x;
    }

    static const num fPartial(const num x, const num y, const num c){
        return -1;
    }

};

struct OpSubRightConstant{

    static const num f(const num x, const num c){
        return x - c;
    }

    static const num fPartial(const num x, const num y, const num c){
        return 1;
    }

};

struct OpDivLeftConstant{

    static const num f(const num x, const num c){
        return c/x;
    }

    static const num fPartial(const num x, const num y, const num c){
        return -y/x;
    }

};

struct OpDivRightConstant{

    static const num f(const num x, const num c){
        return x/c;
    }

    static const num fPartial(const num x, const num y, const num c){
        return 1./c;
    }

};

struct OpPowLeftConstant{

    static const num f(const num x, const num c){
        return pow(c, x);
    }

    static const num fPartial(const num x, const num y, const num c){
        return y*log(c);
    }

};

struct OpPowRightConstant{

    static const num f(const num x, const num c){
        return pow(x,c);
    }

    static const num fPartial(const num x, const num y, const num c){
        return c*y/x;
    }

};

struct OpMaxConstant{

    static const num f(const num x, const num c){
        return max(x,c);
    }

    static const num fPartial(const num x, const num y, const num c){
        return (x>c) ? 1 : 0;
    }

};

struct OpMinConstant{

    static const num f(const num x, const num c){
        return min(c,x);
    }

    static const num fPartial(const num x, const num y, const num c){
        return (x<c) ? 1 : 0;
    }

};

template<class EXPR>
unaryExpr<EXPR, OpNormCDF> NormCDF(const adjointExpr<EXPR> &myexpr){
    return unaryExpr<EXPR, OpNormCDF>(myexpr);
}

template<class EXPR>
unaryExpr<EXPR, OpExp> exp(const adjointExpr<EXPR> &myexpr){
    return unaryExpr<EXPR, OpExp>(myexpr);
}

template<class EXPR>
unaryExpr<EXPR, OpLog> log(const adjointExpr<EXPR> &myexpr){
    return unaryExpr<EXPR, OpLog>(myexpr);
}

template<class EXPR>
unaryExpr<EXPR, OpSqrt> sqrt(const adjointExpr<EXPR> &myexpr){
    return unaryExpr<EXPR, OpSqrt>(myexpr);
}

template<class EXPR>
unaryExpr<EXPR, OpAbs> abs(const adjointExpr<EXPR> &myexpr){
    return unaryExpr<EXPR, OpAbs>(myexpr);
}

template<class EXPR>
unaryExpr<EXPR, OpRound> round(const adjointExpr<EXPR> &myexpr){
    return unaryExpr<EXPR, OpRound>(myexpr);
}

template<class EXPR>
unaryExpr<EXPR, OpMultConstant> operator*(const adjointExpr<EXPR> &myexpr, const num c){
    return unaryExpr<EXPR, OpMultConstant>(myexpr, c);
}

template<class EXPR>
unaryExpr<EXPR, OpMultConstant> operator*(const num c, const adjointExpr<EXPR> &myexpr){
    return unaryExpr<EXPR, OpMultConstant>(myexpr, c);
}

template<class EXPR>
unaryExpr<EXPR, OpAddConstant> operator+(const adjointExpr<EXPR> &myexpr, const num c){
    return unaryExpr<EXPR, OpAddConstant>(myexpr, c);
}

template<class EXPR>
unaryExpr<EXPR, OpAddConstant> operator+(const num c, const adjointExpr<EXPR> &myexpr){
    return unaryExpr<EXPR, OpAddConstant>(myexpr, c);
}

template<class EXPR>
unaryExpr<EXPR, OpSubLeftConstant> operator-(const num &c, const adjointExpr<EXPR> &myexpr){
    return unaryExpr<EXPR, OpSubLeftConstant>(myexpr, c);
}

template<class EXPR>
unaryExpr<EXPR, OpSubRightConstant> operator-(const adjointExpr<EXPR> &myexpr, const num &c){
    return unaryExpr<EXPR, OpSubRightConstant>(myexpr, c);
}

template<class EXPR>
unaryExpr<EXPR, OpDivLeftConstant> operator/(const num c, const adjointExpr<EXPR> &myexpr){
    return unaryExpr<EXPR, OpDivLeftConstant>(myexpr, c);
}

template<class EXPR>
unaryExpr<EXPR, OpDivRightConstant> operator/(const adjointExpr<EXPR> &myexpr, const num c){
    return unaryExpr<EXPR, OpDivRightConstant>(myexpr, c);
}

template<class EXPR>
unaryExpr<EXPR, OpPowLeftConstant> pow(const num &c, const adjointExpr<EXPR> &myexpr){
    return unaryExpr<EXPR, OpPowLeftConstant>(myexpr, c);
}

template<class EXPR>
unaryExpr<EXPR, OpPowRightConstant> pow(const adjointExpr<EXPR> &myexpr, const num c){
    return unaryExpr<EXPR, OpPowRightConstant>(myexpr, c);
}

template<class EXPR>
unaryExpr<EXPR, OpSubLeftConstant> operator-(const adjointExpr<EXPR> &rhs){
    return 0. -rhs;
}

// The leaf
class adjointIntegral: public adjointExpr<adjointIntegral>{
private:

    num val;
    adjointNode *node;

    static thread_local adjointTape* localTape;

    // Creat node from expression
    template<class ARGTYPE>
    void makeNode(const adjointExpr<ARGTYPE>& myExpr){

        node = localTape->record<ARGTYPE::nActives>();

        // Find expression partial with AAD (this does not get saved to the numbers adjoint, but to the partial)
        static_cast<const ARGTYPE&>(myExpr).template propagatePartial<0>(node, 1.);

    }

public:

    // This has to be compile time and counts how many active constants we have left,
    // so that we do not differentiate code that is passive.
    enum {nActives= 1};

    adjointIntegral(){
        val = 0.;
        node = localTape->record<0>();
        node->setIsLeaf(true);
    };

    adjointIntegral(bool &putontape){
        if (putontape){
            val = 0.;
            node = localTape->record<0>();
            node->setIsLeaf(true);
        }
    };

    explicit adjointIntegral(const num &myVal){
        val = myVal;
        node = localTape->record<0>();
        node->setIsLeaf(true);
    }

    void setupLeaf(){
        node = localTape->record<0>();
        node->setIsLeaf(true);
    }

    void setupLeaf(const num &myVal){
        val = myVal;
        node = localTape->record<0>();
        node->setIsLeaf(true);
    }

    // Check tape position
    const adjointNode* getNode(){
        return node;
    }

    void propagateAll(){

        adjointNode* currNode = node;
        localTape->poolNodes.deallocateToInc(currNode);
        while(currNode != nullptr){
            // Push adjoints
            currNode->pushAdjoints();

            // Clear the tape
            localTape->poolPartials.deallocateToInc(currNode->getPartials());
            localTape->poolArgAdjoints.deallocateToInc(currNode->getArgAdjoint());

            // Get prior Node in graph
            currNode = localTape->poolNodes.deallocateToInc(currNode);

        }
        if (currNode != nullptr){
            // Push adjoints
            currNode->pushAdjoints();
        }

        // Get prior Node in graph
        //currNode = localTape->poolNodes.deallocateToInc(currNode);
        cout << "Propagated"<< endl;

    }

    void setPause(){
        node->setPause(true);
    }

    void propagateToPause(){

        adjointNode* currNode = node;
        localTape->poolNodes.deallocateToInc(currNode);
        while(currNode != nullptr and !currNode->getPause()){
            // Push adjoints
            currNode->pushAdjoints();

            // Clear the tape
            localTape->poolPartials.deallocateToInc(currNode->getPartials());
            localTape->poolArgAdjoints.deallocateToInc(currNode->getArgAdjoint());

            // Get prior Node in graph
            currNode = localTape->poolNodes.deallocateToInc(currNode);

        }
        if (currNode != nullptr){
            // Push adjoints
            currNode->pushAdjoints();
        }

        // Get prior Node in graph
        //currNode = localTape->poolNodes.deallocateToInc(currNode);
        //cout << "Propagated"<< endl;

    }

    void deleteNodesToPause(){

        adjointNode* currNode = node;
        localTape->poolNodes.deallocateToInc(currNode);
        while(currNode != nullptr and !currNode->getPause()){
            // Push adjoints
            //currNode->pushAdjoints();

            // Clear the tape
            localTape->poolPartials.deallocateToInc(currNode->getPartials());
            localTape->poolArgAdjoints.deallocateToInc(currNode->getArgAdjoint());

            // Get prior Node in graph
            currNode = localTape->poolNodes.deallocateToInc(currNode);

        }
    }

    // Includes the pause
    void propagateToPauseInc(){

        adjointNode* currNode = node;
        localTape->poolNodes.deallocateToInc(currNode);
        while(currNode != nullptr and !currNode->getPause()){
            // Push adjoints
            currNode->pushAdjoints();

            // Clear the tape
            //localTape->poolPartials.deallocateToInc(currNode->getPartials());
            //localTape->poolArgAdjoints.deallocateToInc(currNode->getArgAdjoint());

            // Get prior Node in graph
            currNode = localTape->poolNodes.deallocateToInc(currNode);

        }
        if (currNode != nullptr){
            // Push adjoints
            currNode->pushAdjoints();

            // Clear the tape
            localTape->poolPartials.deallocateToInc(currNode->getPartials());
            localTape->poolArgAdjoints.deallocateToInc(currNode->getArgAdjoint());
            localTape->poolNodes.deallocateToInc(currNode);
        } else {
            localTape->poolPartials.deallocateToStart();
            localTape->poolArgAdjoints.deallocateToStart();
        }

        // Get prior Node in graph
        //currNode = localTape->poolNodes.deallocateToInc(currNode);
        cout << "Propagated"<< endl;

    }

    // When we try to propagate the expressions Partial derivative, then the end desitination is the
    // adjointIntegral class at thos point we should store the partial derivatives, which are currently stored
    // as adjoint down in to our partials.
    template<size_t idxActives>
    void propagatePartial(adjointNode* exprNode, const num adjoint) const{
        // Adjoint is still what it was before this flattend expresion, since we take the whole expresion
        // in one step
        exprNode->setArgAdjoint<idxActives>(&node->getAdjoint());
        //exprNode->setAdjoint<idxActives>(node->getAdjoint());

        // The partial is thus the adjoint we have calculated so far by setting t
        exprNode->setPartial<idxActives>(adjoint);
    };

    adjointIntegral& operator=(const num &myVal) {
        val = myVal;
        node = localTape->record<0>();
        node->setIsLeaf(true);
        return *this;
    };

    template<class EXPR>
    adjointIntegral(const adjointExpr<EXPR>& expr): val(expr.getVal()){
        makeNode<EXPR>(static_cast<const EXPR&>(expr));
    }

    template<class EXPR>
    adjointIntegral& operator=(const adjointExpr<EXPR>&expr){
        val = expr.getVal();
        makeNode<EXPR>(static_cast<const EXPR&>(expr));
        return *this;
    }

    explicit operator num () {return val;};
    explicit operator num& () {return val;};
    explicit operator string () {return to_string(val)+(getAdjoint()>=0 ? "+d":"-d")+ to_string(getAdjoint());};

    const num& getVal() const{return val;};
    void setVal(const num myVal){val = myVal;}

    const num& getAdjoint()const {return node->getAdjoint();};
    void setAdjoint(const num adj){node->setAdjoint(adj);};

    const num* getPartials() const{return node->getPartials();};
    template<size_t idx>
    void setPartial(const num& partial) {node->setPartial<idx>(partial);};

    /*adjointIntegral operator-() const{
        adjointIntegral temp(-this->getVal());
        temp.setAdjoint(-this->getAdjoint());
        return temp;
    }*/

    template<class ARGTYPE>
    adjointIntegral& operator+=(const adjointExpr<ARGTYPE> &expr){
        return *this = *this + expr;
    }

    template<class ARGTYPE>
    adjointIntegral& operator-=(const adjointExpr<ARGTYPE> &expr){
        return *this = *this - expr;
    }

    template<class ARGTYPE>
    adjointIntegral& operator*=(const adjointExpr<ARGTYPE> &expr){
        return *this = *this * expr;
    }

    template<class ARGTYPE>
    adjointIntegral& operator/=(const adjointExpr<ARGTYPE> &expr){
        return *this = *this / expr;
    }

    adjointIntegral& operator+=(const num& myVal){
        return *this = *this + myVal;
    }

    adjointIntegral& operator-=(const num& myVal){
        return *this = *this - myVal;
    }

    adjointIntegral& operator*=(const num& myVal){
        return *this = *this * myVal;
    }

    adjointIntegral& operator/=(const num& myVal){
        return *this = *this / myVal;
    }

    inline friend bool operator<(const adjointIntegral& LHS, const num& myVal){
        return LHS.getVal()<myVal;
    }

    bool operator<=(const num& myVal) const{
        return this->getVal()<myVal;
    }

    inline friend bool operator>(const adjointIntegral& LHS, const num& myVal){
        return LHS.getVal()>myVal;
    }

    bool operator>=(const num& myVal) const{
        return this->getVal()>myVal;
    }

    bool operator==(const num& myVal) const{
        return this->getVal()==myVal;
    }

    inline friend bool operator<(const adjointIntegral& LHS, const adjointIntegral& RHS){
        return LHS.getVal()<RHS.getVal();
    }

    bool operator<=(const adjointIntegral& myVal) const{
        return this->getVal()<myVal.getVal();
    }

    bool operator>(const adjointIntegral& myVal) const{
        return this->getVal()>myVal.getVal();
    }

    bool operator>=(const adjointIntegral& myVal) const{
        return this->getVal()>myVal.getVal();
    }

    bool operator==(const adjointIntegral& myVal) const{
        return this->getVal()==myVal.getVal();
    }


};

inline num to_num(const adjointIntegral& rhs){
    return rhs.getVal();
}

template<class fp, typename  = enable_if<is_floating_point_v<fp>>>
inline num to_num(const fp& rhs){
    return rhs;
}

// output stream operator
inline ostream& operator<<(ostream& os, const adjointIntegral& myAdj) {

    os << myAdj.getVal() << ((myAdj.getAdjoint() >= 0) ?  "+d" : "-d") << abs(myAdj.getAdjoint());
    return os;
}

inline string to_string(const adjointIntegral &rhs) {
    string res = to_string(rhs.getVal()) + ((rhs.getAdjoint() >= 0) ?  "+d" : "-d") + to_string(abs(rhs.getAdjoint()));
    return res;
}


#endif //CODE_ADJOINTNUMERIC_H
