#pragma once
 

#include "./mean.h"
#include "Eigen/Dense"

namespace Ponca
{

    /*!
        \brief Point Set Surfaces fitting algorithm

        This class implements the Point Set Surfaces algorithm described in part 1 of \cite Alexa:2009:Hermite
    */
 
template < class DataPoint, class _WFunctor, typename T >
class PSSPrim : public T
{
PONCA_FITTING_DECLARE_DEFAULT_TYPES
PONCA_FITTING_DECLARE_MATRIX_TYPE
 
protected:
    enum { Check = Base::PROVIDES_MEAN_POSITION && Base::PROVIDES_MEAN_NORMAL }; // && Base::MEAN_POSITION_DERIVATIVE ?
 
public:
    {
     
    PONCA_MULTIARCH inline Scalar potential (const VectorType& _x) const
    {
        // The potential is the distance from the point to the surface : f(x) = n(x)^T(x − c(x)) = 0 
        // as equation 4 given in \cite Alexa:2009:Hermite    
        return (((Base::m_sumN/Base::m_sumW).transpose()).dot( (_x - Base::barycenter())));
    }   

 
    PONCA_MULTIARCH inline VectorType gradient (const VectorType& _x) const
    {
        // The gradient of the potential is the normal at the point : \nabla f (x) = Jn(x)(c(x) − x) + (Jc(x) − I)n(x) with J the Jacobian matrix 
        // as equation 5 given in \cite Alexa:2009:Hermite

        Jn = // \todo
        Jc = // \todo Potentially barycenterDerivatives() ?
        return  Jn*(Base::barycenter() - _x) + (Jc - MatrixType::Identity())*(Base::m_sumN/Base::m_sumW);
    }
    }
};
} //namespace Ponca