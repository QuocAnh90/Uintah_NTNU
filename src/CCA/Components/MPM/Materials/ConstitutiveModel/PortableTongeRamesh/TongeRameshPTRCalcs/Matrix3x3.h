/*
 * This project constitutes a work of the United States Government and is not
 * subject to domestic copyright protection under 17 USC § 105.
 *
 * However, because the project utilizes code licensed from contributors and other
 * third parties, it therefore is licensed under the MIT License.
 * http://opensource.org/licenses/mit-license.php.
 *
 * Under that license, permission is granted free of charge, to any
 * person obtaining a copy of this software and associated
 * documentation files (the "Software"), to deal in the Software
 * without restriction, including without limitation the rights to
 * use, copy, modify, merge, publish, distribute, sublicense, and/or
 * sell copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the conditions that any
 * appropriate copyright notices and this permission notice are
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
 * IN THE SOFTWARE.
 * 
 * Matrix3x3.h -
 *
 * This file is based on a file by a similar name in the Uintah 
 * project (uintah.utah.edu). 
 * 
 */


#ifndef MATRIX3X3_H
#define MATRIX3X3_H

#ifndef EQN_EPS
#define EQN_EPS 1e-12
#endif

#ifndef IsZero
#define IsZero( x ) ((x) > -EQN_EPS && (x) < EQN_EPS)
#endif

#include <ostream>
#include <cmath>                // sqrt

class Matrix3x3
{
        public:
                Matrix3x3();
                Matrix3x3(const double val);
                Matrix3x3(const bool isIdentity);
                Matrix3x3(      const double v0,
                                        const double v1,
                                        const double v2,
                                        const double v3,
                                        const double v4,
                                        const double v5,
                                        const double v6,
                                        const double v7,
                                        const double v8);

                virtual ~Matrix3x3();
                
                void identity();
                Matrix3x3 transpose() const;

                double determinant() const;
                double trace() const;
                double normSquared() const;
                inline double norm() const
                {
                        return sqrt(normSquared()); 
                }

                const Matrix3x3 inverse() const;
                double get(const int i, const int j) const;
                void set(const int i, const int j, const double val);
                
                void swap(Matrix3x3 *rhs);
                
                Matrix3x3 operator+= (const Matrix3x3 rhs);
                const Matrix3x3 operator+(const Matrix3x3 rhs) const;
                Matrix3x3 operator-= (const Matrix3x3 rhs);
                const Matrix3x3 operator-(const Matrix3x3 rhs) const;
                Matrix3x3 operator*= (const double rhs);
                const Matrix3x3 operator*(const double rhs) const;
                Matrix3x3 operator/= (const double rhs);
                const Matrix3x3 operator/(const double rhs) const;
                const Matrix3x3 operator*(const Matrix3x3 rhs) const;
                int getEigenValues(double *e1, double *e2, double *e3) const;
  
                inline double Contract(const Matrix3x3 mat) const
                {
                        // Return the contraction of this matrix with another 
                        double contract = 0.0;
                        for (int i = 0; i< 3; i++) {
                                for(int j=0;j<3;j++){
                                        contract += get(i, j)*(mat.get(i, j));
                                }
                        }
                        return contract;
                }

        void polarRotationRMB(Matrix3x3 *R) const;
        void polarDecompositionRMB(Matrix3x3 *U, Matrix3x3 *R) const;


        template <class T> void singularValueDecomposition(Matrix3x3 *U, Matrix3x3 *Singma, Matrix3x3 *V) const;
        template <class T> class GivensRotation {
        public:
            int rowi;
            int rowk;
            T c;
            T s;

            inline GivensRotation(int rowi_in, int rowk_in)
                : rowi(rowi_in), rowk(rowk_in), c(1), s(0) {}

            inline GivensRotation(T a, T b, int rowi_in, int rowk_in)
                : rowi(rowi_in), rowk(rowk_in) {
                compute(a, b);
            }

            ~GivensRotation() {}

            inline void setIdentity() {
                c = 1;
                s = 0;
            }

            inline void transposeInPlace() { s = -s; }

            /**
                Compute c and s from a and b so that
                ( c -s ) ( a )  =  ( * )
                s  c     b       ( 0 )
                */
            inline void compute(const T a, const T b) {
                using std::sqrt;

                T d = a * a + b * b;
                c = 1;
                s = 0;
                T sqrtd = sqrt(d);
                // T t = MATH_TOOLS::rsqrt(d);
                if (sqrtd) {
                    T t = 1 / sqrtd;
                    c = a * t;
                    s = -b * t;
                }
            }

            /**
                This function computes c and s so that
                ( c -s ) ( a )  =  ( 0 )
                s  c     b       ( * )
                */
            inline void computeUnconventional(const T a, const T b) {
                using std::sqrt;

                T d = a * a + b * b;
                c = 0;
                s = 1;
                T sqrtd = sqrt(d);
                // T t = MATH_TOOLS::rsqrt(d);
                if (sqrtd) {
                    T t = 1 / sqrtd;
                    s = a * t;
                    c = b * t;
                }
            }

            /**
              Fill the R with the entries of this rotation
                */
            template <class Matrix3x3>
            inline void fill(const Matrix3x3& R) const {
                Matrix3x3& A = const_cast<Matrix3x3&>(R);
                A.identity();
                A(rowi, rowi) = c;
                A(rowk, rowi) = -s;
                A(rowi, rowk) = s;
                A(rowk, rowk) = c;
            }

            /**
                This function does something like Q^T A -> A
                [ c -s  0 ]
                [ s  c  0 ] A -> A
                [ 0  0  1 ]
                It only affects row i and row k of A.
                */
            template <class Matrix3x3>
            inline void rowRotation(Matrix3x3& A) const {
                for (int j = 0; j < A.cols(); j++) {
                    T tau1 = A(rowi, j);
                    T tau2 = A(rowk, j);
                    A(rowi, j) = c * tau1 - s * tau2;
                    A(rowk, j) = s * tau1 + c * tau2;
                }
                // not type safe :/
            }

            /**
                This function does something like A Q -> A
                   [ c  s  0 ]
                A  [-s  c  0 ]  -> A
                   [ 0  0  1 ]
                It only affects column i and column k of A.
                */
            template <class Matrix3x3>
            inline void columnRotation(Matrix3x3& A) const {
                for (int j = 0; j < A.rows(); j++) {
                    T tau1 = A(j, rowi);
                    T tau2 = A(j, rowk);
                    A(j, rowi) = c * tau1 - s * tau2;
                    A(j, rowk) = s * tau1 + c * tau2;
                }
                // not type safe :/
            }

            /**
              Multiply givens must be for same row and column
              **/
            inline void operator*=(const GivensRotation<T>& A) {
                T new_c = c * A.c - s * A.s;
                T new_s = s * A.c + c * A.s;
                c = new_c;
                s = new_s;
            }

            /**
              Multiply givens must be for same row and column
              **/
            inline GivensRotation<T>
                operator*(const GivensRotation<T>& A) const {
                GivensRotation<T> r(*this);
                r *= A;
                return r;
            }
        };

        /**
            \brief zero chasing the 3X3 matrix to bidiagonal form
            original form of H:
            x x 0
            x x x
            0 0 x
            after zero chase:
            x x 0
            0 x x
            0 0 x
            */


        template <class T>
        inline void zeroChase(Matrix3x3& H, Matrix3x3& U,
            Matrix3x3& V) {

            /**
               Reduce H to of form
               x x +
               0 x x
               0 0 x
               */
            GivensRotation<T> r1(H._values[0], H._values[3], 0, 1);
            /**
                Reduce H to of form
                x x 0
                0 x x
                0 + x
                Can calculate r2 without multiplying by r1 since both entries are in first
               two rows thus no need to divide by sqrt(a^2+b^2)
                */
            GivensRotation<T> r2(1, 2);
            if (H._values[3] != 0)
                r2.compute(H._values[0] * H._values[1] + H._values[3] * H._values[4],
                    H._values[0] * H._values[2] + H._values[3] * H._values[5]);
            else
                r2.compute(H._values[1], H._values[2]);

            r1.rowRotation(H);

            /* GivensRotation<T> r2(H._values[1], H._values[2], 1, 2); */
            r2.columnRotation(H);
            r2.columnRotation(V);

            /**
                Reduce H to of form
                x x 0
                0 x x
                0 0 x
                */
            GivensRotation<T> r3(H._values[4], H._values[7], 1, 2);
            r3.rowRotation(H);

            // Save this till end for better cache coherency
            // r1.rowRotation(u_transpose);
            // r3.rowRotation(u_transpose);
            r1.columnRotation(U);
            r3.columnRotation(U);

        }

        /**
             \brief make a 3X3 matrix to upper bidiagonal form
             original form of H:   x x x
                                   x x x
                                   x x x
             after zero chase:
                                   x x 0
                                   0 x x
                                   0 0 x
          */
        template <class T>
        inline void
            makeUpperBidiag(Matrix3x3& H, Matrix3x3& U, Matrix3x3& V) {
            U.identity();
            V.identity();

            /**
              Reduce H to of form
                                  x x x
                                  x x x
                                  0 x x
            */

            GivensRotation<T> r(H._values[3], H._values[6], 1, 2);
            r.rowRotation(H);
            // r.rowRotation(u_transpose);
            r.columnRotation(U);
            // zeroChase(H, u_transpose, V);
            //zeroChase(H, U, V);

            // Copy zeroChase down here

             /**
                Reduce H to of form
                x x +
                0 x x
                0 0 x
                */
            GivensRotation<T> r1(H._values[0], H._values[3], 0, 1);
            /**
                Reduce H to of form
                x x 0
                0 x x
                0 + x
                Can calculate r2 without multiplying by r1 since both entries are in first
               two rows thus no need to divide by sqrt(a^2+b^2)
                */
            GivensRotation<T> r2(1, 2);
            if (H._values[3] != 0)
                r2.compute(H._values[0] * H._values[1] + H._values[3] * H._values[4],
                    H._values[0] * H._values[2] + H._values[3] * H._values[5]);
            else
                r2.compute(H._values[1], H._values[2]);

            r1.rowRotation(H);

            /* GivensRotation<T> r2(H._values[1], H._values[2], 1, 2); */
            r2.columnRotation(H);
            r2.columnRotation(V);

            /**
                Reduce H to of form
                x x 0
                0 x x
                0 0 x
                */
            GivensRotation<T> r3(H._values[4], H._values[7], 1, 2);
            r3.rowRotation(H);

            // Save this till end for better cache coherency
            // r1.rowRotation(u_transpose);
            // r3.rowRotation(u_transpose);
            r1.columnRotation(U);
            r3.columnRotation(U);
        }


        /**
             \brief make a 3X3 matrix to lambda shape
             original form of H:   x x x
             *                     x x x
             *                     x x x
             after :
             *                     x 0 0
             *                     x x 0
             *                     x 0 x
          */
        template <class T>
        inline void
            makeLambdaShape(Matrix3x3& H, Matrix3x3& U, Matrix3x3& V) {
            U.identity();
            V.identity();

            /**
              Reduce H to of form
              *                    x x 0
              *                    x x x
              *                    x x x
              */

            GivensRotation<T> r1(H._values[1], H._values[2], 1, 2);
            r1.columnRotation(H);
            r1.columnRotation(V);

            /**
              Reduce H to of form
              *                    x x 0
              *                    x x 0
              *                    x x x
              */

            r1.computeUnconventional(H._values[5], H._values[8]);
            r1.rowRotation(H);
            r1.columnRotation(U);

            /**
              Reduce H to of form
              *                    x x 0
              *                    x x 0
              *                    x 0 x
              */

            GivensRotation<T> r2(H._values[6], H._values[7], 0, 1);
            r2.columnRotation(H);
            r2.columnRotation(V);

            /**
              Reduce H to of form
              *                    x 0 0
              *                    x x 0
              *                    x 0 x
              */
            r2.computeUnconventional(H._values[1], H._values[4]);
            r2.rowRotation(H);
            r2.columnRotation(U);
        }


        private:
                double _values[9];
                /* add your private declarations */
};

std::ostream& operator<<(std::ostream& out, const Matrix3x3 rhs);

int SolveLinear(double *c, double *s);
int SolveQuadratic(double *c, double *s);
int SolveCubic  (       double c[4],
                                        double s[3]      
                                );

#endif /* MATRIX3X3_H */ 
