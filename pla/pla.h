
#ifdef __cplusplus
extern "C" {
#endif //__cplusplus

#ifndef PORTABLE_LINEAR_ALGEBRA_H
#define PORTABLE_LINEAR_ALGEBRA_H


#include <math.h>


#define LINMATH_H_DEFINE_VEC(n) \
\
typedef float vec##n[n]; \
typedef const float const_vec##n[n]; \
\
static inline void vec##n##_zero(vec##n a) \
{ \
    int i; \
    for(i=0; i<(n); ++i) \
        a[i] = 0; \
} \
\
static inline void vec##n##_add_n(vec##n r, vec##n a, vec##n b) \
{ \
    int i; \
    for(i=0; i<(n); ++i) \
        r[i] = a[i] + b[i]; \
} \
\
static inline void vec##n##_add_nc(vec##n r, const_vec##n a, const_vec##n b) \
{ \
    int i; \
    for(i=0; i<(n); ++i) \
        r[i] = a[i] + b[i]; \
} \
\
static inline float* vec##n##_add_rn(vec##n r, vec##n a, vec##n b) \
{ \
    int i; \
    for(i=0; i<(n); ++i) \
        r[i] = a[i] + b[i]; \
    return r; \
} \
\
static inline float* vec##n##_add_rnc(vec##n r, const_vec##n a, const_vec##n b) \
{ \
    int i; \
    for(i=0; i<(n); ++i) \
        r[i] = a[i] + b[i]; \
    return r; \
} \
\
static inline float* vec##n##_add_r(vec##n a, vec##n b) \
{ \
    int i; \
    for(i=0; i<(n); ++i) \
        a[i] += b[i]; \
    return a; \
} \
\
static inline float* vec##n##_add_rc(vec##n a, const_vec##n b) \
{ \
    int i; \
    for(i=0; i<(n); ++i) \
        a[i] += b[i]; \
    return a; \
} \
\
static inline void vec##n##_add(vec##n a, vec##n b) \
{ \
    int i; \
    for(i=0; i<(n); ++i) \
        a[i] += b[i]; \
} \
\
static inline void vec##n##_add_c(vec##n a, const_vec##n b) \
{ \
    int i; \
    for(i=0; i<(n); ++i) \
        a[i] += b[i]; \
} \
\
static inline void vec##n##_sub_n(vec##n r, vec##n a, vec##n b) \
{ \
    int i; \
    for(i=0; i<(n); ++i) \
        r[i] = a[i] - b[i]; \
} \
\
static inline void vec##n##_sub_nc(vec##n r, const_vec##n a, const_vec##n b) \
{ \
    int i; \
    for(i=0; i<(n); ++i) \
        r[i] = a[i] - b[i]; \
} \
\
static inline float* vec##n##_sub_r(vec##n a, vec##n b) \
{ \
    int i; \
    for(i=0; i<(n); ++i) \
        a[i] -= b[i]; \
    return a; \
} \
\
static inline float* vec##n##_sub_rc(vec##n a, const_vec##n b) \
{ \
    int i; \
    for(i=0; i<(n); ++i) \
        a[i] -= b[i]; \
    return a; \
} \
\
static inline void vec##n##_sub(vec##n a, vec##n b) \
{ \
    int i; \
    for(i=0; i<(n); ++i) \
        a[i] -= b[i]; \
} \
\
static inline void vec##n##_sub_c(vec##n a, const_vec##n b) \
{ \
    int i; \
    for(i=0; i<(n); ++i) \
        a[i] -= b[i]; \
} \
\
static inline void vec##n##_scale_n(vec##n r, vec##n v, float s) \
{ \
    int i; \
    for(i=0; i<(n); ++i) \
        r[i] = v[i] * s; \
} \
\
static inline void vec##n##_scale_nc(vec##n r, const_vec##n v, float s) \
{ \
    int i; \
    for(i=0; i<(n); ++i) \
        r[i] = v[i] * s; \
} \
\
static inline float* vec##n##_scale_rn(vec##n r, vec##n v, float s) \
{ \
    int i; \
    for(i=0; i<(n); ++i) \
        r[i] = v[i] * s; \
    return r; \
} \
\
static inline float* vec##n##_scale_rnc(vec##n r, const_vec##n v, float s) \
{ \
    int i; \
    for(i=0; i<(n); ++i) \
        r[i] = v[i] * s; \
    return r; \
} \
\
static inline float* vec##n##_scale_r(vec##n v, float s) \
{ \
    int i; \
    for(i=0; i<(n); ++i) \
        v[i] *= s; \
    return v; \
} \
\
static inline void vec##n##_scale(vec##n v, float s) \
{ \
    int i; \
    for(i=0; i<(n); ++i) \
        v[i] *= s; \
} \
\
static inline float vec##n##_mul_inner(vec##n a, vec##n b) \
{ \
    float p = 0.0f; \
    int i; \
    for(i=0; i<(n); ++i) \
        p += b[i]*a[i]; \
    return p; \
} \
\
static inline float vec##n##_mul_inner_c(const_vec##n a, const_vec##n b) \
{ \
    float p = 0.0f; \
    int i; \
    for(i=0; i<(n); ++i) \
        p += b[i]*a[i]; \
    return p; \
} \
\
static inline float vec##n##_norm(vec##n v) \
{ \
    return sqrtf(vec##n##_mul_inner(v,v)); \
} \
\
static inline float vec##n##_norm_c(const_vec##n v) \
{ \
    return sqrtf(vec##n##_mul_inner_c(v,v)); \
} \
\
static inline void vec##n##_normalize_n(vec##n r, vec##n v) \
{ \
    float k = 1.0f / vec##n##_norm(v); \
    vec##n##_scale_n(r, v, k); \
} \
\
static inline void vec##n##_normalize_nc(vec##n r, const_vec##n v) \
{ \
    float k = 1.0f / vec##n##_norm_c(v); \
    vec##n##_scale_nc(r, v, k); \
} \
\
static inline void vec##n##_normalize(vec##n v) \
{ \
    float k = 1.0f / vec##n##_norm(v); \
    vec##n##_scale(v, k); \
} \
\
static inline float* vec##n##_normalize_r(vec##n v) \
{ \
    float k = 1.0f / vec##n##_norm(v); \
    return vec##n##_scale_r(v, k); \
} \
\
static inline void vec##n##_min(vec##n r, vec##n a, vec##n b) \
{ \
    int i; \
    for(i=0; i<(n); ++i) \
        r[i] = a[i]<b[i] ? a[i] : b[i]; \
} \
\
static inline void vec##n##_min_c(vec##n r, const_vec##n a, const_vec##n b) \
{ \
    int i; \
    for(i=0; i<(n); ++i) \
        r[i] = a[i]<b[i] ? a[i] : b[i]; \
} \
\
static inline void vec##n##_max(vec##n r, vec##n a, vec##n b) \
{ \
    int i; \
    for(i=0; i<(n); ++i) \
        r[i] = a[i]>b[i] ? a[i] : b[i]; \
} \
\
static inline void vec##n##_max_c(vec##n r, const_vec##n a, const_vec##n b) \
{ \
    int i; \
    for(i=0; i<(n); ++i) \
        r[i] = a[i]>b[i] ? a[i] : b[i]; \
} \
\
static inline void vec##n##_copy(vec##n a, vec##n b) \
{ \
    int i; \
    for(i=0; i<(n); ++i) \
        a[i] = b[i]; \
} \
\
static inline void vec##n##_copy_c(vec##n a, const_vec##n b) \
{ \
    int i; \
    for(i=0; i<(n); ++i) \
        a[i] = b[i]; \
}

LINMATH_H_DEFINE_VEC(2)
LINMATH_H_DEFINE_VEC(3)
LINMATH_H_DEFINE_VEC(4)
LINMATH_H_DEFINE_VEC(5)
LINMATH_H_DEFINE_VEC(6)


#define ZERO2  {0.0f,0.0f}
#define ZERO3  {0.0f,0.0f,0.0f}
#define ZERO4  {0.0f,0.0f,0.0f,0.0f}
#define ZERO5  {0.0f,0.0f,0.0f,0.0f,0.0f}
#define ZERO6  {0.0f,0.0f,0.0f,0.0f,0.0f,0.0f}


static inline void vec3_mul_cross_n(vec3 r, vec3 a, vec3 b)
{
    r[0] = a[1] * b[2] - a[2] * b[1];
    r[1] = a[2] * b[0] - a[0] * b[2];
    r[2] = a[0] * b[1] - a[1] * b[0];
}

static inline void vec3_mul_cross_nc(vec3 r, const_vec3 a, const_vec3 b)
{
    r[0] = a[1] * b[2] - a[2] * b[1];
    r[1] = a[2] * b[0] - a[0] * b[2];
    r[2] = a[0] * b[1] - a[1] * b[0];
}

static inline float * vec3_mul_cross_rn(vec3 r, vec3 a, vec3 b)
{
    r[0] = a[1] * b[2] - a[2] * b[1];
    r[1] = a[2] * b[0] - a[0] * b[2];
    r[2] = a[0] * b[1] - a[1] * b[0];
    return r;
}

static inline float * vec3_mul_cross_rnc(vec3 r, const_vec3 a, const_vec3 b)
{
    r[0] = a[1] * b[2] - a[2] * b[1];
    r[1] = a[2] * b[0] - a[0] * b[2];
    r[2] = a[0] * b[1] - a[1] * b[0];
    return r;
}

static inline float * vec3_mul_cross_r(vec3 a, vec3 b)
{
    vec3 result = ZERO3;
    result[0] = a[1] * b[2] - a[2] * b[1];
    result[1] = a[2] * b[0] - a[0] * b[2];
    result[2] = a[0] * b[1] - a[1] * b[0];
    vec3_copy(a, result);
    return a;
}

static inline float * vec3_mul_cross_rc(vec3 a, const_vec3 b)
{
    vec3 result = ZERO3;
    result[0] = a[1] * b[2] - a[2] * b[1];
    result[1] = a[2] * b[0] - a[0] * b[2];
    result[2] = a[0] * b[1] - a[1] * b[0];
    vec3_copy(a, result);
    return a;
}

static inline void vec3_mul_cross(vec3 a, vec3 b)
{
    vec3 result = ZERO3;
    result[0] = a[1] * b[2] - a[2] * b[1];
    result[1] = a[2] * b[0] - a[0] * b[2];
    result[2] = a[0] * b[1] - a[1] * b[0];
    vec3_copy(a, result);
}

static inline void vec3_mul_cross_c(vec3 a, const_vec3 b)
{
    vec3 result = ZERO3;
    result[0] = a[1] * b[2] - a[2] * b[1];
    result[1] = a[2] * b[0] - a[0] * b[2];
    result[2] = a[0] * b[1] - a[1] * b[0];
    vec3_copy(a, result);
}

static inline void vec3_reflect_n(vec3 r, vec3 v, vec3 n)
{
    float p = 2.0f * vec3_mul_inner(v, n);
    int i;
    for (i = 0; i < 3; ++ i) {
        r[i] = v[i] - p * n[i];
    }
}

static inline float * vec3_reflect_r(vec3 v, vec3 n) {
    float p = 2.f * vec3_mul_inner(v, n);
    vec3 result = {0, 0, 0};
    int i;
    for (i = 0; i < 3; ++i) {
        result[i] = v[i] - p * n[i];
    }
    vec3_copy(v, result);
    return v;
}

static inline void vec3_reflect(vec3 v, vec3 n) {
    float p = 2.f * vec3_mul_inner(v, n);
    vec3 result = {0, 0, 0};
    int i;
    for (i = 0; i < 3; ++i) {
        result[i] = v[i] - p * n[i];
    }
    vec3_copy(v, result);
}

static inline void vec4_mul_cross_n(vec4 result, vec4 a, vec4 b) {
    result[0] = a[1] * b[2] - a[2] * b[1];
    result[1] = a[2] * b[0] - a[0] * b[2];
    result[2] = a[0] * b[1] - a[1] * b[0];
    result[3] = 1.0f;
}

static inline void vec4_mul_cross(vec4 a, vec4 b) {
    vec4 result = {0, 0, 0, 0};
    result[0] = a[1] * b[2] - a[2] * b[1];
    result[1] = a[2] * b[0] - a[0] * b[2];
    result[2] = a[0] * b[1] - a[1] * b[0];
    result[3] = 1.0f;
    vec4_copy(a, result);
}

static inline float * vec4_mul_cross_r(vec4 a, vec4 b) {
    vec4 result = {0, 0, 0, 0};
    result[0] = a[1] * b[2] - a[2] * b[1];
    result[1] = a[2] * b[0] - a[0] * b[2];
    result[2] = a[0] * b[1] - a[1] * b[0];
    result[3] = 1.0f;
    vec4_copy(a, result);
    return a;
}

static inline void vec4_reflect_n(vec4 result, vec4 v, vec4 n) {
    float p = 2.f * vec4_mul_inner(v, n);
    int i;
    for (i = 0; i < 4; ++ i) {
        result[i] = v[i] - p * n[i];
    }
}

static inline void vec4_reflect(vec4 v, vec4 n) {
    vec4 result = {0, 0, 0, 0};
    float p = 2.f * vec4_mul_inner(v, n);
    int i;
    for (i = 0; i < 4; ++ i) {
        result[i] = v[i] - p * n[i];
    }
    vec4_copy(v, result);
}

static inline float * vec4_reflect_r(vec4 v, vec4 n) {
    vec4 result = {0, 0, 0, 0};
    float p = 2.f * vec4_mul_inner(v, n);
    int i;
    for (i = 0; i < 4; ++ i) {
        result[i] = v[i] - p * n[i];
    }
    vec4_copy(v, result);
    return v;
}


#define LINMATH_H_DEFINE_SQUARE_MAT(n) \
\
typedef vec##n mat##n##x##n[n]; \
typedef const_vec##n const_mat##n##x##n[n]; \
\
static inline void mat##n##x##n##_copy(mat##n##x##n A, mat##n##x##n B) { \
    int row, col; \
    for (row=0; row<(n); ++row) { \
        for (col=0; col<(n); ++col) { \
            A[row][col] = B[row][col]; \
        } \
    } \
} \
static inline void mat##n##x##n##_identity(mat##n##x##n A) { \
    int row, col; \
    for (row=0; row<(n); ++row) { \
        for (col=0; col<(n); ++col) { \
            A[row][col] = (row == col) ? 1.0f : 0.0f; \
        } \
    } \
} \
static inline vec##n * mat##n##x##n##_identity_r(mat##n##x##n A) { \
    int row, col; \
    for (row=0; row<(n); ++row) { \
        for (col=0; col<(n); ++col) { \
            A[row][col] = (row == col) ? 1.0f : 0.0f; \
        } \
    } \
    return A; \
} \
static inline void mat##n##x##n##_row(vec##n row, mat##n##x##n A, unsigned int row_index) { \
    unsigned int col; \
    for (col=0; col<(n); ++col) { \
        row[col] = A[row_index][col]; \
    } \
} \
static inline float * mat##n##x##n##_row_r(vec##n row, mat##n##x##n A, unsigned int row_index) { \
    unsigned int col; \
    for (col=0; col<(n); ++col) { \
        row[col] = A[row_index][col]; \
    } \
    return row; \
} \
static inline void mat##n##x##n##_col(vec##n col, mat##n##x##n A, unsigned int col_index) { \
    unsigned int row; \
    for (row=0; row<(n); ++row) { \
        col[row] = A[row][col_index]; \
    } \
} \
static inline float * mat##n##x##n##_col_r(vec##n col, mat##n##x##n A, unsigned int col_index) { \
    unsigned int row; \
    for (row=0; row<(n); ++row) { \
        col[row] = A[row][col_index]; \
    } \
    return col; \
} \
static inline void mat##n##x##n##_transpose(mat##n##x##n A) { \
    mat##n##x##n temp; \
    int row, col; \
    for (row=0; row<(n); ++row) { \
        for (col=0; col<(n); ++col) { \
            temp[row][col] = A[col][row]; \
        } \
    } \
    mat##n##x##n##_copy(A,temp); \
} \
static inline vec##n * mat##n##x##n##_transpose_r(mat##n##x##n A) { \
    mat##n##x##n temp; \
    int row, col; \
    for (row=0; row<(n); ++row) { \
        for (col=0; col<(n); ++col) { \
            temp[row][col] = A[col][row]; \
        } \
    } \
    mat##n##x##n##_copy(A,temp); \
    return A; \
} \
static inline void mat##n##x##n##_transpose_n(mat##n##x##n M, mat##n##x##n A) { \
    int row, col; \
    for (row=0; row<(n); ++row) { \
        for (col=0; col<(n); ++col) { \
            M[row][col] = A[col][row]; \
        } \
    } \
} \
static inline void mat##n##x##n##_add_n(mat##n##x##n result, mat##n##x##n A, mat##n##x##n B) \
{ \
    int row, col; \
    for (row=0; row<(n); ++row) { \
        for (col=0; col<(n); ++col) { \
            result[row][col] = A[row][col]+B[row][col]; \
        } \
    } \
}\
static inline void mat##n##x##n##_add(mat##n##x##n A, mat##n##x##n B) \
{ \
    int row, col; \
    for (row=0; row<(n); ++row) { \
        for (col=0; col<(n); ++col) { \
            A[row][col] += B[row][col]; \
        } \
    } \
} \
static inline vec##n * mat##n##x##n##_add_r(mat##n##x##n A, mat##n##x##n B) \
{ \
    int row, col; \
    for (row=0; row<(n); ++row) { \
        for (col=0; col<(n); ++col) { \
            A[row][col] += B[row][col]; \
        } \
    } \
    return A; \
} \
static inline void mat##n##x##n##_sub_n(mat##n##x##n result, mat##n##x##n A, mat##n##x##n B) \
{ \
    int row, col; \
    for (row=0; row<(n); ++row) { \
        for (col=0; col<(n); ++col) { \
            result[row][col] = A[row][col]-B[row][col]; \
        } \
    } \
}\
static inline void mat##n##x##n##_sub(mat##n##x##n A, mat##n##x##n B) \
{ \
    int row, col; \
    for (row=0; row<(n); ++row) { \
        for (col=0; col<(n); ++col) { \
            A[row][col] -= B[row][col]; \
        } \
    } \
} \
static inline vec##n * mat##n##x##n##_sub_r(mat##n##x##n A, mat##n##x##n B) \
{ \
    int row, col; \
    for (row=0; row<(n); ++row) { \
        for (col=0; col<(n); ++col) { \
            A[row][col] -= B[row][col]; \
        } \
    } \
    return A; \
} \
static inline void mat##n##x##n##_scale_n(mat##n##x##n result, mat##n##x##n A, float scale) \
{ \
    int row, col; \
    for (row=0; row<(n); ++row) { \
        for (col=0; col<(n); ++col) { \
            result[row][col] = A[row][col] * scale; \
        } \
    } \
}\
static inline void mat##n##x##n##_scale(mat##n##x##n A, float scale) \
{ \
    int row, col; \
    for (row=0; row<(n); ++row) { \
        for (col=0; col<(n); ++col) { \
            A[row][col] *= scale; \
        } \
    } \
} \
static inline vec##n * mat##n##x##n##_scale_r(mat##n##x##n A, float scale) \
{ \
    int row, col; \
    for (row=0; row<(n); ++row) { \
        for (col=0; col<(n); ++col) { \
            A[row][col] *= scale; \
        } \
    } \
    return A; \
} \
static inline void mat##n##x##n##_mul_vec##n##_n(vec##n result, mat##n##x##n A, vec##n B) { \
    int row, col; \
    for(row=0; row<(n); ++row){ \
        result[row] = 0.0f; \
        for(col=0; col<(n); ++col){ \
            result[row] += A[row][col] * B[col]; \
        } \
    } \
} \
static inline float * mat##n##x##n##_mul_vec##n##_r(mat##n##x##n A, vec##n B) { \
    vec##n result; \
    int row, col; \
    for(row=0; row<(n); ++row){ \
        result[row] = 0.0f; \
        for(col=0; col<(n); ++col){ \
            result[row] += A[row][col] * B[col]; \
        } \
    } \
    vec##n##_copy(B,result); \
    return B; \
} \
static inline float * mat##n##x##n##_mul_vec##n##_rn(vec##n result, mat##n##x##n A, vec##n B) { \
    int row, col; \
    for(row=0; row<(n); ++row){ \
        result[row] = 0.0f; \
        for(col=0; col<(n); ++col){ \
            result[row] += A[row][col] * B[col]; \
        } \
    } \
    return result; \
}\
static inline void mat##n##x##n##_mul(mat##n##x##n a, mat##n##x##n b) { \
    mat##n##x##n temp;\
    unsigned int row, col, keep; \
    for (row=0; row<(n); ++row) { \
        for (col = 0; col<(n); ++col) { \
            temp[row][col] = 0.0f; \
            for (keep = 0; keep<(n); ++keep) { \
                temp[row][col] += a[row][keep] * b[keep][col]; \
            } \
        } \
    } \
    mat##n##x##n##_copy(a, temp); \
} \
static inline void mat##n##x##n##_mul_n(mat##n##x##n result, mat##n##x##n a, mat##n##x##n b) { \
    mat##n##x##n temp;\
    unsigned int row, col, keep; \
    for (row=0; row<(n); ++row) { \
        for (col = 0; col<(n); ++col) { \
            temp[row][col] = 0.0f; \
            for (keep = 0; keep<(n); ++keep) { \
                temp[row][col] += a[row][keep] * b[keep][col]; \
            } \
        } \
    } \
    mat##n##x##n##_copy(result, temp); \
} \
static inline vec##n * mat##n##x##n##_mul_nr(mat##n##x##n result, mat##n##x##n a, mat##n##x##n b) { \
    unsigned int row, col, keep; \
    for (row=0; row<(n); ++row) { \
        for (col = 0; col<(n); ++col) { \
            result[row][col] = 0.0f; \
            for (keep = 0; keep<(n); ++keep) { \
                result[row][col] += a[row][keep] * b[keep][col]; \
            } \
        } \
    } \
    return result; \
} \
static inline vec##n * mat##n##x##n##_mul_r(mat##n##x##n a, mat##n##x##n b) { \
    mat##n##x##n temp;\
    unsigned int row, col, keep; \
    for (row=0; row<(n); ++row) { \
        for (col = 0; col<(n); ++col) { \
            temp[row][col] = 0.0f; \
            for (keep = 0; keep<(n); ++keep) { \
                temp[row][col] += a[row][keep] * b[keep][col]; \
            } \
        } \
    } \
    mat##n##x##n##_copy(a, temp); \
    return a; \
}

LINMATH_H_DEFINE_SQUARE_MAT(2)
LINMATH_H_DEFINE_SQUARE_MAT(3)
LINMATH_H_DEFINE_SQUARE_MAT(4)
LINMATH_H_DEFINE_SQUARE_MAT(5)
LINMATH_H_DEFINE_SQUARE_MAT(6)


#define IDENTITY2x2 {{1,0},{0,1}}
#define IDENTITY3x3 {{1,0,0},{0,1,0},{0,0,1}}
#define IDENTITY4x4 {{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}}
#define IDENTITY5x5 {{1,0,0,0,0},{0,1,0,0,0},{0,0,1,0,0},{0,0,0,1,0},{0,0,0,0,1}}
#define IDENTITY6x6 {{1,0,0,0,0,0},{0,1,0,0,0,0},{0,0,1,0,0,0},{0,0,0,1,0,0},{0,0,0,0,1,0},{0,0,0,0,0,1}}


#define ZERO2x2 {{0,0},{0,0}}
#define ZERO3x3 {{0,0,0},{0,0,0},{0,0,0}}
#define ZERO4x4 {{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}}
#define ZERO5x5 {{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0}}
#define ZERO6x6 {{0,0,0,0,0,0},{0,0,0,0,0,0},{0,0,0,0,0,0},{0,0,0,0,0,0},{0,0,0,0,0,0},{0,0,0,0,0,0}}


static inline void mat4x4_get_rotational(mat3x3 M, mat4x4 A) {
    int row, col;
    for (row=0; row<3; ++row) {
        for(col=0; col<3; ++col) {
            M[row][col] = A[row][col];
        }
    }
}

static inline vec3 * mat4x4_get_rotational_rn(mat3x3 M, mat4x4 A) {
    int row, col;
    for (row=0; row<3; ++row) {
        for(col=0; col<3; ++col) {
            M[row][col] = A[row][col];
        }
    }
    return M;
}

static inline void mat4x4_translate(mat4x4 T, float x, float y, float z) {
    T[0][3] += x;
    T[1][3] += y;
    T[2][3] += z;
}


// This is generically written but its intent is to create a
// rotation matrix from an axis-angle representation of the rotation
static inline void mat4x4_from_vec3_mul_outer(mat4x4 M, vec3 a, vec3 b) {
    int row, col;
    for (row=0; row<4; ++row) {
        for (col=0; col<4; ++col) {
            M[row][col] = (row < 3 && col < 3) ? a[row]*b[col] : 0.0f;
        }
    }
}


// Axis Angle Rotation with Radian
static inline void mat4x4_rotate(mat4x4 R, mat4x4 M,
                                 float x, float y, float z, float angle) {
    float sin_angle     = sinf(angle);
    float cos_angle     = cosf(angle);
    vec3  rotation_axis = {x, y, z};

    if (vec3_norm(rotation_axis) > 1e-4) {
        vec3_normalize(rotation_axis);
        mat4x4 T;
        mat4x4_from_vec3_mul_outer(T, rotation_axis, rotation_axis);

        mat4x4 skew_sym = {
            {                0, -rotation_axis[2],  rotation_axis[1], 0},
            { rotation_axis[2],                 0, -rotation_axis[0], 0},
            {-rotation_axis[1],  rotation_axis[0],                 0, 0},
            {                0,                 0,                 0, 0}
        };
        mat4x4_scale(skew_sym, sin_angle);

        mat4x4 C;
        mat4x4_identity(C);
        mat4x4_sub(C, T);

        mat4x4_scale(C, cos_angle);

        mat4x4_add(T, C);
        mat4x4_add(T, skew_sym);

        T[3][3] = 1.0f;
        mat4x4_mul_n(R, M, T);
    } else {
        mat4x4_copy(R, M);
    }
}


static inline void mat4x4_rotate_about_x(mat4x4 Q, mat4x4 M, float angle) {
    float s = sinf(angle);
    float c = cosf(angle);
    mat4x4 R = {
        { 1.0f, 0.0f, 0.0f, 0.0f},
        { 0.0f,    c,   -s, 0.0f},
        { 0.0f,    s,    c, 0.0f},
        { 0.0f, 0.0f, 0.0f, 1.0f}
    };
    mat4x4_mul_n(Q, M, R);
}


static inline void mat4x4_rotate_about_y(mat4x4 Q, mat4x4 M, float angle) {
    float s = sinf(angle);
    float c = cosf(angle);
    mat4x4 R = {
        {    c, 0.0f,    s, 0.0f},
        { 0.0f, 1.0f, 0.0f, 0.0f},
        {   -s, 0.0f,    c, 0.0f},
        { 0.0f, 0.0f, 0.0f, 1.0f}
    };
    mat4x4_mul_n(Q, M, R);
}


static inline void mat4x4_rotate_about_z(mat4x4 Q, mat4x4 M, float angle) {
    float s = sinf(angle);
    float c = cosf(angle);
    mat4x4 R = {
        {    c,   -s, 0.0f, 0.0f},
        {    s,    c, 0.0f, 0.0f},
        { 0.0f, 0.0f, 1.0f, 0.0f},
        { 0.0f, 0.0f, 0.0f, 1.0f}
    };
    mat4x4_mul_n(Q, M, R);
}


static inline void mat6x6_get_sub_mat3x3(mat3x3 sub, mat6x6 mat, unsigned int start_row,
                                         unsigned int end_row, unsigned int start_col, unsigned int end_col)
{
    unsigned int mat_row, mat_col, sub_row, sub_col;
    for (mat_row=start_row, sub_row=0; mat_row<end_row; ++mat_row, ++sub_row) {
        for (mat_col=start_col, sub_col=0; mat_col<end_col; ++mat_col, ++sub_col) {
            sub[sub_row][sub_col] = mat[mat_row][mat_col];
        }
    }
}


static inline void combine_mat3x3_to_mat6x6(mat6x6 mat, mat3x3 a, mat3x3 b, mat3x3 c, mat3x3 d)
{
    unsigned int row, col;
    for (row=0; row<3; ++row) {
        for (col=0; col<3; ++col) {
            mat[  row][  col] = a[row][col];
            mat[  row][col+3] = b[row][col];
            mat[row+3][  col] = c[row][col];
            mat[row+3][col+3] = d[row][col];
        }
    }
}


static inline void mat3x3_invert(mat3x3 T, mat3x3 M)
{
    float a0 = M[0][0];
    float a1 = M[1][0];
    float a2 = M[2][0];
    float b0 = M[0][1];
    float b1 = M[1][1];
    float b2 = M[2][1];
    float c0 = M[0][2];
    float c1 = M[1][2];
    float c2 = M[2][2];
    float t2 = a0*b1*c2;
    float t3 = a1*b2*c0;
    float t4 = a2*b0*c1;
    float t7 = a0*b2*c1;
    float t8 = a1*b0*c2;
    float t9 = a2*b1*c0;
    float t5 = t2+t3+t4-t7-t8-t9;
    float t6 = 1.0f/t5;

    T[0][0] =  t6*(b1*c2-b2*c1);
    T[1][0] = -t6*(a1*c2-a2*c1);
    T[2][0] =  t6*(a1*b2-a2*b1);
    T[0][1] = -t6*(b0*c2-b2*c0);
    T[1][1] =  t6*(a0*c2-a2*c0);
    T[2][1] = -t6*(a0*b2-a2*b0);
    T[0][2] =  t6*(b0*c1-b1*c0);
    T[1][2] = -t6*(a0*c1-a1*c0);
    T[2][2] =  t6*(a0*b1-a1*b0);
}


static inline void mat4x4_invert(mat4x4 T, mat4x4 M) {
    float s[6];
    float c[6];

    s[0] = M[0][0] * M[1][1] - M[1][0] * M[0][1];
    s[1] = M[0][0] * M[1][2] - M[1][0] * M[0][2];
    s[2] = M[0][0] * M[1][3] - M[1][0] * M[0][3];
    s[3] = M[0][1] * M[1][2] - M[1][1] * M[0][2];
    s[4] = M[0][1] * M[1][3] - M[1][1] * M[0][3];
    s[5] = M[0][2] * M[1][3] - M[1][2] * M[0][3];

    c[0] = M[2][0] * M[3][1] - M[3][0] * M[2][1];
    c[1] = M[2][0] * M[3][2] - M[3][0] * M[2][2];
    c[2] = M[2][0] * M[3][3] - M[3][0] * M[2][3];
    c[3] = M[2][1] * M[3][2] - M[3][1] * M[2][2];
    c[4] = M[2][1] * M[3][3] - M[3][1] * M[2][3];
    c[5] = M[2][2] * M[3][3] - M[3][2] * M[2][3];

    float idet = 1.0f /
        (s[0] * c[5] - s[1] * c[4] + s[2] * c[3] + s[3] * c[2] - s[4] * c[1] + s[5] * c[0]);

    T[0][0] = (  M[1][1] * c[5] - M[1][2] * c[4] + M[1][3] * c[3]) * idet;
    T[0][1] = (- M[0][1] * c[5] + M[0][2] * c[4] - M[0][3] * c[3]) * idet;
    T[0][2] = (  M[3][1] * s[5] - M[3][2] * s[4] + M[3][3] * s[3]) * idet;
    T[0][3] = (- M[2][1] * s[5] + M[2][2] * s[4] - M[2][3] * s[3]) * idet;

    T[1][0] = (- M[1][0] * c[5] + M[1][2] * c[2] - M[1][3] * c[1]) * idet;
    T[1][1] = (  M[0][0] * c[5] - M[0][2] * c[2] + M[0][3] * c[1]) * idet;
    T[1][2] = (- M[3][0] * s[5] + M[3][2] * s[2] - M[3][3] * s[1]) * idet;
    T[1][3] = (  M[2][0] * s[5] - M[2][2] * s[2] + M[2][3] * s[1]) * idet;

    T[2][0] = (  M[1][0] * c[4] - M[1][1] * c[2] + M[1][3] * c[0]) * idet;
    T[2][1] = (- M[0][0] * c[4] + M[0][1] * c[2] - M[0][3] * c[0]) * idet;
    T[2][2] = (  M[3][0] * s[4] - M[3][1] * s[2] + M[3][3] * s[0]) * idet;
    T[2][3] = (- M[2][0] * s[4] + M[2][1] * s[2] - M[2][3] * s[0]) * idet;

    T[3][0] = (- M[1][0] * c[3] + M[1][1] * c[1] - M[1][2] * c[0]) * idet;
    T[3][1] = (  M[0][0] * c[3] - M[0][1] * c[1] + M[0][2] * c[0]) * idet;
    T[3][2] = (- M[3][0] * s[3] + M[3][1] * s[1] - M[3][2] * s[0]) * idet;
    T[3][3] = (  M[2][0] * s[3] - M[2][1] * s[1] + M[2][2] * s[0]) * idet;
}


// Quaternion is defined as {x, y, z, w}
typedef float quat[4];

#define IDENTITY_QUAT {0.0f, 0.0f, 0.0f, 1.0f};

#define quat_copy vec4_copy

static inline void set_quat_identity(quat q) {
    q[0] = q[1] = q[2] = 0.0f;
    q[3] = 1.0f;
}


static inline void quat_add(quat r, quat a, quat b) {
    int i;
    for (i = 0; i < 4; ++ i) {
        r[i] = a[i] + b[i];
    }
}


static inline void quat_sub(quat r, quat a, quat b) {
    int i;
    for (i = 0; i < 4; ++ i) {
        r[i] = a[i] - b[i];
    }
}


static inline void quat_mul(quat r, quat p, quat q) {
    vec3 w;
    vec3_mul_cross_n(r, p, q);
    vec3_scale_n(w, p, q[3]);
    vec3_add_n(r, r, w);
    vec3_scale_n(w, q, p[3]);
    vec3_add_n(r, r, w);
    r[3] = p[3] * q[3] - vec3_mul_inner(p, q);
}


static inline void quat_scale(quat r, quat v, float s) {
    int i;
    for (i = 0; i < 4; ++ i) {
        r[i] = v[i] * s;
    }
}


static inline float quat_inner_product(quat a, quat b) {
    float p = 0.0f;
    int i;
    for (i = 0; i < 4; ++ i) {
        p += a[i] * b[i];
    }
    return p;
}


static inline void quat_conj(quat r, quat q) {
    int i;
    for (i = 0; i < 3; ++ i) {
        r[i] = - q[i];
    }
    r[3] = q[3];
}


static inline void quat_from_angle_axis(quat r, float angle, vec3 axis) {
    vec3 v = ZERO3;
    vec3_scale_n(v, axis, sinf(angle / 2));
    int i;
    for (i = 0; i < 3; ++ i) {
        r[i] = v[i];
    }
    r[3] = cosf(angle / 2);
}


#define quat_norm        vec4_norm
#define quat_normalize_n vec4_normalize_n
#define quat_normalize   vec4_normalize


static inline void quat_mul_vec3(vec3 r, quat q, vec3 v) {
/*
 * Method by Fabian 'ryg' Giessen (of Farbrausch)
t = 2 * cross(q.xyz, v)
v' = v + q.w * t + cross(q.xyz, t)
 */
    vec3 t;
    vec3 q_xyz = {q[0], q[1], q[2]};
    vec3 u = {q[0], q[1], q[2]};

    vec3_mul_cross_n(t, q_xyz, v);
    vec3_scale_n(t, t, 2);

    vec3_mul_cross_n(u, q_xyz, t);
    vec3_scale_n(t, t, q[3]);

    vec3_add_n(r, v, t);
    vec3_add_n(r, r, u);
}


// WARNING! Quaternion must be normalized!
static inline void mat4x4_from_quat(mat4x4 M, quat q) {
    float xx = q[0]*q[0];
    float xy = q[0]*q[1];
    float xz = q[0]*q[2];
    float xw = q[0]*q[3];
    float yy = q[1]*q[1];
    float yz = q[1]*q[2];
    float yw = q[1]*q[3];
    float zz = q[2]*q[2];
    float zw = q[2]*q[3];

    M[0][0] = 1.0f - 2.0f * (yy + zz);
    M[0][1] = 2.0f * (xy - zw);
    M[0][2] = 2.0f * (xz + yw);

    M[1][0] = 2.0f * (xy + zw);
    M[1][1] = 1.0f - 2.0f * (xx + zz);
    M[1][2] = 2.0f * (yz - xw);

    M[2][0] = 2.0f * (xz - yw);
    M[2][1] = 2.0f * (yz + xw);
    M[2][2] = 1.0f - 2.0f * (xx + yy);

    M[0][3] = M[1][3] = M[2][3] = M[3][0] = M[3][1] = M[3][2] = 0.0f;
    M[3][3] = 1.0f;
}


/*
 * Warning! The rotation matrix must be orhtonormalized!
 * Therefore, it must represent a rotational matrix.
 * Calc from http://www.euclideanspace.com/
 */
static inline void quat_from_mat4x4(quat q, mat4x4 M) {
    float trace = M[0][0] + M[1][1] + M[2][2];

    if (trace > 0.0f) {
        float trace_sqrt = sqrtf(trace+1.0f) * 2.0f;
        q[3] = 0.25f * trace_sqrt;
        q[0] = (M[2][1] - M[1][2]) / trace_sqrt;
        q[1] = (M[0][2] - M[2][0]) / trace_sqrt;
        q[2] = (M[1][0] - M[0][1]) / trace_sqrt;
    } else if ((M[0][0] > M[1][1]) && (M[0][0] > M[2][2])) {
        float trace_sqrt = sqrtf(1.0f + M[0][0] - M[1][1] - M[2][2]) * 2.0f;
        q[3] = (M[2][1] - M[1][2]) / trace_sqrt;
        q[0] = 0.25f * trace_sqrt;
        q[1] = (M[0][1] + M[1][0]) / trace_sqrt;
        q[2] = (M[0][2] + M[2][0]) / trace_sqrt;
    } else if (M[1][1] > M[2][2]) {
        float trace_sqrt = sqrtf(1.0f + M[1][1] - M[0][0] - M[2][2]) * 2.0f;
        q[3] = (M[0][2] - M[2][0]) / trace_sqrt;
        q[0] = (M[0][1] + M[1][0]) / trace_sqrt;
        q[1] = 0.25f * trace_sqrt;
        q[2] = (M[2][1] + M[1][2]) / trace_sqrt;
    } else {
        float trace_sqrt = sqrtf(1.0f + M[2][2] - M[0][0] - M[1][1]) * 2.0f;
        q[3] = (M[1][0] - M[0][1]) / trace_sqrt;
        q[0] = (M[0][2] + M[2][0]) / trace_sqrt;
        q[1] = (M[1][2] + M[2][1]) / trace_sqrt;
        q[2] = 0.25f * trace_sqrt;
    }
}


static inline void mat6x6_invert(mat6x6 T, mat6x6 M)
{
    // 6 square matrix inversion with blockwise inversion
    // TODO: no safety check yet for singularities
    mat3x3 a, b, c, d, a_inv, b_inv, c_inv, d_inv, a_final, b_final, c_final, d_final;
    mat6x6_get_sub_mat3x3(a, M, 0, 3, 0, 3);
    mat6x6_get_sub_mat3x3(b, M, 0, 3, 3, 6);
    mat6x6_get_sub_mat3x3(c, M, 3, 6, 0, 3);
    mat6x6_get_sub_mat3x3(d, M, 3, 6, 3, 6);
    mat3x3_invert(a_inv, a);
    mat3x3_invert(b_inv, b);
    mat3x3_invert(c_inv, c);
    mat3x3_invert(d_inv, d);

    mat3x3 inter1 = IDENTITY3x3;
    mat3x3 inter2 = IDENTITY3x3;
    mat3x3 inter3 = IDENTITY3x3;

    mat3x3 dcab = ZERO3x3;
    mat3x3_mul_n(inter1, a_inv, b);
    mat3x3_mul_n(inter2, c, inter1);
    mat3x3_sub_n(inter3, d, inter2);
    mat3x3_invert(dcab, inter3);

    // compute new A
    mat3x3_mul_n(inter1, c, a_inv);

    mat3x3_mul_n(inter2, dcab, inter1);
    mat3x3_mul_n(inter3, b, inter2);
    mat3x3_mul_n(inter2, a_inv, inter3);

    mat3x3_add_n(a_final, a_inv, inter2);

    // compute new B
    mat3x3_mul_n(inter2, b, dcab);
    mat3x3_mul_n(inter1, a_inv, inter2);
    mat3x3_scale_n(b_final, inter1, -1.0f);

    // compute new C
    mat3x3_mul_n(inter1, c, a_inv);

    mat3x3_scale_n(inter2, dcab, -1.0f);

    mat3x3_mul_n(c_final, inter2, inter1);

    // compute new D
    mat3x3_copy(d_final, dcab);

    combine_mat3x3_to_mat6x6(T, a_final, b_final, c_final, d_final);
}


#endif //PORTABLE_LINEAR_ALGEBRA_H

#ifdef __cplusplus
}
#endif //__cplusplus
