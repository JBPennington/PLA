
#ifdef __cplusplus
extern "C" {
#endif //__cplusplus

#ifndef TORQUECALCULATOR_LINEAR_ALGEBRA_H
#define TORQUECALCULATOR_LINEAR_ALGEBRA_H


#include <math.h>


#define LINMATH_H_DEFINE_VEC(n) \
typedef float vec##n[n]; \
static inline void vec##n##_zero(vec##n a) \
{ \
    int i; \
    for(i=0; i<(n); ++i) \
        a[i] = 0; \
} \
static inline void vec##n##_add_n(vec##n r, vec##n const a, vec##n const b) \
{ \
    int i; \
    for(i=0; i<(n); ++i) \
        r[i] = a[i] + b[i]; \
} \
static inline float* vec##n##_add_rn(vec##n r, vec##n const a, vec##n const b) \
{ \
    int i; \
    for(i=0; i<(n); ++i) \
        r[i] = a[i] + b[i]; \
    return r; \
} \
static inline float* vec##n##_add_r(vec##n a, vec##n const b) \
{ \
    int i; \
    for(i=0; i<(n); ++i) \
        a[i] += b[i]; \
    return a; \
} \
static inline void vec##n##_add(vec##n a, vec##n const b) \
{ \
    int i; \
    for(i=0; i<(n); ++i) \
        a[i] += b[i]; \
} \
static inline void vec##n##_sub_n(vec##n r, vec##n const a, vec##n const b) \
{ \
    int i; \
    for(i=0; i<(n); ++i) \
        r[i] = a[i] - b[i]; \
} \
static inline float* vec##n##_sub_r(vec##n a, vec##n const b) \
{ \
    int i; \
    for(i=0; i<(n); ++i) \
        a[i] -= b[i]; \
    return a; \
} \
static inline void vec##n##_sub(vec##n a, vec##n const b) \
{ \
    int i; \
    for(i=0; i<(n); ++i) \
        a[i] -= b[i]; \
} \
static inline void vec##n##_scale_n(vec##n r, vec##n const v, float const s) \
{ \
    int i; \
    for(i=0; i<(n); ++i) \
        r[i] = v[i] * s; \
} \
static inline float* vec##n##_scale_rn(vec##n r, vec##n const v, float const s) \
{ \
    int i; \
    for(i=0; i<(n); ++i) \
        r[i] = v[i] * s; \
    return r; \
} \
static inline float* vec##n##_scale_r(vec##n v, float const s) \
{ \
    int i; \
    for(i=0; i<(n); ++i) \
        v[i] *= s; \
    return v; \
} \
static inline void vec##n##_scale(vec##n v, float const s) \
{ \
    int i; \
    for(i=0; i<(n); ++i) \
        v[i] *= s; \
} \
static inline float vec##n##_mul_inner(vec##n const a, vec##n const b) \
{ \
    float p = 0.0f; \
    int i; \
    for(i=0; i<(n); ++i) \
        p += b[i]*a[i]; \
    return p; \
} \
static inline float vec##n##_norm(vec##n const v) \
{ \
    return sqrtf(vec##n##_mul_inner(v,v)); \
} \
static inline void vec##n##_normalize_n(vec##n r, vec##n const v) \
{ \
    float k = 1.0 / vec##n##_norm(v); \
    vec##n##_scale_n(r, v, k); \
} \
static inline void vec##n##_normalize(vec##n v) \
{ \
    float k = 1.0 / vec##n##_norm(v); \
    vec##n##_scale(v, k); \
} \
static inline float* vec##n##_normalize_r(vec##n v) \
{ \
    float k = 1.0 / vec##n##_norm(v); \
    return vec##n##_scale_r(v, k); \
} \
static inline void vec##n##_min(vec##n r, vec##n a, vec##n b) \
{ \
    int i; \
    for(i=0; i<(n); ++i) \
        r[i] = a[i]<b[i] ? a[i] : b[i]; \
} \
static inline void vec##n##_max(vec##n r, vec##n a, vec##n b) \
{ \
    int i; \
    for(i=0; i<(n); ++i) \
        r[i] = a[i]>b[i] ? a[i] : b[i]; \
} \
static inline void vec##n##_copy(vec##n a, vec##n b) \
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

#define EMPTYVEC3 {0,0,0}


static inline void vec3_mul_cross_n(vec3 r, vec3 const a, vec3 const b) {
    r[0] = a[1] * b[2] - a[2] * b[1];
    r[1] = a[2] * b[0] - a[0] * b[2];
    r[2] = a[0] * b[1] - a[1] * b[0];
}

static inline float * vec3_mul_cross_rn(vec3 r, vec3 const a, vec3 const b) {
    r[0] = a[1] * b[2] - a[2] * b[1];
    r[1] = a[2] * b[0] - a[0] * b[2];
    r[2] = a[0] * b[1] - a[1] * b[0];
    return r;
}

static inline float * vec3_mul_cross_r(vec3 a, vec3 const b) {
    vec3 result = {0,0,0};
    result[0] = a[1] * b[2] - a[2] * b[1];
    result[1] = a[2] * b[0] - a[0] * b[2];
    result[2] = a[0] * b[1] - a[1] * b[0];
    vec3_copy(a, result);
    return a;
}

static inline void vec3_mul_cross(vec3 a, vec3 const b) {
    vec3 result = {0,0,0};
    result[0] = a[1] * b[2] - a[2] * b[1];
    result[1] = a[2] * b[0] - a[0] * b[2];
    result[2] = a[0] * b[1] - a[1] * b[0];
    vec3_copy(a, result);
}

static inline void vec3_reflect_n(vec3 r, vec3 const v, vec3 const n) {
    float p = 2.f * vec3_mul_inner(v, n);
    int i;
    for (i = 0; i < 3; ++ i) {
        r[i] = v[i] - p * n[i];
    }
}

static inline float * vec3_reflect_r(vec3 v, vec3 const n) {
    float p = 2.f * vec3_mul_inner(v, n);
    vec3 result = {0, 0, 0};
    int i;
    for (i = 0; i < 3; ++i) {
        result[i] = v[i] - p * n[i];
    }
    vec3_copy(v, result);
    return v;
}

static inline void vec3_reflect(vec3 v, vec3 const n) {
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
typedef vec##n mat##n##x##n[n]; \
static inline void mat##n##x##n##_copy(mat##n##x##n A, mat##n##x##n const B) { \
    int row, col; \
    for (row=0; row<n; ++row) { \
        for (col=0; col<n; ++col) { \
            A[row][col] = B[row][col]; \
        } \
    } \
} \
static inline void mat##n##x##n##_identity(mat##n##x##n A) { \
    int row, col; \
    for (row=0; row<n; ++row) { \
        for (col=0; col<n; ++col) { \
            A[row][col] = (row == col) ? 1.0f : 0.0f; \
        } \
    } \
} \
static inline vec##n * mat##n##x##n##_identity_r(mat##n##x##n A) { \
    int row, col; \
    for (row=0; row<n; ++row) { \
        for (col=0; col<n; ++col) { \
            A[row][col] = (row == col) ? 1.0f : 0.0f; \
        } \
    } \
    return A; \
} \
static inline void mat##n##x##n##_row(vec##n row, mat##n##x##n const A, unsigned int row_index) { \
    unsigned int col; \
    for (col=0; col<n; ++col) { \
        row[col] = A[row_index][col]; \
    } \
} \
static inline float * mat##n##x##n##_row_r(vec##n row, mat##n##x##n const A, unsigned int row_index) { \
    unsigned int col; \
    for (col=0; col<n; ++col) { \
        row[col] = A[row_index][col]; \
    } \
    return row; \
} \
static inline void mat##n##x##n##_col(vec##n col, mat##n##x##n const A, unsigned int col_index) { \
    unsigned int row; \
    for (row=0; row<n; ++row) { \
        col[row] = A[row][col_index]; \
    } \
} \
static inline float * mat##n##x##n##_col_r(vec##n col, mat##n##x##n const A, unsigned int col_index) { \
    unsigned int row; \
    for (row=0; row<n; ++row) { \
        col[row] = A[row][col_index]; \
    } \
    return col; \
} \
static inline void mat##n##x##n##_transpose(mat##n##x##n A) { \
    mat##n##x##n temp; \
    int row, col; \
    for (row=0; row<n; ++row) { \
        for (col=0; col<n; ++col) { \
            temp[row][col] = A[col][row]; \
        } \
    } \
    mat##n##x##n##_copy(A,temp); \
} \
static inline vec##n * mat##n##x##n##_transpose_r(mat##n##x##n A) { \
    mat##n##x##n temp; \
    int row, col; \
    for (row=0; row<n; ++row) { \
        for (col=0; col<n; ++col) { \
            temp[row][col] = A[col][row]; \
        } \
    } \
    mat##n##x##n##_copy(A,temp); \
    return A; \
} \
static inline void mat##n##x##n##_transpose_n(mat##n##x##n M, mat##n##x##n A) { \
    int row, col; \
    for (row=0; row<n; ++row) { \
        for (col=0; col<n; ++col) { \
            M[row][col] = A[col][row]; \
        } \
    } \
} \
static inline void mat##n##x##n##_add_n(mat##n##x##n result, mat##n##x##n const A, mat##n##x##n const B) \
{ \
    int row, col; \
    for (row=0; row<n; ++row) { \
        for (col=0; col<n; ++col) { \
            result[row][col] = A[row][col]+B[row][col]; \
        } \
    } \
}\
static inline void mat##n##x##n##_add(mat##n##x##n A, mat##n##x##n const B) \
{ \
    int row, col; \
    for (row=0; row<n; ++row) { \
        for (col=0; col<n; ++col) { \
            A[row][col] += B[row][col]; \
        } \
    } \
} \
static inline vec##n * mat##n##x##n##_add_r(mat##n##x##n A, mat##n##x##n const B) \
{ \
    int row, col; \
    for (row=0; row<n; ++row) { \
        for (col=0; col<n; ++col) { \
            A[row][col] += B[row][col]; \
        } \
    } \
    return A; \
} \
static inline void mat##n##x##n##_sub_n(mat##n##x##n result, mat##n##x##n const A, mat##n##x##n const B) \
{ \
    int row, col; \
    for (row=0; row<n; ++row) { \
        for (col=0; col<n; ++col) { \
            result[row][col] = A[row][col]-B[row][col]; \
        } \
    } \
}\
static inline void mat##n##x##n##_sub(mat##n##x##n A, mat##n##x##n const B) \
{ \
    int row, col; \
    for (row=0; row<n; ++row) { \
        for (col=0; col<n; ++col) { \
            A[row][col] -= B[row][col]; \
        } \
    } \
} \
static inline vec##n * mat##n##x##n##_sub_r(mat##n##x##n A, mat##n##x##n const B) \
{ \
    int row, col; \
    for (row=0; row<n; ++row) { \
        for (col=0; col<n; ++col) { \
            A[row][col] -= B[row][col]; \
        } \
    } \
    return A; \
} \
static inline void mat##n##x##n##_scale_n(mat##n##x##n result, mat##n##x##n const A, const float scale) \
{ \
    int row, col; \
    for (row=0; row<n; ++row) { \
        for (col=0; col<n; ++col) { \
            result[row][col] = A[row][col] * scale; \
        } \
    } \
}\
static inline void mat##n##x##n##_scale(mat##n##x##n A, const float scale) \
{ \
    int row, col; \
    for (row=0; row<n; ++row) { \
        for (col=0; col<n; ++col) { \
            A[row][col] *= scale; \
        } \
    } \
} \
static inline vec##n * mat##n##x##n##_scale_r(mat##n##x##n A, const float scale) \
{ \
    int row, col; \
    for (row=0; row<n; ++row) { \
        for (col=0; col<n; ++col) { \
            A[row][col] *= scale; \
        } \
    } \
    return A; \
} \
static inline void mat##n##x##n##_mul_vec##n##_n(vec##n result, mat##n##x##n const A, vec##n const B) { \
    int row, col; \
    for(row=0; row<n; ++row){ \
        result[row] = 0.0f; \
        for(col=0; col<n; ++col){ \
            result[row] += A[row][col] * B[col]; \
        } \
    } \
} \
static inline float * mat##n##x##n##_mul_vec##n##_r(mat##n##x##n const A, vec##n B) { \
    vec##n result; \
    int row, col; \
    for(row=0; row<n; ++row){ \
        result[row] = 0.0f; \
        for(col=0; col<n; ++col){ \
            result[row] += A[row][col] * B[col]; \
        } \
    } \
    vec##n##_copy(B,result); \
    return B; \
} \
static inline float * mat##n##x##n##_mul_vec##n##_rn(vec##n result, mat##n##x##n const A, vec##n const B) { \
    int row, col; \
    for(row=0; row<n; ++row){ \
        result[row] = 0.0f; \
        for(col=0; col<n; ++col){ \
            result[row] += A[row][col] * B[col]; \
        } \
    } \
    return result; \
}\
static inline void mat##n##x##n##_mul(mat##n##x##n a, mat##n##x##n const b) { \
    mat##n##x##n temp;\
    unsigned int row, col, keep; \
    for (row=0; row<n; ++row) { \
        for (col = 0; col<n; ++col) { \
            temp[row][col] = 0.0f; \
            for (keep = 0; keep<n; ++keep) { \
                temp[row][col] += a[row][keep] * b[keep][col]; \
            } \
        } \
    } \
    mat##n##x##n##_copy(a, temp); \
} \
static inline void mat##n##x##n##_mul_n(mat##n##x##n result, mat##n##x##n const a, mat##n##x##n const b) { \
    mat##n##x##n temp;\
    unsigned int row, col, keep; \
    for (row=0; row<n; ++row) { \
        for (col = 0; col<n; ++col) { \
            temp[row][col] = 0.0f; \
            for (keep = 0; keep<n; ++keep) { \
                temp[row][col] += a[row][keep] * b[keep][col]; \
            } \
        } \
    } \
    mat##n##x##n##_copy(result, temp); \
} \
static inline vec##n * mat##n##x##n##_mul_nr(mat##n##x##n result, mat##n##x##n const a, mat##n##x##n const b) { \
    unsigned int row, col, keep; \
    for (row=0; row<n; ++row) { \
        for (col = 0; col<n; ++col) { \
            result[row][col] = 0.0f; \
            for (keep = 0; keep<n; ++keep) { \
                result[row][col] += a[row][keep] * b[keep][col]; \
            } \
        } \
    } \
    return result; \
} \
static inline vec##n * mat##n##x##n##_mul_r(mat##n##x##n a, mat##n##x##n const b) { \
    mat##n##x##n temp;\
    unsigned int row, col, keep; \
    for (row=0; row<n; ++row) { \
        for (col = 0; col<n; ++col) { \
            temp[row][col] = 0.0f; \
            for (keep = 0; keep<n; ++keep) { \
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


#define IDENTITY3x3 {{1,0,0},{0,1,0},{0,0,1}}
#define IDENTITY4x4 {{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}}


static inline void mat4x4_get_rotational(mat3x3 M, mat4x4 const A) {
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


//static inline void mat4x4_from_vec3_mul_outer(mat4x4 M, vec3 a, vec3 b) {
//    int i, j;
//    for (i = 0; i < 4; ++ i)
//        for (j = 0; j < 4; ++ j)
//            M[i][j] = i < 3 && j < 3 ? a[i] * b[j] : 0.f;
//}
//
//static inline void mat4x4_rotate(mat4x4 R, mat4x4 M, float x, float y, float z, float angle) {
//    float s = sinf(angle);
//    float c = cosf(angle);
//    vec3 u = {x, y, z};
//
//    if (vec3_norm(u) > 1e-4) {
//        vec3_normalize(u);
//        mat4x4 T;
//        mat4x4_from_vec3_mul_outer(T, u, u);
//
//        mat4x4 S = {
//            {0, u[2], - u[1], 0},
//            {- u[2], 0, u[0], 0},
//            {u[1], - u[0], 0, 0},
//            {0, 0, 0, 0}
//        };
//        mat4x4_scale(S, S, s);
//
//        mat4x4 C;
//        mat4x4_identity(C);
//        mat4x4_sub(C, C, T);
//
//        mat4x4_scale(C, C, c);
//
//        mat4x4_add(T, T, C);
//        mat4x4_add(T, T, S);
//
//        T[3][3] = 1.;
//        mat4x4_mul(R, M, T);
//    } else {
//        mat4x4_copy(R, M);
//    }
//}
//
//static inline void mat4x4_rotate_X(mat4x4 Q, mat4x4 M, float angle) {
//    float s = sinf(angle);
//    float c = cosf(angle);
//    mat4x4 R = {
//        {1.f, 0.f, 0.f, 0.f},
//        {0.f, c, s, 0.f},
//        {0.f, - s, c, 0.f},
//        {0.f, 0.f, 0.f, 1.f}
//    };
//    mat4x4_mul(Q, M, R);
//}
//
//static inline void mat4x4_rotate_Y(mat4x4 Q, mat4x4 M, float angle) {
//    float s = sinf(angle);
//    float c = cosf(angle);
//    mat4x4 R = {
//        {c, 0.f, s, 0.f},
//        {0.f, 1.f, 0.f, 0.f},
//        {- s, 0.f, c, 0.f},
//        {0.f, 0.f, 0.f, 1.f}
//    };
//    mat4x4_mul(Q, M, R);
//}
//
//static inline void mat4x4_rotate_Z(mat4x4 Q, mat4x4 M, float angle) {
//    float s = sinf(angle);
//    float c = cosf(angle);
//    mat4x4 R = {
//        {c, s, 0.f, 0.f},
//        {- s, c, 0.f, 0.f},
//        {0.f, 0.f, 1.f, 0.f},
//        {0.f, 0.f, 0.f, 1.f}
//    };
//    mat4x4_mul(Q, M, R);
//}
//
//static inline void mat4x4_invert(mat4x4 T, mat4x4 M) {
//    float s[6];
//    float c[6];
//    s[0] = M[0][0] * M[1][1] - M[1][0] * M[0][1];
//    s[1] = M[0][0] * M[1][2] - M[1][0] * M[0][2];
//    s[2] = M[0][0] * M[1][3] - M[1][0] * M[0][3];
//    s[3] = M[0][1] * M[1][2] - M[1][1] * M[0][2];
//    s[4] = M[0][1] * M[1][3] - M[1][1] * M[0][3];
//    s[5] = M[0][2] * M[1][3] - M[1][2] * M[0][3];
//
//    c[0] = M[2][0] * M[3][1] - M[3][0] * M[2][1];
//    c[1] = M[2][0] * M[3][2] - M[3][0] * M[2][2];
//    c[2] = M[2][0] * M[3][3] - M[3][0] * M[2][3];
//    c[3] = M[2][1] * M[3][2] - M[3][1] * M[2][2];
//    c[4] = M[2][1] * M[3][3] - M[3][1] * M[2][3];
//    c[5] = M[2][2] * M[3][3] - M[3][2] * M[2][3];
//
//    /* Assumes it is invertible */
//    float idet = 1.0f / (s[0] * c[5] - s[1] * c[4] + s[2] * c[3] + s[3] * c[2] - s[4] * c[1] + s[5] * c[0]);
//
//    T[0][0] = (M[1][1] * c[5] - M[1][2] * c[4] + M[1][3] * c[3]) * idet;
//    T[0][1] = (- M[0][1] * c[5] + M[0][2] * c[4] - M[0][3] * c[3]) * idet;
//    T[0][2] = (M[3][1] * s[5] - M[3][2] * s[4] + M[3][3] * s[3]) * idet;
//    T[0][3] = (- M[2][1] * s[5] + M[2][2] * s[4] - M[2][3] * s[3]) * idet;
//
//    T[1][0] = (- M[1][0] * c[5] + M[1][2] * c[2] - M[1][3] * c[1]) * idet;
//    T[1][1] = (M[0][0] * c[5] - M[0][2] * c[2] + M[0][3] * c[1]) * idet;
//    T[1][2] = (- M[3][0] * s[5] + M[3][2] * s[2] - M[3][3] * s[1]) * idet;
//    T[1][3] = (M[2][0] * s[5] - M[2][2] * s[2] + M[2][3] * s[1]) * idet;
//
//    T[2][0] = (M[1][0] * c[4] - M[1][1] * c[2] + M[1][3] * c[0]) * idet;
//    T[2][1] = (- M[0][0] * c[4] + M[0][1] * c[2] - M[0][3] * c[0]) * idet;
//    T[2][2] = (M[3][0] * s[4] - M[3][1] * s[2] + M[3][3] * s[0]) * idet;
//    T[2][3] = (- M[2][0] * s[4] + M[2][1] * s[2] - M[2][3] * s[0]) * idet;
//
//    T[3][0] = (- M[1][0] * c[3] + M[1][1] * c[1] - M[1][2] * c[0]) * idet;
//    T[3][1] = (M[0][0] * c[3] - M[0][1] * c[1] + M[0][2] * c[0]) * idet;
//    T[3][2] = (- M[3][0] * s[3] + M[3][1] * s[1] - M[3][2] * s[0]) * idet;
//    T[3][3] = (M[2][0] * s[3] - M[2][1] * s[1] + M[2][2] * s[0]) * idet;
//}
//
//static inline void mat4x4_orthonormalize(mat4x4 R, mat4x4 M) {
//    mat4x4_dup(R, M);
//    float s = 1.0f;
//    vec3 h;
//
//    vec3_normalize(R[2]);
//
//    s = vec3_mul_inner(R[1], R[2]);
//    vec3_scale_n(h, R[2], s);
//    vec3_sub_n(R[1], R[1], h);
//    vec3_normalize(R[2]);
//
//    s = vec3_mul_inner(R[1], R[2]);
//    vec3_scale_n(h, R[2], s);
//    vec3_sub_n(R[1], R[1], h);
//    vec3_normalize(R[1]);
//
//    s = vec3_mul_inner(R[0], R[1]);
//    vec3_scale_n(h, R[1], s);
//    vec3_sub_n(R[0], R[0], h);
//    vec3_normalize(R[0]);
//}

//static inline void mat4x4_frustum(mat4x4 M, float l, float r, float b, float t, float n, float f) {
//    M[0][0] = 2.0f * n / (r - l);
//    M[0][1] = M[0][2] = M[0][3] = 0.f;
//
//    M[1][1] = 2.0f * n / (t - b);
//    M[1][0] = M[1][2] = M[1][3] = 0.f;
//
//    M[2][0] = (r + l) / (r - l);
//    M[2][1] = (t + b) / (t - b);
//    M[2][2] = - (f + n) / (f - n);
//    M[2][3] = - 1.f;
//
//    M[3][2] = - 2.f * (f * n) / (f - n);
//    M[3][0] = M[3][1] = M[3][3] = 0.f;
//}
//
//static inline void mat4x4_ortho(mat4x4 M, float l, float r, float b, float t, float n, float f) {
//    M[0][0] = 2.f / (r - l);
//    M[0][1] = M[0][2] = M[0][3] = 0.f;
//
//    M[1][1] = 2.f / (t - b);
//    M[1][0] = M[1][2] = M[1][3] = 0.f;
//
//    M[2][2] = - 2.f / (f - n);
//    M[2][0] = M[2][1] = M[2][3] = 0.f;
//
//    M[3][0] = - (r + l) / (r - l);
//    M[3][1] = - (t + b) / (t - b);
//    M[3][2] = - (f + n) / (f - n);
//    M[3][3] = 1.f;
//}
//
//static inline void mat4x4_perspective(mat4x4 m, float y_fov, float aspect, float n, float f) {
//    /* NOTE: Degrees are an unhandy unit to work with.
//     * linmath.h uses radians for everything! */
//    float const a = 1.f / tan(y_fov / 2.f);
//
//    m[0][0] = a / aspect;
//    m[0][1] = 0.f;
//    m[0][2] = 0.f;
//    m[0][3] = 0.f;
//
//    m[1][0] = 0.f;
//    m[1][1] = a;
//    m[1][2] = 0.f;
//    m[1][3] = 0.f;
//
//    m[2][0] = 0.f;
//    m[2][1] = 0.f;
//    m[2][2] = - ((f + n) / (f - n));
//    m[2][3] = - 1.f;
//
//    m[3][0] = 0.f;
//    m[3][1] = 0.f;
//    m[3][2] = - ((2.f * f * n) / (f - n));
//    m[3][3] = 0.f;
//}
//
//static inline void mat4x4_look_at(mat4x4 m, vec3 eye, vec3 center, vec3 up) {
//    /* Adapted from Android's OpenGL Matrix.java.                        */
//    /* See the OpenGL GLUT documentation for gluLookAt for a description */
//    /* of the algorithm. We implement it in a straightforward way:       */
//
//    /* TODO: The negation of of can be spared by swapping the order of
//     *       operands in the following cross products in the right way. */
//    vec3 f;
//    vec3_sub_n(f, center, eye);
//    vec3_normalize(f);
//
//    vec3 s;
//    vec3_mul_cross_n(s, f, up);
//    vec3_normalize(s);
//
//    vec3 t;
//    vec3_mul_cross_n(t, s, f);
//
//    m[0][0] = s[0];
//    m[0][1] = t[0];
//    m[0][2] = - f[0];
//    m[0][3] = 0.f;
//
//    m[1][0] = s[1];
//    m[1][1] = t[1];
//    m[1][2] = - f[1];
//    m[1][3] = 0.f;
//
//    m[2][0] = s[2];
//    m[2][1] = t[2];
//    m[2][2] = - f[2];
//    m[2][3] = 0.f;
//
//    m[3][0] = 0.f;
//    m[3][1] = 0.f;
//    m[3][2] = 0.f;
//    m[3][3] = 1.f;
//
//    mat4x4_translate_in_place(m, - eye[0], - eye[1], - eye[2]);
//}
//
//typedef float quat[4];
//
//static inline void quat_identity(quat q) {
//    q[0] = q[1] = q[2] = 0.f;
//    q[3] = 1.f;
//}
//
//static inline void quat_add(quat r, quat a, quat b) {
//    int i;
//    for (i = 0; i < 4; ++ i)
//        r[i] = a[i] + b[i];
//}
//
//static inline void quat_sub(quat r, quat a, quat b) {
//    int i;
//    for (i = 0; i < 4; ++ i)
//        r[i] = a[i] - b[i];
//}
//
//static inline void quat_mul(quat r, quat p, quat q) {
//    vec3 w;
//    vec3_mul_cross_n(r, p, q);
//    vec3_scale_n(w, p, q[3]);
//    vec3_add_n(r, r, w);
//    vec3_scale_n(w, q, p[3]);
//    vec3_add_n(r, r, w);
//    r[3] = p[3] * q[3] - vec3_mul_inner(p, q);
//}
//
//static inline void quat_scale(quat r, quat v, float s) {
//    int i;
//    for (i = 0; i < 4; ++ i)
//        r[i] = v[i] * s;
//}
//
//static inline float quat_inner_product(quat a, quat b) {
//    float p = 0.f;
//    int i;
//    for (i = 0; i < 4; ++ i)
//        p += b[i] * a[i];
//    return p;
//}
//
//static inline void quat_conj(quat r, quat q) {
//    int i;
//    for (i = 0; i < 3; ++ i)
//        r[i] = - q[i];
//    r[3] = q[3];
//}
//
//static inline void quat_rotate(quat r, float angle, vec3 axis) {
//    vec3 v;
//    vec3_scale_n(v, axis, sinf(angle / 2));
//    int i;
//    for (i = 0; i < 3; ++ i)
//        r[i] = v[i];
//    r[3] = cosf(angle / 2);
//}
//
//#define quat_norm vec4_norm
//
//static inline void quat_mul_vec3(vec3 r, quat q, vec3 v) {
///*
// * Method by Fabian 'ryg' Giessen (of Farbrausch)
//t = 2 * cross(q.xyz, v)
//v' = v + q.w * t + cross(q.xyz, t)
// */
//    vec3 t;
//    vec3 q_xyz = {q[0], q[1], q[2]};
//    vec3 u = {q[0], q[1], q[2]};
//
//    vec3_mul_cross_n(t, q_xyz, v);
//    vec3_scale_n(t, t, 2);
//
//    vec3_mul_cross_n(u, q_xyz, t);
//    vec3_scale_n(t, t, q[3]);
//
//    vec3_add_n(r, v, t);
//    vec3_add_n(r, r, u);
//}
//
//static inline void mat4x4_from_quat(mat4x4 M, quat q) {
//    float a = q[3];
//    float b = q[0];
//    float c = q[1];
//    float d = q[2];
//    float a2 = a * a;
//    float b2 = b * b;
//    float c2 = c * c;
//    float d2 = d * d;
//
//    M[0][0] = a2 + b2 - c2 - d2;
//    M[0][1] = 2.f * (b * c + a * d);
//    M[0][2] = 2.f * (b * d - a * c);
//    M[0][3] = 0.f;
//
//    M[1][0] = 2 * (b * c - a * d);
//    M[1][1] = a2 - b2 + c2 - d2;
//    M[1][2] = 2.f * (c * d + a * b);
//    M[1][3] = 0.f;
//
//    M[2][0] = 2.f * (b * d + a * c);
//    M[2][1] = 2.f * (c * d - a * b);
//    M[2][2] = a2 - b2 - c2 + d2;
//    M[2][3] = 0.f;
//
//    M[3][0] = M[3][1] = M[3][2] = 0.f;
//    M[3][3] = 1.f;
//}
//
//static inline void mat4x4o_mul_quat(mat4x4 R, mat4x4 M, quat q) {
///*  XXX: The way this is written only works for othogonal matrices. */
///* TODO: Take care of non-orthogonal case. */
//    quat_mul_vec3(R[0], q, M[0]);
//    quat_mul_vec3(R[1], q, M[1]);
//    quat_mul_vec3(R[2], q, M[2]);
//
//    R[3][0] = R[3][1] = R[3][2] = 0.f;
//    R[3][3] = 1.f;
//}
//
//static inline void quat_from_mat4x4(quat q, mat4x4 M) {
//    float r = 0.f;
//    int i;
//
//    int perm[] = {0, 1, 2, 0, 1};
//    int *p = perm;
//
//    for (i = 0; i < 3; i ++) {
//        float m = M[i][i];
//        if (m < r)
//            continue;
//        m = r;
//        p = &perm[i];
//    }
//
//    r = sqrtf(1.f + M[p[0]][p[0]] - M[p[1]][p[1]] - M[p[2]][p[2]]);
//
//    if (r < 1e-6) {
//        q[0] = 1.f;
//        q[1] = q[2] = q[3] = 0.f;
//        return;
//    }
//
//    q[0] = r / 2.f;
//    q[1] = (M[p[0]][p[1]] - M[p[1]][p[0]]) / (2.f * r);
//    q[2] = (M[p[2]][p[0]] - M[p[0]][p[2]]) / (2.f * r);
//    q[3] = (M[p[2]][p[1]] - M[p[1]][p[2]]) / (2.f * r);
//}
//

#endif //TORQUECALCULATOR_LINEAR_ALGEBRA_H

#ifdef __cplusplus
}
#endif //__cplusplus
