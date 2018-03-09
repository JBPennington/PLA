
#ifdef __cplusplus
extern "C" {
#endif //__cplusplus

#ifndef PORTABLE_LINEAR_ALGEBRA_H
#define PORTABLE_LINEAR_ALGEBRA_H


#include <math.h>


#define LINMATH_H_DEFINE_VEC(n) \
typedef float vec##n[n]; \
static inline void vec##n##_zero(vec##n a) \
{ \
    int i; \
    for(i=0; i<(n); ++i) \
        a[i] = 0; \
} \
static inline void vec##n##_add_n(vec##n r, vec##n a, vec##n b) \
{ \
    int i; \
    for(i=0; i<(n); ++i) \
        r[i] = a[i] + b[i]; \
} \
static inline float* vec##n##_add_rn(vec##n r, vec##n a, vec##n b) \
{ \
    int i; \
    for(i=0; i<(n); ++i) \
        r[i] = a[i] + b[i]; \
    return r; \
} \
static inline float* vec##n##_add_r(vec##n a, vec##n b) \
{ \
    int i; \
    for(i=0; i<(n); ++i) \
        a[i] += b[i]; \
    return a; \
} \
static inline void vec##n##_add(vec##n a, vec##n b) \
{ \
    int i; \
    for(i=0; i<(n); ++i) \
        a[i] += b[i]; \
} \
static inline void vec##n##_sub_n(vec##n r, vec##n a, vec##n b) \
{ \
    int i; \
    for(i=0; i<(n); ++i) \
        r[i] = a[i] - b[i]; \
} \
static inline float* vec##n##_sub_r(vec##n a, vec##n b) \
{ \
    int i; \
    for(i=0; i<(n); ++i) \
        a[i] -= b[i]; \
    return a; \
} \
static inline void vec##n##_sub(vec##n a, vec##n b) \
{ \
    int i; \
    for(i=0; i<(n); ++i) \
        a[i] -= b[i]; \
} \
static inline void vec##n##_scale_n(vec##n r, vec##n v, float s) \
{ \
    int i; \
    for(i=0; i<(n); ++i) \
        r[i] = v[i] * s; \
} \
static inline float* vec##n##_scale_rn(vec##n r, vec##n v, float s) \
{ \
    int i; \
    for(i=0; i<(n); ++i) \
        r[i] = v[i] * s; \
    return r; \
} \
static inline float* vec##n##_scale_r(vec##n v, float s) \
{ \
    int i; \
    for(i=0; i<(n); ++i) \
        v[i] *= s; \
    return v; \
} \
static inline void vec##n##_scale(vec##n v, float s) \
{ \
    int i; \
    for(i=0; i<(n); ++i) \
        v[i] *= s; \
} \
static inline float vec##n##_mul_inner(vec##n a, vec##n b) \
{ \
    float p = 0.0f; \
    int i; \
    for(i=0; i<(n); ++i) \
        p += b[i]*a[i]; \
    return p; \
} \
static inline float vec##n##_norm(vec##n v) \
{ \
    return sqrtf(vec##n##_mul_inner(v,v)); \
} \
static inline void vec##n##_normalize_n(vec##n r, vec##n v) \
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
#define ZEROVEC3  {0.0f,0.0f,0.0f}

static inline void vec3_mul_cross_n(vec3 r, vec3 a, vec3 b) {
    r[0] = a[1] * b[2] - a[2] * b[1];
    r[1] = a[2] * b[0] - a[0] * b[2];
    r[2] = a[0] * b[1] - a[1] * b[0];
}

static inline float * vec3_mul_cross_rn(vec3 r, vec3 a, vec3 b) {
    r[0] = a[1] * b[2] - a[2] * b[1];
    r[1] = a[2] * b[0] - a[0] * b[2];
    r[2] = a[0] * b[1] - a[1] * b[0];
    return r;
}

static inline float * vec3_mul_cross_r(vec3 a, vec3 b) {
    vec3 result = {0,0,0};
    result[0] = a[1] * b[2] - a[2] * b[1];
    result[1] = a[2] * b[0] - a[0] * b[2];
    result[2] = a[0] * b[1] - a[1] * b[0];
    vec3_copy(a, result);
    return a;
}

static inline void vec3_mul_cross(vec3 a, vec3 b) {
    vec3 result = {0,0,0};
    result[0] = a[1] * b[2] - a[2] * b[1];
    result[1] = a[2] * b[0] - a[0] * b[2];
    result[2] = a[0] * b[1] - a[1] * b[0];
    vec3_copy(a, result);
}

static inline void vec3_reflect_n(vec3 r, vec3 v, vec3 n) {
    float p = 2.f * vec3_mul_inner(v, n);
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
typedef vec##n mat##n##x##n[n]; \
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
    vec3 v = EMPTYVEC3;
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
    float a0 = M[0][0];
    float a1 = M[1][0];
    float a2 = M[2][0];
    float a3 = M[3][0];
    float a4 = M[4][0];
    float a5 = M[5][0];
    float b0 = M[0][1];
    float b1 = M[1][1];
    float b2 = M[2][1];
    float b3 = M[3][1];
    float b4 = M[4][1];
    float b5 = M[5][1];
    float c0 = M[0][2];
    float c1 = M[1][2];
    float c2 = M[2][2];
    float c3 = M[3][2];
    float c4 = M[4][2];
    float c5 = M[5][2];
    float d0 = M[0][3];
    float d1 = M[1][3];
    float d2 = M[2][3];
    float d3 = M[3][3];
    float d4 = M[4][3];
    float d5 = M[5][3];
    float e0 = M[0][4];
    float e1 = M[1][4];
    float e2 = M[2][4];
    float e3 = M[3][4];
    float e4 = M[4][4];
    float e5 = M[5][4];
    float f0 = M[0][5];
    float f1 = M[1][5];
    float f2 = M[2][5];
    float f3 = M[3][5];
    float f4 = M[4][5];
    float f5 = M[5][5];
    float t2 = a0*b1*c2*d3*e4*f5;
    float t3 = a0*b1*c2*d4*e5*f3;
    float t4 = a0*b1*c2*d5*e3*f4;
    float t5 = a0*b1*c3*d2*e5*f4;
    float t6 = a0*b1*c3*d4*e2*f5;
    float t7 = a0*b1*c3*d5*e4*f2;
    float t8 = a0*b1*c4*d2*e3*f5;
    float t9 = a0*b1*c4*d3*e5*f2;
    float t10 = a0*b1*c4*d5*e2*f3;
    float t11 = a0*b1*c5*d2*e4*f3;
    float t12 = a0*b1*c5*d3*e2*f4;
    float t13 = a0*b1*c5*d4*e3*f2;
    float t14 = a0*b2*c1*d3*e5*f4;
    float t15 = a0*b2*c1*d4*e3*f5;
    float t16 = a0*b2*c1*d5*e4*f3;
    float t17 = a0*b2*c3*d1*e4*f5;
    float t18 = a0*b2*c3*d4*e5*f1;
    float t19 = a0*b2*c3*d5*e1*f4;
    float t20 = a0*b2*c4*d1*e5*f3;
    float t21 = a0*b2*c4*d3*e1*f5;
    float t22 = a0*b2*c4*d5*e3*f1;
    float t23 = a0*b2*c5*d1*e3*f4;
    float t24 = a0*b2*c5*d3*e4*f1;
    float t25 = a0*b2*c5*d4*e1*f3;
    float t26 = a0*b3*c1*d2*e4*f5;
    float t27 = a0*b3*c1*d4*e5*f2;
    float t28 = a0*b3*c1*d5*e2*f4;
    float t29 = a0*b3*c2*d1*e5*f4;
    float t30 = a0*b3*c2*d4*e1*f5;
    float t31 = a0*b3*c2*d5*e4*f1;
    float t32 = a0*b3*c4*d1*e2*f5;
    float t33 = a0*b3*c4*d2*e5*f1;
    float t34 = a0*b3*c4*d5*e1*f2;
    float t35 = a0*b3*c5*d1*e4*f2;
    float t36 = a0*b3*c5*d2*e1*f4;
    float t37 = a0*b3*c5*d4*e2*f1;
    float t38 = a0*b4*c1*d2*e5*f3;
    float t39 = a0*b4*c1*d3*e2*f5;
    float t40 = a0*b4*c1*d5*e3*f2;
    float t41 = a0*b4*c2*d1*e3*f5;
    float t42 = a0*b4*c2*d3*e5*f1;
    float t43 = a0*b4*c2*d5*e1*f3;
    float t44 = a0*b4*c3*d1*e5*f2;
    float t45 = a0*b4*c3*d2*e1*f5;
    float t46 = a0*b4*c3*d5*e2*f1;
    float t47 = a0*b4*c5*d1*e2*f3;
    float t48 = a0*b4*c5*d2*e3*f1;
    float t49 = a0*b4*c5*d3*e1*f2;
    float t50 = a0*b5*c1*d2*e3*f4;
    float t51 = a0*b5*c1*d3*e4*f2;
    float t52 = a0*b5*c1*d4*e2*f3;
    float t53 = a0*b5*c2*d1*e4*f3;
    float t54 = a0*b5*c2*d3*e1*f4;
    float t55 = a0*b5*c2*d4*e3*f1;
    float t56 = a0*b5*c3*d1*e2*f4;
    float t57 = a0*b5*c3*d2*e4*f1;
    float t58 = a0*b5*c3*d4*e1*f2;
    float t59 = a0*b5*c4*d1*e3*f2;
    float t60 = a0*b5*c4*d2*e1*f3;
    float t61 = a0*b5*c4*d3*e2*f1;
    float t62 = a1*b0*c2*d3*e5*f4;
    float t63 = a1*b0*c2*d4*e3*f5;
    float t64 = a1*b0*c2*d5*e4*f3;
    float t65 = a1*b0*c3*d2*e4*f5;
    float t66 = a1*b0*c3*d4*e5*f2;
    float t67 = a1*b0*c3*d5*e2*f4;
    float t68 = a1*b0*c4*d2*e5*f3;
    float t69 = a1*b0*c4*d3*e2*f5;
    float t70 = a1*b0*c4*d5*e3*f2;
    float t71 = a1*b0*c5*d2*e3*f4;
    float t72 = a1*b0*c5*d3*e4*f2;
    float t73 = a1*b0*c5*d4*e2*f3;
    float t74 = a1*b2*c0*d3*e4*f5;
    float t75 = a1*b2*c0*d4*e5*f3;
    float t76 = a1*b2*c0*d5*e3*f4;
    float t77 = a1*b2*c3*d0*e5*f4;
    float t78 = a1*b2*c3*d4*e0*f5;
    float t79 = a1*b2*c3*d5*e4*f0;
    float t80 = a1*b2*c4*d0*e3*f5;
    float t81 = a1*b2*c4*d3*e5*f0;
    float t82 = a1*b2*c4*d5*e0*f3;
    float t83 = a1*b2*c5*d0*e4*f3;
    float t84 = a1*b2*c5*d3*e0*f4;
    float t85 = a1*b2*c5*d4*e3*f0;
    float t86 = a1*b3*c0*d2*e5*f4;
    float t87 = a1*b3*c0*d4*e2*f5;
    float t88 = a1*b3*c0*d5*e4*f2;
    float t89 = a1*b3*c2*d0*e4*f5;
    float t90 = a1*b3*c2*d4*e5*f0;
    float t91 = a1*b3*c2*d5*e0*f4;
    float t92 = a1*b3*c4*d0*e5*f2;
    float t93 = a1*b3*c4*d2*e0*f5;
    float t94 = a1*b3*c4*d5*e2*f0;
    float t95 = a1*b3*c5*d0*e2*f4;
    float t96 = a1*b3*c5*d2*e4*f0;
    float t97 = a1*b3*c5*d4*e0*f2;
    float t98 = a1*b4*c0*d2*e3*f5;
    float t99 = a1*b4*c0*d3*e5*f2;
    float t100 = a1*b4*c0*d5*e2*f3;
    float t101 = a1*b4*c2*d0*e5*f3;
    float t102 = a1*b4*c2*d3*e0*f5;
    float t103 = a1*b4*c2*d5*e3*f0;
    float t104 = a1*b4*c3*d0*e2*f5;
    float t105 = a1*b4*c3*d2*e5*f0;
    float t106 = a1*b4*c3*d5*e0*f2;
    float t107 = a1*b4*c5*d0*e3*f2;
    float t108 = a1*b4*c5*d2*e0*f3;
    float t109 = a1*b4*c5*d3*e2*f0;
    float t110 = a1*b5*c0*d2*e4*f3;
    float t111 = a1*b5*c0*d3*e2*f4;
    float t112 = a1*b5*c0*d4*e3*f2;
    float t113 = a1*b5*c2*d0*e3*f4;
    float t114 = a1*b5*c2*d3*e4*f0;
    float t115 = a1*b5*c2*d4*e0*f3;
    float t116 = a1*b5*c3*d0*e4*f2;
    float t117 = a1*b5*c3*d2*e0*f4;
    float t118 = a1*b5*c3*d4*e2*f0;
    float t119 = a1*b5*c4*d0*e2*f3;
    float t120 = a1*b5*c4*d2*e3*f0;
    float t121 = a1*b5*c4*d3*e0*f2;
    float t122 = a2*b0*c1*d3*e4*f5;
    float t123 = a2*b0*c1*d4*e5*f3;
    float t124 = a2*b0*c1*d5*e3*f4;
    float t125 = a2*b0*c3*d1*e5*f4;
    float t126 = a2*b0*c3*d4*e1*f5;
    float t127 = a2*b0*c3*d5*e4*f1;
    float t128 = a2*b0*c4*d1*e3*f5;
    float t129 = a2*b0*c4*d3*e5*f1;
    float t130 = a2*b0*c4*d5*e1*f3;
    float t131 = a2*b0*c5*d1*e4*f3;
    float t132 = a2*b0*c5*d3*e1*f4;
    float t133 = a2*b0*c5*d4*e3*f1;
    float t134 = a2*b1*c0*d3*e5*f4;
    float t135 = a2*b1*c0*d4*e3*f5;
    float t136 = a2*b1*c0*d5*e4*f3;
    float t137 = a2*b1*c3*d0*e4*f5;
    float t138 = a2*b1*c3*d4*e5*f0;
    float t139 = a2*b1*c3*d5*e0*f4;
    float t140 = a2*b1*c4*d0*e5*f3;
    float t141 = a2*b1*c4*d3*e0*f5;
    float t142 = a2*b1*c4*d5*e3*f0;
    float t143 = a2*b1*c5*d0*e3*f4;
    float t144 = a2*b1*c5*d3*e4*f0;
    float t145 = a2*b1*c5*d4*e0*f3;
    float t146 = a2*b3*c0*d1*e4*f5;
    float t147 = a2*b3*c0*d4*e5*f1;
    float t148 = a2*b3*c0*d5*e1*f4;
    float t149 = a2*b3*c1*d0*e5*f4;
    float t150 = a2*b3*c1*d4*e0*f5;
    float t151 = a2*b3*c1*d5*e4*f0;
    float t152 = a2*b3*c4*d0*e1*f5;
    float t153 = a2*b3*c4*d1*e5*f0;
    float t154 = a2*b3*c4*d5*e0*f1;
    float t155 = a2*b3*c5*d0*e4*f1;
    float t156 = a2*b3*c5*d1*e0*f4;
    float t157 = a2*b3*c5*d4*e1*f0;
    float t158 = a2*b4*c0*d1*e5*f3;
    float t159 = a2*b4*c0*d3*e1*f5;
    float t160 = a2*b4*c0*d5*e3*f1;
    float t161 = a2*b4*c1*d0*e3*f5;
    float t162 = a2*b4*c1*d3*e5*f0;
    float t163 = a2*b4*c1*d5*e0*f3;
    float t164 = a2*b4*c3*d0*e5*f1;
    float t165 = a2*b4*c3*d1*e0*f5;
    float t166 = a2*b4*c3*d5*e1*f0;
    float t167 = a2*b4*c5*d0*e1*f3;
    float t168 = a2*b4*c5*d1*e3*f0;
    float t169 = a2*b4*c5*d3*e0*f1;
    float t170 = a2*b5*c0*d1*e3*f4;
    float t171 = a2*b5*c0*d3*e4*f1;
    float t172 = a2*b5*c0*d4*e1*f3;
    float t173 = a2*b5*c1*d0*e4*f3;
    float t174 = a2*b5*c1*d3*e0*f4;
    float t175 = a2*b5*c1*d4*e3*f0;
    float t176 = a2*b5*c3*d0*e1*f4;
    float t177 = a2*b5*c3*d1*e4*f0;
    float t178 = a2*b5*c3*d4*e0*f1;
    float t179 = a2*b5*c4*d0*e3*f1;
    float t180 = a2*b5*c4*d1*e0*f3;
    float t181 = a2*b5*c4*d3*e1*f0;
    float t182 = a3*b0*c1*d2*e5*f4;
    float t183 = a3*b0*c1*d4*e2*f5;
    float t184 = a3*b0*c1*d5*e4*f2;
    float t185 = a3*b0*c2*d1*e4*f5;
    float t186 = a3*b0*c2*d4*e5*f1;
    float t187 = a3*b0*c2*d5*e1*f4;
    float t188 = a3*b0*c4*d1*e5*f2;
    float t189 = a3*b0*c4*d2*e1*f5;
    float t190 = a3*b0*c4*d5*e2*f1;
    float t191 = a3*b0*c5*d1*e2*f4;
    float t192 = a3*b0*c5*d2*e4*f1;
    float t193 = a3*b0*c5*d4*e1*f2;
    float t194 = a3*b1*c0*d2*e4*f5;
    float t195 = a3*b1*c0*d4*e5*f2;
    float t196 = a3*b1*c0*d5*e2*f4;
    float t197 = a3*b1*c2*d0*e5*f4;
    float t198 = a3*b1*c2*d4*e0*f5;
    float t199 = a3*b1*c2*d5*e4*f0;
    float t200 = a3*b1*c4*d0*e2*f5;
    float t201 = a3*b1*c4*d2*e5*f0;
    float t202 = a3*b1*c4*d5*e0*f2;
    float t203 = a3*b1*c5*d0*e4*f2;
    float t204 = a3*b1*c5*d2*e0*f4;
    float t205 = a3*b1*c5*d4*e2*f0;
    float t206 = a3*b2*c0*d1*e5*f4;
    float t207 = a3*b2*c0*d4*e1*f5;
    float t208 = a3*b2*c0*d5*e4*f1;
    float t209 = a3*b2*c1*d0*e4*f5;
    float t210 = a3*b2*c1*d4*e5*f0;
    float t211 = a3*b2*c1*d5*e0*f4;
    float t212 = a3*b2*c4*d0*e5*f1;
    float t213 = a3*b2*c4*d1*e0*f5;
    float t214 = a3*b2*c4*d5*e1*f0;
    float t215 = a3*b2*c5*d0*e1*f4;
    float t216 = a3*b2*c5*d1*e4*f0;
    float t217 = a3*b2*c5*d4*e0*f1;
    float t218 = a3*b4*c0*d1*e2*f5;
    float t219 = a3*b4*c0*d2*e5*f1;
    float t220 = a3*b4*c0*d5*e1*f2;
    float t221 = a3*b4*c1*d0*e5*f2;
    float t222 = a3*b4*c1*d2*e0*f5;
    float t223 = a3*b4*c1*d5*e2*f0;
    float t224 = a3*b4*c2*d0*e1*f5;
    float t225 = a3*b4*c2*d1*e5*f0;
    float t226 = a3*b4*c2*d5*e0*f1;
    float t227 = a3*b4*c5*d0*e2*f1;
    float t228 = a3*b4*c5*d1*e0*f2;
    float t229 = a3*b4*c5*d2*e1*f0;
    float t230 = a3*b5*c0*d1*e4*f2;
    float t231 = a3*b5*c0*d2*e1*f4;
    float t232 = a3*b5*c0*d4*e2*f1;
    float t233 = a3*b5*c1*d0*e2*f4;
    float t234 = a3*b5*c1*d2*e4*f0;
    float t235 = a3*b5*c1*d4*e0*f2;
    float t236 = a3*b5*c2*d0*e4*f1;
    float t237 = a3*b5*c2*d1*e0*f4;
    float t238 = a3*b5*c2*d4*e1*f0;
    float t239 = a3*b5*c4*d0*e1*f2;
    float t240 = a3*b5*c4*d1*e2*f0;
    float t241 = a3*b5*c4*d2*e0*f1;
    float t242 = a4*b0*c1*d2*e3*f5;
    float t243 = a4*b0*c1*d3*e5*f2;
    float t244 = a4*b0*c1*d5*e2*f3;
    float t245 = a4*b0*c2*d1*e5*f3;
    float t246 = a4*b0*c2*d3*e1*f5;
    float t247 = a4*b0*c2*d5*e3*f1;
    float t248 = a4*b0*c3*d1*e2*f5;
    float t249 = a4*b0*c3*d2*e5*f1;
    float t250 = a4*b0*c3*d5*e1*f2;
    float t251 = a4*b0*c5*d1*e3*f2;
    float t252 = a4*b0*c5*d2*e1*f3;
    float t253 = a4*b0*c5*d3*e2*f1;
    float t254 = a4*b1*c0*d2*e5*f3;
    float t255 = a4*b1*c0*d3*e2*f5;
    float t256 = a4*b1*c0*d5*e3*f2;
    float t257 = a4*b1*c2*d0*e3*f5;
    float t258 = a4*b1*c2*d3*e5*f0;
    float t259 = a4*b1*c2*d5*e0*f3;
    float t260 = a4*b1*c3*d0*e5*f2;
    float t261 = a4*b1*c3*d2*e0*f5;
    float t262 = a4*b1*c3*d5*e2*f0;
    float t263 = a4*b1*c5*d0*e2*f3;
    float t264 = a4*b1*c5*d2*e3*f0;
    float t265 = a4*b1*c5*d3*e0*f2;
    float t266 = a4*b2*c0*d1*e3*f5;
    float t267 = a4*b2*c0*d3*e5*f1;
    float t268 = a4*b2*c0*d5*e1*f3;
    float t269 = a4*b2*c1*d0*e5*f3;
    float t270 = a4*b2*c1*d3*e0*f5;
    float t271 = a4*b2*c1*d5*e3*f0;
    float t272 = a4*b2*c3*d0*e1*f5;
    float t273 = a4*b2*c3*d1*e5*f0;
    float t274 = a4*b2*c3*d5*e0*f1;
    float t275 = a4*b2*c5*d0*e3*f1;
    float t276 = a4*b2*c5*d1*e0*f3;
    float t277 = a4*b2*c5*d3*e1*f0;
    float t278 = a4*b3*c0*d1*e5*f2;
    float t279 = a4*b3*c0*d2*e1*f5;
    float t280 = a4*b3*c0*d5*e2*f1;
    float t281 = a4*b3*c1*d0*e2*f5;
    float t282 = a4*b3*c1*d2*e5*f0;
    float t283 = a4*b3*c1*d5*e0*f2;
    float t284 = a4*b3*c2*d0*e5*f1;
    float t285 = a4*b3*c2*d1*e0*f5;
    float t286 = a4*b3*c2*d5*e1*f0;
    float t287 = a4*b3*c5*d0*e1*f2;
    float t288 = a4*b3*c5*d1*e2*f0;
    float t289 = a4*b3*c5*d2*e0*f1;
    float t290 = a4*b5*c0*d1*e2*f3;
    float t291 = a4*b5*c0*d2*e3*f1;
    float t292 = a4*b5*c0*d3*e1*f2;
    float t293 = a4*b5*c1*d0*e3*f2;
    float t294 = a4*b5*c1*d2*e0*f3;
    float t295 = a4*b5*c1*d3*e2*f0;
    float t296 = a4*b5*c2*d0*e1*f3;
    float t297 = a4*b5*c2*d1*e3*f0;
    float t298 = a4*b5*c2*d3*e0*f1;
    float t299 = a4*b5*c3*d0*e2*f1;
    float t300 = a4*b5*c3*d1*e0*f2;
    float t301 = a4*b5*c3*d2*e1*f0;
    float t302 = a5*b0*c1*d2*e4*f3;
    float t303 = a5*b0*c1*d3*e2*f4;
    float t304 = a5*b0*c1*d4*e3*f2;
    float t305 = a5*b0*c2*d1*e3*f4;
    float t306 = a5*b0*c2*d3*e4*f1;
    float t307 = a5*b0*c2*d4*e1*f3;
    float t308 = a5*b0*c3*d1*e4*f2;
    float t309 = a5*b0*c3*d2*e1*f4;
    float t310 = a5*b0*c3*d4*e2*f1;
    float t311 = a5*b0*c4*d1*e2*f3;
    float t312 = a5*b0*c4*d2*e3*f1;
    float t313 = a5*b0*c4*d3*e1*f2;
    float t314 = a5*b1*c0*d2*e3*f4;
    float t315 = a5*b1*c0*d3*e4*f2;
    float t316 = a5*b1*c0*d4*e2*f3;
    float t317 = a5*b1*c2*d0*e4*f3;
    float t318 = a5*b1*c2*d3*e0*f4;
    float t319 = a5*b1*c2*d4*e3*f0;
    float t320 = a5*b1*c3*d0*e2*f4;
    float t321 = a5*b1*c3*d2*e4*f0;
    float t322 = a5*b1*c3*d4*e0*f2;
    float t323 = a5*b1*c4*d0*e3*f2;
    float t324 = a5*b1*c4*d2*e0*f3;
    float t325 = a5*b1*c4*d3*e2*f0;
    float t326 = a5*b2*c0*d1*e4*f3;
    float t327 = a5*b2*c0*d3*e1*f4;
    float t328 = a5*b2*c0*d4*e3*f1;
    float t329 = a5*b2*c1*d0*e3*f4;
    float t330 = a5*b2*c1*d3*e4*f0;
    float t331 = a5*b2*c1*d4*e0*f3;
    float t332 = a5*b2*c3*d0*e4*f1;
    float t333 = a5*b2*c3*d1*e0*f4;
    float t334 = a5*b2*c3*d4*e1*f0;
    float t335 = a5*b2*c4*d0*e1*f3;
    float t336 = a5*b2*c4*d1*e3*f0;
    float t337 = a5*b2*c4*d3*e0*f1;
    float t338 = a5*b3*c0*d1*e2*f4;
    float t339 = a5*b3*c0*d2*e4*f1;
    float t340 = a5*b3*c0*d4*e1*f2;
    float t341 = a5*b3*c1*d0*e4*f2;
    float t342 = a5*b3*c1*d2*e0*f4;
    float t343 = a5*b3*c1*d4*e2*f0;
    float t344 = a5*b3*c2*d0*e1*f4;
    float t345 = a5*b3*c2*d1*e4*f0;
    float t346 = a5*b3*c2*d4*e0*f1;
    float t347 = a5*b3*c4*d0*e2*f1;
    float t348 = a5*b3*c4*d1*e0*f2;
    float t349 = a5*b3*c4*d2*e1*f0;
    float t350 = a5*b4*c0*d1*e3*f2;
    float t351 = a5*b4*c0*d2*e1*f3;
    float t352 = a5*b4*c0*d3*e2*f1;
    float t353 = a5*b4*c1*d0*e2*f3;
    float t354 = a5*b4*c1*d2*e3*f0;
    float t355 = a5*b4*c1*d3*e0*f2;
    float t356 = a5*b4*c2*d0*e3*f1;
    float t357 = a5*b4*c2*d1*e0*f3;
    float t358 = a5*b4*c2*d3*e1*f0;
    float t359 = a5*b4*c3*d0*e1*f2;
    float t360 = a5*b4*c3*d1*e2*f0;
    float t361 = a5*b4*c3*d2*e0*f1;
    float t364 = a0*b1*c2*d3*e5*f4;
    float t365 = a0*b1*c2*d4*e3*f5;
    float t366 = a0*b1*c2*d5*e4*f3;
    float t367 = a0*b1*c3*d2*e4*f5;
    float t368 = a0*b1*c3*d4*e5*f2;
    float t369 = a0*b1*c3*d5*e2*f4;
    float t370 = a0*b1*c4*d2*e5*f3;
    float t371 = a0*b1*c4*d3*e2*f5;
    float t372 = a0*b1*c4*d5*e3*f2;
    float t373 = a0*b1*c5*d2*e3*f4;
    float t374 = a0*b1*c5*d3*e4*f2;
    float t375 = a0*b1*c5*d4*e2*f3;
    float t376 = a0*b2*c1*d3*e4*f5;
    float t377 = a0*b2*c1*d4*e5*f3;
    float t378 = a0*b2*c1*d5*e3*f4;
    float t379 = a0*b2*c3*d1*e5*f4;
    float t380 = a0*b2*c3*d4*e1*f5;
    float t381 = a0*b2*c3*d5*e4*f1;
    float t382 = a0*b2*c4*d1*e3*f5;
    float t383 = a0*b2*c4*d3*e5*f1;
    float t384 = a0*b2*c4*d5*e1*f3;
    float t385 = a0*b2*c5*d1*e4*f3;
    float t386 = a0*b2*c5*d3*e1*f4;
    float t387 = a0*b2*c5*d4*e3*f1;
    float t388 = a0*b3*c1*d2*e5*f4;
    float t389 = a0*b3*c1*d4*e2*f5;
    float t390 = a0*b3*c1*d5*e4*f2;
    float t391 = a0*b3*c2*d1*e4*f5;
    float t392 = a0*b3*c2*d4*e5*f1;
    float t393 = a0*b3*c2*d5*e1*f4;
    float t394 = a0*b3*c4*d1*e5*f2;
    float t395 = a0*b3*c4*d2*e1*f5;
    float t396 = a0*b3*c4*d5*e2*f1;
    float t397 = a0*b3*c5*d1*e2*f4;
    float t398 = a0*b3*c5*d2*e4*f1;
    float t399 = a0*b3*c5*d4*e1*f2;
    float t400 = a0*b4*c1*d2*e3*f5;
    float t401 = a0*b4*c1*d3*e5*f2;
    float t402 = a0*b4*c1*d5*e2*f3;
    float t403 = a0*b4*c2*d1*e5*f3;
    float t404 = a0*b4*c2*d3*e1*f5;
    float t405 = a0*b4*c2*d5*e3*f1;
    float t406 = a0*b4*c3*d1*e2*f5;
    float t407 = a0*b4*c3*d2*e5*f1;
    float t408 = a0*b4*c3*d5*e1*f2;
    float t409 = a0*b4*c5*d1*e3*f2;
    float t410 = a0*b4*c5*d2*e1*f3;
    float t411 = a0*b4*c5*d3*e2*f1;
    float t412 = a0*b5*c1*d2*e4*f3;
    float t413 = a0*b5*c1*d3*e2*f4;
    float t414 = a0*b5*c1*d4*e3*f2;
    float t415 = a0*b5*c2*d1*e3*f4;
    float t416 = a0*b5*c2*d3*e4*f1;
    float t417 = a0*b5*c2*d4*e1*f3;
    float t418 = a0*b5*c3*d1*e4*f2;
    float t419 = a0*b5*c3*d2*e1*f4;
    float t420 = a0*b5*c3*d4*e2*f1;
    float t421 = a0*b5*c4*d1*e2*f3;
    float t422 = a0*b5*c4*d2*e3*f1;
    float t423 = a0*b5*c4*d3*e1*f2;
    float t424 = a1*b0*c2*d3*e4*f5;
    float t425 = a1*b0*c2*d4*e5*f3;
    float t426 = a1*b0*c2*d5*e3*f4;
    float t427 = a1*b0*c3*d2*e5*f4;
    float t428 = a1*b0*c3*d4*e2*f5;
    float t429 = a1*b0*c3*d5*e4*f2;
    float t430 = a1*b0*c4*d2*e3*f5;
    float t431 = a1*b0*c4*d3*e5*f2;
    float t432 = a1*b0*c4*d5*e2*f3;
    float t433 = a1*b0*c5*d2*e4*f3;
    float t434 = a1*b0*c5*d3*e2*f4;
    float t435 = a1*b0*c5*d4*e3*f2;
    float t436 = a1*b2*c0*d3*e5*f4;
    float t437 = a1*b2*c0*d4*e3*f5;
    float t438 = a1*b2*c0*d5*e4*f3;
    float t439 = a1*b2*c3*d0*e4*f5;
    float t440 = a1*b2*c3*d4*e5*f0;
    float t441 = a1*b2*c3*d5*e0*f4;
    float t442 = a1*b2*c4*d0*e5*f3;
    float t443 = a1*b2*c4*d3*e0*f5;
    float t444 = a1*b2*c4*d5*e3*f0;
    float t445 = a1*b2*c5*d0*e3*f4;
    float t446 = a1*b2*c5*d3*e4*f0;
    float t447 = a1*b2*c5*d4*e0*f3;
    float t448 = a1*b3*c0*d2*e4*f5;
    float t449 = a1*b3*c0*d4*e5*f2;
    float t450 = a1*b3*c0*d5*e2*f4;
    float t451 = a1*b3*c2*d0*e5*f4;
    float t452 = a1*b3*c2*d4*e0*f5;
    float t453 = a1*b3*c2*d5*e4*f0;
    float t454 = a1*b3*c4*d0*e2*f5;
    float t455 = a1*b3*c4*d2*e5*f0;
    float t456 = a1*b3*c4*d5*e0*f2;
    float t457 = a1*b3*c5*d0*e4*f2;
    float t458 = a1*b3*c5*d2*e0*f4;
    float t459 = a1*b3*c5*d4*e2*f0;
    float t460 = a1*b4*c0*d2*e5*f3;
    float t461 = a1*b4*c0*d3*e2*f5;
    float t462 = a1*b4*c0*d5*e3*f2;
    float t463 = a1*b4*c2*d0*e3*f5;
    float t464 = a1*b4*c2*d3*e5*f0;
    float t465 = a1*b4*c2*d5*e0*f3;
    float t466 = a1*b4*c3*d0*e5*f2;
    float t467 = a1*b4*c3*d2*e0*f5;
    float t468 = a1*b4*c3*d5*e2*f0;
    float t469 = a1*b4*c5*d0*e2*f3;
    float t470 = a1*b4*c5*d2*e3*f0;
    float t471 = a1*b4*c5*d3*e0*f2;
    float t472 = a1*b5*c0*d2*e3*f4;
    float t473 = a1*b5*c0*d3*e4*f2;
    float t474 = a1*b5*c0*d4*e2*f3;
    float t475 = a1*b5*c2*d0*e4*f3;
    float t476 = a1*b5*c2*d3*e0*f4;
    float t477 = a1*b5*c2*d4*e3*f0;
    float t478 = a1*b5*c3*d0*e2*f4;
    float t479 = a1*b5*c3*d2*e4*f0;
    float t480 = a1*b5*c3*d4*e0*f2;
    float t481 = a1*b5*c4*d0*e3*f2;
    float t482 = a1*b5*c4*d2*e0*f3;
    float t483 = a1*b5*c4*d3*e2*f0;
    float t484 = a2*b0*c1*d3*e5*f4;
    float t485 = a2*b0*c1*d4*e3*f5;
    float t486 = a2*b0*c1*d5*e4*f3;
    float t487 = a2*b0*c3*d1*e4*f5;
    float t488 = a2*b0*c3*d4*e5*f1;
    float t489 = a2*b0*c3*d5*e1*f4;
    float t490 = a2*b0*c4*d1*e5*f3;
    float t491 = a2*b0*c4*d3*e1*f5;
    float t492 = a2*b0*c4*d5*e3*f1;
    float t493 = a2*b0*c5*d1*e3*f4;
    float t494 = a2*b0*c5*d3*e4*f1;
    float t495 = a2*b0*c5*d4*e1*f3;
    float t496 = a2*b1*c0*d3*e4*f5;
    float t497 = a2*b1*c0*d4*e5*f3;
    float t498 = a2*b1*c0*d5*e3*f4;
    float t499 = a2*b1*c3*d0*e5*f4;
    float t500 = a2*b1*c3*d4*e0*f5;
    float t501 = a2*b1*c3*d5*e4*f0;
    float t502 = a2*b1*c4*d0*e3*f5;
    float t503 = a2*b1*c4*d3*e5*f0;
    float t504 = a2*b1*c4*d5*e0*f3;
    float t505 = a2*b1*c5*d0*e4*f3;
    float t506 = a2*b1*c5*d3*e0*f4;
    float t507 = a2*b1*c5*d4*e3*f0;
    float t508 = a2*b3*c0*d1*e5*f4;
    float t509 = a2*b3*c0*d4*e1*f5;
    float t510 = a2*b3*c0*d5*e4*f1;
    float t511 = a2*b3*c1*d0*e4*f5;
    float t512 = a2*b3*c1*d4*e5*f0;
    float t513 = a2*b3*c1*d5*e0*f4;
    float t514 = a2*b3*c4*d0*e5*f1;
    float t515 = a2*b3*c4*d1*e0*f5;
    float t516 = a2*b3*c4*d5*e1*f0;
    float t517 = a2*b3*c5*d0*e1*f4;
    float t518 = a2*b3*c5*d1*e4*f0;
    float t519 = a2*b3*c5*d4*e0*f1;
    float t520 = a2*b4*c0*d1*e3*f5;
    float t521 = a2*b4*c0*d3*e5*f1;
    float t522 = a2*b4*c0*d5*e1*f3;
    float t523 = a2*b4*c1*d0*e5*f3;
    float t524 = a2*b4*c1*d3*e0*f5;
    float t525 = a2*b4*c1*d5*e3*f0;
    float t526 = a2*b4*c3*d0*e1*f5;
    float t527 = a2*b4*c3*d1*e5*f0;
    float t528 = a2*b4*c3*d5*e0*f1;
    float t529 = a2*b4*c5*d0*e3*f1;
    float t530 = a2*b4*c5*d1*e0*f3;
    float t531 = a2*b4*c5*d3*e1*f0;
    float t532 = a2*b5*c0*d1*e4*f3;
    float t533 = a2*b5*c0*d3*e1*f4;
    float t534 = a2*b5*c0*d4*e3*f1;
    float t535 = a2*b5*c1*d0*e3*f4;
    float t536 = a2*b5*c1*d3*e4*f0;
    float t537 = a2*b5*c1*d4*e0*f3;
    float t538 = a2*b5*c3*d0*e4*f1;
    float t539 = a2*b5*c3*d1*e0*f4;
    float t540 = a2*b5*c3*d4*e1*f0;
    float t541 = a2*b5*c4*d0*e1*f3;
    float t542 = a2*b5*c4*d1*e3*f0;
    float t543 = a2*b5*c4*d3*e0*f1;
    float t544 = a3*b0*c1*d2*e4*f5;
    float t545 = a3*b0*c1*d4*e5*f2;
    float t546 = a3*b0*c1*d5*e2*f4;
    float t547 = a3*b0*c2*d1*e5*f4;
    float t548 = a3*b0*c2*d4*e1*f5;
    float t549 = a3*b0*c2*d5*e4*f1;
    float t550 = a3*b0*c4*d1*e2*f5;
    float t551 = a3*b0*c4*d2*e5*f1;
    float t552 = a3*b0*c4*d5*e1*f2;
    float t553 = a3*b0*c5*d1*e4*f2;
    float t554 = a3*b0*c5*d2*e1*f4;
    float t555 = a3*b0*c5*d4*e2*f1;
    float t556 = a3*b1*c0*d2*e5*f4;
    float t557 = a3*b1*c0*d4*e2*f5;
    float t558 = a3*b1*c0*d5*e4*f2;
    float t559 = a3*b1*c2*d0*e4*f5;
    float t560 = a3*b1*c2*d4*e5*f0;
    float t561 = a3*b1*c2*d5*e0*f4;
    float t562 = a3*b1*c4*d0*e5*f2;
    float t563 = a3*b1*c4*d2*e0*f5;
    float t564 = a3*b1*c4*d5*e2*f0;
    float t565 = a3*b1*c5*d0*e2*f4;
    float t566 = a3*b1*c5*d2*e4*f0;
    float t567 = a3*b1*c5*d4*e0*f2;
    float t568 = a3*b2*c0*d1*e4*f5;
    float t569 = a3*b2*c0*d4*e5*f1;
    float t570 = a3*b2*c0*d5*e1*f4;
    float t571 = a3*b2*c1*d0*e5*f4;
    float t572 = a3*b2*c1*d4*e0*f5;
    float t573 = a3*b2*c1*d5*e4*f0;
    float t574 = a3*b2*c4*d0*e1*f5;
    float t575 = a3*b2*c4*d1*e5*f0;
    float t576 = a3*b2*c4*d5*e0*f1;
    float t577 = a3*b2*c5*d0*e4*f1;
    float t578 = a3*b2*c5*d1*e0*f4;
    float t579 = a3*b2*c5*d4*e1*f0;
    float t580 = a3*b4*c0*d1*e5*f2;
    float t581 = a3*b4*c0*d2*e1*f5;
    float t582 = a3*b4*c0*d5*e2*f1;
    float t583 = a3*b4*c1*d0*e2*f5;
    float t584 = a3*b4*c1*d2*e5*f0;
    float t585 = a3*b4*c1*d5*e0*f2;
    float t586 = a3*b4*c2*d0*e5*f1;
    float t587 = a3*b4*c2*d1*e0*f5;
    float t588 = a3*b4*c2*d5*e1*f0;
    float t589 = a3*b4*c5*d0*e1*f2;
    float t590 = a3*b4*c5*d1*e2*f0;
    float t591 = a3*b4*c5*d2*e0*f1;
    float t592 = a3*b5*c0*d1*e2*f4;
    float t593 = a3*b5*c0*d2*e4*f1;
    float t594 = a3*b5*c0*d4*e1*f2;
    float t595 = a3*b5*c1*d0*e4*f2;
    float t596 = a3*b5*c1*d2*e0*f4;
    float t597 = a3*b5*c1*d4*e2*f0;
    float t598 = a3*b5*c2*d0*e1*f4;
    float t599 = a3*b5*c2*d1*e4*f0;
    float t600 = a3*b5*c2*d4*e0*f1;
    float t601 = a3*b5*c4*d0*e2*f1;
    float t602 = a3*b5*c4*d1*e0*f2;
    float t603 = a3*b5*c4*d2*e1*f0;
    float t604 = a4*b0*c1*d2*e5*f3;
    float t605 = a4*b0*c1*d3*e2*f5;
    float t606 = a4*b0*c1*d5*e3*f2;
    float t607 = a4*b0*c2*d1*e3*f5;
    float t608 = a4*b0*c2*d3*e5*f1;
    float t609 = a4*b0*c2*d5*e1*f3;
    float t610 = a4*b0*c3*d1*e5*f2;
    float t611 = a4*b0*c3*d2*e1*f5;
    float t612 = a4*b0*c3*d5*e2*f1;
    float t613 = a4*b0*c5*d1*e2*f3;
    float t614 = a4*b0*c5*d2*e3*f1;
    float t615 = a4*b0*c5*d3*e1*f2;
    float t616 = a4*b1*c0*d2*e3*f5;
    float t617 = a4*b1*c0*d3*e5*f2;
    float t618 = a4*b1*c0*d5*e2*f3;
    float t619 = a4*b1*c2*d0*e5*f3;
    float t620 = a4*b1*c2*d3*e0*f5;
    float t621 = a4*b1*c2*d5*e3*f0;
    float t622 = a4*b1*c3*d0*e2*f5;
    float t623 = a4*b1*c3*d2*e5*f0;
    float t624 = a4*b1*c3*d5*e0*f2;
    float t625 = a4*b1*c5*d0*e3*f2;
    float t626 = a4*b1*c5*d2*e0*f3;
    float t627 = a4*b1*c5*d3*e2*f0;
    float t628 = a4*b2*c0*d1*e5*f3;
    float t629 = a4*b2*c0*d3*e1*f5;
    float t630 = a4*b2*c0*d5*e3*f1;
    float t631 = a4*b2*c1*d0*e3*f5;
    float t632 = a4*b2*c1*d3*e5*f0;
    float t633 = a4*b2*c1*d5*e0*f3;
    float t634 = a4*b2*c3*d0*e5*f1;
    float t635 = a4*b2*c3*d1*e0*f5;
    float t636 = a4*b2*c3*d5*e1*f0;
    float t637 = a4*b2*c5*d0*e1*f3;
    float t638 = a4*b2*c5*d1*e3*f0;
    float t639 = a4*b2*c5*d3*e0*f1;
    float t640 = a4*b3*c0*d1*e2*f5;
    float t641 = a4*b3*c0*d2*e5*f1;
    float t642 = a4*b3*c0*d5*e1*f2;
    float t643 = a4*b3*c1*d0*e5*f2;
    float t644 = a4*b3*c1*d2*e0*f5;
    float t645 = a4*b3*c1*d5*e2*f0;
    float t646 = a4*b3*c2*d0*e1*f5;
    float t647 = a4*b3*c2*d1*e5*f0;
    float t648 = a4*b3*c2*d5*e0*f1;
    float t649 = a4*b3*c5*d0*e2*f1;
    float t650 = a4*b3*c5*d1*e0*f2;
    float t651 = a4*b3*c5*d2*e1*f0;
    float t652 = a4*b5*c0*d1*e3*f2;
    float t653 = a4*b5*c0*d2*e1*f3;
    float t654 = a4*b5*c0*d3*e2*f1;
    float t655 = a4*b5*c1*d0*e2*f3;
    float t656 = a4*b5*c1*d2*e3*f0;
    float t657 = a4*b5*c1*d3*e0*f2;
    float t658 = a4*b5*c2*d0*e3*f1;
    float t659 = a4*b5*c2*d1*e0*f3;
    float t660 = a4*b5*c2*d3*e1*f0;
    float t661 = a4*b5*c3*d0*e1*f2;
    float t662 = a4*b5*c3*d1*e2*f0;
    float t663 = a4*b5*c3*d2*e0*f1;
    float t664 = a5*b0*c1*d2*e3*f4;
    float t665 = a5*b0*c1*d3*e4*f2;
    float t666 = a5*b0*c1*d4*e2*f3;
    float t667 = a5*b0*c2*d1*e4*f3;
    float t668 = a5*b0*c2*d3*e1*f4;
    float t669 = a5*b0*c2*d4*e3*f1;
    float t670 = a5*b0*c3*d1*e2*f4;
    float t671 = a5*b0*c3*d2*e4*f1;
    float t672 = a5*b0*c3*d4*e1*f2;
    float t673 = a5*b0*c4*d1*e3*f2;
    float t674 = a5*b0*c4*d2*e1*f3;
    float t675 = a5*b0*c4*d3*e2*f1;
    float t676 = a5*b1*c0*d2*e4*f3;
    float t677 = a5*b1*c0*d3*e2*f4;
    float t678 = a5*b1*c0*d4*e3*f2;
    float t679 = a5*b1*c2*d0*e3*f4;
    float t680 = a5*b1*c2*d3*e4*f0;
    float t681 = a5*b1*c2*d4*e0*f3;
    float t682 = a5*b1*c3*d0*e4*f2;
    float t683 = a5*b1*c3*d2*e0*f4;
    float t684 = a5*b1*c3*d4*e2*f0;
    float t685 = a5*b1*c4*d0*e2*f3;
    float t686 = a5*b1*c4*d2*e3*f0;
    float t687 = a5*b1*c4*d3*e0*f2;
    float t688 = a5*b2*c0*d1*e3*f4;
    float t689 = a5*b2*c0*d3*e4*f1;
    float t690 = a5*b2*c0*d4*e1*f3;
    float t691 = a5*b2*c1*d0*e4*f3;
    float t692 = a5*b2*c1*d3*e0*f4;
    float t693 = a5*b2*c1*d4*e3*f0;
    float t694 = a5*b2*c3*d0*e1*f4;
    float t695 = a5*b2*c3*d1*e4*f0;
    float t696 = a5*b2*c3*d4*e0*f1;
    float t697 = a5*b2*c4*d0*e3*f1;
    float t698 = a5*b2*c4*d1*e0*f3;
    float t699 = a5*b2*c4*d3*e1*f0;
    float t700 = a5*b3*c0*d1*e4*f2;
    float t701 = a5*b3*c0*d2*e1*f4;
    float t702 = a5*b3*c0*d4*e2*f1;
    float t703 = a5*b3*c1*d0*e2*f4;
    float t704 = a5*b3*c1*d2*e4*f0;
    float t705 = a5*b3*c1*d4*e0*f2;
    float t706 = a5*b3*c2*d0*e4*f1;
    float t707 = a5*b3*c2*d1*e0*f4;
    float t708 = a5*b3*c2*d4*e1*f0;
    float t709 = a5*b3*c4*d0*e1*f2;
    float t710 = a5*b3*c4*d1*e2*f0;
    float t711 = a5*b3*c4*d2*e0*f1;
    float t712 = a5*b4*c0*d1*e2*f3;
    float t713 = a5*b4*c0*d2*e3*f1;
    float t714 = a5*b4*c0*d3*e1*f2;
    float t715 = a5*b4*c1*d0*e3*f2;
    float t716 = a5*b4*c1*d2*e0*f3;
    float t717 = a5*b4*c1*d3*e2*f0;
    float t718 = a5*b4*c2*d0*e1*f3;
    float t719 = a5*b4*c2*d1*e3*f0;
    float t720 = a5*b4*c2*d3*e0*f1;
    float t721 = a5*b4*c3*d0*e2*f1;
    float t722 = a5*b4*c3*d1*e0*f2;
    float t723 = a5*b4*c3*d2*e1*f0;
    float t362 = t2+t3+t4+t5+t6+t7+t8+t9+t10+t11+t12+t13+t14+t15+t16+t17+t18+t19+t20+t21+t22+t23+t24+t25+t26+t27+t28
                 +t29+t30+t31+t32+t33+t34+t35+t36+t37+t38+t39+t40+t41+t42+t43+t44+t45+t46+t47+t48+t49+t50+t51+t52+t53
                 +t54+t55+t56+t57+t58+t59+t60+t61+t62+t63+t64+t65+t66+t67+t68+t69+t70+t71+t72+t73+t74+t75+t76+t77+t78
                 +t79+t80+t81+t82+t83+t84+t85+t86+t87+t88+t89+t90+t91+t92+t93+t94+t95+t96+t97+t98+t99+t100+t101+t102
                 +t103+t104+t105+t106+t107+t108+t109+t110+t111+t112+t113+t114+t115+t116+t117+t118+t119+t120+t121+t122
                 +t123+t124+t125+t126+t127+t128+t129+t130+t131+t132+t133+t134+t135+t136+t137+t138+t139+t140+t141+t142
                 +t143+t144+t145+t146+t147+t148+t149+t150+t151+t152+t153+t154+t155+t156+t157+t158+t159+t160+t161+t162
                 +t163+t164+t165+t166+t167+t168+t169+t170+t171+t172+t173+t174+t175+t176+t177+t178+t179+t180+t181+t182
                 +t183+t184+t185+t186+t187+t188+t189+t190+t191+t192+t193+t194+t195+t196+t197+t198+t199+t200+t201+t202
                 +t203+t204+t205+t206+t207+t208+t209+t210+t211+t212+t213+t214+t215+t216+t217+t218+t219+t220+t221+t222
                 +t223+t224+t225+t226+t227+t228+t229+t230+t231+t232+t233+t234+t235+t236+t237+t238+t239+t240+t241+t242
                 +t243+t244+t245+t246+t247+t248+t249+t250+t251+t252+t253+t254+t255+t256+t257+t258+t259+t260+t261+t262
                 +t263+t264+t265+t266+t267+t268+t269+t270+t271+t272+t273+t274+t275+t276+t277+t278+t279+t280+t281+t282
                 +t283+t284+t285+t286+t287+t288+t289+t290+t291+t292+t293+t294+t295+t296+t297+t298+t299+t300+t301+t302
                 +t303+t304+t305+t306+t307+t308+t309+t310+t311+t312+t313+t314+t315+t316+t317+t318+t319+t320+t321+t322
                 +t323+t324+t325+t326+t327+t328+t329+t330+t331+t332+t333+t334+t335+t336+t337+t338+t339+t340+t341+t342
                 +t343+t344+t345+t346+t347+t348+t349+t350+t351+t352+t353+t354+t355+t356+t357+t358+t359+t360+t361-t364
                 -t365-t366-t367-t368-t369-t370-t371-t372-t373-t374-t375-t376-t377-t378-t379-t380-t381-t382-t383-t384
                 -t385-t386-t387-t388-t389-t390-t391-t392-t393-t394-t395-t396-t397-t398-t399-t400-t401-t402-t403-t404
                 -t405-t406-t407-t408-t409-t410-t411-t412-t413-t414-t415-t416-t417-t418-t419-t420-t421-t422-t423-t424
                 -t425-t426-t427-t428-t429-t430-t431-t432-t433-t434-t435-t436-t437-t438-t439-t440-t441-t442-t443-t444
                 -t445-t446-t447-t448-t449-t450-t451-t452-t453-t454-t455-t456-t457-t458-t459-t460-t461-t462-t463-t464
                 -t465-t466-t467-t468-t469-t470-t471-t472-t473-t474-t475-t476-t477-t478-t479-t480-t481-t482-t483-t484
                 -t485-t486-t487-t488-t489-t490-t491-t492-t493-t494-t495-t496-t497-t498-t499-t500-t501-t502-t503-t504
                 -t505-t506-t507-t508-t509-t510-t511-t512-t513-t514-t515-t516-t517-t518-t519-t520-t521-t522-t523-t524
                 -t525-t526-t527-t528-t529-t530-t531-t532-t533-t534-t535-t536-t537-t538-t539-t540-t541-t542-t543-t544
                 -t545-t546-t547-t548-t549-t550-t551-t552-t553-t554-t555-t556-t557-t558-t559-t560-t561-t562-t563-t564
                 -t565-t566-t567-t568-t569-t570-t571-t572-t573-t574-t575-t576-t577-t578-t579-t580-t581-t582-t583-t584
                 -t585-t586-t587-t588-t589-t590-t591-t592-t593-t594-t595-t596-t597-t598-t599-t600-t601-t602-t603-t604
                 -t605-t606-t607-t608-t609-t610-t611-t612-t613-t614-t615-t616-t617-t618-t619-t620-t621-t622-t623-t624
                 -t625-t626-t627-t628-t629-t630-t631-t632-t633-t634-t635-t636-t637-t638-t639-t640-t641-t642-t643-t644
                 -t645-t646-t647-t648-t649-t650-t651-t652-t653-t654-t655-t656-t657-t658-t659-t660-t661-t662-t663-t664
                 -t665-t666-t667-t668-t669-t670-t671-t672-t673-t674-t675-t676-t677-t678-t679-t680-t681-t682-t683-t684
                 -t685-t686-t687-t688-t689-t690-t691-t692-t693-t694-t695-t696-t697-t698-t699-t700-t701-t702-t703-t704
                 -t705-t706-t707-t708-t709-t710-t711-t712-t713-t714-t715-t716-t717-t718-t719-t720-t721-t722-t723;
    float t363 = 1.0f/t362;
    if (t363) {
        T[0][0] = t363 *
                  (b1 * c2 * d3 * e4 * f5 - b1 * c2 * d3 * e5 * f4 - b1 * c2 * d4 * e3 * f5 + b1 * c2 * d4 * e5 * f3 +
                   b1 * c2 * d5 * e3 * f4 - b1 * c2 * d5 * e4 * f3 - b1 * c3 * d2 * e4 * f5 + b1 * c3 * d2 * e5 * f4 +
                   b1 * c3 * d4 * e2 * f5 - b1 * c3 * d4 * e5 * f2 - b1 * c3 * d5 * e2 * f4 + b1 * c3 * d5 * e4 * f2 +
                   b1 * c4 * d2 * e3 * f5 - b1 * c4 * d2 * e5 * f3 - b1 * c4 * d3 * e2 * f5 + b1 * c4 * d3 * e5 * f2 +
                   b1 * c4 * d5 * e2 * f3 - b1 * c4 * d5 * e3 * f2 - b1 * c5 * d2 * e3 * f4 + b1 * c5 * d2 * e4 * f3 +
                   b1 * c5 * d3 * e2 * f4 - b1 * c5 * d3 * e4 * f2 - b1 * c5 * d4 * e2 * f3 + b1 * c5 * d4 * e3 * f2 -
                   b2 * c1 * d3 * e4 * f5 + b2 * c1 * d3 * e5 * f4 + b2 * c1 * d4 * e3 * f5 - b2 * c1 * d4 * e5 * f3 -
                   b2 * c1 * d5 * e3 * f4 + b2 * c1 * d5 * e4 * f3 + b2 * c3 * d1 * e4 * f5 - b2 * c3 * d1 * e5 * f4 -
                   b2 * c3 * d4 * e1 * f5 + b2 * c3 * d4 * e5 * f1 + b2 * c3 * d5 * e1 * f4 - b2 * c3 * d5 * e4 * f1 -
                   b2 * c4 * d1 * e3 * f5 + b2 * c4 * d1 * e5 * f3 + b2 * c4 * d3 * e1 * f5 - b2 * c4 * d3 * e5 * f1 -
                   b2 * c4 * d5 * e1 * f3 + b2 * c4 * d5 * e3 * f1 + b2 * c5 * d1 * e3 * f4 - b2 * c5 * d1 * e4 * f3 -
                   b2 * c5 * d3 * e1 * f4 + b2 * c5 * d3 * e4 * f1 + b2 * c5 * d4 * e1 * f3 - b2 * c5 * d4 * e3 * f1 +
                   b3 * c1 * d2 * e4 * f5 - b3 * c1 * d2 * e5 * f4 - b3 * c1 * d4 * e2 * f5 + b3 * c1 * d4 * e5 * f2 +
                   b3 * c1 * d5 * e2 * f4 - b3 * c1 * d5 * e4 * f2 - b3 * c2 * d1 * e4 * f5 + b3 * c2 * d1 * e5 * f4 +
                   b3 * c2 * d4 * e1 * f5 - b3 * c2 * d4 * e5 * f1 - b3 * c2 * d5 * e1 * f4 + b3 * c2 * d5 * e4 * f1 +
                   b3 * c4 * d1 * e2 * f5 - b3 * c4 * d1 * e5 * f2 - b3 * c4 * d2 * e1 * f5 + b3 * c4 * d2 * e5 * f1 +
                   b3 * c4 * d5 * e1 * f2 - b3 * c4 * d5 * e2 * f1 - b3 * c5 * d1 * e2 * f4 + b3 * c5 * d1 * e4 * f2 +
                   b3 * c5 * d2 * e1 * f4 - b3 * c5 * d2 * e4 * f1 - b3 * c5 * d4 * e1 * f2 + b3 * c5 * d4 * e2 * f1 -
                   b4 * c1 * d2 * e3 * f5 + b4 * c1 * d2 * e5 * f3 + b4 * c1 * d3 * e2 * f5 - b4 * c1 * d3 * e5 * f2 -
                   b4 * c1 * d5 * e2 * f3 + b4 * c1 * d5 * e3 * f2 + b4 * c2 * d1 * e3 * f5 - b4 * c2 * d1 * e5 * f3 -
                   b4 * c2 * d3 * e1 * f5 + b4 * c2 * d3 * e5 * f1 + b4 * c2 * d5 * e1 * f3 - b4 * c2 * d5 * e3 * f1 -
                   b4 * c3 * d1 * e2 * f5 + b4 * c3 * d1 * e5 * f2 + b4 * c3 * d2 * e1 * f5 - b4 * c3 * d2 * e5 * f1 -
                   b4 * c3 * d5 * e1 * f2 + b4 * c3 * d5 * e2 * f1 + b4 * c5 * d1 * e2 * f3 - b4 * c5 * d1 * e3 * f2 -
                   b4 * c5 * d2 * e1 * f3 + b4 * c5 * d2 * e3 * f1 + b4 * c5 * d3 * e1 * f2 - b4 * c5 * d3 * e2 * f1 +
                   b5 * c1 * d2 * e3 * f4 - b5 * c1 * d2 * e4 * f3 - b5 * c1 * d3 * e2 * f4 + b5 * c1 * d3 * e4 * f2 +
                   b5 * c1 * d4 * e2 * f3 - b5 * c1 * d4 * e3 * f2 - b5 * c2 * d1 * e3 * f4 + b5 * c2 * d1 * e4 * f3 +
                   b5 * c2 * d3 * e1 * f4 - b5 * c2 * d3 * e4 * f1 - b5 * c2 * d4 * e1 * f3 + b5 * c2 * d4 * e3 * f1 +
                   b5 * c3 * d1 * e2 * f4 - b5 * c3 * d1 * e4 * f2 - b5 * c3 * d2 * e1 * f4 + b5 * c3 * d2 * e4 * f1 +
                   b5 * c3 * d4 * e1 * f2 - b5 * c3 * d4 * e2 * f1 - b5 * c4 * d1 * e2 * f3 + b5 * c4 * d1 * e3 * f2 +
                   b5 * c4 * d2 * e1 * f3 - b5 * c4 * d2 * e3 * f1 - b5 * c4 * d3 * e1 * f2 + b5 * c4 * d3 * e2 * f1);
        T[1][0] = -t363 *
                  (a1 * c2 * d3 * e4 * f5 - a1 * c2 * d3 * e5 * f4 - a1 * c2 * d4 * e3 * f5 + a1 * c2 * d4 * e5 * f3 +
                   a1 * c2 * d5 * e3 * f4 - a1 * c2 * d5 * e4 * f3 - a1 * c3 * d2 * e4 * f5 + a1 * c3 * d2 * e5 * f4 +
                   a1 * c3 * d4 * e2 * f5 - a1 * c3 * d4 * e5 * f2 - a1 * c3 * d5 * e2 * f4 + a1 * c3 * d5 * e4 * f2 +
                   a1 * c4 * d2 * e3 * f5 - a1 * c4 * d2 * e5 * f3 - a1 * c4 * d3 * e2 * f5 + a1 * c4 * d3 * e5 * f2 +
                   a1 * c4 * d5 * e2 * f3 - a1 * c4 * d5 * e3 * f2 - a1 * c5 * d2 * e3 * f4 + a1 * c5 * d2 * e4 * f3 +
                   a1 * c5 * d3 * e2 * f4 - a1 * c5 * d3 * e4 * f2 - a1 * c5 * d4 * e2 * f3 + a1 * c5 * d4 * e3 * f2 -
                   a2 * c1 * d3 * e4 * f5 + a2 * c1 * d3 * e5 * f4 + a2 * c1 * d4 * e3 * f5 - a2 * c1 * d4 * e5 * f3 -
                   a2 * c1 * d5 * e3 * f4 + a2 * c1 * d5 * e4 * f3 + a2 * c3 * d1 * e4 * f5 - a2 * c3 * d1 * e5 * f4 -
                   a2 * c3 * d4 * e1 * f5 + a2 * c3 * d4 * e5 * f1 + a2 * c3 * d5 * e1 * f4 - a2 * c3 * d5 * e4 * f1 -
                   a2 * c4 * d1 * e3 * f5 + a2 * c4 * d1 * e5 * f3 + a2 * c4 * d3 * e1 * f5 - a2 * c4 * d3 * e5 * f1 -
                   a2 * c4 * d5 * e1 * f3 + a2 * c4 * d5 * e3 * f1 + a2 * c5 * d1 * e3 * f4 - a2 * c5 * d1 * e4 * f3 -
                   a2 * c5 * d3 * e1 * f4 + a2 * c5 * d3 * e4 * f1 + a2 * c5 * d4 * e1 * f3 - a2 * c5 * d4 * e3 * f1 +
                   a3 * c1 * d2 * e4 * f5 - a3 * c1 * d2 * e5 * f4 - a3 * c1 * d4 * e2 * f5 + a3 * c1 * d4 * e5 * f2 +
                   a3 * c1 * d5 * e2 * f4 - a3 * c1 * d5 * e4 * f2 - a3 * c2 * d1 * e4 * f5 + a3 * c2 * d1 * e5 * f4 +
                   a3 * c2 * d4 * e1 * f5 - a3 * c2 * d4 * e5 * f1 - a3 * c2 * d5 * e1 * f4 + a3 * c2 * d5 * e4 * f1 +
                   a3 * c4 * d1 * e2 * f5 - a3 * c4 * d1 * e5 * f2 - a3 * c4 * d2 * e1 * f5 + a3 * c4 * d2 * e5 * f1 +
                   a3 * c4 * d5 * e1 * f2 - a3 * c4 * d5 * e2 * f1 - a3 * c5 * d1 * e2 * f4 + a3 * c5 * d1 * e4 * f2 +
                   a3 * c5 * d2 * e1 * f4 - a3 * c5 * d2 * e4 * f1 - a3 * c5 * d4 * e1 * f2 + a3 * c5 * d4 * e2 * f1 -
                   a4 * c1 * d2 * e3 * f5 + a4 * c1 * d2 * e5 * f3 + a4 * c1 * d3 * e2 * f5 - a4 * c1 * d3 * e5 * f2 -
                   a4 * c1 * d5 * e2 * f3 + a4 * c1 * d5 * e3 * f2 + a4 * c2 * d1 * e3 * f5 - a4 * c2 * d1 * e5 * f3 -
                   a4 * c2 * d3 * e1 * f5 + a4 * c2 * d3 * e5 * f1 + a4 * c2 * d5 * e1 * f3 - a4 * c2 * d5 * e3 * f1 -
                   a4 * c3 * d1 * e2 * f5 + a4 * c3 * d1 * e5 * f2 + a4 * c3 * d2 * e1 * f5 - a4 * c3 * d2 * e5 * f1 -
                   a4 * c3 * d5 * e1 * f2 + a4 * c3 * d5 * e2 * f1 + a4 * c5 * d1 * e2 * f3 - a4 * c5 * d1 * e3 * f2 -
                   a4 * c5 * d2 * e1 * f3 + a4 * c5 * d2 * e3 * f1 + a4 * c5 * d3 * e1 * f2 - a4 * c5 * d3 * e2 * f1 +
                   a5 * c1 * d2 * e3 * f4 - a5 * c1 * d2 * e4 * f3 - a5 * c1 * d3 * e2 * f4 + a5 * c1 * d3 * e4 * f2 +
                   a5 * c1 * d4 * e2 * f3 - a5 * c1 * d4 * e3 * f2 - a5 * c2 * d1 * e3 * f4 + a5 * c2 * d1 * e4 * f3 +
                   a5 * c2 * d3 * e1 * f4 - a5 * c2 * d3 * e4 * f1 - a5 * c2 * d4 * e1 * f3 + a5 * c2 * d4 * e3 * f1 +
                   a5 * c3 * d1 * e2 * f4 - a5 * c3 * d1 * e4 * f2 - a5 * c3 * d2 * e1 * f4 + a5 * c3 * d2 * e4 * f1 +
                   a5 * c3 * d4 * e1 * f2 - a5 * c3 * d4 * e2 * f1 - a5 * c4 * d1 * e2 * f3 + a5 * c4 * d1 * e3 * f2 +
                   a5 * c4 * d2 * e1 * f3 - a5 * c4 * d2 * e3 * f1 - a5 * c4 * d3 * e1 * f2 + a5 * c4 * d3 * e2 * f1);
        T[2][0] = t363 *
                  (a1 * b2 * d3 * e4 * f5 - a1 * b2 * d3 * e5 * f4 - a1 * b2 * d4 * e3 * f5 + a1 * b2 * d4 * e5 * f3 +
                   a1 * b2 * d5 * e3 * f4 - a1 * b2 * d5 * e4 * f3 - a1 * b3 * d2 * e4 * f5 + a1 * b3 * d2 * e5 * f4 +
                   a1 * b3 * d4 * e2 * f5 - a1 * b3 * d4 * e5 * f2 - a1 * b3 * d5 * e2 * f4 + a1 * b3 * d5 * e4 * f2 +
                   a1 * b4 * d2 * e3 * f5 - a1 * b4 * d2 * e5 * f3 - a1 * b4 * d3 * e2 * f5 + a1 * b4 * d3 * e5 * f2 +
                   a1 * b4 * d5 * e2 * f3 - a1 * b4 * d5 * e3 * f2 - a1 * b5 * d2 * e3 * f4 + a1 * b5 * d2 * e4 * f3 +
                   a1 * b5 * d3 * e2 * f4 - a1 * b5 * d3 * e4 * f2 - a1 * b5 * d4 * e2 * f3 + a1 * b5 * d4 * e3 * f2 -
                   a2 * b1 * d3 * e4 * f5 + a2 * b1 * d3 * e5 * f4 + a2 * b1 * d4 * e3 * f5 - a2 * b1 * d4 * e5 * f3 -
                   a2 * b1 * d5 * e3 * f4 + a2 * b1 * d5 * e4 * f3 + a2 * b3 * d1 * e4 * f5 - a2 * b3 * d1 * e5 * f4 -
                   a2 * b3 * d4 * e1 * f5 + a2 * b3 * d4 * e5 * f1 + a2 * b3 * d5 * e1 * f4 - a2 * b3 * d5 * e4 * f1 -
                   a2 * b4 * d1 * e3 * f5 + a2 * b4 * d1 * e5 * f3 + a2 * b4 * d3 * e1 * f5 - a2 * b4 * d3 * e5 * f1 -
                   a2 * b4 * d5 * e1 * f3 + a2 * b4 * d5 * e3 * f1 + a2 * b5 * d1 * e3 * f4 - a2 * b5 * d1 * e4 * f3 -
                   a2 * b5 * d3 * e1 * f4 + a2 * b5 * d3 * e4 * f1 + a2 * b5 * d4 * e1 * f3 - a2 * b5 * d4 * e3 * f1 +
                   a3 * b1 * d2 * e4 * f5 - a3 * b1 * d2 * e5 * f4 - a3 * b1 * d4 * e2 * f5 + a3 * b1 * d4 * e5 * f2 +
                   a3 * b1 * d5 * e2 * f4 - a3 * b1 * d5 * e4 * f2 - a3 * b2 * d1 * e4 * f5 + a3 * b2 * d1 * e5 * f4 +
                   a3 * b2 * d4 * e1 * f5 - a3 * b2 * d4 * e5 * f1 - a3 * b2 * d5 * e1 * f4 + a3 * b2 * d5 * e4 * f1 +
                   a3 * b4 * d1 * e2 * f5 - a3 * b4 * d1 * e5 * f2 - a3 * b4 * d2 * e1 * f5 + a3 * b4 * d2 * e5 * f1 +
                   a3 * b4 * d5 * e1 * f2 - a3 * b4 * d5 * e2 * f1 - a3 * b5 * d1 * e2 * f4 + a3 * b5 * d1 * e4 * f2 +
                   a3 * b5 * d2 * e1 * f4 - a3 * b5 * d2 * e4 * f1 - a3 * b5 * d4 * e1 * f2 + a3 * b5 * d4 * e2 * f1 -
                   a4 * b1 * d2 * e3 * f5 + a4 * b1 * d2 * e5 * f3 + a4 * b1 * d3 * e2 * f5 - a4 * b1 * d3 * e5 * f2 -
                   a4 * b1 * d5 * e2 * f3 + a4 * b1 * d5 * e3 * f2 + a4 * b2 * d1 * e3 * f5 - a4 * b2 * d1 * e5 * f3 -
                   a4 * b2 * d3 * e1 * f5 + a4 * b2 * d3 * e5 * f1 + a4 * b2 * d5 * e1 * f3 - a4 * b2 * d5 * e3 * f1 -
                   a4 * b3 * d1 * e2 * f5 + a4 * b3 * d1 * e5 * f2 + a4 * b3 * d2 * e1 * f5 - a4 * b3 * d2 * e5 * f1 -
                   a4 * b3 * d5 * e1 * f2 + a4 * b3 * d5 * e2 * f1 + a4 * b5 * d1 * e2 * f3 - a4 * b5 * d1 * e3 * f2 -
                   a4 * b5 * d2 * e1 * f3 + a4 * b5 * d2 * e3 * f1 + a4 * b5 * d3 * e1 * f2 - a4 * b5 * d3 * e2 * f1 +
                   a5 * b1 * d2 * e3 * f4 - a5 * b1 * d2 * e4 * f3 - a5 * b1 * d3 * e2 * f4 + a5 * b1 * d3 * e4 * f2 +
                   a5 * b1 * d4 * e2 * f3 - a5 * b1 * d4 * e3 * f2 - a5 * b2 * d1 * e3 * f4 + a5 * b2 * d1 * e4 * f3 +
                   a5 * b2 * d3 * e1 * f4 - a5 * b2 * d3 * e4 * f1 - a5 * b2 * d4 * e1 * f3 + a5 * b2 * d4 * e3 * f1 +
                   a5 * b3 * d1 * e2 * f4 - a5 * b3 * d1 * e4 * f2 - a5 * b3 * d2 * e1 * f4 + a5 * b3 * d2 * e4 * f1 +
                   a5 * b3 * d4 * e1 * f2 - a5 * b3 * d4 * e2 * f1 - a5 * b4 * d1 * e2 * f3 + a5 * b4 * d1 * e3 * f2 +
                   a5 * b4 * d2 * e1 * f3 - a5 * b4 * d2 * e3 * f1 - a5 * b4 * d3 * e1 * f2 + a5 * b4 * d3 * e2 * f1);
        T[3][0] = -t363 *
                  (a1 * b2 * c3 * e4 * f5 - a1 * b2 * c3 * e5 * f4 - a1 * b2 * c4 * e3 * f5 + a1 * b2 * c4 * e5 * f3 +
                   a1 * b2 * c5 * e3 * f4 - a1 * b2 * c5 * e4 * f3 - a1 * b3 * c2 * e4 * f5 + a1 * b3 * c2 * e5 * f4 +
                   a1 * b3 * c4 * e2 * f5 - a1 * b3 * c4 * e5 * f2 - a1 * b3 * c5 * e2 * f4 + a1 * b3 * c5 * e4 * f2 +
                   a1 * b4 * c2 * e3 * f5 - a1 * b4 * c2 * e5 * f3 - a1 * b4 * c3 * e2 * f5 + a1 * b4 * c3 * e5 * f2 +
                   a1 * b4 * c5 * e2 * f3 - a1 * b4 * c5 * e3 * f2 - a1 * b5 * c2 * e3 * f4 + a1 * b5 * c2 * e4 * f3 +
                   a1 * b5 * c3 * e2 * f4 - a1 * b5 * c3 * e4 * f2 - a1 * b5 * c4 * e2 * f3 + a1 * b5 * c4 * e3 * f2 -
                   a2 * b1 * c3 * e4 * f5 + a2 * b1 * c3 * e5 * f4 + a2 * b1 * c4 * e3 * f5 - a2 * b1 * c4 * e5 * f3 -
                   a2 * b1 * c5 * e3 * f4 + a2 * b1 * c5 * e4 * f3 + a2 * b3 * c1 * e4 * f5 - a2 * b3 * c1 * e5 * f4 -
                   a2 * b3 * c4 * e1 * f5 + a2 * b3 * c4 * e5 * f1 + a2 * b3 * c5 * e1 * f4 - a2 * b3 * c5 * e4 * f1 -
                   a2 * b4 * c1 * e3 * f5 + a2 * b4 * c1 * e5 * f3 + a2 * b4 * c3 * e1 * f5 - a2 * b4 * c3 * e5 * f1 -
                   a2 * b4 * c5 * e1 * f3 + a2 * b4 * c5 * e3 * f1 + a2 * b5 * c1 * e3 * f4 - a2 * b5 * c1 * e4 * f3 -
                   a2 * b5 * c3 * e1 * f4 + a2 * b5 * c3 * e4 * f1 + a2 * b5 * c4 * e1 * f3 - a2 * b5 * c4 * e3 * f1 +
                   a3 * b1 * c2 * e4 * f5 - a3 * b1 * c2 * e5 * f4 - a3 * b1 * c4 * e2 * f5 + a3 * b1 * c4 * e5 * f2 +
                   a3 * b1 * c5 * e2 * f4 - a3 * b1 * c5 * e4 * f2 - a3 * b2 * c1 * e4 * f5 + a3 * b2 * c1 * e5 * f4 +
                   a3 * b2 * c4 * e1 * f5 - a3 * b2 * c4 * e5 * f1 - a3 * b2 * c5 * e1 * f4 + a3 * b2 * c5 * e4 * f1 +
                   a3 * b4 * c1 * e2 * f5 - a3 * b4 * c1 * e5 * f2 - a3 * b4 * c2 * e1 * f5 + a3 * b4 * c2 * e5 * f1 +
                   a3 * b4 * c5 * e1 * f2 - a3 * b4 * c5 * e2 * f1 - a3 * b5 * c1 * e2 * f4 + a3 * b5 * c1 * e4 * f2 +
                   a3 * b5 * c2 * e1 * f4 - a3 * b5 * c2 * e4 * f1 - a3 * b5 * c4 * e1 * f2 + a3 * b5 * c4 * e2 * f1 -
                   a4 * b1 * c2 * e3 * f5 + a4 * b1 * c2 * e5 * f3 + a4 * b1 * c3 * e2 * f5 - a4 * b1 * c3 * e5 * f2 -
                   a4 * b1 * c5 * e2 * f3 + a4 * b1 * c5 * e3 * f2 + a4 * b2 * c1 * e3 * f5 - a4 * b2 * c1 * e5 * f3 -
                   a4 * b2 * c3 * e1 * f5 + a4 * b2 * c3 * e5 * f1 + a4 * b2 * c5 * e1 * f3 - a4 * b2 * c5 * e3 * f1 -
                   a4 * b3 * c1 * e2 * f5 + a4 * b3 * c1 * e5 * f2 + a4 * b3 * c2 * e1 * f5 - a4 * b3 * c2 * e5 * f1 -
                   a4 * b3 * c5 * e1 * f2 + a4 * b3 * c5 * e2 * f1 + a4 * b5 * c1 * e2 * f3 - a4 * b5 * c1 * e3 * f2 -
                   a4 * b5 * c2 * e1 * f3 + a4 * b5 * c2 * e3 * f1 + a4 * b5 * c3 * e1 * f2 - a4 * b5 * c3 * e2 * f1 +
                   a5 * b1 * c2 * e3 * f4 - a5 * b1 * c2 * e4 * f3 - a5 * b1 * c3 * e2 * f4 + a5 * b1 * c3 * e4 * f2 +
                   a5 * b1 * c4 * e2 * f3 - a5 * b1 * c4 * e3 * f2 - a5 * b2 * c1 * e3 * f4 + a5 * b2 * c1 * e4 * f3 +
                   a5 * b2 * c3 * e1 * f4 - a5 * b2 * c3 * e4 * f1 - a5 * b2 * c4 * e1 * f3 + a5 * b2 * c4 * e3 * f1 +
                   a5 * b3 * c1 * e2 * f4 - a5 * b3 * c1 * e4 * f2 - a5 * b3 * c2 * e1 * f4 + a5 * b3 * c2 * e4 * f1 +
                   a5 * b3 * c4 * e1 * f2 - a5 * b3 * c4 * e2 * f1 - a5 * b4 * c1 * e2 * f3 + a5 * b4 * c1 * e3 * f2 +
                   a5 * b4 * c2 * e1 * f3 - a5 * b4 * c2 * e3 * f1 - a5 * b4 * c3 * e1 * f2 + a5 * b4 * c3 * e2 * f1);
        T[4][0] = t363 *
                  (a1 * b2 * c3 * d4 * f5 - a1 * b2 * c3 * d5 * f4 - a1 * b2 * c4 * d3 * f5 + a1 * b2 * c4 * d5 * f3 +
                   a1 * b2 * c5 * d3 * f4 - a1 * b2 * c5 * d4 * f3 - a1 * b3 * c2 * d4 * f5 + a1 * b3 * c2 * d5 * f4 +
                   a1 * b3 * c4 * d2 * f5 - a1 * b3 * c4 * d5 * f2 - a1 * b3 * c5 * d2 * f4 + a1 * b3 * c5 * d4 * f2 +
                   a1 * b4 * c2 * d3 * f5 - a1 * b4 * c2 * d5 * f3 - a1 * b4 * c3 * d2 * f5 + a1 * b4 * c3 * d5 * f2 +
                   a1 * b4 * c5 * d2 * f3 - a1 * b4 * c5 * d3 * f2 - a1 * b5 * c2 * d3 * f4 + a1 * b5 * c2 * d4 * f3 +
                   a1 * b5 * c3 * d2 * f4 - a1 * b5 * c3 * d4 * f2 - a1 * b5 * c4 * d2 * f3 + a1 * b5 * c4 * d3 * f2 -
                   a2 * b1 * c3 * d4 * f5 + a2 * b1 * c3 * d5 * f4 + a2 * b1 * c4 * d3 * f5 - a2 * b1 * c4 * d5 * f3 -
                   a2 * b1 * c5 * d3 * f4 + a2 * b1 * c5 * d4 * f3 + a2 * b3 * c1 * d4 * f5 - a2 * b3 * c1 * d5 * f4 -
                   a2 * b3 * c4 * d1 * f5 + a2 * b3 * c4 * d5 * f1 + a2 * b3 * c5 * d1 * f4 - a2 * b3 * c5 * d4 * f1 -
                   a2 * b4 * c1 * d3 * f5 + a2 * b4 * c1 * d5 * f3 + a2 * b4 * c3 * d1 * f5 - a2 * b4 * c3 * d5 * f1 -
                   a2 * b4 * c5 * d1 * f3 + a2 * b4 * c5 * d3 * f1 + a2 * b5 * c1 * d3 * f4 - a2 * b5 * c1 * d4 * f3 -
                   a2 * b5 * c3 * d1 * f4 + a2 * b5 * c3 * d4 * f1 + a2 * b5 * c4 * d1 * f3 - a2 * b5 * c4 * d3 * f1 +
                   a3 * b1 * c2 * d4 * f5 - a3 * b1 * c2 * d5 * f4 - a3 * b1 * c4 * d2 * f5 + a3 * b1 * c4 * d5 * f2 +
                   a3 * b1 * c5 * d2 * f4 - a3 * b1 * c5 * d4 * f2 - a3 * b2 * c1 * d4 * f5 + a3 * b2 * c1 * d5 * f4 +
                   a3 * b2 * c4 * d1 * f5 - a3 * b2 * c4 * d5 * f1 - a3 * b2 * c5 * d1 * f4 + a3 * b2 * c5 * d4 * f1 +
                   a3 * b4 * c1 * d2 * f5 - a3 * b4 * c1 * d5 * f2 - a3 * b4 * c2 * d1 * f5 + a3 * b4 * c2 * d5 * f1 +
                   a3 * b4 * c5 * d1 * f2 - a3 * b4 * c5 * d2 * f1 - a3 * b5 * c1 * d2 * f4 + a3 * b5 * c1 * d4 * f2 +
                   a3 * b5 * c2 * d1 * f4 - a3 * b5 * c2 * d4 * f1 - a3 * b5 * c4 * d1 * f2 + a3 * b5 * c4 * d2 * f1 -
                   a4 * b1 * c2 * d3 * f5 + a4 * b1 * c2 * d5 * f3 + a4 * b1 * c3 * d2 * f5 - a4 * b1 * c3 * d5 * f2 -
                   a4 * b1 * c5 * d2 * f3 + a4 * b1 * c5 * d3 * f2 + a4 * b2 * c1 * d3 * f5 - a4 * b2 * c1 * d5 * f3 -
                   a4 * b2 * c3 * d1 * f5 + a4 * b2 * c3 * d5 * f1 + a4 * b2 * c5 * d1 * f3 - a4 * b2 * c5 * d3 * f1 -
                   a4 * b3 * c1 * d2 * f5 + a4 * b3 * c1 * d5 * f2 + a4 * b3 * c2 * d1 * f5 - a4 * b3 * c2 * d5 * f1 -
                   a4 * b3 * c5 * d1 * f2 + a4 * b3 * c5 * d2 * f1 + a4 * b5 * c1 * d2 * f3 - a4 * b5 * c1 * d3 * f2 -
                   a4 * b5 * c2 * d1 * f3 + a4 * b5 * c2 * d3 * f1 + a4 * b5 * c3 * d1 * f2 - a4 * b5 * c3 * d2 * f1 +
                   a5 * b1 * c2 * d3 * f4 - a5 * b1 * c2 * d4 * f3 - a5 * b1 * c3 * d2 * f4 + a5 * b1 * c3 * d4 * f2 +
                   a5 * b1 * c4 * d2 * f3 - a5 * b1 * c4 * d3 * f2 - a5 * b2 * c1 * d3 * f4 + a5 * b2 * c1 * d4 * f3 +
                   a5 * b2 * c3 * d1 * f4 - a5 * b2 * c3 * d4 * f1 - a5 * b2 * c4 * d1 * f3 + a5 * b2 * c4 * d3 * f1 +
                   a5 * b3 * c1 * d2 * f4 - a5 * b3 * c1 * d4 * f2 - a5 * b3 * c2 * d1 * f4 + a5 * b3 * c2 * d4 * f1 +
                   a5 * b3 * c4 * d1 * f2 - a5 * b3 * c4 * d2 * f1 - a5 * b4 * c1 * d2 * f3 + a5 * b4 * c1 * d3 * f2 +
                   a5 * b4 * c2 * d1 * f3 - a5 * b4 * c2 * d3 * f1 - a5 * b4 * c3 * d1 * f2 + a5 * b4 * c3 * d2 * f1);
        T[5][0] = -t363 *
                  (a1 * b2 * c3 * d4 * e5 - a1 * b2 * c3 * d5 * e4 - a1 * b2 * c4 * d3 * e5 + a1 * b2 * c4 * d5 * e3 +
                   a1 * b2 * c5 * d3 * e4 - a1 * b2 * c5 * d4 * e3 - a1 * b3 * c2 * d4 * e5 + a1 * b3 * c2 * d5 * e4 +
                   a1 * b3 * c4 * d2 * e5 - a1 * b3 * c4 * d5 * e2 - a1 * b3 * c5 * d2 * e4 + a1 * b3 * c5 * d4 * e2 +
                   a1 * b4 * c2 * d3 * e5 - a1 * b4 * c2 * d5 * e3 - a1 * b4 * c3 * d2 * e5 + a1 * b4 * c3 * d5 * e2 +
                   a1 * b4 * c5 * d2 * e3 - a1 * b4 * c5 * d3 * e2 - a1 * b5 * c2 * d3 * e4 + a1 * b5 * c2 * d4 * e3 +
                   a1 * b5 * c3 * d2 * e4 - a1 * b5 * c3 * d4 * e2 - a1 * b5 * c4 * d2 * e3 + a1 * b5 * c4 * d3 * e2 -
                   a2 * b1 * c3 * d4 * e5 + a2 * b1 * c3 * d5 * e4 + a2 * b1 * c4 * d3 * e5 - a2 * b1 * c4 * d5 * e3 -
                   a2 * b1 * c5 * d3 * e4 + a2 * b1 * c5 * d4 * e3 + a2 * b3 * c1 * d4 * e5 - a2 * b3 * c1 * d5 * e4 -
                   a2 * b3 * c4 * d1 * e5 + a2 * b3 * c4 * d5 * e1 + a2 * b3 * c5 * d1 * e4 - a2 * b3 * c5 * d4 * e1 -
                   a2 * b4 * c1 * d3 * e5 + a2 * b4 * c1 * d5 * e3 + a2 * b4 * c3 * d1 * e5 - a2 * b4 * c3 * d5 * e1 -
                   a2 * b4 * c5 * d1 * e3 + a2 * b4 * c5 * d3 * e1 + a2 * b5 * c1 * d3 * e4 - a2 * b5 * c1 * d4 * e3 -
                   a2 * b5 * c3 * d1 * e4 + a2 * b5 * c3 * d4 * e1 + a2 * b5 * c4 * d1 * e3 - a2 * b5 * c4 * d3 * e1 +
                   a3 * b1 * c2 * d4 * e5 - a3 * b1 * c2 * d5 * e4 - a3 * b1 * c4 * d2 * e5 + a3 * b1 * c4 * d5 * e2 +
                   a3 * b1 * c5 * d2 * e4 - a3 * b1 * c5 * d4 * e2 - a3 * b2 * c1 * d4 * e5 + a3 * b2 * c1 * d5 * e4 +
                   a3 * b2 * c4 * d1 * e5 - a3 * b2 * c4 * d5 * e1 - a3 * b2 * c5 * d1 * e4 + a3 * b2 * c5 * d4 * e1 +
                   a3 * b4 * c1 * d2 * e5 - a3 * b4 * c1 * d5 * e2 - a3 * b4 * c2 * d1 * e5 + a3 * b4 * c2 * d5 * e1 +
                   a3 * b4 * c5 * d1 * e2 - a3 * b4 * c5 * d2 * e1 - a3 * b5 * c1 * d2 * e4 + a3 * b5 * c1 * d4 * e2 +
                   a3 * b5 * c2 * d1 * e4 - a3 * b5 * c2 * d4 * e1 - a3 * b5 * c4 * d1 * e2 + a3 * b5 * c4 * d2 * e1 -
                   a4 * b1 * c2 * d3 * e5 + a4 * b1 * c2 * d5 * e3 + a4 * b1 * c3 * d2 * e5 - a4 * b1 * c3 * d5 * e2 -
                   a4 * b1 * c5 * d2 * e3 + a4 * b1 * c5 * d3 * e2 + a4 * b2 * c1 * d3 * e5 - a4 * b2 * c1 * d5 * e3 -
                   a4 * b2 * c3 * d1 * e5 + a4 * b2 * c3 * d5 * e1 + a4 * b2 * c5 * d1 * e3 - a4 * b2 * c5 * d3 * e1 -
                   a4 * b3 * c1 * d2 * e5 + a4 * b3 * c1 * d5 * e2 + a4 * b3 * c2 * d1 * e5 - a4 * b3 * c2 * d5 * e1 -
                   a4 * b3 * c5 * d1 * e2 + a4 * b3 * c5 * d2 * e1 + a4 * b5 * c1 * d2 * e3 - a4 * b5 * c1 * d3 * e2 -
                   a4 * b5 * c2 * d1 * e3 + a4 * b5 * c2 * d3 * e1 + a4 * b5 * c3 * d1 * e2 - a4 * b5 * c3 * d2 * e1 +
                   a5 * b1 * c2 * d3 * e4 - a5 * b1 * c2 * d4 * e3 - a5 * b1 * c3 * d2 * e4 + a5 * b1 * c3 * d4 * e2 +
                   a5 * b1 * c4 * d2 * e3 - a5 * b1 * c4 * d3 * e2 - a5 * b2 * c1 * d3 * e4 + a5 * b2 * c1 * d4 * e3 +
                   a5 * b2 * c3 * d1 * e4 - a5 * b2 * c3 * d4 * e1 - a5 * b2 * c4 * d1 * e3 + a5 * b2 * c4 * d3 * e1 +
                   a5 * b3 * c1 * d2 * e4 - a5 * b3 * c1 * d4 * e2 - a5 * b3 * c2 * d1 * e4 + a5 * b3 * c2 * d4 * e1 +
                   a5 * b3 * c4 * d1 * e2 - a5 * b3 * c4 * d2 * e1 - a5 * b4 * c1 * d2 * e3 + a5 * b4 * c1 * d3 * e2 +
                   a5 * b4 * c2 * d1 * e3 - a5 * b4 * c2 * d3 * e1 - a5 * b4 * c3 * d1 * e2 + a5 * b4 * c3 * d2 * e1);
        T[0][1] = -t363 *
                  (b0 * c2 * d3 * e4 * f5 - b0 * c2 * d3 * e5 * f4 - b0 * c2 * d4 * e3 * f5 + b0 * c2 * d4 * e5 * f3 +
                   b0 * c2 * d5 * e3 * f4 - b0 * c2 * d5 * e4 * f3 - b0 * c3 * d2 * e4 * f5 + b0 * c3 * d2 * e5 * f4 +
                   b0 * c3 * d4 * e2 * f5 - b0 * c3 * d4 * e5 * f2 - b0 * c3 * d5 * e2 * f4 + b0 * c3 * d5 * e4 * f2 +
                   b0 * c4 * d2 * e3 * f5 - b0 * c4 * d2 * e5 * f3 - b0 * c4 * d3 * e2 * f5 + b0 * c4 * d3 * e5 * f2 +
                   b0 * c4 * d5 * e2 * f3 - b0 * c4 * d5 * e3 * f2 - b0 * c5 * d2 * e3 * f4 + b0 * c5 * d2 * e4 * f3 +
                   b0 * c5 * d3 * e2 * f4 - b0 * c5 * d3 * e4 * f2 - b0 * c5 * d4 * e2 * f3 + b0 * c5 * d4 * e3 * f2 -
                   b2 * c0 * d3 * e4 * f5 + b2 * c0 * d3 * e5 * f4 + b2 * c0 * d4 * e3 * f5 - b2 * c0 * d4 * e5 * f3 -
                   b2 * c0 * d5 * e3 * f4 + b2 * c0 * d5 * e4 * f3 + b2 * c3 * d0 * e4 * f5 - b2 * c3 * d0 * e5 * f4 -
                   b2 * c3 * d4 * e0 * f5 + b2 * c3 * d4 * e5 * f0 + b2 * c3 * d5 * e0 * f4 - b2 * c3 * d5 * e4 * f0 -
                   b2 * c4 * d0 * e3 * f5 + b2 * c4 * d0 * e5 * f3 + b2 * c4 * d3 * e0 * f5 - b2 * c4 * d3 * e5 * f0 -
                   b2 * c4 * d5 * e0 * f3 + b2 * c4 * d5 * e3 * f0 + b2 * c5 * d0 * e3 * f4 - b2 * c5 * d0 * e4 * f3 -
                   b2 * c5 * d3 * e0 * f4 + b2 * c5 * d3 * e4 * f0 + b2 * c5 * d4 * e0 * f3 - b2 * c5 * d4 * e3 * f0 +
                   b3 * c0 * d2 * e4 * f5 - b3 * c0 * d2 * e5 * f4 - b3 * c0 * d4 * e2 * f5 + b3 * c0 * d4 * e5 * f2 +
                   b3 * c0 * d5 * e2 * f4 - b3 * c0 * d5 * e4 * f2 - b3 * c2 * d0 * e4 * f5 + b3 * c2 * d0 * e5 * f4 +
                   b3 * c2 * d4 * e0 * f5 - b3 * c2 * d4 * e5 * f0 - b3 * c2 * d5 * e0 * f4 + b3 * c2 * d5 * e4 * f0 +
                   b3 * c4 * d0 * e2 * f5 - b3 * c4 * d0 * e5 * f2 - b3 * c4 * d2 * e0 * f5 + b3 * c4 * d2 * e5 * f0 +
                   b3 * c4 * d5 * e0 * f2 - b3 * c4 * d5 * e2 * f0 - b3 * c5 * d0 * e2 * f4 + b3 * c5 * d0 * e4 * f2 +
                   b3 * c5 * d2 * e0 * f4 - b3 * c5 * d2 * e4 * f0 - b3 * c5 * d4 * e0 * f2 + b3 * c5 * d4 * e2 * f0 -
                   b4 * c0 * d2 * e3 * f5 + b4 * c0 * d2 * e5 * f3 + b4 * c0 * d3 * e2 * f5 - b4 * c0 * d3 * e5 * f2 -
                   b4 * c0 * d5 * e2 * f3 + b4 * c0 * d5 * e3 * f2 + b4 * c2 * d0 * e3 * f5 - b4 * c2 * d0 * e5 * f3 -
                   b4 * c2 * d3 * e0 * f5 + b4 * c2 * d3 * e5 * f0 + b4 * c2 * d5 * e0 * f3 - b4 * c2 * d5 * e3 * f0 -
                   b4 * c3 * d0 * e2 * f5 + b4 * c3 * d0 * e5 * f2 + b4 * c3 * d2 * e0 * f5 - b4 * c3 * d2 * e5 * f0 -
                   b4 * c3 * d5 * e0 * f2 + b4 * c3 * d5 * e2 * f0 + b4 * c5 * d0 * e2 * f3 - b4 * c5 * d0 * e3 * f2 -
                   b4 * c5 * d2 * e0 * f3 + b4 * c5 * d2 * e3 * f0 + b4 * c5 * d3 * e0 * f2 - b4 * c5 * d3 * e2 * f0 +
                   b5 * c0 * d2 * e3 * f4 - b5 * c0 * d2 * e4 * f3 - b5 * c0 * d3 * e2 * f4 + b5 * c0 * d3 * e4 * f2 +
                   b5 * c0 * d4 * e2 * f3 - b5 * c0 * d4 * e3 * f2 - b5 * c2 * d0 * e3 * f4 + b5 * c2 * d0 * e4 * f3 +
                   b5 * c2 * d3 * e0 * f4 - b5 * c2 * d3 * e4 * f0 - b5 * c2 * d4 * e0 * f3 + b5 * c2 * d4 * e3 * f0 +
                   b5 * c3 * d0 * e2 * f4 - b5 * c3 * d0 * e4 * f2 - b5 * c3 * d2 * e0 * f4 + b5 * c3 * d2 * e4 * f0 +
                   b5 * c3 * d4 * e0 * f2 - b5 * c3 * d4 * e2 * f0 - b5 * c4 * d0 * e2 * f3 + b5 * c4 * d0 * e3 * f2 +
                   b5 * c4 * d2 * e0 * f3 - b5 * c4 * d2 * e3 * f0 - b5 * c4 * d3 * e0 * f2 + b5 * c4 * d3 * e2 * f0);
        T[1][1] = t363 *
                  (a0 * c2 * d3 * e4 * f5 - a0 * c2 * d3 * e5 * f4 - a0 * c2 * d4 * e3 * f5 + a0 * c2 * d4 * e5 * f3 +
                   a0 * c2 * d5 * e3 * f4 - a0 * c2 * d5 * e4 * f3 - a0 * c3 * d2 * e4 * f5 + a0 * c3 * d2 * e5 * f4 +
                   a0 * c3 * d4 * e2 * f5 - a0 * c3 * d4 * e5 * f2 - a0 * c3 * d5 * e2 * f4 + a0 * c3 * d5 * e4 * f2 +
                   a0 * c4 * d2 * e3 * f5 - a0 * c4 * d2 * e5 * f3 - a0 * c4 * d3 * e2 * f5 + a0 * c4 * d3 * e5 * f2 +
                   a0 * c4 * d5 * e2 * f3 - a0 * c4 * d5 * e3 * f2 - a0 * c5 * d2 * e3 * f4 + a0 * c5 * d2 * e4 * f3 +
                   a0 * c5 * d3 * e2 * f4 - a0 * c5 * d3 * e4 * f2 - a0 * c5 * d4 * e2 * f3 + a0 * c5 * d4 * e3 * f2 -
                   a2 * c0 * d3 * e4 * f5 + a2 * c0 * d3 * e5 * f4 + a2 * c0 * d4 * e3 * f5 - a2 * c0 * d4 * e5 * f3 -
                   a2 * c0 * d5 * e3 * f4 + a2 * c0 * d5 * e4 * f3 + a2 * c3 * d0 * e4 * f5 - a2 * c3 * d0 * e5 * f4 -
                   a2 * c3 * d4 * e0 * f5 + a2 * c3 * d4 * e5 * f0 + a2 * c3 * d5 * e0 * f4 - a2 * c3 * d5 * e4 * f0 -
                   a2 * c4 * d0 * e3 * f5 + a2 * c4 * d0 * e5 * f3 + a2 * c4 * d3 * e0 * f5 - a2 * c4 * d3 * e5 * f0 -
                   a2 * c4 * d5 * e0 * f3 + a2 * c4 * d5 * e3 * f0 + a2 * c5 * d0 * e3 * f4 - a2 * c5 * d0 * e4 * f3 -
                   a2 * c5 * d3 * e0 * f4 + a2 * c5 * d3 * e4 * f0 + a2 * c5 * d4 * e0 * f3 - a2 * c5 * d4 * e3 * f0 +
                   a3 * c0 * d2 * e4 * f5 - a3 * c0 * d2 * e5 * f4 - a3 * c0 * d4 * e2 * f5 + a3 * c0 * d4 * e5 * f2 +
                   a3 * c0 * d5 * e2 * f4 - a3 * c0 * d5 * e4 * f2 - a3 * c2 * d0 * e4 * f5 + a3 * c2 * d0 * e5 * f4 +
                   a3 * c2 * d4 * e0 * f5 - a3 * c2 * d4 * e5 * f0 - a3 * c2 * d5 * e0 * f4 + a3 * c2 * d5 * e4 * f0 +
                   a3 * c4 * d0 * e2 * f5 - a3 * c4 * d0 * e5 * f2 - a3 * c4 * d2 * e0 * f5 + a3 * c4 * d2 * e5 * f0 +
                   a3 * c4 * d5 * e0 * f2 - a3 * c4 * d5 * e2 * f0 - a3 * c5 * d0 * e2 * f4 + a3 * c5 * d0 * e4 * f2 +
                   a3 * c5 * d2 * e0 * f4 - a3 * c5 * d2 * e4 * f0 - a3 * c5 * d4 * e0 * f2 + a3 * c5 * d4 * e2 * f0 -
                   a4 * c0 * d2 * e3 * f5 + a4 * c0 * d2 * e5 * f3 + a4 * c0 * d3 * e2 * f5 - a4 * c0 * d3 * e5 * f2 -
                   a4 * c0 * d5 * e2 * f3 + a4 * c0 * d5 * e3 * f2 + a4 * c2 * d0 * e3 * f5 - a4 * c2 * d0 * e5 * f3 -
                   a4 * c2 * d3 * e0 * f5 + a4 * c2 * d3 * e5 * f0 + a4 * c2 * d5 * e0 * f3 - a4 * c2 * d5 * e3 * f0 -
                   a4 * c3 * d0 * e2 * f5 + a4 * c3 * d0 * e5 * f2 + a4 * c3 * d2 * e0 * f5 - a4 * c3 * d2 * e5 * f0 -
                   a4 * c3 * d5 * e0 * f2 + a4 * c3 * d5 * e2 * f0 + a4 * c5 * d0 * e2 * f3 - a4 * c5 * d0 * e3 * f2 -
                   a4 * c5 * d2 * e0 * f3 + a4 * c5 * d2 * e3 * f0 + a4 * c5 * d3 * e0 * f2 - a4 * c5 * d3 * e2 * f0 +
                   a5 * c0 * d2 * e3 * f4 - a5 * c0 * d2 * e4 * f3 - a5 * c0 * d3 * e2 * f4 + a5 * c0 * d3 * e4 * f2 +
                   a5 * c0 * d4 * e2 * f3 - a5 * c0 * d4 * e3 * f2 - a5 * c2 * d0 * e3 * f4 + a5 * c2 * d0 * e4 * f3 +
                   a5 * c2 * d3 * e0 * f4 - a5 * c2 * d3 * e4 * f0 - a5 * c2 * d4 * e0 * f3 + a5 * c2 * d4 * e3 * f0 +
                   a5 * c3 * d0 * e2 * f4 - a5 * c3 * d0 * e4 * f2 - a5 * c3 * d2 * e0 * f4 + a5 * c3 * d2 * e4 * f0 +
                   a5 * c3 * d4 * e0 * f2 - a5 * c3 * d4 * e2 * f0 - a5 * c4 * d0 * e2 * f3 + a5 * c4 * d0 * e3 * f2 +
                   a5 * c4 * d2 * e0 * f3 - a5 * c4 * d2 * e3 * f0 - a5 * c4 * d3 * e0 * f2 + a5 * c4 * d3 * e2 * f0);
        T[2][1] = -t363 *
                  (a0 * b2 * d3 * e4 * f5 - a0 * b2 * d3 * e5 * f4 - a0 * b2 * d4 * e3 * f5 + a0 * b2 * d4 * e5 * f3 +
                   a0 * b2 * d5 * e3 * f4 - a0 * b2 * d5 * e4 * f3 - a0 * b3 * d2 * e4 * f5 + a0 * b3 * d2 * e5 * f4 +
                   a0 * b3 * d4 * e2 * f5 - a0 * b3 * d4 * e5 * f2 - a0 * b3 * d5 * e2 * f4 + a0 * b3 * d5 * e4 * f2 +
                   a0 * b4 * d2 * e3 * f5 - a0 * b4 * d2 * e5 * f3 - a0 * b4 * d3 * e2 * f5 + a0 * b4 * d3 * e5 * f2 +
                   a0 * b4 * d5 * e2 * f3 - a0 * b4 * d5 * e3 * f2 - a0 * b5 * d2 * e3 * f4 + a0 * b5 * d2 * e4 * f3 +
                   a0 * b5 * d3 * e2 * f4 - a0 * b5 * d3 * e4 * f2 - a0 * b5 * d4 * e2 * f3 + a0 * b5 * d4 * e3 * f2 -
                   a2 * b0 * d3 * e4 * f5 + a2 * b0 * d3 * e5 * f4 + a2 * b0 * d4 * e3 * f5 - a2 * b0 * d4 * e5 * f3 -
                   a2 * b0 * d5 * e3 * f4 + a2 * b0 * d5 * e4 * f3 + a2 * b3 * d0 * e4 * f5 - a2 * b3 * d0 * e5 * f4 -
                   a2 * b3 * d4 * e0 * f5 + a2 * b3 * d4 * e5 * f0 + a2 * b3 * d5 * e0 * f4 - a2 * b3 * d5 * e4 * f0 -
                   a2 * b4 * d0 * e3 * f5 + a2 * b4 * d0 * e5 * f3 + a2 * b4 * d3 * e0 * f5 - a2 * b4 * d3 * e5 * f0 -
                   a2 * b4 * d5 * e0 * f3 + a2 * b4 * d5 * e3 * f0 + a2 * b5 * d0 * e3 * f4 - a2 * b5 * d0 * e4 * f3 -
                   a2 * b5 * d3 * e0 * f4 + a2 * b5 * d3 * e4 * f0 + a2 * b5 * d4 * e0 * f3 - a2 * b5 * d4 * e3 * f0 +
                   a3 * b0 * d2 * e4 * f5 - a3 * b0 * d2 * e5 * f4 - a3 * b0 * d4 * e2 * f5 + a3 * b0 * d4 * e5 * f2 +
                   a3 * b0 * d5 * e2 * f4 - a3 * b0 * d5 * e4 * f2 - a3 * b2 * d0 * e4 * f5 + a3 * b2 * d0 * e5 * f4 +
                   a3 * b2 * d4 * e0 * f5 - a3 * b2 * d4 * e5 * f0 - a3 * b2 * d5 * e0 * f4 + a3 * b2 * d5 * e4 * f0 +
                   a3 * b4 * d0 * e2 * f5 - a3 * b4 * d0 * e5 * f2 - a3 * b4 * d2 * e0 * f5 + a3 * b4 * d2 * e5 * f0 +
                   a3 * b4 * d5 * e0 * f2 - a3 * b4 * d5 * e2 * f0 - a3 * b5 * d0 * e2 * f4 + a3 * b5 * d0 * e4 * f2 +
                   a3 * b5 * d2 * e0 * f4 - a3 * b5 * d2 * e4 * f0 - a3 * b5 * d4 * e0 * f2 + a3 * b5 * d4 * e2 * f0 -
                   a4 * b0 * d2 * e3 * f5 + a4 * b0 * d2 * e5 * f3 + a4 * b0 * d3 * e2 * f5 - a4 * b0 * d3 * e5 * f2 -
                   a4 * b0 * d5 * e2 * f3 + a4 * b0 * d5 * e3 * f2 + a4 * b2 * d0 * e3 * f5 - a4 * b2 * d0 * e5 * f3 -
                   a4 * b2 * d3 * e0 * f5 + a4 * b2 * d3 * e5 * f0 + a4 * b2 * d5 * e0 * f3 - a4 * b2 * d5 * e3 * f0 -
                   a4 * b3 * d0 * e2 * f5 + a4 * b3 * d0 * e5 * f2 + a4 * b3 * d2 * e0 * f5 - a4 * b3 * d2 * e5 * f0 -
                   a4 * b3 * d5 * e0 * f2 + a4 * b3 * d5 * e2 * f0 + a4 * b5 * d0 * e2 * f3 - a4 * b5 * d0 * e3 * f2 -
                   a4 * b5 * d2 * e0 * f3 + a4 * b5 * d2 * e3 * f0 + a4 * b5 * d3 * e0 * f2 - a4 * b5 * d3 * e2 * f0 +
                   a5 * b0 * d2 * e3 * f4 - a5 * b0 * d2 * e4 * f3 - a5 * b0 * d3 * e2 * f4 + a5 * b0 * d3 * e4 * f2 +
                   a5 * b0 * d4 * e2 * f3 - a5 * b0 * d4 * e3 * f2 - a5 * b2 * d0 * e3 * f4 + a5 * b2 * d0 * e4 * f3 +
                   a5 * b2 * d3 * e0 * f4 - a5 * b2 * d3 * e4 * f0 - a5 * b2 * d4 * e0 * f3 + a5 * b2 * d4 * e3 * f0 +
                   a5 * b3 * d0 * e2 * f4 - a5 * b3 * d0 * e4 * f2 - a5 * b3 * d2 * e0 * f4 + a5 * b3 * d2 * e4 * f0 +
                   a5 * b3 * d4 * e0 * f2 - a5 * b3 * d4 * e2 * f0 - a5 * b4 * d0 * e2 * f3 + a5 * b4 * d0 * e3 * f2 +
                   a5 * b4 * d2 * e0 * f3 - a5 * b4 * d2 * e3 * f0 - a5 * b4 * d3 * e0 * f2 + a5 * b4 * d3 * e2 * f0);
        T[3][1] = t363 *
                  (a0 * b2 * c3 * e4 * f5 - a0 * b2 * c3 * e5 * f4 - a0 * b2 * c4 * e3 * f5 + a0 * b2 * c4 * e5 * f3 +
                   a0 * b2 * c5 * e3 * f4 - a0 * b2 * c5 * e4 * f3 - a0 * b3 * c2 * e4 * f5 + a0 * b3 * c2 * e5 * f4 +
                   a0 * b3 * c4 * e2 * f5 - a0 * b3 * c4 * e5 * f2 - a0 * b3 * c5 * e2 * f4 + a0 * b3 * c5 * e4 * f2 +
                   a0 * b4 * c2 * e3 * f5 - a0 * b4 * c2 * e5 * f3 - a0 * b4 * c3 * e2 * f5 + a0 * b4 * c3 * e5 * f2 +
                   a0 * b4 * c5 * e2 * f3 - a0 * b4 * c5 * e3 * f2 - a0 * b5 * c2 * e3 * f4 + a0 * b5 * c2 * e4 * f3 +
                   a0 * b5 * c3 * e2 * f4 - a0 * b5 * c3 * e4 * f2 - a0 * b5 * c4 * e2 * f3 + a0 * b5 * c4 * e3 * f2 -
                   a2 * b0 * c3 * e4 * f5 + a2 * b0 * c3 * e5 * f4 + a2 * b0 * c4 * e3 * f5 - a2 * b0 * c4 * e5 * f3 -
                   a2 * b0 * c5 * e3 * f4 + a2 * b0 * c5 * e4 * f3 + a2 * b3 * c0 * e4 * f5 - a2 * b3 * c0 * e5 * f4 -
                   a2 * b3 * c4 * e0 * f5 + a2 * b3 * c4 * e5 * f0 + a2 * b3 * c5 * e0 * f4 - a2 * b3 * c5 * e4 * f0 -
                   a2 * b4 * c0 * e3 * f5 + a2 * b4 * c0 * e5 * f3 + a2 * b4 * c3 * e0 * f5 - a2 * b4 * c3 * e5 * f0 -
                   a2 * b4 * c5 * e0 * f3 + a2 * b4 * c5 * e3 * f0 + a2 * b5 * c0 * e3 * f4 - a2 * b5 * c0 * e4 * f3 -
                   a2 * b5 * c3 * e0 * f4 + a2 * b5 * c3 * e4 * f0 + a2 * b5 * c4 * e0 * f3 - a2 * b5 * c4 * e3 * f0 +
                   a3 * b0 * c2 * e4 * f5 - a3 * b0 * c2 * e5 * f4 - a3 * b0 * c4 * e2 * f5 + a3 * b0 * c4 * e5 * f2 +
                   a3 * b0 * c5 * e2 * f4 - a3 * b0 * c5 * e4 * f2 - a3 * b2 * c0 * e4 * f5 + a3 * b2 * c0 * e5 * f4 +
                   a3 * b2 * c4 * e0 * f5 - a3 * b2 * c4 * e5 * f0 - a3 * b2 * c5 * e0 * f4 + a3 * b2 * c5 * e4 * f0 +
                   a3 * b4 * c0 * e2 * f5 - a3 * b4 * c0 * e5 * f2 - a3 * b4 * c2 * e0 * f5 + a3 * b4 * c2 * e5 * f0 +
                   a3 * b4 * c5 * e0 * f2 - a3 * b4 * c5 * e2 * f0 - a3 * b5 * c0 * e2 * f4 + a3 * b5 * c0 * e4 * f2 +
                   a3 * b5 * c2 * e0 * f4 - a3 * b5 * c2 * e4 * f0 - a3 * b5 * c4 * e0 * f2 + a3 * b5 * c4 * e2 * f0 -
                   a4 * b0 * c2 * e3 * f5 + a4 * b0 * c2 * e5 * f3 + a4 * b0 * c3 * e2 * f5 - a4 * b0 * c3 * e5 * f2 -
                   a4 * b0 * c5 * e2 * f3 + a4 * b0 * c5 * e3 * f2 + a4 * b2 * c0 * e3 * f5 - a4 * b2 * c0 * e5 * f3 -
                   a4 * b2 * c3 * e0 * f5 + a4 * b2 * c3 * e5 * f0 + a4 * b2 * c5 * e0 * f3 - a4 * b2 * c5 * e3 * f0 -
                   a4 * b3 * c0 * e2 * f5 + a4 * b3 * c0 * e5 * f2 + a4 * b3 * c2 * e0 * f5 - a4 * b3 * c2 * e5 * f0 -
                   a4 * b3 * c5 * e0 * f2 + a4 * b3 * c5 * e2 * f0 + a4 * b5 * c0 * e2 * f3 - a4 * b5 * c0 * e3 * f2 -
                   a4 * b5 * c2 * e0 * f3 + a4 * b5 * c2 * e3 * f0 + a4 * b5 * c3 * e0 * f2 - a4 * b5 * c3 * e2 * f0 +
                   a5 * b0 * c2 * e3 * f4 - a5 * b0 * c2 * e4 * f3 - a5 * b0 * c3 * e2 * f4 + a5 * b0 * c3 * e4 * f2 +
                   a5 * b0 * c4 * e2 * f3 - a5 * b0 * c4 * e3 * f2 - a5 * b2 * c0 * e3 * f4 + a5 * b2 * c0 * e4 * f3 +
                   a5 * b2 * c3 * e0 * f4 - a5 * b2 * c3 * e4 * f0 - a5 * b2 * c4 * e0 * f3 + a5 * b2 * c4 * e3 * f0 +
                   a5 * b3 * c0 * e2 * f4 - a5 * b3 * c0 * e4 * f2 - a5 * b3 * c2 * e0 * f4 + a5 * b3 * c2 * e4 * f0 +
                   a5 * b3 * c4 * e0 * f2 - a5 * b3 * c4 * e2 * f0 - a5 * b4 * c0 * e2 * f3 + a5 * b4 * c0 * e3 * f2 +
                   a5 * b4 * c2 * e0 * f3 - a5 * b4 * c2 * e3 * f0 - a5 * b4 * c3 * e0 * f2 + a5 * b4 * c3 * e2 * f0);
        T[4][1] = -t363 *
                  (a0 * b2 * c3 * d4 * f5 - a0 * b2 * c3 * d5 * f4 - a0 * b2 * c4 * d3 * f5 + a0 * b2 * c4 * d5 * f3 +
                   a0 * b2 * c5 * d3 * f4 - a0 * b2 * c5 * d4 * f3 - a0 * b3 * c2 * d4 * f5 + a0 * b3 * c2 * d5 * f4 +
                   a0 * b3 * c4 * d2 * f5 - a0 * b3 * c4 * d5 * f2 - a0 * b3 * c5 * d2 * f4 + a0 * b3 * c5 * d4 * f2 +
                   a0 * b4 * c2 * d3 * f5 - a0 * b4 * c2 * d5 * f3 - a0 * b4 * c3 * d2 * f5 + a0 * b4 * c3 * d5 * f2 +
                   a0 * b4 * c5 * d2 * f3 - a0 * b4 * c5 * d3 * f2 - a0 * b5 * c2 * d3 * f4 + a0 * b5 * c2 * d4 * f3 +
                   a0 * b5 * c3 * d2 * f4 - a0 * b5 * c3 * d4 * f2 - a0 * b5 * c4 * d2 * f3 + a0 * b5 * c4 * d3 * f2 -
                   a2 * b0 * c3 * d4 * f5 + a2 * b0 * c3 * d5 * f4 + a2 * b0 * c4 * d3 * f5 - a2 * b0 * c4 * d5 * f3 -
                   a2 * b0 * c5 * d3 * f4 + a2 * b0 * c5 * d4 * f3 + a2 * b3 * c0 * d4 * f5 - a2 * b3 * c0 * d5 * f4 -
                   a2 * b3 * c4 * d0 * f5 + a2 * b3 * c4 * d5 * f0 + a2 * b3 * c5 * d0 * f4 - a2 * b3 * c5 * d4 * f0 -
                   a2 * b4 * c0 * d3 * f5 + a2 * b4 * c0 * d5 * f3 + a2 * b4 * c3 * d0 * f5 - a2 * b4 * c3 * d5 * f0 -
                   a2 * b4 * c5 * d0 * f3 + a2 * b4 * c5 * d3 * f0 + a2 * b5 * c0 * d3 * f4 - a2 * b5 * c0 * d4 * f3 -
                   a2 * b5 * c3 * d0 * f4 + a2 * b5 * c3 * d4 * f0 + a2 * b5 * c4 * d0 * f3 - a2 * b5 * c4 * d3 * f0 +
                   a3 * b0 * c2 * d4 * f5 - a3 * b0 * c2 * d5 * f4 - a3 * b0 * c4 * d2 * f5 + a3 * b0 * c4 * d5 * f2 +
                   a3 * b0 * c5 * d2 * f4 - a3 * b0 * c5 * d4 * f2 - a3 * b2 * c0 * d4 * f5 + a3 * b2 * c0 * d5 * f4 +
                   a3 * b2 * c4 * d0 * f5 - a3 * b2 * c4 * d5 * f0 - a3 * b2 * c5 * d0 * f4 + a3 * b2 * c5 * d4 * f0 +
                   a3 * b4 * c0 * d2 * f5 - a3 * b4 * c0 * d5 * f2 - a3 * b4 * c2 * d0 * f5 + a3 * b4 * c2 * d5 * f0 +
                   a3 * b4 * c5 * d0 * f2 - a3 * b4 * c5 * d2 * f0 - a3 * b5 * c0 * d2 * f4 + a3 * b5 * c0 * d4 * f2 +
                   a3 * b5 * c2 * d0 * f4 - a3 * b5 * c2 * d4 * f0 - a3 * b5 * c4 * d0 * f2 + a3 * b5 * c4 * d2 * f0 -
                   a4 * b0 * c2 * d3 * f5 + a4 * b0 * c2 * d5 * f3 + a4 * b0 * c3 * d2 * f5 - a4 * b0 * c3 * d5 * f2 -
                   a4 * b0 * c5 * d2 * f3 + a4 * b0 * c5 * d3 * f2 + a4 * b2 * c0 * d3 * f5 - a4 * b2 * c0 * d5 * f3 -
                   a4 * b2 * c3 * d0 * f5 + a4 * b2 * c3 * d5 * f0 + a4 * b2 * c5 * d0 * f3 - a4 * b2 * c5 * d3 * f0 -
                   a4 * b3 * c0 * d2 * f5 + a4 * b3 * c0 * d5 * f2 + a4 * b3 * c2 * d0 * f5 - a4 * b3 * c2 * d5 * f0 -
                   a4 * b3 * c5 * d0 * f2 + a4 * b3 * c5 * d2 * f0 + a4 * b5 * c0 * d2 * f3 - a4 * b5 * c0 * d3 * f2 -
                   a4 * b5 * c2 * d0 * f3 + a4 * b5 * c2 * d3 * f0 + a4 * b5 * c3 * d0 * f2 - a4 * b5 * c3 * d2 * f0 +
                   a5 * b0 * c2 * d3 * f4 - a5 * b0 * c2 * d4 * f3 - a5 * b0 * c3 * d2 * f4 + a5 * b0 * c3 * d4 * f2 +
                   a5 * b0 * c4 * d2 * f3 - a5 * b0 * c4 * d3 * f2 - a5 * b2 * c0 * d3 * f4 + a5 * b2 * c0 * d4 * f3 +
                   a5 * b2 * c3 * d0 * f4 - a5 * b2 * c3 * d4 * f0 - a5 * b2 * c4 * d0 * f3 + a5 * b2 * c4 * d3 * f0 +
                   a5 * b3 * c0 * d2 * f4 - a5 * b3 * c0 * d4 * f2 - a5 * b3 * c2 * d0 * f4 + a5 * b3 * c2 * d4 * f0 +
                   a5 * b3 * c4 * d0 * f2 - a5 * b3 * c4 * d2 * f0 - a5 * b4 * c0 * d2 * f3 + a5 * b4 * c0 * d3 * f2 +
                   a5 * b4 * c2 * d0 * f3 - a5 * b4 * c2 * d3 * f0 - a5 * b4 * c3 * d0 * f2 + a5 * b4 * c3 * d2 * f0);
        T[5][1] = t363 *
                  (a0 * b2 * c3 * d4 * e5 - a0 * b2 * c3 * d5 * e4 - a0 * b2 * c4 * d3 * e5 + a0 * b2 * c4 * d5 * e3 +
                   a0 * b2 * c5 * d3 * e4 - a0 * b2 * c5 * d4 * e3 - a0 * b3 * c2 * d4 * e5 + a0 * b3 * c2 * d5 * e4 +
                   a0 * b3 * c4 * d2 * e5 - a0 * b3 * c4 * d5 * e2 - a0 * b3 * c5 * d2 * e4 + a0 * b3 * c5 * d4 * e2 +
                   a0 * b4 * c2 * d3 * e5 - a0 * b4 * c2 * d5 * e3 - a0 * b4 * c3 * d2 * e5 + a0 * b4 * c3 * d5 * e2 +
                   a0 * b4 * c5 * d2 * e3 - a0 * b4 * c5 * d3 * e2 - a0 * b5 * c2 * d3 * e4 + a0 * b5 * c2 * d4 * e3 +
                   a0 * b5 * c3 * d2 * e4 - a0 * b5 * c3 * d4 * e2 - a0 * b5 * c4 * d2 * e3 + a0 * b5 * c4 * d3 * e2 -
                   a2 * b0 * c3 * d4 * e5 + a2 * b0 * c3 * d5 * e4 + a2 * b0 * c4 * d3 * e5 - a2 * b0 * c4 * d5 * e3 -
                   a2 * b0 * c5 * d3 * e4 + a2 * b0 * c5 * d4 * e3 + a2 * b3 * c0 * d4 * e5 - a2 * b3 * c0 * d5 * e4 -
                   a2 * b3 * c4 * d0 * e5 + a2 * b3 * c4 * d5 * e0 + a2 * b3 * c5 * d0 * e4 - a2 * b3 * c5 * d4 * e0 -
                   a2 * b4 * c0 * d3 * e5 + a2 * b4 * c0 * d5 * e3 + a2 * b4 * c3 * d0 * e5 - a2 * b4 * c3 * d5 * e0 -
                   a2 * b4 * c5 * d0 * e3 + a2 * b4 * c5 * d3 * e0 + a2 * b5 * c0 * d3 * e4 - a2 * b5 * c0 * d4 * e3 -
                   a2 * b5 * c3 * d0 * e4 + a2 * b5 * c3 * d4 * e0 + a2 * b5 * c4 * d0 * e3 - a2 * b5 * c4 * d3 * e0 +
                   a3 * b0 * c2 * d4 * e5 - a3 * b0 * c2 * d5 * e4 - a3 * b0 * c4 * d2 * e5 + a3 * b0 * c4 * d5 * e2 +
                   a3 * b0 * c5 * d2 * e4 - a3 * b0 * c5 * d4 * e2 - a3 * b2 * c0 * d4 * e5 + a3 * b2 * c0 * d5 * e4 +
                   a3 * b2 * c4 * d0 * e5 - a3 * b2 * c4 * d5 * e0 - a3 * b2 * c5 * d0 * e4 + a3 * b2 * c5 * d4 * e0 +
                   a3 * b4 * c0 * d2 * e5 - a3 * b4 * c0 * d5 * e2 - a3 * b4 * c2 * d0 * e5 + a3 * b4 * c2 * d5 * e0 +
                   a3 * b4 * c5 * d0 * e2 - a3 * b4 * c5 * d2 * e0 - a3 * b5 * c0 * d2 * e4 + a3 * b5 * c0 * d4 * e2 +
                   a3 * b5 * c2 * d0 * e4 - a3 * b5 * c2 * d4 * e0 - a3 * b5 * c4 * d0 * e2 + a3 * b5 * c4 * d2 * e0 -
                   a4 * b0 * c2 * d3 * e5 + a4 * b0 * c2 * d5 * e3 + a4 * b0 * c3 * d2 * e5 - a4 * b0 * c3 * d5 * e2 -
                   a4 * b0 * c5 * d2 * e3 + a4 * b0 * c5 * d3 * e2 + a4 * b2 * c0 * d3 * e5 - a4 * b2 * c0 * d5 * e3 -
                   a4 * b2 * c3 * d0 * e5 + a4 * b2 * c3 * d5 * e0 + a4 * b2 * c5 * d0 * e3 - a4 * b2 * c5 * d3 * e0 -
                   a4 * b3 * c0 * d2 * e5 + a4 * b3 * c0 * d5 * e2 + a4 * b3 * c2 * d0 * e5 - a4 * b3 * c2 * d5 * e0 -
                   a4 * b3 * c5 * d0 * e2 + a4 * b3 * c5 * d2 * e0 + a4 * b5 * c0 * d2 * e3 - a4 * b5 * c0 * d3 * e2 -
                   a4 * b5 * c2 * d0 * e3 + a4 * b5 * c2 * d3 * e0 + a4 * b5 * c3 * d0 * e2 - a4 * b5 * c3 * d2 * e0 +
                   a5 * b0 * c2 * d3 * e4 - a5 * b0 * c2 * d4 * e3 - a5 * b0 * c3 * d2 * e4 + a5 * b0 * c3 * d4 * e2 +
                   a5 * b0 * c4 * d2 * e3 - a5 * b0 * c4 * d3 * e2 - a5 * b2 * c0 * d3 * e4 + a5 * b2 * c0 * d4 * e3 +
                   a5 * b2 * c3 * d0 * e4 - a5 * b2 * c3 * d4 * e0 - a5 * b2 * c4 * d0 * e3 + a5 * b2 * c4 * d3 * e0 +
                   a5 * b3 * c0 * d2 * e4 - a5 * b3 * c0 * d4 * e2 - a5 * b3 * c2 * d0 * e4 + a5 * b3 * c2 * d4 * e0 +
                   a5 * b3 * c4 * d0 * e2 - a5 * b3 * c4 * d2 * e0 - a5 * b4 * c0 * d2 * e3 + a5 * b4 * c0 * d3 * e2 +
                   a5 * b4 * c2 * d0 * e3 - a5 * b4 * c2 * d3 * e0 - a5 * b4 * c3 * d0 * e2 + a5 * b4 * c3 * d2 * e0);
        T[0][2] = t363 *
                  (b0 * c1 * d3 * e4 * f5 - b0 * c1 * d3 * e5 * f4 - b0 * c1 * d4 * e3 * f5 + b0 * c1 * d4 * e5 * f3 +
                   b0 * c1 * d5 * e3 * f4 - b0 * c1 * d5 * e4 * f3 - b0 * c3 * d1 * e4 * f5 + b0 * c3 * d1 * e5 * f4 +
                   b0 * c3 * d4 * e1 * f5 - b0 * c3 * d4 * e5 * f1 - b0 * c3 * d5 * e1 * f4 + b0 * c3 * d5 * e4 * f1 +
                   b0 * c4 * d1 * e3 * f5 - b0 * c4 * d1 * e5 * f3 - b0 * c4 * d3 * e1 * f5 + b0 * c4 * d3 * e5 * f1 +
                   b0 * c4 * d5 * e1 * f3 - b0 * c4 * d5 * e3 * f1 - b0 * c5 * d1 * e3 * f4 + b0 * c5 * d1 * e4 * f3 +
                   b0 * c5 * d3 * e1 * f4 - b0 * c5 * d3 * e4 * f1 - b0 * c5 * d4 * e1 * f3 + b0 * c5 * d4 * e3 * f1 -
                   b1 * c0 * d3 * e4 * f5 + b1 * c0 * d3 * e5 * f4 + b1 * c0 * d4 * e3 * f5 - b1 * c0 * d4 * e5 * f3 -
                   b1 * c0 * d5 * e3 * f4 + b1 * c0 * d5 * e4 * f3 + b1 * c3 * d0 * e4 * f5 - b1 * c3 * d0 * e5 * f4 -
                   b1 * c3 * d4 * e0 * f5 + b1 * c3 * d4 * e5 * f0 + b1 * c3 * d5 * e0 * f4 - b1 * c3 * d5 * e4 * f0 -
                   b1 * c4 * d0 * e3 * f5 + b1 * c4 * d0 * e5 * f3 + b1 * c4 * d3 * e0 * f5 - b1 * c4 * d3 * e5 * f0 -
                   b1 * c4 * d5 * e0 * f3 + b1 * c4 * d5 * e3 * f0 + b1 * c5 * d0 * e3 * f4 - b1 * c5 * d0 * e4 * f3 -
                   b1 * c5 * d3 * e0 * f4 + b1 * c5 * d3 * e4 * f0 + b1 * c5 * d4 * e0 * f3 - b1 * c5 * d4 * e3 * f0 +
                   b3 * c0 * d1 * e4 * f5 - b3 * c0 * d1 * e5 * f4 - b3 * c0 * d4 * e1 * f5 + b3 * c0 * d4 * e5 * f1 +
                   b3 * c0 * d5 * e1 * f4 - b3 * c0 * d5 * e4 * f1 - b3 * c1 * d0 * e4 * f5 + b3 * c1 * d0 * e5 * f4 +
                   b3 * c1 * d4 * e0 * f5 - b3 * c1 * d4 * e5 * f0 - b3 * c1 * d5 * e0 * f4 + b3 * c1 * d5 * e4 * f0 +
                   b3 * c4 * d0 * e1 * f5 - b3 * c4 * d0 * e5 * f1 - b3 * c4 * d1 * e0 * f5 + b3 * c4 * d1 * e5 * f0 +
                   b3 * c4 * d5 * e0 * f1 - b3 * c4 * d5 * e1 * f0 - b3 * c5 * d0 * e1 * f4 + b3 * c5 * d0 * e4 * f1 +
                   b3 * c5 * d1 * e0 * f4 - b3 * c5 * d1 * e4 * f0 - b3 * c5 * d4 * e0 * f1 + b3 * c5 * d4 * e1 * f0 -
                   b4 * c0 * d1 * e3 * f5 + b4 * c0 * d1 * e5 * f3 + b4 * c0 * d3 * e1 * f5 - b4 * c0 * d3 * e5 * f1 -
                   b4 * c0 * d5 * e1 * f3 + b4 * c0 * d5 * e3 * f1 + b4 * c1 * d0 * e3 * f5 - b4 * c1 * d0 * e5 * f3 -
                   b4 * c1 * d3 * e0 * f5 + b4 * c1 * d3 * e5 * f0 + b4 * c1 * d5 * e0 * f3 - b4 * c1 * d5 * e3 * f0 -
                   b4 * c3 * d0 * e1 * f5 + b4 * c3 * d0 * e5 * f1 + b4 * c3 * d1 * e0 * f5 - b4 * c3 * d1 * e5 * f0 -
                   b4 * c3 * d5 * e0 * f1 + b4 * c3 * d5 * e1 * f0 + b4 * c5 * d0 * e1 * f3 - b4 * c5 * d0 * e3 * f1 -
                   b4 * c5 * d1 * e0 * f3 + b4 * c5 * d1 * e3 * f0 + b4 * c5 * d3 * e0 * f1 - b4 * c5 * d3 * e1 * f0 +
                   b5 * c0 * d1 * e3 * f4 - b5 * c0 * d1 * e4 * f3 - b5 * c0 * d3 * e1 * f4 + b5 * c0 * d3 * e4 * f1 +
                   b5 * c0 * d4 * e1 * f3 - b5 * c0 * d4 * e3 * f1 - b5 * c1 * d0 * e3 * f4 + b5 * c1 * d0 * e4 * f3 +
                   b5 * c1 * d3 * e0 * f4 - b5 * c1 * d3 * e4 * f0 - b5 * c1 * d4 * e0 * f3 + b5 * c1 * d4 * e3 * f0 +
                   b5 * c3 * d0 * e1 * f4 - b5 * c3 * d0 * e4 * f1 - b5 * c3 * d1 * e0 * f4 + b5 * c3 * d1 * e4 * f0 +
                   b5 * c3 * d4 * e0 * f1 - b5 * c3 * d4 * e1 * f0 - b5 * c4 * d0 * e1 * f3 + b5 * c4 * d0 * e3 * f1 +
                   b5 * c4 * d1 * e0 * f3 - b5 * c4 * d1 * e3 * f0 - b5 * c4 * d3 * e0 * f1 + b5 * c4 * d3 * e1 * f0);
        T[1][2] = -t363 *
                  (a0 * c1 * d3 * e4 * f5 - a0 * c1 * d3 * e5 * f4 - a0 * c1 * d4 * e3 * f5 + a0 * c1 * d4 * e5 * f3 +
                   a0 * c1 * d5 * e3 * f4 - a0 * c1 * d5 * e4 * f3 - a0 * c3 * d1 * e4 * f5 + a0 * c3 * d1 * e5 * f4 +
                   a0 * c3 * d4 * e1 * f5 - a0 * c3 * d4 * e5 * f1 - a0 * c3 * d5 * e1 * f4 + a0 * c3 * d5 * e4 * f1 +
                   a0 * c4 * d1 * e3 * f5 - a0 * c4 * d1 * e5 * f3 - a0 * c4 * d3 * e1 * f5 + a0 * c4 * d3 * e5 * f1 +
                   a0 * c4 * d5 * e1 * f3 - a0 * c4 * d5 * e3 * f1 - a0 * c5 * d1 * e3 * f4 + a0 * c5 * d1 * e4 * f3 +
                   a0 * c5 * d3 * e1 * f4 - a0 * c5 * d3 * e4 * f1 - a0 * c5 * d4 * e1 * f3 + a0 * c5 * d4 * e3 * f1 -
                   a1 * c0 * d3 * e4 * f5 + a1 * c0 * d3 * e5 * f4 + a1 * c0 * d4 * e3 * f5 - a1 * c0 * d4 * e5 * f3 -
                   a1 * c0 * d5 * e3 * f4 + a1 * c0 * d5 * e4 * f3 + a1 * c3 * d0 * e4 * f5 - a1 * c3 * d0 * e5 * f4 -
                   a1 * c3 * d4 * e0 * f5 + a1 * c3 * d4 * e5 * f0 + a1 * c3 * d5 * e0 * f4 - a1 * c3 * d5 * e4 * f0 -
                   a1 * c4 * d0 * e3 * f5 + a1 * c4 * d0 * e5 * f3 + a1 * c4 * d3 * e0 * f5 - a1 * c4 * d3 * e5 * f0 -
                   a1 * c4 * d5 * e0 * f3 + a1 * c4 * d5 * e3 * f0 + a1 * c5 * d0 * e3 * f4 - a1 * c5 * d0 * e4 * f3 -
                   a1 * c5 * d3 * e0 * f4 + a1 * c5 * d3 * e4 * f0 + a1 * c5 * d4 * e0 * f3 - a1 * c5 * d4 * e3 * f0 +
                   a3 * c0 * d1 * e4 * f5 - a3 * c0 * d1 * e5 * f4 - a3 * c0 * d4 * e1 * f5 + a3 * c0 * d4 * e5 * f1 +
                   a3 * c0 * d5 * e1 * f4 - a3 * c0 * d5 * e4 * f1 - a3 * c1 * d0 * e4 * f5 + a3 * c1 * d0 * e5 * f4 +
                   a3 * c1 * d4 * e0 * f5 - a3 * c1 * d4 * e5 * f0 - a3 * c1 * d5 * e0 * f4 + a3 * c1 * d5 * e4 * f0 +
                   a3 * c4 * d0 * e1 * f5 - a3 * c4 * d0 * e5 * f1 - a3 * c4 * d1 * e0 * f5 + a3 * c4 * d1 * e5 * f0 +
                   a3 * c4 * d5 * e0 * f1 - a3 * c4 * d5 * e1 * f0 - a3 * c5 * d0 * e1 * f4 + a3 * c5 * d0 * e4 * f1 +
                   a3 * c5 * d1 * e0 * f4 - a3 * c5 * d1 * e4 * f0 - a3 * c5 * d4 * e0 * f1 + a3 * c5 * d4 * e1 * f0 -
                   a4 * c0 * d1 * e3 * f5 + a4 * c0 * d1 * e5 * f3 + a4 * c0 * d3 * e1 * f5 - a4 * c0 * d3 * e5 * f1 -
                   a4 * c0 * d5 * e1 * f3 + a4 * c0 * d5 * e3 * f1 + a4 * c1 * d0 * e3 * f5 - a4 * c1 * d0 * e5 * f3 -
                   a4 * c1 * d3 * e0 * f5 + a4 * c1 * d3 * e5 * f0 + a4 * c1 * d5 * e0 * f3 - a4 * c1 * d5 * e3 * f0 -
                   a4 * c3 * d0 * e1 * f5 + a4 * c3 * d0 * e5 * f1 + a4 * c3 * d1 * e0 * f5 - a4 * c3 * d1 * e5 * f0 -
                   a4 * c3 * d5 * e0 * f1 + a4 * c3 * d5 * e1 * f0 + a4 * c5 * d0 * e1 * f3 - a4 * c5 * d0 * e3 * f1 -
                   a4 * c5 * d1 * e0 * f3 + a4 * c5 * d1 * e3 * f0 + a4 * c5 * d3 * e0 * f1 - a4 * c5 * d3 * e1 * f0 +
                   a5 * c0 * d1 * e3 * f4 - a5 * c0 * d1 * e4 * f3 - a5 * c0 * d3 * e1 * f4 + a5 * c0 * d3 * e4 * f1 +
                   a5 * c0 * d4 * e1 * f3 - a5 * c0 * d4 * e3 * f1 - a5 * c1 * d0 * e3 * f4 + a5 * c1 * d0 * e4 * f3 +
                   a5 * c1 * d3 * e0 * f4 - a5 * c1 * d3 * e4 * f0 - a5 * c1 * d4 * e0 * f3 + a5 * c1 * d4 * e3 * f0 +
                   a5 * c3 * d0 * e1 * f4 - a5 * c3 * d0 * e4 * f1 - a5 * c3 * d1 * e0 * f4 + a5 * c3 * d1 * e4 * f0 +
                   a5 * c3 * d4 * e0 * f1 - a5 * c3 * d4 * e1 * f0 - a5 * c4 * d0 * e1 * f3 + a5 * c4 * d0 * e3 * f1 +
                   a5 * c4 * d1 * e0 * f3 - a5 * c4 * d1 * e3 * f0 - a5 * c4 * d3 * e0 * f1 + a5 * c4 * d3 * e1 * f0);
        T[2][2] = t363 *
                  (a0 * b1 * d3 * e4 * f5 - a0 * b1 * d3 * e5 * f4 - a0 * b1 * d4 * e3 * f5 + a0 * b1 * d4 * e5 * f3 +
                   a0 * b1 * d5 * e3 * f4 - a0 * b1 * d5 * e4 * f3 - a0 * b3 * d1 * e4 * f5 + a0 * b3 * d1 * e5 * f4 +
                   a0 * b3 * d4 * e1 * f5 - a0 * b3 * d4 * e5 * f1 - a0 * b3 * d5 * e1 * f4 + a0 * b3 * d5 * e4 * f1 +
                   a0 * b4 * d1 * e3 * f5 - a0 * b4 * d1 * e5 * f3 - a0 * b4 * d3 * e1 * f5 + a0 * b4 * d3 * e5 * f1 +
                   a0 * b4 * d5 * e1 * f3 - a0 * b4 * d5 * e3 * f1 - a0 * b5 * d1 * e3 * f4 + a0 * b5 * d1 * e4 * f3 +
                   a0 * b5 * d3 * e1 * f4 - a0 * b5 * d3 * e4 * f1 - a0 * b5 * d4 * e1 * f3 + a0 * b5 * d4 * e3 * f1 -
                   a1 * b0 * d3 * e4 * f5 + a1 * b0 * d3 * e5 * f4 + a1 * b0 * d4 * e3 * f5 - a1 * b0 * d4 * e5 * f3 -
                   a1 * b0 * d5 * e3 * f4 + a1 * b0 * d5 * e4 * f3 + a1 * b3 * d0 * e4 * f5 - a1 * b3 * d0 * e5 * f4 -
                   a1 * b3 * d4 * e0 * f5 + a1 * b3 * d4 * e5 * f0 + a1 * b3 * d5 * e0 * f4 - a1 * b3 * d5 * e4 * f0 -
                   a1 * b4 * d0 * e3 * f5 + a1 * b4 * d0 * e5 * f3 + a1 * b4 * d3 * e0 * f5 - a1 * b4 * d3 * e5 * f0 -
                   a1 * b4 * d5 * e0 * f3 + a1 * b4 * d5 * e3 * f0 + a1 * b5 * d0 * e3 * f4 - a1 * b5 * d0 * e4 * f3 -
                   a1 * b5 * d3 * e0 * f4 + a1 * b5 * d3 * e4 * f0 + a1 * b5 * d4 * e0 * f3 - a1 * b5 * d4 * e3 * f0 +
                   a3 * b0 * d1 * e4 * f5 - a3 * b0 * d1 * e5 * f4 - a3 * b0 * d4 * e1 * f5 + a3 * b0 * d4 * e5 * f1 +
                   a3 * b0 * d5 * e1 * f4 - a3 * b0 * d5 * e4 * f1 - a3 * b1 * d0 * e4 * f5 + a3 * b1 * d0 * e5 * f4 +
                   a3 * b1 * d4 * e0 * f5 - a3 * b1 * d4 * e5 * f0 - a3 * b1 * d5 * e0 * f4 + a3 * b1 * d5 * e4 * f0 +
                   a3 * b4 * d0 * e1 * f5 - a3 * b4 * d0 * e5 * f1 - a3 * b4 * d1 * e0 * f5 + a3 * b4 * d1 * e5 * f0 +
                   a3 * b4 * d5 * e0 * f1 - a3 * b4 * d5 * e1 * f0 - a3 * b5 * d0 * e1 * f4 + a3 * b5 * d0 * e4 * f1 +
                   a3 * b5 * d1 * e0 * f4 - a3 * b5 * d1 * e4 * f0 - a3 * b5 * d4 * e0 * f1 + a3 * b5 * d4 * e1 * f0 -
                   a4 * b0 * d1 * e3 * f5 + a4 * b0 * d1 * e5 * f3 + a4 * b0 * d3 * e1 * f5 - a4 * b0 * d3 * e5 * f1 -
                   a4 * b0 * d5 * e1 * f3 + a4 * b0 * d5 * e3 * f1 + a4 * b1 * d0 * e3 * f5 - a4 * b1 * d0 * e5 * f3 -
                   a4 * b1 * d3 * e0 * f5 + a4 * b1 * d3 * e5 * f0 + a4 * b1 * d5 * e0 * f3 - a4 * b1 * d5 * e3 * f0 -
                   a4 * b3 * d0 * e1 * f5 + a4 * b3 * d0 * e5 * f1 + a4 * b3 * d1 * e0 * f5 - a4 * b3 * d1 * e5 * f0 -
                   a4 * b3 * d5 * e0 * f1 + a4 * b3 * d5 * e1 * f0 + a4 * b5 * d0 * e1 * f3 - a4 * b5 * d0 * e3 * f1 -
                   a4 * b5 * d1 * e0 * f3 + a4 * b5 * d1 * e3 * f0 + a4 * b5 * d3 * e0 * f1 - a4 * b5 * d3 * e1 * f0 +
                   a5 * b0 * d1 * e3 * f4 - a5 * b0 * d1 * e4 * f3 - a5 * b0 * d3 * e1 * f4 + a5 * b0 * d3 * e4 * f1 +
                   a5 * b0 * d4 * e1 * f3 - a5 * b0 * d4 * e3 * f1 - a5 * b1 * d0 * e3 * f4 + a5 * b1 * d0 * e4 * f3 +
                   a5 * b1 * d3 * e0 * f4 - a5 * b1 * d3 * e4 * f0 - a5 * b1 * d4 * e0 * f3 + a5 * b1 * d4 * e3 * f0 +
                   a5 * b3 * d0 * e1 * f4 - a5 * b3 * d0 * e4 * f1 - a5 * b3 * d1 * e0 * f4 + a5 * b3 * d1 * e4 * f0 +
                   a5 * b3 * d4 * e0 * f1 - a5 * b3 * d4 * e1 * f0 - a5 * b4 * d0 * e1 * f3 + a5 * b4 * d0 * e3 * f1 +
                   a5 * b4 * d1 * e0 * f3 - a5 * b4 * d1 * e3 * f0 - a5 * b4 * d3 * e0 * f1 + a5 * b4 * d3 * e1 * f0);
        T[3][2] = -t363 *
                  (a0 * b1 * c3 * e4 * f5 - a0 * b1 * c3 * e5 * f4 - a0 * b1 * c4 * e3 * f5 + a0 * b1 * c4 * e5 * f3 +
                   a0 * b1 * c5 * e3 * f4 - a0 * b1 * c5 * e4 * f3 - a0 * b3 * c1 * e4 * f5 + a0 * b3 * c1 * e5 * f4 +
                   a0 * b3 * c4 * e1 * f5 - a0 * b3 * c4 * e5 * f1 - a0 * b3 * c5 * e1 * f4 + a0 * b3 * c5 * e4 * f1 +
                   a0 * b4 * c1 * e3 * f5 - a0 * b4 * c1 * e5 * f3 - a0 * b4 * c3 * e1 * f5 + a0 * b4 * c3 * e5 * f1 +
                   a0 * b4 * c5 * e1 * f3 - a0 * b4 * c5 * e3 * f1 - a0 * b5 * c1 * e3 * f4 + a0 * b5 * c1 * e4 * f3 +
                   a0 * b5 * c3 * e1 * f4 - a0 * b5 * c3 * e4 * f1 - a0 * b5 * c4 * e1 * f3 + a0 * b5 * c4 * e3 * f1 -
                   a1 * b0 * c3 * e4 * f5 + a1 * b0 * c3 * e5 * f4 + a1 * b0 * c4 * e3 * f5 - a1 * b0 * c4 * e5 * f3 -
                   a1 * b0 * c5 * e3 * f4 + a1 * b0 * c5 * e4 * f3 + a1 * b3 * c0 * e4 * f5 - a1 * b3 * c0 * e5 * f4 -
                   a1 * b3 * c4 * e0 * f5 + a1 * b3 * c4 * e5 * f0 + a1 * b3 * c5 * e0 * f4 - a1 * b3 * c5 * e4 * f0 -
                   a1 * b4 * c0 * e3 * f5 + a1 * b4 * c0 * e5 * f3 + a1 * b4 * c3 * e0 * f5 - a1 * b4 * c3 * e5 * f0 -
                   a1 * b4 * c5 * e0 * f3 + a1 * b4 * c5 * e3 * f0 + a1 * b5 * c0 * e3 * f4 - a1 * b5 * c0 * e4 * f3 -
                   a1 * b5 * c3 * e0 * f4 + a1 * b5 * c3 * e4 * f0 + a1 * b5 * c4 * e0 * f3 - a1 * b5 * c4 * e3 * f0 +
                   a3 * b0 * c1 * e4 * f5 - a3 * b0 * c1 * e5 * f4 - a3 * b0 * c4 * e1 * f5 + a3 * b0 * c4 * e5 * f1 +
                   a3 * b0 * c5 * e1 * f4 - a3 * b0 * c5 * e4 * f1 - a3 * b1 * c0 * e4 * f5 + a3 * b1 * c0 * e5 * f4 +
                   a3 * b1 * c4 * e0 * f5 - a3 * b1 * c4 * e5 * f0 - a3 * b1 * c5 * e0 * f4 + a3 * b1 * c5 * e4 * f0 +
                   a3 * b4 * c0 * e1 * f5 - a3 * b4 * c0 * e5 * f1 - a3 * b4 * c1 * e0 * f5 + a3 * b4 * c1 * e5 * f0 +
                   a3 * b4 * c5 * e0 * f1 - a3 * b4 * c5 * e1 * f0 - a3 * b5 * c0 * e1 * f4 + a3 * b5 * c0 * e4 * f1 +
                   a3 * b5 * c1 * e0 * f4 - a3 * b5 * c1 * e4 * f0 - a3 * b5 * c4 * e0 * f1 + a3 * b5 * c4 * e1 * f0 -
                   a4 * b0 * c1 * e3 * f5 + a4 * b0 * c1 * e5 * f3 + a4 * b0 * c3 * e1 * f5 - a4 * b0 * c3 * e5 * f1 -
                   a4 * b0 * c5 * e1 * f3 + a4 * b0 * c5 * e3 * f1 + a4 * b1 * c0 * e3 * f5 - a4 * b1 * c0 * e5 * f3 -
                   a4 * b1 * c3 * e0 * f5 + a4 * b1 * c3 * e5 * f0 + a4 * b1 * c5 * e0 * f3 - a4 * b1 * c5 * e3 * f0 -
                   a4 * b3 * c0 * e1 * f5 + a4 * b3 * c0 * e5 * f1 + a4 * b3 * c1 * e0 * f5 - a4 * b3 * c1 * e5 * f0 -
                   a4 * b3 * c5 * e0 * f1 + a4 * b3 * c5 * e1 * f0 + a4 * b5 * c0 * e1 * f3 - a4 * b5 * c0 * e3 * f1 -
                   a4 * b5 * c1 * e0 * f3 + a4 * b5 * c1 * e3 * f0 + a4 * b5 * c3 * e0 * f1 - a4 * b5 * c3 * e1 * f0 +
                   a5 * b0 * c1 * e3 * f4 - a5 * b0 * c1 * e4 * f3 - a5 * b0 * c3 * e1 * f4 + a5 * b0 * c3 * e4 * f1 +
                   a5 * b0 * c4 * e1 * f3 - a5 * b0 * c4 * e3 * f1 - a5 * b1 * c0 * e3 * f4 + a5 * b1 * c0 * e4 * f3 +
                   a5 * b1 * c3 * e0 * f4 - a5 * b1 * c3 * e4 * f0 - a5 * b1 * c4 * e0 * f3 + a5 * b1 * c4 * e3 * f0 +
                   a5 * b3 * c0 * e1 * f4 - a5 * b3 * c0 * e4 * f1 - a5 * b3 * c1 * e0 * f4 + a5 * b3 * c1 * e4 * f0 +
                   a5 * b3 * c4 * e0 * f1 - a5 * b3 * c4 * e1 * f0 - a5 * b4 * c0 * e1 * f3 + a5 * b4 * c0 * e3 * f1 +
                   a5 * b4 * c1 * e0 * f3 - a5 * b4 * c1 * e3 * f0 - a5 * b4 * c3 * e0 * f1 + a5 * b4 * c3 * e1 * f0);
        T[4][2] = t363 *
                  (a0 * b1 * c3 * d4 * f5 - a0 * b1 * c3 * d5 * f4 - a0 * b1 * c4 * d3 * f5 + a0 * b1 * c4 * d5 * f3 +
                   a0 * b1 * c5 * d3 * f4 - a0 * b1 * c5 * d4 * f3 - a0 * b3 * c1 * d4 * f5 + a0 * b3 * c1 * d5 * f4 +
                   a0 * b3 * c4 * d1 * f5 - a0 * b3 * c4 * d5 * f1 - a0 * b3 * c5 * d1 * f4 + a0 * b3 * c5 * d4 * f1 +
                   a0 * b4 * c1 * d3 * f5 - a0 * b4 * c1 * d5 * f3 - a0 * b4 * c3 * d1 * f5 + a0 * b4 * c3 * d5 * f1 +
                   a0 * b4 * c5 * d1 * f3 - a0 * b4 * c5 * d3 * f1 - a0 * b5 * c1 * d3 * f4 + a0 * b5 * c1 * d4 * f3 +
                   a0 * b5 * c3 * d1 * f4 - a0 * b5 * c3 * d4 * f1 - a0 * b5 * c4 * d1 * f3 + a0 * b5 * c4 * d3 * f1 -
                   a1 * b0 * c3 * d4 * f5 + a1 * b0 * c3 * d5 * f4 + a1 * b0 * c4 * d3 * f5 - a1 * b0 * c4 * d5 * f3 -
                   a1 * b0 * c5 * d3 * f4 + a1 * b0 * c5 * d4 * f3 + a1 * b3 * c0 * d4 * f5 - a1 * b3 * c0 * d5 * f4 -
                   a1 * b3 * c4 * d0 * f5 + a1 * b3 * c4 * d5 * f0 + a1 * b3 * c5 * d0 * f4 - a1 * b3 * c5 * d4 * f0 -
                   a1 * b4 * c0 * d3 * f5 + a1 * b4 * c0 * d5 * f3 + a1 * b4 * c3 * d0 * f5 - a1 * b4 * c3 * d5 * f0 -
                   a1 * b4 * c5 * d0 * f3 + a1 * b4 * c5 * d3 * f0 + a1 * b5 * c0 * d3 * f4 - a1 * b5 * c0 * d4 * f3 -
                   a1 * b5 * c3 * d0 * f4 + a1 * b5 * c3 * d4 * f0 + a1 * b5 * c4 * d0 * f3 - a1 * b5 * c4 * d3 * f0 +
                   a3 * b0 * c1 * d4 * f5 - a3 * b0 * c1 * d5 * f4 - a3 * b0 * c4 * d1 * f5 + a3 * b0 * c4 * d5 * f1 +
                   a3 * b0 * c5 * d1 * f4 - a3 * b0 * c5 * d4 * f1 - a3 * b1 * c0 * d4 * f5 + a3 * b1 * c0 * d5 * f4 +
                   a3 * b1 * c4 * d0 * f5 - a3 * b1 * c4 * d5 * f0 - a3 * b1 * c5 * d0 * f4 + a3 * b1 * c5 * d4 * f0 +
                   a3 * b4 * c0 * d1 * f5 - a3 * b4 * c0 * d5 * f1 - a3 * b4 * c1 * d0 * f5 + a3 * b4 * c1 * d5 * f0 +
                   a3 * b4 * c5 * d0 * f1 - a3 * b4 * c5 * d1 * f0 - a3 * b5 * c0 * d1 * f4 + a3 * b5 * c0 * d4 * f1 +
                   a3 * b5 * c1 * d0 * f4 - a3 * b5 * c1 * d4 * f0 - a3 * b5 * c4 * d0 * f1 + a3 * b5 * c4 * d1 * f0 -
                   a4 * b0 * c1 * d3 * f5 + a4 * b0 * c1 * d5 * f3 + a4 * b0 * c3 * d1 * f5 - a4 * b0 * c3 * d5 * f1 -
                   a4 * b0 * c5 * d1 * f3 + a4 * b0 * c5 * d3 * f1 + a4 * b1 * c0 * d3 * f5 - a4 * b1 * c0 * d5 * f3 -
                   a4 * b1 * c3 * d0 * f5 + a4 * b1 * c3 * d5 * f0 + a4 * b1 * c5 * d0 * f3 - a4 * b1 * c5 * d3 * f0 -
                   a4 * b3 * c0 * d1 * f5 + a4 * b3 * c0 * d5 * f1 + a4 * b3 * c1 * d0 * f5 - a4 * b3 * c1 * d5 * f0 -
                   a4 * b3 * c5 * d0 * f1 + a4 * b3 * c5 * d1 * f0 + a4 * b5 * c0 * d1 * f3 - a4 * b5 * c0 * d3 * f1 -
                   a4 * b5 * c1 * d0 * f3 + a4 * b5 * c1 * d3 * f0 + a4 * b5 * c3 * d0 * f1 - a4 * b5 * c3 * d1 * f0 +
                   a5 * b0 * c1 * d3 * f4 - a5 * b0 * c1 * d4 * f3 - a5 * b0 * c3 * d1 * f4 + a5 * b0 * c3 * d4 * f1 +
                   a5 * b0 * c4 * d1 * f3 - a5 * b0 * c4 * d3 * f1 - a5 * b1 * c0 * d3 * f4 + a5 * b1 * c0 * d4 * f3 +
                   a5 * b1 * c3 * d0 * f4 - a5 * b1 * c3 * d4 * f0 - a5 * b1 * c4 * d0 * f3 + a5 * b1 * c4 * d3 * f0 +
                   a5 * b3 * c0 * d1 * f4 - a5 * b3 * c0 * d4 * f1 - a5 * b3 * c1 * d0 * f4 + a5 * b3 * c1 * d4 * f0 +
                   a5 * b3 * c4 * d0 * f1 - a5 * b3 * c4 * d1 * f0 - a5 * b4 * c0 * d1 * f3 + a5 * b4 * c0 * d3 * f1 +
                   a5 * b4 * c1 * d0 * f3 - a5 * b4 * c1 * d3 * f0 - a5 * b4 * c3 * d0 * f1 + a5 * b4 * c3 * d1 * f0);
        T[5][2] = -t363 *
                  (a0 * b1 * c3 * d4 * e5 - a0 * b1 * c3 * d5 * e4 - a0 * b1 * c4 * d3 * e5 + a0 * b1 * c4 * d5 * e3 +
                   a0 * b1 * c5 * d3 * e4 - a0 * b1 * c5 * d4 * e3 - a0 * b3 * c1 * d4 * e5 + a0 * b3 * c1 * d5 * e4 +
                   a0 * b3 * c4 * d1 * e5 - a0 * b3 * c4 * d5 * e1 - a0 * b3 * c5 * d1 * e4 + a0 * b3 * c5 * d4 * e1 +
                   a0 * b4 * c1 * d3 * e5 - a0 * b4 * c1 * d5 * e3 - a0 * b4 * c3 * d1 * e5 + a0 * b4 * c3 * d5 * e1 +
                   a0 * b4 * c5 * d1 * e3 - a0 * b4 * c5 * d3 * e1 - a0 * b5 * c1 * d3 * e4 + a0 * b5 * c1 * d4 * e3 +
                   a0 * b5 * c3 * d1 * e4 - a0 * b5 * c3 * d4 * e1 - a0 * b5 * c4 * d1 * e3 + a0 * b5 * c4 * d3 * e1 -
                   a1 * b0 * c3 * d4 * e5 + a1 * b0 * c3 * d5 * e4 + a1 * b0 * c4 * d3 * e5 - a1 * b0 * c4 * d5 * e3 -
                   a1 * b0 * c5 * d3 * e4 + a1 * b0 * c5 * d4 * e3 + a1 * b3 * c0 * d4 * e5 - a1 * b3 * c0 * d5 * e4 -
                   a1 * b3 * c4 * d0 * e5 + a1 * b3 * c4 * d5 * e0 + a1 * b3 * c5 * d0 * e4 - a1 * b3 * c5 * d4 * e0 -
                   a1 * b4 * c0 * d3 * e5 + a1 * b4 * c0 * d5 * e3 + a1 * b4 * c3 * d0 * e5 - a1 * b4 * c3 * d5 * e0 -
                   a1 * b4 * c5 * d0 * e3 + a1 * b4 * c5 * d3 * e0 + a1 * b5 * c0 * d3 * e4 - a1 * b5 * c0 * d4 * e3 -
                   a1 * b5 * c3 * d0 * e4 + a1 * b5 * c3 * d4 * e0 + a1 * b5 * c4 * d0 * e3 - a1 * b5 * c4 * d3 * e0 +
                   a3 * b0 * c1 * d4 * e5 - a3 * b0 * c1 * d5 * e4 - a3 * b0 * c4 * d1 * e5 + a3 * b0 * c4 * d5 * e1 +
                   a3 * b0 * c5 * d1 * e4 - a3 * b0 * c5 * d4 * e1 - a3 * b1 * c0 * d4 * e5 + a3 * b1 * c0 * d5 * e4 +
                   a3 * b1 * c4 * d0 * e5 - a3 * b1 * c4 * d5 * e0 - a3 * b1 * c5 * d0 * e4 + a3 * b1 * c5 * d4 * e0 +
                   a3 * b4 * c0 * d1 * e5 - a3 * b4 * c0 * d5 * e1 - a3 * b4 * c1 * d0 * e5 + a3 * b4 * c1 * d5 * e0 +
                   a3 * b4 * c5 * d0 * e1 - a3 * b4 * c5 * d1 * e0 - a3 * b5 * c0 * d1 * e4 + a3 * b5 * c0 * d4 * e1 +
                   a3 * b5 * c1 * d0 * e4 - a3 * b5 * c1 * d4 * e0 - a3 * b5 * c4 * d0 * e1 + a3 * b5 * c4 * d1 * e0 -
                   a4 * b0 * c1 * d3 * e5 + a4 * b0 * c1 * d5 * e3 + a4 * b0 * c3 * d1 * e5 - a4 * b0 * c3 * d5 * e1 -
                   a4 * b0 * c5 * d1 * e3 + a4 * b0 * c5 * d3 * e1 + a4 * b1 * c0 * d3 * e5 - a4 * b1 * c0 * d5 * e3 -
                   a4 * b1 * c3 * d0 * e5 + a4 * b1 * c3 * d5 * e0 + a4 * b1 * c5 * d0 * e3 - a4 * b1 * c5 * d3 * e0 -
                   a4 * b3 * c0 * d1 * e5 + a4 * b3 * c0 * d5 * e1 + a4 * b3 * c1 * d0 * e5 - a4 * b3 * c1 * d5 * e0 -
                   a4 * b3 * c5 * d0 * e1 + a4 * b3 * c5 * d1 * e0 + a4 * b5 * c0 * d1 * e3 - a4 * b5 * c0 * d3 * e1 -
                   a4 * b5 * c1 * d0 * e3 + a4 * b5 * c1 * d3 * e0 + a4 * b5 * c3 * d0 * e1 - a4 * b5 * c3 * d1 * e0 +
                   a5 * b0 * c1 * d3 * e4 - a5 * b0 * c1 * d4 * e3 - a5 * b0 * c3 * d1 * e4 + a5 * b0 * c3 * d4 * e1 +
                   a5 * b0 * c4 * d1 * e3 - a5 * b0 * c4 * d3 * e1 - a5 * b1 * c0 * d3 * e4 + a5 * b1 * c0 * d4 * e3 +
                   a5 * b1 * c3 * d0 * e4 - a5 * b1 * c3 * d4 * e0 - a5 * b1 * c4 * d0 * e3 + a5 * b1 * c4 * d3 * e0 +
                   a5 * b3 * c0 * d1 * e4 - a5 * b3 * c0 * d4 * e1 - a5 * b3 * c1 * d0 * e4 + a5 * b3 * c1 * d4 * e0 +
                   a5 * b3 * c4 * d0 * e1 - a5 * b3 * c4 * d1 * e0 - a5 * b4 * c0 * d1 * e3 + a5 * b4 * c0 * d3 * e1 +
                   a5 * b4 * c1 * d0 * e3 - a5 * b4 * c1 * d3 * e0 - a5 * b4 * c3 * d0 * e1 + a5 * b4 * c3 * d1 * e0);
        T[0][3] = -t363 *
                  (b0 * c1 * d2 * e4 * f5 - b0 * c1 * d2 * e5 * f4 - b0 * c1 * d4 * e2 * f5 + b0 * c1 * d4 * e5 * f2 +
                   b0 * c1 * d5 * e2 * f4 - b0 * c1 * d5 * e4 * f2 - b0 * c2 * d1 * e4 * f5 + b0 * c2 * d1 * e5 * f4 +
                   b0 * c2 * d4 * e1 * f5 - b0 * c2 * d4 * e5 * f1 - b0 * c2 * d5 * e1 * f4 + b0 * c2 * d5 * e4 * f1 +
                   b0 * c4 * d1 * e2 * f5 - b0 * c4 * d1 * e5 * f2 - b0 * c4 * d2 * e1 * f5 + b0 * c4 * d2 * e5 * f1 +
                   b0 * c4 * d5 * e1 * f2 - b0 * c4 * d5 * e2 * f1 - b0 * c5 * d1 * e2 * f4 + b0 * c5 * d1 * e4 * f2 +
                   b0 * c5 * d2 * e1 * f4 - b0 * c5 * d2 * e4 * f1 - b0 * c5 * d4 * e1 * f2 + b0 * c5 * d4 * e2 * f1 -
                   b1 * c0 * d2 * e4 * f5 + b1 * c0 * d2 * e5 * f4 + b1 * c0 * d4 * e2 * f5 - b1 * c0 * d4 * e5 * f2 -
                   b1 * c0 * d5 * e2 * f4 + b1 * c0 * d5 * e4 * f2 + b1 * c2 * d0 * e4 * f5 - b1 * c2 * d0 * e5 * f4 -
                   b1 * c2 * d4 * e0 * f5 + b1 * c2 * d4 * e5 * f0 + b1 * c2 * d5 * e0 * f4 - b1 * c2 * d5 * e4 * f0 -
                   b1 * c4 * d0 * e2 * f5 + b1 * c4 * d0 * e5 * f2 + b1 * c4 * d2 * e0 * f5 - b1 * c4 * d2 * e5 * f0 -
                   b1 * c4 * d5 * e0 * f2 + b1 * c4 * d5 * e2 * f0 + b1 * c5 * d0 * e2 * f4 - b1 * c5 * d0 * e4 * f2 -
                   b1 * c5 * d2 * e0 * f4 + b1 * c5 * d2 * e4 * f0 + b1 * c5 * d4 * e0 * f2 - b1 * c5 * d4 * e2 * f0 +
                   b2 * c0 * d1 * e4 * f5 - b2 * c0 * d1 * e5 * f4 - b2 * c0 * d4 * e1 * f5 + b2 * c0 * d4 * e5 * f1 +
                   b2 * c0 * d5 * e1 * f4 - b2 * c0 * d5 * e4 * f1 - b2 * c1 * d0 * e4 * f5 + b2 * c1 * d0 * e5 * f4 +
                   b2 * c1 * d4 * e0 * f5 - b2 * c1 * d4 * e5 * f0 - b2 * c1 * d5 * e0 * f4 + b2 * c1 * d5 * e4 * f0 +
                   b2 * c4 * d0 * e1 * f5 - b2 * c4 * d0 * e5 * f1 - b2 * c4 * d1 * e0 * f5 + b2 * c4 * d1 * e5 * f0 +
                   b2 * c4 * d5 * e0 * f1 - b2 * c4 * d5 * e1 * f0 - b2 * c5 * d0 * e1 * f4 + b2 * c5 * d0 * e4 * f1 +
                   b2 * c5 * d1 * e0 * f4 - b2 * c5 * d1 * e4 * f0 - b2 * c5 * d4 * e0 * f1 + b2 * c5 * d4 * e1 * f0 -
                   b4 * c0 * d1 * e2 * f5 + b4 * c0 * d1 * e5 * f2 + b4 * c0 * d2 * e1 * f5 - b4 * c0 * d2 * e5 * f1 -
                   b4 * c0 * d5 * e1 * f2 + b4 * c0 * d5 * e2 * f1 + b4 * c1 * d0 * e2 * f5 - b4 * c1 * d0 * e5 * f2 -
                   b4 * c1 * d2 * e0 * f5 + b4 * c1 * d2 * e5 * f0 + b4 * c1 * d5 * e0 * f2 - b4 * c1 * d5 * e2 * f0 -
                   b4 * c2 * d0 * e1 * f5 + b4 * c2 * d0 * e5 * f1 + b4 * c2 * d1 * e0 * f5 - b4 * c2 * d1 * e5 * f0 -
                   b4 * c2 * d5 * e0 * f1 + b4 * c2 * d5 * e1 * f0 + b4 * c5 * d0 * e1 * f2 - b4 * c5 * d0 * e2 * f1 -
                   b4 * c5 * d1 * e0 * f2 + b4 * c5 * d1 * e2 * f0 + b4 * c5 * d2 * e0 * f1 - b4 * c5 * d2 * e1 * f0 +
                   b5 * c0 * d1 * e2 * f4 - b5 * c0 * d1 * e4 * f2 - b5 * c0 * d2 * e1 * f4 + b5 * c0 * d2 * e4 * f1 +
                   b5 * c0 * d4 * e1 * f2 - b5 * c0 * d4 * e2 * f1 - b5 * c1 * d0 * e2 * f4 + b5 * c1 * d0 * e4 * f2 +
                   b5 * c1 * d2 * e0 * f4 - b5 * c1 * d2 * e4 * f0 - b5 * c1 * d4 * e0 * f2 + b5 * c1 * d4 * e2 * f0 +
                   b5 * c2 * d0 * e1 * f4 - b5 * c2 * d0 * e4 * f1 - b5 * c2 * d1 * e0 * f4 + b5 * c2 * d1 * e4 * f0 +
                   b5 * c2 * d4 * e0 * f1 - b5 * c2 * d4 * e1 * f0 - b5 * c4 * d0 * e1 * f2 + b5 * c4 * d0 * e2 * f1 +
                   b5 * c4 * d1 * e0 * f2 - b5 * c4 * d1 * e2 * f0 - b5 * c4 * d2 * e0 * f1 + b5 * c4 * d2 * e1 * f0);
        T[1][3] = t363 *
                  (a0 * c1 * d2 * e4 * f5 - a0 * c1 * d2 * e5 * f4 - a0 * c1 * d4 * e2 * f5 + a0 * c1 * d4 * e5 * f2 +
                   a0 * c1 * d5 * e2 * f4 - a0 * c1 * d5 * e4 * f2 - a0 * c2 * d1 * e4 * f5 + a0 * c2 * d1 * e5 * f4 +
                   a0 * c2 * d4 * e1 * f5 - a0 * c2 * d4 * e5 * f1 - a0 * c2 * d5 * e1 * f4 + a0 * c2 * d5 * e4 * f1 +
                   a0 * c4 * d1 * e2 * f5 - a0 * c4 * d1 * e5 * f2 - a0 * c4 * d2 * e1 * f5 + a0 * c4 * d2 * e5 * f1 +
                   a0 * c4 * d5 * e1 * f2 - a0 * c4 * d5 * e2 * f1 - a0 * c5 * d1 * e2 * f4 + a0 * c5 * d1 * e4 * f2 +
                   a0 * c5 * d2 * e1 * f4 - a0 * c5 * d2 * e4 * f1 - a0 * c5 * d4 * e1 * f2 + a0 * c5 * d4 * e2 * f1 -
                   a1 * c0 * d2 * e4 * f5 + a1 * c0 * d2 * e5 * f4 + a1 * c0 * d4 * e2 * f5 - a1 * c0 * d4 * e5 * f2 -
                   a1 * c0 * d5 * e2 * f4 + a1 * c0 * d5 * e4 * f2 + a1 * c2 * d0 * e4 * f5 - a1 * c2 * d0 * e5 * f4 -
                   a1 * c2 * d4 * e0 * f5 + a1 * c2 * d4 * e5 * f0 + a1 * c2 * d5 * e0 * f4 - a1 * c2 * d5 * e4 * f0 -
                   a1 * c4 * d0 * e2 * f5 + a1 * c4 * d0 * e5 * f2 + a1 * c4 * d2 * e0 * f5 - a1 * c4 * d2 * e5 * f0 -
                   a1 * c4 * d5 * e0 * f2 + a1 * c4 * d5 * e2 * f0 + a1 * c5 * d0 * e2 * f4 - a1 * c5 * d0 * e4 * f2 -
                   a1 * c5 * d2 * e0 * f4 + a1 * c5 * d2 * e4 * f0 + a1 * c5 * d4 * e0 * f2 - a1 * c5 * d4 * e2 * f0 +
                   a2 * c0 * d1 * e4 * f5 - a2 * c0 * d1 * e5 * f4 - a2 * c0 * d4 * e1 * f5 + a2 * c0 * d4 * e5 * f1 +
                   a2 * c0 * d5 * e1 * f4 - a2 * c0 * d5 * e4 * f1 - a2 * c1 * d0 * e4 * f5 + a2 * c1 * d0 * e5 * f4 +
                   a2 * c1 * d4 * e0 * f5 - a2 * c1 * d4 * e5 * f0 - a2 * c1 * d5 * e0 * f4 + a2 * c1 * d5 * e4 * f0 +
                   a2 * c4 * d0 * e1 * f5 - a2 * c4 * d0 * e5 * f1 - a2 * c4 * d1 * e0 * f5 + a2 * c4 * d1 * e5 * f0 +
                   a2 * c4 * d5 * e0 * f1 - a2 * c4 * d5 * e1 * f0 - a2 * c5 * d0 * e1 * f4 + a2 * c5 * d0 * e4 * f1 +
                   a2 * c5 * d1 * e0 * f4 - a2 * c5 * d1 * e4 * f0 - a2 * c5 * d4 * e0 * f1 + a2 * c5 * d4 * e1 * f0 -
                   a4 * c0 * d1 * e2 * f5 + a4 * c0 * d1 * e5 * f2 + a4 * c0 * d2 * e1 * f5 - a4 * c0 * d2 * e5 * f1 -
                   a4 * c0 * d5 * e1 * f2 + a4 * c0 * d5 * e2 * f1 + a4 * c1 * d0 * e2 * f5 - a4 * c1 * d0 * e5 * f2 -
                   a4 * c1 * d2 * e0 * f5 + a4 * c1 * d2 * e5 * f0 + a4 * c1 * d5 * e0 * f2 - a4 * c1 * d5 * e2 * f0 -
                   a4 * c2 * d0 * e1 * f5 + a4 * c2 * d0 * e5 * f1 + a4 * c2 * d1 * e0 * f5 - a4 * c2 * d1 * e5 * f0 -
                   a4 * c2 * d5 * e0 * f1 + a4 * c2 * d5 * e1 * f0 + a4 * c5 * d0 * e1 * f2 - a4 * c5 * d0 * e2 * f1 -
                   a4 * c5 * d1 * e0 * f2 + a4 * c5 * d1 * e2 * f0 + a4 * c5 * d2 * e0 * f1 - a4 * c5 * d2 * e1 * f0 +
                   a5 * c0 * d1 * e2 * f4 - a5 * c0 * d1 * e4 * f2 - a5 * c0 * d2 * e1 * f4 + a5 * c0 * d2 * e4 * f1 +
                   a5 * c0 * d4 * e1 * f2 - a5 * c0 * d4 * e2 * f1 - a5 * c1 * d0 * e2 * f4 + a5 * c1 * d0 * e4 * f2 +
                   a5 * c1 * d2 * e0 * f4 - a5 * c1 * d2 * e4 * f0 - a5 * c1 * d4 * e0 * f2 + a5 * c1 * d4 * e2 * f0 +
                   a5 * c2 * d0 * e1 * f4 - a5 * c2 * d0 * e4 * f1 - a5 * c2 * d1 * e0 * f4 + a5 * c2 * d1 * e4 * f0 +
                   a5 * c2 * d4 * e0 * f1 - a5 * c2 * d4 * e1 * f0 - a5 * c4 * d0 * e1 * f2 + a5 * c4 * d0 * e2 * f1 +
                   a5 * c4 * d1 * e0 * f2 - a5 * c4 * d1 * e2 * f0 - a5 * c4 * d2 * e0 * f1 + a5 * c4 * d2 * e1 * f0);
        T[2][3] = -t363 *
                  (a0 * b1 * d2 * e4 * f5 - a0 * b1 * d2 * e5 * f4 - a0 * b1 * d4 * e2 * f5 + a0 * b1 * d4 * e5 * f2 +
                   a0 * b1 * d5 * e2 * f4 - a0 * b1 * d5 * e4 * f2 - a0 * b2 * d1 * e4 * f5 + a0 * b2 * d1 * e5 * f4 +
                   a0 * b2 * d4 * e1 * f5 - a0 * b2 * d4 * e5 * f1 - a0 * b2 * d5 * e1 * f4 + a0 * b2 * d5 * e4 * f1 +
                   a0 * b4 * d1 * e2 * f5 - a0 * b4 * d1 * e5 * f2 - a0 * b4 * d2 * e1 * f5 + a0 * b4 * d2 * e5 * f1 +
                   a0 * b4 * d5 * e1 * f2 - a0 * b4 * d5 * e2 * f1 - a0 * b5 * d1 * e2 * f4 + a0 * b5 * d1 * e4 * f2 +
                   a0 * b5 * d2 * e1 * f4 - a0 * b5 * d2 * e4 * f1 - a0 * b5 * d4 * e1 * f2 + a0 * b5 * d4 * e2 * f1 -
                   a1 * b0 * d2 * e4 * f5 + a1 * b0 * d2 * e5 * f4 + a1 * b0 * d4 * e2 * f5 - a1 * b0 * d4 * e5 * f2 -
                   a1 * b0 * d5 * e2 * f4 + a1 * b0 * d5 * e4 * f2 + a1 * b2 * d0 * e4 * f5 - a1 * b2 * d0 * e5 * f4 -
                   a1 * b2 * d4 * e0 * f5 + a1 * b2 * d4 * e5 * f0 + a1 * b2 * d5 * e0 * f4 - a1 * b2 * d5 * e4 * f0 -
                   a1 * b4 * d0 * e2 * f5 + a1 * b4 * d0 * e5 * f2 + a1 * b4 * d2 * e0 * f5 - a1 * b4 * d2 * e5 * f0 -
                   a1 * b4 * d5 * e0 * f2 + a1 * b4 * d5 * e2 * f0 + a1 * b5 * d0 * e2 * f4 - a1 * b5 * d0 * e4 * f2 -
                   a1 * b5 * d2 * e0 * f4 + a1 * b5 * d2 * e4 * f0 + a1 * b5 * d4 * e0 * f2 - a1 * b5 * d4 * e2 * f0 +
                   a2 * b0 * d1 * e4 * f5 - a2 * b0 * d1 * e5 * f4 - a2 * b0 * d4 * e1 * f5 + a2 * b0 * d4 * e5 * f1 +
                   a2 * b0 * d5 * e1 * f4 - a2 * b0 * d5 * e4 * f1 - a2 * b1 * d0 * e4 * f5 + a2 * b1 * d0 * e5 * f4 +
                   a2 * b1 * d4 * e0 * f5 - a2 * b1 * d4 * e5 * f0 - a2 * b1 * d5 * e0 * f4 + a2 * b1 * d5 * e4 * f0 +
                   a2 * b4 * d0 * e1 * f5 - a2 * b4 * d0 * e5 * f1 - a2 * b4 * d1 * e0 * f5 + a2 * b4 * d1 * e5 * f0 +
                   a2 * b4 * d5 * e0 * f1 - a2 * b4 * d5 * e1 * f0 - a2 * b5 * d0 * e1 * f4 + a2 * b5 * d0 * e4 * f1 +
                   a2 * b5 * d1 * e0 * f4 - a2 * b5 * d1 * e4 * f0 - a2 * b5 * d4 * e0 * f1 + a2 * b5 * d4 * e1 * f0 -
                   a4 * b0 * d1 * e2 * f5 + a4 * b0 * d1 * e5 * f2 + a4 * b0 * d2 * e1 * f5 - a4 * b0 * d2 * e5 * f1 -
                   a4 * b0 * d5 * e1 * f2 + a4 * b0 * d5 * e2 * f1 + a4 * b1 * d0 * e2 * f5 - a4 * b1 * d0 * e5 * f2 -
                   a4 * b1 * d2 * e0 * f5 + a4 * b1 * d2 * e5 * f0 + a4 * b1 * d5 * e0 * f2 - a4 * b1 * d5 * e2 * f0 -
                   a4 * b2 * d0 * e1 * f5 + a4 * b2 * d0 * e5 * f1 + a4 * b2 * d1 * e0 * f5 - a4 * b2 * d1 * e5 * f0 -
                   a4 * b2 * d5 * e0 * f1 + a4 * b2 * d5 * e1 * f0 + a4 * b5 * d0 * e1 * f2 - a4 * b5 * d0 * e2 * f1 -
                   a4 * b5 * d1 * e0 * f2 + a4 * b5 * d1 * e2 * f0 + a4 * b5 * d2 * e0 * f1 - a4 * b5 * d2 * e1 * f0 +
                   a5 * b0 * d1 * e2 * f4 - a5 * b0 * d1 * e4 * f2 - a5 * b0 * d2 * e1 * f4 + a5 * b0 * d2 * e4 * f1 +
                   a5 * b0 * d4 * e1 * f2 - a5 * b0 * d4 * e2 * f1 - a5 * b1 * d0 * e2 * f4 + a5 * b1 * d0 * e4 * f2 +
                   a5 * b1 * d2 * e0 * f4 - a5 * b1 * d2 * e4 * f0 - a5 * b1 * d4 * e0 * f2 + a5 * b1 * d4 * e2 * f0 +
                   a5 * b2 * d0 * e1 * f4 - a5 * b2 * d0 * e4 * f1 - a5 * b2 * d1 * e0 * f4 + a5 * b2 * d1 * e4 * f0 +
                   a5 * b2 * d4 * e0 * f1 - a5 * b2 * d4 * e1 * f0 - a5 * b4 * d0 * e1 * f2 + a5 * b4 * d0 * e2 * f1 +
                   a5 * b4 * d1 * e0 * f2 - a5 * b4 * d1 * e2 * f0 - a5 * b4 * d2 * e0 * f1 + a5 * b4 * d2 * e1 * f0);
        T[3][3] = t363 *
                  (a0 * b1 * c2 * e4 * f5 - a0 * b1 * c2 * e5 * f4 - a0 * b1 * c4 * e2 * f5 + a0 * b1 * c4 * e5 * f2 +
                   a0 * b1 * c5 * e2 * f4 - a0 * b1 * c5 * e4 * f2 - a0 * b2 * c1 * e4 * f5 + a0 * b2 * c1 * e5 * f4 +
                   a0 * b2 * c4 * e1 * f5 - a0 * b2 * c4 * e5 * f1 - a0 * b2 * c5 * e1 * f4 + a0 * b2 * c5 * e4 * f1 +
                   a0 * b4 * c1 * e2 * f5 - a0 * b4 * c1 * e5 * f2 - a0 * b4 * c2 * e1 * f5 + a0 * b4 * c2 * e5 * f1 +
                   a0 * b4 * c5 * e1 * f2 - a0 * b4 * c5 * e2 * f1 - a0 * b5 * c1 * e2 * f4 + a0 * b5 * c1 * e4 * f2 +
                   a0 * b5 * c2 * e1 * f4 - a0 * b5 * c2 * e4 * f1 - a0 * b5 * c4 * e1 * f2 + a0 * b5 * c4 * e2 * f1 -
                   a1 * b0 * c2 * e4 * f5 + a1 * b0 * c2 * e5 * f4 + a1 * b0 * c4 * e2 * f5 - a1 * b0 * c4 * e5 * f2 -
                   a1 * b0 * c5 * e2 * f4 + a1 * b0 * c5 * e4 * f2 + a1 * b2 * c0 * e4 * f5 - a1 * b2 * c0 * e5 * f4 -
                   a1 * b2 * c4 * e0 * f5 + a1 * b2 * c4 * e5 * f0 + a1 * b2 * c5 * e0 * f4 - a1 * b2 * c5 * e4 * f0 -
                   a1 * b4 * c0 * e2 * f5 + a1 * b4 * c0 * e5 * f2 + a1 * b4 * c2 * e0 * f5 - a1 * b4 * c2 * e5 * f0 -
                   a1 * b4 * c5 * e0 * f2 + a1 * b4 * c5 * e2 * f0 + a1 * b5 * c0 * e2 * f4 - a1 * b5 * c0 * e4 * f2 -
                   a1 * b5 * c2 * e0 * f4 + a1 * b5 * c2 * e4 * f0 + a1 * b5 * c4 * e0 * f2 - a1 * b5 * c4 * e2 * f0 +
                   a2 * b0 * c1 * e4 * f5 - a2 * b0 * c1 * e5 * f4 - a2 * b0 * c4 * e1 * f5 + a2 * b0 * c4 * e5 * f1 +
                   a2 * b0 * c5 * e1 * f4 - a2 * b0 * c5 * e4 * f1 - a2 * b1 * c0 * e4 * f5 + a2 * b1 * c0 * e5 * f4 +
                   a2 * b1 * c4 * e0 * f5 - a2 * b1 * c4 * e5 * f0 - a2 * b1 * c5 * e0 * f4 + a2 * b1 * c5 * e4 * f0 +
                   a2 * b4 * c0 * e1 * f5 - a2 * b4 * c0 * e5 * f1 - a2 * b4 * c1 * e0 * f5 + a2 * b4 * c1 * e5 * f0 +
                   a2 * b4 * c5 * e0 * f1 - a2 * b4 * c5 * e1 * f0 - a2 * b5 * c0 * e1 * f4 + a2 * b5 * c0 * e4 * f1 +
                   a2 * b5 * c1 * e0 * f4 - a2 * b5 * c1 * e4 * f0 - a2 * b5 * c4 * e0 * f1 + a2 * b5 * c4 * e1 * f0 -
                   a4 * b0 * c1 * e2 * f5 + a4 * b0 * c1 * e5 * f2 + a4 * b0 * c2 * e1 * f5 - a4 * b0 * c2 * e5 * f1 -
                   a4 * b0 * c5 * e1 * f2 + a4 * b0 * c5 * e2 * f1 + a4 * b1 * c0 * e2 * f5 - a4 * b1 * c0 * e5 * f2 -
                   a4 * b1 * c2 * e0 * f5 + a4 * b1 * c2 * e5 * f0 + a4 * b1 * c5 * e0 * f2 - a4 * b1 * c5 * e2 * f0 -
                   a4 * b2 * c0 * e1 * f5 + a4 * b2 * c0 * e5 * f1 + a4 * b2 * c1 * e0 * f5 - a4 * b2 * c1 * e5 * f0 -
                   a4 * b2 * c5 * e0 * f1 + a4 * b2 * c5 * e1 * f0 + a4 * b5 * c0 * e1 * f2 - a4 * b5 * c0 * e2 * f1 -
                   a4 * b5 * c1 * e0 * f2 + a4 * b5 * c1 * e2 * f0 + a4 * b5 * c2 * e0 * f1 - a4 * b5 * c2 * e1 * f0 +
                   a5 * b0 * c1 * e2 * f4 - a5 * b0 * c1 * e4 * f2 - a5 * b0 * c2 * e1 * f4 + a5 * b0 * c2 * e4 * f1 +
                   a5 * b0 * c4 * e1 * f2 - a5 * b0 * c4 * e2 * f1 - a5 * b1 * c0 * e2 * f4 + a5 * b1 * c0 * e4 * f2 +
                   a5 * b1 * c2 * e0 * f4 - a5 * b1 * c2 * e4 * f0 - a5 * b1 * c4 * e0 * f2 + a5 * b1 * c4 * e2 * f0 +
                   a5 * b2 * c0 * e1 * f4 - a5 * b2 * c0 * e4 * f1 - a5 * b2 * c1 * e0 * f4 + a5 * b2 * c1 * e4 * f0 +
                   a5 * b2 * c4 * e0 * f1 - a5 * b2 * c4 * e1 * f0 - a5 * b4 * c0 * e1 * f2 + a5 * b4 * c0 * e2 * f1 +
                   a5 * b4 * c1 * e0 * f2 - a5 * b4 * c1 * e2 * f0 - a5 * b4 * c2 * e0 * f1 + a5 * b4 * c2 * e1 * f0);
        T[4][3] = -t363 *
                  (a0 * b1 * c2 * d4 * f5 - a0 * b1 * c2 * d5 * f4 - a0 * b1 * c4 * d2 * f5 + a0 * b1 * c4 * d5 * f2 +
                   a0 * b1 * c5 * d2 * f4 - a0 * b1 * c5 * d4 * f2 - a0 * b2 * c1 * d4 * f5 + a0 * b2 * c1 * d5 * f4 +
                   a0 * b2 * c4 * d1 * f5 - a0 * b2 * c4 * d5 * f1 - a0 * b2 * c5 * d1 * f4 + a0 * b2 * c5 * d4 * f1 +
                   a0 * b4 * c1 * d2 * f5 - a0 * b4 * c1 * d5 * f2 - a0 * b4 * c2 * d1 * f5 + a0 * b4 * c2 * d5 * f1 +
                   a0 * b4 * c5 * d1 * f2 - a0 * b4 * c5 * d2 * f1 - a0 * b5 * c1 * d2 * f4 + a0 * b5 * c1 * d4 * f2 +
                   a0 * b5 * c2 * d1 * f4 - a0 * b5 * c2 * d4 * f1 - a0 * b5 * c4 * d1 * f2 + a0 * b5 * c4 * d2 * f1 -
                   a1 * b0 * c2 * d4 * f5 + a1 * b0 * c2 * d5 * f4 + a1 * b0 * c4 * d2 * f5 - a1 * b0 * c4 * d5 * f2 -
                   a1 * b0 * c5 * d2 * f4 + a1 * b0 * c5 * d4 * f2 + a1 * b2 * c0 * d4 * f5 - a1 * b2 * c0 * d5 * f4 -
                   a1 * b2 * c4 * d0 * f5 + a1 * b2 * c4 * d5 * f0 + a1 * b2 * c5 * d0 * f4 - a1 * b2 * c5 * d4 * f0 -
                   a1 * b4 * c0 * d2 * f5 + a1 * b4 * c0 * d5 * f2 + a1 * b4 * c2 * d0 * f5 - a1 * b4 * c2 * d5 * f0 -
                   a1 * b4 * c5 * d0 * f2 + a1 * b4 * c5 * d2 * f0 + a1 * b5 * c0 * d2 * f4 - a1 * b5 * c0 * d4 * f2 -
                   a1 * b5 * c2 * d0 * f4 + a1 * b5 * c2 * d4 * f0 + a1 * b5 * c4 * d0 * f2 - a1 * b5 * c4 * d2 * f0 +
                   a2 * b0 * c1 * d4 * f5 - a2 * b0 * c1 * d5 * f4 - a2 * b0 * c4 * d1 * f5 + a2 * b0 * c4 * d5 * f1 +
                   a2 * b0 * c5 * d1 * f4 - a2 * b0 * c5 * d4 * f1 - a2 * b1 * c0 * d4 * f5 + a2 * b1 * c0 * d5 * f4 +
                   a2 * b1 * c4 * d0 * f5 - a2 * b1 * c4 * d5 * f0 - a2 * b1 * c5 * d0 * f4 + a2 * b1 * c5 * d4 * f0 +
                   a2 * b4 * c0 * d1 * f5 - a2 * b4 * c0 * d5 * f1 - a2 * b4 * c1 * d0 * f5 + a2 * b4 * c1 * d5 * f0 +
                   a2 * b4 * c5 * d0 * f1 - a2 * b4 * c5 * d1 * f0 - a2 * b5 * c0 * d1 * f4 + a2 * b5 * c0 * d4 * f1 +
                   a2 * b5 * c1 * d0 * f4 - a2 * b5 * c1 * d4 * f0 - a2 * b5 * c4 * d0 * f1 + a2 * b5 * c4 * d1 * f0 -
                   a4 * b0 * c1 * d2 * f5 + a4 * b0 * c1 * d5 * f2 + a4 * b0 * c2 * d1 * f5 - a4 * b0 * c2 * d5 * f1 -
                   a4 * b0 * c5 * d1 * f2 + a4 * b0 * c5 * d2 * f1 + a4 * b1 * c0 * d2 * f5 - a4 * b1 * c0 * d5 * f2 -
                   a4 * b1 * c2 * d0 * f5 + a4 * b1 * c2 * d5 * f0 + a4 * b1 * c5 * d0 * f2 - a4 * b1 * c5 * d2 * f0 -
                   a4 * b2 * c0 * d1 * f5 + a4 * b2 * c0 * d5 * f1 + a4 * b2 * c1 * d0 * f5 - a4 * b2 * c1 * d5 * f0 -
                   a4 * b2 * c5 * d0 * f1 + a4 * b2 * c5 * d1 * f0 + a4 * b5 * c0 * d1 * f2 - a4 * b5 * c0 * d2 * f1 -
                   a4 * b5 * c1 * d0 * f2 + a4 * b5 * c1 * d2 * f0 + a4 * b5 * c2 * d0 * f1 - a4 * b5 * c2 * d1 * f0 +
                   a5 * b0 * c1 * d2 * f4 - a5 * b0 * c1 * d4 * f2 - a5 * b0 * c2 * d1 * f4 + a5 * b0 * c2 * d4 * f1 +
                   a5 * b0 * c4 * d1 * f2 - a5 * b0 * c4 * d2 * f1 - a5 * b1 * c0 * d2 * f4 + a5 * b1 * c0 * d4 * f2 +
                   a5 * b1 * c2 * d0 * f4 - a5 * b1 * c2 * d4 * f0 - a5 * b1 * c4 * d0 * f2 + a5 * b1 * c4 * d2 * f0 +
                   a5 * b2 * c0 * d1 * f4 - a5 * b2 * c0 * d4 * f1 - a5 * b2 * c1 * d0 * f4 + a5 * b2 * c1 * d4 * f0 +
                   a5 * b2 * c4 * d0 * f1 - a5 * b2 * c4 * d1 * f0 - a5 * b4 * c0 * d1 * f2 + a5 * b4 * c0 * d2 * f1 +
                   a5 * b4 * c1 * d0 * f2 - a5 * b4 * c1 * d2 * f0 - a5 * b4 * c2 * d0 * f1 + a5 * b4 * c2 * d1 * f0);
        T[5][3] = t363 *
                  (a0 * b1 * c2 * d4 * e5 - a0 * b1 * c2 * d5 * e4 - a0 * b1 * c4 * d2 * e5 + a0 * b1 * c4 * d5 * e2 +
                   a0 * b1 * c5 * d2 * e4 - a0 * b1 * c5 * d4 * e2 - a0 * b2 * c1 * d4 * e5 + a0 * b2 * c1 * d5 * e4 +
                   a0 * b2 * c4 * d1 * e5 - a0 * b2 * c4 * d5 * e1 - a0 * b2 * c5 * d1 * e4 + a0 * b2 * c5 * d4 * e1 +
                   a0 * b4 * c1 * d2 * e5 - a0 * b4 * c1 * d5 * e2 - a0 * b4 * c2 * d1 * e5 + a0 * b4 * c2 * d5 * e1 +
                   a0 * b4 * c5 * d1 * e2 - a0 * b4 * c5 * d2 * e1 - a0 * b5 * c1 * d2 * e4 + a0 * b5 * c1 * d4 * e2 +
                   a0 * b5 * c2 * d1 * e4 - a0 * b5 * c2 * d4 * e1 - a0 * b5 * c4 * d1 * e2 + a0 * b5 * c4 * d2 * e1 -
                   a1 * b0 * c2 * d4 * e5 + a1 * b0 * c2 * d5 * e4 + a1 * b0 * c4 * d2 * e5 - a1 * b0 * c4 * d5 * e2 -
                   a1 * b0 * c5 * d2 * e4 + a1 * b0 * c5 * d4 * e2 + a1 * b2 * c0 * d4 * e5 - a1 * b2 * c0 * d5 * e4 -
                   a1 * b2 * c4 * d0 * e5 + a1 * b2 * c4 * d5 * e0 + a1 * b2 * c5 * d0 * e4 - a1 * b2 * c5 * d4 * e0 -
                   a1 * b4 * c0 * d2 * e5 + a1 * b4 * c0 * d5 * e2 + a1 * b4 * c2 * d0 * e5 - a1 * b4 * c2 * d5 * e0 -
                   a1 * b4 * c5 * d0 * e2 + a1 * b4 * c5 * d2 * e0 + a1 * b5 * c0 * d2 * e4 - a1 * b5 * c0 * d4 * e2 -
                   a1 * b5 * c2 * d0 * e4 + a1 * b5 * c2 * d4 * e0 + a1 * b5 * c4 * d0 * e2 - a1 * b5 * c4 * d2 * e0 +
                   a2 * b0 * c1 * d4 * e5 - a2 * b0 * c1 * d5 * e4 - a2 * b0 * c4 * d1 * e5 + a2 * b0 * c4 * d5 * e1 +
                   a2 * b0 * c5 * d1 * e4 - a2 * b0 * c5 * d4 * e1 - a2 * b1 * c0 * d4 * e5 + a2 * b1 * c0 * d5 * e4 +
                   a2 * b1 * c4 * d0 * e5 - a2 * b1 * c4 * d5 * e0 - a2 * b1 * c5 * d0 * e4 + a2 * b1 * c5 * d4 * e0 +
                   a2 * b4 * c0 * d1 * e5 - a2 * b4 * c0 * d5 * e1 - a2 * b4 * c1 * d0 * e5 + a2 * b4 * c1 * d5 * e0 +
                   a2 * b4 * c5 * d0 * e1 - a2 * b4 * c5 * d1 * e0 - a2 * b5 * c0 * d1 * e4 + a2 * b5 * c0 * d4 * e1 +
                   a2 * b5 * c1 * d0 * e4 - a2 * b5 * c1 * d4 * e0 - a2 * b5 * c4 * d0 * e1 + a2 * b5 * c4 * d1 * e0 -
                   a4 * b0 * c1 * d2 * e5 + a4 * b0 * c1 * d5 * e2 + a4 * b0 * c2 * d1 * e5 - a4 * b0 * c2 * d5 * e1 -
                   a4 * b0 * c5 * d1 * e2 + a4 * b0 * c5 * d2 * e1 + a4 * b1 * c0 * d2 * e5 - a4 * b1 * c0 * d5 * e2 -
                   a4 * b1 * c2 * d0 * e5 + a4 * b1 * c2 * d5 * e0 + a4 * b1 * c5 * d0 * e2 - a4 * b1 * c5 * d2 * e0 -
                   a4 * b2 * c0 * d1 * e5 + a4 * b2 * c0 * d5 * e1 + a4 * b2 * c1 * d0 * e5 - a4 * b2 * c1 * d5 * e0 -
                   a4 * b2 * c5 * d0 * e1 + a4 * b2 * c5 * d1 * e0 + a4 * b5 * c0 * d1 * e2 - a4 * b5 * c0 * d2 * e1 -
                   a4 * b5 * c1 * d0 * e2 + a4 * b5 * c1 * d2 * e0 + a4 * b5 * c2 * d0 * e1 - a4 * b5 * c2 * d1 * e0 +
                   a5 * b0 * c1 * d2 * e4 - a5 * b0 * c1 * d4 * e2 - a5 * b0 * c2 * d1 * e4 + a5 * b0 * c2 * d4 * e1 +
                   a5 * b0 * c4 * d1 * e2 - a5 * b0 * c4 * d2 * e1 - a5 * b1 * c0 * d2 * e4 + a5 * b1 * c0 * d4 * e2 +
                   a5 * b1 * c2 * d0 * e4 - a5 * b1 * c2 * d4 * e0 - a5 * b1 * c4 * d0 * e2 + a5 * b1 * c4 * d2 * e0 +
                   a5 * b2 * c0 * d1 * e4 - a5 * b2 * c0 * d4 * e1 - a5 * b2 * c1 * d0 * e4 + a5 * b2 * c1 * d4 * e0 +
                   a5 * b2 * c4 * d0 * e1 - a5 * b2 * c4 * d1 * e0 - a5 * b4 * c0 * d1 * e2 + a5 * b4 * c0 * d2 * e1 +
                   a5 * b4 * c1 * d0 * e2 - a5 * b4 * c1 * d2 * e0 - a5 * b4 * c2 * d0 * e1 + a5 * b4 * c2 * d1 * e0);
        T[0][4] = t363 *
                  (b0 * c1 * d2 * e3 * f5 - b0 * c1 * d2 * e5 * f3 - b0 * c1 * d3 * e2 * f5 + b0 * c1 * d3 * e5 * f2 +
                   b0 * c1 * d5 * e2 * f3 - b0 * c1 * d5 * e3 * f2 - b0 * c2 * d1 * e3 * f5 + b0 * c2 * d1 * e5 * f3 +
                   b0 * c2 * d3 * e1 * f5 - b0 * c2 * d3 * e5 * f1 - b0 * c2 * d5 * e1 * f3 + b0 * c2 * d5 * e3 * f1 +
                   b0 * c3 * d1 * e2 * f5 - b0 * c3 * d1 * e5 * f2 - b0 * c3 * d2 * e1 * f5 + b0 * c3 * d2 * e5 * f1 +
                   b0 * c3 * d5 * e1 * f2 - b0 * c3 * d5 * e2 * f1 - b0 * c5 * d1 * e2 * f3 + b0 * c5 * d1 * e3 * f2 +
                   b0 * c5 * d2 * e1 * f3 - b0 * c5 * d2 * e3 * f1 - b0 * c5 * d3 * e1 * f2 + b0 * c5 * d3 * e2 * f1 -
                   b1 * c0 * d2 * e3 * f5 + b1 * c0 * d2 * e5 * f3 + b1 * c0 * d3 * e2 * f5 - b1 * c0 * d3 * e5 * f2 -
                   b1 * c0 * d5 * e2 * f3 + b1 * c0 * d5 * e3 * f2 + b1 * c2 * d0 * e3 * f5 - b1 * c2 * d0 * e5 * f3 -
                   b1 * c2 * d3 * e0 * f5 + b1 * c2 * d3 * e5 * f0 + b1 * c2 * d5 * e0 * f3 - b1 * c2 * d5 * e3 * f0 -
                   b1 * c3 * d0 * e2 * f5 + b1 * c3 * d0 * e5 * f2 + b1 * c3 * d2 * e0 * f5 - b1 * c3 * d2 * e5 * f0 -
                   b1 * c3 * d5 * e0 * f2 + b1 * c3 * d5 * e2 * f0 + b1 * c5 * d0 * e2 * f3 - b1 * c5 * d0 * e3 * f2 -
                   b1 * c5 * d2 * e0 * f3 + b1 * c5 * d2 * e3 * f0 + b1 * c5 * d3 * e0 * f2 - b1 * c5 * d3 * e2 * f0 +
                   b2 * c0 * d1 * e3 * f5 - b2 * c0 * d1 * e5 * f3 - b2 * c0 * d3 * e1 * f5 + b2 * c0 * d3 * e5 * f1 +
                   b2 * c0 * d5 * e1 * f3 - b2 * c0 * d5 * e3 * f1 - b2 * c1 * d0 * e3 * f5 + b2 * c1 * d0 * e5 * f3 +
                   b2 * c1 * d3 * e0 * f5 - b2 * c1 * d3 * e5 * f0 - b2 * c1 * d5 * e0 * f3 + b2 * c1 * d5 * e3 * f0 +
                   b2 * c3 * d0 * e1 * f5 - b2 * c3 * d0 * e5 * f1 - b2 * c3 * d1 * e0 * f5 + b2 * c3 * d1 * e5 * f0 +
                   b2 * c3 * d5 * e0 * f1 - b2 * c3 * d5 * e1 * f0 - b2 * c5 * d0 * e1 * f3 + b2 * c5 * d0 * e3 * f1 +
                   b2 * c5 * d1 * e0 * f3 - b2 * c5 * d1 * e3 * f0 - b2 * c5 * d3 * e0 * f1 + b2 * c5 * d3 * e1 * f0 -
                   b3 * c0 * d1 * e2 * f5 + b3 * c0 * d1 * e5 * f2 + b3 * c0 * d2 * e1 * f5 - b3 * c0 * d2 * e5 * f1 -
                   b3 * c0 * d5 * e1 * f2 + b3 * c0 * d5 * e2 * f1 + b3 * c1 * d0 * e2 * f5 - b3 * c1 * d0 * e5 * f2 -
                   b3 * c1 * d2 * e0 * f5 + b3 * c1 * d2 * e5 * f0 + b3 * c1 * d5 * e0 * f2 - b3 * c1 * d5 * e2 * f0 -
                   b3 * c2 * d0 * e1 * f5 + b3 * c2 * d0 * e5 * f1 + b3 * c2 * d1 * e0 * f5 - b3 * c2 * d1 * e5 * f0 -
                   b3 * c2 * d5 * e0 * f1 + b3 * c2 * d5 * e1 * f0 + b3 * c5 * d0 * e1 * f2 - b3 * c5 * d0 * e2 * f1 -
                   b3 * c5 * d1 * e0 * f2 + b3 * c5 * d1 * e2 * f0 + b3 * c5 * d2 * e0 * f1 - b3 * c5 * d2 * e1 * f0 +
                   b5 * c0 * d1 * e2 * f3 - b5 * c0 * d1 * e3 * f2 - b5 * c0 * d2 * e1 * f3 + b5 * c0 * d2 * e3 * f1 +
                   b5 * c0 * d3 * e1 * f2 - b5 * c0 * d3 * e2 * f1 - b5 * c1 * d0 * e2 * f3 + b5 * c1 * d0 * e3 * f2 +
                   b5 * c1 * d2 * e0 * f3 - b5 * c1 * d2 * e3 * f0 - b5 * c1 * d3 * e0 * f2 + b5 * c1 * d3 * e2 * f0 +
                   b5 * c2 * d0 * e1 * f3 - b5 * c2 * d0 * e3 * f1 - b5 * c2 * d1 * e0 * f3 + b5 * c2 * d1 * e3 * f0 +
                   b5 * c2 * d3 * e0 * f1 - b5 * c2 * d3 * e1 * f0 - b5 * c3 * d0 * e1 * f2 + b5 * c3 * d0 * e2 * f1 +
                   b5 * c3 * d1 * e0 * f2 - b5 * c3 * d1 * e2 * f0 - b5 * c3 * d2 * e0 * f1 + b5 * c3 * d2 * e1 * f0);
        T[1][4] = -t363 *
                  (a0 * c1 * d2 * e3 * f5 - a0 * c1 * d2 * e5 * f3 - a0 * c1 * d3 * e2 * f5 + a0 * c1 * d3 * e5 * f2 +
                   a0 * c1 * d5 * e2 * f3 - a0 * c1 * d5 * e3 * f2 - a0 * c2 * d1 * e3 * f5 + a0 * c2 * d1 * e5 * f3 +
                   a0 * c2 * d3 * e1 * f5 - a0 * c2 * d3 * e5 * f1 - a0 * c2 * d5 * e1 * f3 + a0 * c2 * d5 * e3 * f1 +
                   a0 * c3 * d1 * e2 * f5 - a0 * c3 * d1 * e5 * f2 - a0 * c3 * d2 * e1 * f5 + a0 * c3 * d2 * e5 * f1 +
                   a0 * c3 * d5 * e1 * f2 - a0 * c3 * d5 * e2 * f1 - a0 * c5 * d1 * e2 * f3 + a0 * c5 * d1 * e3 * f2 +
                   a0 * c5 * d2 * e1 * f3 - a0 * c5 * d2 * e3 * f1 - a0 * c5 * d3 * e1 * f2 + a0 * c5 * d3 * e2 * f1 -
                   a1 * c0 * d2 * e3 * f5 + a1 * c0 * d2 * e5 * f3 + a1 * c0 * d3 * e2 * f5 - a1 * c0 * d3 * e5 * f2 -
                   a1 * c0 * d5 * e2 * f3 + a1 * c0 * d5 * e3 * f2 + a1 * c2 * d0 * e3 * f5 - a1 * c2 * d0 * e5 * f3 -
                   a1 * c2 * d3 * e0 * f5 + a1 * c2 * d3 * e5 * f0 + a1 * c2 * d5 * e0 * f3 - a1 * c2 * d5 * e3 * f0 -
                   a1 * c3 * d0 * e2 * f5 + a1 * c3 * d0 * e5 * f2 + a1 * c3 * d2 * e0 * f5 - a1 * c3 * d2 * e5 * f0 -
                   a1 * c3 * d5 * e0 * f2 + a1 * c3 * d5 * e2 * f0 + a1 * c5 * d0 * e2 * f3 - a1 * c5 * d0 * e3 * f2 -
                   a1 * c5 * d2 * e0 * f3 + a1 * c5 * d2 * e3 * f0 + a1 * c5 * d3 * e0 * f2 - a1 * c5 * d3 * e2 * f0 +
                   a2 * c0 * d1 * e3 * f5 - a2 * c0 * d1 * e5 * f3 - a2 * c0 * d3 * e1 * f5 + a2 * c0 * d3 * e5 * f1 +
                   a2 * c0 * d5 * e1 * f3 - a2 * c0 * d5 * e3 * f1 - a2 * c1 * d0 * e3 * f5 + a2 * c1 * d0 * e5 * f3 +
                   a2 * c1 * d3 * e0 * f5 - a2 * c1 * d3 * e5 * f0 - a2 * c1 * d5 * e0 * f3 + a2 * c1 * d5 * e3 * f0 +
                   a2 * c3 * d0 * e1 * f5 - a2 * c3 * d0 * e5 * f1 - a2 * c3 * d1 * e0 * f5 + a2 * c3 * d1 * e5 * f0 +
                   a2 * c3 * d5 * e0 * f1 - a2 * c3 * d5 * e1 * f0 - a2 * c5 * d0 * e1 * f3 + a2 * c5 * d0 * e3 * f1 +
                   a2 * c5 * d1 * e0 * f3 - a2 * c5 * d1 * e3 * f0 - a2 * c5 * d3 * e0 * f1 + a2 * c5 * d3 * e1 * f0 -
                   a3 * c0 * d1 * e2 * f5 + a3 * c0 * d1 * e5 * f2 + a3 * c0 * d2 * e1 * f5 - a3 * c0 * d2 * e5 * f1 -
                   a3 * c0 * d5 * e1 * f2 + a3 * c0 * d5 * e2 * f1 + a3 * c1 * d0 * e2 * f5 - a3 * c1 * d0 * e5 * f2 -
                   a3 * c1 * d2 * e0 * f5 + a3 * c1 * d2 * e5 * f0 + a3 * c1 * d5 * e0 * f2 - a3 * c1 * d5 * e2 * f0 -
                   a3 * c2 * d0 * e1 * f5 + a3 * c2 * d0 * e5 * f1 + a3 * c2 * d1 * e0 * f5 - a3 * c2 * d1 * e5 * f0 -
                   a3 * c2 * d5 * e0 * f1 + a3 * c2 * d5 * e1 * f0 + a3 * c5 * d0 * e1 * f2 - a3 * c5 * d0 * e2 * f1 -
                   a3 * c5 * d1 * e0 * f2 + a3 * c5 * d1 * e2 * f0 + a3 * c5 * d2 * e0 * f1 - a3 * c5 * d2 * e1 * f0 +
                   a5 * c0 * d1 * e2 * f3 - a5 * c0 * d1 * e3 * f2 - a5 * c0 * d2 * e1 * f3 + a5 * c0 * d2 * e3 * f1 +
                   a5 * c0 * d3 * e1 * f2 - a5 * c0 * d3 * e2 * f1 - a5 * c1 * d0 * e2 * f3 + a5 * c1 * d0 * e3 * f2 +
                   a5 * c1 * d2 * e0 * f3 - a5 * c1 * d2 * e3 * f0 - a5 * c1 * d3 * e0 * f2 + a5 * c1 * d3 * e2 * f0 +
                   a5 * c2 * d0 * e1 * f3 - a5 * c2 * d0 * e3 * f1 - a5 * c2 * d1 * e0 * f3 + a5 * c2 * d1 * e3 * f0 +
                   a5 * c2 * d3 * e0 * f1 - a5 * c2 * d3 * e1 * f0 - a5 * c3 * d0 * e1 * f2 + a5 * c3 * d0 * e2 * f1 +
                   a5 * c3 * d1 * e0 * f2 - a5 * c3 * d1 * e2 * f0 - a5 * c3 * d2 * e0 * f1 + a5 * c3 * d2 * e1 * f0);
        T[2][4] = t363 *
                  (a0 * b1 * d2 * e3 * f5 - a0 * b1 * d2 * e5 * f3 - a0 * b1 * d3 * e2 * f5 + a0 * b1 * d3 * e5 * f2 +
                   a0 * b1 * d5 * e2 * f3 - a0 * b1 * d5 * e3 * f2 - a0 * b2 * d1 * e3 * f5 + a0 * b2 * d1 * e5 * f3 +
                   a0 * b2 * d3 * e1 * f5 - a0 * b2 * d3 * e5 * f1 - a0 * b2 * d5 * e1 * f3 + a0 * b2 * d5 * e3 * f1 +
                   a0 * b3 * d1 * e2 * f5 - a0 * b3 * d1 * e5 * f2 - a0 * b3 * d2 * e1 * f5 + a0 * b3 * d2 * e5 * f1 +
                   a0 * b3 * d5 * e1 * f2 - a0 * b3 * d5 * e2 * f1 - a0 * b5 * d1 * e2 * f3 + a0 * b5 * d1 * e3 * f2 +
                   a0 * b5 * d2 * e1 * f3 - a0 * b5 * d2 * e3 * f1 - a0 * b5 * d3 * e1 * f2 + a0 * b5 * d3 * e2 * f1 -
                   a1 * b0 * d2 * e3 * f5 + a1 * b0 * d2 * e5 * f3 + a1 * b0 * d3 * e2 * f5 - a1 * b0 * d3 * e5 * f2 -
                   a1 * b0 * d5 * e2 * f3 + a1 * b0 * d5 * e3 * f2 + a1 * b2 * d0 * e3 * f5 - a1 * b2 * d0 * e5 * f3 -
                   a1 * b2 * d3 * e0 * f5 + a1 * b2 * d3 * e5 * f0 + a1 * b2 * d5 * e0 * f3 - a1 * b2 * d5 * e3 * f0 -
                   a1 * b3 * d0 * e2 * f5 + a1 * b3 * d0 * e5 * f2 + a1 * b3 * d2 * e0 * f5 - a1 * b3 * d2 * e5 * f0 -
                   a1 * b3 * d5 * e0 * f2 + a1 * b3 * d5 * e2 * f0 + a1 * b5 * d0 * e2 * f3 - a1 * b5 * d0 * e3 * f2 -
                   a1 * b5 * d2 * e0 * f3 + a1 * b5 * d2 * e3 * f0 + a1 * b5 * d3 * e0 * f2 - a1 * b5 * d3 * e2 * f0 +
                   a2 * b0 * d1 * e3 * f5 - a2 * b0 * d1 * e5 * f3 - a2 * b0 * d3 * e1 * f5 + a2 * b0 * d3 * e5 * f1 +
                   a2 * b0 * d5 * e1 * f3 - a2 * b0 * d5 * e3 * f1 - a2 * b1 * d0 * e3 * f5 + a2 * b1 * d0 * e5 * f3 +
                   a2 * b1 * d3 * e0 * f5 - a2 * b1 * d3 * e5 * f0 - a2 * b1 * d5 * e0 * f3 + a2 * b1 * d5 * e3 * f0 +
                   a2 * b3 * d0 * e1 * f5 - a2 * b3 * d0 * e5 * f1 - a2 * b3 * d1 * e0 * f5 + a2 * b3 * d1 * e5 * f0 +
                   a2 * b3 * d5 * e0 * f1 - a2 * b3 * d5 * e1 * f0 - a2 * b5 * d0 * e1 * f3 + a2 * b5 * d0 * e3 * f1 +
                   a2 * b5 * d1 * e0 * f3 - a2 * b5 * d1 * e3 * f0 - a2 * b5 * d3 * e0 * f1 + a2 * b5 * d3 * e1 * f0 -
                   a3 * b0 * d1 * e2 * f5 + a3 * b0 * d1 * e5 * f2 + a3 * b0 * d2 * e1 * f5 - a3 * b0 * d2 * e5 * f1 -
                   a3 * b0 * d5 * e1 * f2 + a3 * b0 * d5 * e2 * f1 + a3 * b1 * d0 * e2 * f5 - a3 * b1 * d0 * e5 * f2 -
                   a3 * b1 * d2 * e0 * f5 + a3 * b1 * d2 * e5 * f0 + a3 * b1 * d5 * e0 * f2 - a3 * b1 * d5 * e2 * f0 -
                   a3 * b2 * d0 * e1 * f5 + a3 * b2 * d0 * e5 * f1 + a3 * b2 * d1 * e0 * f5 - a3 * b2 * d1 * e5 * f0 -
                   a3 * b2 * d5 * e0 * f1 + a3 * b2 * d5 * e1 * f0 + a3 * b5 * d0 * e1 * f2 - a3 * b5 * d0 * e2 * f1 -
                   a3 * b5 * d1 * e0 * f2 + a3 * b5 * d1 * e2 * f0 + a3 * b5 * d2 * e0 * f1 - a3 * b5 * d2 * e1 * f0 +
                   a5 * b0 * d1 * e2 * f3 - a5 * b0 * d1 * e3 * f2 - a5 * b0 * d2 * e1 * f3 + a5 * b0 * d2 * e3 * f1 +
                   a5 * b0 * d3 * e1 * f2 - a5 * b0 * d3 * e2 * f1 - a5 * b1 * d0 * e2 * f3 + a5 * b1 * d0 * e3 * f2 +
                   a5 * b1 * d2 * e0 * f3 - a5 * b1 * d2 * e3 * f0 - a5 * b1 * d3 * e0 * f2 + a5 * b1 * d3 * e2 * f0 +
                   a5 * b2 * d0 * e1 * f3 - a5 * b2 * d0 * e3 * f1 - a5 * b2 * d1 * e0 * f3 + a5 * b2 * d1 * e3 * f0 +
                   a5 * b2 * d3 * e0 * f1 - a5 * b2 * d3 * e1 * f0 - a5 * b3 * d0 * e1 * f2 + a5 * b3 * d0 * e2 * f1 +
                   a5 * b3 * d1 * e0 * f2 - a5 * b3 * d1 * e2 * f0 - a5 * b3 * d2 * e0 * f1 + a5 * b3 * d2 * e1 * f0);
        T[3][4] = -t363 *
                  (a0 * b1 * c2 * e3 * f5 - a0 * b1 * c2 * e5 * f3 - a0 * b1 * c3 * e2 * f5 + a0 * b1 * c3 * e5 * f2 +
                   a0 * b1 * c5 * e2 * f3 - a0 * b1 * c5 * e3 * f2 - a0 * b2 * c1 * e3 * f5 + a0 * b2 * c1 * e5 * f3 +
                   a0 * b2 * c3 * e1 * f5 - a0 * b2 * c3 * e5 * f1 - a0 * b2 * c5 * e1 * f3 + a0 * b2 * c5 * e3 * f1 +
                   a0 * b3 * c1 * e2 * f5 - a0 * b3 * c1 * e5 * f2 - a0 * b3 * c2 * e1 * f5 + a0 * b3 * c2 * e5 * f1 +
                   a0 * b3 * c5 * e1 * f2 - a0 * b3 * c5 * e2 * f1 - a0 * b5 * c1 * e2 * f3 + a0 * b5 * c1 * e3 * f2 +
                   a0 * b5 * c2 * e1 * f3 - a0 * b5 * c2 * e3 * f1 - a0 * b5 * c3 * e1 * f2 + a0 * b5 * c3 * e2 * f1 -
                   a1 * b0 * c2 * e3 * f5 + a1 * b0 * c2 * e5 * f3 + a1 * b0 * c3 * e2 * f5 - a1 * b0 * c3 * e5 * f2 -
                   a1 * b0 * c5 * e2 * f3 + a1 * b0 * c5 * e3 * f2 + a1 * b2 * c0 * e3 * f5 - a1 * b2 * c0 * e5 * f3 -
                   a1 * b2 * c3 * e0 * f5 + a1 * b2 * c3 * e5 * f0 + a1 * b2 * c5 * e0 * f3 - a1 * b2 * c5 * e3 * f0 -
                   a1 * b3 * c0 * e2 * f5 + a1 * b3 * c0 * e5 * f2 + a1 * b3 * c2 * e0 * f5 - a1 * b3 * c2 * e5 * f0 -
                   a1 * b3 * c5 * e0 * f2 + a1 * b3 * c5 * e2 * f0 + a1 * b5 * c0 * e2 * f3 - a1 * b5 * c0 * e3 * f2 -
                   a1 * b5 * c2 * e0 * f3 + a1 * b5 * c2 * e3 * f0 + a1 * b5 * c3 * e0 * f2 - a1 * b5 * c3 * e2 * f0 +
                   a2 * b0 * c1 * e3 * f5 - a2 * b0 * c1 * e5 * f3 - a2 * b0 * c3 * e1 * f5 + a2 * b0 * c3 * e5 * f1 +
                   a2 * b0 * c5 * e1 * f3 - a2 * b0 * c5 * e3 * f1 - a2 * b1 * c0 * e3 * f5 + a2 * b1 * c0 * e5 * f3 +
                   a2 * b1 * c3 * e0 * f5 - a2 * b1 * c3 * e5 * f0 - a2 * b1 * c5 * e0 * f3 + a2 * b1 * c5 * e3 * f0 +
                   a2 * b3 * c0 * e1 * f5 - a2 * b3 * c0 * e5 * f1 - a2 * b3 * c1 * e0 * f5 + a2 * b3 * c1 * e5 * f0 +
                   a2 * b3 * c5 * e0 * f1 - a2 * b3 * c5 * e1 * f0 - a2 * b5 * c0 * e1 * f3 + a2 * b5 * c0 * e3 * f1 +
                   a2 * b5 * c1 * e0 * f3 - a2 * b5 * c1 * e3 * f0 - a2 * b5 * c3 * e0 * f1 + a2 * b5 * c3 * e1 * f0 -
                   a3 * b0 * c1 * e2 * f5 + a3 * b0 * c1 * e5 * f2 + a3 * b0 * c2 * e1 * f5 - a3 * b0 * c2 * e5 * f1 -
                   a3 * b0 * c5 * e1 * f2 + a3 * b0 * c5 * e2 * f1 + a3 * b1 * c0 * e2 * f5 - a3 * b1 * c0 * e5 * f2 -
                   a3 * b1 * c2 * e0 * f5 + a3 * b1 * c2 * e5 * f0 + a3 * b1 * c5 * e0 * f2 - a3 * b1 * c5 * e2 * f0 -
                   a3 * b2 * c0 * e1 * f5 + a3 * b2 * c0 * e5 * f1 + a3 * b2 * c1 * e0 * f5 - a3 * b2 * c1 * e5 * f0 -
                   a3 * b2 * c5 * e0 * f1 + a3 * b2 * c5 * e1 * f0 + a3 * b5 * c0 * e1 * f2 - a3 * b5 * c0 * e2 * f1 -
                   a3 * b5 * c1 * e0 * f2 + a3 * b5 * c1 * e2 * f0 + a3 * b5 * c2 * e0 * f1 - a3 * b5 * c2 * e1 * f0 +
                   a5 * b0 * c1 * e2 * f3 - a5 * b0 * c1 * e3 * f2 - a5 * b0 * c2 * e1 * f3 + a5 * b0 * c2 * e3 * f1 +
                   a5 * b0 * c3 * e1 * f2 - a5 * b0 * c3 * e2 * f1 - a5 * b1 * c0 * e2 * f3 + a5 * b1 * c0 * e3 * f2 +
                   a5 * b1 * c2 * e0 * f3 - a5 * b1 * c2 * e3 * f0 - a5 * b1 * c3 * e0 * f2 + a5 * b1 * c3 * e2 * f0 +
                   a5 * b2 * c0 * e1 * f3 - a5 * b2 * c0 * e3 * f1 - a5 * b2 * c1 * e0 * f3 + a5 * b2 * c1 * e3 * f0 +
                   a5 * b2 * c3 * e0 * f1 - a5 * b2 * c3 * e1 * f0 - a5 * b3 * c0 * e1 * f2 + a5 * b3 * c0 * e2 * f1 +
                   a5 * b3 * c1 * e0 * f2 - a5 * b3 * c1 * e2 * f0 - a5 * b3 * c2 * e0 * f1 + a5 * b3 * c2 * e1 * f0);
        T[4][4] = t363 *
                  (a0 * b1 * c2 * d3 * f5 - a0 * b1 * c2 * d5 * f3 - a0 * b1 * c3 * d2 * f5 + a0 * b1 * c3 * d5 * f2 +
                   a0 * b1 * c5 * d2 * f3 - a0 * b1 * c5 * d3 * f2 - a0 * b2 * c1 * d3 * f5 + a0 * b2 * c1 * d5 * f3 +
                   a0 * b2 * c3 * d1 * f5 - a0 * b2 * c3 * d5 * f1 - a0 * b2 * c5 * d1 * f3 + a0 * b2 * c5 * d3 * f1 +
                   a0 * b3 * c1 * d2 * f5 - a0 * b3 * c1 * d5 * f2 - a0 * b3 * c2 * d1 * f5 + a0 * b3 * c2 * d5 * f1 +
                   a0 * b3 * c5 * d1 * f2 - a0 * b3 * c5 * d2 * f1 - a0 * b5 * c1 * d2 * f3 + a0 * b5 * c1 * d3 * f2 +
                   a0 * b5 * c2 * d1 * f3 - a0 * b5 * c2 * d3 * f1 - a0 * b5 * c3 * d1 * f2 + a0 * b5 * c3 * d2 * f1 -
                   a1 * b0 * c2 * d3 * f5 + a1 * b0 * c2 * d5 * f3 + a1 * b0 * c3 * d2 * f5 - a1 * b0 * c3 * d5 * f2 -
                   a1 * b0 * c5 * d2 * f3 + a1 * b0 * c5 * d3 * f2 + a1 * b2 * c0 * d3 * f5 - a1 * b2 * c0 * d5 * f3 -
                   a1 * b2 * c3 * d0 * f5 + a1 * b2 * c3 * d5 * f0 + a1 * b2 * c5 * d0 * f3 - a1 * b2 * c5 * d3 * f0 -
                   a1 * b3 * c0 * d2 * f5 + a1 * b3 * c0 * d5 * f2 + a1 * b3 * c2 * d0 * f5 - a1 * b3 * c2 * d5 * f0 -
                   a1 * b3 * c5 * d0 * f2 + a1 * b3 * c5 * d2 * f0 + a1 * b5 * c0 * d2 * f3 - a1 * b5 * c0 * d3 * f2 -
                   a1 * b5 * c2 * d0 * f3 + a1 * b5 * c2 * d3 * f0 + a1 * b5 * c3 * d0 * f2 - a1 * b5 * c3 * d2 * f0 +
                   a2 * b0 * c1 * d3 * f5 - a2 * b0 * c1 * d5 * f3 - a2 * b0 * c3 * d1 * f5 + a2 * b0 * c3 * d5 * f1 +
                   a2 * b0 * c5 * d1 * f3 - a2 * b0 * c5 * d3 * f1 - a2 * b1 * c0 * d3 * f5 + a2 * b1 * c0 * d5 * f3 +
                   a2 * b1 * c3 * d0 * f5 - a2 * b1 * c3 * d5 * f0 - a2 * b1 * c5 * d0 * f3 + a2 * b1 * c5 * d3 * f0 +
                   a2 * b3 * c0 * d1 * f5 - a2 * b3 * c0 * d5 * f1 - a2 * b3 * c1 * d0 * f5 + a2 * b3 * c1 * d5 * f0 +
                   a2 * b3 * c5 * d0 * f1 - a2 * b3 * c5 * d1 * f0 - a2 * b5 * c0 * d1 * f3 + a2 * b5 * c0 * d3 * f1 +
                   a2 * b5 * c1 * d0 * f3 - a2 * b5 * c1 * d3 * f0 - a2 * b5 * c3 * d0 * f1 + a2 * b5 * c3 * d1 * f0 -
                   a3 * b0 * c1 * d2 * f5 + a3 * b0 * c1 * d5 * f2 + a3 * b0 * c2 * d1 * f5 - a3 * b0 * c2 * d5 * f1 -
                   a3 * b0 * c5 * d1 * f2 + a3 * b0 * c5 * d2 * f1 + a3 * b1 * c0 * d2 * f5 - a3 * b1 * c0 * d5 * f2 -
                   a3 * b1 * c2 * d0 * f5 + a3 * b1 * c2 * d5 * f0 + a3 * b1 * c5 * d0 * f2 - a3 * b1 * c5 * d2 * f0 -
                   a3 * b2 * c0 * d1 * f5 + a3 * b2 * c0 * d5 * f1 + a3 * b2 * c1 * d0 * f5 - a3 * b2 * c1 * d5 * f0 -
                   a3 * b2 * c5 * d0 * f1 + a3 * b2 * c5 * d1 * f0 + a3 * b5 * c0 * d1 * f2 - a3 * b5 * c0 * d2 * f1 -
                   a3 * b5 * c1 * d0 * f2 + a3 * b5 * c1 * d2 * f0 + a3 * b5 * c2 * d0 * f1 - a3 * b5 * c2 * d1 * f0 +
                   a5 * b0 * c1 * d2 * f3 - a5 * b0 * c1 * d3 * f2 - a5 * b0 * c2 * d1 * f3 + a5 * b0 * c2 * d3 * f1 +
                   a5 * b0 * c3 * d1 * f2 - a5 * b0 * c3 * d2 * f1 - a5 * b1 * c0 * d2 * f3 + a5 * b1 * c0 * d3 * f2 +
                   a5 * b1 * c2 * d0 * f3 - a5 * b1 * c2 * d3 * f0 - a5 * b1 * c3 * d0 * f2 + a5 * b1 * c3 * d2 * f0 +
                   a5 * b2 * c0 * d1 * f3 - a5 * b2 * c0 * d3 * f1 - a5 * b2 * c1 * d0 * f3 + a5 * b2 * c1 * d3 * f0 +
                   a5 * b2 * c3 * d0 * f1 - a5 * b2 * c3 * d1 * f0 - a5 * b3 * c0 * d1 * f2 + a5 * b3 * c0 * d2 * f1 +
                   a5 * b3 * c1 * d0 * f2 - a5 * b3 * c1 * d2 * f0 - a5 * b3 * c2 * d0 * f1 + a5 * b3 * c2 * d1 * f0);
        T[5][4] = -t363 *
                  (a0 * b1 * c2 * d3 * e5 - a0 * b1 * c2 * d5 * e3 - a0 * b1 * c3 * d2 * e5 + a0 * b1 * c3 * d5 * e2 +
                   a0 * b1 * c5 * d2 * e3 - a0 * b1 * c5 * d3 * e2 - a0 * b2 * c1 * d3 * e5 + a0 * b2 * c1 * d5 * e3 +
                   a0 * b2 * c3 * d1 * e5 - a0 * b2 * c3 * d5 * e1 - a0 * b2 * c5 * d1 * e3 + a0 * b2 * c5 * d3 * e1 +
                   a0 * b3 * c1 * d2 * e5 - a0 * b3 * c1 * d5 * e2 - a0 * b3 * c2 * d1 * e5 + a0 * b3 * c2 * d5 * e1 +
                   a0 * b3 * c5 * d1 * e2 - a0 * b3 * c5 * d2 * e1 - a0 * b5 * c1 * d2 * e3 + a0 * b5 * c1 * d3 * e2 +
                   a0 * b5 * c2 * d1 * e3 - a0 * b5 * c2 * d3 * e1 - a0 * b5 * c3 * d1 * e2 + a0 * b5 * c3 * d2 * e1 -
                   a1 * b0 * c2 * d3 * e5 + a1 * b0 * c2 * d5 * e3 + a1 * b0 * c3 * d2 * e5 - a1 * b0 * c3 * d5 * e2 -
                   a1 * b0 * c5 * d2 * e3 + a1 * b0 * c5 * d3 * e2 + a1 * b2 * c0 * d3 * e5 - a1 * b2 * c0 * d5 * e3 -
                   a1 * b2 * c3 * d0 * e5 + a1 * b2 * c3 * d5 * e0 + a1 * b2 * c5 * d0 * e3 - a1 * b2 * c5 * d3 * e0 -
                   a1 * b3 * c0 * d2 * e5 + a1 * b3 * c0 * d5 * e2 + a1 * b3 * c2 * d0 * e5 - a1 * b3 * c2 * d5 * e0 -
                   a1 * b3 * c5 * d0 * e2 + a1 * b3 * c5 * d2 * e0 + a1 * b5 * c0 * d2 * e3 - a1 * b5 * c0 * d3 * e2 -
                   a1 * b5 * c2 * d0 * e3 + a1 * b5 * c2 * d3 * e0 + a1 * b5 * c3 * d0 * e2 - a1 * b5 * c3 * d2 * e0 +
                   a2 * b0 * c1 * d3 * e5 - a2 * b0 * c1 * d5 * e3 - a2 * b0 * c3 * d1 * e5 + a2 * b0 * c3 * d5 * e1 +
                   a2 * b0 * c5 * d1 * e3 - a2 * b0 * c5 * d3 * e1 - a2 * b1 * c0 * d3 * e5 + a2 * b1 * c0 * d5 * e3 +
                   a2 * b1 * c3 * d0 * e5 - a2 * b1 * c3 * d5 * e0 - a2 * b1 * c5 * d0 * e3 + a2 * b1 * c5 * d3 * e0 +
                   a2 * b3 * c0 * d1 * e5 - a2 * b3 * c0 * d5 * e1 - a2 * b3 * c1 * d0 * e5 + a2 * b3 * c1 * d5 * e0 +
                   a2 * b3 * c5 * d0 * e1 - a2 * b3 * c5 * d1 * e0 - a2 * b5 * c0 * d1 * e3 + a2 * b5 * c0 * d3 * e1 +
                   a2 * b5 * c1 * d0 * e3 - a2 * b5 * c1 * d3 * e0 - a2 * b5 * c3 * d0 * e1 + a2 * b5 * c3 * d1 * e0 -
                   a3 * b0 * c1 * d2 * e5 + a3 * b0 * c1 * d5 * e2 + a3 * b0 * c2 * d1 * e5 - a3 * b0 * c2 * d5 * e1 -
                   a3 * b0 * c5 * d1 * e2 + a3 * b0 * c5 * d2 * e1 + a3 * b1 * c0 * d2 * e5 - a3 * b1 * c0 * d5 * e2 -
                   a3 * b1 * c2 * d0 * e5 + a3 * b1 * c2 * d5 * e0 + a3 * b1 * c5 * d0 * e2 - a3 * b1 * c5 * d2 * e0 -
                   a3 * b2 * c0 * d1 * e5 + a3 * b2 * c0 * d5 * e1 + a3 * b2 * c1 * d0 * e5 - a3 * b2 * c1 * d5 * e0 -
                   a3 * b2 * c5 * d0 * e1 + a3 * b2 * c5 * d1 * e0 + a3 * b5 * c0 * d1 * e2 - a3 * b5 * c0 * d2 * e1 -
                   a3 * b5 * c1 * d0 * e2 + a3 * b5 * c1 * d2 * e0 + a3 * b5 * c2 * d0 * e1 - a3 * b5 * c2 * d1 * e0 +
                   a5 * b0 * c1 * d2 * e3 - a5 * b0 * c1 * d3 * e2 - a5 * b0 * c2 * d1 * e3 + a5 * b0 * c2 * d3 * e1 +
                   a5 * b0 * c3 * d1 * e2 - a5 * b0 * c3 * d2 * e1 - a5 * b1 * c0 * d2 * e3 + a5 * b1 * c0 * d3 * e2 +
                   a5 * b1 * c2 * d0 * e3 - a5 * b1 * c2 * d3 * e0 - a5 * b1 * c3 * d0 * e2 + a5 * b1 * c3 * d2 * e0 +
                   a5 * b2 * c0 * d1 * e3 - a5 * b2 * c0 * d3 * e1 - a5 * b2 * c1 * d0 * e3 + a5 * b2 * c1 * d3 * e0 +
                   a5 * b2 * c3 * d0 * e1 - a5 * b2 * c3 * d1 * e0 - a5 * b3 * c0 * d1 * e2 + a5 * b3 * c0 * d2 * e1 +
                   a5 * b3 * c1 * d0 * e2 - a5 * b3 * c1 * d2 * e0 - a5 * b3 * c2 * d0 * e1 + a5 * b3 * c2 * d1 * e0);
        T[0][5] = -t363 *
                  (b0 * c1 * d2 * e3 * f4 - b0 * c1 * d2 * e4 * f3 - b0 * c1 * d3 * e2 * f4 + b0 * c1 * d3 * e4 * f2 +
                   b0 * c1 * d4 * e2 * f3 - b0 * c1 * d4 * e3 * f2 - b0 * c2 * d1 * e3 * f4 + b0 * c2 * d1 * e4 * f3 +
                   b0 * c2 * d3 * e1 * f4 - b0 * c2 * d3 * e4 * f1 - b0 * c2 * d4 * e1 * f3 + b0 * c2 * d4 * e3 * f1 +
                   b0 * c3 * d1 * e2 * f4 - b0 * c3 * d1 * e4 * f2 - b0 * c3 * d2 * e1 * f4 + b0 * c3 * d2 * e4 * f1 +
                   b0 * c3 * d4 * e1 * f2 - b0 * c3 * d4 * e2 * f1 - b0 * c4 * d1 * e2 * f3 + b0 * c4 * d1 * e3 * f2 +
                   b0 * c4 * d2 * e1 * f3 - b0 * c4 * d2 * e3 * f1 - b0 * c4 * d3 * e1 * f2 + b0 * c4 * d3 * e2 * f1 -
                   b1 * c0 * d2 * e3 * f4 + b1 * c0 * d2 * e4 * f3 + b1 * c0 * d3 * e2 * f4 - b1 * c0 * d3 * e4 * f2 -
                   b1 * c0 * d4 * e2 * f3 + b1 * c0 * d4 * e3 * f2 + b1 * c2 * d0 * e3 * f4 - b1 * c2 * d0 * e4 * f3 -
                   b1 * c2 * d3 * e0 * f4 + b1 * c2 * d3 * e4 * f0 + b1 * c2 * d4 * e0 * f3 - b1 * c2 * d4 * e3 * f0 -
                   b1 * c3 * d0 * e2 * f4 + b1 * c3 * d0 * e4 * f2 + b1 * c3 * d2 * e0 * f4 - b1 * c3 * d2 * e4 * f0 -
                   b1 * c3 * d4 * e0 * f2 + b1 * c3 * d4 * e2 * f0 + b1 * c4 * d0 * e2 * f3 - b1 * c4 * d0 * e3 * f2 -
                   b1 * c4 * d2 * e0 * f3 + b1 * c4 * d2 * e3 * f0 + b1 * c4 * d3 * e0 * f2 - b1 * c4 * d3 * e2 * f0 +
                   b2 * c0 * d1 * e3 * f4 - b2 * c0 * d1 * e4 * f3 - b2 * c0 * d3 * e1 * f4 + b2 * c0 * d3 * e4 * f1 +
                   b2 * c0 * d4 * e1 * f3 - b2 * c0 * d4 * e3 * f1 - b2 * c1 * d0 * e3 * f4 + b2 * c1 * d0 * e4 * f3 +
                   b2 * c1 * d3 * e0 * f4 - b2 * c1 * d3 * e4 * f0 - b2 * c1 * d4 * e0 * f3 + b2 * c1 * d4 * e3 * f0 +
                   b2 * c3 * d0 * e1 * f4 - b2 * c3 * d0 * e4 * f1 - b2 * c3 * d1 * e0 * f4 + b2 * c3 * d1 * e4 * f0 +
                   b2 * c3 * d4 * e0 * f1 - b2 * c3 * d4 * e1 * f0 - b2 * c4 * d0 * e1 * f3 + b2 * c4 * d0 * e3 * f1 +
                   b2 * c4 * d1 * e0 * f3 - b2 * c4 * d1 * e3 * f0 - b2 * c4 * d3 * e0 * f1 + b2 * c4 * d3 * e1 * f0 -
                   b3 * c0 * d1 * e2 * f4 + b3 * c0 * d1 * e4 * f2 + b3 * c0 * d2 * e1 * f4 - b3 * c0 * d2 * e4 * f1 -
                   b3 * c0 * d4 * e1 * f2 + b3 * c0 * d4 * e2 * f1 + b3 * c1 * d0 * e2 * f4 - b3 * c1 * d0 * e4 * f2 -
                   b3 * c1 * d2 * e0 * f4 + b3 * c1 * d2 * e4 * f0 + b3 * c1 * d4 * e0 * f2 - b3 * c1 * d4 * e2 * f0 -
                   b3 * c2 * d0 * e1 * f4 + b3 * c2 * d0 * e4 * f1 + b3 * c2 * d1 * e0 * f4 - b3 * c2 * d1 * e4 * f0 -
                   b3 * c2 * d4 * e0 * f1 + b3 * c2 * d4 * e1 * f0 + b3 * c4 * d0 * e1 * f2 - b3 * c4 * d0 * e2 * f1 -
                   b3 * c4 * d1 * e0 * f2 + b3 * c4 * d1 * e2 * f0 + b3 * c4 * d2 * e0 * f1 - b3 * c4 * d2 * e1 * f0 +
                   b4 * c0 * d1 * e2 * f3 - b4 * c0 * d1 * e3 * f2 - b4 * c0 * d2 * e1 * f3 + b4 * c0 * d2 * e3 * f1 +
                   b4 * c0 * d3 * e1 * f2 - b4 * c0 * d3 * e2 * f1 - b4 * c1 * d0 * e2 * f3 + b4 * c1 * d0 * e3 * f2 +
                   b4 * c1 * d2 * e0 * f3 - b4 * c1 * d2 * e3 * f0 - b4 * c1 * d3 * e0 * f2 + b4 * c1 * d3 * e2 * f0 +
                   b4 * c2 * d0 * e1 * f3 - b4 * c2 * d0 * e3 * f1 - b4 * c2 * d1 * e0 * f3 + b4 * c2 * d1 * e3 * f0 +
                   b4 * c2 * d3 * e0 * f1 - b4 * c2 * d3 * e1 * f0 - b4 * c3 * d0 * e1 * f2 + b4 * c3 * d0 * e2 * f1 +
                   b4 * c3 * d1 * e0 * f2 - b4 * c3 * d1 * e2 * f0 - b4 * c3 * d2 * e0 * f1 + b4 * c3 * d2 * e1 * f0);
        T[1][5] = t363 *
                  (a0 * c1 * d2 * e3 * f4 - a0 * c1 * d2 * e4 * f3 - a0 * c1 * d3 * e2 * f4 + a0 * c1 * d3 * e4 * f2 +
                   a0 * c1 * d4 * e2 * f3 - a0 * c1 * d4 * e3 * f2 - a0 * c2 * d1 * e3 * f4 + a0 * c2 * d1 * e4 * f3 +
                   a0 * c2 * d3 * e1 * f4 - a0 * c2 * d3 * e4 * f1 - a0 * c2 * d4 * e1 * f3 + a0 * c2 * d4 * e3 * f1 +
                   a0 * c3 * d1 * e2 * f4 - a0 * c3 * d1 * e4 * f2 - a0 * c3 * d2 * e1 * f4 + a0 * c3 * d2 * e4 * f1 +
                   a0 * c3 * d4 * e1 * f2 - a0 * c3 * d4 * e2 * f1 - a0 * c4 * d1 * e2 * f3 + a0 * c4 * d1 * e3 * f2 +
                   a0 * c4 * d2 * e1 * f3 - a0 * c4 * d2 * e3 * f1 - a0 * c4 * d3 * e1 * f2 + a0 * c4 * d3 * e2 * f1 -
                   a1 * c0 * d2 * e3 * f4 + a1 * c0 * d2 * e4 * f3 + a1 * c0 * d3 * e2 * f4 - a1 * c0 * d3 * e4 * f2 -
                   a1 * c0 * d4 * e2 * f3 + a1 * c0 * d4 * e3 * f2 + a1 * c2 * d0 * e3 * f4 - a1 * c2 * d0 * e4 * f3 -
                   a1 * c2 * d3 * e0 * f4 + a1 * c2 * d3 * e4 * f0 + a1 * c2 * d4 * e0 * f3 - a1 * c2 * d4 * e3 * f0 -
                   a1 * c3 * d0 * e2 * f4 + a1 * c3 * d0 * e4 * f2 + a1 * c3 * d2 * e0 * f4 - a1 * c3 * d2 * e4 * f0 -
                   a1 * c3 * d4 * e0 * f2 + a1 * c3 * d4 * e2 * f0 + a1 * c4 * d0 * e2 * f3 - a1 * c4 * d0 * e3 * f2 -
                   a1 * c4 * d2 * e0 * f3 + a1 * c4 * d2 * e3 * f0 + a1 * c4 * d3 * e0 * f2 - a1 * c4 * d3 * e2 * f0 +
                   a2 * c0 * d1 * e3 * f4 - a2 * c0 * d1 * e4 * f3 - a2 * c0 * d3 * e1 * f4 + a2 * c0 * d3 * e4 * f1 +
                   a2 * c0 * d4 * e1 * f3 - a2 * c0 * d4 * e3 * f1 - a2 * c1 * d0 * e3 * f4 + a2 * c1 * d0 * e4 * f3 +
                   a2 * c1 * d3 * e0 * f4 - a2 * c1 * d3 * e4 * f0 - a2 * c1 * d4 * e0 * f3 + a2 * c1 * d4 * e3 * f0 +
                   a2 * c3 * d0 * e1 * f4 - a2 * c3 * d0 * e4 * f1 - a2 * c3 * d1 * e0 * f4 + a2 * c3 * d1 * e4 * f0 +
                   a2 * c3 * d4 * e0 * f1 - a2 * c3 * d4 * e1 * f0 - a2 * c4 * d0 * e1 * f3 + a2 * c4 * d0 * e3 * f1 +
                   a2 * c4 * d1 * e0 * f3 - a2 * c4 * d1 * e3 * f0 - a2 * c4 * d3 * e0 * f1 + a2 * c4 * d3 * e1 * f0 -
                   a3 * c0 * d1 * e2 * f4 + a3 * c0 * d1 * e4 * f2 + a3 * c0 * d2 * e1 * f4 - a3 * c0 * d2 * e4 * f1 -
                   a3 * c0 * d4 * e1 * f2 + a3 * c0 * d4 * e2 * f1 + a3 * c1 * d0 * e2 * f4 - a3 * c1 * d0 * e4 * f2 -
                   a3 * c1 * d2 * e0 * f4 + a3 * c1 * d2 * e4 * f0 + a3 * c1 * d4 * e0 * f2 - a3 * c1 * d4 * e2 * f0 -
                   a3 * c2 * d0 * e1 * f4 + a3 * c2 * d0 * e4 * f1 + a3 * c2 * d1 * e0 * f4 - a3 * c2 * d1 * e4 * f0 -
                   a3 * c2 * d4 * e0 * f1 + a3 * c2 * d4 * e1 * f0 + a3 * c4 * d0 * e1 * f2 - a3 * c4 * d0 * e2 * f1 -
                   a3 * c4 * d1 * e0 * f2 + a3 * c4 * d1 * e2 * f0 + a3 * c4 * d2 * e0 * f1 - a3 * c4 * d2 * e1 * f0 +
                   a4 * c0 * d1 * e2 * f3 - a4 * c0 * d1 * e3 * f2 - a4 * c0 * d2 * e1 * f3 + a4 * c0 * d2 * e3 * f1 +
                   a4 * c0 * d3 * e1 * f2 - a4 * c0 * d3 * e2 * f1 - a4 * c1 * d0 * e2 * f3 + a4 * c1 * d0 * e3 * f2 +
                   a4 * c1 * d2 * e0 * f3 - a4 * c1 * d2 * e3 * f0 - a4 * c1 * d3 * e0 * f2 + a4 * c1 * d3 * e2 * f0 +
                   a4 * c2 * d0 * e1 * f3 - a4 * c2 * d0 * e3 * f1 - a4 * c2 * d1 * e0 * f3 + a4 * c2 * d1 * e3 * f0 +
                   a4 * c2 * d3 * e0 * f1 - a4 * c2 * d3 * e1 * f0 - a4 * c3 * d0 * e1 * f2 + a4 * c3 * d0 * e2 * f1 +
                   a4 * c3 * d1 * e0 * f2 - a4 * c3 * d1 * e2 * f0 - a4 * c3 * d2 * e0 * f1 + a4 * c3 * d2 * e1 * f0);
        T[2][5] = -t363 *
                  (a0 * b1 * d2 * e3 * f4 - a0 * b1 * d2 * e4 * f3 - a0 * b1 * d3 * e2 * f4 + a0 * b1 * d3 * e4 * f2 +
                   a0 * b1 * d4 * e2 * f3 - a0 * b1 * d4 * e3 * f2 - a0 * b2 * d1 * e3 * f4 + a0 * b2 * d1 * e4 * f3 +
                   a0 * b2 * d3 * e1 * f4 - a0 * b2 * d3 * e4 * f1 - a0 * b2 * d4 * e1 * f3 + a0 * b2 * d4 * e3 * f1 +
                   a0 * b3 * d1 * e2 * f4 - a0 * b3 * d1 * e4 * f2 - a0 * b3 * d2 * e1 * f4 + a0 * b3 * d2 * e4 * f1 +
                   a0 * b3 * d4 * e1 * f2 - a0 * b3 * d4 * e2 * f1 - a0 * b4 * d1 * e2 * f3 + a0 * b4 * d1 * e3 * f2 +
                   a0 * b4 * d2 * e1 * f3 - a0 * b4 * d2 * e3 * f1 - a0 * b4 * d3 * e1 * f2 + a0 * b4 * d3 * e2 * f1 -
                   a1 * b0 * d2 * e3 * f4 + a1 * b0 * d2 * e4 * f3 + a1 * b0 * d3 * e2 * f4 - a1 * b0 * d3 * e4 * f2 -
                   a1 * b0 * d4 * e2 * f3 + a1 * b0 * d4 * e3 * f2 + a1 * b2 * d0 * e3 * f4 - a1 * b2 * d0 * e4 * f3 -
                   a1 * b2 * d3 * e0 * f4 + a1 * b2 * d3 * e4 * f0 + a1 * b2 * d4 * e0 * f3 - a1 * b2 * d4 * e3 * f0 -
                   a1 * b3 * d0 * e2 * f4 + a1 * b3 * d0 * e4 * f2 + a1 * b3 * d2 * e0 * f4 - a1 * b3 * d2 * e4 * f0 -
                   a1 * b3 * d4 * e0 * f2 + a1 * b3 * d4 * e2 * f0 + a1 * b4 * d0 * e2 * f3 - a1 * b4 * d0 * e3 * f2 -
                   a1 * b4 * d2 * e0 * f3 + a1 * b4 * d2 * e3 * f0 + a1 * b4 * d3 * e0 * f2 - a1 * b4 * d3 * e2 * f0 +
                   a2 * b0 * d1 * e3 * f4 - a2 * b0 * d1 * e4 * f3 - a2 * b0 * d3 * e1 * f4 + a2 * b0 * d3 * e4 * f1 +
                   a2 * b0 * d4 * e1 * f3 - a2 * b0 * d4 * e3 * f1 - a2 * b1 * d0 * e3 * f4 + a2 * b1 * d0 * e4 * f3 +
                   a2 * b1 * d3 * e0 * f4 - a2 * b1 * d3 * e4 * f0 - a2 * b1 * d4 * e0 * f3 + a2 * b1 * d4 * e3 * f0 +
                   a2 * b3 * d0 * e1 * f4 - a2 * b3 * d0 * e4 * f1 - a2 * b3 * d1 * e0 * f4 + a2 * b3 * d1 * e4 * f0 +
                   a2 * b3 * d4 * e0 * f1 - a2 * b3 * d4 * e1 * f0 - a2 * b4 * d0 * e1 * f3 + a2 * b4 * d0 * e3 * f1 +
                   a2 * b4 * d1 * e0 * f3 - a2 * b4 * d1 * e3 * f0 - a2 * b4 * d3 * e0 * f1 + a2 * b4 * d3 * e1 * f0 -
                   a3 * b0 * d1 * e2 * f4 + a3 * b0 * d1 * e4 * f2 + a3 * b0 * d2 * e1 * f4 - a3 * b0 * d2 * e4 * f1 -
                   a3 * b0 * d4 * e1 * f2 + a3 * b0 * d4 * e2 * f1 + a3 * b1 * d0 * e2 * f4 - a3 * b1 * d0 * e4 * f2 -
                   a3 * b1 * d2 * e0 * f4 + a3 * b1 * d2 * e4 * f0 + a3 * b1 * d4 * e0 * f2 - a3 * b1 * d4 * e2 * f0 -
                   a3 * b2 * d0 * e1 * f4 + a3 * b2 * d0 * e4 * f1 + a3 * b2 * d1 * e0 * f4 - a3 * b2 * d1 * e4 * f0 -
                   a3 * b2 * d4 * e0 * f1 + a3 * b2 * d4 * e1 * f0 + a3 * b4 * d0 * e1 * f2 - a3 * b4 * d0 * e2 * f1 -
                   a3 * b4 * d1 * e0 * f2 + a3 * b4 * d1 * e2 * f0 + a3 * b4 * d2 * e0 * f1 - a3 * b4 * d2 * e1 * f0 +
                   a4 * b0 * d1 * e2 * f3 - a4 * b0 * d1 * e3 * f2 - a4 * b0 * d2 * e1 * f3 + a4 * b0 * d2 * e3 * f1 +
                   a4 * b0 * d3 * e1 * f2 - a4 * b0 * d3 * e2 * f1 - a4 * b1 * d0 * e2 * f3 + a4 * b1 * d0 * e3 * f2 +
                   a4 * b1 * d2 * e0 * f3 - a4 * b1 * d2 * e3 * f0 - a4 * b1 * d3 * e0 * f2 + a4 * b1 * d3 * e2 * f0 +
                   a4 * b2 * d0 * e1 * f3 - a4 * b2 * d0 * e3 * f1 - a4 * b2 * d1 * e0 * f3 + a4 * b2 * d1 * e3 * f0 +
                   a4 * b2 * d3 * e0 * f1 - a4 * b2 * d3 * e1 * f0 - a4 * b3 * d0 * e1 * f2 + a4 * b3 * d0 * e2 * f1 +
                   a4 * b3 * d1 * e0 * f2 - a4 * b3 * d1 * e2 * f0 - a4 * b3 * d2 * e0 * f1 + a4 * b3 * d2 * e1 * f0);
        T[3][5] = t363 *
                  (a0 * b1 * c2 * e3 * f4 - a0 * b1 * c2 * e4 * f3 - a0 * b1 * c3 * e2 * f4 + a0 * b1 * c3 * e4 * f2 +
                   a0 * b1 * c4 * e2 * f3 - a0 * b1 * c4 * e3 * f2 - a0 * b2 * c1 * e3 * f4 + a0 * b2 * c1 * e4 * f3 +
                   a0 * b2 * c3 * e1 * f4 - a0 * b2 * c3 * e4 * f1 - a0 * b2 * c4 * e1 * f3 + a0 * b2 * c4 * e3 * f1 +
                   a0 * b3 * c1 * e2 * f4 - a0 * b3 * c1 * e4 * f2 - a0 * b3 * c2 * e1 * f4 + a0 * b3 * c2 * e4 * f1 +
                   a0 * b3 * c4 * e1 * f2 - a0 * b3 * c4 * e2 * f1 - a0 * b4 * c1 * e2 * f3 + a0 * b4 * c1 * e3 * f2 +
                   a0 * b4 * c2 * e1 * f3 - a0 * b4 * c2 * e3 * f1 - a0 * b4 * c3 * e1 * f2 + a0 * b4 * c3 * e2 * f1 -
                   a1 * b0 * c2 * e3 * f4 + a1 * b0 * c2 * e4 * f3 + a1 * b0 * c3 * e2 * f4 - a1 * b0 * c3 * e4 * f2 -
                   a1 * b0 * c4 * e2 * f3 + a1 * b0 * c4 * e3 * f2 + a1 * b2 * c0 * e3 * f4 - a1 * b2 * c0 * e4 * f3 -
                   a1 * b2 * c3 * e0 * f4 + a1 * b2 * c3 * e4 * f0 + a1 * b2 * c4 * e0 * f3 - a1 * b2 * c4 * e3 * f0 -
                   a1 * b3 * c0 * e2 * f4 + a1 * b3 * c0 * e4 * f2 + a1 * b3 * c2 * e0 * f4 - a1 * b3 * c2 * e4 * f0 -
                   a1 * b3 * c4 * e0 * f2 + a1 * b3 * c4 * e2 * f0 + a1 * b4 * c0 * e2 * f3 - a1 * b4 * c0 * e3 * f2 -
                   a1 * b4 * c2 * e0 * f3 + a1 * b4 * c2 * e3 * f0 + a1 * b4 * c3 * e0 * f2 - a1 * b4 * c3 * e2 * f0 +
                   a2 * b0 * c1 * e3 * f4 - a2 * b0 * c1 * e4 * f3 - a2 * b0 * c3 * e1 * f4 + a2 * b0 * c3 * e4 * f1 +
                   a2 * b0 * c4 * e1 * f3 - a2 * b0 * c4 * e3 * f1 - a2 * b1 * c0 * e3 * f4 + a2 * b1 * c0 * e4 * f3 +
                   a2 * b1 * c3 * e0 * f4 - a2 * b1 * c3 * e4 * f0 - a2 * b1 * c4 * e0 * f3 + a2 * b1 * c4 * e3 * f0 +
                   a2 * b3 * c0 * e1 * f4 - a2 * b3 * c0 * e4 * f1 - a2 * b3 * c1 * e0 * f4 + a2 * b3 * c1 * e4 * f0 +
                   a2 * b3 * c4 * e0 * f1 - a2 * b3 * c4 * e1 * f0 - a2 * b4 * c0 * e1 * f3 + a2 * b4 * c0 * e3 * f1 +
                   a2 * b4 * c1 * e0 * f3 - a2 * b4 * c1 * e3 * f0 - a2 * b4 * c3 * e0 * f1 + a2 * b4 * c3 * e1 * f0 -
                   a3 * b0 * c1 * e2 * f4 + a3 * b0 * c1 * e4 * f2 + a3 * b0 * c2 * e1 * f4 - a3 * b0 * c2 * e4 * f1 -
                   a3 * b0 * c4 * e1 * f2 + a3 * b0 * c4 * e2 * f1 + a3 * b1 * c0 * e2 * f4 - a3 * b1 * c0 * e4 * f2 -
                   a3 * b1 * c2 * e0 * f4 + a3 * b1 * c2 * e4 * f0 + a3 * b1 * c4 * e0 * f2 - a3 * b1 * c4 * e2 * f0 -
                   a3 * b2 * c0 * e1 * f4 + a3 * b2 * c0 * e4 * f1 + a3 * b2 * c1 * e0 * f4 - a3 * b2 * c1 * e4 * f0 -
                   a3 * b2 * c4 * e0 * f1 + a3 * b2 * c4 * e1 * f0 + a3 * b4 * c0 * e1 * f2 - a3 * b4 * c0 * e2 * f1 -
                   a3 * b4 * c1 * e0 * f2 + a3 * b4 * c1 * e2 * f0 + a3 * b4 * c2 * e0 * f1 - a3 * b4 * c2 * e1 * f0 +
                   a4 * b0 * c1 * e2 * f3 - a4 * b0 * c1 * e3 * f2 - a4 * b0 * c2 * e1 * f3 + a4 * b0 * c2 * e3 * f1 +
                   a4 * b0 * c3 * e1 * f2 - a4 * b0 * c3 * e2 * f1 - a4 * b1 * c0 * e2 * f3 + a4 * b1 * c0 * e3 * f2 +
                   a4 * b1 * c2 * e0 * f3 - a4 * b1 * c2 * e3 * f0 - a4 * b1 * c3 * e0 * f2 + a4 * b1 * c3 * e2 * f0 +
                   a4 * b2 * c0 * e1 * f3 - a4 * b2 * c0 * e3 * f1 - a4 * b2 * c1 * e0 * f3 + a4 * b2 * c1 * e3 * f0 +
                   a4 * b2 * c3 * e0 * f1 - a4 * b2 * c3 * e1 * f0 - a4 * b3 * c0 * e1 * f2 + a4 * b3 * c0 * e2 * f1 +
                   a4 * b3 * c1 * e0 * f2 - a4 * b3 * c1 * e2 * f0 - a4 * b3 * c2 * e0 * f1 + a4 * b3 * c2 * e1 * f0);
        T[4][5] = -t363 *
                  (a0 * b1 * c2 * d3 * f4 - a0 * b1 * c2 * d4 * f3 - a0 * b1 * c3 * d2 * f4 + a0 * b1 * c3 * d4 * f2 +
                   a0 * b1 * c4 * d2 * f3 - a0 * b1 * c4 * d3 * f2 - a0 * b2 * c1 * d3 * f4 + a0 * b2 * c1 * d4 * f3 +
                   a0 * b2 * c3 * d1 * f4 - a0 * b2 * c3 * d4 * f1 - a0 * b2 * c4 * d1 * f3 + a0 * b2 * c4 * d3 * f1 +
                   a0 * b3 * c1 * d2 * f4 - a0 * b3 * c1 * d4 * f2 - a0 * b3 * c2 * d1 * f4 + a0 * b3 * c2 * d4 * f1 +
                   a0 * b3 * c4 * d1 * f2 - a0 * b3 * c4 * d2 * f1 - a0 * b4 * c1 * d2 * f3 + a0 * b4 * c1 * d3 * f2 +
                   a0 * b4 * c2 * d1 * f3 - a0 * b4 * c2 * d3 * f1 - a0 * b4 * c3 * d1 * f2 + a0 * b4 * c3 * d2 * f1 -
                   a1 * b0 * c2 * d3 * f4 + a1 * b0 * c2 * d4 * f3 + a1 * b0 * c3 * d2 * f4 - a1 * b0 * c3 * d4 * f2 -
                   a1 * b0 * c4 * d2 * f3 + a1 * b0 * c4 * d3 * f2 + a1 * b2 * c0 * d3 * f4 - a1 * b2 * c0 * d4 * f3 -
                   a1 * b2 * c3 * d0 * f4 + a1 * b2 * c3 * d4 * f0 + a1 * b2 * c4 * d0 * f3 - a1 * b2 * c4 * d3 * f0 -
                   a1 * b3 * c0 * d2 * f4 + a1 * b3 * c0 * d4 * f2 + a1 * b3 * c2 * d0 * f4 - a1 * b3 * c2 * d4 * f0 -
                   a1 * b3 * c4 * d0 * f2 + a1 * b3 * c4 * d2 * f0 + a1 * b4 * c0 * d2 * f3 - a1 * b4 * c0 * d3 * f2 -
                   a1 * b4 * c2 * d0 * f3 + a1 * b4 * c2 * d3 * f0 + a1 * b4 * c3 * d0 * f2 - a1 * b4 * c3 * d2 * f0 +
                   a2 * b0 * c1 * d3 * f4 - a2 * b0 * c1 * d4 * f3 - a2 * b0 * c3 * d1 * f4 + a2 * b0 * c3 * d4 * f1 +
                   a2 * b0 * c4 * d1 * f3 - a2 * b0 * c4 * d3 * f1 - a2 * b1 * c0 * d3 * f4 + a2 * b1 * c0 * d4 * f3 +
                   a2 * b1 * c3 * d0 * f4 - a2 * b1 * c3 * d4 * f0 - a2 * b1 * c4 * d0 * f3 + a2 * b1 * c4 * d3 * f0 +
                   a2 * b3 * c0 * d1 * f4 - a2 * b3 * c0 * d4 * f1 - a2 * b3 * c1 * d0 * f4 + a2 * b3 * c1 * d4 * f0 +
                   a2 * b3 * c4 * d0 * f1 - a2 * b3 * c4 * d1 * f0 - a2 * b4 * c0 * d1 * f3 + a2 * b4 * c0 * d3 * f1 +
                   a2 * b4 * c1 * d0 * f3 - a2 * b4 * c1 * d3 * f0 - a2 * b4 * c3 * d0 * f1 + a2 * b4 * c3 * d1 * f0 -
                   a3 * b0 * c1 * d2 * f4 + a3 * b0 * c1 * d4 * f2 + a3 * b0 * c2 * d1 * f4 - a3 * b0 * c2 * d4 * f1 -
                   a3 * b0 * c4 * d1 * f2 + a3 * b0 * c4 * d2 * f1 + a3 * b1 * c0 * d2 * f4 - a3 * b1 * c0 * d4 * f2 -
                   a3 * b1 * c2 * d0 * f4 + a3 * b1 * c2 * d4 * f0 + a3 * b1 * c4 * d0 * f2 - a3 * b1 * c4 * d2 * f0 -
                   a3 * b2 * c0 * d1 * f4 + a3 * b2 * c0 * d4 * f1 + a3 * b2 * c1 * d0 * f4 - a3 * b2 * c1 * d4 * f0 -
                   a3 * b2 * c4 * d0 * f1 + a3 * b2 * c4 * d1 * f0 + a3 * b4 * c0 * d1 * f2 - a3 * b4 * c0 * d2 * f1 -
                   a3 * b4 * c1 * d0 * f2 + a3 * b4 * c1 * d2 * f0 + a3 * b4 * c2 * d0 * f1 - a3 * b4 * c2 * d1 * f0 +
                   a4 * b0 * c1 * d2 * f3 - a4 * b0 * c1 * d3 * f2 - a4 * b0 * c2 * d1 * f3 + a4 * b0 * c2 * d3 * f1 +
                   a4 * b0 * c3 * d1 * f2 - a4 * b0 * c3 * d2 * f1 - a4 * b1 * c0 * d2 * f3 + a4 * b1 * c0 * d3 * f2 +
                   a4 * b1 * c2 * d0 * f3 - a4 * b1 * c2 * d3 * f0 - a4 * b1 * c3 * d0 * f2 + a4 * b1 * c3 * d2 * f0 +
                   a4 * b2 * c0 * d1 * f3 - a4 * b2 * c0 * d3 * f1 - a4 * b2 * c1 * d0 * f3 + a4 * b2 * c1 * d3 * f0 +
                   a4 * b2 * c3 * d0 * f1 - a4 * b2 * c3 * d1 * f0 - a4 * b3 * c0 * d1 * f2 + a4 * b3 * c0 * d2 * f1 +
                   a4 * b3 * c1 * d0 * f2 - a4 * b3 * c1 * d2 * f0 - a4 * b3 * c2 * d0 * f1 + a4 * b3 * c2 * d1 * f0);
        T[5][5] = t363 *
                  (a0 * b1 * c2 * d3 * e4 - a0 * b1 * c2 * d4 * e3 - a0 * b1 * c3 * d2 * e4 + a0 * b1 * c3 * d4 * e2 +
                   a0 * b1 * c4 * d2 * e3 - a0 * b1 * c4 * d3 * e2 - a0 * b2 * c1 * d3 * e4 + a0 * b2 * c1 * d4 * e3 +
                   a0 * b2 * c3 * d1 * e4 - a0 * b2 * c3 * d4 * e1 - a0 * b2 * c4 * d1 * e3 + a0 * b2 * c4 * d3 * e1 +
                   a0 * b3 * c1 * d2 * e4 - a0 * b3 * c1 * d4 * e2 - a0 * b3 * c2 * d1 * e4 + a0 * b3 * c2 * d4 * e1 +
                   a0 * b3 * c4 * d1 * e2 - a0 * b3 * c4 * d2 * e1 - a0 * b4 * c1 * d2 * e3 + a0 * b4 * c1 * d3 * e2 +
                   a0 * b4 * c2 * d1 * e3 - a0 * b4 * c2 * d3 * e1 - a0 * b4 * c3 * d1 * e2 + a0 * b4 * c3 * d2 * e1 -
                   a1 * b0 * c2 * d3 * e4 + a1 * b0 * c2 * d4 * e3 + a1 * b0 * c3 * d2 * e4 - a1 * b0 * c3 * d4 * e2 -
                   a1 * b0 * c4 * d2 * e3 + a1 * b0 * c4 * d3 * e2 + a1 * b2 * c0 * d3 * e4 - a1 * b2 * c0 * d4 * e3 -
                   a1 * b2 * c3 * d0 * e4 + a1 * b2 * c3 * d4 * e0 + a1 * b2 * c4 * d0 * e3 - a1 * b2 * c4 * d3 * e0 -
                   a1 * b3 * c0 * d2 * e4 + a1 * b3 * c0 * d4 * e2 + a1 * b3 * c2 * d0 * e4 - a1 * b3 * c2 * d4 * e0 -
                   a1 * b3 * c4 * d0 * e2 + a1 * b3 * c4 * d2 * e0 + a1 * b4 * c0 * d2 * e3 - a1 * b4 * c0 * d3 * e2 -
                   a1 * b4 * c2 * d0 * e3 + a1 * b4 * c2 * d3 * e0 + a1 * b4 * c3 * d0 * e2 - a1 * b4 * c3 * d2 * e0 +
                   a2 * b0 * c1 * d3 * e4 - a2 * b0 * c1 * d4 * e3 - a2 * b0 * c3 * d1 * e4 + a2 * b0 * c3 * d4 * e1 +
                   a2 * b0 * c4 * d1 * e3 - a2 * b0 * c4 * d3 * e1 - a2 * b1 * c0 * d3 * e4 + a2 * b1 * c0 * d4 * e3 +
                   a2 * b1 * c3 * d0 * e4 - a2 * b1 * c3 * d4 * e0 - a2 * b1 * c4 * d0 * e3 + a2 * b1 * c4 * d3 * e0 +
                   a2 * b3 * c0 * d1 * e4 - a2 * b3 * c0 * d4 * e1 - a2 * b3 * c1 * d0 * e4 + a2 * b3 * c1 * d4 * e0 +
                   a2 * b3 * c4 * d0 * e1 - a2 * b3 * c4 * d1 * e0 - a2 * b4 * c0 * d1 * e3 + a2 * b4 * c0 * d3 * e1 +
                   a2 * b4 * c1 * d0 * e3 - a2 * b4 * c1 * d3 * e0 - a2 * b4 * c3 * d0 * e1 + a2 * b4 * c3 * d1 * e0 -
                   a3 * b0 * c1 * d2 * e4 + a3 * b0 * c1 * d4 * e2 + a3 * b0 * c2 * d1 * e4 - a3 * b0 * c2 * d4 * e1 -
                   a3 * b0 * c4 * d1 * e2 + a3 * b0 * c4 * d2 * e1 + a3 * b1 * c0 * d2 * e4 - a3 * b1 * c0 * d4 * e2 -
                   a3 * b1 * c2 * d0 * e4 + a3 * b1 * c2 * d4 * e0 + a3 * b1 * c4 * d0 * e2 - a3 * b1 * c4 * d2 * e0 -
                   a3 * b2 * c0 * d1 * e4 + a3 * b2 * c0 * d4 * e1 + a3 * b2 * c1 * d0 * e4 - a3 * b2 * c1 * d4 * e0 -
                   a3 * b2 * c4 * d0 * e1 + a3 * b2 * c4 * d1 * e0 + a3 * b4 * c0 * d1 * e2 - a3 * b4 * c0 * d2 * e1 -
                   a3 * b4 * c1 * d0 * e2 + a3 * b4 * c1 * d2 * e0 + a3 * b4 * c2 * d0 * e1 - a3 * b4 * c2 * d1 * e0 +
                   a4 * b0 * c1 * d2 * e3 - a4 * b0 * c1 * d3 * e2 - a4 * b0 * c2 * d1 * e3 + a4 * b0 * c2 * d3 * e1 +
                   a4 * b0 * c3 * d1 * e2 - a4 * b0 * c3 * d2 * e1 - a4 * b1 * c0 * d2 * e3 + a4 * b1 * c0 * d3 * e2 +
                   a4 * b1 * c2 * d0 * e3 - a4 * b1 * c2 * d3 * e0 - a4 * b1 * c3 * d0 * e2 + a4 * b1 * c3 * d2 * e0 +
                   a4 * b2 * c0 * d1 * e3 - a4 * b2 * c0 * d3 * e1 - a4 * b2 * c1 * d0 * e3 + a4 * b2 * c1 * d3 * e0 +
                   a4 * b2 * c3 * d0 * e1 - a4 * b2 * c3 * d1 * e0 - a4 * b3 * c0 * d1 * e2 + a4 * b3 * c0 * d2 * e1 +
                   a4 * b3 * c1 * d0 * e2 - a4 * b3 * c1 * d2 * e0 - a4 * b3 * c2 * d0 * e1 + a4 * b3 * c2 * d1 * e0);
    }
}


#endif //PORTABLE_LINEAR_ALGEBRA_H

#ifdef __cplusplus
}
#endif //__cplusplus
