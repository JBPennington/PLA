
extern "C" {
#include "../pla.h"
}

#include <gtest/gtest.h>
#include <iostream>
#include <map>


class Linear_Algebra_Quat_ : public testing::Test {
    void SetUp() {

    }
    void TearDown() {

    }
public:
    const float tolerance = 0.00001f;
    mat2x2 Mat2_zero {{0,0},{0,0}};
    mat3x3 Mat3_zero {{0,0,0},{0,0,0},{0,0,0}};
    mat4x4 Mat4_zero {{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}};
    mat5x5 Mat5_zero {{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0}};
    mat6x6 Mat6_zero {{0,0,0,0,0,0},{0,0,0,0,0,0},{0,0,0,0,0,0},{0,0,0,0,0,0},{0,0,0,0,0,0},{0,0,0,0,0,0}};
    mat2x2 Mat2A     {{0,1},{2,3}};
    mat3x3 Mat3A     {{0,1,2},{3,4,5},{6,7,8}};
    mat4x4 Mat4A     {{0,1,2,3},{4,5,6,7},{8,9,10,11},{12,13,14,15}};
    mat5x5 Mat5A     {{0,1,2,3,4},{5,6,7,8,9},{10,11,12,13,14},{15,16,17,18,19},{20,21,22,23,24}};
    mat6x6 Mat6A     {{0,1,2,3,4,5},{6,7,8,9,10,11},{12,13,14,15,16,17},{18,19,20,21,22,23},{24,25,26,27,28,29},{30,31,32,33,34,35}};
    mat2x2 Mat2B     {{0,1},{2,3}};
    mat3x3 Mat3B     {{0,1,2},{3,4,5},{6,7,8}};
    mat4x4 Mat4B     {{0,1,2,3},{4,5,6,7},{8,9,10,11},{12,13,14,15}};
    mat5x5 Mat5B     {{0,1,2,3,4},{5,6,7,8,9},{10,11,12,13,14},{15,16,17,18,19},{20,21,22,23,24}};
    mat6x6 Mat6B     {{ 0, 1, 2, 3, 4, 5},
                      { 6, 7, 8, 9,10,11},
                      {12,13,14,15,16,17},
                      {18,19,20,21,22,23},
                      {24,25,26,27,28,29},
                      {30,31,32,33,34,35}};
    quat QuatA        {0.7071f, 0.0f, 0.0f, 0.7071f}; // Represents a rotation of PI/2 about x
    quat QuatB        {0.7071f, 0.0f, 0.0f, 0.7071f};
};


TEST_F(Linear_Algebra_Quat_, Quat_Identity) {
    set_quat_identity(QuatA);

    quat expect = IDENTITY_QUAT;

    uint32_t iter;
    for (iter=0; iter<4; ++iter) {
        EXPECT_EQ(expect[iter], QuatA[iter]) << "Index: " << iter;
    }
}


TEST_F(Linear_Algebra_Quat_, Quat_Add_New) {
    quat result = IDENTITY_QUAT;

    quat_add(result, QuatA, QuatB);

    quat expect = {1.4142f, 0.0f, 0.0f, 1.4142f};

    uint32_t iter;
    for (iter=0; iter<4; ++iter) {
        EXPECT_EQ(expect[iter], result[iter]) << "Index: " << iter;
    }
}


TEST_F(Linear_Algebra_Quat_, Quat_Sub_New) {
    quat result = IDENTITY_QUAT;

    quat_sub(result, QuatA, QuatB);

    quat expect = {0.0f, 0.0f, 0.0f, 0.0f};

    uint32_t iter;
    for (iter=0; iter<4; ++iter) {
        EXPECT_EQ(expect[iter], result[iter]) << "Index: " << iter;
    }
}


TEST_F(Linear_Algebra_Quat_, Quat_Mul_New) {
    quat result = IDENTITY_QUAT;

    quat_mul(result, QuatA, QuatB);

    quat expect = {0.99999f, 0.0f, 0.0f, 0.0f};

    uint32_t iter;
    for (iter=0; iter<4; ++iter) {
        EXPECT_NEAR(expect[iter], result[iter], tolerance) << "Index: " << iter;
    }
}


TEST_F(Linear_Algebra_Quat_, Quat_Scale_New) {
    quat result = IDENTITY_QUAT;
    float scale = 0.5f;

    quat_scale(result, QuatA, scale);

    quat expect = {0.35355f, 0.0f, 0.0f, 0.35355f};

    uint32_t iter;
    for (iter=0; iter<4; ++iter) {
        EXPECT_NEAR(expect[iter], result[iter], tolerance) << "Index: " << iter;
    }
}


TEST_F(Linear_Algebra_Quat_, Quat_Mult_Inner_New) {
    float result = quat_inner_product(QuatA, QuatB);
    float expect = 0.99999f;

    EXPECT_NEAR(expect, result, tolerance);
}


TEST_F(Linear_Algebra_Quat_, Quat_Conjugate_New) {
    quat result = IDENTITY_QUAT;

    quat_conj(result, QuatA);

    quat expect = {-0.7071f, 0.0f, 0.0f, 0.7071f};

    uint32_t iter;
    for (iter=0; iter<4; ++iter) {
        EXPECT_NEAR(expect[iter], result[iter], tolerance) << "Index: " << iter;
    }
}


TEST_F(Linear_Algebra_Quat_, Quat_From_Angle_Axis_New) {
    quat result = IDENTITY_QUAT;
    quat_copy(result, QuatA);

    float rotation_angle = static_cast<float>(M_PI_4);
    vec3  rotation_axis  = {0.0f,1.0f,0.0f};

    // Test Rotation about Y axis
    quat_from_angle_axis(result, rotation_angle, rotation_axis);

    quat expect = {0.0f, 0.38268f, 0.0f, 0.9238795};

    uint32_t iter;
    for (iter=0; iter<4; ++iter) {
        EXPECT_NEAR(expect[iter], result[iter], tolerance) << "Index: " << iter;
    }
}


TEST_F(Linear_Algebra_Quat_, Quat_Norm) {
    float result = quat_norm(QuatA);

    float expect = 0.99999f;

    EXPECT_NEAR(expect, result, tolerance);
}


TEST_F(Linear_Algebra_Quat_, Quat_Normalize_New) {
    quat result = IDENTITY_QUAT;

    quat_normalize_n(result, QuatA);

    quat expect = {0.7071f, 0.0f, 0.0f, 0.7071f};

    uint32_t iter;
    for (iter=0; iter<4; ++iter) {
        EXPECT_NEAR(expect[iter], result[iter], tolerance) << "Index: " << iter;
    }
}


TEST_F(Linear_Algebra_Quat_, Quat_Normalize) {
    quat_normalize(QuatA);

    quat expect = {0.7071f, 0.0f, 0.0f, 0.7071f};

    uint32_t iter;
    for (iter=0; iter<4; ++iter) {
        EXPECT_NEAR(expect[iter], QuatA[iter], tolerance) << "Index: " << iter;
    }
}


TEST_F(Linear_Algebra_Quat_, Quat_Mul_Vec3_New) {
    vec3 result = EMPTYVEC3;
    vec3 axis = {0.0f, 0.5f, 0.0f};

    quat_mul_vec3(result, QuatA, axis);

    vec3 expect = {0.0f, 0.0f, 0.5f};

    uint32_t iter;
    for (iter=0; iter<3; ++iter) {
        EXPECT_NEAR(expect[iter], result[iter], tolerance) << "Index: " << iter;
    }
}


TEST_F(Linear_Algebra_Quat_, Matrix_From_Quat) {
    mat4x4 result = IDENTITY4x4;

    mat4x4_from_quat(result, QuatA);

    mat4x4 expect = {
        {  1.0f,  0.0f,  0.0f,  0.0f,},
        {  0.0f,  0.0f, -1.0f,  0.0f,},
        {  0.0f,  1.0f,  0.0f,  0.0f,},
        {  0.0f,  0.0f,  0.0f,  1.0f,},
    };

    uint32_t row, col;
    for (row=0; row<4; ++row) {
        for (col=0; col<4; ++col) {
            EXPECT_NEAR(expect[row][col], result[row][col], 1e2*tolerance)
                            << "Row: " << row << " Col: " << col;
        }
    }
}


TEST_F(Linear_Algebra_Quat_, Quat_From_Matrix) {
    quat result = IDENTITY_QUAT;
    mat4x4 test = {
        {  0.0f,  0.0f, -1.0f,  0.0f,},
        {  0.0f,  1.0f,  0.0f,  0.0f,},
        {  1.0f,  0.0f,  0.0f,  0.0f,},
        {  0.0f,  0.0f,  0.0f,  1.0f,},
    };

    quat_from_mat4x4(result, test);

    quat expect = {0, -0.7071068f, 0, 0.7071068f};

    uint32_t iter;
    for (iter=0; iter<4; ++iter) {
        EXPECT_NEAR(expect[iter], result[iter], tolerance) << "Index: " << iter;
    }
}
