
extern "C" {
#include "../linear_algebra.h"
}

#include <gtest/gtest.h>
#include <iostream>


class Linear_Algebra_Scale_ : public testing::Test {
    void SetUp() {}
    void TearDown() {}

public:
    const float tolerance = 0.00000001f;
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
    mat6x6 Mat6B     {{0,1,2,3,4,5},{6,7,8,9,10,11},{12,13,14,15,16,17},{18,19,20,21,22,23},{24,25,26,27,28,29},{30,31,32,33,34,35}};
};

TEST_F(Linear_Algebra_Scale_, Scale_Two_Vec2_Test) {
    vec2    a = {0, 1};
    float   b = 7;

    vec2 result = {0, 0};
    vec2_scale_n(result, a, b);

    vec2 expect = {0, 7};

    uint32_t iter;
    for (iter=0; iter<2; ++iter) {
        EXPECT_EQ(expect[iter], result[iter]);
    }
}

TEST_F(Linear_Algebra_Scale_, Scale_Two_Vec3_Test) {
    vec3    a = {0, 1, 2};
    float   b = 7;

    vec3 result = {0, 0, 0};
    vec3_scale_n(result, a, b);

    vec3 expect = {0, 7, 14};

    uint32_t iter;
    for (iter=0; iter<3; ++iter) {
        EXPECT_EQ(expect[iter], result[iter]);
    }
}

TEST_F(Linear_Algebra_Scale_, Scale_Two_Vec4_Test) {
    vec4    a = {0, 1, 2, 3};
    float   b = 7;

    vec4 result = {0, 0, 0, 0};
    vec4_scale_n(result, a, b);

    vec4 expect = {0, 7, 14, 21};

    uint32_t iter;
    for (iter=0; iter<4; ++iter) {
        EXPECT_EQ(expect[iter], result[iter]);
    }
}

TEST_F(Linear_Algebra_Scale_, Scale_Two_Vec5_Test) {
    vec5    a = {0, 1, 2, 3, 4};
    float   b = 7;

    vec5 result = {0, 0, 0, 0, 0};
    vec5_scale_n(result, a, b);

    vec5 expect = {0, 7, 14, 21, 28};

    uint32_t iter;
    for (iter=0; iter<5; ++iter) {
        EXPECT_EQ(expect[iter], result[iter]);
    }
}

TEST_F(Linear_Algebra_Scale_, Scale_Two_Vec6_Test) {
    vec6    a = {0, 1, 2, 3, 4, 5};
    float   b = 7;

    vec6 result = {0, 0, 0, 0, 0, 0};
    vec6_scale_n(result, a, b);

    vec6 expect = {0, 7, 14, 21, 28, 35};

    uint32_t iter;
    for (iter=0; iter<6; ++iter) {
        EXPECT_EQ(expect[iter], result[iter]);
    }
}

TEST_F(Linear_Algebra_Scale_, Scale_Two_Vec2_In_Place_Test) {
    vec2    a = {0, 1};
    float   b = 7;

    vec2_scale(a, b);

    vec2 expect = {0, 7};
    uint32_t iter;
    for (iter=0; iter<2; ++iter) {
        EXPECT_EQ(expect[iter], a[iter]);
    }
}

TEST_F(Linear_Algebra_Scale_, Scale_Two_Vec3_In_Place_Test) {
    vec3    a = {0, 1, 2};
    float   b = 7;

    vec3_scale(a, b);

    vec3 expect = {0, 7, 14};
    uint32_t iter;
    for (iter=0; iter<3; ++iter) {
        EXPECT_EQ(expect[iter], a[iter]);
    }
}

TEST_F(Linear_Algebra_Scale_, Scale_Two_Vec4_In_Place_Test) {
    vec4    a = {0, 1, 2, 3};
    float   b = 7;

    vec4_scale(a, b);

    vec4 expect = {0, 7, 14, 21};
    uint32_t iter;
    for (iter=0; iter<4; ++iter) {
        EXPECT_EQ(expect[iter], a[iter]);
    }
}

TEST_F(Linear_Algebra_Scale_, Scale_Two_Vec5_In_Place_Test) {
    vec5    a = {0, 1, 2, 3, 4};
    float   b = 7;

    vec5_scale(a, b);

    vec5 expect = {0, 7, 14, 21, 28};
    uint32_t iter;
    for (iter=0; iter<5; ++iter) {
        EXPECT_EQ(expect[iter], a[iter]);
    }
}

TEST_F(Linear_Algebra_Scale_, Scale_Two_Vec6_In_Place_Test) {
    vec6    a = {0, 1, 2, 3, 4, 5};
    float   b = 7;

    vec6_scale(a, b);

    vec6 expect = {0, 7, 14, 21, 28, 35};
    uint32_t iter;
    for (iter=0; iter<6; ++iter) {
        EXPECT_EQ(expect[iter], a[iter]);
    }
}

TEST_F(Linear_Algebra_Scale_, Scale_Return_Two_Vec2_Test) {
    vec2    a = {0, 1};
    float   b = 7;

    float * result = vec2_scale_r(a, b);

    vec2 expect = {0, 7};

    uint32_t iter;
    for (iter=0; iter<2; ++iter) {
        EXPECT_EQ(expect[iter], result[iter]);
    }
}

TEST_F(Linear_Algebra_Scale_, Scale_Return_Two_Vec3_Test) {
    vec3    a = {0, 1, 2};
    float   b = 7;

    float * result = vec3_scale_r(a, b);

    vec3 expect = {0, 7, 14};

    uint32_t iter;
    for (iter=0; iter<3; ++iter) {
        EXPECT_EQ(expect[iter], result[iter]);
    }
}

TEST_F(Linear_Algebra_Scale_, Scale_Return_Two_Vec4_Test) {
    vec4    a = {0, 1, 2, 3};
    float   b = 7;

    float * result = vec4_scale_r(a, b);

    vec4 expect = {0, 7, 14, 21};

    uint32_t iter;
    for (iter=0; iter<4; ++iter) {
        EXPECT_EQ(expect[iter], result[iter]);
    }
}

TEST_F(Linear_Algebra_Scale_, Scale_Return_Two_Vec5_Test) {
    vec5    a = {0, 1, 2, 3, 4};
    float   b = 7;

    float * result = vec5_scale_r(a, b);

    vec5 expect = {0, 7, 14, 21, 28};

    uint32_t iter;
    for (iter=0; iter<5; ++iter) {
        EXPECT_EQ(expect[iter], result[iter]);
    }
}

TEST_F(Linear_Algebra_Scale_, Scale_Return_Two_Vec6_Test) {
    vec6    a = {0, 1, 2, 3, 4, 5};
    float   b = 7;

    float * result = vec6_scale_r(a, b);

    vec6 expect = {0, 7, 14, 21, 28, 35};

    uint32_t iter;
    for (iter=0; iter<6; ++iter) {
        EXPECT_EQ(expect[iter], result[iter]);
    }
}


TEST_F(Linear_Algebra_Scale_, Two_Through_Six_Square_Mat_Scale_to_New_Mat_Test) {

    mat2x2 Mat2_expected {{0,.5},{1,1.5}};
    mat3x3 Mat3_expected {{0,.5,1},{1.5,2,2.5},{3,3.5,4}};
    mat4x4 Mat4_expected {{0,.5,1,1.5},{2,2.5,3,3.5},{4,4.5,5,5.5},{6,6.5,7,7.5}};
    mat5x5 Mat5_expected {{0,.5,1,1.5,2},{2.5,3,3.5,4,4.5},{5,5.5,6,6.5,7},{7.5,8,8.5,9,9.5},{10,10.5,11,11.5,12}};
    mat6x6 Mat6_expected {{0,.5,1,1.5,2,2.5},{3,3.5,4,4.5,5,5.5},{6,6.5,7,7.5,8,8.5},{9,9.5,10,10.5,11,11.5},{12,12.5,13,13.5,14,14.5},{15,15.5,16,16.5,17,17.5}};

    mat2x2 Mat2_result; mat2x2_copy(Mat2_result, Mat2_zero);
    mat3x3 Mat3_result; mat3x3_copy(Mat3_result, Mat3_zero);
    mat4x4 Mat4_result; mat4x4_copy(Mat4_result, Mat4_zero);
    mat5x5 Mat5_result; mat5x5_copy(Mat5_result, Mat5_zero);
    mat6x6 Mat6_result; mat6x6_copy(Mat6_result, Mat6_zero);

    for (std::size_t size = 2; size<7; ++size) {
        if (size == 2) mat2x2_scale_n(Mat2_result, Mat2A, 0.5f);
        if (size == 3) mat3x3_scale_n(Mat3_result, Mat3A, 0.5f);
        if (size == 4) mat4x4_scale_n(Mat4_result, Mat4A, 0.5f);
        if (size == 5) mat5x5_scale_n(Mat5_result, Mat5A, 0.5f);
        if (size == 6) mat6x6_scale_n(Mat6_result, Mat6A, 0.5f);
        for (std::size_t row = 0; row < size; ++row) {
            for (std::size_t col = 0; col < size; ++col) {
                if (size == 2) EXPECT_NEAR(Mat2_expected[row][col], Mat2_result[row][col], tolerance) << "Matrix: " << size << " Row: " << row << " Col: " << col;
                if (size == 3) EXPECT_NEAR(Mat3_expected[row][col], Mat3_result[row][col], tolerance) << "Matrix: " << size << " Row: " << row << " Col: " << col;
                if (size == 4) EXPECT_NEAR(Mat4_expected[row][col], Mat4_result[row][col], tolerance) << "Matrix: " << size << " Row: " << row << " Col: " << col;
                if (size == 5) EXPECT_NEAR(Mat5_expected[row][col], Mat5_result[row][col], tolerance) << "Matrix: " << size << " Row: " << row << " Col: " << col;
                if (size == 6) EXPECT_NEAR(Mat6_expected[row][col], Mat6_result[row][col], tolerance) << "Matrix: " << size << " Row: " << row << " Col: " << col;
            }
        }
    }
}

TEST_F(Linear_Algebra_Scale_, Two_Through_Six_Square_Mat_Scale_In_Place_Test) {

    mat2x2 Mat2_expected {{0,.5},{1,1.5}};
    mat3x3 Mat3_expected {{0,.5,1},{1.5,2,2.5},{3,3.5,4}};
    mat4x4 Mat4_expected {{0,.5,1,1.5},{2,2.5,3,3.5},{4,4.5,5,5.5},{6,6.5,7,7.5}};
    mat5x5 Mat5_expected {{0,.5,1,1.5,2},{2.5,3,3.5,4,4.5},{5,5.5,6,6.5,7},{7.5,8,8.5,9,9.5},{10,10.5,11,11.5,12}};
    mat6x6 Mat6_expected {{0,.5,1,1.5,2,2.5},{3,3.5,4,4.5,5,5.5},{6,6.5,7,7.5,8,8.5},{9,9.5,10,10.5,11,11.5},{12,12.5,13,13.5,14,14.5},{15,15.5,16,16.5,17,17.5}};

    mat2x2 Mat2_result; mat2x2_copy(Mat2_result, Mat2A);
    mat3x3 Mat3_result; mat3x3_copy(Mat3_result, Mat3A);
    mat4x4 Mat4_result; mat4x4_copy(Mat4_result, Mat4A);
    mat5x5 Mat5_result; mat5x5_copy(Mat5_result, Mat5A);
    mat6x6 Mat6_result; mat6x6_copy(Mat6_result, Mat6A);

    for (std::size_t size = 2; size<7; ++size) {
        if (size == 2) mat2x2_scale(Mat2_result, 0.5f);
        if (size == 3) mat3x3_scale(Mat3_result, 0.5f);
        if (size == 4) mat4x4_scale(Mat4_result, 0.5f);
        if (size == 5) mat5x5_scale(Mat5_result, 0.5f);
        if (size == 6) mat6x6_scale(Mat6_result, 0.5f);
        for (std::size_t row = 0; row < size; ++row) {
            for (std::size_t col = 0; col < size; ++col) {
                if (size == 2) EXPECT_NEAR(Mat2_expected[row][col], Mat2_result[row][col], tolerance) << "Matrix: " << size << " Row: " << row << " Col: " << col;
                if (size == 3) EXPECT_NEAR(Mat3_expected[row][col], Mat3_result[row][col], tolerance) << "Matrix: " << size << " Row: " << row << " Col: " << col;
                if (size == 4) EXPECT_NEAR(Mat4_expected[row][col], Mat4_result[row][col], tolerance) << "Matrix: " << size << " Row: " << row << " Col: " << col;
                if (size == 5) EXPECT_NEAR(Mat5_expected[row][col], Mat5_result[row][col], tolerance) << "Matrix: " << size << " Row: " << row << " Col: " << col;
                if (size == 6) EXPECT_NEAR(Mat6_expected[row][col], Mat6_result[row][col], tolerance) << "Matrix: " << size << " Row: " << row << " Col: " << col;
            }
        }
    }
}

TEST_F(Linear_Algebra_Scale_, Two_Through_Six_Square_Mat_Scale_In_Place_And_Return_Test) {

    mat2x2 Mat2_expected {{0,.5},{1,1.5}};
    mat3x3 Mat3_expected {{0,.5,1},{1.5,2,2.5},{3,3.5,4}};
    mat4x4 Mat4_expected {{0,.5,1,1.5},{2,2.5,3,3.5},{4,4.5,5,5.5},{6,6.5,7,7.5}};
    mat5x5 Mat5_expected {{0,.5,1,1.5,2},{2.5,3,3.5,4,4.5},{5,5.5,6,6.5,7},{7.5,8,8.5,9,9.5},{10,10.5,11,11.5,12}};
    mat6x6 Mat6_expected {{0,.5,1,1.5,2,2.5},{3,3.5,4,4.5,5,5.5},{6,6.5,7,7.5,8,8.5},{9,9.5,10,10.5,11,11.5},{12,12.5,13,13.5,14,14.5},{15,15.5,16,16.5,17,17.5}};

    vec2 * Mat2_result = mat2x2_scale_r(Mat2A, 0.5f);
    vec3 * Mat3_result = mat3x3_scale_r(Mat3A, 0.5f);
    vec4 * Mat4_result = mat4x4_scale_r(Mat4A, 0.5f);
    vec5 * Mat5_result = mat5x5_scale_r(Mat5A, 0.5f);
    vec6 * Mat6_result = mat6x6_scale_r(Mat6A, 0.5f);

    for (std::size_t size = 2; size<7; ++size) {
        for (std::size_t row = 0; row < size; ++row) {
            for (std::size_t col = 0; col < size; ++col) {
                if (size == 2) EXPECT_NEAR(Mat2_expected[row][col], Mat2_result[row][col], tolerance) << "Matrix: " << size << " Row: " << row << " Col: " << col;
                if (size == 3) EXPECT_NEAR(Mat3_expected[row][col], Mat3_result[row][col], tolerance) << "Matrix: " << size << " Row: " << row << " Col: " << col;
                if (size == 4) EXPECT_NEAR(Mat4_expected[row][col], Mat4_result[row][col], tolerance) << "Matrix: " << size << " Row: " << row << " Col: " << col;
                if (size == 5) EXPECT_NEAR(Mat5_expected[row][col], Mat5_result[row][col], tolerance) << "Matrix: " << size << " Row: " << row << " Col: " << col;
                if (size == 6) EXPECT_NEAR(Mat6_expected[row][col], Mat6_result[row][col], tolerance) << "Matrix: " << size << " Row: " << row << " Col: " << col;
            }
        }
    }
}

