
extern "C" {
#include "../pla.h"
}

#include <gtest/gtest.h>
#include <iostream>
#include <map>


static const float PI_2 = static_cast<float>(M_PI_2);


class Linear_Algebra_Misc_ : public testing::Test {
    void SetUp() {

    }
    void TearDown() {

    }
public:
    const float tolerance = 0.0000001f;
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
};

TEST_F(Linear_Algebra_Misc_, Vec3_Cross_Multiplication_Set_New_Value_Test) {
    vec3 a = {0, 1, 2};
    vec3 b = {3, 4, 5};

    vec3 result = {0,0,0};
    vec3_mul_cross_n(result, a, b);

    vec3 expect = {-3, 6, -3};

    uint32_t iter;
    for (iter=0; iter<3; ++iter) {
        EXPECT_EQ(expect[iter], result[iter]);
    }
}

TEST_F(Linear_Algebra_Misc_, Vec3_Cross_Multiplication_Set_First_Test) {
    vec3 a = {0, 1, 2};
    vec3 b = {3, 4, 5};

    vec3_mul_cross(a, b);

    vec3 expect = {-3, 6, -3};

    uint32_t iter;
    for (iter=0; iter<3; ++iter) {
        EXPECT_EQ(expect[iter], a[iter]);
    }
}

TEST_F(Linear_Algebra_Misc_, Vec3_Cross_Multiplication_Set_First_And_Return_Test) {
    vec3 a = {0, 1, 2};
    vec3 b = {3, 4, 5};

    float * result = vec3_mul_cross_r(a, b);

    vec3 expect = {-3, 6, -3};

    uint32_t iter;
    for (iter=0; iter<3; ++iter) {
        EXPECT_EQ(expect[iter], result[iter]);
    }
}

TEST_F(Linear_Algebra_Misc_, Vec3_Reflect_Set_New_Value_Test) {
    vec3 a = {1, 0, 0};
    vec3 n = {1, 1, 1};

    vec3 result = {0, 0, 0};
    vec3_reflect_n(result, a, n);

    vec3 expect = {-1, -2, -2};

    uint32_t iter;
    for (iter=0; iter<3; ++iter) {
        EXPECT_NEAR(expect[iter], result[iter], tolerance);
    }
}

TEST_F(Linear_Algebra_Misc_, Vec3_Reflect_Set_First_Test) {
    vec3 a = {1, 0, 0};
    vec3 n = {1, 1, 1};

    vec3_reflect(a, n);

    vec3 expect = {-1, -2, -2};

    uint32_t iter;
    for (iter=0; iter<3; ++iter) {
        EXPECT_NEAR(expect[iter], a[iter], tolerance);
    }
}

TEST_F(Linear_Algebra_Misc_, Vec3_Reflect_Set_First_And_Return_Test) {
    vec3 a = {1, 0, 0};
    vec3 n = {1, 1, 1};

    float * result = vec3_reflect_r(a, n);

    vec3 expect = {-1, -2, -2};

    uint32_t iter;
    for (iter=0; iter<3; ++iter) {
        EXPECT_NEAR(expect[iter], result[iter], tolerance);
    }
}

TEST_F(Linear_Algebra_Misc_, Vec4_Cross_Multiplication_Set_New_Value_Test) {
    vec4 a = {0, 1, 2, 3};
    vec4 b = {4, 5, 6, 7};

    vec4 result = {0,0,0,0};
    vec4_mul_cross_n(result, a, b);

    vec4 expect = {-4, 8, -4, 1};

    uint32_t iter;
    for (iter=0; iter<4; ++iter) {
        EXPECT_EQ(expect[iter], result[iter]);
    }
}

TEST_F(Linear_Algebra_Misc_, Vec4_Cross_Multiplication_Set_First_Test) {
    vec4 a = {0, 1, 2, 3};
    vec4 b = {4, 5, 6, 7};

    vec4_mul_cross(a, b);

    vec4 expect = {-4, 8, -4, 1};

    uint32_t iter;
    for (iter=0; iter<4; ++iter) {
        EXPECT_EQ(expect[iter], a[iter]);
    }
}

TEST_F(Linear_Algebra_Misc_, Vec4_Cross_Multiplication_Set_First_And_Return_Test) {
    vec4 a = {0, 1, 2, 3};
    vec4 b = {4, 5, 6, 7};

    float * result = vec4_mul_cross_r(a, b);

    vec4 expect = {-4, 8, -4, 1};

    uint32_t iter;
    for (iter=0; iter<4; ++iter) {
        EXPECT_EQ(expect[iter], result[iter]);
    }
}

TEST_F(Linear_Algebra_Misc_, Vec4_Reflect_Set_New_Value_Test) {
    vec4 a = {0, 1, 2, 3};
    vec4 n = {1, 1, 1, 1};

    vec4 result = {0, 0, 0, 0};
    vec4_reflect_n(result, a, n);

    vec4 expect = {-12, -11, -10, -9};

    uint32_t iter;
    for (iter=0; iter<4; ++iter) {
        EXPECT_NEAR(expect[iter], result[iter], tolerance);
    }
}

TEST_F(Linear_Algebra_Misc_, Vec4_Reflect_Set_First_Test) {
    vec4 a = {0, 1, 2, 3};
    vec4 n = {1, 1, 1, 1};

    vec4_reflect(a, n);

    vec4 expect = {-12, -11, -10, -9};

    uint32_t iter;
    for (iter=0; iter<4; ++iter) {
        EXPECT_NEAR(expect[iter], a[iter], tolerance);
    }
}

TEST_F(Linear_Algebra_Misc_, Vec4_Reflect_Set_First_And_Return_Test) {
    vec4 a = {0, 1, 2, 3};
    vec4 n = {1, 1, 1, 1};

    float * result = vec4_reflect_r(a, n);

    vec4 expect = {-12, -11, -10, -9};

    uint32_t iter;
    for (iter=0; iter<4; ++iter) {
        EXPECT_NEAR(expect[iter], result[iter], tolerance);
    }
}


TEST_F(Linear_Algebra_Misc_, Vec4_Add_In_Place_And_Return_Cross_Product_Test) {
    vec4 a = { 0,  1,  2,  3};
    vec4 b = { 4,  5,  6,  7};
    vec4 c = { 8,  9, 10, 11};
    vec4 d = {12, 13, 14, 15};

    float * result = vec4_mul_cross_r(vec4_add_r(a, b), vec4_sub_r(d,c));

    vec4 expect = {-8, 16, -8, 1};
    uint32_t iter;
    for (iter=0; iter<4; ++iter) {
        EXPECT_NEAR(expect[iter], result[iter], tolerance);
    }
}

TEST_F(Linear_Algebra_Misc_, Vec4_Add_And_Cross_Product_Test) {
    vec4 a = { 0,  1,  2,  3};
    vec4 b = { 4,  5,  6,  7};
    vec4 c = { 8,  9, 10, 11};
    vec4 d = {12, 13, 14, 15};

    vec4 temp1 = {0, 0, 0, 0};
    vec4 temp2 = {0, 0, 0, 0};
    vec4 temp3 = {0, 0, 0, 0};
    vec4_add_n(temp1, a, b);
    vec4_sub_n(temp2, d, c);
    vec4_mul_cross_n(temp3, temp1, temp2);

    vec4 expect = {-8, 16, -8, 1};
    uint32_t iter;
    for (iter=0; iter<4; ++iter) {
        EXPECT_NEAR(expect[iter], temp3[iter], tolerance);
    }
}

TEST_F(Linear_Algebra_Misc_, Two_Through_Six_Square_Mat_Initialization_Test) {

    // Creation Pushed to Setup now, not sure how I feel about that

    for (std::size_t size = 2; size<7; ++size) {
        for (std::size_t row = 0; row < size; ++row) {
            for (std::size_t col = 0; col < size; ++col) {
                if (size == 2) EXPECT_NEAR(0, Mat2_zero[row][col], tolerance) << "Matrix: " << size << " Row: " << row << " Col: " << col;
                if (size == 3) EXPECT_NEAR(0, Mat3_zero[row][col], tolerance) << "Matrix: " << size << " Row: " << row << " Col: " << col;
                if (size == 4) EXPECT_NEAR(0, Mat4_zero[row][col], tolerance) << "Matrix: " << size << " Row: " << row << " Col: " << col;
                if (size == 5) EXPECT_NEAR(0, Mat5_zero[row][col], tolerance) << "Matrix: " << size << " Row: " << row << " Col: " << col;
                if (size == 6) EXPECT_NEAR(0, Mat6_zero[row][col], tolerance) << "Matrix: " << size << " Row: " << row << " Col: " << col;
            }
        }
    }
}


TEST_F(Linear_Algebra_Misc_, Two_Through_Six_Square_Mat_Identity_Test) {

    mat2x2 Mat2_expected {{1,0},{0,1}};
    mat3x3 Mat3_expected {{1,0,0},{0,1,0},{0,0,1}};
    mat4x4 Mat4_expected {{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}};
    mat5x5 Mat5_expected {{1,0,0,0,0},{0,1,0,0,0},{0,0,1,0,0},{0,0,0,1,0},{0,0,0,0,1}};
    mat6x6 Mat6_expected {{1,0,0,0,0,0},{0,1,0,0,0,0},{0,0,1,0,0,0},{0,0,0,1,0,0},{0,0,0,0,1,0},{0,0,0,0,0,1}};

    mat2x2 Mat2_result; mat2x2_identity(Mat2_result);
    mat3x3 Mat3_result; mat3x3_identity(Mat3_result);
    mat4x4 Mat4_result; mat4x4_identity(Mat4_result);
    mat5x5 Mat5_result; mat5x5_identity(Mat5_result);
    mat6x6 Mat6_result; mat6x6_identity(Mat6_result);

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

TEST_F(Linear_Algebra_Misc_, Two_Through_Six_Square_Mat_Identity_Return_Test) {

    mat2x2 Mat2_expected {{1,0},{0,1}};
    mat3x3 Mat3_expected {{1,0,0},{0,1,0},{0,0,1}};
    mat4x4 Mat4_expected {{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}};
    mat5x5 Mat5_expected {{1,0,0,0,0},{0,1,0,0,0},{0,0,1,0,0},{0,0,0,1,0},{0,0,0,0,1}};
    mat6x6 Mat6_expected {{1,0,0,0,0,0},{0,1,0,0,0,0},{0,0,1,0,0,0},{0,0,0,1,0,0},{0,0,0,0,1,0},{0,0,0,0,0,1}};

    vec2* Mat2_result = mat2x2_identity_r(Mat2_zero);
    vec3* Mat3_result = mat3x3_identity_r(Mat3_zero);
    vec4* Mat4_result = mat4x4_identity_r(Mat4_zero);
    vec5* Mat5_result = mat5x5_identity_r(Mat5_zero);
    vec6* Mat6_result = mat6x6_identity_r(Mat6_zero);

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

TEST_F(Linear_Algebra_Misc_, Two_Through_Six_Square_Mat_Row_Test) {

    mat2x2 Mat2_expected {{1,0},{0,1}};
    mat3x3 Mat3_expected {{1,0,0},{0,1,0},{0,0,1}};
    mat4x4 Mat4_expected {{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}};
    mat5x5 Mat5_expected {{1,0,0,0,0},{0,1,0,0,0},{0,0,1,0,0},{0,0,0,1,0},{0,0,0,0,1}};
    mat6x6 Mat6_expected {{1,0,0,0,0,0},{0,1,0,0,0,0},{0,0,1,0,0,0},{0,0,0,1,0,0},{0,0,0,0,1,0},{0,0,0,0,0,1}};

    mat2x2 Mat2_result; mat2x2_identity(Mat2_result);
    mat3x3 Mat3_result; mat3x3_identity(Mat3_result);
    mat4x4 Mat4_result; mat4x4_identity(Mat4_result);
    mat5x5 Mat5_result; mat5x5_identity(Mat5_result);
    mat6x6 Mat6_result; mat6x6_identity(Mat6_result);

    vec2 vec2_result; mat2x2_row(vec2_result, Mat2_result, 1);
    vec3 vec3_result; mat3x3_row(vec3_result, Mat3_result, 1);
    vec4 vec4_result; mat4x4_row(vec4_result, Mat4_result, 1);
    vec5 vec5_result; mat5x5_row(vec5_result, Mat5_result, 1);
    vec6 vec6_result; mat6x6_row(vec6_result, Mat6_result, 1);

    for (std::size_t size = 2; size<7; ++size) {
        for (std::size_t col = 0; col < size; ++col) {
            if (size == 2) EXPECT_NEAR(Mat2_expected[1][col], vec2_result[col], tolerance) << "Matrix: " << size << " Col: " << col;
            if (size == 3) EXPECT_NEAR(Mat3_expected[1][col], vec3_result[col], tolerance) << "Matrix: " << size << " Col: " << col;
            if (size == 4) EXPECT_NEAR(Mat4_expected[1][col], vec4_result[col], tolerance) << "Matrix: " << size << " Col: " << col;
            if (size == 5) EXPECT_NEAR(Mat5_expected[1][col], vec5_result[col], tolerance) << "Matrix: " << size << " Col: " << col;
            if (size == 6) EXPECT_NEAR(Mat6_expected[1][col], vec6_result[col], tolerance) << "Matrix: " << size << " Col: " << col;
        }
    }
}

TEST_F(Linear_Algebra_Misc_, Two_Through_Six_Square_Mat_Row_Return_Test) {

    mat2x2 Mat2_expected {{1,0},{0,1}};
    mat3x3 Mat3_expected {{1,0,0},{0,1,0},{0,0,1}};
    mat4x4 Mat4_expected {{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}};
    mat5x5 Mat5_expected {{1,0,0,0,0},{0,1,0,0,0},{0,0,1,0,0},{0,0,0,1,0},{0,0,0,0,1}};
    mat6x6 Mat6_expected {{1,0,0,0,0,0},{0,1,0,0,0,0},{0,0,1,0,0,0},{0,0,0,1,0,0},{0,0,0,0,1,0},{0,0,0,0,0,1}};

    mat2x2 Mat2_result; mat2x2_identity(Mat2_result); vec2 vec2_A;
    mat3x3 Mat3_result; mat3x3_identity(Mat3_result); vec3 vec3_A;
    mat4x4 Mat4_result; mat4x4_identity(Mat4_result); vec4 vec4_A;
    mat5x5 Mat5_result; mat5x5_identity(Mat5_result); vec5 vec5_A;
    mat6x6 Mat6_result; mat6x6_identity(Mat6_result); vec6 vec6_A;

    float* vec2_result = mat2x2_row_r(vec2_A, Mat2_result, 1);
    float* vec3_result = mat3x3_row_r(vec3_A, Mat3_result, 1);
    float* vec4_result = mat4x4_row_r(vec4_A, Mat4_result, 1);
    float* vec5_result = mat5x5_row_r(vec5_A, Mat5_result, 1);
    float* vec6_result = mat6x6_row_r(vec6_A, Mat6_result, 1);

    for (std::size_t size = 2; size<7; ++size) {
        for (std::size_t col = 0; col < size; ++col) {
            if (size == 2) EXPECT_NEAR(Mat2_expected[1][col], vec2_result[col], tolerance) << "Matrix: " << size << " Col: " << col;
            if (size == 3) EXPECT_NEAR(Mat3_expected[1][col], vec3_result[col], tolerance) << "Matrix: " << size << " Col: " << col;
            if (size == 4) EXPECT_NEAR(Mat4_expected[1][col], vec4_result[col], tolerance) << "Matrix: " << size << " Col: " << col;
            if (size == 5) EXPECT_NEAR(Mat5_expected[1][col], vec5_result[col], tolerance) << "Matrix: " << size << " Col: " << col;
            if (size == 6) EXPECT_NEAR(Mat6_expected[1][col], vec6_result[col], tolerance) << "Matrix: " << size << " Col: " << col;
        }
    }
}

TEST_F(Linear_Algebra_Misc_, Two_Through_Six_Square_Mat_Col_Test) {

    mat2x2 Mat2_expected {{1,0},{0,1}};
    mat3x3 Mat3_expected {{1,0,0},{0,1,0},{0,0,1}};
    mat4x4 Mat4_expected {{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}};
    mat5x5 Mat5_expected {{1,0,0,0,0},{0,1,0,0,0},{0,0,1,0,0},{0,0,0,1,0},{0,0,0,0,1}};
    mat6x6 Mat6_expected {{1,0,0,0,0,0},{0,1,0,0,0,0},{0,0,1,0,0,0},{0,0,0,1,0,0},{0,0,0,0,1,0},{0,0,0,0,0,1}};

    mat2x2 Mat2_result; mat2x2_identity(Mat2_result);
    mat3x3 Mat3_result; mat3x3_identity(Mat3_result);
    mat4x4 Mat4_result; mat4x4_identity(Mat4_result);
    mat5x5 Mat5_result; mat5x5_identity(Mat5_result);
    mat6x6 Mat6_result; mat6x6_identity(Mat6_result);

    vec2 vec2_result; mat2x2_col(vec2_result, Mat2_result, 1);
    vec3 vec3_result; mat3x3_col(vec3_result, Mat3_result, 1);
    vec4 vec4_result; mat4x4_col(vec4_result, Mat4_result, 1);
    vec5 vec5_result; mat5x5_col(vec5_result, Mat5_result, 1);
    vec6 vec6_result; mat6x6_col(vec6_result, Mat6_result, 1);

    for (std::size_t size = 2; size<7; ++size) {
        for (std::size_t row = 0; row < size; ++row) {
            if (size == 2) EXPECT_NEAR(Mat2_expected[row][1], vec2_result[row], tolerance) << "Matrix: " << size << " Row: " << row;
            if (size == 3) EXPECT_NEAR(Mat3_expected[row][1], vec3_result[row], tolerance) << "Matrix: " << size << " Row: " << row;
            if (size == 4) EXPECT_NEAR(Mat4_expected[row][1], vec4_result[row], tolerance) << "Matrix: " << size << " Row: " << row;
            if (size == 5) EXPECT_NEAR(Mat5_expected[row][1], vec5_result[row], tolerance) << "Matrix: " << size << " Row: " << row;
            if (size == 6) EXPECT_NEAR(Mat6_expected[row][1], vec6_result[row], tolerance) << "Matrix: " << size << " Row: " << row;
        }
    }
}

TEST_F(Linear_Algebra_Misc_, Two_Through_Six_Square_Mat_Col_Return_Test) {

    mat2x2 Mat2_expected {{1,0},{0,1}};
    mat3x3 Mat3_expected {{1,0,0},{0,1,0},{0,0,1}};
    mat4x4 Mat4_expected {{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}};
    mat5x5 Mat5_expected {{1,0,0,0,0},{0,1,0,0,0},{0,0,1,0,0},{0,0,0,1,0},{0,0,0,0,1}};
    mat6x6 Mat6_expected {{1,0,0,0,0,0},{0,1,0,0,0,0},{0,0,1,0,0,0},{0,0,0,1,0,0},{0,0,0,0,1,0},{0,0,0,0,0,1}};

    mat2x2 Mat2_result; mat2x2_identity(Mat2_result); vec2 vec2_A;
    mat3x3 Mat3_result; mat3x3_identity(Mat3_result); vec3 vec3_A;
    mat4x4 Mat4_result; mat4x4_identity(Mat4_result); vec4 vec4_A;
    mat5x5 Mat5_result; mat5x5_identity(Mat5_result); vec5 vec5_A;
    mat6x6 Mat6_result; mat6x6_identity(Mat6_result); vec6 vec6_A;

    float* vec2_result = mat2x2_col_r(vec2_A, Mat2_result, 1);
    float* vec3_result = mat3x3_col_r(vec3_A, Mat3_result, 1);
    float* vec4_result = mat4x4_col_r(vec4_A, Mat4_result, 1);
    float* vec5_result = mat5x5_col_r(vec5_A, Mat5_result, 1);
    float* vec6_result = mat6x6_col_r(vec6_A, Mat6_result, 1);

    for (std::size_t size = 2; size<7; ++size) {
        for (std::size_t row = 0; row < size; ++row) {
            if (size == 2) EXPECT_NEAR(Mat2_expected[row][1], vec2_result[row], tolerance) << "Matrix: " << size << " Row: " << row;
            if (size == 3) EXPECT_NEAR(Mat3_expected[row][1], vec3_result[row], tolerance) << "Matrix: " << size << " Row: " << row;
            if (size == 4) EXPECT_NEAR(Mat4_expected[row][1], vec4_result[row], tolerance) << "Matrix: " << size << " Row: " << row;
            if (size == 5) EXPECT_NEAR(Mat5_expected[row][1], vec5_result[row], tolerance) << "Matrix: " << size << " Row: " << row;
            if (size == 6) EXPECT_NEAR(Mat6_expected[row][1], vec6_result[row], tolerance) << "Matrix: " << size << " Row: " << row;
        }
    }
}

TEST_F(Linear_Algebra_Misc_, Two_Through_Six_Square_Mat_Transpose_Test) {

    mat2x2 Mat2_expected {{0,2},{1,3}};
    mat3x3 Mat3_expected {{0,3,6},{1,4,7},{2,5,8}};
    mat4x4 Mat4_expected {{0,4,8,12},{1,5,9,13},{2,6,10,14},{3,7,11,15}};
    mat5x5 Mat5_expected {{0,5,10,15,20},{1,6,11,16,21},{2,7,12,17,22},{3,8,13,18,23},{4,9,14,19,24}};
    mat6x6 Mat6_expected {{0,6,12,18,24,30},{1,7,13,19,25,31},{2,8,14,20,26,32},{3,9,15,21,27,33},{4,10,16,22,28,34},{5,11,17,23,29,35}};

    mat2x2 Mat2_result; mat2x2_copy(Mat2_result, Mat2A);
    mat3x3 Mat3_result; mat3x3_copy(Mat3_result, Mat3A);
    mat4x4 Mat4_result; mat4x4_copy(Mat4_result, Mat4A);
    mat5x5 Mat5_result; mat5x5_copy(Mat5_result, Mat5A);
    mat6x6 Mat6_result; mat6x6_copy(Mat6_result, Mat6A);

    mat2x2_transpose(Mat2_result);
    mat3x3_transpose(Mat3_result);
    mat4x4_transpose(Mat4_result);
    mat5x5_transpose(Mat5_result);
    mat6x6_transpose(Mat6_result);

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

TEST_F(Linear_Algebra_Misc_, Two_Through_Six_Square_Mat_Transpose_Return_Test) {

    mat2x2 Mat2_expected {{0,2},{1,3}};
    mat3x3 Mat3_expected {{0,3,6},{1,4,7},{2,5,8}};
    mat4x4 Mat4_expected {{0,4,8,12},{1,5,9,13},{2,6,10,14},{3,7,11,15}};
    mat5x5 Mat5_expected {{0,5,10,15,20},{1,6,11,16,21},{2,7,12,17,22},{3,8,13,18,23},{4,9,14,19,24}};
    mat6x6 Mat6_expected {{0,6,12,18,24,30},{1,7,13,19,25,31},{2,8,14,20,26,32},{3,9,15,21,27,33},{4,10,16,22,28,34},{5,11,17,23,29,35}};

    vec2 * Mat2_result = mat2x2_transpose_r(Mat2A);
    vec3 * Mat3_result = mat3x3_transpose_r(Mat3A);
    vec4 * Mat4_result = mat4x4_transpose_r(Mat4A);
    vec5 * Mat5_result = mat5x5_transpose_r(Mat5A);
    vec6 * Mat6_result = mat6x6_transpose_r(Mat6A);

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

TEST_F(Linear_Algebra_Misc_, Four_Square_Mat_Get_Rotational_Mat_Test) {
    mat3x3 mat3_expected {{0,1,2},{4,5,6},{8,9,10}};

    mat3x3 mat3_result; mat3x3_copy(mat3_result, Mat3_zero);

    mat4x4_get_rotational(mat3_result, Mat4A);

    for (std::size_t row = 0; row < 3; ++row) {
        for (std::size_t col = 0; col < 3; ++col) {
            EXPECT_NEAR(mat3_expected[row][col], mat3_result[row][col], tolerance) << " Row: " << row << " Col: " << col;
        }
    }
}

TEST_F(Linear_Algebra_Misc_, Four_Square_Mat_Return_Rotational_Mat_Test) {
    mat3x3 mat3_expected {{0,1,2},{4,5,6},{8,9,10}};

    vec3 * mat3_result = mat4x4_get_rotational_rn(Mat3_zero, Mat4A);

    for (std::size_t row = 0; row < 3; ++row) {
        for (std::size_t col = 0; col < 3; ++col) {
            EXPECT_NEAR(mat3_expected[row][col], mat3_result[row][col], tolerance) << " Row: " << row << " Col: " << col;
        }
    }
}

TEST_F(Linear_Algebra_Misc_, Four_Square_Mat_Translate_Mat_Test) {
    mat4x4 mat4_expected {{0,1,2,4},{4,5,6,8},{8,9,10,12},{12,13,14,15}};

    mat4x4_translate(Mat4A, 1.0f, 1.0f, 1.0f);

    for (std::size_t row = 0; row < 4; ++row) {
        for (std::size_t col = 0; col < 4; ++col) {
            EXPECT_NEAR(mat4_expected[row][col], Mat4A[row][col], tolerance) << " Row: " << row << " Col: " << col;
        }
    }
}


TEST_F(Linear_Algebra_Misc_, Four_Square_Mat_From_Vec_3_by_Vec3) {
    vec3 a = {1, 2, 3};
    vec3 b = {1, 2, 3};

    mat4x4 mat4_expected {{ 1, 2, 3, 0},
                          { 2, 4, 6, 0},
                          { 3, 6, 9, 0},
                          { 0, 0, 0, 0}};

    mat4x4_from_vec3_mul_outer(Mat4A, a, b);

    for (std::size_t row = 0; row < 4; ++row) {
        for (std::size_t col = 0; col < 4; ++col) {
            EXPECT_NEAR(mat4_expected[row][col], Mat4A[row][col], tolerance) << " Row: " << row << " Col: " << col;
        }
    }
}


TEST_F(Linear_Algebra_Misc_, Four_Square_Mat_Rotate_With_Angel_Axis) {
    mat4x4 identity = IDENTITY4x4;
    mat4x4 result   = IDENTITY4x4;
    mat4x4 expected = IDENTITY4x4;

    mat4x4_rotate(result, identity, 0, 0, 0, 0);

    for (std::size_t row = 0; row < 4; ++row) {
        for (std::size_t col = 0; col < 4; ++col) {
            EXPECT_NEAR(expected[row][col], result[row][col], tolerance) << " Row: " << row << " Col: " << col;
        }
    }

    mat4x4 expect_x_rotation = {{  1, 0, 0, 0},
                                {  0, 0,-1, 0},
                                {  0, 1, 0, 0},
                                {  0, 0, 0, 1}};

    mat4x4_rotate(result, identity, 1.0f, 0, 0, PI_2);

    for (std::size_t row = 0; row < 4; ++row) {
        for (std::size_t col = 0; col < 4; ++col) {
            EXPECT_NEAR(expect_x_rotation[row][col], result[row][col], tolerance) << " Row: " << row << " Col: " << col;
        }
    }

    mat4x4 expect_y_rotation = {{  0, 0, 1, 0},
                                {  0, 1, 0, 0},
                                { -1, 0, 0, 0},
                                {  0, 0, 0, 1}};

    mat4x4_rotate(result, identity, 0, 1, 0, PI_2);

    for (std::size_t row = 0; row < 4; ++row) {
        for (std::size_t col = 0; col < 4; ++col) {
            EXPECT_NEAR(expect_y_rotation[row][col], result[row][col], tolerance) << " Row: " << row << " Col: " << col;
        }
    }

    mat4x4 expect_z_rotation = {{  0,-1, 0, 0},
                                {  1, 0, 0, 0},
                                {  0, 0, 1, 0},
                                {  0, 0, 0, 1}};

    mat4x4_rotate(result, identity, 0, 0, 1, PI_2);

    for (std::size_t row = 0; row < 4; ++row) {
        for (std::size_t col = 0; col < 4; ++col) {
            EXPECT_NEAR(expect_z_rotation[row][col], result[row][col], tolerance) << " Row: " << row << " Col: " << col;
        }
    }
}


TEST_F(Linear_Algebra_Misc_, Four_Square_Mat_Rotate_With_Direct_Rotation) {
    mat4x4 identity = IDENTITY4x4;
    mat4x4 result   = IDENTITY4x4;

    mat4x4 expect_x_rotation = {{  1, 0, 0, 0},
                                {  0, 0,-1, 0},
                                {  0, 1, 0, 0},
                                {  0, 0, 0, 1}};

    mat4x4_rotate_about_x(result, identity, PI_2);

    for (std::size_t row = 0; row < 4; ++row) {
        for (std::size_t col = 0; col < 4; ++col) {
            EXPECT_NEAR(expect_x_rotation[row][col], result[row][col], tolerance) << " Row: " << row << " Col: " << col;
        }
    }

    mat4x4 expect_y_rotation = {{  0, 0, 1, 0},
                                {  0, 1, 0, 0},
                                { -1, 0, 0, 0},
                                {  0, 0, 0, 1}};

    mat4x4_rotate_about_y(result, identity, PI_2);

    for (std::size_t row = 0; row < 4; ++row) {
        for (std::size_t col = 0; col < 4; ++col) {
            EXPECT_NEAR(expect_y_rotation[row][col], result[row][col], tolerance) << " Row: " << row << " Col: " << col;
        }
    }

    mat4x4 expect_z_rotation = {{  0,-1, 0, 0},
                                {  1, 0, 0, 0},
                                {  0, 0, 1, 0},
                                {  0, 0, 0, 1}};

    mat4x4_rotate_about_z(result, identity, PI_2);

    for (std::size_t row = 0; row < 4; ++row) {
        for (std::size_t col = 0; col < 4; ++col) {
            EXPECT_NEAR(expect_z_rotation[row][col], result[row][col], tolerance) << " Row: " << row << " Col: " << col;
        }
    }
}


TEST_F(Linear_Algebra_Misc_, Four_Square_Mat_Inversion) {
    float local_tolerance = 0.001;

    mat4x4 result   = IDENTITY4x4;
    mat4x4 test     = {
        {0.4218f, 0.6557f, 0.6787f, 0.6555f},
        {0.9157f, 0.0357f, 0.7577f, 0.1712f},
        {0.7922f, 0.8491f, 0.7431f, 0.7060f},
        {0.9595f, 0.9340f, 0.3922f, 0.0318f},
    };

    mat4x4_invert(result, test);

    mat4x4 expected = {{ -4.2510f,  0.2757f,  3.9060f, -0.5772f},
                       {  2.0134f, -0.7508f, -1.7446f,  1.2719f},
                       {  5.9198f,  1.1647f, -5.8282f,  1.0980f},
                       { -3.8823f, -0.6323f,  5.2662f, -2.0378f}};

    for (std::size_t row = 0; row < 4; ++row) {
        for (std::size_t col = 0; col < 4; ++col) {
            EXPECT_NEAR(expected[row][col], result[row][col], local_tolerance) << " Row: " << row << " Col: " << col;
        }
    }
}

