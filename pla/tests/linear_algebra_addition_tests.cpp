
extern "C" {
#include "../pla.h"
}

#include <gtest/gtest.h>
#include <iostream>


class Linear_Algebra_Add_ : public testing::Test {
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

TEST_F(Linear_Algebra_Add_, Add_Two_Vec2_Test) {
    vec2 a = {0, 1};
    vec2 b = {2, 3};

    vec2 result = {0, 0};
    vec2_add_n(result, a, b);

    vec2 expect = {2, 4};

    uint32_t iter;
    for (iter=0; iter<2; ++iter) {
        EXPECT_EQ(expect[iter], result[iter]);
    }
}

TEST_F(Linear_Algebra_Add_, Add_Two_Vec3_Test) {
    vec3 a = {0, 1, 2};
    vec3 b = {3, 4, 5};

    vec3 result = {0, 0, 0};
    vec3_add_n(result, a, b);

    vec3 expect = {3, 5, 7};

    uint32_t iter;
    for (iter=0; iter<3; ++iter) {
        EXPECT_EQ(expect[iter], result[iter]);
    }
}

TEST_F(Linear_Algebra_Add_, Add_Two_Vec4_Test) {
    vec4 a = {0, 1, 2, 3};
    vec4 b = {4, 5, 6, 7};

    vec4 result = {0, 0, 0, 0};
    vec4_add_n(result, a, b);

    vec4 expect = {4, 6, 8, 10};

    uint32_t iter;
    for (iter=0; iter<4; ++iter) {
        EXPECT_EQ(expect[iter], result[iter]);
    }
}

TEST_F(Linear_Algebra_Add_, Add_Two_Vec5_Test) {
    vec5 a = {0, 1, 2, 3, 4};
    vec5 b = {5, 6, 7, 8, 9};

    vec5 result = {0, 0, 0, 0, 0};
    vec5_add_n(result, a, b);

    vec5 expect = {5, 7, 9, 11, 13};

    uint32_t iter;
    for (iter=0; iter<5; ++iter) {
        EXPECT_EQ(expect[iter], result[iter]);
    }
}

TEST_F(Linear_Algebra_Add_, Add_Two_Vec6_Test) {
    vec6 a = {0, 1, 2, 3, 4, 5};
    vec6 b = {6, 7, 8, 9, 10, 11};

    vec6 result = {0, 0, 0, 0, 0, 0};
    vec6_add_n(result, a, b);

    vec6 expect = {6, 8, 10, 12, 14, 16};

    uint32_t iter;
    for (iter=0; iter<6; ++iter) {
        EXPECT_EQ(expect[iter], result[iter]);
    }
}

TEST_F(Linear_Algebra_Add_, Add_Two_Vec2_In_Place_Test) {
    vec2 a = {0, 1};
    vec2 b = {2, 3};

    vec2_add(a, b);

    vec2 expect = {2, 4};
    uint32_t iter;
    for (iter=0; iter<2; ++iter) {
        EXPECT_EQ(expect[iter], a[iter]);
    }
}

TEST_F(Linear_Algebra_Add_, Add_Two_Vec3_In_Place_Test) {
    vec3 a = {0, 1, 2};
    vec3 b = {3, 4, 5};

    vec3_add(a, b);

    vec3 expect = {3, 5, 7};
    uint32_t iter;
    for (iter=0; iter<3; ++iter) {
        EXPECT_EQ(expect[iter], a[iter]);
    }
}

TEST_F(Linear_Algebra_Add_, Add_Two_Vec4_In_Place_Test) {
    vec4 a = {0, 1, 2, 3};
    vec4 b = {4, 5, 6, 7};

    vec4_add(a, b);

    vec4 expect = {4, 6, 8, 10};
    uint32_t iter;
    for (iter=0; iter<4; ++iter) {
        EXPECT_EQ(expect[iter], a[iter]);
    }
}

TEST_F(Linear_Algebra_Add_, Add_Two_Vec5_In_Place_Test) {
    vec5 a = {0, 1, 2, 3, 4};
    vec5 b = {5, 6, 7, 8, 9};

    vec5_add(a, b);

    vec5 expect = {5, 7, 9, 11, 13};
    uint32_t iter;
    for (iter=0; iter<5; ++iter) {
        EXPECT_EQ(expect[iter], a[iter]);
    }
}

TEST_F(Linear_Algebra_Add_, Add_Two_Vec6_In_Place_Test) {
    vec6 a = {0, 1, 2, 3, 4, 5};
    vec6 b = {6, 7, 8, 9, 10, 11};

    vec6_add(a, b);

    vec6 expect = {6, 8, 10, 12, 14, 16};
    uint32_t iter;
    for (iter=0; iter<6; ++iter) {
        EXPECT_EQ(expect[iter], a[iter]);
    }
}

TEST_F(Linear_Algebra_Add_, Add_Return_Two_Vec2_Test) {
    vec2 a = {0, 1};
    vec2 b = {2, 3};

    float * result = vec2_add_r(a, b);

    vec2 expect = {2, 4};

    uint32_t iter;
    for (iter=0; iter<2; ++iter) {
        EXPECT_EQ(expect[iter], result[iter]);
    }
}

TEST_F(Linear_Algebra_Add_, Add_Return_Two_Vec3_Test) {
    vec3 a = {0, 1, 2};
    vec3 b = {3, 4, 5};

    float * result = vec3_add_r(a, b);

    vec3 expect = {3, 5, 7};

    uint32_t iter;
    for (iter=0; iter<3; ++iter) {
        EXPECT_EQ(expect[iter], result[iter]);
    }
}

TEST_F(Linear_Algebra_Add_, Add_Return_Two_Vec4_Test) {
    vec4 a = {0, 1, 2, 3};
    vec4 b = {4, 5, 6, 7};

    float * result = vec4_add_r(a, b);

    vec4 expect = {4, 6, 8, 10};

    uint32_t iter;
    for (iter=0; iter<4; ++iter) {
        EXPECT_EQ(expect[iter], result[iter]);
    }
}

TEST_F(Linear_Algebra_Add_, Add_Return_Two_Vec5_Test) {
    vec5 a = {0, 1, 2, 3, 4};
    vec5 b = {5, 6, 7, 8, 9};

    float * result = vec5_add_r(a, b);

    vec5 expect = {5, 7, 9, 11, 13};

    uint32_t iter;
    for (iter=0; iter<5; ++iter) {
        EXPECT_EQ(expect[iter], result[iter]);
    }
}

TEST_F(Linear_Algebra_Add_, Add_Return_Two_Vec6_Test) {
    vec6 a = {0, 1, 2, 3, 4, 5};
    vec6 b = {6, 7, 8, 9, 10, 11};

    float * result = vec6_add_r(a, b);

    vec6 expect = {6, 8, 10, 12, 14, 16};

    uint32_t iter;
    for (iter=0; iter<6; ++iter) {
        EXPECT_EQ(expect[iter], result[iter]);
    }
}

TEST_F(Linear_Algebra_Add_, Two_Through_Six_Square_Mat_Add_to_New_Mat_Test) {

    mat2x2 Mat2_expected {{0,2},{4,6}};
    mat3x3 Mat3_expected {{0,2,4},{6,8,10},{12,14,16}};
    mat4x4 Mat4_expected {{0,2,4,6},{8,10,12,14},{16,18,20,22},{24,26,28,30}};
    mat5x5 Mat5_expected {{0,2,4,6,8},{10,12,14,16,18},{20,22,24,26,28},{30,32,34,36,38},{40,42,44,46,48}};
    mat6x6 Mat6_expected {{0,2,4,6,8,10},{12,14,16,18,20,22},{24,26,28,30,32,34},{36,38,40,42,44,46},{48,50,52,54,56,58},{60,62,64,66,68,70}};

    mat2x2 Mat2_result; mat2x2_copy(Mat2_result, Mat2_zero);
    mat3x3 Mat3_result; mat3x3_copy(Mat3_result, Mat3_zero);
    mat4x4 Mat4_result; mat4x4_copy(Mat4_result, Mat4_zero);
    mat5x5 Mat5_result; mat5x5_copy(Mat5_result, Mat5_zero);
    mat6x6 Mat6_result; mat6x6_copy(Mat6_result, Mat6_zero);

    for (std::size_t size = 2; size<7; ++size) {
        if (size == 2) mat2x2_add_n(Mat2_result, Mat2A, Mat2B);;
        if (size == 3) mat3x3_add_n(Mat3_result, Mat3A, Mat3B);;
        if (size == 4) mat4x4_add_n(Mat4_result, Mat4A, Mat4B);;
        if (size == 5) mat5x5_add_n(Mat5_result, Mat5A, Mat5B);;
        if (size == 6) mat6x6_add_n(Mat6_result, Mat6A, Mat6B);;
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

TEST_F(Linear_Algebra_Add_, Two_Through_Six_Square_Mat_Add_In_Place_Test) {

    mat2x2 Mat2_expected {{0,2},{4,6}};
    mat3x3 Mat3_expected {{0,2,4},{6,8,10},{12,14,16}};
    mat4x4 Mat4_expected {{0,2,4,6},{8,10,12,14},{16,18,20,22},{24,26,28,30}};
    mat5x5 Mat5_expected {{0,2,4,6,8},{10,12,14,16,18},{20,22,24,26,28},{30,32,34,36,38},{40,42,44,46,48}};
    mat6x6 Mat6_expected {{0,2,4,6,8,10},{12,14,16,18,20,22},{24,26,28,30,32,34},{36,38,40,42,44,46},{48,50,52,54,56,58},{60,62,64,66,68,70}};

    mat2x2 Mat2_result; mat2x2_copy(Mat2_result, Mat2A);
    mat3x3 Mat3_result; mat3x3_copy(Mat3_result, Mat3A);
    mat4x4 Mat4_result; mat4x4_copy(Mat4_result, Mat4A);
    mat5x5 Mat5_result; mat5x5_copy(Mat5_result, Mat5A);
    mat6x6 Mat6_result; mat6x6_copy(Mat6_result, Mat6A);

    for (std::size_t size = 2; size<7; ++size) {
        if (size == 2) mat2x2_add(Mat2_result, Mat2B);
        if (size == 3) mat3x3_add(Mat3_result, Mat3B);
        if (size == 4) mat4x4_add(Mat4_result, Mat4B);
        if (size == 5) mat5x5_add(Mat5_result, Mat5B);
        if (size == 6) mat6x6_add(Mat6_result, Mat6B);
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

TEST_F(Linear_Algebra_Add_, Two_Through_Six_Square_Mat_Add_In_Place_And_Return_Test) {

    mat2x2 Mat2_expected {{0,2},{4,6}};
    mat3x3 Mat3_expected {{0,2,4},{6,8,10},{12,14,16}};
    mat4x4 Mat4_expected {{0,2,4,6},{8,10,12,14},{16,18,20,22},{24,26,28,30}};
    mat5x5 Mat5_expected {{0,2,4,6,8},{10,12,14,16,18},{20,22,24,26,28},{30,32,34,36,38},{40,42,44,46,48}};
    mat6x6 Mat6_expected {{0,2,4,6,8,10},{12,14,16,18,20,22},{24,26,28,30,32,34},{36,38,40,42,44,46},{48,50,52,54,56,58},{60,62,64,66,68,70}};

    vec2 * Mat2_result = mat2x2_add_r(Mat2A, Mat2B);
    vec3 * Mat3_result = mat3x3_add_r(Mat3A, Mat3B);
    vec4 * Mat4_result = mat4x4_add_r(Mat4A, Mat4B);
    vec5 * Mat5_result = mat5x5_add_r(Mat5A, Mat5B);
    vec6 * Mat6_result = mat6x6_add_r(Mat6A, Mat6B);

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
