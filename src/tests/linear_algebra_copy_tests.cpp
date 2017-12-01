
extern "C" {
#include "../linear_algebra.h"
}

#include <gtest/gtest.h>
#include <iostream>


class Linear_Algebra_Copy_ : public testing::Test {
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

TEST_F(Linear_Algebra_Copy_, Copy_Two_Vec2_Test) {
    vec2 a = {0, 1};
    vec2 b = {2, 3};

    vec2_copy(a, b);

    vec2 expect = {2, 3};

    uint32_t iter;
    for (iter=0; iter<2; ++iter) {
        EXPECT_EQ(expect[iter], a[iter]);
    }
}

TEST_F(Linear_Algebra_Copy_, Copy_Two_Vec3_Test) {
    vec3 a = {0, 1, 2};
    vec3 b = {3, 4, 5};

    vec3_copy(a, b);

    vec3 expect = {3, 4, 5};

    uint32_t iter;
    for (iter=0; iter<3; ++iter) {
        EXPECT_EQ(expect[iter], a[iter]);
    }
}

TEST_F(Linear_Algebra_Copy_, Copy_Two_Vec4_Test) {
    vec4 a = {0, 1, 2, 3};
    vec4 b = {4, 5, 6, 7};

    vec4_copy(a, b);

    vec4 expect = {4, 5, 6, 7};

    uint32_t iter;
    for (iter=0; iter<4; ++iter) {
        EXPECT_EQ(expect[iter], a[iter]);
    }
}

TEST_F(Linear_Algebra_Copy_, Copy_Two_Vec5_Test) {
    vec5 a = {0, 1, 2, 3, 4};
    vec5 b = {5, 6, 7, 8, 9};

    vec5_copy(a, b);

    vec5 expect = {5, 6, 7, 8, 9};

    uint32_t iter;
    for (iter=0; iter<5; ++iter) {
        EXPECT_EQ(expect[iter], a[iter]);
    }
}

TEST_F(Linear_Algebra_Copy_, Copy_Two_Vec6_Test) {
    vec6 a = {0, 1, 2, 3, 4, 5};
    vec6 b = {6, 7, 8, 9, 10, 11};

    vec6_copy(a, b);

    vec6 expect = {6, 7, 8, 9, 10, 11};

    uint32_t iter;
    for (iter=0; iter<6; ++iter) {
        EXPECT_EQ(expect[iter], a[iter]);
    }
}

TEST_F(Linear_Algebra_Copy_, Two_Through_Six_Square_Mat_Copy_Test) {

    // Creation Pushed to Setup now, not sure how I feel about that

    for (std::size_t size = 2; size<7; ++size) {
        switch(size) {
            case 2: mat2x2_copy(Mat2_zero, Mat2A);
            case 3: mat3x3_copy(Mat3_zero, Mat3A);
            case 4: mat4x4_copy(Mat4_zero, Mat4A);
            case 5: mat5x5_copy(Mat5_zero, Mat5A);
            case 6: mat6x6_copy(Mat6_zero, Mat6A);
            default:;
        }
        for (std::size_t row = 0; row < size; ++row) {
            for (std::size_t col = 0; col < size; ++col) {
                switch(size) {
                    case 2: EXPECT_NEAR(Mat2A[row][col], Mat2_zero[row][col], tolerance);
                    case 3: EXPECT_NEAR(Mat3A[row][col], Mat3_zero[row][col], tolerance);
                    case 4: EXPECT_NEAR(Mat4A[row][col], Mat4_zero[row][col], tolerance);
                    case 5: EXPECT_NEAR(Mat5A[row][col], Mat5_zero[row][col], tolerance);
                    case 6: EXPECT_NEAR(Mat6A[row][col], Mat6_zero[row][col], tolerance);
                    default:;
                }
            }
        }
    }
}

