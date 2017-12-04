
extern "C" {
#include "../linear_algebra.h"
}

#include <gtest/gtest.h>
#include <iostream>


class Linear_Algebra_Mat_Mul_Vec_ : public testing::Test {
    void SetUp() {

    }
    void TearDown() {

    }
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
    vec2 Vec2_zero   {0,0};
    vec3 Vec3_zero   {0,0,0};
    vec4 Vec4_zero   {0,0,0,0};
    vec5 Vec5_zero   {0,0,0,0,0};
    vec6 Vec6_zero   {0,0,0,0,0,0};
    vec2 Vec2A       {0,1};
    vec3 Vec3A       {0,1,2};
    vec4 Vec4A       {0,1,2,3};
    vec5 Vec5A       {0,1,2,3,4};
    vec6 Vec6A       {0,1,2,3,4,5};
    vec2 Vec2B       {0,1};
    vec3 Vec3B       {0,1,2};
    vec4 Vec4B       {0,1,2,3};
    vec5 Vec5B       {0,1,2,3,4};
    vec6 Vec6B       {0,1,2,3,4,5};
};

TEST_F(Linear_Algebra_Mat_Mul_Vec_, Two_Through_Six_Square_Mat_New_Vec_Mat_Cross_Vec_Test) {
    vec2 vec2_expected {1,3};
    vec3 vec3_expected {5,14,23};
    vec4 vec4_expected {14,38,62,86};
    vec5 vec5_expected {30,80,130,180,230};
    vec6 vec6_expected {55,145,235,325,415,505};

    vec2 vec2_result; vec2_copy(vec2_result, Vec2_zero);
    vec3 vec3_result; vec3_copy(vec3_result, Vec3_zero);
    vec4 vec4_result; vec4_copy(vec4_result, Vec4_zero);
    vec5 vec5_result; vec5_copy(vec5_result, Vec5_zero);
    vec6 vec6_result; vec6_copy(vec6_result, Vec6_zero);

    for (std::size_t size = 2; size<7; ++size) {
        if (size == 2) mat2x2_mul_vec2_n(vec2_result, Mat2A, Vec2A);
        if (size == 3) mat3x3_mul_vec3_n(vec3_result, Mat3A, Vec3A);
        if (size == 4) mat4x4_mul_vec4_n(vec4_result, Mat4A, Vec4A);
        if (size == 5) mat5x5_mul_vec5_n(vec5_result, Mat5A, Vec5A);
        if (size == 6) mat6x6_mul_vec6_n(vec6_result, Mat6A, Vec6A);
        for (std::size_t row = 0; row < size; ++row) {
            if (size == 2) EXPECT_NEAR(vec2_expected[row], vec2_result[row], tolerance) << "Matrix: " << size << " Row: " << row;
            if (size == 3) EXPECT_NEAR(vec3_expected[row], vec3_result[row], tolerance) << "Matrix: " << size << " Row: " << row;
            if (size == 4) EXPECT_NEAR(vec4_expected[row], vec4_result[row], tolerance) << "Matrix: " << size << " Row: " << row;
            if (size == 5) EXPECT_NEAR(vec5_expected[row], vec5_result[row], tolerance) << "Matrix: " << size << " Row: " << row;
            if (size == 6) EXPECT_NEAR(vec6_expected[row], vec6_result[row], tolerance) << "Matrix: " << size << " Row: " << row;
        }
    }
}

TEST_F(Linear_Algebra_Mat_Mul_Vec_, Two_Through_Six_Square_Mat_New_Vec_Mat_Cross_Vec_Return_Test) {
    vec2 vec2_expected {1,3};
    vec3 vec3_expected {5,14,23};
    vec4 vec4_expected {14,38,62,86};
    vec5 vec5_expected {30,80,130,180,230};
    vec6 vec6_expected {55,145,235,325,415,505};

    float * vec2_result = mat2x2_mul_vec2_r(Mat2A, Vec2A);
    float * vec3_result = mat3x3_mul_vec3_r(Mat3A, Vec3A);
    float * vec4_result = mat4x4_mul_vec4_r(Mat4A, Vec4A);
    float * vec5_result = mat5x5_mul_vec5_r(Mat5A, Vec5A);
    float * vec6_result = mat6x6_mul_vec6_r(Mat6A, Vec6A);

    for (std::size_t size = 2; size<7; ++size) {
        for (std::size_t row = 0; row < size; ++row) {
            if (size == 2) EXPECT_NEAR(vec2_expected[row], vec2_result[row], tolerance) << "Matrix: " << size << " Row: " << row;
            if (size == 3) EXPECT_NEAR(vec3_expected[row], vec3_result[row], tolerance) << "Matrix: " << size << " Row: " << row;
            if (size == 4) EXPECT_NEAR(vec4_expected[row], vec4_result[row], tolerance) << "Matrix: " << size << " Row: " << row;
            if (size == 5) EXPECT_NEAR(vec5_expected[row], vec5_result[row], tolerance) << "Matrix: " << size << " Row: " << row;
            if (size == 6) EXPECT_NEAR(vec6_expected[row], vec6_result[row], tolerance) << "Matrix: " << size << " Row: " << row;
        }
    }
}

TEST_F(Linear_Algebra_Mat_Mul_Vec_, Two_Through_Six_Square_Mat_New_Vec_Mat_Cross_Vec_Return_New_Test) {
    vec2 vec2_expected {1,3};
    vec3 vec3_expected {5,14,23};
    vec4 vec4_expected {14,38,62,86};
    vec5 vec5_expected {30,80,130,180,230};
    vec6 vec6_expected {55,145,235,325,415,505};

    vec2 vec2_temp; vec2_copy(vec2_temp, Vec2_zero);
    vec3 vec3_temp; vec3_copy(vec3_temp, Vec3_zero);
    vec4 vec4_temp; vec4_copy(vec4_temp, Vec4_zero);
    vec5 vec5_temp; vec5_copy(vec5_temp, Vec5_zero);
    vec6 vec6_temp; vec6_copy(vec6_temp, Vec6_zero);

    float * vec2_result = mat2x2_mul_vec2_rn(vec2_temp, Mat2A, Vec2A);
    float * vec3_result = mat3x3_mul_vec3_rn(vec3_temp, Mat3A, Vec3A);
    float * vec4_result = mat4x4_mul_vec4_rn(vec4_temp, Mat4A, Vec4A);
    float * vec5_result = mat5x5_mul_vec5_rn(vec5_temp, Mat5A, Vec5A);
    float * vec6_result = mat6x6_mul_vec6_rn(vec6_temp, Mat6A, Vec6A);

    for (std::size_t size = 2; size<7; ++size) {
        for (std::size_t row = 0; row < size; ++row) {
            if (size == 2) EXPECT_NEAR(vec2_expected[row], vec2_result[row], tolerance) << "Matrix: " << size << " Row: " << row;
            if (size == 3) EXPECT_NEAR(vec3_expected[row], vec3_result[row], tolerance) << "Matrix: " << size << " Row: " << row;
            if (size == 4) EXPECT_NEAR(vec4_expected[row], vec4_result[row], tolerance) << "Matrix: " << size << " Row: " << row;
            if (size == 5) EXPECT_NEAR(vec5_expected[row], vec5_result[row], tolerance) << "Matrix: " << size << " Row: " << row;
            if (size == 6) EXPECT_NEAR(vec6_expected[row], vec6_result[row], tolerance) << "Matrix: " << size << " Row: " << row;
        }
    }
}
