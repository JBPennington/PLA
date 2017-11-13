
extern "C" {
#include "../linear_algebra.h"
}

#include <gtest/gtest.h>
#include <iostream>


class Linear_Algebra_Mul_Inner : public testing::Test {
    void SetUp() {}
    void TearDown() {}
};

TEST_F(Linear_Algebra_Mul_Inner, Mul_Inner_Two_Vec2_Test) {
    vec2 a = {0, 1};
    vec2 b = {2, 3};

    float result = vec2_mul_inner(a, b);

    float expect = 3;
    EXPECT_EQ(expect, result);
}

TEST_F(Linear_Algebra_Mul_Inner, Mul_Inner_Two_Vec3_Test) {
    vec3 a = {0, 1, 2};
    vec3 b = {3, 4, 5};

    float result = vec3_mul_inner(a, b);

    float expect = 14;
    EXPECT_EQ(expect, result);
}

TEST_F(Linear_Algebra_Mul_Inner, Mul_Inner_Two_Vec4_Test) {
    vec4 a = {0, 1, 2, 3};
    vec4 b = {4, 5, 6, 7};

    float result = vec4_mul_inner(a, b);

    float expect = 38;
    EXPECT_EQ(expect, result);
}

TEST_F(Linear_Algebra_Mul_Inner, Mul_Inner_Two_Vec5_Test) {
    vec5 a = {0, 1, 2, 3, 4};
    vec5 b = {5, 6, 7, 8, 9};

    float result = vec5_mul_inner(a, b);

    float expect = 80;
    EXPECT_EQ(expect, result);
}

TEST_F(Linear_Algebra_Mul_Inner, Mul_Inner_Two_Vec6_Test) {
    vec6 a = {0, 1, 2, 3,  4,  5};
    vec6 b = {6, 7, 8, 9, 10, 11};

    float result = vec6_mul_inner(a, b);

    float expect = 145;
    EXPECT_EQ(expect, result);
}
