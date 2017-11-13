
extern "C" {
#include "../linear_algebra.h"
}

#include <gtest/gtest.h>
#include <iostream>


class Linear_Algebra_Mag_ : public testing::Test {
    void SetUp() {}
    void TearDown() {}
};

TEST_F(Linear_Algebra_Mag_, Mag_Vec2_Test) {
    vec2 a = {0, 1};

    float result = vec2_mag(a);

    float expect = 1;
    EXPECT_NEAR(expect, result, 0.001);
}

TEST_F(Linear_Algebra_Mag_, Mag_Vec3_Test) {
    vec3 a = {0, 1, 2};

    float result = vec3_mag(a);

    float expect = 2.23607;
    EXPECT_NEAR(expect, result, 0.001);
}

TEST_F(Linear_Algebra_Mag_, Mag_Vec4_Test) {
    vec4 a = {0, 1, 2, 3};

    float result = vec4_mag(a);

    float expect = 3.74166;
    EXPECT_NEAR(expect, result, 0.001);
}

TEST_F(Linear_Algebra_Mag_, Mag_Vec5_Test) {
    vec5    a = {0, 1, 2, 3, 4};

    float result = vec5_mag(a);

    float expect = 5.47723;
    EXPECT_NEAR(expect, result, 0.001);
}

TEST_F(Linear_Algebra_Mag_, Mag_Vec6_Test) {
    vec6    a = {0, 1, 2, 3, 4, 5};

    float result = vec6_mag(a);

    float expect = 7.4162;
    EXPECT_NEAR(expect, result, 0.001);
}
