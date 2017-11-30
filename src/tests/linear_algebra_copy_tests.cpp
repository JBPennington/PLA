
extern "C" {
#include "../linear_algebra.h"
}

#include <gtest/gtest.h>
#include <iostream>


class Linear_Algebra_Copy_ : public testing::Test {
    void SetUp() {}
    void TearDown() {}
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

