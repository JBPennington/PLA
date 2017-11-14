
extern "C" {
#include "../linear_algebra.h"
}

#include <gtest/gtest.h>
#include <iostream>


class Linear_Algebra_Scale_ : public testing::Test {
    void SetUp() {}
    void TearDown() {}
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

