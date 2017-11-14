
extern "C" {
#include "../linear_algebra.h"
}

#include <gtest/gtest.h>
#include <iostream>


class Linear_Algebra_Add_ : public testing::Test {
    void SetUp() {}
    void TearDown() {}
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
