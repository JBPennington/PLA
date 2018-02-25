
extern "C" {
#include "../pla.h"
}

#include <gtest/gtest.h>
#include <iostream>


class Linear_Algebra_Min_ : public testing::Test {
    void SetUp() {}
    void TearDown() {}
};

class Linear_Algebra_Max_ : public testing::Test {
    void SetUp() {}
    void TearDown() {}
};

TEST_F(Linear_Algebra_Min_, Min_Two_Vec2_Test) {
    vec2    a = { 0, 1};
    vec2    b = {-1, 2};

    vec2 result = {0, 0};
    vec2_min(result, a, b);

    vec2 expect = {-1, 1};

    uint32_t iter;
    for (iter=0; iter<2; ++iter) {
        EXPECT_EQ(expect[iter], result[iter]);
    }
}

TEST_F(Linear_Algebra_Min_, Min_Two_Vec3_Test) {
    vec3    a = {0, 1, 2};
    vec3    b = {1, 0, 1};

    vec3 result = {0, 0, 0};
    vec3_min(result, a, b);

    vec3 expect = {0, 0, 1};

    uint32_t iter;
    for (iter=0; iter<3; ++iter) {
        EXPECT_EQ(expect[iter], result[iter]);
    }
}

TEST_F(Linear_Algebra_Min_, Min_Two_Vec4_Test) {
    vec4    a = {0, 1, 2, 3};
    vec4    b = {3, 2, 1, 0};

    vec4 result = {0, 0, 0, 0};
    vec4_min(result, a, b);

    vec4 expect = {0, 1, 1, 0};

    uint32_t iter;
    for (iter=0; iter<4; ++iter) {
        EXPECT_EQ(expect[iter], result[iter]);
    }
}

TEST_F(Linear_Algebra_Min_, Min_Two_Vec5_Test) {
    vec5    a = {0, 1, 2, 3, 4};
    vec5    b = {4, 3, 2, 1, 0};

    vec5 result = {0, 0, 0, 0, 0};
    vec5_min(result, a, b);

    vec5 expect = {0, 1, 2, 1, 0};

    uint32_t iter;
    for (iter=0; iter<5; ++iter) {
        EXPECT_EQ(expect[iter], result[iter]);
    }
}

TEST_F(Linear_Algebra_Min_, Min_Two_Vec6_Test) {
    vec6    a = {0, 1, 2, 3, 4, 5};
    vec6    b = {5, 4, 3, 2, 1, 0};

    vec6 result = {0, 0, 0, 0, 0, 0};
    vec6_min(result, a, b);

    vec6 expect = {0, 1, 2, 2, 1, 0};

    uint32_t iter;
    for (iter=0; iter<6; ++iter) {
        EXPECT_EQ(expect[iter], result[iter]);
    }
}

TEST_F(Linear_Algebra_Max_, Max_Two_Vec2_Test) {
    vec2    a = { 0, 1};
    vec2    b = {-1, 2};

    vec2 result = {0, 0};
    vec2_max(result, a, b);

    vec2 expect = {0,2};

    uint32_t iter;
    for (iter=0; iter<2; ++iter) {
        EXPECT_EQ(expect[iter], result[iter]);
    }
}

TEST_F(Linear_Algebra_Max_, Max_Two_Vec3_Test) {
    vec3    a = {0, 1, 2};
    vec3    b = {1, 0, 1};

    vec3 result = {0, 0, 0};
    vec3_max(result, a, b);

    vec3 expect = {1,1,2};

    uint32_t iter;
    for (iter=0; iter<3; ++iter) {
        EXPECT_EQ(expect[iter], result[iter]);
    }
}

TEST_F(Linear_Algebra_Max_, Max_Two_Vec4_Test) {
    vec4    a = {0, 1, 2, 3};
    vec4    b = {3, 2, 1, 0};

    vec4 result = {0, 0, 0, 0};
    vec4_max(result, a, b);

    vec4 expect = {3,2,2,3};

    uint32_t iter;
    for (iter=0; iter<4; ++iter) {
        EXPECT_EQ(expect[iter], result[iter]);
    }
}

TEST_F(Linear_Algebra_Max_, Max_Two_Vec5_Test) {
    vec5    a = {0, 1, 2, 3, 4};
    vec5    b = {4, 3, 2, 1, 0};

    vec5 result = {0, 0, 0, 0, 0};
    vec5_max(result, a, b);

    vec5 expect = {4,3,2,3,4};

    uint32_t iter;
    for (iter=0; iter<5; ++iter) {
        EXPECT_EQ(expect[iter], result[iter]);
    }
}

TEST_F(Linear_Algebra_Max_, Max_Two_Vec6_Test) {
    vec6    a = {0, 1, 2, 3, 4, 5};
    vec6    b = {5, 4, 3, 2, 1, 0};

    vec6 result = {0, 0, 0, 0, 0, 0};
    vec6_max(result, a, b);

    vec6 expect = {5,4,3,3,4,5};

    uint32_t iter;
    for (iter=0; iter<6; ++iter) {
        EXPECT_EQ(expect[iter], result[iter]);
    }
}

