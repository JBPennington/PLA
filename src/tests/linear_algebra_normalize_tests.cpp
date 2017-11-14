
extern "C" {
#include "../linear_algebra.h"
}

#include <gtest/gtest.h>
#include <iostream>


class Linear_Algebra_Norm_ : public testing::Test {
    void SetUp() {}
    void TearDown() {}
};

TEST_F(Linear_Algebra_Norm_, Normalize_Vec2_Test) {
    vec2 a = {0, 1};

    vec2 result = {0, 0};
    vec2_normalize_n(result, a);

    vec2 expect = {0, 1};
    uint32_t iter;
    for (iter=0; iter<2; ++iter) {
        EXPECT_NEAR(expect[iter], result[iter], 0.0001);
    }
}

TEST_F(Linear_Algebra_Norm_, Normalize_Vec3_Test) {
    vec3 a = {0, 1, 2};

    vec3 result = {0, 0, 0};
    vec3_normalize_n(result, a);

    vec3 expect = {0, 0.4472, 0.8944};
    uint32_t iter;
    for (iter=0; iter<3; ++iter) {
        EXPECT_NEAR(expect[iter], result[iter], 0.0001);
    }
}

TEST_F(Linear_Algebra_Norm_, Normalize_Vec4_Test) {
    vec4 a = {0, 1, 2, 3};

    vec4 result = {0, 0, 0, 0};
    vec4_normalize_n(result, a);

    vec4 expect = {0, 0.2673, 0.5345, 0.8018};
    uint32_t iter;
    for (iter=0; iter<4; ++iter) {
        EXPECT_NEAR(expect[iter], result[iter], 0.0001);
    }
}

TEST_F(Linear_Algebra_Norm_, Normalize_Vec5_Test) {
    vec5 a = {0, 1, 2, 3, 4};

    vec5 result = {0, 0, 0, 0, 0};
    vec5_normalize_n(result, a);

    vec5 expect = {0, 0.1826, 0.3651, 0.5477, 0.7303};
    uint32_t iter;
    for (iter=0; iter<5; ++iter) {
        EXPECT_NEAR(expect[iter], result[iter], 0.0001);
    }
}

TEST_F(Linear_Algebra_Norm_, Normalize_Vec6_Test) {
    vec6 a = {0, 1, 2, 3, 4, 5};

    vec6 result = {0, 0, 0, 0, 0, 0};
    vec6_normalize_n(result, a);

    vec6 expect = {0, 0.1348, 0.2697, 0.4045, 0.5394, 0.6742};
    uint32_t iter;
    for (iter=0; iter<6; ++iter) {
        EXPECT_NEAR(expect[iter], result[iter], 0.0001);
    }
}

TEST_F(Linear_Algebra_Norm_, Normalize_In_Place_Vec2_Test) {
    vec2 a = {0, 1};

    vec2_normalize(a);

    vec2 expect = {0, 1};
    uint32_t iter;
    for (iter=0; iter<2; ++iter) {
        EXPECT_NEAR(expect[iter], a[iter], 0.0001);
    }
}

TEST_F(Linear_Algebra_Norm_, Normalize_In_Place_Vec3_Test) {
    vec3 a = {0, 1, 2};

    vec3_normalize(a);

    vec3 expect = {0, 0.4472, 0.8944};
    uint32_t iter;
    for (iter=0; iter<3; ++iter) {
        EXPECT_NEAR(expect[iter], a[iter], 0.0001);
    }
}

TEST_F(Linear_Algebra_Norm_, Normalize_In_Place_Vec4_Test) {
    vec4 a = {0, 1, 2, 3};

    vec4_normalize(a);

    vec4 expect = {0, 0.2673, 0.5345, 0.8018};
    uint32_t iter;
    for (iter=0; iter<4; ++iter) {
        EXPECT_NEAR(expect[iter], a[iter], 0.0001);
    }
}

TEST_F(Linear_Algebra_Norm_, Normalize_In_Place_Vec5_Test) {
    vec5 a = {0, 1, 2, 3, 4};

    vec5_normalize(a);

    vec5 expect = {0, 0.1826, 0.3651, 0.5477, 0.7303};
    uint32_t iter;
    for (iter=0; iter<5; ++iter) {
        EXPECT_NEAR(expect[iter], a[iter], 0.0001);
    }
}

TEST_F(Linear_Algebra_Norm_, Normalize_In_Place_Vec6_Test) {
    vec6 a = {0, 1, 2, 3, 4, 5};

    vec6_normalize(a);

    vec6 expect = {0, 0.1348, 0.2697, 0.4045, 0.5394, 0.6742};
    uint32_t iter;
    for (iter=0; iter<6; ++iter) {
        EXPECT_NEAR(expect[iter], a[iter], 0.0001);
    }
}

TEST_F(Linear_Algebra_Norm_, Normalize_Return_Vec2_Test) {
    vec2 a = {0, 1};

    float * result = vec2_normalize_r(a);

    vec2 expect = {0, 1};
    uint32_t iter;
    for (iter=0; iter<2; ++iter) {
        EXPECT_NEAR(expect[iter], result[iter], 0.0001);
    }
}

TEST_F(Linear_Algebra_Norm_, Normalize_Return_Vec3_Test) {
    vec3 a = {0, 1, 2};

    float * result = vec3_normalize_r(a);

    vec3 expect = {0, 0.4472, 0.8944};
    uint32_t iter;
    for (iter=0; iter<3; ++iter) {
        EXPECT_NEAR(expect[iter], result[iter], 0.0001);
    }
}

TEST_F(Linear_Algebra_Norm_, Normalize_Return_Vec4_Test) {
    vec4 a = {0, 1, 2, 3};

    float * result = vec4_normalize_r(a);

    vec4 expect = {0, 0.2673, 0.5345, 0.8018};
    uint32_t iter;
    for (iter=0; iter<4; ++iter) {
        EXPECT_NEAR(expect[iter], result[iter], 0.0001);
    }
}

TEST_F(Linear_Algebra_Norm_, Normalize_Return_Vec5_Test) {
    vec5 a = {0, 1, 2, 3, 4};

    float * result = vec5_normalize_r(a);

    vec5 expect = {0, 0.1826, 0.3651, 0.5477, 0.7303};
    uint32_t iter;
    for (iter=0; iter<5; ++iter) {
        EXPECT_NEAR(expect[iter], result[iter], 0.0001);
    }
}

TEST_F(Linear_Algebra_Norm_, Normalize_Return_Vec6_Test) {
    vec6 a = {0, 1, 2, 3, 4, 5};

    float * result = vec6_normalize_r(a);

    vec6 expect = {0, 0.1348, 0.2697, 0.4045, 0.5394, 0.6742};
    uint32_t iter;
    for (iter=0; iter<6; ++iter) {
        EXPECT_NEAR(expect[iter], result[iter], 0.0001);
    }
}
