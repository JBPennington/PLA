
extern "C" {
#include "../linear_algebra.h"
}

#include <gtest/gtest.h>
#include <iostream>


class Linear_Algebra_Misc_ : public testing::Test {
    void SetUp() {}
    void TearDown() {}
public:
    const float tolerance = 0.00000001f;
};

TEST_F(Linear_Algebra_Misc_, Vec3_Cross_Multiplication_Set_New_Value_Test) {
    vec3 a = {0, 1, 2};
    vec3 b = {3, 4, 5};

    vec3 result = {0,0,0};
    vec3_mul_cross_n(result, a, b);

    vec3 expect = {-3, 6, -3};

    uint32_t iter;
    for (iter=0; iter<3; ++iter) {
        EXPECT_EQ(expect[iter], result[iter]);
    }
}

TEST_F(Linear_Algebra_Misc_, Vec3_Cross_Multiplication_Set_First_Test) {
    vec3 a = {0, 1, 2};
    vec3 b = {3, 4, 5};

    vec3_mul_cross(a, b);

    vec3 expect = {-3, 6, -3};

    uint32_t iter;
    for (iter=0; iter<3; ++iter) {
        EXPECT_EQ(expect[iter], a[iter]);
    }
}

TEST_F(Linear_Algebra_Misc_, Vec3_Cross_Multiplication_Set_First_And_Return_Test) {
    vec3 a = {0, 1, 2};
    vec3 b = {3, 4, 5};

    float * result = vec3_mul_cross_r(a, b);

    vec3 expect = {-3, 6, -3};

    uint32_t iter;
    for (iter=0; iter<3; ++iter) {
        EXPECT_EQ(expect[iter], result[iter]);
    }
}

TEST_F(Linear_Algebra_Misc_, Vec3_Reflect_Set_New_Value_Test) {
    vec3 a = {1, 0, 0};
    vec3 n = {1, 1, 1};

    vec3 result = {0, 0, 0};
    vec3_reflect_n(result, a, n);

    vec3 expect = {-1, -2, -2};

    uint32_t iter;
    for (iter=0; iter<3; ++iter) {
        EXPECT_NEAR(expect[iter], result[iter], tolerance);
    }
}

TEST_F(Linear_Algebra_Misc_, Vec3_Reflect_Set_First_Test) {
    vec3 a = {1, 0, 0};
    vec3 n = {1, 1, 1};

    vec3_reflect(a, n);

    vec3 expect = {-1, -2, -2};

    uint32_t iter;
    for (iter=0; iter<3; ++iter) {
        EXPECT_NEAR(expect[iter], a[iter], tolerance);
    }
}

TEST_F(Linear_Algebra_Misc_, Vec3_Reflect_Set_First_And_Return_Test) {
    vec3 a = {1, 0, 0};
    vec3 n = {1, 1, 1};

    float * result = vec3_reflect_r(a, n);

    vec3 expect = {-1, -2, -2};

    uint32_t iter;
    for (iter=0; iter<3; ++iter) {
        EXPECT_NEAR(expect[iter], result[iter], tolerance);
    }
}

TEST_F(Linear_Algebra_Misc_, Vec4_Cross_Multiplication_Set_New_Value_Test) {
    vec4 a = {0, 1, 2, 3};
    vec4 b = {4, 5, 6, 7};

    vec4 result = {0,0,0,0};
    vec4_mul_cross_n(result, a, b);

    vec4 expect = {-4, 8, -4, 1};

    uint32_t iter;
    for (iter=0; iter<4; ++iter) {
        EXPECT_EQ(expect[iter], result[iter]);
    }
}

TEST_F(Linear_Algebra_Misc_, Vec4_Cross_Multiplication_Set_First_Test) {
    vec4 a = {0, 1, 2, 3};
    vec4 b = {4, 5, 6, 7};

    vec4_mul_cross(a, b);

    vec4 expect = {-4, 8, -4, 1};

    uint32_t iter;
    for (iter=0; iter<4; ++iter) {
        EXPECT_EQ(expect[iter], a[iter]);
    }
}

TEST_F(Linear_Algebra_Misc_, Vec4_Cross_Multiplication_Set_First_And_Return_Test) {
    vec4 a = {0, 1, 2, 3};
    vec4 b = {4, 5, 6, 7};

    float * result = vec4_mul_cross_r(a, b);

    vec4 expect = {-4, 8, -4, 1};

    uint32_t iter;
    for (iter=0; iter<4; ++iter) {
        EXPECT_EQ(expect[iter], result[iter]);
    }
}

TEST_F(Linear_Algebra_Misc_, Vec4_Reflect_Set_New_Value_Test) {
    vec4 a = {0, 1, 2, 3};
    vec4 n = {1, 1, 1, 1};

    vec4 result = {0, 0, 0, 0};
    vec4_reflect_n(result, a, n);

    vec4 expect = {-12, -11, -10, -9};

    uint32_t iter;
    for (iter=0; iter<4; ++iter) {
        EXPECT_NEAR(expect[iter], result[iter], tolerance);
    }
}

TEST_F(Linear_Algebra_Misc_, Vec4_Reflect_Set_First_Test) {
    vec4 a = {0, 1, 2, 3};
    vec4 n = {1, 1, 1, 1};

    vec4_reflect(a, n);

    vec4 expect = {-12, -11, -10, -9};

    uint32_t iter;
    for (iter=0; iter<4; ++iter) {
        EXPECT_NEAR(expect[iter], a[iter], tolerance);
    }
}

TEST_F(Linear_Algebra_Misc_, Vec4_Reflect_Set_First_And_Return_Test) {
    vec4 a = {0, 1, 2, 3};
    vec4 n = {1, 1, 1, 1};

    float * result = vec4_reflect_r(a, n);

    vec4 expect = {-12, -11, -10, -9};

    uint32_t iter;
    for (iter=0; iter<4; ++iter) {
        EXPECT_NEAR(expect[iter], result[iter], tolerance);
    }
}


TEST_F(Linear_Algebra_Misc_, Vec4_Add_In_Place_And_Return_Cross_Product_Test) {
    vec4 a = { 0,  1,  2,  3};
    vec4 b = { 4,  5,  6,  7};
    vec4 c = { 8,  9, 10, 11};
    vec4 d = {12, 13, 14, 15};

    float * result = vec4_mul_cross_r(vec4_add_r(a, b), vec4_sub_r(d,c));

    vec4 expect = {-8, 16, -8, 1};
    uint32_t iter;
    for (iter=0; iter<4; ++iter) {
        EXPECT_NEAR(expect[iter], result[iter], tolerance);
    }
}

TEST_F(Linear_Algebra_Misc_, Vec4_Add_And_Cross_Product_Test) {
    vec4 a = { 0,  1,  2,  3};
    vec4 b = { 4,  5,  6,  7};
    vec4 c = { 8,  9, 10, 11};
    vec4 d = {12, 13, 14, 15};

    vec4 temp1 = {0, 0, 0, 0};
    vec4 temp2 = {0, 0, 0, 0};
    vec4 temp3 = {0, 0, 0, 0};
    vec4_add_n(temp1, a, b);
    vec4_sub_n(temp2, d, c);
    vec4_mul_cross_n(temp3, temp1, temp2);

    vec4 expect = {-8, 16, -8, 1};
    uint32_t iter;
    for (iter=0; iter<4; ++iter) {
        EXPECT_NEAR(expect[iter], temp3[iter], tolerance);
    }
}
