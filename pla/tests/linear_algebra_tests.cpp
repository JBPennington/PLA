
extern "C" {
#include "../pla.h"
}

#include <gtest/gtest.h>
#include <iostream>


class Linear_Algebra_ : public testing::Test {

    void SetUp() {}

    void TearDown() {}
};


//TEST_F(Linear_Algebra_, Multiply_4_By_4_Matrix_Test) {
//    mat4x4 A = {{0, 1, 2, 3}, {4, 5, 6, 7}, {8, 9, 10, 11}, {12, 13, 14, 15}};
//    mat4x4 B = {{0, 1, 2, 3}, {4, 5, 6, 7}, {8, 9, 10, 11}, {12, 13, 14, 15}};
//
//    mat4x4 M;
//
//    mat4x4_mul(M, A, B);
//
//    mat4x4 expected = {{56, 62, 68, 74}, {152, 174, 196, 218}, {248,286, 324, 362}, {344, 398, 452, 506}};
//    uint32_t i = 0;
//    uint32_t j = 0;
//    for ( i=0; i<4; i++) {
//        for ( j=0; j<4; j++) {
//            EXPECT_EQ(expected[i][j], M[i][j]);
//        }
//    }
//}
//
//
//TEST_F(Linear_Algebra_, Multiply_3_By_3_Matrix_Test) {
//    mat3x3 A = {{0, 1, 2}, {3, 4, 5}, {6, 7, 8}};
//    mat3x3 B = {{0, 1, 2}, {3, 4, 5}, {6, 7, 8}};
//
//    mat3x3 M;
//
//    mat3x3_mul(M, A, B);
//
//    mat3x3 expected = {{15, 18, 21}, {42, 54, 66}, {69, 90, 111}};
//    uint32_t i = 0;
//    uint32_t j = 0;
//    for ( i=0; i<3; i++) {
//        for ( j=0; j<3; j++) {
//            EXPECT_EQ(expected[i][j], M[i][j]);
//        }
//    }
//}
//
//
//TEST_F(Linear_Algebra_, Multiply_6_By_6_Matrix_by_6_Vec_Test) {
//    mat6x6 A = {{  0,  1,  2,  3,  4,  5},
//                {  6,  7,  8,  9, 10, 11},
//                { 12, 13, 14, 15, 16, 17},
//                { 18, 19, 20, 21, 22, 23},
//                { 24, 25, 26, 27, 28, 29},
//                { 30, 31, 32, 33, 34, 35}};
//
//    vec6 B = {0, 1, 2, 3, 4, 5};
//
//    vec6 answer;
//
//    mat6x6_mul_vec6(answer, A, B);
//
//    vec6 expected = {55, 145, 235, 325, 415, 505};
//    uint32_t i = 0;
//    for ( i=0; i<6; i++) {
//        EXPECT_EQ(expected[i], answer[i]);
//    }
//}
//
//
//TEST_F(Linear_Algebra_, Get_Rotation_Matrix_Test) {
//    mat4x4 A = {{0, 1, 2, 3}, {4, 5, 6, 7}, {8, 9, 10, 11}, {12, 13, 14, 15}};
//
//    mat3x3 M;
//
//    mat4x4_get_rotational(M, A);
//
//    mat4x4 expected = {{0, 1, 2}, {4, 5, 6}, {8, 9, 10}};
//    uint32_t i = 0;
//    uint32_t j = 0;
//    for ( i=0; i<3; i++) {
//        for ( j=0; j<3; j++) {
//            EXPECT_EQ(expected[i][j], M[i][j]);
//        }
//    }
//}
//
//TEST_F(Linear_Algebra_, Multiply_3_By_3_Matrix_By_Vec3_Test) {
//    mat3x3  A = {{0, 1, 2}, {3, 4, 5}, {6, 7, 8}};
//    vec3    B = {0,1,2};
//    vec3    M;
//
//    mat3x3_mul_vec3(M, A, B);
//
//    vec3    expected = {5, 14, 23};
//    uint32_t i = 0;
//    for ( i=0; i<3; i++) {
//        EXPECT_EQ(expected[i], M[i]);
//    }
//}
//
//
