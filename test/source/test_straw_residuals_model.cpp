#include <cmath>
#include <random>

#include <gtest/gtest.h>

#include "../../source/straw_residual_model.hpp"

static double accuracy = 1e-12;
static size_t n_repeats = 100000;

using ROOT::Math::Rotation3D;
using ROOT::Math::Transform3D;
using ROOT::Math::XYZPoint;
using ROOT::Math::XYZVector;

TEST(StrawResiduals, GlobalLocalTransformation)
{
    // tuple:
    //  { global + corrections, local}
    //  { local_points, global_points }
    //  { local_vectors, global_vectors }
    std::vector<std::tuple<std::pair<std::array<double, 12>, std::array<double, 6>>,
                           std::vector<std::pair<XYZPoint, XYZPoint>>,
                           std::vector<std::pair<XYZVector, XYZVector>>>>
        transform_data = {// no rotations
                          {{{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0}},
                           {{{0, 0, 0}, {0, 0, 0}}, {{0, 0, 0}, {0, 0, 0}}, {{0, 0, 0}, {0, 0, 0}}, {{0, 0, 0}, {0, 0, 0}}},
                           {{{0, 0, 0}, {0, 0, 0}}, {{0, 0, 0}, {0, 0, 0}}, {{0, 0, 0}, {0, 0, 0}}, {{0, 0, 0}, {0, 0, 0}}}},
                          // transform global by x=1
                          {{{1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0}},
                           {{{0, 0, 0}, {1, 0, 0}}, {{0, 1, 0}, {1, 1, 0}}, {{0, 0, 1}, {1, 0, 1}}, {{-5, 0, 0}, {-4, 0, 0}}},
                           {{{0, 0, 0}, {0, 0, 0}}, {{1, 0, 0}, {1, 0, 0}}, {{0, 1, 0}, {0, 1, 0}}, {{1, -1, 2}, {1, -1, 2}}}}};

    for (const auto& d : transform_data) {
        auto model = hsa::sts_residual<double, promille::euler::zyz>(std::get<0>(d).first[0],
                                                                     std::get<0>(d).first[1],
                                                                     std::get<0>(d).first[2],
                                                                     std::get<0>(d).first[3],
                                                                     std::get<0>(d).first[4],
                                                                     std::get<0>(d).first[5],
                                                                     std::get<0>(d).first[6],
                                                                     std::get<0>(d).first[7],
                                                                     std::get<0>(d).first[8],
                                                                     std::get<0>(d).first[9],
                                                                     std::get<0>(d).first[10],
                                                                     std::get<0>(d).first[11]);
        model.set_local_params(std::get<0>(d).second[0],
                               std::get<0>(d).second[1],
                               std::get<0>(d).second[2],
                               std::get<0>(d).second[3],
                               std::get<0>(d).second[4],
                               std::get<0>(d).second[5]);

        for (const auto& t_p : std::get<1>(d)) {
            EXPECT_EQ(t_p.second, model.to_global(t_p.first));
            EXPECT_EQ(t_p.first, model.to_local(t_p.second));
        }

        for (const auto& t_v : std::get<2>(d)) {
            EXPECT_EQ(t_v.second, model.to_global(t_v.first));
            EXPECT_EQ(t_v.first, model.to_local(t_v.second));
        }
    }
}

TEST(StrawResiduals, ResidualSign)
{
    // tuple:
    //  { global + corrections, local}
    //  { track_base_g, track_vector_g }
    //  { straw_base_L, sign_-1_or_1 }
    std::vector<std::tuple<std::pair<std::array<double, 12>, std::array<double, 6>>,
                           std::vector<std::pair<std::tuple<XYZPoint, XYZVector, XYZPoint>, double>>>>
        test_data = {
            // no rotations
            {{{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0}},
             {{{{0, 0, 0}, {0, 0, 0}, {0, 0, 0}}, 1}, {{{0, 0, 0}, {0, 0, 0}, {1, 0, 0}}, -1}, {{{0, 0, 0}, {0, 0, 0}, {-1, 0, 0}}, 1}}}};

    int test_set_counter = 0;
    for (const auto& d : test_data) {
        auto model = hsa::sts_residual<double, promille::euler::zyz>(std::get<0>(d).first[0],
                                                                     std::get<0>(d).first[1],
                                                                     std::get<0>(d).first[2],
                                                                     std::get<0>(d).first[3],
                                                                     std::get<0>(d).first[4],
                                                                     std::get<0>(d).first[5],
                                                                     std::get<0>(d).first[6],
                                                                     std::get<0>(d).first[7],
                                                                     std::get<0>(d).first[8],
                                                                     std::get<0>(d).first[9],
                                                                     std::get<0>(d).first[10],
                                                                     std::get<0>(d).first[11]);
        model.set_local_params(std::get<0>(d).second[0],
                               std::get<0>(d).second[1],
                               std::get<0>(d).second[2],
                               std::get<0>(d).second[3],
                               std::get<0>(d).second[4],
                               std::get<0>(d).second[5]);

        int point_set_counter = 0;
        for (const auto& tries : std::get<1>(d)) {
            EXPECT_EQ(model.calc_residual_sign(std::get<0>(tries.first), std::get<1>(tries.first), std::get<2>(tries.first)), tries.second)
                << " for test_set_counter=" << test_set_counter << " point_set_counter=" << point_set_counter;

            point_set_counter++;
        }

        test_set_counter++;
    }
}

/*
TEST(StrawResiduals, StsDerivatives)
{
    std::vector<std::tuple<std::array<double, 18>, std::array<double, 7>, double, std::array<double, 12>, std::array<double, 4>>> data = {
        {{1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
         {0, 0, 0, 0, 0, 0, 0},
         1,
         {1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0},
         {-1, 0, 0, 0}},
        {{0, 1, 2, 0, 0, 0, 6, 7, 8, 0, 0, 0, 12, 13, 14, 0, 0, 0},
         {18, 19, 20, 21, 22, 23, 24},
         43.3069,
         {-0.0434372, 0.0, 0.999056, -18.5975, 17.1143, 0.808588, -0.0434372, 0.0, 0.999056, -18.5975, 17.1143, 0.808588},
         {0.0434372, 0.0, 0.0318813, 0.0}},
        {{0, 1, 2, 0.03, 0.04, 0.05, 6, 7, 8, 0, 0, 0, 12, 13, 14, 0.015, 0.016, 0.017},
         {18, 19, 20, 21, 22, 23, 24},
         43.404,
         {-0.0847495, 0.0397345, 0.99561, -21.243, 16.1736, 1.16278, -0.0847495, 0.0397345, 0.99561, -21.7111, 16.9432, 0.784717},
         {0.0847495, -0.0397345, 0.0674815, -0.0316385}}};

    int test_counter = 0;
    for (const auto& d : data) {
        auto derivs = hsa::sts_residual<double, promille::euler::zyz>(std::get<0>(d)[0],
                                                                      std::get<0>(d)[1],
                                                                      std::get<0>(d)[2],
                                                                      std::get<0>(d)[3],
                                                                      std::get<0>(d)[4],
                                                                      std::get<0>(d)[5],
                                                                      std::get<0>(d)[12],
                                                                      std::get<0>(d)[13],
                                                                      std::get<0>(d)[14],
                                                                      std::get<0>(d)[15],
                                                                      std::get<0>(d)[16],
                                                                      std::get<0>(d)[17]);

        derivs.set_local_params(
            std::get<0>(d)[6], std::get<0>(d)[7], std::get<0>(d)[8], std::get<0>(d)[9], std::get<0>(d)[10], std::get<0>(d)[11]);

        derivs.recalculate(XYZPoint(std::get<1>(d)[3], std::get<1>(d)[4], 0),
                           XYZVector(std::get<1>(d)[5], std::get<1>(d)[6], 1),
                           XYZPoint(std::get<1>(d)[0], std::get<1>(d)[1], std::get<1>(d)[2]),
                           XYZVector(0, 1, 0),
                           0);

        auto residua = derivs.residual();

        EXPECT_NEAR(residua, std::get<2>(d), accuracy) << "residua at test " << test_counter << derivs;

        EXPECT_NEAR(derivs.global_derivative(0), std::get<3>(d)[0], accuracy) << "dr/dxg at test " << test_counter;
        EXPECT_NEAR(derivs.global_derivative(1), std::get<3>(d)[1], accuracy) << "dr/dyg at test " << test_counter;
        EXPECT_NEAR(derivs.global_derivative(2), std::get<3>(d)[2], accuracy) << "dr/dzg at test " << test_counter;

        EXPECT_NEAR(derivs.global_derivative(3), std::get<3>(d)[3], accuracy) << "dr/dag at test " << test_counter;
        EXPECT_NEAR(derivs.global_derivative(4), std::get<3>(d)[4], accuracy) << "dr/dbg at test " << test_counter;
        EXPECT_NEAR(derivs.global_derivative(5), std::get<3>(d)[5], accuracy) << "dr/dcg at test " << test_counter;

        EXPECT_NEAR(derivs.global_derivative(6), std::get<3>(d)[6], accuracy) << "dr/dxc at test " << test_counter;
        EXPECT_NEAR(derivs.global_derivative(7), std::get<3>(d)[7], accuracy) << "dr/dyc at test " << test_counter;
        EXPECT_NEAR(derivs.global_derivative(8), std::get<3>(d)[8], accuracy) << "dr/dzc at test " << test_counter;

        EXPECT_NEAR(derivs.global_derivative(9), std::get<3>(d)[9], accuracy) << "dr/dac at test " << test_counter;
        EXPECT_NEAR(derivs.global_derivative(10), std::get<3>(d)[10], accuracy) << "dr/dbc at test " << test_counter;
        EXPECT_NEAR(derivs.global_derivative(11), std::get<3>(d)[11], accuracy) << "dr/dcc at test " << test_counter;

        EXPECT_NEAR(derivs.local_derivative(0), std::get<4>(d)[0], accuracy) << "dr/dbx at test " << test_counter;
        EXPECT_NEAR(derivs.local_derivative(1), std::get<4>(d)[1], accuracy) << "dr/dby at test " << test_counter;

        EXPECT_NEAR(derivs.local_derivative(2), std::get<4>(d)[2], accuracy) << "dr/dtx at test " << test_counter;
        EXPECT_NEAR(derivs.local_derivative(3), std::get<4>(d)[3], accuracy) << "dr/dtx at test " << test_counter;

        test_counter++;
    }
}
*/

TEST(StrawResiduals, PlanesResiduasX)
{
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis_straw_u(-200.0, 200.0);
    std::uniform_real_distribution<> dis_straw_z(2000.0, 5000.0);

    std::uniform_real_distribution<> dis_track_base(-20.0, 20.0);
    std::uniform_real_distribution<> dis_track_dir(-0.05, 0.05);

    auto the_model = hsa::sts_residual<double, promille::euler::zyz>(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);

    const auto su = 0.0;
    const auto sz = 0.0;
    const auto xl = 0.0;
    const auto yl = 0.0;
    const auto zl = 0.0;
    const auto al = 0.0;
    const auto bl = 0.0;
    const auto cl = 0.0;
    const auto bx = 0.0;
    const auto by = 0.0;
    const auto tx = 0.0;
    const auto ty = 0.0;

    the_model.set_local_params(xl, yl, zl, al, bl, cl);
    the_model.recalculate({bx, by, 0}, {tx, ty, 1.0}, {su, 0, sz}, {0, 1, 0}, 0.0);

    const auto expected_value = fabs(su + xl - bx - tx * (sz + zl)) / sqrt(1.0 + tx * tx);

    EXPECT_EQ(the_model.residual(), expected_value);

    for (size_t i = 0; i < n_repeats; ++i) {
        const auto su = dis_straw_u(gen);
        const auto sz = dis_straw_z(gen);
        const auto xl = 0.0;
        const auto yl = 0.0;
        const auto zl = 0.0;
        const auto al = 0.0;
        const auto bl = 0.0;
        const auto cl = 0.0;
        const auto bx = dis_track_base(gen);
        const auto by = 0.0;
        const auto tx = dis_track_dir(gen);
        const auto ty = 0.0;

        the_model.set_local_params(xl, yl, zl, al, bl, cl);
        the_model.recalculate({bx, by, 0}, {tx, ty, 1.0}, {su, 0, sz}, {0, 1, 0}, 0.0);

        const auto nom = -bx + su + xl - tx * (sz + zl);
        const auto denom_1 = sqrt(1.0 + tx * tx);

        const auto expected_doca = fabs(nom) / denom_1;

        EXPECT_NEAR(fabs(the_model.residual()), expected_doca, accuracy) << i << the_model;

        const auto expected_drdbx = -nom / (denom_1 * fabs(nom));

        EXPECT_NEAR(the_model.local_derivative(0), expected_drdbx, accuracy) << i << the_model;
        EXPECT_NEAR(the_model.local_derivative(1), 0, accuracy) << i << the_model;

        const auto expected_drdtx = (-nom * (sz + zl)) / (denom_1 * fabs(nom)) - (tx * fabs(nom)) / pow(1 + tx * tx, 3.0 / 2.0);

        EXPECT_NEAR(the_model.local_derivative(2), expected_drdtx, accuracy) << i << the_model;
        EXPECT_NEAR(the_model.local_derivative(3), 0, accuracy) << i << the_model;

        const auto expected_drdxg = nom / (denom_1 * fabs(nom));
        const auto expected_drdzg = -(tx * nom) / (denom_1 * fabs(nom));

        EXPECT_NEAR(the_model.global_derivative(0), expected_drdxg, accuracy) << i << the_model;
        EXPECT_NEAR(the_model.global_derivative(1), 0, accuracy) << i << the_model;
        EXPECT_NEAR(the_model.global_derivative(2), expected_drdzg, accuracy) << i << the_model;

        const auto expected_drdbg = (sz - su * tx) * nom / (denom_1 * fabs(nom));

        EXPECT_NEAR(the_model.global_derivative(3), 0, accuracy) << i << the_model;
        EXPECT_NEAR(the_model.global_derivative(4), expected_drdbg, accuracy) << i << the_model;
        EXPECT_NEAR(the_model.global_derivative(5), 0, accuracy) << i << the_model;

        const auto expected_drdxc = nom / (denom_1 * fabs(nom));
        const auto expected_drdzc = -(tx * nom) / (denom_1 * fabs(nom));

        EXPECT_NEAR(the_model.global_derivative(6), expected_drdxc, accuracy) << i << the_model;
        EXPECT_NEAR(the_model.global_derivative(7), 0, accuracy) << i << the_model;
        EXPECT_NEAR(the_model.global_derivative(8), expected_drdzc, accuracy) << i << the_model;

        const auto expected_drdbc = (sz - su * tx) * nom / (denom_1 * fabs(nom));

        EXPECT_NEAR(the_model.global_derivative(9), 0, accuracy) << i << the_model;
        EXPECT_NEAR(the_model.global_derivative(10), expected_drdbc, accuracy) << i << the_model;
        EXPECT_NEAR(the_model.global_derivative(11), 0, accuracy) << i << the_model;
    }
}

TEST(StrawResiduals, PlanesResiduasY)
{
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis_straw_u(-200.0, 200.0);
    std::uniform_real_distribution<> dis_straw_z(2000.0, 5000.0);

    std::uniform_real_distribution<> dis_track_base(-20.0, 20.0);
    std::uniform_real_distribution<> dis_track_dir(-0.05, 0.05);

    auto the_model = hsa::sts_residual<double, promille::euler::zyz>(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);

    const auto su = 0.0;
    const auto sz = 0.0;
    const auto xl = 0.0;
    const auto yl = 0.0;
    const auto zl = 0.0;
    const auto al = 0.0;
    const auto bl = 0.0;
    const auto cl = M_PI_2;
    const auto bx = 0.0;
    const auto by = 0.0;
    const auto tx = 0.0;
    const auto ty = 0.0;

    the_model.set_local_params(xl, yl, zl, al, bl, cl);
    the_model.recalculate({bx, by, 0}, {tx, ty, 1.0}, {su, 0, sz}, {0, 1, 0}, 0.0);

    const auto expected_value = fabs(su + yl - by - ty * (sz + zl)) / sqrt(1.0 + ty * ty);

    EXPECT_EQ(the_model.residual(), expected_value);

    for (size_t i = 0; i < n_repeats; ++i) {
        const auto su = dis_straw_u(gen);
        const auto sz = dis_straw_z(gen);
        const auto xl = 0.0;
        const auto yl = 0.0;
        const auto zl = 0.0;
        const auto al = 0.0;
        const auto bl = 0.0;
        const auto cl = M_PI_2;
        const auto bx = 0.0;
        const auto by = dis_track_base(gen);
        const auto tx = 0.0;
        const auto ty = dis_track_dir(gen);

        the_model.set_local_params(xl, yl, zl, al, bl, cl);
        the_model.recalculate({bx, by, 0}, {tx, ty, 1.0}, {su, 0, sz}, {0, 1, 0}, 0.0);

        const auto nom = -by + su + yl - ty * (sz + zl);
        const auto denom_1 = sqrt(1.0 + ty * ty);

        const auto expected_doca = fabs(nom) / denom_1;

        EXPECT_NEAR(fabs(the_model.residual()), expected_doca, accuracy) << i << the_model;

        const auto expected_drdby = -nom / (denom_1 * fabs(nom));

        EXPECT_NEAR(the_model.local_derivative(0), 0, accuracy) << i << the_model;
        EXPECT_NEAR(the_model.local_derivative(1), expected_drdby, accuracy) << i << the_model;

        const auto expected_drdty = (-nom * (sz + zl)) / (denom_1 * fabs(nom)) - (ty * fabs(nom)) / pow(1 + ty * ty, 3.0 / 2.0);

        EXPECT_NEAR(the_model.local_derivative(2), 0, accuracy) << i << the_model;
        EXPECT_NEAR(the_model.local_derivative(3), expected_drdty, accuracy) << i << the_model;

        const auto expected_drdyg = nom / (denom_1 * fabs(nom));
        const auto expected_drdzg = -(ty * nom) / (denom_1 * fabs(nom));

        EXPECT_NEAR(the_model.global_derivative(0), 0, accuracy) << i << the_model;
        EXPECT_NEAR(the_model.global_derivative(1), expected_drdyg, accuracy) << i << the_model;
        EXPECT_NEAR(the_model.global_derivative(2), expected_drdzg, accuracy) << i << the_model;

        const auto expected_drdag = (sz + su * ty) * nom / (denom_1 * fabs(nom));

        EXPECT_NEAR(the_model.global_derivative(3), expected_drdag, accuracy) << i << the_model;
        EXPECT_NEAR(the_model.global_derivative(4), 0, accuracy) << i << the_model;
        EXPECT_NEAR(the_model.global_derivative(5), 0, accuracy) << i << the_model;

        const auto expected_drdyc = nom / (denom_1 * fabs(nom));
        const auto expected_drdzc = -(ty * nom) / (denom_1 * fabs(nom));

        EXPECT_NEAR(the_model.global_derivative(6), 0, accuracy) << i << the_model;
        EXPECT_NEAR(the_model.global_derivative(7), expected_drdyc, accuracy) << i << the_model;
        EXPECT_NEAR(the_model.global_derivative(8), expected_drdzc, accuracy) << i << the_model;

        const auto expected_drdbc = (sz - su * ty) * nom / (denom_1 * fabs(nom));

        EXPECT_NEAR(the_model.global_derivative(9), 0, accuracy) << i << the_model;
        EXPECT_NEAR(the_model.global_derivative(10), expected_drdbc, accuracy) << i << the_model;
        EXPECT_NEAR(the_model.global_derivative(11), 0, accuracy) << i << the_model;
    }
}
