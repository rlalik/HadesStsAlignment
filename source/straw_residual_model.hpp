#pragma once

// mille builder
#include <promille/euler_angles.hpp>
#include <promille/promille.hpp>

// ROOT
#include <Math/Point3D.h>
#include <Math/Rotation3D.h>
#include <Math/Transform3D.h>
#include <Math/Vector3D.h>

/*
 * Rotation matrix, see https://mathworld.wolfram.com/EulerAngles.html (48..50)
 * for definitions.
 */

namespace hsa
{

using ROOT::Math::Rotation3D;
using ROOT::Math::Transform3D;
using ROOT::Math::XYZPoint;
using ROOT::Math::XYZVector;

template<typename T, template<class> class R, size_t Nlocals = 4, size_t Nglobals = 12>
struct sts_residual final : promille::residual_model_base<T, Nlocals, Nglobals>
{
    // global translational corrections
    T x_g {0};
    T y_g {0};
    T z_g {0};

    // global rotational corrections
    T a_g {0};
    T b_g {0};
    T c_g {0};

    // global translational alignment
    T x_l {0};
    T y_l {0};
    T z_l {0};

    // global rotational alignment
    T a_l {0};
    T b_l {0};
    T c_l {0};

    // translational corrections
    T x_c {0};
    T y_c {0};
    T z_c {0};

    // rotational corrections
    T a_c {0};
    T b_c {0};
    T c_c {0};

    R<T> wm;

    Transform3D global_transform;
    Transform3D local_transform;
    Transform3D correction_transform;

    Transform3D total_transform;
    Transform3D inverse_total_transform;

    XYZPoint straw_base_loc;
    XYZVector straw_dir_loc;

    XYZPoint straw_base_lab;
    XYZVector straw_dir_lab;

    XYZPoint track_base;
    XYZVector track_dir;

    T drift_radius {0};
    T track_wire_distance {0};
    T res_sign {0};

    sts_residual(T gx, T gy, T gz, T ga, T gb, T gc, T cx, T cy, T cz, T ca, T cb, T cc)
        : x_g(gx)
        , y_g(gy)
        , z_g(gz)
        , a_g(ga)
        , b_g(gb)
        , c_g(gc)
        , x_c(cx)
        , y_c(cy)
        , z_c(cz)
        , a_c(ca)
        , b_c(cb)
        , c_c(cc)
        , wm(0, 0, 0)
        , global_transform(Rotation3D(1, -c_g, b_g, c_g, 1, -a_g, -b_g, a_g, 1), XYZVector(x_g, y_g, z_g))
        , correction_transform(Rotation3D(1, -c_c, b_c, c_c, 1, -a_c, -b_c, a_c, 1), XYZVector(x_c, y_c, z_c))
    {
    }

    auto set_local_params(T lx, T ly, T lz, T la, T lb, T lc) -> void
    {
        x_l = lx;
        y_l = ly;
        z_l = lz;
        a_l = la;
        b_l = lb;
        c_l = lc;

        wm = R<T>(a_l, b_l, c_l);
        local_transform = Transform3D(promille::euler::make_rotation_matrix(wm), XYZVector(x_l, y_l, z_l));
        // local_rotation = promille::euler::make_rotation_matrix(wm);
        // local_translation = {x_l, y_l, z_l};

        // total_rotation = global_rotation * local_rotation * correction_rotation;
        // total_translation = global_translation + local_translation + correction_translation;

        total_transform = global_transform * local_transform * correction_transform;
        inverse_total_transform = total_transform.Inverse();
        // total_translation = global_translation + local_translation + correction_translation;
    }

    auto print_params() const -> void
    {
        // std::cout << "=====\nGLOBAL: Tg" << global_translation << "  Rg:" << global_rotation << "        Tl" << local_translation
        //           << "  Rl: " << local_rotation << "        Tc" << correction_translation << "  Rc:" << correction_rotation
        //           << "\nTOTAL : Tt" << total_translation << "  Rt" << total_rotation << "\nLOCAL : Sb" << straw_base_loc << "->"
        //           << straw_base_lab << "  Sd" << straw_dir_loc << "->" << straw_dir_lab
        //           << "\n        Tb" << track_base << "  Td" << track_dir << "  dr=" << drift_radius << "  |T-W|=" << track_wire_distance
        //           << "\n=====\n";

        std::cout << "=====\nGT: " << global_transform << "\nLT: " << local_transform << "\nCT: " << correction_transform
                  << "\nTT:" << total_transform << "\nIT:" << inverse_total_transform << "\nLOCAL : Sb" << straw_base_loc << "->"
                  << straw_base_lab << "  Sd" << straw_dir_loc << "->" << straw_dir_lab << "\n        Tb" << track_base << "  Td"
                  << track_dir << "  dr=" << drift_radius << "  |T-W|=" << track_wire_distance << "\n=====\n";

        this->print();

        std::cout << "\n=====\n";
    }

    template<typename Type>
    auto to_global(Type t) const -> Type
    {
        return total_transform * t;
    }

    template<typename Type>
    auto to_local(Type t) const -> Type
    {
        return inverse_total_transform * t;
    }

    /** We need to determine sign of residual. By convention if track is on the left from wire, the residul are positive, if on the right,
     * then negative.
     *
     * One need to rotate track into module local coordination system, and then compare the x of track intersection at plane=z with c of
     * straw.
     * @param track_base base of track in global coordinate system
     * @param track_dir direction of track in global coordinate system, requires z()==1
     * @param straw_bsae base of track in local coordinate system
     * @return -1 or 1 depends on the relative coordinate
     */
    auto calc_residual_sign(XYZPoint track_base, XYZVector track_dir, XYZPoint straw_base) -> T
    {
        const auto l_track_base = to_local(track_base);
        const auto l_track_dir = to_local(track_dir);

        const auto intersection = l_track_base + l_track_dir * (straw_base.z() - l_track_base.z());

        if (intersection.x() < straw_base.x())
            return -1;
        else
            return 1;
    }

    static auto calc_track_wire_distance(XYZPoint track_base, XYZVector track_dir, XYZPoint straw_base_lab, XYZVector straw_dir_lab) -> T
    {
        return std::fabs((track_base - straw_base_lab).Dot(track_dir.Cross(straw_dir_lab)) / track_dir.Cross(straw_dir_lab).R());
    }

    auto residual() const -> T override { return res_sign * (track_wire_distance - drift_radius); }

    auto recalculate(XYZPoint track_base, XYZVector track_dir, XYZPoint straw_base, XYZVector straw_dir, T dr) -> void
    {
        straw_base_loc = straw_base;
        straw_dir_loc = straw_dir;

        this->track_base = track_base;
        this->track_dir = track_dir;
        this->drift_radius = drift_radius;

        // straw_base_lab = total_rotation * straw_base + total_translation;
        // straw_dir_lab = total_rotation * straw_dir;

        straw_base_lab = total_transform * straw_base;
        straw_dir_lab = total_transform * straw_dir;

        drift_radius = dr;
        track_wire_distance = calc_track_wire_distance(track_base, track_dir, straw_base_lab, straw_dir_lab);
        res_sign = calc_residual_sign(track_base, track_dir, straw_base);

        auto su = straw_base_loc.X();
        auto sv = straw_base_loc.Y();
        auto sz = straw_base_loc.Z();
        auto bx = track_base.X();
        auto by = track_base.Y();
        auto bz = track_base.Z();
        auto tx = track_dir.X();
        auto ty = track_dir.Y();
        auto tz = track_dir.Z();

        // auto a_c__R11 = a_c * wm.R11;
        // auto a_c__R12 = a_c * wm.R12;
        // auto a_c__R13 = a_c * wm.R13;
        // auto a_c__R21 = a_c * wm.R21;
        // auto a_c__R22 = a_c * wm.R22;
        // auto a_c__R23 = a_c * wm.R23;
        // auto a_c__R31 = a_c * wm.R31;
        // auto a_c__R32 = a_c * wm.R32;
        // auto a_c__R33 = a_c * wm.R33;
        // auto b_c__R11 = b_c * wm.R11;
        // auto b_c__R12 = b_c * wm.R12;
        // auto b_c__R13 = b_c * wm.R13;
        // auto b_c__R21 = b_c * wm.R21;
        // auto b_c__R22 = b_c * wm.R22;
        // auto b_c__R23 = b_c * wm.R23;
        // auto b_c__R31 = b_c * wm.R31;
        // auto b_c__R32 = b_c * wm.R32;
        // auto b_c__R33 = b_c * wm.R33;
        // auto c_c__R11 = c_c * wm.R11;
        // auto c_c__R12 = c_c * wm.R12;
        // auto c_c__R13 = c_c * wm.R13;
        // auto c_c__R21 = c_c * wm.R21;
        // auto c_c__R22 = c_c * wm.R22;
        // auto c_c__R23 = c_c * wm.R23;
        // auto c_c__R31 = c_c * wm.R31;
        // auto c_c__R32 = c_c * wm.R32;
        // auto c_c__R33 = c_c * wm.R33;
        // auto x_c__R11 = x_c * wm.R11;
        // auto x_c__R12 = x_c * wm.R12;
        // auto x_c__R13 = x_c * wm.R13;
        // auto x_c__R21 = x_c * wm.R21;
        // auto x_c__R22 = x_c * wm.R22;
        // auto x_c__R23 = x_c * wm.R23;
        // auto x_c__R31 = x_c * wm.R31;
        // auto x_c__R32 = x_c * wm.R32;
        // auto x_c__R33 = x_c * wm.R33;
        // auto y_c__R11 = y_c * wm.R11;
        // auto y_c__R12 = y_c * wm.R12;
        // auto y_c__R13 = y_c * wm.R13;
        // auto y_c__R21 = y_c * wm.R21;
        // auto y_c__R22 = y_c * wm.R22;
        // auto y_c__R23 = y_c * wm.R23;
        // auto y_c__R31 = y_c * wm.R31;
        // auto y_c__R32 = y_c * wm.R32;
        // auto y_c__R33 = y_c * wm.R33;
        // auto z_c__R11 = z_c * wm.R11;
        // auto z_c__R12 = z_c * wm.R12;
        // auto z_c__R13 = z_c * wm.R13;
        // auto z_c__R21 = z_c * wm.R21;
        // auto z_c__R22 = z_c * wm.R22;
        // auto z_c__R23 = z_c * wm.R23;
        // auto z_c__R31 = z_c * wm.R31;
        // auto z_c__R32 = z_c * wm.R32;
        // auto z_c__R33 = z_c * wm.R33;
        // auto a_g__R11 = a_g * wm.R11;
        // auto a_g__R12 = a_g * wm.R12;
        // auto a_g__R13 = a_g * wm.R13;
        // auto a_g__R21 = a_g * wm.R21;
        // auto a_g__R22 = a_g * wm.R22;
        // auto a_g__R23 = a_g * wm.R23;
        // auto a_g__R31 = a_g * wm.R31;
        // auto a_g__R32 = a_g * wm.R32;
        // auto a_g__R33 = a_g * wm.R33;
        // auto b_g__R11 = b_g * wm.R11;
        // auto b_g__R12 = b_g * wm.R12;
        // auto b_g__R13 = b_g * wm.R13;
        // auto b_g__R21 = b_g * wm.R21;
        // auto b_g__R22 = b_g * wm.R22;
        // auto b_g__R23 = b_g * wm.R23;
        // auto b_g__R31 = b_g * wm.R31;
        // auto b_g__R32 = b_g * wm.R32;
        // auto b_g__R33 = b_g * wm.R33;
        // auto c_g__R11 = c_g * wm.R11;
        // auto c_g__R12 = c_g * wm.R12;
        // auto c_g__R13 = c_g * wm.R13;
        // auto c_g__R21 = c_g * wm.R21;
        // auto c_g__R22 = c_g * wm.R22;
        // auto c_g__R23 = c_g * wm.R23;
        // auto c_g__R31 = c_g * wm.R31;
        // auto c_g__R32 = c_g * wm.R32;
        // auto c_g__R33 = c_g * wm.R33;
        // auto x_g__R11 = x_g * wm.R11;
        // auto x_g__R12 = x_g * wm.R12;
        // auto x_g__R13 = x_g * wm.R13;
        // auto x_g__R21 = x_g * wm.R21;
        // auto x_g__R22 = x_g * wm.R22;
        // auto x_g__R23 = x_g * wm.R23;
        // auto x_g__R31 = x_g * wm.R31;
        // auto x_g__R32 = x_g * wm.R32;
        // auto x_g__R33 = x_g * wm.R33;
        // auto y_g__R11 = y_g * wm.R11;
        // auto y_g__R12 = y_g * wm.R12;
        // auto y_g__R13 = y_g * wm.R13;
        // auto y_g__R21 = y_g * wm.R21;
        // auto y_g__R22 = y_g * wm.R22;
        // auto y_g__R23 = y_g * wm.R23;
        // auto y_g__R31 = y_g * wm.R31;
        // auto y_g__R32 = y_g * wm.R32;
        // auto y_g__R33 = y_g * wm.R33;
        // auto z_g__R11 = z_g * wm.R11;
        // auto z_g__R12 = z_g * wm.R12;
        // auto z_g__R13 = z_g * wm.R13;
        // auto z_g__R21 = z_g * wm.R21;
        // auto z_g__R22 = z_g * wm.R22;
        // auto z_g__R23 = z_g * wm.R23;
        // auto z_g__R31 = z_g * wm.R31;
        // auto z_g__R32 = z_g * wm.R32;
        // auto z_g__R33 = z_g * wm.R33;
        // auto tx__R11 = tx * wm.R11;
        // auto tx__R12 = tx * wm.R12;
        // auto tx__R13 = tx * wm.R13;
        // auto tx__R21 = tx * wm.R21;
        // auto tx__R22 = tx * wm.R22;
        // auto tx__R23 = tx * wm.R23;
        // auto tx__R31 = tx * wm.R31;
        // auto tx__R32 = tx * wm.R32;
        // auto tx__R33 = tx * wm.R33;
        // auto ty__R11 = ty * wm.R11;
        // auto ty__R12 = ty * wm.R12;
        // auto ty__R13 = ty * wm.R13;
        // auto ty__R21 = ty * wm.R21;
        // auto ty__R22 = ty * wm.R22;
        // auto ty__R23 = ty * wm.R23;
        // auto ty__R31 = ty * wm.R31;
        // auto ty__R32 = ty * wm.R32;
        // auto ty__R33 = ty * wm.R33;
        // auto tz__R11 = tz * wm.R11;
        // auto tz__R12 = tz * wm.R12;
        // auto tz__R13 = tz * wm.R13;
        // auto tz__R21 = tz * wm.R21;
        // auto tz__R22 = tz * wm.R22;
        // auto tz__R23 = tz * wm.R23;
        // auto tz__R31 = tz * wm.R31;
        // auto tz__R32 = tz * wm.R32;
        // auto tz__R33 = tz * wm.R33;
        // auto bx__R11 = bx * wm.R11;
        // auto bx__R12 = bx * wm.R12;
        // auto bx__R13 = bx * wm.R13;
        // auto bx__R21 = bx * wm.R21;
        // auto bx__R22 = bx * wm.R22;
        // auto bx__R23 = bx * wm.R23;
        // auto bx__R31 = bx * wm.R31;
        // auto bx__R32 = bx * wm.R32;
        // auto bx__R33 = bx * wm.R33;
        // auto by__R11 = by * wm.R11;
        // auto by__R12 = by * wm.R12;
        // auto by__R13 = by * wm.R13;
        // auto by__R21 = by * wm.R21;
        // auto by__R22 = by * wm.R22;
        // auto by__R23 = by * wm.R23;
        // auto by__R31 = by * wm.R31;
        // auto by__R32 = by * wm.R32;
        // auto by__R33 = by * wm.R33;
        // auto bz__R11 = bz * wm.R11;
        // auto bz__R12 = bz * wm.R12;
        // auto bz__R13 = bz * wm.R13;
        // auto bz__R21 = bz * wm.R21;
        // auto bz__R22 = bz * wm.R22;
        // auto bz__R23 = bz * wm.R23;
        // auto bz__R31 = bz * wm.R31;
        // auto bz__R32 = bz * wm.R32;
        // auto bz__R33 = bz * wm.R33;
        // auto a_g__a_c__R11 = a_g * a_c__R11;
        // auto a_g__a_c__R12 = a_g * a_c__R12;
        // auto a_g__a_c__R13 = a_g * a_c__R13;
        // auto a_g__a_c__R21 = a_g * a_c__R21;
        // auto a_g__a_c__R22 = a_g * a_c__R22;
        // auto a_g__a_c__R23 = a_g * a_c__R23;
        // auto a_g__a_c__R31 = a_g * a_c__R31;
        // auto a_g__a_c__R32 = a_g * a_c__R32;
        // auto a_g__a_c__R33 = a_g * a_c__R33;
        // auto a_g__b_c__R11 = a_g * b_c__R11;
        // auto a_g__b_c__R12 = a_g * b_c__R12;
        // auto a_g__b_c__R13 = a_g * b_c__R13;
        // auto a_g__b_c__R21 = a_g * b_c__R21;
        // auto a_g__b_c__R22 = a_g * b_c__R22;
        // auto a_g__b_c__R23 = a_g * b_c__R23;
        // auto a_g__b_c__R31 = a_g * b_c__R31;
        // auto a_g__b_c__R32 = a_g * b_c__R32;
        // auto a_g__b_c__R33 = a_g * b_c__R33;
        // auto a_g__c_c__R11 = a_g * c_c__R11;
        // auto a_g__c_c__R12 = a_g * c_c__R12;
        // auto a_g__c_c__R13 = a_g * c_c__R13;
        // auto a_g__c_c__R21 = a_g * c_c__R21;
        // auto a_g__c_c__R22 = a_g * c_c__R22;
        // auto a_g__c_c__R23 = a_g * c_c__R23;
        // auto a_g__c_c__R31 = a_g * c_c__R31;
        // auto a_g__c_c__R32 = a_g * c_c__R32;
        // auto a_g__c_c__R33 = a_g * c_c__R33;
        // auto a_g__x_c__R11 = a_g * x_c__R11;
        // auto a_g__x_c__R12 = a_g * x_c__R12;
        // auto a_g__x_c__R13 = a_g * x_c__R13;
        // auto a_g__x_c__R21 = a_g * x_c__R21;
        // auto a_g__x_c__R22 = a_g * x_c__R22;
        // auto a_g__x_c__R23 = a_g * x_c__R23;
        // auto a_g__x_c__R31 = a_g * x_c__R31;
        // auto a_g__x_c__R32 = a_g * x_c__R32;
        // auto a_g__x_c__R33 = a_g * x_c__R33;
        // auto a_g__y_c__R11 = a_g * y_c__R11;
        // auto a_g__y_c__R12 = a_g * y_c__R12;
        // auto a_g__y_c__R13 = a_g * y_c__R13;
        // auto a_g__y_c__R21 = a_g * y_c__R21;
        // auto a_g__y_c__R22 = a_g * y_c__R22;
        // auto a_g__y_c__R23 = a_g * y_c__R23;
        // auto a_g__y_c__R31 = a_g * y_c__R31;
        // auto a_g__y_c__R32 = a_g * y_c__R32;
        // auto a_g__y_c__R33 = a_g * y_c__R33;
        // auto a_g__z_c__R11 = a_g * z_c__R11;
        // auto a_g__z_c__R12 = a_g * z_c__R12;
        // auto a_g__z_c__R13 = a_g * z_c__R13;
        // auto a_g__z_c__R21 = a_g * z_c__R21;
        // auto a_g__z_c__R22 = a_g * z_c__R22;
        // auto a_g__z_c__R23 = a_g * z_c__R23;
        // auto a_g__z_c__R31 = a_g * z_c__R31;
        // auto a_g__z_c__R32 = a_g * z_c__R32;
        // auto a_g__z_c__R33 = a_g * z_c__R33;
        // auto b_g__a_c__R11 = b_g * a_c__R11;
        // auto b_g__a_c__R12 = b_g * a_c__R12;
        // auto b_g__a_c__R13 = b_g * a_c__R13;
        // auto b_g__a_c__R21 = b_g * a_c__R21;
        // auto b_g__a_c__R22 = b_g * a_c__R22;
        // auto b_g__a_c__R23 = b_g * a_c__R23;
        // auto b_g__a_c__R31 = b_g * a_c__R31;
        // auto b_g__a_c__R32 = b_g * a_c__R32;
        // auto b_g__a_c__R33 = b_g * a_c__R33;
        // auto b_g__b_c__R11 = b_g * b_c__R11;
        // auto b_g__b_c__R12 = b_g * b_c__R12;
        // auto b_g__b_c__R13 = b_g * b_c__R13;
        // auto b_g__b_c__R21 = b_g * b_c__R21;
        // auto b_g__b_c__R22 = b_g * b_c__R22;
        // auto b_g__b_c__R23 = b_g * b_c__R23;
        // auto b_g__b_c__R31 = b_g * b_c__R31;
        // auto b_g__b_c__R32 = b_g * b_c__R32;
        // auto b_g__b_c__R33 = b_g * b_c__R33;
        // auto b_g__c_c__R11 = b_g * c_c__R11;
        // auto b_g__c_c__R12 = b_g * c_c__R12;
        // auto b_g__c_c__R13 = b_g * c_c__R13;
        // auto b_g__c_c__R21 = b_g * c_c__R21;
        // auto b_g__c_c__R22 = b_g * c_c__R22;
        // auto b_g__c_c__R23 = b_g * c_c__R23;
        // auto b_g__c_c__R31 = b_g * c_c__R31;
        // auto b_g__c_c__R32 = b_g * c_c__R32;
        // auto b_g__c_c__R33 = b_g * c_c__R33;
        // auto b_g__x_c__R11 = b_g * x_c__R11;
        // auto b_g__x_c__R12 = b_g * x_c__R12;
        // auto b_g__x_c__R13 = b_g * x_c__R13;
        // auto b_g__x_c__R21 = b_g * x_c__R21;
        // auto b_g__x_c__R22 = b_g * x_c__R22;
        // auto b_g__x_c__R23 = b_g * x_c__R23;
        // auto b_g__x_c__R31 = b_g * x_c__R31;
        // auto b_g__x_c__R32 = b_g * x_c__R32;
        // auto b_g__x_c__R33 = b_g * x_c__R33;
        // auto b_g__y_c__R11 = b_g * y_c__R11;
        // auto b_g__y_c__R12 = b_g * y_c__R12;
        // auto b_g__y_c__R13 = b_g * y_c__R13;
        // auto b_g__y_c__R21 = b_g * y_c__R21;
        // auto b_g__y_c__R22 = b_g * y_c__R22;
        // auto b_g__y_c__R23 = b_g * y_c__R23;
        // auto b_g__y_c__R31 = b_g * y_c__R31;
        // auto b_g__y_c__R32 = b_g * y_c__R32;
        // auto b_g__y_c__R33 = b_g * y_c__R33;
        // auto b_g__z_c__R11 = b_g * z_c__R11;
        // auto b_g__z_c__R12 = b_g * z_c__R12;
        // auto b_g__z_c__R13 = b_g * z_c__R13;
        // auto b_g__z_c__R21 = b_g * z_c__R21;
        // auto b_g__z_c__R22 = b_g * z_c__R22;
        // auto b_g__z_c__R23 = b_g * z_c__R23;
        // auto b_g__z_c__R31 = b_g * z_c__R31;
        // auto b_g__z_c__R32 = b_g * z_c__R32;
        // auto b_g__z_c__R33 = b_g * z_c__R33;
        // auto c_g__a_c__R11 = c_g * a_c__R11;
        // auto c_g__a_c__R12 = c_g * a_c__R12;
        // auto c_g__a_c__R13 = c_g * a_c__R13;
        // auto c_g__a_c__R21 = c_g * a_c__R21;
        // auto c_g__a_c__R22 = c_g * a_c__R22;
        // auto c_g__a_c__R23 = c_g * a_c__R23;
        // auto c_g__a_c__R31 = c_g * a_c__R31;
        // auto c_g__a_c__R32 = c_g * a_c__R32;
        // auto c_g__a_c__R33 = c_g * a_c__R33;
        // auto c_g__b_c__R11 = c_g * b_c__R11;
        // auto c_g__b_c__R12 = c_g * b_c__R12;
        // auto c_g__b_c__R13 = c_g * b_c__R13;
        // auto c_g__b_c__R21 = c_g * b_c__R21;
        // auto c_g__b_c__R22 = c_g * b_c__R22;
        // auto c_g__b_c__R23 = c_g * b_c__R23;
        // auto c_g__b_c__R31 = c_g * b_c__R31;
        // auto c_g__b_c__R32 = c_g * b_c__R32;
        // auto c_g__b_c__R33 = c_g * b_c__R33;
        // auto c_g__c_c__R11 = c_g * c_c__R11;
        // auto c_g__c_c__R12 = c_g * c_c__R12;
        // auto c_g__c_c__R13 = c_g * c_c__R13;
        // auto c_g__c_c__R21 = c_g * c_c__R21;
        // auto c_g__c_c__R22 = c_g * c_c__R22;
        // auto c_g__c_c__R23 = c_g * c_c__R23;
        // auto c_g__c_c__R31 = c_g * c_c__R31;
        // auto c_g__c_c__R32 = c_g * c_c__R32;
        // auto c_g__c_c__R33 = c_g * c_c__R33;
        // auto c_g__x_c__R11 = c_g * x_c__R11;
        // auto c_g__x_c__R12 = c_g * x_c__R12;
        // auto c_g__x_c__R13 = c_g * x_c__R13;
        // auto c_g__x_c__R21 = c_g * x_c__R21;
        // auto c_g__x_c__R22 = c_g * x_c__R22;
        // auto c_g__x_c__R23 = c_g * x_c__R23;
        // auto c_g__x_c__R31 = c_g * x_c__R31;
        // auto c_g__x_c__R32 = c_g * x_c__R32;
        // auto c_g__x_c__R33 = c_g * x_c__R33;
        // auto c_g__y_c__R11 = c_g * y_c__R11;
        // auto c_g__y_c__R12 = c_g * y_c__R12;
        // auto c_g__y_c__R13 = c_g * y_c__R13;
        // auto c_g__y_c__R21 = c_g * y_c__R21;
        // auto c_g__y_c__R22 = c_g * y_c__R22;
        // auto c_g__y_c__R23 = c_g * y_c__R23;
        // auto c_g__y_c__R31 = c_g * y_c__R31;
        // auto c_g__y_c__R32 = c_g * y_c__R32;
        // auto c_g__y_c__R33 = c_g * y_c__R33;
        // auto c_g__z_c__R11 = c_g * z_c__R11;
        // auto c_g__z_c__R12 = c_g * z_c__R12;
        // auto c_g__z_c__R13 = c_g * z_c__R13;
        // auto c_g__z_c__R21 = c_g * z_c__R21;
        // auto c_g__z_c__R22 = c_g * z_c__R22;
        // auto c_g__z_c__R23 = c_g * z_c__R23;
        // auto c_g__z_c__R31 = c_g * z_c__R31;
        // auto c_g__z_c__R32 = c_g * z_c__R32;
        // auto c_g__z_c__R33 = c_g * z_c__R33;
        // auto x_g__a_c__R11 = x_g * a_c__R11;
        // auto x_g__a_c__R12 = x_g * a_c__R12;
        // auto x_g__a_c__R13 = x_g * a_c__R13;
        // auto x_g__a_c__R21 = x_g * a_c__R21;
        // auto x_g__a_c__R22 = x_g * a_c__R22;
        // auto x_g__a_c__R23 = x_g * a_c__R23;
        // auto x_g__a_c__R31 = x_g * a_c__R31;
        // auto x_g__a_c__R32 = x_g * a_c__R32;
        // auto x_g__a_c__R33 = x_g * a_c__R33;
        // auto x_g__b_c__R11 = x_g * b_c__R11;
        // auto x_g__b_c__R12 = x_g * b_c__R12;
        // auto x_g__b_c__R13 = x_g * b_c__R13;
        // auto x_g__b_c__R21 = x_g * b_c__R21;
        // auto x_g__b_c__R22 = x_g * b_c__R22;
        // auto x_g__b_c__R23 = x_g * b_c__R23;
        // auto x_g__b_c__R31 = x_g * b_c__R31;
        // auto x_g__b_c__R32 = x_g * b_c__R32;
        // auto x_g__b_c__R33 = x_g * b_c__R33;
        // auto x_g__c_c__R11 = x_g * c_c__R11;
        // auto x_g__c_c__R12 = x_g * c_c__R12;
        // auto x_g__c_c__R13 = x_g * c_c__R13;
        // auto x_g__c_c__R21 = x_g * c_c__R21;
        // auto x_g__c_c__R22 = x_g * c_c__R22;
        // auto x_g__c_c__R23 = x_g * c_c__R23;
        // auto x_g__c_c__R31 = x_g * c_c__R31;
        // auto x_g__c_c__R32 = x_g * c_c__R32;
        // auto x_g__c_c__R33 = x_g * c_c__R33;
        // auto x_g__x_c__R11 = x_g * x_c__R11;
        // auto x_g__x_c__R12 = x_g * x_c__R12;
        // auto x_g__x_c__R13 = x_g * x_c__R13;
        // auto x_g__x_c__R21 = x_g * x_c__R21;
        // auto x_g__x_c__R22 = x_g * x_c__R22;
        // auto x_g__x_c__R23 = x_g * x_c__R23;
        // auto x_g__x_c__R31 = x_g * x_c__R31;
        // auto x_g__x_c__R32 = x_g * x_c__R32;
        // auto x_g__x_c__R33 = x_g * x_c__R33;
        // auto x_g__y_c__R11 = x_g * y_c__R11;
        // auto x_g__y_c__R12 = x_g * y_c__R12;
        // auto x_g__y_c__R13 = x_g * y_c__R13;
        // auto x_g__y_c__R21 = x_g * y_c__R21;
        // auto x_g__y_c__R22 = x_g * y_c__R22;
        // auto x_g__y_c__R23 = x_g * y_c__R23;
        // auto x_g__y_c__R31 = x_g * y_c__R31;
        // auto x_g__y_c__R32 = x_g * y_c__R32;
        // auto x_g__y_c__R33 = x_g * y_c__R33;
        // auto x_g__z_c__R11 = x_g * z_c__R11;
        // auto x_g__z_c__R12 = x_g * z_c__R12;
        // auto x_g__z_c__R13 = x_g * z_c__R13;
        // auto x_g__z_c__R21 = x_g * z_c__R21;
        // auto x_g__z_c__R22 = x_g * z_c__R22;
        // auto x_g__z_c__R23 = x_g * z_c__R23;
        // auto x_g__z_c__R31 = x_g * z_c__R31;
        // auto x_g__z_c__R32 = x_g * z_c__R32;
        // auto x_g__z_c__R33 = x_g * z_c__R33;
        // auto y_g__a_c__R11 = y_g * a_c__R11;
        // auto y_g__a_c__R12 = y_g * a_c__R12;
        // auto y_g__a_c__R13 = y_g * a_c__R13;
        // auto y_g__a_c__R21 = y_g * a_c__R21;
        // auto y_g__a_c__R22 = y_g * a_c__R22;
        // auto y_g__a_c__R23 = y_g * a_c__R23;
        // auto y_g__a_c__R31 = y_g * a_c__R31;
        // auto y_g__a_c__R32 = y_g * a_c__R32;
        // auto y_g__a_c__R33 = y_g * a_c__R33;
        // auto y_g__b_c__R11 = y_g * b_c__R11;
        // auto y_g__b_c__R12 = y_g * b_c__R12;
        // auto y_g__b_c__R13 = y_g * b_c__R13;
        // auto y_g__b_c__R21 = y_g * b_c__R21;
        // auto y_g__b_c__R22 = y_g * b_c__R22;
        // auto y_g__b_c__R23 = y_g * b_c__R23;
        // auto y_g__b_c__R31 = y_g * b_c__R31;
        // auto y_g__b_c__R32 = y_g * b_c__R32;
        // auto y_g__b_c__R33 = y_g * b_c__R33;
        // auto y_g__c_c__R11 = y_g * c_c__R11;
        // auto y_g__c_c__R12 = y_g * c_c__R12;
        // auto y_g__c_c__R13 = y_g * c_c__R13;
        // auto y_g__c_c__R21 = y_g * c_c__R21;
        // auto y_g__c_c__R22 = y_g * c_c__R22;
        // auto y_g__c_c__R23 = y_g * c_c__R23;
        // auto y_g__c_c__R31 = y_g * c_c__R31;
        // auto y_g__c_c__R32 = y_g * c_c__R32;
        // auto y_g__c_c__R33 = y_g * c_c__R33;
        // auto y_g__x_c__R11 = y_g * x_c__R11;
        // auto y_g__x_c__R12 = y_g * x_c__R12;
        // auto y_g__x_c__R13 = y_g * x_c__R13;
        // auto y_g__x_c__R21 = y_g * x_c__R21;
        // auto y_g__x_c__R22 = y_g * x_c__R22;
        // auto y_g__x_c__R23 = y_g * x_c__R23;
        // auto y_g__x_c__R31 = y_g * x_c__R31;
        // auto y_g__x_c__R32 = y_g * x_c__R32;
        // auto y_g__x_c__R33 = y_g * x_c__R33;
        // auto y_g__y_c__R11 = y_g * y_c__R11;
        // auto y_g__y_c__R12 = y_g * y_c__R12;
        // auto y_g__y_c__R13 = y_g * y_c__R13;
        // auto y_g__y_c__R21 = y_g * y_c__R21;
        // auto y_g__y_c__R22 = y_g * y_c__R22;
        // auto y_g__y_c__R23 = y_g * y_c__R23;
        // auto y_g__y_c__R31 = y_g * y_c__R31;
        // auto y_g__y_c__R32 = y_g * y_c__R32;
        // auto y_g__y_c__R33 = y_g * y_c__R33;
        // auto y_g__z_c__R11 = y_g * z_c__R11;
        // auto y_g__z_c__R12 = y_g * z_c__R12;
        // auto y_g__z_c__R13 = y_g * z_c__R13;
        // auto y_g__z_c__R21 = y_g * z_c__R21;
        // auto y_g__z_c__R22 = y_g * z_c__R22;
        // auto y_g__z_c__R23 = y_g * z_c__R23;
        // auto y_g__z_c__R31 = y_g * z_c__R31;
        // auto y_g__z_c__R32 = y_g * z_c__R32;
        // auto y_g__z_c__R33 = y_g * z_c__R33;
        // auto z_g__a_c__R11 = z_g * a_c__R11;
        // auto z_g__a_c__R12 = z_g * a_c__R12;
        // auto z_g__a_c__R13 = z_g * a_c__R13;
        // auto z_g__a_c__R21 = z_g * a_c__R21;
        // auto z_g__a_c__R22 = z_g * a_c__R22;
        // auto z_g__a_c__R23 = z_g * a_c__R23;
        // auto z_g__a_c__R31 = z_g * a_c__R31;
        // auto z_g__a_c__R32 = z_g * a_c__R32;
        // auto z_g__a_c__R33 = z_g * a_c__R33;
        // auto z_g__b_c__R11 = z_g * b_c__R11;
        // auto z_g__b_c__R12 = z_g * b_c__R12;
        // auto z_g__b_c__R13 = z_g * b_c__R13;
        // auto z_g__b_c__R21 = z_g * b_c__R21;
        // auto z_g__b_c__R22 = z_g * b_c__R22;
        // auto z_g__b_c__R23 = z_g * b_c__R23;
        // auto z_g__b_c__R31 = z_g * b_c__R31;
        // auto z_g__b_c__R32 = z_g * b_c__R32;
        // auto z_g__b_c__R33 = z_g * b_c__R33;
        // auto z_g__c_c__R11 = z_g * c_c__R11;
        // auto z_g__c_c__R12 = z_g * c_c__R12;
        // auto z_g__c_c__R13 = z_g * c_c__R13;
        // auto z_g__c_c__R21 = z_g * c_c__R21;
        // auto z_g__c_c__R22 = z_g * c_c__R22;
        // auto z_g__c_c__R23 = z_g * c_c__R23;
        // auto z_g__c_c__R31 = z_g * c_c__R31;
        // auto z_g__c_c__R32 = z_g * c_c__R32;
        // auto z_g__c_c__R33 = z_g * c_c__R33;
        // auto z_g__x_c__R11 = z_g * x_c__R11;
        // auto z_g__x_c__R12 = z_g * x_c__R12;
        // auto z_g__x_c__R13 = z_g * x_c__R13;
        // auto z_g__x_c__R21 = z_g * x_c__R21;
        // auto z_g__x_c__R22 = z_g * x_c__R22;
        // auto z_g__x_c__R23 = z_g * x_c__R23;
        // auto z_g__x_c__R31 = z_g * x_c__R31;
        // auto z_g__x_c__R32 = z_g * x_c__R32;
        // auto z_g__x_c__R33 = z_g * x_c__R33;
        // auto z_g__y_c__R11 = z_g * y_c__R11;
        // auto z_g__y_c__R12 = z_g * y_c__R12;
        // auto z_g__y_c__R13 = z_g * y_c__R13;
        // auto z_g__y_c__R21 = z_g * y_c__R21;
        // auto z_g__y_c__R22 = z_g * y_c__R22;
        // auto z_g__y_c__R23 = z_g * y_c__R23;
        // auto z_g__y_c__R31 = z_g * y_c__R31;
        // auto z_g__y_c__R32 = z_g * y_c__R32;
        // auto z_g__y_c__R33 = z_g * y_c__R33;
        // auto z_g__z_c__R11 = z_g * z_c__R11;
        // auto z_g__z_c__R12 = z_g * z_c__R12;
        // auto z_g__z_c__R13 = z_g * z_c__R13;
        // auto z_g__z_c__R21 = z_g * z_c__R21;
        // auto z_g__z_c__R22 = z_g * z_c__R22;
        // auto z_g__z_c__R23 = z_g * z_c__R23;
        // auto z_g__z_c__R31 = z_g * z_c__R31;
        // auto z_g__z_c__R32 = z_g * z_c__R32;
        // auto z_g__z_c__R33 = z_g * z_c__R33;
        // auto tx__a_c__R11 = tx * a_c__R11;
        // auto tx__a_c__R12 = tx * a_c__R12;
        // auto tx__a_c__R13 = tx * a_c__R13;
        // auto tx__a_c__R21 = tx * a_c__R21;
        // auto tx__a_c__R22 = tx * a_c__R22;
        // auto tx__a_c__R23 = tx * a_c__R23;
        // auto tx__a_c__R31 = tx * a_c__R31;
        // auto tx__a_c__R32 = tx * a_c__R32;
        // auto tx__a_c__R33 = tx * a_c__R33;
        // auto tx__b_c__R11 = tx * b_c__R11;
        // auto tx__b_c__R12 = tx * b_c__R12;
        // auto tx__b_c__R13 = tx * b_c__R13;
        // auto tx__b_c__R21 = tx * b_c__R21;
        // auto tx__b_c__R22 = tx * b_c__R22;
        // auto tx__b_c__R23 = tx * b_c__R23;
        // auto tx__b_c__R31 = tx * b_c__R31;
        // auto tx__b_c__R32 = tx * b_c__R32;
        // auto tx__b_c__R33 = tx * b_c__R33;
        // auto tx__c_c__R11 = tx * c_c__R11;
        // auto tx__c_c__R12 = tx * c_c__R12;
        // auto tx__c_c__R13 = tx * c_c__R13;
        // auto tx__c_c__R21 = tx * c_c__R21;
        // auto tx__c_c__R22 = tx * c_c__R22;
        // auto tx__c_c__R23 = tx * c_c__R23;
        // auto tx__c_c__R31 = tx * c_c__R31;
        // auto tx__c_c__R32 = tx * c_c__R32;
        // auto tx__c_c__R33 = tx * c_c__R33;
        // auto tx__x_c__R11 = tx * x_c__R11;
        // auto tx__x_c__R12 = tx * x_c__R12;
        // auto tx__x_c__R13 = tx * x_c__R13;
        // auto tx__x_c__R21 = tx * x_c__R21;
        // auto tx__x_c__R22 = tx * x_c__R22;
        // auto tx__x_c__R23 = tx * x_c__R23;
        // auto tx__x_c__R31 = tx * x_c__R31;
        // auto tx__x_c__R32 = tx * x_c__R32;
        // auto tx__x_c__R33 = tx * x_c__R33;
        // auto tx__y_c__R11 = tx * y_c__R11;
        // auto tx__y_c__R12 = tx * y_c__R12;
        // auto tx__y_c__R13 = tx * y_c__R13;
        // auto tx__y_c__R21 = tx * y_c__R21;
        // auto tx__y_c__R22 = tx * y_c__R22;
        // auto tx__y_c__R23 = tx * y_c__R23;
        // auto tx__y_c__R31 = tx * y_c__R31;
        // auto tx__y_c__R32 = tx * y_c__R32;
        // auto tx__y_c__R33 = tx * y_c__R33;
        // auto tx__z_c__R11 = tx * z_c__R11;
        // auto tx__z_c__R12 = tx * z_c__R12;
        // auto tx__z_c__R13 = tx * z_c__R13;
        // auto tx__z_c__R21 = tx * z_c__R21;
        // auto tx__z_c__R22 = tx * z_c__R22;
        // auto tx__z_c__R23 = tx * z_c__R23;
        // auto tx__z_c__R31 = tx * z_c__R31;
        // auto tx__z_c__R32 = tx * z_c__R32;
        // auto tx__z_c__R33 = tx * z_c__R33;
        // auto ty__a_c__R11 = ty * a_c__R11;
        // auto ty__a_c__R12 = ty * a_c__R12;
        // auto ty__a_c__R13 = ty * a_c__R13;
        // auto ty__a_c__R21 = ty * a_c__R21;
        // auto ty__a_c__R22 = ty * a_c__R22;
        // auto ty__a_c__R23 = ty * a_c__R23;
        // auto ty__a_c__R31 = ty * a_c__R31;
        // auto ty__a_c__R32 = ty * a_c__R32;
        // auto ty__a_c__R33 = ty * a_c__R33;
        // auto ty__b_c__R11 = ty * b_c__R11;
        // auto ty__b_c__R12 = ty * b_c__R12;
        // auto ty__b_c__R13 = ty * b_c__R13;
        // auto ty__b_c__R21 = ty * b_c__R21;
        // auto ty__b_c__R22 = ty * b_c__R22;
        // auto ty__b_c__R23 = ty * b_c__R23;
        // auto ty__b_c__R31 = ty * b_c__R31;
        // auto ty__b_c__R32 = ty * b_c__R32;
        // auto ty__b_c__R33 = ty * b_c__R33;
        // auto ty__c_c__R11 = ty * c_c__R11;
        // auto ty__c_c__R12 = ty * c_c__R12;
        // auto ty__c_c__R13 = ty * c_c__R13;
        // auto ty__c_c__R21 = ty * c_c__R21;
        // auto ty__c_c__R22 = ty * c_c__R22;
        // auto ty__c_c__R23 = ty * c_c__R23;
        // auto ty__c_c__R31 = ty * c_c__R31;
        // auto ty__c_c__R32 = ty * c_c__R32;
        // auto ty__c_c__R33 = ty * c_c__R33;
        // auto ty__x_c__R11 = ty * x_c__R11;
        // auto ty__x_c__R12 = ty * x_c__R12;
        // auto ty__x_c__R13 = ty * x_c__R13;
        // auto ty__x_c__R21 = ty * x_c__R21;
        // auto ty__x_c__R22 = ty * x_c__R22;
        // auto ty__x_c__R23 = ty * x_c__R23;
        // auto ty__x_c__R31 = ty * x_c__R31;
        // auto ty__x_c__R32 = ty * x_c__R32;
        // auto ty__x_c__R33 = ty * x_c__R33;
        // auto ty__y_c__R11 = ty * y_c__R11;
        // auto ty__y_c__R12 = ty * y_c__R12;
        // auto ty__y_c__R13 = ty * y_c__R13;
        // auto ty__y_c__R21 = ty * y_c__R21;
        // auto ty__y_c__R22 = ty * y_c__R22;
        // auto ty__y_c__R23 = ty * y_c__R23;
        // auto ty__y_c__R31 = ty * y_c__R31;
        // auto ty__y_c__R32 = ty * y_c__R32;
        // auto ty__y_c__R33 = ty * y_c__R33;
        // auto ty__z_c__R11 = ty * z_c__R11;
        // auto ty__z_c__R12 = ty * z_c__R12;
        // auto ty__z_c__R13 = ty * z_c__R13;
        // auto ty__z_c__R21 = ty * z_c__R21;
        // auto ty__z_c__R22 = ty * z_c__R22;
        // auto ty__z_c__R23 = ty * z_c__R23;
        // auto ty__z_c__R31 = ty * z_c__R31;
        // auto ty__z_c__R32 = ty * z_c__R32;
        // auto ty__z_c__R33 = ty * z_c__R33;
        // auto tz__a_c__R11 = tz * a_c__R11;
        // auto tz__a_c__R12 = tz * a_c__R12;
        // auto tz__a_c__R13 = tz * a_c__R13;
        // auto tz__a_c__R21 = tz * a_c__R21;
        // auto tz__a_c__R22 = tz * a_c__R22;
        // auto tz__a_c__R23 = tz * a_c__R23;
        // auto tz__a_c__R31 = tz * a_c__R31;
        // auto tz__a_c__R32 = tz * a_c__R32;
        // auto tz__a_c__R33 = tz * a_c__R33;
        // auto tz__b_c__R11 = tz * b_c__R11;
        // auto tz__b_c__R12 = tz * b_c__R12;
        // auto tz__b_c__R13 = tz * b_c__R13;
        // auto tz__b_c__R21 = tz * b_c__R21;
        // auto tz__b_c__R22 = tz * b_c__R22;
        // auto tz__b_c__R23 = tz * b_c__R23;
        // auto tz__b_c__R31 = tz * b_c__R31;
        // auto tz__b_c__R32 = tz * b_c__R32;
        // auto tz__b_c__R33 = tz * b_c__R33;
        // auto tz__c_c__R11 = tz * c_c__R11;
        // auto tz__c_c__R12 = tz * c_c__R12;
        // auto tz__c_c__R13 = tz * c_c__R13;
        // auto tz__c_c__R21 = tz * c_c__R21;
        // auto tz__c_c__R22 = tz * c_c__R22;
        // auto tz__c_c__R23 = tz * c_c__R23;
        // auto tz__c_c__R31 = tz * c_c__R31;
        // auto tz__c_c__R32 = tz * c_c__R32;
        // auto tz__c_c__R33 = tz * c_c__R33;
        // auto tz__x_c__R11 = tz * x_c__R11;
        // auto tz__x_c__R12 = tz * x_c__R12;
        // auto tz__x_c__R13 = tz * x_c__R13;
        // auto tz__x_c__R21 = tz * x_c__R21;
        // auto tz__x_c__R22 = tz * x_c__R22;
        // auto tz__x_c__R23 = tz * x_c__R23;
        // auto tz__x_c__R31 = tz * x_c__R31;
        // auto tz__x_c__R32 = tz * x_c__R32;
        // auto tz__x_c__R33 = tz * x_c__R33;
        // auto tz__y_c__R11 = tz * y_c__R11;
        // auto tz__y_c__R12 = tz * y_c__R12;
        // auto tz__y_c__R13 = tz * y_c__R13;
        // auto tz__y_c__R21 = tz * y_c__R21;
        // auto tz__y_c__R22 = tz * y_c__R22;
        // auto tz__y_c__R23 = tz * y_c__R23;
        // auto tz__y_c__R31 = tz * y_c__R31;
        // auto tz__y_c__R32 = tz * y_c__R32;
        // auto tz__y_c__R33 = tz * y_c__R33;
        // auto tz__z_c__R11 = tz * z_c__R11;
        // auto tz__z_c__R12 = tz * z_c__R12;
        // auto tz__z_c__R13 = tz * z_c__R13;
        // auto tz__z_c__R21 = tz * z_c__R21;
        // auto tz__z_c__R22 = tz * z_c__R22;
        // auto tz__z_c__R23 = tz * z_c__R23;
        // auto tz__z_c__R31 = tz * z_c__R31;
        // auto tz__z_c__R32 = tz * z_c__R32;
        // auto tz__z_c__R33 = tz * z_c__R33;
        // auto bx__a_c__R11 = bx * a_c__R11;
        // auto bx__a_c__R12 = bx * a_c__R12;
        // auto bx__a_c__R13 = bx * a_c__R13;
        // auto bx__a_c__R21 = bx * a_c__R21;
        // auto bx__a_c__R22 = bx * a_c__R22;
        // auto bx__a_c__R23 = bx * a_c__R23;
        // auto bx__a_c__R31 = bx * a_c__R31;
        // auto bx__a_c__R32 = bx * a_c__R32;
        // auto bx__a_c__R33 = bx * a_c__R33;
        // auto bx__b_c__R11 = bx * b_c__R11;
        // auto bx__b_c__R12 = bx * b_c__R12;
        // auto bx__b_c__R13 = bx * b_c__R13;
        // auto bx__b_c__R21 = bx * b_c__R21;
        // auto bx__b_c__R22 = bx * b_c__R22;
        // auto bx__b_c__R23 = bx * b_c__R23;
        // auto bx__b_c__R31 = bx * b_c__R31;
        // auto bx__b_c__R32 = bx * b_c__R32;
        // auto bx__b_c__R33 = bx * b_c__R33;
        // auto bx__c_c__R11 = bx * c_c__R11;
        // auto bx__c_c__R12 = bx * c_c__R12;
        // auto bx__c_c__R13 = bx * c_c__R13;
        // auto bx__c_c__R21 = bx * c_c__R21;
        // auto bx__c_c__R22 = bx * c_c__R22;
        // auto bx__c_c__R23 = bx * c_c__R23;
        // auto bx__c_c__R31 = bx * c_c__R31;
        // auto bx__c_c__R32 = bx * c_c__R32;
        // auto bx__c_c__R33 = bx * c_c__R33;
        // auto bx__x_c__R11 = bx * x_c__R11;
        // auto bx__x_c__R12 = bx * x_c__R12;
        // auto bx__x_c__R13 = bx * x_c__R13;
        // auto bx__x_c__R21 = bx * x_c__R21;
        // auto bx__x_c__R22 = bx * x_c__R22;
        // auto bx__x_c__R23 = bx * x_c__R23;
        // auto bx__x_c__R31 = bx * x_c__R31;
        // auto bx__x_c__R32 = bx * x_c__R32;
        // auto bx__x_c__R33 = bx * x_c__R33;
        // auto bx__y_c__R11 = bx * y_c__R11;
        // auto bx__y_c__R12 = bx * y_c__R12;
        // auto bx__y_c__R13 = bx * y_c__R13;
        // auto bx__y_c__R21 = bx * y_c__R21;
        // auto bx__y_c__R22 = bx * y_c__R22;
        // auto bx__y_c__R23 = bx * y_c__R23;
        // auto bx__y_c__R31 = bx * y_c__R31;
        // auto bx__y_c__R32 = bx * y_c__R32;
        // auto bx__y_c__R33 = bx * y_c__R33;
        // auto bx__z_c__R11 = bx * z_c__R11;
        // auto bx__z_c__R12 = bx * z_c__R12;
        // auto bx__z_c__R13 = bx * z_c__R13;
        // auto bx__z_c__R21 = bx * z_c__R21;
        // auto bx__z_c__R22 = bx * z_c__R22;
        // auto bx__z_c__R23 = bx * z_c__R23;
        // auto bx__z_c__R31 = bx * z_c__R31;
        // auto bx__z_c__R32 = bx * z_c__R32;
        // auto bx__z_c__R33 = bx * z_c__R33;
        // auto by__a_c__R11 = by * a_c__R11;
        // auto by__a_c__R12 = by * a_c__R12;
        // auto by__a_c__R13 = by * a_c__R13;
        // auto by__a_c__R21 = by * a_c__R21;
        // auto by__a_c__R22 = by * a_c__R22;
        // auto by__a_c__R23 = by * a_c__R23;
        // auto by__a_c__R31 = by * a_c__R31;
        // auto by__a_c__R32 = by * a_c__R32;
        // auto by__a_c__R33 = by * a_c__R33;
        // auto by__b_c__R11 = by * b_c__R11;
        // auto by__b_c__R12 = by * b_c__R12;
        // auto by__b_c__R13 = by * b_c__R13;
        // auto by__b_c__R21 = by * b_c__R21;
        // auto by__b_c__R22 = by * b_c__R22;
        // auto by__b_c__R23 = by * b_c__R23;
        // auto by__b_c__R31 = by * b_c__R31;
        // auto by__b_c__R32 = by * b_c__R32;
        // auto by__b_c__R33 = by * b_c__R33;
        // auto by__c_c__R11 = by * c_c__R11;
        // auto by__c_c__R12 = by * c_c__R12;
        // auto by__c_c__R13 = by * c_c__R13;
        // auto by__c_c__R21 = by * c_c__R21;
        // auto by__c_c__R22 = by * c_c__R22;
        // auto by__c_c__R23 = by * c_c__R23;
        // auto by__c_c__R31 = by * c_c__R31;
        // auto by__c_c__R32 = by * c_c__R32;
        // auto by__c_c__R33 = by * c_c__R33;
        // auto by__x_c__R11 = by * x_c__R11;
        // auto by__x_c__R12 = by * x_c__R12;
        // auto by__x_c__R13 = by * x_c__R13;
        // auto by__x_c__R21 = by * x_c__R21;
        // auto by__x_c__R22 = by * x_c__R22;
        // auto by__x_c__R23 = by * x_c__R23;
        // auto by__x_c__R31 = by * x_c__R31;
        // auto by__x_c__R32 = by * x_c__R32;
        // auto by__x_c__R33 = by * x_c__R33;
        // auto by__y_c__R11 = by * y_c__R11;
        // auto by__y_c__R12 = by * y_c__R12;
        // auto by__y_c__R13 = by * y_c__R13;
        // auto by__y_c__R21 = by * y_c__R21;
        // auto by__y_c__R22 = by * y_c__R22;
        // auto by__y_c__R23 = by * y_c__R23;
        // auto by__y_c__R31 = by * y_c__R31;
        // auto by__y_c__R32 = by * y_c__R32;
        // auto by__y_c__R33 = by * y_c__R33;
        // auto by__z_c__R11 = by * z_c__R11;
        // auto by__z_c__R12 = by * z_c__R12;
        // auto by__z_c__R13 = by * z_c__R13;
        // auto by__z_c__R21 = by * z_c__R21;
        // auto by__z_c__R22 = by * z_c__R22;
        // auto by__z_c__R23 = by * z_c__R23;
        // auto by__z_c__R31 = by * z_c__R31;
        // auto by__z_c__R32 = by * z_c__R32;
        // auto by__z_c__R33 = by * z_c__R33;
        // auto bz__a_c__R11 = bz * a_c__R11;
        // auto bz__a_c__R12 = bz * a_c__R12;
        // auto bz__a_c__R13 = bz * a_c__R13;
        // auto bz__a_c__R21 = bz * a_c__R21;
        // auto bz__a_c__R22 = bz * a_c__R22;
        // auto bz__a_c__R23 = bz * a_c__R23;
        // auto bz__a_c__R31 = bz * a_c__R31;
        // auto bz__a_c__R32 = bz * a_c__R32;
        // auto bz__a_c__R33 = bz * a_c__R33;
        // auto bz__b_c__R11 = bz * b_c__R11;
        // auto bz__b_c__R12 = bz * b_c__R12;
        // auto bz__b_c__R13 = bz * b_c__R13;
        // auto bz__b_c__R21 = bz * b_c__R21;
        // auto bz__b_c__R22 = bz * b_c__R22;
        // auto bz__b_c__R23 = bz * b_c__R23;
        // auto bz__b_c__R31 = bz * b_c__R31;
        // auto bz__b_c__R32 = bz * b_c__R32;
        // auto bz__b_c__R33 = bz * b_c__R33;
        // auto bz__c_c__R11 = bz * c_c__R11;
        // auto bz__c_c__R12 = bz * c_c__R12;
        // auto bz__c_c__R13 = bz * c_c__R13;
        // auto bz__c_c__R21 = bz * c_c__R21;
        // auto bz__c_c__R22 = bz * c_c__R22;
        // auto bz__c_c__R23 = bz * c_c__R23;
        // auto bz__c_c__R31 = bz * c_c__R31;
        // auto bz__c_c__R32 = bz * c_c__R32;
        // auto bz__c_c__R33 = bz * c_c__R33;
        // auto bz__x_c__R11 = bz * x_c__R11;
        // auto bz__x_c__R12 = bz * x_c__R12;
        // auto bz__x_c__R13 = bz * x_c__R13;
        // auto bz__x_c__R21 = bz * x_c__R21;
        // auto bz__x_c__R22 = bz * x_c__R22;
        // auto bz__x_c__R23 = bz * x_c__R23;
        // auto bz__x_c__R31 = bz * x_c__R31;
        // auto bz__x_c__R32 = bz * x_c__R32;
        // auto bz__x_c__R33 = bz * x_c__R33;
        // auto bz__y_c__R11 = bz * y_c__R11;
        // auto bz__y_c__R12 = bz * y_c__R12;
        // auto bz__y_c__R13 = bz * y_c__R13;
        // auto bz__y_c__R21 = bz * y_c__R21;
        // auto bz__y_c__R22 = bz * y_c__R22;
        // auto bz__y_c__R23 = bz * y_c__R23;
        // auto bz__y_c__R31 = bz * y_c__R31;
        // auto bz__y_c__R32 = bz * y_c__R32;
        // auto bz__y_c__R33 = bz * y_c__R33;
        // auto bz__z_c__R11 = bz * z_c__R11;
        // auto bz__z_c__R12 = bz * z_c__R12;
        // auto bz__z_c__R13 = bz * z_c__R13;
        // auto bz__z_c__R21 = bz * z_c__R21;
        // auto bz__z_c__R22 = bz * z_c__R22;
        // auto bz__z_c__R23 = bz * z_c__R23;
        // auto bz__z_c__R31 = bz * z_c__R31;
        // auto bz__z_c__R32 = bz * z_c__R32;
        // auto bz__z_c__R33 = bz * z_c__R33;
        //
        // auto tx__a_g__R11 = tx * a_g__R11;
        // auto tx__a_g__R12 = tx * a_g__R12;
        // auto tx__a_g__R13 = tx * a_g__R13;
        // auto tx__a_g__R21 = tx * a_g__R21;
        // auto tx__a_g__R22 = tx * a_g__R22;
        // auto tx__a_g__R23 = tx * a_g__R23;
        // auto tx__a_g__R31 = tx * a_g__R31;
        // auto tx__a_g__R32 = tx * a_g__R32;
        // auto tx__a_g__R33 = tx * a_g__R33;
        // auto tx__b_g__R11 = tx * b_g__R11;
        // auto tx__b_g__R12 = tx * b_g__R12;
        // auto tx__b_g__R13 = tx * b_g__R13;
        // auto tx__b_g__R21 = tx * b_g__R21;
        // auto tx__b_g__R22 = tx * b_g__R22;
        // auto tx__b_g__R23 = tx * b_g__R23;
        // auto tx__b_g__R31 = tx * b_g__R31;
        // auto tx__b_g__R32 = tx * b_g__R32;
        // auto tx__b_g__R33 = tx * b_g__R33;
        // auto tx__c_g__R11 = tx * c_g__R11;
        // auto tx__c_g__R12 = tx * c_g__R12;
        // auto tx__c_g__R13 = tx * c_g__R13;
        // auto tx__c_g__R21 = tx * c_g__R21;
        // auto tx__c_g__R22 = tx * c_g__R22;
        // auto tx__c_g__R23 = tx * c_g__R23;
        // auto tx__c_g__R31 = tx * c_g__R31;
        // auto tx__c_g__R32 = tx * c_g__R32;
        // auto tx__c_g__R33 = tx * c_g__R33;
        // auto tx__x_g__R11 = tx * x_g__R11;
        // auto tx__x_g__R12 = tx * x_g__R12;
        // auto tx__x_g__R13 = tx * x_g__R13;
        // auto tx__x_g__R21 = tx * x_g__R21;
        // auto tx__x_g__R22 = tx * x_g__R22;
        // auto tx__x_g__R23 = tx * x_g__R23;
        // auto tx__x_g__R31 = tx * x_g__R31;
        // auto tx__x_g__R32 = tx * x_g__R32;
        // auto tx__x_g__R33 = tx * x_g__R33;
        // auto tx__y_g__R11 = tx * y_g__R11;
        // auto tx__y_g__R12 = tx * y_g__R12;
        // auto tx__y_g__R13 = tx * y_g__R13;
        // auto tx__y_g__R21 = tx * y_g__R21;
        // auto tx__y_g__R22 = tx * y_g__R22;
        // auto tx__y_g__R23 = tx * y_g__R23;
        // auto tx__y_g__R31 = tx * y_g__R31;
        // auto tx__y_g__R32 = tx * y_g__R32;
        // auto tx__y_g__R33 = tx * y_g__R33;
        // auto tx__z_g__R11 = tx * z_g__R11;
        // auto tx__z_g__R12 = tx * z_g__R12;
        // auto tx__z_g__R13 = tx * z_g__R13;
        // auto tx__z_g__R21 = tx * z_g__R21;
        // auto tx__z_g__R22 = tx * z_g__R22;
        // auto tx__z_g__R23 = tx * z_g__R23;
        // auto tx__z_g__R31 = tx * z_g__R31;
        // auto tx__z_g__R32 = tx * z_g__R32;
        // auto tx__z_g__R33 = tx * z_g__R33;
        // auto ty__a_g__R11 = ty * a_g__R11;
        // auto ty__a_g__R12 = ty * a_g__R12;
        // auto ty__a_g__R13 = ty * a_g__R13;
        // auto ty__a_g__R21 = ty * a_g__R21;
        // auto ty__a_g__R22 = ty * a_g__R22;
        // auto ty__a_g__R23 = ty * a_g__R23;
        // auto ty__a_g__R31 = ty * a_g__R31;
        // auto ty__a_g__R32 = ty * a_g__R32;
        // auto ty__a_g__R33 = ty * a_g__R33;
        // auto ty__b_g__R11 = ty * b_g__R11;
        // auto ty__b_g__R12 = ty * b_g__R12;
        // auto ty__b_g__R13 = ty * b_g__R13;
        // auto ty__b_g__R21 = ty * b_g__R21;
        // auto ty__b_g__R22 = ty * b_g__R22;
        // auto ty__b_g__R23 = ty * b_g__R23;
        // auto ty__b_g__R31 = ty * b_g__R31;
        // auto ty__b_g__R32 = ty * b_g__R32;
        // auto ty__b_g__R33 = ty * b_g__R33;
        // auto ty__c_g__R11 = ty * c_g__R11;
        // auto ty__c_g__R12 = ty * c_g__R12;
        // auto ty__c_g__R13 = ty * c_g__R13;
        // auto ty__c_g__R21 = ty * c_g__R21;
        // auto ty__c_g__R22 = ty * c_g__R22;
        // auto ty__c_g__R23 = ty * c_g__R23;
        // auto ty__c_g__R31 = ty * c_g__R31;
        // auto ty__c_g__R32 = ty * c_g__R32;
        // auto ty__c_g__R33 = ty * c_g__R33;
        // auto ty__x_g__R11 = ty * x_g__R11;
        // auto ty__x_g__R12 = ty * x_g__R12;
        // auto ty__x_g__R13 = ty * x_g__R13;
        // auto ty__x_g__R21 = ty * x_g__R21;
        // auto ty__x_g__R22 = ty * x_g__R22;
        // auto ty__x_g__R23 = ty * x_g__R23;
        // auto ty__x_g__R31 = ty * x_g__R31;
        // auto ty__x_g__R32 = ty * x_g__R32;
        // auto ty__x_g__R33 = ty * x_g__R33;
        // auto ty__y_g__R11 = ty * y_g__R11;
        // auto ty__y_g__R12 = ty * y_g__R12;
        // auto ty__y_g__R13 = ty * y_g__R13;
        // auto ty__y_g__R21 = ty * y_g__R21;
        // auto ty__y_g__R22 = ty * y_g__R22;
        // auto ty__y_g__R23 = ty * y_g__R23;
        // auto ty__y_g__R31 = ty * y_g__R31;
        // auto ty__y_g__R32 = ty * y_g__R32;
        // auto ty__y_g__R33 = ty * y_g__R33;
        // auto ty__z_g__R11 = ty * z_g__R11;
        // auto ty__z_g__R12 = ty * z_g__R12;
        // auto ty__z_g__R13 = ty * z_g__R13;
        // auto ty__z_g__R21 = ty * z_g__R21;
        // auto ty__z_g__R22 = ty * z_g__R22;
        // auto ty__z_g__R23 = ty * z_g__R23;
        // auto ty__z_g__R31 = ty * z_g__R31;
        // auto ty__z_g__R32 = ty * z_g__R32;
        // auto ty__z_g__R33 = ty * z_g__R33;
        // auto tz__a_g__R11 = tz * a_g__R11;
        // auto tz__a_g__R12 = tz * a_g__R12;
        // auto tz__a_g__R13 = tz * a_g__R13;
        // auto tz__a_g__R21 = tz * a_g__R21;
        // auto tz__a_g__R22 = tz * a_g__R22;
        // auto tz__a_g__R23 = tz * a_g__R23;
        // auto tz__a_g__R31 = tz * a_g__R31;
        // auto tz__a_g__R32 = tz * a_g__R32;
        // auto tz__a_g__R33 = tz * a_g__R33;
        // auto tz__b_g__R11 = tz * b_g__R11;
        // auto tz__b_g__R12 = tz * b_g__R12;
        // auto tz__b_g__R13 = tz * b_g__R13;
        // auto tz__b_g__R21 = tz * b_g__R21;
        // auto tz__b_g__R22 = tz * b_g__R22;
        // auto tz__b_g__R23 = tz * b_g__R23;
        // auto tz__b_g__R31 = tz * b_g__R31;
        // auto tz__b_g__R32 = tz * b_g__R32;
        // auto tz__b_g__R33 = tz * b_g__R33;
        // auto tz__c_g__R11 = tz * c_g__R11;
        // auto tz__c_g__R12 = tz * c_g__R12;
        // auto tz__c_g__R13 = tz * c_g__R13;
        // auto tz__c_g__R21 = tz * c_g__R21;
        // auto tz__c_g__R22 = tz * c_g__R22;
        // auto tz__c_g__R23 = tz * c_g__R23;
        // auto tz__c_g__R31 = tz * c_g__R31;
        // auto tz__c_g__R32 = tz * c_g__R32;
        // auto tz__c_g__R33 = tz * c_g__R33;
        // auto tz__x_g__R11 = tz * x_g__R11;
        // auto tz__x_g__R12 = tz * x_g__R12;
        // auto tz__x_g__R13 = tz * x_g__R13;
        // auto tz__x_g__R21 = tz * x_g__R21;
        // auto tz__x_g__R22 = tz * x_g__R22;
        // auto tz__x_g__R23 = tz * x_g__R23;
        // auto tz__x_g__R31 = tz * x_g__R31;
        // auto tz__x_g__R32 = tz * x_g__R32;
        // auto tz__x_g__R33 = tz * x_g__R33;
        // auto tz__y_g__R11 = tz * y_g__R11;
        // auto tz__y_g__R12 = tz * y_g__R12;
        // auto tz__y_g__R13 = tz * y_g__R13;
        // auto tz__y_g__R21 = tz * y_g__R21;
        // auto tz__y_g__R22 = tz * y_g__R22;
        // auto tz__y_g__R23 = tz * y_g__R23;
        // auto tz__y_g__R31 = tz * y_g__R31;
        // auto tz__y_g__R32 = tz * y_g__R32;
        // auto tz__y_g__R33 = tz * y_g__R33;
        // auto tz__z_g__R11 = tz * z_g__R11;
        // auto tz__z_g__R12 = tz * z_g__R12;
        // auto tz__z_g__R13 = tz * z_g__R13;
        // auto tz__z_g__R21 = tz * z_g__R21;
        // auto tz__z_g__R22 = tz * z_g__R22;
        // auto tz__z_g__R23 = tz * z_g__R23;
        // auto tz__z_g__R31 = tz * z_g__R31;
        // auto tz__z_g__R32 = tz * z_g__R32;
        // auto tz__z_g__R33 = tz * z_g__R33;
        // auto bx__a_g__R11 = bx * a_g__R11;
        // auto bx__a_g__R12 = bx * a_g__R12;
        // auto bx__a_g__R13 = bx * a_g__R13;
        // auto bx__a_g__R21 = bx * a_g__R21;
        // auto bx__a_g__R22 = bx * a_g__R22;
        // auto bx__a_g__R23 = bx * a_g__R23;
        // auto bx__a_g__R31 = bx * a_g__R31;
        // auto bx__a_g__R32 = bx * a_g__R32;
        // auto bx__a_g__R33 = bx * a_g__R33;
        // auto bx__b_g__R11 = bx * b_g__R11;
        // auto bx__b_g__R12 = bx * b_g__R12;
        // auto bx__b_g__R13 = bx * b_g__R13;
        // auto bx__b_g__R21 = bx * b_g__R21;
        // auto bx__b_g__R22 = bx * b_g__R22;
        // auto bx__b_g__R23 = bx * b_g__R23;
        // auto bx__b_g__R31 = bx * b_g__R31;
        // auto bx__b_g__R32 = bx * b_g__R32;
        // auto bx__b_g__R33 = bx * b_g__R33;
        // auto bx__c_g__R11 = bx * c_g__R11;
        // auto bx__c_g__R12 = bx * c_g__R12;
        // auto bx__c_g__R13 = bx * c_g__R13;
        // auto bx__c_g__R21 = bx * c_g__R21;
        // auto bx__c_g__R22 = bx * c_g__R22;
        // auto bx__c_g__R23 = bx * c_g__R23;
        // auto bx__c_g__R31 = bx * c_g__R31;
        // auto bx__c_g__R32 = bx * c_g__R32;
        // auto bx__c_g__R33 = bx * c_g__R33;
        // auto bx__x_g__R11 = bx * x_g__R11;
        // auto bx__x_g__R12 = bx * x_g__R12;
        // auto bx__x_g__R13 = bx * x_g__R13;
        // auto bx__x_g__R21 = bx * x_g__R21;
        // auto bx__x_g__R22 = bx * x_g__R22;
        // auto bx__x_g__R23 = bx * x_g__R23;
        // auto bx__x_g__R31 = bx * x_g__R31;
        // auto bx__x_g__R32 = bx * x_g__R32;
        // auto bx__x_g__R33 = bx * x_g__R33;
        // auto bx__y_g__R11 = bx * y_g__R11;
        // auto bx__y_g__R12 = bx * y_g__R12;
        // auto bx__y_g__R13 = bx * y_g__R13;
        // auto bx__y_g__R21 = bx * y_g__R21;
        // auto bx__y_g__R22 = bx * y_g__R22;
        // auto bx__y_g__R23 = bx * y_g__R23;
        // auto bx__y_g__R31 = bx * y_g__R31;
        // auto bx__y_g__R32 = bx * y_g__R32;
        // auto bx__y_g__R33 = bx * y_g__R33;
        // auto bx__z_g__R11 = bx * z_g__R11;
        // auto bx__z_g__R12 = bx * z_g__R12;
        // auto bx__z_g__R13 = bx * z_g__R13;
        // auto bx__z_g__R21 = bx * z_g__R21;
        // auto bx__z_g__R22 = bx * z_g__R22;
        // auto bx__z_g__R23 = bx * z_g__R23;
        // auto bx__z_g__R31 = bx * z_g__R31;
        // auto bx__z_g__R32 = bx * z_g__R32;
        // auto bx__z_g__R33 = bx * z_g__R33;
        // auto by__a_g__R11 = by * a_g__R11;
        // auto by__a_g__R12 = by * a_g__R12;
        // auto by__a_g__R13 = by * a_g__R13;
        // auto by__a_g__R21 = by * a_g__R21;
        // auto by__a_g__R22 = by * a_g__R22;
        // auto by__a_g__R23 = by * a_g__R23;
        // auto by__a_g__R31 = by * a_g__R31;
        // auto by__a_g__R32 = by * a_g__R32;
        // auto by__a_g__R33 = by * a_g__R33;
        // auto by__b_g__R11 = by * b_g__R11;
        // auto by__b_g__R12 = by * b_g__R12;
        // auto by__b_g__R13 = by * b_g__R13;
        // auto by__b_g__R21 = by * b_g__R21;
        // auto by__b_g__R22 = by * b_g__R22;
        // auto by__b_g__R23 = by * b_g__R23;
        // auto by__b_g__R31 = by * b_g__R31;
        // auto by__b_g__R32 = by * b_g__R32;
        // auto by__b_g__R33 = by * b_g__R33;
        // auto by__c_g__R11 = by * c_g__R11;
        // auto by__c_g__R12 = by * c_g__R12;
        // auto by__c_g__R13 = by * c_g__R13;
        // auto by__c_g__R21 = by * c_g__R21;
        // auto by__c_g__R22 = by * c_g__R22;
        // auto by__c_g__R23 = by * c_g__R23;
        // auto by__c_g__R31 = by * c_g__R31;
        // auto by__c_g__R32 = by * c_g__R32;
        // auto by__c_g__R33 = by * c_g__R33;
        // auto by__x_g__R11 = by * x_g__R11;
        // auto by__x_g__R12 = by * x_g__R12;
        // auto by__x_g__R13 = by * x_g__R13;
        // auto by__x_g__R21 = by * x_g__R21;
        // auto by__x_g__R22 = by * x_g__R22;
        // auto by__x_g__R23 = by * x_g__R23;
        // auto by__x_g__R31 = by * x_g__R31;
        // auto by__x_g__R32 = by * x_g__R32;
        // auto by__x_g__R33 = by * x_g__R33;
        // auto by__y_g__R11 = by * y_g__R11;
        // auto by__y_g__R12 = by * y_g__R12;
        // auto by__y_g__R13 = by * y_g__R13;
        // auto by__y_g__R21 = by * y_g__R21;
        // auto by__y_g__R22 = by * y_g__R22;
        // auto by__y_g__R23 = by * y_g__R23;
        // auto by__y_g__R31 = by * y_g__R31;
        // auto by__y_g__R32 = by * y_g__R32;
        // auto by__y_g__R33 = by * y_g__R33;
        // auto by__z_g__R11 = by * z_g__R11;
        // auto by__z_g__R12 = by * z_g__R12;
        // auto by__z_g__R13 = by * z_g__R13;
        // auto by__z_g__R21 = by * z_g__R21;
        // auto by__z_g__R22 = by * z_g__R22;
        // auto by__z_g__R23 = by * z_g__R23;
        // auto by__z_g__R31 = by * z_g__R31;
        // auto by__z_g__R32 = by * z_g__R32;
        // auto by__z_g__R33 = by * z_g__R33;
        // auto bz__a_g__R11 = bz * a_g__R11;
        // auto bz__a_g__R12 = bz * a_g__R12;
        // auto bz__a_g__R13 = bz * a_g__R13;
        // auto bz__a_g__R21 = bz * a_g__R21;
        // auto bz__a_g__R22 = bz * a_g__R22;
        // auto bz__a_g__R23 = bz * a_g__R23;
        // auto bz__a_g__R31 = bz * a_g__R31;
        // auto bz__a_g__R32 = bz * a_g__R32;
        // auto bz__a_g__R33 = bz * a_g__R33;
        // auto bz__b_g__R11 = bz * b_g__R11;
        // auto bz__b_g__R12 = bz * b_g__R12;
        // auto bz__b_g__R13 = bz * b_g__R13;
        // auto bz__b_g__R21 = bz * b_g__R21;
        // auto bz__b_g__R22 = bz * b_g__R22;
        // auto bz__b_g__R23 = bz * b_g__R23;
        // auto bz__b_g__R31 = bz * b_g__R31;
        // auto bz__b_g__R32 = bz * b_g__R32;
        // auto bz__b_g__R33 = bz * b_g__R33;
        // auto bz__c_g__R11 = bz * c_g__R11;
        // auto bz__c_g__R12 = bz * c_g__R12;
        // auto bz__c_g__R13 = bz * c_g__R13;
        // auto bz__c_g__R21 = bz * c_g__R21;
        // auto bz__c_g__R22 = bz * c_g__R22;
        // auto bz__c_g__R23 = bz * c_g__R23;
        // auto bz__c_g__R31 = bz * c_g__R31;
        // auto bz__c_g__R32 = bz * c_g__R32;
        // auto bz__c_g__R33 = bz * c_g__R33;
        // auto bz__x_g__R11 = bz * x_g__R11;
        // auto bz__x_g__R12 = bz * x_g__R12;
        // auto bz__x_g__R13 = bz * x_g__R13;
        // auto bz__x_g__R21 = bz * x_g__R21;
        // auto bz__x_g__R22 = bz * x_g__R22;
        // auto bz__x_g__R23 = bz * x_g__R23;
        // auto bz__x_g__R31 = bz * x_g__R31;
        // auto bz__x_g__R32 = bz * x_g__R32;
        // auto bz__x_g__R33 = bz * x_g__R33;
        // auto bz__y_g__R11 = bz * y_g__R11;
        // auto bz__y_g__R12 = bz * y_g__R12;
        // auto bz__y_g__R13 = bz * y_g__R13;
        // auto bz__y_g__R21 = bz * y_g__R21;
        // auto bz__y_g__R22 = bz * y_g__R22;
        // auto bz__y_g__R23 = bz * y_g__R23;
        // auto bz__y_g__R31 = bz * y_g__R31;
        // auto bz__y_g__R32 = bz * y_g__R32;
        // auto bz__y_g__R33 = bz * y_g__R33;
        // auto bz__z_g__R11 = bz * z_g__R11;
        // auto bz__z_g__R12 = bz * z_g__R12;
        // auto bz__z_g__R13 = bz * z_g__R13;
        // auto bz__z_g__R21 = bz * z_g__R21;
        // auto bz__z_g__R22 = bz * z_g__R22;
        // auto bz__z_g__R23 = bz * z_g__R23;
        // auto bz__z_g__R31 = bz * z_g__R31;
        // auto bz__z_g__R32 = bz * z_g__R32;
        // auto bz__z_g__R33 = bz * z_g__R33;
        //
        // auto tx__a_g__a_c__R11 = tx * a_g__a_c__R11;
        // auto tx__a_g__a_c__R12 = tx * a_g__a_c__R12;
        // auto tx__a_g__a_c__R13 = tx * a_g__a_c__R13;
        // auto tx__a_g__a_c__R21 = tx * a_g__a_c__R21;
        // auto tx__a_g__a_c__R22 = tx * a_g__a_c__R22;
        // auto tx__a_g__a_c__R23 = tx * a_g__a_c__R23;
        // auto tx__a_g__a_c__R31 = tx * a_g__a_c__R31;
        // auto tx__a_g__a_c__R32 = tx * a_g__a_c__R32;
        // auto tx__a_g__a_c__R33 = tx * a_g__a_c__R33;
        // auto tx__a_g__b_c__R11 = tx * a_g__b_c__R11;
        // auto tx__a_g__b_c__R12 = tx * a_g__b_c__R12;
        // auto tx__a_g__b_c__R13 = tx * a_g__b_c__R13;
        // auto tx__a_g__b_c__R21 = tx * a_g__b_c__R21;
        // auto tx__a_g__b_c__R22 = tx * a_g__b_c__R22;
        // auto tx__a_g__b_c__R23 = tx * a_g__b_c__R23;
        // auto tx__a_g__b_c__R31 = tx * a_g__b_c__R31;
        // auto tx__a_g__b_c__R32 = tx * a_g__b_c__R32;
        // auto tx__a_g__b_c__R33 = tx * a_g__b_c__R33;
        // auto tx__a_g__c_c__R11 = tx * a_g__c_c__R11;
        // auto tx__a_g__c_c__R12 = tx * a_g__c_c__R12;
        // auto tx__a_g__c_c__R13 = tx * a_g__c_c__R13;
        // auto tx__a_g__c_c__R21 = tx * a_g__c_c__R21;
        // auto tx__a_g__c_c__R22 = tx * a_g__c_c__R22;
        // auto tx__a_g__c_c__R23 = tx * a_g__c_c__R23;
        // auto tx__a_g__c_c__R31 = tx * a_g__c_c__R31;
        // auto tx__a_g__c_c__R32 = tx * a_g__c_c__R32;
        // auto tx__a_g__c_c__R33 = tx * a_g__c_c__R33;
        // auto tx__a_g__x_c__R11 = tx * a_g__x_c__R11;
        // auto tx__a_g__x_c__R12 = tx * a_g__x_c__R12;
        // auto tx__a_g__x_c__R13 = tx * a_g__x_c__R13;
        // auto tx__a_g__x_c__R21 = tx * a_g__x_c__R21;
        // auto tx__a_g__x_c__R22 = tx * a_g__x_c__R22;
        // auto tx__a_g__x_c__R23 = tx * a_g__x_c__R23;
        // auto tx__a_g__x_c__R31 = tx * a_g__x_c__R31;
        // auto tx__a_g__x_c__R32 = tx * a_g__x_c__R32;
        // auto tx__a_g__x_c__R33 = tx * a_g__x_c__R33;
        // auto tx__a_g__y_c__R11 = tx * a_g__y_c__R11;
        // auto tx__a_g__y_c__R12 = tx * a_g__y_c__R12;
        // auto tx__a_g__y_c__R13 = tx * a_g__y_c__R13;
        // auto tx__a_g__y_c__R21 = tx * a_g__y_c__R21;
        // auto tx__a_g__y_c__R22 = tx * a_g__y_c__R22;
        // auto tx__a_g__y_c__R23 = tx * a_g__y_c__R23;
        // auto tx__a_g__y_c__R31 = tx * a_g__y_c__R31;
        // auto tx__a_g__y_c__R32 = tx * a_g__y_c__R32;
        // auto tx__a_g__y_c__R33 = tx * a_g__y_c__R33;
        // auto tx__a_g__z_c__R11 = tx * a_g__z_c__R11;
        // auto tx__a_g__z_c__R12 = tx * a_g__z_c__R12;
        // auto tx__a_g__z_c__R13 = tx * a_g__z_c__R13;
        // auto tx__a_g__z_c__R21 = tx * a_g__z_c__R21;
        // auto tx__a_g__z_c__R22 = tx * a_g__z_c__R22;
        // auto tx__a_g__z_c__R23 = tx * a_g__z_c__R23;
        // auto tx__a_g__z_c__R31 = tx * a_g__z_c__R31;
        // auto tx__a_g__z_c__R32 = tx * a_g__z_c__R32;
        // auto tx__a_g__z_c__R33 = tx * a_g__z_c__R33;
        // auto tx__b_g__a_c__R11 = tx * b_g__a_c__R11;
        // auto tx__b_g__a_c__R12 = tx * b_g__a_c__R12;
        // auto tx__b_g__a_c__R13 = tx * b_g__a_c__R13;
        // auto tx__b_g__a_c__R21 = tx * b_g__a_c__R21;
        // auto tx__b_g__a_c__R22 = tx * b_g__a_c__R22;
        // auto tx__b_g__a_c__R23 = tx * b_g__a_c__R23;
        // auto tx__b_g__a_c__R31 = tx * b_g__a_c__R31;
        // auto tx__b_g__a_c__R32 = tx * b_g__a_c__R32;
        // auto tx__b_g__a_c__R33 = tx * b_g__a_c__R33;
        // auto tx__b_g__b_c__R11 = tx * b_g__b_c__R11;
        // auto tx__b_g__b_c__R12 = tx * b_g__b_c__R12;
        // auto tx__b_g__b_c__R13 = tx * b_g__b_c__R13;
        // auto tx__b_g__b_c__R21 = tx * b_g__b_c__R21;
        // auto tx__b_g__b_c__R22 = tx * b_g__b_c__R22;
        // auto tx__b_g__b_c__R23 = tx * b_g__b_c__R23;
        // auto tx__b_g__b_c__R31 = tx * b_g__b_c__R31;
        // auto tx__b_g__b_c__R32 = tx * b_g__b_c__R32;
        // auto tx__b_g__b_c__R33 = tx * b_g__b_c__R33;
        // auto tx__b_g__c_c__R11 = tx * b_g__c_c__R11;
        // auto tx__b_g__c_c__R12 = tx * b_g__c_c__R12;
        // auto tx__b_g__c_c__R13 = tx * b_g__c_c__R13;
        // auto tx__b_g__c_c__R21 = tx * b_g__c_c__R21;
        // auto tx__b_g__c_c__R22 = tx * b_g__c_c__R22;
        // auto tx__b_g__c_c__R23 = tx * b_g__c_c__R23;
        // auto tx__b_g__c_c__R31 = tx * b_g__c_c__R31;
        // auto tx__b_g__c_c__R32 = tx * b_g__c_c__R32;
        // auto tx__b_g__c_c__R33 = tx * b_g__c_c__R33;
        // auto tx__b_g__x_c__R11 = tx * b_g__x_c__R11;
        // auto tx__b_g__x_c__R12 = tx * b_g__x_c__R12;
        // auto tx__b_g__x_c__R13 = tx * b_g__x_c__R13;
        // auto tx__b_g__x_c__R21 = tx * b_g__x_c__R21;
        // auto tx__b_g__x_c__R22 = tx * b_g__x_c__R22;
        // auto tx__b_g__x_c__R23 = tx * b_g__x_c__R23;
        // auto tx__b_g__x_c__R31 = tx * b_g__x_c__R31;
        // auto tx__b_g__x_c__R32 = tx * b_g__x_c__R32;
        // auto tx__b_g__x_c__R33 = tx * b_g__x_c__R33;
        // auto tx__b_g__y_c__R11 = tx * b_g__y_c__R11;
        // auto tx__b_g__y_c__R12 = tx * b_g__y_c__R12;
        // auto tx__b_g__y_c__R13 = tx * b_g__y_c__R13;
        // auto tx__b_g__y_c__R21 = tx * b_g__y_c__R21;
        // auto tx__b_g__y_c__R22 = tx * b_g__y_c__R22;
        // auto tx__b_g__y_c__R23 = tx * b_g__y_c__R23;
        // auto tx__b_g__y_c__R31 = tx * b_g__y_c__R31;
        // auto tx__b_g__y_c__R32 = tx * b_g__y_c__R32;
        // auto tx__b_g__y_c__R33 = tx * b_g__y_c__R33;
        // auto tx__b_g__z_c__R11 = tx * b_g__z_c__R11;
        // auto tx__b_g__z_c__R12 = tx * b_g__z_c__R12;
        // auto tx__b_g__z_c__R13 = tx * b_g__z_c__R13;
        // auto tx__b_g__z_c__R21 = tx * b_g__z_c__R21;
        // auto tx__b_g__z_c__R22 = tx * b_g__z_c__R22;
        // auto tx__b_g__z_c__R23 = tx * b_g__z_c__R23;
        // auto tx__b_g__z_c__R31 = tx * b_g__z_c__R31;
        // auto tx__b_g__z_c__R32 = tx * b_g__z_c__R32;
        // auto tx__b_g__z_c__R33 = tx * b_g__z_c__R33;
        // auto tx__c_g__a_c__R11 = tx * c_g__a_c__R11;
        // auto tx__c_g__a_c__R12 = tx * c_g__a_c__R12;
        // auto tx__c_g__a_c__R13 = tx * c_g__a_c__R13;
        // auto tx__c_g__a_c__R21 = tx * c_g__a_c__R21;
        // auto tx__c_g__a_c__R22 = tx * c_g__a_c__R22;
        // auto tx__c_g__a_c__R23 = tx * c_g__a_c__R23;
        // auto tx__c_g__a_c__R31 = tx * c_g__a_c__R31;
        // auto tx__c_g__a_c__R32 = tx * c_g__a_c__R32;
        // auto tx__c_g__a_c__R33 = tx * c_g__a_c__R33;
        // auto tx__c_g__b_c__R11 = tx * c_g__b_c__R11;
        // auto tx__c_g__b_c__R12 = tx * c_g__b_c__R12;
        // auto tx__c_g__b_c__R13 = tx * c_g__b_c__R13;
        // auto tx__c_g__b_c__R21 = tx * c_g__b_c__R21;
        // auto tx__c_g__b_c__R22 = tx * c_g__b_c__R22;
        // auto tx__c_g__b_c__R23 = tx * c_g__b_c__R23;
        // auto tx__c_g__b_c__R31 = tx * c_g__b_c__R31;
        // auto tx__c_g__b_c__R32 = tx * c_g__b_c__R32;
        // auto tx__c_g__b_c__R33 = tx * c_g__b_c__R33;
        // auto tx__c_g__c_c__R11 = tx * c_g__c_c__R11;
        // auto tx__c_g__c_c__R12 = tx * c_g__c_c__R12;
        // auto tx__c_g__c_c__R13 = tx * c_g__c_c__R13;
        // auto tx__c_g__c_c__R21 = tx * c_g__c_c__R21;
        // auto tx__c_g__c_c__R22 = tx * c_g__c_c__R22;
        // auto tx__c_g__c_c__R23 = tx * c_g__c_c__R23;
        // auto tx__c_g__c_c__R31 = tx * c_g__c_c__R31;
        // auto tx__c_g__c_c__R32 = tx * c_g__c_c__R32;
        // auto tx__c_g__c_c__R33 = tx * c_g__c_c__R33;
        // auto tx__c_g__x_c__R11 = tx * c_g__x_c__R11;
        // auto tx__c_g__x_c__R12 = tx * c_g__x_c__R12;
        // auto tx__c_g__x_c__R13 = tx * c_g__x_c__R13;
        // auto tx__c_g__x_c__R21 = tx * c_g__x_c__R21;
        // auto tx__c_g__x_c__R22 = tx * c_g__x_c__R22;
        // auto tx__c_g__x_c__R23 = tx * c_g__x_c__R23;
        // auto tx__c_g__x_c__R31 = tx * c_g__x_c__R31;
        // auto tx__c_g__x_c__R32 = tx * c_g__x_c__R32;
        // auto tx__c_g__x_c__R33 = tx * c_g__x_c__R33;
        // auto tx__c_g__y_c__R11 = tx * c_g__y_c__R11;
        // auto tx__c_g__y_c__R12 = tx * c_g__y_c__R12;
        // auto tx__c_g__y_c__R13 = tx * c_g__y_c__R13;
        // auto tx__c_g__y_c__R21 = tx * c_g__y_c__R21;
        // auto tx__c_g__y_c__R22 = tx * c_g__y_c__R22;
        // auto tx__c_g__y_c__R23 = tx * c_g__y_c__R23;
        // auto tx__c_g__y_c__R31 = tx * c_g__y_c__R31;
        // auto tx__c_g__y_c__R32 = tx * c_g__y_c__R32;
        // auto tx__c_g__y_c__R33 = tx * c_g__y_c__R33;
        // auto tx__c_g__z_c__R11 = tx * c_g__z_c__R11;
        // auto tx__c_g__z_c__R12 = tx * c_g__z_c__R12;
        // auto tx__c_g__z_c__R13 = tx * c_g__z_c__R13;
        // auto tx__c_g__z_c__R21 = tx * c_g__z_c__R21;
        // auto tx__c_g__z_c__R22 = tx * c_g__z_c__R22;
        // auto tx__c_g__z_c__R23 = tx * c_g__z_c__R23;
        // auto tx__c_g__z_c__R31 = tx * c_g__z_c__R31;
        // auto tx__c_g__z_c__R32 = tx * c_g__z_c__R32;
        // auto tx__c_g__z_c__R33 = tx * c_g__z_c__R33;
        // auto tx__x_g__a_c__R11 = tx * x_g__a_c__R11;
        // auto tx__x_g__a_c__R12 = tx * x_g__a_c__R12;
        // auto tx__x_g__a_c__R13 = tx * x_g__a_c__R13;
        // auto tx__x_g__a_c__R21 = tx * x_g__a_c__R21;
        // auto tx__x_g__a_c__R22 = tx * x_g__a_c__R22;
        // auto tx__x_g__a_c__R23 = tx * x_g__a_c__R23;
        // auto tx__x_g__a_c__R31 = tx * x_g__a_c__R31;
        // auto tx__x_g__a_c__R32 = tx * x_g__a_c__R32;
        // auto tx__x_g__a_c__R33 = tx * x_g__a_c__R33;
        // auto tx__x_g__b_c__R11 = tx * x_g__b_c__R11;
        // auto tx__x_g__b_c__R12 = tx * x_g__b_c__R12;
        // auto tx__x_g__b_c__R13 = tx * x_g__b_c__R13;
        // auto tx__x_g__b_c__R21 = tx * x_g__b_c__R21;
        // auto tx__x_g__b_c__R22 = tx * x_g__b_c__R22;
        // auto tx__x_g__b_c__R23 = tx * x_g__b_c__R23;
        // auto tx__x_g__b_c__R31 = tx * x_g__b_c__R31;
        // auto tx__x_g__b_c__R32 = tx * x_g__b_c__R32;
        // auto tx__x_g__b_c__R33 = tx * x_g__b_c__R33;
        // auto tx__x_g__c_c__R11 = tx * x_g__c_c__R11;
        // auto tx__x_g__c_c__R12 = tx * x_g__c_c__R12;
        // auto tx__x_g__c_c__R13 = tx * x_g__c_c__R13;
        // auto tx__x_g__c_c__R21 = tx * x_g__c_c__R21;
        // auto tx__x_g__c_c__R22 = tx * x_g__c_c__R22;
        // auto tx__x_g__c_c__R23 = tx * x_g__c_c__R23;
        // auto tx__x_g__c_c__R31 = tx * x_g__c_c__R31;
        // auto tx__x_g__c_c__R32 = tx * x_g__c_c__R32;
        // auto tx__x_g__c_c__R33 = tx * x_g__c_c__R33;
        // auto tx__x_g__x_c__R11 = tx * x_g__x_c__R11;
        // auto tx__x_g__x_c__R12 = tx * x_g__x_c__R12;
        // auto tx__x_g__x_c__R13 = tx * x_g__x_c__R13;
        // auto tx__x_g__x_c__R21 = tx * x_g__x_c__R21;
        // auto tx__x_g__x_c__R22 = tx * x_g__x_c__R22;
        // auto tx__x_g__x_c__R23 = tx * x_g__x_c__R23;
        // auto tx__x_g__x_c__R31 = tx * x_g__x_c__R31;
        // auto tx__x_g__x_c__R32 = tx * x_g__x_c__R32;
        // auto tx__x_g__x_c__R33 = tx * x_g__x_c__R33;
        // auto tx__x_g__y_c__R11 = tx * x_g__y_c__R11;
        // auto tx__x_g__y_c__R12 = tx * x_g__y_c__R12;
        // auto tx__x_g__y_c__R13 = tx * x_g__y_c__R13;
        // auto tx__x_g__y_c__R21 = tx * x_g__y_c__R21;
        // auto tx__x_g__y_c__R22 = tx * x_g__y_c__R22;
        // auto tx__x_g__y_c__R23 = tx * x_g__y_c__R23;
        // auto tx__x_g__y_c__R31 = tx * x_g__y_c__R31;
        // auto tx__x_g__y_c__R32 = tx * x_g__y_c__R32;
        // auto tx__x_g__y_c__R33 = tx * x_g__y_c__R33;
        // auto tx__x_g__z_c__R11 = tx * x_g__z_c__R11;
        // auto tx__x_g__z_c__R12 = tx * x_g__z_c__R12;
        // auto tx__x_g__z_c__R13 = tx * x_g__z_c__R13;
        // auto tx__x_g__z_c__R21 = tx * x_g__z_c__R21;
        // auto tx__x_g__z_c__R22 = tx * x_g__z_c__R22;
        // auto tx__x_g__z_c__R23 = tx * x_g__z_c__R23;
        // auto tx__x_g__z_c__R31 = tx * x_g__z_c__R31;
        // auto tx__x_g__z_c__R32 = tx * x_g__z_c__R32;
        // auto tx__x_g__z_c__R33 = tx * x_g__z_c__R33;
        // auto tx__y_g__a_c__R11 = tx * y_g__a_c__R11;
        // auto tx__y_g__a_c__R12 = tx * y_g__a_c__R12;
        // auto tx__y_g__a_c__R13 = tx * y_g__a_c__R13;
        // auto tx__y_g__a_c__R21 = tx * y_g__a_c__R21;
        // auto tx__y_g__a_c__R22 = tx * y_g__a_c__R22;
        // auto tx__y_g__a_c__R23 = tx * y_g__a_c__R23;
        // auto tx__y_g__a_c__R31 = tx * y_g__a_c__R31;
        // auto tx__y_g__a_c__R32 = tx * y_g__a_c__R32;
        // auto tx__y_g__a_c__R33 = tx * y_g__a_c__R33;
        // auto tx__y_g__b_c__R11 = tx * y_g__b_c__R11;
        // auto tx__y_g__b_c__R12 = tx * y_g__b_c__R12;
        // auto tx__y_g__b_c__R13 = tx * y_g__b_c__R13;
        // auto tx__y_g__b_c__R21 = tx * y_g__b_c__R21;
        // auto tx__y_g__b_c__R22 = tx * y_g__b_c__R22;
        // auto tx__y_g__b_c__R23 = tx * y_g__b_c__R23;
        // auto tx__y_g__b_c__R31 = tx * y_g__b_c__R31;
        // auto tx__y_g__b_c__R32 = tx * y_g__b_c__R32;
        // auto tx__y_g__b_c__R33 = tx * y_g__b_c__R33;
        // auto tx__y_g__c_c__R11 = tx * y_g__c_c__R11;
        // auto tx__y_g__c_c__R12 = tx * y_g__c_c__R12;
        // auto tx__y_g__c_c__R13 = tx * y_g__c_c__R13;
        // auto tx__y_g__c_c__R21 = tx * y_g__c_c__R21;
        // auto tx__y_g__c_c__R22 = tx * y_g__c_c__R22;
        // auto tx__y_g__c_c__R23 = tx * y_g__c_c__R23;
        // auto tx__y_g__c_c__R31 = tx * y_g__c_c__R31;
        // auto tx__y_g__c_c__R32 = tx * y_g__c_c__R32;
        // auto tx__y_g__c_c__R33 = tx * y_g__c_c__R33;
        // auto tx__y_g__x_c__R11 = tx * y_g__x_c__R11;
        // auto tx__y_g__x_c__R12 = tx * y_g__x_c__R12;
        // auto tx__y_g__x_c__R13 = tx * y_g__x_c__R13;
        // auto tx__y_g__x_c__R21 = tx * y_g__x_c__R21;
        // auto tx__y_g__x_c__R22 = tx * y_g__x_c__R22;
        // auto tx__y_g__x_c__R23 = tx * y_g__x_c__R23;
        // auto tx__y_g__x_c__R31 = tx * y_g__x_c__R31;
        // auto tx__y_g__x_c__R32 = tx * y_g__x_c__R32;
        // auto tx__y_g__x_c__R33 = tx * y_g__x_c__R33;
        // auto tx__y_g__y_c__R11 = tx * y_g__y_c__R11;
        // auto tx__y_g__y_c__R12 = tx * y_g__y_c__R12;
        // auto tx__y_g__y_c__R13 = tx * y_g__y_c__R13;
        // auto tx__y_g__y_c__R21 = tx * y_g__y_c__R21;
        // auto tx__y_g__y_c__R22 = tx * y_g__y_c__R22;
        // auto tx__y_g__y_c__R23 = tx * y_g__y_c__R23;
        // auto tx__y_g__y_c__R31 = tx * y_g__y_c__R31;
        // auto tx__y_g__y_c__R32 = tx * y_g__y_c__R32;
        // auto tx__y_g__y_c__R33 = tx * y_g__y_c__R33;
        // auto tx__y_g__z_c__R11 = tx * y_g__z_c__R11;
        // auto tx__y_g__z_c__R12 = tx * y_g__z_c__R12;
        // auto tx__y_g__z_c__R13 = tx * y_g__z_c__R13;
        // auto tx__y_g__z_c__R21 = tx * y_g__z_c__R21;
        // auto tx__y_g__z_c__R22 = tx * y_g__z_c__R22;
        // auto tx__y_g__z_c__R23 = tx * y_g__z_c__R23;
        // auto tx__y_g__z_c__R31 = tx * y_g__z_c__R31;
        // auto tx__y_g__z_c__R32 = tx * y_g__z_c__R32;
        // auto tx__y_g__z_c__R33 = tx * y_g__z_c__R33;
        // auto tx__z_g__a_c__R11 = tx * z_g__a_c__R11;
        // auto tx__z_g__a_c__R12 = tx * z_g__a_c__R12;
        // auto tx__z_g__a_c__R13 = tx * z_g__a_c__R13;
        // auto tx__z_g__a_c__R21 = tx * z_g__a_c__R21;
        // auto tx__z_g__a_c__R22 = tx * z_g__a_c__R22;
        // auto tx__z_g__a_c__R23 = tx * z_g__a_c__R23;
        // auto tx__z_g__a_c__R31 = tx * z_g__a_c__R31;
        // auto tx__z_g__a_c__R32 = tx * z_g__a_c__R32;
        // auto tx__z_g__a_c__R33 = tx * z_g__a_c__R33;
        // auto tx__z_g__b_c__R11 = tx * z_g__b_c__R11;
        // auto tx__z_g__b_c__R12 = tx * z_g__b_c__R12;
        // auto tx__z_g__b_c__R13 = tx * z_g__b_c__R13;
        // auto tx__z_g__b_c__R21 = tx * z_g__b_c__R21;
        // auto tx__z_g__b_c__R22 = tx * z_g__b_c__R22;
        // auto tx__z_g__b_c__R23 = tx * z_g__b_c__R23;
        // auto tx__z_g__b_c__R31 = tx * z_g__b_c__R31;
        // auto tx__z_g__b_c__R32 = tx * z_g__b_c__R32;
        // auto tx__z_g__b_c__R33 = tx * z_g__b_c__R33;
        // auto tx__z_g__c_c__R11 = tx * z_g__c_c__R11;
        // auto tx__z_g__c_c__R12 = tx * z_g__c_c__R12;
        // auto tx__z_g__c_c__R13 = tx * z_g__c_c__R13;
        // auto tx__z_g__c_c__R21 = tx * z_g__c_c__R21;
        // auto tx__z_g__c_c__R22 = tx * z_g__c_c__R22;
        // auto tx__z_g__c_c__R23 = tx * z_g__c_c__R23;
        // auto tx__z_g__c_c__R31 = tx * z_g__c_c__R31;
        // auto tx__z_g__c_c__R32 = tx * z_g__c_c__R32;
        // auto tx__z_g__c_c__R33 = tx * z_g__c_c__R33;
        // auto tx__z_g__x_c__R11 = tx * z_g__x_c__R11;
        // auto tx__z_g__x_c__R12 = tx * z_g__x_c__R12;
        // auto tx__z_g__x_c__R13 = tx * z_g__x_c__R13;
        // auto tx__z_g__x_c__R21 = tx * z_g__x_c__R21;
        // auto tx__z_g__x_c__R22 = tx * z_g__x_c__R22;
        // auto tx__z_g__x_c__R23 = tx * z_g__x_c__R23;
        // auto tx__z_g__x_c__R31 = tx * z_g__x_c__R31;
        // auto tx__z_g__x_c__R32 = tx * z_g__x_c__R32;
        // auto tx__z_g__x_c__R33 = tx * z_g__x_c__R33;
        // auto tx__z_g__y_c__R11 = tx * z_g__y_c__R11;
        // auto tx__z_g__y_c__R12 = tx * z_g__y_c__R12;
        // auto tx__z_g__y_c__R13 = tx * z_g__y_c__R13;
        // auto tx__z_g__y_c__R21 = tx * z_g__y_c__R21;
        // auto tx__z_g__y_c__R22 = tx * z_g__y_c__R22;
        // auto tx__z_g__y_c__R23 = tx * z_g__y_c__R23;
        // auto tx__z_g__y_c__R31 = tx * z_g__y_c__R31;
        // auto tx__z_g__y_c__R32 = tx * z_g__y_c__R32;
        // auto tx__z_g__y_c__R33 = tx * z_g__y_c__R33;
        // auto tx__z_g__z_c__R11 = tx * z_g__z_c__R11;
        // auto tx__z_g__z_c__R12 = tx * z_g__z_c__R12;
        // auto tx__z_g__z_c__R13 = tx * z_g__z_c__R13;
        // auto tx__z_g__z_c__R21 = tx * z_g__z_c__R21;
        // auto tx__z_g__z_c__R22 = tx * z_g__z_c__R22;
        // auto tx__z_g__z_c__R23 = tx * z_g__z_c__R23;
        // auto tx__z_g__z_c__R31 = tx * z_g__z_c__R31;
        // auto tx__z_g__z_c__R32 = tx * z_g__z_c__R32;
        // auto tx__z_g__z_c__R33 = tx * z_g__z_c__R33;
        // auto ty__a_g__a_c__R11 = ty * a_g__a_c__R11;
        // auto ty__a_g__a_c__R12 = ty * a_g__a_c__R12;
        // auto ty__a_g__a_c__R13 = ty * a_g__a_c__R13;
        // auto ty__a_g__a_c__R21 = ty * a_g__a_c__R21;
        // auto ty__a_g__a_c__R22 = ty * a_g__a_c__R22;
        // auto ty__a_g__a_c__R23 = ty * a_g__a_c__R23;
        // auto ty__a_g__a_c__R31 = ty * a_g__a_c__R31;
        // auto ty__a_g__a_c__R32 = ty * a_g__a_c__R32;
        // auto ty__a_g__a_c__R33 = ty * a_g__a_c__R33;
        // auto ty__a_g__b_c__R11 = ty * a_g__b_c__R11;
        // auto ty__a_g__b_c__R12 = ty * a_g__b_c__R12;
        // auto ty__a_g__b_c__R13 = ty * a_g__b_c__R13;
        // auto ty__a_g__b_c__R21 = ty * a_g__b_c__R21;
        // auto ty__a_g__b_c__R22 = ty * a_g__b_c__R22;
        // auto ty__a_g__b_c__R23 = ty * a_g__b_c__R23;
        // auto ty__a_g__b_c__R31 = ty * a_g__b_c__R31;
        // auto ty__a_g__b_c__R32 = ty * a_g__b_c__R32;
        // auto ty__a_g__b_c__R33 = ty * a_g__b_c__R33;
        // auto ty__a_g__c_c__R11 = ty * a_g__c_c__R11;
        // auto ty__a_g__c_c__R12 = ty * a_g__c_c__R12;
        // auto ty__a_g__c_c__R13 = ty * a_g__c_c__R13;
        // auto ty__a_g__c_c__R21 = ty * a_g__c_c__R21;
        // auto ty__a_g__c_c__R22 = ty * a_g__c_c__R22;
        // auto ty__a_g__c_c__R23 = ty * a_g__c_c__R23;
        // auto ty__a_g__c_c__R31 = ty * a_g__c_c__R31;
        // auto ty__a_g__c_c__R32 = ty * a_g__c_c__R32;
        // auto ty__a_g__c_c__R33 = ty * a_g__c_c__R33;
        // auto ty__a_g__x_c__R11 = ty * a_g__x_c__R11;
        // auto ty__a_g__x_c__R12 = ty * a_g__x_c__R12;
        // auto ty__a_g__x_c__R13 = ty * a_g__x_c__R13;
        // auto ty__a_g__x_c__R21 = ty * a_g__x_c__R21;
        // auto ty__a_g__x_c__R22 = ty * a_g__x_c__R22;
        // auto ty__a_g__x_c__R23 = ty * a_g__x_c__R23;
        // auto ty__a_g__x_c__R31 = ty * a_g__x_c__R31;
        // auto ty__a_g__x_c__R32 = ty * a_g__x_c__R32;
        // auto ty__a_g__x_c__R33 = ty * a_g__x_c__R33;
        // auto ty__a_g__y_c__R11 = ty * a_g__y_c__R11;
        // auto ty__a_g__y_c__R12 = ty * a_g__y_c__R12;
        // auto ty__a_g__y_c__R13 = ty * a_g__y_c__R13;
        // auto ty__a_g__y_c__R21 = ty * a_g__y_c__R21;
        // auto ty__a_g__y_c__R22 = ty * a_g__y_c__R22;
        // auto ty__a_g__y_c__R23 = ty * a_g__y_c__R23;
        // auto ty__a_g__y_c__R31 = ty * a_g__y_c__R31;
        // auto ty__a_g__y_c__R32 = ty * a_g__y_c__R32;
        // auto ty__a_g__y_c__R33 = ty * a_g__y_c__R33;
        // auto ty__a_g__z_c__R11 = ty * a_g__z_c__R11;
        // auto ty__a_g__z_c__R12 = ty * a_g__z_c__R12;
        // auto ty__a_g__z_c__R13 = ty * a_g__z_c__R13;
        // auto ty__a_g__z_c__R21 = ty * a_g__z_c__R21;
        // auto ty__a_g__z_c__R22 = ty * a_g__z_c__R22;
        // auto ty__a_g__z_c__R23 = ty * a_g__z_c__R23;
        // auto ty__a_g__z_c__R31 = ty * a_g__z_c__R31;
        // auto ty__a_g__z_c__R32 = ty * a_g__z_c__R32;
        // auto ty__a_g__z_c__R33 = ty * a_g__z_c__R33;
        // auto ty__b_g__a_c__R11 = ty * b_g__a_c__R11;
        // auto ty__b_g__a_c__R12 = ty * b_g__a_c__R12;
        // auto ty__b_g__a_c__R13 = ty * b_g__a_c__R13;
        // auto ty__b_g__a_c__R21 = ty * b_g__a_c__R21;
        // auto ty__b_g__a_c__R22 = ty * b_g__a_c__R22;
        // auto ty__b_g__a_c__R23 = ty * b_g__a_c__R23;
        // auto ty__b_g__a_c__R31 = ty * b_g__a_c__R31;
        // auto ty__b_g__a_c__R32 = ty * b_g__a_c__R32;
        // auto ty__b_g__a_c__R33 = ty * b_g__a_c__R33;
        // auto ty__b_g__b_c__R11 = ty * b_g__b_c__R11;
        // auto ty__b_g__b_c__R12 = ty * b_g__b_c__R12;
        // auto ty__b_g__b_c__R13 = ty * b_g__b_c__R13;
        // auto ty__b_g__b_c__R21 = ty * b_g__b_c__R21;
        // auto ty__b_g__b_c__R22 = ty * b_g__b_c__R22;
        // auto ty__b_g__b_c__R23 = ty * b_g__b_c__R23;
        // auto ty__b_g__b_c__R31 = ty * b_g__b_c__R31;
        // auto ty__b_g__b_c__R32 = ty * b_g__b_c__R32;
        // auto ty__b_g__b_c__R33 = ty * b_g__b_c__R33;
        // auto ty__b_g__c_c__R11 = ty * b_g__c_c__R11;
        // auto ty__b_g__c_c__R12 = ty * b_g__c_c__R12;
        // auto ty__b_g__c_c__R13 = ty * b_g__c_c__R13;
        // auto ty__b_g__c_c__R21 = ty * b_g__c_c__R21;
        // auto ty__b_g__c_c__R22 = ty * b_g__c_c__R22;
        // auto ty__b_g__c_c__R23 = ty * b_g__c_c__R23;
        // auto ty__b_g__c_c__R31 = ty * b_g__c_c__R31;
        // auto ty__b_g__c_c__R32 = ty * b_g__c_c__R32;
        // auto ty__b_g__c_c__R33 = ty * b_g__c_c__R33;
        // auto ty__b_g__x_c__R11 = ty * b_g__x_c__R11;
        // auto ty__b_g__x_c__R12 = ty * b_g__x_c__R12;
        // auto ty__b_g__x_c__R13 = ty * b_g__x_c__R13;
        // auto ty__b_g__x_c__R21 = ty * b_g__x_c__R21;
        // auto ty__b_g__x_c__R22 = ty * b_g__x_c__R22;
        // auto ty__b_g__x_c__R23 = ty * b_g__x_c__R23;
        // auto ty__b_g__x_c__R31 = ty * b_g__x_c__R31;
        // auto ty__b_g__x_c__R32 = ty * b_g__x_c__R32;
        // auto ty__b_g__x_c__R33 = ty * b_g__x_c__R33;
        // auto ty__b_g__y_c__R11 = ty * b_g__y_c__R11;
        // auto ty__b_g__y_c__R12 = ty * b_g__y_c__R12;
        // auto ty__b_g__y_c__R13 = ty * b_g__y_c__R13;
        // auto ty__b_g__y_c__R21 = ty * b_g__y_c__R21;
        // auto ty__b_g__y_c__R22 = ty * b_g__y_c__R22;
        // auto ty__b_g__y_c__R23 = ty * b_g__y_c__R23;
        // auto ty__b_g__y_c__R31 = ty * b_g__y_c__R31;
        // auto ty__b_g__y_c__R32 = ty * b_g__y_c__R32;
        // auto ty__b_g__y_c__R33 = ty * b_g__y_c__R33;
        // auto ty__b_g__z_c__R11 = ty * b_g__z_c__R11;
        // auto ty__b_g__z_c__R12 = ty * b_g__z_c__R12;
        // auto ty__b_g__z_c__R13 = ty * b_g__z_c__R13;
        // auto ty__b_g__z_c__R21 = ty * b_g__z_c__R21;
        // auto ty__b_g__z_c__R22 = ty * b_g__z_c__R22;
        // auto ty__b_g__z_c__R23 = ty * b_g__z_c__R23;
        // auto ty__b_g__z_c__R31 = ty * b_g__z_c__R31;
        // auto ty__b_g__z_c__R32 = ty * b_g__z_c__R32;
        // auto ty__b_g__z_c__R33 = ty * b_g__z_c__R33;
        // auto ty__c_g__a_c__R11 = ty * c_g__a_c__R11;
        // auto ty__c_g__a_c__R12 = ty * c_g__a_c__R12;
        // auto ty__c_g__a_c__R13 = ty * c_g__a_c__R13;
        // auto ty__c_g__a_c__R21 = ty * c_g__a_c__R21;
        // auto ty__c_g__a_c__R22 = ty * c_g__a_c__R22;
        // auto ty__c_g__a_c__R23 = ty * c_g__a_c__R23;
        // auto ty__c_g__a_c__R31 = ty * c_g__a_c__R31;
        // auto ty__c_g__a_c__R32 = ty * c_g__a_c__R32;
        // auto ty__c_g__a_c__R33 = ty * c_g__a_c__R33;
        // auto ty__c_g__b_c__R11 = ty * c_g__b_c__R11;
        // auto ty__c_g__b_c__R12 = ty * c_g__b_c__R12;
        // auto ty__c_g__b_c__R13 = ty * c_g__b_c__R13;
        // auto ty__c_g__b_c__R21 = ty * c_g__b_c__R21;
        // auto ty__c_g__b_c__R22 = ty * c_g__b_c__R22;
        // auto ty__c_g__b_c__R23 = ty * c_g__b_c__R23;
        // auto ty__c_g__b_c__R31 = ty * c_g__b_c__R31;
        // auto ty__c_g__b_c__R32 = ty * c_g__b_c__R32;
        // auto ty__c_g__b_c__R33 = ty * c_g__b_c__R33;
        // auto ty__c_g__c_c__R11 = ty * c_g__c_c__R11;
        // auto ty__c_g__c_c__R12 = ty * c_g__c_c__R12;
        // auto ty__c_g__c_c__R13 = ty * c_g__c_c__R13;
        // auto ty__c_g__c_c__R21 = ty * c_g__c_c__R21;
        // auto ty__c_g__c_c__R22 = ty * c_g__c_c__R22;
        // auto ty__c_g__c_c__R23 = ty * c_g__c_c__R23;
        // auto ty__c_g__c_c__R31 = ty * c_g__c_c__R31;
        // auto ty__c_g__c_c__R32 = ty * c_g__c_c__R32;
        // auto ty__c_g__c_c__R33 = ty * c_g__c_c__R33;
        // auto ty__c_g__x_c__R11 = ty * c_g__x_c__R11;
        // auto ty__c_g__x_c__R12 = ty * c_g__x_c__R12;
        // auto ty__c_g__x_c__R13 = ty * c_g__x_c__R13;
        // auto ty__c_g__x_c__R21 = ty * c_g__x_c__R21;
        // auto ty__c_g__x_c__R22 = ty * c_g__x_c__R22;
        // auto ty__c_g__x_c__R23 = ty * c_g__x_c__R23;
        // auto ty__c_g__x_c__R31 = ty * c_g__x_c__R31;
        // auto ty__c_g__x_c__R32 = ty * c_g__x_c__R32;
        // auto ty__c_g__x_c__R33 = ty * c_g__x_c__R33;
        // auto ty__c_g__y_c__R11 = ty * c_g__y_c__R11;
        // auto ty__c_g__y_c__R12 = ty * c_g__y_c__R12;
        // auto ty__c_g__y_c__R13 = ty * c_g__y_c__R13;
        // auto ty__c_g__y_c__R21 = ty * c_g__y_c__R21;
        // auto ty__c_g__y_c__R22 = ty * c_g__y_c__R22;
        // auto ty__c_g__y_c__R23 = ty * c_g__y_c__R23;
        // auto ty__c_g__y_c__R31 = ty * c_g__y_c__R31;
        // auto ty__c_g__y_c__R32 = ty * c_g__y_c__R32;
        // auto ty__c_g__y_c__R33 = ty * c_g__y_c__R33;
        // auto ty__c_g__z_c__R11 = ty * c_g__z_c__R11;
        // auto ty__c_g__z_c__R12 = ty * c_g__z_c__R12;
        // auto ty__c_g__z_c__R13 = ty * c_g__z_c__R13;
        // auto ty__c_g__z_c__R21 = ty * c_g__z_c__R21;
        // auto ty__c_g__z_c__R22 = ty * c_g__z_c__R22;
        // auto ty__c_g__z_c__R23 = ty * c_g__z_c__R23;
        // auto ty__c_g__z_c__R31 = ty * c_g__z_c__R31;
        // auto ty__c_g__z_c__R32 = ty * c_g__z_c__R32;
        // auto ty__c_g__z_c__R33 = ty * c_g__z_c__R33;
        // auto ty__x_g__a_c__R11 = ty * x_g__a_c__R11;
        // auto ty__x_g__a_c__R12 = ty * x_g__a_c__R12;
        // auto ty__x_g__a_c__R13 = ty * x_g__a_c__R13;
        // auto ty__x_g__a_c__R21 = ty * x_g__a_c__R21;
        // auto ty__x_g__a_c__R22 = ty * x_g__a_c__R22;
        // auto ty__x_g__a_c__R23 = ty * x_g__a_c__R23;
        // auto ty__x_g__a_c__R31 = ty * x_g__a_c__R31;
        // auto ty__x_g__a_c__R32 = ty * x_g__a_c__R32;
        // auto ty__x_g__a_c__R33 = ty * x_g__a_c__R33;
        // auto ty__x_g__b_c__R11 = ty * x_g__b_c__R11;
        // auto ty__x_g__b_c__R12 = ty * x_g__b_c__R12;
        // auto ty__x_g__b_c__R13 = ty * x_g__b_c__R13;
        // auto ty__x_g__b_c__R21 = ty * x_g__b_c__R21;
        // auto ty__x_g__b_c__R22 = ty * x_g__b_c__R22;
        // auto ty__x_g__b_c__R23 = ty * x_g__b_c__R23;
        // auto ty__x_g__b_c__R31 = ty * x_g__b_c__R31;
        // auto ty__x_g__b_c__R32 = ty * x_g__b_c__R32;
        // auto ty__x_g__b_c__R33 = ty * x_g__b_c__R33;
        // auto ty__x_g__c_c__R11 = ty * x_g__c_c__R11;
        // auto ty__x_g__c_c__R12 = ty * x_g__c_c__R12;
        // auto ty__x_g__c_c__R13 = ty * x_g__c_c__R13;
        // auto ty__x_g__c_c__R21 = ty * x_g__c_c__R21;
        // auto ty__x_g__c_c__R22 = ty * x_g__c_c__R22;
        // auto ty__x_g__c_c__R23 = ty * x_g__c_c__R23;
        // auto ty__x_g__c_c__R31 = ty * x_g__c_c__R31;
        // auto ty__x_g__c_c__R32 = ty * x_g__c_c__R32;
        // auto ty__x_g__c_c__R33 = ty * x_g__c_c__R33;
        // auto ty__x_g__x_c__R11 = ty * x_g__x_c__R11;
        // auto ty__x_g__x_c__R12 = ty * x_g__x_c__R12;
        // auto ty__x_g__x_c__R13 = ty * x_g__x_c__R13;
        // auto ty__x_g__x_c__R21 = ty * x_g__x_c__R21;
        // auto ty__x_g__x_c__R22 = ty * x_g__x_c__R22;
        // auto ty__x_g__x_c__R23 = ty * x_g__x_c__R23;
        // auto ty__x_g__x_c__R31 = ty * x_g__x_c__R31;
        // auto ty__x_g__x_c__R32 = ty * x_g__x_c__R32;
        // auto ty__x_g__x_c__R33 = ty * x_g__x_c__R33;
        // auto ty__x_g__y_c__R11 = ty * x_g__y_c__R11;
        // auto ty__x_g__y_c__R12 = ty * x_g__y_c__R12;
        // auto ty__x_g__y_c__R13 = ty * x_g__y_c__R13;
        // auto ty__x_g__y_c__R21 = ty * x_g__y_c__R21;
        // auto ty__x_g__y_c__R22 = ty * x_g__y_c__R22;
        // auto ty__x_g__y_c__R23 = ty * x_g__y_c__R23;
        // auto ty__x_g__y_c__R31 = ty * x_g__y_c__R31;
        // auto ty__x_g__y_c__R32 = ty * x_g__y_c__R32;
        // auto ty__x_g__y_c__R33 = ty * x_g__y_c__R33;
        // auto ty__x_g__z_c__R11 = ty * x_g__z_c__R11;
        // auto ty__x_g__z_c__R12 = ty * x_g__z_c__R12;
        // auto ty__x_g__z_c__R13 = ty * x_g__z_c__R13;
        // auto ty__x_g__z_c__R21 = ty * x_g__z_c__R21;
        // auto ty__x_g__z_c__R22 = ty * x_g__z_c__R22;
        // auto ty__x_g__z_c__R23 = ty * x_g__z_c__R23;
        // auto ty__x_g__z_c__R31 = ty * x_g__z_c__R31;
        // auto ty__x_g__z_c__R32 = ty * x_g__z_c__R32;
        // auto ty__x_g__z_c__R33 = ty * x_g__z_c__R33;
        // auto ty__y_g__a_c__R11 = ty * y_g__a_c__R11;
        // auto ty__y_g__a_c__R12 = ty * y_g__a_c__R12;
        // auto ty__y_g__a_c__R13 = ty * y_g__a_c__R13;
        // auto ty__y_g__a_c__R21 = ty * y_g__a_c__R21;
        // auto ty__y_g__a_c__R22 = ty * y_g__a_c__R22;
        // auto ty__y_g__a_c__R23 = ty * y_g__a_c__R23;
        // auto ty__y_g__a_c__R31 = ty * y_g__a_c__R31;
        // auto ty__y_g__a_c__R32 = ty * y_g__a_c__R32;
        // auto ty__y_g__a_c__R33 = ty * y_g__a_c__R33;
        // auto ty__y_g__b_c__R11 = ty * y_g__b_c__R11;
        // auto ty__y_g__b_c__R12 = ty * y_g__b_c__R12;
        // auto ty__y_g__b_c__R13 = ty * y_g__b_c__R13;
        // auto ty__y_g__b_c__R21 = ty * y_g__b_c__R21;
        // auto ty__y_g__b_c__R22 = ty * y_g__b_c__R22;
        // auto ty__y_g__b_c__R23 = ty * y_g__b_c__R23;
        // auto ty__y_g__b_c__R31 = ty * y_g__b_c__R31;
        // auto ty__y_g__b_c__R32 = ty * y_g__b_c__R32;
        // auto ty__y_g__b_c__R33 = ty * y_g__b_c__R33;
        // auto ty__y_g__c_c__R11 = ty * y_g__c_c__R11;
        // auto ty__y_g__c_c__R12 = ty * y_g__c_c__R12;
        // auto ty__y_g__c_c__R13 = ty * y_g__c_c__R13;
        // auto ty__y_g__c_c__R21 = ty * y_g__c_c__R21;
        // auto ty__y_g__c_c__R22 = ty * y_g__c_c__R22;
        // auto ty__y_g__c_c__R23 = ty * y_g__c_c__R23;
        // auto ty__y_g__c_c__R31 = ty * y_g__c_c__R31;
        // auto ty__y_g__c_c__R32 = ty * y_g__c_c__R32;
        // auto ty__y_g__c_c__R33 = ty * y_g__c_c__R33;
        // auto ty__y_g__x_c__R11 = ty * y_g__x_c__R11;
        // auto ty__y_g__x_c__R12 = ty * y_g__x_c__R12;
        // auto ty__y_g__x_c__R13 = ty * y_g__x_c__R13;
        // auto ty__y_g__x_c__R21 = ty * y_g__x_c__R21;
        // auto ty__y_g__x_c__R22 = ty * y_g__x_c__R22;
        // auto ty__y_g__x_c__R23 = ty * y_g__x_c__R23;
        // auto ty__y_g__x_c__R31 = ty * y_g__x_c__R31;
        // auto ty__y_g__x_c__R32 = ty * y_g__x_c__R32;
        // auto ty__y_g__x_c__R33 = ty * y_g__x_c__R33;
        // auto ty__y_g__y_c__R11 = ty * y_g__y_c__R11;
        // auto ty__y_g__y_c__R12 = ty * y_g__y_c__R12;
        // auto ty__y_g__y_c__R13 = ty * y_g__y_c__R13;
        // auto ty__y_g__y_c__R21 = ty * y_g__y_c__R21;
        // auto ty__y_g__y_c__R22 = ty * y_g__y_c__R22;
        // auto ty__y_g__y_c__R23 = ty * y_g__y_c__R23;
        // auto ty__y_g__y_c__R31 = ty * y_g__y_c__R31;
        // auto ty__y_g__y_c__R32 = ty * y_g__y_c__R32;
        // auto ty__y_g__y_c__R33 = ty * y_g__y_c__R33;
        // auto ty__y_g__z_c__R11 = ty * y_g__z_c__R11;
        // auto ty__y_g__z_c__R12 = ty * y_g__z_c__R12;
        // auto ty__y_g__z_c__R13 = ty * y_g__z_c__R13;
        // auto ty__y_g__z_c__R21 = ty * y_g__z_c__R21;
        // auto ty__y_g__z_c__R22 = ty * y_g__z_c__R22;
        // auto ty__y_g__z_c__R23 = ty * y_g__z_c__R23;
        // auto ty__y_g__z_c__R31 = ty * y_g__z_c__R31;
        // auto ty__y_g__z_c__R32 = ty * y_g__z_c__R32;
        // auto ty__y_g__z_c__R33 = ty * y_g__z_c__R33;
        // auto ty__z_g__a_c__R11 = ty * z_g__a_c__R11;
        // auto ty__z_g__a_c__R12 = ty * z_g__a_c__R12;
        // auto ty__z_g__a_c__R13 = ty * z_g__a_c__R13;
        // auto ty__z_g__a_c__R21 = ty * z_g__a_c__R21;
        // auto ty__z_g__a_c__R22 = ty * z_g__a_c__R22;
        // auto ty__z_g__a_c__R23 = ty * z_g__a_c__R23;
        // auto ty__z_g__a_c__R31 = ty * z_g__a_c__R31;
        // auto ty__z_g__a_c__R32 = ty * z_g__a_c__R32;
        // auto ty__z_g__a_c__R33 = ty * z_g__a_c__R33;
        // auto ty__z_g__b_c__R11 = ty * z_g__b_c__R11;
        // auto ty__z_g__b_c__R12 = ty * z_g__b_c__R12;
        // auto ty__z_g__b_c__R13 = ty * z_g__b_c__R13;
        // auto ty__z_g__b_c__R21 = ty * z_g__b_c__R21;
        // auto ty__z_g__b_c__R22 = ty * z_g__b_c__R22;
        // auto ty__z_g__b_c__R23 = ty * z_g__b_c__R23;
        // auto ty__z_g__b_c__R31 = ty * z_g__b_c__R31;
        // auto ty__z_g__b_c__R32 = ty * z_g__b_c__R32;
        // auto ty__z_g__b_c__R33 = ty * z_g__b_c__R33;
        // auto ty__z_g__c_c__R11 = ty * z_g__c_c__R11;
        // auto ty__z_g__c_c__R12 = ty * z_g__c_c__R12;
        // auto ty__z_g__c_c__R13 = ty * z_g__c_c__R13;
        // auto ty__z_g__c_c__R21 = ty * z_g__c_c__R21;
        // auto ty__z_g__c_c__R22 = ty * z_g__c_c__R22;
        // auto ty__z_g__c_c__R23 = ty * z_g__c_c__R23;
        // auto ty__z_g__c_c__R31 = ty * z_g__c_c__R31;
        // auto ty__z_g__c_c__R32 = ty * z_g__c_c__R32;
        // auto ty__z_g__c_c__R33 = ty * z_g__c_c__R33;
        // auto ty__z_g__x_c__R11 = ty * z_g__x_c__R11;
        // auto ty__z_g__x_c__R12 = ty * z_g__x_c__R12;
        // auto ty__z_g__x_c__R13 = ty * z_g__x_c__R13;
        // auto ty__z_g__x_c__R21 = ty * z_g__x_c__R21;
        // auto ty__z_g__x_c__R22 = ty * z_g__x_c__R22;
        // auto ty__z_g__x_c__R23 = ty * z_g__x_c__R23;
        // auto ty__z_g__x_c__R31 = ty * z_g__x_c__R31;
        // auto ty__z_g__x_c__R32 = ty * z_g__x_c__R32;
        // auto ty__z_g__x_c__R33 = ty * z_g__x_c__R33;
        // auto ty__z_g__y_c__R11 = ty * z_g__y_c__R11;
        // auto ty__z_g__y_c__R12 = ty * z_g__y_c__R12;
        // auto ty__z_g__y_c__R13 = ty * z_g__y_c__R13;
        // auto ty__z_g__y_c__R21 = ty * z_g__y_c__R21;
        // auto ty__z_g__y_c__R22 = ty * z_g__y_c__R22;
        // auto ty__z_g__y_c__R23 = ty * z_g__y_c__R23;
        // auto ty__z_g__y_c__R31 = ty * z_g__y_c__R31;
        // auto ty__z_g__y_c__R32 = ty * z_g__y_c__R32;
        // auto ty__z_g__y_c__R33 = ty * z_g__y_c__R33;
        // auto ty__z_g__z_c__R11 = ty * z_g__z_c__R11;
        // auto ty__z_g__z_c__R12 = ty * z_g__z_c__R12;
        // auto ty__z_g__z_c__R13 = ty * z_g__z_c__R13;
        // auto ty__z_g__z_c__R21 = ty * z_g__z_c__R21;
        // auto ty__z_g__z_c__R22 = ty * z_g__z_c__R22;
        // auto ty__z_g__z_c__R23 = ty * z_g__z_c__R23;
        // auto ty__z_g__z_c__R31 = ty * z_g__z_c__R31;
        // auto ty__z_g__z_c__R32 = ty * z_g__z_c__R32;
        // auto ty__z_g__z_c__R33 = ty * z_g__z_c__R33;
        // auto tz__a_g__a_c__R11 = tz * a_g__a_c__R11;
        // auto tz__a_g__a_c__R12 = tz * a_g__a_c__R12;
        // auto tz__a_g__a_c__R13 = tz * a_g__a_c__R13;
        // auto tz__a_g__a_c__R21 = tz * a_g__a_c__R21;
        // auto tz__a_g__a_c__R22 = tz * a_g__a_c__R22;
        // auto tz__a_g__a_c__R23 = tz * a_g__a_c__R23;
        // auto tz__a_g__a_c__R31 = tz * a_g__a_c__R31;
        // auto tz__a_g__a_c__R32 = tz * a_g__a_c__R32;
        // auto tz__a_g__a_c__R33 = tz * a_g__a_c__R33;
        // auto tz__a_g__b_c__R11 = tz * a_g__b_c__R11;
        // auto tz__a_g__b_c__R12 = tz * a_g__b_c__R12;
        // auto tz__a_g__b_c__R13 = tz * a_g__b_c__R13;
        // auto tz__a_g__b_c__R21 = tz * a_g__b_c__R21;
        // auto tz__a_g__b_c__R22 = tz * a_g__b_c__R22;
        // auto tz__a_g__b_c__R23 = tz * a_g__b_c__R23;
        // auto tz__a_g__b_c__R31 = tz * a_g__b_c__R31;
        // auto tz__a_g__b_c__R32 = tz * a_g__b_c__R32;
        // auto tz__a_g__b_c__R33 = tz * a_g__b_c__R33;
        // auto tz__a_g__c_c__R11 = tz * a_g__c_c__R11;
        // auto tz__a_g__c_c__R12 = tz * a_g__c_c__R12;
        // auto tz__a_g__c_c__R13 = tz * a_g__c_c__R13;
        // auto tz__a_g__c_c__R21 = tz * a_g__c_c__R21;
        // auto tz__a_g__c_c__R22 = tz * a_g__c_c__R22;
        // auto tz__a_g__c_c__R23 = tz * a_g__c_c__R23;
        // auto tz__a_g__c_c__R31 = tz * a_g__c_c__R31;
        // auto tz__a_g__c_c__R32 = tz * a_g__c_c__R32;
        // auto tz__a_g__c_c__R33 = tz * a_g__c_c__R33;
        // auto tz__a_g__x_c__R11 = tz * a_g__x_c__R11;
        // auto tz__a_g__x_c__R12 = tz * a_g__x_c__R12;
        // auto tz__a_g__x_c__R13 = tz * a_g__x_c__R13;
        // auto tz__a_g__x_c__R21 = tz * a_g__x_c__R21;
        // auto tz__a_g__x_c__R22 = tz * a_g__x_c__R22;
        // auto tz__a_g__x_c__R23 = tz * a_g__x_c__R23;
        // auto tz__a_g__x_c__R31 = tz * a_g__x_c__R31;
        // auto tz__a_g__x_c__R32 = tz * a_g__x_c__R32;
        // auto tz__a_g__x_c__R33 = tz * a_g__x_c__R33;
        // auto tz__a_g__y_c__R11 = tz * a_g__y_c__R11;
        // auto tz__a_g__y_c__R12 = tz * a_g__y_c__R12;
        // auto tz__a_g__y_c__R13 = tz * a_g__y_c__R13;
        // auto tz__a_g__y_c__R21 = tz * a_g__y_c__R21;
        // auto tz__a_g__y_c__R22 = tz * a_g__y_c__R22;
        // auto tz__a_g__y_c__R23 = tz * a_g__y_c__R23;
        // auto tz__a_g__y_c__R31 = tz * a_g__y_c__R31;
        // auto tz__a_g__y_c__R32 = tz * a_g__y_c__R32;
        // auto tz__a_g__y_c__R33 = tz * a_g__y_c__R33;
        // auto tz__a_g__z_c__R11 = tz * a_g__z_c__R11;
        // auto tz__a_g__z_c__R12 = tz * a_g__z_c__R12;
        // auto tz__a_g__z_c__R13 = tz * a_g__z_c__R13;
        // auto tz__a_g__z_c__R21 = tz * a_g__z_c__R21;
        // auto tz__a_g__z_c__R22 = tz * a_g__z_c__R22;
        // auto tz__a_g__z_c__R23 = tz * a_g__z_c__R23;
        // auto tz__a_g__z_c__R31 = tz * a_g__z_c__R31;
        // auto tz__a_g__z_c__R32 = tz * a_g__z_c__R32;
        // auto tz__a_g__z_c__R33 = tz * a_g__z_c__R33;
        // auto tz__b_g__a_c__R11 = tz * b_g__a_c__R11;
        // auto tz__b_g__a_c__R12 = tz * b_g__a_c__R12;
        // auto tz__b_g__a_c__R13 = tz * b_g__a_c__R13;
        // auto tz__b_g__a_c__R21 = tz * b_g__a_c__R21;
        // auto tz__b_g__a_c__R22 = tz * b_g__a_c__R22;
        // auto tz__b_g__a_c__R23 = tz * b_g__a_c__R23;
        // auto tz__b_g__a_c__R31 = tz * b_g__a_c__R31;
        // auto tz__b_g__a_c__R32 = tz * b_g__a_c__R32;
        // auto tz__b_g__a_c__R33 = tz * b_g__a_c__R33;
        // auto tz__b_g__b_c__R11 = tz * b_g__b_c__R11;
        // auto tz__b_g__b_c__R12 = tz * b_g__b_c__R12;
        // auto tz__b_g__b_c__R13 = tz * b_g__b_c__R13;
        // auto tz__b_g__b_c__R21 = tz * b_g__b_c__R21;
        // auto tz__b_g__b_c__R22 = tz * b_g__b_c__R22;
        // auto tz__b_g__b_c__R23 = tz * b_g__b_c__R23;
        // auto tz__b_g__b_c__R31 = tz * b_g__b_c__R31;
        // auto tz__b_g__b_c__R32 = tz * b_g__b_c__R32;
        // auto tz__b_g__b_c__R33 = tz * b_g__b_c__R33;
        // auto tz__b_g__c_c__R11 = tz * b_g__c_c__R11;
        // auto tz__b_g__c_c__R12 = tz * b_g__c_c__R12;
        // auto tz__b_g__c_c__R13 = tz * b_g__c_c__R13;
        // auto tz__b_g__c_c__R21 = tz * b_g__c_c__R21;
        // auto tz__b_g__c_c__R22 = tz * b_g__c_c__R22;
        // auto tz__b_g__c_c__R23 = tz * b_g__c_c__R23;
        // auto tz__b_g__c_c__R31 = tz * b_g__c_c__R31;
        // auto tz__b_g__c_c__R32 = tz * b_g__c_c__R32;
        // auto tz__b_g__c_c__R33 = tz * b_g__c_c__R33;
        // auto tz__b_g__x_c__R11 = tz * b_g__x_c__R11;
        // auto tz__b_g__x_c__R12 = tz * b_g__x_c__R12;
        // auto tz__b_g__x_c__R13 = tz * b_g__x_c__R13;
        // auto tz__b_g__x_c__R21 = tz * b_g__x_c__R21;
        // auto tz__b_g__x_c__R22 = tz * b_g__x_c__R22;
        // auto tz__b_g__x_c__R23 = tz * b_g__x_c__R23;
        // auto tz__b_g__x_c__R31 = tz * b_g__x_c__R31;
        // auto tz__b_g__x_c__R32 = tz * b_g__x_c__R32;
        // auto tz__b_g__x_c__R33 = tz * b_g__x_c__R33;
        // auto tz__b_g__y_c__R11 = tz * b_g__y_c__R11;
        // auto tz__b_g__y_c__R12 = tz * b_g__y_c__R12;
        // auto tz__b_g__y_c__R13 = tz * b_g__y_c__R13;
        // auto tz__b_g__y_c__R21 = tz * b_g__y_c__R21;
        // auto tz__b_g__y_c__R22 = tz * b_g__y_c__R22;
        // auto tz__b_g__y_c__R23 = tz * b_g__y_c__R23;
        // auto tz__b_g__y_c__R31 = tz * b_g__y_c__R31;
        // auto tz__b_g__y_c__R32 = tz * b_g__y_c__R32;
        // auto tz__b_g__y_c__R33 = tz * b_g__y_c__R33;
        // auto tz__b_g__z_c__R11 = tz * b_g__z_c__R11;
        // auto tz__b_g__z_c__R12 = tz * b_g__z_c__R12;
        // auto tz__b_g__z_c__R13 = tz * b_g__z_c__R13;
        // auto tz__b_g__z_c__R21 = tz * b_g__z_c__R21;
        // auto tz__b_g__z_c__R22 = tz * b_g__z_c__R22;
        // auto tz__b_g__z_c__R23 = tz * b_g__z_c__R23;
        // auto tz__b_g__z_c__R31 = tz * b_g__z_c__R31;
        // auto tz__b_g__z_c__R32 = tz * b_g__z_c__R32;
        // auto tz__b_g__z_c__R33 = tz * b_g__z_c__R33;
        // auto tz__c_g__a_c__R11 = tz * c_g__a_c__R11;
        // auto tz__c_g__a_c__R12 = tz * c_g__a_c__R12;
        // auto tz__c_g__a_c__R13 = tz * c_g__a_c__R13;
        // auto tz__c_g__a_c__R21 = tz * c_g__a_c__R21;
        // auto tz__c_g__a_c__R22 = tz * c_g__a_c__R22;
        // auto tz__c_g__a_c__R23 = tz * c_g__a_c__R23;
        // auto tz__c_g__a_c__R31 = tz * c_g__a_c__R31;
        // auto tz__c_g__a_c__R32 = tz * c_g__a_c__R32;
        // auto tz__c_g__a_c__R33 = tz * c_g__a_c__R33;
        // auto tz__c_g__b_c__R11 = tz * c_g__b_c__R11;
        // auto tz__c_g__b_c__R12 = tz * c_g__b_c__R12;
        // auto tz__c_g__b_c__R13 = tz * c_g__b_c__R13;
        // auto tz__c_g__b_c__R21 = tz * c_g__b_c__R21;
        // auto tz__c_g__b_c__R22 = tz * c_g__b_c__R22;
        // auto tz__c_g__b_c__R23 = tz * c_g__b_c__R23;
        // auto tz__c_g__b_c__R31 = tz * c_g__b_c__R31;
        // auto tz__c_g__b_c__R32 = tz * c_g__b_c__R32;
        // auto tz__c_g__b_c__R33 = tz * c_g__b_c__R33;
        // auto tz__c_g__c_c__R11 = tz * c_g__c_c__R11;
        // auto tz__c_g__c_c__R12 = tz * c_g__c_c__R12;
        // auto tz__c_g__c_c__R13 = tz * c_g__c_c__R13;
        // auto tz__c_g__c_c__R21 = tz * c_g__c_c__R21;
        // auto tz__c_g__c_c__R22 = tz * c_g__c_c__R22;
        // auto tz__c_g__c_c__R23 = tz * c_g__c_c__R23;
        // auto tz__c_g__c_c__R31 = tz * c_g__c_c__R31;
        // auto tz__c_g__c_c__R32 = tz * c_g__c_c__R32;
        // auto tz__c_g__c_c__R33 = tz * c_g__c_c__R33;
        // auto tz__c_g__x_c__R11 = tz * c_g__x_c__R11;
        // auto tz__c_g__x_c__R12 = tz * c_g__x_c__R12;
        // auto tz__c_g__x_c__R13 = tz * c_g__x_c__R13;
        // auto tz__c_g__x_c__R21 = tz * c_g__x_c__R21;
        // auto tz__c_g__x_c__R22 = tz * c_g__x_c__R22;
        // auto tz__c_g__x_c__R23 = tz * c_g__x_c__R23;
        // auto tz__c_g__x_c__R31 = tz * c_g__x_c__R31;
        // auto tz__c_g__x_c__R32 = tz * c_g__x_c__R32;
        // auto tz__c_g__x_c__R33 = tz * c_g__x_c__R33;
        // auto tz__c_g__y_c__R11 = tz * c_g__y_c__R11;
        // auto tz__c_g__y_c__R12 = tz * c_g__y_c__R12;
        // auto tz__c_g__y_c__R13 = tz * c_g__y_c__R13;
        // auto tz__c_g__y_c__R21 = tz * c_g__y_c__R21;
        // auto tz__c_g__y_c__R22 = tz * c_g__y_c__R22;
        // auto tz__c_g__y_c__R23 = tz * c_g__y_c__R23;
        // auto tz__c_g__y_c__R31 = tz * c_g__y_c__R31;
        // auto tz__c_g__y_c__R32 = tz * c_g__y_c__R32;
        // auto tz__c_g__y_c__R33 = tz * c_g__y_c__R33;
        // auto tz__c_g__z_c__R11 = tz * c_g__z_c__R11;
        // auto tz__c_g__z_c__R12 = tz * c_g__z_c__R12;
        // auto tz__c_g__z_c__R13 = tz * c_g__z_c__R13;
        // auto tz__c_g__z_c__R21 = tz * c_g__z_c__R21;
        // auto tz__c_g__z_c__R22 = tz * c_g__z_c__R22;
        // auto tz__c_g__z_c__R23 = tz * c_g__z_c__R23;
        // auto tz__c_g__z_c__R31 = tz * c_g__z_c__R31;
        // auto tz__c_g__z_c__R32 = tz * c_g__z_c__R32;
        // auto tz__c_g__z_c__R33 = tz * c_g__z_c__R33;
        // auto tz__x_g__a_c__R11 = tz * x_g__a_c__R11;
        // auto tz__x_g__a_c__R12 = tz * x_g__a_c__R12;
        // auto tz__x_g__a_c__R13 = tz * x_g__a_c__R13;
        // auto tz__x_g__a_c__R21 = tz * x_g__a_c__R21;
        // auto tz__x_g__a_c__R22 = tz * x_g__a_c__R22;
        // auto tz__x_g__a_c__R23 = tz * x_g__a_c__R23;
        // auto tz__x_g__a_c__R31 = tz * x_g__a_c__R31;
        // auto tz__x_g__a_c__R32 = tz * x_g__a_c__R32;
        // auto tz__x_g__a_c__R33 = tz * x_g__a_c__R33;
        // auto tz__x_g__b_c__R11 = tz * x_g__b_c__R11;
        // auto tz__x_g__b_c__R12 = tz * x_g__b_c__R12;
        // auto tz__x_g__b_c__R13 = tz * x_g__b_c__R13;
        // auto tz__x_g__b_c__R21 = tz * x_g__b_c__R21;
        // auto tz__x_g__b_c__R22 = tz * x_g__b_c__R22;
        // auto tz__x_g__b_c__R23 = tz * x_g__b_c__R23;
        // auto tz__x_g__b_c__R31 = tz * x_g__b_c__R31;
        // auto tz__x_g__b_c__R32 = tz * x_g__b_c__R32;
        // auto tz__x_g__b_c__R33 = tz * x_g__b_c__R33;
        // auto tz__x_g__c_c__R11 = tz * x_g__c_c__R11;
        // auto tz__x_g__c_c__R12 = tz * x_g__c_c__R12;
        // auto tz__x_g__c_c__R13 = tz * x_g__c_c__R13;
        // auto tz__x_g__c_c__R21 = tz * x_g__c_c__R21;
        // auto tz__x_g__c_c__R22 = tz * x_g__c_c__R22;
        // auto tz__x_g__c_c__R23 = tz * x_g__c_c__R23;
        // auto tz__x_g__c_c__R31 = tz * x_g__c_c__R31;
        // auto tz__x_g__c_c__R32 = tz * x_g__c_c__R32;
        // auto tz__x_g__c_c__R33 = tz * x_g__c_c__R33;
        // auto tz__x_g__x_c__R11 = tz * x_g__x_c__R11;
        // auto tz__x_g__x_c__R12 = tz * x_g__x_c__R12;
        // auto tz__x_g__x_c__R13 = tz * x_g__x_c__R13;
        // auto tz__x_g__x_c__R21 = tz * x_g__x_c__R21;
        // auto tz__x_g__x_c__R22 = tz * x_g__x_c__R22;
        // auto tz__x_g__x_c__R23 = tz * x_g__x_c__R23;
        // auto tz__x_g__x_c__R31 = tz * x_g__x_c__R31;
        // auto tz__x_g__x_c__R32 = tz * x_g__x_c__R32;
        // auto tz__x_g__x_c__R33 = tz * x_g__x_c__R33;
        // auto tz__x_g__y_c__R11 = tz * x_g__y_c__R11;
        // auto tz__x_g__y_c__R12 = tz * x_g__y_c__R12;
        // auto tz__x_g__y_c__R13 = tz * x_g__y_c__R13;
        // auto tz__x_g__y_c__R21 = tz * x_g__y_c__R21;
        // auto tz__x_g__y_c__R22 = tz * x_g__y_c__R22;
        // auto tz__x_g__y_c__R23 = tz * x_g__y_c__R23;
        // auto tz__x_g__y_c__R31 = tz * x_g__y_c__R31;
        // auto tz__x_g__y_c__R32 = tz * x_g__y_c__R32;
        // auto tz__x_g__y_c__R33 = tz * x_g__y_c__R33;
        // auto tz__x_g__z_c__R11 = tz * x_g__z_c__R11;
        // auto tz__x_g__z_c__R12 = tz * x_g__z_c__R12;
        // auto tz__x_g__z_c__R13 = tz * x_g__z_c__R13;
        // auto tz__x_g__z_c__R21 = tz * x_g__z_c__R21;
        // auto tz__x_g__z_c__R22 = tz * x_g__z_c__R22;
        // auto tz__x_g__z_c__R23 = tz * x_g__z_c__R23;
        // auto tz__x_g__z_c__R31 = tz * x_g__z_c__R31;
        // auto tz__x_g__z_c__R32 = tz * x_g__z_c__R32;
        // auto tz__x_g__z_c__R33 = tz * x_g__z_c__R33;
        // auto tz__y_g__a_c__R11 = tz * y_g__a_c__R11;
        // auto tz__y_g__a_c__R12 = tz * y_g__a_c__R12;
        // auto tz__y_g__a_c__R13 = tz * y_g__a_c__R13;
        // auto tz__y_g__a_c__R21 = tz * y_g__a_c__R21;
        // auto tz__y_g__a_c__R22 = tz * y_g__a_c__R22;
        // auto tz__y_g__a_c__R23 = tz * y_g__a_c__R23;
        // auto tz__y_g__a_c__R31 = tz * y_g__a_c__R31;
        // auto tz__y_g__a_c__R32 = tz * y_g__a_c__R32;
        // auto tz__y_g__a_c__R33 = tz * y_g__a_c__R33;
        // auto tz__y_g__b_c__R11 = tz * y_g__b_c__R11;
        // auto tz__y_g__b_c__R12 = tz * y_g__b_c__R12;
        // auto tz__y_g__b_c__R13 = tz * y_g__b_c__R13;
        // auto tz__y_g__b_c__R21 = tz * y_g__b_c__R21;
        // auto tz__y_g__b_c__R22 = tz * y_g__b_c__R22;
        // auto tz__y_g__b_c__R23 = tz * y_g__b_c__R23;
        // auto tz__y_g__b_c__R31 = tz * y_g__b_c__R31;
        // auto tz__y_g__b_c__R32 = tz * y_g__b_c__R32;
        // auto tz__y_g__b_c__R33 = tz * y_g__b_c__R33;
        // auto tz__y_g__c_c__R11 = tz * y_g__c_c__R11;
        // auto tz__y_g__c_c__R12 = tz * y_g__c_c__R12;
        // auto tz__y_g__c_c__R13 = tz * y_g__c_c__R13;
        // auto tz__y_g__c_c__R21 = tz * y_g__c_c__R21;
        // auto tz__y_g__c_c__R22 = tz * y_g__c_c__R22;
        // auto tz__y_g__c_c__R23 = tz * y_g__c_c__R23;
        // auto tz__y_g__c_c__R31 = tz * y_g__c_c__R31;
        // auto tz__y_g__c_c__R32 = tz * y_g__c_c__R32;
        // auto tz__y_g__c_c__R33 = tz * y_g__c_c__R33;
        // auto tz__y_g__x_c__R11 = tz * y_g__x_c__R11;
        // auto tz__y_g__x_c__R12 = tz * y_g__x_c__R12;
        // auto tz__y_g__x_c__R13 = tz * y_g__x_c__R13;
        // auto tz__y_g__x_c__R21 = tz * y_g__x_c__R21;
        // auto tz__y_g__x_c__R22 = tz * y_g__x_c__R22;
        // auto tz__y_g__x_c__R23 = tz * y_g__x_c__R23;
        // auto tz__y_g__x_c__R31 = tz * y_g__x_c__R31;
        // auto tz__y_g__x_c__R32 = tz * y_g__x_c__R32;
        // auto tz__y_g__x_c__R33 = tz * y_g__x_c__R33;
        // auto tz__y_g__y_c__R11 = tz * y_g__y_c__R11;
        // auto tz__y_g__y_c__R12 = tz * y_g__y_c__R12;
        // auto tz__y_g__y_c__R13 = tz * y_g__y_c__R13;
        // auto tz__y_g__y_c__R21 = tz * y_g__y_c__R21;
        // auto tz__y_g__y_c__R22 = tz * y_g__y_c__R22;
        // auto tz__y_g__y_c__R23 = tz * y_g__y_c__R23;
        // auto tz__y_g__y_c__R31 = tz * y_g__y_c__R31;
        // auto tz__y_g__y_c__R32 = tz * y_g__y_c__R32;
        // auto tz__y_g__y_c__R33 = tz * y_g__y_c__R33;
        // auto tz__y_g__z_c__R11 = tz * y_g__z_c__R11;
        // auto tz__y_g__z_c__R12 = tz * y_g__z_c__R12;
        // auto tz__y_g__z_c__R13 = tz * y_g__z_c__R13;
        // auto tz__y_g__z_c__R21 = tz * y_g__z_c__R21;
        // auto tz__y_g__z_c__R22 = tz * y_g__z_c__R22;
        // auto tz__y_g__z_c__R23 = tz * y_g__z_c__R23;
        // auto tz__y_g__z_c__R31 = tz * y_g__z_c__R31;
        // auto tz__y_g__z_c__R32 = tz * y_g__z_c__R32;
        // auto tz__y_g__z_c__R33 = tz * y_g__z_c__R33;
        // auto tz__z_g__a_c__R11 = tz * z_g__a_c__R11;
        // auto tz__z_g__a_c__R12 = tz * z_g__a_c__R12;
        // auto tz__z_g__a_c__R13 = tz * z_g__a_c__R13;
        // auto tz__z_g__a_c__R21 = tz * z_g__a_c__R21;
        // auto tz__z_g__a_c__R22 = tz * z_g__a_c__R22;
        // auto tz__z_g__a_c__R23 = tz * z_g__a_c__R23;
        // auto tz__z_g__a_c__R31 = tz * z_g__a_c__R31;
        // auto tz__z_g__a_c__R32 = tz * z_g__a_c__R32;
        // auto tz__z_g__a_c__R33 = tz * z_g__a_c__R33;
        // auto tz__z_g__b_c__R11 = tz * z_g__b_c__R11;
        // auto tz__z_g__b_c__R12 = tz * z_g__b_c__R12;
        // auto tz__z_g__b_c__R13 = tz * z_g__b_c__R13;
        // auto tz__z_g__b_c__R21 = tz * z_g__b_c__R21;
        // auto tz__z_g__b_c__R22 = tz * z_g__b_c__R22;
        // auto tz__z_g__b_c__R23 = tz * z_g__b_c__R23;
        // auto tz__z_g__b_c__R31 = tz * z_g__b_c__R31;
        // auto tz__z_g__b_c__R32 = tz * z_g__b_c__R32;
        // auto tz__z_g__b_c__R33 = tz * z_g__b_c__R33;
        // auto tz__z_g__c_c__R11 = tz * z_g__c_c__R11;
        // auto tz__z_g__c_c__R12 = tz * z_g__c_c__R12;
        // auto tz__z_g__c_c__R13 = tz * z_g__c_c__R13;
        // auto tz__z_g__c_c__R21 = tz * z_g__c_c__R21;
        // auto tz__z_g__c_c__R22 = tz * z_g__c_c__R22;
        // auto tz__z_g__c_c__R23 = tz * z_g__c_c__R23;
        // auto tz__z_g__c_c__R31 = tz * z_g__c_c__R31;
        // auto tz__z_g__c_c__R32 = tz * z_g__c_c__R32;
        // auto tz__z_g__c_c__R33 = tz * z_g__c_c__R33;
        // auto tz__z_g__x_c__R11 = tz * z_g__x_c__R11;
        // auto tz__z_g__x_c__R12 = tz * z_g__x_c__R12;
        // auto tz__z_g__x_c__R13 = tz * z_g__x_c__R13;
        // auto tz__z_g__x_c__R21 = tz * z_g__x_c__R21;
        // auto tz__z_g__x_c__R22 = tz * z_g__x_c__R22;
        // auto tz__z_g__x_c__R23 = tz * z_g__x_c__R23;
        // auto tz__z_g__x_c__R31 = tz * z_g__x_c__R31;
        // auto tz__z_g__x_c__R32 = tz * z_g__x_c__R32;
        // auto tz__z_g__x_c__R33 = tz * z_g__x_c__R33;
        // auto tz__z_g__y_c__R11 = tz * z_g__y_c__R11;
        // auto tz__z_g__y_c__R12 = tz * z_g__y_c__R12;
        // auto tz__z_g__y_c__R13 = tz * z_g__y_c__R13;
        // auto tz__z_g__y_c__R21 = tz * z_g__y_c__R21;
        // auto tz__z_g__y_c__R22 = tz * z_g__y_c__R22;
        // auto tz__z_g__y_c__R23 = tz * z_g__y_c__R23;
        // auto tz__z_g__y_c__R31 = tz * z_g__y_c__R31;
        // auto tz__z_g__y_c__R32 = tz * z_g__y_c__R32;
        // auto tz__z_g__y_c__R33 = tz * z_g__y_c__R33;
        // auto tz__z_g__z_c__R11 = tz * z_g__z_c__R11;
        // auto tz__z_g__z_c__R12 = tz * z_g__z_c__R12;
        // auto tz__z_g__z_c__R13 = tz * z_g__z_c__R13;
        // auto tz__z_g__z_c__R21 = tz * z_g__z_c__R21;
        // auto tz__z_g__z_c__R22 = tz * z_g__z_c__R22;
        // auto tz__z_g__z_c__R23 = tz * z_g__z_c__R23;
        // auto tz__z_g__z_c__R31 = tz * z_g__z_c__R31;
        // auto tz__z_g__z_c__R32 = tz * z_g__z_c__R32;
        // auto tz__z_g__z_c__R33 = tz * z_g__z_c__R33;
        // auto bx__a_g__a_c__R11 = bx * a_g__a_c__R11;
        // auto bx__a_g__a_c__R12 = bx * a_g__a_c__R12;
        // auto bx__a_g__a_c__R13 = bx * a_g__a_c__R13;
        // auto bx__a_g__a_c__R21 = bx * a_g__a_c__R21;
        // auto bx__a_g__a_c__R22 = bx * a_g__a_c__R22;
        // auto bx__a_g__a_c__R23 = bx * a_g__a_c__R23;
        // auto bx__a_g__a_c__R31 = bx * a_g__a_c__R31;
        // auto bx__a_g__a_c__R32 = bx * a_g__a_c__R32;
        // auto bx__a_g__a_c__R33 = bx * a_g__a_c__R33;
        // auto bx__a_g__b_c__R11 = bx * a_g__b_c__R11;
        // auto bx__a_g__b_c__R12 = bx * a_g__b_c__R12;
        // auto bx__a_g__b_c__R13 = bx * a_g__b_c__R13;
        // auto bx__a_g__b_c__R21 = bx * a_g__b_c__R21;
        // auto bx__a_g__b_c__R22 = bx * a_g__b_c__R22;
        // auto bx__a_g__b_c__R23 = bx * a_g__b_c__R23;
        // auto bx__a_g__b_c__R31 = bx * a_g__b_c__R31;
        // auto bx__a_g__b_c__R32 = bx * a_g__b_c__R32;
        // auto bx__a_g__b_c__R33 = bx * a_g__b_c__R33;
        // auto bx__a_g__c_c__R11 = bx * a_g__c_c__R11;
        // auto bx__a_g__c_c__R12 = bx * a_g__c_c__R12;
        // auto bx__a_g__c_c__R13 = bx * a_g__c_c__R13;
        // auto bx__a_g__c_c__R21 = bx * a_g__c_c__R21;
        // auto bx__a_g__c_c__R22 = bx * a_g__c_c__R22;
        // auto bx__a_g__c_c__R23 = bx * a_g__c_c__R23;
        // auto bx__a_g__c_c__R31 = bx * a_g__c_c__R31;
        // auto bx__a_g__c_c__R32 = bx * a_g__c_c__R32;
        // auto bx__a_g__c_c__R33 = bx * a_g__c_c__R33;
        // auto bx__a_g__x_c__R11 = bx * a_g__x_c__R11;
        // auto bx__a_g__x_c__R12 = bx * a_g__x_c__R12;
        // auto bx__a_g__x_c__R13 = bx * a_g__x_c__R13;
        // auto bx__a_g__x_c__R21 = bx * a_g__x_c__R21;
        // auto bx__a_g__x_c__R22 = bx * a_g__x_c__R22;
        // auto bx__a_g__x_c__R23 = bx * a_g__x_c__R23;
        // auto bx__a_g__x_c__R31 = bx * a_g__x_c__R31;
        // auto bx__a_g__x_c__R32 = bx * a_g__x_c__R32;
        // auto bx__a_g__x_c__R33 = bx * a_g__x_c__R33;
        // auto bx__a_g__y_c__R11 = bx * a_g__y_c__R11;
        // auto bx__a_g__y_c__R12 = bx * a_g__y_c__R12;
        // auto bx__a_g__y_c__R13 = bx * a_g__y_c__R13;
        // auto bx__a_g__y_c__R21 = bx * a_g__y_c__R21;
        // auto bx__a_g__y_c__R22 = bx * a_g__y_c__R22;
        // auto bx__a_g__y_c__R23 = bx * a_g__y_c__R23;
        // auto bx__a_g__y_c__R31 = bx * a_g__y_c__R31;
        // auto bx__a_g__y_c__R32 = bx * a_g__y_c__R32;
        // auto bx__a_g__y_c__R33 = bx * a_g__y_c__R33;
        // auto bx__a_g__z_c__R11 = bx * a_g__z_c__R11;
        // auto bx__a_g__z_c__R12 = bx * a_g__z_c__R12;
        // auto bx__a_g__z_c__R13 = bx * a_g__z_c__R13;
        // auto bx__a_g__z_c__R21 = bx * a_g__z_c__R21;
        // auto bx__a_g__z_c__R22 = bx * a_g__z_c__R22;
        // auto bx__a_g__z_c__R23 = bx * a_g__z_c__R23;
        // auto bx__a_g__z_c__R31 = bx * a_g__z_c__R31;
        // auto bx__a_g__z_c__R32 = bx * a_g__z_c__R32;
        // auto bx__a_g__z_c__R33 = bx * a_g__z_c__R33;
        // auto bx__b_g__a_c__R11 = bx * b_g__a_c__R11;
        // auto bx__b_g__a_c__R12 = bx * b_g__a_c__R12;
        // auto bx__b_g__a_c__R13 = bx * b_g__a_c__R13;
        // auto bx__b_g__a_c__R21 = bx * b_g__a_c__R21;
        // auto bx__b_g__a_c__R22 = bx * b_g__a_c__R22;
        // auto bx__b_g__a_c__R23 = bx * b_g__a_c__R23;
        // auto bx__b_g__a_c__R31 = bx * b_g__a_c__R31;
        // auto bx__b_g__a_c__R32 = bx * b_g__a_c__R32;
        // auto bx__b_g__a_c__R33 = bx * b_g__a_c__R33;
        // auto bx__b_g__b_c__R11 = bx * b_g__b_c__R11;
        // auto bx__b_g__b_c__R12 = bx * b_g__b_c__R12;
        // auto bx__b_g__b_c__R13 = bx * b_g__b_c__R13;
        // auto bx__b_g__b_c__R21 = bx * b_g__b_c__R21;
        // auto bx__b_g__b_c__R22 = bx * b_g__b_c__R22;
        // auto bx__b_g__b_c__R23 = bx * b_g__b_c__R23;
        // auto bx__b_g__b_c__R31 = bx * b_g__b_c__R31;
        // auto bx__b_g__b_c__R32 = bx * b_g__b_c__R32;
        // auto bx__b_g__b_c__R33 = bx * b_g__b_c__R33;
        // auto bx__b_g__c_c__R11 = bx * b_g__c_c__R11;
        // auto bx__b_g__c_c__R12 = bx * b_g__c_c__R12;
        // auto bx__b_g__c_c__R13 = bx * b_g__c_c__R13;
        // auto bx__b_g__c_c__R21 = bx * b_g__c_c__R21;
        // auto bx__b_g__c_c__R22 = bx * b_g__c_c__R22;
        // auto bx__b_g__c_c__R23 = bx * b_g__c_c__R23;
        // auto bx__b_g__c_c__R31 = bx * b_g__c_c__R31;
        // auto bx__b_g__c_c__R32 = bx * b_g__c_c__R32;
        // auto bx__b_g__c_c__R33 = bx * b_g__c_c__R33;
        // auto bx__b_g__x_c__R11 = bx * b_g__x_c__R11;
        // auto bx__b_g__x_c__R12 = bx * b_g__x_c__R12;
        // auto bx__b_g__x_c__R13 = bx * b_g__x_c__R13;
        // auto bx__b_g__x_c__R21 = bx * b_g__x_c__R21;
        // auto bx__b_g__x_c__R22 = bx * b_g__x_c__R22;
        // auto bx__b_g__x_c__R23 = bx * b_g__x_c__R23;
        // auto bx__b_g__x_c__R31 = bx * b_g__x_c__R31;
        // auto bx__b_g__x_c__R32 = bx * b_g__x_c__R32;
        // auto bx__b_g__x_c__R33 = bx * b_g__x_c__R33;
        // auto bx__b_g__y_c__R11 = bx * b_g__y_c__R11;
        // auto bx__b_g__y_c__R12 = bx * b_g__y_c__R12;
        // auto bx__b_g__y_c__R13 = bx * b_g__y_c__R13;
        // auto bx__b_g__y_c__R21 = bx * b_g__y_c__R21;
        // auto bx__b_g__y_c__R22 = bx * b_g__y_c__R22;
        // auto bx__b_g__y_c__R23 = bx * b_g__y_c__R23;
        // auto bx__b_g__y_c__R31 = bx * b_g__y_c__R31;
        // auto bx__b_g__y_c__R32 = bx * b_g__y_c__R32;
        // auto bx__b_g__y_c__R33 = bx * b_g__y_c__R33;
        // auto bx__b_g__z_c__R11 = bx * b_g__z_c__R11;
        // auto bx__b_g__z_c__R12 = bx * b_g__z_c__R12;
        // auto bx__b_g__z_c__R13 = bx * b_g__z_c__R13;
        // auto bx__b_g__z_c__R21 = bx * b_g__z_c__R21;
        // auto bx__b_g__z_c__R22 = bx * b_g__z_c__R22;
        // auto bx__b_g__z_c__R23 = bx * b_g__z_c__R23;
        // auto bx__b_g__z_c__R31 = bx * b_g__z_c__R31;
        // auto bx__b_g__z_c__R32 = bx * b_g__z_c__R32;
        // auto bx__b_g__z_c__R33 = bx * b_g__z_c__R33;
        // auto bx__c_g__a_c__R11 = bx * c_g__a_c__R11;
        // auto bx__c_g__a_c__R12 = bx * c_g__a_c__R12;
        // auto bx__c_g__a_c__R13 = bx * c_g__a_c__R13;
        // auto bx__c_g__a_c__R21 = bx * c_g__a_c__R21;
        // auto bx__c_g__a_c__R22 = bx * c_g__a_c__R22;
        // auto bx__c_g__a_c__R23 = bx * c_g__a_c__R23;
        // auto bx__c_g__a_c__R31 = bx * c_g__a_c__R31;
        // auto bx__c_g__a_c__R32 = bx * c_g__a_c__R32;
        // auto bx__c_g__a_c__R33 = bx * c_g__a_c__R33;
        // auto bx__c_g__b_c__R11 = bx * c_g__b_c__R11;
        // auto bx__c_g__b_c__R12 = bx * c_g__b_c__R12;
        // auto bx__c_g__b_c__R13 = bx * c_g__b_c__R13;
        // auto bx__c_g__b_c__R21 = bx * c_g__b_c__R21;
        // auto bx__c_g__b_c__R22 = bx * c_g__b_c__R22;
        // auto bx__c_g__b_c__R23 = bx * c_g__b_c__R23;
        // auto bx__c_g__b_c__R31 = bx * c_g__b_c__R31;
        // auto bx__c_g__b_c__R32 = bx * c_g__b_c__R32;
        // auto bx__c_g__b_c__R33 = bx * c_g__b_c__R33;
        // auto bx__c_g__c_c__R11 = bx * c_g__c_c__R11;
        // auto bx__c_g__c_c__R12 = bx * c_g__c_c__R12;
        // auto bx__c_g__c_c__R13 = bx * c_g__c_c__R13;
        // auto bx__c_g__c_c__R21 = bx * c_g__c_c__R21;
        // auto bx__c_g__c_c__R22 = bx * c_g__c_c__R22;
        // auto bx__c_g__c_c__R23 = bx * c_g__c_c__R23;
        // auto bx__c_g__c_c__R31 = bx * c_g__c_c__R31;
        // auto bx__c_g__c_c__R32 = bx * c_g__c_c__R32;
        // auto bx__c_g__c_c__R33 = bx * c_g__c_c__R33;
        // auto bx__c_g__x_c__R11 = bx * c_g__x_c__R11;
        // auto bx__c_g__x_c__R12 = bx * c_g__x_c__R12;
        // auto bx__c_g__x_c__R13 = bx * c_g__x_c__R13;
        // auto bx__c_g__x_c__R21 = bx * c_g__x_c__R21;
        // auto bx__c_g__x_c__R22 = bx * c_g__x_c__R22;
        // auto bx__c_g__x_c__R23 = bx * c_g__x_c__R23;
        // auto bx__c_g__x_c__R31 = bx * c_g__x_c__R31;
        // auto bx__c_g__x_c__R32 = bx * c_g__x_c__R32;
        // auto bx__c_g__x_c__R33 = bx * c_g__x_c__R33;
        // auto bx__c_g__y_c__R11 = bx * c_g__y_c__R11;
        // auto bx__c_g__y_c__R12 = bx * c_g__y_c__R12;
        // auto bx__c_g__y_c__R13 = bx * c_g__y_c__R13;
        // auto bx__c_g__y_c__R21 = bx * c_g__y_c__R21;
        // auto bx__c_g__y_c__R22 = bx * c_g__y_c__R22;
        // auto bx__c_g__y_c__R23 = bx * c_g__y_c__R23;
        // auto bx__c_g__y_c__R31 = bx * c_g__y_c__R31;
        // auto bx__c_g__y_c__R32 = bx * c_g__y_c__R32;
        // auto bx__c_g__y_c__R33 = bx * c_g__y_c__R33;
        // auto bx__c_g__z_c__R11 = bx * c_g__z_c__R11;
        // auto bx__c_g__z_c__R12 = bx * c_g__z_c__R12;
        // auto bx__c_g__z_c__R13 = bx * c_g__z_c__R13;
        // auto bx__c_g__z_c__R21 = bx * c_g__z_c__R21;
        // auto bx__c_g__z_c__R22 = bx * c_g__z_c__R22;
        // auto bx__c_g__z_c__R23 = bx * c_g__z_c__R23;
        // auto bx__c_g__z_c__R31 = bx * c_g__z_c__R31;
        // auto bx__c_g__z_c__R32 = bx * c_g__z_c__R32;
        // auto bx__c_g__z_c__R33 = bx * c_g__z_c__R33;
        // auto bx__x_g__a_c__R11 = bx * x_g__a_c__R11;
        // auto bx__x_g__a_c__R12 = bx * x_g__a_c__R12;
        // auto bx__x_g__a_c__R13 = bx * x_g__a_c__R13;
        // auto bx__x_g__a_c__R21 = bx * x_g__a_c__R21;
        // auto bx__x_g__a_c__R22 = bx * x_g__a_c__R22;
        // auto bx__x_g__a_c__R23 = bx * x_g__a_c__R23;
        // auto bx__x_g__a_c__R31 = bx * x_g__a_c__R31;
        // auto bx__x_g__a_c__R32 = bx * x_g__a_c__R32;
        // auto bx__x_g__a_c__R33 = bx * x_g__a_c__R33;
        // auto bx__x_g__b_c__R11 = bx * x_g__b_c__R11;
        // auto bx__x_g__b_c__R12 = bx * x_g__b_c__R12;
        // auto bx__x_g__b_c__R13 = bx * x_g__b_c__R13;
        // auto bx__x_g__b_c__R21 = bx * x_g__b_c__R21;
        // auto bx__x_g__b_c__R22 = bx * x_g__b_c__R22;
        // auto bx__x_g__b_c__R23 = bx * x_g__b_c__R23;
        // auto bx__x_g__b_c__R31 = bx * x_g__b_c__R31;
        // auto bx__x_g__b_c__R32 = bx * x_g__b_c__R32;
        // auto bx__x_g__b_c__R33 = bx * x_g__b_c__R33;
        // auto bx__x_g__c_c__R11 = bx * x_g__c_c__R11;
        // auto bx__x_g__c_c__R12 = bx * x_g__c_c__R12;
        // auto bx__x_g__c_c__R13 = bx * x_g__c_c__R13;
        // auto bx__x_g__c_c__R21 = bx * x_g__c_c__R21;
        // auto bx__x_g__c_c__R22 = bx * x_g__c_c__R22;
        // auto bx__x_g__c_c__R23 = bx * x_g__c_c__R23;
        // auto bx__x_g__c_c__R31 = bx * x_g__c_c__R31;
        // auto bx__x_g__c_c__R32 = bx * x_g__c_c__R32;
        // auto bx__x_g__c_c__R33 = bx * x_g__c_c__R33;
        // auto bx__x_g__x_c__R11 = bx * x_g__x_c__R11;
        // auto bx__x_g__x_c__R12 = bx * x_g__x_c__R12;
        // auto bx__x_g__x_c__R13 = bx * x_g__x_c__R13;
        // auto bx__x_g__x_c__R21 = bx * x_g__x_c__R21;
        // auto bx__x_g__x_c__R22 = bx * x_g__x_c__R22;
        // auto bx__x_g__x_c__R23 = bx * x_g__x_c__R23;
        // auto bx__x_g__x_c__R31 = bx * x_g__x_c__R31;
        // auto bx__x_g__x_c__R32 = bx * x_g__x_c__R32;
        // auto bx__x_g__x_c__R33 = bx * x_g__x_c__R33;
        // auto bx__x_g__y_c__R11 = bx * x_g__y_c__R11;
        // auto bx__x_g__y_c__R12 = bx * x_g__y_c__R12;
        // auto bx__x_g__y_c__R13 = bx * x_g__y_c__R13;
        // auto bx__x_g__y_c__R21 = bx * x_g__y_c__R21;
        // auto bx__x_g__y_c__R22 = bx * x_g__y_c__R22;
        // auto bx__x_g__y_c__R23 = bx * x_g__y_c__R23;
        // auto bx__x_g__y_c__R31 = bx * x_g__y_c__R31;
        // auto bx__x_g__y_c__R32 = bx * x_g__y_c__R32;
        // auto bx__x_g__y_c__R33 = bx * x_g__y_c__R33;
        // auto bx__x_g__z_c__R11 = bx * x_g__z_c__R11;
        // auto bx__x_g__z_c__R12 = bx * x_g__z_c__R12;
        // auto bx__x_g__z_c__R13 = bx * x_g__z_c__R13;
        // auto bx__x_g__z_c__R21 = bx * x_g__z_c__R21;
        // auto bx__x_g__z_c__R22 = bx * x_g__z_c__R22;
        // auto bx__x_g__z_c__R23 = bx * x_g__z_c__R23;
        // auto bx__x_g__z_c__R31 = bx * x_g__z_c__R31;
        // auto bx__x_g__z_c__R32 = bx * x_g__z_c__R32;
        // auto bx__x_g__z_c__R33 = bx * x_g__z_c__R33;
        // auto bx__y_g__a_c__R11 = bx * y_g__a_c__R11;
        // auto bx__y_g__a_c__R12 = bx * y_g__a_c__R12;
        // auto bx__y_g__a_c__R13 = bx * y_g__a_c__R13;
        // auto bx__y_g__a_c__R21 = bx * y_g__a_c__R21;
        // auto bx__y_g__a_c__R22 = bx * y_g__a_c__R22;
        // auto bx__y_g__a_c__R23 = bx * y_g__a_c__R23;
        // auto bx__y_g__a_c__R31 = bx * y_g__a_c__R31;
        // auto bx__y_g__a_c__R32 = bx * y_g__a_c__R32;
        // auto bx__y_g__a_c__R33 = bx * y_g__a_c__R33;
        // auto bx__y_g__b_c__R11 = bx * y_g__b_c__R11;
        // auto bx__y_g__b_c__R12 = bx * y_g__b_c__R12;
        // auto bx__y_g__b_c__R13 = bx * y_g__b_c__R13;
        // auto bx__y_g__b_c__R21 = bx * y_g__b_c__R21;
        // auto bx__y_g__b_c__R22 = bx * y_g__b_c__R22;
        // auto bx__y_g__b_c__R23 = bx * y_g__b_c__R23;
        // auto bx__y_g__b_c__R31 = bx * y_g__b_c__R31;
        // auto bx__y_g__b_c__R32 = bx * y_g__b_c__R32;
        // auto bx__y_g__b_c__R33 = bx * y_g__b_c__R33;
        // auto bx__y_g__c_c__R11 = bx * y_g__c_c__R11;
        // auto bx__y_g__c_c__R12 = bx * y_g__c_c__R12;
        // auto bx__y_g__c_c__R13 = bx * y_g__c_c__R13;
        // auto bx__y_g__c_c__R21 = bx * y_g__c_c__R21;
        // auto bx__y_g__c_c__R22 = bx * y_g__c_c__R22;
        // auto bx__y_g__c_c__R23 = bx * y_g__c_c__R23;
        // auto bx__y_g__c_c__R31 = bx * y_g__c_c__R31;
        // auto bx__y_g__c_c__R32 = bx * y_g__c_c__R32;
        // auto bx__y_g__c_c__R33 = bx * y_g__c_c__R33;
        // auto bx__y_g__x_c__R11 = bx * y_g__x_c__R11;
        // auto bx__y_g__x_c__R12 = bx * y_g__x_c__R12;
        // auto bx__y_g__x_c__R13 = bx * y_g__x_c__R13;
        // auto bx__y_g__x_c__R21 = bx * y_g__x_c__R21;
        // auto bx__y_g__x_c__R22 = bx * y_g__x_c__R22;
        // auto bx__y_g__x_c__R23 = bx * y_g__x_c__R23;
        // auto bx__y_g__x_c__R31 = bx * y_g__x_c__R31;
        // auto bx__y_g__x_c__R32 = bx * y_g__x_c__R32;
        // auto bx__y_g__x_c__R33 = bx * y_g__x_c__R33;
        // auto bx__y_g__y_c__R11 = bx * y_g__y_c__R11;
        // auto bx__y_g__y_c__R12 = bx * y_g__y_c__R12;
        // auto bx__y_g__y_c__R13 = bx * y_g__y_c__R13;
        // auto bx__y_g__y_c__R21 = bx * y_g__y_c__R21;
        // auto bx__y_g__y_c__R22 = bx * y_g__y_c__R22;
        // auto bx__y_g__y_c__R23 = bx * y_g__y_c__R23;
        // auto bx__y_g__y_c__R31 = bx * y_g__y_c__R31;
        // auto bx__y_g__y_c__R32 = bx * y_g__y_c__R32;
        // auto bx__y_g__y_c__R33 = bx * y_g__y_c__R33;
        // auto bx__y_g__z_c__R11 = bx * y_g__z_c__R11;
        // auto bx__y_g__z_c__R12 = bx * y_g__z_c__R12;
        // auto bx__y_g__z_c__R13 = bx * y_g__z_c__R13;
        // auto bx__y_g__z_c__R21 = bx * y_g__z_c__R21;
        // auto bx__y_g__z_c__R22 = bx * y_g__z_c__R22;
        // auto bx__y_g__z_c__R23 = bx * y_g__z_c__R23;
        // auto bx__y_g__z_c__R31 = bx * y_g__z_c__R31;
        // auto bx__y_g__z_c__R32 = bx * y_g__z_c__R32;
        // auto bx__y_g__z_c__R33 = bx * y_g__z_c__R33;
        // auto bx__z_g__a_c__R11 = bx * z_g__a_c__R11;
        // auto bx__z_g__a_c__R12 = bx * z_g__a_c__R12;
        // auto bx__z_g__a_c__R13 = bx * z_g__a_c__R13;
        // auto bx__z_g__a_c__R21 = bx * z_g__a_c__R21;
        // auto bx__z_g__a_c__R22 = bx * z_g__a_c__R22;
        // auto bx__z_g__a_c__R23 = bx * z_g__a_c__R23;
        // auto bx__z_g__a_c__R31 = bx * z_g__a_c__R31;
        // auto bx__z_g__a_c__R32 = bx * z_g__a_c__R32;
        // auto bx__z_g__a_c__R33 = bx * z_g__a_c__R33;
        // auto bx__z_g__b_c__R11 = bx * z_g__b_c__R11;
        // auto bx__z_g__b_c__R12 = bx * z_g__b_c__R12;
        // auto bx__z_g__b_c__R13 = bx * z_g__b_c__R13;
        // auto bx__z_g__b_c__R21 = bx * z_g__b_c__R21;
        // auto bx__z_g__b_c__R22 = bx * z_g__b_c__R22;
        // auto bx__z_g__b_c__R23 = bx * z_g__b_c__R23;
        // auto bx__z_g__b_c__R31 = bx * z_g__b_c__R31;
        // auto bx__z_g__b_c__R32 = bx * z_g__b_c__R32;
        // auto bx__z_g__b_c__R33 = bx * z_g__b_c__R33;
        // auto bx__z_g__c_c__R11 = bx * z_g__c_c__R11;
        // auto bx__z_g__c_c__R12 = bx * z_g__c_c__R12;
        // auto bx__z_g__c_c__R13 = bx * z_g__c_c__R13;
        // auto bx__z_g__c_c__R21 = bx * z_g__c_c__R21;
        // auto bx__z_g__c_c__R22 = bx * z_g__c_c__R22;
        // auto bx__z_g__c_c__R23 = bx * z_g__c_c__R23;
        // auto bx__z_g__c_c__R31 = bx * z_g__c_c__R31;
        // auto bx__z_g__c_c__R32 = bx * z_g__c_c__R32;
        // auto bx__z_g__c_c__R33 = bx * z_g__c_c__R33;
        // auto bx__z_g__x_c__R11 = bx * z_g__x_c__R11;
        // auto bx__z_g__x_c__R12 = bx * z_g__x_c__R12;
        // auto bx__z_g__x_c__R13 = bx * z_g__x_c__R13;
        // auto bx__z_g__x_c__R21 = bx * z_g__x_c__R21;
        // auto bx__z_g__x_c__R22 = bx * z_g__x_c__R22;
        // auto bx__z_g__x_c__R23 = bx * z_g__x_c__R23;
        // auto bx__z_g__x_c__R31 = bx * z_g__x_c__R31;
        // auto bx__z_g__x_c__R32 = bx * z_g__x_c__R32;
        // auto bx__z_g__x_c__R33 = bx * z_g__x_c__R33;
        // auto bx__z_g__y_c__R11 = bx * z_g__y_c__R11;
        // auto bx__z_g__y_c__R12 = bx * z_g__y_c__R12;
        // auto bx__z_g__y_c__R13 = bx * z_g__y_c__R13;
        // auto bx__z_g__y_c__R21 = bx * z_g__y_c__R21;
        // auto bx__z_g__y_c__R22 = bx * z_g__y_c__R22;
        // auto bx__z_g__y_c__R23 = bx * z_g__y_c__R23;
        // auto bx__z_g__y_c__R31 = bx * z_g__y_c__R31;
        // auto bx__z_g__y_c__R32 = bx * z_g__y_c__R32;
        // auto bx__z_g__y_c__R33 = bx * z_g__y_c__R33;
        // auto bx__z_g__z_c__R11 = bx * z_g__z_c__R11;
        // auto bx__z_g__z_c__R12 = bx * z_g__z_c__R12;
        // auto bx__z_g__z_c__R13 = bx * z_g__z_c__R13;
        // auto bx__z_g__z_c__R21 = bx * z_g__z_c__R21;
        // auto bx__z_g__z_c__R22 = bx * z_g__z_c__R22;
        // auto bx__z_g__z_c__R23 = bx * z_g__z_c__R23;
        // auto bx__z_g__z_c__R31 = bx * z_g__z_c__R31;
        // auto bx__z_g__z_c__R32 = bx * z_g__z_c__R32;
        // auto bx__z_g__z_c__R33 = bx * z_g__z_c__R33;
        // auto by__a_g__a_c__R11 = by * a_g__a_c__R11;
        // auto by__a_g__a_c__R12 = by * a_g__a_c__R12;
        // auto by__a_g__a_c__R13 = by * a_g__a_c__R13;
        // auto by__a_g__a_c__R21 = by * a_g__a_c__R21;
        // auto by__a_g__a_c__R22 = by * a_g__a_c__R22;
        // auto by__a_g__a_c__R23 = by * a_g__a_c__R23;
        // auto by__a_g__a_c__R31 = by * a_g__a_c__R31;
        // auto by__a_g__a_c__R32 = by * a_g__a_c__R32;
        // auto by__a_g__a_c__R33 = by * a_g__a_c__R33;
        // auto by__a_g__b_c__R11 = by * a_g__b_c__R11;
        // auto by__a_g__b_c__R12 = by * a_g__b_c__R12;
        // auto by__a_g__b_c__R13 = by * a_g__b_c__R13;
        // auto by__a_g__b_c__R21 = by * a_g__b_c__R21;
        // auto by__a_g__b_c__R22 = by * a_g__b_c__R22;
        // auto by__a_g__b_c__R23 = by * a_g__b_c__R23;
        // auto by__a_g__b_c__R31 = by * a_g__b_c__R31;
        // auto by__a_g__b_c__R32 = by * a_g__b_c__R32;
        // auto by__a_g__b_c__R33 = by * a_g__b_c__R33;
        // auto by__a_g__c_c__R11 = by * a_g__c_c__R11;
        // auto by__a_g__c_c__R12 = by * a_g__c_c__R12;
        // auto by__a_g__c_c__R13 = by * a_g__c_c__R13;
        // auto by__a_g__c_c__R21 = by * a_g__c_c__R21;
        // auto by__a_g__c_c__R22 = by * a_g__c_c__R22;
        // auto by__a_g__c_c__R23 = by * a_g__c_c__R23;
        // auto by__a_g__c_c__R31 = by * a_g__c_c__R31;
        // auto by__a_g__c_c__R32 = by * a_g__c_c__R32;
        // auto by__a_g__c_c__R33 = by * a_g__c_c__R33;
        // auto by__a_g__x_c__R11 = by * a_g__x_c__R11;
        // auto by__a_g__x_c__R12 = by * a_g__x_c__R12;
        // auto by__a_g__x_c__R13 = by * a_g__x_c__R13;
        // auto by__a_g__x_c__R21 = by * a_g__x_c__R21;
        // auto by__a_g__x_c__R22 = by * a_g__x_c__R22;
        // auto by__a_g__x_c__R23 = by * a_g__x_c__R23;
        // auto by__a_g__x_c__R31 = by * a_g__x_c__R31;
        // auto by__a_g__x_c__R32 = by * a_g__x_c__R32;
        // auto by__a_g__x_c__R33 = by * a_g__x_c__R33;
        // auto by__a_g__y_c__R11 = by * a_g__y_c__R11;
        // auto by__a_g__y_c__R12 = by * a_g__y_c__R12;
        // auto by__a_g__y_c__R13 = by * a_g__y_c__R13;
        // auto by__a_g__y_c__R21 = by * a_g__y_c__R21;
        // auto by__a_g__y_c__R22 = by * a_g__y_c__R22;
        // auto by__a_g__y_c__R23 = by * a_g__y_c__R23;
        // auto by__a_g__y_c__R31 = by * a_g__y_c__R31;
        // auto by__a_g__y_c__R32 = by * a_g__y_c__R32;
        // auto by__a_g__y_c__R33 = by * a_g__y_c__R33;
        // auto by__a_g__z_c__R11 = by * a_g__z_c__R11;
        // auto by__a_g__z_c__R12 = by * a_g__z_c__R12;
        // auto by__a_g__z_c__R13 = by * a_g__z_c__R13;
        // auto by__a_g__z_c__R21 = by * a_g__z_c__R21;
        // auto by__a_g__z_c__R22 = by * a_g__z_c__R22;
        // auto by__a_g__z_c__R23 = by * a_g__z_c__R23;
        // auto by__a_g__z_c__R31 = by * a_g__z_c__R31;
        // auto by__a_g__z_c__R32 = by * a_g__z_c__R32;
        // auto by__a_g__z_c__R33 = by * a_g__z_c__R33;
        // auto by__b_g__a_c__R11 = by * b_g__a_c__R11;
        // auto by__b_g__a_c__R12 = by * b_g__a_c__R12;
        // auto by__b_g__a_c__R13 = by * b_g__a_c__R13;
        // auto by__b_g__a_c__R21 = by * b_g__a_c__R21;
        // auto by__b_g__a_c__R22 = by * b_g__a_c__R22;
        // auto by__b_g__a_c__R23 = by * b_g__a_c__R23;
        // auto by__b_g__a_c__R31 = by * b_g__a_c__R31;
        // auto by__b_g__a_c__R32 = by * b_g__a_c__R32;
        // auto by__b_g__a_c__R33 = by * b_g__a_c__R33;
        // auto by__b_g__b_c__R11 = by * b_g__b_c__R11;
        // auto by__b_g__b_c__R12 = by * b_g__b_c__R12;
        // auto by__b_g__b_c__R13 = by * b_g__b_c__R13;
        // auto by__b_g__b_c__R21 = by * b_g__b_c__R21;
        // auto by__b_g__b_c__R22 = by * b_g__b_c__R22;
        // auto by__b_g__b_c__R23 = by * b_g__b_c__R23;
        // auto by__b_g__b_c__R31 = by * b_g__b_c__R31;
        // auto by__b_g__b_c__R32 = by * b_g__b_c__R32;
        // auto by__b_g__b_c__R33 = by * b_g__b_c__R33;
        // auto by__b_g__c_c__R11 = by * b_g__c_c__R11;
        // auto by__b_g__c_c__R12 = by * b_g__c_c__R12;
        // auto by__b_g__c_c__R13 = by * b_g__c_c__R13;
        // auto by__b_g__c_c__R21 = by * b_g__c_c__R21;
        // auto by__b_g__c_c__R22 = by * b_g__c_c__R22;
        // auto by__b_g__c_c__R23 = by * b_g__c_c__R23;
        // auto by__b_g__c_c__R31 = by * b_g__c_c__R31;
        // auto by__b_g__c_c__R32 = by * b_g__c_c__R32;
        // auto by__b_g__c_c__R33 = by * b_g__c_c__R33;
        // auto by__b_g__x_c__R11 = by * b_g__x_c__R11;
        // auto by__b_g__x_c__R12 = by * b_g__x_c__R12;
        // auto by__b_g__x_c__R13 = by * b_g__x_c__R13;
        // auto by__b_g__x_c__R21 = by * b_g__x_c__R21;
        // auto by__b_g__x_c__R22 = by * b_g__x_c__R22;
        // auto by__b_g__x_c__R23 = by * b_g__x_c__R23;
        // auto by__b_g__x_c__R31 = by * b_g__x_c__R31;
        // auto by__b_g__x_c__R32 = by * b_g__x_c__R32;
        // auto by__b_g__x_c__R33 = by * b_g__x_c__R33;
        // auto by__b_g__y_c__R11 = by * b_g__y_c__R11;
        // auto by__b_g__y_c__R12 = by * b_g__y_c__R12;
        // auto by__b_g__y_c__R13 = by * b_g__y_c__R13;
        // auto by__b_g__y_c__R21 = by * b_g__y_c__R21;
        // auto by__b_g__y_c__R22 = by * b_g__y_c__R22;
        // auto by__b_g__y_c__R23 = by * b_g__y_c__R23;
        // auto by__b_g__y_c__R31 = by * b_g__y_c__R31;
        // auto by__b_g__y_c__R32 = by * b_g__y_c__R32;
        // auto by__b_g__y_c__R33 = by * b_g__y_c__R33;
        // auto by__b_g__z_c__R11 = by * b_g__z_c__R11;
        // auto by__b_g__z_c__R12 = by * b_g__z_c__R12;
        // auto by__b_g__z_c__R13 = by * b_g__z_c__R13;
        // auto by__b_g__z_c__R21 = by * b_g__z_c__R21;
        // auto by__b_g__z_c__R22 = by * b_g__z_c__R22;
        // auto by__b_g__z_c__R23 = by * b_g__z_c__R23;
        // auto by__b_g__z_c__R31 = by * b_g__z_c__R31;
        // auto by__b_g__z_c__R32 = by * b_g__z_c__R32;
        // auto by__b_g__z_c__R33 = by * b_g__z_c__R33;
        // auto by__c_g__a_c__R11 = by * c_g__a_c__R11;
        // auto by__c_g__a_c__R12 = by * c_g__a_c__R12;
        // auto by__c_g__a_c__R13 = by * c_g__a_c__R13;
        // auto by__c_g__a_c__R21 = by * c_g__a_c__R21;
        // auto by__c_g__a_c__R22 = by * c_g__a_c__R22;
        // auto by__c_g__a_c__R23 = by * c_g__a_c__R23;
        // auto by__c_g__a_c__R31 = by * c_g__a_c__R31;
        // auto by__c_g__a_c__R32 = by * c_g__a_c__R32;
        // auto by__c_g__a_c__R33 = by * c_g__a_c__R33;
        // auto by__c_g__b_c__R11 = by * c_g__b_c__R11;
        // auto by__c_g__b_c__R12 = by * c_g__b_c__R12;
        // auto by__c_g__b_c__R13 = by * c_g__b_c__R13;
        // auto by__c_g__b_c__R21 = by * c_g__b_c__R21;
        // auto by__c_g__b_c__R22 = by * c_g__b_c__R22;
        // auto by__c_g__b_c__R23 = by * c_g__b_c__R23;
        // auto by__c_g__b_c__R31 = by * c_g__b_c__R31;
        // auto by__c_g__b_c__R32 = by * c_g__b_c__R32;
        // auto by__c_g__b_c__R33 = by * c_g__b_c__R33;
        // auto by__c_g__c_c__R11 = by * c_g__c_c__R11;
        // auto by__c_g__c_c__R12 = by * c_g__c_c__R12;
        // auto by__c_g__c_c__R13 = by * c_g__c_c__R13;
        // auto by__c_g__c_c__R21 = by * c_g__c_c__R21;
        // auto by__c_g__c_c__R22 = by * c_g__c_c__R22;
        // auto by__c_g__c_c__R23 = by * c_g__c_c__R23;
        // auto by__c_g__c_c__R31 = by * c_g__c_c__R31;
        // auto by__c_g__c_c__R32 = by * c_g__c_c__R32;
        // auto by__c_g__c_c__R33 = by * c_g__c_c__R33;
        // auto by__c_g__x_c__R11 = by * c_g__x_c__R11;
        // auto by__c_g__x_c__R12 = by * c_g__x_c__R12;
        // auto by__c_g__x_c__R13 = by * c_g__x_c__R13;
        // auto by__c_g__x_c__R21 = by * c_g__x_c__R21;
        // auto by__c_g__x_c__R22 = by * c_g__x_c__R22;
        // auto by__c_g__x_c__R23 = by * c_g__x_c__R23;
        // auto by__c_g__x_c__R31 = by * c_g__x_c__R31;
        // auto by__c_g__x_c__R32 = by * c_g__x_c__R32;
        // auto by__c_g__x_c__R33 = by * c_g__x_c__R33;
        // auto by__c_g__y_c__R11 = by * c_g__y_c__R11;
        // auto by__c_g__y_c__R12 = by * c_g__y_c__R12;
        // auto by__c_g__y_c__R13 = by * c_g__y_c__R13;
        // auto by__c_g__y_c__R21 = by * c_g__y_c__R21;
        // auto by__c_g__y_c__R22 = by * c_g__y_c__R22;
        // auto by__c_g__y_c__R23 = by * c_g__y_c__R23;
        // auto by__c_g__y_c__R31 = by * c_g__y_c__R31;
        // auto by__c_g__y_c__R32 = by * c_g__y_c__R32;
        // auto by__c_g__y_c__R33 = by * c_g__y_c__R33;
        // auto by__c_g__z_c__R11 = by * c_g__z_c__R11;
        // auto by__c_g__z_c__R12 = by * c_g__z_c__R12;
        // auto by__c_g__z_c__R13 = by * c_g__z_c__R13;
        // auto by__c_g__z_c__R21 = by * c_g__z_c__R21;
        // auto by__c_g__z_c__R22 = by * c_g__z_c__R22;
        // auto by__c_g__z_c__R23 = by * c_g__z_c__R23;
        // auto by__c_g__z_c__R31 = by * c_g__z_c__R31;
        // auto by__c_g__z_c__R32 = by * c_g__z_c__R32;
        // auto by__c_g__z_c__R33 = by * c_g__z_c__R33;
        // auto by__x_g__a_c__R11 = by * x_g__a_c__R11;
        // auto by__x_g__a_c__R12 = by * x_g__a_c__R12;
        // auto by__x_g__a_c__R13 = by * x_g__a_c__R13;
        // auto by__x_g__a_c__R21 = by * x_g__a_c__R21;
        // auto by__x_g__a_c__R22 = by * x_g__a_c__R22;
        // auto by__x_g__a_c__R23 = by * x_g__a_c__R23;
        // auto by__x_g__a_c__R31 = by * x_g__a_c__R31;
        // auto by__x_g__a_c__R32 = by * x_g__a_c__R32;
        // auto by__x_g__a_c__R33 = by * x_g__a_c__R33;
        // auto by__x_g__b_c__R11 = by * x_g__b_c__R11;
        // auto by__x_g__b_c__R12 = by * x_g__b_c__R12;
        // auto by__x_g__b_c__R13 = by * x_g__b_c__R13;
        // auto by__x_g__b_c__R21 = by * x_g__b_c__R21;
        // auto by__x_g__b_c__R22 = by * x_g__b_c__R22;
        // auto by__x_g__b_c__R23 = by * x_g__b_c__R23;
        // auto by__x_g__b_c__R31 = by * x_g__b_c__R31;
        // auto by__x_g__b_c__R32 = by * x_g__b_c__R32;
        // auto by__x_g__b_c__R33 = by * x_g__b_c__R33;
        // auto by__x_g__c_c__R11 = by * x_g__c_c__R11;
        // auto by__x_g__c_c__R12 = by * x_g__c_c__R12;
        // auto by__x_g__c_c__R13 = by * x_g__c_c__R13;
        // auto by__x_g__c_c__R21 = by * x_g__c_c__R21;
        // auto by__x_g__c_c__R22 = by * x_g__c_c__R22;
        // auto by__x_g__c_c__R23 = by * x_g__c_c__R23;
        // auto by__x_g__c_c__R31 = by * x_g__c_c__R31;
        // auto by__x_g__c_c__R32 = by * x_g__c_c__R32;
        // auto by__x_g__c_c__R33 = by * x_g__c_c__R33;
        // auto by__x_g__x_c__R11 = by * x_g__x_c__R11;
        // auto by__x_g__x_c__R12 = by * x_g__x_c__R12;
        // auto by__x_g__x_c__R13 = by * x_g__x_c__R13;
        // auto by__x_g__x_c__R21 = by * x_g__x_c__R21;
        // auto by__x_g__x_c__R22 = by * x_g__x_c__R22;
        // auto by__x_g__x_c__R23 = by * x_g__x_c__R23;
        // auto by__x_g__x_c__R31 = by * x_g__x_c__R31;
        // auto by__x_g__x_c__R32 = by * x_g__x_c__R32;
        // auto by__x_g__x_c__R33 = by * x_g__x_c__R33;
        // auto by__x_g__y_c__R11 = by * x_g__y_c__R11;
        // auto by__x_g__y_c__R12 = by * x_g__y_c__R12;
        // auto by__x_g__y_c__R13 = by * x_g__y_c__R13;
        // auto by__x_g__y_c__R21 = by * x_g__y_c__R21;
        // auto by__x_g__y_c__R22 = by * x_g__y_c__R22;
        // auto by__x_g__y_c__R23 = by * x_g__y_c__R23;
        // auto by__x_g__y_c__R31 = by * x_g__y_c__R31;
        // auto by__x_g__y_c__R32 = by * x_g__y_c__R32;
        // auto by__x_g__y_c__R33 = by * x_g__y_c__R33;
        // auto by__x_g__z_c__R11 = by * x_g__z_c__R11;
        // auto by__x_g__z_c__R12 = by * x_g__z_c__R12;
        // auto by__x_g__z_c__R13 = by * x_g__z_c__R13;
        // auto by__x_g__z_c__R21 = by * x_g__z_c__R21;
        // auto by__x_g__z_c__R22 = by * x_g__z_c__R22;
        // auto by__x_g__z_c__R23 = by * x_g__z_c__R23;
        // auto by__x_g__z_c__R31 = by * x_g__z_c__R31;
        // auto by__x_g__z_c__R32 = by * x_g__z_c__R32;
        // auto by__x_g__z_c__R33 = by * x_g__z_c__R33;
        // auto by__y_g__a_c__R11 = by * y_g__a_c__R11;
        // auto by__y_g__a_c__R12 = by * y_g__a_c__R12;
        // auto by__y_g__a_c__R13 = by * y_g__a_c__R13;
        // auto by__y_g__a_c__R21 = by * y_g__a_c__R21;
        // auto by__y_g__a_c__R22 = by * y_g__a_c__R22;
        // auto by__y_g__a_c__R23 = by * y_g__a_c__R23;
        // auto by__y_g__a_c__R31 = by * y_g__a_c__R31;
        // auto by__y_g__a_c__R32 = by * y_g__a_c__R32;
        // auto by__y_g__a_c__R33 = by * y_g__a_c__R33;
        // auto by__y_g__b_c__R11 = by * y_g__b_c__R11;
        // auto by__y_g__b_c__R12 = by * y_g__b_c__R12;
        // auto by__y_g__b_c__R13 = by * y_g__b_c__R13;
        // auto by__y_g__b_c__R21 = by * y_g__b_c__R21;
        // auto by__y_g__b_c__R22 = by * y_g__b_c__R22;
        // auto by__y_g__b_c__R23 = by * y_g__b_c__R23;
        // auto by__y_g__b_c__R31 = by * y_g__b_c__R31;
        // auto by__y_g__b_c__R32 = by * y_g__b_c__R32;
        // auto by__y_g__b_c__R33 = by * y_g__b_c__R33;
        // auto by__y_g__c_c__R11 = by * y_g__c_c__R11;
        // auto by__y_g__c_c__R12 = by * y_g__c_c__R12;
        // auto by__y_g__c_c__R13 = by * y_g__c_c__R13;
        // auto by__y_g__c_c__R21 = by * y_g__c_c__R21;
        // auto by__y_g__c_c__R22 = by * y_g__c_c__R22;
        // auto by__y_g__c_c__R23 = by * y_g__c_c__R23;
        // auto by__y_g__c_c__R31 = by * y_g__c_c__R31;
        // auto by__y_g__c_c__R32 = by * y_g__c_c__R32;
        // auto by__y_g__c_c__R33 = by * y_g__c_c__R33;
        // auto by__y_g__x_c__R11 = by * y_g__x_c__R11;
        // auto by__y_g__x_c__R12 = by * y_g__x_c__R12;
        // auto by__y_g__x_c__R13 = by * y_g__x_c__R13;
        // auto by__y_g__x_c__R21 = by * y_g__x_c__R21;
        // auto by__y_g__x_c__R22 = by * y_g__x_c__R22;
        // auto by__y_g__x_c__R23 = by * y_g__x_c__R23;
        // auto by__y_g__x_c__R31 = by * y_g__x_c__R31;
        // auto by__y_g__x_c__R32 = by * y_g__x_c__R32;
        // auto by__y_g__x_c__R33 = by * y_g__x_c__R33;
        // auto by__y_g__y_c__R11 = by * y_g__y_c__R11;
        // auto by__y_g__y_c__R12 = by * y_g__y_c__R12;
        // auto by__y_g__y_c__R13 = by * y_g__y_c__R13;
        // auto by__y_g__y_c__R21 = by * y_g__y_c__R21;
        // auto by__y_g__y_c__R22 = by * y_g__y_c__R22;
        // auto by__y_g__y_c__R23 = by * y_g__y_c__R23;
        // auto by__y_g__y_c__R31 = by * y_g__y_c__R31;
        // auto by__y_g__y_c__R32 = by * y_g__y_c__R32;
        // auto by__y_g__y_c__R33 = by * y_g__y_c__R33;
        // auto by__y_g__z_c__R11 = by * y_g__z_c__R11;
        // auto by__y_g__z_c__R12 = by * y_g__z_c__R12;
        // auto by__y_g__z_c__R13 = by * y_g__z_c__R13;
        // auto by__y_g__z_c__R21 = by * y_g__z_c__R21;
        // auto by__y_g__z_c__R22 = by * y_g__z_c__R22;
        // auto by__y_g__z_c__R23 = by * y_g__z_c__R23;
        // auto by__y_g__z_c__R31 = by * y_g__z_c__R31;
        // auto by__y_g__z_c__R32 = by * y_g__z_c__R32;
        // auto by__y_g__z_c__R33 = by * y_g__z_c__R33;
        // auto by__z_g__a_c__R11 = by * z_g__a_c__R11;
        // auto by__z_g__a_c__R12 = by * z_g__a_c__R12;
        // auto by__z_g__a_c__R13 = by * z_g__a_c__R13;
        // auto by__z_g__a_c__R21 = by * z_g__a_c__R21;
        // auto by__z_g__a_c__R22 = by * z_g__a_c__R22;
        // auto by__z_g__a_c__R23 = by * z_g__a_c__R23;
        // auto by__z_g__a_c__R31 = by * z_g__a_c__R31;
        // auto by__z_g__a_c__R32 = by * z_g__a_c__R32;
        // auto by__z_g__a_c__R33 = by * z_g__a_c__R33;
        // auto by__z_g__b_c__R11 = by * z_g__b_c__R11;
        // auto by__z_g__b_c__R12 = by * z_g__b_c__R12;
        // auto by__z_g__b_c__R13 = by * z_g__b_c__R13;
        // auto by__z_g__b_c__R21 = by * z_g__b_c__R21;
        // auto by__z_g__b_c__R22 = by * z_g__b_c__R22;
        // auto by__z_g__b_c__R23 = by * z_g__b_c__R23;
        // auto by__z_g__b_c__R31 = by * z_g__b_c__R31;
        // auto by__z_g__b_c__R32 = by * z_g__b_c__R32;
        // auto by__z_g__b_c__R33 = by * z_g__b_c__R33;
        // auto by__z_g__c_c__R11 = by * z_g__c_c__R11;
        // auto by__z_g__c_c__R12 = by * z_g__c_c__R12;
        // auto by__z_g__c_c__R13 = by * z_g__c_c__R13;
        // auto by__z_g__c_c__R21 = by * z_g__c_c__R21;
        // auto by__z_g__c_c__R22 = by * z_g__c_c__R22;
        // auto by__z_g__c_c__R23 = by * z_g__c_c__R23;
        // auto by__z_g__c_c__R31 = by * z_g__c_c__R31;
        // auto by__z_g__c_c__R32 = by * z_g__c_c__R32;
        // auto by__z_g__c_c__R33 = by * z_g__c_c__R33;
        // auto by__z_g__x_c__R11 = by * z_g__x_c__R11;
        // auto by__z_g__x_c__R12 = by * z_g__x_c__R12;
        // auto by__z_g__x_c__R13 = by * z_g__x_c__R13;
        // auto by__z_g__x_c__R21 = by * z_g__x_c__R21;
        // auto by__z_g__x_c__R22 = by * z_g__x_c__R22;
        // auto by__z_g__x_c__R23 = by * z_g__x_c__R23;
        // auto by__z_g__x_c__R31 = by * z_g__x_c__R31;
        // auto by__z_g__x_c__R32 = by * z_g__x_c__R32;
        // auto by__z_g__x_c__R33 = by * z_g__x_c__R33;
        // auto by__z_g__y_c__R11 = by * z_g__y_c__R11;
        // auto by__z_g__y_c__R12 = by * z_g__y_c__R12;
        // auto by__z_g__y_c__R13 = by * z_g__y_c__R13;
        // auto by__z_g__y_c__R21 = by * z_g__y_c__R21;
        // auto by__z_g__y_c__R22 = by * z_g__y_c__R22;
        // auto by__z_g__y_c__R23 = by * z_g__y_c__R23;
        // auto by__z_g__y_c__R31 = by * z_g__y_c__R31;
        // auto by__z_g__y_c__R32 = by * z_g__y_c__R32;
        // auto by__z_g__y_c__R33 = by * z_g__y_c__R33;
        // auto by__z_g__z_c__R11 = by * z_g__z_c__R11;
        // auto by__z_g__z_c__R12 = by * z_g__z_c__R12;
        // auto by__z_g__z_c__R13 = by * z_g__z_c__R13;
        // auto by__z_g__z_c__R21 = by * z_g__z_c__R21;
        // auto by__z_g__z_c__R22 = by * z_g__z_c__R22;
        // auto by__z_g__z_c__R23 = by * z_g__z_c__R23;
        // auto by__z_g__z_c__R31 = by * z_g__z_c__R31;
        // auto by__z_g__z_c__R32 = by * z_g__z_c__R32;
        // auto by__z_g__z_c__R33 = by * z_g__z_c__R33;
        // auto bz__a_g__a_c__R11 = bz * a_g__a_c__R11;
        // auto bz__a_g__a_c__R12 = bz * a_g__a_c__R12;
        // auto bz__a_g__a_c__R13 = bz * a_g__a_c__R13;
        // auto bz__a_g__a_c__R21 = bz * a_g__a_c__R21;
        // auto bz__a_g__a_c__R22 = bz * a_g__a_c__R22;
        // auto bz__a_g__a_c__R23 = bz * a_g__a_c__R23;
        // auto bz__a_g__a_c__R31 = bz * a_g__a_c__R31;
        // auto bz__a_g__a_c__R32 = bz * a_g__a_c__R32;
        // auto bz__a_g__a_c__R33 = bz * a_g__a_c__R33;
        // auto bz__a_g__b_c__R11 = bz * a_g__b_c__R11;
        // auto bz__a_g__b_c__R12 = bz * a_g__b_c__R12;
        // auto bz__a_g__b_c__R13 = bz * a_g__b_c__R13;
        // auto bz__a_g__b_c__R21 = bz * a_g__b_c__R21;
        // auto bz__a_g__b_c__R22 = bz * a_g__b_c__R22;
        // auto bz__a_g__b_c__R23 = bz * a_g__b_c__R23;
        // auto bz__a_g__b_c__R31 = bz * a_g__b_c__R31;
        // auto bz__a_g__b_c__R32 = bz * a_g__b_c__R32;
        // auto bz__a_g__b_c__R33 = bz * a_g__b_c__R33;
        // auto bz__a_g__c_c__R11 = bz * a_g__c_c__R11;
        // auto bz__a_g__c_c__R12 = bz * a_g__c_c__R12;
        // auto bz__a_g__c_c__R13 = bz * a_g__c_c__R13;
        // auto bz__a_g__c_c__R21 = bz * a_g__c_c__R21;
        // auto bz__a_g__c_c__R22 = bz * a_g__c_c__R22;
        // auto bz__a_g__c_c__R23 = bz * a_g__c_c__R23;
        // auto bz__a_g__c_c__R31 = bz * a_g__c_c__R31;
        // auto bz__a_g__c_c__R32 = bz * a_g__c_c__R32;
        // auto bz__a_g__c_c__R33 = bz * a_g__c_c__R33;
        // auto bz__a_g__x_c__R11 = bz * a_g__x_c__R11;
        // auto bz__a_g__x_c__R12 = bz * a_g__x_c__R12;
        // auto bz__a_g__x_c__R13 = bz * a_g__x_c__R13;
        // auto bz__a_g__x_c__R21 = bz * a_g__x_c__R21;
        // auto bz__a_g__x_c__R22 = bz * a_g__x_c__R22;
        // auto bz__a_g__x_c__R23 = bz * a_g__x_c__R23;
        // auto bz__a_g__x_c__R31 = bz * a_g__x_c__R31;
        // auto bz__a_g__x_c__R32 = bz * a_g__x_c__R32;
        // auto bz__a_g__x_c__R33 = bz * a_g__x_c__R33;
        // auto bz__a_g__y_c__R11 = bz * a_g__y_c__R11;
        // auto bz__a_g__y_c__R12 = bz * a_g__y_c__R12;
        // auto bz__a_g__y_c__R13 = bz * a_g__y_c__R13;
        // auto bz__a_g__y_c__R21 = bz * a_g__y_c__R21;
        // auto bz__a_g__y_c__R22 = bz * a_g__y_c__R22;
        // auto bz__a_g__y_c__R23 = bz * a_g__y_c__R23;
        // auto bz__a_g__y_c__R31 = bz * a_g__y_c__R31;
        // auto bz__a_g__y_c__R32 = bz * a_g__y_c__R32;
        // auto bz__a_g__y_c__R33 = bz * a_g__y_c__R33;
        // auto bz__a_g__z_c__R11 = bz * a_g__z_c__R11;
        // auto bz__a_g__z_c__R12 = bz * a_g__z_c__R12;
        // auto bz__a_g__z_c__R13 = bz * a_g__z_c__R13;
        // auto bz__a_g__z_c__R21 = bz * a_g__z_c__R21;
        // auto bz__a_g__z_c__R22 = bz * a_g__z_c__R22;
        // auto bz__a_g__z_c__R23 = bz * a_g__z_c__R23;
        // auto bz__a_g__z_c__R31 = bz * a_g__z_c__R31;
        // auto bz__a_g__z_c__R32 = bz * a_g__z_c__R32;
        // auto bz__a_g__z_c__R33 = bz * a_g__z_c__R33;
        // auto bz__b_g__a_c__R11 = bz * b_g__a_c__R11;
        // auto bz__b_g__a_c__R12 = bz * b_g__a_c__R12;
        // auto bz__b_g__a_c__R13 = bz * b_g__a_c__R13;
        // auto bz__b_g__a_c__R21 = bz * b_g__a_c__R21;
        // auto bz__b_g__a_c__R22 = bz * b_g__a_c__R22;
        // auto bz__b_g__a_c__R23 = bz * b_g__a_c__R23;
        // auto bz__b_g__a_c__R31 = bz * b_g__a_c__R31;
        // auto bz__b_g__a_c__R32 = bz * b_g__a_c__R32;
        // auto bz__b_g__a_c__R33 = bz * b_g__a_c__R33;
        // auto bz__b_g__b_c__R11 = bz * b_g__b_c__R11;
        // auto bz__b_g__b_c__R12 = bz * b_g__b_c__R12;
        // auto bz__b_g__b_c__R13 = bz * b_g__b_c__R13;
        // auto bz__b_g__b_c__R21 = bz * b_g__b_c__R21;
        // auto bz__b_g__b_c__R22 = bz * b_g__b_c__R22;
        // auto bz__b_g__b_c__R23 = bz * b_g__b_c__R23;
        // auto bz__b_g__b_c__R31 = bz * b_g__b_c__R31;
        // auto bz__b_g__b_c__R32 = bz * b_g__b_c__R32;
        // auto bz__b_g__b_c__R33 = bz * b_g__b_c__R33;
        // auto bz__b_g__c_c__R11 = bz * b_g__c_c__R11;
        // auto bz__b_g__c_c__R12 = bz * b_g__c_c__R12;
        // auto bz__b_g__c_c__R13 = bz * b_g__c_c__R13;
        // auto bz__b_g__c_c__R21 = bz * b_g__c_c__R21;
        // auto bz__b_g__c_c__R22 = bz * b_g__c_c__R22;
        // auto bz__b_g__c_c__R23 = bz * b_g__c_c__R23;
        // auto bz__b_g__c_c__R31 = bz * b_g__c_c__R31;
        // auto bz__b_g__c_c__R32 = bz * b_g__c_c__R32;
        // auto bz__b_g__c_c__R33 = bz * b_g__c_c__R33;
        // auto bz__b_g__x_c__R11 = bz * b_g__x_c__R11;
        // auto bz__b_g__x_c__R12 = bz * b_g__x_c__R12;
        // auto bz__b_g__x_c__R13 = bz * b_g__x_c__R13;
        // auto bz__b_g__x_c__R21 = bz * b_g__x_c__R21;
        // auto bz__b_g__x_c__R22 = bz * b_g__x_c__R22;
        // auto bz__b_g__x_c__R23 = bz * b_g__x_c__R23;
        // auto bz__b_g__x_c__R31 = bz * b_g__x_c__R31;
        // auto bz__b_g__x_c__R32 = bz * b_g__x_c__R32;
        // auto bz__b_g__x_c__R33 = bz * b_g__x_c__R33;
        // auto bz__b_g__y_c__R11 = bz * b_g__y_c__R11;
        // auto bz__b_g__y_c__R12 = bz * b_g__y_c__R12;
        // auto bz__b_g__y_c__R13 = bz * b_g__y_c__R13;
        // auto bz__b_g__y_c__R21 = bz * b_g__y_c__R21;
        // auto bz__b_g__y_c__R22 = bz * b_g__y_c__R22;
        // auto bz__b_g__y_c__R23 = bz * b_g__y_c__R23;
        // auto bz__b_g__y_c__R31 = bz * b_g__y_c__R31;
        // auto bz__b_g__y_c__R32 = bz * b_g__y_c__R32;
        // auto bz__b_g__y_c__R33 = bz * b_g__y_c__R33;
        // auto bz__b_g__z_c__R11 = bz * b_g__z_c__R11;
        // auto bz__b_g__z_c__R12 = bz * b_g__z_c__R12;
        // auto bz__b_g__z_c__R13 = bz * b_g__z_c__R13;
        // auto bz__b_g__z_c__R21 = bz * b_g__z_c__R21;
        // auto bz__b_g__z_c__R22 = bz * b_g__z_c__R22;
        // auto bz__b_g__z_c__R23 = bz * b_g__z_c__R23;
        // auto bz__b_g__z_c__R31 = bz * b_g__z_c__R31;
        // auto bz__b_g__z_c__R32 = bz * b_g__z_c__R32;
        // auto bz__b_g__z_c__R33 = bz * b_g__z_c__R33;
        // auto bz__c_g__a_c__R11 = bz * c_g__a_c__R11;
        // auto bz__c_g__a_c__R12 = bz * c_g__a_c__R12;
        // auto bz__c_g__a_c__R13 = bz * c_g__a_c__R13;
        // auto bz__c_g__a_c__R21 = bz * c_g__a_c__R21;
        // auto bz__c_g__a_c__R22 = bz * c_g__a_c__R22;
        // auto bz__c_g__a_c__R23 = bz * c_g__a_c__R23;
        // auto bz__c_g__a_c__R31 = bz * c_g__a_c__R31;
        // auto bz__c_g__a_c__R32 = bz * c_g__a_c__R32;
        // auto bz__c_g__a_c__R33 = bz * c_g__a_c__R33;
        // auto bz__c_g__b_c__R11 = bz * c_g__b_c__R11;
        // auto bz__c_g__b_c__R12 = bz * c_g__b_c__R12;
        // auto bz__c_g__b_c__R13 = bz * c_g__b_c__R13;
        // auto bz__c_g__b_c__R21 = bz * c_g__b_c__R21;
        // auto bz__c_g__b_c__R22 = bz * c_g__b_c__R22;
        // auto bz__c_g__b_c__R23 = bz * c_g__b_c__R23;
        // auto bz__c_g__b_c__R31 = bz * c_g__b_c__R31;
        // auto bz__c_g__b_c__R32 = bz * c_g__b_c__R32;
        // auto bz__c_g__b_c__R33 = bz * c_g__b_c__R33;
        // auto bz__c_g__c_c__R11 = bz * c_g__c_c__R11;
        // auto bz__c_g__c_c__R12 = bz * c_g__c_c__R12;
        // auto bz__c_g__c_c__R13 = bz * c_g__c_c__R13;
        // auto bz__c_g__c_c__R21 = bz * c_g__c_c__R21;
        // auto bz__c_g__c_c__R22 = bz * c_g__c_c__R22;
        // auto bz__c_g__c_c__R23 = bz * c_g__c_c__R23;
        // auto bz__c_g__c_c__R31 = bz * c_g__c_c__R31;
        // auto bz__c_g__c_c__R32 = bz * c_g__c_c__R32;
        // auto bz__c_g__c_c__R33 = bz * c_g__c_c__R33;
        // auto bz__c_g__x_c__R11 = bz * c_g__x_c__R11;
        // auto bz__c_g__x_c__R12 = bz * c_g__x_c__R12;
        // auto bz__c_g__x_c__R13 = bz * c_g__x_c__R13;
        // auto bz__c_g__x_c__R21 = bz * c_g__x_c__R21;
        // auto bz__c_g__x_c__R22 = bz * c_g__x_c__R22;
        // auto bz__c_g__x_c__R23 = bz * c_g__x_c__R23;
        // auto bz__c_g__x_c__R31 = bz * c_g__x_c__R31;
        // auto bz__c_g__x_c__R32 = bz * c_g__x_c__R32;
        // auto bz__c_g__x_c__R33 = bz * c_g__x_c__R33;
        // auto bz__c_g__y_c__R11 = bz * c_g__y_c__R11;
        // auto bz__c_g__y_c__R12 = bz * c_g__y_c__R12;
        // auto bz__c_g__y_c__R13 = bz * c_g__y_c__R13;
        // auto bz__c_g__y_c__R21 = bz * c_g__y_c__R21;
        // auto bz__c_g__y_c__R22 = bz * c_g__y_c__R22;
        // auto bz__c_g__y_c__R23 = bz * c_g__y_c__R23;
        // auto bz__c_g__y_c__R31 = bz * c_g__y_c__R31;
        // auto bz__c_g__y_c__R32 = bz * c_g__y_c__R32;
        // auto bz__c_g__y_c__R33 = bz * c_g__y_c__R33;
        // auto bz__c_g__z_c__R11 = bz * c_g__z_c__R11;
        // auto bz__c_g__z_c__R12 = bz * c_g__z_c__R12;
        // auto bz__c_g__z_c__R13 = bz * c_g__z_c__R13;
        // auto bz__c_g__z_c__R21 = bz * c_g__z_c__R21;
        // auto bz__c_g__z_c__R22 = bz * c_g__z_c__R22;
        // auto bz__c_g__z_c__R23 = bz * c_g__z_c__R23;
        // auto bz__c_g__z_c__R31 = bz * c_g__z_c__R31;
        // auto bz__c_g__z_c__R32 = bz * c_g__z_c__R32;
        // auto bz__c_g__z_c__R33 = bz * c_g__z_c__R33;
        // auto bz__x_g__a_c__R11 = bz * x_g__a_c__R11;
        // auto bz__x_g__a_c__R12 = bz * x_g__a_c__R12;
        // auto bz__x_g__a_c__R13 = bz * x_g__a_c__R13;
        // auto bz__x_g__a_c__R21 = bz * x_g__a_c__R21;
        // auto bz__x_g__a_c__R22 = bz * x_g__a_c__R22;
        // auto bz__x_g__a_c__R23 = bz * x_g__a_c__R23;
        // auto bz__x_g__a_c__R31 = bz * x_g__a_c__R31;
        // auto bz__x_g__a_c__R32 = bz * x_g__a_c__R32;
        // auto bz__x_g__a_c__R33 = bz * x_g__a_c__R33;
        // auto bz__x_g__b_c__R11 = bz * x_g__b_c__R11;
        // auto bz__x_g__b_c__R12 = bz * x_g__b_c__R12;
        // auto bz__x_g__b_c__R13 = bz * x_g__b_c__R13;
        // auto bz__x_g__b_c__R21 = bz * x_g__b_c__R21;
        // auto bz__x_g__b_c__R22 = bz * x_g__b_c__R22;
        // auto bz__x_g__b_c__R23 = bz * x_g__b_c__R23;
        // auto bz__x_g__b_c__R31 = bz * x_g__b_c__R31;
        // auto bz__x_g__b_c__R32 = bz * x_g__b_c__R32;
        // auto bz__x_g__b_c__R33 = bz * x_g__b_c__R33;
        // auto bz__x_g__c_c__R11 = bz * x_g__c_c__R11;
        // auto bz__x_g__c_c__R12 = bz * x_g__c_c__R12;
        // auto bz__x_g__c_c__R13 = bz * x_g__c_c__R13;
        // auto bz__x_g__c_c__R21 = bz * x_g__c_c__R21;
        // auto bz__x_g__c_c__R22 = bz * x_g__c_c__R22;
        // auto bz__x_g__c_c__R23 = bz * x_g__c_c__R23;
        // auto bz__x_g__c_c__R31 = bz * x_g__c_c__R31;
        // auto bz__x_g__c_c__R32 = bz * x_g__c_c__R32;
        // auto bz__x_g__c_c__R33 = bz * x_g__c_c__R33;
        // auto bz__x_g__x_c__R11 = bz * x_g__x_c__R11;
        // auto bz__x_g__x_c__R12 = bz * x_g__x_c__R12;
        // auto bz__x_g__x_c__R13 = bz * x_g__x_c__R13;
        // auto bz__x_g__x_c__R21 = bz * x_g__x_c__R21;
        // auto bz__x_g__x_c__R22 = bz * x_g__x_c__R22;
        // auto bz__x_g__x_c__R23 = bz * x_g__x_c__R23;
        // auto bz__x_g__x_c__R31 = bz * x_g__x_c__R31;
        // auto bz__x_g__x_c__R32 = bz * x_g__x_c__R32;
        // auto bz__x_g__x_c__R33 = bz * x_g__x_c__R33;
        // auto bz__x_g__y_c__R11 = bz * x_g__y_c__R11;
        // auto bz__x_g__y_c__R12 = bz * x_g__y_c__R12;
        // auto bz__x_g__y_c__R13 = bz * x_g__y_c__R13;
        // auto bz__x_g__y_c__R21 = bz * x_g__y_c__R21;
        // auto bz__x_g__y_c__R22 = bz * x_g__y_c__R22;
        // auto bz__x_g__y_c__R23 = bz * x_g__y_c__R23;
        // auto bz__x_g__y_c__R31 = bz * x_g__y_c__R31;
        // auto bz__x_g__y_c__R32 = bz * x_g__y_c__R32;
        // auto bz__x_g__y_c__R33 = bz * x_g__y_c__R33;
        // auto bz__x_g__z_c__R11 = bz * x_g__z_c__R11;
        // auto bz__x_g__z_c__R12 = bz * x_g__z_c__R12;
        // auto bz__x_g__z_c__R13 = bz * x_g__z_c__R13;
        // auto bz__x_g__z_c__R21 = bz * x_g__z_c__R21;
        // auto bz__x_g__z_c__R22 = bz * x_g__z_c__R22;
        // auto bz__x_g__z_c__R23 = bz * x_g__z_c__R23;
        // auto bz__x_g__z_c__R31 = bz * x_g__z_c__R31;
        // auto bz__x_g__z_c__R32 = bz * x_g__z_c__R32;
        // auto bz__x_g__z_c__R33 = bz * x_g__z_c__R33;
        // auto bz__y_g__a_c__R11 = bz * y_g__a_c__R11;
        // auto bz__y_g__a_c__R12 = bz * y_g__a_c__R12;
        // auto bz__y_g__a_c__R13 = bz * y_g__a_c__R13;
        // auto bz__y_g__a_c__R21 = bz * y_g__a_c__R21;
        // auto bz__y_g__a_c__R22 = bz * y_g__a_c__R22;
        // auto bz__y_g__a_c__R23 = bz * y_g__a_c__R23;
        // auto bz__y_g__a_c__R31 = bz * y_g__a_c__R31;
        // auto bz__y_g__a_c__R32 = bz * y_g__a_c__R32;
        // auto bz__y_g__a_c__R33 = bz * y_g__a_c__R33;
        // auto bz__y_g__b_c__R11 = bz * y_g__b_c__R11;
        // auto bz__y_g__b_c__R12 = bz * y_g__b_c__R12;
        // auto bz__y_g__b_c__R13 = bz * y_g__b_c__R13;
        // auto bz__y_g__b_c__R21 = bz * y_g__b_c__R21;
        // auto bz__y_g__b_c__R22 = bz * y_g__b_c__R22;
        // auto bz__y_g__b_c__R23 = bz * y_g__b_c__R23;
        // auto bz__y_g__b_c__R31 = bz * y_g__b_c__R31;
        // auto bz__y_g__b_c__R32 = bz * y_g__b_c__R32;
        // auto bz__y_g__b_c__R33 = bz * y_g__b_c__R33;
        // auto bz__y_g__c_c__R11 = bz * y_g__c_c__R11;
        // auto bz__y_g__c_c__R12 = bz * y_g__c_c__R12;
        // auto bz__y_g__c_c__R13 = bz * y_g__c_c__R13;
        // auto bz__y_g__c_c__R21 = bz * y_g__c_c__R21;
        // auto bz__y_g__c_c__R22 = bz * y_g__c_c__R22;
        // auto bz__y_g__c_c__R23 = bz * y_g__c_c__R23;
        // auto bz__y_g__c_c__R31 = bz * y_g__c_c__R31;
        // auto bz__y_g__c_c__R32 = bz * y_g__c_c__R32;
        // auto bz__y_g__c_c__R33 = bz * y_g__c_c__R33;
        // auto bz__y_g__x_c__R11 = bz * y_g__x_c__R11;
        // auto bz__y_g__x_c__R12 = bz * y_g__x_c__R12;
        // auto bz__y_g__x_c__R13 = bz * y_g__x_c__R13;
        // auto bz__y_g__x_c__R21 = bz * y_g__x_c__R21;
        // auto bz__y_g__x_c__R22 = bz * y_g__x_c__R22;
        // auto bz__y_g__x_c__R23 = bz * y_g__x_c__R23;
        // auto bz__y_g__x_c__R31 = bz * y_g__x_c__R31;
        // auto bz__y_g__x_c__R32 = bz * y_g__x_c__R32;
        // auto bz__y_g__x_c__R33 = bz * y_g__x_c__R33;
        // auto bz__y_g__y_c__R11 = bz * y_g__y_c__R11;
        // auto bz__y_g__y_c__R12 = bz * y_g__y_c__R12;
        // auto bz__y_g__y_c__R13 = bz * y_g__y_c__R13;
        // auto bz__y_g__y_c__R21 = bz * y_g__y_c__R21;
        // auto bz__y_g__y_c__R22 = bz * y_g__y_c__R22;
        // auto bz__y_g__y_c__R23 = bz * y_g__y_c__R23;
        // auto bz__y_g__y_c__R31 = bz * y_g__y_c__R31;
        // auto bz__y_g__y_c__R32 = bz * y_g__y_c__R32;
        // auto bz__y_g__y_c__R33 = bz * y_g__y_c__R33;
        // auto bz__y_g__z_c__R11 = bz * y_g__z_c__R11;
        // auto bz__y_g__z_c__R12 = bz * y_g__z_c__R12;
        // auto bz__y_g__z_c__R13 = bz * y_g__z_c__R13;
        // auto bz__y_g__z_c__R21 = bz * y_g__z_c__R21;
        // auto bz__y_g__z_c__R22 = bz * y_g__z_c__R22;
        // auto bz__y_g__z_c__R23 = bz * y_g__z_c__R23;
        // auto bz__y_g__z_c__R31 = bz * y_g__z_c__R31;
        // auto bz__y_g__z_c__R32 = bz * y_g__z_c__R32;
        // auto bz__y_g__z_c__R33 = bz * y_g__z_c__R33;
        // auto bz__z_g__a_c__R11 = bz * z_g__a_c__R11;
        // auto bz__z_g__a_c__R12 = bz * z_g__a_c__R12;
        // auto bz__z_g__a_c__R13 = bz * z_g__a_c__R13;
        // auto bz__z_g__a_c__R21 = bz * z_g__a_c__R21;
        // auto bz__z_g__a_c__R22 = bz * z_g__a_c__R22;
        // auto bz__z_g__a_c__R23 = bz * z_g__a_c__R23;
        // auto bz__z_g__a_c__R31 = bz * z_g__a_c__R31;
        // auto bz__z_g__a_c__R32 = bz * z_g__a_c__R32;
        // auto bz__z_g__a_c__R33 = bz * z_g__a_c__R33;
        // auto bz__z_g__b_c__R11 = bz * z_g__b_c__R11;
        // auto bz__z_g__b_c__R12 = bz * z_g__b_c__R12;
        // auto bz__z_g__b_c__R13 = bz * z_g__b_c__R13;
        // auto bz__z_g__b_c__R21 = bz * z_g__b_c__R21;
        // auto bz__z_g__b_c__R22 = bz * z_g__b_c__R22;
        // auto bz__z_g__b_c__R23 = bz * z_g__b_c__R23;
        // auto bz__z_g__b_c__R31 = bz * z_g__b_c__R31;
        // auto bz__z_g__b_c__R32 = bz * z_g__b_c__R32;
        // auto bz__z_g__b_c__R33 = bz * z_g__b_c__R33;
        // auto bz__z_g__c_c__R11 = bz * z_g__c_c__R11;
        // auto bz__z_g__c_c__R12 = bz * z_g__c_c__R12;
        // auto bz__z_g__c_c__R13 = bz * z_g__c_c__R13;
        // auto bz__z_g__c_c__R21 = bz * z_g__c_c__R21;
        // auto bz__z_g__c_c__R22 = bz * z_g__c_c__R22;
        // auto bz__z_g__c_c__R23 = bz * z_g__c_c__R23;
        // auto bz__z_g__c_c__R31 = bz * z_g__c_c__R31;
        // auto bz__z_g__c_c__R32 = bz * z_g__c_c__R32;
        // auto bz__z_g__c_c__R33 = bz * z_g__c_c__R33;
        // auto bz__z_g__x_c__R11 = bz * z_g__x_c__R11;
        // auto bz__z_g__x_c__R12 = bz * z_g__x_c__R12;
        // auto bz__z_g__x_c__R13 = bz * z_g__x_c__R13;
        // auto bz__z_g__x_c__R21 = bz * z_g__x_c__R21;
        // auto bz__z_g__x_c__R22 = bz * z_g__x_c__R22;
        // auto bz__z_g__x_c__R23 = bz * z_g__x_c__R23;
        // auto bz__z_g__x_c__R31 = bz * z_g__x_c__R31;
        // auto bz__z_g__x_c__R32 = bz * z_g__x_c__R32;
        // auto bz__z_g__x_c__R33 = bz * z_g__x_c__R33;
        // auto bz__z_g__y_c__R11 = bz * z_g__y_c__R11;
        // auto bz__z_g__y_c__R12 = bz * z_g__y_c__R12;
        // auto bz__z_g__y_c__R13 = bz * z_g__y_c__R13;
        // auto bz__z_g__y_c__R21 = bz * z_g__y_c__R21;
        // auto bz__z_g__y_c__R22 = bz * z_g__y_c__R22;
        // auto bz__z_g__y_c__R23 = bz * z_g__y_c__R23;
        // auto bz__z_g__y_c__R31 = bz * z_g__y_c__R31;
        // auto bz__z_g__y_c__R32 = bz * z_g__y_c__R32;
        // auto bz__z_g__y_c__R33 = bz * z_g__y_c__R33;
        // auto bz__z_g__z_c__R11 = bz * z_g__z_c__R11;
        // auto bz__z_g__z_c__R12 = bz * z_g__z_c__R12;
        // auto bz__z_g__z_c__R13 = bz * z_g__z_c__R13;
        // auto bz__z_g__z_c__R21 = bz * z_g__z_c__R21;
        // auto bz__z_g__z_c__R22 = bz * z_g__z_c__R22;
        // auto bz__z_g__z_c__R23 = bz * z_g__z_c__R23;
        // auto bz__z_g__z_c__R31 = bz * z_g__z_c__R31;
        // auto bz__z_g__z_c__R32 = bz * z_g__z_c__R32;
        // auto bz__z_g__z_c__R33 = bz * z_g__z_c__R33;

        // Manually optimized shortcuts for frequently appearing expressions.
        // Touch it on your own risk.
        // auto common_000 = ;
        // auto common_001 = -(c_c__R11) + tx__b_g__c_c__R11 + wm.R12 - tx__b_g__R12 - a_c__R13 + tx__b_g__a_c__R13 + c_g__c_c__R21
        //     - tx__a_g__c_c__R21 - c_g__R22 + tx__a_g__R22 + c_g__a_c__R23 - tx__a_g__a_c__R23 - b_g__c_c__R31 + tx__c_c__R31 + b_g__R32
        //     - tx__R32 - b_g__a_c__R33 + tx__a_c__R33;
        // auto common_002 = -(c_g__c_c__R11)-ty__b_g__c_c__R11 + c_g__R12 + ty__b_g__R12 - c_g__a_c__R13 - ty__b_g__a_c__R13 + c_c__R21
        //     + ty__a_g__c_c__R21 - wm.R22 - ty__a_g__R22 + a_c__R23 + ty__a_g__a_c__R23 + a_g__c_c__R31 - ty__c_c__R31 - a_g__R32 +
        //     ty__R32
        //     + a_g__a_c__R33 - ty__a_c__R33;
        // auto common_003 = tx__c_g__c_c__R11 + ty__c_c__R11 - tx__c_g__R12 - ty__R12 + tx__c_g__a_c__R13 + ty__a_c__R13 - tx__c_c__R21
        //     - ty__c_g__c_c__R21 + tx__R22 + ty__c_g__R22 - tx__a_c__R23 - ty__c_g__a_c__R23 - tx__a_g__c_c__R31 + ty__b_g__c_c__R31
        //     + tx__a_g__R32 - ty__b_g__R32 - tx__a_g__a_c__R33 + ty__b_g__a_c__R33;
        //
        // auto common_010 = su * (b_g__R11 - a_g__R21 + wm.R31 - c_c * (b_g__R12 - a_g__R22 + wm.R32) + b_c * (b_g__R13 - a_g__R23 +
        // wm.R33)); auto common_011 = sy * (b_g__R12 - a_g__R22 - c_c * (b_g__R11 - a_g__R21 + wm.R31) + wm.R32 - a_c * (b_g__R13 -
        // a_g__R23 + wm.R33)); auto common_012 = sz * (b_g__R13 - a_g__R23 + b_c * (b_g__R11 - a_g__R21 + wm.R31) + a_c * (b_g__R12 -
        // a_g__R22 + wm.R32) + wm.R33);
        //
        // auto common_021 = wm.R23 - c_g__R13;
        // auto common_022 = wm.R13 - c_g__R23;
        // auto common_023 = wm.R21 - c_g__R11;
        // auto common_024 = wm.R12 - c_g__R22;
        //
        // auto common_031 = common_022 + b_g__R33;
        // auto common_032 = common_024 + b_g__R32;
        // auto common_033 = wm.R11 - c_g__R21 + b_g__R31;
        //
        // auto common_034 = common_021 + a_g__R33;
        // auto common_035 = -(c_g__R12) + wm.R22 + a_g__R32;
        // auto common_036 = common_023 + a_g__R31;
        //
        // // auto common_031 = ;
        // // auto common_032 = ;
        // // auto common_033 = ;
        // //
        // auto common_051 = b_c * common_031;
        // auto common_052 = a_c * common_031;
        // // auto common_031 = ;
        //
        // auto common_100 = pow(common_001, 2) + pow(common_002, 2) + pow(common_003, 2);
        // // auto common_011 = common_2 * common_5 + common_0 * common_4 + common_1 * common_3;
        //
        // auto common_200 = sqrt(common_100);
        // auto common_201 = pow(common_100, 3. / 2.);
        // // auto common_032 = fabs(common_11);

        // drdxg
        this->template set_global_derivative<0>() =
            ((c_c * c_g * wm.R11 + b_g * c_c * ty * wm.R11 - c_g * wm.R12 - b_g * ty * wm.R12 + a_c * c_g * wm.R13 + a_c * b_g * ty * wm.R13
              - c_c * wm.R21 - a_g * c_c * ty * wm.R21 + wm.R22 + a_g * ty * wm.R22 - a_c * wm.R23 - a_c * a_g * ty * wm.R23
              - a_g * c_c * wm.R31 + c_c * ty * wm.R31 + a_g * wm.R32 - ty * wm.R32 - a_c * a_g * wm.R33 + a_c * ty * wm.R33)
             * ((c_c * c_g * tx * wm.R11 + c_c * ty * wm.R11 - c_g * tx * wm.R12 - ty * wm.R12 + a_c * c_g * tx * wm.R13 + a_c * ty * wm.R13
                 - c_c * tx * wm.R21 - c_c * c_g * ty * wm.R21 + tx * wm.R22 + c_g * ty * wm.R22 - a_c * tx * wm.R23
                 - a_c * c_g * ty * wm.R23 - a_g * c_c * tx * wm.R31 + b_g * c_c * ty * wm.R31 + a_g * tx * wm.R32 - b_g * ty * wm.R32
                 - a_c * a_g * tx * wm.R33 + a_c * b_g * ty * wm.R33)
                    * (-z_c - z_g - z_l
                       - sz
                           * (b_g * wm.R13 - a_g * wm.R23 + b_c * (b_g * wm.R11 - a_g * wm.R21 + wm.R31)
                              + a_c * (b_g * wm.R12 - a_g * wm.R22 + wm.R32) + wm.R33)
                       - sv
                           * (b_g * wm.R12 - a_g * wm.R22 - c_c * (b_g * wm.R11 - a_g * wm.R21 + wm.R31) + wm.R32
                              - a_c * (b_g * wm.R13 - a_g * wm.R23 + wm.R33))
                       - su
                           * (b_g * wm.R11 - a_g * wm.R21 + wm.R31 - c_c * (b_g * wm.R12 - a_g * wm.R22 + wm.R32)
                              + b_c * (b_g * wm.R13 - a_g * wm.R23 + wm.R33)))
                + (-(c_c * wm.R11) + b_g * c_c * tx * wm.R11 + wm.R12 - b_g * tx * wm.R12 - a_c * wm.R13 + a_c * b_g * tx * wm.R13
                   + c_c * c_g * wm.R21 - a_g * c_c * tx * wm.R21 - c_g * wm.R22 + a_g * tx * wm.R22 + a_c * c_g * wm.R23
                   - a_c * a_g * tx * wm.R23 - b_g * c_c * wm.R31 + c_c * tx * wm.R31 + b_g * wm.R32 - tx * wm.R32 - a_c * b_g * wm.R33
                   + a_c * tx * wm.R33)
                    * (by - y_c - y_g - y_l
                       - sz
                           * (-(c_g * wm.R13) + wm.R23 + b_c * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31)
                              + a_c * (-(c_g * wm.R12) + wm.R22 + a_g * wm.R32) + a_g * wm.R33)
                       - sv
                           * (-(c_g * wm.R12) + wm.R22 - c_c * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31) + a_g * wm.R32
                              - a_c * (-(c_g * wm.R13) + wm.R23 + a_g * wm.R33))
                       - su
                           * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31 - c_c * (-(c_g * wm.R12) + wm.R22 + a_g * wm.R32)
                              + b_c * (-(c_g * wm.R13) + wm.R23 + a_g * wm.R33)))
                + (-(c_c * c_g * wm.R11) - b_g * c_c * ty * wm.R11 + c_g * wm.R12 + b_g * ty * wm.R12 - a_c * c_g * wm.R13
                   - a_c * b_g * ty * wm.R13 + c_c * wm.R21 + a_g * c_c * ty * wm.R21 - wm.R22 - a_g * ty * wm.R22 + a_c * wm.R23
                   + a_c * a_g * ty * wm.R23 + a_g * c_c * wm.R31 - c_c * ty * wm.R31 - a_g * wm.R32 + ty * wm.R32 + a_c * a_g * wm.R33
                   - a_c * ty * wm.R33)
                    * (bx - x_c - x_g - x_l
                       - sz
                           * (wm.R13 - c_g * wm.R23 + b_c * (wm.R11 - c_g * wm.R21 + b_g * wm.R31)
                              + a_c * (wm.R12 - c_g * wm.R22 + b_g * wm.R32) + b_g * wm.R33)
                       - sv
                           * (wm.R12 - c_g * wm.R22 - c_c * (wm.R11 - c_g * wm.R21 + b_g * wm.R31) + b_g * wm.R32
                              - a_c * (wm.R13 - c_g * wm.R23 + b_g * wm.R33))
                       - su
                           * (wm.R11 - c_g * wm.R21 + b_g * wm.R31 - c_c * (wm.R12 - c_g * wm.R22 + b_g * wm.R32)
                              + b_c * (wm.R13 - c_g * wm.R23 + b_g * wm.R33)))))
            / (sqrt(pow(-(c_c * wm.R11) + b_g * c_c * tx * wm.R11 + wm.R12 - b_g * tx * wm.R12 - a_c * wm.R13 + a_c * b_g * tx * wm.R13
                            + c_c * c_g * wm.R21 - a_g * c_c * tx * wm.R21 - c_g * wm.R22 + a_g * tx * wm.R22 + a_c * c_g * wm.R23
                            - a_c * a_g * tx * wm.R23 - b_g * c_c * wm.R31 + c_c * tx * wm.R31 + b_g * wm.R32 - tx * wm.R32
                            - a_c * b_g * wm.R33 + a_c * tx * wm.R33,
                        2)
                    + pow(-(c_c * c_g * wm.R11) - b_g * c_c * ty * wm.R11 + c_g * wm.R12 + b_g * ty * wm.R12 - a_c * c_g * wm.R13
                              - a_c * b_g * ty * wm.R13 + c_c * wm.R21 + a_g * c_c * ty * wm.R21 - wm.R22 - a_g * ty * wm.R22 + a_c * wm.R23
                              + a_c * a_g * ty * wm.R23 + a_g * c_c * wm.R31 - c_c * ty * wm.R31 - a_g * wm.R32 + ty * wm.R32
                              + a_c * a_g * wm.R33 - a_c * ty * wm.R33,
                          2)
                    + pow(c_c * c_g * tx * wm.R11 + c_c * ty * wm.R11 - c_g * tx * wm.R12 - ty * wm.R12 + a_c * c_g * tx * wm.R13
                              + a_c * ty * wm.R13 - c_c * tx * wm.R21 - c_c * c_g * ty * wm.R21 + tx * wm.R22 + c_g * ty * wm.R22
                              - a_c * tx * wm.R23 - a_c * c_g * ty * wm.R23 - a_g * c_c * tx * wm.R31 + b_g * c_c * ty * wm.R31
                              + a_g * tx * wm.R32 - b_g * ty * wm.R32 - a_c * a_g * tx * wm.R33 + a_c * b_g * ty * wm.R33,
                          2))
               * fabs(((c_c * c_g * tx * wm.R11 + c_c * ty * wm.R11 - c_g * tx * wm.R12 - ty * wm.R12 + a_c * c_g * tx * wm.R13
                        + a_c * ty * wm.R13 - c_c * tx * wm.R21 - c_c * c_g * ty * wm.R21 + tx * wm.R22 + c_g * ty * wm.R22
                        - a_c * tx * wm.R23 - a_c * c_g * ty * wm.R23 - a_g * c_c * tx * wm.R31 + b_g * c_c * ty * wm.R31
                        + a_g * tx * wm.R32 - b_g * ty * wm.R32 - a_c * a_g * tx * wm.R33 + a_c * b_g * ty * wm.R33)
                           * (-z_c - z_g - z_l
                              - sz
                                  * (b_g * wm.R13 - a_g * wm.R23 + b_c * (b_g * wm.R11 - a_g * wm.R21 + wm.R31)
                                     + a_c * (b_g * wm.R12 - a_g * wm.R22 + wm.R32) + wm.R33)
                              - sv
                                  * (b_g * wm.R12 - a_g * wm.R22 - c_c * (b_g * wm.R11 - a_g * wm.R21 + wm.R31) + wm.R32
                                     - a_c * (b_g * wm.R13 - a_g * wm.R23 + wm.R33))
                              - su
                                  * (b_g * wm.R11 - a_g * wm.R21 + wm.R31 - c_c * (b_g * wm.R12 - a_g * wm.R22 + wm.R32)
                                     + b_c * (b_g * wm.R13 - a_g * wm.R23 + wm.R33)))
                       + (-(c_c * wm.R11) + b_g * c_c * tx * wm.R11 + wm.R12 - b_g * tx * wm.R12 - a_c * wm.R13 + a_c * b_g * tx * wm.R13
                          + c_c * c_g * wm.R21 - a_g * c_c * tx * wm.R21 - c_g * wm.R22 + a_g * tx * wm.R22 + a_c * c_g * wm.R23
                          - a_c * a_g * tx * wm.R23 - b_g * c_c * wm.R31 + c_c * tx * wm.R31 + b_g * wm.R32 - tx * wm.R32
                          - a_c * b_g * wm.R33 + a_c * tx * wm.R33)
                           * (by - y_c - y_g - y_l
                              - sz
                                  * (-(c_g * wm.R13) + wm.R23 + b_c * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31)
                                     + a_c * (-(c_g * wm.R12) + wm.R22 + a_g * wm.R32) + a_g * wm.R33)
                              - sv
                                  * (-(c_g * wm.R12) + wm.R22 - c_c * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31) + a_g * wm.R32
                                     - a_c * (-(c_g * wm.R13) + wm.R23 + a_g * wm.R33))
                              - su
                                  * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31 - c_c * (-(c_g * wm.R12) + wm.R22 + a_g * wm.R32)
                                     + b_c * (-(c_g * wm.R13) + wm.R23 + a_g * wm.R33)))
                       + (-(c_c * c_g * wm.R11) - b_g * c_c * ty * wm.R11 + c_g * wm.R12 + b_g * ty * wm.R12 - a_c * c_g * wm.R13
                          - a_c * b_g * ty * wm.R13 + c_c * wm.R21 + a_g * c_c * ty * wm.R21 - wm.R22 - a_g * ty * wm.R22 + a_c * wm.R23
                          + a_c * a_g * ty * wm.R23 + a_g * c_c * wm.R31 - c_c * ty * wm.R31 - a_g * wm.R32 + ty * wm.R32
                          + a_c * a_g * wm.R33 - a_c * ty * wm.R33)
                           * (bx - x_c - x_g - x_l
                              - sz
                                  * (wm.R13 - c_g * wm.R23 + b_c * (wm.R11 - c_g * wm.R21 + b_g * wm.R31)
                                     + a_c * (wm.R12 - c_g * wm.R22 + b_g * wm.R32) + b_g * wm.R33)
                              - sv
                                  * (wm.R12 - c_g * wm.R22 - c_c * (wm.R11 - c_g * wm.R21 + b_g * wm.R31) + b_g * wm.R32
                                     - a_c * (wm.R13 - c_g * wm.R23 + b_g * wm.R33))
                              - su
                                  * (wm.R11 - c_g * wm.R21 + b_g * wm.R31 - c_c * (wm.R12 - c_g * wm.R22 + b_g * wm.R32)
                                     + b_c * (wm.R13 - c_g * wm.R23 + b_g * wm.R33))))));
        // drdyg
        this->template set_global_derivative<1>() =
            ((c_c * wm.R11 - b_g * c_c * tx * wm.R11 - wm.R12 + b_g * tx * wm.R12 + a_c * wm.R13 - a_c * b_g * tx * wm.R13
              - c_c * c_g * wm.R21 + a_g * c_c * tx * wm.R21 + c_g * wm.R22 - a_g * tx * wm.R22 - a_c * c_g * wm.R23
              + a_c * a_g * tx * wm.R23 + b_g * c_c * wm.R31 - c_c * tx * wm.R31 - b_g * wm.R32 + tx * wm.R32 + a_c * b_g * wm.R33
              - a_c * tx * wm.R33)
             * ((c_c * c_g * tx * wm.R11 + c_c * ty * wm.R11 - c_g * tx * wm.R12 - ty * wm.R12 + a_c * c_g * tx * wm.R13 + a_c * ty * wm.R13
                 - c_c * tx * wm.R21 - c_c * c_g * ty * wm.R21 + tx * wm.R22 + c_g * ty * wm.R22 - a_c * tx * wm.R23
                 - a_c * c_g * ty * wm.R23 - a_g * c_c * tx * wm.R31 + b_g * c_c * ty * wm.R31 + a_g * tx * wm.R32 - b_g * ty * wm.R32
                 - a_c * a_g * tx * wm.R33 + a_c * b_g * ty * wm.R33)
                    * (-z_c - z_g - z_l
                       - sz
                           * (b_g * wm.R13 - a_g * wm.R23 + b_c * (b_g * wm.R11 - a_g * wm.R21 + wm.R31)
                              + a_c * (b_g * wm.R12 - a_g * wm.R22 + wm.R32) + wm.R33)
                       - sv
                           * (b_g * wm.R12 - a_g * wm.R22 - c_c * (b_g * wm.R11 - a_g * wm.R21 + wm.R31) + wm.R32
                              - a_c * (b_g * wm.R13 - a_g * wm.R23 + wm.R33))
                       - su
                           * (b_g * wm.R11 - a_g * wm.R21 + wm.R31 - c_c * (b_g * wm.R12 - a_g * wm.R22 + wm.R32)
                              + b_c * (b_g * wm.R13 - a_g * wm.R23 + wm.R33)))
                + (-(c_c * wm.R11) + b_g * c_c * tx * wm.R11 + wm.R12 - b_g * tx * wm.R12 - a_c * wm.R13 + a_c * b_g * tx * wm.R13
                   + c_c * c_g * wm.R21 - a_g * c_c * tx * wm.R21 - c_g * wm.R22 + a_g * tx * wm.R22 + a_c * c_g * wm.R23
                   - a_c * a_g * tx * wm.R23 - b_g * c_c * wm.R31 + c_c * tx * wm.R31 + b_g * wm.R32 - tx * wm.R32 - a_c * b_g * wm.R33
                   + a_c * tx * wm.R33)
                    * (by - y_c - y_g - y_l
                       - sz
                           * (-(c_g * wm.R13) + wm.R23 + b_c * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31)
                              + a_c * (-(c_g * wm.R12) + wm.R22 + a_g * wm.R32) + a_g * wm.R33)
                       - sv
                           * (-(c_g * wm.R12) + wm.R22 - c_c * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31) + a_g * wm.R32
                              - a_c * (-(c_g * wm.R13) + wm.R23 + a_g * wm.R33))
                       - su
                           * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31 - c_c * (-(c_g * wm.R12) + wm.R22 + a_g * wm.R32)
                              + b_c * (-(c_g * wm.R13) + wm.R23 + a_g * wm.R33)))
                + (-(c_c * c_g * wm.R11) - b_g * c_c * ty * wm.R11 + c_g * wm.R12 + b_g * ty * wm.R12 - a_c * c_g * wm.R13
                   - a_c * b_g * ty * wm.R13 + c_c * wm.R21 + a_g * c_c * ty * wm.R21 - wm.R22 - a_g * ty * wm.R22 + a_c * wm.R23
                   + a_c * a_g * ty * wm.R23 + a_g * c_c * wm.R31 - c_c * ty * wm.R31 - a_g * wm.R32 + ty * wm.R32 + a_c * a_g * wm.R33
                   - a_c * ty * wm.R33)
                    * (bx - x_c - x_g - x_l
                       - sz
                           * (wm.R13 - c_g * wm.R23 + b_c * (wm.R11 - c_g * wm.R21 + b_g * wm.R31)
                              + a_c * (wm.R12 - c_g * wm.R22 + b_g * wm.R32) + b_g * wm.R33)
                       - sv
                           * (wm.R12 - c_g * wm.R22 - c_c * (wm.R11 - c_g * wm.R21 + b_g * wm.R31) + b_g * wm.R32
                              - a_c * (wm.R13 - c_g * wm.R23 + b_g * wm.R33))
                       - su
                           * (wm.R11 - c_g * wm.R21 + b_g * wm.R31 - c_c * (wm.R12 - c_g * wm.R22 + b_g * wm.R32)
                              + b_c * (wm.R13 - c_g * wm.R23 + b_g * wm.R33)))))
            / (sqrt(pow(-(c_c * wm.R11) + b_g * c_c * tx * wm.R11 + wm.R12 - b_g * tx * wm.R12 - a_c * wm.R13 + a_c * b_g * tx * wm.R13
                            + c_c * c_g * wm.R21 - a_g * c_c * tx * wm.R21 - c_g * wm.R22 + a_g * tx * wm.R22 + a_c * c_g * wm.R23
                            - a_c * a_g * tx * wm.R23 - b_g * c_c * wm.R31 + c_c * tx * wm.R31 + b_g * wm.R32 - tx * wm.R32
                            - a_c * b_g * wm.R33 + a_c * tx * wm.R33,
                        2)
                    + pow(-(c_c * c_g * wm.R11) - b_g * c_c * ty * wm.R11 + c_g * wm.R12 + b_g * ty * wm.R12 - a_c * c_g * wm.R13
                              - a_c * b_g * ty * wm.R13 + c_c * wm.R21 + a_g * c_c * ty * wm.R21 - wm.R22 - a_g * ty * wm.R22 + a_c * wm.R23
                              + a_c * a_g * ty * wm.R23 + a_g * c_c * wm.R31 - c_c * ty * wm.R31 - a_g * wm.R32 + ty * wm.R32
                              + a_c * a_g * wm.R33 - a_c * ty * wm.R33,
                          2)
                    + pow(c_c * c_g * tx * wm.R11 + c_c * ty * wm.R11 - c_g * tx * wm.R12 - ty * wm.R12 + a_c * c_g * tx * wm.R13
                              + a_c * ty * wm.R13 - c_c * tx * wm.R21 - c_c * c_g * ty * wm.R21 + tx * wm.R22 + c_g * ty * wm.R22
                              - a_c * tx * wm.R23 - a_c * c_g * ty * wm.R23 - a_g * c_c * tx * wm.R31 + b_g * c_c * ty * wm.R31
                              + a_g * tx * wm.R32 - b_g * ty * wm.R32 - a_c * a_g * tx * wm.R33 + a_c * b_g * ty * wm.R33,
                          2))
               * fabs(((c_c * c_g * tx * wm.R11 + c_c * ty * wm.R11 - c_g * tx * wm.R12 - ty * wm.R12 + a_c * c_g * tx * wm.R13
                        + a_c * ty * wm.R13 - c_c * tx * wm.R21 - c_c * c_g * ty * wm.R21 + tx * wm.R22 + c_g * ty * wm.R22
                        - a_c * tx * wm.R23 - a_c * c_g * ty * wm.R23 - a_g * c_c * tx * wm.R31 + b_g * c_c * ty * wm.R31
                        + a_g * tx * wm.R32 - b_g * ty * wm.R32 - a_c * a_g * tx * wm.R33 + a_c * b_g * ty * wm.R33)
                           * (-z_c - z_g - z_l
                              - sz
                                  * (b_g * wm.R13 - a_g * wm.R23 + b_c * (b_g * wm.R11 - a_g * wm.R21 + wm.R31)
                                     + a_c * (b_g * wm.R12 - a_g * wm.R22 + wm.R32) + wm.R33)
                              - sv
                                  * (b_g * wm.R12 - a_g * wm.R22 - c_c * (b_g * wm.R11 - a_g * wm.R21 + wm.R31) + wm.R32
                                     - a_c * (b_g * wm.R13 - a_g * wm.R23 + wm.R33))
                              - su
                                  * (b_g * wm.R11 - a_g * wm.R21 + wm.R31 - c_c * (b_g * wm.R12 - a_g * wm.R22 + wm.R32)
                                     + b_c * (b_g * wm.R13 - a_g * wm.R23 + wm.R33)))
                       + (-(c_c * wm.R11) + b_g * c_c * tx * wm.R11 + wm.R12 - b_g * tx * wm.R12 - a_c * wm.R13 + a_c * b_g * tx * wm.R13
                          + c_c * c_g * wm.R21 - a_g * c_c * tx * wm.R21 - c_g * wm.R22 + a_g * tx * wm.R22 + a_c * c_g * wm.R23
                          - a_c * a_g * tx * wm.R23 - b_g * c_c * wm.R31 + c_c * tx * wm.R31 + b_g * wm.R32 - tx * wm.R32
                          - a_c * b_g * wm.R33 + a_c * tx * wm.R33)
                           * (by - y_c - y_g - y_l
                              - sz
                                  * (-(c_g * wm.R13) + wm.R23 + b_c * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31)
                                     + a_c * (-(c_g * wm.R12) + wm.R22 + a_g * wm.R32) + a_g * wm.R33)
                              - sv
                                  * (-(c_g * wm.R12) + wm.R22 - c_c * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31) + a_g * wm.R32
                                     - a_c * (-(c_g * wm.R13) + wm.R23 + a_g * wm.R33))
                              - su
                                  * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31 - c_c * (-(c_g * wm.R12) + wm.R22 + a_g * wm.R32)
                                     + b_c * (-(c_g * wm.R13) + wm.R23 + a_g * wm.R33)))
                       + (-(c_c * c_g * wm.R11) - b_g * c_c * ty * wm.R11 + c_g * wm.R12 + b_g * ty * wm.R12 - a_c * c_g * wm.R13
                          - a_c * b_g * ty * wm.R13 + c_c * wm.R21 + a_g * c_c * ty * wm.R21 - wm.R22 - a_g * ty * wm.R22 + a_c * wm.R23
                          + a_c * a_g * ty * wm.R23 + a_g * c_c * wm.R31 - c_c * ty * wm.R31 - a_g * wm.R32 + ty * wm.R32
                          + a_c * a_g * wm.R33 - a_c * ty * wm.R33)
                           * (bx - x_c - x_g - x_l
                              - sz
                                  * (wm.R13 - c_g * wm.R23 + b_c * (wm.R11 - c_g * wm.R21 + b_g * wm.R31)
                                     + a_c * (wm.R12 - c_g * wm.R22 + b_g * wm.R32) + b_g * wm.R33)
                              - sv
                                  * (wm.R12 - c_g * wm.R22 - c_c * (wm.R11 - c_g * wm.R21 + b_g * wm.R31) + b_g * wm.R32
                                     - a_c * (wm.R13 - c_g * wm.R23 + b_g * wm.R33))
                              - su
                                  * (wm.R11 - c_g * wm.R21 + b_g * wm.R31 - c_c * (wm.R12 - c_g * wm.R22 + b_g * wm.R32)
                                     + b_c * (wm.R13 - c_g * wm.R23 + b_g * wm.R33))))));
        // drdzg
        this->template set_global_derivative<2>() =
            ((-(c_c * c_g * tx * wm.R11) - c_c * ty * wm.R11 + c_g * tx * wm.R12 + ty * wm.R12 - a_c * c_g * tx * wm.R13 - a_c * ty * wm.R13
              + c_c * tx * wm.R21 + c_c * c_g * ty * wm.R21 - tx * wm.R22 - c_g * ty * wm.R22 + a_c * tx * wm.R23 + a_c * c_g * ty * wm.R23
              + a_g * c_c * tx * wm.R31 - b_g * c_c * ty * wm.R31 - a_g * tx * wm.R32 + b_g * ty * wm.R32 + a_c * a_g * tx * wm.R33
              - a_c * b_g * ty * wm.R33)
             * ((c_c * c_g * tx * wm.R11 + c_c * ty * wm.R11 - c_g * tx * wm.R12 - ty * wm.R12 + a_c * c_g * tx * wm.R13 + a_c * ty * wm.R13
                 - c_c * tx * wm.R21 - c_c * c_g * ty * wm.R21 + tx * wm.R22 + c_g * ty * wm.R22 - a_c * tx * wm.R23
                 - a_c * c_g * ty * wm.R23 - a_g * c_c * tx * wm.R31 + b_g * c_c * ty * wm.R31 + a_g * tx * wm.R32 - b_g * ty * wm.R32
                 - a_c * a_g * tx * wm.R33 + a_c * b_g * ty * wm.R33)
                    * (-z_c - z_g - z_l
                       - sz
                           * (b_g * wm.R13 - a_g * wm.R23 + b_c * (b_g * wm.R11 - a_g * wm.R21 + wm.R31)
                              + a_c * (b_g * wm.R12 - a_g * wm.R22 + wm.R32) + wm.R33)
                       - sv
                           * (b_g * wm.R12 - a_g * wm.R22 - c_c * (b_g * wm.R11 - a_g * wm.R21 + wm.R31) + wm.R32
                              - a_c * (b_g * wm.R13 - a_g * wm.R23 + wm.R33))
                       - su
                           * (b_g * wm.R11 - a_g * wm.R21 + wm.R31 - c_c * (b_g * wm.R12 - a_g * wm.R22 + wm.R32)
                              + b_c * (b_g * wm.R13 - a_g * wm.R23 + wm.R33)))
                + (-(c_c * wm.R11) + b_g * c_c * tx * wm.R11 + wm.R12 - b_g * tx * wm.R12 - a_c * wm.R13 + a_c * b_g * tx * wm.R13
                   + c_c * c_g * wm.R21 - a_g * c_c * tx * wm.R21 - c_g * wm.R22 + a_g * tx * wm.R22 + a_c * c_g * wm.R23
                   - a_c * a_g * tx * wm.R23 - b_g * c_c * wm.R31 + c_c * tx * wm.R31 + b_g * wm.R32 - tx * wm.R32 - a_c * b_g * wm.R33
                   + a_c * tx * wm.R33)
                    * (by - y_c - y_g - y_l
                       - sz
                           * (-(c_g * wm.R13) + wm.R23 + b_c * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31)
                              + a_c * (-(c_g * wm.R12) + wm.R22 + a_g * wm.R32) + a_g * wm.R33)
                       - sv
                           * (-(c_g * wm.R12) + wm.R22 - c_c * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31) + a_g * wm.R32
                              - a_c * (-(c_g * wm.R13) + wm.R23 + a_g * wm.R33))
                       - su
                           * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31 - c_c * (-(c_g * wm.R12) + wm.R22 + a_g * wm.R32)
                              + b_c * (-(c_g * wm.R13) + wm.R23 + a_g * wm.R33)))
                + (-(c_c * c_g * wm.R11) - b_g * c_c * ty * wm.R11 + c_g * wm.R12 + b_g * ty * wm.R12 - a_c * c_g * wm.R13
                   - a_c * b_g * ty * wm.R13 + c_c * wm.R21 + a_g * c_c * ty * wm.R21 - wm.R22 - a_g * ty * wm.R22 + a_c * wm.R23
                   + a_c * a_g * ty * wm.R23 + a_g * c_c * wm.R31 - c_c * ty * wm.R31 - a_g * wm.R32 + ty * wm.R32 + a_c * a_g * wm.R33
                   - a_c * ty * wm.R33)
                    * (bx - x_c - x_g - x_l
                       - sz
                           * (wm.R13 - c_g * wm.R23 + b_c * (wm.R11 - c_g * wm.R21 + b_g * wm.R31)
                              + a_c * (wm.R12 - c_g * wm.R22 + b_g * wm.R32) + b_g * wm.R33)
                       - sv
                           * (wm.R12 - c_g * wm.R22 - c_c * (wm.R11 - c_g * wm.R21 + b_g * wm.R31) + b_g * wm.R32
                              - a_c * (wm.R13 - c_g * wm.R23 + b_g * wm.R33))
                       - su
                           * (wm.R11 - c_g * wm.R21 + b_g * wm.R31 - c_c * (wm.R12 - c_g * wm.R22 + b_g * wm.R32)
                              + b_c * (wm.R13 - c_g * wm.R23 + b_g * wm.R33)))))
            / (sqrt(pow(-(c_c * wm.R11) + b_g * c_c * tx * wm.R11 + wm.R12 - b_g * tx * wm.R12 - a_c * wm.R13 + a_c * b_g * tx * wm.R13
                            + c_c * c_g * wm.R21 - a_g * c_c * tx * wm.R21 - c_g * wm.R22 + a_g * tx * wm.R22 + a_c * c_g * wm.R23
                            - a_c * a_g * tx * wm.R23 - b_g * c_c * wm.R31 + c_c * tx * wm.R31 + b_g * wm.R32 - tx * wm.R32
                            - a_c * b_g * wm.R33 + a_c * tx * wm.R33,
                        2)
                    + pow(-(c_c * c_g * wm.R11) - b_g * c_c * ty * wm.R11 + c_g * wm.R12 + b_g * ty * wm.R12 - a_c * c_g * wm.R13
                              - a_c * b_g * ty * wm.R13 + c_c * wm.R21 + a_g * c_c * ty * wm.R21 - wm.R22 - a_g * ty * wm.R22 + a_c * wm.R23
                              + a_c * a_g * ty * wm.R23 + a_g * c_c * wm.R31 - c_c * ty * wm.R31 - a_g * wm.R32 + ty * wm.R32
                              + a_c * a_g * wm.R33 - a_c * ty * wm.R33,
                          2)
                    + pow(c_c * c_g * tx * wm.R11 + c_c * ty * wm.R11 - c_g * tx * wm.R12 - ty * wm.R12 + a_c * c_g * tx * wm.R13
                              + a_c * ty * wm.R13 - c_c * tx * wm.R21 - c_c * c_g * ty * wm.R21 + tx * wm.R22 + c_g * ty * wm.R22
                              - a_c * tx * wm.R23 - a_c * c_g * ty * wm.R23 - a_g * c_c * tx * wm.R31 + b_g * c_c * ty * wm.R31
                              + a_g * tx * wm.R32 - b_g * ty * wm.R32 - a_c * a_g * tx * wm.R33 + a_c * b_g * ty * wm.R33,
                          2))
               * fabs(((c_c * c_g * tx * wm.R11 + c_c * ty * wm.R11 - c_g * tx * wm.R12 - ty * wm.R12 + a_c * c_g * tx * wm.R13
                        + a_c * ty * wm.R13 - c_c * tx * wm.R21 - c_c * c_g * ty * wm.R21 + tx * wm.R22 + c_g * ty * wm.R22
                        - a_c * tx * wm.R23 - a_c * c_g * ty * wm.R23 - a_g * c_c * tx * wm.R31 + b_g * c_c * ty * wm.R31
                        + a_g * tx * wm.R32 - b_g * ty * wm.R32 - a_c * a_g * tx * wm.R33 + a_c * b_g * ty * wm.R33)
                           * (-z_c - z_g - z_l
                              - sz
                                  * (b_g * wm.R13 - a_g * wm.R23 + b_c * (b_g * wm.R11 - a_g * wm.R21 + wm.R31)
                                     + a_c * (b_g * wm.R12 - a_g * wm.R22 + wm.R32) + wm.R33)
                              - sv
                                  * (b_g * wm.R12 - a_g * wm.R22 - c_c * (b_g * wm.R11 - a_g * wm.R21 + wm.R31) + wm.R32
                                     - a_c * (b_g * wm.R13 - a_g * wm.R23 + wm.R33))
                              - su
                                  * (b_g * wm.R11 - a_g * wm.R21 + wm.R31 - c_c * (b_g * wm.R12 - a_g * wm.R22 + wm.R32)
                                     + b_c * (b_g * wm.R13 - a_g * wm.R23 + wm.R33)))
                       + (-(c_c * wm.R11) + b_g * c_c * tx * wm.R11 + wm.R12 - b_g * tx * wm.R12 - a_c * wm.R13 + a_c * b_g * tx * wm.R13
                          + c_c * c_g * wm.R21 - a_g * c_c * tx * wm.R21 - c_g * wm.R22 + a_g * tx * wm.R22 + a_c * c_g * wm.R23
                          - a_c * a_g * tx * wm.R23 - b_g * c_c * wm.R31 + c_c * tx * wm.R31 + b_g * wm.R32 - tx * wm.R32
                          - a_c * b_g * wm.R33 + a_c * tx * wm.R33)
                           * (by - y_c - y_g - y_l
                              - sz
                                  * (-(c_g * wm.R13) + wm.R23 + b_c * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31)
                                     + a_c * (-(c_g * wm.R12) + wm.R22 + a_g * wm.R32) + a_g * wm.R33)
                              - sv
                                  * (-(c_g * wm.R12) + wm.R22 - c_c * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31) + a_g * wm.R32
                                     - a_c * (-(c_g * wm.R13) + wm.R23 + a_g * wm.R33))
                              - su
                                  * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31 - c_c * (-(c_g * wm.R12) + wm.R22 + a_g * wm.R32)
                                     + b_c * (-(c_g * wm.R13) + wm.R23 + a_g * wm.R33)))
                       + (-(c_c * c_g * wm.R11) - b_g * c_c * ty * wm.R11 + c_g * wm.R12 + b_g * ty * wm.R12 - a_c * c_g * wm.R13
                          - a_c * b_g * ty * wm.R13 + c_c * wm.R21 + a_g * c_c * ty * wm.R21 - wm.R22 - a_g * ty * wm.R22 + a_c * wm.R23
                          + a_c * a_g * ty * wm.R23 + a_g * c_c * wm.R31 - c_c * ty * wm.R31 - a_g * wm.R32 + ty * wm.R32
                          + a_c * a_g * wm.R33 - a_c * ty * wm.R33)
                           * (bx - x_c - x_g - x_l
                              - sz
                                  * (wm.R13 - c_g * wm.R23 + b_c * (wm.R11 - c_g * wm.R21 + b_g * wm.R31)
                                     + a_c * (wm.R12 - c_g * wm.R22 + b_g * wm.R32) + b_g * wm.R33)
                              - sv
                                  * (wm.R12 - c_g * wm.R22 - c_c * (wm.R11 - c_g * wm.R21 + b_g * wm.R31) + b_g * wm.R32
                                     - a_c * (wm.R13 - c_g * wm.R23 + b_g * wm.R33))
                              - su
                                  * (wm.R11 - c_g * wm.R21 + b_g * wm.R31 - c_c * (wm.R12 - c_g * wm.R22 + b_g * wm.R32)
                                     + b_c * (wm.R13 - c_g * wm.R23 + b_g * wm.R33))))));
        // drdag
        this->template set_global_derivative<3>() =
            (((-(sz * (-(b_c * wm.R21) - a_c * wm.R22 - wm.R23)) - sv * (c_c * wm.R21 - wm.R22 + a_c * wm.R23)
               - su * (-wm.R21 + c_c * wm.R22 - b_c * wm.R23))
                  * (c_c * c_g * tx * wm.R11 + c_c * ty * wm.R11 - c_g * tx * wm.R12 - ty * wm.R12 + a_c * c_g * tx * wm.R13
                     + a_c * ty * wm.R13 - c_c * tx * wm.R21 - c_c * c_g * ty * wm.R21 + tx * wm.R22 + c_g * ty * wm.R22 - a_c * tx * wm.R23
                     - a_c * c_g * ty * wm.R23 - a_g * c_c * tx * wm.R31 + b_g * c_c * ty * wm.R31 + a_g * tx * wm.R32 - b_g * ty * wm.R32
                     - a_c * a_g * tx * wm.R33 + a_c * b_g * ty * wm.R33)
              + (-(c_c * wm.R11) + b_g * c_c * tx * wm.R11 + wm.R12 - b_g * tx * wm.R12 - a_c * wm.R13 + a_c * b_g * tx * wm.R13
                 + c_c * c_g * wm.R21 - a_g * c_c * tx * wm.R21 - c_g * wm.R22 + a_g * tx * wm.R22 + a_c * c_g * wm.R23
                 - a_c * a_g * tx * wm.R23 - b_g * c_c * wm.R31 + c_c * tx * wm.R31 + b_g * wm.R32 - tx * wm.R32 - a_c * b_g * wm.R33
                 + a_c * tx * wm.R33)
                  * (-(sz * (b_c * wm.R31 + a_c * wm.R32 + wm.R33)) - sv * (-(c_c * wm.R31) + wm.R32 - a_c * wm.R33)
                     - su * (wm.R31 - c_c * wm.R32 + b_c * wm.R33))
              + (-(c_c * tx * wm.R31) + tx * wm.R32 - a_c * tx * wm.R33)
                  * (-z_c - z_g - z_l
                     - sz
                         * (b_g * wm.R13 - a_g * wm.R23 + b_c * (b_g * wm.R11 - a_g * wm.R21 + wm.R31)
                            + a_c * (b_g * wm.R12 - a_g * wm.R22 + wm.R32) + wm.R33)
                     - sv
                         * (b_g * wm.R12 - a_g * wm.R22 - c_c * (b_g * wm.R11 - a_g * wm.R21 + wm.R31) + wm.R32
                            - a_c * (b_g * wm.R13 - a_g * wm.R23 + wm.R33))
                     - su
                         * (b_g * wm.R11 - a_g * wm.R21 + wm.R31 - c_c * (b_g * wm.R12 - a_g * wm.R22 + wm.R32)
                            + b_c * (b_g * wm.R13 - a_g * wm.R23 + wm.R33)))
              + (-(c_c * tx * wm.R21) + tx * wm.R22 - a_c * tx * wm.R23)
                  * (by - y_c - y_g - y_l
                     - sz
                         * (-(c_g * wm.R13) + wm.R23 + b_c * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31)
                            + a_c * (-(c_g * wm.R12) + wm.R22 + a_g * wm.R32) + a_g * wm.R33)
                     - sv
                         * (-(c_g * wm.R12) + wm.R22 - c_c * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31) + a_g * wm.R32
                            - a_c * (-(c_g * wm.R13) + wm.R23 + a_g * wm.R33))
                     - su
                         * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31 - c_c * (-(c_g * wm.R12) + wm.R22 + a_g * wm.R32)
                            + b_c * (-(c_g * wm.R13) + wm.R23 + a_g * wm.R33)))
              + (c_c * ty * wm.R21 - ty * wm.R22 + a_c * ty * wm.R23 + c_c * wm.R31 - wm.R32 + a_c * wm.R33)
                  * (bx - x_c - x_g - x_l
                     - sz
                         * (wm.R13 - c_g * wm.R23 + b_c * (wm.R11 - c_g * wm.R21 + b_g * wm.R31)
                            + a_c * (wm.R12 - c_g * wm.R22 + b_g * wm.R32) + b_g * wm.R33)
                     - sv
                         * (wm.R12 - c_g * wm.R22 - c_c * (wm.R11 - c_g * wm.R21 + b_g * wm.R31) + b_g * wm.R32
                            - a_c * (wm.R13 - c_g * wm.R23 + b_g * wm.R33))
                     - su
                         * (wm.R11 - c_g * wm.R21 + b_g * wm.R31 - c_c * (wm.R12 - c_g * wm.R22 + b_g * wm.R32)
                            + b_c * (wm.R13 - c_g * wm.R23 + b_g * wm.R33))))
             * ((c_c * c_g * tx * wm.R11 + c_c * ty * wm.R11 - c_g * tx * wm.R12 - ty * wm.R12 + a_c * c_g * tx * wm.R13 + a_c * ty * wm.R13
                 - c_c * tx * wm.R21 - c_c * c_g * ty * wm.R21 + tx * wm.R22 + c_g * ty * wm.R22 - a_c * tx * wm.R23
                 - a_c * c_g * ty * wm.R23 - a_g * c_c * tx * wm.R31 + b_g * c_c * ty * wm.R31 + a_g * tx * wm.R32 - b_g * ty * wm.R32
                 - a_c * a_g * tx * wm.R33 + a_c * b_g * ty * wm.R33)
                    * (-z_c - z_g - z_l
                       - sz
                           * (b_g * wm.R13 - a_g * wm.R23 + b_c * (b_g * wm.R11 - a_g * wm.R21 + wm.R31)
                              + a_c * (b_g * wm.R12 - a_g * wm.R22 + wm.R32) + wm.R33)
                       - sv
                           * (b_g * wm.R12 - a_g * wm.R22 - c_c * (b_g * wm.R11 - a_g * wm.R21 + wm.R31) + wm.R32
                              - a_c * (b_g * wm.R13 - a_g * wm.R23 + wm.R33))
                       - su
                           * (b_g * wm.R11 - a_g * wm.R21 + wm.R31 - c_c * (b_g * wm.R12 - a_g * wm.R22 + wm.R32)
                              + b_c * (b_g * wm.R13 - a_g * wm.R23 + wm.R33)))
                + (-(c_c * wm.R11) + b_g * c_c * tx * wm.R11 + wm.R12 - b_g * tx * wm.R12 - a_c * wm.R13 + a_c * b_g * tx * wm.R13
                   + c_c * c_g * wm.R21 - a_g * c_c * tx * wm.R21 - c_g * wm.R22 + a_g * tx * wm.R22 + a_c * c_g * wm.R23
                   - a_c * a_g * tx * wm.R23 - b_g * c_c * wm.R31 + c_c * tx * wm.R31 + b_g * wm.R32 - tx * wm.R32 - a_c * b_g * wm.R33
                   + a_c * tx * wm.R33)
                    * (by - y_c - y_g - y_l
                       - sz
                           * (-(c_g * wm.R13) + wm.R23 + b_c * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31)
                              + a_c * (-(c_g * wm.R12) + wm.R22 + a_g * wm.R32) + a_g * wm.R33)
                       - sv
                           * (-(c_g * wm.R12) + wm.R22 - c_c * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31) + a_g * wm.R32
                              - a_c * (-(c_g * wm.R13) + wm.R23 + a_g * wm.R33))
                       - su
                           * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31 - c_c * (-(c_g * wm.R12) + wm.R22 + a_g * wm.R32)
                              + b_c * (-(c_g * wm.R13) + wm.R23 + a_g * wm.R33)))
                + (-(c_c * c_g * wm.R11) - b_g * c_c * ty * wm.R11 + c_g * wm.R12 + b_g * ty * wm.R12 - a_c * c_g * wm.R13
                   - a_c * b_g * ty * wm.R13 + c_c * wm.R21 + a_g * c_c * ty * wm.R21 - wm.R22 - a_g * ty * wm.R22 + a_c * wm.R23
                   + a_c * a_g * ty * wm.R23 + a_g * c_c * wm.R31 - c_c * ty * wm.R31 - a_g * wm.R32 + ty * wm.R32 + a_c * a_g * wm.R33
                   - a_c * ty * wm.R33)
                    * (bx - x_c - x_g - x_l
                       - sz
                           * (wm.R13 - c_g * wm.R23 + b_c * (wm.R11 - c_g * wm.R21 + b_g * wm.R31)
                              + a_c * (wm.R12 - c_g * wm.R22 + b_g * wm.R32) + b_g * wm.R33)
                       - sv
                           * (wm.R12 - c_g * wm.R22 - c_c * (wm.R11 - c_g * wm.R21 + b_g * wm.R31) + b_g * wm.R32
                              - a_c * (wm.R13 - c_g * wm.R23 + b_g * wm.R33))
                       - su
                           * (wm.R11 - c_g * wm.R21 + b_g * wm.R31 - c_c * (wm.R12 - c_g * wm.R22 + b_g * wm.R32)
                              + b_c * (wm.R13 - c_g * wm.R23 + b_g * wm.R33)))))
                / (sqrt(pow(-(c_c * wm.R11) + b_g * c_c * tx * wm.R11 + wm.R12 - b_g * tx * wm.R12 - a_c * wm.R13 + a_c * b_g * tx * wm.R13
                                + c_c * c_g * wm.R21 - a_g * c_c * tx * wm.R21 - c_g * wm.R22 + a_g * tx * wm.R22 + a_c * c_g * wm.R23
                                - a_c * a_g * tx * wm.R23 - b_g * c_c * wm.R31 + c_c * tx * wm.R31 + b_g * wm.R32 - tx * wm.R32
                                - a_c * b_g * wm.R33 + a_c * tx * wm.R33,
                            2)
                        + pow(-(c_c * c_g * wm.R11) - b_g * c_c * ty * wm.R11 + c_g * wm.R12 + b_g * ty * wm.R12 - a_c * c_g * wm.R13
                                  - a_c * b_g * ty * wm.R13 + c_c * wm.R21 + a_g * c_c * ty * wm.R21 - wm.R22 - a_g * ty * wm.R22
                                  + a_c * wm.R23 + a_c * a_g * ty * wm.R23 + a_g * c_c * wm.R31 - c_c * ty * wm.R31 - a_g * wm.R32
                                  + ty * wm.R32 + a_c * a_g * wm.R33 - a_c * ty * wm.R33,
                              2)
                        + pow(c_c * c_g * tx * wm.R11 + c_c * ty * wm.R11 - c_g * tx * wm.R12 - ty * wm.R12 + a_c * c_g * tx * wm.R13
                                  + a_c * ty * wm.R13 - c_c * tx * wm.R21 - c_c * c_g * ty * wm.R21 + tx * wm.R22 + c_g * ty * wm.R22
                                  - a_c * tx * wm.R23 - a_c * c_g * ty * wm.R23 - a_g * c_c * tx * wm.R31 + b_g * c_c * ty * wm.R31
                                  + a_g * tx * wm.R32 - b_g * ty * wm.R32 - a_c * a_g * tx * wm.R33 + a_c * b_g * ty * wm.R33,
                              2))
                   * fabs(((c_c * c_g * tx * wm.R11 + c_c * ty * wm.R11 - c_g * tx * wm.R12 - ty * wm.R12 + a_c * c_g * tx * wm.R13
                            + a_c * ty * wm.R13 - c_c * tx * wm.R21 - c_c * c_g * ty * wm.R21 + tx * wm.R22 + c_g * ty * wm.R22
                            - a_c * tx * wm.R23 - a_c * c_g * ty * wm.R23 - a_g * c_c * tx * wm.R31 + b_g * c_c * ty * wm.R31
                            + a_g * tx * wm.R32 - b_g * ty * wm.R32 - a_c * a_g * tx * wm.R33 + a_c * b_g * ty * wm.R33)
                               * (-z_c - z_g - z_l
                                  - sz
                                      * (b_g * wm.R13 - a_g * wm.R23 + b_c * (b_g * wm.R11 - a_g * wm.R21 + wm.R31)
                                         + a_c * (b_g * wm.R12 - a_g * wm.R22 + wm.R32) + wm.R33)
                                  - sv
                                      * (b_g * wm.R12 - a_g * wm.R22 - c_c * (b_g * wm.R11 - a_g * wm.R21 + wm.R31) + wm.R32
                                         - a_c * (b_g * wm.R13 - a_g * wm.R23 + wm.R33))
                                  - su
                                      * (b_g * wm.R11 - a_g * wm.R21 + wm.R31 - c_c * (b_g * wm.R12 - a_g * wm.R22 + wm.R32)
                                         + b_c * (b_g * wm.R13 - a_g * wm.R23 + wm.R33)))
                           + (-(c_c * wm.R11) + b_g * c_c * tx * wm.R11 + wm.R12 - b_g * tx * wm.R12 - a_c * wm.R13
                              + a_c * b_g * tx * wm.R13 + c_c * c_g * wm.R21 - a_g * c_c * tx * wm.R21 - c_g * wm.R22 + a_g * tx * wm.R22
                              + a_c * c_g * wm.R23 - a_c * a_g * tx * wm.R23 - b_g * c_c * wm.R31 + c_c * tx * wm.R31 + b_g * wm.R32
                              - tx * wm.R32 - a_c * b_g * wm.R33 + a_c * tx * wm.R33)
                               * (by - y_c - y_g - y_l
                                  - sz
                                      * (-(c_g * wm.R13) + wm.R23 + b_c * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31)
                                         + a_c * (-(c_g * wm.R12) + wm.R22 + a_g * wm.R32) + a_g * wm.R33)
                                  - sv
                                      * (-(c_g * wm.R12) + wm.R22 - c_c * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31) + a_g * wm.R32
                                         - a_c * (-(c_g * wm.R13) + wm.R23 + a_g * wm.R33))
                                  - su
                                      * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31 - c_c * (-(c_g * wm.R12) + wm.R22 + a_g * wm.R32)
                                         + b_c * (-(c_g * wm.R13) + wm.R23 + a_g * wm.R33)))
                           + (-(c_c * c_g * wm.R11) - b_g * c_c * ty * wm.R11 + c_g * wm.R12 + b_g * ty * wm.R12 - a_c * c_g * wm.R13
                              - a_c * b_g * ty * wm.R13 + c_c * wm.R21 + a_g * c_c * ty * wm.R21 - wm.R22 - a_g * ty * wm.R22 + a_c * wm.R23
                              + a_c * a_g * ty * wm.R23 + a_g * c_c * wm.R31 - c_c * ty * wm.R31 - a_g * wm.R32 + ty * wm.R32
                              + a_c * a_g * wm.R33 - a_c * ty * wm.R33)
                               * (bx - x_c - x_g - x_l
                                  - sz
                                      * (wm.R13 - c_g * wm.R23 + b_c * (wm.R11 - c_g * wm.R21 + b_g * wm.R31)
                                         + a_c * (wm.R12 - c_g * wm.R22 + b_g * wm.R32) + b_g * wm.R33)
                                  - sv
                                      * (wm.R12 - c_g * wm.R22 - c_c * (wm.R11 - c_g * wm.R21 + b_g * wm.R31) + b_g * wm.R32
                                         - a_c * (wm.R13 - c_g * wm.R23 + b_g * wm.R33))
                                  - su
                                      * (wm.R11 - c_g * wm.R21 + b_g * wm.R31 - c_c * (wm.R12 - c_g * wm.R22 + b_g * wm.R32)
                                         + b_c * (wm.R13 - c_g * wm.R23 + b_g * wm.R33))))))
            - ((2 * (-(c_c * tx * wm.R21) + tx * wm.R22 - a_c * tx * wm.R23)
                    * (-(c_c * wm.R11) + b_g * c_c * tx * wm.R11 + wm.R12 - b_g * tx * wm.R12 - a_c * wm.R13 + a_c * b_g * tx * wm.R13
                       + c_c * c_g * wm.R21 - a_g * c_c * tx * wm.R21 - c_g * wm.R22 + a_g * tx * wm.R22 + a_c * c_g * wm.R23
                       - a_c * a_g * tx * wm.R23 - b_g * c_c * wm.R31 + c_c * tx * wm.R31 + b_g * wm.R32 - tx * wm.R32 - a_c * b_g * wm.R33
                       + a_c * tx * wm.R33)
                + 2 * (c_c * ty * wm.R21 - ty * wm.R22 + a_c * ty * wm.R23 + c_c * wm.R31 - wm.R32 + a_c * wm.R33)
                    * (-(c_c * c_g * wm.R11) - b_g * c_c * ty * wm.R11 + c_g * wm.R12 + b_g * ty * wm.R12 - a_c * c_g * wm.R13
                       - a_c * b_g * ty * wm.R13 + c_c * wm.R21 + a_g * c_c * ty * wm.R21 - wm.R22 - a_g * ty * wm.R22 + a_c * wm.R23
                       + a_c * a_g * ty * wm.R23 + a_g * c_c * wm.R31 - c_c * ty * wm.R31 - a_g * wm.R32 + ty * wm.R32 + a_c * a_g * wm.R33
                       - a_c * ty * wm.R33)
                + 2 * (-(c_c * tx * wm.R31) + tx * wm.R32 - a_c * tx * wm.R33)
                    * (c_c * c_g * tx * wm.R11 + c_c * ty * wm.R11 - c_g * tx * wm.R12 - ty * wm.R12 + a_c * c_g * tx * wm.R13
                       + a_c * ty * wm.R13 - c_c * tx * wm.R21 - c_c * c_g * ty * wm.R21 + tx * wm.R22 + c_g * ty * wm.R22
                       - a_c * tx * wm.R23 - a_c * c_g * ty * wm.R23 - a_g * c_c * tx * wm.R31 + b_g * c_c * ty * wm.R31 + a_g * tx * wm.R32
                       - b_g * ty * wm.R32 - a_c * a_g * tx * wm.R33 + a_c * b_g * ty * wm.R33))
               * fabs(((c_c * c_g * tx * wm.R11 + c_c * ty * wm.R11 - c_g * tx * wm.R12 - ty * wm.R12 + a_c * c_g * tx * wm.R13
                        + a_c * ty * wm.R13 - c_c * tx * wm.R21 - c_c * c_g * ty * wm.R21 + tx * wm.R22 + c_g * ty * wm.R22
                        - a_c * tx * wm.R23 - a_c * c_g * ty * wm.R23 - a_g * c_c * tx * wm.R31 + b_g * c_c * ty * wm.R31
                        + a_g * tx * wm.R32 - b_g * ty * wm.R32 - a_c * a_g * tx * wm.R33 + a_c * b_g * ty * wm.R33)
                           * (-z_c - z_g - z_l
                              - sz
                                  * (b_g * wm.R13 - a_g * wm.R23 + b_c * (b_g * wm.R11 - a_g * wm.R21 + wm.R31)
                                     + a_c * (b_g * wm.R12 - a_g * wm.R22 + wm.R32) + wm.R33)
                              - sv
                                  * (b_g * wm.R12 - a_g * wm.R22 - c_c * (b_g * wm.R11 - a_g * wm.R21 + wm.R31) + wm.R32
                                     - a_c * (b_g * wm.R13 - a_g * wm.R23 + wm.R33))
                              - su
                                  * (b_g * wm.R11 - a_g * wm.R21 + wm.R31 - c_c * (b_g * wm.R12 - a_g * wm.R22 + wm.R32)
                                     + b_c * (b_g * wm.R13 - a_g * wm.R23 + wm.R33)))
                       + (-(c_c * wm.R11) + b_g * c_c * tx * wm.R11 + wm.R12 - b_g * tx * wm.R12 - a_c * wm.R13 + a_c * b_g * tx * wm.R13
                          + c_c * c_g * wm.R21 - a_g * c_c * tx * wm.R21 - c_g * wm.R22 + a_g * tx * wm.R22 + a_c * c_g * wm.R23
                          - a_c * a_g * tx * wm.R23 - b_g * c_c * wm.R31 + c_c * tx * wm.R31 + b_g * wm.R32 - tx * wm.R32
                          - a_c * b_g * wm.R33 + a_c * tx * wm.R33)
                           * (by - y_c - y_g - y_l
                              - sz
                                  * (-(c_g * wm.R13) + wm.R23 + b_c * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31)
                                     + a_c * (-(c_g * wm.R12) + wm.R22 + a_g * wm.R32) + a_g * wm.R33)
                              - sv
                                  * (-(c_g * wm.R12) + wm.R22 - c_c * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31) + a_g * wm.R32
                                     - a_c * (-(c_g * wm.R13) + wm.R23 + a_g * wm.R33))
                              - su
                                  * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31 - c_c * (-(c_g * wm.R12) + wm.R22 + a_g * wm.R32)
                                     + b_c * (-(c_g * wm.R13) + wm.R23 + a_g * wm.R33)))
                       + (-(c_c * c_g * wm.R11) - b_g * c_c * ty * wm.R11 + c_g * wm.R12 + b_g * ty * wm.R12 - a_c * c_g * wm.R13
                          - a_c * b_g * ty * wm.R13 + c_c * wm.R21 + a_g * c_c * ty * wm.R21 - wm.R22 - a_g * ty * wm.R22 + a_c * wm.R23
                          + a_c * a_g * ty * wm.R23 + a_g * c_c * wm.R31 - c_c * ty * wm.R31 - a_g * wm.R32 + ty * wm.R32
                          + a_c * a_g * wm.R33 - a_c * ty * wm.R33)
                           * (bx - x_c - x_g - x_l
                              - sz
                                  * (wm.R13 - c_g * wm.R23 + b_c * (wm.R11 - c_g * wm.R21 + b_g * wm.R31)
                                     + a_c * (wm.R12 - c_g * wm.R22 + b_g * wm.R32) + b_g * wm.R33)
                              - sv
                                  * (wm.R12 - c_g * wm.R22 - c_c * (wm.R11 - c_g * wm.R21 + b_g * wm.R31) + b_g * wm.R32
                                     - a_c * (wm.R13 - c_g * wm.R23 + b_g * wm.R33))
                              - su
                                  * (wm.R11 - c_g * wm.R21 + b_g * wm.R31 - c_c * (wm.R12 - c_g * wm.R22 + b_g * wm.R32)
                                     + b_c * (wm.R13 - c_g * wm.R23 + b_g * wm.R33))))))
                / (2
                   * pow(pow(-(c_c * wm.R11) + b_g * c_c * tx * wm.R11 + wm.R12 - b_g * tx * wm.R12 - a_c * wm.R13 + a_c * b_g * tx * wm.R13
                                 + c_c * c_g * wm.R21 - a_g * c_c * tx * wm.R21 - c_g * wm.R22 + a_g * tx * wm.R22 + a_c * c_g * wm.R23
                                 - a_c * a_g * tx * wm.R23 - b_g * c_c * wm.R31 + c_c * tx * wm.R31 + b_g * wm.R32 - tx * wm.R32
                                 - a_c * b_g * wm.R33 + a_c * tx * wm.R33,
                             2)
                             + pow(-(c_c * c_g * wm.R11) - b_g * c_c * ty * wm.R11 + c_g * wm.R12 + b_g * ty * wm.R12 - a_c * c_g * wm.R13
                                       - a_c * b_g * ty * wm.R13 + c_c * wm.R21 + a_g * c_c * ty * wm.R21 - wm.R22 - a_g * ty * wm.R22
                                       + a_c * wm.R23 + a_c * a_g * ty * wm.R23 + a_g * c_c * wm.R31 - c_c * ty * wm.R31 - a_g * wm.R32
                                       + ty * wm.R32 + a_c * a_g * wm.R33 - a_c * ty * wm.R33,
                                   2)
                             + pow(c_c * c_g * tx * wm.R11 + c_c * ty * wm.R11 - c_g * tx * wm.R12 - ty * wm.R12 + a_c * c_g * tx * wm.R13
                                       + a_c * ty * wm.R13 - c_c * tx * wm.R21 - c_c * c_g * ty * wm.R21 + tx * wm.R22 + c_g * ty * wm.R22
                                       - a_c * tx * wm.R23 - a_c * c_g * ty * wm.R23 - a_g * c_c * tx * wm.R31 + b_g * c_c * ty * wm.R31
                                       + a_g * tx * wm.R32 - b_g * ty * wm.R32 - a_c * a_g * tx * wm.R33 + a_c * b_g * ty * wm.R33,
                                   2),
                         3.0 / 2.0));
        // drdbg
        this->template set_global_derivative<4>() =
            (((-(sz * (b_c * wm.R11 + a_c * wm.R12 + wm.R13)) - sv * (-(c_c * wm.R11) + wm.R12 - a_c * wm.R13)
               - su * (wm.R11 - c_c * wm.R12 + b_c * wm.R13))
                  * (c_c * c_g * tx * wm.R11 + c_c * ty * wm.R11 - c_g * tx * wm.R12 - ty * wm.R12 + a_c * c_g * tx * wm.R13
                     + a_c * ty * wm.R13 - c_c * tx * wm.R21 - c_c * c_g * ty * wm.R21 + tx * wm.R22 + c_g * ty * wm.R22 - a_c * tx * wm.R23
                     - a_c * c_g * ty * wm.R23 - a_g * c_c * tx * wm.R31 + b_g * c_c * ty * wm.R31 + a_g * tx * wm.R32 - b_g * ty * wm.R32
                     - a_c * a_g * tx * wm.R33 + a_c * b_g * ty * wm.R33)
              + (-(c_c * c_g * wm.R11) - b_g * c_c * ty * wm.R11 + c_g * wm.R12 + b_g * ty * wm.R12 - a_c * c_g * wm.R13
                 - a_c * b_g * ty * wm.R13 + c_c * wm.R21 + a_g * c_c * ty * wm.R21 - wm.R22 - a_g * ty * wm.R22 + a_c * wm.R23
                 + a_c * a_g * ty * wm.R23 + a_g * c_c * wm.R31 - c_c * ty * wm.R31 - a_g * wm.R32 + ty * wm.R32 + a_c * a_g * wm.R33
                 - a_c * ty * wm.R33)
                  * (-(sz * (b_c * wm.R31 + a_c * wm.R32 + wm.R33)) - sv * (-(c_c * wm.R31) + wm.R32 - a_c * wm.R33)
                     - su * (wm.R31 - c_c * wm.R32 + b_c * wm.R33))
              + (c_c * ty * wm.R31 - ty * wm.R32 + a_c * ty * wm.R33)
                  * (-z_c - z_g - z_l
                     - sz
                         * (b_g * wm.R13 - a_g * wm.R23 + b_c * (b_g * wm.R11 - a_g * wm.R21 + wm.R31)
                            + a_c * (b_g * wm.R12 - a_g * wm.R22 + wm.R32) + wm.R33)
                     - sv
                         * (b_g * wm.R12 - a_g * wm.R22 - c_c * (b_g * wm.R11 - a_g * wm.R21 + wm.R31) + wm.R32
                            - a_c * (b_g * wm.R13 - a_g * wm.R23 + wm.R33))
                     - su
                         * (b_g * wm.R11 - a_g * wm.R21 + wm.R31 - c_c * (b_g * wm.R12 - a_g * wm.R22 + wm.R32)
                            + b_c * (b_g * wm.R13 - a_g * wm.R23 + wm.R33)))
              + (c_c * tx * wm.R11 - tx * wm.R12 + a_c * tx * wm.R13 - c_c * wm.R31 + wm.R32 - a_c * wm.R33)
                  * (by - y_c - y_g - y_l
                     - sz
                         * (-(c_g * wm.R13) + wm.R23 + b_c * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31)
                            + a_c * (-(c_g * wm.R12) + wm.R22 + a_g * wm.R32) + a_g * wm.R33)
                     - sv
                         * (-(c_g * wm.R12) + wm.R22 - c_c * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31) + a_g * wm.R32
                            - a_c * (-(c_g * wm.R13) + wm.R23 + a_g * wm.R33))
                     - su
                         * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31 - c_c * (-(c_g * wm.R12) + wm.R22 + a_g * wm.R32)
                            + b_c * (-(c_g * wm.R13) + wm.R23 + a_g * wm.R33)))
              + (-(c_c * ty * wm.R11) + ty * wm.R12 - a_c * ty * wm.R13)
                  * (bx - x_c - x_g - x_l
                     - sz
                         * (wm.R13 - c_g * wm.R23 + b_c * (wm.R11 - c_g * wm.R21 + b_g * wm.R31)
                            + a_c * (wm.R12 - c_g * wm.R22 + b_g * wm.R32) + b_g * wm.R33)
                     - sv
                         * (wm.R12 - c_g * wm.R22 - c_c * (wm.R11 - c_g * wm.R21 + b_g * wm.R31) + b_g * wm.R32
                            - a_c * (wm.R13 - c_g * wm.R23 + b_g * wm.R33))
                     - su
                         * (wm.R11 - c_g * wm.R21 + b_g * wm.R31 - c_c * (wm.R12 - c_g * wm.R22 + b_g * wm.R32)
                            + b_c * (wm.R13 - c_g * wm.R23 + b_g * wm.R33))))
             * ((c_c * c_g * tx * wm.R11 + c_c * ty * wm.R11 - c_g * tx * wm.R12 - ty * wm.R12 + a_c * c_g * tx * wm.R13 + a_c * ty * wm.R13
                 - c_c * tx * wm.R21 - c_c * c_g * ty * wm.R21 + tx * wm.R22 + c_g * ty * wm.R22 - a_c * tx * wm.R23
                 - a_c * c_g * ty * wm.R23 - a_g * c_c * tx * wm.R31 + b_g * c_c * ty * wm.R31 + a_g * tx * wm.R32 - b_g * ty * wm.R32
                 - a_c * a_g * tx * wm.R33 + a_c * b_g * ty * wm.R33)
                    * (-z_c - z_g - z_l
                       - sz
                           * (b_g * wm.R13 - a_g * wm.R23 + b_c * (b_g * wm.R11 - a_g * wm.R21 + wm.R31)
                              + a_c * (b_g * wm.R12 - a_g * wm.R22 + wm.R32) + wm.R33)
                       - sv
                           * (b_g * wm.R12 - a_g * wm.R22 - c_c * (b_g * wm.R11 - a_g * wm.R21 + wm.R31) + wm.R32
                              - a_c * (b_g * wm.R13 - a_g * wm.R23 + wm.R33))
                       - su
                           * (b_g * wm.R11 - a_g * wm.R21 + wm.R31 - c_c * (b_g * wm.R12 - a_g * wm.R22 + wm.R32)
                              + b_c * (b_g * wm.R13 - a_g * wm.R23 + wm.R33)))
                + (-(c_c * wm.R11) + b_g * c_c * tx * wm.R11 + wm.R12 - b_g * tx * wm.R12 - a_c * wm.R13 + a_c * b_g * tx * wm.R13
                   + c_c * c_g * wm.R21 - a_g * c_c * tx * wm.R21 - c_g * wm.R22 + a_g * tx * wm.R22 + a_c * c_g * wm.R23
                   - a_c * a_g * tx * wm.R23 - b_g * c_c * wm.R31 + c_c * tx * wm.R31 + b_g * wm.R32 - tx * wm.R32 - a_c * b_g * wm.R33
                   + a_c * tx * wm.R33)
                    * (by - y_c - y_g - y_l
                       - sz
                           * (-(c_g * wm.R13) + wm.R23 + b_c * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31)
                              + a_c * (-(c_g * wm.R12) + wm.R22 + a_g * wm.R32) + a_g * wm.R33)
                       - sv
                           * (-(c_g * wm.R12) + wm.R22 - c_c * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31) + a_g * wm.R32
                              - a_c * (-(c_g * wm.R13) + wm.R23 + a_g * wm.R33))
                       - su
                           * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31 - c_c * (-(c_g * wm.R12) + wm.R22 + a_g * wm.R32)
                              + b_c * (-(c_g * wm.R13) + wm.R23 + a_g * wm.R33)))
                + (-(c_c * c_g * wm.R11) - b_g * c_c * ty * wm.R11 + c_g * wm.R12 + b_g * ty * wm.R12 - a_c * c_g * wm.R13
                   - a_c * b_g * ty * wm.R13 + c_c * wm.R21 + a_g * c_c * ty * wm.R21 - wm.R22 - a_g * ty * wm.R22 + a_c * wm.R23
                   + a_c * a_g * ty * wm.R23 + a_g * c_c * wm.R31 - c_c * ty * wm.R31 - a_g * wm.R32 + ty * wm.R32 + a_c * a_g * wm.R33
                   - a_c * ty * wm.R33)
                    * (bx - x_c - x_g - x_l
                       - sz
                           * (wm.R13 - c_g * wm.R23 + b_c * (wm.R11 - c_g * wm.R21 + b_g * wm.R31)
                              + a_c * (wm.R12 - c_g * wm.R22 + b_g * wm.R32) + b_g * wm.R33)
                       - sv
                           * (wm.R12 - c_g * wm.R22 - c_c * (wm.R11 - c_g * wm.R21 + b_g * wm.R31) + b_g * wm.R32
                              - a_c * (wm.R13 - c_g * wm.R23 + b_g * wm.R33))
                       - su
                           * (wm.R11 - c_g * wm.R21 + b_g * wm.R31 - c_c * (wm.R12 - c_g * wm.R22 + b_g * wm.R32)
                              + b_c * (wm.R13 - c_g * wm.R23 + b_g * wm.R33)))))
                / (sqrt(pow(-(c_c * wm.R11) + b_g * c_c * tx * wm.R11 + wm.R12 - b_g * tx * wm.R12 - a_c * wm.R13 + a_c * b_g * tx * wm.R13
                                + c_c * c_g * wm.R21 - a_g * c_c * tx * wm.R21 - c_g * wm.R22 + a_g * tx * wm.R22 + a_c * c_g * wm.R23
                                - a_c * a_g * tx * wm.R23 - b_g * c_c * wm.R31 + c_c * tx * wm.R31 + b_g * wm.R32 - tx * wm.R32
                                - a_c * b_g * wm.R33 + a_c * tx * wm.R33,
                            2)
                        + pow(-(c_c * c_g * wm.R11) - b_g * c_c * ty * wm.R11 + c_g * wm.R12 + b_g * ty * wm.R12 - a_c * c_g * wm.R13
                                  - a_c * b_g * ty * wm.R13 + c_c * wm.R21 + a_g * c_c * ty * wm.R21 - wm.R22 - a_g * ty * wm.R22
                                  + a_c * wm.R23 + a_c * a_g * ty * wm.R23 + a_g * c_c * wm.R31 - c_c * ty * wm.R31 - a_g * wm.R32
                                  + ty * wm.R32 + a_c * a_g * wm.R33 - a_c * ty * wm.R33,
                              2)
                        + pow(c_c * c_g * tx * wm.R11 + c_c * ty * wm.R11 - c_g * tx * wm.R12 - ty * wm.R12 + a_c * c_g * tx * wm.R13
                                  + a_c * ty * wm.R13 - c_c * tx * wm.R21 - c_c * c_g * ty * wm.R21 + tx * wm.R22 + c_g * ty * wm.R22
                                  - a_c * tx * wm.R23 - a_c * c_g * ty * wm.R23 - a_g * c_c * tx * wm.R31 + b_g * c_c * ty * wm.R31
                                  + a_g * tx * wm.R32 - b_g * ty * wm.R32 - a_c * a_g * tx * wm.R33 + a_c * b_g * ty * wm.R33,
                              2))
                   * fabs(((c_c * c_g * tx * wm.R11 + c_c * ty * wm.R11 - c_g * tx * wm.R12 - ty * wm.R12 + a_c * c_g * tx * wm.R13
                            + a_c * ty * wm.R13 - c_c * tx * wm.R21 - c_c * c_g * ty * wm.R21 + tx * wm.R22 + c_g * ty * wm.R22
                            - a_c * tx * wm.R23 - a_c * c_g * ty * wm.R23 - a_g * c_c * tx * wm.R31 + b_g * c_c * ty * wm.R31
                            + a_g * tx * wm.R32 - b_g * ty * wm.R32 - a_c * a_g * tx * wm.R33 + a_c * b_g * ty * wm.R33)
                               * (-z_c - z_g - z_l
                                  - sz
                                      * (b_g * wm.R13 - a_g * wm.R23 + b_c * (b_g * wm.R11 - a_g * wm.R21 + wm.R31)
                                         + a_c * (b_g * wm.R12 - a_g * wm.R22 + wm.R32) + wm.R33)
                                  - sv
                                      * (b_g * wm.R12 - a_g * wm.R22 - c_c * (b_g * wm.R11 - a_g * wm.R21 + wm.R31) + wm.R32
                                         - a_c * (b_g * wm.R13 - a_g * wm.R23 + wm.R33))
                                  - su
                                      * (b_g * wm.R11 - a_g * wm.R21 + wm.R31 - c_c * (b_g * wm.R12 - a_g * wm.R22 + wm.R32)
                                         + b_c * (b_g * wm.R13 - a_g * wm.R23 + wm.R33)))
                           + (-(c_c * wm.R11) + b_g * c_c * tx * wm.R11 + wm.R12 - b_g * tx * wm.R12 - a_c * wm.R13
                              + a_c * b_g * tx * wm.R13 + c_c * c_g * wm.R21 - a_g * c_c * tx * wm.R21 - c_g * wm.R22 + a_g * tx * wm.R22
                              + a_c * c_g * wm.R23 - a_c * a_g * tx * wm.R23 - b_g * c_c * wm.R31 + c_c * tx * wm.R31 + b_g * wm.R32
                              - tx * wm.R32 - a_c * b_g * wm.R33 + a_c * tx * wm.R33)
                               * (by - y_c - y_g - y_l
                                  - sz
                                      * (-(c_g * wm.R13) + wm.R23 + b_c * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31)
                                         + a_c * (-(c_g * wm.R12) + wm.R22 + a_g * wm.R32) + a_g * wm.R33)
                                  - sv
                                      * (-(c_g * wm.R12) + wm.R22 - c_c * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31) + a_g * wm.R32
                                         - a_c * (-(c_g * wm.R13) + wm.R23 + a_g * wm.R33))
                                  - su
                                      * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31 - c_c * (-(c_g * wm.R12) + wm.R22 + a_g * wm.R32)
                                         + b_c * (-(c_g * wm.R13) + wm.R23 + a_g * wm.R33)))
                           + (-(c_c * c_g * wm.R11) - b_g * c_c * ty * wm.R11 + c_g * wm.R12 + b_g * ty * wm.R12 - a_c * c_g * wm.R13
                              - a_c * b_g * ty * wm.R13 + c_c * wm.R21 + a_g * c_c * ty * wm.R21 - wm.R22 - a_g * ty * wm.R22 + a_c * wm.R23
                              + a_c * a_g * ty * wm.R23 + a_g * c_c * wm.R31 - c_c * ty * wm.R31 - a_g * wm.R32 + ty * wm.R32
                              + a_c * a_g * wm.R33 - a_c * ty * wm.R33)
                               * (bx - x_c - x_g - x_l
                                  - sz
                                      * (wm.R13 - c_g * wm.R23 + b_c * (wm.R11 - c_g * wm.R21 + b_g * wm.R31)
                                         + a_c * (wm.R12 - c_g * wm.R22 + b_g * wm.R32) + b_g * wm.R33)
                                  - sv
                                      * (wm.R12 - c_g * wm.R22 - c_c * (wm.R11 - c_g * wm.R21 + b_g * wm.R31) + b_g * wm.R32
                                         - a_c * (wm.R13 - c_g * wm.R23 + b_g * wm.R33))
                                  - su
                                      * (wm.R11 - c_g * wm.R21 + b_g * wm.R31 - c_c * (wm.R12 - c_g * wm.R22 + b_g * wm.R32)
                                         + b_c * (wm.R13 - c_g * wm.R23 + b_g * wm.R33))))))
            - ((2 * (c_c * tx * wm.R11 - tx * wm.R12 + a_c * tx * wm.R13 - c_c * wm.R31 + wm.R32 - a_c * wm.R33)
                    * (-(c_c * wm.R11) + b_g * c_c * tx * wm.R11 + wm.R12 - b_g * tx * wm.R12 - a_c * wm.R13 + a_c * b_g * tx * wm.R13
                       + c_c * c_g * wm.R21 - a_g * c_c * tx * wm.R21 - c_g * wm.R22 + a_g * tx * wm.R22 + a_c * c_g * wm.R23
                       - a_c * a_g * tx * wm.R23 - b_g * c_c * wm.R31 + c_c * tx * wm.R31 + b_g * wm.R32 - tx * wm.R32 - a_c * b_g * wm.R33
                       + a_c * tx * wm.R33)
                + 2 * (-(c_c * ty * wm.R11) + ty * wm.R12 - a_c * ty * wm.R13)
                    * (-(c_c * c_g * wm.R11) - b_g * c_c * ty * wm.R11 + c_g * wm.R12 + b_g * ty * wm.R12 - a_c * c_g * wm.R13
                       - a_c * b_g * ty * wm.R13 + c_c * wm.R21 + a_g * c_c * ty * wm.R21 - wm.R22 - a_g * ty * wm.R22 + a_c * wm.R23
                       + a_c * a_g * ty * wm.R23 + a_g * c_c * wm.R31 - c_c * ty * wm.R31 - a_g * wm.R32 + ty * wm.R32 + a_c * a_g * wm.R33
                       - a_c * ty * wm.R33)
                + 2 * (c_c * ty * wm.R31 - ty * wm.R32 + a_c * ty * wm.R33)
                    * (c_c * c_g * tx * wm.R11 + c_c * ty * wm.R11 - c_g * tx * wm.R12 - ty * wm.R12 + a_c * c_g * tx * wm.R13
                       + a_c * ty * wm.R13 - c_c * tx * wm.R21 - c_c * c_g * ty * wm.R21 + tx * wm.R22 + c_g * ty * wm.R22
                       - a_c * tx * wm.R23 - a_c * c_g * ty * wm.R23 - a_g * c_c * tx * wm.R31 + b_g * c_c * ty * wm.R31 + a_g * tx * wm.R32
                       - b_g * ty * wm.R32 - a_c * a_g * tx * wm.R33 + a_c * b_g * ty * wm.R33))
               * fabs(((c_c * c_g * tx * wm.R11 + c_c * ty * wm.R11 - c_g * tx * wm.R12 - ty * wm.R12 + a_c * c_g * tx * wm.R13
                        + a_c * ty * wm.R13 - c_c * tx * wm.R21 - c_c * c_g * ty * wm.R21 + tx * wm.R22 + c_g * ty * wm.R22
                        - a_c * tx * wm.R23 - a_c * c_g * ty * wm.R23 - a_g * c_c * tx * wm.R31 + b_g * c_c * ty * wm.R31
                        + a_g * tx * wm.R32 - b_g * ty * wm.R32 - a_c * a_g * tx * wm.R33 + a_c * b_g * ty * wm.R33)
                           * (-z_c - z_g - z_l
                              - sz
                                  * (b_g * wm.R13 - a_g * wm.R23 + b_c * (b_g * wm.R11 - a_g * wm.R21 + wm.R31)
                                     + a_c * (b_g * wm.R12 - a_g * wm.R22 + wm.R32) + wm.R33)
                              - sv
                                  * (b_g * wm.R12 - a_g * wm.R22 - c_c * (b_g * wm.R11 - a_g * wm.R21 + wm.R31) + wm.R32
                                     - a_c * (b_g * wm.R13 - a_g * wm.R23 + wm.R33))
                              - su
                                  * (b_g * wm.R11 - a_g * wm.R21 + wm.R31 - c_c * (b_g * wm.R12 - a_g * wm.R22 + wm.R32)
                                     + b_c * (b_g * wm.R13 - a_g * wm.R23 + wm.R33)))
                       + (-(c_c * wm.R11) + b_g * c_c * tx * wm.R11 + wm.R12 - b_g * tx * wm.R12 - a_c * wm.R13 + a_c * b_g * tx * wm.R13
                          + c_c * c_g * wm.R21 - a_g * c_c * tx * wm.R21 - c_g * wm.R22 + a_g * tx * wm.R22 + a_c * c_g * wm.R23
                          - a_c * a_g * tx * wm.R23 - b_g * c_c * wm.R31 + c_c * tx * wm.R31 + b_g * wm.R32 - tx * wm.R32
                          - a_c * b_g * wm.R33 + a_c * tx * wm.R33)
                           * (by - y_c - y_g - y_l
                              - sz
                                  * (-(c_g * wm.R13) + wm.R23 + b_c * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31)
                                     + a_c * (-(c_g * wm.R12) + wm.R22 + a_g * wm.R32) + a_g * wm.R33)
                              - sv
                                  * (-(c_g * wm.R12) + wm.R22 - c_c * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31) + a_g * wm.R32
                                     - a_c * (-(c_g * wm.R13) + wm.R23 + a_g * wm.R33))
                              - su
                                  * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31 - c_c * (-(c_g * wm.R12) + wm.R22 + a_g * wm.R32)
                                     + b_c * (-(c_g * wm.R13) + wm.R23 + a_g * wm.R33)))
                       + (-(c_c * c_g * wm.R11) - b_g * c_c * ty * wm.R11 + c_g * wm.R12 + b_g * ty * wm.R12 - a_c * c_g * wm.R13
                          - a_c * b_g * ty * wm.R13 + c_c * wm.R21 + a_g * c_c * ty * wm.R21 - wm.R22 - a_g * ty * wm.R22 + a_c * wm.R23
                          + a_c * a_g * ty * wm.R23 + a_g * c_c * wm.R31 - c_c * ty * wm.R31 - a_g * wm.R32 + ty * wm.R32
                          + a_c * a_g * wm.R33 - a_c * ty * wm.R33)
                           * (bx - x_c - x_g - x_l
                              - sz
                                  * (wm.R13 - c_g * wm.R23 + b_c * (wm.R11 - c_g * wm.R21 + b_g * wm.R31)
                                     + a_c * (wm.R12 - c_g * wm.R22 + b_g * wm.R32) + b_g * wm.R33)
                              - sv
                                  * (wm.R12 - c_g * wm.R22 - c_c * (wm.R11 - c_g * wm.R21 + b_g * wm.R31) + b_g * wm.R32
                                     - a_c * (wm.R13 - c_g * wm.R23 + b_g * wm.R33))
                              - su
                                  * (wm.R11 - c_g * wm.R21 + b_g * wm.R31 - c_c * (wm.R12 - c_g * wm.R22 + b_g * wm.R32)
                                     + b_c * (wm.R13 - c_g * wm.R23 + b_g * wm.R33))))))
                / (2
                   * pow(pow(-(c_c * wm.R11) + b_g * c_c * tx * wm.R11 + wm.R12 - b_g * tx * wm.R12 - a_c * wm.R13 + a_c * b_g * tx * wm.R13
                                 + c_c * c_g * wm.R21 - a_g * c_c * tx * wm.R21 - c_g * wm.R22 + a_g * tx * wm.R22 + a_c * c_g * wm.R23
                                 - a_c * a_g * tx * wm.R23 - b_g * c_c * wm.R31 + c_c * tx * wm.R31 + b_g * wm.R32 - tx * wm.R32
                                 - a_c * b_g * wm.R33 + a_c * tx * wm.R33,
                             2)
                             + pow(-(c_c * c_g * wm.R11) - b_g * c_c * ty * wm.R11 + c_g * wm.R12 + b_g * ty * wm.R12 - a_c * c_g * wm.R13
                                       - a_c * b_g * ty * wm.R13 + c_c * wm.R21 + a_g * c_c * ty * wm.R21 - wm.R22 - a_g * ty * wm.R22
                                       + a_c * wm.R23 + a_c * a_g * ty * wm.R23 + a_g * c_c * wm.R31 - c_c * ty * wm.R31 - a_g * wm.R32
                                       + ty * wm.R32 + a_c * a_g * wm.R33 - a_c * ty * wm.R33,
                                   2)
                             + pow(c_c * c_g * tx * wm.R11 + c_c * ty * wm.R11 - c_g * tx * wm.R12 - ty * wm.R12 + a_c * c_g * tx * wm.R13
                                       + a_c * ty * wm.R13 - c_c * tx * wm.R21 - c_c * c_g * ty * wm.R21 + tx * wm.R22 + c_g * ty * wm.R22
                                       - a_c * tx * wm.R23 - a_c * c_g * ty * wm.R23 - a_g * c_c * tx * wm.R31 + b_g * c_c * ty * wm.R31
                                       + a_g * tx * wm.R32 - b_g * ty * wm.R32 - a_c * a_g * tx * wm.R33 + a_c * b_g * ty * wm.R33,
                                   2),
                         3.0 / 2.0));
        // drdcg
        this->template set_global_derivative<5>() =
            (((-(sz * (-(b_c * wm.R11) - a_c * wm.R12 - wm.R13)) - sv * (c_c * wm.R11 - wm.R12 + a_c * wm.R13)
               - su * (-wm.R11 + c_c * wm.R12 - b_c * wm.R13))
                  * (-(c_c * wm.R11) + b_g * c_c * tx * wm.R11 + wm.R12 - b_g * tx * wm.R12 - a_c * wm.R13 + a_c * b_g * tx * wm.R13
                     + c_c * c_g * wm.R21 - a_g * c_c * tx * wm.R21 - c_g * wm.R22 + a_g * tx * wm.R22 + a_c * c_g * wm.R23
                     - a_c * a_g * tx * wm.R23 - b_g * c_c * wm.R31 + c_c * tx * wm.R31 + b_g * wm.R32 - tx * wm.R32 - a_c * b_g * wm.R33
                     + a_c * tx * wm.R33)
              + (-(sz * (-(b_c * wm.R21) - a_c * wm.R22 - wm.R23)) - sv * (c_c * wm.R21 - wm.R22 + a_c * wm.R23)
                 - su * (-wm.R21 + c_c * wm.R22 - b_c * wm.R23))
                  * (-(c_c * c_g * wm.R11) - b_g * c_c * ty * wm.R11 + c_g * wm.R12 + b_g * ty * wm.R12 - a_c * c_g * wm.R13
                     - a_c * b_g * ty * wm.R13 + c_c * wm.R21 + a_g * c_c * ty * wm.R21 - wm.R22 - a_g * ty * wm.R22 + a_c * wm.R23
                     + a_c * a_g * ty * wm.R23 + a_g * c_c * wm.R31 - c_c * ty * wm.R31 - a_g * wm.R32 + ty * wm.R32 + a_c * a_g * wm.R33
                     - a_c * ty * wm.R33)
              + (c_c * tx * wm.R11 - tx * wm.R12 + a_c * tx * wm.R13 - c_c * ty * wm.R21 + ty * wm.R22 - a_c * ty * wm.R23)
                  * (-z_c - z_g - z_l
                     - sz
                         * (b_g * wm.R13 - a_g * wm.R23 + b_c * (b_g * wm.R11 - a_g * wm.R21 + wm.R31)
                            + a_c * (b_g * wm.R12 - a_g * wm.R22 + wm.R32) + wm.R33)
                     - sv
                         * (b_g * wm.R12 - a_g * wm.R22 - c_c * (b_g * wm.R11 - a_g * wm.R21 + wm.R31) + wm.R32
                            - a_c * (b_g * wm.R13 - a_g * wm.R23 + wm.R33))
                     - su
                         * (b_g * wm.R11 - a_g * wm.R21 + wm.R31 - c_c * (b_g * wm.R12 - a_g * wm.R22 + wm.R32)
                            + b_c * (b_g * wm.R13 - a_g * wm.R23 + wm.R33)))
              + (c_c * wm.R21 - wm.R22 + a_c * wm.R23)
                  * (by - y_c - y_g - y_l
                     - sz
                         * (-(c_g * wm.R13) + wm.R23 + b_c * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31)
                            + a_c * (-(c_g * wm.R12) + wm.R22 + a_g * wm.R32) + a_g * wm.R33)
                     - sv
                         * (-(c_g * wm.R12) + wm.R22 - c_c * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31) + a_g * wm.R32
                            - a_c * (-(c_g * wm.R13) + wm.R23 + a_g * wm.R33))
                     - su
                         * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31 - c_c * (-(c_g * wm.R12) + wm.R22 + a_g * wm.R32)
                            + b_c * (-(c_g * wm.R13) + wm.R23 + a_g * wm.R33)))
              + (-(c_c * wm.R11) + wm.R12 - a_c * wm.R13)
                  * (bx - x_c - x_g - x_l
                     - sz
                         * (wm.R13 - c_g * wm.R23 + b_c * (wm.R11 - c_g * wm.R21 + b_g * wm.R31)
                            + a_c * (wm.R12 - c_g * wm.R22 + b_g * wm.R32) + b_g * wm.R33)
                     - sv
                         * (wm.R12 - c_g * wm.R22 - c_c * (wm.R11 - c_g * wm.R21 + b_g * wm.R31) + b_g * wm.R32
                            - a_c * (wm.R13 - c_g * wm.R23 + b_g * wm.R33))
                     - su
                         * (wm.R11 - c_g * wm.R21 + b_g * wm.R31 - c_c * (wm.R12 - c_g * wm.R22 + b_g * wm.R32)
                            + b_c * (wm.R13 - c_g * wm.R23 + b_g * wm.R33))))
             * ((c_c * c_g * tx * wm.R11 + c_c * ty * wm.R11 - c_g * tx * wm.R12 - ty * wm.R12 + a_c * c_g * tx * wm.R13 + a_c * ty * wm.R13
                 - c_c * tx * wm.R21 - c_c * c_g * ty * wm.R21 + tx * wm.R22 + c_g * ty * wm.R22 - a_c * tx * wm.R23
                 - a_c * c_g * ty * wm.R23 - a_g * c_c * tx * wm.R31 + b_g * c_c * ty * wm.R31 + a_g * tx * wm.R32 - b_g * ty * wm.R32
                 - a_c * a_g * tx * wm.R33 + a_c * b_g * ty * wm.R33)
                    * (-z_c - z_g - z_l
                       - sz
                           * (b_g * wm.R13 - a_g * wm.R23 + b_c * (b_g * wm.R11 - a_g * wm.R21 + wm.R31)
                              + a_c * (b_g * wm.R12 - a_g * wm.R22 + wm.R32) + wm.R33)
                       - sv
                           * (b_g * wm.R12 - a_g * wm.R22 - c_c * (b_g * wm.R11 - a_g * wm.R21 + wm.R31) + wm.R32
                              - a_c * (b_g * wm.R13 - a_g * wm.R23 + wm.R33))
                       - su
                           * (b_g * wm.R11 - a_g * wm.R21 + wm.R31 - c_c * (b_g * wm.R12 - a_g * wm.R22 + wm.R32)
                              + b_c * (b_g * wm.R13 - a_g * wm.R23 + wm.R33)))
                + (-(c_c * wm.R11) + b_g * c_c * tx * wm.R11 + wm.R12 - b_g * tx * wm.R12 - a_c * wm.R13 + a_c * b_g * tx * wm.R13
                   + c_c * c_g * wm.R21 - a_g * c_c * tx * wm.R21 - c_g * wm.R22 + a_g * tx * wm.R22 + a_c * c_g * wm.R23
                   - a_c * a_g * tx * wm.R23 - b_g * c_c * wm.R31 + c_c * tx * wm.R31 + b_g * wm.R32 - tx * wm.R32 - a_c * b_g * wm.R33
                   + a_c * tx * wm.R33)
                    * (by - y_c - y_g - y_l
                       - sz
                           * (-(c_g * wm.R13) + wm.R23 + b_c * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31)
                              + a_c * (-(c_g * wm.R12) + wm.R22 + a_g * wm.R32) + a_g * wm.R33)
                       - sv
                           * (-(c_g * wm.R12) + wm.R22 - c_c * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31) + a_g * wm.R32
                              - a_c * (-(c_g * wm.R13) + wm.R23 + a_g * wm.R33))
                       - su
                           * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31 - c_c * (-(c_g * wm.R12) + wm.R22 + a_g * wm.R32)
                              + b_c * (-(c_g * wm.R13) + wm.R23 + a_g * wm.R33)))
                + (-(c_c * c_g * wm.R11) - b_g * c_c * ty * wm.R11 + c_g * wm.R12 + b_g * ty * wm.R12 - a_c * c_g * wm.R13
                   - a_c * b_g * ty * wm.R13 + c_c * wm.R21 + a_g * c_c * ty * wm.R21 - wm.R22 - a_g * ty * wm.R22 + a_c * wm.R23
                   + a_c * a_g * ty * wm.R23 + a_g * c_c * wm.R31 - c_c * ty * wm.R31 - a_g * wm.R32 + ty * wm.R32 + a_c * a_g * wm.R33
                   - a_c * ty * wm.R33)
                    * (bx - x_c - x_g - x_l
                       - sz
                           * (wm.R13 - c_g * wm.R23 + b_c * (wm.R11 - c_g * wm.R21 + b_g * wm.R31)
                              + a_c * (wm.R12 - c_g * wm.R22 + b_g * wm.R32) + b_g * wm.R33)
                       - sv
                           * (wm.R12 - c_g * wm.R22 - c_c * (wm.R11 - c_g * wm.R21 + b_g * wm.R31) + b_g * wm.R32
                              - a_c * (wm.R13 - c_g * wm.R23 + b_g * wm.R33))
                       - su
                           * (wm.R11 - c_g * wm.R21 + b_g * wm.R31 - c_c * (wm.R12 - c_g * wm.R22 + b_g * wm.R32)
                              + b_c * (wm.R13 - c_g * wm.R23 + b_g * wm.R33)))))
                / (sqrt(pow(-(c_c * wm.R11) + b_g * c_c * tx * wm.R11 + wm.R12 - b_g * tx * wm.R12 - a_c * wm.R13 + a_c * b_g * tx * wm.R13
                                + c_c * c_g * wm.R21 - a_g * c_c * tx * wm.R21 - c_g * wm.R22 + a_g * tx * wm.R22 + a_c * c_g * wm.R23
                                - a_c * a_g * tx * wm.R23 - b_g * c_c * wm.R31 + c_c * tx * wm.R31 + b_g * wm.R32 - tx * wm.R32
                                - a_c * b_g * wm.R33 + a_c * tx * wm.R33,
                            2)
                        + pow(-(c_c * c_g * wm.R11) - b_g * c_c * ty * wm.R11 + c_g * wm.R12 + b_g * ty * wm.R12 - a_c * c_g * wm.R13
                                  - a_c * b_g * ty * wm.R13 + c_c * wm.R21 + a_g * c_c * ty * wm.R21 - wm.R22 - a_g * ty * wm.R22
                                  + a_c * wm.R23 + a_c * a_g * ty * wm.R23 + a_g * c_c * wm.R31 - c_c * ty * wm.R31 - a_g * wm.R32
                                  + ty * wm.R32 + a_c * a_g * wm.R33 - a_c * ty * wm.R33,
                              2)
                        + pow(c_c * c_g * tx * wm.R11 + c_c * ty * wm.R11 - c_g * tx * wm.R12 - ty * wm.R12 + a_c * c_g * tx * wm.R13
                                  + a_c * ty * wm.R13 - c_c * tx * wm.R21 - c_c * c_g * ty * wm.R21 + tx * wm.R22 + c_g * ty * wm.R22
                                  - a_c * tx * wm.R23 - a_c * c_g * ty * wm.R23 - a_g * c_c * tx * wm.R31 + b_g * c_c * ty * wm.R31
                                  + a_g * tx * wm.R32 - b_g * ty * wm.R32 - a_c * a_g * tx * wm.R33 + a_c * b_g * ty * wm.R33,
                              2))
                   * fabs(((c_c * c_g * tx * wm.R11 + c_c * ty * wm.R11 - c_g * tx * wm.R12 - ty * wm.R12 + a_c * c_g * tx * wm.R13
                            + a_c * ty * wm.R13 - c_c * tx * wm.R21 - c_c * c_g * ty * wm.R21 + tx * wm.R22 + c_g * ty * wm.R22
                            - a_c * tx * wm.R23 - a_c * c_g * ty * wm.R23 - a_g * c_c * tx * wm.R31 + b_g * c_c * ty * wm.R31
                            + a_g * tx * wm.R32 - b_g * ty * wm.R32 - a_c * a_g * tx * wm.R33 + a_c * b_g * ty * wm.R33)
                               * (-z_c - z_g - z_l
                                  - sz
                                      * (b_g * wm.R13 - a_g * wm.R23 + b_c * (b_g * wm.R11 - a_g * wm.R21 + wm.R31)
                                         + a_c * (b_g * wm.R12 - a_g * wm.R22 + wm.R32) + wm.R33)
                                  - sv
                                      * (b_g * wm.R12 - a_g * wm.R22 - c_c * (b_g * wm.R11 - a_g * wm.R21 + wm.R31) + wm.R32
                                         - a_c * (b_g * wm.R13 - a_g * wm.R23 + wm.R33))
                                  - su
                                      * (b_g * wm.R11 - a_g * wm.R21 + wm.R31 - c_c * (b_g * wm.R12 - a_g * wm.R22 + wm.R32)
                                         + b_c * (b_g * wm.R13 - a_g * wm.R23 + wm.R33)))
                           + (-(c_c * wm.R11) + b_g * c_c * tx * wm.R11 + wm.R12 - b_g * tx * wm.R12 - a_c * wm.R13
                              + a_c * b_g * tx * wm.R13 + c_c * c_g * wm.R21 - a_g * c_c * tx * wm.R21 - c_g * wm.R22 + a_g * tx * wm.R22
                              + a_c * c_g * wm.R23 - a_c * a_g * tx * wm.R23 - b_g * c_c * wm.R31 + c_c * tx * wm.R31 + b_g * wm.R32
                              - tx * wm.R32 - a_c * b_g * wm.R33 + a_c * tx * wm.R33)
                               * (by - y_c - y_g - y_l
                                  - sz
                                      * (-(c_g * wm.R13) + wm.R23 + b_c * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31)
                                         + a_c * (-(c_g * wm.R12) + wm.R22 + a_g * wm.R32) + a_g * wm.R33)
                                  - sv
                                      * (-(c_g * wm.R12) + wm.R22 - c_c * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31) + a_g * wm.R32
                                         - a_c * (-(c_g * wm.R13) + wm.R23 + a_g * wm.R33))
                                  - su
                                      * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31 - c_c * (-(c_g * wm.R12) + wm.R22 + a_g * wm.R32)
                                         + b_c * (-(c_g * wm.R13) + wm.R23 + a_g * wm.R33)))
                           + (-(c_c * c_g * wm.R11) - b_g * c_c * ty * wm.R11 + c_g * wm.R12 + b_g * ty * wm.R12 - a_c * c_g * wm.R13
                              - a_c * b_g * ty * wm.R13 + c_c * wm.R21 + a_g * c_c * ty * wm.R21 - wm.R22 - a_g * ty * wm.R22 + a_c * wm.R23
                              + a_c * a_g * ty * wm.R23 + a_g * c_c * wm.R31 - c_c * ty * wm.R31 - a_g * wm.R32 + ty * wm.R32
                              + a_c * a_g * wm.R33 - a_c * ty * wm.R33)
                               * (bx - x_c - x_g - x_l
                                  - sz
                                      * (wm.R13 - c_g * wm.R23 + b_c * (wm.R11 - c_g * wm.R21 + b_g * wm.R31)
                                         + a_c * (wm.R12 - c_g * wm.R22 + b_g * wm.R32) + b_g * wm.R33)
                                  - sv
                                      * (wm.R12 - c_g * wm.R22 - c_c * (wm.R11 - c_g * wm.R21 + b_g * wm.R31) + b_g * wm.R32
                                         - a_c * (wm.R13 - c_g * wm.R23 + b_g * wm.R33))
                                  - su
                                      * (wm.R11 - c_g * wm.R21 + b_g * wm.R31 - c_c * (wm.R12 - c_g * wm.R22 + b_g * wm.R32)
                                         + b_c * (wm.R13 - c_g * wm.R23 + b_g * wm.R33))))))
            - ((2 * (c_c * wm.R21 - wm.R22 + a_c * wm.R23)
                    * (-(c_c * wm.R11) + b_g * c_c * tx * wm.R11 + wm.R12 - b_g * tx * wm.R12 - a_c * wm.R13 + a_c * b_g * tx * wm.R13
                       + c_c * c_g * wm.R21 - a_g * c_c * tx * wm.R21 - c_g * wm.R22 + a_g * tx * wm.R22 + a_c * c_g * wm.R23
                       - a_c * a_g * tx * wm.R23 - b_g * c_c * wm.R31 + c_c * tx * wm.R31 + b_g * wm.R32 - tx * wm.R32 - a_c * b_g * wm.R33
                       + a_c * tx * wm.R33)
                + 2 * (-(c_c * wm.R11) + wm.R12 - a_c * wm.R13)
                    * (-(c_c * c_g * wm.R11) - b_g * c_c * ty * wm.R11 + c_g * wm.R12 + b_g * ty * wm.R12 - a_c * c_g * wm.R13
                       - a_c * b_g * ty * wm.R13 + c_c * wm.R21 + a_g * c_c * ty * wm.R21 - wm.R22 - a_g * ty * wm.R22 + a_c * wm.R23
                       + a_c * a_g * ty * wm.R23 + a_g * c_c * wm.R31 - c_c * ty * wm.R31 - a_g * wm.R32 + ty * wm.R32 + a_c * a_g * wm.R33
                       - a_c * ty * wm.R33)
                + 2 * (c_c * tx * wm.R11 - tx * wm.R12 + a_c * tx * wm.R13 - c_c * ty * wm.R21 + ty * wm.R22 - a_c * ty * wm.R23)
                    * (c_c * c_g * tx * wm.R11 + c_c * ty * wm.R11 - c_g * tx * wm.R12 - ty * wm.R12 + a_c * c_g * tx * wm.R13
                       + a_c * ty * wm.R13 - c_c * tx * wm.R21 - c_c * c_g * ty * wm.R21 + tx * wm.R22 + c_g * ty * wm.R22
                       - a_c * tx * wm.R23 - a_c * c_g * ty * wm.R23 - a_g * c_c * tx * wm.R31 + b_g * c_c * ty * wm.R31 + a_g * tx * wm.R32
                       - b_g * ty * wm.R32 - a_c * a_g * tx * wm.R33 + a_c * b_g * ty * wm.R33))
               * fabs(((c_c * c_g * tx * wm.R11 + c_c * ty * wm.R11 - c_g * tx * wm.R12 - ty * wm.R12 + a_c * c_g * tx * wm.R13
                        + a_c * ty * wm.R13 - c_c * tx * wm.R21 - c_c * c_g * ty * wm.R21 + tx * wm.R22 + c_g * ty * wm.R22
                        - a_c * tx * wm.R23 - a_c * c_g * ty * wm.R23 - a_g * c_c * tx * wm.R31 + b_g * c_c * ty * wm.R31
                        + a_g * tx * wm.R32 - b_g * ty * wm.R32 - a_c * a_g * tx * wm.R33 + a_c * b_g * ty * wm.R33)
                           * (-z_c - z_g - z_l
                              - sz
                                  * (b_g * wm.R13 - a_g * wm.R23 + b_c * (b_g * wm.R11 - a_g * wm.R21 + wm.R31)
                                     + a_c * (b_g * wm.R12 - a_g * wm.R22 + wm.R32) + wm.R33)
                              - sv
                                  * (b_g * wm.R12 - a_g * wm.R22 - c_c * (b_g * wm.R11 - a_g * wm.R21 + wm.R31) + wm.R32
                                     - a_c * (b_g * wm.R13 - a_g * wm.R23 + wm.R33))
                              - su
                                  * (b_g * wm.R11 - a_g * wm.R21 + wm.R31 - c_c * (b_g * wm.R12 - a_g * wm.R22 + wm.R32)
                                     + b_c * (b_g * wm.R13 - a_g * wm.R23 + wm.R33)))
                       + (-(c_c * wm.R11) + b_g * c_c * tx * wm.R11 + wm.R12 - b_g * tx * wm.R12 - a_c * wm.R13 + a_c * b_g * tx * wm.R13
                          + c_c * c_g * wm.R21 - a_g * c_c * tx * wm.R21 - c_g * wm.R22 + a_g * tx * wm.R22 + a_c * c_g * wm.R23
                          - a_c * a_g * tx * wm.R23 - b_g * c_c * wm.R31 + c_c * tx * wm.R31 + b_g * wm.R32 - tx * wm.R32
                          - a_c * b_g * wm.R33 + a_c * tx * wm.R33)
                           * (by - y_c - y_g - y_l
                              - sz
                                  * (-(c_g * wm.R13) + wm.R23 + b_c * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31)
                                     + a_c * (-(c_g * wm.R12) + wm.R22 + a_g * wm.R32) + a_g * wm.R33)
                              - sv
                                  * (-(c_g * wm.R12) + wm.R22 - c_c * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31) + a_g * wm.R32
                                     - a_c * (-(c_g * wm.R13) + wm.R23 + a_g * wm.R33))
                              - su
                                  * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31 - c_c * (-(c_g * wm.R12) + wm.R22 + a_g * wm.R32)
                                     + b_c * (-(c_g * wm.R13) + wm.R23 + a_g * wm.R33)))
                       + (-(c_c * c_g * wm.R11) - b_g * c_c * ty * wm.R11 + c_g * wm.R12 + b_g * ty * wm.R12 - a_c * c_g * wm.R13
                          - a_c * b_g * ty * wm.R13 + c_c * wm.R21 + a_g * c_c * ty * wm.R21 - wm.R22 - a_g * ty * wm.R22 + a_c * wm.R23
                          + a_c * a_g * ty * wm.R23 + a_g * c_c * wm.R31 - c_c * ty * wm.R31 - a_g * wm.R32 + ty * wm.R32
                          + a_c * a_g * wm.R33 - a_c * ty * wm.R33)
                           * (bx - x_c - x_g - x_l
                              - sz
                                  * (wm.R13 - c_g * wm.R23 + b_c * (wm.R11 - c_g * wm.R21 + b_g * wm.R31)
                                     + a_c * (wm.R12 - c_g * wm.R22 + b_g * wm.R32) + b_g * wm.R33)
                              - sv
                                  * (wm.R12 - c_g * wm.R22 - c_c * (wm.R11 - c_g * wm.R21 + b_g * wm.R31) + b_g * wm.R32
                                     - a_c * (wm.R13 - c_g * wm.R23 + b_g * wm.R33))
                              - su
                                  * (wm.R11 - c_g * wm.R21 + b_g * wm.R31 - c_c * (wm.R12 - c_g * wm.R22 + b_g * wm.R32)
                                     + b_c * (wm.R13 - c_g * wm.R23 + b_g * wm.R33))))))
                / (2
                   * pow(pow(-(c_c * wm.R11) + b_g * c_c * tx * wm.R11 + wm.R12 - b_g * tx * wm.R12 - a_c * wm.R13 + a_c * b_g * tx * wm.R13
                                 + c_c * c_g * wm.R21 - a_g * c_c * tx * wm.R21 - c_g * wm.R22 + a_g * tx * wm.R22 + a_c * c_g * wm.R23
                                 - a_c * a_g * tx * wm.R23 - b_g * c_c * wm.R31 + c_c * tx * wm.R31 + b_g * wm.R32 - tx * wm.R32
                                 - a_c * b_g * wm.R33 + a_c * tx * wm.R33,
                             2)
                             + pow(-(c_c * c_g * wm.R11) - b_g * c_c * ty * wm.R11 + c_g * wm.R12 + b_g * ty * wm.R12 - a_c * c_g * wm.R13
                                       - a_c * b_g * ty * wm.R13 + c_c * wm.R21 + a_g * c_c * ty * wm.R21 - wm.R22 - a_g * ty * wm.R22
                                       + a_c * wm.R23 + a_c * a_g * ty * wm.R23 + a_g * c_c * wm.R31 - c_c * ty * wm.R31 - a_g * wm.R32
                                       + ty * wm.R32 + a_c * a_g * wm.R33 - a_c * ty * wm.R33,
                                   2)
                             + pow(c_c * c_g * tx * wm.R11 + c_c * ty * wm.R11 - c_g * tx * wm.R12 - ty * wm.R12 + a_c * c_g * tx * wm.R13
                                       + a_c * ty * wm.R13 - c_c * tx * wm.R21 - c_c * c_g * ty * wm.R21 + tx * wm.R22 + c_g * ty * wm.R22
                                       - a_c * tx * wm.R23 - a_c * c_g * ty * wm.R23 - a_g * c_c * tx * wm.R31 + b_g * c_c * ty * wm.R31
                                       + a_g * tx * wm.R32 - b_g * ty * wm.R32 - a_c * a_g * tx * wm.R33 + a_c * b_g * ty * wm.R33,
                                   2),
                         3.0 / 2.0));
        // drdxc
        this->template set_global_derivative<6>() =
            ((c_c * c_g * wm.R11 + b_g * c_c * ty * wm.R11 - c_g * wm.R12 - b_g * ty * wm.R12 + a_c * c_g * wm.R13 + a_c * b_g * ty * wm.R13
              - c_c * wm.R21 - a_g * c_c * ty * wm.R21 + wm.R22 + a_g * ty * wm.R22 - a_c * wm.R23 - a_c * a_g * ty * wm.R23
              - a_g * c_c * wm.R31 + c_c * ty * wm.R31 + a_g * wm.R32 - ty * wm.R32 - a_c * a_g * wm.R33 + a_c * ty * wm.R33)
             * ((c_c * c_g * tx * wm.R11 + c_c * ty * wm.R11 - c_g * tx * wm.R12 - ty * wm.R12 + a_c * c_g * tx * wm.R13 + a_c * ty * wm.R13
                 - c_c * tx * wm.R21 - c_c * c_g * ty * wm.R21 + tx * wm.R22 + c_g * ty * wm.R22 - a_c * tx * wm.R23
                 - a_c * c_g * ty * wm.R23 - a_g * c_c * tx * wm.R31 + b_g * c_c * ty * wm.R31 + a_g * tx * wm.R32 - b_g * ty * wm.R32
                 - a_c * a_g * tx * wm.R33 + a_c * b_g * ty * wm.R33)
                    * (-z_c - z_g - z_l
                       - sz
                           * (b_g * wm.R13 - a_g * wm.R23 + b_c * (b_g * wm.R11 - a_g * wm.R21 + wm.R31)
                              + a_c * (b_g * wm.R12 - a_g * wm.R22 + wm.R32) + wm.R33)
                       - sv
                           * (b_g * wm.R12 - a_g * wm.R22 - c_c * (b_g * wm.R11 - a_g * wm.R21 + wm.R31) + wm.R32
                              - a_c * (b_g * wm.R13 - a_g * wm.R23 + wm.R33))
                       - su
                           * (b_g * wm.R11 - a_g * wm.R21 + wm.R31 - c_c * (b_g * wm.R12 - a_g * wm.R22 + wm.R32)
                              + b_c * (b_g * wm.R13 - a_g * wm.R23 + wm.R33)))
                + (-(c_c * wm.R11) + b_g * c_c * tx * wm.R11 + wm.R12 - b_g * tx * wm.R12 - a_c * wm.R13 + a_c * b_g * tx * wm.R13
                   + c_c * c_g * wm.R21 - a_g * c_c * tx * wm.R21 - c_g * wm.R22 + a_g * tx * wm.R22 + a_c * c_g * wm.R23
                   - a_c * a_g * tx * wm.R23 - b_g * c_c * wm.R31 + c_c * tx * wm.R31 + b_g * wm.R32 - tx * wm.R32 - a_c * b_g * wm.R33
                   + a_c * tx * wm.R33)
                    * (by - y_c - y_g - y_l
                       - sz
                           * (-(c_g * wm.R13) + wm.R23 + b_c * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31)
                              + a_c * (-(c_g * wm.R12) + wm.R22 + a_g * wm.R32) + a_g * wm.R33)
                       - sv
                           * (-(c_g * wm.R12) + wm.R22 - c_c * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31) + a_g * wm.R32
                              - a_c * (-(c_g * wm.R13) + wm.R23 + a_g * wm.R33))
                       - su
                           * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31 - c_c * (-(c_g * wm.R12) + wm.R22 + a_g * wm.R32)
                              + b_c * (-(c_g * wm.R13) + wm.R23 + a_g * wm.R33)))
                + (-(c_c * c_g * wm.R11) - b_g * c_c * ty * wm.R11 + c_g * wm.R12 + b_g * ty * wm.R12 - a_c * c_g * wm.R13
                   - a_c * b_g * ty * wm.R13 + c_c * wm.R21 + a_g * c_c * ty * wm.R21 - wm.R22 - a_g * ty * wm.R22 + a_c * wm.R23
                   + a_c * a_g * ty * wm.R23 + a_g * c_c * wm.R31 - c_c * ty * wm.R31 - a_g * wm.R32 + ty * wm.R32 + a_c * a_g * wm.R33
                   - a_c * ty * wm.R33)
                    * (bx - x_c - x_g - x_l
                       - sz
                           * (wm.R13 - c_g * wm.R23 + b_c * (wm.R11 - c_g * wm.R21 + b_g * wm.R31)
                              + a_c * (wm.R12 - c_g * wm.R22 + b_g * wm.R32) + b_g * wm.R33)
                       - sv
                           * (wm.R12 - c_g * wm.R22 - c_c * (wm.R11 - c_g * wm.R21 + b_g * wm.R31) + b_g * wm.R32
                              - a_c * (wm.R13 - c_g * wm.R23 + b_g * wm.R33))
                       - su
                           * (wm.R11 - c_g * wm.R21 + b_g * wm.R31 - c_c * (wm.R12 - c_g * wm.R22 + b_g * wm.R32)
                              + b_c * (wm.R13 - c_g * wm.R23 + b_g * wm.R33)))))
            / (sqrt(pow(-(c_c * wm.R11) + b_g * c_c * tx * wm.R11 + wm.R12 - b_g * tx * wm.R12 - a_c * wm.R13 + a_c * b_g * tx * wm.R13
                            + c_c * c_g * wm.R21 - a_g * c_c * tx * wm.R21 - c_g * wm.R22 + a_g * tx * wm.R22 + a_c * c_g * wm.R23
                            - a_c * a_g * tx * wm.R23 - b_g * c_c * wm.R31 + c_c * tx * wm.R31 + b_g * wm.R32 - tx * wm.R32
                            - a_c * b_g * wm.R33 + a_c * tx * wm.R33,
                        2)
                    + pow(-(c_c * c_g * wm.R11) - b_g * c_c * ty * wm.R11 + c_g * wm.R12 + b_g * ty * wm.R12 - a_c * c_g * wm.R13
                              - a_c * b_g * ty * wm.R13 + c_c * wm.R21 + a_g * c_c * ty * wm.R21 - wm.R22 - a_g * ty * wm.R22 + a_c * wm.R23
                              + a_c * a_g * ty * wm.R23 + a_g * c_c * wm.R31 - c_c * ty * wm.R31 - a_g * wm.R32 + ty * wm.R32
                              + a_c * a_g * wm.R33 - a_c * ty * wm.R33,
                          2)
                    + pow(c_c * c_g * tx * wm.R11 + c_c * ty * wm.R11 - c_g * tx * wm.R12 - ty * wm.R12 + a_c * c_g * tx * wm.R13
                              + a_c * ty * wm.R13 - c_c * tx * wm.R21 - c_c * c_g * ty * wm.R21 + tx * wm.R22 + c_g * ty * wm.R22
                              - a_c * tx * wm.R23 - a_c * c_g * ty * wm.R23 - a_g * c_c * tx * wm.R31 + b_g * c_c * ty * wm.R31
                              + a_g * tx * wm.R32 - b_g * ty * wm.R32 - a_c * a_g * tx * wm.R33 + a_c * b_g * ty * wm.R33,
                          2))
               * fabs(((c_c * c_g * tx * wm.R11 + c_c * ty * wm.R11 - c_g * tx * wm.R12 - ty * wm.R12 + a_c * c_g * tx * wm.R13
                        + a_c * ty * wm.R13 - c_c * tx * wm.R21 - c_c * c_g * ty * wm.R21 + tx * wm.R22 + c_g * ty * wm.R22
                        - a_c * tx * wm.R23 - a_c * c_g * ty * wm.R23 - a_g * c_c * tx * wm.R31 + b_g * c_c * ty * wm.R31
                        + a_g * tx * wm.R32 - b_g * ty * wm.R32 - a_c * a_g * tx * wm.R33 + a_c * b_g * ty * wm.R33)
                           * (-z_c - z_g - z_l
                              - sz
                                  * (b_g * wm.R13 - a_g * wm.R23 + b_c * (b_g * wm.R11 - a_g * wm.R21 + wm.R31)
                                     + a_c * (b_g * wm.R12 - a_g * wm.R22 + wm.R32) + wm.R33)
                              - sv
                                  * (b_g * wm.R12 - a_g * wm.R22 - c_c * (b_g * wm.R11 - a_g * wm.R21 + wm.R31) + wm.R32
                                     - a_c * (b_g * wm.R13 - a_g * wm.R23 + wm.R33))
                              - su
                                  * (b_g * wm.R11 - a_g * wm.R21 + wm.R31 - c_c * (b_g * wm.R12 - a_g * wm.R22 + wm.R32)
                                     + b_c * (b_g * wm.R13 - a_g * wm.R23 + wm.R33)))
                       + (-(c_c * wm.R11) + b_g * c_c * tx * wm.R11 + wm.R12 - b_g * tx * wm.R12 - a_c * wm.R13 + a_c * b_g * tx * wm.R13
                          + c_c * c_g * wm.R21 - a_g * c_c * tx * wm.R21 - c_g * wm.R22 + a_g * tx * wm.R22 + a_c * c_g * wm.R23
                          - a_c * a_g * tx * wm.R23 - b_g * c_c * wm.R31 + c_c * tx * wm.R31 + b_g * wm.R32 - tx * wm.R32
                          - a_c * b_g * wm.R33 + a_c * tx * wm.R33)
                           * (by - y_c - y_g - y_l
                              - sz
                                  * (-(c_g * wm.R13) + wm.R23 + b_c * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31)
                                     + a_c * (-(c_g * wm.R12) + wm.R22 + a_g * wm.R32) + a_g * wm.R33)
                              - sv
                                  * (-(c_g * wm.R12) + wm.R22 - c_c * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31) + a_g * wm.R32
                                     - a_c * (-(c_g * wm.R13) + wm.R23 + a_g * wm.R33))
                              - su
                                  * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31 - c_c * (-(c_g * wm.R12) + wm.R22 + a_g * wm.R32)
                                     + b_c * (-(c_g * wm.R13) + wm.R23 + a_g * wm.R33)))
                       + (-(c_c * c_g * wm.R11) - b_g * c_c * ty * wm.R11 + c_g * wm.R12 + b_g * ty * wm.R12 - a_c * c_g * wm.R13
                          - a_c * b_g * ty * wm.R13 + c_c * wm.R21 + a_g * c_c * ty * wm.R21 - wm.R22 - a_g * ty * wm.R22 + a_c * wm.R23
                          + a_c * a_g * ty * wm.R23 + a_g * c_c * wm.R31 - c_c * ty * wm.R31 - a_g * wm.R32 + ty * wm.R32
                          + a_c * a_g * wm.R33 - a_c * ty * wm.R33)
                           * (bx - x_c - x_g - x_l
                              - sz
                                  * (wm.R13 - c_g * wm.R23 + b_c * (wm.R11 - c_g * wm.R21 + b_g * wm.R31)
                                     + a_c * (wm.R12 - c_g * wm.R22 + b_g * wm.R32) + b_g * wm.R33)
                              - sv
                                  * (wm.R12 - c_g * wm.R22 - c_c * (wm.R11 - c_g * wm.R21 + b_g * wm.R31) + b_g * wm.R32
                                     - a_c * (wm.R13 - c_g * wm.R23 + b_g * wm.R33))
                              - su
                                  * (wm.R11 - c_g * wm.R21 + b_g * wm.R31 - c_c * (wm.R12 - c_g * wm.R22 + b_g * wm.R32)
                                     + b_c * (wm.R13 - c_g * wm.R23 + b_g * wm.R33))))));
        // drdyc
        this->template set_global_derivative<7>() =
            ((c_c * wm.R11 - b_g * c_c * tx * wm.R11 - wm.R12 + b_g * tx * wm.R12 + a_c * wm.R13 - a_c * b_g * tx * wm.R13
              - c_c * c_g * wm.R21 + a_g * c_c * tx * wm.R21 + c_g * wm.R22 - a_g * tx * wm.R22 - a_c * c_g * wm.R23
              + a_c * a_g * tx * wm.R23 + b_g * c_c * wm.R31 - c_c * tx * wm.R31 - b_g * wm.R32 + tx * wm.R32 + a_c * b_g * wm.R33
              - a_c * tx * wm.R33)
             * ((c_c * c_g * tx * wm.R11 + c_c * ty * wm.R11 - c_g * tx * wm.R12 - ty * wm.R12 + a_c * c_g * tx * wm.R13 + a_c * ty * wm.R13
                 - c_c * tx * wm.R21 - c_c * c_g * ty * wm.R21 + tx * wm.R22 + c_g * ty * wm.R22 - a_c * tx * wm.R23
                 - a_c * c_g * ty * wm.R23 - a_g * c_c * tx * wm.R31 + b_g * c_c * ty * wm.R31 + a_g * tx * wm.R32 - b_g * ty * wm.R32
                 - a_c * a_g * tx * wm.R33 + a_c * b_g * ty * wm.R33)
                    * (-z_c - z_g - z_l
                       - sz
                           * (b_g * wm.R13 - a_g * wm.R23 + b_c * (b_g * wm.R11 - a_g * wm.R21 + wm.R31)
                              + a_c * (b_g * wm.R12 - a_g * wm.R22 + wm.R32) + wm.R33)
                       - sv
                           * (b_g * wm.R12 - a_g * wm.R22 - c_c * (b_g * wm.R11 - a_g * wm.R21 + wm.R31) + wm.R32
                              - a_c * (b_g * wm.R13 - a_g * wm.R23 + wm.R33))
                       - su
                           * (b_g * wm.R11 - a_g * wm.R21 + wm.R31 - c_c * (b_g * wm.R12 - a_g * wm.R22 + wm.R32)
                              + b_c * (b_g * wm.R13 - a_g * wm.R23 + wm.R33)))
                + (-(c_c * wm.R11) + b_g * c_c * tx * wm.R11 + wm.R12 - b_g * tx * wm.R12 - a_c * wm.R13 + a_c * b_g * tx * wm.R13
                   + c_c * c_g * wm.R21 - a_g * c_c * tx * wm.R21 - c_g * wm.R22 + a_g * tx * wm.R22 + a_c * c_g * wm.R23
                   - a_c * a_g * tx * wm.R23 - b_g * c_c * wm.R31 + c_c * tx * wm.R31 + b_g * wm.R32 - tx * wm.R32 - a_c * b_g * wm.R33
                   + a_c * tx * wm.R33)
                    * (by - y_c - y_g - y_l
                       - sz
                           * (-(c_g * wm.R13) + wm.R23 + b_c * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31)
                              + a_c * (-(c_g * wm.R12) + wm.R22 + a_g * wm.R32) + a_g * wm.R33)
                       - sv
                           * (-(c_g * wm.R12) + wm.R22 - c_c * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31) + a_g * wm.R32
                              - a_c * (-(c_g * wm.R13) + wm.R23 + a_g * wm.R33))
                       - su
                           * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31 - c_c * (-(c_g * wm.R12) + wm.R22 + a_g * wm.R32)
                              + b_c * (-(c_g * wm.R13) + wm.R23 + a_g * wm.R33)))
                + (-(c_c * c_g * wm.R11) - b_g * c_c * ty * wm.R11 + c_g * wm.R12 + b_g * ty * wm.R12 - a_c * c_g * wm.R13
                   - a_c * b_g * ty * wm.R13 + c_c * wm.R21 + a_g * c_c * ty * wm.R21 - wm.R22 - a_g * ty * wm.R22 + a_c * wm.R23
                   + a_c * a_g * ty * wm.R23 + a_g * c_c * wm.R31 - c_c * ty * wm.R31 - a_g * wm.R32 + ty * wm.R32 + a_c * a_g * wm.R33
                   - a_c * ty * wm.R33)
                    * (bx - x_c - x_g - x_l
                       - sz
                           * (wm.R13 - c_g * wm.R23 + b_c * (wm.R11 - c_g * wm.R21 + b_g * wm.R31)
                              + a_c * (wm.R12 - c_g * wm.R22 + b_g * wm.R32) + b_g * wm.R33)
                       - sv
                           * (wm.R12 - c_g * wm.R22 - c_c * (wm.R11 - c_g * wm.R21 + b_g * wm.R31) + b_g * wm.R32
                              - a_c * (wm.R13 - c_g * wm.R23 + b_g * wm.R33))
                       - su
                           * (wm.R11 - c_g * wm.R21 + b_g * wm.R31 - c_c * (wm.R12 - c_g * wm.R22 + b_g * wm.R32)
                              + b_c * (wm.R13 - c_g * wm.R23 + b_g * wm.R33)))))
            / (sqrt(pow(-(c_c * wm.R11) + b_g * c_c * tx * wm.R11 + wm.R12 - b_g * tx * wm.R12 - a_c * wm.R13 + a_c * b_g * tx * wm.R13
                            + c_c * c_g * wm.R21 - a_g * c_c * tx * wm.R21 - c_g * wm.R22 + a_g * tx * wm.R22 + a_c * c_g * wm.R23
                            - a_c * a_g * tx * wm.R23 - b_g * c_c * wm.R31 + c_c * tx * wm.R31 + b_g * wm.R32 - tx * wm.R32
                            - a_c * b_g * wm.R33 + a_c * tx * wm.R33,
                        2)
                    + pow(-(c_c * c_g * wm.R11) - b_g * c_c * ty * wm.R11 + c_g * wm.R12 + b_g * ty * wm.R12 - a_c * c_g * wm.R13
                              - a_c * b_g * ty * wm.R13 + c_c * wm.R21 + a_g * c_c * ty * wm.R21 - wm.R22 - a_g * ty * wm.R22 + a_c * wm.R23
                              + a_c * a_g * ty * wm.R23 + a_g * c_c * wm.R31 - c_c * ty * wm.R31 - a_g * wm.R32 + ty * wm.R32
                              + a_c * a_g * wm.R33 - a_c * ty * wm.R33,
                          2)
                    + pow(c_c * c_g * tx * wm.R11 + c_c * ty * wm.R11 - c_g * tx * wm.R12 - ty * wm.R12 + a_c * c_g * tx * wm.R13
                              + a_c * ty * wm.R13 - c_c * tx * wm.R21 - c_c * c_g * ty * wm.R21 + tx * wm.R22 + c_g * ty * wm.R22
                              - a_c * tx * wm.R23 - a_c * c_g * ty * wm.R23 - a_g * c_c * tx * wm.R31 + b_g * c_c * ty * wm.R31
                              + a_g * tx * wm.R32 - b_g * ty * wm.R32 - a_c * a_g * tx * wm.R33 + a_c * b_g * ty * wm.R33,
                          2))
               * fabs(((c_c * c_g * tx * wm.R11 + c_c * ty * wm.R11 - c_g * tx * wm.R12 - ty * wm.R12 + a_c * c_g * tx * wm.R13
                        + a_c * ty * wm.R13 - c_c * tx * wm.R21 - c_c * c_g * ty * wm.R21 + tx * wm.R22 + c_g * ty * wm.R22
                        - a_c * tx * wm.R23 - a_c * c_g * ty * wm.R23 - a_g * c_c * tx * wm.R31 + b_g * c_c * ty * wm.R31
                        + a_g * tx * wm.R32 - b_g * ty * wm.R32 - a_c * a_g * tx * wm.R33 + a_c * b_g * ty * wm.R33)
                           * (-z_c - z_g - z_l
                              - sz
                                  * (b_g * wm.R13 - a_g * wm.R23 + b_c * (b_g * wm.R11 - a_g * wm.R21 + wm.R31)
                                     + a_c * (b_g * wm.R12 - a_g * wm.R22 + wm.R32) + wm.R33)
                              - sv
                                  * (b_g * wm.R12 - a_g * wm.R22 - c_c * (b_g * wm.R11 - a_g * wm.R21 + wm.R31) + wm.R32
                                     - a_c * (b_g * wm.R13 - a_g * wm.R23 + wm.R33))
                              - su
                                  * (b_g * wm.R11 - a_g * wm.R21 + wm.R31 - c_c * (b_g * wm.R12 - a_g * wm.R22 + wm.R32)
                                     + b_c * (b_g * wm.R13 - a_g * wm.R23 + wm.R33)))
                       + (-(c_c * wm.R11) + b_g * c_c * tx * wm.R11 + wm.R12 - b_g * tx * wm.R12 - a_c * wm.R13 + a_c * b_g * tx * wm.R13
                          + c_c * c_g * wm.R21 - a_g * c_c * tx * wm.R21 - c_g * wm.R22 + a_g * tx * wm.R22 + a_c * c_g * wm.R23
                          - a_c * a_g * tx * wm.R23 - b_g * c_c * wm.R31 + c_c * tx * wm.R31 + b_g * wm.R32 - tx * wm.R32
                          - a_c * b_g * wm.R33 + a_c * tx * wm.R33)
                           * (by - y_c - y_g - y_l
                              - sz
                                  * (-(c_g * wm.R13) + wm.R23 + b_c * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31)
                                     + a_c * (-(c_g * wm.R12) + wm.R22 + a_g * wm.R32) + a_g * wm.R33)
                              - sv
                                  * (-(c_g * wm.R12) + wm.R22 - c_c * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31) + a_g * wm.R32
                                     - a_c * (-(c_g * wm.R13) + wm.R23 + a_g * wm.R33))
                              - su
                                  * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31 - c_c * (-(c_g * wm.R12) + wm.R22 + a_g * wm.R32)
                                     + b_c * (-(c_g * wm.R13) + wm.R23 + a_g * wm.R33)))
                       + (-(c_c * c_g * wm.R11) - b_g * c_c * ty * wm.R11 + c_g * wm.R12 + b_g * ty * wm.R12 - a_c * c_g * wm.R13
                          - a_c * b_g * ty * wm.R13 + c_c * wm.R21 + a_g * c_c * ty * wm.R21 - wm.R22 - a_g * ty * wm.R22 + a_c * wm.R23
                          + a_c * a_g * ty * wm.R23 + a_g * c_c * wm.R31 - c_c * ty * wm.R31 - a_g * wm.R32 + ty * wm.R32
                          + a_c * a_g * wm.R33 - a_c * ty * wm.R33)
                           * (bx - x_c - x_g - x_l
                              - sz
                                  * (wm.R13 - c_g * wm.R23 + b_c * (wm.R11 - c_g * wm.R21 + b_g * wm.R31)
                                     + a_c * (wm.R12 - c_g * wm.R22 + b_g * wm.R32) + b_g * wm.R33)
                              - sv
                                  * (wm.R12 - c_g * wm.R22 - c_c * (wm.R11 - c_g * wm.R21 + b_g * wm.R31) + b_g * wm.R32
                                     - a_c * (wm.R13 - c_g * wm.R23 + b_g * wm.R33))
                              - su
                                  * (wm.R11 - c_g * wm.R21 + b_g * wm.R31 - c_c * (wm.R12 - c_g * wm.R22 + b_g * wm.R32)
                                     + b_c * (wm.R13 - c_g * wm.R23 + b_g * wm.R33))))));
        // drdzc
        this->template set_global_derivative<8>() =
            ((-(c_c * c_g * tx * wm.R11) - c_c * ty * wm.R11 + c_g * tx * wm.R12 + ty * wm.R12 - a_c * c_g * tx * wm.R13 - a_c * ty * wm.R13
              + c_c * tx * wm.R21 + c_c * c_g * ty * wm.R21 - tx * wm.R22 - c_g * ty * wm.R22 + a_c * tx * wm.R23 + a_c * c_g * ty * wm.R23
              + a_g * c_c * tx * wm.R31 - b_g * c_c * ty * wm.R31 - a_g * tx * wm.R32 + b_g * ty * wm.R32 + a_c * a_g * tx * wm.R33
              - a_c * b_g * ty * wm.R33)
             * ((c_c * c_g * tx * wm.R11 + c_c * ty * wm.R11 - c_g * tx * wm.R12 - ty * wm.R12 + a_c * c_g * tx * wm.R13 + a_c * ty * wm.R13
                 - c_c * tx * wm.R21 - c_c * c_g * ty * wm.R21 + tx * wm.R22 + c_g * ty * wm.R22 - a_c * tx * wm.R23
                 - a_c * c_g * ty * wm.R23 - a_g * c_c * tx * wm.R31 + b_g * c_c * ty * wm.R31 + a_g * tx * wm.R32 - b_g * ty * wm.R32
                 - a_c * a_g * tx * wm.R33 + a_c * b_g * ty * wm.R33)
                    * (-z_c - z_g - z_l
                       - sz
                           * (b_g * wm.R13 - a_g * wm.R23 + b_c * (b_g * wm.R11 - a_g * wm.R21 + wm.R31)
                              + a_c * (b_g * wm.R12 - a_g * wm.R22 + wm.R32) + wm.R33)
                       - sv
                           * (b_g * wm.R12 - a_g * wm.R22 - c_c * (b_g * wm.R11 - a_g * wm.R21 + wm.R31) + wm.R32
                              - a_c * (b_g * wm.R13 - a_g * wm.R23 + wm.R33))
                       - su
                           * (b_g * wm.R11 - a_g * wm.R21 + wm.R31 - c_c * (b_g * wm.R12 - a_g * wm.R22 + wm.R32)
                              + b_c * (b_g * wm.R13 - a_g * wm.R23 + wm.R33)))
                + (-(c_c * wm.R11) + b_g * c_c * tx * wm.R11 + wm.R12 - b_g * tx * wm.R12 - a_c * wm.R13 + a_c * b_g * tx * wm.R13
                   + c_c * c_g * wm.R21 - a_g * c_c * tx * wm.R21 - c_g * wm.R22 + a_g * tx * wm.R22 + a_c * c_g * wm.R23
                   - a_c * a_g * tx * wm.R23 - b_g * c_c * wm.R31 + c_c * tx * wm.R31 + b_g * wm.R32 - tx * wm.R32 - a_c * b_g * wm.R33
                   + a_c * tx * wm.R33)
                    * (by - y_c - y_g - y_l
                       - sz
                           * (-(c_g * wm.R13) + wm.R23 + b_c * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31)
                              + a_c * (-(c_g * wm.R12) + wm.R22 + a_g * wm.R32) + a_g * wm.R33)
                       - sv
                           * (-(c_g * wm.R12) + wm.R22 - c_c * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31) + a_g * wm.R32
                              - a_c * (-(c_g * wm.R13) + wm.R23 + a_g * wm.R33))
                       - su
                           * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31 - c_c * (-(c_g * wm.R12) + wm.R22 + a_g * wm.R32)
                              + b_c * (-(c_g * wm.R13) + wm.R23 + a_g * wm.R33)))
                + (-(c_c * c_g * wm.R11) - b_g * c_c * ty * wm.R11 + c_g * wm.R12 + b_g * ty * wm.R12 - a_c * c_g * wm.R13
                   - a_c * b_g * ty * wm.R13 + c_c * wm.R21 + a_g * c_c * ty * wm.R21 - wm.R22 - a_g * ty * wm.R22 + a_c * wm.R23
                   + a_c * a_g * ty * wm.R23 + a_g * c_c * wm.R31 - c_c * ty * wm.R31 - a_g * wm.R32 + ty * wm.R32 + a_c * a_g * wm.R33
                   - a_c * ty * wm.R33)
                    * (bx - x_c - x_g - x_l
                       - sz
                           * (wm.R13 - c_g * wm.R23 + b_c * (wm.R11 - c_g * wm.R21 + b_g * wm.R31)
                              + a_c * (wm.R12 - c_g * wm.R22 + b_g * wm.R32) + b_g * wm.R33)
                       - sv
                           * (wm.R12 - c_g * wm.R22 - c_c * (wm.R11 - c_g * wm.R21 + b_g * wm.R31) + b_g * wm.R32
                              - a_c * (wm.R13 - c_g * wm.R23 + b_g * wm.R33))
                       - su
                           * (wm.R11 - c_g * wm.R21 + b_g * wm.R31 - c_c * (wm.R12 - c_g * wm.R22 + b_g * wm.R32)
                              + b_c * (wm.R13 - c_g * wm.R23 + b_g * wm.R33)))))
            / (sqrt(pow(-(c_c * wm.R11) + b_g * c_c * tx * wm.R11 + wm.R12 - b_g * tx * wm.R12 - a_c * wm.R13 + a_c * b_g * tx * wm.R13
                            + c_c * c_g * wm.R21 - a_g * c_c * tx * wm.R21 - c_g * wm.R22 + a_g * tx * wm.R22 + a_c * c_g * wm.R23
                            - a_c * a_g * tx * wm.R23 - b_g * c_c * wm.R31 + c_c * tx * wm.R31 + b_g * wm.R32 - tx * wm.R32
                            - a_c * b_g * wm.R33 + a_c * tx * wm.R33,
                        2)
                    + pow(-(c_c * c_g * wm.R11) - b_g * c_c * ty * wm.R11 + c_g * wm.R12 + b_g * ty * wm.R12 - a_c * c_g * wm.R13
                              - a_c * b_g * ty * wm.R13 + c_c * wm.R21 + a_g * c_c * ty * wm.R21 - wm.R22 - a_g * ty * wm.R22 + a_c * wm.R23
                              + a_c * a_g * ty * wm.R23 + a_g * c_c * wm.R31 - c_c * ty * wm.R31 - a_g * wm.R32 + ty * wm.R32
                              + a_c * a_g * wm.R33 - a_c * ty * wm.R33,
                          2)
                    + pow(c_c * c_g * tx * wm.R11 + c_c * ty * wm.R11 - c_g * tx * wm.R12 - ty * wm.R12 + a_c * c_g * tx * wm.R13
                              + a_c * ty * wm.R13 - c_c * tx * wm.R21 - c_c * c_g * ty * wm.R21 + tx * wm.R22 + c_g * ty * wm.R22
                              - a_c * tx * wm.R23 - a_c * c_g * ty * wm.R23 - a_g * c_c * tx * wm.R31 + b_g * c_c * ty * wm.R31
                              + a_g * tx * wm.R32 - b_g * ty * wm.R32 - a_c * a_g * tx * wm.R33 + a_c * b_g * ty * wm.R33,
                          2))
               * fabs(((c_c * c_g * tx * wm.R11 + c_c * ty * wm.R11 - c_g * tx * wm.R12 - ty * wm.R12 + a_c * c_g * tx * wm.R13
                        + a_c * ty * wm.R13 - c_c * tx * wm.R21 - c_c * c_g * ty * wm.R21 + tx * wm.R22 + c_g * ty * wm.R22
                        - a_c * tx * wm.R23 - a_c * c_g * ty * wm.R23 - a_g * c_c * tx * wm.R31 + b_g * c_c * ty * wm.R31
                        + a_g * tx * wm.R32 - b_g * ty * wm.R32 - a_c * a_g * tx * wm.R33 + a_c * b_g * ty * wm.R33)
                           * (-z_c - z_g - z_l
                              - sz
                                  * (b_g * wm.R13 - a_g * wm.R23 + b_c * (b_g * wm.R11 - a_g * wm.R21 + wm.R31)
                                     + a_c * (b_g * wm.R12 - a_g * wm.R22 + wm.R32) + wm.R33)
                              - sv
                                  * (b_g * wm.R12 - a_g * wm.R22 - c_c * (b_g * wm.R11 - a_g * wm.R21 + wm.R31) + wm.R32
                                     - a_c * (b_g * wm.R13 - a_g * wm.R23 + wm.R33))
                              - su
                                  * (b_g * wm.R11 - a_g * wm.R21 + wm.R31 - c_c * (b_g * wm.R12 - a_g * wm.R22 + wm.R32)
                                     + b_c * (b_g * wm.R13 - a_g * wm.R23 + wm.R33)))
                       + (-(c_c * wm.R11) + b_g * c_c * tx * wm.R11 + wm.R12 - b_g * tx * wm.R12 - a_c * wm.R13 + a_c * b_g * tx * wm.R13
                          + c_c * c_g * wm.R21 - a_g * c_c * tx * wm.R21 - c_g * wm.R22 + a_g * tx * wm.R22 + a_c * c_g * wm.R23
                          - a_c * a_g * tx * wm.R23 - b_g * c_c * wm.R31 + c_c * tx * wm.R31 + b_g * wm.R32 - tx * wm.R32
                          - a_c * b_g * wm.R33 + a_c * tx * wm.R33)
                           * (by - y_c - y_g - y_l
                              - sz
                                  * (-(c_g * wm.R13) + wm.R23 + b_c * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31)
                                     + a_c * (-(c_g * wm.R12) + wm.R22 + a_g * wm.R32) + a_g * wm.R33)
                              - sv
                                  * (-(c_g * wm.R12) + wm.R22 - c_c * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31) + a_g * wm.R32
                                     - a_c * (-(c_g * wm.R13) + wm.R23 + a_g * wm.R33))
                              - su
                                  * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31 - c_c * (-(c_g * wm.R12) + wm.R22 + a_g * wm.R32)
                                     + b_c * (-(c_g * wm.R13) + wm.R23 + a_g * wm.R33)))
                       + (-(c_c * c_g * wm.R11) - b_g * c_c * ty * wm.R11 + c_g * wm.R12 + b_g * ty * wm.R12 - a_c * c_g * wm.R13
                          - a_c * b_g * ty * wm.R13 + c_c * wm.R21 + a_g * c_c * ty * wm.R21 - wm.R22 - a_g * ty * wm.R22 + a_c * wm.R23
                          + a_c * a_g * ty * wm.R23 + a_g * c_c * wm.R31 - c_c * ty * wm.R31 - a_g * wm.R32 + ty * wm.R32
                          + a_c * a_g * wm.R33 - a_c * ty * wm.R33)
                           * (bx - x_c - x_g - x_l
                              - sz
                                  * (wm.R13 - c_g * wm.R23 + b_c * (wm.R11 - c_g * wm.R21 + b_g * wm.R31)
                                     + a_c * (wm.R12 - c_g * wm.R22 + b_g * wm.R32) + b_g * wm.R33)
                              - sv
                                  * (wm.R12 - c_g * wm.R22 - c_c * (wm.R11 - c_g * wm.R21 + b_g * wm.R31) + b_g * wm.R32
                                     - a_c * (wm.R13 - c_g * wm.R23 + b_g * wm.R33))
                              - su
                                  * (wm.R11 - c_g * wm.R21 + b_g * wm.R31 - c_c * (wm.R12 - c_g * wm.R22 + b_g * wm.R32)
                                     + b_c * (wm.R13 - c_g * wm.R23 + b_g * wm.R33))))));
        // drdac
        this->template set_global_derivative<9>() =
            (((-(sz * (b_g * wm.R12 - a_g * wm.R22 + wm.R32)) - sv * (-(b_g * wm.R13) + a_g * wm.R23 - wm.R33))
                  * (c_c * c_g * tx * wm.R11 + c_c * ty * wm.R11 - c_g * tx * wm.R12 - ty * wm.R12 + a_c * c_g * tx * wm.R13
                     + a_c * ty * wm.R13 - c_c * tx * wm.R21 - c_c * c_g * ty * wm.R21 + tx * wm.R22 + c_g * ty * wm.R22 - a_c * tx * wm.R23
                     - a_c * c_g * ty * wm.R23 - a_g * c_c * tx * wm.R31 + b_g * c_c * ty * wm.R31 + a_g * tx * wm.R32 - b_g * ty * wm.R32
                     - a_c * a_g * tx * wm.R33 + a_c * b_g * ty * wm.R33)
              + (-(c_c * wm.R11) + b_g * c_c * tx * wm.R11 + wm.R12 - b_g * tx * wm.R12 - a_c * wm.R13 + a_c * b_g * tx * wm.R13
                 + c_c * c_g * wm.R21 - a_g * c_c * tx * wm.R21 - c_g * wm.R22 + a_g * tx * wm.R22 + a_c * c_g * wm.R23
                 - a_c * a_g * tx * wm.R23 - b_g * c_c * wm.R31 + c_c * tx * wm.R31 + b_g * wm.R32 - tx * wm.R32 - a_c * b_g * wm.R33
                 + a_c * tx * wm.R33)
                  * (-(sz * (-(c_g * wm.R12) + wm.R22 + a_g * wm.R32)) - sv * (c_g * wm.R13 - wm.R23 - a_g * wm.R33))
              + (-(c_c * c_g * wm.R11) - b_g * c_c * ty * wm.R11 + c_g * wm.R12 + b_g * ty * wm.R12 - a_c * c_g * wm.R13
                 - a_c * b_g * ty * wm.R13 + c_c * wm.R21 + a_g * c_c * ty * wm.R21 - wm.R22 - a_g * ty * wm.R22 + a_c * wm.R23
                 + a_c * a_g * ty * wm.R23 + a_g * c_c * wm.R31 - c_c * ty * wm.R31 - a_g * wm.R32 + ty * wm.R32 + a_c * a_g * wm.R33
                 - a_c * ty * wm.R33)
                  * (-(sz * (wm.R12 - c_g * wm.R22 + b_g * wm.R32)) - sv * (-wm.R13 + c_g * wm.R23 - b_g * wm.R33))
              + (c_g * tx * wm.R13 + ty * wm.R13 - tx * wm.R23 - c_g * ty * wm.R23 - a_g * tx * wm.R33 + b_g * ty * wm.R33)
                  * (-z_c - z_g - z_l
                     - sz
                         * (b_g * wm.R13 - a_g * wm.R23 + b_c * (b_g * wm.R11 - a_g * wm.R21 + wm.R31)
                            + a_c * (b_g * wm.R12 - a_g * wm.R22 + wm.R32) + wm.R33)
                     - sv
                         * (b_g * wm.R12 - a_g * wm.R22 - c_c * (b_g * wm.R11 - a_g * wm.R21 + wm.R31) + wm.R32
                            - a_c * (b_g * wm.R13 - a_g * wm.R23 + wm.R33))
                     - su
                         * (b_g * wm.R11 - a_g * wm.R21 + wm.R31 - c_c * (b_g * wm.R12 - a_g * wm.R22 + wm.R32)
                            + b_c * (b_g * wm.R13 - a_g * wm.R23 + wm.R33)))
              + (-wm.R13 + b_g * tx * wm.R13 + c_g * wm.R23 - a_g * tx * wm.R23 - b_g * wm.R33 + tx * wm.R33)
                  * (by - y_c - y_g - y_l
                     - sz
                         * (-(c_g * wm.R13) + wm.R23 + b_c * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31)
                            + a_c * (-(c_g * wm.R12) + wm.R22 + a_g * wm.R32) + a_g * wm.R33)
                     - sv
                         * (-(c_g * wm.R12) + wm.R22 - c_c * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31) + a_g * wm.R32
                            - a_c * (-(c_g * wm.R13) + wm.R23 + a_g * wm.R33))
                     - su
                         * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31 - c_c * (-(c_g * wm.R12) + wm.R22 + a_g * wm.R32)
                            + b_c * (-(c_g * wm.R13) + wm.R23 + a_g * wm.R33)))
              + (-(c_g * wm.R13) - b_g * ty * wm.R13 + wm.R23 + a_g * ty * wm.R23 + a_g * wm.R33 - ty * wm.R33)
                  * (bx - x_c - x_g - x_l
                     - sz
                         * (wm.R13 - c_g * wm.R23 + b_c * (wm.R11 - c_g * wm.R21 + b_g * wm.R31)
                            + a_c * (wm.R12 - c_g * wm.R22 + b_g * wm.R32) + b_g * wm.R33)
                     - sv
                         * (wm.R12 - c_g * wm.R22 - c_c * (wm.R11 - c_g * wm.R21 + b_g * wm.R31) + b_g * wm.R32
                            - a_c * (wm.R13 - c_g * wm.R23 + b_g * wm.R33))
                     - su
                         * (wm.R11 - c_g * wm.R21 + b_g * wm.R31 - c_c * (wm.R12 - c_g * wm.R22 + b_g * wm.R32)
                            + b_c * (wm.R13 - c_g * wm.R23 + b_g * wm.R33))))
             * ((c_c * c_g * tx * wm.R11 + c_c * ty * wm.R11 - c_g * tx * wm.R12 - ty * wm.R12 + a_c * c_g * tx * wm.R13 + a_c * ty * wm.R13
                 - c_c * tx * wm.R21 - c_c * c_g * ty * wm.R21 + tx * wm.R22 + c_g * ty * wm.R22 - a_c * tx * wm.R23
                 - a_c * c_g * ty * wm.R23 - a_g * c_c * tx * wm.R31 + b_g * c_c * ty * wm.R31 + a_g * tx * wm.R32 - b_g * ty * wm.R32
                 - a_c * a_g * tx * wm.R33 + a_c * b_g * ty * wm.R33)
                    * (-z_c - z_g - z_l
                       - sz
                           * (b_g * wm.R13 - a_g * wm.R23 + b_c * (b_g * wm.R11 - a_g * wm.R21 + wm.R31)
                              + a_c * (b_g * wm.R12 - a_g * wm.R22 + wm.R32) + wm.R33)
                       - sv
                           * (b_g * wm.R12 - a_g * wm.R22 - c_c * (b_g * wm.R11 - a_g * wm.R21 + wm.R31) + wm.R32
                              - a_c * (b_g * wm.R13 - a_g * wm.R23 + wm.R33))
                       - su
                           * (b_g * wm.R11 - a_g * wm.R21 + wm.R31 - c_c * (b_g * wm.R12 - a_g * wm.R22 + wm.R32)
                              + b_c * (b_g * wm.R13 - a_g * wm.R23 + wm.R33)))
                + (-(c_c * wm.R11) + b_g * c_c * tx * wm.R11 + wm.R12 - b_g * tx * wm.R12 - a_c * wm.R13 + a_c * b_g * tx * wm.R13
                   + c_c * c_g * wm.R21 - a_g * c_c * tx * wm.R21 - c_g * wm.R22 + a_g * tx * wm.R22 + a_c * c_g * wm.R23
                   - a_c * a_g * tx * wm.R23 - b_g * c_c * wm.R31 + c_c * tx * wm.R31 + b_g * wm.R32 - tx * wm.R32 - a_c * b_g * wm.R33
                   + a_c * tx * wm.R33)
                    * (by - y_c - y_g - y_l
                       - sz
                           * (-(c_g * wm.R13) + wm.R23 + b_c * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31)
                              + a_c * (-(c_g * wm.R12) + wm.R22 + a_g * wm.R32) + a_g * wm.R33)
                       - sv
                           * (-(c_g * wm.R12) + wm.R22 - c_c * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31) + a_g * wm.R32
                              - a_c * (-(c_g * wm.R13) + wm.R23 + a_g * wm.R33))
                       - su
                           * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31 - c_c * (-(c_g * wm.R12) + wm.R22 + a_g * wm.R32)
                              + b_c * (-(c_g * wm.R13) + wm.R23 + a_g * wm.R33)))
                + (-(c_c * c_g * wm.R11) - b_g * c_c * ty * wm.R11 + c_g * wm.R12 + b_g * ty * wm.R12 - a_c * c_g * wm.R13
                   - a_c * b_g * ty * wm.R13 + c_c * wm.R21 + a_g * c_c * ty * wm.R21 - wm.R22 - a_g * ty * wm.R22 + a_c * wm.R23
                   + a_c * a_g * ty * wm.R23 + a_g * c_c * wm.R31 - c_c * ty * wm.R31 - a_g * wm.R32 + ty * wm.R32 + a_c * a_g * wm.R33
                   - a_c * ty * wm.R33)
                    * (bx - x_c - x_g - x_l
                       - sz
                           * (wm.R13 - c_g * wm.R23 + b_c * (wm.R11 - c_g * wm.R21 + b_g * wm.R31)
                              + a_c * (wm.R12 - c_g * wm.R22 + b_g * wm.R32) + b_g * wm.R33)
                       - sv
                           * (wm.R12 - c_g * wm.R22 - c_c * (wm.R11 - c_g * wm.R21 + b_g * wm.R31) + b_g * wm.R32
                              - a_c * (wm.R13 - c_g * wm.R23 + b_g * wm.R33))
                       - su
                           * (wm.R11 - c_g * wm.R21 + b_g * wm.R31 - c_c * (wm.R12 - c_g * wm.R22 + b_g * wm.R32)
                              + b_c * (wm.R13 - c_g * wm.R23 + b_g * wm.R33)))))
                / (sqrt(pow(-(c_c * wm.R11) + b_g * c_c * tx * wm.R11 + wm.R12 - b_g * tx * wm.R12 - a_c * wm.R13 + a_c * b_g * tx * wm.R13
                                + c_c * c_g * wm.R21 - a_g * c_c * tx * wm.R21 - c_g * wm.R22 + a_g * tx * wm.R22 + a_c * c_g * wm.R23
                                - a_c * a_g * tx * wm.R23 - b_g * c_c * wm.R31 + c_c * tx * wm.R31 + b_g * wm.R32 - tx * wm.R32
                                - a_c * b_g * wm.R33 + a_c * tx * wm.R33,
                            2)
                        + pow(-(c_c * c_g * wm.R11) - b_g * c_c * ty * wm.R11 + c_g * wm.R12 + b_g * ty * wm.R12 - a_c * c_g * wm.R13
                                  - a_c * b_g * ty * wm.R13 + c_c * wm.R21 + a_g * c_c * ty * wm.R21 - wm.R22 - a_g * ty * wm.R22
                                  + a_c * wm.R23 + a_c * a_g * ty * wm.R23 + a_g * c_c * wm.R31 - c_c * ty * wm.R31 - a_g * wm.R32
                                  + ty * wm.R32 + a_c * a_g * wm.R33 - a_c * ty * wm.R33,
                              2)
                        + pow(c_c * c_g * tx * wm.R11 + c_c * ty * wm.R11 - c_g * tx * wm.R12 - ty * wm.R12 + a_c * c_g * tx * wm.R13
                                  + a_c * ty * wm.R13 - c_c * tx * wm.R21 - c_c * c_g * ty * wm.R21 + tx * wm.R22 + c_g * ty * wm.R22
                                  - a_c * tx * wm.R23 - a_c * c_g * ty * wm.R23 - a_g * c_c * tx * wm.R31 + b_g * c_c * ty * wm.R31
                                  + a_g * tx * wm.R32 - b_g * ty * wm.R32 - a_c * a_g * tx * wm.R33 + a_c * b_g * ty * wm.R33,
                              2))
                   * fabs(((c_c * c_g * tx * wm.R11 + c_c * ty * wm.R11 - c_g * tx * wm.R12 - ty * wm.R12 + a_c * c_g * tx * wm.R13
                            + a_c * ty * wm.R13 - c_c * tx * wm.R21 - c_c * c_g * ty * wm.R21 + tx * wm.R22 + c_g * ty * wm.R22
                            - a_c * tx * wm.R23 - a_c * c_g * ty * wm.R23 - a_g * c_c * tx * wm.R31 + b_g * c_c * ty * wm.R31
                            + a_g * tx * wm.R32 - b_g * ty * wm.R32 - a_c * a_g * tx * wm.R33 + a_c * b_g * ty * wm.R33)
                               * (-z_c - z_g - z_l
                                  - sz
                                      * (b_g * wm.R13 - a_g * wm.R23 + b_c * (b_g * wm.R11 - a_g * wm.R21 + wm.R31)
                                         + a_c * (b_g * wm.R12 - a_g * wm.R22 + wm.R32) + wm.R33)
                                  - sv
                                      * (b_g * wm.R12 - a_g * wm.R22 - c_c * (b_g * wm.R11 - a_g * wm.R21 + wm.R31) + wm.R32
                                         - a_c * (b_g * wm.R13 - a_g * wm.R23 + wm.R33))
                                  - su
                                      * (b_g * wm.R11 - a_g * wm.R21 + wm.R31 - c_c * (b_g * wm.R12 - a_g * wm.R22 + wm.R32)
                                         + b_c * (b_g * wm.R13 - a_g * wm.R23 + wm.R33)))
                           + (-(c_c * wm.R11) + b_g * c_c * tx * wm.R11 + wm.R12 - b_g * tx * wm.R12 - a_c * wm.R13
                              + a_c * b_g * tx * wm.R13 + c_c * c_g * wm.R21 - a_g * c_c * tx * wm.R21 - c_g * wm.R22 + a_g * tx * wm.R22
                              + a_c * c_g * wm.R23 - a_c * a_g * tx * wm.R23 - b_g * c_c * wm.R31 + c_c * tx * wm.R31 + b_g * wm.R32
                              - tx * wm.R32 - a_c * b_g * wm.R33 + a_c * tx * wm.R33)
                               * (by - y_c - y_g - y_l
                                  - sz
                                      * (-(c_g * wm.R13) + wm.R23 + b_c * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31)
                                         + a_c * (-(c_g * wm.R12) + wm.R22 + a_g * wm.R32) + a_g * wm.R33)
                                  - sv
                                      * (-(c_g * wm.R12) + wm.R22 - c_c * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31) + a_g * wm.R32
                                         - a_c * (-(c_g * wm.R13) + wm.R23 + a_g * wm.R33))
                                  - su
                                      * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31 - c_c * (-(c_g * wm.R12) + wm.R22 + a_g * wm.R32)
                                         + b_c * (-(c_g * wm.R13) + wm.R23 + a_g * wm.R33)))
                           + (-(c_c * c_g * wm.R11) - b_g * c_c * ty * wm.R11 + c_g * wm.R12 + b_g * ty * wm.R12 - a_c * c_g * wm.R13
                              - a_c * b_g * ty * wm.R13 + c_c * wm.R21 + a_g * c_c * ty * wm.R21 - wm.R22 - a_g * ty * wm.R22 + a_c * wm.R23
                              + a_c * a_g * ty * wm.R23 + a_g * c_c * wm.R31 - c_c * ty * wm.R31 - a_g * wm.R32 + ty * wm.R32
                              + a_c * a_g * wm.R33 - a_c * ty * wm.R33)
                               * (bx - x_c - x_g - x_l
                                  - sz
                                      * (wm.R13 - c_g * wm.R23 + b_c * (wm.R11 - c_g * wm.R21 + b_g * wm.R31)
                                         + a_c * (wm.R12 - c_g * wm.R22 + b_g * wm.R32) + b_g * wm.R33)
                                  - sv
                                      * (wm.R12 - c_g * wm.R22 - c_c * (wm.R11 - c_g * wm.R21 + b_g * wm.R31) + b_g * wm.R32
                                         - a_c * (wm.R13 - c_g * wm.R23 + b_g * wm.R33))
                                  - su
                                      * (wm.R11 - c_g * wm.R21 + b_g * wm.R31 - c_c * (wm.R12 - c_g * wm.R22 + b_g * wm.R32)
                                         + b_c * (wm.R13 - c_g * wm.R23 + b_g * wm.R33))))))
            - ((2 * (-wm.R13 + b_g * tx * wm.R13 + c_g * wm.R23 - a_g * tx * wm.R23 - b_g * wm.R33 + tx * wm.R33)
                    * (-(c_c * wm.R11) + b_g * c_c * tx * wm.R11 + wm.R12 - b_g * tx * wm.R12 - a_c * wm.R13 + a_c * b_g * tx * wm.R13
                       + c_c * c_g * wm.R21 - a_g * c_c * tx * wm.R21 - c_g * wm.R22 + a_g * tx * wm.R22 + a_c * c_g * wm.R23
                       - a_c * a_g * tx * wm.R23 - b_g * c_c * wm.R31 + c_c * tx * wm.R31 + b_g * wm.R32 - tx * wm.R32 - a_c * b_g * wm.R33
                       + a_c * tx * wm.R33)
                + 2 * (-(c_g * wm.R13) - b_g * ty * wm.R13 + wm.R23 + a_g * ty * wm.R23 + a_g * wm.R33 - ty * wm.R33)
                    * (-(c_c * c_g * wm.R11) - b_g * c_c * ty * wm.R11 + c_g * wm.R12 + b_g * ty * wm.R12 - a_c * c_g * wm.R13
                       - a_c * b_g * ty * wm.R13 + c_c * wm.R21 + a_g * c_c * ty * wm.R21 - wm.R22 - a_g * ty * wm.R22 + a_c * wm.R23
                       + a_c * a_g * ty * wm.R23 + a_g * c_c * wm.R31 - c_c * ty * wm.R31 - a_g * wm.R32 + ty * wm.R32 + a_c * a_g * wm.R33
                       - a_c * ty * wm.R33)
                + 2 * (c_g * tx * wm.R13 + ty * wm.R13 - tx * wm.R23 - c_g * ty * wm.R23 - a_g * tx * wm.R33 + b_g * ty * wm.R33)
                    * (c_c * c_g * tx * wm.R11 + c_c * ty * wm.R11 - c_g * tx * wm.R12 - ty * wm.R12 + a_c * c_g * tx * wm.R13
                       + a_c * ty * wm.R13 - c_c * tx * wm.R21 - c_c * c_g * ty * wm.R21 + tx * wm.R22 + c_g * ty * wm.R22
                       - a_c * tx * wm.R23 - a_c * c_g * ty * wm.R23 - a_g * c_c * tx * wm.R31 + b_g * c_c * ty * wm.R31 + a_g * tx * wm.R32
                       - b_g * ty * wm.R32 - a_c * a_g * tx * wm.R33 + a_c * b_g * ty * wm.R33))
               * fabs(((c_c * c_g * tx * wm.R11 + c_c * ty * wm.R11 - c_g * tx * wm.R12 - ty * wm.R12 + a_c * c_g * tx * wm.R13
                        + a_c * ty * wm.R13 - c_c * tx * wm.R21 - c_c * c_g * ty * wm.R21 + tx * wm.R22 + c_g * ty * wm.R22
                        - a_c * tx * wm.R23 - a_c * c_g * ty * wm.R23 - a_g * c_c * tx * wm.R31 + b_g * c_c * ty * wm.R31
                        + a_g * tx * wm.R32 - b_g * ty * wm.R32 - a_c * a_g * tx * wm.R33 + a_c * b_g * ty * wm.R33)
                           * (-z_c - z_g - z_l
                              - sz
                                  * (b_g * wm.R13 - a_g * wm.R23 + b_c * (b_g * wm.R11 - a_g * wm.R21 + wm.R31)
                                     + a_c * (b_g * wm.R12 - a_g * wm.R22 + wm.R32) + wm.R33)
                              - sv
                                  * (b_g * wm.R12 - a_g * wm.R22 - c_c * (b_g * wm.R11 - a_g * wm.R21 + wm.R31) + wm.R32
                                     - a_c * (b_g * wm.R13 - a_g * wm.R23 + wm.R33))
                              - su
                                  * (b_g * wm.R11 - a_g * wm.R21 + wm.R31 - c_c * (b_g * wm.R12 - a_g * wm.R22 + wm.R32)
                                     + b_c * (b_g * wm.R13 - a_g * wm.R23 + wm.R33)))
                       + (-(c_c * wm.R11) + b_g * c_c * tx * wm.R11 + wm.R12 - b_g * tx * wm.R12 - a_c * wm.R13 + a_c * b_g * tx * wm.R13
                          + c_c * c_g * wm.R21 - a_g * c_c * tx * wm.R21 - c_g * wm.R22 + a_g * tx * wm.R22 + a_c * c_g * wm.R23
                          - a_c * a_g * tx * wm.R23 - b_g * c_c * wm.R31 + c_c * tx * wm.R31 + b_g * wm.R32 - tx * wm.R32
                          - a_c * b_g * wm.R33 + a_c * tx * wm.R33)
                           * (by - y_c - y_g - y_l
                              - sz
                                  * (-(c_g * wm.R13) + wm.R23 + b_c * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31)
                                     + a_c * (-(c_g * wm.R12) + wm.R22 + a_g * wm.R32) + a_g * wm.R33)
                              - sv
                                  * (-(c_g * wm.R12) + wm.R22 - c_c * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31) + a_g * wm.R32
                                     - a_c * (-(c_g * wm.R13) + wm.R23 + a_g * wm.R33))
                              - su
                                  * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31 - c_c * (-(c_g * wm.R12) + wm.R22 + a_g * wm.R32)
                                     + b_c * (-(c_g * wm.R13) + wm.R23 + a_g * wm.R33)))
                       + (-(c_c * c_g * wm.R11) - b_g * c_c * ty * wm.R11 + c_g * wm.R12 + b_g * ty * wm.R12 - a_c * c_g * wm.R13
                          - a_c * b_g * ty * wm.R13 + c_c * wm.R21 + a_g * c_c * ty * wm.R21 - wm.R22 - a_g * ty * wm.R22 + a_c * wm.R23
                          + a_c * a_g * ty * wm.R23 + a_g * c_c * wm.R31 - c_c * ty * wm.R31 - a_g * wm.R32 + ty * wm.R32
                          + a_c * a_g * wm.R33 - a_c * ty * wm.R33)
                           * (bx - x_c - x_g - x_l
                              - sz
                                  * (wm.R13 - c_g * wm.R23 + b_c * (wm.R11 - c_g * wm.R21 + b_g * wm.R31)
                                     + a_c * (wm.R12 - c_g * wm.R22 + b_g * wm.R32) + b_g * wm.R33)
                              - sv
                                  * (wm.R12 - c_g * wm.R22 - c_c * (wm.R11 - c_g * wm.R21 + b_g * wm.R31) + b_g * wm.R32
                                     - a_c * (wm.R13 - c_g * wm.R23 + b_g * wm.R33))
                              - su
                                  * (wm.R11 - c_g * wm.R21 + b_g * wm.R31 - c_c * (wm.R12 - c_g * wm.R22 + b_g * wm.R32)
                                     + b_c * (wm.R13 - c_g * wm.R23 + b_g * wm.R33))))))
                / (2
                   * pow(pow(-(c_c * wm.R11) + b_g * c_c * tx * wm.R11 + wm.R12 - b_g * tx * wm.R12 - a_c * wm.R13 + a_c * b_g * tx * wm.R13
                                 + c_c * c_g * wm.R21 - a_g * c_c * tx * wm.R21 - c_g * wm.R22 + a_g * tx * wm.R22 + a_c * c_g * wm.R23
                                 - a_c * a_g * tx * wm.R23 - b_g * c_c * wm.R31 + c_c * tx * wm.R31 + b_g * wm.R32 - tx * wm.R32
                                 - a_c * b_g * wm.R33 + a_c * tx * wm.R33,
                             2)
                             + pow(-(c_c * c_g * wm.R11) - b_g * c_c * ty * wm.R11 + c_g * wm.R12 + b_g * ty * wm.R12 - a_c * c_g * wm.R13
                                       - a_c * b_g * ty * wm.R13 + c_c * wm.R21 + a_g * c_c * ty * wm.R21 - wm.R22 - a_g * ty * wm.R22
                                       + a_c * wm.R23 + a_c * a_g * ty * wm.R23 + a_g * c_c * wm.R31 - c_c * ty * wm.R31 - a_g * wm.R32
                                       + ty * wm.R32 + a_c * a_g * wm.R33 - a_c * ty * wm.R33,
                                   2)
                             + pow(c_c * c_g * tx * wm.R11 + c_c * ty * wm.R11 - c_g * tx * wm.R12 - ty * wm.R12 + a_c * c_g * tx * wm.R13
                                       + a_c * ty * wm.R13 - c_c * tx * wm.R21 - c_c * c_g * ty * wm.R21 + tx * wm.R22 + c_g * ty * wm.R22
                                       - a_c * tx * wm.R23 - a_c * c_g * ty * wm.R23 - a_g * c_c * tx * wm.R31 + b_g * c_c * ty * wm.R31
                                       + a_g * tx * wm.R32 - b_g * ty * wm.R32 - a_c * a_g * tx * wm.R33 + a_c * b_g * ty * wm.R33,
                                   2),
                         3.0 / 2.0));
        // drdbc
        this->template set_global_derivative<10>() =
            (((c_c * c_g * tx * wm.R11 + c_c * ty * wm.R11 - c_g * tx * wm.R12 - ty * wm.R12 + a_c * c_g * tx * wm.R13 + a_c * ty * wm.R13
               - c_c * tx * wm.R21 - c_c * c_g * ty * wm.R21 + tx * wm.R22 + c_g * ty * wm.R22 - a_c * tx * wm.R23 - a_c * c_g * ty * wm.R23
               - a_g * c_c * tx * wm.R31 + b_g * c_c * ty * wm.R31 + a_g * tx * wm.R32 - b_g * ty * wm.R32 - a_c * a_g * tx * wm.R33
               + a_c * b_g * ty * wm.R33)
                  * (-(sz * (b_g * wm.R11 - a_g * wm.R21 + wm.R31)) - su * (b_g * wm.R13 - a_g * wm.R23 + wm.R33))
              + (-(c_c * wm.R11) + b_g * c_c * tx * wm.R11 + wm.R12 - b_g * tx * wm.R12 - a_c * wm.R13 + a_c * b_g * tx * wm.R13
                 + c_c * c_g * wm.R21 - a_g * c_c * tx * wm.R21 - c_g * wm.R22 + a_g * tx * wm.R22 + a_c * c_g * wm.R23
                 - a_c * a_g * tx * wm.R23 - b_g * c_c * wm.R31 + c_c * tx * wm.R31 + b_g * wm.R32 - tx * wm.R32 - a_c * b_g * wm.R33
                 + a_c * tx * wm.R33)
                  * (-(sz * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31)) - su * (-(c_g * wm.R13) + wm.R23 + a_g * wm.R33))
              + (-(c_c * c_g * wm.R11) - b_g * c_c * ty * wm.R11 + c_g * wm.R12 + b_g * ty * wm.R12 - a_c * c_g * wm.R13
                 - a_c * b_g * ty * wm.R13 + c_c * wm.R21 + a_g * c_c * ty * wm.R21 - wm.R22 - a_g * ty * wm.R22 + a_c * wm.R23
                 + a_c * a_g * ty * wm.R23 + a_g * c_c * wm.R31 - c_c * ty * wm.R31 - a_g * wm.R32 + ty * wm.R32 + a_c * a_g * wm.R33
                 - a_c * ty * wm.R33)
                  * (-(sz * (wm.R11 - c_g * wm.R21 + b_g * wm.R31)) - su * (wm.R13 - c_g * wm.R23 + b_g * wm.R33)))
             * ((c_c * c_g * tx * wm.R11 + c_c * ty * wm.R11 - c_g * tx * wm.R12 - ty * wm.R12 + a_c * c_g * tx * wm.R13 + a_c * ty * wm.R13
                 - c_c * tx * wm.R21 - c_c * c_g * ty * wm.R21 + tx * wm.R22 + c_g * ty * wm.R22 - a_c * tx * wm.R23
                 - a_c * c_g * ty * wm.R23 - a_g * c_c * tx * wm.R31 + b_g * c_c * ty * wm.R31 + a_g * tx * wm.R32 - b_g * ty * wm.R32
                 - a_c * a_g * tx * wm.R33 + a_c * b_g * ty * wm.R33)
                    * (-z_c - z_g - z_l
                       - sz
                           * (b_g * wm.R13 - a_g * wm.R23 + b_c * (b_g * wm.R11 - a_g * wm.R21 + wm.R31)
                              + a_c * (b_g * wm.R12 - a_g * wm.R22 + wm.R32) + wm.R33)
                       - sv
                           * (b_g * wm.R12 - a_g * wm.R22 - c_c * (b_g * wm.R11 - a_g * wm.R21 + wm.R31) + wm.R32
                              - a_c * (b_g * wm.R13 - a_g * wm.R23 + wm.R33))
                       - su
                           * (b_g * wm.R11 - a_g * wm.R21 + wm.R31 - c_c * (b_g * wm.R12 - a_g * wm.R22 + wm.R32)
                              + b_c * (b_g * wm.R13 - a_g * wm.R23 + wm.R33)))
                + (-(c_c * wm.R11) + b_g * c_c * tx * wm.R11 + wm.R12 - b_g * tx * wm.R12 - a_c * wm.R13 + a_c * b_g * tx * wm.R13
                   + c_c * c_g * wm.R21 - a_g * c_c * tx * wm.R21 - c_g * wm.R22 + a_g * tx * wm.R22 + a_c * c_g * wm.R23
                   - a_c * a_g * tx * wm.R23 - b_g * c_c * wm.R31 + c_c * tx * wm.R31 + b_g * wm.R32 - tx * wm.R32 - a_c * b_g * wm.R33
                   + a_c * tx * wm.R33)
                    * (by - y_c - y_g - y_l
                       - sz
                           * (-(c_g * wm.R13) + wm.R23 + b_c * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31)
                              + a_c * (-(c_g * wm.R12) + wm.R22 + a_g * wm.R32) + a_g * wm.R33)
                       - sv
                           * (-(c_g * wm.R12) + wm.R22 - c_c * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31) + a_g * wm.R32
                              - a_c * (-(c_g * wm.R13) + wm.R23 + a_g * wm.R33))
                       - su
                           * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31 - c_c * (-(c_g * wm.R12) + wm.R22 + a_g * wm.R32)
                              + b_c * (-(c_g * wm.R13) + wm.R23 + a_g * wm.R33)))
                + (-(c_c * c_g * wm.R11) - b_g * c_c * ty * wm.R11 + c_g * wm.R12 + b_g * ty * wm.R12 - a_c * c_g * wm.R13
                   - a_c * b_g * ty * wm.R13 + c_c * wm.R21 + a_g * c_c * ty * wm.R21 - wm.R22 - a_g * ty * wm.R22 + a_c * wm.R23
                   + a_c * a_g * ty * wm.R23 + a_g * c_c * wm.R31 - c_c * ty * wm.R31 - a_g * wm.R32 + ty * wm.R32 + a_c * a_g * wm.R33
                   - a_c * ty * wm.R33)
                    * (bx - x_c - x_g - x_l
                       - sz
                           * (wm.R13 - c_g * wm.R23 + b_c * (wm.R11 - c_g * wm.R21 + b_g * wm.R31)
                              + a_c * (wm.R12 - c_g * wm.R22 + b_g * wm.R32) + b_g * wm.R33)
                       - sv
                           * (wm.R12 - c_g * wm.R22 - c_c * (wm.R11 - c_g * wm.R21 + b_g * wm.R31) + b_g * wm.R32
                              - a_c * (wm.R13 - c_g * wm.R23 + b_g * wm.R33))
                       - su
                           * (wm.R11 - c_g * wm.R21 + b_g * wm.R31 - c_c * (wm.R12 - c_g * wm.R22 + b_g * wm.R32)
                              + b_c * (wm.R13 - c_g * wm.R23 + b_g * wm.R33)))))
            / (sqrt(pow(-(c_c * wm.R11) + b_g * c_c * tx * wm.R11 + wm.R12 - b_g * tx * wm.R12 - a_c * wm.R13 + a_c * b_g * tx * wm.R13
                            + c_c * c_g * wm.R21 - a_g * c_c * tx * wm.R21 - c_g * wm.R22 + a_g * tx * wm.R22 + a_c * c_g * wm.R23
                            - a_c * a_g * tx * wm.R23 - b_g * c_c * wm.R31 + c_c * tx * wm.R31 + b_g * wm.R32 - tx * wm.R32
                            - a_c * b_g * wm.R33 + a_c * tx * wm.R33,
                        2)
                    + pow(-(c_c * c_g * wm.R11) - b_g * c_c * ty * wm.R11 + c_g * wm.R12 + b_g * ty * wm.R12 - a_c * c_g * wm.R13
                              - a_c * b_g * ty * wm.R13 + c_c * wm.R21 + a_g * c_c * ty * wm.R21 - wm.R22 - a_g * ty * wm.R22 + a_c * wm.R23
                              + a_c * a_g * ty * wm.R23 + a_g * c_c * wm.R31 - c_c * ty * wm.R31 - a_g * wm.R32 + ty * wm.R32
                              + a_c * a_g * wm.R33 - a_c * ty * wm.R33,
                          2)
                    + pow(c_c * c_g * tx * wm.R11 + c_c * ty * wm.R11 - c_g * tx * wm.R12 - ty * wm.R12 + a_c * c_g * tx * wm.R13
                              + a_c * ty * wm.R13 - c_c * tx * wm.R21 - c_c * c_g * ty * wm.R21 + tx * wm.R22 + c_g * ty * wm.R22
                              - a_c * tx * wm.R23 - a_c * c_g * ty * wm.R23 - a_g * c_c * tx * wm.R31 + b_g * c_c * ty * wm.R31
                              + a_g * tx * wm.R32 - b_g * ty * wm.R32 - a_c * a_g * tx * wm.R33 + a_c * b_g * ty * wm.R33,
                          2))
               * fabs(((c_c * c_g * tx * wm.R11 + c_c * ty * wm.R11 - c_g * tx * wm.R12 - ty * wm.R12 + a_c * c_g * tx * wm.R13
                        + a_c * ty * wm.R13 - c_c * tx * wm.R21 - c_c * c_g * ty * wm.R21 + tx * wm.R22 + c_g * ty * wm.R22
                        - a_c * tx * wm.R23 - a_c * c_g * ty * wm.R23 - a_g * c_c * tx * wm.R31 + b_g * c_c * ty * wm.R31
                        + a_g * tx * wm.R32 - b_g * ty * wm.R32 - a_c * a_g * tx * wm.R33 + a_c * b_g * ty * wm.R33)
                           * (-z_c - z_g - z_l
                              - sz
                                  * (b_g * wm.R13 - a_g * wm.R23 + b_c * (b_g * wm.R11 - a_g * wm.R21 + wm.R31)
                                     + a_c * (b_g * wm.R12 - a_g * wm.R22 + wm.R32) + wm.R33)
                              - sv
                                  * (b_g * wm.R12 - a_g * wm.R22 - c_c * (b_g * wm.R11 - a_g * wm.R21 + wm.R31) + wm.R32
                                     - a_c * (b_g * wm.R13 - a_g * wm.R23 + wm.R33))
                              - su
                                  * (b_g * wm.R11 - a_g * wm.R21 + wm.R31 - c_c * (b_g * wm.R12 - a_g * wm.R22 + wm.R32)
                                     + b_c * (b_g * wm.R13 - a_g * wm.R23 + wm.R33)))
                       + (-(c_c * wm.R11) + b_g * c_c * tx * wm.R11 + wm.R12 - b_g * tx * wm.R12 - a_c * wm.R13 + a_c * b_g * tx * wm.R13
                          + c_c * c_g * wm.R21 - a_g * c_c * tx * wm.R21 - c_g * wm.R22 + a_g * tx * wm.R22 + a_c * c_g * wm.R23
                          - a_c * a_g * tx * wm.R23 - b_g * c_c * wm.R31 + c_c * tx * wm.R31 + b_g * wm.R32 - tx * wm.R32
                          - a_c * b_g * wm.R33 + a_c * tx * wm.R33)
                           * (by - y_c - y_g - y_l
                              - sz
                                  * (-(c_g * wm.R13) + wm.R23 + b_c * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31)
                                     + a_c * (-(c_g * wm.R12) + wm.R22 + a_g * wm.R32) + a_g * wm.R33)
                              - sv
                                  * (-(c_g * wm.R12) + wm.R22 - c_c * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31) + a_g * wm.R32
                                     - a_c * (-(c_g * wm.R13) + wm.R23 + a_g * wm.R33))
                              - su
                                  * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31 - c_c * (-(c_g * wm.R12) + wm.R22 + a_g * wm.R32)
                                     + b_c * (-(c_g * wm.R13) + wm.R23 + a_g * wm.R33)))
                       + (-(c_c * c_g * wm.R11) - b_g * c_c * ty * wm.R11 + c_g * wm.R12 + b_g * ty * wm.R12 - a_c * c_g * wm.R13
                          - a_c * b_g * ty * wm.R13 + c_c * wm.R21 + a_g * c_c * ty * wm.R21 - wm.R22 - a_g * ty * wm.R22 + a_c * wm.R23
                          + a_c * a_g * ty * wm.R23 + a_g * c_c * wm.R31 - c_c * ty * wm.R31 - a_g * wm.R32 + ty * wm.R32
                          + a_c * a_g * wm.R33 - a_c * ty * wm.R33)
                           * (bx - x_c - x_g - x_l
                              - sz
                                  * (wm.R13 - c_g * wm.R23 + b_c * (wm.R11 - c_g * wm.R21 + b_g * wm.R31)
                                     + a_c * (wm.R12 - c_g * wm.R22 + b_g * wm.R32) + b_g * wm.R33)
                              - sv
                                  * (wm.R12 - c_g * wm.R22 - c_c * (wm.R11 - c_g * wm.R21 + b_g * wm.R31) + b_g * wm.R32
                                     - a_c * (wm.R13 - c_g * wm.R23 + b_g * wm.R33))
                              - su
                                  * (wm.R11 - c_g * wm.R21 + b_g * wm.R31 - c_c * (wm.R12 - c_g * wm.R22 + b_g * wm.R32)
                                     + b_c * (wm.R13 - c_g * wm.R23 + b_g * wm.R33))))));
        // drdcc
        this->template set_global_derivative<11>() =
            (((-(sv * (c_g * wm.R11 - wm.R21 - a_g * wm.R31)) - su * (c_g * wm.R12 - wm.R22 - a_g * wm.R32))
                  * (-(c_c * wm.R11) + b_g * c_c * tx * wm.R11 + wm.R12 - b_g * tx * wm.R12 - a_c * wm.R13 + a_c * b_g * tx * wm.R13
                     + c_c * c_g * wm.R21 - a_g * c_c * tx * wm.R21 - c_g * wm.R22 + a_g * tx * wm.R22 + a_c * c_g * wm.R23
                     - a_c * a_g * tx * wm.R23 - b_g * c_c * wm.R31 + c_c * tx * wm.R31 + b_g * wm.R32 - tx * wm.R32 - a_c * b_g * wm.R33
                     + a_c * tx * wm.R33)
              + (-(sv * (-wm.R11 + c_g * wm.R21 - b_g * wm.R31)) - su * (-wm.R12 + c_g * wm.R22 - b_g * wm.R32))
                  * (-(c_c * c_g * wm.R11) - b_g * c_c * ty * wm.R11 + c_g * wm.R12 + b_g * ty * wm.R12 - a_c * c_g * wm.R13
                     - a_c * b_g * ty * wm.R13 + c_c * wm.R21 + a_g * c_c * ty * wm.R21 - wm.R22 - a_g * ty * wm.R22 + a_c * wm.R23
                     + a_c * a_g * ty * wm.R23 + a_g * c_c * wm.R31 - c_c * ty * wm.R31 - a_g * wm.R32 + ty * wm.R32 + a_c * a_g * wm.R33
                     - a_c * ty * wm.R33)
              + (-(sv * (-(b_g * wm.R11) + a_g * wm.R21 - wm.R31)) - su * (-(b_g * wm.R12) + a_g * wm.R22 - wm.R32))
                  * (c_c * c_g * tx * wm.R11 + c_c * ty * wm.R11 - c_g * tx * wm.R12 - ty * wm.R12 + a_c * c_g * tx * wm.R13
                     + a_c * ty * wm.R13 - c_c * tx * wm.R21 - c_c * c_g * ty * wm.R21 + tx * wm.R22 + c_g * ty * wm.R22 - a_c * tx * wm.R23
                     - a_c * c_g * ty * wm.R23 - a_g * c_c * tx * wm.R31 + b_g * c_c * ty * wm.R31 + a_g * tx * wm.R32 - b_g * ty * wm.R32
                     - a_c * a_g * tx * wm.R33 + a_c * b_g * ty * wm.R33)
              + (c_g * tx * wm.R11 + ty * wm.R11 - tx * wm.R21 - c_g * ty * wm.R21 - a_g * tx * wm.R31 + b_g * ty * wm.R31)
                  * (-z_c - z_g - z_l
                     - sz
                         * (b_g * wm.R13 - a_g * wm.R23 + b_c * (b_g * wm.R11 - a_g * wm.R21 + wm.R31)
                            + a_c * (b_g * wm.R12 - a_g * wm.R22 + wm.R32) + wm.R33)
                     - sv
                         * (b_g * wm.R12 - a_g * wm.R22 - c_c * (b_g * wm.R11 - a_g * wm.R21 + wm.R31) + wm.R32
                            - a_c * (b_g * wm.R13 - a_g * wm.R23 + wm.R33))
                     - su
                         * (b_g * wm.R11 - a_g * wm.R21 + wm.R31 - c_c * (b_g * wm.R12 - a_g * wm.R22 + wm.R32)
                            + b_c * (b_g * wm.R13 - a_g * wm.R23 + wm.R33)))
              + (-wm.R11 + b_g * tx * wm.R11 + c_g * wm.R21 - a_g * tx * wm.R21 - b_g * wm.R31 + tx * wm.R31)
                  * (by - y_c - y_g - y_l
                     - sz
                         * (-(c_g * wm.R13) + wm.R23 + b_c * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31)
                            + a_c * (-(c_g * wm.R12) + wm.R22 + a_g * wm.R32) + a_g * wm.R33)
                     - sv
                         * (-(c_g * wm.R12) + wm.R22 - c_c * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31) + a_g * wm.R32
                            - a_c * (-(c_g * wm.R13) + wm.R23 + a_g * wm.R33))
                     - su
                         * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31 - c_c * (-(c_g * wm.R12) + wm.R22 + a_g * wm.R32)
                            + b_c * (-(c_g * wm.R13) + wm.R23 + a_g * wm.R33)))
              + (-(c_g * wm.R11) - b_g * ty * wm.R11 + wm.R21 + a_g * ty * wm.R21 + a_g * wm.R31 - ty * wm.R31)
                  * (bx - x_c - x_g - x_l
                     - sz
                         * (wm.R13 - c_g * wm.R23 + b_c * (wm.R11 - c_g * wm.R21 + b_g * wm.R31)
                            + a_c * (wm.R12 - c_g * wm.R22 + b_g * wm.R32) + b_g * wm.R33)
                     - sv
                         * (wm.R12 - c_g * wm.R22 - c_c * (wm.R11 - c_g * wm.R21 + b_g * wm.R31) + b_g * wm.R32
                            - a_c * (wm.R13 - c_g * wm.R23 + b_g * wm.R33))
                     - su
                         * (wm.R11 - c_g * wm.R21 + b_g * wm.R31 - c_c * (wm.R12 - c_g * wm.R22 + b_g * wm.R32)
                            + b_c * (wm.R13 - c_g * wm.R23 + b_g * wm.R33))))
             * ((c_c * c_g * tx * wm.R11 + c_c * ty * wm.R11 - c_g * tx * wm.R12 - ty * wm.R12 + a_c * c_g * tx * wm.R13 + a_c * ty * wm.R13
                 - c_c * tx * wm.R21 - c_c * c_g * ty * wm.R21 + tx * wm.R22 + c_g * ty * wm.R22 - a_c * tx * wm.R23
                 - a_c * c_g * ty * wm.R23 - a_g * c_c * tx * wm.R31 + b_g * c_c * ty * wm.R31 + a_g * tx * wm.R32 - b_g * ty * wm.R32
                 - a_c * a_g * tx * wm.R33 + a_c * b_g * ty * wm.R33)
                    * (-z_c - z_g - z_l
                       - sz
                           * (b_g * wm.R13 - a_g * wm.R23 + b_c * (b_g * wm.R11 - a_g * wm.R21 + wm.R31)
                              + a_c * (b_g * wm.R12 - a_g * wm.R22 + wm.R32) + wm.R33)
                       - sv
                           * (b_g * wm.R12 - a_g * wm.R22 - c_c * (b_g * wm.R11 - a_g * wm.R21 + wm.R31) + wm.R32
                              - a_c * (b_g * wm.R13 - a_g * wm.R23 + wm.R33))
                       - su
                           * (b_g * wm.R11 - a_g * wm.R21 + wm.R31 - c_c * (b_g * wm.R12 - a_g * wm.R22 + wm.R32)
                              + b_c * (b_g * wm.R13 - a_g * wm.R23 + wm.R33)))
                + (-(c_c * wm.R11) + b_g * c_c * tx * wm.R11 + wm.R12 - b_g * tx * wm.R12 - a_c * wm.R13 + a_c * b_g * tx * wm.R13
                   + c_c * c_g * wm.R21 - a_g * c_c * tx * wm.R21 - c_g * wm.R22 + a_g * tx * wm.R22 + a_c * c_g * wm.R23
                   - a_c * a_g * tx * wm.R23 - b_g * c_c * wm.R31 + c_c * tx * wm.R31 + b_g * wm.R32 - tx * wm.R32 - a_c * b_g * wm.R33
                   + a_c * tx * wm.R33)
                    * (by - y_c - y_g - y_l
                       - sz
                           * (-(c_g * wm.R13) + wm.R23 + b_c * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31)
                              + a_c * (-(c_g * wm.R12) + wm.R22 + a_g * wm.R32) + a_g * wm.R33)
                       - sv
                           * (-(c_g * wm.R12) + wm.R22 - c_c * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31) + a_g * wm.R32
                              - a_c * (-(c_g * wm.R13) + wm.R23 + a_g * wm.R33))
                       - su
                           * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31 - c_c * (-(c_g * wm.R12) + wm.R22 + a_g * wm.R32)
                              + b_c * (-(c_g * wm.R13) + wm.R23 + a_g * wm.R33)))
                + (-(c_c * c_g * wm.R11) - b_g * c_c * ty * wm.R11 + c_g * wm.R12 + b_g * ty * wm.R12 - a_c * c_g * wm.R13
                   - a_c * b_g * ty * wm.R13 + c_c * wm.R21 + a_g * c_c * ty * wm.R21 - wm.R22 - a_g * ty * wm.R22 + a_c * wm.R23
                   + a_c * a_g * ty * wm.R23 + a_g * c_c * wm.R31 - c_c * ty * wm.R31 - a_g * wm.R32 + ty * wm.R32 + a_c * a_g * wm.R33
                   - a_c * ty * wm.R33)
                    * (bx - x_c - x_g - x_l
                       - sz
                           * (wm.R13 - c_g * wm.R23 + b_c * (wm.R11 - c_g * wm.R21 + b_g * wm.R31)
                              + a_c * (wm.R12 - c_g * wm.R22 + b_g * wm.R32) + b_g * wm.R33)
                       - sv
                           * (wm.R12 - c_g * wm.R22 - c_c * (wm.R11 - c_g * wm.R21 + b_g * wm.R31) + b_g * wm.R32
                              - a_c * (wm.R13 - c_g * wm.R23 + b_g * wm.R33))
                       - su
                           * (wm.R11 - c_g * wm.R21 + b_g * wm.R31 - c_c * (wm.R12 - c_g * wm.R22 + b_g * wm.R32)
                              + b_c * (wm.R13 - c_g * wm.R23 + b_g * wm.R33)))))
                / (sqrt(pow(-(c_c * wm.R11) + b_g * c_c * tx * wm.R11 + wm.R12 - b_g * tx * wm.R12 - a_c * wm.R13 + a_c * b_g * tx * wm.R13
                                + c_c * c_g * wm.R21 - a_g * c_c * tx * wm.R21 - c_g * wm.R22 + a_g * tx * wm.R22 + a_c * c_g * wm.R23
                                - a_c * a_g * tx * wm.R23 - b_g * c_c * wm.R31 + c_c * tx * wm.R31 + b_g * wm.R32 - tx * wm.R32
                                - a_c * b_g * wm.R33 + a_c * tx * wm.R33,
                            2)
                        + pow(-(c_c * c_g * wm.R11) - b_g * c_c * ty * wm.R11 + c_g * wm.R12 + b_g * ty * wm.R12 - a_c * c_g * wm.R13
                                  - a_c * b_g * ty * wm.R13 + c_c * wm.R21 + a_g * c_c * ty * wm.R21 - wm.R22 - a_g * ty * wm.R22
                                  + a_c * wm.R23 + a_c * a_g * ty * wm.R23 + a_g * c_c * wm.R31 - c_c * ty * wm.R31 - a_g * wm.R32
                                  + ty * wm.R32 + a_c * a_g * wm.R33 - a_c * ty * wm.R33,
                              2)
                        + pow(c_c * c_g * tx * wm.R11 + c_c * ty * wm.R11 - c_g * tx * wm.R12 - ty * wm.R12 + a_c * c_g * tx * wm.R13
                                  + a_c * ty * wm.R13 - c_c * tx * wm.R21 - c_c * c_g * ty * wm.R21 + tx * wm.R22 + c_g * ty * wm.R22
                                  - a_c * tx * wm.R23 - a_c * c_g * ty * wm.R23 - a_g * c_c * tx * wm.R31 + b_g * c_c * ty * wm.R31
                                  + a_g * tx * wm.R32 - b_g * ty * wm.R32 - a_c * a_g * tx * wm.R33 + a_c * b_g * ty * wm.R33,
                              2))
                   * fabs(((c_c * c_g * tx * wm.R11 + c_c * ty * wm.R11 - c_g * tx * wm.R12 - ty * wm.R12 + a_c * c_g * tx * wm.R13
                            + a_c * ty * wm.R13 - c_c * tx * wm.R21 - c_c * c_g * ty * wm.R21 + tx * wm.R22 + c_g * ty * wm.R22
                            - a_c * tx * wm.R23 - a_c * c_g * ty * wm.R23 - a_g * c_c * tx * wm.R31 + b_g * c_c * ty * wm.R31
                            + a_g * tx * wm.R32 - b_g * ty * wm.R32 - a_c * a_g * tx * wm.R33 + a_c * b_g * ty * wm.R33)
                               * (-z_c - z_g - z_l
                                  - sz
                                      * (b_g * wm.R13 - a_g * wm.R23 + b_c * (b_g * wm.R11 - a_g * wm.R21 + wm.R31)
                                         + a_c * (b_g * wm.R12 - a_g * wm.R22 + wm.R32) + wm.R33)
                                  - sv
                                      * (b_g * wm.R12 - a_g * wm.R22 - c_c * (b_g * wm.R11 - a_g * wm.R21 + wm.R31) + wm.R32
                                         - a_c * (b_g * wm.R13 - a_g * wm.R23 + wm.R33))
                                  - su
                                      * (b_g * wm.R11 - a_g * wm.R21 + wm.R31 - c_c * (b_g * wm.R12 - a_g * wm.R22 + wm.R32)
                                         + b_c * (b_g * wm.R13 - a_g * wm.R23 + wm.R33)))
                           + (-(c_c * wm.R11) + b_g * c_c * tx * wm.R11 + wm.R12 - b_g * tx * wm.R12 - a_c * wm.R13
                              + a_c * b_g * tx * wm.R13 + c_c * c_g * wm.R21 - a_g * c_c * tx * wm.R21 - c_g * wm.R22 + a_g * tx * wm.R22
                              + a_c * c_g * wm.R23 - a_c * a_g * tx * wm.R23 - b_g * c_c * wm.R31 + c_c * tx * wm.R31 + b_g * wm.R32
                              - tx * wm.R32 - a_c * b_g * wm.R33 + a_c * tx * wm.R33)
                               * (by - y_c - y_g - y_l
                                  - sz
                                      * (-(c_g * wm.R13) + wm.R23 + b_c * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31)
                                         + a_c * (-(c_g * wm.R12) + wm.R22 + a_g * wm.R32) + a_g * wm.R33)
                                  - sv
                                      * (-(c_g * wm.R12) + wm.R22 - c_c * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31) + a_g * wm.R32
                                         - a_c * (-(c_g * wm.R13) + wm.R23 + a_g * wm.R33))
                                  - su
                                      * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31 - c_c * (-(c_g * wm.R12) + wm.R22 + a_g * wm.R32)
                                         + b_c * (-(c_g * wm.R13) + wm.R23 + a_g * wm.R33)))
                           + (-(c_c * c_g * wm.R11) - b_g * c_c * ty * wm.R11 + c_g * wm.R12 + b_g * ty * wm.R12 - a_c * c_g * wm.R13
                              - a_c * b_g * ty * wm.R13 + c_c * wm.R21 + a_g * c_c * ty * wm.R21 - wm.R22 - a_g * ty * wm.R22 + a_c * wm.R23
                              + a_c * a_g * ty * wm.R23 + a_g * c_c * wm.R31 - c_c * ty * wm.R31 - a_g * wm.R32 + ty * wm.R32
                              + a_c * a_g * wm.R33 - a_c * ty * wm.R33)
                               * (bx - x_c - x_g - x_l
                                  - sz
                                      * (wm.R13 - c_g * wm.R23 + b_c * (wm.R11 - c_g * wm.R21 + b_g * wm.R31)
                                         + a_c * (wm.R12 - c_g * wm.R22 + b_g * wm.R32) + b_g * wm.R33)
                                  - sv
                                      * (wm.R12 - c_g * wm.R22 - c_c * (wm.R11 - c_g * wm.R21 + b_g * wm.R31) + b_g * wm.R32
                                         - a_c * (wm.R13 - c_g * wm.R23 + b_g * wm.R33))
                                  - su
                                      * (wm.R11 - c_g * wm.R21 + b_g * wm.R31 - c_c * (wm.R12 - c_g * wm.R22 + b_g * wm.R32)
                                         + b_c * (wm.R13 - c_g * wm.R23 + b_g * wm.R33))))))
            - ((2 * (-wm.R11 + b_g * tx * wm.R11 + c_g * wm.R21 - a_g * tx * wm.R21 - b_g * wm.R31 + tx * wm.R31)
                    * (-(c_c * wm.R11) + b_g * c_c * tx * wm.R11 + wm.R12 - b_g * tx * wm.R12 - a_c * wm.R13 + a_c * b_g * tx * wm.R13
                       + c_c * c_g * wm.R21 - a_g * c_c * tx * wm.R21 - c_g * wm.R22 + a_g * tx * wm.R22 + a_c * c_g * wm.R23
                       - a_c * a_g * tx * wm.R23 - b_g * c_c * wm.R31 + c_c * tx * wm.R31 + b_g * wm.R32 - tx * wm.R32 - a_c * b_g * wm.R33
                       + a_c * tx * wm.R33)
                + 2 * (-(c_g * wm.R11) - b_g * ty * wm.R11 + wm.R21 + a_g * ty * wm.R21 + a_g * wm.R31 - ty * wm.R31)
                    * (-(c_c * c_g * wm.R11) - b_g * c_c * ty * wm.R11 + c_g * wm.R12 + b_g * ty * wm.R12 - a_c * c_g * wm.R13
                       - a_c * b_g * ty * wm.R13 + c_c * wm.R21 + a_g * c_c * ty * wm.R21 - wm.R22 - a_g * ty * wm.R22 + a_c * wm.R23
                       + a_c * a_g * ty * wm.R23 + a_g * c_c * wm.R31 - c_c * ty * wm.R31 - a_g * wm.R32 + ty * wm.R32 + a_c * a_g * wm.R33
                       - a_c * ty * wm.R33)
                + 2 * (c_g * tx * wm.R11 + ty * wm.R11 - tx * wm.R21 - c_g * ty * wm.R21 - a_g * tx * wm.R31 + b_g * ty * wm.R31)
                    * (c_c * c_g * tx * wm.R11 + c_c * ty * wm.R11 - c_g * tx * wm.R12 - ty * wm.R12 + a_c * c_g * tx * wm.R13
                       + a_c * ty * wm.R13 - c_c * tx * wm.R21 - c_c * c_g * ty * wm.R21 + tx * wm.R22 + c_g * ty * wm.R22
                       - a_c * tx * wm.R23 - a_c * c_g * ty * wm.R23 - a_g * c_c * tx * wm.R31 + b_g * c_c * ty * wm.R31 + a_g * tx * wm.R32
                       - b_g * ty * wm.R32 - a_c * a_g * tx * wm.R33 + a_c * b_g * ty * wm.R33))
               * fabs(((c_c * c_g * tx * wm.R11 + c_c * ty * wm.R11 - c_g * tx * wm.R12 - ty * wm.R12 + a_c * c_g * tx * wm.R13
                        + a_c * ty * wm.R13 - c_c * tx * wm.R21 - c_c * c_g * ty * wm.R21 + tx * wm.R22 + c_g * ty * wm.R22
                        - a_c * tx * wm.R23 - a_c * c_g * ty * wm.R23 - a_g * c_c * tx * wm.R31 + b_g * c_c * ty * wm.R31
                        + a_g * tx * wm.R32 - b_g * ty * wm.R32 - a_c * a_g * tx * wm.R33 + a_c * b_g * ty * wm.R33)
                           * (-z_c - z_g - z_l
                              - sz
                                  * (b_g * wm.R13 - a_g * wm.R23 + b_c * (b_g * wm.R11 - a_g * wm.R21 + wm.R31)
                                     + a_c * (b_g * wm.R12 - a_g * wm.R22 + wm.R32) + wm.R33)
                              - sv
                                  * (b_g * wm.R12 - a_g * wm.R22 - c_c * (b_g * wm.R11 - a_g * wm.R21 + wm.R31) + wm.R32
                                     - a_c * (b_g * wm.R13 - a_g * wm.R23 + wm.R33))
                              - su
                                  * (b_g * wm.R11 - a_g * wm.R21 + wm.R31 - c_c * (b_g * wm.R12 - a_g * wm.R22 + wm.R32)
                                     + b_c * (b_g * wm.R13 - a_g * wm.R23 + wm.R33)))
                       + (-(c_c * wm.R11) + b_g * c_c * tx * wm.R11 + wm.R12 - b_g * tx * wm.R12 - a_c * wm.R13 + a_c * b_g * tx * wm.R13
                          + c_c * c_g * wm.R21 - a_g * c_c * tx * wm.R21 - c_g * wm.R22 + a_g * tx * wm.R22 + a_c * c_g * wm.R23
                          - a_c * a_g * tx * wm.R23 - b_g * c_c * wm.R31 + c_c * tx * wm.R31 + b_g * wm.R32 - tx * wm.R32
                          - a_c * b_g * wm.R33 + a_c * tx * wm.R33)
                           * (by - y_c - y_g - y_l
                              - sz
                                  * (-(c_g * wm.R13) + wm.R23 + b_c * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31)
                                     + a_c * (-(c_g * wm.R12) + wm.R22 + a_g * wm.R32) + a_g * wm.R33)
                              - sv
                                  * (-(c_g * wm.R12) + wm.R22 - c_c * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31) + a_g * wm.R32
                                     - a_c * (-(c_g * wm.R13) + wm.R23 + a_g * wm.R33))
                              - su
                                  * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31 - c_c * (-(c_g * wm.R12) + wm.R22 + a_g * wm.R32)
                                     + b_c * (-(c_g * wm.R13) + wm.R23 + a_g * wm.R33)))
                       + (-(c_c * c_g * wm.R11) - b_g * c_c * ty * wm.R11 + c_g * wm.R12 + b_g * ty * wm.R12 - a_c * c_g * wm.R13
                          - a_c * b_g * ty * wm.R13 + c_c * wm.R21 + a_g * c_c * ty * wm.R21 - wm.R22 - a_g * ty * wm.R22 + a_c * wm.R23
                          + a_c * a_g * ty * wm.R23 + a_g * c_c * wm.R31 - c_c * ty * wm.R31 - a_g * wm.R32 + ty * wm.R32
                          + a_c * a_g * wm.R33 - a_c * ty * wm.R33)
                           * (bx - x_c - x_g - x_l
                              - sz
                                  * (wm.R13 - c_g * wm.R23 + b_c * (wm.R11 - c_g * wm.R21 + b_g * wm.R31)
                                     + a_c * (wm.R12 - c_g * wm.R22 + b_g * wm.R32) + b_g * wm.R33)
                              - sv
                                  * (wm.R12 - c_g * wm.R22 - c_c * (wm.R11 - c_g * wm.R21 + b_g * wm.R31) + b_g * wm.R32
                                     - a_c * (wm.R13 - c_g * wm.R23 + b_g * wm.R33))
                              - su
                                  * (wm.R11 - c_g * wm.R21 + b_g * wm.R31 - c_c * (wm.R12 - c_g * wm.R22 + b_g * wm.R32)
                                     + b_c * (wm.R13 - c_g * wm.R23 + b_g * wm.R33))))))
                / (2
                   * pow(pow(-(c_c * wm.R11) + b_g * c_c * tx * wm.R11 + wm.R12 - b_g * tx * wm.R12 - a_c * wm.R13 + a_c * b_g * tx * wm.R13
                                 + c_c * c_g * wm.R21 - a_g * c_c * tx * wm.R21 - c_g * wm.R22 + a_g * tx * wm.R22 + a_c * c_g * wm.R23
                                 - a_c * a_g * tx * wm.R23 - b_g * c_c * wm.R31 + c_c * tx * wm.R31 + b_g * wm.R32 - tx * wm.R32
                                 - a_c * b_g * wm.R33 + a_c * tx * wm.R33,
                             2)
                             + pow(-(c_c * c_g * wm.R11) - b_g * c_c * ty * wm.R11 + c_g * wm.R12 + b_g * ty * wm.R12 - a_c * c_g * wm.R13
                                       - a_c * b_g * ty * wm.R13 + c_c * wm.R21 + a_g * c_c * ty * wm.R21 - wm.R22 - a_g * ty * wm.R22
                                       + a_c * wm.R23 + a_c * a_g * ty * wm.R23 + a_g * c_c * wm.R31 - c_c * ty * wm.R31 - a_g * wm.R32
                                       + ty * wm.R32 + a_c * a_g * wm.R33 - a_c * ty * wm.R33,
                                   2)
                             + pow(c_c * c_g * tx * wm.R11 + c_c * ty * wm.R11 - c_g * tx * wm.R12 - ty * wm.R12 + a_c * c_g * tx * wm.R13
                                       + a_c * ty * wm.R13 - c_c * tx * wm.R21 - c_c * c_g * ty * wm.R21 + tx * wm.R22 + c_g * ty * wm.R22
                                       - a_c * tx * wm.R23 - a_c * c_g * ty * wm.R23 - a_g * c_c * tx * wm.R31 + b_g * c_c * ty * wm.R31
                                       + a_g * tx * wm.R32 - b_g * ty * wm.R32 - a_c * a_g * tx * wm.R33 + a_c * b_g * ty * wm.R33,
                                   2),
                         3.0 / 2.0));
        // drdbx
        this->template set_local_derivative<0>() =
            ((-(c_c * c_g * wm.R11) - b_g * c_c * ty * wm.R11 + c_g * wm.R12 + b_g * ty * wm.R12 - a_c * c_g * wm.R13
              - a_c * b_g * ty * wm.R13 + c_c * wm.R21 + a_g * c_c * ty * wm.R21 - wm.R22 - a_g * ty * wm.R22 + a_c * wm.R23
              + a_c * a_g * ty * wm.R23 + a_g * c_c * wm.R31 - c_c * ty * wm.R31 - a_g * wm.R32 + ty * wm.R32 + a_c * a_g * wm.R33
              - a_c * ty * wm.R33)
             * ((c_c * c_g * tx * wm.R11 + c_c * ty * wm.R11 - c_g * tx * wm.R12 - ty * wm.R12 + a_c * c_g * tx * wm.R13 + a_c * ty * wm.R13
                 - c_c * tx * wm.R21 - c_c * c_g * ty * wm.R21 + tx * wm.R22 + c_g * ty * wm.R22 - a_c * tx * wm.R23
                 - a_c * c_g * ty * wm.R23 - a_g * c_c * tx * wm.R31 + b_g * c_c * ty * wm.R31 + a_g * tx * wm.R32 - b_g * ty * wm.R32
                 - a_c * a_g * tx * wm.R33 + a_c * b_g * ty * wm.R33)
                    * (-z_c - z_g - z_l
                       - sz
                           * (b_g * wm.R13 - a_g * wm.R23 + b_c * (b_g * wm.R11 - a_g * wm.R21 + wm.R31)
                              + a_c * (b_g * wm.R12 - a_g * wm.R22 + wm.R32) + wm.R33)
                       - sv
                           * (b_g * wm.R12 - a_g * wm.R22 - c_c * (b_g * wm.R11 - a_g * wm.R21 + wm.R31) + wm.R32
                              - a_c * (b_g * wm.R13 - a_g * wm.R23 + wm.R33))
                       - su
                           * (b_g * wm.R11 - a_g * wm.R21 + wm.R31 - c_c * (b_g * wm.R12 - a_g * wm.R22 + wm.R32)
                              + b_c * (b_g * wm.R13 - a_g * wm.R23 + wm.R33)))
                + (-(c_c * wm.R11) + b_g * c_c * tx * wm.R11 + wm.R12 - b_g * tx * wm.R12 - a_c * wm.R13 + a_c * b_g * tx * wm.R13
                   + c_c * c_g * wm.R21 - a_g * c_c * tx * wm.R21 - c_g * wm.R22 + a_g * tx * wm.R22 + a_c * c_g * wm.R23
                   - a_c * a_g * tx * wm.R23 - b_g * c_c * wm.R31 + c_c * tx * wm.R31 + b_g * wm.R32 - tx * wm.R32 - a_c * b_g * wm.R33
                   + a_c * tx * wm.R33)
                    * (by - y_c - y_g - y_l
                       - sz
                           * (-(c_g * wm.R13) + wm.R23 + b_c * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31)
                              + a_c * (-(c_g * wm.R12) + wm.R22 + a_g * wm.R32) + a_g * wm.R33)
                       - sv
                           * (-(c_g * wm.R12) + wm.R22 - c_c * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31) + a_g * wm.R32
                              - a_c * (-(c_g * wm.R13) + wm.R23 + a_g * wm.R33))
                       - su
                           * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31 - c_c * (-(c_g * wm.R12) + wm.R22 + a_g * wm.R32)
                              + b_c * (-(c_g * wm.R13) + wm.R23 + a_g * wm.R33)))
                + (-(c_c * c_g * wm.R11) - b_g * c_c * ty * wm.R11 + c_g * wm.R12 + b_g * ty * wm.R12 - a_c * c_g * wm.R13
                   - a_c * b_g * ty * wm.R13 + c_c * wm.R21 + a_g * c_c * ty * wm.R21 - wm.R22 - a_g * ty * wm.R22 + a_c * wm.R23
                   + a_c * a_g * ty * wm.R23 + a_g * c_c * wm.R31 - c_c * ty * wm.R31 - a_g * wm.R32 + ty * wm.R32 + a_c * a_g * wm.R33
                   - a_c * ty * wm.R33)
                    * (bx - x_c - x_g - x_l
                       - sz
                           * (wm.R13 - c_g * wm.R23 + b_c * (wm.R11 - c_g * wm.R21 + b_g * wm.R31)
                              + a_c * (wm.R12 - c_g * wm.R22 + b_g * wm.R32) + b_g * wm.R33)
                       - sv
                           * (wm.R12 - c_g * wm.R22 - c_c * (wm.R11 - c_g * wm.R21 + b_g * wm.R31) + b_g * wm.R32
                              - a_c * (wm.R13 - c_g * wm.R23 + b_g * wm.R33))
                       - su
                           * (wm.R11 - c_g * wm.R21 + b_g * wm.R31 - c_c * (wm.R12 - c_g * wm.R22 + b_g * wm.R32)
                              + b_c * (wm.R13 - c_g * wm.R23 + b_g * wm.R33)))))
            / (sqrt(pow(-(c_c * wm.R11) + b_g * c_c * tx * wm.R11 + wm.R12 - b_g * tx * wm.R12 - a_c * wm.R13 + a_c * b_g * tx * wm.R13
                            + c_c * c_g * wm.R21 - a_g * c_c * tx * wm.R21 - c_g * wm.R22 + a_g * tx * wm.R22 + a_c * c_g * wm.R23
                            - a_c * a_g * tx * wm.R23 - b_g * c_c * wm.R31 + c_c * tx * wm.R31 + b_g * wm.R32 - tx * wm.R32
                            - a_c * b_g * wm.R33 + a_c * tx * wm.R33,
                        2)
                    + pow(-(c_c * c_g * wm.R11) - b_g * c_c * ty * wm.R11 + c_g * wm.R12 + b_g * ty * wm.R12 - a_c * c_g * wm.R13
                              - a_c * b_g * ty * wm.R13 + c_c * wm.R21 + a_g * c_c * ty * wm.R21 - wm.R22 - a_g * ty * wm.R22 + a_c * wm.R23
                              + a_c * a_g * ty * wm.R23 + a_g * c_c * wm.R31 - c_c * ty * wm.R31 - a_g * wm.R32 + ty * wm.R32
                              + a_c * a_g * wm.R33 - a_c * ty * wm.R33,
                          2)
                    + pow(c_c * c_g * tx * wm.R11 + c_c * ty * wm.R11 - c_g * tx * wm.R12 - ty * wm.R12 + a_c * c_g * tx * wm.R13
                              + a_c * ty * wm.R13 - c_c * tx * wm.R21 - c_c * c_g * ty * wm.R21 + tx * wm.R22 + c_g * ty * wm.R22
                              - a_c * tx * wm.R23 - a_c * c_g * ty * wm.R23 - a_g * c_c * tx * wm.R31 + b_g * c_c * ty * wm.R31
                              + a_g * tx * wm.R32 - b_g * ty * wm.R32 - a_c * a_g * tx * wm.R33 + a_c * b_g * ty * wm.R33,
                          2))
               * fabs(((c_c * c_g * tx * wm.R11 + c_c * ty * wm.R11 - c_g * tx * wm.R12 - ty * wm.R12 + a_c * c_g * tx * wm.R13
                        + a_c * ty * wm.R13 - c_c * tx * wm.R21 - c_c * c_g * ty * wm.R21 + tx * wm.R22 + c_g * ty * wm.R22
                        - a_c * tx * wm.R23 - a_c * c_g * ty * wm.R23 - a_g * c_c * tx * wm.R31 + b_g * c_c * ty * wm.R31
                        + a_g * tx * wm.R32 - b_g * ty * wm.R32 - a_c * a_g * tx * wm.R33 + a_c * b_g * ty * wm.R33)
                           * (-z_c - z_g - z_l
                              - sz
                                  * (b_g * wm.R13 - a_g * wm.R23 + b_c * (b_g * wm.R11 - a_g * wm.R21 + wm.R31)
                                     + a_c * (b_g * wm.R12 - a_g * wm.R22 + wm.R32) + wm.R33)
                              - sv
                                  * (b_g * wm.R12 - a_g * wm.R22 - c_c * (b_g * wm.R11 - a_g * wm.R21 + wm.R31) + wm.R32
                                     - a_c * (b_g * wm.R13 - a_g * wm.R23 + wm.R33))
                              - su
                                  * (b_g * wm.R11 - a_g * wm.R21 + wm.R31 - c_c * (b_g * wm.R12 - a_g * wm.R22 + wm.R32)
                                     + b_c * (b_g * wm.R13 - a_g * wm.R23 + wm.R33)))
                       + (-(c_c * wm.R11) + b_g * c_c * tx * wm.R11 + wm.R12 - b_g * tx * wm.R12 - a_c * wm.R13 + a_c * b_g * tx * wm.R13
                          + c_c * c_g * wm.R21 - a_g * c_c * tx * wm.R21 - c_g * wm.R22 + a_g * tx * wm.R22 + a_c * c_g * wm.R23
                          - a_c * a_g * tx * wm.R23 - b_g * c_c * wm.R31 + c_c * tx * wm.R31 + b_g * wm.R32 - tx * wm.R32
                          - a_c * b_g * wm.R33 + a_c * tx * wm.R33)
                           * (by - y_c - y_g - y_l
                              - sz
                                  * (-(c_g * wm.R13) + wm.R23 + b_c * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31)
                                     + a_c * (-(c_g * wm.R12) + wm.R22 + a_g * wm.R32) + a_g * wm.R33)
                              - sv
                                  * (-(c_g * wm.R12) + wm.R22 - c_c * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31) + a_g * wm.R32
                                     - a_c * (-(c_g * wm.R13) + wm.R23 + a_g * wm.R33))
                              - su
                                  * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31 - c_c * (-(c_g * wm.R12) + wm.R22 + a_g * wm.R32)
                                     + b_c * (-(c_g * wm.R13) + wm.R23 + a_g * wm.R33)))
                       + (-(c_c * c_g * wm.R11) - b_g * c_c * ty * wm.R11 + c_g * wm.R12 + b_g * ty * wm.R12 - a_c * c_g * wm.R13
                          - a_c * b_g * ty * wm.R13 + c_c * wm.R21 + a_g * c_c * ty * wm.R21 - wm.R22 - a_g * ty * wm.R22 + a_c * wm.R23
                          + a_c * a_g * ty * wm.R23 + a_g * c_c * wm.R31 - c_c * ty * wm.R31 - a_g * wm.R32 + ty * wm.R32
                          + a_c * a_g * wm.R33 - a_c * ty * wm.R33)
                           * (bx - x_c - x_g - x_l
                              - sz
                                  * (wm.R13 - c_g * wm.R23 + b_c * (wm.R11 - c_g * wm.R21 + b_g * wm.R31)
                                     + a_c * (wm.R12 - c_g * wm.R22 + b_g * wm.R32) + b_g * wm.R33)
                              - sv
                                  * (wm.R12 - c_g * wm.R22 - c_c * (wm.R11 - c_g * wm.R21 + b_g * wm.R31) + b_g * wm.R32
                                     - a_c * (wm.R13 - c_g * wm.R23 + b_g * wm.R33))
                              - su
                                  * (wm.R11 - c_g * wm.R21 + b_g * wm.R31 - c_c * (wm.R12 - c_g * wm.R22 + b_g * wm.R32)
                                     + b_c * (wm.R13 - c_g * wm.R23 + b_g * wm.R33))))));
        // drdby
        this->template set_local_derivative<1>() =
            ((-(c_c * wm.R11) + b_g * c_c * tx * wm.R11 + wm.R12 - b_g * tx * wm.R12 - a_c * wm.R13 + a_c * b_g * tx * wm.R13
              + c_c * c_g * wm.R21 - a_g * c_c * tx * wm.R21 - c_g * wm.R22 + a_g * tx * wm.R22 + a_c * c_g * wm.R23
              - a_c * a_g * tx * wm.R23 - b_g * c_c * wm.R31 + c_c * tx * wm.R31 + b_g * wm.R32 - tx * wm.R32 - a_c * b_g * wm.R33
              + a_c * tx * wm.R33)
             * ((c_c * c_g * tx * wm.R11 + c_c * ty * wm.R11 - c_g * tx * wm.R12 - ty * wm.R12 + a_c * c_g * tx * wm.R13 + a_c * ty * wm.R13
                 - c_c * tx * wm.R21 - c_c * c_g * ty * wm.R21 + tx * wm.R22 + c_g * ty * wm.R22 - a_c * tx * wm.R23
                 - a_c * c_g * ty * wm.R23 - a_g * c_c * tx * wm.R31 + b_g * c_c * ty * wm.R31 + a_g * tx * wm.R32 - b_g * ty * wm.R32
                 - a_c * a_g * tx * wm.R33 + a_c * b_g * ty * wm.R33)
                    * (-z_c - z_g - z_l
                       - sz
                           * (b_g * wm.R13 - a_g * wm.R23 + b_c * (b_g * wm.R11 - a_g * wm.R21 + wm.R31)
                              + a_c * (b_g * wm.R12 - a_g * wm.R22 + wm.R32) + wm.R33)
                       - sv
                           * (b_g * wm.R12 - a_g * wm.R22 - c_c * (b_g * wm.R11 - a_g * wm.R21 + wm.R31) + wm.R32
                              - a_c * (b_g * wm.R13 - a_g * wm.R23 + wm.R33))
                       - su
                           * (b_g * wm.R11 - a_g * wm.R21 + wm.R31 - c_c * (b_g * wm.R12 - a_g * wm.R22 + wm.R32)
                              + b_c * (b_g * wm.R13 - a_g * wm.R23 + wm.R33)))
                + (-(c_c * wm.R11) + b_g * c_c * tx * wm.R11 + wm.R12 - b_g * tx * wm.R12 - a_c * wm.R13 + a_c * b_g * tx * wm.R13
                   + c_c * c_g * wm.R21 - a_g * c_c * tx * wm.R21 - c_g * wm.R22 + a_g * tx * wm.R22 + a_c * c_g * wm.R23
                   - a_c * a_g * tx * wm.R23 - b_g * c_c * wm.R31 + c_c * tx * wm.R31 + b_g * wm.R32 - tx * wm.R32 - a_c * b_g * wm.R33
                   + a_c * tx * wm.R33)
                    * (by - y_c - y_g - y_l
                       - sz
                           * (-(c_g * wm.R13) + wm.R23 + b_c * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31)
                              + a_c * (-(c_g * wm.R12) + wm.R22 + a_g * wm.R32) + a_g * wm.R33)
                       - sv
                           * (-(c_g * wm.R12) + wm.R22 - c_c * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31) + a_g * wm.R32
                              - a_c * (-(c_g * wm.R13) + wm.R23 + a_g * wm.R33))
                       - su
                           * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31 - c_c * (-(c_g * wm.R12) + wm.R22 + a_g * wm.R32)
                              + b_c * (-(c_g * wm.R13) + wm.R23 + a_g * wm.R33)))
                + (-(c_c * c_g * wm.R11) - b_g * c_c * ty * wm.R11 + c_g * wm.R12 + b_g * ty * wm.R12 - a_c * c_g * wm.R13
                   - a_c * b_g * ty * wm.R13 + c_c * wm.R21 + a_g * c_c * ty * wm.R21 - wm.R22 - a_g * ty * wm.R22 + a_c * wm.R23
                   + a_c * a_g * ty * wm.R23 + a_g * c_c * wm.R31 - c_c * ty * wm.R31 - a_g * wm.R32 + ty * wm.R32 + a_c * a_g * wm.R33
                   - a_c * ty * wm.R33)
                    * (bx - x_c - x_g - x_l
                       - sz
                           * (wm.R13 - c_g * wm.R23 + b_c * (wm.R11 - c_g * wm.R21 + b_g * wm.R31)
                              + a_c * (wm.R12 - c_g * wm.R22 + b_g * wm.R32) + b_g * wm.R33)
                       - sv
                           * (wm.R12 - c_g * wm.R22 - c_c * (wm.R11 - c_g * wm.R21 + b_g * wm.R31) + b_g * wm.R32
                              - a_c * (wm.R13 - c_g * wm.R23 + b_g * wm.R33))
                       - su
                           * (wm.R11 - c_g * wm.R21 + b_g * wm.R31 - c_c * (wm.R12 - c_g * wm.R22 + b_g * wm.R32)
                              + b_c * (wm.R13 - c_g * wm.R23 + b_g * wm.R33)))))
            / (sqrt(pow(-(c_c * wm.R11) + b_g * c_c * tx * wm.R11 + wm.R12 - b_g * tx * wm.R12 - a_c * wm.R13 + a_c * b_g * tx * wm.R13
                            + c_c * c_g * wm.R21 - a_g * c_c * tx * wm.R21 - c_g * wm.R22 + a_g * tx * wm.R22 + a_c * c_g * wm.R23
                            - a_c * a_g * tx * wm.R23 - b_g * c_c * wm.R31 + c_c * tx * wm.R31 + b_g * wm.R32 - tx * wm.R32
                            - a_c * b_g * wm.R33 + a_c * tx * wm.R33,
                        2)
                    + pow(-(c_c * c_g * wm.R11) - b_g * c_c * ty * wm.R11 + c_g * wm.R12 + b_g * ty * wm.R12 - a_c * c_g * wm.R13
                              - a_c * b_g * ty * wm.R13 + c_c * wm.R21 + a_g * c_c * ty * wm.R21 - wm.R22 - a_g * ty * wm.R22 + a_c * wm.R23
                              + a_c * a_g * ty * wm.R23 + a_g * c_c * wm.R31 - c_c * ty * wm.R31 - a_g * wm.R32 + ty * wm.R32
                              + a_c * a_g * wm.R33 - a_c * ty * wm.R33,
                          2)
                    + pow(c_c * c_g * tx * wm.R11 + c_c * ty * wm.R11 - c_g * tx * wm.R12 - ty * wm.R12 + a_c * c_g * tx * wm.R13
                              + a_c * ty * wm.R13 - c_c * tx * wm.R21 - c_c * c_g * ty * wm.R21 + tx * wm.R22 + c_g * ty * wm.R22
                              - a_c * tx * wm.R23 - a_c * c_g * ty * wm.R23 - a_g * c_c * tx * wm.R31 + b_g * c_c * ty * wm.R31
                              + a_g * tx * wm.R32 - b_g * ty * wm.R32 - a_c * a_g * tx * wm.R33 + a_c * b_g * ty * wm.R33,
                          2))
               * fabs(((c_c * c_g * tx * wm.R11 + c_c * ty * wm.R11 - c_g * tx * wm.R12 - ty * wm.R12 + a_c * c_g * tx * wm.R13
                        + a_c * ty * wm.R13 - c_c * tx * wm.R21 - c_c * c_g * ty * wm.R21 + tx * wm.R22 + c_g * ty * wm.R22
                        - a_c * tx * wm.R23 - a_c * c_g * ty * wm.R23 - a_g * c_c * tx * wm.R31 + b_g * c_c * ty * wm.R31
                        + a_g * tx * wm.R32 - b_g * ty * wm.R32 - a_c * a_g * tx * wm.R33 + a_c * b_g * ty * wm.R33)
                           * (-z_c - z_g - z_l
                              - sz
                                  * (b_g * wm.R13 - a_g * wm.R23 + b_c * (b_g * wm.R11 - a_g * wm.R21 + wm.R31)
                                     + a_c * (b_g * wm.R12 - a_g * wm.R22 + wm.R32) + wm.R33)
                              - sv
                                  * (b_g * wm.R12 - a_g * wm.R22 - c_c * (b_g * wm.R11 - a_g * wm.R21 + wm.R31) + wm.R32
                                     - a_c * (b_g * wm.R13 - a_g * wm.R23 + wm.R33))
                              - su
                                  * (b_g * wm.R11 - a_g * wm.R21 + wm.R31 - c_c * (b_g * wm.R12 - a_g * wm.R22 + wm.R32)
                                     + b_c * (b_g * wm.R13 - a_g * wm.R23 + wm.R33)))
                       + (-(c_c * wm.R11) + b_g * c_c * tx * wm.R11 + wm.R12 - b_g * tx * wm.R12 - a_c * wm.R13 + a_c * b_g * tx * wm.R13
                          + c_c * c_g * wm.R21 - a_g * c_c * tx * wm.R21 - c_g * wm.R22 + a_g * tx * wm.R22 + a_c * c_g * wm.R23
                          - a_c * a_g * tx * wm.R23 - b_g * c_c * wm.R31 + c_c * tx * wm.R31 + b_g * wm.R32 - tx * wm.R32
                          - a_c * b_g * wm.R33 + a_c * tx * wm.R33)
                           * (by - y_c - y_g - y_l
                              - sz
                                  * (-(c_g * wm.R13) + wm.R23 + b_c * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31)
                                     + a_c * (-(c_g * wm.R12) + wm.R22 + a_g * wm.R32) + a_g * wm.R33)
                              - sv
                                  * (-(c_g * wm.R12) + wm.R22 - c_c * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31) + a_g * wm.R32
                                     - a_c * (-(c_g * wm.R13) + wm.R23 + a_g * wm.R33))
                              - su
                                  * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31 - c_c * (-(c_g * wm.R12) + wm.R22 + a_g * wm.R32)
                                     + b_c * (-(c_g * wm.R13) + wm.R23 + a_g * wm.R33)))
                       + (-(c_c * c_g * wm.R11) - b_g * c_c * ty * wm.R11 + c_g * wm.R12 + b_g * ty * wm.R12 - a_c * c_g * wm.R13
                          - a_c * b_g * ty * wm.R13 + c_c * wm.R21 + a_g * c_c * ty * wm.R21 - wm.R22 - a_g * ty * wm.R22 + a_c * wm.R23
                          + a_c * a_g * ty * wm.R23 + a_g * c_c * wm.R31 - c_c * ty * wm.R31 - a_g * wm.R32 + ty * wm.R32
                          + a_c * a_g * wm.R33 - a_c * ty * wm.R33)
                           * (bx - x_c - x_g - x_l
                              - sz
                                  * (wm.R13 - c_g * wm.R23 + b_c * (wm.R11 - c_g * wm.R21 + b_g * wm.R31)
                                     + a_c * (wm.R12 - c_g * wm.R22 + b_g * wm.R32) + b_g * wm.R33)
                              - sv
                                  * (wm.R12 - c_g * wm.R22 - c_c * (wm.R11 - c_g * wm.R21 + b_g * wm.R31) + b_g * wm.R32
                                     - a_c * (wm.R13 - c_g * wm.R23 + b_g * wm.R33))
                              - su
                                  * (wm.R11 - c_g * wm.R21 + b_g * wm.R31 - c_c * (wm.R12 - c_g * wm.R22 + b_g * wm.R32)
                                     + b_c * (wm.R13 - c_g * wm.R23 + b_g * wm.R33))))));
        // drdtx
        this->template set_local_derivative<2>() =
            (((c_c * c_g * wm.R11 - c_g * wm.R12 + a_c * c_g * wm.R13 - c_c * wm.R21 + wm.R22 - a_c * wm.R23 - a_g * c_c * wm.R31
               + a_g * wm.R32 - a_c * a_g * wm.R33)
                  * (-z_c - z_g - z_l
                     - sz
                         * (b_g * wm.R13 - a_g * wm.R23 + b_c * (b_g * wm.R11 - a_g * wm.R21 + wm.R31)
                            + a_c * (b_g * wm.R12 - a_g * wm.R22 + wm.R32) + wm.R33)
                     - sv
                         * (b_g * wm.R12 - a_g * wm.R22 - c_c * (b_g * wm.R11 - a_g * wm.R21 + wm.R31) + wm.R32
                            - a_c * (b_g * wm.R13 - a_g * wm.R23 + wm.R33))
                     - su
                         * (b_g * wm.R11 - a_g * wm.R21 + wm.R31 - c_c * (b_g * wm.R12 - a_g * wm.R22 + wm.R32)
                            + b_c * (b_g * wm.R13 - a_g * wm.R23 + wm.R33)))
              + (b_g * c_c * wm.R11 - b_g * wm.R12 + a_c * b_g * wm.R13 - a_g * c_c * wm.R21 + a_g * wm.R22 - a_c * a_g * wm.R23
                 + c_c * wm.R31 - wm.R32 + a_c * wm.R33)
                  * (by - y_c - y_g - y_l
                     - sz
                         * (-(c_g * wm.R13) + wm.R23 + b_c * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31)
                            + a_c * (-(c_g * wm.R12) + wm.R22 + a_g * wm.R32) + a_g * wm.R33)
                     - sv
                         * (-(c_g * wm.R12) + wm.R22 - c_c * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31) + a_g * wm.R32
                            - a_c * (-(c_g * wm.R13) + wm.R23 + a_g * wm.R33))
                     - su
                         * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31 - c_c * (-(c_g * wm.R12) + wm.R22 + a_g * wm.R32)
                            + b_c * (-(c_g * wm.R13) + wm.R23 + a_g * wm.R33))))
             * ((c_c * c_g * tx * wm.R11 + c_c * ty * wm.R11 - c_g * tx * wm.R12 - ty * wm.R12 + a_c * c_g * tx * wm.R13 + a_c * ty * wm.R13
                 - c_c * tx * wm.R21 - c_c * c_g * ty * wm.R21 + tx * wm.R22 + c_g * ty * wm.R22 - a_c * tx * wm.R23
                 - a_c * c_g * ty * wm.R23 - a_g * c_c * tx * wm.R31 + b_g * c_c * ty * wm.R31 + a_g * tx * wm.R32 - b_g * ty * wm.R32
                 - a_c * a_g * tx * wm.R33 + a_c * b_g * ty * wm.R33)
                    * (-z_c - z_g - z_l
                       - sz
                           * (b_g * wm.R13 - a_g * wm.R23 + b_c * (b_g * wm.R11 - a_g * wm.R21 + wm.R31)
                              + a_c * (b_g * wm.R12 - a_g * wm.R22 + wm.R32) + wm.R33)
                       - sv
                           * (b_g * wm.R12 - a_g * wm.R22 - c_c * (b_g * wm.R11 - a_g * wm.R21 + wm.R31) + wm.R32
                              - a_c * (b_g * wm.R13 - a_g * wm.R23 + wm.R33))
                       - su
                           * (b_g * wm.R11 - a_g * wm.R21 + wm.R31 - c_c * (b_g * wm.R12 - a_g * wm.R22 + wm.R32)
                              + b_c * (b_g * wm.R13 - a_g * wm.R23 + wm.R33)))
                + (-(c_c * wm.R11) + b_g * c_c * tx * wm.R11 + wm.R12 - b_g * tx * wm.R12 - a_c * wm.R13 + a_c * b_g * tx * wm.R13
                   + c_c * c_g * wm.R21 - a_g * c_c * tx * wm.R21 - c_g * wm.R22 + a_g * tx * wm.R22 + a_c * c_g * wm.R23
                   - a_c * a_g * tx * wm.R23 - b_g * c_c * wm.R31 + c_c * tx * wm.R31 + b_g * wm.R32 - tx * wm.R32 - a_c * b_g * wm.R33
                   + a_c * tx * wm.R33)
                    * (by - y_c - y_g - y_l
                       - sz
                           * (-(c_g * wm.R13) + wm.R23 + b_c * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31)
                              + a_c * (-(c_g * wm.R12) + wm.R22 + a_g * wm.R32) + a_g * wm.R33)
                       - sv
                           * (-(c_g * wm.R12) + wm.R22 - c_c * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31) + a_g * wm.R32
                              - a_c * (-(c_g * wm.R13) + wm.R23 + a_g * wm.R33))
                       - su
                           * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31 - c_c * (-(c_g * wm.R12) + wm.R22 + a_g * wm.R32)
                              + b_c * (-(c_g * wm.R13) + wm.R23 + a_g * wm.R33)))
                + (-(c_c * c_g * wm.R11) - b_g * c_c * ty * wm.R11 + c_g * wm.R12 + b_g * ty * wm.R12 - a_c * c_g * wm.R13
                   - a_c * b_g * ty * wm.R13 + c_c * wm.R21 + a_g * c_c * ty * wm.R21 - wm.R22 - a_g * ty * wm.R22 + a_c * wm.R23
                   + a_c * a_g * ty * wm.R23 + a_g * c_c * wm.R31 - c_c * ty * wm.R31 - a_g * wm.R32 + ty * wm.R32 + a_c * a_g * wm.R33
                   - a_c * ty * wm.R33)
                    * (bx - x_c - x_g - x_l
                       - sz
                           * (wm.R13 - c_g * wm.R23 + b_c * (wm.R11 - c_g * wm.R21 + b_g * wm.R31)
                              + a_c * (wm.R12 - c_g * wm.R22 + b_g * wm.R32) + b_g * wm.R33)
                       - sv
                           * (wm.R12 - c_g * wm.R22 - c_c * (wm.R11 - c_g * wm.R21 + b_g * wm.R31) + b_g * wm.R32
                              - a_c * (wm.R13 - c_g * wm.R23 + b_g * wm.R33))
                       - su
                           * (wm.R11 - c_g * wm.R21 + b_g * wm.R31 - c_c * (wm.R12 - c_g * wm.R22 + b_g * wm.R32)
                              + b_c * (wm.R13 - c_g * wm.R23 + b_g * wm.R33)))))
                / (sqrt(pow(-(c_c * wm.R11) + b_g * c_c * tx * wm.R11 + wm.R12 - b_g * tx * wm.R12 - a_c * wm.R13 + a_c * b_g * tx * wm.R13
                                + c_c * c_g * wm.R21 - a_g * c_c * tx * wm.R21 - c_g * wm.R22 + a_g * tx * wm.R22 + a_c * c_g * wm.R23
                                - a_c * a_g * tx * wm.R23 - b_g * c_c * wm.R31 + c_c * tx * wm.R31 + b_g * wm.R32 - tx * wm.R32
                                - a_c * b_g * wm.R33 + a_c * tx * wm.R33,
                            2)
                        + pow(-(c_c * c_g * wm.R11) - b_g * c_c * ty * wm.R11 + c_g * wm.R12 + b_g * ty * wm.R12 - a_c * c_g * wm.R13
                                  - a_c * b_g * ty * wm.R13 + c_c * wm.R21 + a_g * c_c * ty * wm.R21 - wm.R22 - a_g * ty * wm.R22
                                  + a_c * wm.R23 + a_c * a_g * ty * wm.R23 + a_g * c_c * wm.R31 - c_c * ty * wm.R31 - a_g * wm.R32
                                  + ty * wm.R32 + a_c * a_g * wm.R33 - a_c * ty * wm.R33,
                              2)
                        + pow(c_c * c_g * tx * wm.R11 + c_c * ty * wm.R11 - c_g * tx * wm.R12 - ty * wm.R12 + a_c * c_g * tx * wm.R13
                                  + a_c * ty * wm.R13 - c_c * tx * wm.R21 - c_c * c_g * ty * wm.R21 + tx * wm.R22 + c_g * ty * wm.R22
                                  - a_c * tx * wm.R23 - a_c * c_g * ty * wm.R23 - a_g * c_c * tx * wm.R31 + b_g * c_c * ty * wm.R31
                                  + a_g * tx * wm.R32 - b_g * ty * wm.R32 - a_c * a_g * tx * wm.R33 + a_c * b_g * ty * wm.R33,
                              2))
                   * fabs(((c_c * c_g * tx * wm.R11 + c_c * ty * wm.R11 - c_g * tx * wm.R12 - ty * wm.R12 + a_c * c_g * tx * wm.R13
                            + a_c * ty * wm.R13 - c_c * tx * wm.R21 - c_c * c_g * ty * wm.R21 + tx * wm.R22 + c_g * ty * wm.R22
                            - a_c * tx * wm.R23 - a_c * c_g * ty * wm.R23 - a_g * c_c * tx * wm.R31 + b_g * c_c * ty * wm.R31
                            + a_g * tx * wm.R32 - b_g * ty * wm.R32 - a_c * a_g * tx * wm.R33 + a_c * b_g * ty * wm.R33)
                               * (-z_c - z_g - z_l
                                  - sz
                                      * (b_g * wm.R13 - a_g * wm.R23 + b_c * (b_g * wm.R11 - a_g * wm.R21 + wm.R31)
                                         + a_c * (b_g * wm.R12 - a_g * wm.R22 + wm.R32) + wm.R33)
                                  - sv
                                      * (b_g * wm.R12 - a_g * wm.R22 - c_c * (b_g * wm.R11 - a_g * wm.R21 + wm.R31) + wm.R32
                                         - a_c * (b_g * wm.R13 - a_g * wm.R23 + wm.R33))
                                  - su
                                      * (b_g * wm.R11 - a_g * wm.R21 + wm.R31 - c_c * (b_g * wm.R12 - a_g * wm.R22 + wm.R32)
                                         + b_c * (b_g * wm.R13 - a_g * wm.R23 + wm.R33)))
                           + (-(c_c * wm.R11) + b_g * c_c * tx * wm.R11 + wm.R12 - b_g * tx * wm.R12 - a_c * wm.R13
                              + a_c * b_g * tx * wm.R13 + c_c * c_g * wm.R21 - a_g * c_c * tx * wm.R21 - c_g * wm.R22 + a_g * tx * wm.R22
                              + a_c * c_g * wm.R23 - a_c * a_g * tx * wm.R23 - b_g * c_c * wm.R31 + c_c * tx * wm.R31 + b_g * wm.R32
                              - tx * wm.R32 - a_c * b_g * wm.R33 + a_c * tx * wm.R33)
                               * (by - y_c - y_g - y_l
                                  - sz
                                      * (-(c_g * wm.R13) + wm.R23 + b_c * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31)
                                         + a_c * (-(c_g * wm.R12) + wm.R22 + a_g * wm.R32) + a_g * wm.R33)
                                  - sv
                                      * (-(c_g * wm.R12) + wm.R22 - c_c * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31) + a_g * wm.R32
                                         - a_c * (-(c_g * wm.R13) + wm.R23 + a_g * wm.R33))
                                  - su
                                      * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31 - c_c * (-(c_g * wm.R12) + wm.R22 + a_g * wm.R32)
                                         + b_c * (-(c_g * wm.R13) + wm.R23 + a_g * wm.R33)))
                           + (-(c_c * c_g * wm.R11) - b_g * c_c * ty * wm.R11 + c_g * wm.R12 + b_g * ty * wm.R12 - a_c * c_g * wm.R13
                              - a_c * b_g * ty * wm.R13 + c_c * wm.R21 + a_g * c_c * ty * wm.R21 - wm.R22 - a_g * ty * wm.R22 + a_c * wm.R23
                              + a_c * a_g * ty * wm.R23 + a_g * c_c * wm.R31 - c_c * ty * wm.R31 - a_g * wm.R32 + ty * wm.R32
                              + a_c * a_g * wm.R33 - a_c * ty * wm.R33)
                               * (bx - x_c - x_g - x_l
                                  - sz
                                      * (wm.R13 - c_g * wm.R23 + b_c * (wm.R11 - c_g * wm.R21 + b_g * wm.R31)
                                         + a_c * (wm.R12 - c_g * wm.R22 + b_g * wm.R32) + b_g * wm.R33)
                                  - sv
                                      * (wm.R12 - c_g * wm.R22 - c_c * (wm.R11 - c_g * wm.R21 + b_g * wm.R31) + b_g * wm.R32
                                         - a_c * (wm.R13 - c_g * wm.R23 + b_g * wm.R33))
                                  - su
                                      * (wm.R11 - c_g * wm.R21 + b_g * wm.R31 - c_c * (wm.R12 - c_g * wm.R22 + b_g * wm.R32)
                                         + b_c * (wm.R13 - c_g * wm.R23 + b_g * wm.R33))))))
            - ((2
                    * (b_g * c_c * wm.R11 - b_g * wm.R12 + a_c * b_g * wm.R13 - a_g * c_c * wm.R21 + a_g * wm.R22 - a_c * a_g * wm.R23
                       + c_c * wm.R31 - wm.R32 + a_c * wm.R33)
                    * (-(c_c * wm.R11) + b_g * c_c * tx * wm.R11 + wm.R12 - b_g * tx * wm.R12 - a_c * wm.R13 + a_c * b_g * tx * wm.R13
                       + c_c * c_g * wm.R21 - a_g * c_c * tx * wm.R21 - c_g * wm.R22 + a_g * tx * wm.R22 + a_c * c_g * wm.R23
                       - a_c * a_g * tx * wm.R23 - b_g * c_c * wm.R31 + c_c * tx * wm.R31 + b_g * wm.R32 - tx * wm.R32 - a_c * b_g * wm.R33
                       + a_c * tx * wm.R33)
                + 2 * (c_c * c_g * wm.R11 - c_g * wm.R12 + a_c * c_g * wm.R13 - c_c * wm.R21 + wm.R22 - a_c * wm.R23 - a_g * c_c * wm.R31 + a_g * wm.R32 - a_c * a_g * wm.R33)
                    * (c_c * c_g * tx * wm.R11 + c_c * ty * wm.R11 - c_g * tx * wm.R12 - ty * wm.R12 + a_c * c_g * tx * wm.R13
                       + a_c * ty * wm.R13 - c_c * tx * wm.R21 - c_c * c_g * ty * wm.R21 + tx * wm.R22 + c_g * ty * wm.R22
                       - a_c * tx * wm.R23 - a_c * c_g * ty * wm.R23 - a_g * c_c * tx * wm.R31 + b_g * c_c * ty * wm.R31 + a_g * tx * wm.R32
                       - b_g * ty * wm.R32 - a_c * a_g * tx * wm.R33 + a_c * b_g * ty * wm.R33))
               * fabs(((c_c * c_g * tx * wm.R11 + c_c * ty * wm.R11 - c_g * tx * wm.R12 - ty * wm.R12 + a_c * c_g * tx * wm.R13
                        + a_c * ty * wm.R13 - c_c * tx * wm.R21 - c_c * c_g * ty * wm.R21 + tx * wm.R22 + c_g * ty * wm.R22
                        - a_c * tx * wm.R23 - a_c * c_g * ty * wm.R23 - a_g * c_c * tx * wm.R31 + b_g * c_c * ty * wm.R31
                        + a_g * tx * wm.R32 - b_g * ty * wm.R32 - a_c * a_g * tx * wm.R33 + a_c * b_g * ty * wm.R33)
                           * (-z_c - z_g - z_l
                              - sz
                                  * (b_g * wm.R13 - a_g * wm.R23 + b_c * (b_g * wm.R11 - a_g * wm.R21 + wm.R31)
                                     + a_c * (b_g * wm.R12 - a_g * wm.R22 + wm.R32) + wm.R33)
                              - sv
                                  * (b_g * wm.R12 - a_g * wm.R22 - c_c * (b_g * wm.R11 - a_g * wm.R21 + wm.R31) + wm.R32
                                     - a_c * (b_g * wm.R13 - a_g * wm.R23 + wm.R33))
                              - su
                                  * (b_g * wm.R11 - a_g * wm.R21 + wm.R31 - c_c * (b_g * wm.R12 - a_g * wm.R22 + wm.R32)
                                     + b_c * (b_g * wm.R13 - a_g * wm.R23 + wm.R33)))
                       + (-(c_c * wm.R11) + b_g * c_c * tx * wm.R11 + wm.R12 - b_g * tx * wm.R12 - a_c * wm.R13 + a_c * b_g * tx * wm.R13
                          + c_c * c_g * wm.R21 - a_g * c_c * tx * wm.R21 - c_g * wm.R22 + a_g * tx * wm.R22 + a_c * c_g * wm.R23
                          - a_c * a_g * tx * wm.R23 - b_g * c_c * wm.R31 + c_c * tx * wm.R31 + b_g * wm.R32 - tx * wm.R32
                          - a_c * b_g * wm.R33 + a_c * tx * wm.R33)
                           * (by - y_c - y_g - y_l
                              - sz
                                  * (-(c_g * wm.R13) + wm.R23 + b_c * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31)
                                     + a_c * (-(c_g * wm.R12) + wm.R22 + a_g * wm.R32) + a_g * wm.R33)
                              - sv
                                  * (-(c_g * wm.R12) + wm.R22 - c_c * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31) + a_g * wm.R32
                                     - a_c * (-(c_g * wm.R13) + wm.R23 + a_g * wm.R33))
                              - su
                                  * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31 - c_c * (-(c_g * wm.R12) + wm.R22 + a_g * wm.R32)
                                     + b_c * (-(c_g * wm.R13) + wm.R23 + a_g * wm.R33)))
                       + (-(c_c * c_g * wm.R11) - b_g * c_c * ty * wm.R11 + c_g * wm.R12 + b_g * ty * wm.R12 - a_c * c_g * wm.R13
                          - a_c * b_g * ty * wm.R13 + c_c * wm.R21 + a_g * c_c * ty * wm.R21 - wm.R22 - a_g * ty * wm.R22 + a_c * wm.R23
                          + a_c * a_g * ty * wm.R23 + a_g * c_c * wm.R31 - c_c * ty * wm.R31 - a_g * wm.R32 + ty * wm.R32
                          + a_c * a_g * wm.R33 - a_c * ty * wm.R33)
                           * (bx - x_c - x_g - x_l
                              - sz
                                  * (wm.R13 - c_g * wm.R23 + b_c * (wm.R11 - c_g * wm.R21 + b_g * wm.R31)
                                     + a_c * (wm.R12 - c_g * wm.R22 + b_g * wm.R32) + b_g * wm.R33)
                              - sv
                                  * (wm.R12 - c_g * wm.R22 - c_c * (wm.R11 - c_g * wm.R21 + b_g * wm.R31) + b_g * wm.R32
                                     - a_c * (wm.R13 - c_g * wm.R23 + b_g * wm.R33))
                              - su
                                  * (wm.R11 - c_g * wm.R21 + b_g * wm.R31 - c_c * (wm.R12 - c_g * wm.R22 + b_g * wm.R32)
                                     + b_c * (wm.R13 - c_g * wm.R23 + b_g * wm.R33))))))
                / (2
                   * pow(pow(-(c_c * wm.R11) + b_g * c_c * tx * wm.R11 + wm.R12 - b_g * tx * wm.R12 - a_c * wm.R13 + a_c * b_g * tx * wm.R13
                                 + c_c * c_g * wm.R21 - a_g * c_c * tx * wm.R21 - c_g * wm.R22 + a_g * tx * wm.R22 + a_c * c_g * wm.R23
                                 - a_c * a_g * tx * wm.R23 - b_g * c_c * wm.R31 + c_c * tx * wm.R31 + b_g * wm.R32 - tx * wm.R32
                                 - a_c * b_g * wm.R33 + a_c * tx * wm.R33,
                             2)
                             + pow(-(c_c * c_g * wm.R11) - b_g * c_c * ty * wm.R11 + c_g * wm.R12 + b_g * ty * wm.R12 - a_c * c_g * wm.R13
                                       - a_c * b_g * ty * wm.R13 + c_c * wm.R21 + a_g * c_c * ty * wm.R21 - wm.R22 - a_g * ty * wm.R22
                                       + a_c * wm.R23 + a_c * a_g * ty * wm.R23 + a_g * c_c * wm.R31 - c_c * ty * wm.R31 - a_g * wm.R32
                                       + ty * wm.R32 + a_c * a_g * wm.R33 - a_c * ty * wm.R33,
                                   2)
                             + pow(c_c * c_g * tx * wm.R11 + c_c * ty * wm.R11 - c_g * tx * wm.R12 - ty * wm.R12 + a_c * c_g * tx * wm.R13
                                       + a_c * ty * wm.R13 - c_c * tx * wm.R21 - c_c * c_g * ty * wm.R21 + tx * wm.R22 + c_g * ty * wm.R22
                                       - a_c * tx * wm.R23 - a_c * c_g * ty * wm.R23 - a_g * c_c * tx * wm.R31 + b_g * c_c * ty * wm.R31
                                       + a_g * tx * wm.R32 - b_g * ty * wm.R32 - a_c * a_g * tx * wm.R33 + a_c * b_g * ty * wm.R33,
                                   2),
                         3.0 / 2.0));
        // drdty
        this->template set_local_derivative<3>() =
            (((c_c * wm.R11 - wm.R12 + a_c * wm.R13 - c_c * c_g * wm.R21 + c_g * wm.R22 - a_c * c_g * wm.R23 + b_g * c_c * wm.R31
               - b_g * wm.R32 + a_c * b_g * wm.R33)
                  * (-z_c - z_g - z_l
                     - sz
                         * (b_g * wm.R13 - a_g * wm.R23 + b_c * (b_g * wm.R11 - a_g * wm.R21 + wm.R31)
                            + a_c * (b_g * wm.R12 - a_g * wm.R22 + wm.R32) + wm.R33)
                     - sv
                         * (b_g * wm.R12 - a_g * wm.R22 - c_c * (b_g * wm.R11 - a_g * wm.R21 + wm.R31) + wm.R32
                            - a_c * (b_g * wm.R13 - a_g * wm.R23 + wm.R33))
                     - su
                         * (b_g * wm.R11 - a_g * wm.R21 + wm.R31 - c_c * (b_g * wm.R12 - a_g * wm.R22 + wm.R32)
                            + b_c * (b_g * wm.R13 - a_g * wm.R23 + wm.R33)))
              + (-(b_g * c_c * wm.R11) + b_g * wm.R12 - a_c * b_g * wm.R13 + a_g * c_c * wm.R21 - a_g * wm.R22 + a_c * a_g * wm.R23
                 - c_c * wm.R31 + wm.R32 - a_c * wm.R33)
                  * (bx - x_c - x_g - x_l
                     - sz
                         * (wm.R13 - c_g * wm.R23 + b_c * (wm.R11 - c_g * wm.R21 + b_g * wm.R31)
                            + a_c * (wm.R12 - c_g * wm.R22 + b_g * wm.R32) + b_g * wm.R33)
                     - sv
                         * (wm.R12 - c_g * wm.R22 - c_c * (wm.R11 - c_g * wm.R21 + b_g * wm.R31) + b_g * wm.R32
                            - a_c * (wm.R13 - c_g * wm.R23 + b_g * wm.R33))
                     - su
                         * (wm.R11 - c_g * wm.R21 + b_g * wm.R31 - c_c * (wm.R12 - c_g * wm.R22 + b_g * wm.R32)
                            + b_c * (wm.R13 - c_g * wm.R23 + b_g * wm.R33))))
             * ((c_c * c_g * tx * wm.R11 + c_c * ty * wm.R11 - c_g * tx * wm.R12 - ty * wm.R12 + a_c * c_g * tx * wm.R13 + a_c * ty * wm.R13
                 - c_c * tx * wm.R21 - c_c * c_g * ty * wm.R21 + tx * wm.R22 + c_g * ty * wm.R22 - a_c * tx * wm.R23
                 - a_c * c_g * ty * wm.R23 - a_g * c_c * tx * wm.R31 + b_g * c_c * ty * wm.R31 + a_g * tx * wm.R32 - b_g * ty * wm.R32
                 - a_c * a_g * tx * wm.R33 + a_c * b_g * ty * wm.R33)
                    * (-z_c - z_g - z_l
                       - sz
                           * (b_g * wm.R13 - a_g * wm.R23 + b_c * (b_g * wm.R11 - a_g * wm.R21 + wm.R31)
                              + a_c * (b_g * wm.R12 - a_g * wm.R22 + wm.R32) + wm.R33)
                       - sv
                           * (b_g * wm.R12 - a_g * wm.R22 - c_c * (b_g * wm.R11 - a_g * wm.R21 + wm.R31) + wm.R32
                              - a_c * (b_g * wm.R13 - a_g * wm.R23 + wm.R33))
                       - su
                           * (b_g * wm.R11 - a_g * wm.R21 + wm.R31 - c_c * (b_g * wm.R12 - a_g * wm.R22 + wm.R32)
                              + b_c * (b_g * wm.R13 - a_g * wm.R23 + wm.R33)))
                + (-(c_c * wm.R11) + b_g * c_c * tx * wm.R11 + wm.R12 - b_g * tx * wm.R12 - a_c * wm.R13 + a_c * b_g * tx * wm.R13
                   + c_c * c_g * wm.R21 - a_g * c_c * tx * wm.R21 - c_g * wm.R22 + a_g * tx * wm.R22 + a_c * c_g * wm.R23
                   - a_c * a_g * tx * wm.R23 - b_g * c_c * wm.R31 + c_c * tx * wm.R31 + b_g * wm.R32 - tx * wm.R32 - a_c * b_g * wm.R33
                   + a_c * tx * wm.R33)
                    * (by - y_c - y_g - y_l
                       - sz
                           * (-(c_g * wm.R13) + wm.R23 + b_c * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31)
                              + a_c * (-(c_g * wm.R12) + wm.R22 + a_g * wm.R32) + a_g * wm.R33)
                       - sv
                           * (-(c_g * wm.R12) + wm.R22 - c_c * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31) + a_g * wm.R32
                              - a_c * (-(c_g * wm.R13) + wm.R23 + a_g * wm.R33))
                       - su
                           * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31 - c_c * (-(c_g * wm.R12) + wm.R22 + a_g * wm.R32)
                              + b_c * (-(c_g * wm.R13) + wm.R23 + a_g * wm.R33)))
                + (-(c_c * c_g * wm.R11) - b_g * c_c * ty * wm.R11 + c_g * wm.R12 + b_g * ty * wm.R12 - a_c * c_g * wm.R13
                   - a_c * b_g * ty * wm.R13 + c_c * wm.R21 + a_g * c_c * ty * wm.R21 - wm.R22 - a_g * ty * wm.R22 + a_c * wm.R23
                   + a_c * a_g * ty * wm.R23 + a_g * c_c * wm.R31 - c_c * ty * wm.R31 - a_g * wm.R32 + ty * wm.R32 + a_c * a_g * wm.R33
                   - a_c * ty * wm.R33)
                    * (bx - x_c - x_g - x_l
                       - sz
                           * (wm.R13 - c_g * wm.R23 + b_c * (wm.R11 - c_g * wm.R21 + b_g * wm.R31)
                              + a_c * (wm.R12 - c_g * wm.R22 + b_g * wm.R32) + b_g * wm.R33)
                       - sv
                           * (wm.R12 - c_g * wm.R22 - c_c * (wm.R11 - c_g * wm.R21 + b_g * wm.R31) + b_g * wm.R32
                              - a_c * (wm.R13 - c_g * wm.R23 + b_g * wm.R33))
                       - su
                           * (wm.R11 - c_g * wm.R21 + b_g * wm.R31 - c_c * (wm.R12 - c_g * wm.R22 + b_g * wm.R32)
                              + b_c * (wm.R13 - c_g * wm.R23 + b_g * wm.R33)))))
                / (sqrt(pow(-(c_c * wm.R11) + b_g * c_c * tx * wm.R11 + wm.R12 - b_g * tx * wm.R12 - a_c * wm.R13 + a_c * b_g * tx * wm.R13
                                + c_c * c_g * wm.R21 - a_g * c_c * tx * wm.R21 - c_g * wm.R22 + a_g * tx * wm.R22 + a_c * c_g * wm.R23
                                - a_c * a_g * tx * wm.R23 - b_g * c_c * wm.R31 + c_c * tx * wm.R31 + b_g * wm.R32 - tx * wm.R32
                                - a_c * b_g * wm.R33 + a_c * tx * wm.R33,
                            2)
                        + pow(-(c_c * c_g * wm.R11) - b_g * c_c * ty * wm.R11 + c_g * wm.R12 + b_g * ty * wm.R12 - a_c * c_g * wm.R13
                                  - a_c * b_g * ty * wm.R13 + c_c * wm.R21 + a_g * c_c * ty * wm.R21 - wm.R22 - a_g * ty * wm.R22
                                  + a_c * wm.R23 + a_c * a_g * ty * wm.R23 + a_g * c_c * wm.R31 - c_c * ty * wm.R31 - a_g * wm.R32
                                  + ty * wm.R32 + a_c * a_g * wm.R33 - a_c * ty * wm.R33,
                              2)
                        + pow(c_c * c_g * tx * wm.R11 + c_c * ty * wm.R11 - c_g * tx * wm.R12 - ty * wm.R12 + a_c * c_g * tx * wm.R13
                                  + a_c * ty * wm.R13 - c_c * tx * wm.R21 - c_c * c_g * ty * wm.R21 + tx * wm.R22 + c_g * ty * wm.R22
                                  - a_c * tx * wm.R23 - a_c * c_g * ty * wm.R23 - a_g * c_c * tx * wm.R31 + b_g * c_c * ty * wm.R31
                                  + a_g * tx * wm.R32 - b_g * ty * wm.R32 - a_c * a_g * tx * wm.R33 + a_c * b_g * ty * wm.R33,
                              2))
                   * fabs(((c_c * c_g * tx * wm.R11 + c_c * ty * wm.R11 - c_g * tx * wm.R12 - ty * wm.R12 + a_c * c_g * tx * wm.R13
                            + a_c * ty * wm.R13 - c_c * tx * wm.R21 - c_c * c_g * ty * wm.R21 + tx * wm.R22 + c_g * ty * wm.R22
                            - a_c * tx * wm.R23 - a_c * c_g * ty * wm.R23 - a_g * c_c * tx * wm.R31 + b_g * c_c * ty * wm.R31
                            + a_g * tx * wm.R32 - b_g * ty * wm.R32 - a_c * a_g * tx * wm.R33 + a_c * b_g * ty * wm.R33)
                               * (-z_c - z_g - z_l
                                  - sz
                                      * (b_g * wm.R13 - a_g * wm.R23 + b_c * (b_g * wm.R11 - a_g * wm.R21 + wm.R31)
                                         + a_c * (b_g * wm.R12 - a_g * wm.R22 + wm.R32) + wm.R33)
                                  - sv
                                      * (b_g * wm.R12 - a_g * wm.R22 - c_c * (b_g * wm.R11 - a_g * wm.R21 + wm.R31) + wm.R32
                                         - a_c * (b_g * wm.R13 - a_g * wm.R23 + wm.R33))
                                  - su
                                      * (b_g * wm.R11 - a_g * wm.R21 + wm.R31 - c_c * (b_g * wm.R12 - a_g * wm.R22 + wm.R32)
                                         + b_c * (b_g * wm.R13 - a_g * wm.R23 + wm.R33)))
                           + (-(c_c * wm.R11) + b_g * c_c * tx * wm.R11 + wm.R12 - b_g * tx * wm.R12 - a_c * wm.R13
                              + a_c * b_g * tx * wm.R13 + c_c * c_g * wm.R21 - a_g * c_c * tx * wm.R21 - c_g * wm.R22 + a_g * tx * wm.R22
                              + a_c * c_g * wm.R23 - a_c * a_g * tx * wm.R23 - b_g * c_c * wm.R31 + c_c * tx * wm.R31 + b_g * wm.R32
                              - tx * wm.R32 - a_c * b_g * wm.R33 + a_c * tx * wm.R33)
                               * (by - y_c - y_g - y_l
                                  - sz
                                      * (-(c_g * wm.R13) + wm.R23 + b_c * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31)
                                         + a_c * (-(c_g * wm.R12) + wm.R22 + a_g * wm.R32) + a_g * wm.R33)
                                  - sv
                                      * (-(c_g * wm.R12) + wm.R22 - c_c * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31) + a_g * wm.R32
                                         - a_c * (-(c_g * wm.R13) + wm.R23 + a_g * wm.R33))
                                  - su
                                      * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31 - c_c * (-(c_g * wm.R12) + wm.R22 + a_g * wm.R32)
                                         + b_c * (-(c_g * wm.R13) + wm.R23 + a_g * wm.R33)))
                           + (-(c_c * c_g * wm.R11) - b_g * c_c * ty * wm.R11 + c_g * wm.R12 + b_g * ty * wm.R12 - a_c * c_g * wm.R13
                              - a_c * b_g * ty * wm.R13 + c_c * wm.R21 + a_g * c_c * ty * wm.R21 - wm.R22 - a_g * ty * wm.R22 + a_c * wm.R23
                              + a_c * a_g * ty * wm.R23 + a_g * c_c * wm.R31 - c_c * ty * wm.R31 - a_g * wm.R32 + ty * wm.R32
                              + a_c * a_g * wm.R33 - a_c * ty * wm.R33)
                               * (bx - x_c - x_g - x_l
                                  - sz
                                      * (wm.R13 - c_g * wm.R23 + b_c * (wm.R11 - c_g * wm.R21 + b_g * wm.R31)
                                         + a_c * (wm.R12 - c_g * wm.R22 + b_g * wm.R32) + b_g * wm.R33)
                                  - sv
                                      * (wm.R12 - c_g * wm.R22 - c_c * (wm.R11 - c_g * wm.R21 + b_g * wm.R31) + b_g * wm.R32
                                         - a_c * (wm.R13 - c_g * wm.R23 + b_g * wm.R33))
                                  - su
                                      * (wm.R11 - c_g * wm.R21 + b_g * wm.R31 - c_c * (wm.R12 - c_g * wm.R22 + b_g * wm.R32)
                                         + b_c * (wm.R13 - c_g * wm.R23 + b_g * wm.R33))))))
            - ((2
                    * (-(b_g * c_c * wm.R11) + b_g * wm.R12 - a_c * b_g * wm.R13 + a_g * c_c * wm.R21 - a_g * wm.R22 + a_c * a_g * wm.R23
                       - c_c * wm.R31 + wm.R32 - a_c * wm.R33)
                    * (-(c_c * c_g * wm.R11) - b_g * c_c * ty * wm.R11 + c_g * wm.R12 + b_g * ty * wm.R12 - a_c * c_g * wm.R13
                       - a_c * b_g * ty * wm.R13 + c_c * wm.R21 + a_g * c_c * ty * wm.R21 - wm.R22 - a_g * ty * wm.R22 + a_c * wm.R23
                       + a_c * a_g * ty * wm.R23 + a_g * c_c * wm.R31 - c_c * ty * wm.R31 - a_g * wm.R32 + ty * wm.R32 + a_c * a_g * wm.R33
                       - a_c * ty * wm.R33)
                + 2 * (c_c * wm.R11 - wm.R12 + a_c * wm.R13 - c_c * c_g * wm.R21 + c_g * wm.R22 - a_c * c_g * wm.R23 + b_g * c_c * wm.R31 - b_g * wm.R32 + a_c * b_g * wm.R33)
                    * (c_c * c_g * tx * wm.R11 + c_c * ty * wm.R11 - c_g * tx * wm.R12 - ty * wm.R12 + a_c * c_g * tx * wm.R13
                       + a_c * ty * wm.R13 - c_c * tx * wm.R21 - c_c * c_g * ty * wm.R21 + tx * wm.R22 + c_g * ty * wm.R22
                       - a_c * tx * wm.R23 - a_c * c_g * ty * wm.R23 - a_g * c_c * tx * wm.R31 + b_g * c_c * ty * wm.R31 + a_g * tx * wm.R32
                       - b_g * ty * wm.R32 - a_c * a_g * tx * wm.R33 + a_c * b_g * ty * wm.R33))
               * fabs(((c_c * c_g * tx * wm.R11 + c_c * ty * wm.R11 - c_g * tx * wm.R12 - ty * wm.R12 + a_c * c_g * tx * wm.R13
                        + a_c * ty * wm.R13 - c_c * tx * wm.R21 - c_c * c_g * ty * wm.R21 + tx * wm.R22 + c_g * ty * wm.R22
                        - a_c * tx * wm.R23 - a_c * c_g * ty * wm.R23 - a_g * c_c * tx * wm.R31 + b_g * c_c * ty * wm.R31
                        + a_g * tx * wm.R32 - b_g * ty * wm.R32 - a_c * a_g * tx * wm.R33 + a_c * b_g * ty * wm.R33)
                           * (-z_c - z_g - z_l
                              - sz
                                  * (b_g * wm.R13 - a_g * wm.R23 + b_c * (b_g * wm.R11 - a_g * wm.R21 + wm.R31)
                                     + a_c * (b_g * wm.R12 - a_g * wm.R22 + wm.R32) + wm.R33)
                              - sv
                                  * (b_g * wm.R12 - a_g * wm.R22 - c_c * (b_g * wm.R11 - a_g * wm.R21 + wm.R31) + wm.R32
                                     - a_c * (b_g * wm.R13 - a_g * wm.R23 + wm.R33))
                              - su
                                  * (b_g * wm.R11 - a_g * wm.R21 + wm.R31 - c_c * (b_g * wm.R12 - a_g * wm.R22 + wm.R32)
                                     + b_c * (b_g * wm.R13 - a_g * wm.R23 + wm.R33)))
                       + (-(c_c * wm.R11) + b_g * c_c * tx * wm.R11 + wm.R12 - b_g * tx * wm.R12 - a_c * wm.R13 + a_c * b_g * tx * wm.R13
                          + c_c * c_g * wm.R21 - a_g * c_c * tx * wm.R21 - c_g * wm.R22 + a_g * tx * wm.R22 + a_c * c_g * wm.R23
                          - a_c * a_g * tx * wm.R23 - b_g * c_c * wm.R31 + c_c * tx * wm.R31 + b_g * wm.R32 - tx * wm.R32
                          - a_c * b_g * wm.R33 + a_c * tx * wm.R33)
                           * (by - y_c - y_g - y_l
                              - sz
                                  * (-(c_g * wm.R13) + wm.R23 + b_c * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31)
                                     + a_c * (-(c_g * wm.R12) + wm.R22 + a_g * wm.R32) + a_g * wm.R33)
                              - sv
                                  * (-(c_g * wm.R12) + wm.R22 - c_c * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31) + a_g * wm.R32
                                     - a_c * (-(c_g * wm.R13) + wm.R23 + a_g * wm.R33))
                              - su
                                  * (-(c_g * wm.R11) + wm.R21 + a_g * wm.R31 - c_c * (-(c_g * wm.R12) + wm.R22 + a_g * wm.R32)
                                     + b_c * (-(c_g * wm.R13) + wm.R23 + a_g * wm.R33)))
                       + (-(c_c * c_g * wm.R11) - b_g * c_c * ty * wm.R11 + c_g * wm.R12 + b_g * ty * wm.R12 - a_c * c_g * wm.R13
                          - a_c * b_g * ty * wm.R13 + c_c * wm.R21 + a_g * c_c * ty * wm.R21 - wm.R22 - a_g * ty * wm.R22 + a_c * wm.R23
                          + a_c * a_g * ty * wm.R23 + a_g * c_c * wm.R31 - c_c * ty * wm.R31 - a_g * wm.R32 + ty * wm.R32
                          + a_c * a_g * wm.R33 - a_c * ty * wm.R33)
                           * (bx - x_c - x_g - x_l
                              - sz
                                  * (wm.R13 - c_g * wm.R23 + b_c * (wm.R11 - c_g * wm.R21 + b_g * wm.R31)
                                     + a_c * (wm.R12 - c_g * wm.R22 + b_g * wm.R32) + b_g * wm.R33)
                              - sv
                                  * (wm.R12 - c_g * wm.R22 - c_c * (wm.R11 - c_g * wm.R21 + b_g * wm.R31) + b_g * wm.R32
                                     - a_c * (wm.R13 - c_g * wm.R23 + b_g * wm.R33))
                              - su
                                  * (wm.R11 - c_g * wm.R21 + b_g * wm.R31 - c_c * (wm.R12 - c_g * wm.R22 + b_g * wm.R32)
                                     + b_c * (wm.R13 - c_g * wm.R23 + b_g * wm.R33))))))
                / (2
                   * pow(pow(-(c_c * wm.R11) + b_g * c_c * tx * wm.R11 + wm.R12 - b_g * tx * wm.R12 - a_c * wm.R13 + a_c * b_g * tx * wm.R13
                                 + c_c * c_g * wm.R21 - a_g * c_c * tx * wm.R21 - c_g * wm.R22 + a_g * tx * wm.R22 + a_c * c_g * wm.R23
                                 - a_c * a_g * tx * wm.R23 - b_g * c_c * wm.R31 + c_c * tx * wm.R31 + b_g * wm.R32 - tx * wm.R32
                                 - a_c * b_g * wm.R33 + a_c * tx * wm.R33,
                             2)
                             + pow(-(c_c * c_g * wm.R11) - b_g * c_c * ty * wm.R11 + c_g * wm.R12 + b_g * ty * wm.R12 - a_c * c_g * wm.R13
                                       - a_c * b_g * ty * wm.R13 + c_c * wm.R21 + a_g * c_c * ty * wm.R21 - wm.R22 - a_g * ty * wm.R22
                                       + a_c * wm.R23 + a_c * a_g * ty * wm.R23 + a_g * c_c * wm.R31 - c_c * ty * wm.R31 - a_g * wm.R32
                                       + ty * wm.R32 + a_c * a_g * wm.R33 - a_c * ty * wm.R33,
                                   2)
                             + pow(c_c * c_g * tx * wm.R11 + c_c * ty * wm.R11 - c_g * tx * wm.R12 - ty * wm.R12 + a_c * c_g * tx * wm.R13
                                       + a_c * ty * wm.R13 - c_c * tx * wm.R21 - c_c * c_g * ty * wm.R21 + tx * wm.R22 + c_g * ty * wm.R22
                                       - a_c * tx * wm.R23 - a_c * c_g * ty * wm.R23 - a_g * c_c * tx * wm.R31 + b_g * c_c * ty * wm.R31
                                       + a_g * tx * wm.R32 - b_g * ty * wm.R32 - a_c * a_g * tx * wm.R33 + a_c * b_g * ty * wm.R33,
                                   2),
                         3.0 / 2.0));
    }
};

}  // namespace hsa
