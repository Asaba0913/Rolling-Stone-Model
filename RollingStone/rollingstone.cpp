// ============================================================================
//  岩境 / ROCK DYNAMICS LAB  rollingstone3d.cpp
//  ---------------------------------------------------------------------------
//  严格刚体动力学（竖直平面内） + 三维透视渲染
//    · 圆形/球体：纯滚动 (5/7) 系数、打滑判据、滚动阻力、空气阻力
//    · 任意凸多边形：冲量法碰撞、绕棱翻转、能量闸门，惯量由多边形精确积分
//    · 多块体同时运动，块体间碰撞用等效圆盘冲量
//    · 三维：地形高程场 + 透视投影 + 画家算法深度排序 + 朗伯光照 + 大气雾
//  编译：见同目录 vcxproj（VS 直接 F5），或
//        cl -O2 -EHsc -utf-8 -DUNICODE -D_UNICODE -I<EasyX inc> rollingstone3d.cpp ...
//  运行：正常启动为三维演示；命令行加 --test 运行物理自检并输出 selftest.txt
// ============================================================================
#define _USE_MATH_DEFINES
#define NOMINMAX
#include <graphics.h>
#include <conio.h>
#include <tchar.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <vector>
#include <algorithm>

// ---------------------------------------------------------------- 全局常量
static const double G = 9.81;      // 重力加速度 m/s^2
static const int    WIN_W = std::min(1360, GetSystemMetrics(SM_CXSCREEN) - 60);
static const int    WIN_H = std::min(820, GetSystemMetrics(SM_CYSCREEN) - 100);
static const double DT = 0.0005;    // 物理步长 s（经收敛性验证）
static bool g_rain = false, g_snow = false, g_wind = false;

static const double RHO_AIR = 1.20;      // 空气密度 kg/m^3

// ================================================================ 二维向量
struct V2 {
    double x, z;
    V2() : x(0), z(0) {}
    V2(double a, double b) : x(a), z(b) {}
    V2 operator+(const V2& o) const { return V2(x + o.x, z + o.z); }
    V2 operator-(const V2& o) const { return V2(x - o.x, z - o.z); }
    V2 operator*(double k) const { return V2(x * k, z * k); }
    V2 operator/(double k) const { return V2(x / k, z / k); }
    V2& operator+=(const V2& o) { x += o.x; z += o.z; return *this; }
    V2& operator-=(const V2& o) { x -= o.x; z -= o.z; return *this; }
};
static double dot(const V2& a, const V2& b) { return a.x * b.x + a.z * b.z; }
static double crs(const V2& a, const V2& b) { return a.x * b.z - a.z * b.x; }
static double len(const V2& a) { return sqrt(a.x * a.x + a.z * a.z); }
static V2 nrm(const V2& a) { double L = len(a); return L < 1e-12 ? V2(0, 0) : V2(a.x / L, a.z / L); }
// 二维叉积：标量 ω 与矢量 r 作用，得到速度 (ω ẑ × r)
static V2 crossW(double w, const V2& r) { return V2(-w * r.z, w * r.x); }

// ================================================================ 三维向量
struct V3 {
    double x, y, z;
    V3() : x(0), y(0), z(0) {}
    V3(double a, double b, double c) : x(a), y(b), z(c) {}
    V3 operator+(const V3& o) const { return V3(x + o.x, y + o.y, z + o.z); }
    V3 operator-(const V3& o) const { return V3(x - o.x, y - o.y, z - o.z); }
    V3 operator*(double k) const { return V3(x * k, y * k, z * k); }
    V3& operator+=(const V3& o) { x += o.x; y += o.y; z += o.z; return *this; }
};
static double dot(const V3& a, const V3& b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
static V3 crs3(const V3& a, const V3& b) {
    return V3(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x);
}
static double len3(const V3& a) { return sqrt(a.x * a.x + a.y * a.y + a.z * a.z); }
static V3 nrm3(const V3& a) { double L = len3(a); return L < 1e-12 ? V3(0, 0, 0) : V3(a.x / L, a.y / L, a.z / L); }

// ================================================================ 随机数
struct Rng {
    unsigned long long s;
    explicit Rng(unsigned long long seed = 88172645463325252ULL) : s(seed ? seed : 1) {}
    double uni() {
        s ^= s << 13; s ^= s >> 7; s ^= s << 17;
        return (double)((s >> 11) % 9007199254740992ULL) / 9007199254740992.0;
    }
    double range(double a, double b) { return a + (b - a) * uni(); }
    // Box-Muller 标准正态
    double norm() {
        double u1 = uni(); if (u1 < 1e-12) u1 = 1e-12;
        double u2 = uni();
        return sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
    }
};
static Rng g_rng(20241201ULL);

// ============================================================== 地形高程场
//  用保形三次 Hermite（PCHIP，Fritsch–Carlson）插值控制点，
//  保证 C1 连续 —— 这是原折线地形"凹折点当撞击"问题的根本解法：
//  曲率连续后，块体经过坡折处不再产生虚假的法向速度突变。
struct Terrain {
    std::vector<double> kx, kz, km;   // 节点、高程、节点斜率

    void build() {
        // 控制点：山脊（源区）→ 陡坡（45°）→ 缓台（10°）→ 二段陡坡（34°）
        //        → 坡脚 → 堆积缓坡
        //  前 6 个点是山脊背面的下坡，只为视觉上封住地形近端断面：
        //  相机站在山脊后方的高处，背面一路下坡消失在地平线下，
        //  画面上就不会出现"地形忽然被切断"的悬空边。
        static const double X[] = { -150, -118, -92, -70, -52, -36, -22, -12, 0,   10,  24,  38,  52,  64,  76,  88, 100, 112, 126, 146, 176, 220, 280 };
        static const double Z[] = { 0.4,  2.6,  6.5, 14,  26,  40,  52,  62,  70,  60,  47,  35,  30.5, 28.5, 22,  14,  8.0, 4.6, 3.2, 2.4, 1.7, 1.0, 0.4 };
        int n = (int)(sizeof(X) / sizeof(X[0]));
        kx.assign(X, X + n);
        kz.assign(Z, Z + n);
        km.assign(n, 0.0);
        // Fritsch–Carlson 保形斜率
        std::vector<double> d(n - 1);
        for (int i = 0; i + 1 < n; i++) d[i] = (kz[i + 1] - kz[i]) / (kx[i + 1] - kx[i]);
        km[0] = d[0];
        km[n - 1] = d[n - 2];
        for (int i = 1; i + 1 < n; i++) {
            if (d[i - 1] * d[i] <= 0.0) km[i] = 0.0;
            else {
                double w1 = 2.0 * (kx[i + 1] - kx[i]) + (kx[i] - kx[i - 1]);
                double w2 = (kx[i + 1] - kx[i]) + 2.0 * (kx[i] - kx[i - 1]);
                km[i] = (w1 + w2) / (w1 / d[i - 1] + w2 / d[i]);
            }
        }
    }
    int seg(double x) const {
        int n = (int)kx.size();
        if (x <= kx[0]) return 0;
        if (x >= kx[n - 1]) return n - 2;
        int lo = 0, hi = n - 1;
        while (hi - lo > 1) { int mid = (lo + hi) >> 1; if (kx[mid] <= x) lo = mid; else hi = mid; }
        return lo;
    }
    double z(double x) const {
        int i = seg(x);
        double h = kx[i + 1] - kx[i];
        if (x <= kx.front()) return kz.front();
        if (x >= kx.back())  return kz.back();
        double t = (x - kx[i]) / h;
        double t2 = t * t, t3 = t2 * t;
        double h00 = 2 * t3 - 3 * t2 + 1, h10 = t3 - 2 * t2 + t;
        double h01 = -2 * t3 + 3 * t2, h11 = t3 - t2;
        return h00 * kz[i] + h10 * h * km[i] + h01 * kz[i + 1] + h11 * h * km[i + 1];
    }
    double dz(double x) const {
        int i = seg(x);
        double h = kx[i + 1] - kx[i];
        if (x <= kx.front()) return 0.0;
        if (x >= kx.back())  return 0.0;
        double t = (x - kx[i]) / h;
        double t2 = t * t;
        double dh00 = 6 * t2 - 6 * t, dh10 = 3 * t2 - 4 * t + 1;
        double dh01 = -6 * t2 + 6 * t, dh11 = 3 * t2 - 2 * t;
        return (dh00 * kz[i] + dh01 * kz[i + 1]) / h + dh10 * km[i] + dh11 * km[i + 1];
    }
    // 坡面单位法向（指向空中）
    V2 normal(double x) const { return nrm(V2(-dz(x), 1.0)); }
    // 坡角（度）
    double angleDeg(double x) const { return atan(fabs(dz(x))) * 180.0 / M_PI; }
};
static Terrain g_T;

// 把半径 R 的球心放到坡面上：沿坡面法向偏移 R（不是竖直偏移 R）
static V2 placeOnTerrain(double x, double R) {
    V2 n = g_T.normal(x);
    return V2(x, g_T.z(x)) + n * R;
}

// ============================================================== 块体形状库
//  每个形状给出：局部凸多边形顶点（逆时针）、等效半径、三维挤出厚度、
//  面积、质心（归一化到原点）、绕质心转动惯量
enum ShapeId {
    SH_SPHERE = 0, SH_TRI, SH_SQUARE, SH_PENTA, SH_HEXA, SH_OCTA,
    SH_DODECA, SH_IRREG, SH_ELONG, SH_SLAB, SH_BOULDER, SH_COUNT
};
static const TCHAR* SHAPE_NAME[SH_COUNT] = {
    _T("球体"), _T("三棱柱"), _T("正方体"), _T("五棱柱"), _T("六棱柱"), _T("八棱柱"),
    _T("十二棱柱"), _T("不规则块石"), _T("扁长块石"), _T("板状块石"), _T("巨砾")
};

struct Shape {
    std::vector<V2> v;     // 顶点（质心在原点，逆时针）
    double area, I_per_m;  // 面积、单位质量转动惯量
    double thickness;      // 三维挤出厚度
    double rEff;           // 等效圆半径（块体间碰撞用）
    bool   isSphere;
};

// 由顶点表构造形状：把质心平移到原点，并计算面积与惯量
static void finalizeShape(Shape& s) {
    int n = (int)s.v.size();
    double A = 0, cx = 0, cz = 0;
    for (int i = 0; i < n; i++) {
        const V2& p = s.v[i]; const V2& q = s.v[(i + 1) % n];
        double c = crs(p, q);
        A += c;
        cx += (p.x + q.x) * c;
        cz += (p.z + q.z) * c;
    }
    A *= 0.5;
    cx /= (6.0 * A); cz /= (6.0 * A);
    for (int i = 0; i < n; i++) { s.v[i].x -= cx; s.v[i].z -= cz; }
    // 面积与绕质心惯量（三角形分解精确积分）
    double area = 0, I = 0;
    for (int i = 0; i < n; i++) {
        const V2& p = s.v[i]; const V2& q = s.v[(i + 1) % n];
        double c = crs(p, q);
        area += c * 0.5;
        I += c * (p.x * p.x + p.x * q.x + q.x * q.x + p.z * p.z + p.z * q.z + q.z * q.z) / 12.0;
    }
    if (area < 0) { area = -area; I = -I; }
    s.area = area;
    s.I_per_m = I / area;            // 除以面密度 → 单位质量惯量 (m^2)
    double rmax = 0;
    for (int i = 0; i < n; i++) rmax = std::max(rmax, len(s.v[i]));
    s.rEff = s.isSphere ? rmax : rmax * 0.86;
}

static Shape makeShape(int id, double R, double aspect, Rng& rng) {
    Shape s; s.isSphere = false; s.thickness = 2.0 * R;
    int n = 0;
    switch (id) {
    case SH_SPHERE:  n = 3; s.isSphere = true; break;
    case SH_TRI:     n = 3; break;
    case SH_SQUARE:  n = 4; break;
    case SH_PENTA:   n = 5; break;
    case SH_HEXA:    n = 6; break;
    case SH_OCTA:    n = 8; break;
    case SH_DODECA:  n = 12; break;
    case SH_IRREG:   n = 7; break;
    case SH_ELONG:   n = 6; break;
    case SH_SLAB:    n = 5; break;
    case SH_BOULDER: n = 11; break;
    default:         n = 6; break;
    }
    if (s.isSphere) {
        // 球体：用 36 边形近似（渲染走球面着色通道，物理仍是圆）
        n = 36;
        for (int i = 0; i < n; i++) {
            double a = 2 * M_PI * i / n;
            s.v.push_back(V2(R * cos(a), R * sin(a)));
        }
        s.thickness = 2 * R;
    }
    else {
        for (int i = 0; i < n; i++) {
            double a = 2 * M_PI * (i + 0.5) / n;
            double rr = R;
            if (id == SH_IRREG || id == SH_ELONG || id == SH_SLAB || id == SH_BOULDER)
                rr = R * rng.range(0.74, 1.06);          // 天然块石的棱长不均匀
            double ax = rr * cos(a), az = rr * sin(a);
            if (id == SH_ELONG) ax *= 1.62;               // 扁长
            if (id == SH_SLAB) { ax *= 1.95; az *= 0.62; } // 板状
            s.v.push_back(V2(ax, az));
        }
        if (aspect != 1.0) for (size_t i = 0; i < s.v.size(); i++) s.v[i].x *= aspect;
        s.thickness = (id == SH_SLAB) ? R * 0.55 : (id == SH_ELONG ? R * 1.15 : R * 1.5);
    }
    finalizeShape(s);
    return s;
}

// ================================================================== 刚体
struct Body {
    int      shapeId;
    Shape    sh;
    double   mass, inertia;
    V2       pos;        // 质心（x 水平，z 竖直向上）
    double   ang, omega;
    V2       vel;
    double   y;          // 横向位置（车道），m
    double   mu, e, crr, cd;
    COLORREF col, colDark;
    bool     alive, reached;
    double   tReach, vToe;
    std::vector<V2> trail;
    double   trailAcc;
    int      bounces;
    double   maxClear;   // 最大离地高度
    int      airFrames;
    double   energy0;

    double radius() const { return sh.rEff; }
    double kinetic() const { return 0.5 * mass * dot(vel, vel) + 0.5 * inertia * omega * omega; }
    double energy() const { return kinetic() + mass * G * pos.z; }
};

// 把多边形顶点变换到世界系
static void worldVerts(const Body& b, std::vector<V2>& out) {
    out.resize(b.sh.v.size());
    double c = cos(b.ang), s = sin(b.ang);
    for (size_t i = 0; i < b.sh.v.size(); i++) {
        const V2& p = b.sh.v[i];
        out[i] = V2(b.pos.x + c * p.x - s * p.z, b.pos.z + s * p.x + c * p.z);
    }
}

// ---------------------------------------------------------------- 接触求解
//  在接触点 p、法向 n 上施加冲量（法向 + 库仑摩擦），返回法向冲量大小
static double resolveContact(Body& b, const V2& p, const V2& n, double e, double mu) {
    V2 r = p - b.pos;
    V2 vc = b.vel + crossW(b.omega, r);
    double vn = dot(vc, n);
    if (vn > 0) return 0.0;                       // 已在分离
    double rn = crs(r, n);
    double denom = 1.0 / b.mass + rn * rn / b.inertia;
    if (denom < 1e-12) return 0.0;
    double j = -(1.0 + e) * vn / denom;
    if (j < 0) j = 0;
    V2 P = n * j;
    b.vel += P / b.mass;
    b.omega += crs(r, P) / b.inertia;
    // 切向：库仑摩擦，限幅 mu*j
    V2 t(-n.z, n.x);
    vc = b.vel + crossW(b.omega, r);
    double vt = dot(vc, t);
    double rt = crs(r, t);
    double dt2 = 1.0 / b.mass + rt * rt / b.inertia;
    double jt = (dt2 > 1e-12) ? -vt / dt2 : 0.0;
    double lim = mu * j;
    if (jt > lim) jt = lim;
    if (jt < -lim) jt = -lim;
    V2 P2 = t * jt;
    b.vel += P2 / b.mass;
    b.omega += crs(r, P2) / b.inertia;
    return j;
}

// ------------------------------------------- 地形最近点（圆形块体接触用）
//  竖直距离 below = z - Z(x) 只在水平面成立；斜坡上真实垂距是 below·cosθ。
//  这里直接在 [x-1.6R, x+1.6R] 内搜索最近点并三轮细化，
//  得到精确的接触点与法向，消除"过推 1/cosθ"造成的伪能量注入。
struct Hit { V2 p, n; double dist; bool ok; };
static Hit nearestTerrain(double px, double pz, double R) {
    Hit h;
    double a = px - R * 1.7, b = px + R * 1.7;
    int M = 64;
    double best = 1e30, bx = px;
    for (int i = 0; i <= M; i++) {
        double x = a + (b - a) * i / M;
        double dx = x - px, dz = g_T.z(x) - pz;
        double d2 = dx * dx + dz * dz;
        if (d2 < best) { best = d2; bx = x; }
    }
    double w = (b - a) / M;
    for (int r = 0; r < 3; r++) {
        double lo = bx - w, hi = bx + w, nb = bx, nd = best;
        for (int i = -8; i <= 8; i++) {
            double x = lo + (hi - lo) * (i + 8) / 16.0;
            double dx = x - px, dz = g_T.z(x) - pz;
            double d2 = dx * dx + dz * dz;
            if (d2 < nd) { nd = d2; nb = x; }
        }
        bx = nb; best = nd; w /= 8.0;
    }
    h.p = V2(bx, g_T.z(bx));
    h.dist = sqrt(best);
    V2 dir = V2(px - bx, pz - h.p.z);
    h.n = (len(dir) > 1e-9) ? nrm(dir) : g_T.normal(bx);
    if (h.n.z < 0.0) h.n = g_T.normal(bx);
    h.ok = (h.dist < R);
    return h;
}

static long g_nContact = 0, g_nStep = 0;
// ------------------------------------------------------------ 圆形块体步进
static void stepCircle(Body& b, double dt) {
    g_nStep++;
    // 重力 + 空气阻力
    V2 acc(0.0, -G);
    double sp = len(b.vel);
    if (sp > 1e-9 && b.cd > 0) {
        double A = M_PI * b.radius() * b.radius();
        double k = 0.5 * RHO_AIR * b.cd * A / b.mass;
        V2 rel = b.vel - V2(g_wind ? 18.0 : 0.0, 0.0); acc -= rel * (k * len(rel));
    }
    b.vel += acc * dt;
    b.pos += b.vel * dt;
    b.ang += b.omega * dt;

    // 接触检测：搜索真实最近点
    Hit hh = nearestTerrain(b.pos.x, b.pos.z, b.radius());
    double zx = g_T.z(b.pos.x);
    V2 n = hh.n;
    bool contact = hh.ok;

    if (contact) {
        g_nContact++;
        if (b.airFrames * dt > 0.03 && b.vel.z < -1.0) b.bounces++;
        b.airFrames = 0;
        // 位置修正：沿法向推出真实穿透量
        double pen = b.radius() - hh.dist;
        if (pen > 0) b.pos += n * pen;
        V2 cp = hh.p;
        double vnBefore = dot(b.vel, n);
        // 掠射判定：法向相对速度很小 → 属于滚动接触，恢复系数取 0，
        // 仅用法向冲量把法向速度归零、切向由库仑摩擦约束。
        // 这样接触点法向力不做功，机械能损失与 vn^2 同阶，可忽略。
        double eEff = (vnBefore > -0.30) ? 0.0 : b.e;
        resolveContact(b, cp, n, eEff, b.mu * (g_snow ? 0.22 : g_rain ? 0.65 : 1.0));
        // 滚动阻力：作为接触点切向力施加（能量一致）
        V2 t(n.z, -n.x);
        double vt2 = dot(b.vel, t);
        double Nforce = b.mass * G * fabs(n.z);
        double jr = -b.crr * Nforce * dt;
        if (vt2 < 0) jr = -jr;
        V2 Pr = t * jr;
        b.vel += Pr / b.mass;
        b.omega += crs(cp - b.pos, Pr) / b.inertia;
        // 纯滚动约束：ω = -v_t / R（接触点瞬时速度为零，即纯滚动）
        double vt3 = dot(b.vel, t);
        // No forced rolling: Coulomb impulses determine angular velocity.
    }
    else {
        b.airFrames++;
    }
    double clear = hh.dist - b.radius();
    if (clear > b.maxClear) b.maxClear = clear;
}

// -------------------------------------------------------- 凸多边形块体步进
static void stepPoly(Body& b, double dt) {
    V2 acc(0.0, -G);
    double sp = len(b.vel);
    if (sp > 1e-9 && b.cd > 0) {
        double A = b.sh.area;
        double CdA = b.cd * A;
        double k = 0.5 * RHO_AIR * CdA / b.mass;
        V2 rel = b.vel - V2(g_wind ? 18.0 : 0.0, 0.0); acc -= rel * (k * len(rel));
    }
    b.vel += acc * dt;
    b.pos += b.vel * dt;
    b.ang += b.omega * dt;

    std::vector<V2> wv;
    bool touched = false;
    // 少量迭代处理多点接触（面接触稳定）
    for (int it = 0; it < 4; it++) {
        worldVerts(b, wv);
        int deepIdx = -1; double deep = 0.0;
        for (size_t i = 0; i < wv.size(); i++) {
            double h = wv[i].z - g_T.z(wv[i].x);
            if (h < deep) { deep = h; deepIdx = (int)i; }
        }
        if (deepIdx < 0) break;
        touched = true;
        V2 cp = wv[deepIdx];
        V2 n = g_T.normal(cp.x);
        // 位置修正
        double pen = -deep;
        if (pen > 1e-5) b.pos += n * (pen * 0.85);
        resolveContact(b, cp, n, (it == 0 ? b.e : 0.0), b.mu * (g_snow ? 0.22 : g_rain ? 0.65 : 1.0));
        if (it == 0) {
            // 滚动阻力（切向接触力）
            V2 t(n.z, -n.x);
            double vt = dot(b.vel, t);
            double Nforce = b.mass * G * fabs(n.z);
            double jr = -b.crr * Nforce * dt * 1.6;
            if (vt < 0) jr = -jr;
            V2 Pr = t * jr;
            b.vel += Pr / b.mass;
            b.omega += crs(cp - b.pos, Pr) / b.inertia;
        }
    }
    // 角速度阻尼（微观塑性变形耗能，取小值）
    b.omega *= (1.0 - 0.10 * dt);
    if (touched) {
        if (b.airFrames * dt > 0.03 && b.vel.z < -0.8) b.bounces++;
        b.airFrames = 0;
    }
    else {
        b.airFrames++;
    }
    double zx = g_T.z(b.pos.x);
    double clear = b.pos.z - zx - b.sh.rEff;
    if (clear > b.maxClear) b.maxClear = clear;
}

static void stepBody(Body& b, double dt) {
    if (b.sh.isSphere) stepCircle(b, dt);
    else               stepPoly(b, dt);
}

// ============================================================ 三维渲染管线
//  画家算法：所有面（地形、块体、树木、石头）统一进缓冲，按相机空间深度
//  由远及近绘制，因此块体与地形的遮挡关系天然正确。
struct Face {
    V3 p[26];        // 相机空间顶点
    int n;
    double depth;    // 质心深度
    int kind;        // 0=多边形 1=球体 2=椭圆阴影
    COLORREF col;
    double rad;      // 球体屏幕半径
    V3 cen;          // 球心/椭圆中心（相机空间）
};

static std::vector<Face> g_faces;
static Rng g_srng(998244353ULL);

struct Camera {
    V3 pos, fwd, right, up;
    double focal, ox, oy;
    double yaw, pitch, dist;
    V3 target;
} g_cam;

static void updateCamera() {
    V3 dir(cos(g_cam.pitch) * cos(g_cam.yaw), sin(g_cam.yaw) * cos(g_cam.pitch), sin(g_cam.pitch));
    g_cam.pos = g_cam.target + dir * g_cam.dist;
    g_cam.fwd = nrm3(g_cam.target - g_cam.pos);
    V3 worldUp(0, 0, 1);
    g_cam.right = nrm3(crs3(g_cam.fwd, worldUp));
    g_cam.up = crs3(g_cam.right, g_cam.fwd);
    g_cam.focal = WIN_H * 0.95;
    g_cam.ox = WIN_W * 0.5;
    g_cam.oy = WIN_H * 0.55;
}

static bool toCam(const V3& p, V3& o) {
    V3 d = p - g_cam.pos;
    double X = dot(d, g_cam.right);
    double Y = dot(d, g_cam.up);
    double Z = dot(d, g_cam.fwd);
    o = V3(X, Y, Z);
    return Z > 0.6;
}

static COLORREF mixc(COLORREF a, COLORREF b, double t) {
    if (t < 0) t = 0; if (t > 1) t = 1;
    int r = (int)(GetRValue(a) * (1 - t) + GetRValue(b) * t);
    int g = (int)(GetGValue(a) * (1 - t) + GetGValue(b) * t);
    int bl = (int)(GetBValue(a) * (1 - t) + GetBValue(b) * t);
    return RGB(r, g, bl);
}

static const COLORREF C_SKY_TOP = RGB(126, 176, 214);
static const COLORREF C_SKY_HOR = RGB(206, 224, 232);
static V3 g_sunDir;      // 光线传播方向

static COLORREF shade(COLORREF base, double lam, double fog) {
    double k = 0.38 + 0.80 * lam;
    int r = (int)(GetRValue(base) * k); if (r > 255) r = 255;
    int g = (int)(GetGValue(base) * k); if (g > 255) g = 255;
    int b = (int)(GetBValue(base) * k); if (b > 255) b = 255;
    return mixc(RGB(r, g, b), C_SKY_HOR, fog);
}

// 加入一个世界空间多边形面（n ≤ 26）
static void addFace(const V3* wp, int n, COLORREF base, const V3& nrmW, double extraFog = 0.0) {
    // 背面剔除
    V3 cen(0, 0, 0);
    for (int i = 0; i < n; i++) cen += wp[i];
    cen = cen * (1.0 / n);
    if (dot(nrmW, g_cam.pos - cen) <= 0.0) return;
    // ---- 近平面裁剪 ----
    // 不能简单丢弃"含相机后方顶点"的面：那会在近处拉出巨大三角形长条
    // （典型症状：屏幕上出现贯穿画面的绿色/灰色细带）。这里用
    // Sutherland–Hodgman 对 z >= NZ 平面做一次裁剪。
    const double NZ = 0.55;
    V3 cv[26];
    for (int i = 0; i < n; i++) toCam(wp[i], cv[i]);
    V3 cl[30];
    int m = 0;
    for (int i = 0; i < n; i++) {
        const V3& a = cv[i];
        const V3& b = cv[(i + 1) % n];
        bool ina = (a.z >= NZ), inb = (b.z >= NZ);
        if (ina && m < 26) cl[m++] = a;
        if (ina != inb && m < 26) {
            double t = (NZ - a.z) / (b.z - a.z);
            V3 q = a + (b - a) * t;
            q.z = NZ;
            cl[m++] = q;
        }
    }
    if (m < 3) return;
    Face f; f.n = m; f.kind = 0;
    for (int i = 0; i < m; i++) f.p[i] = cl[i];
    double zsum = 0;
    for (int i = 0; i < m; i++) zsum += f.p[i].z;
    f.depth = zsum / m;
    double lam = dot(nrmW, V3(-g_sunDir.x, -g_sunDir.y, -g_sunDir.z));
    if (lam < 0) lam = 0;
    double fog = (f.depth - 150.0) / 1150.0 + extraFog;
    f.col = shade(base, lam, fog);
    g_faces.push_back(f);
}

// 由顶点自动求法向，并翻向相机一侧（用于零厚度片状物，如树冠）
static void addFaceAuto(const V3* wp, int n, COLORREF base) {
    V3 nw(0, 0, 0);
    for (int i = 0; i < n; i++) {
        const V3& a = wp[i];
        const V3& b = wp[(i + 1) % n];
        nw.x += (a.y - b.y) * (a.z + b.z);
        nw.y += (a.z - b.z) * (a.x + b.x);
        nw.z += (a.x - b.x) * (a.y + b.y);
    }
    nw = nrm3(nw);
    V3 cen(0, 0, 0);
    for (int i = 0; i < n; i++) cen += wp[i];
    cen = cen * (1.0 / n);
    if (dot(nw, g_cam.pos - cen) < 0) nw = nw * -1.0;
    addFace(wp, n, base, nw);
}

static void addSphere(const V3& c, double R, COLORREF base) {
    V3 cc;
    if (!toCam(c, cc)) return;
    if (cc.z < R) return;
    Face f; f.n = 0; f.kind = 1;
    f.cen = cc; f.depth = cc.z;
    f.rad = g_cam.focal * R / cc.z;
    f.col = base;
    g_faces.push_back(f);
}

static void addShadow(const V3& c, double rx, double rz, COLORREF base) {
    V3 cc;
    if (!toCam(c, cc)) return;
    Face f; f.n = 0; f.kind = 2;
    f.cen = cc; f.depth = cc.z + 0.05;
    f.rad = g_cam.focal * rx / cc.z;
    f.col = base;
    // 用 rad 传长轴，另存短轴到 p[0].x
    f.p[0].x = g_cam.focal * rz / cc.z;
    g_faces.push_back(f);
}

// ---------------------------------------------------------------- 地形颜色
static double noise2(double x, double y) {
    double s = sin(x * 0.37 + y * 0.11) * 0.5
        + sin(x * 0.13 - y * 0.29) * 0.3
        + sin(x * 1.07 + y * 0.83) * 0.2;
    return s;
}

static COLORREF terrainColor(double x, double y, const V2& n) {
    double slope = fabs(g_T.dz(x));
    double rock = std::min(1.0, slope / 0.85);
    COLORREF grass = RGB(108, 130, 78);
    COLORREF rockC = RGB(136, 126, 113);
    COLORREF gravel = RGB(163, 152, 134);
    COLORREF c;
    if (x < 126) c = mixc(grass, rockC, rock * 0.92 + 0.05);
    else {
        double t = std::min(1.0, (x - 126) / 60.0);
        c = mixc(rockC, gravel, t * 0.75);
    }
    // 细颗粒噪声
    double nz = noise2(x * 1.0, y * 1.0) * 0.5 + noise2(x * 3.1, y * 2.7) * 0.28;
    int dr = (int)(nz * 20.0);
    int r = GetRValue(c) + dr, g = GetGValue(c) + dr, b = GetBValue(c) + dr;
    if (r < 0) r = 0; if (r > 255) r = 255;
    if (g < 0) g = 0; if (g > 255) g = 255;
    if (b < 0) b = 0; if (b > 255) b = 255;
    (void)n;
    return RGB(r, g, b);
}

// --------------------------------------------------------------- 景深辅助
//  按距相机的水平距离生成自适应采样的 x 序列（近密远疏，LOD）
static void buildXGrid(std::vector<double>& xs, double cx) {
    xs.clear(); for (int i = 0;i <= 240;i++) xs.push_back(-30.0 + i);
    (void)cx;
}

// ------------------------------------------------ 谷地横向剖面（纯渲染用）
//  物理是严格二维（沿 x 方向滚落），但这不代表画面只能是"一根沟槽"。
//  这里给地形叠一个横向剖面 h(x, y)：
//    · |y| ≤ 半宽 halfW(x)  → h = 0，与物理模型完全一致，不干扰任何计算
//    · 再往外 → 两侧山壁平滑抬起，到谷肩处收平
//  半宽沿程变大、山壁沿程变矮，于是上游是狭窄深切沟、下游张开成扇形堆积体，
//  地形就从"浮空的长条"变成真正的山谷。
static double valleyHalfW(double x) {
    if (x <= 20.0) return 12.0;
    if (x >= 190.0) return 64.0;
    double t = (x - 20.0) / 170.0;
    return 12.0 + 52.0 * (0.5 - 0.5 * cos(M_PI * t));
}
static double valleyAmp(double x) {
    if (x <= 20.0) return 21.0;
    if (x >= 190.0) return 3.5;
    double t = (x - 20.0) / 170.0;
    return 21.0 - 17.5 * (0.5 - 0.5 * cos(M_PI * t));
}
static double wallH(double x, double y) {
    double d = fabs(y) - valleyHalfW(x);
    if (d <= 0.0) return 0.0;                 // 谷底：与物理模型严格一致
    double t = d / 60.0;
    if (t > 1.0) t = 1.0;
    double h = valleyAmp(x) * (0.5 - 0.5 * cos(M_PI * t));
    // 谷肩以上继续以固定坡率抬升成两侧山体，于是相机是"站在山谷里"，
    // 而不是"看着一块放在地上的板"。
    if (d > 60.0) h += 0.06 * (d - 60.0);
    // 沿沟起伏：打掉那条又直又长的"肩线"，山脊线才有自然感。
    // 振幅随着离谷底的距离增大 —— 远处是起伏的山体，不是一块水平的台地。
    double grow = std::min(1.0, d / 45.0) * (1.0 + d / 110.0);
    h += (2.6 * sin(x * 0.031 + y * 0.042)
        + 1.8 * sin(x * 0.055 - y * 0.075)
        + 1.0 * sin(x * 0.11 + y * 0.19)) * grow;
    return h;
}

// ---------------------------------------------------------------- 绘制地形
//  横向网格按"半宽的倍数"取列，于是网格天然贴着谷形走：
//  u=±1 恰好在谷底与谷壁的交界，往外是山壁，往内是沟底。
static const double VAL_U[23] = { -4.0, -3.3, -2.7, -2.2, -1.8, -1.5, -1.28, -1.10, -0.95, -0.62, -0.30, 0.0,
                                  0.30, 0.62, 0.95, 1.10, 1.28, 1.50, 1.80, 2.2, 2.7, 3.3, 4.0 };
static void buildTerrain(double camX) {
    std::vector<double> xs;
    buildXGrid(xs, camX);
    const int NU = (int)(sizeof(VAL_U) / sizeof(VAL_U[0]));
    for (size_t i = 0; i + 1 < xs.size(); i++) {
        double x0 = xs[i], x1 = xs[i + 1];
        double hw0 = valleyHalfW(x0), hw1 = valleyHalfW(x1);
        for (int j = 0; j + 1 < NU; j++) {
            double ua = VAL_U[j], ub = VAL_U[j + 1];
            double ya0 = ua * 24.0, ya1 = ua * 24.0;
            double yb0 = ub * 24.0, yb1 = ub * 24.0;
            V3 a(x0, ya0, g_T.z(x0) + wallH(x0, ya0));
            V3 b(x1, ya1, g_T.z(x1) + wallH(x1, ya1));
            V3 c(x1, yb1, g_T.z(x1) + wallH(x1, yb1));
            V3 d(x0, yb0, g_T.z(x0) + wallH(x0, yb0));
            V3 n = nrm3(crs3(b - a, d - a));
            if (n.z < 0) n = n * -1.0;
            double xc = (x0 + x1) * 0.5, yc = (ya0 + yb0) * 0.5;
            COLORREF col = terrainColor(xc, yc, V2(n.x, n.z));
            // 山壁越陡、离谷底越远，越偏岩色（植被向谷底聚集）
            double over = fabs(yc) / valleyHalfW(xc);
            double rock = std::min(1.0, (over - 0.85) / 0.9);
            if (rock > 0) col = mixc(col, RGB(128, 120, 108), rock * 0.62);
            // 谷底冲积带：偏砂砾
            if (over < 0.55) col = mixc(col, RGB(150, 142, 126), (0.55 - over) * 0.55);
            // 横向越远、雾越重：让网格外缘自然溶进远山与天空，
            // 而不是像一块板那样留一条清楚的边
            double extra = (fabs(yc) - 110.0) / 620.0;
            if (extra < 0.0) extra = 0.0;
            if (extra > 0.45) extra = 0.45;
            V3 q[4] = { a, b, c, d };
            addFace(q, 4, col, n, extra);
        }
    }
}

// ---------------------------------------------------------------- 绘制块体
static void drawBodyFaces(const Body& b, int idx) {
    COLORREF base = b.col;
    if (b.sh.isSphere) {
        V3 c(b.pos.x, b.y, b.pos.z);
        // 接触阴影
        double zx = g_T.z(b.pos.x);
        addShadow(V3(b.pos.x, b.y, zx + 0.02), b.sh.rEff * 0.95, b.sh.rEff * 0.34, RGB(96, 88, 72));
        addSphere(c, b.sh.rEff, base);
        return;
    }
    // 挤出棱柱：两个端面 + 侧面
    std::vector<V2> wv;
    worldVerts(b, wv);
    int n = (int)wv.size();
    double half = b.sh.thickness * 0.5;
    // 端面（y = b.y ± half）
    for (int side = 0; side < 2; side++) {
        double yy = b.y + (side ? half : -half);
        V3* poly = new V3[n];
        for (int i = 0; i < n; i++) poly[i] = V3(wv[i].x, yy, wv[i].z);
        V3 nrm(0, side ? 1.0 : -1.0, 0);
        // 材质：铝灰岩，端面稍亮
        COLORREF c = (side ? mixc(base, RGB(255, 255, 255), 0.08)
            : mixc(base, RGB(0, 0, 0), 0.10));
        addFace(poly, n > 26 ? 26 : n, c, nrm);
        delete[] poly;
    }
    // 侧面
    for (int i = 0; i < n; i++) {
        int j = (i + 1) % n;
        V2 e = wv[j] - wv[i];
        V2 outN(e.z, -e.x);
        outN = nrm(outN);
        V3 p0(wv[i].x, b.y - half, wv[i].z);
        V3 p1(wv[j].x, b.y - half, wv[j].z);
        V3 p2(wv[j].x, b.y + half, wv[j].z);
        V3 p3(wv[i].x, b.y + half, wv[i].z);
        V3 q[4] = { p0, p1, p2, p3 };
        V3 nrmW(outN.x, 0, outN.z);
        // 侧面法向必须指向外侧
        V3 cen = (p0 + p1 + p2 + p3) * 0.25;
        V3 radial = cen - V3(b.pos.x, b.y, b.pos.z);
        if (dot(nrmW, radial) < 0) nrmW = nrmW * -1.0;
        addFace(q, 4, base, nrmW);
    }
    (void)idx;
}

// -------------------------------------------------------------- 散布景物
struct Prop { V3 p; double r; int kind; COLORREF c; };
static std::vector<Prop> g_props;

static void buildProps() {
    if (!g_props.empty()) return;
    Rng r(13579ULL);
    for (int i = 0; i < 420; i++) {
        double x = r.range(-24.0, 204.0);
        double hw = valleyHalfW(x);
        double side = (r.uni() < 0.5) ? -1.0 : 1.0;
        double off = r.range(2.0, 62.0);
        double y = side * (hw + off);
        if (fabs(y) > 87.0) continue;
        Prop p;
        // 关键：必须落在"含横向剖面"的真实地表上，否则石块会浮在空中 / 陷进山壁
        p.p = V3(x, y, g_T.z(x) + wallH(x, y));
        // 谷壁以上的坡面才有树（谷底是滚石通道，不长树）
        bool tree = (off > 12.0) && (r.uni() < 0.60);
        if (tree) { p.kind = 1; p.r = r.range(2.2, 4.2); }
        else {
            p.kind = 0;
            p.r = r.range(0.35, 1.6) * (off > 10 ? 1.2 : 0.85);
            if (x > 120 && r.uni() < 0.6) p.r *= 0.7;
        }
        int t = (int)r.range(0, 3);
        p.c = (t == 0) ? RGB(132, 126, 116) : (t == 1 ? RGB(112, 104, 96) : RGB(148, 142, 130));
        g_props.push_back(p);
    }
}

static void buildSection() {
    const double base = -18.0;
    for (int side = -1;side <= 1;side += 2) {
        double y = side * 96.0;
        for (int i = 0;i < 240;i++) {
            double x = -30.0 + i, z0 = g_T.z(x) + wallH(x, y), z1 = g_T.z(x + 1) + wallH(x + 1, y);
            for (int layer = 0;layer < 4;layer++) {
                double t0 = layer / 4.0, t1 = (layer + 1) / 4.0;
                V3 q[4] = { V3(x,y,base + (z0 - base) * t0),V3(x + 1,y,base + (z1 - base) * t0),V3(x + 1,y,base + (z1 - base) * t1),V3(x,y,base + (z0 - base) * t1) };
                COLORREF colors[4] = { RGB(68,57,52),RGB(110,85,65),RGB(149,117,81),RGB(174,145,102) };
                addFace(q, 4, colors[layer], V3(0, side, 0));
            }
        }
    }
    for (int end = 0;end < 2;end++) {
        double x = end ? 210 : -30;
        for (int j = 0;j < 96;j++) {
            double y = -96 + j * 2, z0 = g_T.z(x) + wallH(x, y), z1 = g_T.z(x) + wallH(x, y + 2);
            V3 q[4] = { V3(x,y,base),V3(x,y + 2,base),V3(x,y + 2,z1),V3(x,y,z0) };
            addFace(q, 4, RGB(130, 103, 77), V3(end ? 1 : -1, 0, 0));
        }
    }
    V3 bottom[4] = { V3(-30,-96,base),V3(210,-96,base),V3(210,96,base),V3(-30,96,base) };
    addFace(bottom, 4, RGB(60, 53, 50), V3(0, 0, -1));
}
static void drawProps(double camX) {
    for (size_t i = 0; i < g_props.size(); i++) {
        const Prop& p = g_props[i];
        if (p.kind == 0) continue; // Hide decorative box-shaped rocks.


        if (p.kind == 0) {
            // 散落石块：带侧面的小方块，看起来有体积
            double R = p.r * 1.15;
            double hz = p.p.z + R * 0.55;
            double zbase = p.p.z - 1.2;
            V3 top[4] = {
                V3(p.p.x - R, p.p.y - R * 0.8, hz),
                V3(p.p.x + R, p.p.y - R * 0.8, hz),
                V3(p.p.x + R, p.p.y + R * 0.8, hz),
                V3(p.p.x - R, p.p.y + R * 0.8, hz) };
            V3 bse[4] = {
                V3(p.p.x - R, p.p.y - R * 0.8, zbase),
                V3(p.p.x + R, p.p.y - R * 0.8, zbase),
                V3(p.p.x + R, p.p.y + R * 0.8, zbase),
                V3(p.p.x - R, p.p.y + R * 0.8, zbase) };
            for (int e = 0; e < 4; e++) {
                int f2 = (e + 1) % 4;
                V3 q[4] = { bse[e], bse[f2], top[f2], top[e] };
                addFaceAuto(q, 4, mixc(p.c, RGB(0, 0, 0), 0.18));
            }
            addFaceAuto(top, 4, p.c);
        }
        else {
            // 针叶树：三层面向相机的"广告牌"三角形 —— 从任何角度看都像树，
            // 不会像固定朝向的薄片那样在某些视角变成一条细线。
            // 往下压一点：网格在列与列之间是线性插值，而道具是按解析高度放的，
            // 两者会差几米，压下去才能保证树根"扎"在地里而不是浮在半空。
            V3 b = p.p; b.z -= 3.0;
            double H = p.r;
            for (int k = 0; k < 3; k++) {
                double hb = H * 0.30 * k;
                double hh = H * (0.60 - 0.09 * k);
                double w = H * (0.30 - 0.055 * k);
                V3 l = b - g_cam.right * w; l.z = b.z + hb;
                V3 rr = b + g_cam.right * w; rr.z = b.z + hb;
                V3 ap(b.x, b.y, b.z + hb + hh);
                V3 tri[3] = { l, rr, ap };
                addFaceAuto(tri, 3, k == 2 ? RGB(64, 104, 68)
                    : (k == 1 ? RGB(52, 90, 58) : RGB(44, 80, 50)));
            }
        }
    }
}

// -------------------------------------------------------------- 天空与远山
//  画面构成（远 → 近）：天空渐变 → 远山剪影 → 地平线以下的谷地雾霭 → 三维地形。
//  "地平线以下的谷地"这一层是关键：没有它，地形边缘之外会直接露出天空色，
//   视觉上就变成"一块浮空的板"；补上之后地形像是从大地里长出来的。
static void drawSky() {
    for (int y = 0; y < WIN_H; y++) {
        double t = (double)y / WIN_H;
        double k = pow(t, 0.75);
        COLORREF c = mixc(C_SKY_TOP, C_SKY_HOR, k);
        setlinecolor(c);
        line(0, y, WIN_W, y);
    }
    // 地平线屏幕高度：由相机俯仰推出来，而不是写死
    double fz = g_cam.fwd.z;
    if (fz > 0.96) fz = 0.96;
    if (fz < -0.96) fz = -0.96;
    double hy = g_cam.oy + g_cam.focal * fz / sqrt(1.0 - fz * fz);

    // ---- 远山：三层越远越淡 ----
    for (int layer = 0; layer < 3; layer++) {
        double baseY = hy + 6.0 + layer * 9.0;
        double amp = 96.0 - layer * 26.0;
        COLORREF c = mixc(RGB(158, 174, 178), RGB(202, 212, 212), 0.30 + layer * 0.30);
        POINT pts[64];
        int n = 0;
        for (int i = 0; i <= 40; i++) {
            double t = (double)i / 40.0;
            pts[n].x = (LONG)(t * WIN_W);
            double yy = baseY
                - amp * (0.50 + 0.50 * sin(t * 5.9 + layer * 2.3))
                - amp * 0.30 * sin(t * 15.3 + layer * 4.1);
            pts[n].y = (LONG)yy; n++;
        }
        pts[n].x = WIN_W; pts[n].y = WIN_H; n++;
        pts[n].x = 0;     pts[n].y = WIN_H; n++;
        setfillcolor(c);
        solidpolygon(pts, n);
    }
    // ---- 地平线以下的谷地雾霭（由近及远变淡） ----
    int y0 = (int)hy;
    if (y0 < 0) y0 = 0;
    for (int y = y0; y < WIN_H; y++) {
        double t = (double)(y - y0) / (WIN_H - y0 + 1.0);   // 0=地平线 1=画面底
        double k = pow(1.0 - t, 1.45);
        COLORREF c = mixc(RGB(112, 126, 100), RGB(170, 182, 178), k);
        setlinecolor(c);
        line(0, y, WIN_W, y);
    }
}

// ---------------------------------------------------------------- 面绘制
static void paintFaces() {
    std::sort(g_faces.begin(), g_faces.end(),
        [](const Face& a, const Face& b) { return a.depth > b.depth; });
    for (size_t i = 0; i < g_faces.size(); i++) {
        const Face& f = g_faces[i];
        if (f.kind == 0) {
            POINT pts[26];
            double cx = 0, cy = 0;
            for (int j = 0; j < f.n; j++) {
                double sx = g_cam.ox + g_cam.focal * f.p[j].x / f.p[j].z;
                double sy = g_cam.oy - g_cam.focal * f.p[j].y / f.p[j].z;
                pts[j].x = (LONG)sx; pts[j].y = (LONG)sy;
                cx += sx; cy += sy;
            }
            cx /= f.n; cy /= f.n;
            for (int j = 0; j < f.n; j++) {     // 轻微外扩，消除面缝
                pts[j].x = (LONG)(cx + (pts[j].x - cx) * 1.035);
                pts[j].y = (LONG)(cy + (pts[j].y - cy) * 1.035);
            }
            setfillcolor(f.col);
            solidpolygon(pts, f.n);
        }
        else if (f.kind == 1) {
            double sx = g_cam.ox + g_cam.focal * f.cen.x / f.cen.z;
            double sy = g_cam.oy - g_cam.focal * f.cen.y / f.cen.z;
            double R = f.rad;
            int steps = (int)(R * 1.6); if (steps < 6) steps = 6; if (steps > 42) steps = 42;
            // 由外向内画同心圆模拟球面光照
            for (int k = steps; k >= 1; k--) {
                double t = (double)k / steps;
                double rr = R * t;
                double lam = 0.30 + 0.95 * pow(1.0 - t, 0.55);
                double ox = -R * 0.30 * (1.0 - t);
                double oy = -R * 0.34 * (1.0 - t);
                int r = (int)(GetRValue(f.col) * lam); if (r > 255) r = 255;
                int g = (int)(GetGValue(f.col) * lam); if (g > 255) g = 255;
                int b = (int)(GetBValue(f.col) * lam); if (b > 255) b = 255;
                setfillcolor(RGB(r, g, b));
                solidcircle((int)(sx + ox), (int)(sy + oy), (int)(rr + 0.5));
            }
        }
        else {
            double sx = g_cam.ox + g_cam.focal * f.cen.x / f.cen.z;
            double sy = g_cam.oy - g_cam.focal * f.cen.y / f.cen.z;
            double rx = f.rad, rz = f.p[0].x;
            if (rx < 1.2 || rz < 0.8) continue;
            if (f.cen.z < 6.0 || rx > WIN_W * 0.32) continue;
            setfillcolor(f.col);
            solidellipse((int)(sx - rx), (int)(sy - rz), (int)(sx + rx), (int)(sy + rz));
        }
    }
}

// ==================================================================== 世界
struct World {
    std::vector<Body> bodies;
    double t;
    bool   paused;
    double speed;        // 时间倍率
    bool   follow, trails;
    double toeX;         // 坡脚参考位置
    int    frames;
    double fps;
} g_world;

static const COLORREF PALETTE[12] = {
    RGB(196, 118, 92), RGB(116, 146, 176), RGB(178, 156, 96),
    RGB(132, 158, 134), RGB(172, 128, 168), RGB(120, 172, 178),
    RGB(198, 150, 108), RGB(148, 132, 172), RGB(200, 176, 132),
    RGB(140, 176, 140), RGB(184, 132, 132), RGB(126, 150, 170)
};

static bool g_shapeEnabled[SH_COUNT] = { false,false,false,false,false,false,false,true,false,false,false };
static Body makeBody(int shapeId, double R, double laneY, int colorIdx, Rng& rng) {
    Body b;
    b.shapeId = shapeId;
    b.sh = makeShape(shapeId, R, 1.0, rng);
    double rho = 2650.0;                       // 花岗岩密度 kg/m^3
    if (b.sh.isSphere) {
        double vol = 4.0 / 3.0 * M_PI * R * R * R;
        b.mass = rho * vol;
        b.inertia = 0.4 * b.mass * R * R;      // 实心球 I = 2/5 mR²
    }
    else {
        b.mass = rho * b.sh.area * b.sh.thickness;
        b.inertia = b.mass * b.sh.I_per_m;
    }
    b.pos = placeOnTerrain(1.0, R);
    b.ang = rng.range(0.0, 2 * M_PI);
    b.omega = 0;
    b.vel = V2(0.0, 0.0);
    b.y = laneY;
    b.mu = 0.55; b.e = 0.35; b.crr = 0.045; b.cd = 0.55;
    b.col = PALETTE[colorIdx % 12];
    b.colDark = mixc(b.col, RGB(0, 0, 0), 0.3);
    b.alive = true; b.reached = false;
    b.tReach = -1; b.vToe = 0;
    b.trailAcc = 0; b.bounces = 0; b.maxClear = 0; b.airFrames = 0;
    b.energy0 = b.energy();
    return b;
}

static void initWorld() {
    g_T.build();
    g_world.bodies.clear();
    Rng rng(20241201ULL);
    // 11 种形状，全部同时起跑，各自一条车道
    double lanes[11] = { -9.0, -5.4, -1.9, 1.9, 5.4, 9.0, -7.2, -3.6, 3.6, 7.2, 0.0 };
    double radii[11] = { 1.15, 1.05, 1.00, 0.95, 0.90, 0.85, 0.80, 0.95, 1.05, 1.10, 1.30 };
    for (int i = 0; i < SH_COUNT; i++) {
        if (!g_shapeEnabled[i]) continue;
        g_world.bodies.push_back(makeBody(i, radii[i], lanes[i], i, rng));
    }
    g_world.t = 0;
    g_world.paused = true;
    g_world.speed = 1.0;
    g_world.follow = false;
    g_world.trails = false;
    g_world.toeX = 112.0;
    g_world.frames = 0;
    g_world.fps = 0;
}

// --------------------------------------------------------- 块体间碰撞
static void collideBodies() {
    int n = (int)g_world.bodies.size();
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            Body& a = g_world.bodies[i];
            Body& b = g_world.bodies[j];
            if (fabs(a.y - b.y) > (a.radius() + b.radius()) * 1.05) continue;
            V2 d = b.pos - a.pos;
            double dist = len(d);
            double sum = a.radius() + b.radius();
            if (dist > sum || dist < 1e-9) continue;
            V2 n = d / dist;
            V2 rel = b.vel - a.vel;
            double vn = dot(rel, n);
            double e = 0.35;
            if (vn < 0) {
                double j = -(1.0 + e) * vn / (1.0 / a.mass + 1.0 / b.mass);
                a.vel -= n * (j / a.mass);
                b.vel += n * (j / b.mass);
                // 摩擦带来的自转
                a.omega -= 0.12 * j / (a.mass * a.radius());
                b.omega += 0.12 * j / (b.mass * b.radius());
            }
            // 分离
            double pen = sum - dist;
            if (pen > 0) {
                double w = a.mass / (a.mass + b.mass);
                a.pos -= n * (pen * w * 0.8);
                b.pos += n * (pen * (1 - w) * 0.8);
            }
        }
    }
}

static void worldStep(double dt) {
    for (size_t i = 0; i < g_world.bodies.size(); i++) {
        Body& b = g_world.bodies[i];
        if (!b.alive) continue;
        stepBody(b, dt);
        if (!b.reached && b.pos.x >= g_world.toeX) {
            b.reached = true;
            b.tReach = g_world.t;
            b.vToe = len(b.vel);
        }
        // 轨迹
        if (g_world.trails) {
            b.trailAcc += dt;
            if (b.trailAcc > 0.06) {
                b.trailAcc = 0;
                b.trail.push_back(V2(b.pos.x, b.pos.z));
                if (b.trail.size() > 900) b.trail.erase(b.trail.begin());
            }
        }
    }
    collideBodies();
    // 出域处理
    for (size_t i = 0; i < g_world.bodies.size(); i++) {
        Body& b = g_world.bodies[i];
        if (b.pos.x > 380.0 || b.pos.z < -30.0) { b.vel = V2(0, 0); b.omega = 0; }
    }
}

// ================================================================== 界面
static void txt(int x, int y, const TCHAR* s, COLORREF c, int size, bool bold = false) {
    LOGFONT lf;
    gettextstyle(&lf);
    lf.lfHeight = size;
    lf.lfWidth = 0;
    lf.lfWeight = bold ? FW_BOLD : FW_NORMAL;
    lf.lfQuality = ANTIALIASED_QUALITY;
    lf.lfItalic = 0;
    _tcscpy_s(lf.lfFaceName, LF_FACESIZE, _T("微软雅黑"));
    setbkmode(TRANSPARENT);
    settextstyle(&lf);
    settextcolor(c);
    outtextxy(x, y, s);
}

static void panel(int x, int y, int w, int h, COLORREF fill, double alpha) {
    (void)alpha;
    setfillcolor(fill);
    setlinecolor(mixc(fill, RGB(255, 255, 255), 0.25));
    solidrectangle(x, y, x + w, y + h);
}

// Mouse-operated shape and weather switches. New shapes spawn at the source.
static bool hitBox(int x, int y, int l, int t, int w, int h) {
    return x >= l && x < l + w && y >= t && y < t + h;
}
static void choiceInput() {
    // Drain EasyX mouse events every frame, preserving the existing polling controls.
    ExMessage event;
    while (peekmessage(&event, EX_MOUSE)) {
        if (event.message == WM_MOUSEWHEEL) {
            g_cam.dist *= exp(-double(event.wheel) / 120.0 * 0.12);
            g_cam.dist = std::max(12.0, std::min(400.0, g_cam.dist));
        }
    }
    static bool wasDown = false;
    bool down = (GetAsyncKeyState(VK_LBUTTON) & 0x8000) != 0;
    if (down && !wasDown) {
        POINT p; GetCursorPos(&p); ScreenToClient(GetHWnd(), &p);
        if (hitBox(p.x, p.y, WIN_W - 290, 380, 128, 34))g_world.paused = false;
        if (hitBox(p.x, p.y, WIN_W - 154, 380, 128, 34))g_world.paused = true;
        for (int i = 0;i < SH_COUNT;i++) {
            int x = WIN_W - 290 + (i % 2) * 136, y = 76 + (i / 2) * 35;
            if (hitBox(p.x, p.y, x, y, 128, 29)) {
                g_shapeEnabled[i] = !g_shapeEnabled[i];
                if (!g_shapeEnabled[i]) {
                    for (size_t j = g_world.bodies.size();j > 0;--j)
                        if (g_world.bodies[j - 1].shapeId == i) g_world.bodies.erase(g_world.bodies.begin() + j - 1);
                }
                else {
                    Rng rng(20261115ULL + i);
                    g_world.bodies.push_back(makeBody(i, 1.0, (i - 5) * 1.7, i, rng));
                }
            }
        }
        if (hitBox(p.x, p.y, WIN_W - 290, 334, 80, 32))g_rain = !g_rain;
        if (hitBox(p.x, p.y, WIN_W - 202, 334, 80, 32))g_snow = !g_snow;
        if (hitBox(p.x, p.y, WIN_W - 114, 334, 80, 32))g_wind = !g_wind;
    }
    static POINT last = { 0,0 }; static bool dragging = false;
    POINT now; GetCursorPos(&now); ScreenToClient(GetHWnd(), &now);
    if (down && !wasDown)dragging = now.x > 320 && now.x < WIN_W - 310 && now.y>60 && now.y < WIN_H - 40;
    if (down && wasDown && dragging) {
        g_cam.yaw -= (now.x - last.x) * 0.007;
        g_cam.pitch = std::max(0.12, std::min(1.35, g_cam.pitch + (now.y - last.y) * 0.005));
        g_world.follow = false;
    }
    if (!down)dragging = false;
    last = now;
    wasDown = down;
}
static void drawWeather() {
    static LARGE_INTEGER freq = { 0 }, last = { 0 }; static double t = 0;
    LARGE_INTEGER now; QueryPerformanceCounter(&now);
    if (!freq.QuadPart) { QueryPerformanceFrequency(&freq);last = now; }
    double dt = double(now.QuadPart - last.QuadPart) / freq.QuadPart; last = now;
    if (!g_world.paused)t += std::min(dt, 0.05);
    if (g_rain) {
        setlinecolor(RGB(125, 195, 235)); setlinestyle(PS_SOLID, 1);
        for (int i = 0;i < 240;i++) {
            double speed = 380 + (i % 7) * 48;
            int y = int(fmod(i * 71.73 + t * speed, WIN_H + 40)) - 20;
            int x = int(fmod(i * 137.31 + t * (g_wind ? 210 : 35), WIN_W + 100)) - 50;
            line(x, y, x + (g_wind ? 10 : 2), y + 10 + i % 9);
        }
    }
    if (g_snow) {
        setfillcolor(RGB(230, 243, 255));
        for (int i = 0;i < 150;i++) {
            int y = int(fmod(i * 91.13 + t * (28 + i % 35), WIN_H + 20)) - 10;
            int x = int(fmod(i * 117.31 + sin(t + i) * 14 + t * (g_wind ? 110 : 12) + 100, WIN_W + 40)) - 20;
            solidcircle(x, y, 1 + i % 3);
        }
    }
    if (g_wind) {
        setlinecolor(RGB(150, 207, 218)); setlinestyle(PS_SOLID, 1);
        for (int i = 0;i < 26;i++) {
            int x = int(fmod(i * 153.7 + t * 320, WIN_W + 150)) - 150;
            int y = 80 + (i * 83) % 680 + int(sin(t * 2 + i) * 12);
            line(x, y, x + 42 + i % 30, y - 5);
        }
    }
}
static void drawChoices() {
    panel(WIN_W - 302, 62, 290, 366, RGB(8, 18, 30), 1);
    for (int i = 0;i < SH_COUNT;i++) {
        int x = WIN_W - 290 + (i % 2) * 136, y = 76 + (i / 2) * 35;
        bool on = g_shapeEnabled[i];
        panel(x, y, 128, 29, on ? RGB(12, 76, 80) : RGB(24, 32, 47), 1);
        TCHAR text[80];_stprintf_s(text, 80, _T("%s %s"), on ? _T("●") : _T("○"), SHAPE_NAME[i]);
        txt(x + 7, y + 5, text, on ? RGB(80, 255, 225) : RGB(160, 175, 195), 13);
    }
    txt(WIN_W - 290, 300, _T("天气 / 可叠加 · 暂停时动画冻结"), RGB(80, 230, 220), 14);
    panel(WIN_W - 290, 380, 128, 34, RGB(12, 90, 80), 1);
    txt(WIN_W - 253, 388, _T("开始"), RGB(100, 255, 225), 15, true);
    panel(WIN_W - 154, 380, 128, 34, RGB(70, 44, 30), 1);
    txt(WIN_W - 117, 388, _T("暂停"), RGB(255, 200, 125), 15, true);
    const TCHAR* names[3] = { _T("雨滴"),_T("飘雪"),_T("大风") };
    bool states[3] = { g_rain,g_snow,g_wind };
    for (int i = 0;i < 3;i++) {
        int x = WIN_W - 290 + i * 88;
        panel(x, 334, 80, 32, states[i] ? RGB(12, 76, 80) : RGB(24, 32, 47), 1);
        txt(x + 16, 341, names[i], states[i] ? RGB(80, 255, 225) : RGB(160, 175, 195), 14);
    }
}
static void drawHUD() {
    TCHAR buf[256];
    // 顶部标题条
    panel(0, 0, WIN_W, 54, RGB(8, 12, 25), 1.0);
    txt(20, 12, _T("岩境 / ROCK DYNAMICS LAB"), RGB(238, 244, 250), 22, true);
    _stprintf_s(buf, 200, _T("t = %6.2f s     时间倍率 x%.1f     %2.0f FPS"),
        g_world.t, g_world.speed, g_world.fps);
    txt(WIN_W - 400, 16, buf, RGB(61, 230, 215), 16);

    // 左侧块体信息
    int panelW = 300;
    int listTop = 74;
    int rowH = 30;
    int rows = (int)g_world.bodies.size();
    panel(12, listTop, panelW, 40 + rows * rowH, RGB(10, 19, 33), 0.92);
    txt(24, listTop + 10, _T("块体状态（同一山坡、同一时刻）"), RGB(180, 200, 218), 14, true);
    for (int i = 0; i < rows; i++) {
        const Body& b = g_world.bodies[i];
        int yy = listTop + 36 + i * rowH;
        setfillcolor(b.col);
        solidrectangle(24, yy + 6, 40, yy + 22);
        setlinecolor(mixc(b.col, RGB(255, 255, 255), 0.4));
        rectangle(24, yy + 6, 40, yy + 22);
        _stprintf_s(buf, 200, _T("%-7s"), SHAPE_NAME[b.shapeId]);
        txt(48, yy + 5, buf, RGB(226, 236, 246), 14);
        _stprintf_s(buf, 200, _T("%5.1f m/s  x=%5.1f m  离地%.2f m"),
            len(b.vel), b.pos.x, b.maxClear);
        txt(120, yy + 5, buf, RGB(146, 170, 192), 13);
    }

    // 底部按键提示
    panel(0, WIN_H - 34, WIN_W, 34, RGB(20, 28, 38), 1.0);
    txt(14, WIN_H - 29,
        _T("空格 暂停/继续   R 重新开始   +/- 速度   F 跟随镜头   T 轨迹   ")
        _T("A/D 旋转  W/S 俯仰  Q/E 远近   N 加一块体   C 清空   Esc 退出"),
        RGB(158, 182, 204), 13);

    // 坡脚标注
    if (g_world.bodies.size() > 0) {
        V3 toe(g_world.toeX, -20.0, g_T.z(g_world.toeX) + 0.2);
        V3 cam;
        if (toCam(toe, cam)) {
            double sx = g_cam.ox + g_cam.focal * cam.x / cam.z;
            double sy = g_cam.oy - g_cam.focal * cam.y / cam.z;
            if (sx > -50 && sx < WIN_W + 50) {
                setlinecolor(RGB(226, 132, 96));
                setlinestyle(PS_DASH, 2);
                line((int)sx, (int)sy - 8, (int)sx, (int)sy - 62);
                setlinestyle(PS_SOLID, 1);
                setfillcolor(RGB(226, 132, 96));
                POINT tri[3] = { { (LONG)sx, (LONG)sy },{ (LONG)sx - 6, (LONG)sy - 12 },{ (LONG)sx + 6, (LONG)sy - 12 } };
                solidpolygon(tri, 3);
                txt((int)sx + 9, (int)sy - 62, _T("坡脚"), RGB(226, 132, 96), 14, true);
            }
        }
    }
}

// ================================================================ 物理自检
//  输出统一用 UTF-8 写盘，避免控制台/文本编辑器编码混乱
static FILE* g_sf = NULL;
static void tlog(const TCHAR* fmt, ...) {
    if (!g_sf) return;
    TCHAR wbuf[1024];
    va_list ap;
    va_start(ap, fmt);
    _vsnwprintf_s(wbuf, 1024, _TRUNCATE, fmt, ap);
    va_end(ap);
    char u8[4096];
    int n = WideCharToMultiByte(CP_UTF8, 0, wbuf, -1, u8, 4096, NULL, NULL);
    if (n > 1) fwrite(u8, 1, n - 1, g_sf);
}

static void selftest() {
    fopen_s(&g_sf, "selftest.txt", "wb");
    if (!g_sf) return;
    // UTF-8 BOM，便于记事本识别
    fwrite("\xEF\xBB\xBF", 1, 3, g_sf);

    tlog(_T("=== 山区滚石三维模拟器 v3.0  物理自检 ===\n"));
    tlog(_T("模型：半隐式（辛）欧拉积分，Δt = %.4f s，重力 g = %.2f m/s²\n\n"), DT, G);

    // ---- [1] 单一均质斜面纯滚动加速度 ----
    {
        double theta = 30.0 * M_PI / 180.0;
        Terrain saved = g_T;
        g_T.kx.clear(); g_T.kz.clear(); g_T.km.clear();
        g_T.kx.push_back(-20000.0); g_T.kz.push_back(20000.0 * tan(theta));
        g_T.km.push_back(-tan(theta));
        g_T.kx.push_back(20000.0); g_T.kz.push_back(-20000.0 * tan(theta));
        g_T.km.push_back(-tan(theta));
        double mu = 0.9, crr = 0.02;
        Body b;
        b.sh = makeShape(SH_SPHERE, 1.0, 1.0, g_rng);
        b.mass = 1000.0; b.inertia = 0.4 * b.mass;
        b.pos = placeOnTerrain(0.0, 1.0);
        b.ang = 0; b.omega = 0;
        b.vel = V2(0.0, 0.0);
        b.y = 0; b.mu = mu; b.e = 0.3; b.crr = crr; b.cd = 0.0;
        b.alive = true; b.reached = false; b.tReach = 0; b.vToe = 0;
        b.trailAcc = 0; b.bounces = 0; b.maxClear = 0; b.airFrames = 0;
        b.energy0 = b.energy();
        double x0 = b.pos.x;               // 起点水平坐标（含法向偏移）
        int steps = (int)(2.0 / DT);
        for (int i = 0; i < steps; i++) stepBody(b, DT);
        double ds = b.pos.x - x0;
        double path = ds / cos(theta);
        double v = len(b.vel);
        double aNum = (v * v) / (2 * path);
        double aTh = (5.0 / 7.0) * G * sin(theta) - crr * G * cos(theta);
        tlog(_T("[1] 30° 均质斜面，球体纯滚动\n"));
        tlog(_T("    理论 a = (5/7)g·sinθ − c_rr·g·cosθ = %.4f m/s²\n"), aTh);
        tlog(_T("    模拟 a = v²/(2s)                      = %.4f m/s²\n"), aNum);
        tlog(_T("    相对误差 %.2f %%   （判据 < 2%%）\n\n"), 100 * fabs(aNum - aTh) / aTh);
        g_T = saved;
    }

    // ---- [2] 抛体运动 ----
    {
        Body b;
        b.sh = makeShape(SH_SPHERE, 0.6, 1.0, g_rng);
        b.mass = 500; b.inertia = 0.4 * 500 * 0.36;
        b.pos = V2(0, 5000.0); b.vel = V2(12.0, 15.0);
        b.ang = 0; b.omega = 0; b.y = 0;
        b.mu = 0.5; b.e = 0.3; b.crr = 0; b.cd = 0;
        b.alive = true; b.reached = false; b.tReach = 0; b.vToe = 0;
        b.trailAcc = 0; b.bounces = 0; b.maxClear = 0; b.airFrames = 0;
        b.energy0 = b.energy();
        double z0 = b.pos.z, vz0 = b.vel.z;
        double zMax = z0, tApex = 0;
        int steps = (int)(3.0 / DT);
        for (int i = 0; i < steps; i++) {
            stepBody(b, DT);
            if (b.pos.z > zMax) { zMax = b.pos.z; tApex = (i + 1) * DT; }
        }
        double tTh = vz0 / G, hTh = vz0 * vz0 / (2 * G);
        tlog(_T("[2] 斜抛运动（无空气阻力）\n"));
        tlog(_T("    理论 顶点高度 %.3f m，上升时间 %.4f s\n"), hTh, tTh);
        tlog(_T("    模拟 顶点高度 %.3f m，上升时间 %.4f s\n"), zMax - z0, tApex);
        tlog(_T("    高度误差 %.3f %%，时间误差 %.3f %%\n\n"),
            100 * fabs((zMax - z0) - hTh) / hTh, 100 * fabs(tApex - tTh) / tTh);
    }

    // ---- [3] 无耗散能量守恒（均质斜面，积分器自洽性） ----
    {
        double theta = 30.0 * M_PI / 180.0;
        Terrain saved = g_T;
        g_T.kx.clear(); g_T.kz.clear(); g_T.km.clear();
        g_T.kx.push_back(-20000.0); g_T.kz.push_back(20000.0 * tan(theta));
        g_T.km.push_back(-tan(theta));
        g_T.kx.push_back(20000.0); g_T.kz.push_back(-20000.0 * tan(theta));
        g_T.km.push_back(-tan(theta));
        Body b;
        b.sh = makeShape(SH_SPHERE, 1.0, 1.0, g_rng);
        b.mass = 1000; b.inertia = 400;
        b.pos = placeOnTerrain(0.0, 1.0);
        b.vel = V2(0.0, 0.0); b.ang = 0; b.omega = 0;
        b.y = 0; b.mu = 0.9; b.e = 0.4; b.crr = 0.0; b.cd = 0.0;
        b.alive = true; b.reached = false; b.tReach = 0; b.vToe = 0;
        b.trailAcc = 0; b.bounces = 0; b.maxClear = 0; b.airFrames = 0;
        b.energy0 = b.energy();
        double E0 = b.energy();
        g_nContact = 0; g_nStep = 0;
        int steps = (int)(10.0 / DT);
        for (int i = 0; i < steps; i++) stepBody(b, DT);
        double E1 = b.energy();
        double aTh = (5.0 / 7.0) * G * sin(theta);
        double vTh = aTh * 10.0;
        tlog(_T("[3] 无耗散（c_rr = 0, c_d = 0）30° 斜面滚动 10 s 的能量守恒\n"));
        tlog(_T("    初始 E = %.0f J，终止 E = %.0f J，漂移 %.4f %%（判据 < 1%%）\n"),
            E0, E1, 100 * (E1 - E0) / E0);
        tlog(_T("    末速度 模拟 %.3f m/s，理论 (5/7)g·sinθ·t = %.3f m/s，偏差 %.3f %%\n"),
            len(b.vel), vTh, 100 * fabs(len(b.vel) - vTh) / vTh);
        tlog(_T("    诊断：接触步数占比 %.2f %%，末位置 x=%.2f m, z=%.2f m, |ω|R=%.3f m/s, v_t=%.3f m/s\n\n"),
            100.0 * g_nContact / (g_nStep > 0 ? g_nStep : 1),
            b.pos.x, b.pos.z, fabs(b.omega) * b.radius(), len(b.vel));
        g_T = saved;
    }

    // ---- [3b] 真实地形上的坡脚速度（含折角耗散，信息性指标） ----
    {
        Body b;
        b.sh = makeShape(SH_SPHERE, 1.0, 1.0, g_rng);
        b.mass = 1000; b.inertia = 400;
        b.pos = placeOnTerrain(2.0, 1.0);
        V2 tg = nrm(V2(1.0, g_T.dz(2.0)));
        b.vel = tg * 1.0; b.ang = 0; b.omega = -1.0;
        b.y = 0; b.mu = 0.95; b.e = 0.5; b.crr = 0.0; b.cd = 0.0;
        b.alive = true; b.reached = false; b.tReach = 0; b.vToe = 0;
        b.trailAcc = 0; b.bounces = 0; b.maxClear = 0; b.airFrames = 0;
        b.energy0 = b.energy();
        double vToe = 0;
        int steps = (int)(20.0 / DT);
        for (int i = 0; i < steps; i++) {
            stepBody(b, DT);
            if (b.pos.x >= 112.0 && vToe == 0) vToe = len(b.vel);
        }
        double vIdeal = sqrt(10.0 * G * (g_T.z(2.0) - g_T.z(112.0)) / 7.0 + 1.0);
        tlog(_T("[3b] 真实地形、零耗散参数下的坡脚速度\n"));
        tlog(_T("    模拟 %.2f m/s，理想 √(10gh/7 + v0²) = %.2f m/s，比值 %.1f %%\n"),
            vToe, vIdeal, 100 * vToe / vIdeal);
        tlog(_T("    差额来自缓台棱线处的腾空—落地（真实地形几何使然，非数值误差）\n\n"));
    }

    // ---- [4] 多边形静止稳定性 ----
    {
        Body b;
        b.sh = makeShape(SH_SQUARE, 1.0, 1.0, g_rng);
        b.mass = 1000; b.inertia = b.mass * b.sh.I_per_m;
        b.pos = placeOnTerrain(52.0, b.sh.rEff + 0.02);
        b.ang = 0.12; b.omega = 0; b.vel = V2(0, 0);
        b.y = 0; b.mu = 0.6; b.e = 0.3; b.crr = 0.05; b.cd = 0;
        b.alive = true; b.reached = false; b.tReach = 0; b.vToe = 0;
        b.trailAcc = 0; b.bounces = 0; b.maxClear = 0; b.airFrames = 0;
        b.energy0 = b.energy();
        int steps = (int)(4.0 / DT);
        for (int i = 0; i < steps; i++) stepBody(b, DT);
        tlog(_T("[4] 缓台上正方体（初始倾斜 6.9°）放置 4 s 后的稳定性\n"));
        tlog(_T("    |v| = %.5f m/s，|ω| = %.5f rad/s（应趋近 0，判据 < 0.02）\n\n"),
            len(b.vel), fabs(b.omega));
    }

    // ---- [5] 时间步长收敛性 ----
    {
        double res[2] = { 0, 0 };
        for (int pass = 0; pass < 2; pass++) {
            double dt = pass == 0 ? 0.002 : 0.0005;
            Body b;
            b.sh = makeShape(SH_SPHERE, 1.1, 1.0, g_rng);
            b.mass = 1000; b.inertia = 484;
            b.pos = placeOnTerrain(1.0, 1.1);
            V2 tg = nrm(V2(1.0, g_T.dz(1.0)));
            b.vel = tg * 0.8; b.ang = 0; b.omega = -0.8 / 1.1;
            b.y = 0; b.mu = 0.55; b.e = 0.35; b.crr = 0.045; b.cd = 0.55;
            b.alive = true; b.reached = false; b.tReach = 0; b.vToe = 0;
            b.trailAcc = 0; b.bounces = 0; b.maxClear = 0; b.airFrames = 0;
            b.energy0 = b.energy();
            int steps = (int)(30.0 / dt);
            for (int i = 0; i < steps; i++) {
                stepBody(b, dt);
                if (b.pos.x >= 112.0) break;
            }
            res[pass] = len(b.vel);
        }
        tlog(_T("[5] 时间步长收敛性\n"));
        tlog(_T("    Δt = 2.0 ms  坡脚速度 %.3f m/s\n"), res[0]);
        tlog(_T("    Δt = 0.5 ms  坡脚速度 %.3f m/s\n"), res[1]);
        tlog(_T("    相对差 %.2f %%（步长缩小 4 倍，判据 < 1%%）\n\n"),
            100 * fabs(res[0] - res[1]) / res[1]);
    }

    // ---- [6] 各形状同条件对比 ----
    {
        tlog(_T("[6] 十一种形状同条件对比（同一起点、同一材质参数）\n"));
        tlog(_T("    %-10s %8s %8s %9s %10s\n"),
            _T("形状"), _T("坡脚速度"), _T("到坡脚"), _T("停止位置"), _T("弹跳次数"));
        for (int s2 = 0; s2 < SH_COUNT; s2++) {
            Body b;
            b.sh = makeShape(s2, 1.0, 1.0, g_rng);
            double rho = 2650.0;
            if (b.sh.isSphere) { b.mass = rho * 4.0 / 3.0 * M_PI; b.inertia = 0.4 * b.mass; }
            else { b.mass = rho * b.sh.area * b.sh.thickness; b.inertia = b.mass * b.sh.I_per_m; }
            b.pos = placeOnTerrain(1.0, b.sh.rEff + 0.02);
            V2 tg = nrm(V2(1.0, g_T.dz(1.0)));
            b.vel = tg * 0.8; b.ang = 0; b.omega = 0;
            b.y = 0; b.mu = 0.55; b.e = 0.35; b.crr = 0.045; b.cd = 0.55;
            b.alive = true; b.reached = false; b.tReach = 0; b.vToe = 0;
            b.trailAcc = 0; b.bounces = 0; b.maxClear = 0; b.airFrames = 0;
            b.energy0 = b.energy();
            int steps = (int)(60.0 / DT);
            double tReach = -1, vToe = 0;
            for (int i = 0; i < steps; i++) {
                stepBody(b, DT);
                if (b.pos.x >= 112.0 && tReach < 0) { tReach = i * DT; vToe = len(b.vel); }
                if (len(b.vel) < 0.05 && b.pos.x > 20 && i * DT > 3.0) break;
            }
            if (tReach < 0)
                tlog(_T("    %-10s %8s %8s %9.1f %10d\n"), SHAPE_NAME[s2],
                    _T("未到达"), _T("—"), b.pos.x, b.bounces);
            else
                tlog(_T("    %-10s %7.2f %7.2f %9.1f %10d\n"), SHAPE_NAME[s2],
                    vToe, tReach, b.pos.x, b.bounces);
        }
        tlog(_T("\n    说明：坡度 45°/40° 段后接 10° 缓台，tan10° = 0.176 < μ = 0.55，\n"));
        tlog(_T("    棱角块体在缓台上被摩擦锁死而停住，球体靠滚动阻力（c_rr = 0.045）\n"));
        tlog(_T("    仍可继续滚动——这正是形状决定性作用的物理来源。\n"));
    }
    fclose(g_sf);
    g_sf = NULL;
}

// ==================================================================== 主程序
static bool g_paused = false;

static void restart() {
    initWorld();
}

static void renderScene();

static void runLoop() {
    initgraph(WIN_W, WIN_H);
    setbkcolor(C_SKY_HOR);
    cleardevice();
    SetWindowText(GetHWnd(), _T("岩境 / ROCK DYNAMICS LAB  v3.0"));
    BeginBatchDraw();

    g_T.build();
    g_sunDir = nrm3(V3(0.55, 0.32, -0.77));
    g_cam.yaw = -0.581;
    g_cam.pitch = 0.65;
    g_cam.dist = 380.0;
    g_cam.target = V3(85, 0, 22);
    buildProps();
    initWorld();

    LARGE_INTEGER freq, last;
    QueryPerformanceFrequency(&freq);
    QueryPerformanceCounter(&last);
    double accFrame = 0;

    while (true) {
        choiceInput();
        // ---- 输入
        if (GetAsyncKeyState(VK_ESCAPE) & 0x8000) break;
        static bool prevSpace = false;
        bool spaceDown = (GetAsyncKeyState(VK_SPACE) & 0x8000) != 0;
        if (spaceDown && !prevSpace)g_world.paused = !g_world.paused;
        prevSpace = spaceDown;
        if (GetAsyncKeyState('R') & 0x8000) { restart(); Sleep(180); }
        if (GetAsyncKeyState('F') & 0x8000) { g_world.follow = !g_world.follow; Sleep(180); }
        if (GetAsyncKeyState('T') & 0x8000) { g_world.trails = !g_world.trails; Sleep(180); }
        if (GetAsyncKeyState('C') & 0x8000) {
            g_world.bodies.clear(); for (int i = 0;i < SH_COUNT;i++)g_shapeEnabled[i] = false;
            Sleep(180);
        }
        if (GetAsyncKeyState('N') & 0x8000) {
            static int addIdx = 0;
            int id = addIdx % SH_COUNT;
            addIdx++; g_shapeEnabled[id] = true;
            Rng r(20250101ULL + addIdx);
            g_world.bodies.push_back(makeBody(id, 0.95, ((addIdx % 5) - 2) * 1.6, addIdx + 3, r));
            Sleep(180);
        }
        if (GetAsyncKeyState(VK_OEM_PLUS) & 0x8000) { g_world.speed *= 1.02; if (g_world.speed > 6) g_world.speed = 6; }
        if (GetAsyncKeyState(VK_OEM_MINUS) & 0x8000) { g_world.speed *= 0.98; if (g_world.speed < 0.08) g_world.speed = 0.08; }
        if (GetAsyncKeyState('A') & 0x8000) g_cam.yaw -= 0.012;
        if (GetAsyncKeyState('D') & 0x8000) g_cam.yaw += 0.012;
        if (GetAsyncKeyState('W') & 0x8000) g_cam.pitch += 0.008;
        if (GetAsyncKeyState('S') & 0x8000) g_cam.pitch -= 0.008;
        if (GetAsyncKeyState('Q') & 0x8000) g_cam.dist *= 0.988;
        if (GetAsyncKeyState('E') & 0x8000) g_cam.dist *= 1.012;
        if (g_cam.pitch > 1.35) g_cam.pitch = 1.35;
        if (g_cam.pitch < -0.10) g_cam.pitch = -0.10;
        if (g_cam.dist < 12) g_cam.dist = 12;
        if (g_cam.dist > 400) g_cam.dist = 400;

        // ---- 时间
        LARGE_INTEGER now;
        QueryPerformanceCounter(&now);
        double dtReal = (double)(now.QuadPart - last.QuadPart) / freq.QuadPart;
        last = now;
        if (dtReal > 0.1) dtReal = 0.1;
        g_world.frames++;
        accFrame += dtReal;
        if (accFrame > 0.35) { g_world.fps = g_world.frames / accFrame; accFrame = 0; g_world.frames = 0; }

        if (!g_world.paused) {
            double simDt = dtReal * g_world.speed;
            int nsteps = (int)(simDt / DT);
            if (nsteps > 80) nsteps = 80;
            for (int i = 0; i < nsteps; i++) {
                worldStep(DT);
                g_world.t += DT;
            }
        }

        // ---- 相机跟随（跟随领先块体，兼顾整体）
        if (g_world.follow && !g_world.bodies.empty()) {
            double bx = -1e9, bz = 0;
            double sumx = 0; int cnt = 0;
            for (size_t i = 0; i < g_world.bodies.size(); i++) {
                const Body& b = g_world.bodies[i];
                if (b.pos.x > bx) { bx = b.pos.x; bz = b.pos.z; }
                sumx += b.pos.x; cnt++;
            }
            double lead = bx;
            double tgt = std::min(lead, sumx / std::max(1, cnt) + 45.0);
            g_cam.target.x += (tgt - g_cam.target.x) * 0.045;
            g_cam.target.z += ((g_T.z(g_cam.target.x) + 12.0) - g_cam.target.z) * 0.05;
            (void)bz;
        }

        // ---- 渲染
        renderScene();
        FlushBatchDraw();
        Sleep(1);
        continue;
    }
    EndBatchDraw();
    closegraph();
}

// ---------------------------------------------------------------- 渲染一帧
static void renderScene() {
    updateCamera();
    g_faces.clear();
    g_faces.reserve(2600);
    drawSky();
    buildTerrain(g_cam.target.x);
    buildSection();
    drawProps(g_cam.target.x);
    // 轨迹（贴着地面的细条）
    if (g_world.trails) {
        for (size_t i = 0; i < g_world.bodies.size(); i++) {
            const Body& b = g_world.bodies[i];
            for (size_t k = 1; k < b.trail.size(); k++) {
                V2 p0 = b.trail[k - 1], p1 = b.trail[k];
                if (p0.x < g_cam.target.x - 40 || p0.x > g_cam.target.x + 150) continue;
                double z0 = g_T.z(p0.x) + 0.06, z1 = g_T.z(p1.x) + 0.06;
                double w = 0.16;
                V3 q[4] = {
                    V3(p0.x, b.y - w, z0), V3(p1.x, b.y - w, z1),
                    V3(p1.x, b.y + w, z1), V3(p0.x, b.y + w, z0) };
                addFace(q, 4, mixc(b.col, RGB(255, 255, 255), 0.12), V3(0, 0, 1), 0.30);
            }
        }
    }
    for (size_t i = 0; i < g_world.bodies.size(); i++)
        drawBodyFaces(g_world.bodies[i], (int)i);
    paintFaces();
    drawWeather();
    drawHUD();
    drawChoices();
}

// ------------------------------------------------------------ 离屏截图模式
//  用于生成报告插图：跑到指定时刻，渲染一帧，存成 PNG，不做交互
static void shotMode(double secs, const TCHAR* out) {
    initgraph(WIN_W, WIN_H);
    setbkcolor(C_SKY_HOR);
    cleardevice();
    SetWindowText(GetHWnd(), _T("岩境 / ROCK DYNAMICS LAB  v3.0"));
    g_T.build();
    g_sunDir = nrm3(V3(0.55, 0.32, -0.77));
    g_cam.yaw = -0.581; g_cam.pitch = 0.65; g_cam.dist = 380.0;
    g_cam.target = V3(85, 0, 22);
    buildProps();
    initWorld();
    int steps = (int)(secs / DT);
    for (int i = 0; i < steps; i++) { worldStep(DT); g_world.t += DT; }
    // 相机对准领先块体
    double bx = -1e9, sumx = 0; int cnt = 0;
    for (size_t i = 0; i < g_world.bodies.size(); i++) {
        const Body& b = g_world.bodies[i];
        if (b.pos.x > bx) bx = b.pos.x;
        sumx += b.pos.x; cnt++;
    }
    g_cam.target.x = std::min(bx, sumx / std::max(1, cnt) + 45.0);
    g_cam.target.z = g_T.z(g_cam.target.x) + 12.0;
    renderScene();
    saveimage(out);
    closegraph();
}

static bool cmdArg(const char* flag) {
    const char* c = GetCommandLineA();
    return c && strstr(c, flag) != NULL;
}

int main(int argc, char** argv) {
    (void)argc; (void)argv;
    if (cmdArg("--test")) {
        g_T.build();
        selftest();
        return 0;
    }
    if (cmdArg("--shot")) {
        double secs = 5.0;
        const char* c = GetCommandLineA();
        const char* q = c ? strstr(c, "--shot=") : NULL;
        if (q) secs = atof(q + 7);
        shotMode(secs, _T("preview.png"));
        return 0;
    }
    runLoop();
    return 0;
}
