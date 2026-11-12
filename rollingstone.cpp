#include <graphics.h>
#include <conio.h>
#include <math.h>
#include <cstring>
#include <cmath>
#include <tchar.h>
#include <windows.h>
#include <cstdlib>
#include <cstdio>

// ============================================================================
//  山区滚石模拟（米制版）
//  ---------------------------------------------------------------
//  物理量全部用国际单位制：位置 m、速度 m/s、时间 s、质量 kg、能量 J
//  只有在“画图”这一步才换算成像素，换算比例 PX_PER_M
//  物理模型：
//    · 触地滚动：a = (5/7) g sinθ − c_rr g cosθ − 空气阻力      （实心球纯滚动）
//    · 陡坡打滑：tanθ > μ 时摩擦不足以维持纯滚动，改用 a = g(sinθ − μcosθ)
//    · 空中飞行：竖直重力 + 二次空气阻力 ½ρCdAv²
//    · 碰撞弹跳：法向恢复系数 e，切向按库仑冲量损失
//    · 停止：速度足够小且坡度小于摩擦角（不规则块体被地面抓住）
// ============================================================================

// ------------------------------------------------------------------ 画面
#define WIDTH 1000          // 窗口宽度（像素）
#define HEIGHT 600          // 窗口高度（像素）
#define PX_PER_M 8.0        // 显示比例：1 米 = 8 像素（窗口即 125 m × 75 m）
#define ORIGIN_Y 490.0      // 高程 z=0 对应的屏幕 y 坐标
#define GROUND_Y 490        // 坡脚基准线（屏幕 y）

// ------------------------------------------------------------------ 物理常数
#define GRAVITY 9.81        // 重力加速度 (m/s²)
#define DT 0.002            // 物理步长 (s)
#define STEPS_PER_FRAME 8   // 每帧 8 步 × 2 ms = 16 ms ≈ 实时
#define ROCK_DENSITY 2600.0 // 岩石密度 (kg/m³)
#define AIR_DENSITY 1.225   // 空气密度 (kg/m³)
#define MIN_RADIUS 0.2      // 半径可调下限 (m)，只用来防止误操作，不是默认半径
#define MAX_RADIUS 3.0      // 半径可调上限 (m)
#define BOUNCE_MIN 0.15     // 法向速度小于此值不再弹跳，按接触处理 (m/s)
#define REST_SPEED 0.05     // 判定“停住”的速度阈值 (m/s)

// 可调参数（程序运行时通过菜单修改）
double g_rollResist = 0.02;      // 滚动阻力系数 c_rr（0.02~0.20）
double g_slideFriction = 0.45;   // 地面摩擦系数 μ（tan 24°，坡面摩擦角约 24°）
double g_restitution = 0.45;     // 法向恢复系数 e
double g_dragCoef = 1.0;         // 空气阻力系数 Cd（0 表示真空）
int g_playbackFactor = 1;        // 1x,2x,4x,8x

// ------------------------------------------------------------------ 坐标换算
double toPxX(double x) { return x * PX_PER_M; }
double toPxY(double z) { return ORIGIN_Y - z * PX_PER_M; }

// 地形点（x 水平距离 m，z 高程 m，向上为正）
struct Point {
    double x, z;
    Point(double x = 0, double z = 0) : x(x), z(z) {}
};

// 时间-速度数据结构（替代 std::pair<double,double>）
struct TimeSpeed {
    double t;
    double s;
    TimeSpeed(double tt = 0.0, double ss = 0.0) : t(tt), s(ss) {}
};

// 简单动态数组模板（替代 std::vector）
template<typename T>
class DynamicArray {
private:
    T* data;
    size_t cnt;
    size_t cap;

    void grow() {
        size_t newCap = cap == 0 ? 4 : cap * 2;
        T* newData = new T[newCap];
        if (data) {
            for (size_t i = 0; i < cnt; ++i) newData[i] = data[i];
            delete[] data;
        }
        data = newData;
        cap = newCap;
    }

public:
    DynamicArray() : data(nullptr), cnt(0), cap(0) {}
    ~DynamicArray() { delete[] data; }
    DynamicArray(const DynamicArray& other) : data(nullptr), cnt(0), cap(0) {
        if (other.cnt > 0) {
            data = new T[other.cnt];
            for (size_t i = 0; i < other.cnt; ++i) data[i] = other.data[i];
            cnt = other.cnt; cap = other.cnt;
        }
    }
    DynamicArray& operator=(const DynamicArray& other) {
        if (this == &other) return *this;
        delete[] data;
        data = nullptr; cnt = 0; cap = 0;
        if (other.cnt > 0) {
            data = new T[other.cnt];
            for (size_t i = 0; i < other.cnt; ++i) data[i] = other.data[i];
            cnt = other.cnt; cap = other.cnt;
        }
        return *this;
    }

    void push_back(const T& v) {
        if (cnt >= cap) grow();
        data[cnt++] = v;
    }

    void erase_front() {
        if (cnt == 0) return;
        for (size_t i = 1; i < cnt; ++i) data[i - 1] = data[i];
        --cnt;
    }

    size_t size() const { return cnt; }
    bool empty() const { return cnt == 0; }
    T& operator[](size_t i) { return data[i]; }
    const T& operator[](size_t i) const { return data[i]; }
    void clear() { cnt = 0; }
};

// ------------------------------------------------------------------ 山坡地形
// 13 个控制点，x 从 0 到 120 m，高程从 55.6 m 降到 0 m（坡脚）
class HillSlope {
private:
    DynamicArray<Point> points;  // 地形控制点（米）
    int startX, endX;            // 只用于快速判断，单位仍是米（取整）

public:
    HillSlope() {
        // 一条长 120 m、落差 55.6 m 的山区坡面：上部陡、中部有缓台、坡脚有一段平台
        points.push_back(Point(0.0, 55.6));    // 坡顶（危岩位置）
        points.push_back(Point(10.0, 51.9));
        points.push_back(Point(11.25, 50.6));
        points.push_back(Point(30.0, 35.6));   // 陡坡段
        points.push_back(Point(40.0, 34.4));   // 缓台（约 7°）
        points.push_back(Point(50.0, 23.1));   // 陡坡段
        points.push_back(Point(60.0, 20.6));
        points.push_back(Point(70.0, 15.6));
        points.push_back(Point(80.0, 10.0));
        points.push_back(Point(90.0, 8.75));   // 缓台
        points.push_back(Point(100.0, 2.88));
        points.push_back(Point(120.0, 0.0));   // 坡脚平台（道路位置）

        startX = (int)points[0].x;
        endX = (int)points[points.size() - 1].x;
    }

    // 指定 x 处的地面高程 z (m)
    double getHeight(double x) const {
        if (x <= points[0].x) return points[0].z;
        if (x >= points[points.size() - 1].x) return points[points.size() - 1].z;

        for (size_t i = 0; i < points.size() - 1; i++) {
            if (x >= points[i].x && x <= points[i + 1].x) {
                double t = (x - points[i].x) / (points[i + 1].x - points[i].x);
                return points[i].z + t * (points[i + 1].z - points[i].z);
            }
        }
        return points[0].z;
    }

    // 指定 x 处的坡角 θ（弧度，>0 表示顺坡向下）
    double getSlopeAngle(double x) const {
        if (x <= points[0].x) {
            double dz = points[1].z - points[0].z;
            double dx = points[1].x - points[0].x;
            return -atan2(dz, dx);
        }
        if (x >= points[points.size() - 1].x) {
            double dz = points[points.size() - 1].z - points[points.size() - 2].z;
            double dx = points[points.size() - 1].x - points[points.size() - 2].x;
            return -atan2(dz, dx);
        }

        for (size_t i = 0; i < points.size() - 1; i++) {
            if (x >= points[i].x && x <= points[i + 1].x) {
                double dz = points[i + 1].z - points[i].z;
                double dx = points[i + 1].x - points[i].x;
                return -atan2(dz, dx);
            }
        }
        return 0;
    }

    // 绘制山坡（唯一需要把米换成像素的地方）
    void draw() const {
        for (size_t i = 0; i < points.size() - 1; i++) {
            int x1 = (int)toPxX(points[i].x);
            int y1 = (int)toPxY(points[i].z);
            int x2 = (int)toPxX(points[i + 1].x);
            int y2 = (int)toPxY(points[i + 1].z);

            setcolor(RGB(34, 139, 34));
            setlinestyle(PS_SOLID, 3);
            line(x1, y1, x2, y2);

            setfillcolor(RGB(139, 90, 43));
            int pts[8] = { x1, y1, x2, y2, x2, HEIGHT, x1, HEIGHT };
            fillpoly(4, pts);
        }

        // 坡脚基准线（z = 0）
        setcolor(RGB(160, 82, 45));
        setlinestyle(PS_SOLID, 1);
        line(0, (int)toPxY(0), WIDTH, (int)toPxY(0));
    }

    int getStartX() const { return startX; }
    int getEndX() const { return endX; }

    int getPeakX() const {
        if (points.empty()) return startX;
        double bestZ = points[0].z;
        int bestX = (int)points[0].x;
        for (size_t i = 1; i < points.size(); ++i) {
            if (points[i].z > bestZ) {   // 高程最大者即坡顶
                bestZ = points[i].z;
                bestX = (int)points[i].x;
            }
        }
        return bestX;
    }

    int getValleyX() const {
        if (points.empty()) return startX;
        double lowZ = points[0].z;
        int lowX = (int)points[0].x;
        for (size_t i = 1; i < points.size(); ++i) {
            if (points[i].z < lowZ) {    // 高程最小者即坡脚
                lowZ = points[i].z;
                lowX = (int)points[i].x;
            }
        }
        return lowX;
    }
};

// ------------------------------------------------------------------ 滚石
class RollingStone {
private:
    double x, z;        // 球心位置（m）
    double vx, vz;      // 速度（m/s）
    double radius;      // 半径（m）
    bool isRolling;     // 是否还在运动
    double energy;      // 动能（J）
    double omega;       // 角速度（rad/s），纯滚动时 ω = v/r
    double mass;        // 质量（kg），由半径和密度算出
    DynamicArray<Point> trail;
    int trailMax;

    void updateMass() {
        double vol = 4.0 / 3.0 * 3.14159265358979 * radius * radius * radius;
        mass = ROCK_DENSITY * vol;
    }

public:
    // 半径没有默认值：必须由调用者（main）给定，对象只负责保存和使用
    RollingStone(double startX, double startZ, double r, double initVx)
        : x(startX), z(startZ), vx(initVx), vz(0), radius(r), isRolling(true),
          energy(0), omega(0), mass(0), trailMax(40) {
        if (radius < MIN_RADIUS) radius = MIN_RADIUS;
        if (radius > MAX_RADIUS) radius = MAX_RADIUS;
        updateMass();
        // 让球心正好贴在坡面上：沿法向抬高一个半径（真正的坡角在 adjustToGround 里再修正）
        z = startZ + radius;
    }

    // 大气阻力系数：a_drag = kd·v²，kd = ½ρCdA/m
    double dragK() const {
        double A = 3.14159265358979 * radius * radius;
        return 0.5 * AIR_DENSITY * g_dragCoef * A / mass;
    }

    // 物理更新：一步 dt
    void update(const HillSlope& hill, double dt) {
        if (!isRolling) return;

        double zg = hill.getHeight(x);
        double th = hill.getSlopeAngle(x);        // 坡角，>0 顺坡向下
        double ct = cos(th), st = sin(th);
        double nx = st, nz = ct;                  // 坡面外法向
        double dx = ct, dz = -st;                 // 顺坡切向（下坡方向）

        // 球心到坡面的法向距离：竖直高差投影到法向
        double dist = (z - zg) * ct;
        double vt = vx * dx + vz * dz;            // 切向速度（下坡为正）
        double vn = vx * nx + vz * nz;            // 法向速度（离开坡面为正）

        if (dist <= radius + 1e-3) {
            // ---------------- 接触坡面 ----------------
            z = zg + radius / ct;                 // 沿法向推回，恰好在表面

            if (vn < -BOUNCE_MIN) {
                // 弹跳：法向按恢复系数反弹，切向按库仑冲量损失
                double vnOut = -g_restitution * vn;
                double dvt = g_slideFriction * (1.0 + g_restitution) * fabs(vn);
                double sgn = (vt >= 0) ? 1.0 : -1.0;
                double vtOut = fabs(vt) - dvt;
                if (vtOut < 0) vtOut = 0;
                vt = sgn * vtOut;
                vn = vnOut;
                omega = vt / radius;
            }
            else {
                // 贴地运动：纯滚动 or 打滑
                vn = 0;
                double a;
                if (tan(th) > g_slideFriction) {
                    // 摩擦不足以维持纯滚动 → 打滑，转动不再同步
                    a = GRAVITY * (st - g_slideFriction * ct);
                    omega = 0;
                }
                else {
                    // 实心球纯滚动：转动惯量 I = 2/5 mR²，只有 5/7 的重力分量用于平动
                    a = (5.0 / 7.0) * GRAVITY * st - g_rollResist * GRAVITY * ct;
                    omega = vt / radius;
                }
                vt += a * dt;
                // 空气阻力（沿切向，与速度反向）
                double kd = dragK();
                vt -= kd * fabs(vt) * vt * dt;

                // 速度很小且坡度小于摩擦角 → 被地面抓住，停住
                if (fabs(vt) < REST_SPEED && tan(th) <= g_slideFriction) {
                    vt = 0;
                    isRolling = false;
                }
                else if (vt < 0) {
                    vt = 0;
                }
            }
            vx = vt * dx + vn * nx;
            vz = vt * dz + vn * nz;
        }
        else {
            // ---------------- 空中飞行 ----------------
            vz -= GRAVITY * dt;
            double sp = sqrt(vx * vx + vz * vz);
            double kd = dragK();
            if (sp > 1e-9) {
                vx -= kd * sp * vx * dt;
                vz -= kd * sp * vz * dt;
            }
        }

        // 位置更新
        x += vx * dt;
        z += vz * dt;

        // 边界：坡顶以上、坡脚以外
        if (x < hill.getStartX()) {
            x = hill.getStartX();
            vx = 0;
        }
        if (x >= hill.getEndX()) {
            x = hill.getEndX();
            vx = 0;
            isRolling = false;
        }

        // 防止穿过坡面：位置夹回表面，并去掉指向坡内的法向速度
        zg = hill.getHeight(x);
        th = hill.getSlopeAngle(x);
        {
            double ct2 = cos(th), st2 = sin(th);
            if ((z - zg) * ct2 < radius) {
                z = zg + radius / ct2;
                double vn2 = vx * st2 + vz * ct2;
                if (vn2 < 0) {
                    vx -= vn2 * st2;
                    vz -= vn2 * ct2;
                }
            }
        }

        // 动能：滚动时含转动动能 E = ½mv² + ½Iω² = 0.7mv²；其他情况只有平动
        double v2 = vx * vx + vz * vz;
        if (fabs(omega) > 1e-9) energy = 0.5 * mass * v2 + 0.2 * mass * radius * radius * omega * omega;
        else energy = 0.5 * mass * v2;

        // 记录轨迹（世界坐标，米）
        trail.push_back(Point(x, z));
        if ((int)trail.size() > trailMax) trail.erase_front();
    }

    // 绘制石头与轨迹（米 → 像素在这里完成）
    void draw() const {
        // 轨迹
        for (size_t i = 0; i < trail.size(); i++) {
            int alpha = 100 - (int)(i * 100 / trail.size());
            (void)alpha;
            setcolor(RGB(255, 100, 100));
            setfillstyle(BS_SOLID);
            int size = 3 - (int)(i * 2 / trail.size());
            if (size < 1) size = 1;
            solidcircle((int)toPxX(trail[i].x), (int)toPxY(trail[i].z), size);
        }

        // 抗锯齿圆掩码（按屏幕半径缓存）
        struct CircleMask {
            int r;
            int sz;
            unsigned char* alpha;
            CircleMask() :r(0), sz(0), alpha(nullptr) {}
            ~CircleMask() { delete[] alpha; }
            void generate(int radiusPx, int ss = 4) {
                if (radiusPx <= 0) return;
                if (r == radiusPx && alpha) return;
                delete[] alpha; alpha = nullptr;
                r = radiusPx;
                sz = 2 * r + 1;
                alpha = new unsigned char[sz * sz];
                int samples = ss * ss;
                double rr = (double)r;
                for (int iy = -r; iy <= r; ++iy) {
                    for (int ix = -r; ix <= r; ++ix) {
                        int inside = 0;
                        for (int sy = 0; sy < ss; ++sy) {
                            for (int sx = 0; sx < ss; ++sx) {
                                double fx = ix + (sx + 0.5) / ss;
                                double fy = iy + (sy + 0.5) / ss;
                                if (fx * fx + fy * fy <= rr * rr) inside++;
                            }
                        }
                        unsigned char a = (unsigned char)((inside * 255) / samples);
                        alpha[(iy + r) * sz + (ix + r)] = a;
                    }
                }
            }
        };

        int ir = (int)(radius * PX_PER_M);
        if (ir < 4) ir = 4;              // 显示时给一个最小半径，否则 1 m 的石头只有 8 px

        static CircleMask mask;
        mask.generate(ir, 3);

        COLORREF stoneColor = RGB(120, 120, 120);
        int mr = mask.r; int msz = mask.sz;

        // 阴影
        for (int yi = 0; yi < msz; ++yi) {
            for (int xi = 0; xi < msz; ++xi) {
                unsigned char a = mask.alpha[yi * msz + xi];
                if (a == 0) continue;
                int px = (int)toPxX(x) - mr + xi + 2;
                int py = (int)toPxY(z) - mr + yi + 3;
                if (px < 0 || px >= WIDTH || py < 0 || py >= HEIGHT) continue;
                COLORREF bg = getpixel(px, py);
                COLORREF sc = RGB(60, 60, 60);
                int alpha = a;
                int inv = 255 - alpha;
                int br = (GetRValue(sc) * alpha + GetRValue(bg) * inv) / 255;
                int bgc = (GetGValue(sc) * alpha + GetGValue(bg) * inv) / 255;
                int bb = (GetBValue(sc) * alpha + GetBValue(bg) * inv) / 255;
                putpixel(px, py, RGB(br, bgc, bb));
            }
        }

        // 主体
        for (int yi = 0; yi < msz; ++yi) {
            for (int xi = 0; xi < msz; ++xi) {
                unsigned char a = mask.alpha[yi * msz + xi];
                if (a == 0) continue;
                int px = (int)toPxX(x) - mr + xi;
                int py = (int)toPxY(z) - mr + yi;
                if (px < 0 || px >= WIDTH || py < 0 || py >= HEIGHT) continue;
                COLORREF bg = getpixel(px, py);
                int alpha = a;
                int inv = 255 - alpha;
                int rr = (GetRValue(stoneColor) * alpha + GetRValue(bg) * inv) / 255;
                int gg = (GetGValue(stoneColor) * alpha + GetGValue(bg) * inv) / 255;
                int bb = (GetBValue(stoneColor) * alpha + GetBValue(bg) * inv) / 255;
                putpixel(px, py, RGB(rr, gg, bb));
            }
        }

        // 高光
        int hr = max(1, (int)(ir / 3));
        static CircleMask hmask;
        hmask.generate(hr, 3);
        COLORREF highlight = RGB(180, 180, 180);
        int hmr = hmask.r; int hmsz = hmask.sz;
        int hx = (int)(toPxX(x) - ir / 3);
        int hy = (int)(toPxY(z) - ir / 3);
        for (int yi = 0; yi < hmsz; ++yi) {
            for (int xi = 0; xi < hmsz; ++xi) {
                unsigned char a = hmask.alpha[yi * hmsz + xi];
                if (a == 0) continue;
                int px = hx - hmr + xi;
                int py = hy - hmr + yi;
                if (px < 0 || px >= WIDTH || py < 0 || py >= HEIGHT) continue;
                COLORREF bg = getpixel(px, py);
                int alpha = a;
                int inv = 255 - alpha;
                int rr = (GetRValue(highlight) * alpha + GetRValue(bg) * inv) / 255;
                int gg = (GetGValue(highlight) * alpha + GetGValue(bg) * inv) / 255;
                int bb = (GetBValue(highlight) * alpha + GetBValue(bg) * inv) / 255;
                putpixel(px, py, RGB(rr, gg, bb));
            }
        }
    }

    bool isMoving() const { return isRolling; }
    double getX() const { return x; }
    double getZ() const { return z; }
    double getSpeed() const { return sqrt(vx * vx + vz * vz); }
    double getEnergy() const { return energy; }
    double getRadius() const { return radius; }
    double getMass() const { return mass; }
    void setRadius(double r) {
        if (r < MIN_RADIUS) r = MIN_RADIUS;
        if (r > MAX_RADIUS) r = MAX_RADIUS;
        radius = r;
        updateMass();
    }
    void adjustToGround(const HillSlope& hill) {
        double zg = hill.getHeight(x);
        z = zg + radius / cos(hill.getSlopeAngle(x));
        vz = 0;
    }
    // 顺坡释放：初速度沿坡面切向，而不是水平（更接近“石块从坡面脱离”）
    void startAlongSlope(const HillSlope& hill, double v0) {
        double th = hill.getSlopeAngle(x);
        vx = v0 * cos(th);
        vz = -v0 * sin(th);
        omega = v0 / radius;
    }
};

// ------------------------------------------------------------------ 主程序
int main() {
    initgraph(WIDTH, HEIGHT);
    setbkcolor(RGB(135, 206, 235));
    BeginBatchDraw();

    HillSlope hill;
    int initX = hill.getPeakX();
    double initZ = hill.getHeight((double)initX);

    // ---- 滚石参数：在 main 里给定，构造时交给滚石对象（对象内部只保存、不再有默认值）----
    double stoneRadius = 2.0;   // 滚石半径 (m)
    double stoneSpeed = 1.5;    // 顺坡释放的初速度 (m/s)

    RollingStone stone((double)initX, initZ, stoneRadius, 0.0);
    stone.adjustToGround(hill);
    stone.startAlongSlope(hill, stoneSpeed);
    int valleyX = hill.getValleyX();
    bool simEnded = false;

    bool paused = false;
    double stoppedDuration = 0.0;
    bool stoppedTriggered = false;
    bool inMenu = false;
    int menuSel = 0;                       // 0: c_rr 1: μ 2: e 3: Cd 4: 倍速
    const int menuRows = 4;
    const int playbackOptions[4] = { 1, 2, 4, 8 };
    int playbackIndex = 0;
    bool inPauseMenu = false;
    int pauseMenuSel = 0;

    char inputBuf[4][32];
    sprintf_s(inputBuf[0], sizeof(inputBuf[0]), "%.3f", g_rollResist);
    sprintf_s(inputBuf[1], sizeof(inputBuf[1]), "%.3f", g_slideFriction);
    sprintf_s(inputBuf[2], sizeof(inputBuf[2]), "%.3f", g_restitution);
    sprintf_s(inputBuf[3], sizeof(inputBuf[3]), "%.3f", g_dragCoef);

    DynamicArray< TimeSpeed > speedData;
    double simTime = 0.0;
    const double graphWindow = 10.0;

    settextstyle(16, 0, _T("Consolas"));
    setcolor(RGB(50, 50, 50));
    outtextxy(10, 20, _T("Rolling Stone Simulator - SI Units (m, s, kg, J)"));
    outtextxy(10, 45, _T("Controls: SPACE=Pause, M=Parameters, P=Pause menu, R=Reset, ESC=Exit"));

    while (true) {
        if (_kbhit()) {
            int key = _getch();
            if (key == 27) break;
            if (key == 32) {
                paused = !paused;
            }
            if (key == 'p' || key == 'P') {
                inPauseMenu = !inPauseMenu;
                pauseMenuSel = paused ? 1 : 0;
            }
            if (inPauseMenu) {
                if (key == 'w' || key == 'W' || key == 's' || key == 'S') {
                    pauseMenuSel = 1 - pauseMenuSel;
                }
                if (key == 13) {
                    if (pauseMenuSel == 0) paused = true; else paused = false;
                    inPauseMenu = false;
                }
            }
            if (key == 'r' || key == 'R') {
                double keepRadius = stone.getRadius();   // 保留用户用 +/- 改过的半径
                stone = RollingStone((double)initX, initZ, keepRadius, 0.0);
                stone.adjustToGround(hill);
                stone.startAlongSlope(hill, stoneSpeed);
                paused = false;
                simEnded = false;
                stoppedDuration = 0.0;
                stoppedTriggered = false;
                speedData.clear();
                simTime = 0.0;
                cleardevice();
            }
            if (key == '+' || key == '=') {
                stone.setRadius(stone.getRadius() + 0.1);
                stone.adjustToGround(hill);
            }
            if (key == '-' || key == '_') {
                stone.setRadius(stone.getRadius() - 0.1);
                stone.adjustToGround(hill);
            }
            if (key == 'm' || key == 'M') {
                inMenu = !inMenu;
                playbackIndex = 0;
                for (int i = 0; i < 4; i++) if (playbackOptions[i] == g_playbackFactor) playbackIndex = i;
            }
            if (inMenu) {
                if (key == 'w' || key == 'W') { menuSel = (menuSel + menuRows) % (menuRows + 1); }
                if (key == 's' || key == 'S') { menuSel = (menuSel + 1) % (menuRows + 1); }
                if (key == 'a' || key == 'A') {
                    if (menuSel == 0) g_rollResist = max(0.0, g_rollResist - 0.01);
                    if (menuSel == 1) g_slideFriction = max(0.0, g_slideFriction - 0.05);
                    if (menuSel == 2) g_restitution = max(0.0, g_restitution - 0.05);
                    if (menuSel == 3) g_dragCoef = max(0.0, g_dragCoef - 0.1);
                    if (menuSel == 4) { playbackIndex = max(0, playbackIndex - 1); g_playbackFactor = playbackOptions[playbackIndex]; }
                }
                if (key == 'd' || key == 'D') {
                    if (menuSel == 0) g_rollResist = min(0.30, g_rollResist + 0.01);
                    if (menuSel == 1) g_slideFriction = min(1.00, g_slideFriction + 0.05);
                    if (menuSel == 2) g_restitution = min(0.90, g_restitution + 0.05);
                    if (menuSel == 3) g_dragCoef = min(2.00, g_dragCoef + 0.1);
                    if (menuSel == 4) { playbackIndex = min(3, playbackIndex + 1); g_playbackFactor = playbackOptions[playbackIndex]; }
                }
                if (key == 13) { inMenu = false; }
            }
        }

        if (!paused) {
            // 每帧 8 个子步（2 ms × 8 = 16 ms），倍速就是“每帧多跑几步”
            int steps = STEPS_PER_FRAME * g_playbackFactor;
            for (int i = 0; i < steps; ++i) {
                stone.update(hill, DT);
                simTime += DT;
            }
            double curSpeed = stone.getSpeed();
            speedData.push_back(TimeSpeed(simTime, curSpeed));
            while (!speedData.empty() && speedData[0].t < simTime - graphWindow) {
                speedData.erase_front();
            }

            if (!simEnded) {
                if (stone.getX() >= hill.getEndX() - 0.5) {
                    MessageBox(NULL, _T("滚石到达坡脚，模拟结束"), _T("提示"), MB_OK | MB_ICONINFORMATION);
                    EndBatchDraw();
                    closegraph();
                    return 0;
                }
            }

            if (!stoppedTriggered && !simEnded) {
                if (curSpeed < REST_SPEED) stoppedDuration += STEPS_PER_FRAME * DT * g_playbackFactor;
                else stoppedDuration = 0.0;
                if (stoppedDuration >= 3.0) {
                    MessageBox(NULL, _T("滚石停止运动，模拟结束"), _T("提示"), MB_OK | MB_ICONINFORMATION);
                    EndBatchDraw();
                    closegraph();
                    return 0;
                }
            }

            cleardevice();
            setbkcolor(RGB(135, 206, 235));

            if (inMenu) {
                setfillcolor(RGB(240, 240, 240));
                bar(70, 45, 500, 440);
                setcolor(RGB(0, 0, 0));
                settextstyle(16, 0, _T("Consolas"));
                outtextxy(90, 55, _T("Parameters (W/S select, A/D change, Enter apply)"));

                TCHAR bufMenu[160];
                for (int i = 0; i < menuRows; i++) {
                    int y = 100 + i * 55;
                    if (i == menuSel) setcolor(RGB(255, 0, 0)); else setcolor(RGB(0, 0, 0));
                    if (i == 0) _stprintf_s(bufMenu, 160, _T("滚动阻力 c_rr (0.02-0.20):"));
                    else if (i == 1) _stprintf_s(bufMenu, 160, _T("摩擦系数 mu (0.3-0.8):"));
                    else if (i == 2) _stprintf_s(bufMenu, 160, _T("恢复系数 e (0.2-0.6):"));
                    else _stprintf_s(bufMenu, 160, _T("空气阻力 Cd (0=真空):"));
                    outtextxy(90, y, bufMenu);

                    setcolor(RGB(0, 0, 0));
                    TCHAR tbuf[64];
                    size_t outLen = 0;
                    mbstowcs_s(&outLen, tbuf, sizeof(tbuf) / sizeof(TCHAR), inputBuf[i], sizeof(inputBuf[i]));
                    outtextxy(400, y, tbuf);
                }
                int y = 100 + menuRows * 55;
                if (menuSel == menuRows) setcolor(RGB(255, 0, 0)); else setcolor(RGB(0, 0, 0));
                _stprintf_s(bufMenu, 160, _T("播放倍速: x1 x2 x4 x8 (当前 x%d)"), g_playbackFactor);
                outtextxy(90, y, bufMenu);

                setfillcolor(RGB(100, 200, 100));
                bar(250, 380, 380, 415);
                setcolor(RGB(0, 0, 0));
                outtextxy(275, 388, _T("启动模拟"));

                if (GetAsyncKeyState(VK_RETURN) & 0x8000) {
                    double v0 = atof(inputBuf[0]);
                    double v1 = atof(inputBuf[1]);
                    double v2 = atof(inputBuf[2]);
                    double v3 = atof(inputBuf[3]);
                    if (v0 < 0.0) v0 = 0.0; if (v0 > 0.30) v0 = 0.30;
                    if (v1 < 0.0) v1 = 0.0; if (v1 > 1.00) v1 = 1.00;
                    if (v2 < 0.0) v2 = 0.0; if (v2 > 0.90) v2 = 0.90;
                    if (v3 < 0.0) v3 = 0.0; if (v3 > 2.00) v3 = 2.00;
                    g_rollResist = v0; g_slideFriction = v1;
                    g_restitution = v2; g_dragCoef = v3;
                    inMenu = false; paused = false; simEnded = false;
                    stoppedDuration = 0.0; stoppedTriggered = false;
                    cleardevice();
                }
                if (menuSel >= 0 && menuSel <= 3) {
                    for (int c = 0; c < 256; c++) {
                        if (GetAsyncKeyState(c) & 0x8000) {
                            if (c == VK_BACK) {
                                size_t len = strlen(inputBuf[menuSel]);
                                if (len > 0) inputBuf[menuSel][len - 1] = 0;
                            }
                            else if ((c >= '0' && c <= '9') || c == '.' || c == '-') {
                                size_t len = strlen(inputBuf[menuSel]);
                                if (len + 1 < sizeof(inputBuf[menuSel])) {
                                    inputBuf[menuSel][len] = (char)c;
                                    inputBuf[menuSel][len + 1] = 0;
                                }
                            }
                        }
                    }
                }

                setcolor(RGB(80, 80, 80));
                outtextxy(90, 360, _T("按 R 重置模拟，按 M 退出菜单"));
            }

            if (inPauseMenu) {
                setfillcolor(RGB(240, 240, 240));
                bar(WIDTH / 2 - 120, HEIGHT / 2 - 40, WIDTH / 2 + 120, HEIGHT / 2 + 40);
                setcolor(RGB(0, 0, 0));
                settextstyle(16, 0, _T("Consolas"));
                outtextxy(WIDTH / 2 - 60, HEIGHT / 2 - 24, _T("Pause Menu (W/S select, Enter)"));
                if (pauseMenuSel == 0) setcolor(RGB(255, 0, 0)); else setcolor(RGB(0, 0, 0));
                outtextxy(WIDTH / 2 - 40, HEIGHT / 2 - 2, _T("Pause"));
                if (pauseMenuSel == 1) setcolor(RGB(255, 0, 0)); else setcolor(RGB(0, 0, 0));
                outtextxy(WIDTH / 2 - 40, HEIGHT / 2 + 18, _T("Continue"));
            }

            hill.draw();
            stone.draw();

            // -------- HUD（米制） --------
            TCHAR info[200];
            setcolor(RGB(255, 255, 255));
            settextstyle(14, 0, _T("Consolas"));

            _stprintf_s(info, 200, _T("Position: %.1f m"), stone.getX());
            outtextxy(10, HEIGHT - 100, info);

            _stprintf_s(info, 200, _T("Height: %.1f m"), stone.getZ());
            outtextxy(10, HEIGHT - 80, info);

            _stprintf_s(info, 200, _T("Speed: %.1f m/s"), stone.getSpeed());
            outtextxy(10, HEIGHT - 60, info);

            _stprintf_s(info, 200, _T("Energy: %.1f kJ"), stone.getEnergy() / 1000.0);
            outtextxy(10, HEIGHT - 40, info);

            _stprintf_s(info, 200, _T("Diameter: %.2f m   Mass: %.0f kg"), stone.getRadius() * 2.0, stone.getMass());
            outtextxy(10, HEIGHT - 20, info);

            _stprintf_s(info, 200, _T("Status: %s"), stone.isMoving() ? _T("ROLLING") : _T("STOPPED"));
            setcolor(stone.isMoving() ? RGB(255, 200, 0) : RGB(100, 255, 100));
            outtextxy(10, HEIGHT - 120, info);

            double slope = hill.getSlopeAngle(stone.getX()) * 180 / 3.14159;
            setcolor(RGB(200, 200, 200));
            _stprintf_s(info, 200, _T("Slope: %.1f deg"), slope);
            outtextxy(WIDTH - 130, HEIGHT - 60, info);

            TCHAR paramBuf[240];
            _stprintf_s(paramBuf, 240, _T("c_rr:%.2f mu:%.2f e:%.2f Cd:%.1f x%d"),
                        g_rollResist, g_slideFriction, g_restitution, g_dragCoef, g_playbackFactor);
            setcolor(RGB(180, 180, 180));
            settextstyle(12, 0, _T("Consolas"));
            outtextxy(WIDTH - 330, 10, paramBuf);

            // 速度条（0~20 m/s）
            int speedBar = (int)(stone.getSpeed() / 20.0 * 100);
            if (speedBar > 100) speedBar = 100;
            setfillcolor(RGB(255, 100, 100));
            bar(WIDTH - 120, HEIGHT - 40, WIDTH - 120 + speedBar, HEIGHT - 30);
            setcolor(RGB(150, 150, 150));
            rectangle(WIDTH - 120, HEIGHT - 40, WIDTH - 20, HEIGHT - 30);

            // -------- 速度-时间曲线 --------
            const int gx = WIDTH - 260;
            const int gy = 20;
            const int gw = 240;
            const int gh = 120;

            setfillcolor(RGB(135, 206, 235));
            bar(gx, gy, gx + gw, gy + gh);
            setcolor(RGB(200, 200, 200));
            rectangle(gx, gy, gx + gw, gy + gh);

            const int marginL = 36;
            const int marginR = 8;
            const int marginTop = 6;
            const int marginBottom = 20;
            int plotX = gx + marginL;
            int plotY = gy + marginTop;
            int plotW = gw - marginL - marginR;
            int plotH = gh - marginTop - marginBottom;

            double t0 = simTime - graphWindow;
            if (t0 < 0) t0 = 0;

            setcolor(RGB(180, 180, 180));
            line(plotX, plotY, plotX, plotY + plotH);
            line(plotX, plotY + plotH, plotX + plotW, plotY + plotH);

            double gridStart = floor(t0);
            double gridEnd = ceil(simTime);
            setcolor(RGB(220, 220, 220));
            setlinestyle(PS_DASH, 1);
            for (int g = (int)gridStart; g <= (int)gridEnd; ++g) {
                double gxTime = (double)g;
                double gxp = plotX + (gxTime - (simTime - graphWindow)) / graphWindow * plotW;
                if (gxp >= plotX && gxp <= plotX + plotW) {
                    line((int)gxp, plotY, (int)gxp, plotY + plotH);
                }
            }
            setlinestyle(PS_SOLID, 1);

            double yMax = 8.0;   // 纵轴最小量程 8 m/s
            for (size_t _i = 0; _i < speedData.size(); ++_i) {
                auto& p = speedData[_i];
                if (p.s > yMax) yMax = p.s;
            }
            yMax *= 1.1;

            settextstyle(12, 0, _T("Consolas"));
            setcolor(RGB(80, 80, 80));
            TCHAR buf[64];
            _stprintf_s(buf, 64, _T("%.1f"), yMax);
            outtextxy(plotX - 32, plotY - 8, buf);
            _stprintf_s(buf, 64, _T("%.1f"), 0.0);
            outtextxy(plotX - 28, plotY + plotH - 14, buf);

            _stprintf_s(buf, 64, _T("%.1fs"), t0);
            outtextxy(plotX + 2, plotY + plotH + 2, buf);
            _stprintf_s(buf, 64, _T("%.1fs"), simTime);
            outtextxy(plotX + plotW - 40, plotY + plotH + 2, buf);

            _stprintf_s(buf, 64, _T("Speed (m/s)"));
            outtextxy(plotX + 2, plotY - 18, buf);
            _stprintf_s(buf, 64, _T("Time (s)"));
            outtextxy(plotX + plotW / 2 - 22, plotY + plotH + 2, buf);

            if (speedData.size() >= 2) {
                setcolor(RGB(255, 0, 0));
                setlinestyle(PS_SOLID, 2);
                for (size_t i = 1; i < speedData.size(); ++i) {
                    double tx1 = speedData[i - 1].t;
                    double ty1 = speedData[i - 1].s;
                    double tx2 = speedData[i].t;
                    double ty2 = speedData[i].s;

                    double nx1 = plotX + (tx1 - (simTime - graphWindow)) / graphWindow * plotW;
                    double nx2 = plotX + (tx2 - (simTime - graphWindow)) / graphWindow * plotW;
                    double ny1 = plotY + plotH - (ty1 / yMax) * plotH;
                    double ny2 = plotY + plotH - (ty2 / yMax) * plotH;

                    if (nx1 < plotX) nx1 = plotX;
                    if (nx1 > plotX + plotW) nx1 = plotX + plotW;
                    if (nx2 < plotX) nx2 = plotX;
                    if (nx2 > plotX + plotW) nx2 = plotX + plotW;
                    if (ny1 < plotY) ny1 = plotY;
                    if (ny1 > plotY + plotH) ny1 = plotY + plotH;
                    if (ny2 < plotY) ny2 = plotY;
                    if (ny2 > plotY + plotH) ny2 = plotY + plotH;

                    line((int)nx1, (int)ny1, (int)nx2, (int)ny2);
                }
                setlinestyle(PS_SOLID, 1);
            }

            FlushBatchDraw();
            Sleep(16);
        }
        else {
            FlushBatchDraw();
            Sleep(50);
        }

        if (!stone.isMoving() && !paused) {
            setcolor(RGB(0, 255, 0));
            settextstyle(16, 0, _T("Consolas"));
            outtextxy(WIDTH / 2 - 80, HEIGHT - 150, _T("Stone Stopped! Press R to restart"));
        }
    }

    closegraph();
    return 0;
}
