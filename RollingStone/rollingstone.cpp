// Visual Studio 2026 + EasyX: C++17, Console subsystem, Unicode, UTF-8.
#define NOMINMAX
#define WIN32_LEAN_AND_MEAN
#include <graphics.h>
#include <windows.h>
#include <windowsx.h>
#include <commdlg.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>
#pragma comment(lib, "user32.lib")
#pragma comment(lib, "gdi32.lib")
#pragma comment(lib, "comdlg32.lib")
#pragma comment(lib, "shell32.lib")

constexpr double DT = 1.0 / 240.0;
struct Vec
{
    double x = 0, y = 0, z = 0;
    Vec operator+(Vec b) const
    {
        return { x + b.x, y + b.y, z + b.z };
    }
    Vec operator-(Vec b) const
    {
        return { x - b.x, y - b.y, z - b.z };
    }
    Vec operator*(double s) const
    {
        return { x * s, y * s, z * s };
    }
    double dot(Vec b) const
    {
        return x * b.x + y * b.y + z * b.z;
    }
    Vec cross(Vec b) const
    {
        return { y * b.z - z * b.y, z * b.x - x * b.z, x * b.y - y * b.x };
    }
    double length() const
    {
        return std::sqrt(dot(*this));
    }
    Vec unit() const
    {
        return *this * (1.0 / std::max(1e-12, length()));
    }
};
struct Quaternion
{
    double a = 1, x = 0, y = 0, z = 0;
    Vec rotate(Vec v) const
    {
        Vec u{ x, y, z }, t = u.cross(v) * 2;
        return v + t * a + u.cross(t);
    }
    Vec inverse(Vec v) const
    {
        return Quaternion{ a, -x, -y, -z }.rotate(v);
    }
    void integrate(Vec w, double dt)
    {
        Quaternion q = *this;
        a += (-q.x * w.x - q.y * w.y - q.z * w.z) * dt / 2;
        x += (q.a * w.x + q.y * w.z - q.z * w.y) * dt / 2;
        y += (q.a * w.y + q.z * w.x - q.x * w.z) * dt / 2;
        z += (q.a * w.z + q.x * w.y - q.y * w.x) * dt / 2;
        double s = std::sqrt(a * a + x * x + y * y + z * z);
        a /= s;
        x /= s;
        y /= s;
        z /= s;
    }
};
double terrain(double x, double y)
{
    double t = std::clamp(1 - x / 90, 0.0, 1.0);
    return 42 * t + 2 * std::sin(x * .19) * t + .45 * std::sin(y * .24 + x * .11) * t;
}
Vec terrainNormal(double x, double y)
{
    double e = .02;
    return Vec{ -(terrain(x + e, y) - terrain(x - e, y)) / (2 * e),
               -(terrain(x, y + e) - terrain(x, y - e)) / (2 * e), 1 }
    .unit();
}
struct Body
{
    Vec p, v, w, alpha, half{ 1, .75, .85 };
    Quaternion q;
    double mass = 13260, rest = 0;
    int id = 0;
    std::vector<Vec> trail;
    std::array<Vec, 8> vertices() const
    {
        std::array<Vec, 8> out{};
        int i = 0;
        for (int x : {-1, 1})
            for (int y : {-1, 1})
                for (int z : {-1, 1})
                    out[i++] = q.rotate({ x * half.x, y * half.y, z * half.z });
        return out;
    }
    Vec inertia(Vec v, bool inverse) const
    {
        Vec u = q.inverse(v);
        Vec d{ mass * (half.y * half.y + half.z * half.z) / 3,
              mass * (half.x * half.x + half.z * half.z) / 3,
              mass * (half.x * half.x + half.y * half.y) / 3 };
        return q.rotate(inverse ? Vec{ u.x / d.x, u.y / d.y, u.z / d.z }
        : Vec{ u.x * d.x, u.y * d.y, u.z * d.z });
    }
    double energy() const
    {
        return .5 * mass * v.dot(v) + .5 * w.dot(inertia(w, false));
    }
    void impulse(Vec j, Vec r)
    {
        v = v + j * (1 / mass);
        w = w + inertia(r.cross(j), true);
    }
    double effective(Vec r, Vec n) const
    {
        return 1 / mass + n.dot(inertia(r.cross(n), true).cross(r));
    }
    void contact(Vec n, Vec r, double depth, double mu, double e)
    {
        Vec u = v + w.cross(r);
        double vn = u.dot(n);
        if (vn < 0)
        {
            double j = -(1 + (vn < -.8 ? e : 0)) * vn / effective(r, n);
            impulse(n * j, r);
            u = v + w.cross(r);
            Vec t = u - n * u.dot(n);
            double speed = t.length();
            if (speed > 1e-9)
            {
                t = t * (1 / speed);
                impulse(t * (-std::min(mu * j, speed / effective(r, t))), r);
            }
        }
        if (depth > .001)
            p = p + n * ((depth - .001) * .35);
    }
};
struct Row
{
    double time;
    int id;
    Vec p, v, w, a;
    double energy;
};
std::vector<Body> bodies;
std::vector<Row> rows;
bool rain = false, snow = false, wind = false, barrier = true, fracture = false, broken = false,
paused = false;
double initialSpeed = 8, simTime = 0, accumulator = 0, yaw = -.9, pitch = .6, zoom = 8;
int stepCount = 0, rate = 1;
HWND mainWindow;
POINT mouse{};
bool dragging = false;
void reset()
{
    Body b;
    b.p = { 2, 0, terrain(2, 0) + 2 };
    b.v = { initialSpeed, 0, 0 };
    bodies = { b };
    rows.clear();
    broken = false;
    simTime = 0;
    stepCount = 0;
    accumulator = 0;
    paused = false;
}
void step()
{
    simTime += DT;
    double mu = snow ? .12 : rain ? .28 : .55, e = snow ? .08 : .25;
    for (Body& b : bodies)
    {
        Vec old = b.w;
        Vec rel = b.v - (wind ? Vec{ 0, 18, 0 } : Vec{});
        double area = 4 * b.half.x * b.half.y;
        Vec force = Vec{ 0, 0, -9.81 * b.mass } + rel * (-.5 * 1.225 * 1.05 * area * rel.length());
        // Euler's rigid-body equation in world coordinates (full box inertia).
        Vec torque = b.w.cross(b.inertia(b.w, false)) * (-1);
        b.w = b.w + b.inertia(torque, true) * DT;
        b.v = b.v + force * (DT / b.mass);
        b.p = b.p + b.v * DT;
        b.q.integrate(b.w, DT);
        for (int pass = 0; pass < 10; pass++)
            for (Vec r : b.vertices())
            {
                Vec p = b.p + r, n = terrainNormal(p.x, p.y);
                double depth = (terrain(p.x, p.y) - p.z) * n.z;
                if (depth > -.002)
                    b.contact(n, r, depth, mu, e);
                if (barrier && p.x > 98 && p.x < 102 && p.z < 7 && std::abs(p.y) < 28)
                    b.contact({ -1, 0, 0 }, r, p.x - 98, mu, .1);
            }
        b.alpha = (b.w - old) * (1 / DT);
        if (b.v.length() < .12 && b.w.length() < .12)
            b.rest += DT;
        else
            b.rest = 0;
        if (stepCount % 12 == 0)
        {
            b.trail.push_back(b.p);
            if (b.trail.size() > 2000)
                b.trail.erase(b.trail.begin());
            rows.push_back({ simTime, b.id, b.p, b.v, b.w, b.alpha, b.energy() });
        }
    }
    // Optional geometric split of a box across its local y axis. Children exactly
    // partition mass/inertia; inherited rigid velocity preserves P, L and energy.
    if (fracture && !broken && bodies[0].p.x >= 48 &&
        bodies[0].p.z < terrain(48, bodies[0].p.y) + 5)
    {
        Body parent = bodies[0], a = parent, b = parent;
        Vec offset = parent.q.rotate({ 0, parent.half.y / 2, 0 });
        a.mass = b.mass = parent.mass / 2;
        a.half.y = b.half.y = parent.half.y / 2;
        a.p = parent.p - offset;
        b.p = parent.p + offset;
        a.v = parent.v + parent.w.cross(offset * (-1));
        b.v = parent.v + parent.w.cross(offset);
        a.id = 1;
        b.id = 2;
        a.trail.clear();
        b.trail.clear();
        bodies = { a, b };
        broken = true;
    }
    stepCount++;
}
POINT project(Vec p, int width, int height)
{
    double x = p.x - 50, y = p.y, z = p.z - 18;
    double u = x * std::cos(yaw) - y * std::sin(yaw), v = x * std::sin(yaw) + y * std::cos(yaw);
    return { LONG(width * .57 + u * zoom),
            LONG(height * .55 + (v * std::sin(pitch) - z * std::cos(pitch)) * zoom) };
}
double depth(Vec p)
{
    return ((p.x - 50) * std::sin(yaw) + p.y * std::cos(yaw)) * std::cos(pitch) +
        (p.z - 18) * std::sin(pitch);
}
struct Face
{
    std::vector<Vec> points;
    COLORREF color;
    double d = 0;
};
void text(int x, int y, const std::wstring& s)
{
    outtextxy(x, y, s.c_str());
}
void paint(int width, int height)
{
    setbkcolor(RGB(24, 44, 59));
    cleardevice();
    std::vector<Face> fs;
    for (int x = 0; x < 114; x += 3)
        for (int y = -30; y < 30; y += 3)
        {
            double l = std::max(.2, terrainNormal(x, y).dot(Vec{ -.3, -.5, 1 }.unit()));
            int c = int(l * 70);
            fs.push_back(
                { { { double(x), double(y), terrain(x, y) },
                  { double(x + 3), double(y), terrain(x + 3, y) },
                  { double(x + 3), double(y + 3), terrain(x + 3, y + 3) },
                  { double(x), double(y + 3), terrain(x, y + 3) } },
                 snow ? RGB(150 + c, 160 + c, 170 + c) : RGB(45 + c / 2, 65 + c, 50 + c / 2) });
        }
    if (barrier)
        fs.push_back({ { { 98, -28, 0 },{ 98, 28, 0 },{ 98, 28, 7 },{ 98, -28, 7 } }, RGB(210, 160, 90) });
    if (fracture)
        fs.push_back(
            { { { 46, -2, terrain(46, -2) },{ 50, 2, terrain(50, 2) },{ 48, 0, terrain(48, 0) + 3 } },
             RGB(200, 155, 145) });
    const int indices[6][4] = { { 0, 1, 3, 2 },{ 4, 6, 7, 5 },{ 0, 4, 5, 1 },
                               { 2, 3, 7, 6 },{ 0, 2, 6, 4 },{ 1, 5, 7, 3 } };
    for (const Body& b : bodies)
    {
        auto vertices = b.vertices();
        for (int i = 0; i < 6; i++)
        {
            Face f;
            f.color = RGB(140 + i * 15, 95 + i * 12, 75 + i * 10);
            for (int j : indices[i])
                f.points.push_back(b.p + vertices[j]);
            fs.push_back(f);
        }
    }
    for (Face& f : fs)
    {
        for (Vec p : f.points)
            f.d += depth(p);
        f.d /= f.points.size();
    }
    std::sort(fs.begin(), fs.end(), [](const Face& a, const Face& b) { return a.d < b.d; });
    setlinecolor(RGB(70, 90, 85));
    setlinestyle(PS_SOLID, 1);
    for (const Face& f : fs)
    {
        std::vector<POINT> points;
        for (Vec p : f.points)
            points.push_back(project(p, width, height));
        setfillcolor(f.color);
        fillpolygon(points.data(), int(points.size()));
    }
    setlinecolor(RGB(255, 190, 100));
    setlinestyle(PS_SOLID, 2);
    for (auto& b : bodies)
    {
        bool first = true;
        for (Vec p : b.trail)
        {
            POINT a = project(p, width, height);
            if (first)
            {
                moveto(a.x, a.y);
                first = false;
            }
            else
                lineto(a.x, a.y);
        }
    }
    setbkmode(TRANSPARENT);
    settextcolor(RGB(235, 242, 245));
    settextstyle(18, 0, L"Microsoft YaHei");
    text(20, 62, L"三维凸多面体滚石 | 拖动旋转，滚轮缩放");
    wchar_t buf[512];
    const Body& b = bodies[0];
    double energy = 0;
    for (auto& stone : bodies)
        energy += stone.energy();
    swprintf_s(buf, L"时间 %.2f s  初速度 %.1f m/s  倍速 %d×  %s", simTime, initialSpeed, rate,
        paused ? L"已暂停"
        : b.rest > 1 ? L"停稳"
        : L"运动中");
    text(20, height - 150, buf);
    swprintf_s(buf, L"质心位置 (%.2f, %.2f, %.2f) m   速度 %.2f m/s", b.p.x, b.p.y, b.p.z,
        b.v.length());
    text(20, height - 123, buf);
    swprintf_s(buf, L"角速度 %.2f rad/s   角加速度 %.2f rad/s²", b.w.length(), b.alpha.length());
    text(20, height - 96, buf);
    swprintf_s(buf, L"总动能 %.2f kJ   石块数 %zu   %s", energy / 1000, bodies.size(),
        broken ? L"已触发简化破裂" : L"");
    text(20, height - 69, buf);
    text(20, height - 40, L"空格暂停 / R重置 / ↑↓调整释放初速度 / T倍速 / 导出完整三维数据");
    int gx = width - 300, gy = 90, gw = 265, gh = 130;
    setlinecolor(RGB(100, 130, 145));
    rectangle(gx, gy, gx + gw, gy + gh);
    text(gx, gy - 23, L"第一石块速度 / m·s⁻¹，最近20秒");
    double vmax = 5;
    for (const Row& r : rows)
        if (r.time > simTime - 20 && r.id == b.id)
            vmax = std::max(vmax, r.v.length());
    bool first = true;
    setlinecolor(RGB(255, 200, 110));
    setlinestyle(PS_SOLID, 2);
    for (const Row& r : rows)
        if (r.time >= simTime - 20 && r.id == b.id)
        {
            int x = gx + int((r.time - std::max(0.0, simTime - 20)) / 20 * gw),
                y = gy + gh - int(r.v.length() / (vmax * 1.1) * gh);
            if (first)
            {
                moveto(x, y);
                first = false;
            }
            else
                lineto(x, y);
        }
}
void exportCsv()
{
    wchar_t path[MAX_PATH] = L"rollingstone.csv";
    OPENFILENAMEW o{};
    o.lStructSize = sizeof(o);
    o.hwndOwner = mainWindow;
    o.lpstrFilter = L"CSV 文件\0*.csv\0\0";
    o.lpstrFile = path;
    o.nMaxFile = MAX_PATH;
    o.lpstrDefExt = L"csv";
    o.Flags = OFN_OVERWRITEPROMPT | OFN_PATHMUSTEXIST | OFN_NOCHANGEDIR;
    if (!GetSaveFileNameW(&o))
        return;
    std::ofstream out{ std::filesystem::path(path) };
    out << "time_s,id,x_m,y_m,z_m,vx_m_s,vy_m_s,vz_m_s,wx_rad_s,wy_rad_s,wz_rad_s,ax_rad_s2,ay_rad_s2,az_rad_s2,energy_J\n"
        << std::setprecision(10);
    for (const Row& r : rows)
        out << r.time << ',' << r.id << ',' << r.p.x << ',' << r.p.y << ',' << r.p.z << ',' << r.v.x
        << ',' << r.v.y << ',' << r.v.z << ',' << r.w.x << ',' << r.w.y << ',' << r.w.z << ','
        << r.a.x << ',' << r.a.y << ',' << r.a.z << ',' << r.energy << '\n';
    out.close();
    MessageBoxW(mainWindow, out ? L"CSV 已保存" : L"保存失败，请检查路径权限", L"导出", MB_OK);
}

// All controls below are painted and handled through EasyX.
const wchar_t* buttonLabels[] = { L"雨",       L"雪",        L"大风",     L"挡板",
                                 L"简化破裂", L"暂停/继续", L"重置释放", L"导出CSV" };
RECT buttonRect(int i)
{
    return { 12 + i * 112, 12, 116 + i * 112, 44 };
}
void drawButtons()
{
    bool active[] = { rain, snow, wind, barrier, fracture, paused, false, false };
    settextstyle(17, 0, L"Microsoft YaHei");
    setbkmode(TRANSPARENT);
    setlinestyle(PS_SOLID, 1);
    for (int i = 0; i < 8; ++i)
    {
        RECT r = buttonRect(i);
        setfillcolor(active[i] ? RGB(45, 120, 135) : RGB(40, 65, 80));
        setlinecolor(RGB(110, 150, 165));
        fillrectangle(r.left, r.top, r.right, r.bottom);
        settextcolor(WHITE);
        outtextxy(r.left + 10, r.top + 6, buttonLabels[i]);
    }
}
bool clickButton(int x, int y)
{
    for (int i = 0; i < 8; ++i)
    {
        RECT r = buttonRect(i);
        if (x < r.left || x > r.right || y < r.top || y > r.bottom)
            continue;
        switch (i)
        {
        case 0:
            rain = !rain;
            break;
        case 1:
            snow = !snow;
            break;
        case 2:
            wind = !wind;
            break;
        case 3:
            barrier = !barrier;
            break;
        case 4:
            fracture = !fracture;
            break;
        case 5:
            paused = !paused;
            break;
        case 6:
            reset();
            break;
        case 7:
            exportCsv();
            break;
        }
        return true;
    }
    return false;
}
int main()
{
    constexpr int WIDTH = 1280, HEIGHT = 820;
    initgraph(WIDTH, HEIGHT);
    mainWindow = GetHWnd();
    SetWindowTextW(mainWindow, L"山区滚石模拟 · EasyX 三维版本");
    reset();
    BeginBatchDraw();
    LARGE_INTEGER frequency{}, last{};
    QueryPerformanceFrequency(&frequency);
    QueryPerformanceCounter(&last);
    bool running = true;
    while (running && IsWindow(mainWindow))
    {
        ExMessage message;
        while (peekmessage(&message, EX_MOUSE | EX_KEY | EX_WINDOW))
        {
            switch (message.message)
            {
            case WM_CLOSE:
                running = false;
                break;
            case WM_KEYDOWN:
                if (message.prevdown)
                    break;
                if (message.vkcode == VK_ESCAPE)
                    running = false;
                else if (message.vkcode == VK_SPACE)
                    paused = !paused;
                else if (message.vkcode == 'R')
                    reset();
                else if (message.vkcode == VK_UP)
                    initialSpeed = std::min(20.0, initialSpeed + 1);
                else if (message.vkcode == VK_DOWN)
                    initialSpeed = std::max(0.0, initialSpeed - 1);
                else if (message.vkcode == 'T')
                    rate = rate == 4 ? 1 : rate * 2;
                break;
            case WM_LBUTTONDOWN:
                if (!clickButton(message.x, message.y))
                {
                    dragging = true;
                    mouse = { message.x, message.y };
                    SetCapture(mainWindow);
                }
                break;
            case WM_MOUSEMOVE:
                if (dragging && message.lbutton)
                {
                    yaw += (message.x - mouse.x) * .008;
                    pitch = std::clamp(pitch + (message.y - mouse.y) * .006, .12, 1.4);
                    mouse = { message.x, message.y };
                }
                break;
            case WM_LBUTTONUP:
                dragging = false;
                ReleaseCapture();
                break;
            case WM_MOUSEWHEEL:
                zoom = std::clamp(zoom * std::exp(message.wheel * .001), 3.0, 25.0);
                break;
            }
        }
        if (!running || !IsWindow(mainWindow))
            break;
        LARGE_INTEGER now;
        QueryPerformanceCounter(&now);
        double elapsed = std::min(.1, double(now.QuadPart - last.QuadPart) / frequency.QuadPart);
        last = now;
        if (!paused && simTime < 180)
        {
            accumulator += elapsed * rate;
            while (accumulator >= DT && simTime < 180)
            {
                step();
                accumulator -= DT;
            }
            if (simTime >= 180)
                paused = true;
        }
        paint(WIDTH, HEIGHT);
        drawButtons();
        FlushBatchDraw();
        Sleep(8);
    }
    EndBatchDraw();
    closegraph();
    return 0;
}
