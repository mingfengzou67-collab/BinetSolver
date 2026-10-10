%% =========================================================================
%  比耐公式（Binet 公式）求解器 —— 理论力学 · 有心力轨道
%  -------------------------------------------------------------------------
%  【这个脚本干什么】
%     质点在有心力 F(r) 作用下运动。比耐公式把"轨道的形状"和"力"直接连起来：
%
%          d²u/dθ² + u = -F(1/u) / ( m h² u² )       其中 u = 1/r
%
%     于是有两类题（本脚本都做）：
%        用法一：已知力 F(r)  →  解出轨道 r(θ)      …… 例子 1、2、3
%        用法二：已知轨道 r(θ) →  反推出力 F(r)      …… 例子 4
%
%     符号：m 质量；r 到力心的距离；θ 极角；u = 1/r；
%           h = r²·(dθ/dt) 是"单位质量的角动量"，是个常数；
%           F > 0 是斥力（背离力心），F < 0 是引力（指向力心）。
%
%  【怎么用】三步，小白也会
%     第 1 步：在 MATLAB 里打开本文件（把本文件夹拖进 MATLAB 再双击文件）
%     第 2 步：只改下面"参数区"里的 exampleId（1~4），别的地方都不用动
%     第 3 步：按键盘 F5（或点工具栏那个绿色三角"运行"按钮）
%     运行后命令行会列出题目让你按回车选，选完自动弹图。
%
%  【想做自己的题】
%     看文件最后 examples 函数里的"例 1"，照着格式把 F(r) 改掉就行，
%     那里有详细注释。本脚本只用 MATLAB 自带函数，不需要任何工具箱。
% =========================================================================

clc; close all;

%% --------------------------- 参数区（只改这里） ---------------------------
exampleId = 1;      % 选题目：1、2、3、4；               改成 0 = 运行时用键盘选
animate   = true;   % 要不要播放"质点沿轨道跑"的动画：true 要，false 不要
% -------------------------------------------------------------------------


%% ========================= 主程序（下面不用改） =========================

P = examples();                     % 把 4 个现成的题目取出来

if exampleId == 0
    exampleId = chooseExample(P);   % 让使用者用键盘选
end
if ~isscalar(exampleId) || exampleId < 1 || exampleId > numel(P) || exampleId ~= round(exampleId)
    warning('exampleId 只能是 1~%d，已自动改用例子 1。', numel(P));
    exampleId = 1;
end

ex = P{exampleId};
fprintf('\n');
fprintf('==========================================================\n');
fprintf(' 例 %d：%s\n', exampleId, ex.name);
fprintf('==========================================================\n');
fprintf('%s\n\n', ex.note);

if strcmp(ex.mode, 'orbit')
    runOrbitProblem(ex, animate);       % 已知力，求轨道
else
    runForceProblem(ex);                % 已知轨道，求力
end

fprintf('\n提示：换题目就改最上面参数区的 exampleId；\n');
fprintf('      想做自己的题，就改文件末尾 examples() 里例 1 的 Ffun（力）。\n');


%% ===================== 用法一：已知力求轨道 =====================
function runOrbitProblem(ex, animate)
    m = ex.m;
    h = ex.r0 * ex.vt0;             % 单位质量角动量 h = r0 * 横向速度

    fprintf('【初始条件】 m = %g，r0 = %g，径向速度 vr0 = %g，横向速度 vt0 = %g\n', ...
            m, ex.r0, ex.vr0, ex.vt0);
    fprintf('【角动量】   h = r0·vt0 = %g   （比耐公式里 h 必须先用这里算出来）\n', h);
    fprintf('【要解的方程】d²u/dθ² + u = -F(1/u)/(m h² u²)，u = 1/r，未知量 u(θ)\n\n');

    % ---------- 1. 用 ode45 解比耐方程 ----------
    [theta, u, du] = solveBinet(ex.Ffun, m, h, ex.r0, ex.vr0, ex.thetaMax);
    r = 1 ./ u;

    % ---------- 2. 由轨道算出速度、能量（用来检查算得准不准）----------
    vr  = -h * du;                          % 径向速度 dr/dt = -h·(du/dθ)
    vt  =  h * u;                           % 横向速度 r·(dθ/dt) = h·u
    E   = 0.5*m*(vr.^2 + vt.^2) + arrayfun(ex.Vfun, r);
    isBound = mean(E) < 0;                  % 能量为负 = 被绑住 = 闭合轨道

    % ---------- 3. 打印结果 ----------
    fprintf('---- 数值结果 ----\n');
    fprintf('积分到 θ = %.3f rad（约 %.2f 圈）\n', theta(end), theta(end)/(2*pi));
    if isBound
        eNum = (max(r) - min(r)) / (max(r) + min(r));
        fprintf('轨道是闭合轨道（平均能量 E = %.4f < 0）\n', mean(E));
        fprintf('r_min = %.4f（近力心点），r_max = %.4f（远力心点）\n', min(r), max(r));
        fprintf('数值偏心率 e = (r_max - r_min)/(r_max + r_min) = %.4f\n', eNum);
    else
        fprintf('轨道是开放轨道（平均能量 E = %.4f > 0），质点会跑到无穷远\n', mean(E));
        fprintf('r_min = %.4f（最靠近力心的那一点）\n', min(r));
    end
    fprintf('能量相对波动 = %.2e（这么小说明数值解很准）\n', ...
            (max(E) - min(E)) / abs(mean(E)));

    % 平方反比力的题目，可以跟理论公式对答案
    rAnalytic = [];
    if ~isempty(ex.kAn)
        E0     = 0.5*m*(ex.vr0^2 + ex.vt0^2) + ex.Vfun(ex.r0);
        eExact = sqrt(max(0, 1 + 2*E0*m*h^2/ex.kAn^2));   % e = √(1 + 2Emh²/k²)
        pExact = m*h^2/ex.kAn;                            % 半正焦弦 p = mh²/k
        fprintf('\n---- 理论对照（纯平方反比力 F = -k/r²）----\n');
        fprintf('理论偏心率 e = √(1 + 2Emh²/k²) = %.4f  ', eExact);
        if eExact < 1e-8
            fprintf('→ 圆轨道\n');
        elseif eExact < 1
            fprintf('→ 椭圆（0 < e < 1）\n');
        elseif abs(eExact - 1) < 1e-8
            fprintf('→ 抛物线（e = 1）\n');
        else
            fprintf('→ 双曲线（e > 1）\n');
        end
        fprintf('半正焦弦 p = mh²/k = %.4f，理论 r_min = p/(1+e) = %.4f\n', ...
                pExact, pExact/(1+eExact));
        if isBound
            fprintf('数值 e 与理论 e 的差 = %.2e（越接近 0 越好）\n', abs(eNum - eExact));
        end
        if ex.vr0 == 0        % 从最近点/最远点出发时，可以画出理论轨道来对照
            if eExact < 1e-8
                rAnalytic = ex.r0 * ones(size(theta));
            else
                s = (pExact/ex.r0 - 1) / eExact;          % s=+1 从近心点出发，s=-1 从远心点出发
                rAnalytic = pExact ./ (1 + s*eExact*cos(theta));
            end
        end
    end

    % 找近力心点（r 的极小值），进动就是看它每次转到哪儿
    periIdx = find(islocalmin(r));
    fprintf('\n---- 近力心点（r 最小处）的位置 ----\n');
    if numel(periIdx) == 0
        fprintf('本次积分范围内没有出现 r 的极小值（开放轨道，或者还没跑到）\n');
    else
        fprintf('θ = %s rad\n', mat2str(theta(periIdx)', 4));
        dPeri = diff(theta(periIdx));
        if numel(dPeri) >= 1
            fprintf('相邻两次近心点的角度间隔 = %s rad\n', mat2str(dPeri', 4));
            fprintf('每圈进动角 ≈ %s 度（等于 0 就是闭合轨道；不等于 0 就是进动）\n', ...
                    mat2str((dPeri' - 2*pi)*180/pi, 4));
        end
    end

    % ---------- 4. 画图 ----------
    plotOrbit(ex, theta, r, E, rAnalytic, periIdx, animate);
end


%% ===================== 用法二：已知轨道求力 =====================
function runForceProblem(ex)
    m = ex.m;  h = ex.h;
    fprintf('【已知轨道】 r(θ) = %.4g / (1 + %.4g·cosθ)   （圆锥曲线，p = %g，e = %g）\n', ...
            ex.p, ex.ecc, ex.p, ex.ecc);
    fprintf('【角动量】   h = %g（求力这类题里，h 只影响力的整体大小）\n', h);
    fprintf('【思路】     比耐公式反着用：F = -m h² u² (d²u/dθ² + u)，u = 1/r\n\n');

    theta = linspace(0, 2*pi, 4001)';        % 把一圈均匀切成 4000 份
    r     = arrayfun(ex.rfun, theta);
    u     = 1 ./ r;

    % 数值求一阶、二阶导数
    du  = gradient(u, theta);
    d2u = gradient(du, theta);

    % 代入比耐公式，反推出力
    Fnum  = -m * h^2 .* u.^2 .* (d2u + u);
    Ftrue = -(m*h^2/ex.p) ./ r.^2;           % 这种轨道的标准答案：F 正比于 -1/r²

    % 差分的两端不准，比较时去掉两头
    keep = 101:3901;
    err  = abs(Fnum(keep) - Ftrue(keep)) ./ abs(Ftrue(keep));

    fprintf('---- 反推结果 ----\n');
    fprintf('如果反推对了，F·r² 应该是常数：\n');
    fprintf('    数值给出  F·r² = %.6f\n', mean(Fnum(keep).*r(keep).^2));
    fprintf('    理论应为 -mh²/p = %.6f\n', -m*h^2/ex.p);
    fprintf('F(r) 的最大相对误差 = %.3e（很小，说明反推成功）\n', max(err));

    % ---------- 画图 ----------
    newFigure(['比耐公式 - ' ex.name], [80 80 1080 430]);

    subplot(1,2,1);
    loglog(r(keep), abs(Fnum(keep)), 'o', 'MarkerSize', 3, ...
           'MarkerEdgeColor', [0.10 0.35 0.70], 'DisplayName', '比耐公式反推的力');
    hold on;
    loglog(r(keep), abs(Ftrue(keep)), '-', 'LineWidth', 1.8, 'Color', [0.85 0.33 0.10], ...
           'DisplayName', '标准答案 |F| ∝ 1/r²');
    grid on; box on;
    xlabel('r'); ylabel('|F(r)|');
    title('反推出的力（双对数坐标下是直线 → 幂次关系）');
    legend('Location', 'southwest');

    subplot(1,2,2);
    semilogy(r(keep), err, '-', 'LineWidth', 1.2, 'Color', [0.10 0.35 0.70]);
    grid on; box on;
    xlabel('r'); ylabel('相对误差');
    title(sprintf('反推误差（最大 %.1e）', max(err)));
end


%% ===================== 核心：解比耐方程的数值部分 =====================
function [theta, u, du] = solveBinet(Ffun, m, h, r0, vr0, thetaMax)
% 比耐公式是二阶常微分方程：  d²u/dθ² + u = -F(1/u)/(m h² u²)
% 把它拆成两个一阶方程，交给 MATLAB 自带的 ode45，以角度 θ 为自变量求解。
%
%   y(1) = u = 1/r
%   y(2) = du/dθ
%   d(y1)/dθ = y(2)
%   d(y2)/dθ = -y(1) - F(1/y(1)) / (m h² y(1)²)

    u0  = 1 / r0;        % 起点：u = 1/r
    du0 = -vr0 / h;      % 起点斜率：由 dr/dt = -h·(du/dθ) 反推

    odefun = @(th, y) [ y(2);
                        -y(1) - Ffun(1/y(1)) / (m * h^2 * y(1)^2) ];

    % 兜底：r 太大（跑没影了）或太小（撞上力心）就提前停，避免图画得没法看
    uTooSmall = 1 / (1000 * r0);      % r > 1000·r0 时停
    uTooBig   = 1 / (0.001 * r0);     % r < 0.001·r0 时停
    opts = odeset('RelTol', 1e-9, 'AbsTol', 1e-12, ...
                  'Events', @(th, y) stopEvent(th, y, uTooSmall, uTooBig));

    [theta, Y] = ode45(odefun, [0, thetaMax], [u0; du0], opts);
    u  = Y(:, 1);
    du = Y(:, 2);
end

function [value, isterminal, direction] = stopEvent(~, y, uTooSmall, uTooBig)
% 事件函数：u 穿过这两个阈值时，ode45 就结束积分
    value      = [y(1) - uTooSmall;  y(1) - uTooBig];
    isterminal = [1; 1];
    direction  = [-1; 1];
end


%% ===================== 画图 =====================
function f = newFigure(name, pos)
% 统一建图：白底，并尽量固定成浅色主题（R2025a 以后的版本才有 theme 函数，老版本跳过）
    f = figure('Name', name, 'Color', 'w', 'Position', pos);
    try
        theme(f, 'light');
    catch
    end
end

function plotOrbit(ex, theta, r, E, rAnalytic, periIdx, animate)
% 说明：开放轨道（逃逸）会一直跑到无穷远，为了图好看，画的时候在 r = 8·r_min 处截断
    isBound = mean(E) < 0;
    x = r .* cos(theta);
    y = r .* sin(theta);

    if isBound
        rCut = inf;                       % 闭合轨道本来就跑不远，不用截断
    else
        rCut = 8 * min(r);                % 开放轨道：只画到 8 倍最近距离
    end
    show = r <= rCut;                     % 图里要显示的点

    newFigure(['比耐公式 - ' ex.name], [60 60 1060 640]);

    % ---- (1) 轨道形状（左半边，画大一点）----
    ax1 = subplot(2,2,[1 3]);
    plot(x(show), y(show), '-', 'LineWidth', 1.8, 'Color', [0.10 0.35 0.70], ...
         'DisplayName', '数值解');
    hold on;
    if ~isempty(rAnalytic)      % 理论椭圆/双曲线，用来对照
        xa = rAnalytic .* cos(theta);
        ya = rAnalytic .* sin(theta);
        plot(xa(show), ya(show), '--', 'LineWidth', 1.2, ...
             'Color', [0.85 0.33 0.10], 'DisplayName', '理论解（对照）');
    end
    plot(0, 0, 'o', 'MarkerSize', 10, 'MarkerFaceColor', [0.9 0.25 0.15], ...
         'MarkerEdgeColor', 'none', 'DisplayName', '力心');
    if ~isempty(periIdx)
        pv = periIdx(r(periIdx) <= rCut);
        plot(x(pv), y(pv), 'v', 'MarkerSize', 7, 'LineWidth', 1, ...
             'MarkerFaceColor', [0.20 0.65 0.35], 'MarkerEdgeColor', 'k', ...
             'DisplayName', '近力心点');
    end
    grid on; box on; axis equal;
    xlabel('x'); ylabel('y');
    if isBound
        title('轨道（解出来的 r(θ)）');
    else
        title(sprintf('轨道（只画到 r = %.3g，之后一路跑向无穷远）', rCut));
    end
    legend('Location', 'southoutside', 'Orientation', 'horizontal', 'FontSize', 8);
    halfspan = 1.12 * max([abs(x(show)); abs(y(show))]);
    xlim(ax1, [-halfspan, halfspan]);
    ylim(ax1, [-halfspan, halfspan]);

    % ---- (2) 距离随角度变化 ----
    subplot(2,2,2);
    plot(theta, r, '-', 'LineWidth', 1.6, 'Color', [0.10 0.35 0.70]); hold on;
    if ~isempty(periIdx)
        plot(theta(periIdx), r(periIdx), 'v', 'MarkerSize', 7, 'LineWidth', 1, ...
             'MarkerFaceColor', [0.20 0.65 0.35], 'MarkerEdgeColor', 'k');
    end
    if ~isBound
        ylim([0, 1.2*rCut]);      % 不然 y 轴被那个巨大的 r 拉到没法看
    end
    grid on; box on;
    xlabel('\theta（弧度）'); ylabel('r（到力心的距离）');
    title('r 随 θ 的变化（谷底就是近力心点）');

    % ---- (3) 能量守恒检查 ----
    subplot(2,2,4);
    drift = (E - E(1)) / abs(E(1));
    plot(theta, drift, '-', 'LineWidth', 1.2, 'Color', [0.10 0.35 0.70]);
    grid on; box on;
    xlabel('\theta（弧度）'); ylabel('(E - E_0)/|E_0|');
    title(sprintf('能量守恒检查：最大偏差 %.1e', max(abs(drift))));

    % ---- (4) 可选动画：小圆点沿轨道跑 ----
    if animate
        idxShow = find(show);
        hP = plot(ax1, x(idxShow(1)), y(idxShow(1)), 'o', 'MarkerSize', 8, ...
                  'MarkerFaceColor', [0.9 0.25 0.15], 'MarkerEdgeColor', 'none', ...
                  'HandleVisibility', 'off');   % 不进图例，免得图例太长
        step = max(1, round(numel(idxShow)/150));   % 最多画 150 帧，够流畅了
        for k = 1:step:numel(idxShow)
            if ~ishandle(hP), break; end
            set(hP, 'XData', x(idxShow(k)), 'YData', y(idxShow(k)));
            drawnow limitrate;
        end
    end
end


%% ===================== 选题目的小工具 =====================
function id = chooseExample(P)
    fprintf('请选择要算的题目，输入数字后按回车（直接回车 = 例子 1）：\n');
    for i = 1:numel(P)
        fprintf('   %d) %s\n', i, P{i}.name);
    end
    try
        id = input('你的选择：');
    catch
        id = [];        % 非交互环境（比如批处理）里 input 用不了，就用默认值
    end
    if isempty(id) || ~isscalar(id) || id < 1 || id > numel(P) || id ~= round(id)
        id = 1;
    end
end


%% ===================== 题目清单（想做自己的题就改这里） =====================
function P = examples()
% 每个题目是一个 struct。想做自己的题，复制一个例子，改参数就行。
% 写 F(r)、V(r) 时一定要用 ./ 和 .^（点除、点乘方），这样才能一次算很多个 r。
% 字段含义：
%   .name     题目名字（会打印出来）
%   .mode     'orbit' = 已知力求轨道；'force' = 已知轨道求力
%   .m        质点质量
%   .Ffun     力 F(r)，引力写负号，斥力写正号
%   .Vfun     势能 V(r)，满足 F = -dV/dr（只用来检查能量守恒）
%   .r0 .vr0 .vt0   t = 0 时的距离、径向速度、横向（切向）速度
%   .thetaMax 积分到多少角度（2π 就是一圈）
%   .kAn      如果力恰好是 F = -k/r²，这里填 k，脚本会拿理论公式对答案；否则填 []
%   ---- 只有"已知轨道求力"（mode = 'force'）才用得上的字段 ----
%   .h        角动量（求力时它只决定反推出来的力的整体大小）
%   .rfun     轨道函数 r(θ)；.p .ecc 只是轨道 r = p/(1+e·cosθ) 的两个参数

    % ---------- 例 1：平方反比引力 → 椭圆（行星绕太阳就是这个）----------
    e1 = struct();
    e1.name   = '平方反比引力 → 椭圆轨道';
    e1.mode   = 'orbit';
    e1.m      = 1;                          % 质量
    e1.Ffun   = @(r) -3 ./ r.^2;            % 力：F = -k/r²，k = 1，负号=引力
    e1.Vfun   = @(r) 100000 ./ r;               % 势能：V = -k/r
    e1.r0     = 1;                          % 起点离力心 1
    e1.vr0    = 0;                          % 径向速度为 0 → 从最近点或最远点出发
    e1.vt0    = 0.8;                        % 横向速度，决定轨道是椭圆还是双曲线
    e1.thetaMax = 2*pi;                     % 正好积一圈
    e1.kAn    = 1;                          % 这个力是纯平方反比力，可以做理论对照
    e1.note   = '横向速度不太大 → 被引力绑住了，轨道是椭圆，这正对应开普勒第一定律。';

    % ---------- 例 2：同样的力，速度大一点 → 双曲线（逃逸）----------
    e2 = struct();
    e2.name   = '平方反比引力 → 双曲线轨道（逃逸）';
    e2.mode   = 'orbit';
    e2.m      = 1;
    e2.Ffun   = @(r) -1 ./ r.^2;
    e2.Vfun   = @(r) -1 ./ r;
    e2.r0     = 1;
    e2.vr0    = 0;
    e2.vt0    = 1.6;                        % 速度超过了逃逸速度 → 跑掉
    e2.thetaMax = 3;
    e2.kAn    = 1;
    e2.note   = '横向速度超过 √(2k/(m·r0)) 就绑不住了：机械能为正，轨道是双曲线，质点跑到无穷远。';

    % ---------- 例 3：平方反比力 + 一点额外吸引力 → 轨道进动 ----------
    e3 = struct();
    e3.name   = '平方反比力 + 微扰 → 轨道进动';
    e3.mode   = 'orbit';
    e3.m      = 1;
    % 力 = -k/r² 再补一项 -α/r⁴，也就是 F = -(1/r²)·(1 + α/r²)，这里取 α = 0.01
    e3.Ffun   = @(r) -(1 ./ r.^2) .* (1 + 0.01 ./ r.^2);
    e3.Vfun   = @(r) -1 ./ r - 0.01 ./ (3 * r.^3);         % 对应势能 V = -k/r - kα/(3r³)
    e3.r0     = 1;
    e3.vr0    = 0;
    e3.vt0    = 0.8;
    e3.thetaMax = 6*pi;                     % 多转几圈，才看得出进动
    e3.kAn    = [];                         % 不是纯平方反比力，没法解析对照
    e3.note   = ['力比起平方反比力多了一小项（α = 0.01）：轨道不再闭合，' ...
                 '近力心点一圈一圈往前挪，这就是"进动"（水星近日点进动的简化版）。'];

    % ---------- 例 4：已知轨道，反推是什么力 ----------
    e4 = struct();
    e4.name   = '已知圆锥曲线轨道 r(θ) → 反推力 F(r)';
    e4.mode   = 'force';
    e4.m      = 1;
    e4.h      = 1;                          % 求力时角动量要自己给（它只决定力的整体大小）
    e4.p      = 1;                          % 轨道 r = p/(1 + e·cosθ)
    e4.ecc    = 0.6;
    e4.rfun   = @(theta) 1 ./ (1 + 0.6*cos(theta));
    e4.note   = '把轨道的二阶导数代进比耐公式，反推出来的力应该正好是平方反比力 F ∝ -1/r²。';

    P = {e1, e2, e3, e4};       % 用 {} 打包，因为各题的字段不完全一样
end
