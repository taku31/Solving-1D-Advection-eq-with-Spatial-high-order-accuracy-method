% 初期化
clearvars;

% パラメーター
i_max = 100;    % 格子セル数
XL = -1.0;      % 計算領域左端の座標
XR = 1.0;       % 計算領域右端の座標
a = 1.0;        % 線形移流方程式の移流速度
tstop = 1;      % 計算停止時刻
eps = 1;        % 高次精度解法のパラメーター
kappa = 0;    % 高次精度解法のパラメーター

% 配列定義
i = 0;                      % セル番号； (i番目セルの左境界の番号はi、右辺境界の番号はi+1とする)
x = zeros(i_max + 1, 1);    % セル境界の座標
u = zeros(i_max, 1);        % セル平均値 （数値解）
ue = zeros(i_max, 1);       % セル平均値 （厳密解）
ul = zeros(i_max + 1, 1);   % セル境界左側の変数値
ur = zeros(i_max + 1, 1);   % セル境界右側の変数値
f = zeros(i_max + 1, 1);    % セル境界の流束
n = 0;                      % 時間ステップ
t = 0;                      % 計算時間

% 方程式を選択
% 1:線形移流方程式、2:非粘性Burgers方程式
sw1 = 1;

% 初期値の選択
% 1:不連続な分布、2:滑らかな分布
sw2 = 1;

% メッシュの設定
dx = (XR - XL) / (i_max - 4.0);         % 格子間隔 計算領域外に二つずつ余分なセルを準備。周期的境界条件向けの設定
dt = 0.2 * dx;                          % 時間刻み
x(1) = XL - 2.0 * dx;                   % 計算領域外の２セルも考慮した座標を振る。
[x, u] = initc(sw2, i_max, x, dx, u);   % 計算格子，時間刻み，初期条件を設定する

%　厳密解の計算
ue = exact(sw1, sw2, i_max, ue, x, t, dx);

% メインループ
while t <= tstop

    % 時間発展
    n = n + 1;
    t = t + dt;

    % 空間再構築
    [ul, ur] = reconstruction_pc(i_max, u, ul, ur, eps, kappa);
    
    % リーマンソルバー
    f = riemann_roe(i_max, f, ul, ur, sw1);

    % 時間積分
    u = update(i_max, u, dt, dx, f);

    % 境界条件
    u = bc(i_max, u);

    % 厳密解を求める
    ue = exact(sw1, sw2, i_max, ue, x, t, dx);

    % 時間表示
    fprintf("n=%d, t=%f \n", n, t);

    % 動画保存
    if n == 1
        plotconfig(x(1 : end - 1), ue, u, t, eps, kappa)
        filename = ['sw1 = ', num2str(sw1, '%.0f'), ', sw2 = ', num2str(sw2, '%.0f'), ', epsilon = ', num2str(eps, '%.0f'), ', kappa = ', num2str(kappa, '%.3f'),'.mp4'];
        v = VideoWriter(filename,'MPEG-4');
        v.FrameRate = 40;
        open(v);
    else
        plotconfig(x(1 : end - 1), ue, u, t, eps, kappa)
        frame = getframe(gcf);
        writeVideo(v,frame);
    end

end

% 動画ファイルを閉じる
close(v);

%% 以下ローカル関数

% 初期値の設定
function [x, u] = initc(sw2, i_max, x, dx, u)

for  i = 2 : i_max + 1
    x(i) = x(i - 1) + dx;% 格子点の座標
end

switch sw2

    case 1 % 初期の変数値（滑らかな分布）

        for  i = 1 : i_max
            u(i) = 0.1;
        end

        for i = i_max / 2 - 10 : i_max / 2 + 10
            u(i) = 1.0;
        end

    case 2 % 初期の変数値（不連続な分布）

        for  i = 1 : i_max
            u(i) = 0.5*(1.1 + sin(2.* pi *(x(i)-x(3))));% 三番目の要素（XL、XRの座標値）を基準に考える。
        end

end

end

% 厳密解の計算(不連続分布)
function [ue] = exact(sw1, sw2, i_max, ue, x, t, dx)

switch sw2

    case 1 % 不連続分布

        switch sw1

            case 1 % 線形移流方程式

                alpha_12 = 1;% 移流速度
                xc = alpha_12 * t;
                xl = xc - dx * 10;% 中央に対して10要素だけマイナスに移動した位置

                % 周期境界条件　
                if  xl > 1.0
                    xl = -2.0 + xl;
                end

                xr = xc + dx * 10.;

                % 周期境界条件　
                if xr > 1.0
                    xr = -2.0 + xr;
                end

                if  xl <= xr
                    for  i = 1 : i_max
                        ue(i) = 0.1;
                        if dx * (i - i_max / 2) >= xl && dx * (i - i_max / 2) <= xr
                            ue(i) = 1.0;
                        end
                    end
                end

                if xl >= xr
                    for i = 1 : i_max
                        ue(i) = 1.0;
                        if dx * (i - i_max / 2) >= xr && dx * (i - i_max / 2) <= xl
                            ue(i) = 0.1;
                        end
                    end
                end

            case 2 % 非粘性burgers方程式

                xc = - dx * 10 + t;
                xl = - dx * 10 + 0.1 * t;
                xr = dx * 10 + 0.55 * t;

                % 周期境界条件　
                if  xl > 1.0
                    xl = -2.0 + xl;
                end

                % 周期境界条件　
                if xr > 1.0
                    xr = -2.0 + xr;
                end

                % 周期境界条件　
                if xc > 1.0
                    xc = -2.0 + xc;
                end

                for i = 1 : i_max
                    if x(i) <= xl
                        ue(i) = 0.1;
                    end
                    if x(i) >= xl && x(i) <= xc
                        ue(i) =(x(i) - xl) / (xc - xl) * 0.9 + 0.1;
                    end
                    if x(i) >= xc && x(i) <= xr
                        ue(i) = 1.0;
                    end
                    if x(i) >= xr
                        ue(i) = 0.1;
                    end
                end

        end

    case 2 % 連続分布

        switch sw1

            case 1 % 線形移流方程式

                alpha_12 = 1; % 移流速度
                for i = 1 : i_max
                    ue(i) = 0.5 * (1.1 + sin(2. * pi * ((x(i) - x(3)) - alpha_12 * t)));
                end

            case 2  % 非粘性burgers方程式

                for i = 1 : i_max

                    c = 2 * pi;
                    f = ue(i) - 0.5 * (1.1 + sin(c * (x(i)-x(3) - ue(i) * t)));
                    df = 1.0 + 0.5 * c * cos(c *(x(i) - x(3) - ue(i) * t)) * t;
                    count = 0;
                    while abs(f) >= 1.0e-6
                        count = count + 1;
                        ue(i) = ue(i) - f / df; % ニュートン法の計算式*/
                        f = ue(i) - 0.5 * (1.1 + sin(c*((x(i) - x(3) - ue(i) * t))));
                        df = 1.0 + 0.5 * c * cos(c * ((x(i) - x(3) - ue(i) * t))) * t;
                        if count > 10000
                            disp("厳密解が収束しないので反復を打ち切ります。");
                            break
                        end
                    end

                end
        end
end

end

% 厳密解の計算(初期値が滑らかな分布の場合）
function [ue] = exact2(sw1, i_max, ue, x, t, dx)

switch sw1

    case 1 % 線形移流方程式

        alpha_12 = 1; % 移流速度
        for i = 1 : i_max
            ue(i) = 0.5 * (1.1 + sin(2. * pi * ((x(i) - x(3)) - alpha_12 * t)));
        end

    case 2  % 非粘性burgers方程式

        for i = 1 : i_max

            c = 2 * pi;
            f = ue(i) - 0.5 * (1.1 + sin(c * (x(i)-x(3) - ue(i) * t)));
            df = 1.0 + 0.5 * c * cos(c *(x(i) - x(3) - ue(i) * t)) * t;
            count = 0;
            while abs(f) >= 1.0e-6
                count = count + 1;
                ue(i) = ue(i) - f / df; % ニュートン法の計算式*/
                f = ue(i) - 0.5 * (1.1 + sin(c*((x(i) - x(3) - ue(i) * t))));
                df = 1.0 + 0.5 * c * cos(c * ((x(i) - x(3) - ue(i) * t))) * t;
                if count > 10000
                    disp("厳密解が収束しないので反復を打ち切ります。");
                    break
                end
            end

        end
end

end


function [ul, ur] = reconstruction_pc(i_max, u, ul, ur, eps, kappa)

for i = 2 : i_max - 2
    ul(i + 1) = u(i) + eps * (0.25 * (1 - kappa) * (u(i) - u(i - 1))...
        + 0.25 * (1 + kappa) * (u(i + 1) - u(i))); % セル境界(i+1/2)左側の値
    ur(i + 1) = u(i + 1) + eps *(- 0.25 * (1 + kappa) * (u(i + 1) - u(i))...
        - 0.25 * (1 - kappa) * (u(i + 2) - u(i + 1))); % セル境界(i+1/2)右側の値
end

end


function [f] = riemann_roe(i_max, f, ul, ur, sw1) % 流束を計算する

% 移流速度の計算
switch  sw1

    case 1 % 線形移流方程式の数値流束

        for i = 3 : i_max - 1

            alpha_12 = 1;% 移流速度
            f_flux_ul = alpha_12 * ul(i);
            f_flux_ur = alpha_12 * ur(i);
            f(i) = 1.0 / 2.0 * (f_flux_ul + f_flux_ur) - 1.0 / 2.0 * abs(alpha_12) * (ur(i) - ul(i));

        end

    case 2 % 非粘性burgers方程式の数値流束

        for i = 3 : i_max - 1

            alpha_12 = 0.5 * (ur(i) + ul(i));% 移流速度
            f_flux_ul = 0.5 * ul(i) * ul(i);
            f_flux_ur = 0.5 * ur(i) * ur(i);
            f(i) = 1.0 / 2.0 * (f_flux_ul + f_flux_ur) - 1.0 / 2.0 * abs(alpha_12) * (ur(i) - ul(i));

        end

end

end

function [u] = update(i_max, u, dt, dx, f)

for  i = 3 : i_max - 2
    u(i) = u(i) - dt / dx * (f(i + 1) - f(i));% 計算変数を更新する
end

end

function [u] = bc(i_max, u) % 周期境界条件

u(1) = u(i_max - 3); % 計算領域左端の境界条件
u(2) = u(i_max - 2); % 計算領域左端の境界条件
u(i_max - 1) = u(3); % 計算領域右端の境界条件
u(i_max) = u(4); % 計算領域右端の境界条件

end

function [] = plotconfig(x, ue, u, t, eps, kappa)

plot(x, ue, x, u)
%plot(x, u)

title(['time = ', num2str(t, '%.3f')]);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);
axis equal;
axis tight;
axis on;
fig=gcf;
fig.Color='white';
ylim([-0.2 1.5]);
xlabel('x')
ylabel('u')

% 凡例
method = ['\epsilon = ', num2str(eps, '%.0f'), ', \kappa = ', num2str(kappa, '%.3f')];
legend({'exact', method},'Location','southwest','FontSize', 10)

% 新しいプロットの時、軸設定を保持したまま前のグラフィックスオブジェクトを消去
set(gca,'nextplot','replacechildren');

end
