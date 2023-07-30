% ������
clearvars;

% �p�����[�^�[
i_max = 100;    % �i�q�Z����
XL = -1.0;      % �v�Z�̈捶�[�̍��W
XR = 1.0;       % �v�Z�̈�E�[�̍��W
a = 1.0;        % ���`�ڗ��������̈ڗ����x
tstop = 1;      % �v�Z��~����
eps = 1;        % �������x��@�̃p�����[�^�[
kappa = 0;    % �������x��@�̃p�����[�^�[

% �z���`
i = 0;                      % �Z���ԍ��G (i�ԖڃZ���̍����E�̔ԍ���i�A�E�Ӌ��E�̔ԍ���i+1�Ƃ���)
x = zeros(i_max + 1, 1);    % �Z�����E�̍��W
u = zeros(i_max, 1);        % �Z�����ϒl �i���l���j
ue = zeros(i_max, 1);       % �Z�����ϒl �i�������j
ul = zeros(i_max + 1, 1);   % �Z�����E�����̕ϐ��l
ur = zeros(i_max + 1, 1);   % �Z�����E�E���̕ϐ��l
f = zeros(i_max + 1, 1);    % �Z�����E�̗���
n = 0;                      % ���ԃX�e�b�v
t = 0;                      % �v�Z����

% ��������I��
% 1:���`�ڗ��������A2:��S��Burgers������
sw1 = 1;

% �����l�̑I��
% 1:�s�A���ȕ��z�A2:���炩�ȕ��z
sw2 = 1;

% ���b�V���̐ݒ�
dx = (XR - XL) / (i_max - 4.0);         % �i�q�Ԋu �v�Z�̈�O�ɓ���]���ȃZ���������B�����I���E���������̐ݒ�
dt = 0.2 * dx;                          % ���ԍ���
x(1) = XL - 2.0 * dx;                   % �v�Z�̈�O�̂Q�Z�����l���������W��U��B
[x, u] = initc(sw2, i_max, x, dx, u);   % �v�Z�i�q�C���ԍ��݁C����������ݒ肷��

%�@�������̌v�Z
ue = exact(sw1, sw2, i_max, ue, x, t, dx);

% ���C�����[�v
while t <= tstop

    % ���Ԕ��W
    n = n + 1;
    t = t + dt;

    % ��ԍč\�z
    [ul, ur] = reconstruction_pc(i_max, u, ul, ur, eps, kappa);
    
    % ���[�}���\���o�[
    f = riemann_roe(i_max, f, ul, ur, sw1);

    % ���Ԑϕ�
    u = update(i_max, u, dt, dx, f);

    % ���E����
    u = bc(i_max, u);

    % �����������߂�
    ue = exact(sw1, sw2, i_max, ue, x, t, dx);

    % ���ԕ\��
    fprintf("n=%d, t=%f \n", n, t);

    % ����ۑ�
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

% ����t�@�C�������
close(v);

%% �ȉ����[�J���֐�

% �����l�̐ݒ�
function [x, u] = initc(sw2, i_max, x, dx, u)

for  i = 2 : i_max + 1
    x(i) = x(i - 1) + dx;% �i�q�_�̍��W
end

switch sw2

    case 1 % �����̕ϐ��l�i���炩�ȕ��z�j

        for  i = 1 : i_max
            u(i) = 0.1;
        end

        for i = i_max / 2 - 10 : i_max / 2 + 10
            u(i) = 1.0;
        end

    case 2 % �����̕ϐ��l�i�s�A���ȕ��z�j

        for  i = 1 : i_max
            u(i) = 0.5*(1.1 + sin(2.* pi *(x(i)-x(3))));% �O�Ԗڂ̗v�f�iXL�AXR�̍��W�l�j����ɍl����B
        end

end

end

% �������̌v�Z(�s�A�����z)
function [ue] = exact(sw1, sw2, i_max, ue, x, t, dx)

switch sw2

    case 1 % �s�A�����z

        switch sw1

            case 1 % ���`�ڗ�������

                alpha_12 = 1;% �ڗ����x
                xc = alpha_12 * t;
                xl = xc - dx * 10;% �����ɑ΂���10�v�f�����}�C�i�X�Ɉړ������ʒu

                % �������E�����@
                if  xl > 1.0
                    xl = -2.0 + xl;
                end

                xr = xc + dx * 10.;

                % �������E�����@
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

            case 2 % ��S��burgers������

                xc = - dx * 10 + t;
                xl = - dx * 10 + 0.1 * t;
                xr = dx * 10 + 0.55 * t;

                % �������E�����@
                if  xl > 1.0
                    xl = -2.0 + xl;
                end

                % �������E�����@
                if xr > 1.0
                    xr = -2.0 + xr;
                end

                % �������E�����@
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

    case 2 % �A�����z

        switch sw1

            case 1 % ���`�ڗ�������

                alpha_12 = 1; % �ڗ����x
                for i = 1 : i_max
                    ue(i) = 0.5 * (1.1 + sin(2. * pi * ((x(i) - x(3)) - alpha_12 * t)));
                end

            case 2  % ��S��burgers������

                for i = 1 : i_max

                    c = 2 * pi;
                    f = ue(i) - 0.5 * (1.1 + sin(c * (x(i)-x(3) - ue(i) * t)));
                    df = 1.0 + 0.5 * c * cos(c *(x(i) - x(3) - ue(i) * t)) * t;
                    count = 0;
                    while abs(f) >= 1.0e-6
                        count = count + 1;
                        ue(i) = ue(i) - f / df; % �j���[�g���@�̌v�Z��*/
                        f = ue(i) - 0.5 * (1.1 + sin(c*((x(i) - x(3) - ue(i) * t))));
                        df = 1.0 + 0.5 * c * cos(c * ((x(i) - x(3) - ue(i) * t))) * t;
                        if count > 10000
                            disp("���������������Ȃ��̂Ŕ�����ł��؂�܂��B");
                            break
                        end
                    end

                end
        end
end

end

% �������̌v�Z(�����l�����炩�ȕ��z�̏ꍇ�j
function [ue] = exact2(sw1, i_max, ue, x, t, dx)

switch sw1

    case 1 % ���`�ڗ�������

        alpha_12 = 1; % �ڗ����x
        for i = 1 : i_max
            ue(i) = 0.5 * (1.1 + sin(2. * pi * ((x(i) - x(3)) - alpha_12 * t)));
        end

    case 2  % ��S��burgers������

        for i = 1 : i_max

            c = 2 * pi;
            f = ue(i) - 0.5 * (1.1 + sin(c * (x(i)-x(3) - ue(i) * t)));
            df = 1.0 + 0.5 * c * cos(c *(x(i) - x(3) - ue(i) * t)) * t;
            count = 0;
            while abs(f) >= 1.0e-6
                count = count + 1;
                ue(i) = ue(i) - f / df; % �j���[�g���@�̌v�Z��*/
                f = ue(i) - 0.5 * (1.1 + sin(c*((x(i) - x(3) - ue(i) * t))));
                df = 1.0 + 0.5 * c * cos(c * ((x(i) - x(3) - ue(i) * t))) * t;
                if count > 10000
                    disp("���������������Ȃ��̂Ŕ�����ł��؂�܂��B");
                    break
                end
            end

        end
end

end


function [ul, ur] = reconstruction_pc(i_max, u, ul, ur, eps, kappa)

for i = 2 : i_max - 2
    ul(i + 1) = u(i) + eps * (0.25 * (1 - kappa) * (u(i) - u(i - 1))...
        + 0.25 * (1 + kappa) * (u(i + 1) - u(i))); % �Z�����E(i+1/2)�����̒l
    ur(i + 1) = u(i + 1) + eps *(- 0.25 * (1 + kappa) * (u(i + 1) - u(i))...
        - 0.25 * (1 - kappa) * (u(i + 2) - u(i + 1))); % �Z�����E(i+1/2)�E���̒l
end

end


function [f] = riemann_roe(i_max, f, ul, ur, sw1) % �������v�Z����

% �ڗ����x�̌v�Z
switch  sw1

    case 1 % ���`�ڗ��������̐��l����

        for i = 3 : i_max - 1

            alpha_12 = 1;% �ڗ����x
            f_flux_ul = alpha_12 * ul(i);
            f_flux_ur = alpha_12 * ur(i);
            f(i) = 1.0 / 2.0 * (f_flux_ul + f_flux_ur) - 1.0 / 2.0 * abs(alpha_12) * (ur(i) - ul(i));

        end

    case 2 % ��S��burgers�������̐��l����

        for i = 3 : i_max - 1

            alpha_12 = 0.5 * (ur(i) + ul(i));% �ڗ����x
            f_flux_ul = 0.5 * ul(i) * ul(i);
            f_flux_ur = 0.5 * ur(i) * ur(i);
            f(i) = 1.0 / 2.0 * (f_flux_ul + f_flux_ur) - 1.0 / 2.0 * abs(alpha_12) * (ur(i) - ul(i));

        end

end

end

function [u] = update(i_max, u, dt, dx, f)

for  i = 3 : i_max - 2
    u(i) = u(i) - dt / dx * (f(i + 1) - f(i));% �v�Z�ϐ����X�V����
end

end

function [u] = bc(i_max, u) % �������E����

u(1) = u(i_max - 3); % �v�Z�̈捶�[�̋��E����
u(2) = u(i_max - 2); % �v�Z�̈捶�[�̋��E����
u(i_max - 1) = u(3); % �v�Z�̈�E�[�̋��E����
u(i_max) = u(4); % �v�Z�̈�E�[�̋��E����

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

% �}��
method = ['\epsilon = ', num2str(eps, '%.0f'), ', \kappa = ', num2str(kappa, '%.3f')];
legend({'exact', method},'Location','southwest','FontSize', 10)

% �V�����v���b�g�̎��A���ݒ��ێ������܂ܑO�̃O���t�B�b�N�X�I�u�W�F�N�g������
set(gca,'nextplot','replacechildren');

end
