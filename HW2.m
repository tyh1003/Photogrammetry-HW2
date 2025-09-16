% 讀取數據
data_read = readmatrix('output.txt'); % 將 'your_file.txt' 替換為你的檔案名
% 將數據存入表格
data = array2table(data_read, 'VariableNames', {'no', 'x1', 'y1', 'x2', 'y2'});
%輸入基本資料
phie_L = input(' > initial value of phie_L  = ? deg. (e.g. 0.0): ');
kapa_L = input(' > initial value of kapa_L  = ? deg. (e.g. 0.0): ');
omega_R = input(' > initial value of omega_R = ? deg. (e.g. 0.0): ');
phie_R = input(' > initial value of phie_R = ? deg. (e.g. 0.0): ');
kapa_R = input(' > initial value of kapa_R  = ? deg. (e.g. 0.0): ');
B = input(' > A reference baseline B = ? mm. (e.g. 40)');
threshold = input(' > Threshold value adopted in convergence condition = ? (e.g. 0.00000001): ');
f = input('Principal distance = ? mm. (e.g. 152.818)');
%自由度
n=height(data);
u=5;
freedom=n-u;
% 顯示初始值
fprintf('\nInitial values of 5 R.O. parameters (units: degree):\n');
fprintf('phie_L =  %.6f : phie  angle of the left  photo\n', phie_L);
fprintf('kapa_L =  %.6f : kapa  angle of the left  photo\n', kapa_L);
fprintf('omega_R = %.6f : omega angle of the right photo\n', omega_R);
fprintf('phie_R =  %.6f : phie  angle of the right photo\n', phie_R);
fprintf('kapa_R =  %.6f : kapa  angle of the right photo\n', kapa_R);
% 迭代過程
fprintf('\nComputation converges if max. |X| is less than %.8f\n\n', threshold);
fprintf('Iteration Results: (units: radians)\n');
[phie_L, kapa_L, omega_R, phie_R, kapa_R, sum,std_errors] = deda(data, phie_L, kapa_L, omega_R, phie_R, kapa_R, f, B, threshold, freedom);
%自由度
fprintf('\nN = %.i  ; n = %.i ; u = %.i ;\n', height(data),n,u);
fprintf('\nDegree of freedom = %.i\n', freedom);
fprintf('\nPhie_L、Kapa_L、Omega_R、Phie_R、Kapa_R= \n');
% 以度數輸出結果
fprintf('Phie_L =  %.6f +/- %.6f(deg)\n', rad2deg(phie_L), rad2deg(std_errors(1)));
fprintf('Kapa_L =  %.6f +/- %.6f(deg)\n', rad2deg(kapa_L), rad2deg(std_errors(2)));
fprintf('Omega_R = %.6f +/- %.6f(deg)\n', rad2deg(omega_R), rad2deg(std_errors(3)));
fprintf('Phie_R =  %.6f +/- %.6f(deg)\n', rad2deg(phie_R), rad2deg(std_errors(4)));
fprintf('Kapa_R =  %.6f +/- %.6f(deg)\n\n', rad2deg(kapa_R), rad2deg(std_errors(5)));
%輸出data值
fprintf('Input data:\nNO.、xl(mm)、yl(mm)、xr(mm)、yr(mm)\n');
% 遍歷每一行，並格式化輸出
for i = 1:height(data)
    fprintf('%8.0f %11.6f %11.6f %11.6f %11.6f\n', data.no(i), data.x1(i), data.y1(i), data.x2(i), data.y2(i));
end
%fprintf('Correction value:\nNO.、xl correction(mm)、yl correction(mm)、xr correction(mm)、yr correction(mm)\n');
%for i = 1:height(change)
    %fprintf('%8.0f %11.6f %11.6f %11.6f %11.6f\n', data.no(i),change(1,i), change(2,i), change(3,i), change(4,i));
%end
fprintf('\n前四項:點號、U、V、W = model coordinates (units: mm)\n');
fprintf('第五項:F (units: mm^2)；第六項:改正數(F*B) (units: mm^3)\n');
% 輸出 uvw，並顯示到小數點後三位
[uvw,rms,std_fmm] = calculate_uvw(data, phie_L, kapa_L, omega_R, phie_R, kapa_R, f, B );
for i = 1:size(uvw, 1)
    fprintf('%8.0f %8.3f %8.3f %8.3f %15.12f %11.6f\n', uvw(i, 1), uvw(i, 2), uvw(i, 3), uvw(i, 4),uvw(i, 5),uvw(i,6));
end
%其他值
fprintf('\nRMS value of all pseudo-observations of "volumes" = %.6f (mm^3)\n', rms);
fprintf('\nstandard deviation of unit weight = %.8f (mm^2) \n', std_fmm);


% 旋轉矩陣
function R_M = Rot_M(omega, phie, kapa)
    m11 = cos(phie)*cos(kapa);
    m12 = sin(omega)*sin(phie)*cos(kapa) + cos(omega)*sin(kapa);
    m13 = -cos(omega)*sin(phie)*cos(kapa) + sin(omega)*sin(kapa);
    m21 = -cos(phie)*sin(kapa);
    m22 = -sin(omega)*sin(phie)*sin(kapa) + cos(omega)*cos(kapa);
    m23 = cos(omega)*sin(phie)*sin(kapa) + sin(omega)*cos(kapa);
    m31 = sin(phie);
    m32 = -sin(omega)*cos(phie);
    m33 = cos(omega)*cos(phie);
    R_M = [m11, m12, m13; m21, m22, m23; m31, m32, m33];
end

%計算detla
function [V1, V2, V3, det_val] = Detla(x1, y1, x2, y2, f, B, p_1, k_1, o_2, p_2, k_2)
    M1 = Rot_M(0, p_1, k_1)';
    M2 = Rot_M(o_2, p_2, k_2)';
    V1 = [B; 0; 0];
    p1 = [x1; y1; -f];
    p2 = [x2; y2; -f];
    V2 = M1 * p1;
    V3 = M2 * p2;
    det_val = dot(V1, cross(V2, V3));
end

% 計算 b2
function b2_value = b2(x1, y1, f, omega_1, phi_1, kappa_1, V1, V3)
    m21 = -x1*sin(phi_1)*cos(kappa_1) + y1*sin(phi_1)*sin(kappa_1) - f*cos(phi_1);
    m22 = x1*sin(omega_1)*cos(phi_1)*cos(kappa_1) - y1*sin(omega_1)*cos(phi_1)*sin(kappa_1) - f*sin(omega_1)*sin(phi_1);
    m23 = -x1*cos(omega_1)*cos(phi_1)*cos(kappa_1) + y1*cos(omega_1)*cos(phi_1)*sin(kappa_1) + f*cos(omega_1)*sin(phi_1);
    V2 = [m21; m22; m23];
    b2_value = dot(V1, cross(V2, V3));
end

% 計算 b3
function b3_value = b3(x1, y1, V1, V3, M1)
    m21 = x1*M1(2,1) - y1*M1(1,1);
    m22 = x1*M1(2,2) - y1*M1(1,2);
    m23 = x1*M1(2,3) - y1*M1(1,3);
    V2 = [m21; m22; m23];
    b3_value = dot(V1, cross(V2, V3));
end

% 計算 b7
function b7_value = b7(V1, V2, V3)
    V3 = [0; -V3(3); V3(2)];
    b7_value = dot(V1, cross(V2, V3));
end

% 計算 b8
function b8_value = b8(x2, y2, f, omega_2, phi_2, kappa_2, V1, V2)
    m31 = -x2*sin(phi_2)*cos(kappa_2) + y2*sin(phi_2)*sin(kappa_2) - f*cos(phi_2);
    m32 = x2*sin(omega_2)*cos(phi_2)*cos(kappa_2) - y2*sin(omega_2)*cos(phi_2)*sin(kappa_2) - f*sin(omega_2)*sin(phi_2);
    m33 = -x2*cos(omega_2)*cos(phi_2)*cos(kappa_2) + y2*cos(omega_2)*cos(phi_2)*sin(kappa_2) + f*cos(omega_2)*sin(phi_2);
    V3 = [m31; m32; m33];
    b8_value = dot(V1, cross(V2, V3));
end

% 計算 b9
function b9_value = b9(x2, y2, V1, V2, M2)
    m31 = x2*M2(2,1) - y2*M2(1,1);
    m32 = x2*M2(2,2) - y2*M2(1,2);
    m33 = x2*M2(2,3) - y2*M2(1,3);
    V3 = [m31; m32; m33];
    b9_value = dot(V1, cross(V2, V3));
end

% 構建設計矩陣的行
function row = A_r(x1, y1, x2, y2, f, p_1, k_1, o_2, p_2, k_2, V1, V2, V3, M1, M2)
    b2_val = b2(x1, y1, f, 0, p_1, k_1, V1, V3);
    b3_val = b3(x1, y1, V1, V3, M1);
    b7_val = b7(V1, V2, V3);
    b8_val = b8(x2, y2, f, o_2, p_2, k_2, V1, V2);
    b9_val = b9(x2, y2, V1, V2, M2);
    row = [b2_val, b3_val, b7_val, b8_val, b9_val];
end

% 構建設計矩陣 A 和偏差向量 h0
function [A, h0] = cmatrix(data, p_1, k_1, o_2, p_2, k_2, f, B)
    A = [];
    h0 = [];
    for i = 1:height(data)
        x1 = data.x1(i);
        y1 = data.y1(i);
        x2 = data.x2(i);
        y2 = data.y2(i);
        M1 = Rot_M(0, p_1, k_1);
        M2 = Rot_M(o_2, p_2, k_2);
        [V1, V2, V3, det_val] = Detla(x1, y1, x2, y2, f, B, p_1, k_1, o_2, p_2, k_2);
        arow = A_r(x1, y1, x2, y2, f, p_1, k_1, o_2, p_2, k_2, V1, V2, V3, M1, M2);
        A = [A; arow];
        h0 = [h0; -det_val];
        
    end
    
end
% 計算 delta 值
function [del_phi_1, del_kappa_1, del_omega_2, del_phi_2, del_kappa_2, sum] = new(A, h0)
    delta = A \ h0; % 解最小二乘問題
    del_phi_1 = delta(1);
    del_kappa_1 = delta(2);
    del_omega_2 = delta(3);
    del_phi_2 = delta(4);
    del_kappa_2 = delta(5);
    sum = del_phi_1^2 + del_kappa_1^2 + del_omega_2^2 + del_phi_2^2 + del_kappa_2^2;
end
% 迭代過程
function [p_1, k_1, o_2, p_2, k_2, sum,std_errors] = deda(data, p_1, k_1, o_2, p_2, k_2, f, B, threshold,freedom)
    max = 20;
    while max > 0
        [A, h0] = cmatrix(data, p_1, k_1, o_2, p_2, k_2, f, B);
        [delta_phi_1, delta_kappa_1, delta_omega_2, delta_phi_2, delta_kappa_2, sum] = new(A, h0);
        p_1 = p_1 + delta_phi_1;
        k_1 = k_1 + delta_kappa_1;
        o_2 = o_2 + delta_omega_2;
        p_2 = p_2 + delta_phi_2;
        k_2 = k_2 + delta_kappa_2;
        
        deltax=[delta_phi_1,delta_kappa_1,delta_omega_2,delta_phi_2,delta_kappa_2]';
      
        % 打印每次迭代的信息（可選）
        fprintf('Iteration %d: d_phie(L)=%-12.9f, d_kapa(L)=%-12.9f, d_omega(R)=%-12.9f, d_phie(R)=%-12.9f, d_kapa(R)=%-12.9f\n', ...
            21 - max, delta_phi_1, delta_kappa_1, delta_omega_2, delta_phi_2, delta_kappa_2);
        if sum < threshold^2
            v=A*deltax-h0;
            sigma=(v'*v)/freedom;
            cov=sigma*(inv(A'*A));
            std_errors=sqrt(diag(cov));
            return;
        end
        max = max - 1;
    end
    v=A*deltax-h0;
    sigma=(v'*v)/freedom;
    cov=sigma*(inv(A'*A));
    std_errors=sqrt(diag(cov));
end


function [uvw,rms,std_fmm] = calculate_uvw(data, phie_L, kapa_L, omega_R, phie_R, kapa_R, f, B )
    % 初始化模型點矩陣
    squaresum=0;
    uvw = zeros(height(data), 6);
    for i = 1:height(data)
        x1 = data.x1(i);
        y1 = data.y1(i);
        x2 = data.x2(i);
        y2 = data.y2(i);
        no=data.no(i);

        % 計算方向向量
        V1 = [x1; y1; -f];
        V2 = [x2; y2; -f];
        M1 = Rot_M(0, phie_L, kapa_L)';
        M2 = Rot_M(omega_R, phie_R, kapa_R)';
        or1=M1*V1;
        or2=M2*V2;
        
        A = [or1,-or2];
        B_vec = [B; 0;0];
        ab = (A' * A) \ (A' * B_vec);
        a = ab(1);
        b = ab(2);
        
        % 計算 u, v, w
        uvw_vec1 = a * or1; 

        uvw_vec2=b*or2;
        u = uvw_vec1(1);
        v = uvw_vec1(2);
        w = uvw_vec1(3);
        
        fmm=or1(2)*or2(3)-or2(2)*or1(3);
        vol=fmm*B;
        squaresum=squaresum+vol*vol;
        
        % 存入結果
        uvw(i, 1) = no;           % 存放點號
        uvw(i, 2:4) = [u,v,w];
        uvw(i,5)=fmm;
        uvw(i,6)=vol;
    end
    rms=sqrt(squaresum/height(data));
    std_fmm = std(uvw(:, 5));
    
end

