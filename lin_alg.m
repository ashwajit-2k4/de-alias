alp = 0.9;
N = 40;

A = zeros(100,N);

for i = 1:50
    A(i,1) = abs(1/(1-alp^N*exp(-1i*0.00157*i*N)));
    for j = 2:N
        A(i,j) = alp^(j-1)*2*cos(0.00157*i*(j-1))*abs(1/(1-alp^N*exp(-1i*0.00157*i*N)));
    end
end

for i = 51:100
    A(i,1) = abs(1/(1-alp^N*exp(-1i*(0.05495*(i-50)+0.33755)*N)));
    for j = 2:N
        A(i,j) = alp^(j-1)*2*cos((0.05495*(i-50)+0.33755)*(j-1))*abs(1/(1-alp^N*exp(-1i*(0.05495*(i-50)+0.33755)*N)));
    end
end


cvx_begin
    variable d(1,N)
    variable delta_stop nonnegative
    minimize( delta_stop )
    subject to
        for i = 2:50
            0.5 <= sum(A(i,:).*d) <= 1.5;
        end
        sum(A(1,:).*d) == 1;
        for i = 51:100
            0 <= sum(A(i,:).*d) <= delta_stop;
        end
cvx_end

% A = zeros(100,N);
% 
% for i = 1:50
%     A(i,1) = abs(1/(1-alp^N*exp(-1i*0.00157*i*N)));
%     for j = 2:N
%         A(i,j) = alp^(j-1)*2*cos(0.00157*i*(j-1))*abs(1/(1-alp^N*exp(-1i*0.00157*i*N)));
%     end
% end
% 
% for i = 51:100
%     A(i,1) = abs(1/(1-alp^N*exp(-1i*(0.05809*(i-50)+0.17741)*N)));
%     for j = 2:N
%         A(i,j) = alp^(j-1)*2*cos((0.05809*(i-50)+0.17741)*(j-1))*abs(1/(1-alp^N*exp(-1i*(0.05809*(i-50)+0.17741)*N)));
%     end
% end
% 
% 
% cvx_begin
%     variable d1(1,N)
%     variable delta_stop nonnegative
%     minimize( delta_stop )
%     subject to
%         for i = 2:50
%             0.5 <= sum(A(i,:).*d1) <= 1.5;
%         end
%         sum(A(1,:).*d1) == 1;
%         for i = 51:100
%             0 <= sum(A(i,:).*d1) <= delta_stop;
%         end
% cvx_end

% d = wextend('1D','zpd',d,20,'l');
% u = zeros(40,1);
% u(21) = 1;
% sysd = tfest(u,d,10,10,'Ts',1);
% sys = arx(u,d,[20 20 0]);
% sysd
% sys.Structure
d
for i = 1:N
    d(1,i) = d(1,i)*alp^(i-1);
end
d = transpose(d);
d1 = flipud(d);
d1(end) = [];
d2 = [d1; d];
d2 = transpose(d2);
pad = zeros(1,N-1);
denom = [1 pad];
tf(d2, denom)
spectralfact(tf(d2, denom))


% y = fft(d);
% 
n1 = 1000;
f = (0:n1-1)*(1.57/n1);     % frequency range
% power = abs(y).^2/n;    % power of the DFT
% 

hold on

y_array2 = arrayfun(@(f) (20*log10(fin_in(f))), f,'UniformOutput',false);
y_array2 = cell2mat(y_array2);
plot(f/0.157, (y_array2));

y_array = arrayfun(@(f) (20*log10(sqrt(fin_out(f,d1)))), f,'UniformOutput',false);
y_array = cell2mat(y_array);
plot(f/0.157, (y_array));

y_array = arrayfun(@(f) (20*log10(sqrt(fin_out(f,d)))), f,'UniformOutput',false);
y_array = cell2mat(y_array);
plot(f/0.157, (y_array));
xlabel('Frequency (xFs)')
ylabel('Power (dB)')

legend('Single pole filter','FA with Fstop = 1.5Fs','FA with Fstop = 2.5Fs')

hold off

function outp = fin_out(w,d)
    alp = 0.9;
    N = 40;
    B = zeros(1,N);
    B(1,1) = abs(1/(1-alp^N*exp(-1i*w*N)));
    for j = 2:N
        B(1,j) = alp^(j-1)*2*cos(w*(j-1))*abs(1/(1-alp^N*exp(-1i*w*N)));
    end
    outp = abs(sum(B(1,:).*d));
end

function outp = fin_in(w)
    alp = 0.9;
    outp = abs(1/(1-alp*exp(-1i*w)))/10;
end


