f = [0 0.0001 0.0002 0.05 0.15 1];
a = [1 1 1 0.5 0 0];
w = [10000000 1 1000000];
b = firpm(20, f, a, w);
n1 = 512;

alp = 0.9;
N = 40;

n1 = 1000;
f = (0:n1-1)*(1.57/n1);
A = zeros(100,N);

for i = 1:50
    A(i,1) = abs(1/(1-alp^N*exp(-1i*0.00157*i*N)));
    for j = 2:N
        A(i,j) = alp^(j-1)*2*cos(0.00157*i*(j-1))*abs(1/(1-alp^N*exp(-1i*0.00157*i*N)));
    end
end

for i = 51:100
    A(i,1) = abs(1/(1-alp^N*exp(-1i*(0.059275*(i-51)+0.2355)*N)));
    for j = 2:N
        A(i,j) = alp^(j-1)*2*cos((0.059275*(i-51)+0.2355)*(j-1))*abs(1/(1-alp^N*exp(-1i*(0.059275*(i-51)+0.2355)*N)));
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
d

hold on
y_array = arrayfun(@(f) (20*log10(sqrt(fin_out(f,d)))), f,'UniformOutput',false);
y_array = cell2mat(y_array);
length(f)
length(y_array)
plot(f*2/0.157, y_array);

[h,w] = freqz(b,1,512);
plot(w/pi*20,20*log10(abs(h)))
x = zeros(512);
for i = 1:  
    x(i) = 1;
end
for i = 26:77
    x(i) = 1-(i-26)/51;
end
plot(w/pi*20,x)

f = [0 0.0001 0.0002 0.05 0.25 1];
a = [1 1 1 0.5 0 0];
w = [100000 1 100000];
b = firpm(20, f, a, w);

% [h,w] = freqz(b,1,512);
% plot(w/pi*20,20*log10(abs(h)))

xlabel('Frequency (xFs)')
ylabel('Amplitude')

legend('Using convex optimisation', 'Using Parks-McClellan','FA with Fstop = 2.5Fs')

hold off
for i = 1:length(b)
    x(i) = b(i) / exp(-1*(i-1)/10);
end

% [h,w] = freqz(x,1,512);
% plot(w/pi,log10(abs(h)))
% x = fliplr(x)
x;
% y = fft(x);                               % Compute DFT of x
% m = abs(y);                               % Magnitude
% p = unwrap(angle(y));                     % Phase
% m/m(1);
% p*180/pi;
% f = (0:length(y)-1)*100/length(y);      % Frequency vector

% subplot(2,1,1)
% plot(f,m)
% title('Magnitude')
% ax = gca;
% ax.XTick = [15 40 60 85];
% 
% subplot(2,1,2)
% plot(f,p*180/pi)
% title('Phase')
% ax = gca;
% ax.XTick = [15 40 60 85];
function outp = fin_in(w)
    alp = 0.8;
    outp = abs(1/(1-alp*exp(-1i*w)))/5;
end

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

