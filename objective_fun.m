function [fun, info] = objective_fun(func)
%OBJECTIVE_FUN Return a benchmark objective from the QPSO paper.
%   FUN = OBJECTIVE_FUN('f1') returns a function handle.
%   [FUN, INFO] = OBJECTIVE_FUN('f1') also returns benchmark metadata.

defs = benchmark_definitions();
ids = {defs.Id};
names = {defs.Name};
idx = find(strcmpi(func, ids) | strcmpi(func, names), 1);

if isempty(idx)
    error('objective_fun:UnknownFunction', ...
        'Unknown benchmark "%s". Use f1 ... f24.', func);
end

info = defs(idx);
fun = @(x) evaluate_benchmark(x, info.Id);
end

function y = evaluate_benchmark(x, funcId)
x = x(:)';
n = numel(x);

switch lower(funcId)
    case 'f1' % Sphere
        y = sum(x.^2);

    case 'f2' % Schwefel's No. 2.22
        y = sum(abs(x)) + prod(abs(x));

    case 'f3' % Schwefel's No. 1.2
        y = sum(cumsum(x).^2);

    case 'f4' % Schwefel's No. 2.21
        y = max(abs(x));

    case 'f5' % Rosenbrock
        y = sum(100 .* (x(2:end) - x(1:end-1).^2).^2 + ...
            (x(1:end-1) - 1).^2);

    case 'f6' % Step
        y = sum(floor(x + 0.5).^2);

    case 'f7' % Quartic with noise
        y = sum((1:n) .* (x.^4)) + rand();

    case 'f8' % X.S. Yang No. 7
        epsilon = rand(1, n);
        y = sum(epsilon .* abs(x) .^ (1 ./ (1:n)));

    case 'f9' % Schwefel's No. 2.26
        y = -sum(x .* sin(sqrt(abs(x))));

    case 'f10' % Rastrigin
        y = sum(x.^2 - 10 .* cos(2 .* pi .* x) + 10);

    case 'f11' % Ackley
        y = -20 .* exp(-0.2 .* sqrt(sum(x.^2) ./ n)) - ...
            exp(sum(cos(2 .* pi .* x)) ./ n) + 20 + exp(1);

    case 'f12' % Griewank
        y = sum(x.^2) ./ 4000 - prod(cos(x ./ sqrt(1:n))) + 1;

    case 'f13' % Penalized No. 1
        z = 1 + (x + 1) ./ 4;
        y = (pi ./ n) .* (10 .* sin(pi .* z(1)).^2 + ...
            sum((z(1:end-1) - 1).^2 .* ...
            (1 + 10 .* sin(pi .* z(2:end)).^2)) + ...
            (z(end) - 1).^2) + sum(penalty_u(x, 10, 100, 4));

    case 'f14' % Penalized No. 2
        y = 0.1 .* (sin(3 .* pi .* x(1)).^2 + ...
            sum((x(1:end-1) - 1).^2 .* ...
            (1 + sin(3 .* pi .* x(2:end)).^2)) + ...
            (x(end) - 1).^2 .* (1 + sin(2 .* pi .* x(end)).^2)) + ...
            sum(penalty_u(x, 5, 100, 4));

    case 'f15' % Alpine No. 1
        y = sum(abs(x .* sin(x) + 0.1 .* x));

    case 'f16' % Shekel's Foxholes
        gridValues = [-32 -16 0 16 32];
        [a1, a2] = meshgrid(gridValues, gridValues);
        A = [a1(:)'; a2(:)'];
        xx = repmat(x(1:2)', 1, 25);
        denominator = (1:25) + sum((xx - A).^6, 1);
        y = (1 ./ 500 + sum(1 ./ denominator)) .^ (-1);

    case 'f17' % Kowalik
        a = [0.1957 0.1947 0.1735 0.1600 0.0844 0.0627 ...
             0.0456 0.0342 0.0323 0.0235 0.0246];
        b = [0.25 0.5 1 2 4 6 8 10 12 14 16];
        y = sum((a - (x(1) .* (b.^2 + x(2) .* b)) ./ ...
            (b.^2 + x(3) .* b + x(4))).^2);

    case 'f18' % Six-Hump Camel-Back
        x1 = x(1); x2 = x(2);
        y = 4 .* x1.^2 - 2.1 .* x1.^4 + (1 ./ 3) .* x1.^6 + ...
            x1 .* x2 - 4 .* x2.^2 + 4 .* x2.^4;

    case 'f19' % Branin's RCOS No. 1
        x1 = x(1); x2 = x(2);
        y = (x2 - (5.1 .* x1.^2) ./ (4 .* pi.^2) + ...
            (5 .* x1) ./ pi - 6).^2 + ...
            10 .* (1 - 1 ./ (8 .* pi)) .* cos(x1) + 10;

    case 'f20' % Goldstein-Price
        x1 = x(1); x2 = x(2);
        y = (1 + (x1 + x2 + 1).^2 .* ...
            (19 - 14 .* x1 + 3 .* x1.^2 - 14 .* x2 + ...
            6 .* x1 .* x2 + 3 .* x2.^2)) .* ...
            (30 + (2 .* x1 - 3 .* x2).^2 .* ...
            (18 - 32 .* x1 + 12 .* x1.^2 + 48 .* x2 - ...
            36 .* x1 .* x2 + 27 .* x2.^2));

    case 'f21' % Hartman's No. 3
        alpha = [1 1.2 3 3.2];
        A = [3 10 30; 0.1 10 35; 3 10 30; 0.1 10 35];
        P = 1e-4 .* [3689 1170 2673; 4699 4387 7470; ...
            1091 8732 5547; 381 5743 8828];
        y = hartman_value(x, alpha, A, P);

    case 'f22' % Hartman's No. 6
        alpha = [1 1.2 3 3.2];
        A = [10 3 17 3.5 1.7 8; ...
            0.05 10 17 0.1 8 14; ...
            3 3.5 1.7 10 17 8; ...
            17 8 0.05 10 0.1 14];
        P = 1e-4 .* [1312 1696 5569 124 8283 5886; ...
            2329 4135 8307 3736 1004 9991; ...
            2348 1451 3522 2883 3047 6650; ...
            4047 8828 8732 5743 1091 381];
        y = hartman_value(x, alpha, A, P);

    case 'f23' % Shekel No. 7
        y = shekel_value(x, 7);

    case 'f24' % Shekel No. 10
        y = shekel_value(x, 10);

    otherwise
        error('objective_fun:UnknownFunction', ...
            'Unknown benchmark "%s".', funcId);
end
end

function y = penalty_u(x, a, k, m)
y = zeros(size(x));
y(x > a) = k .* (x(x > a) - a).^m;
y(x < -a) = k .* (-x(x < -a) - a).^m;
end

function y = hartman_value(x, alpha, A, P)
outer = 0;
for i = 1:numel(alpha)
    inner = sum(A(i, :) .* (x - P(i, :)).^2);
    outer = outer + alpha(i) .* exp(-inner);
end
y = -outer;
end

function y = shekel_value(x, m)
A = [4 4 4 4; ...
     1 1 1 1; ...
     8 8 8 8; ...
     6 6 6 6; ...
     3 7 3 7; ...
     2 9 2 9; ...
     5 5 3 3; ...
     8 1 8 1; ...
     6 2 6 2; ...
     7 3.6 7 3.6];
c = [0.1 0.2 0.2 0.4 0.4 0.6 0.3 0.7 0.5 0.5];

y = 0;
for i = 1:m
    y = y + 1 ./ (sum((x - A(i, :)).^2) + c(i));
end
y = -y;
end
