function defs = benchmark_definitions()
%BENCHMARK_DEFINITIONS Metadata for the 24 functions in Tables 7-9.

template = struct('Id', '', 'Name', '', 'Type', '', 'Dimension', [], ...
    'LowerBound', [], 'UpperBound', [], 'FMin', []);
defs = repmat(template, 24, 1);

defs(1)  = def('f1',  'Sphere',                 'Unimodal',          30, -100, 100, 0);
defs(2)  = def('f2',  'Schwefel No. 2.22',      'Unimodal',          30, -10,  10,  0);
defs(3)  = def('f3',  'Schwefel No. 1.2',       'Unimodal',          30, -100, 100, 0);
defs(4)  = def('f4',  'Schwefel No. 2.21',      'Unimodal',          30, -100, 100, 0);
defs(5)  = def('f5',  'Rosenbrock',             'Unimodal',          30, -30,  30,  0);
defs(6)  = def('f6',  'Step',                   'Unimodal',          30, -100, 100, 0);
defs(7)  = def('f7',  'Quartic',                'Unimodal',          30, -1.28, 1.28, 0);
defs(8)  = def('f8',  'X.S. Yang No. 7',        'Unimodal',          30, -5,   5,   0);

defs(9)  = def('f9',  'Schwefel No. 2.26',      'Multimodal',        30, -500, 500, -12569.5);
defs(10) = def('f10', 'Rastrigin',              'Multimodal',        30, -5.12, 5.12, 0);
defs(11) = def('f11', 'Ackley',                 'Multimodal',        30, -32,  32,  0);
defs(12) = def('f12', 'Griewank',               'Multimodal',        30, -600, 600, 0);
defs(13) = def('f13', 'Penalized No. 1',        'Multimodal',        30, -50,  50,  0);
defs(14) = def('f14', 'Penalized No. 2',        'Multimodal',        30, -50,  50,  0);
defs(15) = def('f15', 'Alpine No. 1',           'Multimodal',        30, -10,  10,  0);

defs(16) = def('f16', 'Shekel Foxholes',        'FixedMultimodal',    2, -65,  65,  1);
defs(17) = def('f17', 'Kowalik',                'FixedMultimodal',    4, -5,   5,   0.0003075);
defs(18) = def('f18', 'Six-Hump Camel-Back',    'FixedMultimodal',    2, -5,   5,  -1.03163);
defs(19) = def('f19', 'Branin RCOS No. 1',      'FixedMultimodal',    2, -5,   15,  0.39800);
defs(20) = def('f20', 'Goldstein-Price',        'FixedMultimodal',    2, -2,   2,   3);
defs(21) = def('f21', 'Hartman No. 3',          'FixedMultimodal',    3, 0,    1,  -3.8628);
defs(22) = def('f22', 'Hartman No. 6',          'FixedMultimodal',    6, 0,    1,  -3.32);
defs(23) = def('f23', 'Shekel No. 7',           'FixedMultimodal',    4, 0,    10, -10.4028);
defs(24) = def('f24', 'Shekel No. 10',          'FixedMultimodal',    4, 0,    10, -10.5363);
end

function s = def(id, name, typeName, dimension, lb, ub, fmin)
s = struct('Id', id, 'Name', name, 'Type', typeName, ...
    'Dimension', dimension, 'LowerBound', lb, 'UpperBound', ub, ...
    'FMin', fmin);
end
