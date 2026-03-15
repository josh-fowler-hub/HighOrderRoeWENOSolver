function filePre = makeFilePrefix(Schemes, test_num, n, cfl, riemann_solver, CD_Term_Order)
% Construct a prefix for output filenames based on configuration
if strcmp(riemann_solver, 'Roe')
    solver = 'Roe';
elseif strcmp(riemann_solver, 'Roe-Pike')
    solver = 'RoePike';
elseif strcmp(riemann_solver, 'HLL')
    solver = 'HLL';
elseif strcmp(riemann_solver, 'HLLC')
    solver = 'HLLC';
elseif strcmp(riemann_solver, 'Osher')
    solver = 'Osher';
else
    solver = '';
end

Scheme_Path = [];
for i = 1:length(Schemes)
    if i ~= length(Schemes)
        curr = [char(Schemes(i)) 'v'];
    else
        curr = [char(Schemes(i))];
    end
    Scheme_Path = [Scheme_Path curr];
end

if CD_Term_Order == 1
    order = '';
else
    order = ['CD' num2str(CD_Term_Order)];
end

file_pre = ['Test' num2str(test_num) order Scheme_Path 'CFL'...
            num2str(cfl) 'Nodes' num2str(n) solver];
end
