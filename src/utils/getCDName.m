function name = getCDName(n)
% Convert central differencing order to string label
switch n
    case 1
        name = 'Roe';
    case 2
        name = 'CD2';
    case 4
        name = 'CD4';
    case 6
        name = 'CD6';
    otherwise
        name = ['CD' num2str(n)];
end
end
