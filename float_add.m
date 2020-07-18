function [base, rem] = float_add(base,rem, add)

bef = base;
rem = rem + add;
base = base + rem;
delta = base - bef;
rem = rem - delta;

end

