function output = TakeDerivative(T, StepSize)
T_temp = diff(T )./StepSize;
output  =   [ T_temp  T_temp(end)];
end