function [np] = Order_of_Accuracy(Lkp1, Lk, r)

np = log(Lkp1/Lk)/log(r);

end