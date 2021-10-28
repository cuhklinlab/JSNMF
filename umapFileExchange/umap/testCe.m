function C=testCe(reset)
if nargin<1 || reset
close all force;
try
BasicMap.Global(true);
CytoGate.Get(true)
catch ex
end
end
C=ColorsEditor([],[],[],true,false);
end