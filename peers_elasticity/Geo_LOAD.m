function [c,e,m] = Geo_LOAD(Type)

switch case
    case 1
        load coordinates.txt; load element.txt; load mc_corr.txt
        c = coordinates; e = element; m = mc_corr;
    case 2
        c=0; e=0; m=0;

return