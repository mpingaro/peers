function [NX,NY,TYPE] = inputcook()
fprintf ('==============================================================\n')
fprintf ('                C O O K    C A N T I L E V E R                \n')
fprintf ('==============================================================\n\n')
NX = input('Number of partitions in x-direction: ');
NY = input('Number of partitions in y-direction: ');
fprintf ('\nSelect the type of mesh: \n')
fprintf ('  1) Triangular mesh\n')
fprintf ('  2) Quadrilateral mesh\n')
TYPE = input('Insert the number relative to the choose (1 or 2): ');
end

