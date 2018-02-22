function [linIndsOfEdges] = GetLinearIndicesOfImageEdge(sizeOfImage)

edgePixels1 = [1:sizeOfImage(1)-1; ones(1,sizeOfImage(1)-1)];
edgePixels2 =[ sizeOfImage(1)*ones(1,sizeOfImage(2)-1); 1:sizeOfImage(2)-1 ];   
edgePixels3 = [1:sizeOfImage(1); sizeOfImage(2)*ones(1,sizeOfImage(1))];
edgePixels4 =[ ones(1,sizeOfImage(2)-1); 1:sizeOfImage(2)-1 ];   
edgePixels = [edgePixels1,edgePixels2,edgePixels3,edgePixels4];
linIndsOfEdges = sub2ind(sizeOfImage,edgePixels(1,:),edgePixels(2,:))'; 
end