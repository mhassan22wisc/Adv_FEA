%% PlotVtk from Text

ii = 0;
while true
    
    NodalCoords = readcell( strcat ('NodalCoords',num2str(ii)));
    %if( isempty( NodalCoords )  == 0 )
    %    break;
    %end
    ElemConnect = readcell( strcat ('ElemConnect',num2str(ii)));
    varargin = readcell( strcat ('varargin',num2str(ii)));
    OutputVTKFileName = strcat ('vtkfile',num2str(ii));
    fem_to_vtk_Vector (OutputVTKFileName, NodalCoords, ElemConnect, varargin);
    ii = ii+1;

end