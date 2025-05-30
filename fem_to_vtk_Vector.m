function fem_to_vtk_Vector (OutputVTKFileName, NodalCoords, ElemConnect, varargin)

%*****************************************************************************80
%
% FEM_TO_VTK_Vector converts an FEM model into a VTK ASCII file.
%
%  Usage:
%
%    fem_to_vtk_Vector (OutputVTKFileName, NodalCoords, ElemConnect, VectorValues(optional) )
%
%    * OutputVTKFileName: the prefix of the VTK unstructured mesh format file OutputVTKFileName.vtu which will be created by FEM_TO_VTK
%    * NodalCoords: the node coordinates array (Size:numNodes X num_dims). num_dims is 2/3.
%    * ElemConnect: the element connectivity array (Size:num_Elements X numNodesPerElement). .
%    * VectorValues[optional]:  the vector values defined at each node (Size:numNodes X num_components of the vector). 
%
% Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
% Author:
%   Siva Shankar Rudraraju @ comp-phys.engin.umich.edu
% Created:
%    6 November 2010
% Modified:
%    6 November 2010
  
  fprintf ('Reading data arrays\n' );
%  Get the file prefix.
%
  if (nargin < 3)
    fprintf ('Atleast 3 arguments (OutputVTKFileName, NodalCoords, ElemConnect) required.\n' ); return;
  elseif (~ischar(OutputVTKFileName))
    fprintf ('First argument (OutputVTKFileName) must be a string corresponding to the output filename desired.\n' ); return;  
  elseif ~(isnumeric(NodalCoords) && isnumeric(ElemConnect))
    fprintf ('Second and Third arguments (NodalCoords, ElemConnect) must be a numeric arrays generated by the meshing programs.\n' ); return;    
  end
  
  if (nargin==3)
      values=zeros(length(NodalCoords(:,1)),3); value_num=3; 
  elseif (nargin>3)    
      if (~isnumeric(varargin{1})), 
          fprintf ('Fourth argument (Values) must be a numeric array representing the nodal field values.\n' ); return; 
      end   
      values=transpose(varargin{1}); value_num=length(values(:,1)); 
      if (value_num==1)
          fprintf ('The vector field array should have atleast two components at each node. For scalar fields use fem_to_vtk.m\n' ); return; 
      elseif(value_num==2)
          values=[values(1,:);values(2,:);values(1,:).*0.0]; value_num=3; 
      elseif(value_num>3)    
          fprintf ('The vector field array should have atmost 3 components at each node.\n' ); return; 
      end
  end
%
%  Create the filenames.
%
  vtk_filename = strcat (OutputVTKFileName, '.vtu');
%
%  Read the FEM data.
%
  node_coord=NodalCoords';
  element_node=ElemConnect';
  dim_num=length(node_coord(:,1)); node_num=length(node_coord(1,:));
  element_num=length(element_node(1,:)); element_order=length(element_node(:,1));
%
%  Write the VTK data.
%
  fprintf ('Writing VTK file\n' );
  vtk_write ( vtk_filename, dim_num, node_num, element_num, ...
    element_order, value_num, node_coord, element_node, values );

%
%  Terminate.
%
  fprintf ('\nFEM_TO_VTK\n' );
  fprintf ('  Normal end of execution.\n' );
  return
end
%
%**************************************************************************
%
function vtk_write ( vtk_filename, dim_num, node_num, element_num, ...
  element_order, value_num, node_coord, element_node, values )

%Identifying vtk cell type
element_type=0;
if ((dim_num==2) && (element_order==3))
    element_type=5;
elseif ((dim_num==2) && (element_order==4))    
    element_type=9;
elseif ((dim_num==3) && (element_order==4))
    element_type=10;
elseif ((dim_num==3) && (element_order==8))    
    element_type=12;
else
    ME = MException('ElementType:Error', ...
       'Element type with %i nodes per element in %iD not supported by FEM_TO_VTK\n', element_order, dim_num);
    throw(ME)
end
vtkFile = fopen(vtk_filename, 'wt' );
if (dim_num==2)
    dim_num=3; node_coord(3,:)=zeros(1,node_num);
end
%
%  Write the vtk header.
%
fprintf (vtkFile, '<?xml version="1.0"?>\n');
fprintf (vtkFile, '<VTKFile type="UnstructuredGrid"  version="0.1"  >\n');
fprintf (vtkFile, '<UnstructuredGrid>\n');
%Writing nodal coordinates array
fprintf (vtkFile, '<Piece  NumberOfPoints="%i" NumberOfCells="%i">\n', node_num, element_num);
fprintf (vtkFile, '<Points>\n');
fprintf (vtkFile, '<DataArray  type="Float32"  NumberOfComponents="%i"  format="ascii">\n', dim_num);
for node=1:node_num
    for dim=1:dim_num
        fprintf (vtkFile, '%f ', node_coord(dim,node));
    end
    fprintf (vtkFile, '\n');
end
fprintf (vtkFile, '</DataArray>\n');
fprintf (vtkFile, '</Points>\n');
%Writing element connectivity array
fprintf (vtkFile, '<Cells>\n');
fprintf (vtkFile, '<DataArray  type="UInt32"  Name="connectivity"  format="ascii">\n');
for elem=1:element_num
    for node=1:element_order
        fprintf (vtkFile, '%i ', element_node(node, elem)-1);
    end
    fprintf (vtkFile, '\n');
end
fprintf (vtkFile, '</DataArray>\n');
%Writing element offsets vector(required in VTK format)
offsets=0;
fprintf (vtkFile, '<DataArray  type="UInt32"  Name="offsets"  format="ascii">\n');
for elem=1:element_num
    offsets=offsets+length(element_node(:, elem));
    fprintf (vtkFile, '%i\n', offsets);
end
fprintf (vtkFile, '</DataArray>\n');
%Writing element types vector(required in VTK format)
fprintf (vtkFile, '<DataArray  type="UInt8"  Name="types"  format="ascii">\n');
for elem=1:element_num
    fprintf (vtkFile, '%i\n', element_type);
end
fprintf (vtkFile, '</DataArray>\n');
fprintf (vtkFile, '</Cells>\n');
%Writing field values (if any) array
fprintf (vtkFile, '<PointData  Scalars="u">\n');
fprintf (vtkFile, '<DataArray  type="Float32"  Name="VectorField%i"  NumberOfComponents="%u" format="ascii">\n', 1, value_num);
for node=1:node_num
    for field=1:value_num
        fprintf (vtkFile, '%f ', values(field, node));
    end
    fprintf (vtkFile, '\n');
end
fprintf (vtkFile, '</DataArray> \n');
%VTK file footer
fprintf (vtkFile, '</PointData> \n');
fprintf (vtkFile, '</Piece> \n');
fprintf (vtkFile, '</UnstructuredGrid> \n');
fprintf (vtkFile, '</VTKFile> \n');
%closing VTK file IO
fprintf('Data written in VTK unstructured format to "%s".\n', vtk_filename);
fclose(vtkFile);
end




