function [ output_args ] = outputDEMVTK( file_string,n_dem,q,r,v,material_parameters )
fp = fopen(file_string,'w');

    fprintf(fp,'<?xml version="1.0"?>\n');
    fprintf(fp,'<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">\n');

    fprintf(fp,' <UnstructuredGrid>\n');
    fprintf(fp,'  <Piece NumberOfPoints="%d" NumberOfCells="0">\n',n_dem);
    fprintf(fp,'   <Points>\n');
    fprintf(fp,'    <DataArray name="Position" type="Float32" NumberOfComponents="3" format="ascii">\n');
    for ii = 1:n_dem
%         string_write = strcat('    ',num2str(q(ii,:)));
%         fprintf(fp,strcat(string_write,'\n'));
        fprintf(fp,'    %.16g %.16g %.16g\n',q(ii,1),q(ii,2),q(ii,3));
    end
    fprintf(fp,"    </DataArray>\n");
    fprintf(fp,"   </Points>\n");
    fprintf(fp,'   <PointData Vectors="vectors">\n');
    
%   Prints out velocity components
    fprintf(fp,'    <DataArray type="Float32" Name="Velocity" NumberOfComponents="3" format="ascii">\n');
    for ii = 1:n_dem
%         string_write = strcat('    ',num2str(v(ii,:)));
%         fprintf(fp,strcat(string_write,'\n'));
        fprintf(fp,'    %.16g %.16g %.16g\n',v(ii,1),v(ii,2),v(ii,3));
    end
    fprintf(fp,"    </DataArray>\n");
    
    %   Prints out radius
    fprintf(fp,'    <DataArray type="Float32" Name="Radius" NumberOfComponents="1" format="ascii">\n');
    for ii = 1:n_dem
%         string_write = strcat('    ',num2str(r(ii)));
%         fprintf(fp,strcat(string_write,'\n'));
        fprintf(fp,'    %.16g\n',r(ii));
    end
    fprintf(fp,"    </DataArray>\n");
    
%   Prints out coordinates so don't need to use calculator in paraview
    fprintf(fp,'    <DataArray type="Float32" Name="Coordinates" NumberOfComponents="3" format="ascii">\n');
    for ii = 1:n_dem
%         string_write = strcat('    ',num2str(q(ii,:)));
%         fprintf(fp,strcat(string_write,'\n'));
        fprintf(fp,'    %.16g %.16g %.16g\n',q(ii,1),q(ii,2),q(ii,3));
    end
    fprintf(fp,"    </DataArray>\n");

    fprintf(fp,'   </PointData>\n');
    fprintf(fp,'   <Cells>\n');
    fprintf(fp,'    <DataArray type="Int32" Name="connectivity" format="ascii">\n');
    fprintf(fp,'    </DataArray>\n');
    fprintf(fp,'    <DataArray type="Int32" Name="offsets" format="ascii">\n');
    fprintf(fp,'    </DataArray>\n');
    fprintf(fp,'    <DataArray type="Int32" Name="types" format="ascii">\n');
    fprintf(fp,'    </DataArray>\n');
    fprintf(fp,'   </Cells>\n');
    fprintf(fp,'  </Piece>\n');
    fprintf(fp,' </UnstructuredGrid>\n');
    fprintf(fp,' <Field>\n');
    fprintf(fp,' <!--k_n(normal spring coeff) k_t gamma_n(normal viscous damp) gamma_t rho(density) mu(friction coefficient)-->\n');
    fprintf(fp,'  <DataArray type="Float32" Name="Material_parameters" NumberOfComponents="6" format="ascii">\n');
    fprintf(fp,'  %.16g %.16g %.16g %.16g %.16g %.16g\n',material_parameters);
    fprintf(fp,'  </DataArray>\n');
    fprintf(fp,' </Field>\n');
    fprintf(fp,'</VTKFile>\n');





fclose(fp);

end

