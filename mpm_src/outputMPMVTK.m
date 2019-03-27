function [ output_args ] = outputMPMVTK( file_string,n_dem,q,v,sigma,gamma_bar_dot_p )
sigma_voigt(:,1) = sigma(:,1,1);
sigma_voigt(:,2) = sigma(:,2,2);
sigma_voigt(:,3) = sigma(:,3,3);
sigma_voigt(:,4) = sigma(:,2,3);
sigma_voigt(:,5) = sigma(:,3,1);
sigma_voigt(:,6) = sigma(:,1,2);

fp = fopen(file_string,'w');

    fprintf(fp,'<?xml version="1.0"?>\n');
    fprintf(fp,'<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">\n');

    fprintf(fp,' <UnstructuredGrid>\n');
    fprintf(fp,'  <Piece NumberOfPoints="%d" NumberOfCells="0">\n',n_dem);
    fprintf(fp,'   <Points>\n');
    fprintf(fp,'    <DataArray name="Position" type="Float32" NumberOfComponents="3" format="ascii">\n');
    for ii = 1:n_dem
        string_write = strcat('    ',num2str(q(ii,:)));
        fprintf(fp,strcat(string_write,'\n'));
    end
    fprintf(fp,"    </DataArray>\n");
    fprintf(fp,"   </Points>\n");
    fprintf(fp,'   <PointData Vectors="vectors">\n');
% 
%     // Prints out sigmaxx, sigmayy, sigmaxy, in that order
%     fprintf(fp,"    <DataArray type="Float32" Name="Cauchy Stress" NumberOfComponents="3" format="ascii">n");
%     for (int ii = 0); ii < bodies.m_numBodies); ii++)
%         for (int p = 0); p < bodies.m_pointObjects[ii]->m_nPoints); p++)
%             fprintf(fp,"     "<<bodies.m_pointObjects[ii]->m_pointStressNew[p](0,0)<<" "<<bodies.m_pointObjects[ii]->m_pointStressNew[p](0,1)<<" "<<bodies.m_pointObjects[ii]->m_pointStressNew[p](1,1)<<"n");
%     fprintf(fp,"    </DataArray>n");
% 
%     // Prints out velocity components
    fprintf(fp,'    <DataArray type="Float32" Name="Velocity" NumberOfComponents="3" format="ascii">\n');
    for ii = 1:n_dem
%         v(ii,:)
        string_write = strcat('    ',num2str(v(ii,:)));
        fprintf(fp,strcat(string_write,'\n'));
    end
    fprintf(fp,"    </DataArray>\n");
% 
%     // Prints out coordinates so don't need to use calculator in paraview
    fprintf(fp,'    <DataArray type="Float32" Name="Coordinates" NumberOfComponents="3" format="ascii">\n');
    for ii = 1:n_dem
        string_write = strcat('    ',num2str(q(ii,:)));
        fprintf(fp,strcat(string_write,'\n'));
    end
    fprintf(fp,"    </DataArray>\n");
    
    %     // Prints out sigma in voigt notation
    fprintf(fp,'    <DataArray type="Float32" Name="Stress" NumberOfComponents="6" format="ascii">\n');
    for ii = 1:n_dem
        string_write = strcat('    ',num2str(sigma_voigt(ii,:)));
        fprintf(fp,strcat(string_write,'\n'));
    end
    fprintf(fp,"    </DataArray>\n");

    %     // Prints out gamma bar dot p
    fprintf(fp,'    <DataArray type="Float32" Name="Gamma_bar_dot_p" NumberOfComponents="1" format="ascii">\n');
    for ii = 1:n_dem
        string_write = strcat('    ',num2str(gamma_bar_dot_p(ii)));
        fprintf(fp,strcat(string_write,'\n'));
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
    fprintf(fp,'</VTKFile>\n');





fclose(fp);

end

