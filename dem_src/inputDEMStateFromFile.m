function [] = inputDEMStateFromFile(state_filename,contact_filename)
state_fp = fopen(state_filename);
tline = fgetl(state_fp);
num_points = 0;
q_read_flag = 0;
ctr = 1;
while ischar(tline)
    disp(tline)
    
    if(q_read_flag == 1)
        
    end
    
    if(contains(tline,'NumberOfPoints'))
        tline = split(tline,'"');
        num_points = str2num(tline{2});
    end
    
    if(contains(tline,'Position'))
        q_read_flag = 1;
    end
    tline = fgetl(state_fp);
end

