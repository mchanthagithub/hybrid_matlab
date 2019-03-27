function [ state_file_string, contact_file_string ] = makeDEMOutputFileStrings( num_saves_total,save_ctr, state_pre_string, contact_pre_string)
num_digits = numel(num2str(num_saves_total));
zero_string = '';
for ii = 1:(num_digits - numel(num2str(save_ctr)))
    zero_string = strcat(zero_string,'0');
end
temp_string = strcat(state_pre_string,zero_string);
temp_string2 = strcat(temp_string,num2str(save_ctr));
state_file_string = strcat(temp_string2,'.vtu');

temp_string_contacts = strcat(contact_pre_string,zero_string);
temp_string_contacts2 = strcat(temp_string_contacts,num2str(save_ctr));
contact_file_string = strcat(temp_string_contacts2,'.txt');

end

