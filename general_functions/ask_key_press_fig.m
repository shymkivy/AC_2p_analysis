function [answer] = ask_key_press_fig(accept, reject)

accept_num = int8(accept);
reject_num = int8(reject);

temp_prompt = 0;
temp_finish = 0;
while temp_finish == 0
    k = waitforbuttonpress;  % check if key or mouse 
    if k == 1               % if key was pressed
        temp_value = double(get(gcf,'CurrentCharacter'));
        if isempty(temp_value)
            fprintf('Use keyboard\n');
        elseif temp_value == accept_num || temp_value == (accept_num - 32)        % y = 121 Y = 89
            temp_finish = 1;
            answer = 1; % answer yes
        elseif temp_value == reject_num || temp_value == (reject_num - 32)    % n = 110 N = 78
            answer = 0; % answer no
            temp_finish = 1;
        else
            temp_prompt = temp_prompt + 1;
            if mod(temp_prompt,5)==0
                fprintf('COME ON!!!!\n');
            elseif mod(temp_prompt,7)==0
                fprintf('Stop messing around!\n');
            elseif temp_prompt == 16
                fprintf('Do it one more time and I will tell Rafa\n');
            elseif temp_prompt == 17
                fprintf('Email will be sent out soon...\n');
            else
                fprintf('It is a Yes or No question... [Y/N]\n');
            end
        end
    end
end