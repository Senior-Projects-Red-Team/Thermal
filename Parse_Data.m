

function[datafile] = Parse_Data(filename)

fid = fopen(filename);
lines = textscan(fid,'%s','Delimiter','\n');
fclose(fid);

lines = lines{1};

% Preallocate (estimate half the lines are data rows)
n = floor(length(lines)/2);

TC0 = zeros(n,1);
TC1 = zeros(n,1);
TC2 = zeros(n,1);
TC3 = zeros(n,1);
TMP1 = zeros(n,1);
TMP2 = zeros(n,1);

Pin9  = zeros(n,1);
Pin10 = zeros(n,1);
Pin11 = zeros(n,1);
Pin6  = zeros(n,1);

row = 0;

for i = 1:length(lines)
    
    % ----- Temperature Line -----
    if contains(lines{i},'TC 0:')
        row = row + 1;
        
        nums = regexp(lines{i},'[-+]?\d*\.?\d+','match');
        nums = str2double(nums);
        
        TC0(row) = nums(1);
        TC1(row) = nums(2);
        TC2(row) = nums(3);
        TC3(row) = nums(4);
        TMP1(row)= nums(5);
        TMP2(row)= nums(6);
    end
    
    % ----- Relay Line -----
    if contains(lines{i},'Relay Duty Cycles')
        nums = regexp(lines{i},'\d+','match');
        nums = str2double(nums);
        
        Pin9(row)  = nums(1);
        Pin10(row) = nums(2);
        Pin11(row) = nums(3);
        Pin6(row)  = nums(end);
    end
end

% Trim in case preallocation was slightly large
TC0 = TC0(1:row);
TC1 = TC1(1:row);
TC2 = TC2(1:row);
TC3 = TC3(1:row);
TMP1 = TMP1(1:row);
TMP2 = TMP2(1:row);
Pin9 = Pin9(1:row);
Pin10 = Pin10(1:row);
Pin11 = Pin11(1:row);
Pin6 = Pin6(1:row);

% Create clean table
datafile = table(TC0,TC1,TC2,TC3,TMP1,TMP2,Pin9,Pin10,Pin11,Pin6);
end